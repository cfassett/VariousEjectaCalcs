import os, optparse, subprocess, sys, tempfile, ogr, osr, gdal
import pdb
from math import *
import numpy as np
import matplotlib.pyplot as plt
import copy
import random
from pynverse import inversefunc

def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
This program generates craters over the entire surface of the Moon for ejecta hazard assessment.
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Crater:
    def __init__(self, x=0.0, y=0.0, diam=0.0, dist=0.0, ejm=0.0, ejke=0.0):
        # crater properties for each crater
        self.x = x
        self.y = y
        self.diam=diam
        self.dist=dist
        self.ejm=ejm
        self.ejke=ejke


    def distance(self):
        """ Compute distance with Haversnie formula
        TO DO: This would be improved using a real geodesic answer; see pyproj"""
        
        lon1, lat1, lon2, lat2 = map(radians, [0.0, 0.0, self.x, self.y])
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2.0) ** 2.0 + cos(lat1) * cos(lat2) * sin(dlon / 2.0) ** 2.0
        hdist= 2.0 * 1747400.0 * asin(sqrt(a))
        
        return hdist
        
    def locate(self):
        """ Finds a random lon lat on the lunar sphere, correctly handling spherical geometry"""
        esriwktgeographic=['GEOGCS["Moon 2000",DATUM["D_Moon_2000",SPHEROID["Moon_2000_IAU_IAG",1737400.0,0.0]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]']
        esriwktcea=['PROJCS["Moon_CEA",GEOGCS["GCS_Moon_2000",DATUM["D_Moon_2000",SPHEROID["Moon_2000_IAU_IAG",1737400.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",0.0],UNIT["Meter",1.0]]']
        spatialRefgeographic = osr.SpatialReference()
        spatialRefgeographic.ImportFromESRI(esriwktgeographic)
        spatialRefcea = osr.SpatialReference()
        spatialRefcea.ImportFromESRI(esriwktcea)
        newx=random.uniform(-5458203.07634690677,5458203.07634690677)
        newy=random.uniform(-1737400.0,1737400.0)
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(newx, newy)
        coordTransform = osr.CoordinateTransformation(spatialRefcea, spatialRefgeographic)
        point.Transform(coordTransform)
        endx=point.GetX()
        endy=point.GetY()
        
        return endx, endy
    
        
    def ivanovimpactorsize(self):
        # This is Schmidt/Housen scaling from Ivanov 2001
        # It is what I commonly use, because of simplicity, and is the basis for some of the NPF calcs.  
        # However it may make logical sense to use Holsapple/Housen instead.  Answers are close in planetary science terms, at least with the right H2.
        
        gravity=1.62        # m/s2
        alpha=45.0          # impact angle degrees
        targdensity=2700.0  #kg/m3 (rho)
        impdensity=2700.0   #kg/m3 (delta)
        velocity=16000.0    # m/s    COULD choose a distribution?
        
        Dt=0.85*self.diam
        Dp=(Dt*((1.16*((impdensity/targdensity)**(1.0/3.0))*((velocity*sin(radians(alpha)))**0.43)*(gravity**-0.22))**-1.0))**(1.0/0.78)
        
        return Dp
    
    def himpactorsize(self):
        # This is Holsapple 1993 style
        # This particular scaling I think might have a typo in it since I don't think k1 should be so small.
        # not using for anything for now.
        
        gravity=1.62        # m/s2
        #strength=0.327e6      # Pa
        strength=0.3e6      # Pa
        targdensity=2700.0  #kg/m3 (rho)
        impdensity=2700.0   #kg/m3 (delta)
        velocity=16000.0    # m/s    COULD choose a distribution?   
        alpha=45            # impact angle degrees
        effvelocity=velocity*sin(radians(alpha))
        nu=0.4              # ~1/3 to 0.4
        mu=0.45             # ~0.4 to 0.55
        k1=0.004            # this is currently exremely dubious.  Probably this plus c on line 134 have two errors that cancel out.
        k2=0.8
        Kr=1.1              # R=Kr*V**0.333 from HH
        
        # a is the radius of the impactor
        
        Dt=0.85*self.diam
        #a=0.01    .... if gravity is important and we are solving for a, we need to do something iterative.  Assume strength
        #pi2=(gravity*a)/(effvelocity**2.0)     
        pi3=strength/(targdensity*(effvelocity**2.0))       
        densityratio=(targdensity/impdensity)
        
        exp1=(6.0*nu-2.0-mu)/(3.0*mu)
        exp2=(6.0*nu-2.0)/(3.0*mu)
        exp3=((2.0+mu)/2.0)
        exp4=(-3.0*mu)/(2.0+mu)
        
        #pi_v=k1*(pi2*(densityratio**exp1)+((k2*pi3*(densityratio**(exp2)))**exp3))**exp4
        pi_v=k1*(((k2*pi3*(densityratio**(exp2)))**exp3))**exp4
        #print pi2, pi3, pi_v
        
        
        ###### TO DO: Verify...
        
        # pi_v is the Crater Mass / impactor mass ratio
        # Crater Mass ratio calculations = ((Kr**3.0 * R**3.0 * targdensity)/((4.0/3.0)*pi*  a**3.0 * projdensity))*0.33333
        
        #Kr, sphere constants
        c=1.465
        cratersizeprojsizeratio=c*(densityratio**(-1.0/3.0))*pi_v
        #not totally sure these constants and calculation end up right... 
       
        
        Dp=Dt/cratersizeprojsizeratio  
        
        
        return Dp
        
    def hohoimpactorsize(self):
        # This is Housen Holsapple 2011 style
        
        gravity=1.62        # m/s2
        strength=5e3      # Pa
        targdensity=2700.0  #kg/m3 (rho)
        impdensity=2700.0   #kg/m3 (delta)
        velocity=16000.0    # m/s    COULD choose a distribution?   
        alpha=45            # impact angle degrees
        effvelocity=velocity*sin(radians(alpha))
        nu=0.4              # ~1/3 to 0.4
        mu=0.41             # ~0.4 to 0.55
        H2=1.0        
        Kr=1.1              # R=Kr*V**0.333 from HH
        
        # a is the radius of the impactor
        
        Dt=0.85*self.diam
        #a=0.01    .... if gravity is important and we are solving for a, we need to do something.  Assuming strength regime only matters for now.  This is a very good assumption.
        # TO DO: Add a check to verify strength regime
        #pi2=(gravity*a)/(effvelocity**2.0)     
        pi3=strength/(targdensity*(effvelocity**2.0))       
        densityratio=(targdensity/impdensity)
        
        exp1=(1.0-3.0*nu)
        exp2=(-3.0*mu)/(2.0)
        
        
        
        rhs=(H2**3.0)*(densityratio**(exp1))*(pi3**exp2)        
        acubed=densityratio*((Dt/2.0)**3.0)/(((4.0*pi)/3.0)*rhs)
        a=acubed**(1.0/3.0)
                    
        
        Dp=a*2.0
        
        
        return Dp
    
    def TotalEjectaMassGtEqVel(self, velocity=1.0):
         # This is Housen Holsapple 2011: POWER LAW PART ONLY -- quick appx
        gravity=1.62        # m/s2
        strength=5e3        # Pa
        targdensity=2700.0  #kg/m3 (rho)
        impdensity=2700.0   #kg/m3 (delta)
        nu=0.4              # ~1/3 to 0.4
        mu=0.41             # ~0.4 to 0.55
        n1=1.2
        H2=1.0   
        k=0.3               # range in table is 0.2-0.5
        C1=0.55             # range in table is 0.18-1.5
        C4=(3.0/4.0)*(1/pi)*(k)*C1**(3.0*mu)
        C6=C4*(H2**-3.0)
        
        # Assuming strength regime, add a check to verify strength regime.
        # TO DO: Add a check to verify strength regime
        
        # TO DO: error check, no zero velocity mofv gets has vel~^-3*mu
        
        Dt=0.85*self.diam
        R=0.5*Dt
        mofv=(targdensity*(R**3.0))*C6*(velocity*((targdensity/strength)**2.0))**(-3.0**mu)
        
        
        # inside n1 a, no ejecta
        a=self.hohoimpactorsize()/2.0
        secondterm=n1*(a**3.0)*targdensity*k
        
        maxm=k*targdensity*(R**3.0)-secondterm
        if mofv >= maxm:
            mofv = maxm
            
        #print maxm
        return mofv

    def TotalEjectaMassGtEqVelTable(self):
         # This is Housen Holsapple 2011
        gravity=1.62        # m/s2
        strength=5e3        # Pa
        targdensity=2700.0  #kg/m3 (rho)
        impdensity=2700.0   #kg/m3 (delta)
        impvelocity=16000.0 # m/s    COULD choose a distribution?   CAREFUL, not the ejecta velocity
        nu=0.4              # ~1/3 to 0.4
        mu=0.41             # ~0.4 to 0.55
        n1=1.2
        n2=1.0
        H2=1.0   
        p=0.3
        k=0.3               # range in table is 0.2-0.5
        C1=0.55             # range in table is 0.18-1.5
        C4=(3.0/4.0)*(1/pi)*(k)*C1**(3.0*mu)
        C6=C4*(H2**-3.0)
        densityratio=(targdensity/impdensity)

        
        # Assuming strength regime, add a check to verify strength regime.
        # TO DO: Add a check to verify strength regime
        
        Dt=0.85*self.diam
        R=0.5*Dt
        
        a=self.hohoimpactorsize()/2.0
        minlaunchlocovera=n1
        maxlaunchlocovera=(n2*R/a)
        intervals=200                                                        # Sets the velocity resolution of the table.  Higher is better, but slower
        spacing=(maxlaunchlocovera-minlaunchlocovera)/intervals
        xp=[]
        fp=[]
        for f in np.arange(minlaunchlocovera,maxlaunchlocovera,spacing):
            #launch positions in sequence
            #velocity at this launch position:
            vhere=impvelocity*C1*((f*(densityratio**nu))**(-1.0/mu))*((1.0-((f*a)/(n2*R)))**p)      #eq 14 from Housen and Holsapple 2011
            mhere=targdensity*(a**3.0)*k*((f**3.0)-(n1**3.0))
            #print vhere, mhere
            xp.insert(0,vhere)
            fp.insert(0,mhere)
        #print "0 "+str(targdensity*(a**3.0)*k*((maxlaunchlocovera**3.0)-(n1**3.0)))   # last bit of total excavated mass
        xp.insert(0,0)
        fp.insert(0,targdensity*(a**3.0)*k*((maxlaunchlocovera**3.0)-(n1**3.0)))   
        return xp,fp
        
    def TotalEjectaMassGtEqVelforV(self,velocity=1.0):                  
        if velocity<=0: velocity=0
        xp,fp=self.TotalEjectaMassGtEqVelTable()
        
        return np.interp(velocity, xp, fp)
        
    def ejectamasspersqm(self,distance,zonewidth=10):
        #distance is range to rim.  For nearfield, it would be important that range comes from the appropriate spot on interior of source crater.  can capture that with extra.
        craterR=self.diam/2.0
        extra=craterR*0.965    # note that extra should be at least n1*0.5*self.hohoimpactorsize() and might be as large as R for the nearfield.
        #to do this right for crater ejecta thickness mapping, we'd have to carefully correct for the changing extra with distance.  
        #currently extra is empirical to roughly match old rim ejecta thickness rules of thumb, which is fine.
        
        #in far field, distance >> extra so this correction is irrelevant.

        plusminus=zonewidth/2.0
        inner=distance-plusminus-craterR
        outer=distance+plusminus-craterR
        ejm=self.TotalEjectaMassGtEqVelforV(eqvel(distance-plusminus-extra))-self.TotalEjectaMassGtEqVelforV(eqvel(distance+plusminus-extra))     #this is TOTAL in a zonewidth-wide ring at that distance from the crater
        areatoouter=pi*((distance+plusminus))**2.0
        areatoinner=pi*((distance-plusminus))**2.0
        ejmpersquarem=ejm/(areatoouter-areatoinner)
        return ejmpersquarem

        
    
        
def npfcalc(size):
    '''this finds the NPF ==> N(>=size) for crater of "size"
    number for whole Moon (hard wired in, easy to change, see surfacearea)
    sizes must be in km''' 
    

    #lunar surface area
    surfacearea=37932328.1
    correctionfrom1millionyrsto1yr=1e-6
       
    log10freq=0.0
   
    npfcoeffs=[-6.076755981,-3.557528,0.781027,1.021521,-0.156012,-0.444058,0.019977,0.08685,-0.005874,-0.006809,0.000825,0.0000554]
    for a in range(12):
        log10freq=log10freq+npfcoeffs[a]*((log10(size))**a)
    
    globalfreq=correctionfrom1millionyrsto1yr*surfacearea*(10**log10freq)
        
    return globalfreq
    
def neukum(timestep=1.0):
    '''This takes the NPF for the Moon and generates craters'''
    '''units of timestep are years'''
    '''craters returned are in m, calculations in km and log10 km'''     

    
    # this threshold could be changed upward at some cost to performance.  A 1 km crater forms somewhere on the Moon every ~30000 yrs.
    # This is in km but needs to be a reasonable value in log10 land
    uppersizethreshold=1.0
    
    # NPF is not defined smaller than 10m, so if max hazard is at 10 m size, then we may want to extrapolate.
    # This is in km but needs to be a reasonable value in log10 land
    lowersizethreshold=0.01
    loglower=log10(lowersizethreshold)
   
   
    diams=[]  #start with an empty diameter list
    logsizestep=0.1 #start with this for large craters, since the odds that they form is laughably small   
    upperbinsize=log10(uppersizethreshold) #start at largest craters
    
    while (upperbinsize>loglower):   
        lowerbinsize=upperbinsize-logsizestep
        
        expectednumberinthisbin=timestep*(npfcalc(10**lowerbinsize)-npfcalc(10**upperbinsize))
        
        foundnumber=np.random.poisson(expectednumberinthisbin)
        
        
        if foundnumber>0.0:
            choiceset=((10**lowerbinsize)+(10**upperbinsize-(10**lowerbinsize))*np.random.random_sample(1000))
            wtssum=np.sum((choiceset)**-3.1)   #inbin distribution
            wts=((choiceset)**-3.1)/wtssum
            diams.extend((1000.0*np.random.choice(choiceset,foundnumber,p=wts)).tolist())            
               
        if lowerbinsize<-1.1: logsizestep=0.01 
        upperbinsize=copy.copy(lowerbinsize)
        
        
    return diams
        
        
def eqvel(dist=0):
    """ Finds the equivalent velocity needed to reach the distance of the crater away from the arbitrary point at the origin
       
        opposite of:
    
        def calcrange(velocity):
        gravity=1.62    # m/s2
        alpha=45.0      #degrees, ejecta launch angle
        Rp=1737400.0    #radius

        thetatop=0.5*(velocity**2)*sin(radians(alpha*2))
        thetabottom=(gravity*Rp)-(velocity*cos(radians(alpha)))**2.0
        theta=thetatop/thetabottom
        range=2.0*Rp*atan(theta)

        return range"""
    g=1.62    # m/s2
    alpha=45.0      #degrees, ejecta launch angle
    Rp=1737400.0    #radius
    
    z=tan(0.5*dist/Rp)
    top=z*g*Rp
    bottom=(0.5*sin(radians(alpha*2)))+z*(cos(radians(alpha))**2.0)
    vsq=top/bottom
    
    vel=vsq**0.5

    return vel 

def computeCDFejm(listofcrats):
    x=np.empty([len(listofcrats)])
    y=np.empty([len(listofcrats)])
    i=0
    sumofejm=0.0
    for crat in listofcrats:
        sumofejm=sumofejm+crat.ejm
        x[i]=crat.diam
        y[i]=sumofejm   
        i=i+1        
    return x,y
    
def computeCDFejke(listofcrats):
    x=np.empty([len(listofcrats)])
    y=np.empty([len(listofcrats)])
    i=0
    sumofejke=0.0
    for crat in listofcrats:
        sumofejke=sumofejke+crat.ejke
        x[i]=crat.diam
        y[i]=sumofejke
        i=i+1
    return x,y


def ethicknessmap(Crater):
    R=np.linspace(Crater.diam/2.0+0.1,5.0*Crater.diam,450)
    for curr in R:
        #print curr, Crater.ejectamasspersqm(curr,zonewidth=0.2)/2700.0  #2700 to convert to ht in m
        #convert curr to range from rim by subtracting a radius
        print str(curr-(Crater.diam/2.0))+","+str(Crater.ejectamasspersqm(curr,zonewidth=0.2)/2700.0) #2700 to convert to ht in m
    
    
    


def main():
    try:
        try:
            usage = "usage: ejecta.py\n"
            parser = optparse.OptionParser(usage=usage)
            (options, inargs) = parser.parse_args()
            #if not inargs: parser.error("Example text")    Currently setup for no arguments so this is pointless
            #firstarg=inargs[0]
            
        except optparse.OptionError, msg:
            raise Usage(msg)            
        

            # diams=neukum(timestep=1)
            # diams.sort(reverse=True)
            
            # maxm=0.0
            # cratlist=[]
            # for crats in diams:
                # newcrat=Crater(diam=crats)
                # newcrat.x,newcrat.y=newcrat.locate()
                # newcrat.dist=newcrat.distance()
                # newcrat.ejm=newcrat.ejectamasspersqm(newcrat.dist)
                # newcrat.ejke=0.5*newcrat.ejectamasspersqm(newcrat.dist)*(eqvel(newcrat.dist)**2.0)
                # if newcrat.ejm>maxm:
                    # maxm=newcrat.ejm
                    # maxd=newcrat.diam
                # cratlist.append(newcrat)
               
            # #x,y=computeCDFejm(cratlist)
            # #x,y=computeCDFejke(cratlist)
            
            # print maxd,maxm
            # plt.scatter(maxd,maxm,2)      
            # if a==0:
                # plt.xscale('log')
                # plt.yscale('log')
                # plt.xlim(10.0,1000.0)
                # plt.ylim(1.0e-13,1.0e-7)     

        # plt.show()
        
            
        newcrat=Crater(diam=100.0)
        print "range-from-rim, thickness"
        ethicknessmap(newcrat)
        #print newcrat.TotalEjectaMassGtEqVelTable()
      
    except Usage, err:
        print >>sys.stderr, err.msg
        # print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())