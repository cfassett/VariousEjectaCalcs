import os, optparse, sys, ogr, osr
from math import *
import numpy as np
import copy
import random
from scalingparams import *

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
        """ Compute great circle distances on sphere with Haversine formula, 
        This could be improved (perhaps unnecessarily) using a more complicated geodesic calculation"""
        
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
    
    def locate_inrange(self):
        """ Finds a random pt inside a particular distance on the lunar sphere, correctly handling spherical geometry"""
        
        distance=range_cutoff(self.diam)
        esriwktgeographic=['GEOGCS["Moon 2000",DATUM["D_Moon_2000",SPHEROID["Moon_2000_IAU_IAG",1737400.0,0.0]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]']
        esriwktaed=['PROJCS["Moon_AED",GEOGCS["GCS_Moon_2000",DATUM["D_Moon_2000",SPHEROID["Moon_2000_IAU_IAG",1737400.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Azimuthal_Equidistant"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]']
        esriwktlaea=['PROJCS["Moon_LAEA",GEOGCS["GCS_Moon_2000",DATUM["D_Moon_2000",SPHEROID["Moon_2000_IAU_IAG",1737400.0,0.0]],PRIMEM["Reference_Meridian",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]']
        spatialRefgeographic = osr.SpatialReference()
        spatialRefgeographic.ImportFromESRI(esriwktgeographic)
        spatialRefaed = osr.SpatialReference()
        spatialRefaed.ImportFromESRI(esriwktaed)
        spatialReflaea = osr.SpatialReference()
        spatialReflaea.ImportFromESRI(esriwktlaea)
        
        if distance>(1737400.0*pi):  #If there is no range limit, use the old code
            endx, endy = self.locate()
        else:
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(0.0, 0.0)
            poly=point.Buffer(distance)
            coordTransform = osr.CoordinateTransformation(spatialRefaed, spatialReflaea)
            poly.Transform(coordTransform)
            envelopeofpoly=poly.GetEnvelope()
            pointnotyet=True
            while pointnotyet:
                newx=random.uniform(envelopeofpoly[0],envelopeofpoly[1])  
                newy=random.uniform(envelopeofpoly[2],envelopeofpoly[3])
                testpoint = ogr.Geometry(ogr.wkbPoint)
                testpoint.AddPoint(newx,newy)
                if poly.Contains(testpoint):
                    pointnotyet=False
                    coordTransformgeo = osr.CoordinateTransformation(spatialReflaea, spatialRefgeographic)
                    testpoint.Transform(coordTransformgeo)
                    endx=testpoint.GetX()
                    endy=testpoint.GetY()
        return endx, endy
    
        
    def ivanovimpactorsize(self):
        # This is Schmidt/Housen scaling from Ivanov 2001
        # It is what I commonly use because of simplicity, and is the basis for some of the NPF calcs.  
        # However it may make logical sense to use Holsapple/Housen instead.  Answers are close in planetary science terms, at least with the right H2.
        
        #  ----- this is not currently used.  ---- 

        Dt=0.85*self.diam
        Dp=(Dt*((1.16*((impdensity/targdensity)**(1.0/3.0))*((velocity*sin(radians(alpha)))**0.43)*(gravity**-0.22))**-1.0))**(1.0/0.78)
        
        return Dp
        
    def overridestrength(self):
        # this takes a crater event size, and then turns it back to a strength estimate
        # it is roughly based on Fig. 9.26 of the lunar source book, although the uppper limit is a little arbitrary 
        # and the regolith thickness is kludgey.  It's probably more realistic than having a constant (low) strength, though.
        # 
        
        if ((self.diam/10.0)<10.0):
            stmodel=2500.0+(self.diam/10.0)*targdensity*tan(radians(40.0))  #Pa.  This is the regolith branch
        else:
            stmodel=5.0e6+(self.diam/10.0)*targdensity*tan(radians(60.0))  #Pa.  This is the fractured bedrock branch
    
        return stmodel
        
    def hohoimpactorsize(self):
        # This is Housen Holsapple 2011 style  (CRATER SIZE -> IMPACTORSIZE)
        # RELIES ON COMMON SCALING PARAMETERS FROM SCALINGPARAMS.PY   
        #
        # Solving for a, the radius of the impactor
        
        Dt=0.85*self.diam
        strength=self.overridestrength()
        pi3=strength/(targdensity*(effvelocity**2.0))       
                
        exp1=(1.0-3.0*nu)
        exp2=(-3.0*mu)/(2.0)     
                
        rhs=(H2**3.0)*(densityratio**(exp1))*(pi3**exp2)        
        acubed=densityratio*((Dt/2.0)**3.0)/(((4.0*pi)/3.0)*rhs)
        a=acubed**(1.0/3.0)
        Dp=a*2.0
        pi2=(gravity*a)/(effvelocity**2.0)     
        
        exp3=(2.0+mu)/2.0
        strengthgravitytransition=((H1/H2)**(exp3))*(densityratio**nu)*(pi3**(exp3))
        
        
        
        gravitycheck=False 
        if pi2>=strengthgravitytransition: gravitycheck=True
        
        if gravitycheck:
            z=(((4.0*pi)/3.0)**(1.0/3.0))*H1*(densityratio**((-2.0*nu)/(2+mu)))
            b=((Dt/(2.0*z))**(exp3))*((gravity/(effvelocity**2.0))**(mu/2.0))
            Dp=2.0*b
            pi2=(gravity*b)/(effvelocity**2.0)  
            exp4=((-2.0*nu)/(2.0+mu)) 
            exp5=(-mu/(2.0+mu))
            ratio=(((4.0/3.0)*pi)*(1.0/3.0))*H1*(densityratio**(exp4))*(pi2**exp5)            
            Dp=0.85*self.diam/ratio            
        #old     Dp=(((4.0*pi)/3.0)**(-(2.0+mu)/6.0))*((Dt/2.0)**(exp3))*(H1**(-exp3))*(densityratio**nu)*((gravity/(effvelocity**2.0))**(0.5*mu))*2.0
        else:
            Dp=2.0*a
        
        #Useful commands for checking the scaling
        #print self.diam, Dp, gravitycheck
        
        
        return Dp
    
    def TotalEjectaMassGtEqVelTable(self):
         # This is Housen Holsapple 2011
         # RELIES ON COMMON SCALING PARAMETERS FROM SCALINGPARAMS.PY  
         
         # This is the heart of the matter.
        
        
        strength=self.overridestrength()
        a=self.hohoimpactorsize()/2.0   #impactorradios
        pi2=(gravity*a)/(effvelocity**2.0)     
        pi3=strength/(targdensity*(effvelocity**2.0))       
        exp3=(2.0+mu)/2.0
        strengthgravitytransition=((H1/H2)**(exp3))*(densityratio**nu)*(pi3**(exp3))
        gravitycheck=False    
        if pi2>strengthgravitytransition: 
            n2=n2_g
        else:
            n2=n2_s
        
       
        Dt=0.85*self.diam
        R=0.5*Dt
        
        minlaunchlocovera=n1
        maxlaunchlocovera=(n2*R/a)
        intervals=10000                                                        # Sets the velocity resolution of the table.  Higher is better, but slower
        spacing=(maxlaunchlocovera-minlaunchlocovera)/intervals
        xp=[]
        fp=[]
        for f in np.arange(minlaunchlocovera,maxlaunchlocovera,spacing):
            #launch positions in sequence
            #velocity at this launch position:
            np.seterr(all='ignore')
            vhere=effvelocity*C1*((f*(densityratio**nu))**(-1.0/mu))*((1.0-((f*a)/(n2*R)))**p)      #eq 14 from Housen and Holsapple 2011
            mhere=targdensity*(a**3.0)*k*((f**3.0)-(n1**3.0))
            if np.isnan(vhere):
                pass
            else:
                xp.insert(0,vhere)
                fp.insert(0,mhere)
            np.seterr(all='warn')
        #print "0 "+str(targdensity*(a**3.0)*k*((maxlaunchlocovera**3.0)-(n1**3.0)))   # last bit of total excavated mass
        xp.insert(0,0)
        fp.insert(0,targdensity*(a**3.0)*k*((maxlaunchlocovera**3.0)-(n1**3.0)))   
        return xp,fp
        
    def TotalEjectaMassGtEqVelforV(self,velocity=1.0):                  
        assert velocity>=0
        xp,fp=self.TotalEjectaMassGtEqVelTable()
        
        return np.interp(velocity, xp, fp)
        
    def ejectamasspersqm(self,zonewidth=10):
        plusminus=zonewidth/2.0
        extra=0.0    # note that extra should be at least n1*0.5*self.hohoimpactorsize() and might be as large as R for the nearfield.
        ejm=self.TotalEjectaMassGtEqVelforV(eqvel(self.dist-plusminus-extra))-self.TotalEjectaMassGtEqVelforV(eqvel(self.dist+plusminus-extra))     #this is TOTAL in a zonewidth-wide ring at that distance from the crater
        areatoouter=pi*((self.dist+plusminus))**2.0
        areatoinner=pi*((self.dist-plusminus))**2.0
        ejmpersquarem=ejm/(areatoouter-areatoinner)
        return ejmpersquarem

def hohocratersize(impactorsize):
    # This is Housen Holsapple 2011 style  (IMPACTOR DIAMETER->CRATER SIZE in DIAM)
    # RELIES ON COMMON SCALING PARAMETERS FROM SCALINGPARAMS.PY  
    
    empcorr=1.024743105 # for some reason, there is an issue getting gravity scaling to match .  Minor, 2%.  This is a kludge

    
    impradius=impactorsize/2.0
    
    strength=self.overridestrength()
    pi2=(gravity*impradius)/(effvelocity**2.0)     
    pi3=strength/(targdensity*(effvelocity**2.0))       
    exp3=(2.0+mu)/2.0
    strengthgravitytransition=((H1/H2)**(exp3))*(densityratio**nu)*(pi3**(exp3))
    gravitycheck=False    
    if pi2>strengthgravitytransition: gravitycheck=True
    if gravitycheck:
        exp4=((-2.0*nu)/(2.0+mu))           
        exp5=(-mu/(2.0+mu))
        ratio=(((4.0/3.0)*pi)*(1.0/3.0))*H1*(densityratio**(exp4))*(pi2**exp5)
        cratersize=ratio*2.0*impradius*(1.0/0.85)*empcorr
    else:        
        exp1=(1.0-3.0*nu)
        exp2=(-3.0*mu)/(2.0)               
        rhs=(H2**3.0)*(densityratio**(exp1))*(pi3**exp2)        
        cratersize=2.0*impradius*(1.0/0.85)*(((4.0/3.0)*pi*(densityratio**-1.0)*rhs)**(1.0/3.0))
    
    return cratersize
    
    
def areaforadistance(distance):
    """ if distance is greater than piR, then just returns lunar surface area
    this version uses the spherical cap formula and some trig.  There is another
    version using map-projections that works in afdtest2 which might be useful for random point generation"""
    
    lunarradius=1737400.0
    lunarsurfacearea=4*pi*lunarradius**2.0
    
    if distance>(lunarradius*pi):
        area=lunarsurfacearea
    else:
        theta=distance/lunarradius
        h=lunarradius*(1-cos(theta))
        area=2*pi*lunarradius*h
    return area

def gruncalc(size):
    '''this finds the Grun ==> N(>=size) for crater with -crater size- size
    number for whole Moon (hard wired in, easy to change, see surfacearea)
    input sizes must be in m (important to track carefully, I know)
    
    >=10m just calls the NPF    
    
    There was originally an error in this calculation in (in the grundfluxeslog) that has been corrected.
    ''' 
    
    #theseareparticleradii in cm
    #grunsizelist=[4.57078E-07, 9.84745E-07, 2.12157E-06, 4.57078E-06, 9.84745E-06, 2.12157E-05, 4.57078E-05, 9.84745E-05, 0.000212157,\
    #0.000457078, 0.000984745, 0.002121569, 0.004570781, 0.00984745, 0.021215688, 0.045707815, 0.098474502, 0.212156884, 0.45707815, 0.984745022,\
    #2.121568836, 4.570781497, 9.847450218, 21.21568836]      
    #for i in grunsizelist:
    #    print hohocratersize(i/100.0)   #/100 is because of meters
    
    
    #these are in m, assuming hohoscaling -- calculated above.  Can switch back to radii in cm with newcrat.hohoimpactorsize()/0.020   if desired.
    grundiameters=[6.91E-07,1.49E-06,3.21E-06,6.91E-06,1.49E-05,3.21E-05,6.91E-05,0.000148787,0.000320552,0.000690607,0.001487869,0.003205516,0.006906074,0.014878687,0.03205516,\
    0.069060751,0.148786877,0.32055161,0.690607508,1.487868771,3.205516094,6.906075072,14.87868771,32.05516094]

  
    grundiameterslog=[-6.160768847, -5.827435381, -5.4941018, -5.160768847, -4.827435381, -4.4941018, -4.160768847, -3.827435381, -3.4941018, -3.160768847, -2.827435381, -2.494102005, -2.160768752, -1.827435381, -1.494102046, -1.160768705, -0.827435372, -0.494102037, -0.160768705, 0.172564628, 0.505897962, 0.839231295, 1.172564628, 1.505897962]
    
    grundfluxeslog=[14.51350132, 14.33044304, 13.49557686, 12.53602555, 11.57428174, 10.61423746, 9.668560138, 8.818172003, 8.266651358, 7.940392878, 7.553930208, 6.975730092, 6.171552705, 5.16532, 4.014713602, 2.775673296,1.487500212, 0.173514165, -1.153632971, -2.487387613, -3.824433773, -4.835279357, -5.940706041, -6.96962463]

    surfacearea=37932328.1
          
    if size<=10.0:
        slog=log10(size)
        fluxintlog=np.interp(slog, grundiameterslog, grundfluxeslog)
        fluxint=10.0**fluxintlog
        globalfreq=surfacearea*fluxint
    else:
        globalfreq=npfcalc(size/1000.0)
        
    return globalfreq
    
    

def npfcalc(size):
    '''this finds the NPF ==> N(>=size) for crater with -crater size- size
    number for whole Moon (hard wired in, easy to change, see surfacearea)
    sizes must be in km''' 
    

    #lunar surface area
    surfacearea=37932328.1
    correctionfrom1millionyrsto1yr=1.0e-6
       
    log10freq=0.0
   
    npfcoeffs=[-6.076755981,-3.557528,0.781027,1.021521,-0.156012,-0.444058,0.019977,0.08685,-0.005874,-0.006809,0.000825,0.0000554]
    for a in range(12):
        log10freq=log10freq+npfcoeffs[a]*((log10(size))**a)
    
    globalfreq=correctionfrom1millionyrsto1yr*surfacearea*(10.0**log10freq)
        
    return globalfreq
    
def neukum(timestep=1.0, lowersizethreshold=0.01, uppersizethreshold=1.0):
    '''This takes the NPF for the Moon and generates craters'''
    '''units of timestep are years'''
    '''craters returned are in m, calculations in km and log10 km'''     

    
    # this threshold could be changed upward at some cost to performance.  A 1 km crater forms somewhere on the Moon every ~30000 yrs.
    # This is in km but needs to be a reasonable value in log10 land
    #uppersizethreshold=1.0  DEFAULT
    
    # NPF is not defined smaller than 10m, so if max hazard is at 10 m size, then we may want to extrapolate.
    # This is in km but needs to be a reasonable value in log10 land
    # lowersizethreshold=0.01  DEFAULT
    
    
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
    
    
def grun(timestep=10.0*(1/86400.0)/365.25, lowersizethreshold=0.1, uppersizethreshold=1.0):
    '''This takes the grun extending into the NPF and  generates craters'''
    '''units of timestep are years'''
    '''craters returned are in m, calculations in m and log10 m'''     

    
    # this threshold could be changed upward at some cost to performance.  A 1 km crater forms somewhere on the Moon every ~30000 yrs.
    # This is in km but needs to be a reasonable value in log10 land
    #uppersizethreshold=1.0  DEFAULT
    
    # NPF is not defined smaller than 10m, so if max hazard is at 10 m size, then we may want to extrapolate.
    # This is in km but needs to be a reasonable value in log10 land
    # lowersizethreshold=0.01  DEFAULT
    
    
    loglower=log10(lowersizethreshold)     
    diams=[]  #start with an empty diameter list
    logsizestep=0.1 #start with this 
    upperbinsize=log10(uppersizethreshold) #start at largest craters
    
    eta=1.0e-6   # this is to avoid floating point nastiness that might result in heading into the next bin,arg 
    while (upperbinsize>(loglower+eta)):   
        lowerbinsize=upperbinsize-logsizestep        
        expectednumberinthisbin=timestep*(gruncalc(10**lowerbinsize)-gruncalc(10**upperbinsize))        
        #print (10**lowerbinsize), (10**upperbinsize), expectednumberinthisbin
        foundnumber=np.random.poisson(expectednumberinthisbin)
        
        
        if foundnumber>0.0:
            choiceset=((10**lowerbinsize)+(10**upperbinsize-(10**lowerbinsize))*np.random.random_sample(foundnumber*10))


            if lowerbinsize<=1:
                alpha=-2.3
            else:
                alpha=-3.1
            wtssum=np.sum((choiceset)**alpha)   #inbin distribution
            wts=((choiceset)**alpha)/wtssum
            diams.extend((np.random.choice(choiceset,foundnumber,p=wts)).tolist())              # could error check these to make sure that nothing < 10**lowerbinsize are made.  
               
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
 
def mass_to_numberofparticles_table(mass):
    # mass bins in kg
    # equivalent sizes are diameters (in mm) of 1.5, 1, 0.4, 0.3, 0.2, 0.1, 0.075, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005 0.001

    massbins=[4.41786E-06, 1.309E-06, 8.37758E-08, 3.53429E-08, 1.0472E-08, 1.309E-09, 5.52233E-10, 1.63625E-10, 8.37758E-11, 3.53429E-11, 1.0472E-11, 1.309E-12,1.63625E-13,1.309E-15]
    
    # This is from the Carrier standard lunar regolith PSD 
    particlesperkg=[11317.68484, 22918.31181, 1074295.866, 1414710.605, 5729577.951, 114591559, 126758070.2, 488923985.2, 596831036.6, 1697652726, 8594366927, 61115498147, 3.66693E+11, 5.34761E+13]

    particles_scaled=[i*mass for i in particlesperkg]   #particles for this mass

    return massbins, particles_scaled

def numbers_given_bintable(mass,limitingmass=1.309E-09):   #default limiting mass is eq. to >100micron particles
    massbins,particles=mass_to_numberofparticles_table(mass)
    assert limitingmass<massbins[0]
    assert limitingmass>massbins[-1]
    cumnum=0
    partialbin=False
    interplastmass=[]
    interplastnum=[]
    for i in range(len(massbins)):
        if limitingmass < massbins[i]:
            cumnum=cumnum+particles[i]
        else:
            if not partialbin:
                partialbin=True
                interplastmass.append(log10(massbins[i]))
                interplastmass.append(log10(massbins[i-1]))
                interplastnum.append(log10(particles[i]))
                interplastnum.append(log10(particles[i-1]))
                lastnumlog=np.interp(log10(limitingmass), interplastmass, interplastnum)
                cumnum=cumnum+10.0**lastnumlog
    
    return cumnum

def numbers_given_limiting_mass_formula(mass,limitingmass=1.309E-09):   #default limiting mass is eq. to >100micron particles
    # these power laws are approximations and might be a bit imperfect
    if limitingmass<1.309e-6:   #if limiting mass is smaller than a mm
        cumnumperkg=0.1604*(limitingmass**-0.981)
    else:                       #if limiting mass is larger than a mm, power law is slightly less steep
        cumnumperkg=0.1507*(limitingmass**-0.91)
    cumnum=mass*cumnumperkg   
        
    return cumnum

def range_cutoff(cratsize):   #Ejecta flux relevant craters are only nearby for small sizes.  But generating a fixed number artificially puts a bunch far away.
                              #using a range cutoff solves this problem.  Here, the range cutoff is chosen to be approximately where the kg/m2 falls below 1e-8.
                                  
                              #This formula cutoffs were calculated with ejecta_montecarlo_masspervelocitythresholdrangecalc.py and approximated in excel.  
    if cratsize>1.259:
        rangecutoff=24000.0*(cratsize**1.3996)
    else:
        rangecutoff=13000.0*(cratsize**1.1615)
                               
                                                           
        
    return rangecutoff


def main():
    try:
        try:
            usage = "usage: ejecta_montecarlo.py\n"
            parser = optparse.OptionParser(usage=usage)
            (options, inargs) = parser.parse_args()
            #if not inargs: parser.error("Example text")    Currently setup for no arguments so this is pointless.
            #firstarg=inargs[0]
            
        except optparse.OptionError, msg:
            raise Usage(msg)            
            
        lunarsurfacearea=4*pi*1737400.0**2.0
        mlist=[]
        
        
        # for <2.5 m, monte carlo a fixed number of events and then rescale.  See fixed number below for more.  for >2.5 m, just use poisson to generate impacts randomly correctly     
        timesample=1.0 # default is 1 == 1 year. Increasing this increases the sampling timesample.        
        logmin=-2.0     #May make sense to drop to -2.4  Lower mass delivery threshold makes this a bit smaller than canonical recent model.
        logmax=3.0      #maximum crater diameter =10^3 = 1 km.  (Very big, uncommon)
        logstep=0.1
       
        ranger=np.arange(logmin,logmax,logstep).tolist()
        
        minvelcutoff=5.0 # Very slow velocities are ignored. (<minvelcutoff)  With a 5 m/s cutoff, and the current rangecutoff scheme, ejecta from craters smaller than 3 cm are irrelevant.  This doesn't seem crazy.
        
        
                               
        trials=100       #Number of monte carlo trials.  A big number here is definitely better as well, but the median is not sensitive to this being super large.  Capturing the tail risk might be!  Bigger is definitely slower.
        
        
        
              
        
          
        
        listofxpfppairs=[]
        for loglower in ranger:
            newcrat=Crater(diam=(10**(0.5*logstep+loglower)))
            xp,fp=newcrat.TotalEjectaMassGtEqVelTable()   #generate the ejecta mass table at a given size only once for efficiency
            listofxpfppairs.append([xp,fp])
        
        
        for monte in range(trials):                       #do monte carlo trials.
            print monte
            totalmassperunitarea=0.0            
            ct=0  #increment so the appropriate massvelocity table is pulled
            for loglower in ranger: 
                massinsizebin=0.0
                logupper=loglower+logstep               
                rc=range_cutoff((10**(0.5*logstep+loglower)))  #how far to look fro craters
                expectednumberinthisbinonMoon=timesample*(gruncalc(10**loglower)-gruncalc(10**logupper))   # per timesample, wholeMoon
                expectednumberinrange=expectednumberinthisbinonMoon*(areaforadistance(rc)/lunarsurfacearea)
                #print loglower, expectednumberinrange, areaforadistance(rc)/lunarsurfacearea
                newcrat=Crater(diam=(10**(0.5*logstep+loglower)))           
                making=np.random.poisson(expectednumberinrange)
                for a in range(making):
                    newcrat=Crater(diam=(10**(0.5*logstep+loglower)))
                    newcrat.x,newcrat.y=newcrat.locate_inrange()  
                    newcrat.dist=newcrat.distance()
                    vel=eqvel(newcrat.dist)
                    if vel>minvelcutoff:
                        plusminus=newcrat.diam/100.0                        
                        vela=eqvel(newcrat.dist-plusminus)
                        velb=eqvel(newcrat.dist+plusminus)
                        xp=listofxpfppairs[ct][0]
                        fp=listofxpfppairs[ct][1]
                        ma=np.interp(vela, xp, fp)
                        mb=np.interp(velb, xp, fp)
                        areatoouter=pi*((newcrat.dist+plusminus))**2.0
                        if newcrat.dist>plusminus:
                            areatoinner=pi*((newcrat.dist-plusminus))**2.0
                            ejmpersquarem=(ma-mb)/(areatoouter-areatoinner)
                        else:
                            ma=np.interp(0.0, xp, fp)
                            ejmpersquarem=(ma-mb)/areatoouter
                        massinsizebin=massinsizebin+ejmpersquarem                         
                totalmassperunitarea=totalmassperunitarea+massinsizebin  
                ct=ct+1
            mlist.append(totalmassperunitarea/timesample)                                                            #This list is recording each trials total mass reaching an arbitrary sq m per year (normalizing out timesample)
            
        masarray=np.asarray(mlist)
        averagemass=np.average(masarray)                                                                            # what is the average flux?        
        lowermasskglog=-14.0                                                                                             
        uppermasskglog=-5.0 
        mranger=np.arange(lowermasskglog,uppermasskglog,1.0).tolist()
        for limitingmasslog in mranger:                                                                                 # This partitions the median flux into an # of particles of given particle mass or greater.
           print "avg",limitingmasslog, numbers_given_limiting_mass_formula(averagemass,10.0**limitingmasslog)           # -14 is 10^-14 kg = 10^-11 g 
                                                                                                                        #-8 is 10^-8 kg = 10^-5 g or roughly a 76 micron particle.

        
    except Usage, err:
        print >>sys.stderr, err.msg
        # print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())