####!/home/yuan/source/yt-conda/bin/python
import yt
yt.enable_parallelism()
import matplotlib
matplotlib.use('Agg')
import numpy as numpy
import matplotlib.pylab as pylab
import os as os
import math

homedir = os.environ['HOME']
name_pattern = 'DD%04i/stest_%04i'
GravConst = 6.67e-8
SolarMass = 1.989e33
kpc=3.086e21
keV=1.1604e7
conc=6.81
M_vir=8.5e14
Rs=0.0224   #SphereCoreRadius=0.0224
TimeUnits=1.78152e+17
LengthUnits=4.93708012753e+25
VelocityUnits=LengthUnits/TimeUnits
kboltz=1.38e-16
m_p=1.67e-24

pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
pylab.rc('lines', linewidth=2)
font = {'weight' : 'bold', 'size'   : 18}
pylab.rc('font', **font)
pylab.rc('text', usetex=False)
#pylab.rc('text', usetex=True)
pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
pylab.rcParams['mathtext.fontset'] = 'custom'
pylab.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
pylab.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
pylab.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

def profile_more(imin=0,imax=None,rmin=1, rmax=1000.0, mod=1,CoolingTimeOn=1,xbins=128):
 if imax is None:
    imax=imin+1
 #cm=pylab.cm.Blues_r
 #cm=pylab.cm.autumn
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
# pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 #cm=pylab.cm.gist_heat
 cm=pylab.cm.gist_rainbow
 conc=6.81
 M_vir=8.5e14
 Rs=0.0224   #SphereCoreRadius=0.0224
 loop=range(imin,imax,mod)
# loop.insert(-1,0)
 for i in loop:
    fn = name_pattern %(i,i)
    pf = yt.load(fn)
#    slc = yt.SlicePlot(pf, "z", "temperature")
#    slc.save()
    sp = pf.sphere([0.5, 0.5, 0.5], (rmax, 'kpc'))
    prof = yt.Profile1D(sp, "radius", xbins, rmin*kpc, rmax*kpc, True, weight_field="cell_mass")
    prof.add_fields(["temperature","density","pressure"])
    r=numpy.array(prof.x)/kpc    # from cm to kpc
    Temperature=numpy.array(prof["temperature"])
    Density=numpy.array(prof["density"])
    Pressure=numpy.array(prof["pressure"])*6.242e8   # from ergs to keV
    pylab.figure(0,figsize=(8, 8))
    if i==0:
       pylab.loglog(r,Temperature,color="k",linewidth=3.0)
    else:
       pylab.loglog(r,Temperature,color=cm(float(i-imin)/(imax-imin)))
    #pylab.loglog(np.array(prof.x), np.array(prof["temperature"]),color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('Radius $(g/cm^3)$')
    pylab.ylabel('Temperature $(K)$')
    pylab.xlim([rmin,rmax])
    pylab.ylim([4e6,2e8])
    pylab.xlabel('r (kpc)')
    pylab.ylabel('T (K)')
    pylab.tight_layout()
    pylab.savefig('Temperature_rainbow.png')

    pylab.figure(1,figsize=(8, 8))
    if i==0:
       pylab.loglog(r,Pressure,color="k",linewidth=3.0)
    else:
       pylab.loglog(r,Pressure,color=cm(float(i-imin)/(imax-imin)))
    pylab.xlim([rmin,rmax])
    pylab.ylim([4e-3,10])
    pylab.xlabel('r (kpc)')
    pylab.ylabel('P '+r'$(\rm{keV}\ /\ \rm{cm}^3)$')
    pylab.tight_layout()
    pylab.savefig('Pressure_rainbow.png')

    pylab.figure(2,figsize=(8, 8))
    if i==0:
       pylab.loglog(r,Density,color="k",linewidth=3.0)
    else:
       pylab.loglog(r,Density,color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('r (kpc)')
    pylab.ylabel('Density '+r'$(\rm{g}\ /\ \rm{cm}^3)$')
    pylab.xlim([rmin,rmax])
    pylab.ylim([1e-27,3e-23])
    pylab.tight_layout()
    pylab.savefig('Density_rainbow.png')

    if CoolingTimeOn==1:
      prof.add_fields(["cooling_time"])
      CoolingTime=numpy.array(prof["cooling_time"])/3.15e7 # from s to year
      x1=r/(Rs*kpc)
      M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass #DM+BCG + BH , cgs
      t_dyn=numpy.pi*(r*kpc)**1.5/(2.0*(GravConst*M_dyn)**0.5)/3.15e13   #from s to Myr
      fig = pylab.figure(figsize=(8, 8))
      pylab.figure(3)
      pylab.loglog(r,CoolingTime/1.0e6,color=cm(float(i-imin)/(imax-imin)))
      pylab.loglog(r,t_dyn,color="green",linewidth=2,linestyle='--')
      pylab.xlim([rmin,rmax])
      pylab.ylim([4,7e4])
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Cooling Time (Myr)')
      pylab.tight_layout()
      pylab.savefig('CoolingTime_rainbow.png')
 pylab.clf()
 return



def slice(imin=0,imax=None,widthkpc=None,mod=1,cc=None):
 if imax is None:
    imax=imin+1
 for i in range(imin,imax,mod):
    fn = name_pattern %(i,i)
    pf = yt.load(fn)
    tt=yt.SlicePlot(pf, 'y', ["temperature","density","cooling_time","SN_Colour","Color_Fraction","pressure","CoolingRate"],fontsize=40)
    tt._setup_plots()
    ax = tt.plots['temperature'].axes
    ax.set_xticks([-7.5, 7.5])
#    tt.annotate_grids()
    if widthkpc is not None:
       tt.set_width(widthkpc, 'kpc')
    if cc is not None:
       tt.set_center(cc,unit="code_length")
    tt.set_zlim("temperature",1e4,4e6)
    tt.set_zlim("SN_Colour",1e-4,100)
    tt.set_cmap("SN_Colour","coolwarm")
    tt.set_zlim("pressure",1e-12,1e-11)
    tt.set_cmap("Color_Fraction","coolwarm")
    tt.set_zlim("Color_Fraction",1e-6,2)
    tt.set_colorbar_label("SN_Colour", "Tracer Density")
    tt._setup_plots()
    ax = tt.plots["pressure"].axes
    ax.set_xticks([-1, 0, 1])
    ax = tt.plots["temperature"].axes
    ax.set_xticks([-1, 0, 1])
    ax = tt.plots["SN_Colour"].axes
    ax.set_xticks([-1, 0, 1])

    tt.save()

    nn=yt.SlicePlot(pf, 'y', "density",fontsize=40)
    nn.set_zlim("density", 4e-26, 1e-22)
    if widthkpc is not None:
       nn.set_width(widthkpc, 'kpc')
    if cc is not None:
       nn.set_center(cc,unit="code_length")
#    nn.annotate_velocity(factor = 16)
    nn._setup_plots()
    ax = nn.plots["density"].axes
    ax.set_xticks([-1, 0, 1])
    nn.save()
 return

def Mira_x_slice(imin=0,imax=1,nslice=20):
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   stest = 'stest_%04i' %i
   pf = yt.load(fn)
   for j in range(0,nslice):
    xx=4.0*j/nslice
    cc=[xx,0.5,0.5]
    #tt=yt.SlicePlot(pf, 'x', ["temperature","density","cooling_time","SN_Colour","Color_Fraction","pressure"],center=cc)
    tt=yt.SlicePlot(pf, 'x', "temperature",center=cc)
#    print(cc)
#    tt.set_center(cc,unit="code_length")
    tt.set_zlim("temperature",1e4,4e6)
#    tt.set_zlim("SN_Colour",1e-4,100)
#    tt.set_cmap("SN_Colour","coolwarm")
#    tt.set_zlim("pressure",1e-12,1e-11)
#    tt.set_cmap("Color_Fraction","coolwarm")
#    tt.set_zlim("Color_Fraction",0.01,2)
    tt.save(stest+"_x_slice_%i.png"%j)



from yt import derived_field
@derived_field(name="ZeusHeatingRate", units="erg/s/cm**3")
def _ZeusHeatingRate(field, data):
    ret = 2.0*data['dx'].in_cgs()**2*data['density'].in_cgs()*data['velocity_divergence'].in_cgs()**3
    ret[data['velocity_divergence'] > 0] = 0
    return -ret

@derived_field(name="ZeusHeatingCell",units = "erg/s")
def _ZeusHeatingCell(field,data):
#    return numpy.array(data["ZeusHeatingRate"])*numpy.array(data["index", "cell_volume"].in_cgs()) ##error
    return data["ZeusHeatingRate"].in_cgs()*data["index", "cell_volume"].in_cgs()

@derived_field(name="Cell_ThermalEnergy", units= "erg")
def _Cell_ThermalEnergy(field,data):
  return data["thermal_energy"]*data["gas","cell_mass"]

@derived_field(name="logP",units="")
def _logP(field,data):
    return numpy.log(numpy.array(data["gas","pressure"]))

@derived_field(name="CoolingDistance",units="cm")
def _CoolingDistance(field,data):
    CD=numpy.abs(data["x-velocity"].in_cgs())*data["cooling_time"].in_cgs()
    return CD

#thermal_energy gives strange results (off by two orders of mag) what is thermal_energy? (units indicate it is energy per mass)
#@derived_field(name="CoolingRateCell", units="erg/s")
#def _CoolingRateCell(field,data):    # this is the rate per cell
#   return data["gas","thermal_energy"]*data["gas","cell_mass"]/(data["cooling_time"])

#@derived_field(name="CoolingRate", units="erg/s/cm**3")
#def _CoolingRate(field,data):    # this is the rate per volume
#   return data["gas","thermal_energy"]*data["gas","cell_mass"]/(data["cooling_time"]*data["index","cell_volume"].in_cgs())

def ReadCoolingCurve():
    FileIn = open(homedir+"/code/cool_rates.in")
    data = numpy.loadtxt(FileIn)
    FileIn.close()
    return (data[:,0],data[:,2])

@derived_field(name="CoolingRate", units="erg/s/cm**3")
def _CoolingRate(field,data):    # this is the rate per volume
   (LogT, LogCoolRate) = ReadCoolingCurve()
##   Lambda = 10.0**numpy.interp(numpy.log10(numpy.array(data["Temperature"])), LogT, LogCoolRate)  #this is lambda, but I used to call it "CoolingRate" in allyt
   Lambda = 10.0**numpy.interp(numpy.log10(numpy.array(data["gas","temperature"])), LogT, LogCoolRate)  #this is lambda, but I used to call it "CoolingRate" in allyt
   CoolingRate=(numpy.array(data["gas","density"])/(0.59*m_p*2.0))**2*Lambda  #ne is half of n
   return data.ds.arr(CoolingRate,"erg/s/cm**3")

#@derived_field(name="CoolingRateCell", units="erg/s")
#def _CoolingRateCell(field,data):    # this is the rate per cell
#   return data["CoolingRate"]*data["index","cell_volume"].in_cgs()


def Heating(imin=0,imax=None,rmax=100,rmin=0.5,mod=1,n_bins=40,Tmin=1e4,Tmax=1e9,useSam=True):
 if imax is None:
   imax=imin+1
 r=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 CoolingRate=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 HeatingRate=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 ShockHeatingRate=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 Cooling=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 Heating=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 ShockHeating=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 CoolingRate_T=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 HeatingRate_T=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 ShockHeatingRate_T=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 density_T=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 volume=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 shock_volume=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 strong_shock_volume=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 T=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 ShellThermal=numpy.zeros((math.ceil(float(imax-imin)/mod),n_bins))
 time=[]
 Edot=[]
 TotalCooling=[]
 TotalCoolingCell=[]
 j=0
 for i in range(imin,imax,mod):
   fn = name_pattern %(i,i)
   ds = yt.load(fn)
   time.append(ds.current_time)
   Edot.append(ds["ClusterSMBHJetEdot"]*2.0*1.0e44)
   ds.add_gradient_fields(("gas","temperature"))
   ds.add_gradient_fields(("gas","entropy"))
   ds.add_gradient_fields(("gas","logP"))

   def _ShockMach(field,data):  #value is shock Mach number
     ep=0.002  # Mach of 1.001
#     ep=0.223  # Mach of 1.1
     new_field = numpy.zeros(data["pressure"].shape, dtype='float64')
     dog=numpy.maximum(numpy.array(data["logP_gradient_x"])*numpy.array(data["index", "dx"].in_cgs()),numpy.array(data["logP_gradient_y"])*numpy.array(data["index", "dy"].in_cgs()))
     dlogP=numpy.maximum(dog,numpy.array(data["logP_gradient_z"])*numpy.array(data["index", "dz"].in_cgs()))
     too_big=dlogP>2  # if the shock is too strong (mach > 5ish), mouse will cause overflow problem
     dlogP[too_big]=2
     if useSam==True:
       TSx=numpy.array(data['gas', 'temperature_gradient_x'])*numpy.array(data['gas', 'entropy_gradient_x'])
       TSy=numpy.array(data['gas', 'temperature_gradient_y'])*numpy.array(data['gas', 'entropy_gradient_y'])
       TSz=numpy.array(data['gas', 'temperature_gradient_z'])*numpy.array(data['gas', 'entropy_gradient_z'])
       TS1=numpy.logical_or((TSx>0), (TSy>0))
       TS=numpy.logical_or(TS1,  (TSz>0))
       good=(TS & (numpy.array(data["gas", "velocity_divergence"])<0)) & (dlogP > ep)
     else:
       cat=numpy.array(data["gas","thermal_energy"]*data["gas","cell_mass"])/numpy.array(data["gas","kinetic_energy"]*data["index","cell_volume"].in_cgs()) > (1.0/9.0)
       good=(cat & (numpy.array(data["gas", "velocity_divergence"])<0))& (dlogP > ep) 
#     good=dlogP > ep
     mouse = (numpy.exp(dlogP)-1.0)*0.4 + 1
     new_field[good]=mouse[good]
     return new_field
   ds.add_field(('gas','ShockMach'), function=_ShockMach, units="", take_log=False,display_name='Shock Mach')


   def _ZeusShockHeatingRate(field,data):
     ZeusShockHeatingRate=numpy.array(data["ZeusHeatingRate"])
     bad = numpy.array(data["gas","ShockMach"])<1
     ZeusShockHeatingRate[bad]=0.0
#   return data.ds.arr(ZeusShockHeatingRate,data["ZeusHeatingRate"].units)
     return data.ds.arr(ZeusShockHeatingRate,data["ZeusHeatingRate"].units)
   ds.add_field(('gas','ZeusShockHeatingRate'), function=_ZeusShockHeatingRate, units="erg/s/cm**3")

   def _ZeusShockHeatingCell(field,data):
#    return numpy.array(data["ZeusShockHeatingRate"])*numpy.array(data["index", "cell_volume"].in_cgs())
      return data["ZeusShockHeatingRate"]*data["index", "cell_volume"]
   ds.add_field(('gas','ZeusShockHeatingCell'), function=_ZeusShockHeatingCell, units="erg/s")

   sp_all = ds.sphere([0.5, 0.5, 0.5], (rmax,'kpc'))
#   sp = sp_all.cut_region(["obj['index','x']>0.5","obj['index','y']>0.5","obj['index','z']>0.5"])
#   sp=sp_all
#   slc = yt.ProjectionPlot(sp, "z", ["ZeusHeatingRate","ShockMach"], width=(rmax/2, 'kpc'))
  # sp_hot=sp.cut_region(["obj['temperature'] > 1e6","obj['ShockMach'] < 2"])
#   sp_hot=sp.cut_region(["obj['temperature'] > 1e6"])
   sp=sp_all.cut_region(["obj['temperature'] > 1e6","obj['density']<1.0e-23"])  #cut out SF region
   #slc = yt.SlicePlot(ds, "x", ["ZeusHeatingRate","ShockMach","pressure",("gas","entropy"),"mach_number"], width=(rmax*2, 'kpc'),data_source=sp)
   slc = yt.SlicePlot(ds, "x", ["ZeusHeatingRate","ShockMach","pressure",("gas","entropy"),"mach_number"], width=(rmax*2, 'kpc'))
   slc.set_zlim("pressure", 1e-12,1e-7)   
   slc.set_zlim(("gas","entropy"), 1e-13,1e-4)   
   slc.set_zlim("ShockMach", 1, 1.2)   
   slc.set_zlim("ZeusHeatingRate", 1e-40,1e-20)
#   slc.set_zlim("Mach", 1, 1.2) # could not find Mach
   slc.set_zlim("mach_number", 1, 1.2)
   slc.save()

   stest = 'stest_%04i' %i
   #prof=yt.create_profile(sp,"radius",["CoolingRateCell","ZeusHeatingCell","ZeusShockHeatingCell"], units={"radius": "kpc"}, logs = {"radius": True},accumulation=True,n_bins=n_bins,weight_field=None)
   prof=yt.create_profile(sp,"radius",["CoolingRateCell","ZeusHeatingCell","ZeusShockHeatingCell","Cell_ThermalEnergy"],units={"radius": "kpc"}, logs = {"radius": True},accumulation=False,n_bins=n_bins,weight_field=None,extrema={'radius':[rmin,rmax]})
   r[j,:]=numpy.array(prof.x.value)
   Cooling[j,:]=numpy.array(prof["CoolingRateCell"].value)
   Heating[j,:]=numpy.array(prof["ZeusHeatingCell"].value)
   ShockHeating[j,:]=numpy.array(prof["ZeusShockHeatingCell"].value)
   ShellThermal[j,:]=numpy.array(prof["Cell_ThermalEnergy"].value)
   #prof2=yt.create_profile(sp,"radius",["CoolingRate","ZeusHeatingRate","ZeusShockHeatingRate"], units={"radius": "kpc"}, logs = {"radius": True},accumulation=False, n_bins=n_bins,weight_field="cell_mass")
   prof2=yt.create_profile(sp,"radius",["CoolingRate","ZeusHeatingRate","ZeusShockHeatingRate"], units={"radius": "kpc"}, logs = {"radius": True},accumulation=False, n_bins=n_bins,weight_field="cell_volume",extrema={'radius':[rmin,rmax]})
   CoolingRate[j,:]=numpy.array(prof2["CoolingRate"].value)
   HeatingRate[j,:]=numpy.array(prof2["ZeusHeatingRate"].value)
   ShockHeatingRate[j,:]=numpy.array(prof2["ZeusShockHeatingRate"].value)

   #prof3=yt.create_profile(sp,"temperature",["CoolingRate","ZeusHeatingRate","ZeusShockHeatingRate"], accumulation=False,n_bins=n_bins,weight_field="ones")
   prof3=yt.create_profile(sp,"temperature",["CoolingRate","ZeusHeatingRate","ZeusShockHeatingRate","density"], accumulation=False,n_bins=n_bins,weight_field=("gas", "cell_mass"),extrema={'temperature':[Tmin,Tmax]}) #this should be done in a better way: first get lambda/n^2 field
   T[j,:]=numpy.array(prof3.x.value)
   CoolingRate_T[j,:]=numpy.array(prof3["CoolingRate"].value)
   HeatingRate_T[j,:]=numpy.array(prof3["ZeusHeatingRate"].value)
   ShockHeatingRate_T[j,:]=numpy.array(prof3["ZeusShockHeatingRate"].value)
   density_T[j,:]=numpy.array(prof3["density"].value)
#   TotalThermal.append(sp.quantities.total_quantity(["Cell_ThermalEnergy"]))
#   print "right", sp.quantities.total_quantity(["Cell_ThermalEnergy"])
#   print "wrong", sp.quantities.total_quantity(["thermal_energy"])
   TotalCooling.append(sp.quantities.total_quantity(["CoolingRate"]))
   TotalCoolingCell.append(sp.quantities.total_quantity(["CoolingRateCell"]))
#   plot = yt.PhasePlot(sp, "temperature", "CoolingRate", ["cell_mass"], weight_field=None)
#   plot = yt.PhasePlot(sp, "density", "temperature", ["cell_mass"],weight_field=None)
#   plot.set_xlim([1e-27,1e-24])
#   plot.set_ylim([1e7,1e9])
#   plot.save(stest+"rho_T_phase.png")
#   plot = yt.PhasePlot(sp, "density", "temperature", ["ones"], weight_field=None)
#   plot.save(stest+"rho_T_phase_ones.png")
#   plot=yt.PhasePlot(sp,"temperature","Cooling_Time",["cell_mass"], weight_field=None)
#   plot.save(stest+"T_tcool_phase.png")
   
#   sp_shock=sp.cut_region(["obj['ShockMach']>1.001"]) #cut_region used twice will give strange results 
   sp_shock=sp_all.cut_region(["obj['temperature'] > 1e6","obj['density']<1.0e-23","obj['ShockMach']>1.001"])
   sp_strong_shock=sp_all.cut_region(["obj['temperature'] > 1e6","obj['density']<1.0e-23","obj['ShockMach']>1.1"])
   prof_sp=yt.create_profile(sp,"radius",["cell_volume"], units={"radius": "kpc"}, logs = {"radius": True},accumulation=True,n_bins=n_bins,weight_field=None,extrema={'radius':[rmin,rmax]})
   volume[j,:]=numpy.array(prof_sp["cell_volume"].value)
   prof_shock=yt.create_profile(sp_shock,"radius",["cell_volume"], units={"radius": "kpc"}, logs = {"radius": True},accumulation=True,n_bins=n_bins,weight_field=None,extrema={'radius':[rmin,rmax]})
   shock_volume[j,:]=numpy.array(prof_shock["cell_volume"].value)
   prof_sshock=yt.create_profile(sp_strong_shock,"radius",["cell_volume"], units={"radius": "kpc"}, logs = {"radius": True},accumulation=True,n_bins=n_bins,weight_field=None,extrema={'radius':[rmin,rmax]})
   strong_shock_volume[j,:]=numpy.array(prof_sshock["cell_volume"].value)

#   print sp.quantities.total_quantity('cell_volume')
#   print "j=",j,"shock_volume", sp_shock.quantities.total_quantity('cell_volume')
#   print "strong_shock_volume", sp_strong_shock.quantities.total_quantity('cell_volume')
#   plot = yt.PhasePlot(sp_strong_shock, "density", "temperature", ["cell_mass"], weight_field=None)
#   plot.set_xlim([1e-27,1e-24])
#   plot.set_ylim([1e7,1e9])
#   plot.save(stest+"rho_T_phase_strong_shock.png")
   j=j+1

 numpy.savez("yt3_Heating_Zeus_%04i_%04i_%ikpc" %(imin,imax,rmax), r=r,Cooling=Cooling, Heating=Heating, ShockHeating=ShockHeating,ShellThermal=ShellThermal,CoolingRate=CoolingRate,HeatingRate=HeatingRate,ShockHeatingRate=ShockHeatingRate,CoolingRate_T=CoolingRate_T,HeatingRate_T=HeatingRate_T,density_T=density_T, ShockHeatingRate_T=ShockHeatingRate_T,T=T,time=time,TotalCooling=TotalCooling,TotalCoolingCell=TotalCoolingCell, Edot=Edot,volume=volume,shock_volume=shock_volume,strong_shock_volume=strong_shock_volume)
 return

from yt import derived_field
@derived_field(name="Vx_squared", units="cm**2/s**2")
def _Vx_squared(field, data):
    return data[('enzo', 'x-velocity')].in_cgs()**2


def Vmap(imin,imax=None,mod=1,rmax=100):
 if imax is None:
   imax=imin+1
 for i in range(imin,imax,mod):
   fn = name_pattern %(i,i)
   stest = 'stest_%04i' %i
   ds = yt.load(fn)
   sp_all=ds.sphere([0.5, 0.5, 0.5], (rmax, 'kpc'))
   sp1=sp_all.cut_region(["obj[('enzo', 'x-velocity')] > 0","obj['temperature']<1e5"])
   #proj = yt.ProjectionPlot(ds, 'x', ('enzo', 'x-velocity'), weight_field="density", data_source=sp1)
   proj = yt.ProjectionPlot(ds, 'x', [('enzo', 'x-velocity'),'Vx_squared'], weight_field="density", data_source=sp1)
   proj_frb=proj.data_source.to_frb((rmax, "kpc"), 512)
   proj_vx1=numpy.array(proj_frb[('enzo', 'x-velocity')])
   proj_Vsquared1=numpy.array(proj_frb['Vx_squared'])

   sp2=sp_all.cut_region(["obj[('enzo', 'x-velocity')] < 0","obj['temperature']<1e5"])
   proj = yt.ProjectionPlot(ds, 'x', [('enzo', 'x-velocity'),'Vx_squared'], weight_field="density",data_source=sp2)
   proj_frb=proj.data_source.to_frb((rmax, "kpc"), 512)
   proj_vx2=numpy.array(proj_frb[('enzo', 'x-velocity')])
   proj_Vsquared2=numpy.array(proj_frb['Vx_squared'])

   slc=yt.SlicePlot(ds,'x',[('enzo', 'y-velocity'),('enzo', 'z-velocity')], data_source=sp_all)
   slc_frb=slc.data_source.to_frb((rmax, "kpc"), 512)
   slc_vy=numpy.array(slc_frb[('enzo', 'y-velocity')])
   slc_vz=numpy.array(slc_frb[('enzo', 'z-velocity')])
   p=yt.SlicePlot(ds, 'x', "density", width = (rmax*2, 'kpc'))
   p.annotate_velocity(factor = 16)
   p.save()
 numpy.savez(stest+"test_v", proj_vx1=proj_vx1, proj_vx2=proj_vx2,proj_Vsquared1=proj_Vsquared1,proj_Vsquared2=proj_Vsquared2,slc_vy=slc_vy,slc_vz=slc_vz)
 return

def Mark_K_log(imin=0,imax=None,rmin=1,rmax=100,mod=1,y_min=1,y_max=100):
 kev_to_erg=6.2415e8
 codetime=0.048108   ##0036
 #codetime=0.1132
 if imax is None:
    imax=imin+1
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = yt.load(fn)
    sp = pf.sphere([0.5, 0.5, 0.5], (rmax, 'kpc'))
    hot_shell=sp.cut_region(["obj[('index', 'radius')].in_units('kpc') > 1","obj[('gas','temperature')]>1e5"])
    #prof2D=yt.create_profile(data_source=hot_shell,bin_fields=["radius","entropy"],fields=["cell_mass"],weight_field=None,n_bins=128,extrema=dict(radius=(1.0,rmax),entropy=(y_min/kev_to_erg,y_max/kev_to_erg)),units=dict(radius="kpc",cell_mass="Msun"),logs=dict(radius=True,entropy=True))
    prof2D=yt.create_profile(data_source=hot_shell,bin_fields=["radius","entropy"],fields=["cell_mass"],weight_field=None,n_bins=128,extrema=dict(radius=(1.0,rmax),entropy=(y_min,y_max)),units=dict(radius="kpc",cell_mass="Msun"),logs=dict(radius=True,entropy=True))
#    prof2D=yt.create_profile(data_source=hot_shell, bin_fields=[('index', 'radius'),('gas','entropy')], fields=[('gas','cell_mass')],weight_field=None,n_bins=128,extrema=dict(radius=(rmin,rmax), entropy=(y_min/kev_to_erg,y_max/kev_to_erg)), units={'radius':'kpc','cell_mass':'Msun'},logs={'radius': True,'entropy',True})
    mm=prof2D[('gas','cell_mass')]
    logmm=numpy.log10(mm)
    TT=logmm.T
##begin percentage calculation
    mm_cumsum=numpy.cumsum(mm,axis=1)
    mm_percentage=numpy.zeros(mm_cumsum.shape)
    lines=numpy.empty([128,3])
    y_locations=numpy.linspace(numpy.log10(y_min),numpy.log10(y_max),num=128)
    for j in xrange(128):
       mm_percentage[j,:]=mm_cumsum[j,:]/mm_cumsum[j,-2]
       for k in xrange(3):
           index = (numpy.abs(mm_percentage[j,:] - percentiles[k])).argmin()
           lines[j,k] = y_locations[index]
           k+=1
       j+=1
##end percentage calculation. Now we have "lines" to plot.
    pylab.figure(i+1000)
    #pylab.imshow(TT,interpolation='nearest',extent=[rmin,rmax,y_min,y_max],aspect='auto',origin='lower',cmap='Oranges')
    pylab.imshow(TT,interpolation='nearest',extent=[numpy.log10(rmin),numpy.log10(rmax),numpy.log10(y_min),numpy.log10(y_max)],aspect='auto',origin='lower',cmap='Oranges')
    #pylab.imshow(TT,interpolation='nearest',aspect='auto',origin='lower',cmap='Blues')
##plot percentile lines:
    x_locations=numpy.linspace(numpy.log10(rmin),numpy.log10(rmax),num=128)
    pylab.xlim([numpy.log10(rmin),numpy.log10(rmax)])
    #pylab.ylim([y_min,y_max])
    pylab.ylim([numpy.log10(y_min),numpy.log10(y_max)])
    y_locs=[0,1,2]
    y_labels=[0,1,2]
    pylab.yticks(y_locs, y_labels)
    x_locs=[1,2]
    x_labels=[1,2]
    pylab.xticks(x_locs, x_labels)
    time=codetime*TimeUnits/3.16e16   ## code unit to Gyr
    pylab.text(1.9, 0.3, "t = %.2f Gyr"%time, fontdict={'size':24, 'color':'black'})
    pylab.xlabel('log r (kpc)')
    #pylab.ylabel("Entropy K (keV "+r"$cm^2$"+")")
    pylab.ylabel("log K (keV "+r"$\rm{cm}^2$"+")")
    pylab.tight_layout()
    pylab.savefig(stest+'Mark_K_2D.png')
    #pylab.savefig(stest+'Mark_K_2D.pdf',format='pdf')
    for k in xrange(3):
        good=lines[:,k]>1e-6
        bad = lines[:,k]<1e-6
        lines[bad,k]=numpy.interp(x_locations[bad],x_locations[good],lines[good,k])
        pylab.plot(x_locations,lines[:,k],color="k",linestyle='--',linewidth=1)
        #pylab.plot(x_locations,lines_good,color="k",linestyle='--')
    pylab.plot(x_locations,lines[:,1],color="k",linewidth=1)
    pylab.savefig(stest+'Mark_K_2D_lines.png')
    pylab.savefig(stest+'Mark_K_2D_lines.pdf',format='pdf')
    codetime+=0.0017
 pylab.clf()
 return


def Mark_rho(imin=0,imax=None,rmin=1,rmax=100,mod=1,y_min=2e-27,y_max=1e-23,percentiles=[0.2,0.5,0.8]):
 kev_to_erg=6.2415e8
 codetime=0.048108   ##0036
 if imax is None:
    imax=imin+1
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = yt.load(fn)
    sp = pf.sphere([0.5, 0.5, 0.5], (rmax, 'kpc'))
    hot_shell=sp.cut_region(["obj[('index', 'radius')].in_units('kpc') > 1","obj[('gas','temperature')]>1e5"])
    prof2D=yt.create_profile(data_source=hot_shell,bin_fields=["radius","density"],fields=["cell_mass"],weight_field=None,n_bins=128,extrema=dict(radius=(1.0,rmax),density=(y_min,y_max)),units=dict(radius="kpc",cell_mass="Msun"),logs=dict(radius=True,density=True))
    mm=prof2D[('gas','cell_mass')]
    logmm=numpy.log10(mm)
    TT=logmm.T
##begin percentage calculation
    mm_cumsum=numpy.cumsum(mm,axis=1)
    mm_percentage=numpy.zeros(mm_cumsum.shape)
    lines=numpy.empty([128,3])
    y_locations=numpy.linspace(numpy.log10(y_min),numpy.log10(y_max),num=128)
    for j in xrange(128):
       mm_percentage[j,:]=mm_cumsum[j,:]/mm_cumsum[j,-2]
       for k in xrange(3):
           index = (numpy.abs(mm_percentage[j,:] - percentiles[k])).argmin()
           lines[j,k] = y_locations[index]
##end percentage calculation. Now we have "lines" to plot.
    pylab.figure(i+1000)
    pylab.imshow(TT,interpolation='nearest',extent=[numpy.log10(rmin),numpy.log10(rmax),numpy.log10(y_min),numpy.log10(y_max)],aspect='auto',origin='lower',cmap='Oranges')
##plot percentile lines:
    x_locations=numpy.linspace(numpy.log10(rmin),numpy.log10(rmax),num=128)
    pylab.xlim([numpy.log10(rmin),numpy.log10(rmax)])
    pylab.ylim([numpy.log10(y_min),numpy.log10(y_max)])
    y_locs=[-26,-25,-24]
    y_labels=[-26,-25,-24]
    pylab.yticks(y_locs, y_labels)
    x_locs=[1,2]
    x_labels=[1,2]
    pylab.xticks(x_locs, x_labels)
    time=codetime*TimeUnits/3.16e16   ## code unit to Gyr
    pylab.text(1.9, 0.3, "t = %.2f Gyr"%time, fontdict={'size':24, 'color':'black'})
    pylab.xlabel('log r (kpc)')
    pylab.ylabel('log'+r'$\rm{\rho}(\rm{g}\ /\ \rm{cm}^3)$')
    pylab.tight_layout()
#    pylab.savefig(stest+'Mark_rho_2D.png')
#this is cheating but it is ok:
    lines[-1,0]=lines[-1,1]
    for k in xrange(3):
        good=lines[:,k]>numpy.log10(y_min)+1.0e-3
        bad = lines[:,k]<numpy.log10(y_min)+1.0e-3
        lines[bad,k]=numpy.interp(x_locations[bad],x_locations[good],lines[good,k])
        pylab.plot(x_locations,lines[:,k],color="k",linestyle='--',linewidth=1)
    pylab.plot(x_locations,lines[:,1],color="k",linewidth=1)
#    pylab.savefig(stest+'Mark_rho_2D_lines.png')
    codetime+=0.0017
##begin all time average
    if i==imin:
      mm_all=mm
    else:
      mm_all=numpy.add(mm,mm_all)
    if i==(imax-1):
      logmm_all=numpy.log10(mm_all)
      TT_all=logmm_all.T
      mm_cumsum_all=numpy.cumsum(mm_all,axis=1)
      for j in xrange(128):
         mm_percentage[j,:]=mm_cumsum_all[j,:]/mm_cumsum_all[j,-2]  #new mm_percentage
         for k in xrange(3):
             index = (numpy.abs(mm_percentage[j,:] - percentiles[k])).argmin()
             lines[j,k] = y_locations[index] #new lines for all time
      pylab.clf()
      pylab.imshow(TT_all,interpolation='nearest',extent=[numpy.log10(rmin),numpy.log10(rmax),numpy.log10(y_min),numpy.log10(y_max)],aspect='auto',origin='lower',cmap='Oranges')
      pylab.xlim([numpy.log10(rmin),numpy.log10(rmax)])
      pylab.ylim([numpy.log10(y_min),numpy.log10(y_max)])
      y_locs=[-26,-25,-24]
      y_labels=[-26,-25,-24]
      pylab.yticks(y_locs, y_labels)
      x_locs=[1,2]
      x_labels=[1,2]
      pylab.xticks(x_locs, x_labels)
      pylab.xlabel('log r (kpc)')
      pylab.ylabel('log'+r'$\mathbf{\rm{\rho}(\rm{g}\ /\ \rm{cm}^3)}$')
      pylab.tight_layout()
      pylab.savefig('All_time_Mark_rho_2D.png')
      lines[-1,0]=lines[-1,1]
      for k in xrange(3):
        good=lines[:,k]>numpy.log10(y_min)+1.0e-3
        bad = lines[:,k]<numpy.log10(y_min)+1.0e-3
        lines[bad,k]=numpy.interp(x_locations[bad],x_locations[good],lines[good,k])
        pylab.plot(x_locations,lines[:,k],color="k",linestyle='--',linewidth=1)
      pylab.plot(x_locations,lines[:,1],color="k",linewidth=1)
      pylab.savefig('All_time_Mark_rho_2D_lines.png')
      pylab.clf()
#electron density
      pylab.imshow(TT_all,interpolation='nearest',extent=[numpy.log10(rmin),numpy.log10(rmax),numpy.log10(y_min),numpy.log10(y_max)],aspect='auto',origin='lower',cmap='Oranges')
      pylab.xlim([numpy.log10(rmin),numpy.log10(rmax)])
      pylab.ylim([numpy.log10(y_min),numpy.log10(y_max)])
      y_locs=numpy.array([-2,-1,0])+numpy.log10(0.59*m_p)
      y_labels=[-2,-1,0]
      pylab.yticks(y_locs, y_labels)
      x_locs=[1,2]
      x_labels=[1,2]
      pylab.xticks(x_locs, x_labels)
      pylab.xlabel('log r (kpc)')
      pylab.ylabel('log'+r'$\mathbf{\rm{\rho}(\rm{g}\ /\ \rm{cm}^3)}$')
      pylab.tight_layout()
      pylab.savefig('All_time_Mark_rho_2D.png')
      lines[-1,0]=lines[-1,1]
      for k in xrange(3):
        good=lines[:,k]>numpy.log10(y_min)+1.0e-3
        bad = lines[:,k]<numpy.log10(y_min)+1.0e-3
        lines[bad,k]=numpy.interp(x_locations[bad],x_locations[good],lines[good,k])
        pylab.plot(x_locations,lines[:,k],color="k",linestyle='--',linewidth=1)
      pylab.plot(x_locations,lines[:,1],color="k",linewidth=1)
      pylab.savefig('All_time_Mark_ne_2D_lines.png')
 pylab.clf()
 return

#this is to calculate the change in total thermal energy and compare with Zeus
@derived_field(name="Cell_KineticEnergy", units= "erg")
def _Cell_KineticEnergy(field,data):
   return data["kinetic_energy"]*data["cell_volume"]

def Zeus_test(imin=0,imax=1):
 name_pattern='DD%04i/data%04i'
 Total_th=[]
 Total_k=[]
 Time=[]
 Zeus_Rate=[]
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   pf = yt.load(fn)
   pf.periodicity=(True,True,True)
   box = pf.all_data()
#   box=pf.region([0.5, 0.5, 0.5], [0, 0, 0], [1,1,1])
   Time.append(pf.current_time)
   Zeus_Rate.append(box.quantities.total_quantity(["ZeusHeatingCell"]))
   Total_th.append(box.quantities.total_quantity(["Cell_ThermalEnergy"]))
   Total_k.append(box.quantities.total_quantity(["Cell_KineticEnergy"]))
 numpy.savez("Zeus_Test",Time=Time,Zeus_Rate=Zeus_Rate,Total_th=Total_th,Total_k=Total_k)
 return



def Level_profile(i=0):
 fn = name_pattern %(i,i)
 ds=yt.load(fn)
 sp=ds.sphere("c",(100,"kpc"))
 yt.ProfilePlot(sp,"radius","grid_level",weight_field="cell_mass").set_unit("radius","kpc").save()
 return



def simple_projections(imin=0,imax=1,widthkpc=None,p_center=[0.2,0.5,0.5]):
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   ds = yt.load(fn)
   P=yt.ProjectionPlot(ds, "y", "temperature", weight_field="density")
   if widthkpc is not None:
    P.set_width(widthkpc,"kpc")
#   P.set_zlim("temperature",2e7,6e7).save()
   P.set_zlim("temperature",1e4,2e6).save()
 return

def Marc_projection(imin=212,imax=213,widthkpc=50,center=[0.5,0.5,0.5]):
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   ds = yt.load(fn)
   P=yt.OffAxisProjectionPlot(ds, [1,0,0], "density", weight_field="density",center=center,north_vector=[0,0,-1])
   P.set_width(widthkpc,"kpc").save()
#   P.set_zlim("temperature",2e7,6e7).save()
#   P.set_zlim("temperature",1e4,2e6).save()
 return

def Mira_ray(imin=0,imax=None,start=[0.2,0.5,0.5],end=[0.2,0.5,1]):
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   ds = yt.load(fn)
   #ray=ds.r[start:end] #first a ray along z
   ray=ds.ray(start, end)
   ray_t=numpy.array(ray["t"])
   print(ray["t"])
   ray_density=numpy.array(ray["density"])
   ray_pressure=numpy.array(ray["pressure"])
   ray_temperature=numpy.array(ray["temperature"])
   ray_vz=numpy.array(ray["z-velocity"])
   ray_z=numpy.array(ray["z"])
   numpy.savez(stest+"_ray_z.npz",ray_density=ray_density,ray_pressure=ray_pressure,ray_temperature=ray_temperature,ray_z=ray_z,ray_vz=ray_vz,ray_t=ray_t)

   end=[0,0.5,0.5] #now a ray along -x
   ray=ds.r[start:end]
   ray_t=numpy.array(ray["t"])
   ray_density=numpy.array(ray["density"])
   ray_pressure=numpy.array(ray["pressure"])
   ray_temperature=numpy.array(ray["temperature"])
   ray_vx=numpy.array(ray["x-velocity"])
   ray_x=numpy.array(ray["x"])
   numpy.savez(stest+"_ray_x.npz",ray_density=ray_density,ray_pressure=ray_pressure,ray_temperature=ray_temperature,ray_x=ray_x,ray_vx=ray_vx,ray_t=ray_t)

 return




@derived_field(name="Color_Fraction", units= "")
def _Color_Fraction(field,data):
   return data["SN_Colour"]*0.496/data["density"]

@derived_field(name="cold_cell_mass", units="g")
def _cold_cell_mass(field,data):
   ret=data[("gas","cell_mass")]
   ret[data["temperature"].in_cgs()>2e4]=0
   return ret



def Mira_phase(imin=0,imax=1,p_center=[0.2,0.5,0.5]):
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   stest = 'stest_%04i' %i
   ds=yt.load(fn)
   box = ds.box([0,0,0],[4,1,1])
#   box1=ds.box([0.,0.25,0.25],[3,0.75,0.75])
#   box=box1.cut_region(["obj['SN_Colour'] > 0.1"])
#   box=box1
   T_phase=yt.PhasePlot(box,"x","temperature","cell_mass",weight_field=None)
   T_phase.set_log("x",False)
   T_phase.set_unit("x","pc")
   T_phase.set_ylim(1e3,2e6)
   T_phase.set_xlim(0,8)
   T_phase.save()
   tcool_color=yt.PhasePlot(box,"SN_Colour","cooling_time","cell_mass",weight_field=None)
   tcool_color.set_xlim(1e-2,100)
   tcool_color.set_ylim(5e9,3e14)
#   tcool_color.set_unit("cooling_time","Myr")
   tcool_color.save()
   tcool_Fcolor=yt.PhasePlot(box,"Color_Fraction","cooling_time","cell_mass",weight_field=None)
   tcool_Fcolor.set_xlim(1e-2,2)
   tcool_Fcolor.set_ylim(5e9,3e14)
   tcool_Fcolor.save()
   T_color=yt.PhasePlot(box,"SN_Colour","temperature","cell_mass",weight_field=None)
   T_color.set_xlim(1e-2,100)
   T_color.set_ylim(300,1.3e6)
   T_color.save()
   rho_color=yt.PhasePlot(box,"SN_Colour","density","cell_mass",weight_field=None)
   rho_color.set_xlim(1e-2,100)
   #rho_color.set_ylim(300,1.3e6)
   rho_color.save()

   x_v=yt.PhasePlot(box,"x",("enzo", "x-velocity"),"cooling_time",weight_field="cell_mass")
   x_v.set_log("x",False)
   x_v.set_log("x-velocity",False)
   x_v.set_unit("x","pc")
   x_v.set_unit("x-velocity","km/s")
   x_v.set_unit("cooling_time","Myr")
   x_v.set_xlim(0,8)
   x_v.set_ylim(-50,200)
   x_v.save()
   x_v2=yt.PhasePlot(box,"x",("enzo", "x-velocity"),"cell_mass",weight_field=None)
   x_v2.set_log("x",False)
   x_v2.set_log("x-velocity",False)
   x_v2.set_unit("x","pc")
   x_v2.set_unit("x-velocity","km/s")
   x_v2.set_xlim(0,8)
   x_v2.set_ylim(-50,200)
   x_v2.save()
#   small_box=ds.box([0.25,0.25,0.25],[3,0.75,0.75])
   small_box=box
   rho_v=yt.PhasePlot(small_box,("enzo", "x-velocity"),"density","cell_mass",weight_field=None)
#   rho_v.set_log("x-velocity",False)
   rho_v.set_unit("x-velocity","km/s")
   rho_v.set_xlim(1,200)
   rho_v.set_ylim(3e-26,3e-22)
   rho_v.save()   
   rho_T=yt.PhasePlot(small_box,"temperature","density","cell_mass",weight_field=None) 
   rho_T.set_ylim(3e-26,3e-22)
   rho_T.set_xlim(300,1.3e6)
   rho_T.save()
   P_K=yt.PhasePlot(small_box,"pressure","entropy","cell_mass",weight_field=None)
#   P_K.set_xlim(
   P_K.save()
   entropy_tcool=yt.PhasePlot(small_box,"cooling_time","entropy","cell_mass",weight_field=None)
   entropy_tcool.save()   
   T_Fcolor=yt.PhasePlot(box,"Color_Fraction","temperature","cell_mass",weight_field=None)
   T_Fcolor.set_xlim(0.01,3)
   T_Fcolor.save()

   rho_Fcolor=yt.PhasePlot(box,"Color_Fraction","density","cell_mass",weight_field=None)
   rho_Fcolor.set_xlim(0.01,3)
   rho_Fcolor.save()

   x_v3=yt.PhasePlot(box,"x",("enzo", "x-velocity"),"CoolingDistance",weight_field="cell_mass").set_ylim(-50,200)
   x_v3.set_log("x",False)
   x_v3.set_log("x-velocity",False)
   x_v3.set_unit("x","pc")
   x_v3.set_unit("x-velocity","km/s")
   x_v3.set_unit("CoolingDistance","pc")
   x_v3.set_xlim(0,8)
   x_v3.set_ylim(-50,200)
   x_v3.save()


   prof_b=yt.create_profile(box,["Color_Fraction","density"],('gas', 'cell_mass'), weight_field=None,extrema={"Color_Fraction":[1e-2,2.5], "density":[2e-26,1e-21]})
   prof_b_x=prof_b.x
   prof_b_y=prof_b.y
   prof_b_prof=prof_b['gas', 'cell_mass']



   #prof2d=yt.create_profile(box, [("enzo", "x-velocity"),("gas", "density")], ('gas', 'cell_mass'), weight_field=None,units={"x-velocity": "km/s"},extrema={"x-velocity":[1,170]}) 
   prof2d=yt.create_profile(box, [("enzo", "x-velocity"),("gas", "density")], ('gas', 'cell_mass'), weight_field=None,logs = {"x-velocity":False},units={"x-velocity": "km/s"},extrema={("enzo","x-velocity"):[1,170], "density":[2e-26,1e-21]}) 
#   prof2d=yt.create_profile(box, [("enzo", "x-velocity"),("gas", "density")], ('gas', 'cell_mass'), weight_field=None,units={"x-velocity": "km/s"},extrema=dict(x-velocity=(1,170)))
   v=prof2d.x
   rho=prof2d.y
   rho_v_profile=prof2d['gas', 'cell_mass']

   prof_c=yt.create_profile(box,["x",("enzo", "x-velocity")],('gas', 'cold_cell_mass'), weight_field=None,logs = {"x": False, "x-velocity":False},units={"x-velocity": "km/s","x":"pc"},extrema={("enzo","x-velocity"):[1,170],"x":[0,8]})
   prof_c_x=prof_c.x
   prof_c_y=prof_c.y
   prof_c_prof=prof_c['gas', 'cold_cell_mass']
   prof_d=yt.create_profile(box,["x",("gas", "temperature")],('gas', 'cell_mass'), weight_field=None,logs = {"x": False, "temperature":True},units={"temperature": "K","x":"pc"},extrema={"x":[0,8],"temperature":[100,1.5e6]})
   prof_d_x=prof_d.x
   prof_d_y=prof_d.y
   prof_d_prof=prof_d['gas', 'cell_mass']

   prof_e=yt.create_profile(box,[("gas","density"),("gas", "temperature")],('gas', 'cell_mass'), weight_field=None)
   prof_e_x=prof_e.x
   prof_e_y=prof_e.y
   prof_e_prof=prof_e['gas', 'cell_mass']

   numpy.savez(stest+"phase",v=v,rho=rho,rho_v_profile=rho_v_profile,prof_b_x=prof_b_x,prof_b_y=prof_b_y,prof_b_prof=prof_b_prof,prof_c_x=prof_c_x,prof_c_y=prof_c_y,prof_c_prof=prof_c_prof,prof_d_x=prof_d_x,prof_d_y=prof_d_y,prof_d_prof=prof_d_prof,prof_e_x=prof_e_x,prof_e_y=prof_e_y,prof_e_prof=prof_e_prof)
 return

#define flux (per cell) 
@derived_field(name="Mira_hot_flux",units="g/s")
def _Mira_hot_flux(field, data):
    ret = data['dx'].in_cgs()**2*data['density'].in_cgs()*data['x-velocity'].in_cgs()
    ret[data['temperature'] < 1e6] = 0
    return ret

@derived_field(name="Mira_cold_flux",units="g/s")
def _Mira_hot_flux(field, data):
    ret = data['dx'].in_cgs()**2*data['density'].in_cgs()*data['x-velocity'].in_cgs()
    ret[data['temperature'] > 2e4] = 0
    return ret

@derived_field(name="Mira_all_flux",units="g/s")
def _Mira_hot_flux(field, data):
    ret = data['dx'].in_cgs()**2*data['density'].in_cgs()*data['x-velocity'].in_cgs()
#    ret[data['temperature'] < 9e5] = 0
    return ret

# try another way of computing flux: do cell_area weighted rho_v and then multiply by total area
@derived_field(name="cell_area", units="cm**2")
def _cell_area(field,data):
    return data['dx'].in_cgs()**2

@derived_field(name="Mira_rho_v",units="g/s/cm**2")
def _Mira_rho_v(field,data):
    return data['density'].in_cgs()*data['x-velocity'].in_cgs()

@derived_field(name="Mira_cold_rho_v",units="g/s/cm**2")
def _Mira_rho_v(field,data):
    ret = data['density'].in_cgs()*data['x-velocity'].in_cgs()
    ret[data['temperature'] > 2e4] = 0
    return ret

def Mira_flux_b(imin=0,imax=1,do_all=False,nbins=64):
 LengthUnits  = 6.1700002089743e+18
 all_time_flux=numpy.zeros((imax-imin,nbins))
 all_time_mass=numpy.zeros((imax-imin,nbins))
 for i in range(imin,imax):
    fn = name_pattern %(i,i)
    stest = 'stest_%04i' %i
    ds = yt.load(fn)
    ds.periodicity=(True,True,True)
    box = ds.box([0,0,0],[4,1,1])
    prof=yt.create_profile(box,"x",["Mira_rho_v","Mira_cold_rho_v"],logs = {"x": False}, accumulation=False, weight_field="cell_area", units={"x":"pc"},n_bins=nbins,extrema={'x':[0,8]})
    all_flux=numpy.array(prof["Mira_rho_v"])*LengthUnits**2/(2.0e33/3.15e13)
    cold_flux=numpy.array(prof["Mira_cold_rho_v"])*LengthUnits**2/(2.0e33/3.15e13)
    r=numpy.array(prof.x)
    box_cold=box.cut_region(["obj['temperature'] < 2e4"])
    prof_cold=yt.create_profile(box,"x",["cold_cell_mass","Mira_rho_v"],logs = {"x": False}, accumulation=False, weight_field=None, units={"x":"pc"},n_bins=nbins,extrema={'x':[0,8]})
    cold_mass=numpy.array(prof_cold["cold_cell_mass"])
    #prof_cold=yt.create_profile(box_cold,"x",[("gas","cell_mass"),"Mira_rho_v"],logs = {"x": False}, accumulation=False, weight_field=None, units={"x":"pc"},n_bins=nbins,extrema={'x':[0,8]})
    #cold_mass=numpy.array(prof_cold["cell_mass"])
    pylab.figure(i)
    pylab.plot(r,all_flux,color="red")
    pylab.plot(r,cold_flux,color="blue")
    pylab.xlabel("x (pc)")
    pylab.ylabel("Flux (M_sun/Myr)")
    pylab.legend(("total","cold"),loc="best")
    pylab.grid()
    pylab.tight_layout()
    pylab.savefig(stest+"_flux_b.png")
    if do_all==True:
       all_time_flux[i-imin,:]=cold_flux
       all_time_mass[i-imin,:]=cold_mass
    else:
       numpy.savez(stest+"flux.npz", r=r,all_flux=all_flux,cold_flux=cold_flux,cold_mass=cold_mass)
#    sp = ds.sphere([0.2,0.5,0.5],(0.1,"code_length")) #0.02 is r_w
#    surf = ds.surface(sp,"radius",(0.01,"code_length"))
#    flux = surf.calculate_flux("x-velocity","y-velocity","z-velocity","density")
#    flux_print=flux*ds.mass_unit/(2.0e33/3.15e13)
#    print("flux out of the star is ", flux,flux_print)
#    slc = yt.SlicePlot(ds, "x", ["density","x-velocity","Mira_rho_v"],center=[0.5,0.5,0.5],data_source=box_cold)
#    #slc = yt.SlicePlot(ds, "x", ["density","x-velocity","Mira_cold_rho_v"],center=[0.5,0.5,0.5],data_source=box)
#    width=(0.002,"kpc")
#    res=[100,100]
#    slc_frb=slc.data_source.to_frb(width, res, center=[0.5,0.5,0.5])
#    slc_vx=numpy.array(slc_frb["x-velocity"])
#    slc_density=numpy.array(slc_frb["density"])
#    slc_cold_rho_v=numpy.array(slc_frb["Mira_rho_v"])
#    numpy.savez(stest+"slice_frb.npz",slc_vs=slc_vx,slc_density=slc_density,slc_cold_rho_v=slc_cold_rho_v)

#    box = ds.all_data() 
#    surf_x=ds.surface(box,"x",(0.5-numpy.spacing(0.5),"code_length"))
#    print(numpy.sum(surf_x["density"]),"is the sum of density")
#    print(numpy.sum(surf_x["cell_area"]),"is the sum of cell area")
#    flux_x=surf_x.calculate_flux("x-velocity","y-velocity","z-velocity","density")
#    flux_x_print=flux_x #*ds.mass_unit/(2.0e33/3.15e13)
#    print("flux in the box is ", flux_x_print)
 if do_all == True:
    numpy.savez(stest+"all_time_flux.npz",all_time_flux=all_time_flux,all_time_mass=all_time_mass,r=r)
 return



def Mira_flux(imin=0,imax=1):  #this is wrong (normalization is off). Use flux_b instead
 LengthUnits  = 6.1700002089743e+18
 for i in range(imin,imax):
    fn = name_pattern %(i,i)
    stest = 'stest_%04i' %i
    ds = yt.load(fn)
    box = ds.box([0,0,0],[4,1,1])
#    #prof=yt.Profile1D(box,"x",100, 0.3, 4,x_log=False,weight_field=None)  # this does not work. Must be a bug 
#    prof.add_fields(["Mira_hot_flux","Mira_cold_flux","Mira_all_flux",("gas","cell_mass")])
#    prof.set_x_unit("pc")
    prof=yt.create_profile(box,"x",["Mira_hot_flux","Mira_cold_flux","Mira_all_flux",("gas","cell_mass")],logs = {"x": False}, accumulation=False, weight_field=None, units={"x":"pc"})
    box_cold=box.cut_region(["obj['temperature'] < 2e4"])
 #   box_cold=box.cut_region("temperature<2.0e4")
    prof_cold=yt.create_profile(box_cold,"x",("gas","cell_mass"),logs = {"x": False}, accumulation=False, weight_field=None, units={"x":"pc"})
    r=numpy.array(prof.x)
    hot_flux=numpy.array(prof["Mira_hot_flux"])/(2.0e33/3.15e13)
    cold_flux=numpy.array(prof["Mira_cold_flux"])/(2.0e33/3.15e13)
    all_flux=numpy.array(prof["Mira_all_flux"])/(2.0e33/3.15e13)
    mass=numpy.array(prof["cell_mass"])
    cold_mass=numpy.array(prof_cold["cell_mass"])
    pylab.figure(i)
    print("all_flux",all_flux)
    print(("gas","cell_mass"),mass)
    #pylab.semilogy(r,all_flux,color="red")
    pylab.plot(r,all_flux,color="k")
#    pylab.plot(r,cold_flux+hot_flux)
    pylab.plot(r,cold_flux,color="blue")
#    pylab.plot(r,hot_flux,color="red")
    pylab.xlabel("x (pc)")
    pylab.ylabel("Flux (M_sun/Myr)")
#    pylab.fill_between(r,
    pylab.legend(("total","cold"),loc="best")
    pylab.tight_layout()
    #pylab.show()
    pylab.savefig(stest+"_flux.png")
    pylab.clf()
    pylab.plot(r,cold_mass)
    pylab.xlabel("x (pc)")
    pylab.ylabel("Mass")
    pylab.savefig(stest+"mass_r.png")
    numpy.savez(stest+"flux.npz",r=r,cold_mass=cold_mass)

 return


def Mira_3D(imin=0,imax=1):
 for i in range(imin,imax):
   fn=name_pattern %(i,i)
   stest= 'stest_%04i' %i
   ds=yt.load(fn)
   ds.periodicity=(True,True,True)
   small_box=ds.box([0.,0,0],[3,1,1])
   ad = ds.all_data()
   surf=ds.surface(small_box,"density",2e-24)
   surf.export_obj("Mira_3D",color_field="temperature",color_map="hot",color_log=True)
#   upload_id = surf.export_sketchfab(title="Mira",description="Mira",color_field="temperature",color_map="hot",color_log=True)
 return



def Mira_Shock(imin=0,imax=1,useSam=False):
 time=[] 
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   ds = yt.load(fn)
   ds.periodicity=(True,True,True)
   time.append(ds.current_time)
   ds.add_gradient_fields(("gas","temperature"))
   ds.add_gradient_fields(("gas","entropy"))
   ds.add_gradient_fields(("gas","logP"))

   def _ShockMach(field,data):  #value is shock Mach number
     ep=0.002  # Mach of 1.001
#     ep=0.223  # Mach of 1.1
     new_field = numpy.zeros(data["pressure"].shape, dtype='float64')
     dog=numpy.maximum(numpy.array(data["logP_gradient_x"])*numpy.array(data["index", "dx"].in_cgs()),numpy.array(data["logP_gradient_y"])*numpy.array(data["index", "dy"].in_cgs()))
     dlogP=numpy.maximum(dog,numpy.array(data["logP_gradient_z"])*numpy.array(data["index", "dz"].in_cgs()))
     too_big=dlogP>2  # if the shock is too strong (mach > 5ish), mouse will cause overflow problem
     dlogP[too_big]=2
     if useSam==True:
       TSx=numpy.array(data['gas', 'temperature_gradient_x'])*numpy.array(data['gas', 'entropy_gradient_x'])
       TSy=numpy.array(data['gas', 'temperature_gradient_y'])*numpy.array(data['gas', 'entropy_gradient_y'])
       TSz=numpy.array(data['gas', 'temperature_gradient_z'])*numpy.array(data['gas', 'entropy_gradient_z'])
       TS1=numpy.logical_or((TSx>0), (TSy>0))
       TS=numpy.logical_or(TS1,  (TSz>0))
       good=(TS & (numpy.array(data["gas", "velocity_divergence"])<0)) & (dlogP > ep)
     else:
       cat=numpy.array(data["gas","thermal_energy"]*data["gas","cell_mass"])/numpy.array(data["gas","kinetic_energy"]*data["index","cell_volume"].in_cgs()) > (1.0/9.0)
       good=(cat & (numpy.array(data["gas", "velocity_divergence"])<0))& (dlogP > ep)
#     good=dlogP > ep
     mouse = (numpy.exp(dlogP)-1.0)*0.4 + 1
     new_field[good]=mouse[good]
     return new_field
   ds.add_field(('gas','ShockMach'), function=_ShockMach, units="", take_log=False,display_name='Shock Mach')
   slc=yt.SlicePlot(ds, "y",["ShockMach","entropy"])
   slc.set_zlim("ShockMach", 1, 1.05)
   slc.save()
 return

def Mira_projection(imin=0,imax=1):
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   stest = 'stest_%04i' %i
   ds = yt.load(fn)
   ds.periodicity=(True,True,True)
   box=ds.all_data()
   small_box = ds.box([0,0,0],[1,1,1])
   width=(0.002,"kpc")
   res=[1000,1000]

   proj = yt.ProjectionPlot(ds,"y",["density","SN_Colour"],weight_field=None, data_source=box)
   proj.save()
   proj_frb=proj.data_source.to_frb(width, res, center=[0.5,0.5,0.5])
   proj_Colour=numpy.array(proj_frb["SN_Colour"])
   proj_density=numpy.array(proj_frb["density"])
   proj_frb2=proj.data_source.to_frb(width, res, center=[1.5,1.5,1.5])
   proj_Colour=numpy.append(proj_Colour,numpy.array(proj_frb2["SN_Colour"]),axis=0)
   proj_density=numpy.append(proj_density,numpy.array(proj_frb2["density"]), axis=0)
   numpy.savez(stest+"proj_SN_Colour.npz", proj_Colour=proj_Colour,proj_density=proj_density)

 return


#from yt.units import kpc
@derived_field(name="rrvir",units = "")
def rrvir(field,data):
#    return data['radius'].in_units('kpc')/(218.,'kpc')
    return data['radius'].in_units('kpc')/data.ds.quan(700, 'kpc')

def AGN_Phase(imin=0,imax=1,rmax=244.0):
  for i in range(imin,imax):
    fn = name_pattern %(i,i)
    ds = yt.load(fn)
    stest = 'stest_%04i' %i
    sphere=ds.sphere([0.5, 0.5, 0.5], (rmax,'kpc'))
#    plot = yt.PhasePlot(sphere, "density", "temperature", "cell_mass", weight_field=None)
    plot = yt.PhasePlot(sphere, "temperature", "density", "cell_mass", weight_field=None, x_bins=512, y_bins=512).set_unit("density","g/cm**3").set_unit("temperature","K").set_unit("cell_mass","Msun").set_xlim(300,1e10).set_ylim(4e-28,4e-22)
    plot.set_unit('cell_mass', 'Msun')
    plot.set_zlim("cell_mass",1e3,1e10)
    plot.save()
#    X1=plot.profile.x.in_units("K").v
#    Y1=plot.profile.y.in_units("g/cm**3").v
#    Z1=plot.profile["cell_mass"].in_units("Msun").v
#    numpy.savez("phase_xyz.npz",X1=X1,Y1=Y1,Z1=Z1)
#    pylab.pcolormesh(X1,Y1,numpy.transpose(Z1), cmap="Oranges")#,norm=LogNorm(1e4, 1e10))
#    pylab.savefig(stest+"_rhoT.png",dpi=300)
    rho_r=yt.PhasePlot(sphere, "radius","density","cell_mass",weight_field=None, x_bins=512, y_bins=512).set_unit("density","g/cm**3").set_unit("cell_mass","Msun").set_unit("radius","kpc").set_xlim(1,rmax).set_ylim(4e-28,4e-22).set_zlim("cell_mass",1e3,1e10)
    rho=rho_r.profile.y.in_units("g/cm**3").v
    rr_rho=rho_r.profile["cell_mass"].in_units("Msun").v
    rho_r.save()
    T_r=yt.PhasePlot(sphere, "radius","temperature","cell_mass",weight_field=None, x_bins=512, y_bins=512).set_unit("cell_mass","Msun").set_unit("radius","kpc").set_xlim(1,rmax).set_ylim(300,1e10).set_zlim("cell_mass",1e3,1e10)
    rr2=T_r.profile.x.in_units("kpc").v
    T=T_r.profile.y.in_units("K").v
    rr_T=T_r.profile["cell_mass"].in_units("Msun").v
#    numpy.savez(stest+"phase.npz",rr=rr,rr2=rr2,rho=rho,T=T,rr_rho=rr_rho,rr_T=rr_T)
    T_r.save()
    dd=ds.all_data()
    rho_v=yt.PhasePlot(sphere,("gas","radial_velocity"),"density","cell_mass",weight_field=None, x_bins=512, y_bins=512).set_unit("density","g/cm**3").set_unit("cell_mass","Msun").set_log("radial_velocity",False)
    rho_v.set_xlim(-6e8,1e9)
    rho_v.set_zlim("cell_mass",1e3,1e10)
    rho_v.save()
    P_K=yt.PhasePlot(sphere,("gas","entropy"),"pressure","cell_mass",weight_field=None, x_bins=512, y_bins=512).set_unit("cell_mass","Msun").set_zlim("cell_mass",1e3,1e10)
    P_K.save()
    fP_r=yt.PhasePlot(sphere,"radius","fP","cell_mass",weight_field=None, x_bins=256, y_bins=256).set_unit("cell_mass","Msun").set_zlim("cell_mass",1e3,1e10).set_log("fP",False).set_unit("radius","kpc").set_xlim(1,rmax).set_ylim(0,1)
    fP_r.save()
    fP1=fP_r.profile["cell_mass"].in_units("Msun").v
    if i==imin:
       fP_all=fP1
    else:
       fP_all=numpy.add(fP1,fP_all)
    vr_r=yt.PhasePlot(sphere,"radius",("gas","radial_velocity"),"cell_mass",weight_field=None, x_bins=256, y_bins=256).set_unit("cell_mass","Msun").set_zlim("cell_mass",1e3,1e10).set_unit("radius","kpc").set_xlim(1,rmax).set_log("radial_velocity",False)
    vr_r.set_unit("radial_velocity","km/s")
    vr_r.set_ylim(-6e3,1e4)
    vr_r.save()
    vt_r=yt.PhasePlot(sphere,"radius",("gas","tangential_velocity"),"cell_mass",weight_field=None, x_bins=256, y_bins=256).set_unit("cell_mass","Msun").set_zlim("cell_mass",1e3,1e10).set_unit("radius","kpc").set_xlim(1,rmax).set_log("tangential_velocity",False)
    rr=vt_r.profile.x.in_units("kpc").v
    vt_r.save()
    vt=vt_r.profile["cell_mass"].in_units("Msun").v
    if i==imin:
       vt_all=numpy.median(vt,axis=0)
       vt_all=numpy.add(numpy.median(vt,axis=0),vt_all)
  numpy.savez("fP.npz",rr=rr,fP_all=fP_all,vt_all=vt_all)
  return


#fraction of non-thermal pressure
@derived_field(name="fP",units = "")
def rrvir(field,data):
    return data['pressure']/(data['density']*data['velocity_magnitude']**2+data['pressure'])



def Clump_Finding(imin=0,imax=None,mod=1,rmax=50):
 for i in range(imin,imax,mod):
   mass=[]
   fn = name_pattern %(i,i)
   stest = 'stest_%04i' %i
   pf = yt.load(fn)
   print("time=", pf.current_time)
   field = "density" # this is the field we look for contours over -- we could do
   step = 5.0 # This is the multiplicative interval between contours.
#   center=[0.5, 0.5, 0.5]
   #data_source=pf.sphere(center,20/pf['kpc'])
   sp_all = pf.sphere([0.5, 0.5, 0.5], (rmax, 'kpc'))
   sp = sp_all.cut_region(["obj['index','x']>0.5","obj['index','y']>0.5","obj['index','z']>0.5"])
   data_source = sp #use only 1/8 of the whole domain
#   c_min = 10**numpy.floor(numpy.log10(data_source[field]).min()  )
#   c_max = 10**numpy.floor(numpy.log10(data_source[field]).max()+1)
   master_clump = Clump(data_source, field)
   c_min = data_source["gas", "density"].min()
   c_max = data_source["gas", "density"].max()
   master_clump.add_validator("min_cells", 8)
   find_clumps(master_clump, c_min, c_max, step)
   master_clump.add_info_item("cell_mass")
   leaf_clumps = master_clump.leaves
   for child in leaf_clumps:
      child_mass=child.quantities.total_mass()
      mass.append(child_mass.v[0])
   numpy.savez(stest+"_clump_mass",mass=mass)
 return


@derived_field(name="Halpha_Emissivity",units = "erg/s/cm**3")
def Halpha_Emissivity(field,data):
   hden_n_bins,hden_min, hden_max = 11,-1,4
   T_n_bins, T_min, T_max = 51,3,8
   patt=homedir+"/code/cloudy_cooling_tools/h_emissivity1/h_emissivity_run%i.dat"
#   patt=Britton+"/cloudy_cooling_tools/h_emissivity1/h_emissivity_run%i.dat"
#   hden_n_bins,hden_min, hden_max = 5,0,4
#   T_n_bins, T_min, T_max = 101,4,6
#   patt=homedir+"/yuan/Britton/cloudy_cooling_tools/h_emissivity4/h_emissivity_run%i.dat"
   from scipy import interpolate
   hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
   T=numpy.linspace(T_min,T_max, T_n_bins)
   table = numpy.zeros((hden_n_bins, T_n_bins))
   for i in range(hden_n_bins):
       table[i,:]=[float(l.split()[-1]) for l in open(patt%(i+1)) if l[0] != "#"]
   sp=interpolate.RectBivariateSpline(hden,T,table)
   good=data["temperature"].shape
   H_N=numpy.log10(numpy.array(data["gas","number_density"]))*3.0/7.0 #convert to H_NumberDensity
   Temperature=numpy.log10(numpy.array(data["temperature"]))
   H_N=H_N.reshape(H_N.size)
   Temperature=Temperature.reshape(Temperature.size)
   dia=sp.ev(H_N,Temperature)
   dia=dia.reshape(good)
   Halpha=(10.0**dia)*((data["gas","number_density"].v*3.0/7.0)**2.0)
   return YTArray(Halpha,"erg/s/cm**3")   #*1.87e-12   # not yet
#add_field("Halpha_Emissivity",units=r"\rm{ergs s^{-1}cm^{-3}",function=_Halpha_Emissivity)

@derived_field(name="Halpha_Luminosity",units = "erg/s")
def Halpha_Luminosity(field,data):
   return data["Halpha_Emissivity"]*data["index", "cell_volume"].in_cgs()

def Halpha_Fits(imin=0,imax=None,res=1000):
 if imax is None:
    imax=imin+1
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   ds = yt.load(fn)
   proj = ds.proj("Halpha_Emissivity", axis='x',weight_field=None)
   proj2 = ds.proj(('enzo', 'x-velocity'), axis='x',weight_field="Halpha_Emissivity")
   frb = proj.to_frb((100.,"kpc"), res)
   frb2=proj2.to_frb((100.,"kpc"), res)
   prj_fits = yt.FITSImageData(frb, fields=["Halpha_Emissivity"],units="kpc")
   prj_fits2 = yt.FITSImageData(frb2, fields=[('enzo', 'x-velocity')],units="kpc")
   prj_fits3 = yt.FITSImageData.from_images([prj_fits, prj_fits2])
   prj_fits3.writeto(stest+"_Halpha.fits", clobber=True)
 return

