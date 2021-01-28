import matplotlib
matplotlib.use('Agg')
from yt.mods import *
from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"
import numpy as numpy
import matplotlib.pylab as pylab
import os as os
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
kboltz=1.38e-16

pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
pylab.rc('lines', linewidth=2)
font = {'weight' : 'bold', 'size'   : 18}
pylab.rc('font', **font)
#pylab.rc('text', usetex=False)
pylab.rc('text', usetex=True)
pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
def slice(imin=None,imax=None,width=None, upload=False,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if width is None:
    width = 20
 for i in range(imin,imax,mod):
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_slice("Density", 0).modify["velocity"]()
#   pc.add_slice("Density", 1).modify["velocity"]()
   pc.add_slice("ZeusHeatingRate", 0)
   pc.add_slice("Temperature", 0)
   pc.add_slice("Pressure", 0)
   pc.add_slice("Entropy", 0)
   pc.add_slice("Entropy", 2)
   pc.add_slice("GridLevel", 0)
   pc.add_slice("CoolingTime",0)
   pc.add_slice("CoolingTime",1)
   pc.add_slice("CoolingTime",2)
   pc.set_width(width, 'kpc')
   pc.save()
   if upload==True:
      pc.hub_upload()
 return

def slice_paper(i_list=[225,255],width=40):
 for i in i_list:
   fn = name_pattern %(i,i)
   stest = 'stest_%04i' %i
   pf = load(fn)
   field_list=["Density","Temperature","Pressure","GridLevel","Entropy"]
   for name in field_list:
     slc = SlicePlot(pf,0,name, center='c', width=(width,'kpc'))
     slc.set_font({'size':24})
     slc.save(stest+name+"slice.png")
     pylab.clf()
 return

Britton="/code/Britton/"
def _Halpha_Emissivity(field,data):
   hden_n_bins,hden_min, hden_max = 11,-1,4
   T_n_bins, T_min, T_max = 51,3,8
#   patt=homedir+"/source/cloudy_cooling_tools/h_emissivity1/h_emissivity_run%i.dat"
   patt=homedir+Britton+"/cloudy_cooling_tools/h_emissivity1/h_emissivity_run%i.dat"
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
   good=data["Temperature"].shape
   H_N=numpy.log10(numpy.array(data["H_NumberDensity"]))
   Temperature=numpy.log10(numpy.array(data["Temperature"]))
   H_N=H_N.reshape(H_N.size)
   Temperature=Temperature.reshape(Temperature.size)
   dia=sp.ev(H_N,Temperature)
   dia=dia.reshape(good)
   Halpha=(10.0**dia)*(data["H_NumberDensity"]**2.0)
   return Halpha   #*1.87e-12   # not yet
add_field("Halpha_Emissivity",units=r"\rm{ergs s^{-1}cm^{-3}",function=_Halpha_Emissivity)
#add_field("Halpha_Emissivity",units=r"\rm{ergs s^{-1}cm^{-3}arcsec^{-2}",function=_Halpha_Emissivity)

def _Halpha_Luminosity(field,data):
   return data["Halpha_Emissivity"]*data["CellVolume"]
add_field("Halpha_Luminosity", units=r"\rm{ergs s^{-1}", function=_Halpha_Luminosity)


def emission_projection(imin=None,imax=None,width=100,Halpha=True, Xray=True, upload=False):    #Sam
 import yt.analysis_modules.spectral_integrator.spectral_frequency_integrator as SI
 SI.add_xray_emissivity_field(0.5,3.0,with_metals=False,constant_metallicity=0.5)
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   if Halpha:
      pc.add_projection("Halpha_Emissivity", 0)
   if Xray:
      pc.add_projection("Xray_Emissivity_0.5_3.0keV", 0)
   pc.set_width(width, 'kpc')
   pc.save()
   if upload==True:
      pc.hub_upload()
 return

def off_projection(imin=None,imax=None,width=40,thicknesskpc=4000,upload=False,n_pics=10,theta=0,Halpha=True, Xray=True, Hardness=True, off_Velocity=False, Hard=False, Soft=False):
 import yt.analysis_modules.spectral_integrator.spectral_frequency_integrator as SI
 SI.add_xray_emissivity_field(0.5,3.0,with_metals=False,constant_metallicity=0.5)
 SI.add_xray_emissivity_field(0.3,1.5,with_metals=False,constant_metallicity=0.5)
 SI.add_xray_emissivity_field(1.5,7.5,with_metals=False,constant_metallicity=0.5)
 def _Hardness_Ratio(field,data):
   return data["Xray_Emissivity_1.5_7.5keV"]/data["Xray_Emissivity_0.3_1.5keV"]
 add_field("Hardness_Ratio",function=_Hardness_Ratio,take_log=False)
 def _off_Velocity(field,data):
   return (data["y-velocity"]*numpy.sin(theta)+data["z-velocity"]*numpy.cos(theta))/1.0e5   #km/s
 add_field("off_Velocity",units=r"km/s",function=_off_Velocity,take_log=False)
 pylab.rc('text', fontsize=16)
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   pf = load(fn)
   center=[0.5,0.5,0.5]
   thickness = thicknesskpc/pf['kpc']
   width=width/pf['kpc']
   wu=(width,width,thickness)
   print wu
   N = 512   #resolution
   dtheta = numpy.pi*0.5/n_pics
   for j in range(n_pics):
     L = [0,numpy.sin(theta), numpy.cos(theta)]
     if Xray:
        image = off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_0.5_3.0keV")
        #write_image(na.log10(image), "%s_offaxis_projection%.1f_Xray.png" % (pf, theta))
        write_projection(image, "%s_offaxis_projection%.1f_Xray.png" % (pf, theta), colorbar_label="X ray Surface Brightness "+r"($ergs \, s^{-1} \, cm^{-2}$)")
#        write_projection(image, "%s_offaxis_projection%.1f_Xray.png" % (pf, theta), colorbar_label=r"X-ray "+r"($ergs \, s^{-1} \, cm^{-2}$)")
     if Halpha:
        image2=off_axis_projection(pf, center, L, (width,width,thickness), N, "Halpha_Emissivity")
#        image2=numpy.rot90(image2)
        write_projection(image2*1.87e-12, "%s_offaxis_projection%.1f_Halpha.png" % (pf, theta), colorbar_label=r"$H_{\alpha} \,(ergs \, s^{-1} \, cm^{-2} \, arcsec^{-2})$", limits=[1e-18,5e-15],cmap_name='hot')
        #p = OffAxisProjectionPlot(pf, L, "Halpha_Emissivity", center='c', width=width, weight_field=None).set_zlim("Halpha_Emissivity",1e-16,1e-11)   ##set_zlim does not work
        #p.save()
     if Hardness:
        image3=off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_1.5_7.5keV")
        image4=off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_0.3_1.5keV")
        image5=image3/image4
        write_projection(image5, "%s_offaxis_projection%.1f_Hardness.png" % (pf, theta),take_log=False, colorbar_label=r"Hardness Ratio")
     if off_Velocity:
        image6=off_axis_projection(pf, center, L, (width,width,thickness), N, "off_Velocity", "Halpha_Emissivity")
        write_projection(image6, "%s_offaxis_projection%.1f_off_Velocity.png" % (pf, theta),take_log=False, colorbar_label=r"Line-of-sight Velocity " + r"($km/s$)")  
     if Hard:
        image3=off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_1.5_7.5keV")
        image3=numpy.rot90(numpy.array(image3))
        write_projection(image3, "%s_offaxis_projection%.1f_HardX.png" % (pf, theta),colorbar_label=r"Hard X-ray",limits=[1e-3,0.2])
     if Soft:
        image4=off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_0.3_1.5keV")
        image4=numpy.rot90(image4)
        write_projection(image4, "%s_offaxis_projection%.1f_SoftX.png" % (pf, theta),colorbar_label=r"Soft X-ray",limits=[1e-3,0.1])

     theta += dtheta
 return


#def L_Halpha(list=["TP_188/DD0188/stest_0188","TP_188/DD0400/stest_0400","Stampede/DD0255/stest_0255"],rmax=40):
# for name in list:
#    pf=load(homedir+"/big/"+name) 
#    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
#    L_Halpha=sp.quantities["TotalQuantity"](["Halpha_Luminosity"])[0]
#    print name, L_Halpha
# return

def thin_projection(imin=None,imax=None,width=None, thicknesskpc=None, upload=False, px=True, py=False, pz=False,mod=1, stars=True):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if width is None:
    width = 80   ##kpc
 if thicknesskpc is None:
    thicknesskpc = 80  ##kpc
 from yt.analysis_modules.spectral_integrator.api import add_xray_emissivity_field
 add_xray_emissivity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
 cmap_UV=matplotlib.cm.Blues_r
 cmap_UV.set_bad('k')
 cmap_star=matplotlib.cm.Greys
 cmap_star.set_bad('w')
 cmap_X=matplotlib.cm.Purples_r
 for i in range(imin,imax,mod):
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   thickness = thicknesskpc/pf['kpc']
   if px == True:
#     pc.add_thin_projection("Density", 0, thickness).set_zlim(8e-3,0.5)
#     pc.add_thin_projection("Temperature", 0, thickness, weight_field="Density").set_zlim(1e6,7e7)
#     pc.add_thin_projection("Entropy", 0, thickness, weight_field="Density").set_zlim(1e-9,3e-8)
#     pc.add_thin_projection("Pressure", 0, thickness, weight_field="Density").set_zlim(1e-10,1e-8)
     pc.add_thin_projection("ZeusHeatingRate", 0, thickness)
#     p=pc.add_thin_projection("Xray_Emissivity_0.5_9.9keV", 0, thickness, weight_field=None)
#     p.set_cmap(cmap_X)
#     p.set_zlim(9.0e-4,6.0e-2)
     if stars == True:
        try:
          pp=pc.add_thin_projection("star_density", 0, thickness, weight_field=None).set_cmap(cmap_star)
          pp.set_zlim(5.0e-5,1.0e-2)
          q=pc.add_thin_projection("young_star_density", 0, thickness, weight_field=None).set_cmap(cmap_UV)
          q.set_zlim(5.0e-6,1.0e-2)
#          pc.add_thin_projection("star_density", 0, thickness, weight_field=None)
          #pc.add_thin_projection("young_star_density", 0, thickness, weight_field=None).set_cmap(cmap_UV)
          #pc.add_thin_projection("young_star_density", 0, thickness, weight_field="star_creation_time")
#          pc.add_thin_projection("StarAgeYears", 0, thickness, weight_field="young_star_density").set_zlim(1e6,5e7)
        except:
          print "no stars formed yet in ", pf
   if py == True:
     pc.add_thin_projection("Density", 1, thickness).set_zlim(8e-3,0.5)
     pc.add_thin_projection("Temperature", 1, thickness, weight_field="Density").set_zlim(1e6,7e7)
     pc.add_thin_projection("Entropy", 1, thickness, weight_field="Density").set_zlim(1e-9,3e-8)
     pc.add_thin_projection("Pressure", 1, thickness, weight_field="Density").set_zlim(1e-10,1e-8)
     if stars == True:
        try:
          pc.add_thin_projection("star_density", 1, thickness, weight_field=None)
          pc.add_thin_projection("young_star_density", 1, thickness, weight_field=None).set_cmap(cmap_UV)
          #pc.add_thin_projection("young_star_density", 1, thickness, weight_field=None)
        except:
          print "no stars formed yet in ", pf
   if pz == True:
     pc.add_thin_projection("Density", 2, thickness).set_zlim(8e-3,0.5)
     pc.add_thin_projection("Temperature", 2, thickness, weight_field="Density").set_zlim(1e6,7e7)
     pc.add_thin_projection("Entropy", 2, thickness, weight_field="Density").set_zlim(1e-9,3e-8)
     pc.add_thin_projection("Pressure", 2, thickness, weight_field="Density").set_zlim(1e-10,1e-8)
     if stars == True:
        try:
          pc.add_thin_projection("star_density", 2, thickness, weight_field=None)
          pc.add_thin_projection("young_star_density", 2, thickness, weight_field=None)
        except:
          print "no stars formed yet in ", pf
   pc.set_width(width, 'kpc')
   pc.save()
   if upload==True:
      pc.hub_upload()
 return

def thin_projection_paper(rmax=20, i_list=[225,255], imin=None, imax=None,mod=1):
# for i in i_list:
 for i in range(imin,imax,mod): 
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   t=pf.current_time
   t*=(TimeUnits/3.16e16)   ## code unit to Gyr
   center= numpy.array([0.5, 0.5, 0.5])
   LE=center-rmax/pf['kpc']
   RE=center+rmax/pf['kpc']
   reg = pf.h.region(center,LE,RE)
   proj = pf.h.proj(0, 'Density', source=reg)
   plot = proj.to_pw(width=(rmax*2,'kpc'))
   plot.set_font({'size':24})
   plot.set_zlim("Density",8e-3,0.5)
   plot.annotate_text((0.75, 0.92), "t = %.2f Gyr"%t, text_args={'size':22, 'color':'w'})
   plot.save(stest+"Density_proj.png")
   pylab.clf()
   field_list=["Temperature","Entropy","Pressure"]
   z_list=[(8e5,9e7),(9e-10,2e-7),(8e-11,1e-8)]
   for name,zz in zip(field_list,z_list):
      proj = pf.h.proj(0, name, weight_field="Density", source=reg)
      plot = proj.to_pw(width=(rmax*2,'kpc'))
      plot.set_font({'size':24})
      plot.set_zlim(name,zz[0],zz[1])
      plot.annotate_text((0.75, 0.92), "t = %.2f Gyr"%t, text_args={'size':22, 'color':'w'})
      plot.save(stest+name+"_proj.png")
      pylab.clf()
 return


def thin_projection_particles(imin=None,imax=None,width=None, thicknesskpc=None, upload=False):
 import h5py
 if imin is None:
    imin=300
 if imax is None:
    imax=345
 if width is None:
    width = 20   ##kpc
 if thicknesskpc is None:
    thicknesskpc = 20  ##kpc
 for i in range(imin,imax):
   j=i-imin
   f = h5py.File("TracerParticles_%08i.h5" % j, "r")
#   p_pos=[f["position_x"],f["position_y"],f["position_z"]]
   xpos,ypos,zpos=f["position_x"],f["position_y"],f["position_z"]
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   thickness = thicknesskpc/pf['kpc']
   p=pc.add_thin_projection("Density", 0, thickness).modify["particles"](width=0.05, p_size=1.0, col='k', marker='o',ptype=3)
#   p.modify["marker"](pos=[xpos[j],ypos[j],zpos[j]],marker='x')
#   p.modify["particles"](ptype=3)
   pc.set_width(width, 'kpc')
   pc.save()
 return

def particles(imin=None,nfiles=10, ncheck=5, nTP=None, case=1, deadcat=40, hot=1, width=20, ColdTemperature=1.0e6, HotTemperature=5e6):
 import h5py
 pylab.rc('text', fontsize=16)
 Rs = 8.31e7
 lu = 4.93708012753e+25   # length unit
 xpos = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 ypos = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 zpos = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 T = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 rho = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 Vx = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 Vy = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 Vz = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 P = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 radius = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 Vr = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 Vmag = numpy.zeros(shape=(nTP,nfiles),dtype='float32')
 time = numpy.zeros(nfiles)
 for j in range(nfiles):
   f = h5py.File("TracerParticles_%08i.h5" % j, "r")
   time[j] = f["/"].attrs["time"]
   xpos[:,j],ypos[:,j],zpos[:,j]=f["position_x"],f["position_y"],f["position_z"]
   Vx[:,j],Vy[:,j],Vz[:,j]=f["velocity_x"],f["velocity_y"],f["velocity_z"]
   radius[:,j]=((xpos[:,j]-0.5)**2.0+(ypos[:,j]-0.5)**2.0+(zpos[:,j]-0.5)**2.0)**0.5 ##in codeunit
   Vr[:,j]=(((xpos[:,j]-0.5)*Vx[:,j]+(ypos[:,j]-0.5)*Vy[:,j]+(zpos[:,j]-0.5)*Vz[:,j])/radius[:,j])/1.0e5  ## from cm/s to km/s
######   Vr[:,j]=f["RadialVelocity"]/1.0e5  does not work
   radius[:,j]*= (lu/kpc) ## now in kpc 
   T[:,j]=f["temperature"]
   rho[:,j]=f["density"]
   P[:, j] = rho[:, j]*Rs*T[:,j]
##   Vmag=f["VelocityMagnitude"]/1.0e5   does not work
   Vmag[:,j]=(Vx[:,j]**2.0+Vy[:,j]**2.0+Vz[:,j]**2.0)**0.5/1.0e5
 if hot == 0:    #hot before deadcat
   cold=(T[:,ncheck-1] < ColdTemperature)
   if (deadcat>0):   #we do not want the particles that have always been cold
      for cat in range(deadcat):
        cold = cold & (T[:,cat] > HotTemperature)     ### if we require the particles to be "hot" at first
        cold = cold & (radius[:,cat] > 0.5)   # initially larger than 0.5 kpc from the center
 if hot == 1:    # hot at deadcat 
   cold=(T[:,ncheck-1] < ColdTemperature)
   if (deadcat>0):   #we do not want the particles that have always been cold
      cold = cold & (T[:,nfiles-deadcat] > HotTemperature)     ### only the hot->cold particles
      cold = cold & (rho[:,nfiles-deadcat] < 1e-24)   ##and not dense
      cold = cold & (radius[:,nfiles-deadcat] > 1)   # larger than 1 kpc from the center
 if hot == 2:   ##select the particles that were cold but then became hot
   cold=(T[:,ncheck-1] > HotTemperature)
   if (deadcat>0):   #we do not want the particles that have always been cold
      cold = cold & (T[:,0] < ColdTemperature)     ### only the hot->cold particles
#      cold = cold & (rho[:,nfiles-deadcat] < 1e-24)   ##and not dense
 time=np.arange(nfiles)/30.0   ##1/30 Myr is the timestep; redefine time here -- or not
 radius_cold=radius[cold,:]
 xpos_cold=xpos[cold,:]
 ypos_cold=ypos[cold,:]
 zpos_cold=zpos[cold,:]
 T_cold=T[cold,:]
 rho_cold=rho[cold,:]
 Vr_cold=Vr[cold,:]
 Vx_cold=Vx[cold,:]
 Vy_cold=Vy[cold,:]
 Vz_cold=Vz[cold,:]
 Vmag_cold=Vmag[cold,:]
 P_cold=P[cold,:]
 n_cold=radius_cold[:,0].size # number of particles that end up being cold
 cool_step=numpy.zeros(n_cold, dtype=numpy.int) # the step when the particle is about to cool
 for k in range(n_cold):
   j = ncheck-1
   if hot == 2:
      while  ((T_cold[k,j]>HotTemperature) and (j >0)):
         j=j-1
   else:
#      while ((T_cold[k,j]<ColdTemperature) and (j >0)):    ###????which one????
      while ((T_cold[k,j]<HotTemperature) and (j >0)):
         j=j-1
   cool_step[k]=j   #the k particle cools at step j (cool_step[k])
 print "n_cold, cool_step=", n_cold, cool_step 
# print "T_cold last=", T_cold[:,nfiles-1]
 nplot=deadcat+nfiles-ncheck
 rho_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 Vr_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 T_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 r_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 xpos_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 ypos_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 zpos_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 Vmag_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 P_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 Vx_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 Vy_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
 Vz_plot=numpy.zeros(shape=(n_cold,nplot),dtype='float32')
# for i in range(deadcat):
 for i in range(nplot):
    for k in range(n_cold):
      #tt=numpy.maximum(0,cool_step[k]+1-i)   # +1 is the step where the gas has just cooled
      tt=numpy.maximum(0,cool_step[k]+1+nfiles-ncheck-i)   # +1 is the step where the gas has just cooled
      #print "tt", tt
      rho_plot[k,i]=rho_cold[k,tt]   #rho_plot[k,0] is the last step, and rho_plot[k,(nfiles-ncheck)] is the step when the gas cools
      Vr_plot[k,i]=Vr_cold[k,tt]
      T_plot[k,i]=T_cold[k,tt]
      r_plot[k,i]=radius_cold[k,tt]
      xpos_plot[k,i]=xpos_cold[k,tt]
      ypos_plot[k,i]=ypos_cold[k,tt]
      zpos_plot[k,i]=zpos_cold[k,tt]
      Vmag_plot[k,i]=Vmag_cold[k,tt]
      P_plot[k,i]=P_cold[k,tt]
      Vx_plot[k,i]=Vx_cold[k,tt]
      Vy_plot[k,i]=Vy_cold[k,tt]
      Vz_plot[k,i]=Vz_cold[k,tt]
    pylab.figure(i) 
    if case==1:
       pylab.plot(numpy.log10(rho_plot[:,i]),Vr_plot[:,i],'.')
       pylab.xlim(-26,-22)
       pylab.ylim(-1.0e3,1.0e3)
       pylab.savefig('rho_Vr_cool%i.png' %i)
    if case == 11:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho_plot[:,i]), Vr_plot[:,i], bins=100,range=([-25,-22],[-1.0e3,1.0e3]))
       #H, xedges, yedges = numpy.histogram2d(numpy.log10(rho_plot[:,i]), numpy.log10(abs(Vr_plot[:,i]))*numpy.sign(Vr_plot[:,i]), bins=100,range=([-26,-22],[-8,8]))    # don't you ever take the log of a negative number!!
       #extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]   # never trust the internet!
       #pylab.imshow(H, interpolation='nearest',vmax=135)
       pylab.imshow(H, interpolation='nearest')  ##,norm=???) ## how do I normalize it?
       pylab.xlabel(r'$v_r\, (km/s)$')
       pylab.ylabel('log '+r'$\rho \, (g/cm^3)$')
       y_locs=numpy.arange(0,100,10)
       y_labels=numpy.arange(-25,-22,3.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,100,20)
       x_labels=numpy.arange(-1e3,1e3,4.0e2)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.tight_layout()
       pylab.savefig('rho_Vr_cool_hist%i.png' %i)
    if case == 2:
       pylab.plot(numpy.log10(rho_plot[:,i]),numpy.log10(T_plot[:,i]),'.')
       pylab.xlim(-26,-22)
       pylab.ylim(4,8)
       pylab.savefig('rho_T_cool%i.png' %i)
    if case == 22:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho_plot[:,i]), numpy.log10(P_plot[:,i]), bins=300,range=([-26,-22],[-11,-8]))
       pylab.imshow(H, interpolation='nearest')
       pylab.xlabel('P')
       pylab.ylabel('rho')
       y_locs=numpy.arange(0,300,30)
       y_labels=numpy.arange(-26,-22,4.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,300,30)
       x_labels=numpy.arange(-11,-8,3.0/10)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.savefig('rho_P_cool_hist%i.png' %i)
    if case== 33:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho_plot[:,i]), r_plot[:,i], bins=300,range=([-26,-22],[0,width/2.0]))
       pylab.imshow(H, interpolation='nearest')
       pylab.xlabel('r')
       pylab.ylabel('rho')
       y_locs=numpy.arange(0,300,30)
       y_labels=numpy.arange(-26,-22,4.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,300,30)
       x_labels=numpy.arange(0,width/2.0,width/20.0)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.savefig('rho_r_cool_hist%i.png' %i)
    if case== 44:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho_plot[:,i]), xpos_plot[:,i], bins=300,range=([-26,-22],[0,width/2.0]))
       pylab.imshow(H, interpolation='nearest')
       pylab.xlabel('x')
       pylab.ylabel('rho')
       y_locs=numpy.arange(0,300,30)
       y_labels=numpy.arange(-26,-22,4.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,300,30)
       x_labels=numpy.arange(0,width/2.0,width/20.0)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.savefig('rho_xpos_cool_hist%i.png' %i)       
    if case== 55:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho_plot[:,i]), zpos_plot[:,i], bins=300,range=([-26,-22],[0,width/2.0]))
       pylab.imshow(H, interpolation='nearest')
       pylab.xlabel('z')
       pylab.ylabel('rho')
       y_locs=numpy.arange(0,300,30)
       y_labels=numpy.arange(-26,-22,4.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,300,30)
       x_labels=numpy.arange(0,width/2.0,width/20.0)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.savefig('rho_zpos_cool_hist%i.png' %i)       
# pylab.xlabel('Time (Myr)')
# pylab.ylabel('Radius ()')
# for k in range(nTP):  # k is particle index
#   pylab.semilogy(time,radius[k,:])
# pylab.figure(1)
# for k in range(n_cold):
#   pylab.semilogy(time,radius_cold[k,:])
#   pylab.savefig('TP_cold_radius_time.png')
#   pylab.semilogy(time,T_cold[k,:], marker='x')
#   pylab.savefig('TP_cold_T_time.png')
#   pylab.semilogy(time,rho_cold[k,:])
#   pylab.savefig('TP_cold_rho_time.png')
 pylab.clf()
 numpy.savez("cool_particles", rho_plot=rho_plot, Vr_plot=Vr_plot, T_plot=T_plot, r_plot=r_plot, xpos_plot=xpos_plot, ypos_plot=ypos_plot, zpos_plot=zpos_plot, Vmag_plot=Vmag_plot, P_plot=P_plot, Vx_plot=Vx_plot, Vy_plot=Vy_plot, Vz_plot=Vz_plot)
 numpy.savez("all_particles", rho=rho, Vr=Vr, T=T, radius=radius, xpos=xpos, ypos=ypos, zpos=zpos, Vmag=Vmag, P=P, Vx=Vx, Vy=Vy, Vz=Vz)
 for i in range(nfiles):
    pylab.figure(i+1000)
    if case == 11:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho[:,i]), Vr[:,i], bins=100,range=([-25,-22],[-1.0e3,1.0e3]))
       #pylab.imshow(H, interpolation='nearest',vmax=135)
       pylab.imshow(H, interpolation='nearest')
       pylab.xlabel(r'$v_r \, (km/s)$')
       pylab.ylabel('log '+r'$\rho\, (g/cm^{3})$')
       y_locs=numpy.arange(0,100,10)
       y_labels=numpy.arange(-25,-22,3.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,100,20)
       x_labels=numpy.arange(-1e3,1e3,4.0e2)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.tight_layout()
       pylab.savefig('rho_Vr_cool_hist_all%i.png' %i)
    if case == 22:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho[:,i]), numpy.log10(P[:,i]), bins=300,range=([-26,-22],[-11,-8]))
       pylab.imshow(H, interpolation='nearest')
       pylab.xlabel('P')
       pylab.ylabel('rho')
       y_locs=numpy.arange(0,300,30)
       y_labels=numpy.arange(-26,-22,4.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,300,30)
       x_labels=numpy.arange(-11,-8,3.0/10)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.savefig('rho_P_cool_hist_all%i.png' %i)
    if case== 33:
       H, xedges, yedges = numpy.histogram2d(numpy.log10(rho[:,i]), radius[:,i], bins=300,range=([-26,-22],[0,width/2.0]))
       pylab.imshow(H, interpolation='nearest')
       pylab.xlabel('r')
       pylab.ylabel('rho')
       y_locs=numpy.arange(0,300,30)
       y_labels=numpy.arange(-26,-22,4.0/10)
       pylab.yticks(y_locs, y_labels)
       x_locs=numpy.arange(0,300,30)
       x_labels=numpy.arange(0,width/2.0,width/20.0)
       pylab.xticks(x_locs, x_labels)
       pylab.colorbar()
       pylab.savefig('rho_r_cool_hist_all%i.png' %i)
 pylab.clf()

 if case == 11:
    for i in numpy.arange(-3,3,1):  #range is totally random
       pylab.figure(i+2013) 
       pylab.hist(Vr_plot[:,nfiles-ncheck+i], bins=30,normed=True,hatch='/',histtype='step') #the nth step before cooling
       pylab.hist(Vr[:,nfiles-ncheck],bins=100,normed=True,hatch='/',histtype='step',color='green')  #all of them
       pylab.xlim([-400,1200])
       pylab.xlabel(r'$v_r\,(km/s)$')
       pylab.legend(("Cold","All"),loc="upper right")
       pylab.tight_layout()
       pylab.savefig('Vr_cool_1Dhist_all%i.png' %i)
 if case == 111:
    for i in numpy.arange(-3,3,1):  #range is totally random
       pylab.figure(i+2013)
       pylab.hist(r_plot[:,nfiles-ncheck+i], bins=30,normed=True,hatch='/',histtype='step') #the nth step before cooling
       pylab.hist(radius[:,nfiles-ncheck],bins=60,normed=True,hatch='/',histtype='step',color='green')  #all of them
       pylab.xlabel('Radius (kpc)')
       pylab.savefig('r_cool_1Dhist_all%i.png' %i)
 if case== 1111:
    for i in numpy.arange(-3,3,1):  #range is totally random
       pylab.figure(i+3000)
       pylab.hist(Vmag_plot[:,nfiles-ncheck+i], bins=30,normed=True,hatch='/',histtype='step') #the nth step before cooling
       pylab.hist(Vmag[:,nfiles-ncheck],bins=60,normed=True,hatch='/',histtype='step',color='green')  #all of them
#       pylab.ylim([0,1])
       pylab.xlabel('Vmag')
       pylab.savefig('Vmag_cool_1Dhist_all%i.png' %i)
 if case== 1984:
    t_kyr=numpy.arange(nfiles-ncheck)*100.0
    for i in numpy.arange(100):
       pylab.figure(1984+i)
       pylab.plot(t_kyr,Vr_plot[i,-(nfiles-ncheck):])
       pylab.xlabel("time (thousand yr)")
       pylab.ylabel("Vr")
       pylab.savefig("Vr_of_cold_TP%i.png" %i)
 pylab.clf()

 if case == 250:
    pylab.figure(0)
    for i in range(n_cold-1):
       pylab.plot(time,radius_cold[i,:])
    pylab.savefig("r_t_cold.png")
    pylab.figure(1)
    for i in range(nfiles-1):
      if (i%10==0):
         pylab.plot(time,radius[i,:])
    pylab.savefig("r_t_all.png")
 pylab.clf()
 return


def projection_plot(imin=None,imax=None,width=None, thicknesskpc=None, upload=False):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if width is None:
    width = 40   ##kpc
 if thicknesskpc is None:
    thicknesskpc = 40  ##kpc
 for i in range(imin,imax):
    fn = name_pattern %(i,i)
    pf = load(fn)
    p=ProjectionPlot(pf, "z", "Density", center=[0.5, 0.5, 0.5], width=(width,'kpc'))
    p.annotate_particles(width=0.05, p_size=5.0, col='k', marker='o')
    p.save()
 return

def v_projection_plot(imin=None,imax=None,width=None, thicknesskpc=None, upload=False):   #for Kevin's Alma proposal
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if width is None:
    width = 40   ##kpc
 if thicknesskpc is None:
    thicknesskpc = 40  ##kpc
 for i in range(imin,imax):
    fn = name_pattern %(i,i)
    pf = load(fn)
#    RE= numpy.array([0.5, 0.5, 0.5])
    rmax=width/4.0
#    center=RE-rmax/pf['kpc']
#    LE=center-rmax/pf['kpc']
    center = numpy.array([0.5, 0.5, 0.5]) - numpy.array([rmax/pf['kpc'],0,rmax/pf['kpc']])
    LE=center-rmax/pf['kpc']
    RE=center+rmax/pf['kpc']
    reg = pf.h.region(center,LE,RE)
    reg=reg.cut_region(["grid['Temperature']<1.0e4"])
    p=ProjectionPlot(pf, "x", ["x-velocity"], center=center, width=(2*rmax,'kpc'),weight_field="Halpha_Emissivity", data_source=reg)
    p.set_zlim("x-velocity",-7e7,7e7)
    #p.annotate_particles(width=0.05, p_size=5.0, col='k', marker='o')
    p.annotate_particles(width=0.05, p_size=20.0, col='k', marker='o')
    p.save()
 return



def rotate_projection(imin=0,imax=None,width=20, thicknesskpc=20, upload=False):
 if imax is None:
    imax=imin + 1
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   pf = load(fn)
   thickness = thicknesskpc/pf['kpc']
   c = [0.5, 0.5, 0.5]
   N = 512
   theta = 0
   n_pics = 30  #total number of projections
   dtheta = numpy.pi*2/n_pics
   for j in range(n_pics):
     theta += dtheta
     L = [numpy.cos(theta), numpy.sin(theta), 0]
     image = off_axis_projection(pf, c, L, thickness, N, "Density", no_ghost=False)
     write_image(na.log10(image), "%s_offaxis_projection%.1f.png" % (pf, theta))
 return

def move_camera(imin=0,imax=None,thicknesskpc=40, upload=False,phi=0,i_list=[66,177,188,281,483], N=512, n_round=40, theta=0):
 if imax is None:
   imax=imin + 1
 from yt.analysis_modules.spectral_integrator.api import add_xray_emissivity_field
 add_xray_emissivity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
# use_weight=["Xray_Emissivity_0.5_9.9keV","Xray_Emissivity_0.5_9.9keV","Xray_Emissivity_0.5_9.9keV","Xray_Emissivity_0.5_9.9keV",None,None,None]
 use_weight=["Density","Density","Density","Density",None,None,None]
 keys = ['Density', 'Temperature','Entropy','Pressure','star_density', 'young_star_density','Xray_Emissivity_0.5_9.9keV']
 values_titles = [r'$\mathrm{Density}\ (\mathrm{g\ cm^{-3}})$', r'$\mathrm{Temperature}\ (\mathrm{K})$', \
                 r'$\mathrm{Entropy}\ (\mathrm{ergs}\ \mathrm{cm}^{3\gamma-3})$', r'$\mathrm{Pressure}\ (\mathrm{dyne}/\rm{cm}^{2})$',\
                 r'$\mathrm{Stellar\ Density}\ (\mathrm{g\ cm^{-3}})$',r'$\mathrm{Young\ stellar\ density}\ (\mathrm{g\ cm^{-3}})$',
                 r'$\mathrm{L_{Xray}}\ \mathrm{ergs\ s}^{-1}$']
 values_z_list=[(7e-26,1e-23),(3e6,1e8),(9e-10,2e-7),(8e-11,1e-8),(5.0e-5,1.0e-2),(5.0e-6,1.0e-2),(9.0e-4,6.0e-2)]
 titles=dict(zip(keys, values_titles))
 z_list = dict(zip(keys, values_z_list))
 weight_list=dict(zip(keys, use_weight))
 cmap_UV=matplotlib.cm.Blues_r
 cmap_UV.set_bad('k')
 cmap_X=matplotlib.cm.Purples_r
 cmap_star=matplotlib.cm.Greys
 cmap_star.set_bad('w')
 cmap_defaul='algae'
 values_cmap=[cmap_defaul,cmap_defaul,cmap_defaul,cmap_defaul, cmap_star, cmap_UV, cmap_X]
 cmap_list = dict(zip(keys, values_cmap))
 n_pics=0
 d_theta = numpy.pi*2.0/n_round
 d_phi = numpy.pi*2.0/n_round
 center = [0.5, 0.5, 0.5]
 north_vector = [0,0,1]
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   pf = load(fn)
   thickness = thicknesskpc/pf['kpc']
   width=thickness
   time=pf.current_time
   time*=(TimeUnits/3.16e16)   ## code unit to Gyr
   look_field=keys[1]
   if i ==i_list[0]:
     look_field=keys[5]
     for j in range(n_round//2):       
       L = [numpy.cos(theta)*numpy.sin(phi), numpy.sin(theta)*numpy.sin(phi), numpy.cos(phi)]
       image = off_axis_projection(pf, center, L, (width,width,thickness), N, look_field, no_ghost=False,weight=weight_list[look_field],north_vector=north_vector)
       write_projection(image, "camera%04i.png" %n_pics, colorbar_label=titles[look_field], limits=z_list[look_field],cmap_name=cmap_list[look_field],title="t = %.2f Gyr"%time)
       pylab.clf()
       theta = theta + d_theta
       #phi = phi+ d_phi
       n_pics+=1
     look_field=keys[1]
   if i == i_list[2]:
     for look_field in [keys[5], keys[6]]:
       for j in range(n_round//3):
         L = [numpy.cos(theta)*numpy.sin(phi), numpy.sin(theta)*numpy.sin(phi), numpy.cos(phi)]
         image = off_axis_projection(pf, center, L, (width,width,thickness), N, look_field, no_ghost=False,weight=weight_list[look_field],north_vector=north_vector)
         write_projection(image, "camera%04i.png" %n_pics, colorbar_label=titles[look_field], limits=z_list[look_field],cmap_name=cmap_list[look_field],title="t = %.2f Gyr"%time)
         pylab.clf()
         theta = theta + d_theta
         n_pics+=1 
     look_field=keys[1]
   L = [numpy.cos(theta)*numpy.sin(phi), numpy.sin(theta)*numpy.sin(phi), numpy.cos(phi)]
#   L=[1,1,1]
   image = off_axis_projection(pf, center, L, (width,width,thickness), N, look_field, no_ghost=False,weight=weight_list[look_field],north_vector=north_vector)
   write_projection(image, "camera%04i.png" %n_pics, colorbar_label=titles[look_field], limits=z_list[look_field],cmap_name=cmap_list[look_field],title="t = %.2f Gyr"%time)
   pylab.clf()
   n_pics+=1
 return


def _Cell_ThermalEnergy(field,data):
   return data["ThermalEnergy"]*data["CellMass"]
add_field("Cell_ThermalEnergy", function=_Cell_ThermalEnergy)
def _Cell_KineticEnergy(field,data):
   return data["KineticEnergy"]*data["CellVolume"]
add_field("Cell_KineticEnergy", function=_Cell_KineticEnergy)


def profile_energy(imin=None,imax=None,rmax=100, mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 for i in range(imin,imax, mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   dog=BinnedProfile1D(sp, 50, "Radiuskpc", 0.1,rmax,lazy_reader=True,log_space=False)
   r=dog["Radiuskpc"]
   dog.add_fields("Cell_ThermalEnergy",weight=None)
   dog.add_fields("Cell_KineticEnergy",weight=None)
   E_th_total=dog["Cell_ThermalEnergy"]
   E_k_total=dog["Cell_KineticEnergy"]
   if i==imin:
      E_th0=E_th_total
      E_k0=E_k_total
   E_th=E_th_total-E_th0
   E_k=E_k_total-E_k0
   pylab.figure(i)
   pylab.semilogx(r,E_th, color='red')
   pylab.semilogx(r,E_k, color='blue')
   pylab.xlabel('r (kpc)')
   pylab.ylabel('Energy')
   pylab.xlim([0.1,rmax])
   pylab.savefig(stest+'_energy.png')
 return

def profile_one(i=0,rmin=1,rmax=300,xbins=128):
 from yt.analysis_modules.spectral_integrator.api import add_xray_emissivity_field
 add_xray_emissivity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
 fn = name_pattern %(i,i)
 stest = 'stest_%04i' %i
 fn = name_pattern %(i,i)
 pf = load(fn)
 pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
 pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Temperature"],weight="Xray_Emissivity_0.5_9.9keV",x_bins=xbins)
 r=pc.plots[-1].data["Radiuskpc"]
 Temperature=pc.plots[-1].data["Temperature"]
 pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Density"],weight="Xray_Emissivity_0.5_9.9keV",x_bins=xbins)
 Density=pc.plots[-1].data["Density"]
 pylab.loglog(r,Temperature)
 pylab.xlim([rmin,rmax])
 pylab.ylim([4e6,2e8])
 pylab.xlabel('r (kpc)')
 pylab.ylabel('T (K)')
 pylab.tight_layout()
 pylab.savefig(stest+'_Temperature_profile.png')
 pylab.clf()
 pylab.loglog(r,Density)
 pylab.xlabel('r (kpc)')
 pylab.ylabel('Density '+r'$(\rm{g}\ /\ \rm{cm}^3)$')
 pylab.xlim([rmin,rmax])
 pylab.ylim([1e-27,3e-23])
 pylab.tight_layout()
 pylab.savefig(stest+'_Density_profile.png')
 numpy.savez(stest+"profile",r=r,Density=Density,Temperature=Temperature)



def profile_more(imin=0,imax=None,rmin=1, rmax=1000.0, mod=1,CoolingTimeOn=1,xbins=128):
 from yt.analysis_modules.spectral_integrator.api import add_xray_emissivity_field
 add_xray_emissivity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
 if imax is None:
    imax=imin+1
 #cm=pylab.cm.Blues_r
 #cm=pylab.cm.autumn
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 #cm=pylab.cm.gist_heat
 cm=pylab.cm.gist_rainbow
 conc=6.81
 M_vir=8.5e14
 Rs=0.0224   #SphereCoreRadius=0.0224
 loop=range(imin,imax,mod)
 loop.insert(-1,0)
 for i in loop:
    fn = name_pattern %(i,i)
    pf = load(fn)
    pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
    pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Temperature"],weight="Xray_Emissivity_0.5_9.9keV",x_bins=xbins)
    r=pc.plots[-1].data["Radiuskpc"]
    Temperature=pc.plots[-1].data["Temperature"]
    pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Pressure"],weight="Xray_Emissivity_0.5_9.9keV",x_bins=xbins)
    Pressure=pc.plots[-1].data["Pressure"]*6.242e8   # from ergs to keV
    pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Density"],weight="Xray_Emissivity_0.5_9.9keV",x_bins=xbins)
    Density=pc.plots[-1].data["Density"]
    if CoolingTimeOn==1:
      pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Cooling_Time"],weight="Xray_Emissivity_0.5_9.9keV",x_bins=xbins)
      CoolingTime=pc.plots[-1].data["Cooling_Time"]/3.15e7
    pylab.figure(0,figsize=(8, 8))
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
       pylab.loglog(r,Temperature,color="k",linewidth=3.0)
    else:
       pylab.loglog(r,Temperature,color=cm(float(i-imin)/(imax-imin)))
    pylab.xlim([rmin,rmax])
    pylab.ylim([4e6,2e8])
    pylab.xlabel('r (kpc)')
    pylab.ylabel('T (K)')
    pylab.tight_layout()
    pylab.savefig('Temperature_rainbow.png')
    if CoolingTimeOn==1:
      x1=r/(Rs*pf['kpc'])
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


def _HotTemperature(field,data):
   HotTemperature=na.array(data["Temperature"])
   cold=HotTemperature < 1.0e5
   HotTemperature[cold]=0.0
   return HotTemperature
add_field("HotTemperature",function=_HotTemperature,take_log=False)

def _HotDensity(field,data):
   Temperature=na.array(data["Temperature"])
   cold=Temperature < 1.0e5
   HotDensity=na.array(data["Density"])
   HotDensity[cold]=0.0
   return HotDensity
add_field("HotDensity",function=_HotDensity,take_log=False)


def Xray_profile_more(imin=None,imax=None,rmax=None,mod=1):  #using add_profile_object
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=1000.0
 cm=pylab.cm.gist_rainbow
 for i in range(imin,imax,mod):
    fn = name_pattern %(i,i)
    pf = load(fn)
    center = [0.5, 0.5, 0.5]
    sp = pf.h.sphere(center, rmax/pf['kpc'])
    pc = PlotCollection(sp, center)
    p = pc.add_profile_object(sp, ["Radiuskpc", "HotDensity"])
    Radiuskpc = pc.plots[-1].data["Radiuskpc"]
    HotDensity= pc.plots[-1].data["HotDensity"]
    p = pc.add_profile_object(sp, ["Radiuskpc", "HotTemperature"])
    HotTemperature=pc.plots[-1].data["HotTemperature"]/keV
    pylab.figure(0)
    pylab.loglog(Radiuskpc,HotDensity, color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('r (kpc)')
    pylab.ylabel('Density (g/cm^3)')
    pylab.xlim([0.05,rmax])
    pylab.savefig('HotDensity_rainbow.png')
    pylab.figure(1)
    pylab.loglog(Radiuskpc,HotTemperature, color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('r (kpc)')
    pylab.ylabel('Temperature (keV)')
    pylab.xlim([0.05,rmax])
    pylab.savefig('HotTemperature_rainbow.png')
 pylab.clf()
 return



def Vr(imin=None,imax=None,rmax=None):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=100.0
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","RadialVelocity"],weight="CellMassMsun")
   r=pc.plots[-1].data["Radiuskpc"]
   Vr=pc.plots[-1].data["RadialVelocity"]/1.0e5
   pylab.figure(i)
   pylab.semilogx(r,Vr)
   pylab.xlabel('r (kpc)')
   pylab.ylabel('Vr (km/s)')
   pylab.savefig(stest+'_Vr.png')
 return



def Vrotation(imin=None,imax=None,rmax=None):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=100.0
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i) 
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","SpecificAngularMomentumZ"])
   r=pc.plots[-1].data["Radiuskpc"]
   Lz=pc.plots[-1].data["SpecificAngularMomentumZ"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","SpecificAngularMomentumX"])
   Lx=pc.plots[-1].data["SpecificAngularMomentumX"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","SpecificAngularMomentumY"])
   Ly=pc.plots[-1].data["SpecificAngularMomentumY"]
   Vrotation=pow(pow(Lx,2.0)+pow(Ly,2.0)+pow(Lz,2.0),0.5)/(r*kpc)/1.0e5
   Vrotation_z=abs(Lz)/(r*kpc)/1.0e5
   pylab.figure(i)
   pylab.semilogx(r,Vrotation_z)
   pylab.xlim([0.1,rmax])
   pylab.xlabel('r (kpc)')
   pylab.ylabel('Rotational velocity along z (km/s)')
   pylab.savefig(stest+'_RotationZ.png')
   pylab.semilogx(r,Vrotation)
   pylab.xlabel('r (kpc)')
   pylab.ylabel('Rotational velocity (km/s)')
   pylab.savefig(stest+'_Rotation.png')
 return


def Mdot(imin=None,imax=None,rmax=None, mod=1):
# from AddField import *
 if imin is None:
    imin=0
 if imax is None:
    max=imin+1
 if rmax is None:
    rmax=100.0
 M_sun=2.0e33
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","RadialMassFlux"])
   r=pc.plots[-1].data["Radiuskpc"]
   Flux=pc.plots[-1].data["RadialMassFlux"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","NegativeRadialMassFlux"])
   NeFlux=pc.plots[-1].data["NegativeRadialMassFlux"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","PositiveRadialMassFlux"])
   PoFlux=pc.plots[-1].data["PositiveRadialMassFlux"]
   dMdt=(Flux/M_sun)*4.0*3.14*pow(r*kpc,2.0)*3.15e7  ## from cgs to M_sun/yr
   NedMdt=(NeFlux/M_sun)*4.0*3.14*pow(r*kpc,2.0)*3.15e7  ## from cgs to M_sun/yr
   PodMdt=(PoFlux/M_sun)*4.0*3.14*pow(r*kpc,2.0)*3.15e7  ## from cgs to M_sun/yr
   pylab.figure(i)
   pylab.semilogx(r,dMdt)
   pylab.semilogx(r,NedMdt)
   pylab.semilogx(r,PodMdt)
   pylab.xlim([0.1,rmax])
   pylab.xlabel('r (kpc)')
   pylab.ylabel('dMdt (M_sun/yr)')
   pylab.legend(("Total","Inflow", "Outflow"),loc="upper right")
   pylab.savefig(stest+'_Mdot.png')
 return

def M_r(imin=None,imax=None,rmax=10.0, coldtemperature=None, binsize=1):
 import numpy as numpy
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if coldtemperature is None:
    coldtemperature = [1e5,1e6,1e7]
 time=[]
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   current_time = pf.current_time
   time.append(current_time)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   Temperature = sp["Temperature"]
   CellMassMsun = sp["CellMassMsun"]
   Radiuskpc = sp["Radiuskpc"]
   pylab.figure(i)
   print coldtemperature
   for coldT in coldtemperature:
     r=[]
     coldgas=[]
     print coldT
     cold = Temperature < coldT
     for j in range(int(numpy.floor(rmax/binsize))):
       rgood= (Radiuskpc > float(j)*binsize) & (Radiuskpc < float(1+j)*binsize)
       r.append(float(j)*binsize)
       mgood=rgood & cold   
       totalcoldgas=CellMassMsun[mgood]
       coldgas.append(sum(totalcoldgas))
     pylab.plot(r,numpy.array(coldgas)/binsize,drawstyle='steps-pre')
   pylab.xlabel('r (kpc)')
   pylab.ylabel('coldgasmass (M_sun)')
   pylab.xlim([0,rmax-1])
#   pylab.ylim([0,3.5e8])
   pylab.savefig(stest+'M_r_in%ikpc.png' %rmax)
 return

def M_r_more_new(imin=None,imax=None,rmax=None, coldtemperature=1e5,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=10.0   ##kpc
 time=[]
# cm=pylab.cm.gist_rainbow
 cm=pylab.cm.gist_heat
 pylab.figure(0)
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = load(fn)
    current_time = pf.current_time
    time.append(current_time)
    center = [0.5, 0.5, 0.5]
    sp = pf.h.sphere(center, rmax/pf['kpc'])
    pc = PlotCollection(sp, center)
#    p = pc.add_profile_object(sp, ["Radiuskpc", "ColdCellMassMsun"], weight = None,accumulation=True)
    p = pc.add_profile_object(sp, ["Radiuskpc", "ColdCellMassMsun"], accumulation=True,x_bins=rmax)   #binsize is 1 kpc
    ColdCellMassMsun = pc.plots[-1].data["ColdCellMassMsun"]
    bad=ColdCellMassMsun<0.1
    ColdCellMassMsun[bad]=0.1
    Radiuskpc = pc.plots[-1].data["Radiuskpc"]
    pylab.semilogy(Radiuskpc,ColdCellMassMsun, color=cm(float(i-imin)/(imax-imin)))
    #pylab.plot(Radiuskpc,ColdCellMassMsun, color=cm(float(i-imin)/(imax-imin)))
 pylab.xlabel('r (kpc)')
 pylab.ylabel('coldgasmass (M_sun)')
 pylab.xlim([0,rmax-1])
# pylab.ylim([0,3.5e8])
 pylab.savefig('M_r_more_new.png')
 pylab.clf()
 return


def M_r_more(imin=None,imax=None,rmax=10.0, coldtemperature=1e5,mod=1,binsize=0.5,ymax=1e9):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 time=[]
# cm=pylab.cm.gist_rainbow
 cm=pylab.cm.gist_heat
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 pylab.figure(0)
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = load(fn)
    current_time = pf.current_time
    time.append(current_time)
    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
    Temperature = sp["Temperature"]
    CellMassMsun = sp["CellMassMsun"]
    Radiuskpc = sp["Radiuskpc"]
    r=[]
    coldgas=[]
    cold = Temperature < coldtemperature
    for j in range(int(numpy.floor(rmax/binsize))):
       rgood= (Radiuskpc > float(j)*binsize) & (Radiuskpc < float(1+j)*binsize)
       r.append(float(j)*binsize)
       mgood=rgood & cold
       totalcoldgas=CellMassMsun[mgood]
       this_sum=sum(totalcoldgas)
       if this_sum < 0.001:
          this_sum =0.001
       coldgas.append(this_sum)
    #pylab.plot(r,coldgas,drawstyle='steps-post', color=cm(float(i-imin)/(imax-imin)))
    pylab.semilogy(r,coldgas,drawstyle='steps-post', color=cm(float(i-imin)/(imax-imin)))
 pylab.xlabel('r (kpc)')
 pylab.ylabel('Mass of Cold Gas (Solar Mass)')
 pylab.xlim([0,rmax-1])
 pylab.ylim([1e6,2e10])
 pylab.savefig('M_r_more.png')
 pylab.clf()
 return

def M_r_compare(coldtemperature=1e5, binsize=0.5,rmax=20):
 pylab.figure(0)
# pf = load('../Gaspari/'+name_pattern %(25,25))
 for name in ('../Gaspari/DD0025/stest_0025', '../high_epsilon/DD0025/stest_0025','../TP_188/DD0349/stest_0349','../Case2/DD0179/stest_0179'):
    pf=load(name)
    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
    Temperature = sp["Temperature"]
    CellMassMsun = sp["CellMassMsun"]
    Radiuskpc = sp["Radiuskpc"]
    r=[]
    coldgas=[]
    cold = Temperature < coldtemperature
    for j in range(int(numpy.floor(rmax/binsize))):
       rgood= (Radiuskpc > float(j)*binsize) & (Radiuskpc < float(1+j)*binsize)
       r.append(float(j)*binsize)
       mgood=rgood & cold
       totalcoldgas=CellMassMsun[mgood]
       coldgas.append(sum(totalcoldgas))
    pylab.plot(r,coldgas,drawstyle='steps-post')
    pylab.ylim([0,2e8])
    pylab.xlim([0,rmax])
 pylab.legend(("Gaspari","high_epsilon","Case2","Case2f1"),loc="upper right")
 pylab.savefig('M_r_compare.png')

 pylab.figure(1)
# pf = load('../Gaspari/'+name_pattern %(25,25))
 for name in ('../Gaspari/DD0152/stest_0152', '../high_epsilon/DD0192/stest_0192'):
    pf=load(name)
    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
    Temperature = sp["Temperature"]
    CellMassMsun = sp["CellMassMsun"]
    Radiuskpc = sp["Radiuskpc"]
    r=[]
    coldgas=[]
    cold = Temperature < coldtemperature
    for j in range(int(numpy.floor(rmax/binsize))):
       rgood= (Radiuskpc > float(j)*binsize) & (Radiuskpc < float(1+j)*binsize)
       r.append(float(j)*binsize)
       mgood=rgood & cold
       totalcoldgas=CellMassMsun[mgood]
       coldgas.append(sum(totalcoldgas))
    pylab.plot(r,coldgas,drawstyle='steps-post')
    pylab.ylim([0,3e10])
    pylab.xlim([0,rmax])
 pylab.legend(("Gaspari","high_epsilon"),loc="upper right")
 pylab.savefig('M_r_compare_later.png')
 return


def sigma_compare(rmax=100):
 pylab.figure(0)
 figname="sigma"
 name_list=['Gaspari/DD0025/stest_0025', 'high_epsilon/DD0025/stest_0025','TP_188/DD0349/stest_0349','Case2/DD0179/stest_0179']
 for name in name_list:
    pf=load(homedir+"/big/"+name)
    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
#    sp_out=sp.cut_region(["grid['Theta']>0.6"])
    sp_out=sp.cut_region(["grid['Temperature']>1.0e5","grid['Temperature']<5.0e7"])
#    sp_out=sp
    prof_out=BinnedProfile1D(sp_out, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)
    r_prof=prof_out["Radiuskpc"]
    prof_out.add_fields("RadialVelocity",weight="CellMassMsun")
    prof_out.add_fields("CellMassMsun",weight=None)
    Vr_prof=prof_out["RadialVelocity"]
    mass_prof=prof_out["CellMassMsun"]
    CellMassMsun = sp_out["CellMassMsun"]
    Radiuskpc = sp_out["Radiuskpc"]    
    TangentialVelocity = sp_out["TangentialVelocity"]
    RadialVelocity=sp_out["RadialVelocity"]
    Vr=numpy.interp(Radiuskpc, r_prof, Vr_prof)
    sigma=numpy.sqrt((RadialVelocity-Vr)**2.0+TangentialVelocity**2.0)
    sig_prof=[]
    for i in range(len(r_prof)-1):
        rgood=(Radiuskpc > r_prof[i]) & (Radiuskpc < r_prof[i+1])
        this_sig=numpy.multiply(numpy.array(sigma[rgood]),numpy.array(CellMassMsun[rgood]))/mass_prof[i]
        sig_prof.append(numpy.sum(this_sig)/1.0e5)
    pylab.plot(r_prof[:-1],sig_prof)#,drawstyle='steps-post')
    pylab.xlabel("r (kpc)")
    pylab.ylabel(r"$\sigma (km/s)$")
    figname+="_"+str(name[0:5])
 pylab.xlim([1,rmax])
 pylab.ylim([0,1000])
 pylab.savefig(figname+"_ICM.png")
 return


def sigma_simple(imin=None,imax=None,rmax=100,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 cm=pylab.cm.gist_heat
 pylab.figure(0)
 for i in range(imin,imax,mod):
#    pylab.figure(i)
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf=load(fn)
    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
#    sp_out=sp.cut_region(["grid['Theta']>0.6"])
    sp_out=sp.cut_region(["grid['Temperature']<1.0e8"])  # cut out high T
    prof_out=BinnedProfile1D(sp_out, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)
    r_prof=prof_out["Radiuskpc"]
    #prof_out.add_fields("VelocityMagnitude",weight="CellMassMsun")
    prof_out.add_fields("VelocityMagnitude",weight="CellVolume")
    prof_out.add_fields("VorticityMagnitude",weight="CellVolume")
    #prof_out.add_fields("VelocityMagnitude",weight=None)
    sigma_prof=prof_out["VelocityMagnitude"]/1.0e5
    curl=prof_out["VorticityMagnitude"]/1.0e5
#    pylab.plot(r_prof,sigma_prof,color=cm(float(i-imin)/(imax-imin)))#,drawstyle='steps-post')
#    pylab.ylim(0,1e3)
#    pylab.xlim(2,rmax-10)
    pylab.semilogy(r_prof,sigma_prof,color=cm(float(i-imin)/(imax-imin)))#,drawstyle='steps-post')
    pylab.ylim(0,1e3)
    pylab.xlim(2,rmax-10)
    pylab.xlabel('r (kpc)')
    pylab.ylabel('Velocity Magnitude '+r'$(km/s)$')
 pylab.tight_layout()
 pylab.savefig("sigma.png")
 return

def Vorticity_simple(imin=None,imax=None,rmax=100,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 cm=pylab.cm.gist_heat
 pylab.figure(0)
 for i in range(imin,imax,mod):
#    pylab.figure(i)
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf=load(fn)
    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
#    sp_out=sp.cut_region(["grid['Theta']>0.6"])
    sp_out=sp.cut_region(["grid['Temperature']<1.0e8"])  # cut out high T
    prof_out=BinnedProfile1D(sp_out, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)
    r_prof=prof_out["Radiuskpc"]
    #prof_out.add_fields("VelocityMagnitude",weight="CellMassMsun")
    prof_out.add_fields("VelocityMagnitude",weight="CellVolume")
    prof_out.add_fields("VorticityMagnitude",weight="CellVolume")
    #prof_out.add_fields("VelocityMagnitude",weight=None)
    sigma_prof=prof_out["VelocityMagnitude"]/1.0e5
    curl_prof=prof_out["VorticityMagnitude"]/1.0e5
    #pylab.plot(r_prof,curl_prof,color=cm(float(i-imin)/(imax-imin)))#,drawstyle='steps-post')
    pylab.semilogy(r_prof,curl_prof,color=cm(float(i-imin)/(imax-imin)))#,drawstyle='steps-post')
    pylab.ylim(0,0.7e-18)
    pylab.xlim(5,rmax-10)
    pylab.xlabel('r (kpc)')
    pylab.ylabel('Vorticity Magnitude '+r'$(s^{-1})$')
 pylab.tight_layout()
 pylab.savefig("Vorticity.png")
 return

def compressional(imin=None,imax=None,rmax=500,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 pylab.figure(0)
 cm=pylab.cm.gist_rainbow
 for i in range(imin,imax,mod):
#    pylab.figure(i)
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf=load(fn)
    sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
    sp_out=sp.cut_region(["grid['Temperature']>1.0e5","grid['Temperature']<1.0e8"])  # cut out high T and low T
    prof_out=BinnedProfile1D(sp_out, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)
    r_prof=prof_out["Radiuskpc"]
    r=r_prof*kpc   ## from kpc to cm
    prof_out.add_fields("RadialVelocity",weight="CellMassMsun")
    Vr=prof_out["RadialVelocity"]
    prof_out.add_fields("CoolingTime",weight="CellMassMsun")
    t_cool=prof_out["CoolingTime"]   ## in yr   ##*3.16e7  
    nn=len(r)
    del_v=numpy.zeros(nn)
    for j in range(nn-1):
       del_v[j]=(Vr[j+1]-Vr[j])/(r[j+1]-r[j])+(Vr[j+1]+Vr[j])*2.0/(r[j+1]+r[j])
    t_comp=(-(1.0/del_v)/(5.0/3.0-1.0))/3.16e7    #s->yr
#    pylab.plot(r_prof,t_comp/t_cool)
    good=t_comp>0
    pylab.loglog(r_prof[good],t_comp[good]/t_cool[good],color=cm(float(i-imin)/(imax-imin)))
    pylab.xlim(0.5,rmax-10)
    pylab.xlabel('r (kpc)')
    pylab.ylabel('t_comp/t_cool')
# pylab.savefig("t_comp_t_cool.png")
 pylab.savefig("t_comp.png")
 return


def X_composite(i=50,widthkpc=100,thicknesskpc=5000,theta=0.9,scale=[0.1,1.0,1.0], N=512):
 import yt.analysis_modules.spectral_integrator.spectral_frequency_integrator as SI
 SI.add_xray_emissivity_field(0.3,1.2,with_metals=False,constant_metallicity=0.5)   #red
 SI.add_xray_emissivity_field(1.2,2.0,with_metals=False,constant_metallicity=0.5)   #green
 SI.add_xray_emissivity_field(2.0,7.0,with_metals=False,constant_metallicity=0.5)   #blue
 fn = name_pattern %(i,i)
 pf = load(fn)
 center=[0.5,0.5,0.5]
 thickness = thicknesskpc/pf['kpc']
 width=widthkpc/pf['kpc']
 L = [0,numpy.sin(theta), numpy.cos(theta)]
 R = off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_0.3_1.2keV")
 G = off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_1.2_2.0keV")
 B = off_axis_projection(pf, center, L, (width,width,thickness), N, "Xray_Emissivity_2.0_7.0keV")
 write_projection(R, "%s_offaxis_projection%.1f_Xray_R.png" % (pf, theta), colorbar_label="X ray Surface Brightness "+r"($ergs \, s^{-1} \, cm^{-2}$)")
 write_projection(G, "%s_offaxis_projection%.1f_Xray_G.png" % (pf, theta), colorbar_label="X ray Surface Brightness "+r"($ergs \, s^{-1} \, cm^{-2}$)")
 write_projection(B, "%s_offaxis_projection%.1f_Xray_B.png" % (pf, theta), colorbar_label="X ray Surface Brightness "+r"($ergs \, s^{-1} \, cm^{-2}$)")
# im_r=numpy.zeros((R.shape[0],R.shape[1],3),dtype=float)
 im_r=numpy.zeros((R.shape[0],R.shape[1],3),'uint8')
 im_g=numpy.zeros_like(im_r)
 im_b=numpy.zeros_like(im_r)
 RGB=numpy.zeros_like(im_r)
 high=R>scale[0]*numpy.max(R)   ##rescale R
 R[high]=scale[0]*numpy.max(R)
 high=G>scale[1]*numpy.max(G)   ##rescale R
 G[high]=scale[1]*numpy.max(G) 
 high=B>scale[2]*numpy.max(B)   ##rescale R
 B[high]=scale[2]*numpy.max(B)
 im_r[:,:,0]=255.9*(R-numpy.min(R))/(numpy.max(R)-numpy.min(R))
 im_g[:,:,1]=255.9*(G-numpy.min(G))/(numpy.max(G)-numpy.min(G))
 im_b[:,:,2]=255.9*(B-numpy.min(B))/(numpy.max(B)-numpy.min(B))
 RGB[:,:,0]=255.9*(R-numpy.min(R))/(numpy.max(R)-numpy.min(R))
 RGB[:,:,1]=255.9*(G-numpy.min(G))/(numpy.max(G)-numpy.min(G))
 RGB[:,:,2]=255.9*(B-numpy.min(B))/(numpy.max(B)-numpy.min(B))
## im_b[:,:,0]=numpy.log(255.9)*numpy.log(B/numpy.min(B))/numpy.log(numpy.max(B)/numpy.min(B))
 stest = 'stest_%04i' %i
# print im_r,im_g,im_b
 pylab.figure(3)
 pylab.imshow(im_b,origin="lower", interpolation="nearest")
 pylab.savefig("im_b.png")
 pylab.figure(4)
 pylab.imshow(im_g,origin="lower", interpolation="nearest")
 pylab.savefig("im_g.png")
 pylab.figure(5)
 pylab.imshow(im_r,origin="lower", interpolation="nearest")
 pylab.savefig("im_r.png")
 pylab.figure(1)
 pylab.imshow(RGB,origin="lower", interpolation="nearest")
 pylab.tick_params(axis='both',which='both',bottom='off',top='off',left='off',right='off',labelbottom='off',labelleft='off',labelright='off')
 pylab.savefig(stest+"X_composite.png")
 return

def Xray_t(imin=None,imax=None,rmax=None):
 import numpy as numpy
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=numpy.array([50,500])
 time=numpy.zeros(imax-imin)
 Xray=numpy.zeros(((imax-imin),rmax.size))
 pylab.figure(0)
 for j in range(rmax.size):
   for i in range(imin, imax):
     stest = 'stest_%04i' %i
     ii=i-imin
     fn = name_pattern %(i,i)
     pf = load(fn)
     time[ii] = pf.current_time
     sphere = pf.h.sphere([0.5, 0.5, 0.5], rmax[j]/pf['kpc'])
     for k in range(3):
       sphere.set_field_parameter("Ephoton", k+0.5)
       Xray[ii,j]+=sphere.quantities["FreeFree_Luminosity"]()
   print "Xray", Xray[ii,j]
   pylab.semilogy(time,Xray[:,j])
 pylab.xlabel('Time (code unit)')
 pylab.ylabel('FreeFree_Luminosity')
 pylab.savefig(stest+'FreeFree_Luminosity.png')
 return

def Xray_Luminosity(imin=None,imax=None,rmax=500.0):    #Sam
 import yt.analysis_modules.spectral_integrator.spectral_frequency_integrator as SI
 from yt.data_objects.api import add_field
 xray_min = numpy.array([0.5])
 xray_max = numpy.array([3.0])
 patt=homedir+"/code/cloudy_xray_02_48_log/cloudy_xray_02_48_log_run%i.dat"
 table = SI.create_table_from_textfiles(patt, (15,-6,1),(200, 0.202759, 47.3468080), (61,5,8))
# SFI = SI.SpectralFrequencyIntegrator(table,
#        ["NumberDensity", "Temperature"], (-6,1,5,8),
#        (numpy.log10(0.202759),numpy.log10(47.3468080)))
 SFI = SI.SpectralFrequencyIntegrator(table,["H_NumberDensity", "Temperature"], (-6,1,5,8),(numpy.log10(0.202759),numpy.log10(47.3468080)))

 print "We have integrated the table"
 for bin in range(len(xray_min)):
        name = SFI.add_frequency_bin_field(xray_min[bin], xray_max[bin])
        print "We have added the frequency bin called"
        print name
 def _XRayEmission(field,data):
    return data[name]*data["CellVolume"]
 add_field("XRayEmission",units=r"\rm{ergs/s}",function=_XRayEmission)
 if imin is None:
    imin=0
 if imax is None:
    imax=imin + 1
# cm=pylab.cm.gist_rainbow
 time=numpy.zeros(imax-imin)
 Lx=numpy.zeros(imax-imin)
 for i in range(imin,imax):
    ii=i-imin
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = load(fn)
    time[ii] = pf.current_time
    center=[0.5, 0.5, 0.5]
#    sp = pf.h.sphere(center, (rmax, "kpc"))
    sp = pf.h.sphere(center, rmax/pf['kpc'])
    Lx[ii] = sp.quantities["TotalQuantity"](["XRayEmission"])[0]
 print "XRayEmission=", Lx
 pylab.figure(0)
 pylab.semilogy(time, Lx)
 pylab.xlabel('time')
 pylab.ylabel('X-ray Luminosity')
 pylab.savefig(stest+'_Xray_Luminosity%ikpc.png' %rmax)
 pylab.clf()
 return




def ColdGasMass_yt(imin=None,imax=None,rmax=None, coldtemperature=None):
 import numpy as numpy
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=0.5   ##kpc
 if coldtemperature is None:
    coldtemperature = 3e4
 kpc=3.086e21
 time=[]
 coldgas=[]
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   current_time = pf.current_time
   time.append(current_time)
   center= numpy.array([0.5, 0.5, 0.5])
   LE=center-rmax/pf['kpc']
   RE=center+rmax/pf['kpc']
   sp = pf.h.region(center,LE,RE)
   temperature = sp["Temperature"]
   cellmass = sp["CellMassMsun"]
   cold = temperature < coldtemperature
   totalcoldgas=cellmass[cold]
   print "Total cold gas mass in", stest, totalcoldgas
   coldgas.append(sum(totalcoldgas))
 print time, coldgas
 pylab.figure(0)
 pylab.plot(time,coldgas,marker="o")
 pylab.xlabel('Time (code unit)')
 pylab.ylabel('coldgasmass (M_sun)')
 pylab.savefig('coldgasmass'+str(rmax)+'kpc'+str(coldtemperature)+'K.png')
 pylab.clf()
 return


def ColdGasMass_enzo(fn="MT.out"):
 data = numpy.genfromtxt(fn)
 pylab.figure(1)
 Time0 = min(data[:,9])
 Time = (data[:,9] - Time0)/1.7e-4
 Switch = data[:,11]
 ColdGas = data[:,12]
 on = Switch > 0.1
 off = Switch < 0.1
 pylab.plot(Time[on],ColdGas[on],marker='.',linestyle='None')
 pylab.plot(Time[off],ColdGas[off],marker='.',color='red',linestyle='None')
# pylab.plot(time,coldgas,linestyle='None',marker="o")
 pylab.xlabel("Time (Myr)")
# lim_temp=np.array([0.064-0.035,0.072-0.035])/1.7e-4
# pylab.xlim(lim_temp)
 pylab.ylabel('coldgasmass (M_sun)')
 pylab.savefig('coldgasmass_enzo.png')
 return



def ColdGasMass(imin=None,imax=None,rmax=40,mod=20,log=False,Halpha=False,L_Particle=True):
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 coldtemperature = 1e5
 time=[]
 coldgas=[]
 Lx=[]
 Ly=[]
 Lz=[]
 Lx_P=[]
 Ly_P=[]
 Lz_P=[]
 L_P=[]
 colddense=[]
 L_Halpha=[]
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = load(fn)
    current_time = pf.current_time
    time.append(current_time)
    center= numpy.array([0.5, 0.5, 0.5])
    sp_all = pf.h.sphere(center,(rmax, 'kpc'))
#     sp_dense=sp.cut_region(["grid['Density']>1.7e-24"])
#     [totalcoldgas, totalLx, totalLy, totalLz]=sp.quantities["TotalQuantity"](["ColdCellMassMsun","Cold_AngularMomentumX", "Cold_AngularMomentumY", "Cold_AngularMomentumZ"])
#     colddense.append(sp_dense.quantities["TotalQuantity"](["ColdCellMassMsun"]))
    sp=sp_all.cut_region(["grid['Temperature']<1.0e5"])
    sp_dense=sp_all.cut_region(["grid['Temperature']<1.0e5","grid['Density']>1.7e-23"])
    [totalcoldgas, totalLx, totalLy, totalLz]=sp.quantities["TotalQuantity"](["CellMassMsun","AngularMomentumX", "AngularMomentumY", "AngularMomentumZ"])
    colddense.append(sp_dense.quantities["TotalQuantity"](["CellMassMsun"]))
    coldgas.append(totalcoldgas)
    Lx.append(totalLx)
    Ly.append(totalLy)
    Lz.append(totalLz)
    if L_Particle:
      try: 
         [particleLx,particleLy,particleLz,particleL]=sp.quantities["TotalQuantity"](["ParticleAngularMomentumX","ParticleAngularMomentumY","ParticleAngularMomentumZ","ParticleAngularMomentum"])
         Lx_P.append(particleLx)
         Ly_P.append(particleLy)
         Lz_P.append(particleLz)
         L_P.append(particleL)
      except:
         pass
    if Halpha:
      try:
         L_Halpha.append(sp.quantities["TotalQuantity"](["Halpha_Luminosity"])[0])
      except:
         L_Halpha.append(1.0e-10)
 time[:] = [x*TimeUnits/3.16e16 for x in time]  # from code unit to Gyr
 numpy.savez("ColdGasMass_in%ikpc" %rmax, coldgas=coldgas, Lx=Lx, Ly=Ly, Lz=Lz,timeGyr=time,colddense=colddense,L_Halpha=L_Halpha,Lx_P=Lx_P,Ly_P=Ly_P,Lz_P=Lz_P,L_P=L_P)
 L=(numpy.array(Lx)**2.0+numpy.array(Ly)**2.0+numpy.array(Lz)**2.0)**0.5
 pylab.figure(0)
 if log:
   pylab.semilogy(time,coldgas)
   pylab.xlabel("Time (Gyr)")
   pylab.ylabel('Cold Gas Mass (Solar Mass)')
   pylab.tight_layout()
   pylab.savefig('coldgasmass_log_in%ikpc.png' %rmax)
 else:
   pylab.plot(time,coldgas)
   pylab.xlabel("Time (Gyr)")
   pylab.ylabel('Cold Gas Mass (Solar Mass)')
   pylab.tight_layout()
   pylab.savefig('coldgasmass_in%ikpc.png' %rmax)
 pylab.clf()
 pylab.figure(1)
 pylab.plot(time,L)
 pylab.plot(time,Lx)
 pylab.plot(time,Ly)
 pylab.plot(time,Lz)
 pylab.legend(("L","Lx","Ly","Lz"),loc="upper left",prop={'size':"small"})
 pylab.xlabel("Time (Gyr)")
 pylab.ylabel("Angular Momentum"+r"$(\rm{g}\ \rm{cm}^{2} \rm{s}^{-1})$")
 pylab.tight_layout()
 pylab.savefig('coldgasL_in%ikpc.png' %rmax)
 pylab.clf()
 return

def M_T(imin=None,imax=None,rmax=None, height=100,mod=1):
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 cm=pylab.cm.gist_heat
 if imin is None:
    imin=0
 if imax is None:
    imax=1
 if rmax is None:
    rmax=32.0    #kpc
 classic = numpy.genfromtxt(homedir+"/code/massfrac.qdp")
 obs = numpy.genfromtxt(homedir+"/code/massfrac_obs.qdp")
 T_classic=classic[:,0]
 M_classic=classic[:,1]
 [T_obs,T_minus,T_plus,sig1,sig2,sig3]=numpy.array_split(obs,6,axis=1)
 name_pattern = 'DD%04i/stest_%04i'
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   t=pf.current_time
   t*=(TimeUnits/3.16e16)   ## code unit to Gyr
#   sp = pf.h.sphere([0.5, 0.5, 0.5],(rmax, 'kpc'))
#   temperature = sp["Temperature"]/keV
#   cellmass = sp["CellMassMsun"]
   disk=pf.h.disk([0.5, 0.5, 0.5], [0,1,1],rmax/pf['kpc'], height/pf['kpc'])
   temperature = disk["Temperature"]/keV
   cellmass = disk["CellMassMsun"]
   T_min=0
   binsize=0.2
   arraysize=26
   T_max=T_min+binsize*(arraysize-1)
   T_keV=numpy.linspace(T_min+0.5*binsize, T_max+0.5*binsize, arraysize)
   M_Msun=numpy.zeros(arraysize)
   for j in range(arraysize):
      cold = [(temperature > T_keV[j]-0.5*binsize) & (temperature <= T_keV[j]+0.5*binsize)]
      M_Msun[j]=sum(cellmass[cold])
   f=pylab.figure(i)
   pylab.loglog(T_keV, M_Msun,drawstyle='steps-mid',color='black')
   pylab.xlim([0.3,6])
   pylab.ylim([1e7,1e12])
   pylab.loglog(T_classic,M_classic,color='red',linestyle='--')
   T_left=(T_obs+T_minus).ravel()
   T_right=(T_obs+T_plus).ravel()
   sig2=sig2.ravel()
   pylab.loglog([T_left,T_right],[sig2,sig2],'k',color='blue')
   pylab.loglog(T_obs,sig2,linestyle='None',marker='o',color='blue')
   pylab.legend(("Simulation","Classic Cooling Flow","Perseus Obs."),loc="upper left",fontsize=16)
   pylab.xlabel('Temperature (keV)')
   pylab.ylabel('Mass/keV (Sollar Mass/keV)')
   f.text(0.68, 0.82, "t = %.2f Gyr"%t, fontdict={'size':20, 'color':'black'})
   pylab.xticks((0.5,1.0,2.0,5.0),('0.5','1.0','2.0','5.0'))
   pylab.savefig(stest+'_M_T_in%ikpc.png' %rmax)
 return


def M_T_old(imin=None,rmax=32, height=1000, n_bins=10,Xbins=[0.3,1.2,2.0,3.5,5.0,7.0]):  #this is actually the spectrum
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 cm=pylab.cm.gist_heat
 if imin is None:
    imin=0
 classic = numpy.genfromtxt(homedir+"/code/massfrac.qdp")
 obs = numpy.genfromtxt(homedir+"/code/massfrac_obs.qdp")
 T_classic=classic[:,0]
 M_classic=classic[:,1]
 [T_obs,T_minus,T_plus,sig1,sig2,sig3]=numpy.array_split(obs,6,axis=1)
 pylab.figure(0)
 pylab.loglog(T_classic,M_classic,color='red',linestyle='--')
 T_left=(T_obs+T_minus).ravel()
 T_right=(T_obs+T_plus).ravel()
 sig2=sig2.ravel()
 for i in range(len(T_left)):
   pylab.fill_between([T_left[i],T_right[i]],[sig2[i],sig2[i]],1.0e6,alpha=0.2,color='blue')
 pylab.loglog(T_obs,sig2,linestyle='None',marker='o',color='blue')
 pylab.legend(("Classic Cooling Flow","Perseus Obs."),loc="upper left",fontsize=16)
 name_pattern = 'DD%04i/stest_%04i'
 i=imin
 stest = 'stest_%04i' %i
 fn = name_pattern %(i,i)
 pf = load(fn)
 current_time = pf.current_time
 sphere=pf.h.sphere([0.5, 0.5, 0.5],rmax/pf['kpc'])
 prof_out=BinnedProfile1D(sphere,n_bins, "Radiuskpc", 1,rmax,lazy_reader=True)
 r_prof=numpy.array(prof_out["Radiuskpc"])
 print "r_prof=", r_prof
 prof_out.add_fields("Temperature",weight="CoolingRate")
 Tmean=numpy.array(prof_out["Temperature"])
# import yt.analysis_modules.spectral_integrator.spectral_frequency_integrator as SI
 nXbins=len(Xbins)
 from yt.analysis_modules.spectral_integrator.api import add_xray_luminosity_field 
 for k in range(0,nXbins-1):
    add_xray_luminosity_field(Xbins[k], Xbins[k+1], with_metals=False, constant_metallicity=0.5)
 X_ICM=numpy.zeros(nXbins-1)
 X_Clump=numpy.zeros(nXbins-1)
 X_all=numpy.zeros(nXbins-1)
 for j in range(0,n_bins-1):
   sp1 = pf.h.sphere([0.5, 0.5, 0.5], r_prof[j]/pf['kpc'])
   sp2 = pf.h.sphere([0.5, 0.5, 0.5], r_prof[j+1]/pf['kpc'])
   shell = pf.h.boolean([sp2, "NOT", sp1])
#   Tmean=shell.quantities["WeightedAverageQuantity"]("Temperature", "CoolingRate")
   Temp=shell["Temperature"]
   ICM=Temp>=Tmean[j]
   Clump=Temp<Tmean[j]
   for k in range(0,nXbins-1):
       Xray_shell=shell["Xray_Luminosity_%.1f_%.1fkeV" %(Xbins[k], Xbins[k+1])]
       X_ICM[k]+=numpy.sum(Xray_shell[ICM])
       X_Clump[k]+=numpy.sum(Xray_shell[Clump])
   prof1=BinnedProfile1D(shell, 26, "Temperature", 0.5e7, 7.0e7, log_space=False, lazy_reader=True)
   prof1.add_fields("CellMassMsun",weight=None)
   T_plot=prof1["Temperature"]   #in K
   Mass1=prof1["CellMassMsun"]
   low=T_plot<Tmean[j]
   Mass1[low]=0.001
   if j == 0:
      Mass=numpy.array(Mass1)
   else:
      Mass=Mass+numpy.array(Mass1)
# big_sphere = pf.h.sphere([0.5, 0.5, 0.5],rmax/pf['kpc'])
# bin_sphere=big_sphere.cut_region(["Radiuskpc">1])
 disk=pf.h.disk([0.5, 0.5, 0.5], [0,1,1],rmax/pf['kpc'], 1000/pf['kpc'])
# disk=disk.cut_region(["Radiuskpc">1])
 for k in range(0,nXbins-1):
#     X_all[k]=big_sphere.quantities["TotalQuantity"](["Xray_Luminosity_%.1f_%.1fkeV" %(Xbins[k], Xbins[k+1])])[0]
     X_all[k]=disk.quantities["TotalQuantity"](["Xray_Luminosity_%.1f_%.1fkeV" %(Xbins[k], Xbins[k+1])])[0]
 temperature = sphere["Temperature"]/keV
 cellmass = sphere["CellMassMsun"]
 T_min=0
 binsize=0.2
 arraysize=26
 T_max=T_min+binsize*(arraysize-1)
 T_keV=numpy.linspace(T_min+0.5*binsize, T_max+0.5*binsize, arraysize)
 M_Msun=numpy.zeros(arraysize)
 for j in range(arraysize):
      cold = [(temperature > T_keV[j]-0.5*binsize) & (temperature <= T_keV[j]+0.5*binsize)]
      M_Msun[j]=sum(cellmass[cold])
 pylab.loglog(T_keV, M_Msun,drawstyle='steps-mid')
 pylab.loglog(T_plot/keV, Mass,drawstyle='steps-mid',)
 pylab.xlim([0.3,6])
 pylab.ylim([1e7,1e12])
 pylab.xlabel('Temperature (keV)')
 pylab.ylabel('Mass (Sollar Mass)')
 pylab.xticks((0.5,1.0,2.0,5.0),('0.5','1.0','2.0','5.0'))
 pylab.savefig(stest+'M_T_in%ikpc.png' %rmax)
 pylab.clf()
 numpy.savez(stest+"LX", X_ICM=X_ICM, X_Clump=X_Clump,Xbins=Xbins,X_all=X_all)
 return

def M_T_Unsharp(imin=None,imax=None, rmax=32, height=100, n_bins=10,TbinskeV=[0.3,1.2,2.0,3.5,5.0,7.0]):
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 cm=pylab.cm.gist_heat
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 classic = numpy.genfromtxt(homedir+"/code/massfrac.qdp")
 obs = numpy.genfromtxt(homedir+"/code/massfrac_obs.qdp")
 T_classic=classic[:,0]
 M_classic=classic[:,1]
 [T_obs,T_minus,T_plus,sig1,sig2,sig3]=numpy.array_split(obs,6,axis=1)
 pylab.figure(0)
 pylab.loglog(T_classic,M_classic,color='red',linestyle='--')
 T_left=(T_obs+T_minus).ravel()
 T_right=(T_obs+T_plus).ravel()
 sig2=sig2.ravel()
 for i in range(len(T_left)):
   pylab.fill_between([T_left[i],T_right[i]],[sig2[i],sig2[i]],1.0e6,alpha=0.2,color='blue')
 pylab.loglog(T_obs,sig2,linestyle='None',marker='o',color='blue')
 pylab.legend(("Classic Cooling Flow","Perseus Obs."),loc="upper left",fontsize=16)
 name_pattern = 'DD%04i/stest_%04i'
 Tbins=numpy.array(TbinskeV)*keV
 nTbins=len(Tbins)
 X_ICM=numpy.zeros((imax-imin,nTbins-1))
 X_Clump=numpy.zeros((imax-imin,nTbins-1))
 X_all=numpy.zeros((imax-imin,nTbins-1))
 for i in range(imin, imax):
   ii=i-imin
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   current_time = pf.current_time
   sphere=pf.h.sphere([0.5, 0.5, 0.5],rmax/pf['kpc'])
   prof_out=BinnedProfile1D(sphere,n_bins, "Radiuskpc", 1,rmax,lazy_reader=True)
   r_prof=numpy.array(prof_out["Radiuskpc"])
   prof_out.add_fields("Temperature",weight="CoolingRate")
   Tmean=numpy.array(prof_out["Temperature"])
   from yt.analysis_modules.spectral_integrator.api import add_xray_luminosity_field
   add_xray_luminosity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
   for j in range(0,n_bins-1):  # j is r(j)
     sp1 = pf.h.sphere([0.5, 0.5, 0.5], r_prof[j]/pf['kpc'])
     sp2 = pf.h.sphere([0.5, 0.5, 0.5], r_prof[j+1]/pf['kpc'])
     shell = pf.h.boolean([sp2, "NOT", sp1])
     Temp=shell["Temperature"]
     Xray_shell=shell["Xray_Luminosity_0.5_9.9keV"]
     for k in range(0,nTbins-1):  # k is T(k)
        good = (Temp>=Tbins[k]) & (Temp<Tbins[k+1])
        Clump = good & (Temp <Tmean[j])
        ICM = good & (Temp >= Tmean[j])
        X_ICM[ii,k]+=numpy.sum(Xray_shell[ICM])
        X_Clump[ii,k]+=numpy.sum(Xray_shell[Clump])
     prof1=BinnedProfile1D(shell, 26, "Temperature", 0.5e7, 7.0e7, log_space=False, lazy_reader=True)
     prof1.add_fields("CellMassMsun",weight=None)
     T_plot=prof1["Temperature"]   #in K
     Mass1=prof1["CellMassMsun"]
     low=T_plot<Tmean[j]
     Mass1[low]=0.001
     if j == 0:
         Mass=numpy.array(Mass1)
     else:
         Mass=Mass+numpy.array(Mass1)
# big_sphere = pf.h.sphere([0.5, 0.5, 0.5],rmax/pf['kpc'])
# bin_sphere=big_sphere.cut_region(["Radiuskpc">1])
   disk=pf.h.disk([0.5, 0.5, 0.5], [0,1,1],rmax/pf['kpc'], 1000/pf['kpc'])
#   disk=disk.cut_region(["Radiuskpc">1.0])
   T_all=disk["Temperature"]
   L_all=disk["Xray_Luminosity_0.5_9.9keV"]
   for k in range(0,nTbins-1):
      good = (T_all>Tbins[k]) & (T_all<Tbins[k+1])
      X_all[ii,k]=numpy.sum(L_all[good])
   temperature = sphere["Temperature"]/keV
   cellmass = sphere["CellMassMsun"]
   T_min=0
   binsize=0.2
   arraysize=26
   T_max=T_min+binsize*(arraysize-1)
   T_keV=numpy.linspace(T_min+0.5*binsize, T_max+0.5*binsize, arraysize)
   M_Msun=numpy.zeros(arraysize)
   for j in range(arraysize):
      cold = [(temperature > T_keV[j]-0.5*binsize) & (temperature <= T_keV[j]+0.5*binsize)]
      M_Msun[j]=sum(cellmass[cold])
   pylab.loglog(T_keV, M_Msun,drawstyle='steps-mid')
   pylab.loglog(T_plot/keV, Mass,drawstyle='steps-mid',)
   pylab.xlim([0.3,6])
   pylab.ylim([1e7,1e12])
   pylab.xlabel('Temperature (keV)')
   pylab.ylabel('Mass (Sollar Mass)')
   pylab.xticks((0.5,1.0,2.0,5.0),('0.5','1.0','2.0','5.0'))
   pylab.savefig(stest+'M_T_in%ikpc.png' %rmax)
   pylab.clf()
 numpy.savez("LX%04ito%04i"%(imin,imax), X_ICM=X_ICM, X_Clump=X_Clump,Tbins=TbinskeV,X_all=X_all)
 return



def PhasePlots(imin=None,imax=None,rmax=100,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=1
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","RadialVelocity","CellMassMsun"],weight=None,y_log=False)
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","ColdRadialVelocity","Entropy"],y_log=False,y_bounds=(-0.1e9,0.2e9))
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","ColdRadialVelocity","CellMassMsun"],y_log=False,y_bounds=(-0.1e9,0.2e9))
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","RadialVelocity","CellMassMsun"],y_log=False,y_bounds=(-0.1e9,0.2e9))
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","RadialVelocity","CellMassMsun"],y_log=False)
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","ColdRadialVelocity","SpecificEntropy"],y_log=False,y_bounds=(-0.2e9,0.1e9))
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","Temperature","CellMassMsun"],y_bounds=(1e3,1e8))
   pc.add_phase_sphere(rmax,'kpc',["Density","Temperature","CellMassMsun"],x_bounds=(5e-26,1e-22),y_bounds=(1e3,5e8))  #average cell mass
   pc.set_xlim(5e-26,1e-22)
   pc.set_ylim(1e3,5e8)
#   pc.add_phase_sphere(rmax,'kpc',["Density","Temperature","CellMassMsun"],weight=None) # add up
#   pc.add_phase_sphere(rmax,'kpc',["Density","Temperature","ColdCellMassMsun"])
#   pc.add_phase_sphere(rmax,'kpc',["y","z","ColdCellMassMsun"],x_log=False,y_log=False)
#   pc.add_phase_sphere(rmax,'kpc',["y","z","Cold_x_velocity"],x_log=False,y_log=False)
# below is for filaments
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","Theta","ColdCellMassMsun"],weight=None,x_log=False,y_log=False).set_ylim(0,1.57)
#   pc.add_phase_sphere(rmax,'kpc',["y","z","x-velocity"],weight="XRayEmissivity", x_log=False,y_log=False)
#   pc.add_phase_sphere(rmax,'kpc',["y","z","Cold_x_velocity"],weight="XRayEmissivity",x_log=False,y_log=False)
# end filaments
#   pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","Cooling_Time","CellMassMsun"], weight=None)
   #pc.add_phase_sphere(rmax,'kpc',["Radiuskpc","Temperature","CellMassMsun"], weight=None,x_bounds=(1,rmax),y_bounds=(1e3,1e9),lazy_reader=True).set_xlim(1,rmax)
#   pc.add_phase_sphere(rmax,'kpc',["Density","JeansMassMsun","CellMassMsun"], weight=None,lazy_reader=True)
#   pc.add_phase_sphere(rmax,"kpc",["Radiuskpc","CoolingTime","CellMassMsun"], weight=None,x_bounds=(1,rmax),y_bounds=(1e6,1e10),lazy_reader=True).set_xlim(1,rmax)
####ratio vs r
#   pc.add_phase_sphere(rmax,"kpc",["Radiuskpc","tcooltdyn","CellMassMsun"], weight=None,x_bounds=(1,rmax),y_bounds=(0.5,100),lazy_reader=True)
#   pc.set_xlim(1,rmax)
#   pc.set_ylim(0.5,100)
#   pc.set_zlim(2e4,8e11)
####end ratio vs r
#   pc.add_phase_sphere(rmax,"kpc",["Radiuskpc","ColdRadialVelocity","CellMassMsun"],weight=None,y_log=False) 
   pc.save()
 return


def _ColdNegativeRadialMassFlux(field, data):
   Temp=np.array(data["Temperature"])
   not_cold=Temp > 1.0e5
   NegFlux= na.minimum(data["RadialVelocity"],0.0) * data["Density"]
   NegFlux[not_cold]=0
   return NegFlux
add_field("ColdNegativeRadialMassFlux", function=_ColdNegativeRadialMassFlux, units=r"\rm{g}/\rm{cm}^{2}", take_log=False)



def Bondi_BCG(imin=None,imax=None,rmax=10,plot_each=0):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin + 1
 M_BH=3.4e8
 s=0.9
 time=numpy.zeros(imax-imin)
 Mdot_A=numpy.zeros(imax-imin)   #the actual accretion rate at Bondi radius at each time step
 Mdot_B=numpy.zeros(imax-imin)   #the Bondi rate at Bondi radius at each time step
 Mdot_C=numpy.zeros(imax-imin)   #the Bondi rate calculated at r=1kpc
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   ii=i-imin
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","Density"],weight="CellMassMsun")
   r=pc.plots[-1].data["Radiuskpc"]
   rho=pc.plots[-1].data["Density"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","SoundSpeed"],weight="CellMassMsun")
   Cs=pc.plots[-1].data["SoundSpeed"]
   x1=r/(Rs*pf['kpc'])
   M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass #DM+BCG + BH , cgs
   Mdot_Bondi=(4.0*numpy.pi*(GravConst*M_dyn)**2*rho/Cs**3)*3.15e7/SolarMass  #need to change
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","NegativeRadialMassFlux"])
   NeFlux=pc.plots[-1].data["NegativeRadialMassFlux"]
   NedMdt=(NeFlux/SolarMass)*4.0*3.14*pow(r*kpc,2.0)*3.15e7  ## from cgs to M_sun/yr
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","ColdNegativeRadialMassFlux"])
   ColdNeFlux=pc.plots[-1].data["ColdNegativeRadialMassFlux"]
   ColdNedMdt=(ColdNeFlux/SolarMass)*4.0*3.14*pow(r*kpc,2.0)*3.15e7  ## from cgs to M_sun/yr
   print "NedMdt=", NedMdt
   print "ColdNedMdt=", ColdNedMdt
   BCGg=((r**0.5975/3.206e-7)**s+(r**1.849/1.861e-6)**s)**(-1.0/s) 
   SMBHg=GravConst*M_BH*SolarMass/(r*kpc)**2
   v_esc=((BCGg+SMBHg)*r*kpc*2.0)**0.5     #cgs
   time[ii] = pf.current_time
   for j in range(r.size):
      r_B=r[j]
      if v_esc[j]<Cs[j]:
        break
   print "BCG r_B= ", r_B
   Mdot_A[ii] = -NedMdt[j]
   Mdot_B[ii] = Mdot_Bondi[j]
   Mdot_C[ii] = numpy.interp(1.0, r, Mdot_Bondi)
   if plot_each==1:
      pylab.figure(ii)
      pylab.loglog(r,Mdot_Bondi)
      pylab.loglog(r,-NedMdt)
      pylab.loglog(r,-ColdNedMdt)   #test
      pylab.axvline(x=r_B)
      pylab.ylim([0.001,1e4])
      pylab.xlim([0.05,rmax])
      pylab.xlabel("r (kpc)")
      pylab.ylabel('Mdot (M_sun/yr)')
#      pylab.legend(("Bondi Rate", "Inflow"),loc="upper left")
      pylab.legend(("Bondi Rate", "Inflow","Cold Inflow"),loc="upper left")
      pylab.savefig(stest+'_Bondi.png')
      pylab.clf()
 pylab.figure(0)
 pylab.semilogy(time,Mdot_B)
 pylab.semilogy(time,Mdot_A)
 pylab.xlabel("time (code units)")
 pylab.ylabel('Mdot (M_sun/yr)')
 pylab.legend(("Bondi Rate", "Inflow"),loc="upper right")
 pylab.savefig('Bondi_time_BCG.png')
 pylab.figure(1)
 ratio=Mdot_C/Mdot_B
 pylab.semilogy(time,ratio)
 pylab.xlabel("time (code units)")
 pylab.ylabel('ratio')
 pylab.savefig('Bondi_ratio_BCG.png')
 pylab.clf()
 return

def Bondi_BCG_Hot(imin=None,imax=None,rmax=10,plot_each=0):  #use BinnedProfile1D instead of PC, also cut out the jet (Vr>0)
 if imin is None:
    imin=0
 if imax is None:
    imax=imin + 1
 M_BH=3.4e8
 s=0.9
 time=numpy.zeros(imax-imin)
 Mdot_A=numpy.zeros(imax-imin)   #the actual accretion rate at Bondi radius at each time step
 Mdot_B=numpy.zeros(imax-imin)   #the Bondi rate at Bondi radius at each time step
 Mdot_C=numpy.zeros(imax-imin)   #the Bondi rate calculated at r=1kpc
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   ii=i-imin
   pf = load(fn)
   sp=pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
#   sp_in=sp.cut_region(["grid['RadialVelocity']<0"])
   sp_hot=sp.cut_region(["grid['Temperature']>1e5"])
   prof=BinnedProfile1D(sp_hot,300,"Radiuskpc", 1e-2,rmax,lazy_reader=True)
   prof.add_fields("Density",weight="CellMassMsun")
   prof.add_fields("SoundSpeed",weight="CellMassMsun")
   r=prof["Radiuskpc"]
   rho=prof["Density"]
   Cs=prof["SoundSpeed"]
   print r.size,"r=",r
   print Cs.size,"Cs=",Cs
   x1=r/(Rs*pf['kpc'])
   M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass #DM+BCG + BH , cgs
   Mdot_Bondi=(4.0*na.pi*(GravConst*M_dyn)**2*rho/Cs**3)*3.15e7/SolarMass
   prof.add_fields("NegativeRadialMassFlux",weight="CellMassMsun")  #???weight?
   NeFlux=prof["NegativeRadialMassFlux"]
   NedMdt=(NeFlux/SolarMass)*4.0*3.14*pow(r*kpc,2.0)*3.15e7  ## from cgs to M_sun/yr
   BCGg=((r**0.5975/3.206e-7)**s+(r**1.849/1.861e-6)**s)**(-1.0/s)
   SMBHg=GravConst*M_BH*SolarMass/(r*kpc)**2
   v_esc=((BCGg+SMBHg)*r*kpc*2.0)**0.5     #cgs
   print v_esc.size, "v_esc=",v_esc
   time[ii] = pf.current_time
   for j in range(r.size):
      r_B=r[j]
      if v_esc[j]<Cs[j]:
        break 
   print "BCG r_B= ", r_B
   Mdot_A[ii] = -NedMdt[j]
   Mdot_B[ii] = Mdot_Bondi[j]
   Mdot_C[ii] = numpy.interp(1.0, r, Mdot_Bondi)
   if plot_each==1:
      pylab.figure(ii)
      pylab.loglog(r,Mdot_Bondi)
      pylab.loglog(r,-NedMdt)
      pylab.axvline(x=r_B)
      pylab.ylim([0.001,1e4])
      pylab.xlim([0.05,rmax])
      pylab.xlabel("r (kpc)")
      pylab.ylabel('Mdot (M_sun/yr)')
      pylab.legend(("Bondi Rate", "Inflow"),loc="upper left")
      pylab.savefig(stest+'_Bondi_hot.png')
      pylab.clf()
 pylab.figure(0)
 pylab.semilogy(time,Mdot_B)
 pylab.semilogy(time,Mdot_A)
 pylab.xlabel("time (code units)")
 pylab.ylabel('Mdot (M_sun/yr)')
 pylab.legend(("Bondi Rate", "Inflow"),loc="upper right")
 pylab.savefig('Bondi_time_BCG_hot.png')
 pylab.clf()
 return



def VrCs(imin=None,imax=None,rmax=10,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin + 1
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","Density"],weight="CellMassMsun")
   r=pc.plots[-1].data["Radiuskpc"]
   rho=pc.plots[-1].data["Density"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","SoundSpeed"],weight="CellMassMsun")
   Cs=pc.plots[-1].data["SoundSpeed"]   #cm/s
   pylab.figure(i)
   pylab.loglog(r,Cs/1.0e5)   #km/s
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","NegativeRadialMassFlux"])
   NeFlux=pc.plots[-1].data["NegativeRadialMassFlux"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","PositiveRadialMassFlux"])
   PoFlux=pc.plots[-1].data["PositiveRadialMassFlux"]
   PoVr=PoFlux/rho     #cm/s
   NeVr=NeFlux/rho
   for j in reversed(xrange(r.size)):
      r_s=r[j]
      if Cs[j]<-NeVr[j]:
        break
   print "sonic radius=", r_s
   pylab.loglog(r,PoVr/1.0e5)
   pylab.loglog(r,-NeVr/1.0e5)
   pylab.xlim([0.5,rmax])
   pylab.ylim([1.0,1e3])
   pylab.xlabel("r (kpc)")
   pylab.ylabel('V (km/s)')
   pylab.legend(("Sound Speed", "Outflow", "Inflow"),loc="lower right")
   pylab.savefig(stest+'_Cs.png')
   pylab.clf()
 return


def VrCs_Binned(imin=None,imax=None,rmax=40,mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin + 1
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   sp_hot=sp.cut_region(["grid['Temperature']>1.0e5"])
   prof=BinnedProfile1D(sp_hot, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)
   r=prof["Radiuskpc"]
   prof.add_fields("SoundSpeed",weight="CellMassMsun")
   Cs=prof["SoundSpeed"]    #cm/s
   sp_out=sp_hot.cut_region(["grid['RadialVelocity']>0"])
   prof_out=BinnedProfile1D(sp_out, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)
   prof_out.add_fields("RadialVelocity",weight="CellMassMsun")
   Vr_out=prof_out["RadialVelocity"]
   prof_out.add_fields("CellMassMsun", weight=None)
   M_out=prof_out["CellMassMsun"]
   sp_in=sp_hot.cut_region(["grid['RadialVelocity']<0"])
   prof_in=BinnedProfile1D(sp_in, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)
   prof_in.add_fields("RadialVelocity",weight="CellMassMsun")
   Vr_in=prof_in["RadialVelocity"]
   prof_in.add_fields("CellMassMsun", weight=None)
   M_in=prof_in["CellMassMsun"]
   pylab.figure(i)
   pylab.loglog(r,Cs/1.0e5)   #km/s
   pylab.loglog(r,Vr_out/1.0e5)
   pylab.loglog(r,-Vr_in/1.0e5)
   pylab.xlim([0.5,rmax])
   pylab.ylim([1.0,5e3])
   pylab.xlabel("r (kpc)")
   pylab.ylabel('V (km/s)')
   pylab.legend(("Sound Speed", "Outflow", "Inflow"),loc="lower right")
   pylab.savefig(stest+'_Cs_Binned.png')
   pylab.clf()



def Bondi_Radius(imin=None,imax=None,rmax=10):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin + 1
 M_BH=3.4e8
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","Density"],weight="CellMassMsun")
   r=pc.plots[-1].data["Radiuskpc"]
   rho=pc.plots[-1].data["Density"]
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","SoundSpeed"],weight="CellMassMsun")
   Cs=pc.plots[-1].data["SoundSpeed"]
   Bondi_Radius=2.0*GravConst*M_BH*SolarMass/Cs**2  #cgs
   for j in range(r.size):
      r_B=r[j]
      if r_B>(Bondi_Radius[j]/kpc):
        break
   pylab.figure(i)
   pylab.loglog(r,Bondi_Radius/kpc)
   pylab.loglog(r,r)
   pylab.axvline(x=r_B)
   pylab.xlim([0.1,rmax])
   pylab.xlabel("r (kpc)")
   pylab.ylabel('Bondi Radius (pc)')
   pylab.savefig(stest+'_BondiRadius.png')
   pylab.clf()
 return

def t_cool(imin=None,imax=None,rmin=1.0,rmax=10.0, mod=1,nbins=200):
 import glob
 if imin is None:
    imin=0
 if imax is None:
    imax=imin + 1
 cm=pylab.cm.gist_rainbow
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = load(fn)
    center=[0.5, 0.5, 0.5]
    sp1 = pf.h.sphere(center, rmin/pf['kpc'])
    sp2 = pf.h.sphere(center, rmax/pf['kpc'])
    shell = pf.h.boolean([sp2, "NOT", sp1])
    pc = PlotCollection(shell, center)
    p = pc.add_profile_object(shell, ["CoolingTime", "CellMassMsun"], weight = None, x_bins=nbins)
    t_cool=pc.plots[-1].data["CoolingTime"]
    mass_t=pc.plots[-1].data["CellMassMsun"]
#    pc.add_profile_object(shell, ["Density", "CellMassMsun"], weight = None)
    pc.add_profile_object(shell, ["Density", "CellMassMsun"], x_bins=nbins)
    Density=pc.plots[-1].data["Density"]
    mass_rho=pc.plots[-1].data["CellMassMsun"]
    pc.add_profile_object(shell, ["Temperature", "CellMassMsun"], weight = None, x_bins=nbins)
    Temperature=pc.plots[-1].data["Temperature"]
    mass_T=pc.plots[-1].data["CellMassMsun"]
    pc.add_profile_object(shell, ["Entropy", "CellMassMsun"], weight = None, x_bins=nbins)
    Entropy=pc.plots[-1].data["Entropy"]
    mass_s=pc.plots[-1].data["CellMassMsun"]
    pylab.figure(0)
    pylab.loglog(Density,mass_rho, color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('Density')
    pylab.ylabel('mass (Msun)')
    pylab.ylim([1e7,1e11])
    pylab.savefig('Density_in%ikpc.png' %rmax)
    pylab.figure(1)
    pylab.loglog(Temperature,mass_T, color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('Temperature (K)')
    pylab.ylabel('mass (Msun)')
    pylab.ylim([1e7,1e11])
    pylab.xlim([1e4,2e8])
    pylab.savefig('Temperature_in%ikpc.png' %rmax)
    pylab.figure(2)
    pylab.loglog(Entropy,mass_s, color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('Entropy')
    pylab.ylabel('mass (Msun)')
    pylab.ylim([1e7,1e11])
    pylab.savefig('Entropy_in%ikpc.png' %rmax)
    pylab.figure(3)
    pylab.loglog(t_cool,mass_t, color=cm(float(i-imin)/(imax-imin)))
    pylab.xlabel('cooling time (year)')
    pylab.ylabel('mass (Msun)')
    pylab.ylim([1e7,1e11])
    pylab.xlim([1e2,1e11])
    pylab.savefig('t_cool_in%ikpc.png' %rmax)
 pylab.clf()
 return

def timescales(imin=None,imax=None,rmax=None, mod=1, CoolingTimeOn=0):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=500
 conc=6.81 
 M_vir=8.5e14
 Rs=0.0224   #SphereCoreRadius=0.0224
 cm=pylab.cm.gist_rainbow
 min_ratio=[]
 Mdot_BH=[]
 coldgas=[] 
 time=[]
 min_tcool=[]
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   center= numpy.array([0.5, 0.5, 0.5])
   time.append(pf.current_time)
   Mdot_BH.append(pf["ClusterSMBHJetMdot"]*2.0)
   sp = pf.h.sphere(center,(rmax, 'kpc'))
   totalcoldgas=sp.quantities["TotalQuantity"](["ColdCellMassMsun"])[0]
   coldgas.append(totalcoldgas)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Temperature"],weight="CellMassMsun")
   T=pc.plots[-1].data["Temperature"]     #yr
   r=pc.plots[-1].data["Radiuskpc"]
   if CoolingTimeOn==0:
      pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "CoolingTime"],weight="CellMassMsun")
      t_cool=pc.plots[-1].data["CoolingTime"]     #yr
   else:    ## CoolingTime on
      pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "Cooling_Time"],weight="CellMassMsun")
      t_cool=pc.plots[-1].data["Cooling_Time"]     #s
      t_cool/=3.15e7     #yr
   pc.add_profile_sphere(rmax, 'kpc', ["Radiuskpc","SoundSpeed"],weight="CellMassMsun")
   Cs=pc.plots[-1].data["SoundSpeed"]
   x1=r/(Rs*pf['kpc'])
   M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass #DM+BCG + BH , cgs
   t_dyn=numpy.pi*(r*kpc)**1.5/(2.0*(GravConst*M_dyn)**0.5)/3.15e7   #from s to yr
   t_sound=(r*kpc/Cs)/3.15e7   #cgs to yr
   dog= numpy.genfromtxt("cool_rates.in")
   [logT,solar,half_solar]=[dog[:,0],dog[:,1],dog[:,2]]
   T_std=10.0**logT
   fT_std=numpy.zeros(T_std.size)
   for j in range(logT.size-1):
       fT_std[j]=2.0-(T_std[j]/half_solar[j])*(half_solar[j+1]-half_solar[j])/(T_std[j+1]-T_std[j])
   fT_std[-1]=fT_std[-2]
   fT=numpy.interp(T,T_std, fT_std)
   t_TI=(5.0/3.0)*t_cool/fT
   t_ff=8.0**0.5*t_dyn/numpy.pi
   cool_dyn=t_cool/t_dyn
   good=(r>5.0)&(r<rmax-1)
   min_ratio.append(min(cool_dyn[good]))
   min_tcool.append(min(t_cool[good]))
   TI_ff=t_TI/t_ff
   pylab.figure(i)
   pylab.loglog(r,t_cool)
   pylab.loglog(r,t_dyn)
   pylab.loglog(r,t_sound)
   pylab.xlabel('r (kpc)')
   pylab.ylabel('timescales (yr)')
   pylab.savefig(stest+'t_cool_and_t_dyn_and_t_sound.png')
   pylab.figure(1)
   pylab.loglog(r,cool_dyn, color=cm(float(i-imin)/(imax-imin)))
   pylab.xlim(0.08,rmax)
   pylab.ylim(0.9,200)
   pylab.xlabel('r (kpc)')
   pylab.ylabel(r'$t_{cool}/t_{dyn}$')
   pylab.savefig('t_cool_vs_t_dyn_PC.png')
   pylab.figure(2)
   pylab.loglog(r,TI_ff, color=cm(float(i-imin)/(imax-imin)))
   pylab.xlim(0.08,rmax)
   pylab.ylim(0.9,200)
   pylab.xlabel('r (kpc)')
   pylab.ylabel(r'$t_{TI}/t_{ff}$')
   pylab.savefig('t_TI_vs_t_ff_PC.png')
 time[:] = [x*TimeUnits/3.16e16 for x in time]  # from code unit to Gyr
 numpy.savez("min_ratio_etc", min_ratio=min_ratio, coldgas=coldgas,timeGyr=time, min_tcool=min_tcool,Mdot_BH=Mdot_BH)
 pylab.clf()
 return


def ReadCoolingCurve():
    FileIn = open(homedir+"/code/cool_rates.in")
    data = na.loadtxt(FileIn)
    FileIn.close()
    return (data[:,0],data[:,2])

def _CoolingRate_perV(field, data):   #like coolingtime #this is the rate per volume
    (LogT, LogCoolRate) = ReadCoolingCurve()
    return 10.0**na.interp(na.log10(data["Temperature"]), LogT, LogCoolRate)
add_field("CoolingRate_perV", function=_CoolingRate_perV, units = r"ergs s^{-1}cm^{-3}")  

def _FakeLx(field,data):
#    return data["XRayEmissivity"]*data["CellVolume"] 
    #return (data["Density"]**2.0)*data["Temperature"]**0.5*data["CellVolume"]
    return (data["Density"]**2.0)*data["CellVolume"]
add_field("FakeLx", function=_FakeLx, units = r"ergs s^{-1}")

def timescales_Binned(imin=None,imax=None,rmin=5, rmax=500, mod=1, case=1, CoolingTimeOn=0, Particle=0, step=9, Trajectory=0,y_min=6.9, y_max=9,OnlyPhase=0):
##use BinnedProfile1D instead of PC, also cut out the jet (Vr>0)
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 conc=6.81
 M_vir=8.5e14
 Rs=0.0224   #SphereCoreRadius=0.0224
 m_p=1.67e-24
 ks=1.84e-5/37.0
 fc=1.0/3.0
 (LogT, LogCoolRate) = ReadCoolingCurve()
 if Particle==0:
   cm=pylab.cm.gist_rainbow
 else:
   cm=pylab.cm.spring
 pylab.rc('text', fontsize=20)
 line = ["-","--",":","-."]
 from yt.analysis_modules.spectral_integrator.api import add_xray_emissivity_field
 add_xray_emissivity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
 cc=numpy.load(homedir+"/code/conduction_tcool.npz")
 t_cond_plot=cc["CoolingTime"]   #yr
 r_plot=cc["r_plot"]
 loop=range(imin,imax,mod)
# loop.insert(0,0)
 for i in loop:
     stest = 'stest_%04i' %i
     fn = name_pattern %(i,i)
     pf = load(fn)
     sp1 = pf.h.sphere('c', rmin/pf['kpc'])
     sp2 = pf.h.sphere('c', rmax/pf['kpc'])
     shell = pf.h.boolean([sp2, "NOT", sp1])
     hot_shell=shell.cut_region(["grid['Temperature'] > 1.0e5"])   #warmer than 0.1 keV-->not filaments yet
#     prof2D=BinnedProfile2D(data_source=hot_shell,x_n_bins=128,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=128,y_bin_field="CoolingTime",y_lower_bound=1e5,y_upper_bound=1e11,y_log=True,lazy_reader=True)
     prof2D=BinnedProfile2D(data_source=hot_shell,x_n_bins=128,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=128,y_bin_field="CoolingTime",y_lower_bound=10.0**y_min,y_upper_bound=10.0**y_max,y_log=True,lazy_reader=True)
#     prof2D=BinnedProfile2D(data_source=hot_shell,x_n_bins=128,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=128,y_bin_field="Cooling_Time",y_lower_bound=10.0**y_min,y_upper_bound=10.0**y_max,y_log=True,lazy_reader=True)
     prof2D.add_fields("CellMassMsun", weight=None)
     mm=prof2D["CellMassMsun"]
     logmm=numpy.log10(mm)
     TT=logmm.T
     pylab.figure(i+1000)
     pylab.imshow(TT,interpolation='nearest',extent=[rmin,rmax,y_min,y_max],aspect='auto',origin='lower',cmap='Blues')
     sp=pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
#     sp_in=sp.cut_region(["grid['RadialVelocity']<0"])
#     sp_in=sp.cut_region(["grid['Theta']>0.6"])
#     sp_in=sp.cut_region(["grid['Theta']<0.6"])
#     sp_in=sp
     sp_in=sp.cut_region(["grid['Temperature'] > 1.0e6"])   #warmer than 0.1 keV-->not filaments yet
     #prof=BinnedProfile1D(sp_in, 128, "Radiuskpc", 1e-2,rmax,lazy_reader=True)    
     prof=BinnedProfile1D(sp_in, 28, "Radiuskpc", 0.8,rmax,lazy_reader=True)   ### try larger bins for obs !!!!!!!!!
     r=prof["Radiuskpc"]
     #prof.add_fields("Temperature",weight="CellMassMsun")
     prof.add_fields("Temperature",weight="Xray_Emissivity_0.5_9.9keV")
     prof.add_fields("Pressure",weight="CellMassMsun")
     T=numpy.array(prof["Temperature"])
     Pressure=prof["Pressure"]
     if CoolingTimeOn==0:
        prof.add_fields("CoolingTime",weight="CellMassMsun")
        t_cool=prof["CoolingTime"]   #yr
     else:
        prof.add_fields("Cooling_Time",weight="Xray_Emissivity_0.5_9.9keV")
        t_cool=prof["Cooling_Time"]   #yr
        t_cool/=3.15e7     #yr
# here we compute t_cool like the observers
     prof.add_fields("CellVolume", weight=None)  #volume in each shell
     Vshell=numpy.array(prof["CellVolume"])
     prof.add_fields("FakeLx",weight=None)
     FakeLx=numpy.array(prof["FakeLx"])
     #rho_obs=(FakeLx/(Vshell*T**0.5))**0.5
     rho_obs=(FakeLx/Vshell)**0.5
     print "rho_obs", rho_obs
     CoolingRate_obs = 10.0**numpy.interp(numpy.log10(T), LogT, LogCoolRate)
     t_cool_obs=6.0 * kboltz * T / ((rho_obs/(0.59*m_p)) * CoolingRate_obs) / 3.15e7
     print "t_cool_obs", t_cool_obs
     x1=r/(Rs*pf['kpc'])
     M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass #DM + BCG + BH , cgs
     t_dyn=numpy.pi*(r*kpc)**1.5/(2.0*(GravConst*M_dyn)**0.5)/3.15e7   #from s to yr
     dog= numpy.genfromtxt("cool_rates.in")
     [logT,solar,half_solar]=[dog[:,0],dog[:,1],dog[:,2]]
     T_std=10.0**logT
     fT_std=numpy.zeros(T_std.size)
     for j in range(logT.size-1):
       fT_std[j]=2.0-(T_std[j]/half_solar[j])*(half_solar[j+1]-half_solar[j])/(T_std[j+1]-T_std[j])
     fT_std[-1]=fT_std[-2]
     fT=numpy.interp(T,T_std, fT_std)
     t_TI=(5.0/3.0)*t_cool/fT
     t_ff=8.0**0.5*t_dyn/numpy.pi
     cool_dyn=t_cool/t_dyn
     TI_ff=t_TI/t_ff
     pylab.plot(r,numpy.log10(t_dyn),color='green',linewidth=2,linestyle='--')
     pylab.ylabel("log "+r"$t_{cool}$"+"(yr)")
     pylab.xlabel("r (kpc)")
     pylab.xlim([rmin,rmax])
     pylab.ylim([y_min,y_max])
#     t_cond=2.0*3.04*(5.0/3.0/(5.0/3.0-1.0))*Pressure/((2.0*3.14/(r*kpc))**2.0*fc*ks*T**3.5)/3.15e7
     if Particle:
        aa=numpy.load("cool_particles.npz")
        rho_p=aa['rho_plot']
        T_p=aa['T_plot']
        r_p=aa['r_plot']
        (LogT, LogCoolRate) = ReadCoolingCurve()
        CoolingRate = 10.0**numpy.interp(numpy.log10(T_p), LogT, LogCoolRate)
        t_cool_p=5.0 * kboltz * T_p / ((rho_p/(0.59*m_p)) * CoolingRate) / 3.15e7
        H, xedges, yedges = numpy.histogram2d(numpy.log10(t_cool_p[:,step]), r_p[:,step], bins=100,range=([y_min,y_max],[rmin,rmax]))
        X=numpy.linspace(rmin,rmax,num=100)
        Y=numpy.linspace(y_min,y_max,num=100)
        pylab.contour(X,Y,H,4)
        pylab.savefig(stest+"_2D_timescales_step%i.png" %step)
     elif Trajectory:
        aa=numpy.load("cool_particles.npz")
        rho_p=aa['rho_plot']
        T_p=aa['T_plot']
        r_p=aa['r_plot']
        (LogT, LogCoolRate) = ReadCoolingCurve()
        CoolingRate = 10.0**numpy.interp(numpy.log10(T_p), LogT, LogCoolRate)
        t_cool_p=5.0 * kboltz * T_p / ((rho_p/(0.59*m_p)) * CoolingRate) / 3.15e7
        total_steps=len(T_p[0,:])  #total number of steps
        #for pp in [10110, 10910, 11110, 2910]:
        #for pp in [10110, 480, 11110, 2910]
        for pp in [9910, 990, 970,9410]:
#            pylab.plot(r_p[pp,:],numpy.log10(t_cool_p[pp,:]),color="blue")
            for qq in range(total_steps-6,1,-1):
                pylab.plot(r_p[pp,0:qq],numpy.log10(t_cool_p[pp,0:qq]),color=cm(float(qq)/(total_steps-7)),linewidth=2)
        pylab.tight_layout()
        pylab.savefig(stest+"_2D_timescales_trajectory.png")
     elif OnlyPhase:
        pylab.savefig(stest+"_2D_timescales.png")
     else:
        pylab.clf()
     pylab.clf()
     pylab.figure(i)
     pylab.loglog(r,t_cool/1.0e6)
     pylab.loglog(r,t_dyn/1.0e6)
     pylab.xlim([0.4,rmax])
     pylab.ylim([0.5,1.5e5])
     pylab.xlabel('r (kpc)')
     pylab.ylabel('Time scales (Myr)')
     pylab.legend((r"$t_{cool}$",r"$t_{dyn}$"),loc="lower right")
     pylab.tight_layout()
#     pylab.savefig(stest+'t_cool_and_t_dyn.png')
     pylab.clf()
     pylab.figure(i+1000)
     pylab.loglog(r,t_TI/1.0e6)   #Myr
     pylab.loglog(r,t_ff/1.0e6)
     pylab.xlim([0.4,rmax])
     pylab.ylim([0.5,1.5e5])
     pylab.xlabel('r (kpc)')
     pylab.ylabel('Time scales (Myr)')
     pylab.legend((r"$t_{TI}$",r"$t_{ff}$"),loc="lower right")
     pylab.tight_layout()
#     pylab.savefig(stest+'t_TI_and_t_ff.png')
     pylab.clf()
     fig2000=pylab.figure(2000)
     if case ==1:
        pylab.loglog(r,cool_dyn, color=cm(float(i-imin)/(imax-imin)))
     elif case ==2:
        pylab.loglog(r,TI_ff, color=cm(float(i-imin)/(imax-imin)))
     elif case == 3:
        lower_bound=(t_ff/1.0e6)*5.0
        upper_bound=(t_ff/1.0e6)*20.0
        pylab.loglog(r,(t_ff/1.0e6)*10.0,color="black")
        #pylab.loglog(r,(t_cond/1.0e6),color="black",linestyle="--") ##conduction
        pylab.loglog(r_plot,(t_cond_plot/1.0e6),color="black",linestyle="--") ##conduction
        pylab.loglog(r,t_cool/1.0e6,color=cm(float(i-imin)/(imax-imin)))
        pylab.fill_between(r,lower_bound,upper_bound,alpha=0.3,facecolor="grey")
     elif case == 4:   # calculate t_cool like an observer 
        lower_bound=(t_ff/1.0e6)*5.0
        upper_bound=(t_ff/1.0e6)*20.0
        pylab.loglog(r,(t_ff/1.0e6)*10.0,color="black")
        #pylab.loglog(r,(t_cond/1.0e6),color="black",linestyle="--") ##conduction
        pylab.loglog(r_plot,(t_cond_plot/1.0e6),color="black",linestyle="--") ##conduction
        pylab.loglog(r,t_cool_obs/1.0e6,color=cm(float(i-imin)/(imax-imin)))
        pylab.fill_between(r,lower_bound,upper_bound,alpha=0.3,facecolor="grey")
     elif case == 5:   #save each figure to make a movie
        fig2000.set_figheight(6)
        fig2000.set_figwidth(6)
        lower_bound=(t_ff/1.0e6)*5.0
        upper_bound=(t_ff/1.0e6)*20.0
        pylab.loglog(r,(t_ff/1.0e6)*10.0,color="black")
        #pylab.loglog(r,(t_cond/1.0e6),color="black",linestyle="--") ##conduction
        pylab.loglog(r_plot,(t_cond_plot/1.0e6),color="black",linestyle="--") ##conduction
        pylab.loglog(r,t_cool/1.0e6,color=cm(float(i-imin)/(imax-imin)))
        pylab.fill_between(r,lower_bound,upper_bound,alpha=0.3,facecolor="grey")
        pylab.xlim([1,rmax])
        pylab.ylim([20,1.5e5])
        pylab.xlabel('r (kpc)')
        pylab.ylabel('Cooling Time (Myr)')
#     pylab.legend((r'$10 t_{dyn}$'),loc='lower right')
        pylab.tight_layout()
        pylab.savefig(stest+'Mark_t_cool_obs_tff.png')
#        pylab.clf()
 if case == 1:
     pylab.xlabel('r (kpc)')
     pylab.xlim(0.4,rmax)
     pylab.ylim(0.9,200)
     pylab.ylabel('t_cool/t_dyn')
     pylab.savefig('t_cool_vs_t_dyn2.png')
 elif case ==2:
     pylab.xlabel('r (kpc)')
     pylab.xlim(0.4,rmax)
     pylab.ylim(0.9,200)
     pylab.ylabel('t_TI/t_ff')
     pylab.savefig('t_TI_vs_t_ff2.png')
 elif case == 3:
     pylab.xlim([1,rmax])
     pylab.ylim([20,1.5e5])
     pylab.xlabel('r (kpc)')
     pylab.ylabel('Cooling Time (Myr)')
#     pylab.legend((r'$10 t_{dyn}$'),loc='lower right')
     pylab.tight_layout()
     pylab.savefig('Mark_t_cool_tff.png')
 elif case == 4:
     pylab.xlim([1,rmax])
     pylab.ylim([20,1.5e5])
     pylab.xlabel('r (kpc)')
     pylab.ylabel('Cooling Time (Myr)')
#     pylab.legend((r'$10 t_{dyn}$'),loc='lower right')
     pylab.tight_layout()
     pylab.savefig('Mark_t_cool_obs_tff.png')
 pylab.clf()
 return

def Mark_K(imin=0,imax=None,rmin=1,rmax=100,mod=1,y_min=1,y_max=100):
 kev_to_erg=6.2415e8
 if imax is None:
    imax=imin+1
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = load(fn)
    sp1 = pf.h.sphere('c', rmin/pf['kpc'])
    sp2 = pf.h.sphere('c', rmax/pf['kpc'])
    shell = pf.h.boolean([sp2, "NOT", sp1])
    hot_shell=shell.cut_region(["grid['Temperature'] > 1.0e5"])   #warmer than 0.1 keV-->not filaments yet
    #prof2D=BinnedProfile2D(data_source=hot_shell,x_n_bins=128,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=128,y_bin_field="Entropy",y_lower_bound=10.0**y_min,y_upper_bound=10.0**y_max,y_log=True,lazy_reader=True)
    prof2D=BinnedProfile2D(data_source=hot_shell,x_n_bins=128,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=128,y_bin_field="Entropy",y_lower_bound=y_min/kev_to_erg,y_upper_bound=y_max/kev_to_erg,y_log=True,lazy_reader=True)
    prof2D.add_fields("CellMassMsun", weight=None)
    mm=prof2D["CellMassMsun"]
    logmm=numpy.log10(mm)
    TT=logmm.T
    pylab.figure(i+1000)
    #pylab.imshow(TT,interpolation='nearest',extent=[rmin,rmax,y_min,y_max],aspect='auto',origin='lower',cmap='Oranges') 
    pylab.imshow(TT,interpolation='nearest',extent=[rmin,rmax,numpy.log10(y_min),numpy.log10(y_max)],aspect='auto',origin='lower',cmap='Oranges') 
    #pylab.imshow(TT,interpolation='nearest',aspect='auto',origin='lower',cmap='Blues') 
    pylab.xlim([rmin,rmax])
    #pylab.ylim([y_min,y_max])
    pylab.ylim([numpy.log10(y_min),numpy.log10(y_max)])
    y_locs=[0,1,2]
    y_labels=[0,1,2]
    pylab.yticks(y_locs, y_labels)
    pylab.xlabel('r (kpc)')
    #pylab.ylabel("Entropy K (keV "+r"$cm^2$"+")")
    pylab.ylabel("log K (keV "+r"$\rm{cm}^2$"+")")
    pylab.savefig(stest+'Mark_K_2D.png')
 pylab.clf()
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
    pf = load(fn)
    sp1 = pf.h.sphere('c', rmin/pf['kpc'])
    sp2 = pf.h.sphere('c', rmax/pf['kpc'])
    shell = pf.h.boolean([sp2, "NOT", sp1])
    hot_shell=shell.cut_region(["grid['Temperature'] > 1.0e5"])   #warmer than 0.1 keV-->not filaments yet
    #prof2D=BinnedProfile2D(data_source=hot_shell,x_n_bins=128,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=128,y_bin_field="Entropy",y_lower_bound=10.0**y_min,y_upper_bound=10.0**y_max,y_log=True,lazy_reader=True)
    prof2D=BinnedProfile2D(data_source=hot_shell,x_n_bins=128,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=True,y_n_bins=128,y_bin_field="Entropy",y_lower_bound=y_min/kev_to_erg,y_upper_bound=y_max/kev_to_erg,y_log=True,lazy_reader=True)
    prof2D.add_fields("CellMassMsun", weight=None)
    mm=prof2D["CellMassMsun"]
    logmm=numpy.log10(mm)
    TT=logmm.T
##begin percentage calculation
    mm_cumsum=numpy.cumsum(mm,axis=1)
    mm_percentage=numpy.zeros(mm_cumsum.shape)
    percentiles=numpy.array([0.2,0.5,0.8])
    lines=numpy.empty([128,3])
    y_locations=numpy.linspace(numpy.log10(y_min),numpy.log10(y_max),num=128)
    print "y_locations", y_locations
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
#    pylab.text(550,480,"t = %.2f Gyr"%time,fontsize=18)
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


def Shocks(imin=None,imax=None,rmax=None, thetas=None, phis=None):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=40
 if thetas is None:
    thetas=numpy.array([0, 1.0/6.0, 1.0/3.0])*numpy.pi
 if phis is None:
    phis=numpy.array([0, 1.0/2.0])*numpy.pi
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   center = (0.5,0.5,0.5)
   j=0
   for theta in thetas:
    for phi in phis:
      end_x=0.5+rmax*numpy.sin(theta)*numpy.cos(phi)/pf['kpc']
      end_y=0.5+rmax*numpy.sin(theta)*numpy.sin(phi)/pf['kpc']
      end_z=0.5+rmax*numpy.cos(theta)/pf['kpc']
      ray=pf.h.ray(center,(end_x,end_y,end_z))
      rho=ray["Density"]
      Vr=ray["RadialVelocity"]
      T=ray["Temperature"]
      P=ray["Pressure"]
      r=ray['t']*rmax
      Mach=ray["RadialMachNumber"]
      Entropy=ray["Entropy"]
      j=j+1
      pylab.figure(j)
      pylab.subplot(2,2,1)
      pylab.semilogy(r,rho)
      pylab.xlabel('r (kpc)')
      pylab.ylabel(r'$rho (g cm^{-3})$')
      pylab.subplot(2,2,2)
      pylab.semilogy(r,T)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Temperature (K)')
      pylab.subplot(2,2,3)
      pylab.semilogy(r,P)
      pylab.xlabel('r (kpc)')
      pylab.ylabel(r'$Pressure (dyne cm^{-2})$')
      pylab.subplot(2,2,4)
      pylab.semilogy(r,Entropy)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Entropy')
      pylab.tight_layout()
      pylab.savefig(stest+'shock_%.2fpi' %(theta/numpy.pi)+ '_%.2fpi.png' %(phi/numpy.pi))  # add name
      pylab.clf()
 return

def Shocks_paper(i=100, rmin=10,rmax=25,theta=1.0,phi=1.0):
   pylab.rc('axes', linewidth=2,labelsize=12,labelweight='bold')
   pylab.rc('lines', linewidth=2)
   font = {'weight' : 'bold', 'size'   : 14}
   pylab.rc('font', **font)
   pylab.rc('text', usetex=True)
   pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   center = (0.5,0.5,0.5)
   end_x=0.5+rmax*numpy.sin(theta)*numpy.cos(phi)/pf['kpc']
   end_y=0.5+rmax*numpy.sin(theta)*numpy.sin(phi)/pf['kpc']
   end_z=0.5+rmax*numpy.cos(theta)/pf['kpc']
   ray=pf.h.ray(center,(end_x,end_y,end_z))
   rho=ray["Density"]
   Vr=ray["RadialVelocity"]
   T=ray["Temperature"]
   P=ray["Pressure"]
   r=ray['t']*rmax
   Entropy=ray["Entropy"]
   f, (ax1, ax2, ax3) = pylab.subplots(3, sharex=True)
   f.set_size_inches(6,8)
   ax1.plot(r,rho)
   d1=1.3e-25
   d2=1.7e-25
   ax1.set_ylim(3e-26,d2)
   ax1.plot([[14,25,33],[14,25,33]],[[d1,d1,d1],[d2,d2,d2]],'k',color='black')
   ax1.set_yticks((0.5e-25,1.0e-25,1.5e-25))
   ax1.set_yticklabels(labels=('0.5','1.0','1.5'))
#   ax1.set_xlabel('r (kpc)')
   ax1.set_ylabel('Density '+r'$(10^{-25}\ \rm{g}\ \rm{cm}^{-3})$')
   ax2.plot(r,P)
   ax2.set_ylim(2e-10,8.8e-10)
#   ax2.set_xlabel('r (kpc)')
   ax2.set_ylabel('Pressure '+r'$(10^{-10}\ \rm{dyne}\ \rm{cm}^{-2})$')
   ax2.set_yticks((3e-10,5e-10,7e-10))
   ax2.set_yticklabels(labels=('3','5','7'))
   ax3.plot(r,Entropy)
   ax3.set_ylim(2e-8,8e-8)
   ax3.set_yticks((3e-8,5e-8,7e-8))
   ax3.set_yticklabels(labels=('3','5','7'))
   ax3.set_xlabel('r (kpc)')
   ax3.set_ylabel('Entropy ' + r'$(10^{-8}\ \rm{ergs}\ \rm{cm}^{3\gamma-3})$')
   ax1.set_xlim(rmin,rmax)
   ax2.set_xlim(rmin,rmax)
   ax3.set_xlim(rmin,rmax)
   pylab.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
   f.tight_layout()
   f.subplots_adjust(hspace=0)
   pylab.savefig(stest+'shock_%.2fpi' %(theta/numpy.pi)+ '_%.2fpi_paper.png' %(phi/numpy.pi))  # add name
   pylab.clf()
   return

def Shocks_Mach(imin=None,imax=None,rmax=None, theta=0, phi=0.5, rshock=22.0):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=40
 gamma=5.0/3.0
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   center = (0.5,0.5,0.5)
   print "pf['kpc']=", pf['kpc']
   end_x=0.5+rmax*numpy.sin(theta)*numpy.cos(phi)/pf['kpc']
   end_y=0.5+rmax*numpy.sin(theta)*numpy.sin(phi)/pf['kpc']
   end_z=0.5+rmax*numpy.cos(theta)/pf['kpc']
   ray=pf.h.ray(center,(end_x,end_y,end_z))
   rho=ray["Density"]
   T=ray["Temperature"]
   P=ray["Pressure"]
   r=ray['t']*rmax
   shock=(r > rshock-1.0) & (r < rshock+1.0)
   rho2=numpy.amax(rho[shock])
   rho1=numpy.amin(rho[shock])
   P2=numpy.amax(P[shock])
   P1=numpy.amin(P[shock])
   print "rho.size=", rho.size 
   M_rho=(2.0/(rho1*(gamma+1.0)/rho2-(gamma-1.0)))**0.5
   M_P=((P2*(gamma+1.0)/P1+gamma-1.0)/(2.0*gamma))**0.5
   print "M_rho= ", M_rho, " and M_P= ", M_P
 return


def _Xray_Emission(field, data) :

    if data.has_field_parameter("Z") :
        Z = data.get_field_parameter("Z")
    else :
        Z = 1.077 # Primordial H/He plasma
        
    if data.has_field_parameter("mue") :
        mue = data.get_field_parameter("mue")
    else :
        mue = 1./0.875 # Primordial H/He plasma

    if data.has_field_parameter("mui") :
        mui = data.get_field_parameter("mui")
    else :
        mui = 1./0.8125 # Primordial H/He plasma

    n_e = data["Density"]/(mue*mp)
    n_i = data["Density"]/(mui*mp)

    # Assume a constant Gaunt factor g_ff=1.3
    g_ff = 1.3

    #eq(3-56) Spitzer, Free-free emission integrated over all frequency
    L_bol=1.426e-27*Z*Z*n_e*n_i*numpy.sqrt(data["Temperature"])*g_ff   #ergs cm^-3 s^-1
    return L_bol

add_field("Xray_Emission", function=_Xray_Emission)

KtokeV = 8.617e-08 # Convert degrees Kelvin to degrees keV
def _Xray_Emission_Corrected(field, data) :   # correct for metals?
    
    kT = data["Temperature"]*KtokeV

#    if data.has_field_parameter("Z") :
#        Z = data.get_field_parameter("Z")
#    else :
#        Z = 1.077 # Primordial H/He plasma

    if data.has_field_parameter("metallicity") :
        metallicity = data.get_field_parameter("metallicity")
    else :
        metallicity = 0.5    # Perseus is half solar metallicity
    if data.has_field_parameter("hniu1") :
        hniu1 = data.get_field_parameter("hniu1")
    else :
        hniu1 = 0.5    # 0.5 keV
    if data.has_field_parameter("hniu2") :
        hniu2 = data.get_field_parameter("hniu2")
    else :
        hniu2 = 3.0    # 3 keV

    #mutiply L_bol by the correction function (eq. 21 of Bryan and Norman 1998)
    eps= na.zeros(kT.shape)
    eps[kT>2.]=(numpy.exp(-hniu1/kT)-numpy.exp(-hniu2/kT))
    eps[kT<2.]=(numpy.exp(-hniu1/kT)-numpy.exp(-hniu2/kT))*(kT/2.0)**(0.9*numpy.sqrt(metallicity/0.3))   #gamma is 0.9 for 0.5-2.4 keV??! what is this?!!    Is metallicity the same as Z?? Not!?
    L_metal=data["Xray_Emission"]*eps*data["CellVolume"]   # already multiplied by CellVolume!
    return L_metal
add_field("Xray_Emission_Corrected", function=_Xray_Emission_Corrected)



def Xray_Greg(imin=None,imax=None,rmax=None):
 import numpy as numpy
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=numpy.array([50,500])
 time=numpy.zeros(imax-imin)
 Xray=numpy.zeros(((imax-imin),rmax.size))
 pylab.figure(0)
 for j in range(rmax.size):
   for i in range(imin, imax):
     stest = 'stest_%04i' %i
     ii=i-imin
     fn = name_pattern %(i,i)
     pf = load(fn)
#     print "Xray_Luminosity_Greg=", pf.h.all_data().quantities['TotalQuantity']('Xray_Luminosity_Greg')
     time[ii] = pf.current_time
     sphere = pf.h.sphere([0.5, 0.5, 0.5], rmax[j]/pf['kpc'])
     Xray_Emission_Corrected= sphere["Xray_Emission_Corrected"]
     print "Xray_Emission_Corrected=", Xray_Emission_Corrected
     Xray[ii,j]=numpy.sum(Xray_Emission_Corrected)
#     print "Xray=", 
#     Xray[ii,j]=sphere.quantities["TotalQuantity"](["Xray_Luminosity_Greg"])[0]
   print "Xray", Xray[ii,j]
   pylab.semilogy(time,Xray[:,j])
 pylab.xlabel('Time (code unit)')
 pylab.ylabel('FreeFree_Luminosity')
 pylab.savefig(stest+'Xray_Greg.png')
 return


def Extrema(imin=None,imax=None, rmax=None):
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=10
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   center = (0.5,0.5,0.5)
   sp=pf.h.sphere(center, rmax/pf['kpc'])
   (rho_min,rho_max),=sp.quantities["Extrema"]("Density")
   print stest,": rho_min=", rho_min, " rho_max= ",rho_max
 return


def Vr_tcool(imin=None,imax=None, rmin=None, rmax=None, Tempcut=1.0e5, mod=1):
 if imin is None:
    imin=0
 if imax is None:
    imax=1
 if rmin is None:
    rmin=5.0
 if rmax is None:
    rmax=7.0
 for i in range(imin,imax, mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   center=[0.5, 0.5, 0.5]
   sp1 = pf.h.sphere(center, rmin/pf['kpc'])
   sp2 = pf.h.sphere(center, rmax/pf['kpc'])
   shell = pf.h.boolean([sp2, "NOT", sp1])
   hot_shell=shell.cut_region(["grid['Temperature'] > 1.0e5"])   #warmer than 0.1 keV-->not filaments yet
   pc = PlotCollection(hot_shell, center=[0.5, 0.5, 0.5])
   #pc.add_phase_object(hot_shell,["CoolingTime","RadialVelocity","CellMassMsun"],weight=None,y_log=False,x_bounds=[1e5,1e11], y_bounds=[-1e8,6e8])
   pc.add_phase_object(hot_shell,["Radiuskpc","CoolingTime","CellMassMsun"],weight=None,x_log=False,y_bounds=[1e5,1e11], x_bounds=[rmin,rmax])
   pc.save(stest)
 return


def particle_ray(imin=424, rmax=0.2, nTP=None, deadcat=40,ColdTemperature=1.0e5):
 import h5py
# lu = 4.93708012753e+25   # length unit
 xpos = numpy.zeros(nTP)
 ypos = numpy.zeros(nTP)
 zpos = numpy.zeros(nTP)
 T = numpy.zeros(nTP)
 f = h5py.File("TracerParticles_%08i.h5" % (imin-188), "r")
 time = f["/"].attrs["time"]
 xpos[:],ypos[:],zpos[:]=f["position_x"],f["position_y"],f["position_z"]
 T=f["temperature"]
 cold=(T[:] < ColdTemperature)
# if (deadcat>0):   #we do not want the particles that have always been cold
#    for cat in range(deadcat):
#        cold = cold & (T[cat] > HotTemperature)     ### if we require the particles to be "hot" at first
 xpos_cold=xpos[cold]
 ypos_cold=ypos[cold]
 zpos_cold=zpos[cold]
 n_cold=xpos_cold.size
 print "n_cold=", n_cold
 fn = name_pattern %(imin,imin)
 stest = 'stest_%04i' %imin
 pf = load(fn)
 import random as random
 j=0
 for j in range(10):
   i=random.randrange(0,n_cold-1)
   center = (xpos_cold[i],ypos_cold[i],zpos_cold[i])
   phi=random.random()*numpy.pi  #fix phi, random theta
   k=0
   for k in range(3):
      theta=random.random()* numpy.pi/2.0
      print "i, theta, phi=", i, theta, phi
      end_x=xpos_cold[i]+rmax*numpy.sin(theta)*numpy.cos(phi)/pf['kpc']
      end_y=ypos_cold[i]+rmax*numpy.sin(theta)*numpy.sin(phi)/pf['kpc']
      end_z=zpos_cold[i]+rmax*numpy.cos(theta)/pf['kpc']
      ray=pf.h.ray(center,(end_x,end_y,end_z))
      rho=ray["Density"]
      Vr=ray["RadialVelocity"]
      T=ray["Temperature"]
      P=ray["Pressure"]
      r=ray['t']*rmax
      Entropy=ray["Entropy"]
      pylab.figure(j)
      pylab.subplot(2,2,1)
      pylab.semilogy(r,rho)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('rho')
      pylab.subplot(2,2,2)
      pylab.semilogy(r,T)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Temperature')
      pylab.subplot(2,2,3)
      pylab.semilogy(r,P)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Pressure')
      pylab.subplot(2,2,4)
      pylab.semilogy(r,Entropy)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Entropy')
      pylab.savefig(stest+'ray_%.2fpi' %(theta/numpy.pi)+ '_%.2fpi.png' %(phi/numpy.pi))  # add name
      pylab.clf()
   j=j+1
 return

def cool_ray(cool_list=[],dr=1.0):
 all=numpy.load("cool_particles.npz")
 r=all["r_plot"]
 n=len(r[0,:])
 Vr=all["Vr_plot"]
 rho=all["rho_plot"]
 T=all["T_plot"]
 P=all["P_plot"]
 Vx=all["Vx_plot"]
 Vy=all["Vy_plot"]
 Vz=all["Vz_plot"]
 xpos=all["xpos_plot"]
 ypos=all["ypos_plot"]
 zpos=all["ypos_plot"]
 #dr=1.0  #kpc
 #cool_list=[1,2]
 read_step=numpy.genfromtxt("cool_step",skip_footer = 1)
 cool_step=numpy.ndarray.flatten(read_step)
 for j in cool_list:
      imin=int(38+cool_step[j])   ##-1??
      fn=name_pattern %(imin,imin)
      pf=load(fn)
      p_center=numpy.array([xpos[j], ypos[j], zpos[j]])
      begin_point=p_center-dr*numpy.array([Vx[j],Vy[j],Vz[j]])/(Vx[j]**2.0+Vy[j]**2.0+Vz[j]**2.0)**0.5/pf['kpc']
      end_point=p_center+dr*numpy.array([Vx[j],Vy[j],Vz[j]])/(Vx[j]**2.0+Vy[j]**2.0+Vz[j]**2.0)**0.5/pf['kpc']
      ray=pf.h.ray(begin_point,end_point)
      rho=ray["Density"]
      Vr=ray["RadialVelocity"]
      T=ray["Temperature"]
      P=ray["Pressure"]
      r=ray['t']*dr*2.0
      Mach=ray["RadialMachNumber"]
      Entropy=ray["Entropy"]
      pylab.figure(j)
      pylab.subplot(2,2,1)
      pylab.semilogy(r,rho)
      pylab.xlabel('r (kpc)')
      pylab.ylabel(r'$rho (g cm^{-3})$')
      pylab.subplot(2,2,2)
      pylab.semilogy(r,T)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Temperature (K)')
      pylab.subplot(2,2,3)
      pylab.semilogy(r,P)
      pylab.xlabel('r (kpc)')
      pylab.ylabel(r'$Pressure (dyne cm^{-2})$')
      pylab.subplot(2,2,4)
      pylab.semilogy(r,Entropy)
      pylab.xlabel('r (kpc)')
      pylab.ylabel('Entropy')
      pylab.tight_layout()
      pylab.savefig('cool_ray_%i.png' %j)  # add name
      pylab.clf()
 return



def _CoolingRate(field,data):    # this is the rate per cell
   return data["ThermalEnergy"]*data["CellMass"]/(data["CoolingTime"]*3.16e7)
add_field("CoolingRate", units=r"\rm{ergs s^{-1}", function=_CoolingRate)

def Cooling_Heating(imin=None,imax=None,rmax=100, mod=1, heating=True,log=True,r2=300):
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 if imin is None:
    imin=0
 if imax is None:
    max=imin+1
 Cooling=[]
 Cooling2=[]
 time=[]
 Edot=[]
 for i in range(imin,imax,mod):
     stest = 'stest_%04i' %i
     fn = name_pattern %(i,i)
     pf = load(fn)
     current_time=pf.current_time
     time.append(current_time)
     Edot.append(pf["ClusterSMBHJetEdot"]*2.0*1.0e44)
     sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
     Cooling.append(sp.quantities["TotalQuantity"](["CoolingRate"])[0])
     sp2 = pf.h.sphere([0.5, 0.5, 0.5], r2/pf['kpc'])
     Cooling2.append(sp2.quantities["TotalQuantity"](["CoolingRate"])[0])
 if log:
   pylab.semilogy(time,Edot,color='red')
   pylab.semilogy(time,Cooling,color="blue")
   pylab.semilogy(time,Cooling2,color="green")
 else:
   pylab.plot(time,Edot,color='red')
   pylab.plot(time,Cooling,color="blue")
   pylab.plot(time,Cooling2,color="green")
 pylab.ylim([6e42,6e46])
 pylab.legend(("Heating","Cooling in 100 kpc", "Cooling in 300 kpc"),loc="lower right")
 pylab.xlabel("Time (Gyr)")
 pylab.ylabel("Cooling Rate v.s. Jet Heating Rate(ergs/s)")
 pylab.tight_layout()
 pylab.savefig("HeatingRate_in%ikpc.png" %rmax)
 numpy.savez("CoolingHeating_%04i_%04i"%(imin,imax), Cooling=Cooling,Cooling2=Cooling2,time=time,Edot=Edot)
 n_step=imax-imin
 H_all=numpy.zeros(n_step-1)
 C_all=numpy.zeros(n_step-1)
 C_2=numpy.zeros(n_step-1)
 Edot[:]=[x*3.16e16 for x in Edot]   #ergs/Gyr
 Cooling[:] = [x*3.16e16 for x in Cooling]  #ergs/Gyr
 Cooling2[:] = [x*3.16e16 for x in Cooling2]  #ergs/Gyr
 for i in range(n_step-1):
     if i == 0:
       H_all[i]=Edot[i]*(time[i+1]-time[i])
       C_all[i]=Cooling[i]*(time[i+1]-time[i])
       C_2[i]=Cooling2[i]*(time[i+1]-time[i])
     else:
       H_all[i]=H_all[i-1]+Edot[i]*(time[i+1]-time[i])
       C_all[i]=C_all[i-1]+Cooling[i]*(time[i+1]-time[i])
       C_2[i]=C_2[i-1]+Cooling2[i]*(time[i+1]-time[i])
 pylab.figure(0)
 if log:
     pylab.semilogy(time[0:n_step-1],H_all,color='red')
     pylab.semilogy(time[0:n_step-1],C_all,color='blue')
     pylab.semilogy(time[0:n_step-1],C_2,color='green')
 else:
     pylab.plot(time[0:n_step-1],H_all,color='red')
     pylab.plot(time[0:n_step-1],C_all,color='blue')
     pylab.plot(time[0:n_step-1],C_2,color='green')
 pylab.ylim([9e58,6e62])
 pylab.xlabel("Time (Gyr)")
 pylab.ylabel("Total Energy (ergs)")
 pylab.legend(("Heating","Cooling in 100 kpc", "Cooling in 300 kpc"),loc="lower right")
 pylab.tight_layout()
 pylab.savefig("All_heating_cooling.png")
 pylab.clf()
 return

def stars(imin=None,imax=None, rmax=100, mod=5):
 time=[]
 Mstar=[]
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 for i in range(imin,imax,mod):
   fn = name_pattern %(i,i)
   pf = load(fn)
   current_time = pf.current_time
   time.append(current_time)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   try:
      m=sp.quantities["TotalQuantity"](["StarMassMsun"])[0]
   except:
      m=1.0e-5
   Mstar.append(m)
 time[:] = [x*TimeUnits/3.16e16 for x in time]  # from code unit to Gyr
 pylab.figure(1)
 pylab.plot(time,Mstar)
 pylab.xlabel("Time (Gyr)")
 pylab.ylabel("Stellar Mass (Sollar Mass)")
 pylab.savefig('stellar_mass_in%ikpc.png'%rmax)
 pylab.clf()
 return

def _young_star_density(field,data):
   star_time=numpy.array(data["star_creation_time"])
   star_density=numpy.array(data["star_density"])
   #old_time=data.pf.current_time-2.0e8/data.pf.time_units['years'] #200Myr in code unit
   old_time=data.pf.current_time-5.0e7/data.pf.time_units['years'] #50Myr in code unit
   old = star_time < old_time
   star_density[old]=0.
   return star_density
add_field("young_star_density",units=r"\rm{g cm^3", function=_young_star_density)

def _young_star_age(field,data):
   star_time=numpy.array(data["star_creation_time"])
   old_time=data.pf.current_time-2.0e8/data.pf.time_units['years'] #200Myr in code unit
   old = star_time < old_time
   star_time[old]=1.0e-40
   return star_time
add_field("young_star_age", function=_young_star_age)

def _MassDepositionRate(field,data):
    return numpy.array(data["CellMassMsun"])/(numpy.array(data["Cooling_Time"])/3.15e7)   #in Solarmass per year
add_field("MassDepositionRate", units = r"M_\odot/yr",function=_MassDepositionRate)

def _t_dyn(field,data):
   conc=6.81
   M_vir=8.5e14
   Rs=0.0224
   r=numpy.array(data["Radiuskpc"])
   x1=r/(Rs*16000.0)
   M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass
   t_dyn=numpy.pi*(r*kpc)**1.5/(2.0*(GravConst*M_dyn)**0.5)/3.15e7   #from s to yr
   return t_dyn
add_field("t_dyn", units=r"\rm yr", function=_t_dyn)

def _tcooltdyn(field,data):
   ratio=numpy.array(data["CoolingTime"])/numpy.array(data["t_dyn"])
   return ratio
add_field("tcooltdyn", function=_tcooltdyn)



def _MassDepositionRate_unstable(field,data):
    t_dyn=numpy.array(data["t_dyn"])
    t_cool=numpy.array(data["CoolingTime"]) #my coolingtime
    stable= t_cool > t_dyn
    rate=numpy.array(data["MassDepositionRate"])
    rate[stable]=1.0e-10
    return rate
add_field("MassDepositionRate_unstable", units=r"M_\odot/yr",function=_MassDepositionRate_unstable)

def _MassDepositionRate_unstable10(field,data):
    t_dyn=numpy.array(data["t_dyn"])
    t_cool=numpy.array(data["CoolingTime"]) #my coolingtime
    stable= t_cool > t_dyn*10.0
    rate=numpy.array(data["MassDepositionRate"])
    rate[stable]=1.0e-10
    return rate
add_field("MassDepositionRate_unstable10", units=r"M_\odot/yr",function=_MassDepositionRate_unstable10)


def _MassDepositionRate_unstabler(field,data):
    t_dyn=numpy.array(data["t_dyn"])
    t_cool=numpy.array(data["CoolingTime"]) #my coolingtime
    stable= (t_cool > t_dyn) | (numpy.array(data["RadialVelocity"])<0)
    rate=numpy.array(data["MassDepositionRate"])
    rate[stable]=1.0e-10
    return rate
add_field("MassDepositionRate_unstabler", units=r"M_\odot/yr",function=_MassDepositionRate_unstabler)

def SFR(imin=0,imax=None, save=False, rmin=5, rmax=150, mod=5): 
# classic cooling rate: cooling radius is where tcool = 5 Gyr here This is from Mittal et al 2009
# Mark Voit: calculate Mdot in Mark's way--Max(Mdot)
# Brian OShea: summing up Mdot in each cell 
# Mark 2: Mdot within r<r_precipitation
# Mark 3: cell-by-cell (t_cool<t_dyn) gas
# Mark 4: cell-by-cell (t_cool/t_dyn < 10 ) gas
#Brian, Mark 3 and Mark 4 are for hot gas only

 from yt.analysis_modules.spectral_integrator.api import add_xray_emissivity_field
 add_xray_emissivity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
 from yt.analysis_modules.spectral_integrator.api import add_xray_luminosity_field
 add_xray_luminosity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
 time=[]
 SFR=[]
 Mdot_classic=[]  #maximum cooling rate 
 Mdot_BH=[]
 Mdot_Mark=[]
 Mdot_Mark2=[]
 Mdot_Mark3=[]
 Mdot_Mark4=[]
 Mdot_Mark10=[]
 Lx=[]
 Mdot_Brian=[]
 dust_all=[]
 dust_out=[]
 SFR_out=[]
 if imax is None:
   imax=imin+1
 pylab.figure(0)
 for i in range(imin,imax):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   current_time = pf.current_time
   time.append(current_time)
   Mdot_BH.append(pf["ClusterSMBHJetMdot"]*2.0)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   prof=BinnedProfile1D(sp, 128, "Radiuskpc", 1,rmax,lazy_reader=True)
   prof.add_fields("Cooling_Time", weight="Xray_Emissivity_0.5_9.9keV")
   try:
     prof.add_fields("TotalMassMsun", weight=None, accumulation=True)
     mass_msun=numpy.array(prof['TotalMassMsun'][prof['UsedBins']])
   except:
     prof.add_fields("CellMassMsun", weight=None, accumulation=True)
     mass_msun=numpy.array(prof['CellMassMsun'][prof['UsedBins']])
   r=numpy.array(prof['Radiuskpc'][prof['UsedBins']])
   CoolingTime=numpy.array(prof['Cooling_Time'][prof['UsedBins']])/3.15e7
   cool=CoolingTime<5.0*1.0e9 #5Gyr. ok???
   ICM_mass=mass_msun[cool][-1]
   tcool=CoolingTime[cool][-1]
   Mdot_classic.append(ICM_mass/tcool) #solarmass per year
   try:
     star_time=sp["creation_time"]
     star_mass=sp["ParticleMassMsun"]
     dt=pf["dtDataDump"]
     good=(star_time> pf.current_time-dt) & (star_time < pf.current_time)
     SFR.append(sum(star_mass[good])/(dt*TimeUnits/3.16e7))  #solarmass per year
     dust_all.append(sum(star_mass[good])*0.02/(dt*TimeUnits/3.16e7))   # dust creation rate
     star_radius=sp["ParticleRadiuskpc"]
     better=good & (star_radius>rmin)
     dust_out.append(sum(star_mass[better])*0.02/(dt*TimeUnits/3.16e7))
     SFR_out.append(sum(star_mass[better])/(dt*TimeUnits/3.16e7))  #SFR > rmin
   except:
     SFR.append(1.0e-10)  # append a tiny number if stars are not formed yet
     SFR_out.append(1.0e-10)
     dust_all.append(1.0e-10)
     dust_out.append(1.0e-10)
# Mdot Mark
   Mdot_r=mass_msun/CoolingTime
   Mdot_Mark.append(max(Mdot_r))
#compute t_dyn
   x1=r/(Rs*pf['kpc'])
   M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass #DM+BCG + BH , cgs
   t_dyn=numpy.pi*(r*kpc)**1.5/(2.0*(GravConst*M_dyn)**0.5)/3.15e7   #from s to yr
# Mark2
   good=(r<rmax-10)&(CoolingTime/10.0<=t_dyn)
   try:
      Mdot_Mark2.append(Mdot_r[good][-1])
   except:
      Mdot_Mark2.append(1.0e-10)  # some small number when Mdot is 0
# Mark3 cell by cell 
#   print "MassDepositionRate_unstable", sp.quantities["TotalQuantity"](["MassDepositionRate_unstable"])[0]
#   Mdot_Mark3.append(sp.quantities["TotalQuantity"](["MassDepositionRate_unstable"])[0])
#   Mdot_Mark4.append(sp.quantities["TotalQuantity"](["MassDepositionRate_unstable10"])[0])

#Brian
   sp_hot=sp.cut_region(["grid['Temperature']>1e6"])
   Mdot_Brian.append(sp_hot.quantities["TotalQuantity"](["MassDepositionRate"])[0])
# Mark3 cell by cell
   Mdot_Mark3.append(sp_hot.quantities["TotalQuantity"](["MassDepositionRate_unstable"])[0])
   #Mdot_Mark4.append(sp_hot.quantities["TotalQuantity"](["MassDepositionRate_unstable10"])[0])
   Mdot_Mark4.append(sp_hot.quantities["TotalQuantity"](["MassDepositionRate_unstabler"])[0])
   Mdot_Mark10.append(sp_hot.quantities["TotalQuantity"](["MassDepositionRate_unstable10"])[0])
   Lx.append(sp_hot.quantities["TotalQuantity"](["Xray_Luminosity_0.5_9.9keV"])[0])
#      large=(r>8)&(r<rmax-5)
   if (i-imin)%mod==0:
      pylab.plot(r,numpy.array(Mdot_r))
 pylab.xlabel("Radius (kpc)")
 pylab.ylabel("Mdot (Solar Mass/yr)")
 pylab.savefig("Mdot_Mark_r.png")
 pylab.clf()
 time[:] = [x*TimeUnits/3.16e16 for x in time]  # from code unit to Gyr
 pylab.figure(1)
 pylab.plot(time,SFR)
 pylab.plot(time,Mdot_classic)
 pylab.xlabel("Time (Gyr)")
 pylab.ylabel("SFR (Sollar Mass/yr)")
 pylab.savefig('SFR_in%ikpc.png'%rmax)
 pylab.clf()
 pylab.figure(2)
 pylab.plot(time,dust_all)
 pylab.plot(time,dust_out)
 pylab.xlabel("Time (Gyr)")
 pylab.ylabel("dust creation rate (Sollar Mass/yr)")
 pylab.savefig("dust.png")
 pylab.clf()
 if save==True:
    numpy.savez("SFR_classic_in%ikpc" %rmax, SFR=SFR,timeGyr=time,Mdot_classic=Mdot_classic, Mdot_BH=Mdot_BH,Mdot_Mark=Mdot_Mark, Mdot_Mark2=Mdot_Mark2, SFR_out=SFR_out,dust_all=dust_all,dust_out=dust_out, Mdot_Brian=Mdot_Brian,Mdot_Mark3=Mdot_Mark3,Mdot_Mark4=Mdot_Mark4,Mdot_Mark10=Mdot_Mark10,Lx=Lx)
 return


def cold_disk(imin=None,imax=None,rmax=None, mod=1):
 pylab.rc('axes', linewidth=2,labelsize=16,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 16}
 pylab.rc('font', **font)
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax=50
 conc=6.81
 M_vir=8.5e14
 Rs=0.0224   #SphereCoreRadius=0.0224
 cm=pylab.cm.gist_rainbow
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   pc.add_profile_sphere(rmax, "kpc", ["Radiuskpc", "TangentialVelocity"],weight="CellMassMsun")
   V_t=pc.plots[-1].data["TangentialVelocity"]     
   r=pc.plots[-1].data["Radiuskpc"]
   x1=r/(Rs*pf['kpc'])
   M_dyn=(M_vir*SolarMass)*(numpy.log(1.0+x1)-x1/(1.0+x1))/(numpy.log(1.0+conc)-conc/(1.0+conc))+(((r**0.5975)/3.206e-7)**0.9+((r**1.849)/1.861e-6)**0.9)**(-1.0/0.9)*(r*kpc)**2/GravConst+(3.4e8)*SolarMass #DM+BCG + BH , cgs
   V_kep=(GravConst*M_dyn/(r*kpc))**0.5   # kepler velocity
   pylab.figure(i)
   pylab.plot(r, V_t/1.0e5)
   pylab.plot(r, V_kep/1.0e5,linestyle='--')
   pylab.xlabel("r (kpc)")
   pylab.ylabel("Velocity (km/s)")
   pylab.ylim([0,5e2])
   pylab.xlim([0,rmax-1])
   pylab.savefig(stest+'_disk_Velocity.png')
   pylab.clf()
 return   


def delta_density(imin=None,imax=None,rmin=10,rmax=20, mod=1):
 import scipy as scipy
 from scipy.stats import norm
 if imin is None:
    imin=0
 if imax is None:
    imax=imin+1
 if rmax is None:
    rmax = rmin+1
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   cap=sp.cut_region(["grid['Theta']<0.3","grid['Radiuskpc']>10.0","grid['Temperature']>1e5"])   # it is a shell within the jet cone. it is a cap   # also exclude clumps
   rho=cap["Density"]
   mass=cap["CellMass"]
#   hist, bin_edges = np.histogram(rho,weights=mass, bins=50, density=True)
   pylab.figure(i)
   pylab.hist(rho,weights=mass, bins=50, normed = True, histtype='step')
#   pylab.xlim([])
   pylab.savefig(stest+'_rho_distribution_in%i.png'%rmax)
   pylab.clf()
   normed_rho=rho*mass/numpy.sum(mass)
   nloc, nscale = scipy.stats.norm.fit(normed_rho)
   print "nloc, nscale=", nloc, nscale
#   shape, loc, scale = scipy.stats.lognorm.fit(,loc=0)
#  time.append(this_time)    #in code unit
#  Nclumps.append(len(mass))  #number of clumps
#  Mean.append(scale)   #mean mass in Msun
 return


def Voit(imin=0,imax=1,mod=1):
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 m=numpy.load("ColdGasMass_in100kpc.npy")
 ratio=numpy.load("min_ratio.npy")
 j=0
 cm=pylab.cm.gist_rainbow
 for i in range(imin,imax,mod):
   #pylab.semilogy([ratio[j]],[numpy.max(1.0e8,m[i])],marker="o",color=cm(float(i-imin)/(imax-imin)))
   pylab.semilogy([ratio[j]],[m[i]],marker="o",color=cm(float(i-imin)/(imax-imin)))
   j+=1
 pylab.xlabel("Minimum "+r"$t_{cool}/t_{dyn}$")
 pylab.ylabel("Cold Gas Mass (SolarMass)")
 pylab.savefig("Voit.png")
 pylab.clf()
 return

def Cooling_Mdot():
 pylab.rc('axes', linewidth=2,labelsize=18,labelweight='bold')
 pylab.rc('lines', linewidth=2)
 font = {'weight' : 'bold', 'size'   : 18}
 pylab.rc('font', **font)
 pylab.rc('text', usetex=True)
 pylab.rcParams['text.latex.preamble'] = [r'\boldmath']
 coldgas=numpy.load("ColdGasMass_in100kpc.npy")
 #print coldgas, coldgas.shape
 time=numpy.genfromtxt("Time")[:,2]
 Mdot_jet=2.0*numpy.genfromtxt("Mdot")[:,2]   # Msun/yr
 n_step=len(time)
 Mdot=numpy.zeros(n_step)
 Mdot_disk=numpy.zeros(n_step)
 time*=TimeUnits  #in s
 time/=3.16e7    #in yr
 for i in range(n_step-2):
   Mdot_disk[i]=(coldgas[i+1]-coldgas[i])/(time[i+1]-time[i])   ##Msun/yr
 Mdot=Mdot_jet + Mdot_disk
 pylab.figure(0)
 pylab.semilogy(time[:n_step-2]/1.0e9,Mdot[:n_step-2])
 pylab.xlabel("Time (Gyr)")
 pylab.ylabel("Cooling Rate (Solar Mass/yr)")
 pylab.savefig("Mdot.png")
 return

def MultiPlots_Paper(ilist=[0],field_list=["Density","Temperature"], rmax=40):
 import matplotlib.colorbar as cb
 from matplotlib.colors import LogNorm
 from yt.analysis_modules.spectral_integrator.api import add_xray_emissivity_field
 add_xray_emissivity_field(0.5, 9.9, with_metals=False, constant_metallicity=0.5)
 keys = ['Density', 'Temperature','Entropy','Pressure','star_density', 'young_star_density','Xray_Emissivity_0.5_9.9keV']
 values_titles = [r'$\mathrm{Density}\ (\mathrm{g\ cm^{-3}})$', r'$\mathrm{Temperature}\ (\mathrm{K})$', \
                 r'$\mathrm{Entropy}\ (\mathrm{ergs}\ \mathrm{cm}^{3\gamma-3})$', r'$\mathrm{Pressure}\ (\mathrm{dyne}/\rm{cm}^{2})$',\
                 r'$\mathrm{Stellar\ Density}\ (\mathrm{g\ cm^{-3}})$',r'$\mathrm{Young\ stellar\ density}\ (\mathrm{g\ cm^{-3}})$', 
                 r'$\mathrm{L_{Xray}}\ \mathrm{ergs\ s}^{-1}$']
 #values_z_list=[(8e-3,0.5),(1e6,1e8),(9e-10,2e-7),(8e-11,1e-8),(5.0e-5,1.0e-2),(5.0e-6,1.0e-2),(9.0e-4,6.0e-2)]
 values_z_list=[(7e-26,1e-23),(3e6,1e8),(9e-10,2e-7),(8e-11,1e-8),(5.0e-5,1.0e-2),(5.0e-6,1.0e-2),(9.0e-4,6.0e-2)]
 titles=dict(zip(keys, values_titles))
 z_list = dict(zip(keys, values_z_list))
 n_x=len(field_list)
 n_y=len(ilist)
 orient = 'horizontal'
 fig, axes, colorbars = get_multi_plot(n_x, n_y, colorbar=orient, bw = 4)
 center= numpy.array([0.5, 0.5, 0.5])
 ii=0
 res = [1000, 1000]
 plots=[]
 for i in ilist:
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf=load(fn)
   LE=center-rmax/pf['kpc']
   RE=center+rmax/pf['kpc']
   reg = pf.h.region(center,LE,RE)
   width=2*rmax/pf['kpc']
   cmap_UV=matplotlib.cm.Blues_r
   cmap_UV.set_bad('k')
   cmap_X=matplotlib.cm.Purples_r
   cmap_star=matplotlib.cm.Greys
   cmap_star.set_bad('w')
#   cmap_star=matplotlib.cm.Greys_r
#   cmap_star.set_bad('k')
   cmap_defaul='algae'
   values_cmap=[cmap_defaul,cmap_defaul,cmap_defaul,cmap_defaul, cmap_star, cmap_UV, cmap_X]
   cmap_list = dict(zip(keys, values_cmap))
   jj=0
   for j in field_list:
       #if (j=="young_star_density" or j=="Density" or j=="Xray_Emissivity_0.5_9.9keV" or j=="star_density"):
       if (j=="young_star_density" or j=="Xray_Emissivity_0.5_9.9keV" or j=="star_density"):
          proj = pf.h.proj(0, j, source=reg)
       else:
          #proj = pf.h.proj(0, j, source=reg, weight_field="Density")
          proj = pf.h.proj(0, j, source=reg, weight_field="Xray_Emissivity_0.5_9.9keV")
       axes[ii][jj].xaxis.set_visible(False)
       axes[ii][jj].yaxis.set_visible(False)
       frb = proj.to_frb(width,res)
       aa=axes[ii][jj].imshow(frb[j],norm=LogNorm())
       if (j=="star_density"):
         time=(0.048108+0.0017*(i-36))*TimeUnits/3.16e16
         axes[ii][jj].text(680,80,"%.2f Gyr"%time,fontsize=22)
       plots.append(aa)
       plots[-1].set_cmap(cmap_list[j])
       plots[-1].set_clim(z_list[j])
       jj+=1
   ii+=1
 titles_use=[titles[x] for x in field_list]
 for p, cax, t in zip(plots, colorbars,titles_use):
   cbar = fig.colorbar(p, cax=cax, orientation=orient)
   cbar.set_label(t)
 fig.savefig("multi.png")
 return

#proj = pf.h.proj(0, "StarAgeYears", source=sp, weight_field="young_star_density")
#frb = proj.to_frb(width,res)
#StarAge=frb["StarAgeYears"]
#write_projection(StarAge)
#pylab.savefig(frb)

def stars_profile(imin=0,imax=None,rmin=0.1, rmax=200.0, mod=1):
 if imax is None:
    imax=imin+1
 cm=pylab.cm.gist_heat
 pylab.figure(0)
 for i in range(imin,imax,mod):
    stest = 'stest_%04i' %i
    fn = name_pattern %(i,i)
    pf = load(fn)
    center = [0.5, 0.5, 0.5]
    sp = pf.h.sphere(center, rmax/pf['kpc'])
    prof=BinnedProfile1D(sp, 128, "Radiuskpc", rmin,rmax,lazy_reader=True)
    r_prof=prof["Radiuskpc"]
    try:
      prof.add_fields("StarMassMsun",weight=None, accumulation=True)
      star_mass=prof["StarMassMsun"]
      pylab.loglog(r,star_mass,color=cm(float(i-imin)/(imax-imin)))
    except:
      pass
 pylab.xlabel('r (kpc)')
# pylab.ylabel('Stellar Density '+r'$(\rm{g}\ /\ \rm{cm}^3)$')
 pylab.ylabel('Enclosed Stellar Mass '+r'$(M_{\odot})$')
 pylab.xlim([rmin,rmax])
#    pylab.ylim([1e-27,3e-21])
# pylab.tight_layout()
 pylab.savefig('StarMass_rainbow.png')

def _ZeusHeatingRate(field,data):
   q=2.0*(data['dx'].flat[0]*data.convert("cm"))**2.0*numpy.array(data["Density"])*numpy.array(data["DivV"])**2.0
   bad=numpy.array(data["DivV"])>=0
   q[bad]=0.0
   return -q*numpy.array(data["DivV"])
add_field("ZeusHeatingRate",function=_ZeusHeatingRate, units = r"ergs s^{-1}cm^{-3}")

def _ZeusHeatingCell(field,data):
    return numpy.array(data["ZeusHeatingRate"])*numpy.array(data["CellVolume"])
add_field("ZeusHeatingCell",function=_ZeusHeatingCell, units = r"ergs s^{-1}")

#def _ShockOrNot(field,data):
#    return

def _ZeusShockHeatingRate(field,data):
   ZeusShockHeatingRate=numpy.array(data["ZeusHeatingRate"])
   bad = numpy.array(data["MachNumber"])<1
   ZeusShockHeatingRate[bad]=0.0
   return ZeusShockHeatingRate
add_field("ZeusShockHeatingRate", function=_ZeusShockHeatingRate, units = r"ergs s^{-1}cm^{-3}")

def _ZeusShockHeatingCell(field,data):
    return numpy.array(data["ZeusShockHeatingRate"])*numpy.array(data["CellVolume"])
add_field("ZeusShockHeatingCell",function=_ZeusShockHeatingCell, units = r"ergs s^{-1}")

def Heating_Zeus(imin=0,imax=None,rmin=0.5, rmax=300, mod=1,n_bins=40):
 if imax is None:
   imax=imin+1
 r=numpy.zeros(((imax-imin)//mod,n_bins))
 Cooling=numpy.zeros(((imax-imin)//mod,n_bins))
 Heating=numpy.zeros(((imax-imin)//mod,n_bins))
 ShockHeating=numpy.zeros(((imax-imin)//mod,n_bins))
 j=0
 for i in range(imin,imax,mod):
   stest = 'stest_%04i' %i
   fn = name_pattern %(i,i)
   pf = load(fn)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   sp_good=sp.cut_region(["grid['Radiuskpc']>2.0","grid['Temperature']>1e6"])
   prof=BinnedProfile1D(sp_good, n_bins, "Radiuskpc", rmin,rmax,lazy_reader=True)
   prof.add_fields("CoolingRate", weight=None, accumulation=True)
   prof.add_fields("ZeusHeatingCell",weight=None, accumulation=True)
   prof.add_fields("ZeusShockHeatingCell",weight=None, accumulation=True)
   r[j,:]=numpy.array(prof["Radiuskpc"][prof['UsedBins']])
   Cooling[j,:]=numpy.array(prof["CoolingRate"][prof['UsedBins']])
   Heating[j,:]=numpy.array(prof["ZeusHeatingCell"][prof['UsedBins']])
   ShockHeating[j,:]=numpy.array(prof["ZeusShockHeatingCell"][prof['UsedBins']])
   
   j=j+1
 numpy.savez("Heating_Zeus_%04i_%04i" %(imin,imax), r=r,Cooling=Cooling, Heating=Heating, ShockHeating=ShockHeating)
 return


def Zeus_more(imin=0,imax=1,rmax=300):
 Heating=[]
 ShockHeating=[]
 Total_Thermal=[]
 Total_Kinetic=[]
 Time=[]
 dEdt=[]
 dEdt_in=[]
 dKdt_in=[]
 Total_K=[]
 dKdt=[]
 dEk_dt=[]
 for i in range(imin,imax):
   fn = name_pattern %(i,i)
   pf = load(fn)
   print "fn=", fn
   pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
   try:
#     dd=pf.h.region([0.5, 0.5, 0.5], [0.01, 0.01, 0.01], [0.99,0.99,0.99])
     dd=pf.h.sphere([0.5, 0.5, 0.5], (rmax,'kpc'))
     Heating.append(dd.quantities["TotalQuantity"](["ZeusHeatingCell"])[0])
     ShockHeating.append(dd.quantities["TotalQuantity"](["ZeusShockHeatingCell"])[0])
     Total_Thermal.append(dd.quantities["TotalQuantity"](["Cell_ThermalEnergy"])[0])
     Total_Kinetic.append(dd.quantities["TotalQuantity"](["Cell_KineticEnergy"])[0])
     Time.append(pf.current_time*pf["TimeUnits"])
   except:
     print "fail"
   pc.save()
#   TotalThermal.append(ds.quantities.total_quantity(["Cell_ThermalEnergy"]))
#   Heating.append(ds.quantities.total_quantity(["ZeusHeatingCell"]))
#   ShockHeating.append(ds.quantities.total_quantity(["ZeusShockHeatingCell"]))
#   time.append(ds.current_time.in_units('s'))
 TotalThermal=numpy.array(Total_Thermal)
 TotalKinetic=numpy.array(Total_Kinetic)
 time=numpy.array(Time)
#print "time", time
#print "timeUnit", pf["TimeUnits"]

 for i in range(0,len(time)-1):
   dEdt.append((TotalThermal[i+1]-TotalThermal[i])/(time[i+1]-time[i]))
   dEk_dt.append((TotalKinetic[i+1]-TotalKinetic[i])/(time[i+1]-time[i]))
 numpy.savez("Heating_compare",Heating=Heating, TotalThermal=TotalThermal, ShockHeating=ShockHeating, time=time,dEdt=dEdt,TotalKinetic=TotalKinetic,dEk_dt=dEk_dt)
 return


def Filament(imin=0,imax=None,rmin=2,rmax=100,mod=1,x_bins=60,y_bins=100):
 if imax is None:
   imax=imin+1
 x_bins=int(rmax-rmin)
 V_in=numpy.zeros(shape=(imax-imin,x_bins+1,y_bins+1))
 V_out=numpy.zeros(shape=(imax-imin,x_bins+1,y_bins+1))
# V_t_in=numpy.zeros(shape=(imax-imin,x_bins+1))
# V_t_out=numpy.zeros(shape=(imax-imin,x_bins+1))
 V_t_in=numpy.zeros(shape=(imax-imin,x_bins+1,y_bins+1))
 V_t_out=numpy.zeros(shape=(imax-imin,x_bins+1,y_bins+1))
 for i in range(imin,imax,mod):
   ii=i-imin
   fn = name_pattern %(i,i)
   pf = load(fn)
   sp = pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
   sp_out=sp.cut_region(["grid['Radiuskpc']>2.0","grid['Temperature']<1e5","grid['RadialVelocity']>0"])
   sp_in=sp.cut_region(["grid['Radiuskpc']>2.0","grid['Temperature']<1e5","grid['RadialVelocity']<0"])
   prof2D_out=BinnedProfile2D(data_source=sp_out,x_n_bins=x_bins,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=y_bins,y_bin_field="VelocityMagnitude",y_lower_bound=0,y_upper_bound=1e8,y_log=False,lazy_reader=True)
   prof2D_out.add_fields("CellMassMsun", weight=None)
   V_out[ii,:,:]=prof2D_out["CellMassMsun"]
   prof2D_in=BinnedProfile2D(data_source=sp_in,x_n_bins=x_bins,x_bin_field="Radiuskpc",x_lower_bound=rmin,x_upper_bound=rmax,x_log=False,y_n_bins=y_bins,y_bin_field="VelocityMagnitude",y_lower_bound=0,y_upper_bound=1e8,y_log=False,lazy_reader=True)
   prof2D_in.add_fields("CellMassMsun", weight=None)
   V_in[ii,:,:]=prof2D_in["CellMassMsun"]   
#   prof1D_out=BinnedProfile1D(sp_out, x_bins, "Radiuskpc", rmin,rmax,lazy_reader=True)
#   prof1D_out.add_fields("TangentialOverVelocityMagnitude",weight="CellMassMsun")
#   prof1D_in=BinnedProfile1D(sp_in, x_bins, "Radiuskpc", rmin,rmax,lazy_reader=True)
#   prof1D_in.add_fields("TangentialOverVelocityMagnitude",weight="CellMassMsun")
#   V_t_in[ii,:]=prof1D_in["TangentialOverVelocityMagnitude"]
#   V_t_out[ii,:]=prof1D_out["TangentialOverVelocityMagnitude"]
   prof2D_in.add_fields("TangentialOverVelocityMagnitude",weight="CellMassMsun")
   V_t_in[ii,:,:]=prof2D_in["TangentialOverVelocityMagnitude"]
   prof2D_out.add_fields("TangentialOverVelocityMagnitude",weight="CellMassMsun")
   V_t_out[ii,:,:]=prof2D_out["TangentialOverVelocityMagnitude"]
 r_prof=prof2D_in["Radiuskpc"]
 v_prof=prof2D_in["VelocityMagnitude"]
 numpy.savez("Filament_Velocity",r_prof=r_prof, v_prof=v_prof, V_in=V_in,V_out=V_out,V_t_in=V_t_in,V_t_out=V_t_out)
 return


def star_ray(i=50):
 fn = name_pattern %(i,i)
 pf = load(fn)
 stest = 'stest_%04i' %i
 dd=pf.h.all_data()
 xx=dd["particle_position_x"]
 yy=dd["particle_position_y"]
 zz=dd["particle_position_z"]
 sp=pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
 sp_cold=sp.cut_region(["grid['Radiuskpc']>1.0","grid['Temperature']<1e6"])
 numpy.savez(stest+"star_position",xx=xx,yy=yy,zz=zz)
 return 

def rho_cold(imin=0,imax=None,rmax=200):
 if imax is None:
   imax=imin+1
 rho_cold=[]
 for i in range(imin,imax):
  fn = name_pattern %(i,i)
  pf = load(fn)
  stest = 'stest_%04i' %i
  sp=pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
  sp_cold=sp.cut_region(["grid['Radiuskpc']>1.0","grid['Temperature']<1e6"])
  #try:
  #mean_d=sp_cold.quantities["WeightedAverageQuantity"]("Density",None)
  mean_d=sp_cold.quantities["WeightedAverageQuantity"]("Density","Ones")
  print mean_d
  rho_cold.append(mean_d)
  #except:
  # rho_cold.append(0.01)
 numpy.savez("cold_density",rho_cold=rho_cold)
 return

if __name__=="__main__":
  import sys
  profile(sys.argv[1])
