import matplotlib
matplotlib.use('Agg')
from yt.mods import *
from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"

name_pattern = 'DD%04i/stest_%04i'

def AddParticle(i=0,rmax=100.0, skip=2):
 fn = name_pattern %(i,i)
 pf = load(fn)
 sp_all=pf.h.sphere([0.5, 0.5, 0.5], rmax/pf['kpc'])
 sp=sp_all.cut_region(["grid['x']>0.5"])
 sp=sp.cut_region(["grid['y']>0.5"])
 sp=sp.cut_region(["grid['z']>0.5"])
 print "sp.size=", sp_all["x"].size, "   sp_first_octant.size=", sp["x"].size

 start_index = 0

 ndata = {}
 ndata["particle_position_x"] = sp["x"][ : :skip]
 ndata["particle_position_y"] = sp["y"][ : :skip]
 ndata["particle_position_z"] = sp["z"][ : :skip]
 ndata["particle_velocity_x"] = na.zeros(sp["x"][ : :skip].size)
 ndata["particle_velocity_y"] = na.zeros(sp["x"][ : :skip].size)
 ndata["particle_velocity_z"] = na.zeros(sp["x"][ : :skip].size)
 ndata["particle_type"] = na.ones(sp["x"][ : :skip].size, dtype="int64") * 3
 ndata["particle_index"] = na.arange(start_index, start_index + sp["x"][ : :skip].size)
 ndata["particle_mass"] = na.zeros(sp["x"][ : :skip].size)
 print ndata["particle_mass"].size

 import h5py
 f = h5py.File("DD%04i/stest_%04i"%(i,i)+".cpu0000", "a")
 print f.keys()

 g = f["/Grid00000001"]
 for field in ndata:
     g.create_dataset(field, data=ndata[field])

 f.close()

 pf2 = load(fn)
 g2 = pf2.h.grids[0]
 print "particle_masssize=", g2["particle_mass"].size
 return

def convert_TP32():
 import struct
 import numpy
 import glob
 import os
 import h5py
 from collections import defaultdict


 header_struct = '6fq'   # Greg float-longlong  gives header_size=32
 header_size = struct.calcsize(header_struct)
 print "header_size=", header_size

 pcount_struct = 'q'   # Matt longlong
 pcount_size = struct.calcsize(pcount_struct)

 files = [open(f) for f in sorted(glob.glob("TracerOutput*"))]
 ends = [os.stat(f.name).st_size for f in files]

 timestep = 0
 quit = 0
 nattr_vel = 0

 attrs = ["position_x", "position_y", "position_z",
         "density", "temperature"]
 attrs_vel = ["velocity_x", "velocity_y", "velocity_z"]

 while quit == 0:
    data = defaultdict(list)
    for f, e in zip(files, ends):
        header = f.read(header_size)
#        print "header is %s" % (header)
        t0, z, tu, lu, du, vu, nattr = struct.unpack(header_struct, header)
        print "t0, z, tu, lu, du, vu, nattr", t0, z, tu, lu, du, vu, nattr
        total_npart = 0
        while 1:
            npart = struct.unpack(pcount_struct, f.read(pcount_size))[0]
#            print "npart=", npart
            if npart == 0: break
            total_npart += npart

            if (nattr == 8):
                nattr_vel = 3
                nattr = 5

            attr_struct = '%sf' % (nattr)
            #attr_struct = '%sd' % (nattr)    ### Matt double
            attr_size = struct.calcsize(attr_struct)
            print "attr_size= (should be 5)", attr_size
            #attr_data = numpy.fromstring(f.read(npart * attr_size), 'float64')   # Matt
            attr_data = numpy.fromstring(f.read(npart * attr_size), 'float32')  # Greg
            print "attr_data.shape", attr_data.shape, (npart, nattr)
            attr_data.shape = (npart, nattr)
            for i, a in zip(xrange(nattr), attrs):
                data[a].append(attr_data[:,i])

            if (nattr_vel > 0):
                attr_struct = '%sf' % (nattr_vel)
                attr_size = struct.calcsize(attr_struct)
                attr_data = numpy.fromstring(f.read(npart * attr_size), 'float32')  # Greg
                print "attr_vel_data.shape", attr_data.shape, (npart, nattr_vel)
                attr_data.shape = (npart, nattr_vel)
                for i, a in zip(xrange(nattr_vel), attrs_vel):
                    data[a].append(attr_data[:,i])

            data['indices'].append(numpy.fromstring(f.read(npart * pcount_size), 'int64'))
            print data['indices'][-1].min(), data['indices'][-1].max()
        if e == f.tell(): quit = 1
    f = h5py.File("TracerParticles_%08i.h5" % timestep, "w")
    for a in data:
        data[a] = numpy.concatenate(data[a])
        if a == "density":
            data[a] *= du
        if "velocity" in a:
            data[a] *= vu
    ii = numpy.argsort(data["indices"])
#    print "TOTAL: ", data["indices"].size
#    print
    for a in data:
        f.create_dataset("/%s" % a, data = data[a][ii])
    f["/"].attrs["redshift"] = z
    f["/"].attrs["time"] = t0 * tu
    f.close()
    timestep += 1
 return

def particles_npz(nTP=None,nfiles=None): 
 import h5py
 lu = 4.93708012753e+25   # length unit
 kpc=3.086e21
 Rs = 8.31e7
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
 time = numpy.zeros(nfiles)
 for j in range(nfiles):
   f = h5py.File("TracerParticles_%08i.h5" % j, "r")
   time[j] = f["/"].attrs["time"]
   xpos[:,j],ypos[:,j],zpos[:,j]=f["position_x"],f["position_y"],f["position_z"]
   try:
      Vx[:,j],Vy[:,j],Vz[:,j]=f["velocity_x"],f["velocity_y"],f["velocity_z"]
   except:
      pass
   radius[:,j]=((xpos[:,j]-0.5)**2.0+(ypos[:,j]-0.5)**2.0+(zpos[:,j]-0.5)**2.0)**0.5 ##in codeunit
   radius[:,j]*= (lu/kpc) ## now in kpc
   T[:,j]=f["temperature"]
   rho[:,j]=f["density"]
   P[:, j] = rho[:, j]*Rs*T[:,j]
 numpy.savez("particles", rho=rho, T=T, radius=radius, xpos=xpos, ypos=ypos, zpos=zpos, P=P, Vx=Vx, Vy=Vy, Vz=Vz,time=time)
 return
