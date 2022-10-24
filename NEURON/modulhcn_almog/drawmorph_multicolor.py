from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import random

v0 = -80
ca0 = 0.0001
BACdt = 5.0
fs = 8
tstop = 13000.0
icell = 0

h(""" load_file("myrun.hoc") """)

close("all")
f,axarr = subplots(1,1)
plotteds = []
dists = []
mynums0 = []
mynums1 = []
for itree in range(0,3):
  if itree == 0:
    nsec = len(h.dend)
  elif itree == 1:
    nsec = len(h.apic)
  else:
    nsec = 1

  for j in range(nsec-1,-1,-1):
    if itree == 0:
      h("access dend["+str(j)+"]")
    elif itree == 1:
      h("access apic["+str(j)+"]")
    else:
      h("access a_soma")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarz = z3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    h("tmpvarz2 = z3d(n3d()-1)")
    coord1 = [h.tmpvarx,h.tmpvary,h.tmpvarz]
    coord2 = [h.tmpvarx2,h.tmpvary2,h.tmpvarz2]
    thisdist = h.distance(0.5)
    dists.append(h.distance(0.5))

    col = "#000000"
    if itree == 0:
      mynum = int(thisdist/25)
      mychr = hex(15-mynum)[2]
      mychr2 = hex(mynum)[2]
      if mynum>15:
        print("HEP")
      col = "#"+'ff'+mychr2+mychr2+'00'
      mynums0.append(mynum)
    elif itree == 1:
      mynum = int(thisdist/100)
      mychr = hex(mynum)[2]
      mychr2 = hex(15-mynum)[2]
      if mynum>15:
        print("HEP")
      col = "#00"+mychr2+mychr2+"ff"
      mynums1.append(mynum)
    elif itree == 2:
      col = "#000000"

    h("""
myn = n3d()
myx0 = x3d(0)
myy0 = y3d(0)
myz0 = z3d(0)
""")
    oldcoord = [h.myx0, h.myy0, h.myz0]
    for k in range(1,int(h.myn)):
      h("""
myx0 = x3d("""+str(k)+""")
myy0 = y3d("""+str(k)+""")
myz0 = z3d("""+str(k)+""")
mydiam = diam""")
      axarr.plot([oldcoord[0],h.myx0],[oldcoord[1],h.myy0],'k-',linewidth=h.mydiam*0.25,color=col)
      plotteds.append([[oldcoord[0],h.myx0],[oldcoord[1],h.myy0],'k-',h.mydiam*0.25,col])
      oldcoord = [h.myx0, h.myy0, h.myz0]
axis("equal")
f.savefig("morph_multicolor.eps")

file = open('morph_multicolor.sav', 'wb')
pickle.dump(plotteds,file)
file.close()

print(str(max(mynums0)))
print(str(max(mynums1)))
  

