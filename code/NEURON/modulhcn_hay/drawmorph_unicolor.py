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

h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode           
cvode = new CVode()    
cvode.active(1)        
cvode.atol(0.00005)    
load_file("import3d.hoc")
objref L5PC      
load_file(\"models/L5PCbiophys3.hoc\")
load_file(\"models/L5PCtemplate.hoc\")
L5PC = new L5PCtemplate(\"morphologies/cell1.asc\")
access L5PC.soma
""")

close("all")
f,axarr = subplots(1,1)
plotteds = []
for itree in range(0,3):
  if itree == 0:
    nsec = len(h.L5PC.dend)
  elif itree == 1:
    nsec = len(h.L5PC.apic)
  else:
    nsec = 1

  for j in range(nsec-1,-1,-1):
    if itree == 0:
      h("access L5PC.dend["+str(j)+"]")
    elif itree == 1:
      h("access L5PC.apic["+str(j)+"]")
    else:
      h("access L5PC.soma")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarz = z3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    h("tmpvarz2 = z3d(n3d()-1)")
    coord1 = [h.tmpvarx,h.tmpvary,h.tmpvarz]
    coord2 = [h.tmpvarx2,h.tmpvary2,h.tmpvarz2]
    col = "#000000"
    if itree == 0:
      col = "#000000"
    elif 0.5*(coord1[1]+coord2[1]) < 650:
      col = "#000000"
    if itree == 2:
      col = "#000000"

    ## Too sparse. Try the below instead
    #axarr.plot([coord1[0],coord2[0]],[coord1[1],coord2[1]],'k-',linewidth=h.mydiam,color=col)
    
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
      #if itree == 2:
      #    h("mydiam = mydiam/5")
      #axarr.plot([oldcoord[0],h.myx0],[oldcoord[1],h.myy0],'k-',linewidth=h.mydiam*0.67,color=col)
      axarr.plot([oldcoord[0],h.myx0],[oldcoord[1],h.myy0],'k-',linewidth=h.mydiam*0.25,color=col)
      plotteds.append([[oldcoord[0],h.myx0],[oldcoord[1],h.myy0],'k-',h.mydiam*0.25,col])
      oldcoord = [h.myx0, h.myy0, h.myz0]

axis("equal")
f.savefig("morph_unicolor.eps")
  
file = open('morph.sav', 'wb')
pickle.dump(plotteds,file)
file.close()

