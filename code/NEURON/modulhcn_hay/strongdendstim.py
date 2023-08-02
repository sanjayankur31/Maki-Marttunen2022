#cp calcifcurves_dendstim.py strongdendstim.py
#cp ../haymod3e/runcontrol_check.py runcontrol.py
# runcontrols
# A script for determining the control neuron F-I curve and limit cycle.
#
# The input code for the hoc-interface is based on BAC_firing.hoc by Etay Hay (2011)
#
# Tuomo Maki-Marttunen, Oct 2014
# (CC BY)

from neuron import h
import mytools
import pickle
import numpy as np
import sys
from pylab import *

spikfreqsAll = []
timescAll = []
VsomacAll = []
VDerivcAll = []
VDcoeffAll = []
VdendcAll = []
VdDcoeffAll = []
VdDerivcAll = []
CasomacAll = []
CaDerivcAll = []
CaDcoeffAll = []
CadendcAll = []
CadDerivcAll = []
CadDcoeffAll = []
times_controlAll = []
Vsoma_controlAll = []
Vdend_controlAll = []
Casoma_controlAll = []
Cadend_controlAll = []

Ihcoeff = 1.0
denddist = 620
if len(sys.argv) > 1:
  Ihcoeff = float(sys.argv[1])
if len(sys.argv) > 2:
  denddist = float(sys.argv[2])

for icell in range(0,1):
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  v0 = -80
  ca0 = 0.0001

  proximalpoint = 400
  distalpoint = denddist
  BACdt = 5.0

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.0002)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
access L5PC.soma
objref st1

objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma, ihsoma, ihdend
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
ihsoma = new Vector()
ihdend = new Vector()
objref sl,ns,tvec
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "distalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "distalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)
L5PC.apic[siteVec[0]] cvode.record(&ihcn_Ih(0.5),ihdend,tvec)
L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
L5PC.soma cvode.record(&cai(0.5),casoma,tvec)
L5PC.soma cvode.record(&ihcn_Ih(0.5),ihsoma,tvec)
L5PC.apic[siteVec[0]] st1 = new IClamp(siteVec[1])
""")

  #Block or amplify Ih channels:
  h("""forall if(ismembrane("Ih")) { gIhbar_Ih = """+str(Ihcoeff)+"""*gIhbar_Ih } """)

  Is = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 100.0]
  times_all = []
  Vsoma_all = []
  Vdend_all = []
  ihsoma_all = []
  ihdend_all = []
  f,axarr = subplots(2,len(Is))
  for iI in range(0,len(Is)):
    squareAmp = Is[iI]
    squareDur = 0.2
    tstop = 1000
    h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
st1.del = 800
""")
    h.init()
    h.run()

    times=np.array(h.tvec)
    Vsoma=np.array(h.vsoma)
    Vdend=np.array(h.vdend)
    Casoma=np.array(h.casoma)
    Cadend=np.array(h.cadend)
    ihsoma=np.array(h.ihsoma)
    ihdend=np.array(h.ihdend)

    times_all.append(times[:])
    Vsoma_all.append(Vsoma[:])
    Vdend_all.append(Vdend[:])
    ihsoma_all.append(ihsoma[:])
    ihdend_all.append(ihdend[:])
    axarr[0,iI].plot(times,Vdend)
    axarr[1,iI].plot(times,Vsoma)
    axarr[0,iI].set_xlim([750,950])
    axarr[1,iI].set_xlim([750,950])

    picklelist = [times_all,Vsoma_all,Vdend_all,ihsoma_all,ihdend_all]
    file = open('strongdendstim'+str(denddist)+'_Ihcoeff'+str(Ihcoeff)+'.sav', 'wb')
    pickle.dump(picklelist,file)
    file.close()

    f.savefig('strongdendstim'+str(denddist)+'_Ihcoeff'+str(Ihcoeff)+'.eps')




