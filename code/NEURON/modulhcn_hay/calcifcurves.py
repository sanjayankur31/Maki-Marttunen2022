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
if len(sys.argv) > 1:
  Ihcoeff = float(sys.argv[1])

for icell in range(0,1):
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  v0 = -80
  ca0 = 0.0001

  proximalpoint = 400
  distalpoint = 620
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
st1 = new IClamp(0.5)
L5PC.soma st1

objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
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
L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
L5PC.soma cvode.record(&cai(0.5),casoma,tvec)
""")

  #Block or amplify Ih channels:
  h("""forall if(ismembrane("Ih")) { gIhbar_Ih = """+str(Ihcoeff)+"""*gIhbar_Ih } """)

  Is = [0.1*x for x in range(0,16)]
  spikfreqs = len(Is)*[0]
  for iI in range(0,len(Is)):
    squareAmp = Is[iI]
    squareDur = 15800
    tstop = 16000
    h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
st1.del = 200
""")
    h.init()
    h.run()

    times=np.array(h.tvec)
    Vsoma=np.array(h.vsoma)
    Vdend=np.array(h.vdend)
    Casoma=np.array(h.casoma)
    Cadend=np.array(h.cadend)
    spikes = mytools.spike_times(times,Vsoma,-35,100)
    spikfreqs[iI] = sum([1 for x in spikes if x >= 500.0])/15.5

    if abs(Is[iI]-1.0) < 0.0001:
      Vsoma_control = Vsoma
      Casoma_control = Casoma
      Vdend_control = Vdend
      Cadend_control = Cadend
      times_control = times
      spikes_control = spikes

    picklelist = [spikfreqs,spikes,Vsoma,Vdend,Casoma,Cadend]
    file = open('fI_Ihcoeff'+str(Ihcoeff)+'_tmp.sav', 'wb')
    pickle.dump(picklelist,file)
    file.close()

    
  spts = spikes_control[len(spikes_control)-3:len(spikes_control)]
  istart = next((i for i,x in enumerate(times_control) if x > spts[0]))
  iend = next((i for i,x in enumerate(times_control) if x > spts[1]))+4
  nsteps = iend-istart-1
  tdiff = [y-x for x,y in zip(times_control[istart:iend-1],times_control[istart+1:iend])]
  cadiff = [y-x for x,y in zip(Casoma_control[istart:iend-1],Casoma_control[istart+1:iend])]
  caddiff = [y-x for x,y in zip(Cadend_control[istart:iend-1],Cadend_control[istart+1:iend])]
  caderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],cadiff[0:nsteps-1])]
  caderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],cadiff[1:nsteps])]
  caderiv = [(x+y)/2.0 for x,y in zip(caderiv1,caderiv2)]
  cadderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],caddiff[0:nsteps-1])]
  cadderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],caddiff[1:nsteps])]
  cadderiv = [(x+y)/2.0 for x,y in zip(cadderiv1,cadderiv2)]
  vdiff = [y-x for x,y in zip(Vsoma_control[istart:iend-1],Vsoma_control[istart+1:iend])]
  vddiff = [y-x for x,y in zip(Vdend_control[istart:iend-1],Vdend_control[istart+1:iend])]
  vderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],vdiff[0:nsteps-1])]
  vderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],vdiff[1:nsteps])]
  vderiv = [(x+y)/2.0 for x,y in zip(vderiv1,vderiv2)]
  vdderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],vddiff[0:nsteps-1])]
  vdderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],vddiff[1:nsteps])]
  vdderiv = [(x+y)/2.0 for x,y in zip(vdderiv1,vdderiv2)]

  Vsomac = Vsoma_control[istart+1:iend-1]
  Vdendc = Vdend_control[istart+1:iend-1]
  Casomac = Casoma_control[istart+1:iend-1]
  Cadendc = Cadend_control[istart+1:iend-1]
  timesc = times_control[istart+1:iend-1]
  VDerivc = vderiv[:]
  VDcoeff =  mytools.limitcyclescaledv(Vsomac,VDerivc,Vsomac,VDerivc)
  VdDerivc = vdderiv[:]
  VdDcoeff =  mytools.limitcyclescaledv(Vdendc,VdDerivc,Vdendc,VdDerivc)
  CaDerivc = caderiv[:]
  CaDcoeff =  mytools.limitcyclescaledv(Casomac,CaDerivc,Casomac,CaDerivc)
  CadDerivc = cadderiv[:]
  CadDcoeff =  mytools.limitcyclescaledv(Cadendc,CadDerivc,Cadendc,CadDerivc)

  spikfreqsAll.append(spikfreqs[:])
  timescAll.append(timesc[:])
  VsomacAll.append(Vsomac[:])
  VDerivcAll.append(VDerivc[:])
  VDcoeffAll.append(VDcoeff)
  VdendcAll.append(Vdendc[:])
  VdDerivcAll.append(VdDerivc[:])
  VdDcoeffAll.append(VdDcoeff)
  CasomacAll.append(Casomac[:])
  CaDerivcAll.append(CaDerivc[:])
  CaDcoeffAll.append(CaDcoeff)
  CadendcAll.append(Cadendc[:])
  CadDerivcAll.append(CadDerivc[:])
  CadDcoeffAll.append(CadDcoeff)
  times_controlAll.append(times_control[:])
  Vsoma_controlAll.append(Vsoma_control[:])
  Vdend_controlAll.append(Vdend_control[:])
  Casoma_controlAll.append(Casoma_control[:])
  Cadend_controlAll.append(Cadend_control[:])

picklelist = [spikfreqsAll,timescAll,VsomacAll,VDerivcAll,VDcoeffAll,VdendcAll,VdDerivcAll,VdDcoeffAll,CasomacAll,CaDerivcAll,CaDcoeffAll,
              CadendcAll,CadDerivcAll,CadDcoeffAll,times_controlAll,Vsoma_controlAll,Vdend_controlAll,Casoma_controlAll,Cadend_controlAll,Is]
file = open('fI_Ihcoeff'+str(Ihcoeff)+'.sav', 'wb')
pickle.dump(picklelist,file)
file.close()


