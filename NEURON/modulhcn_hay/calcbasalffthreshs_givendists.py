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
myseed = 1
dist1 = 685
dist2 = 885
treename = 'dend'
if len(sys.argv) > 1:
  Ihcoeff = float(sys.argv[1])
if len(sys.argv) > 2:
  dist1 = int(float(sys.argv[2]))
if len(sys.argv) > 3:
  dist2 = int(float(sys.argv[3]))
if len(sys.argv) > 4:
  myseed = int(float(sys.argv[4]))

Econs = [0.0,0.1,0.05]
threshs = [-40,-35,-30,-25,-20]
thresh = -35
Niter = 30
for icell in range(0,1):
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate_withsyns_nogaba.hoc"
  v0 = -80
  ca0 = 0.0001

  proximalpoint = 400
  distalpoint = 620
  BACdt = 5.0
  tstop = 1100
  T_stim = 1000

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
tstop = """+str(tstop)+"""
""")

  #Block or amplify Ih channels:
  h("""forall if(ismembrane("Ih")) { gIhbar_Ih = """+str(Ihcoeff)+"""*gIhbar_Ih } """)

  h("""
  NsynE = 2000
  rdSeed = """+str(myseed)+"""
  L5PC.initRand(rdSeed)
  L5PC.setparameters(0.0,NsynE)
  objref preTrainListLocal
  objref rds1,rds2
  preTrainListLocal = new List()
  
  L5PC.distributeSynGivenDist(\""""+treename+"""\","""+str(dist1)+""","""+str(dist2)+""")
  for(i2=0;i2<(NsynE);i2+=1){
    {preTrainListLocal.append(new Vector())}
    {preTrainListLocal.o[preTrainListLocal.count()-1].append("""+str(T_stim)+""")}
    }
  L5PC.setpretrains(preTrainListLocal)
  L5PC.queuePreTrains()
  """)

  spikes_all = []
  Vmaxes_all = []
  Econs_all = []
  for iiter in range(0,Niter):
    h("""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = 0
st1.dur = 0
  """)

    print("Econ = "+str(Econs[min(2,iiter)]))
    h("""
    L5PC.setparameters("""+str(Econs[min(2,iiter)])+""",NsynE)
    for(syni=0;syni<NsynE;syni+=1) { 
      L5PC.synlist.o[syni].gAMPAmax = """+str(Econs[min(2,iiter)])+"""
      L5PC.synlist.o[syni].gNMDAmax = """+str(Econs[min(2,iiter)])+"""
    }
""")
    h.init()
    h.run()

    times=np.array(h.tvec)
    Vsoma=np.array(h.vsoma)
    Vdend=np.array(h.vdend)
    Casoma=np.array(h.casoma)
    Cadend=np.array(h.cadend)
    spikes = []
    for ithresh in range(0,len(threshs)):
      spikes.append(mytools.spike_times(times,Vsoma,threshs[ithresh],100))
    Vmaxes_all.append(max(Vsoma))
    spikes_all.append(spikes[:])
    Econs_all.append(Econs[min(2,iiter)])
    if iiter > 1 and max(Vsoma) > thresh:
      Econs = [Econs[0],Econs[2],0.5*(Econs[0]+Econs[2])]
    elif iiter > 1:
      Econs = [Econs[2],Econs[1],0.5*(Econs[2]+Econs[1])]

  picklelist = [Econs,spikes_all,Vmaxes_all,Econs_all] #,times,Vsoma,Vdend,Casoma,Cadend]
  file = open('ffthreshs/basalffthreshs_Ihcoeff'+str(Ihcoeff)+'_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav', 'wb')
  pickle.dump(picklelist,file)
  file.close()

    
