#cp ../almogmod2/calcifcurves2.py calcifcurves.py
from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys
from setparams import *
from os.path import exists
import time

v0 = -62
ca0 = 0.0001
BACdt = 5.0
Is = [0.4+0.04*x for x in range(0,20)]
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

myseed = 1
Ihcoeff = 1.0
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

h("""
load_file("myrun_withsyns_nogaba.hoc")
""")
h("""
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma
objref st1,syn1, sl
a_soma st1 = new IClamp(0.5)

double siteVec[2]
sl = new List()
sl=locateSites("apic",620)
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
apic[siteVec[0]] syn1 = new AlphaSynapse(siteVec[1])
//apic[41] syn1 = new AlphaSynapse(0.5)

syn1.onset = 3400
syn1.tau = 3
syn1.gmax = 0.0
syn1.e = 50

objref vsoma, vdend, casoma, cadend, tvec
vsoma = new Vector()
vdend = new Vector()
casoma = new Vector()
cadend = new Vector()
tvec = new Vector()
a_soma cvode.record(&v(0.5),vsoma,tvec)
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
a_soma cvode.record(&cai(0.5),casoma,tvec)
apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)

v_init = -62
dt = 0.025
tstop = 1100
""")

tstop = 1100
T_stim = 1000

styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
#cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  

paramdicts = []
paramdicts.append({})                                                                                                                                         # 4-6 spikes per burst, control
paramdicts.append({'transvec.x(31)': 1.25, 'transvec.x(32)': 1.25})                                                                                           # 4-5 spikes per burst         
paramdicts.append({'transvec.x(31)': 1.5, 'transvec.x(32)': 1.5})                                                                                             # 3-4 spikes per burst         
paramdicts.append({'transvec.x(31)': 2.0, 'transvec.x(32)': 2.0})                                                                                             # 3-4 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0})                                                                                             # 2-3 spikes per burst       
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.3, 'transvec.x(21)': 1.3, 'transvec.x(25)': 1.3, 'transvec.x(26)': 1.3}) # 2 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.6, 'transvec.x(21)': 1.6, 'transvec.x(25)': 1.6, 'transvec.x(26)': 1.6}) # 1-2 spikes per burst     


icell = 0

paramdict = paramdicts[icell]
print("Setting params...")
setparams(paramdict)

h("""forall if(ismembrane("iH")) { gbar_iH = """+str(Ihcoeff)+"""*gbar_iH } """)


#Add noisy inputs to the tip and give activation timings.
h("""
NsynE = 2000
rdSeed = """+str(myseed)+"""
initRand(rdSeed)
setparameters(0.0,NsynE)
objref preTrainListLocal
objref rds1,rds2
preTrainListLocal = new List()
{rds1 = new Random(1000*rdSeed)}

distributeSynGivenDist(\""""+treename+"""\","""+str(dist1)+""","""+str(dist2)+""")
for(i2=0;i2<(NsynE);i2+=1){
  {preTrainListLocal.append(new Vector())}
  {preTrainListLocal.o[preTrainListLocal.count()-1].append("""+str(T_stim)+""")}
}
setpretrains(preTrainListLocal)
queuePreTrains()
""")


#thisCa = h.a_soma.cainf_cad
#cai0_ca_ion = """+str(thisCa)+"""


spTimesThisCoeff = []
spTimesThisCoeff2 = []
ISIs = len(Is)*[0.0]
nSpikes = []
nSpikes1 = []
spikfreqs = len(Is)*[0]
spikfreqs2 = len(Is)*[0]

spikes_all = []
Vmaxes_all = []
Econs_all = []
for iiter in range(0,Niter):
  h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
st1.amp = 0
st1.del = 0
st1.dur = 0
syn1.gmax = 0
syn1.onset = 200 + """+str(BACdt)+""" 
  """)

  print("Econ = "+str(Econs[min(2,iiter)]))
  h("""
  setparameters("""+str(Econs[min(2,iiter)])+""",NsynE)
  for(syni=0;syni<NsynE;syni+=1) { 
    synlist.o[syni].gAMPAmax = """+str(Econs[min(2,iiter)])+"""
    synlist.o[syni].gNMDAmax = """+str(Econs[min(2,iiter)])+"""
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
