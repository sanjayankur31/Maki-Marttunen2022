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
import scipy.io

v0 = -62
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = [0.4+0.01*x for x in range(0,80)]
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

Ihmod = 0.0
if len(sys.argv) > 1:
  Ihmod = float(sys.argv[1])


h("""
load_file("myrun.hoc")
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
tstop = 1000
""")

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

h("""forall if(ismembrane("iH")) { off_iH = off_iH + """+str(Ihmod)+""" } """)

#thisCa = h.a_soma.cainf_cad
#cai0_ca_ion = """+str(thisCa)+"""


spTimesThisCoeff = []
spTimesThisCoeff2 = []
ISIs = len(Is)*[0.0]
nSpikes = []
nSpikes1 = []
spikfreqs = len(Is)*[0]
spikfreqs2 = len(Is)*[0]


Is_this = [0.0,1.0, 0.5]
Niter = 15
Vmaxs_tested = []
Vmaxdends_tested = []
Is_tested = []
for iI in range(0,Niter):
    print(("Running "+str(iI)))
    tstop = 16000.0
    squareAmp = Is_this[min(iI,2)]
    squareDur = 5800.0
    h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 10200
st1.dur = """+str(squareDur)+"""
syn1.gmax = 0
syn1.onset = 10200 + """+str(BACdt)+""" 
  """)
    print("Starting soon...")
    time.sleep(1)
    h.init()
    h.run()
  
    times=np.array(h.tvec)
    Vsoma=np.array(h.vsoma)
    Vdend=np.array(h.vdend)
    Casoma=np.array(h.casoma)
    Cadend=np.array(h.cadend)
    spikes = mytools.spike_times(times,Vsoma,-50,-50)
    spikes2 = mytools.spike_times(times,Vsoma,-50,inf)

    if iI < 2:
      continue
    if len(spikes2) > 0:
      Is_this = [Is_this[0], Is_this[2], (Is_this[0]+Is_this[2])/2]
    else:
      Is_this = [Is_this[2], Is_this[1], (Is_this[2]+Is_this[1])/2]

    Vmaxs_tested.append(max(Vsoma))
    Vmaxdends_tested.append(max(Vdend))
    Is_tested.append(Is_this[min(iI,2)])

    #picklelist = [spikfreqs,spikfreqs2,spikes,spikes2,Vsoma,Vdend,Casoma,Cadend]
    scipy.io.savemat('fI2_wait_Ihmod'+str(Ihmod)+'_thresh.mat',{'Is_tested':Is_tested,'Vmaxs_tested':Vmaxs_tested,'Vmaxdends_tested':Vmaxdends_tested})


  #picklelist = [ISIs_all,spTimesAll,spTimesAll2,MT]
  #file = open('ifcurvesmut.sav', 'w')
  #pickle.dump(picklelist,file)
  #file.close()
  
scipy.io.savemat('fI2_wait_Ihmod'+str(Ihmod)+'_thresh.mat',{'Is_tested':Is_tested,'Vmaxs_tested':Vmaxs_tested,'Vmaxdends_tested':Vmaxdends_tested})
