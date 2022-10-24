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
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

Is = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]

Ihcoeff = 1.0
denddist = 620
if len(sys.argv) > 1:
  Ihcoeff = float(sys.argv[1])
if len(sys.argv) > 2:
  denddist = float(sys.argv[2])

I0 = 0
I1 = 50.0

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

double siteVec[2]
sl = new List()
sl=locateSites("apic","""+str(denddist)+""")
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
syn1.e = 0

objref vsoma, vdend, casoma, cadend, ihsoma, ihdend, tvec
vsoma = new Vector()
vdend = new Vector()
casoma = new Vector()
cadend = new Vector()
ihsoma = new Vector()
ihdend = new Vector()
tvec = new Vector()
a_soma cvode.record(&v(0.5),vsoma,tvec)
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
a_soma cvode.record(&cai(0.5),casoma,tvec)
a_soma cvode.record(&ih_iH(0.5),ihsoma,tvec)
apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)
apic[siteVec[0]] cvode.record(&ih_iH(siteVec[1]),ihdend,tvec)
apic[siteVec[0]] st1 = new IClamp(siteVec[1])

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

h("""forall if(ismembrane("iH")) { gbar_iH = """+str(Ihcoeff)+"""*gbar_iH } """)

#thisCa = h.a_soma.cainf_cad
#cai0_ca_ion = """+str(thisCa)+"""


times_all = []
Vsoma_all = []
Vdend_all = []
ihsoma_all = []
ihdend_all = []
Is_this = [I0, I1, (I0+I1)/2]
Niter = 25
Vmaxs_tested = []
Vmaxdends_tested = []
Is_tested = []
for iI in range(0,Niter):
    print("Running "+str(iI)+": I = "+str(Is_this[min(iI,2)]))
    tstop = 1000.0
    squareAmp = Is_this[min(iI,2)]
    squareDur = 0.2
    h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
st1.amp = 0
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(squareAmp)+"""
syn1.onset = 800 + """+str(BACdt)+""" 
  """)
    print("Starting soon...")
    time.sleep(1)
    h.init()
    h.run()
  
    times=np.array(h.tvec)
    Vsoma=np.array(h.vsoma)
    Vdend=np.array(h.vdend)
    Vmaxs_tested.append(max(Vsoma))
    Vmaxdends_tested.append(max(Vdend))
    Is_tested.append(Is_this[min(iI,2)])

    if iI < 2:
      continue
    if max(Vsoma) > -10:
      Is_this = [Is_this[0], Is_this[2], (Is_this[0]+Is_this[2])/2]
    else:
      Is_this = [Is_this[2], Is_this[1], (Is_this[2]+Is_this[1])/2]


    picklelist = [Is_tested,Vmaxs_tested,Vmaxdends_tested]
    file = open('strongdendcondthresh'+str(denddist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav', 'wb')
    pickle.dump(picklelist,file)
    file.close()

