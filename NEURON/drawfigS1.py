
from pylab import *
import scipy.io
import pickle
from os.path import exists


f,axarr = subplots(12,1)

for iax in range(0,12):
  axarr[iax].set_position([0.1+0.07*iax,0.1,0.07,0.1])
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()

axarr[0].set_position([0.07,0.8,0.2,0.18])
axarr[1].set_position([0.31,0.8,0.2,0.18])
axarr[2].set_position([0.55,0.8,0.2,0.18])
axarr[3].set_position([0.79,0.8,0.2,0.18])

axarr[4].set_position([0.07,0.5,0.2,0.18])
axarr[5].set_position([0.31,0.5,0.2,0.18])
axarr[6].set_position([0.55,0.5,0.2,0.18])
axarr[7].set_position([0.79,0.5,0.2,0.18])

axarr[8].set_position([0.07,0.2,0.2,0.18])
axarr[9].set_position([0.31,0.2,0.2,0.18])
axarr[10].set_position([0.55,0.2,0.2,0.18])
axarr[11].set_position([0.79,0.2,0.2,0.18])

for ix in range(0,2):
  for iy in range(0,2):
    axarr[2*iy+ix].set_position([0.07+0.26*ix,0.4-0.3*iy,0.2,0.18])
  axarr[4+ix].set_position([0.07+0.26*ix,0.7,0.2,0.28])
for i in range(4,12):
  axarr[i].set_visible(False)
  
#Is = [0.69+0.01*x for x in range(0,20)]
#Is_dend = [0.6+0.05*x for x in range(0,20)]
Is1 = [0.4+0.04*x for x in range(0,20)]
Is = [0.4+0.01*x for x in range(0,80)]
Is1_dend = [1.2+0.04*x for x in range(0,20)]
Is_dend = [1.2+0.01*x for x in range(0,80)]

#axarr[0]: fI for Almog, with Ihcoeff 0 and 2
#axarr[2]: fI for Almog, with Ihmod -10 and +10 mV
if exists('fI2_blockdend_Ihcoeff1.0.mat'):
  MAT = scipy.io.loadmat('fI2_blockdend_Ihcoeff1.0.mat')
  axarr[0].plot(Is,MAT['spikfreqs'][0],'k-',lw=0.5)
  axarr[2].plot(Is,MAT['spikfreqs'][0],'k-',lw=0.5)
else:
  print('fI2_blockdend_Ihcoeff1.0.mat does not exist')
  unpicklefile = open('fI_blockdend_Ihcoeff1.0.sav','rb');unpickledlist_control = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
  axarr[0].plot(Is1,unpickledlist_control[0],'k-',lw=0.5)
  axarr[2].plot(Is1,unpickledlist_control[0],'k-',lw=0.5)
MAT = scipy.io.loadmat('fI2_blockdend_Ihcoeff0.5.mat')
axarr[0].plot(Is,MAT['spikfreqs'][0],'b-',lw=0.5)
MAT = scipy.io.loadmat('fI2_blockdend_Ihcoeff1.5.mat')
axarr[0].plot(Is,MAT['spikfreqs'][0],'r-',lw=0.5)

MAT = scipy.io.loadmat('fI2_blockdend_Ihmod-10.0.mat')
axarr[2].plot(Is,MAT['spikfreqs'][0],'b-',lw=0.5)
if exists('fI2_blockdend_Ihmod10.0.mat'):
  MAT = scipy.io.loadmat('fI2_blockdend_Ihmod10.0.mat')
  axarr[2].plot(Is,MAT['spikfreqs'][0],'r-',lw=0.5)
else:
  MAT = scipy.io.loadmat('fI2_blockdend_Ihmod10.0_tmp.mat')
  axarr[2].plot(Is[0:len(MAT['spikfreqs'][0])],MAT['spikfreqs'][0],'r-',lw=0.5)





Is_hay = [0.1*x for x in range(0,16)]

#axarr[1]: fI for Hay, with Ihcoeff 0 and 2
#axarr[3]: fI for Hay, with Ihmod -10 and +10 mV
unpicklefile = open('../modulhcn_hay/fI_blockdend_Ihcoeff1.0.sav','rb');unpickledlist_control_hay = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[1].plot(Is_hay,unpickledlist_control_hay[0][0],'k-',lw=0.5)
unpicklefile = open('../modulhcn_hay/fI_blockdend_Ihcoeff0.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[1].plot(Is_hay,unpickledlist[0][0],'b-',lw=0.5)
unpicklefile = open('../modulhcn_hay/fI_blockdend_Ihcoeff2.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[1].plot(Is_hay,unpickledlist[0][0],'r-',lw=0.5)

Is_hay = [0.1*x for x in range(0,16)]

unpicklefile = open('../modulhcn_hay/fI_blockdend_Ihmod0.0.sav','rb');unpickledlist_control_hay = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[3].plot(Is_hay,unpickledlist_control_hay[0][0],'k-',lw=0.5)
unpicklefile = open('../modulhcn_hay/fI_blockdend_Ihmod-10.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[3].plot(Is_hay,unpickledlist[0][0],'b-',lw=0.5)
unpicklefile = open('../modulhcn_hay/fI_blockdend_Ihmod10.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[3].plot(Is_hay,unpickledlist[0][0],'r-',lw=0.5)


unpicklefile = open('morph.sav','rb');morphdata = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
for iplotted in range(0,len(morphdata)):
  axarr[4].plot(morphdata[iplotted][0],morphdata[iplotted][1],morphdata[iplotted][2],linewidth=morphdata[iplotted][3],color=morphdata[iplotted][4])
axarr[4].axis("equal")

unpicklefile = open('../modulhcn_hay/morph.sav','rb');morphdata = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
for iplotted in range(0,len(morphdata)):
  axarr[5].plot(morphdata[iplotted][0],morphdata[iplotted][1],morphdata[iplotted][2],linewidth=morphdata[iplotted][3],color=morphdata[iplotted][4])
axarr[5].axis("equal")

for i in [0,2]:
  axarr[i].set_ylim([0,60])
for i in range(0,4):
  axarr[i].set_ylabel('$f$ (spikes/s)',fontsize=7)
  axarr[i].set_xlabel('$I$ (nA)',fontsize=7)
                            
for i in range(0,4):
  pos = axarr[i].get_position()
  f.text(pos.x0 - 0.06, pos.y1 - 0.02, chr(ord('A')+i), fontsize=12)
                            



f.savefig("figS1.eps")
