
from pylab import *
import scipy.io
import pickle
from os.path import exists
import mytools

f,axarr = subplots(12,1)

for iax in range(0,12):
  axarr[iax].set_position([0.1+0.07*iax,0.1,0.07,0.1])
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()

for ix in range(0,2):
  for iy in range(0,2):
    axarr[2*iy+ix].set_position([0.07+0.26*ix,0.42-0.31*iy,0.2,0.18])
  axarr[4+ix].set_position([0.07+0.26*ix,0.7,0.2,0.26])
for i in range(6,12):
  axarr[i].set_visible(False)
  
Is1 = [0.4+0.04*x for x in range(0,20)]
Is = [0.4+0.01*x for x in range(0,80)]
Is1_dend = [1.2+0.04*x for x in range(0,20)]
Is_dend = [1.2+0.01*x for x in range(0,80)]

if exists('modulhcn_almog/fI2_Ihcoeff1.0.mat'):
  MAT = scipy.io.loadmat('modulhcn_almog/fI2_Ihcoeff1.0.mat')
  axarr[0].plot(Is,MAT['spikfreqs'][0],'k-',lw=0.5)
  axarr[2].plot(Is,MAT['spikfreqs'][0],'k-',lw=0.5)
else:
  print('modulhcn_almog/fI2_Ihcoeff1.0.mat does not exist')
  unpicklefile = open('modulhcn_almog/fI_Ihcoeff1.0.sav','rb');unpickledlist_control = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
  axarr[0].plot(Is1,unpickledlist_control[0],'k-',lw=0.5)
  axarr[2].plot(Is1,unpickledlist_control[0],'k-',lw=0.5)
MAT = scipy.io.loadmat('modulhcn_almog/fI2_Ihcoeff0.0.mat')
axarr[0].plot(Is,MAT['spikfreqs'][0],'b-',lw=0.5)
MAT = scipy.io.loadmat('modulhcn_almog/fI2_Ihcoeff2.0.mat')
axarr[0].plot(Is,MAT['spikfreqs'][0],'r-',lw=0.5)

MAT = scipy.io.loadmat('modulhcn_almog/fI2_Ihmod-5.0.mat')
axarr[2].plot(Is,MAT['spikfreqs'][0],'b-',lw=0.5)
MAT = scipy.io.loadmat('modulhcn_almog/fI2_Ihmod5.0.mat')
axarr[2].plot(Is,MAT['spikfreqs'][0],'r-',lw=0.5)





Is_hay = [0.1*x for x in range(0,16)]

unpicklefile = open('modulhcn_hay/fI_Ihcoeff1.0.sav','rb');unpickledlist_control_hay = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[1].plot(Is_hay,unpickledlist_control_hay[0][0],'k-',lw=0.5)
unpicklefile = open('modulhcn_hay/fI_Ihcoeff0.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[1].plot(Is_hay,unpickledlist[0][0],'b-',lw=0.5)
unpicklefile = open('modulhcn_hay/fI_Ihcoeff2.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[1].plot(Is_hay,unpickledlist[0][0],'r-',lw=0.5)

Is_hay = [0.1*x for x in range(0,16)]

unpicklefile = open('modulhcn_hay/fI_Ihcoeff1.0.sav','rb');unpickledlist_control_hay = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[3].plot(Is_hay,unpickledlist_control_hay[0][0],'k-',lw=0.5)
unpicklefile = open('modulhcn_hay/fI_Ihmod-10.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[3].plot(Is_hay,unpickledlist[0][0],'b-',lw=0.5)
unpicklefile = open('modulhcn_hay/fI_Ihmod10.0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[3].plot(Is_hay,unpickledlist[0][0],'r-',lw=0.5)


unpicklefile = open('modulhcn_almog/morph.sav','rb');morphdata = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
for iplotted in range(0,len(morphdata)):
  axarr[4].plot(morphdata[iplotted][0],morphdata[iplotted][1],morphdata[iplotted][2],linewidth=morphdata[iplotted][3],color=morphdata[iplotted][4])
axarr[4].axis("equal")

unpicklefile = open('modulhcn_hay/morph.sav','rb');morphdata = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
for iplotted in range(0,len(morphdata)):
  axarr[5].plot(morphdata[iplotted][0],morphdata[iplotted][1],morphdata[iplotted][2],linewidth=morphdata[iplotted][3],color=morphdata[iplotted][4])
axarr[5].axis("equal")

for i in [0,2]:
  axarr[i].set_ylim([0,60])
for i in range(0,4):
  axarr[i].set_ylabel('$f$ (spikes/s)',fontsize=7)
  axarr[i].set_xlabel('$I$ (nA)',fontsize=7)
                            
for i in [4,5]:
  axarr[i].plot([-300-150*(i==4),-300-150*(i==4)],[0,200],'k-')
  axarr[i].text(-300-150*(i==4)-40,100,'200 um',ha='right',va='center',fontsize=6)
  axarr[i].get_xaxis().set_visible(False)
  axarr[i].get_yaxis().set_visible(False)

axnew = []
for i in range(0,4):
  pos = axarr[i].get_position()
  f.text(pos.x0 - 0.06, pos.y1 - 0.02, chr(ord('C')+i), fontsize=12)
  if pos.x0 < 0.2:
    axnew.append(f.add_axes([pos.x0+0.055,pos.y0+0.1,0.02,0.1]))
  else:
    axnew.append(f.add_axes([pos.x0+0.155,pos.y0+0.015,0.02,0.1]))
for i in range(0,2):
  pos = axarr[4+i].get_position()
  f.text(pos.x0 - 0.06, pos.y1 - 0.02, chr(ord('A')+i), fontsize=12)
  
MAT = scipy.io.loadmat('modulhcn_almog/fI2_wait_Ihcoeff1.0_thresh.mat')
axnew[0].bar(0,MAT['Is_tested'][-1],facecolor='#000000')
MAT = scipy.io.loadmat('modulhcn_almog/fI2_wait_Ihcoeff0.0_thresh.mat')
axnew[0].bar(-1,MAT['Is_tested'][-1],facecolor='#0000FF')
MAT = scipy.io.loadmat('modulhcn_almog/fI2_wait_Ihcoeff2.0_thresh.mat')
axnew[0].bar(1,MAT['Is_tested'][-1],facecolor='#FF0000')

unpicklefile = open('modulhcn_hay/fI_wait_Ihcoeff1.0_thresh.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axnew[1].bar(0,unpickledlist[0][-1],facecolor='#000000')
unpicklefile = open('modulhcn_hay/fI_wait_Ihcoeff0.0_thresh.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axnew[1].bar(-1,unpickledlist[0][-1],facecolor='#0000FF')
unpicklefile = open('modulhcn_hay/fI_wait_Ihcoeff2.0_thresh.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axnew[1].bar(1,unpickledlist[0][-1],facecolor='#FF0000')

MAT = scipy.io.loadmat('modulhcn_almog/fI2_wait_Ihcoeff1.0_thresh.mat')
axnew[2].bar(0,MAT['Is_tested'][-1],facecolor='#000000')
MAT = scipy.io.loadmat('modulhcn_almog/fI2_wait_Ihmod-5.0_thresh.mat')
axnew[2].bar(-1,MAT['Is_tested'][-1],facecolor='#0000FF')
MAT = scipy.io.loadmat('modulhcn_almog/fI2_wait_Ihmod5.0_thresh.mat')
axnew[2].bar(1,MAT['Is_tested'][-1],facecolor='#FF0000')

unpicklefile = open('modulhcn_hay/fI_wait_Ihcoeff1.0_thresh.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axnew[3].bar(0,unpickledlist[0][-1],facecolor='#000000')
unpicklefile = open('modulhcn_hay/fI_wait_Ihmod-10.0_thresh.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axnew[3].bar(-1,unpickledlist[0][-1],facecolor='#0000FF')
unpicklefile = open('modulhcn_hay/fI_wait_Ihmod10.0_thresh.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axnew[3].bar(1,unpickledlist[0][-1],facecolor='#FF0000')

axnew[0].set_ylim([0,0.86])
axnew[2].set_ylim([0,0.86])
axnew[1].set_ylim([0,0.43])
axnew[3].set_ylim([0,0.43])

axarr[0].set_xlim([0.5,1.0])
axarr[2].set_xlim([0.5,1.0])
axarr[1].set_xlim([0.25,1.6])
axarr[3].set_xlim([0.25,1.6])

axarr[0].set_ylim([0,33])
axarr[2].set_ylim([0,33])

for iax in range(0,len(axnew)):
  for tick in axnew[iax].xaxis.get_major_ticks() + axnew[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axnew[iax].spines['top'].set_visible(False)
    axnew[iax].spines['right'].set_visible(False)
    axnew[iax].get_xaxis().tick_bottom()
    axnew[iax].get_yaxis().tick_left()
  axnew[iax].set_xticks([])
  axnew[iax].set_ylabel('nA',fontsize=6)

myleg = mytools.mylegend(f,[0.03,0.42+0.18+0.02,0.5,0.04],['b-','k-','r-'],['$I_h$ blocked','Control','$I_h$ over-expressed'],3,2,0.55,0.35,myfontsize=6)
myleg2 = mytools.mylegend(f,[0.03,0.42-0.31+0.18+0.02,0.5,0.04],['b-','k-','r-'],['ACh modulation','Control','DA modulation'],3,2,0.55,0.35,myfontsize=6)
for q in ['top','right','bottom','left']:
  myleg.spines[q].set_visible(False)
  myleg2.spines[q].set_visible(False)

axarr[4].set_title('Almog',fontsize=8)
axarr[5].set_title('Hay',fontsize=8)
f.savefig("fig1.eps")
