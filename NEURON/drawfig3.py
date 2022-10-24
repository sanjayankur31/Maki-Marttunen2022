
from pylab import *
import scipy.io
import pickle
from os.path import exists
from matplotlib.collections import PatchCollection
import mytools

def discontlog(ax,x,y,width=18,height=1.05,col='#000000'):
  ax.plot([x-width/2,x+width/2],[y/sqrt(height),y*sqrt(height)],'k-',lw=2.3,clip_on=False,color=col)
  ax.plot([x-width,x+width],[y/height,y*height],'w-',lw=1.0,clip_on=False,zorder=100)
def discontlogopp(ax,x,y,width=18,height=1.05,col='#000000'):
  ax.plot([x-width/2,x+width/2],[y*sqrt(height),y/sqrt(height)],'k-',lw=2.3,clip_on=False,color=col)
  ax.plot([x-width,x+width],[y*height,y/height],'w-',lw=1.0,clip_on=False,zorder=100)


f,axarr = subplots(6,1)

axarr[0].set_position([0.065,0.75,0.15,0.18])
axarr[1].set_position([0.27,0.75,0.15,0.18])
axarr[2].set_position([0.065,0.45,0.15,0.18])
axarr[3].set_position([0.27,0.45,0.15,0.18])
axarr[4].set_position([0.065,0.15,0.15,0.18])
axarr[5].set_position([0.27,0.15,0.15,0.18])


for iax in range(0,5):
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()
  axarr[iax].set_xlim([0,1000])
  pos = axarr[iax].get_position()
  if iax ==2:
    f.text(pos.x0 - 0.035, pos.y1 - 0.025, chr(ord('A')+iax), fontsize=12)
  else:
    f.text(pos.x0 + 0.01, pos.y1 - 0.025, chr(ord('A')+iax), fontsize=12)


Is1 = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 100.0]
Is = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 100.0]

dists = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
Ihcoeffs = [0.0,1.0]
styles = ['bx-','kx-']
styles_timecourse = ['b-','k-']
styles_timecourse2 = ['b--','k--']


#Thresholds, normal
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_almog/strongdendthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    print(filename+" loaded")
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[4].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0, color='#AAAAFF')
axarr[4].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0, color='#CCCCCC')

#Hay Thresholds, normal
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_hay/strongdendthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[0].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0,color='#AAAAFF')
axarr[0].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0,color='#CCCCCC')
axarr[1].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0,color='#AAAAFF')
axarr[1].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0,color='#CCCCCC')

#Hay: Thresholds, cond-based
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_hay/strongdendcondthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[2].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0,color='#AAAAFF')
axarr[2].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0,color='#CCCCCC')


#Hay Thresholds, no CaLVA
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_hay/strongdendthresh'+str(dist)+'_gCa_LVAstbar_Ca_LVAst_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[0].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[0].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[0].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[0].set_ylabel('$I$ (nA)',fontsize=6)
axarr[0].set_title('Hay, no LVA CCs',fontsize=6)

#Hay Thresholds, no LVA hot zone
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_hay/nohotLVA_strongdendthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[1].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[1].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[1].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[1].set_ylabel('$I$ (nA)',fontsize=6)
axarr[1].set_title('Hay, no hot zone of LVA CCs',fontsize=6)

#Hay: Thresholds, cond-based, no hot zone
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_hay/nohotLVA_strongdendcondthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[2].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[2].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[2].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[2].set_ylabel('$g$ ($\mu$S)',fontsize=6)
axarr[2].set_title('Hay, no hot zone of LVA CCs',fontsize=6)

discontlog(axarr[2],690,30,col='#0000FF')
discontlog(axarr[2],810,30,col='#AAAAFF')
discontlog(axarr[2],760,18,col='#000000')
discontlog(axarr[2],0,30,col='#000000')

#Hay: time course of inactivation
iamp = 1 #0.01 nS
for iIhcoeff in range(0,2):
  filename = 'modulhcn_hay/strongdendstimcond800.0_Ihcoeff'+str(Ihcoeffs[iIhcoeff])+'_recCaLVA_Ih.sav'
  unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
  #picklelist = [times_all,Vsoma_all,Vdend_all,ihsoma_all,ihdend_all,h_ih_all,m_calva_all,h_calva_all]
  axarr[3].plot([x-800 for x in unpickledlist[0][iamp]],unpickledlist[7][iamp],styles_timecourse[iIhcoeff],lw=0.5)
  axarr[3].set_xlim([-5,50])
  axarr[3].set_ylim([0,0.6])
axarr[3].set_xlabel('$t$ (ms)',fontsize=6)
axarr[3].set_ylabel('h',fontsize=6)
axarr[3].set_title('Hay, LVA CC inactivation',fontsize=6)


#Thresholds, Almog, added LVA hot zone
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_almog/strongdendthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_apicalCaLVAHay_0.003_100.0_dists585-985_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    print(filename+" loaded")
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[4].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[4].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[4].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[4].set_ylabel('$I$ (nA)',fontsize=6)
axarr[4].set_title('Almog, with hot zone of LVA CCs',fontsize=6)

axarr[5].set_visible(False)

myleg = mytools.mylegend(f,[0.13,0.968,0.25,0.03],['b-','k-'],['$I_h$ blocked','Control'],2,2,0.55,0.35,myfontsize=7)
for q in ['top','right','bottom','left']:
  myleg.spines[q].set_visible(False)

f.savefig("fig3.eps")
