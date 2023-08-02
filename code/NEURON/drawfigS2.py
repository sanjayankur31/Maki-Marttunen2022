
from pylab import *
import scipy.io
import pickle
from os.path import exists
from matplotlib.collections import PatchCollection

def discontlog(ax,x,y,width=18,height=1.05,col='#000000'):
  ax.plot([x-width/2,x+width/2],[y/sqrt(height),y*sqrt(height)],'k-',lw=2.3,clip_on=False,color=col)
  ax.plot([x-width,x+width],[y/height,y*height],'w-',lw=1.0,clip_on=False,zorder=100)
def discontlogopp(ax,x,y,width=18,height=1.05,col='#000000'):
  ax.plot([x-width/2,x+width/2],[y*sqrt(height),y/sqrt(height)],'k-',lw=2.3,clip_on=False,color=col)
  ax.plot([x-width,x+width],[y*height,y/height],'w-',lw=1.0,clip_on=False,zorder=100)


f,axs = subplots(3,3)
axarr = sum([axs[i].tolist() for i in range(0,len(axs))]+[[]])

for iax in range(0,3):
  for iay in range(0,3):
    axs[iay,iax].set_position([0.07+0.32*iax,0.07+0.32*(2-iay),0.26,0.21])
    for tick in axs[iay,iax].xaxis.get_major_ticks() + axs[iay,iax].yaxis.get_major_ticks():
      tick.label.set_fontsize(3.5)
    axs[iay,iax].spines['top'].set_visible(False)
    axs[iay,iax].spines['right'].set_visible(False)
    axs[iay,iax].get_xaxis().tick_bottom()
    axs[iay,iax].get_yaxis().tick_left()
    axs[iay,iax].set_xlim([0,1000])
    axs[iay,iax].set_ylim([2,110])
    if iax == 0:
      axs[iay,iax].set_ylabel('$I$ (nA)',fontsize=6)
    if iay == 2:
      axs[iay,iax].set_xlabel('distance ($\mu$m)',fontsize=6)

Is1 = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 100.0]
Is = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 100.0]

#dists = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
dists = [100.0,  200.0,  300.0,  400.0,  500.0,  600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0]
Ihcoeffs = [0.0,1.0]
styles = ['bx-','kx-']

additionalblockeds = ['gCa_LVAstbar_Ca_LVAst','gCa_HVAbar_Ca_HVA','gImbar_Im','gK_Pstbar_K_Pst','gK_Tstbar_K_Tst','gNap_Et2bar_Nap_Et2','gSK_E2bar_SK_E2','gSKv3_1bar_SKv3_1','gNaTa_tbar_NaTa_t']
additionalblockedtitles = ['LVA Ca$^{2+}$ channel','HVA Ca$^{2+}$ channel','M-type K$^{+}$ channel','Persistent K$^{+}$ channel','Transient K$^{+}$ channel','Persistent Na$^{+}$ channel','SK channel','Kv3.1 K$^{+}$ channel','Fast Na$^{+}$ channel']
for iblocked in range(0,len(additionalblockeds)):
  #Hay Thresholds, normal
  Ithreshs_all = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    Ithreshs = []
    dists_saved = []
    for idist in range(0,len(dists)):
      dist = dists[idist]
      filename = '../modulhcn_hay/strongdendthresh'+str(dist)+'_'+additionalblockeds[iblocked]+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
      if not exists(filename):
        print(filename+" does not exist")
        continue
      unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
      print(filename+" loaded")
      Ithreshs.append(unpickledlist[0][-1])
      dists_saved.append(dist)
    Ithreshs_all.append(Ithreshs[:])
    axarr[iblocked].semilogy(dists_saved, Ithreshs, styles[iIhcoeff],lw=0.5,mew=0.5,ms=2.0)

  axarr[iblocked].set_title(additionalblockedtitles[iblocked],fontsize=6)

discontlogopp(axarr[8],150,85,width=18,height=1.05,col='#000000')
discontlog(axarr[8],0,85,width=18,height=1.05,col='#000000')
f.savefig("figS2.eps")
