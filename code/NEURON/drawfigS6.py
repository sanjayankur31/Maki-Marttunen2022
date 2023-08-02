from pylab import *
import scipy.io
import pickle
from os.path import exists
from matplotlib.collections import PatchCollection
import scipy.stats

def mystr(x):
  mystr = str(x)
  if '000000' in mystr or '999999' in mystr:
    for q in range(7,0,-1):
      mystr = ('{:.'+str(q)+'f}').format(x)
      ilast = len(mystr)-1
      if 'e' in mystr:
        ilast = mystr.find('e') - 1
      if mystr[ilast] != '0':
        return mystr
  return mystr
def myscistr(x):
  mystr = '{:e}'.format(x)
  if '0000' in mystr or '9999' in mystr:
    for q in range(7,0,-1):
      mystr = ('{:.'+str(q)+'e}').format(x)
      ilast = mystr.find('e') - 1
      if mystr[ilast] != '0':
        return mystr
  return mystr

def discontlog(ax,x,y,width=18,height=1.05,col='#000000'):
  ax.plot([x-width/2,x+width/2],[y/sqrt(height),y*sqrt(height)],'k-',lw=2.3,clip_on=False,color=col)
  ax.plot([x-width,x+width],[y/height,y*height],'w-',lw=1.0,clip_on=False,zorder=100)

f,axarr = subplots(6,1)

for iax in range(0,6):
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()

for ix in range(0,3):
  axarr[ix].set_xlabel('')
for iy in range(0,2):
  for ix in range(0,3):
    axarr[3*iy+ix].set_position([0.09+0.313*ix,0.62-0.37*iy,0.28,0.33])
    axarr[3*iy+ix].set_ylim([1,200])
  

dists = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
styles = ['bx-','kx-']

for imodul in range(0,2):
  txt1stmod = 'DA' if imodul == 0 else 'ACh'
  txt2ndmod = 'ACh' if imodul == 0 else 'DA'
  #Hay Thresholds, proximal DA (+10 mV) or ACh (-10 mV)
  moduldv = 10.0-20*imodul
  Ithreshs_all = []
  for idist in range(0,len(dists)):
    dist = dists[idist]
    Ithreshs_both = []
    for iIhcoeff in range(0,3):
      filename = '../modulhcn_hay/strongdendthresh'+str(dist)+'_Ihcoeff1.0.sav'
      if iIhcoeff == 0:
        filename = '../modulhcn_hay/strongdendthresh'+str(dist)+'_Ihmod2ways'+str(moduldv)+',0.0,500_absbound.sav'
      elif iIhcoeff == 1:
        filename = '../modulhcn_hay/strongdendthresh'+str(dist)+'_Ihmod2ways'+str(moduldv)+','+str(-moduldv)+',500_absbound.sav'
      if not exists(filename):
        print(filename+" does not exist")
        continue
      unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
      print(filename+" loaded")
      Ithreshs_both.append(unpickledlist[0][-1])
    Ithreshs_all.append(Ithreshs_both[:])
  try:
   axarr[3*imodul].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], 'bx-',lw=0.5,mew=0.5,ms=2.0)
   axarr[3*imodul].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], 'mx-',lw=0.5,mew=0.5,ms=2.0)
   axarr[3*imodul].semilogy(dists, [Ithreshs_all[i][2] for i in range(0,len(Ithreshs_all))], 'kx--',lw=0.5,mew=0.5,ms=2.0)
  except:
   print("except")

  #Almog Thresholds, proximal DA (+5 mV) or ACh (-5 mV)
  moduldv = 5.0-10*imodul
  Ithreshs_all = []
  for idist in range(0,len(dists)):
    dist = dists[idist]
    Ithreshs_both = []
    for iIhcoeff in range(0,3):
      filename = 'strongdendthresh'+str(dist)+'_Ihcoeff1.0.sav'
      if iIhcoeff == 0:
        filename = 'strongdendthresh'+str(dist)+'_Ihmod2ways'+str(moduldv)+',0.0,800_absbound.sav'
      elif iIhcoeff == 1:
        filename = 'strongdendthresh'+str(dist)+'_Ihmod2ways'+str(moduldv)+','+str(-moduldv)+',800_absbound.sav'
      if not exists(filename):
        print(filename+" does not exist")
        continue
      unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
      print(filename+" loaded")
      Ithreshs_both.append(unpickledlist[0][-1])
    Ithreshs_all.append(Ithreshs_both[:])

  axarr[3*imodul+1].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], 'bx-',lw=0.5,mew=0.5,ms=2.0)
  axarr[3*imodul+1].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], 'mx-',lw=0.5,mew=0.5,ms=2.0)
  axarr[3*imodul+1].semilogy(dists, [Ithreshs_all[i][2] for i in range(0,len(Ithreshs_all))], 'kx--',lw=0.5,mew=0.5,ms=2.0)
  
  #Almog Thresholds, proximal DA (+5 mV) or ACh (-5 mV)
  moduldv = 5.0-10*imodul
  Ithreshs_all = []
  for idist in range(0,len(dists)):
    dist = dists[idist]
    Ithreshs_both = []
    for iIhcoeff in range(0,3):
      filename = 'strongdendthresh'+str(dist)+'_Ihcoeff1.0_apicalCaLVAHay_0.003_100.0_dists585-985_absbound.sav'
      if iIhcoeff == 0:
        filename = 'strongdendthresh'+str(dist)+'_Ihmod2ways'+str(moduldv)+',0.0,800_apicalCaLVAHay_0.003_100.0_dists585-985_absbound.sav'
      elif iIhcoeff == 1:
        filename = 'strongdendthresh'+str(dist)+'_Ihmod2ways'+str(moduldv)+','+str(-moduldv)+',800_apicalCaLVAHay_0.003_100.0_dists585-985_absbound.sav'
      if not exists(filename):
        print(filename+" does not exist")
        continue
      unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
      print(filename+" loaded")
      Ithreshs_both.append(unpickledlist[0][-1])
    Ithreshs_all.append(Ithreshs_both[:])

  axarr[3*imodul+2].semilogy(dists, [Ithreshs_all[i][2] for i in range(0,len(Ithreshs_all))], 'kx--',lw=0.5,mew=0.5,ms=2.0,label='Control',zorder=5)
  axarr[3*imodul+2].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], 'bx-',lw=0.5,mew=0.5,ms=2.0,label=txt1stmod+' at proximal dendrite')
  axarr[3*imodul+2].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], 'mx-',lw=0.5,mew=0.5,ms=2.0,label=txt1stmod+' at prox., '+txt2ndmod+' at dist. dend.')

for ix in range(0,3):
  axarr[ix].set_xticklabels([])
  axarr[3+ix].set_xlabel('distance ($\mu$m)',fontsize=6)
for iy in range(0,2):
  axarr[3*iy].set_ylabel('$I$ (nA)',fontsize=6)
  for ix in range(1,3):
    axarr[3*iy+ix].set_yticklabels([])

axarr[0].set_title('Hay model',fontsize=6)
axarr[1].set_title('Almog model',fontsize=6)
axarr[2].set_title('Almog model with hot zone of Ca$^{2+}$ channels    ',fontsize=6)

#axarr[0].set_ylabel('DA at prox. dend.',fontsize=6)
#axarr[3].set_ylabel('ACh at prox. dend.',fontsize=6)
axarr[2].legend(fontsize=5,edgecolor='none', facecolor=(1, 1, 1, 0.1),loc='upper left',borderpad=-0.25)
axarr[5].legend(fontsize=5,edgecolor='none', facecolor=(1, 1, 1, 0.1),loc='upper left',borderpad=-0.25)

dist1s = [200,300,400,500,600,700,800,900,1000,1100,1200,1300]
treename = 'apic'
patchxs0 = [[0.5-0.45*sin(0.025*x*2*pi) for x in range(0,21)],[0.5-0.45*sin(0.025*x*2*pi) for x in range(20,41)]]
patchys0 = [[0.5+0.45*cos(0.025*x*2*pi) for x in range(0,21)],[0.5+0.45*cos(0.025*x*2*pi) for x in range(20,41)]]
patchxs = patchxs0[0]+[0.5,0,0,1,1,0.5]+patchxs0[1]
patchys = patchys0[0]+[0,0,1,1,0,0]+patchys0[1]

for iax in range(0,6):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.03, pos.y1 - 0.01, chr(ord('A')+iax), fontsize=9)

f.savefig("figS6.eps")


