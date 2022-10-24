
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


f,axarr = subplots(9,1)

for iax in range(0,9):
  axarr[iax].set_position([0.1+0.07*iax,0.1,0.07,0.1])
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()

axarr[0].set_position([0.055,0.76,0.3,0.18])
axarr[1].set_position([0.415,0.76,0.15,0.18])
axarr[2].set_position([0.625,0.76,0.15,0.18])
axarr[3].set_position([0.835,0.76,0.15,0.18])

axarr[4].set_position([0.055,0.5,0.15,0.18])
axarr[5].set_position([0.205,0.5,0.15,0.18])
axarr[6].set_position([0.415,0.5,0.15,0.18])
axarr[7].set_position([0.625,0.5,0.15,0.18])
axarr[8].set_position([0.835,0.5,0.15,0.18])


for iax in [1,2,3,6,7,8]:
  axarr[iax].set_xlim([0,1000])
for iax in range(0,9):
  pos = axarr[iax].get_position()
  if iax != 5 and iax != 9 and iax != 10:
    f.text(pos.x0 - 0.035, pos.y1 - 0.025, chr(ord('A')+iax), fontsize=12)
  else:
    f.text(pos.x0 + 0.01, pos.y1 - 0.025, chr(ord('A')+iax), fontsize=12)

Is1 = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 100.0]
Is = [1.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 100.0]

dists = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
#dists = [100.0,  200.0,  300.0,  400.0,  500.0,  600.0, 700.0, 800.0]
Ihcoeffs = [0.0,1.0]
styles = ['bx-','kx-']

polygon = Polygon(array([[800,800,800.2,800.2],[140*2+20,-85,-85,140*2+20]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#FFBBBB')
p.set_edgecolor('none')
axarr[0].add_collection(p)

for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  filename = 'modulhcn_almog/strongdendstim500.0_Ihcoeff'+str(Ihcoeff)+'.sav'
  if not exists(filename):
    print(filename+" does not exist")
  else:
    print(filename+" loaded")
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    times_all = unpickledlist[0]
    Vsoma_all = unpickledlist[1]
    iIs = [4,7,8]
    for iiI in range(0,len(iIs)):
      iI = iIs[iiI]
      axarr[0].plot(times_all[iI], [140*iiI+x for x in Vsoma_all[iI]], styles[iIhcoeff].replace('x',''),lw=0.3)
      if iIhcoeff == 0:
        axarr[0].text(820,140*iiI-30,'$I$ = '+str(Is[iI])+' nA',color='#FF0000',fontsize=6)
    axarr[0].set_xlim([796,864])
    axarr[0].set_ylim([-85,140*2+20])
    axarr[0].plot([851,851,856],[140*2-0,140*2-50,140*2-50],'k-',lw=0.7)
    axarr[0].set_xticks([])
    axarr[0].set_yticks([])
axarr[0].set_ylabel('$\mathbf{Almog}$\n\n$V_m$',fontsize=6)

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

axarr[1].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[1].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[1].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[1].set_ylabel('$I$ (nA)',fontsize=6)
#TODO: check this
discontlogopp(axarr[1],925,450,col='#0000FF')
discontlogopp(axarr[1],0,450,col='#000000')
axarr[1].set_ylim([2,550])

#Membrane potential maxima 
for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  Vdendmax_all = []
  for idist in range(0,len(dists)):
    dist = dists[idist]
    filename = 'modulhcn_almog/strongdendthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    Vdend_all = unpickledlist[2]
    Vdendmax_all.append(Vdend_all[-1])
    print(str(std(Vdend_all[-6:])))
    if std(Vdend_all[-6:]) > 0.1:
      print("  "+str(Vdend_all[-6:]))
        
  axarr[2].plot(dists, Vdendmax_all,styles[iIhcoeff],lw=0.5,mew=0.5,ms=2.0)
  axarr[2].set_ylim([-60,140])
axarr[2].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[2].set_ylabel('$V_{\mathrm{max}}$ (mV)',fontsize=6)

#Almog: Thresholds, cond-based
Ithreshs_all = []
for idist in range(0,len(dists)):
  dist = dists[idist]
  Ithreshs_both = []
  for iIhcoeff in range(0,len(Ihcoeffs)):
    Ihcoeff = Ihcoeffs[iIhcoeff]
    filename = 'modulhcn_almog/strongdendcondthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    print(filename+" loaded")
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    Ithreshs_both.append(unpickledlist[0][-1])
  Ithreshs_all.append(Ithreshs_both[:])

axarr[3].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[3].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[3].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[3].set_ylabel('$g$ ($\mu$S)',fontsize=6)
#TODO: check this
discontlog(axarr[3],630,30,col='#0000FF')
discontlog(axarr[3],690,30,col='#000000')
discontlog(axarr[3],0,30,col='#000000')
  

#Hay:
polygon = Polygon(array([[800,800,800.2,800.2],[140*2+60,-95,-95,140*2+60]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#FFBBBB')
p.set_edgecolor('none')
axarr[4].add_collection(p)
polygon = Polygon(array([[800,800,800.2,800.2],[140*2+60,-95,-95,140*2+60]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#FFBBBB')
p.set_edgecolor('none')
axarr[5].add_collection(p)
for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  filename = 'modulhcn_hay/strongdendstim500.0_Ihcoeff'+str(Ihcoeff)+'.sav'
  if not exists(filename):
    print(filename+" does not exist")
  else:
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    times_all = unpickledlist[0]
    Vsoma_all = unpickledlist[1]
    for iiI in range(0,len(iIs)):
      iI = iIs[iiI]
      axarr[4].plot(times_all[iI], [140*iiI+x for x in Vsoma_all[iI]], styles[iIhcoeff].replace('x',''),lw=0.3)
      if iIhcoeff == 0:
        axarr[4].text(816,140*iiI+20,str(Is[iI])+' nA',color='#FF0000',fontsize=6)
    axarr[4].set_xlim([796,830])
    axarr[4].set_ylim([-95,140*2+60])
    #axarr[4].set_xticks([800,820])
    axarr[4].plot([798,798,803],[140*0+20,140*0-30,140*0-30],'k-',lw=0.7)
    axarr[4].set_xticks([])
    axarr[4].set_yticks([])

  filename = 'modulhcn_hay/strongdendstim800.0_Ihcoeff'+str(Ihcoeff)+'.sav'
  if not exists(filename):
    print(filename+" does not exist")
  else:
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    times_all = unpickledlist[0]
    Vsoma_all = unpickledlist[1]
    for iiI in range(0,len(iIs)):
      iI = iIs[iiI]
      axarr[5].plot(times_all[iI], [140*iiI+x for x in Vsoma_all[iI]], styles[iIhcoeff].replace('x',''),lw=0.3)
      if iIhcoeff == 0:
        axarr[5].text(816,140*iiI+20,str(Is[iI])+' nA',color='#FF0000',fontsize=6)
    axarr[5].set_xlim([796,830])
    axarr[5].set_ylim([-95,140*2+60])
    #axarr[5].set_xticks([800,820])
    axarr[5].set_yticklabels([])
    axarr[5].plot([798,798,803],[140*0+20,140*0-30,140*0-30],'k-',lw=0.7)
    axarr[5].set_xticks([])
    axarr[5].set_yticks([])

axarr[4].set_ylabel('$\mathbf{Hay}$\n\n$V_m$',fontsize=6)


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

axarr[6].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[6].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[6].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[6].set_ylabel('$I$ (nA)',fontsize=6)

#Membrane potential maxima 
for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  Vdendmax_all = []
  for idist in range(0,len(dists)):
    dist = dists[idist]
    filename = 'modulhcn_hay/strongdendthresh'+str(dist)+'_Ihcoeff'+str(Ihcoeff)+'_absbound.sav'
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    print(filename+" loaded")
    Vdend_all = unpickledlist[2]
    Vdendmax_all.append(Vdend_all[-1])
    print(str(std(Vdend_all[-6:])))
    if std(Vdend_all[-6:]) > 0.1:
      print("  "+str(Vdend_all[-6:]))
        
  axarr[7].plot(dists, Vdendmax_all,styles[iIhcoeff],lw=0.5,mew=0.5,ms=2.0)
  axarr[7].set_ylim([-60,140])    
axarr[7].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[7].set_ylabel('$V_{\mathrm{max}}$ (mV)',fontsize=6)
      

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

axarr[8].semilogy(dists, [Ithreshs_all[i][0] for i in range(0,len(Ithreshs_all))], styles[0],lw=0.5,mew=0.5,ms=2.0)
axarr[8].semilogy(dists, [Ithreshs_all[i][1] for i in range(0,len(Ithreshs_all))], styles[1],lw=0.5,mew=0.5,ms=2.0)
axarr[8].set_xlabel('distance ($\mu$m)',fontsize=6)
axarr[8].set_ylabel('$g$ ($\mu$S)',fontsize=6)
#TODO: check this
discontlog(axarr[8],690,30,col='#0000FF')
discontlog(axarr[8],810,30,col='#0000FF')
discontlog(axarr[8],760,18,col='#000000')
discontlog(axarr[8],0,30,col='#000000')

myleg = mytools.mylegend(f,[0.33,0.955,0.4,0.04],['b-','k-'],['$I_h$ blocked','Control'],2,2,0.55,0.35,myfontsize=7)
for q in ['top','right','bottom','left']:
  myleg.spines[q].set_visible(False)


f.savefig("fig2.eps")
