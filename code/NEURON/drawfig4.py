
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

f,axarr = subplots(7,1)

for iax in range(0,7):
  axarr[iax].set_position([0.1+0.07*iax,0.1,0.07,0.1])
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()

axarr[0].set_position([0.07,0.78,0.4,0.205])
for ix in range(0,4):
  axarr[1+ix].set_position([0.07+0.11*ix,0.64,0.03,0.06])
  axarr[1+ix].set_xticks([])
for iy in range(0,2):
  axarr[1+4+iy].set_position([0.12,0.34-0.3*iy,0.28,0.28])

#axarr[9].set_position([0.31,0.2,0.2,0.18])
#axarr[10].set_position([0.55,0.2,0.2,0.18])
#axarr[11].set_position([0.79,0.2,0.2,0.18])

#for i in range(9,12):
#  axarr[i].set_visible(False)
  

for icol in range(0,13):
  mychr = hex(15-icol)[2]
  mychr2 = hex(icol)[2]
  #mycol = "#"+mychr+mychr+mychr2+mychr2+mychr2+mychr2
  mycol = "#"+'ff'+mychr2+mychr2+'00'
  mycol2 = "#00"+mychr+mychr+"ff"
  dists = str(25*icol)+'-'+str(25*(icol+1))+' um'
  dists2 = str(100*icol)+'-'+str(100*(icol+1))+' um'

  if icol < 11:
    polygon = Polygon(array([[-670,-620,-620,-670],[-660+71*icol+x for x in [0,0,71,71]]]).T, True)
    p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
    p.set_facecolor(mycol)
    p.set_edgecolor('none')
    axarr[0].add_collection(p)
    axarr[0].text(-610,-660+71*icol+35,dists,va='center',ha='left',fontsize=4)

  polygon = Polygon(array([[1600,1550,1550,1600],[-660+71*icol+x for x in [0,0,71,71]]]).T, True)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  p.set_facecolor(mycol2)
  p.set_edgecolor('none')
  axarr[0].add_collection(p)
  axarr[0].text(1540,-660+71*icol+35,dists2,va='center',ha='right',fontsize=4)

axarr[0].plot([400,600],[200,200],'k-',lw=1.5)
axarr[0].text(500,210,'200 um',va='bottom',ha='center',fontsize=4)
  
unpicklefile = open('modulhcn_almog/morph_multicolor.sav','rb');morphdata = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
for iplotted in range(len(morphdata)-1,-1,-1):
  axarr[0].plot(morphdata[iplotted][1],morphdata[iplotted][0],morphdata[iplotted][2],linewidth=morphdata[iplotted][3],color=morphdata[iplotted][4])

unpicklefile = open('modulhcn_hay/morph_multicolor.sav','rb');morphdata = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
for iplotted in range(len(morphdata)-1,-1,-1):
  axarr[0].plot([0+x for x in morphdata[iplotted][1]], [-500+x for x in morphdata[iplotted][0]],morphdata[iplotted][2],linewidth=morphdata[iplotted][3],color=morphdata[iplotted][4])
axarr[0].axis("equal")
axarr[0].set_ylim([-670,265])
#axarr[0].set_xlabel('$\mu$m',fontsize=6)
#axarr[0].set_ylabel('$\mu$m',fontsize=6,rotation=0)
#axarr[0].spines['top'].set_visible(False)
#axarr[0].spines['right'].set_visible(False)
#axarr[0].get_xaxis().tick_bottom()
#axarr[0].get_yaxis().tick_left()
axarr[0].get_xaxis().set_visible(False)
axarr[0].get_yaxis().set_visible(False)

Ihcoeffs = [0.0,1.0]

#Basal ff, Almog:
cols = ['#0000FF','#000000']
for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  Ithreshs_thisIh = []
  for myseed in range(1,41):
    filename = 'modulhcn_almog/ffthreshs/basalffthreshs_Ihcoeff'+str(Ihcoeff)+'_dend0-356_seed'+str(myseed)+'.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    #print(filename+" loaded")
    Ithreshs_thisIh.append(unpickledlist[0][-1])
  print(str(mean(Ithreshs_thisIh)))

  axarr[1].bar(2-iIhcoeff,1e3*mean(Ithreshs_thisIh),facecolor=cols[iIhcoeff]) #convert from uS to nS
  print('Basal almog Ihcoeff='+str(Ihcoeff)+' thresh = '+str(1e3*mean(Ithreshs_thisIh)))

#Apical ff, Almog:
cols = ['#0000FF','#000000']
for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  Ithreshs_thisIh = []
  for myseed in range(1,41):
    filename = 'modulhcn_almog/ffthreshs/ffthreshs_Ihcoeff'+str(Ihcoeff)+'_apic0-1300_seed'+str(myseed)+'.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    #print(filename+" loaded")
    Ithreshs_thisIh.append(unpickledlist[0][-1])
  print(str(mean(Ithreshs_thisIh)))

  axarr[2].bar(2-iIhcoeff,1e3*mean(Ithreshs_thisIh),facecolor=cols[iIhcoeff]) #convert from uS to nS
  print('Apical almog Ihcoeff='+str(Ihcoeff)+' thresh = '+str(1e3*mean(Ithreshs_thisIh)))


#Basal ff, Hay:
cols = ['#0000FF','#000000']
for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  Ithreshs_thisIh = []
  for myseed in range(1,41):
    filename = 'modulhcn_hay/ffthreshs/basalffthreshs_Ihcoeff'+str(Ihcoeff)+'_dend0-282_seed'+str(myseed)+'.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    #print(filename+" loaded")
    Ithreshs_thisIh.append(unpickledlist[0][-1])
  print(str(mean(Ithreshs_thisIh)))

  axarr[3].bar(2-iIhcoeff,1e3*mean(Ithreshs_thisIh),facecolor=cols[iIhcoeff]) #convert from uS to nS
  print('Basal Hay Ihcoeff='+str(Ihcoeff)+' thresh = '+str(1e3*mean(Ithreshs_thisIh)))

#Apical ff, Hay:
cols = ['#0000FF','#000000']
for iIhcoeff in range(0,len(Ihcoeffs)):
  Ihcoeff = Ihcoeffs[iIhcoeff]
  Ithreshs_thisIh = []
  for myseed in range(1,41):
    filename = 'modulhcn_hay/ffthreshs/ffthreshs_Ihcoeff'+str(Ihcoeff)+'_apic0-1300_seed'+str(myseed)+'.sav'
    if not exists(filename):
      print(filename+" does not exist")
      continue
    unpicklefile = open(filename,'rb');unpickledlist = pickle.load(unpicklefile);unpicklefile.close()
    #print(filename+" loaded")
    Ithreshs_thisIh.append(unpickledlist[0][-1])
  print(str(mean(Ithreshs_thisIh)))

  axarr[4].bar(2-iIhcoeff,1e3*mean(Ithreshs_thisIh),facecolor=cols[iIhcoeff]) #convert from uS to nS
  print('Apical Hay Ihcoeff='+str(Ihcoeff)+' thresh = '+str(1e3*mean(Ithreshs_thisIh)))

titles_bars = ['Hay, basal','Hay, apical','Almog, basal','Almog,apical']
for ix in range(0,4):
  my_xlim = axarr[1+ix].get_xlim()
  my_ylim = axarr[1+ix].get_ylim()
  axarr[1+ix].text(mean(my_xlim),my_ylim[1]*1.29,titles_bars[ix],fontsize=5,ha='center',va='center',clip_on=False)
  axarr[1+ix].set_ylabel('$g$ (nS)',fontsize=5)
  
#Apical, givendists, Almog+Hay:
dist1s = [200,300,400,500,600,700,800,900,1000,1100,1200,1300]
treename = 'apic'

for imodel in range(0,2):
 minmaxes = [inf,-inf]
 coeffs_saved = [[],[]]
 for idist1 in range(0,len(dist1s)):
  dist1 = dist1s[idist1]
  for idist2 in range(0,len(dist1s)):
    dist2 = dist1s[idist2]
    threshEcons_thisdist2 = []
    for iIhcoeff in range(0,2):
      Ihcoeff = Ihcoeffs[iIhcoeff]
      foundOne = 0
      threshEcons = []
      for myseed in range(1,41):
        filename = 'modulhcn_almog/ffthreshs/ffthreshs_Ihcoeff'+str(Ihcoeff)+'_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'
        if imodel == 1:
          filename = 'modulhcn_hay/ffthreshs/ffthreshs_Ihcoeff'+str(Ihcoeff)+'_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'
        if exists(filename):
          unpicklefile = open(filename,'rb')
          unpickledlist = pickle.load(unpicklefile)
          unpicklefile.close()
          threshEcons.append(unpickledlist[0][2])
          foundOne = 1
      #if not foundOne:
      #  print(filename+' does not exist')
      threshEcons_thisdist2.append(threshEcons[:])

      if len(threshEcons) == 0:
        continue
      c = tanh(mean(threshEcons[:])/5e-5)
      minmaxes = [min(minmaxes[0],mean(threshEcons[:])),max(minmaxes[1],(mean(threshEcons[:]) if mean(threshEcons[:]) < 0.0999 else -inf))]

      myhex = hex(255-int(255*c))
      if len(myhex) < 3:
        mycol = '#000000'
      elif len(myhex) < 4:
        mycol = '#'+'0'+myhex[2]+'0'+myhex[2]+'0'+myhex[2]
      elif len(myhex) > 4:
        mycol = '#000000'
      else:
        mycol = '#'+myhex[2:]+myhex[2:]+myhex[2:]

      if iIhcoeff == 0:
        polygon = Polygon(array([[idist1,idist1+1,idist1],[idist2,idist2,idist2+1]]).T, True) #lower left triangle control
      else:
        polygon = Polygon(array([[idist1,idist1+1,idist1+1],[idist2+1,idist2,idist2+1]]).T, True) #upper right triangle Ih blocked
      p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
      p.set_facecolor(mycol); p.set_edgecolor('none')
      axarr[5+imodel].add_collection(p)
    if len(threshEcons_thisdist2) > 1 and len(threshEcons_thisdist2[0]) > 1 and len(threshEcons_thisdist2[1]) > 1:
      s1 = mean(threshEcons_thisdist2[0])
      s2 = mean(threshEcons_thisdist2[1])
      if s1 == s2:
        mycol = '#808080'
      elif s1 > s2: # more spiking with Ih than without
        c = tanh((s1/s2 - 1)/2)
        #print('c = '+str(c))
        myhex = hex(200-int(200*c))
        if len(myhex) < 3:
          mycol = '#000000'
        elif len(myhex) < 4:
          mycol = '#ff'+'0'+myhex[2]+'0'+myhex[2]
        elif len(myhex) > 4:
          mycol = '#ff0000'
        else:
          mycol = '#ff'+myhex[2:]+myhex[2:]
        coeffs_saved[0].append(s1/s2)
      else: # more spiking without Ih than without
        c = tanh((s2/s1 - 1)/2)
        #print('c = '+str(c))
        myhex = hex(200-int(200*c))
        if len(myhex) < 4:
          mycol = '#0'+myhex[2]+'0'+myhex[2]+'ff'
        elif len(myhex) > 4:
          mycol = '#0000ff'
        else:
          mycol = '#'+myhex[2:]+myhex[2:]+'ff'
        coeffs_saved[1].append(s1/s2)
      #print(mycol)
      #print(myhex)
      #print('s1 = '+str(s1)+', s2 = '+str(s2))
      pval = scipy.stats.ranksums(threshEcons_thisdist2[0], threshEcons_thisdist2[1])[1]
      polygon = Polygon(array([[idist2+1,idist2+1,idist2,idist2],[idist1,idist1+1,idist1+1,idist1]]).T, True)
      p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
      p.set_facecolor(mycol); p.set_edgecolor('none')
      axarr[5+imodel].add_collection(p)
      if pval < 0.05/66:
        axarr[5+imodel].plot(idist2+0.5,idist1+0.5,'k*',lw=0.5,mew=0.5,ms=2.0)
 print('minmaxes = '+str(minmaxes))
 print('coeffs_saved = '+str(unique(coeffs_saved[0]))+' and '+str(unique(coeffs_saved[1])))


axnew = f.add_axes([0.026,0.01,0.033,0.59])
axnew2 = f.add_axes([0.407,0.01,0.033,0.59])
zs = [5e-6*i for i in range(0,21)]
for iz in range(0,len(zs)):
  polygon = Polygon(array([[1,2,2,1],[iz,iz,iz+1,iz+1]]).T, True) 
  c = tanh(zs[iz]/5e-5)
  
  myhex = hex(int(255-255*c))
  if len(myhex) < 4:
    mycol = '#'+'0'+myhex[2]+'0'+myhex[2]+'0'+myhex[2]
  elif len(myhex) > 4:
    mycol = '#000000'
  else:
    mycol = '#'+myhex[2:]+myhex[2:]+myhex[2:]
  print('myhex='+myhex)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  p.set_facecolor(mycol); p.set_edgecolor('none')
  axnew.add_collection(p)
  if iz < len(zs)-1:
    axnew.text(0.8,iz+0.5,myscistr(1e3*zs[iz]),fontsize=4,ha='right',va='center') # convert from uS to nS
  else:
    axnew.text(0.8,iz+0.5,'$\geq$0.1',fontsize=4,ha='right',va='center')
axnew.text(0.8,len(zs)+0.1,'(nS)',fontsize=4,ha='right',va='center')
axnew.set_ylim([0,21])
axnew.set_xlim([0,2])
axnew.get_xaxis().set_visible(False)
axnew.get_yaxis().set_visible(False)
for q in ['top','bottom','left','right']:
  axnew.spines[q].set_visible(False)

zs = [1.1]+[1.0+0.5*i for i in range(1,10)]
for iz in range(0,len(zs)):
  polygon = Polygon(array([[1,2,2,1],[iz+1,iz+1,iz+2,iz+2]]).T, True) 
  polygon2 = Polygon(array([[1,2,2,1],[-iz,-iz,-iz-1,-iz-1]]).T, True) 
  c = tanh((zs[iz]-1)/2)
  myhex = hex(200-int(200*c))
  print('myhex='+myhex)
  if len(myhex) < 3:
    mycol = '#000000'
  elif len(myhex) < 4:
    mycol = '#ff'+'0'+myhex[2]+'0'+myhex[2]
  elif len(myhex) > 4:
    mycol = '#ff0000'
  else:
    mycol = '#ff'+myhex[2:]+myhex[2:]

  if len(myhex) < 4:
    mycol2 = '#0'+myhex[2]+'0'+myhex[2]+'ff'
  elif len(myhex) > 4:
    mycol2 = '#0000ff'
  else:
    mycol2 = '#'+myhex[2:]+myhex[2:]+'ff'

  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  p.set_facecolor(mycol); p.set_edgecolor('none')
  axnew2.add_collection(p)
  if iz < len(zs)-1:
    axnew2.text(0.8,iz+1.5,mystr(zs[iz]),fontsize=4,ha='right',va='center')
  else:
    axnew2.text(0.8,iz+1.5,'>5',fontsize=4,ha='right',va='center')

  p = PatchCollection([polygon2], cmap=matplotlib.cm.jet)
  p.set_facecolor(mycol2); p.set_edgecolor('none')
  axnew2.add_collection(p)
  if iz < len(zs)-1:
    axnew2.text(0.8,-iz-0.5,mystr(zs[iz]),fontsize=4,ha='right',va='center')
  else:
    axnew2.text(0.8,-iz-0.5,'>5',fontsize=4,ha='right',va='center')

polygon = Polygon(array([[1,2,2,1],[0,0,1,1]]).T, True) 
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#808080'); p.set_edgecolor('none')
axnew2.add_collection(p)
axnew2.text(0.5,0.5,'n/a',fontsize=4,ha='center',va='center')

axnew2.set_ylim([-11,12])
axnew2.set_xlim([0,2])
axnew2.get_xaxis().set_visible(False)
axnew2.get_yaxis().set_visible(False)
for q in ['top','bottom','left','right']:
  axnew2.spines[q].set_visible(False)

axnew2.text(2.2,2.5,'Threshold larger for $I_h$-blocked',fontsize=4.5,color='#FF0000',ha='left',va='bottom',rotation=90,clip_on=False)
axnew2.text(2.2,-1.0,'Threshold smaller for $I_h$-blocked',fontsize=4.5,color='#0000FF',ha='left',va='top',rotation=90,clip_on=False)

for iax in range(5,7):
  axarr[iax].set_xticks([0.5+x for x in range(0,len(dist1s))])
  axarr[iax].set_xticklabels([str(d) for d in dist1s],fontsize=4,rotation=20)
  axarr[iax].set_yticks([0.5+x for x in range(0,len(dist1s))])
  axarr[iax].set_yticklabels([str(d) for d in dist1s],fontsize=4,rotation=20)

  #for ix in range(0,2):
  #  axarr[0,ix].set_xticks([0.5+x for x in range(0,len(dist1s))])
  #  axarr[0,ix].set_xticklabels([str(d) for d in dist1s],fontsize=4)
  #  axarr[0,ix].set_yticks([0.5+x for x in range(0,len(dist2s))])
  #  axarr[0,ix].set_yticklabels([str(d) for d in dist2s],fontsize=4)
  #  axarr[0,ix].set_xlim([0,len(dist1s)])
  #  axarr[0,ix].set_ylim([0,len(dist2s)])

for iax in range(0,7):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.05, pos.y1 - 0.01 + 0.015*(iax > 0 and iax < 5), chr(ord('A')+iax), fontsize=9)

f.savefig("fig4.eps")
