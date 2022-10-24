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

f,axarr = subplots(4,1)

for iax in range(0,4):
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()

#for ix in range(0,3):
#  axarr[ix].set_xlabel('')
for ix2 in range(0,2):
  for ix in range(0,2):
    axarr[ix+2*ix2].set_position([0.15+0.35*ix2,0.64-0.36*ix,0.29,0.34])


dist1s = [200,300,400,500,600,700,800,900,1000,1100,1200,1300]
treename = 'apic'
patchxs0 = [[0.5-0.45*sin(0.025*x*2*pi) for x in range(0,21)],[0.5-0.45*sin(0.025*x*2*pi) for x in range(20,41)]]
patchys0 = [[0.5+0.45*cos(0.025*x*2*pi) for x in range(0,21)],[0.5+0.45*cos(0.025*x*2*pi) for x in range(20,41)]]
patchxs = patchxs0[0]+[0.5,0,0,1,1,0.5]+patchxs0[1]
patchys = patchys0[0]+[0,0,1,1,0,0]+patchys0[1]


minmaxes = [inf,-inf]
coeffs_saved = [[],[]]
Ihcoeffs = [0.0,1.0]

for imodel in range(0,2):
 for imodul in range(0,2):
  moduldv = 10.0 - 20.0*imodul
  for idist1 in range(0,len(dist1s)):
    dist1 = dist1s[idist1]
    for idist2 in range(idist1+1,len(dist1s)):
      dist2 = dist1s[idist2]
      threshEcons_thisdist2 = []
      for iIhcoeff in range(0,2):
        Ihcoeff = Ihcoeffs[iIhcoeff]
        foundOne = 0
        threshEcons = []
        for myseed in range(1,41):
          if imodel == 0:
            if Ihcoeff == 0:
              filename = '../modulhcn_hay/ffthreshs/ffthreshs_basalthreshcoeff0.5_basaldt0.0_Ihmod'+str(moduldv)+'_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'
            else:
              filename = '../modulhcn_hay/ffthreshs/ffthreshs_basalthreshcoeff0.5_basaldt0.0_Ihcoeff1.0_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'
          else:
            if Ihcoeff == 0:
              filename = '../modulhcn_hay/ffthreshs/ffthreshs_basalgabacond0.0002_unijitter25.0_Ihmod'+str(moduldv)+'_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'
            else:
              filename = '../modulhcn_hay/ffthreshs/ffthreshs_basalgabacond0.0002_unijitter25.0_Ihcoeff1.0_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'

          if not exists(filename):
            #print(filename+' does not exist, trying dist without decimal point')
            filename = filename.replace(',500.0_',',500_')
          if exists(filename):
            unpicklefile = open(filename,'rb')
            unpickledlist = pickle.load(unpicklefile)
            unpicklefile.close()
            threshEcons.append(unpickledlist[0][2])
            foundOne = 1
          else:
            print(filename+' does not exist either')
        if not foundOne:
          print(filename+' does not exist')
        threshEcons_thisdist2.append(threshEcons[:])

        if len(threshEcons) == 0:
          print('Ihcoeff = '+str(Ihcoeff)+', continue')
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
        print('Ihcoeff = '+str(Ihcoeff)+', mycol = '+str(mycol))

        if iIhcoeff == 0:
          polygon = Polygon(array([[idist1,idist1+1,idist1],[idist2,idist2,idist2+1]]).T, True) #lower left triangle Ih blocked
        else:
          polygon = Polygon(array([[idist1,idist1+1,idist1+1],[idist2+1,idist2,idist2+1]]).T, True) #upper right triangle control
        p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
        p.set_facecolor(mycol); p.set_edgecolor('none')
        print('Ihcoeff = '+str(Ihcoeff)+', mycol = '+str(mycol))
        axarr[imodel+2*imodul].add_collection(p)
      if len(threshEcons_thisdist2) > 1 and len(threshEcons_thisdist2[0]) > 1 and len(threshEcons_thisdist2[1]) > 1:
        s1 = mean(threshEcons_thisdist2[0])
        s2 = mean(threshEcons_thisdist2[1])
        print('s1 = '+str(s1)+', s2 = '+str(s2))
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
        axarr[imodel+2*imodul].add_collection(p)
        if pval < 0.05/66:
          axarr[imodel+2*imodul].plot(idist2+0.5,idist1+0.5,'k*',lw=0.5,mew=0.5,ms=2.0)

  print('minmaxes = '+str(minmaxes))

axnew = f.add_axes([0.02,0.27,0.05,0.69])
axnew2 = f.add_axes([0.89,0.27,0.05,0.69])
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

axnew2.text(2.2,1.5,'Threshold larger when $I_h$ modulated',fontsize=4,color='#FF0000',ha='left',va='bottom',rotation=90,clip_on=False)
axnew2.text(2.2,-0.5,'Threshold smaller when $I_h$ blocked',fontsize=4,color='#0000FF',ha='left',va='top',rotation=90,clip_on=False)


for iax in range(0,4):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.03 - (iax > 1)*0.02, pos.y1 - 0.01, chr(ord('A')+iax), fontsize=9)


for iax in range(0,2):
  axarr[iax].set_xticks([0.5+x for x in range(0,len(dist1s))])
  axarr[iax].set_yticks([0.5+x for x in range(0,len(dist1s))])
  axarr[iax].set_xticklabels([str(d) for d in dist1s],fontsize=4,rotation=20)
  axarr[iax].set_yticklabels([str(d) for d in dist1s],fontsize=4,rotation=20)
  axarr[iax].set_ylabel('distance ($\mu$m)',fontsize=6)
  axarr[2*iax+1].set_xlabel('distance ($\mu$m)',fontsize=6)

  axarr[2+iax].set_xticks([0.5+x for x in range(0,len(dist1s))])
  axarr[2+iax].set_yticks([0.5+x for x in range(0,len(dist1s))])
  axarr[2+iax].set_xticklabels([])
  axarr[2+iax].set_yticklabels([])
axarr[0].set_xticklabels([])
axarr[3].set_xticklabels([str(d) for d in dist1s],fontsize=4,rotation=20)
  
f.savefig("figS5.eps")


