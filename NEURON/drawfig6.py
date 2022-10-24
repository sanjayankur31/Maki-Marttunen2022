from pylab import *
import scipy.io
import pickle
from os.path import exists
from matplotlib.collections import PatchCollection
import scipy.stats
import time
import mytools

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

def discont(ax,x,y,width=18,height=0.1,col='#000000'):
  ax.plot([x-width/2,x+width/2],[y-height/2,y+height/2],'k-',lw=2.3,clip_on=False,color=col)
  ax.plot([x-width,x+width],[y-height,y+height],'w-',lw=1.0,clip_on=False,zorder=100)

f,axarr = subplots(4,1)

for iax in range(0,4):
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()

for ix in range(0,2):
  axarr[ix].set_position([0.06+0.50*ix,0.64,0.43,0.26])
  axarr[2+ix].set_position([0.13+0.43*ix,0.13,0.34,0.42])
  

dists = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
styles = ['bx-','kx-']


dist1s = [200,300,400,500,600,700,800,900,1000,1100,1200,1300]
treename = 'apic'
patchxs0 = [[0.5-0.45*sin(0.025*x*2*pi) for x in range(0,21)],[0.5-0.45*sin(0.025*x*2*pi) for x in range(20,41)]]
patchys0 = [[0.5+0.45*cos(0.025*x*2*pi) for x in range(0,21)],[0.5+0.45*cos(0.025*x*2*pi) for x in range(20,41)]]
patchxs = patchxs0[0]+[0.5,0,0,1,1,0.5]+patchxs0[1]
patchys = patchys0[0]+[0,0,1,1,0,0]+patchys0[1]

cols = ['#000000','#0000FF','#FF00FF']
iIhord = [1,0,2]
titles = [['DA at prox. apic. dend.','Control','DA at prox., ACh at distal apic. dend.'],['ACh at prox. apic. dend.','Control','ACh at prox., DA at distal apic. dend.']]

minmaxes = [inf,-inf]
coeffs_saved = [[],[]]
diagvals_all = []
for imodel in range(0,2):
  moduldv = 10.0 - 20*imodel
  moduldv_opposite = -moduldv
  diagvals = []
  for idist1 in range(0,len(dist1s)):
    dist1 = dist1s[idist1]
    for idist2 in range(idist1+1,len(dist1s)):
      dist2 = dist1s[idist2]
      threshEcons_thisdist2 = []
      for iIhcoeff in range(0,3):
        foundOne = 0
        threshEcons = []
        for myseed in range(1,41):
          if iIhcoeff == 0:
            filename = 'modulhcn_hay/ffthreshs/ffthreshs_Ihmod2ways'+str(moduldv)+',0.0,500.0_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'
          elif iIhcoeff == 1:
            filename = 'modulhcn_hay/ffthreshs/ffthreshs_Ihcoeff1.0_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'
          else:
            filename = 'modulhcn_hay/ffthreshs/ffthreshs_Ihmod2ways'+str(moduldv)+','+str(moduldv_opposite)+',500.0_'+treename+str(dist1)+'-'+str(dist2)+'_seed'+str(myseed)+'.sav'

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
          polygon = Polygon(array([[idist1,idist1+1,idist1],[idist2,idist2,idist2+1]]).T, True) #lower left triangle Ih modulated proxim
        elif iIhcoeff == 1:
          polygon = Polygon(array([[idist1,idist1+1,idist1+1],[idist2+1,idist2,idist2+1]]).T, True) #upper right triangle control
        else:
          polygon = Polygon(array([[idist1+x for x in patchxs],[idist2+y for y in patchys]]).T, True) #circle around Ih modulated both
        p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
        p.set_facecolor(mycol); p.set_edgecolor('none')
        axarr[2+imodel].add_collection(p)

        if iIhcoeff != 0:
          if iIhcoeff == 1 and len(threshEcons_thisdist2) > 1 and len(threshEcons_thisdist2[0]) > 1 and len(threshEcons_thisdist2[1]) > 1:
            s1 = mean(threshEcons_thisdist2[0]) #Ih modulated proximally
            s2 = mean(threshEcons_thisdist2[1]) #control
            pval = scipy.stats.ranksums(threshEcons_thisdist2[0], threshEcons_thisdist2[1])[1]
          elif iIhcoeff == 2 and len(threshEcons_thisdist2) > 2 and len(threshEcons_thisdist2[1]) > 1 and len(threshEcons_thisdist2[2]) > 1:
            s2old = mean(threshEcons_thisdist2[1]) #control
            s1 = mean(threshEcons_thisdist2[2]) #Ih modulated both proximally and distally
            s2 = mean(threshEcons_thisdist2[0]) #Ih modulated proximally
            pval = scipy.stats.ranksums(threshEcons_thisdist2[2], threshEcons_thisdist2[0])[1]
            pval2 = scipy.stats.ranksums(threshEcons_thisdist2[2], threshEcons_thisdist2[1])[1]
          else:
            print('Something may be wrong..')
            qwerty
            continue
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
          if iIhcoeff == 2:
            polygon = Polygon(array([[idist2,idist2+1,idist2],[idist1,idist1,idist1+1]]).T, True) #lower left triangle Ih modulated proxim
          elif iIhcoeff == 1:
            polygon = Polygon(array([[idist2,idist2+1,idist2+1],[idist1+1,idist1,idist1+1]]).T, True) #upper right triangle control
          #polygon = Polygon(array([[idist2+1,idist2+1,idist2,idist2],[idist1,idist1+1,idist1+1,idist1]]).T, True)
          p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
          p.set_facecolor(mycol); p.set_edgecolor('none')
          axarr[2+imodel].add_collection(p)
          if pval < 0.05/(3*66):
            if iIhcoeff == 1:
              #axarr[2+imodel].plot(idist2+0.7,idist1+0.7,'k*',lw=0.5,mew=0.5,ms=2.0)
              1 #pass for now
            else:
              if s2old < s2 < s1:
                if pval2 < 0.05/(3*66):
                  #axarr[2+imodel].plot(idist2+0.3,idist1+0.3,'k^',lw=0.5,mew=0.5,ms=2.0)
                  axarr[2+imodel].plot(idist2+0.5,idist1+0.5,'k^',lw=0.5,mew=0.5,ms=2.0)
                else:
                  axarr[2+imodel].plot(idist2+0.3,idist1+0.3,'k^',lw=0.5,mew=0.5,ms=2.0,fillstyle='none')
              elif s2old > s2 > s1:
                if pval2 < 0.05/(3*66):
                  #axarr[2+imodel].plot(idist2+0.3,idist1+0.3,'kv',lw=0.5,mew=0.5,ms=2.0)
                  axarr[2+imodel].plot(idist2+0.5,idist1+0.5,'kv',lw=0.5,mew=0.5,ms=2.0)
                else:
                  axarr[2+imodel].plot(idist2+0.3,idist1+0.3,'kv',lw=0.5,mew=0.5,ms=2.0,fillstyle='none')
              else:
                if pval2 < 0.05/(3*66):
                  #axarr[2+imodel].plot(idist2+0.3,idist1+0.3,'kx',lw=0.5,mew=0.5,ms=2.0)
                  1 #pass for now
                #else:
                #  axarr[2+imodel].plot(idist2+0.3,idist1+0.3,'k4',lw=0.5,mew=0.5,ms=2.0)
      if idist2 == idist1+1:
        diagvals.append(threshEcons_thisdist2[:])
    if len(diagvals) <= idist1:
      continue
    for iiIhcoeff in range(0,len(iIhord)):
      iIhcoeff = iIhord[iiIhcoeff]
      my_y = 1e3*mean(diagvals[idist1][iIhcoeff])
      axarr[imodel].bar(4*idist1+iiIhcoeff,my_y if my_y < 0.035 else 0.036+(my_y-0.036)*0.004,facecolor=cols[iiIhcoeff])
      if my_y > 0.04:
        print('imodel='+str(imodel)+',iIhcoeff='+str(iIhcoeff)+',idist1='+str(idist1)+',dist='+str(dist1s[idist1])+',y='+str(my_y))
      axarr[imodel].plot([4*idist1+iiIhcoeff,4*idist1+iiIhcoeff],[my_y-1e3*std(diagvals[idist1][iIhcoeff]),my_y+1e3*std(diagvals[idist1][iIhcoeff])],'k-',lw=0.3)
  print('minmaxes = '+str(minmaxes))
  diagvals_all.append(diagvals[:])

for iax in range(0,2):
  discont(axarr[iax],-1,0.035,width=0.5,height=0.0006,col='#000000')
  axarr[iax].set_xlim([-1,4*len(diagvals)])
  axarr[iax].set_ylim([0,0.04])
  axarr[iax].set_yticks([0,0.01,0.02,0.03])
  axarr[iax].set_xticks([4*idist1+1 for idist1 in range(0,len(diagvals))])
  axarr[iax].set_xticklabels([str(int(dist1s[idist1]))+'-'+str(int(dist1s[idist1+1])) for idist1 in range(0,len(diagvals))],rotation=20,ha='right',va='top')
axarr[0].set_ylabel('$g$ (nS)',fontsize=6)
for iax in range(0,2):
  diagvals = diagvals_all[iax]
  for ibar in [9,10]:
    for iiIhcoeff in range(0,3):
      iIhcoeff = iIhord[iiIhcoeff]
      my_y = 1e3*mean(diagvals[ibar][iIhcoeff])
      discont(axarr[iax],ibar*4+iiIhcoeff,0.035,width=0.5,height=0.0006,col=cols[iiIhcoeff])
      if ibar == 9:
        axarr[iax].plot([ibar*4+iiIhcoeff,ibar*4+2*iiIhcoeff-(10-ibar)*2],[0.0405,0.044],'k-',lw=0.3,color=cols[iiIhcoeff],clip_on=False,label=titles[iax][iIhcoeff])
      else:
        axarr[iax].plot([ibar*4+iiIhcoeff,ibar*4+2*iiIhcoeff-(10-ibar)*2],[0.0405,0.044],'k-',lw=0.3,color=cols[iiIhcoeff],clip_on=False)        
      if my_y  > 1.0:
        axarr[iax].text(ibar*4+2*iiIhcoeff-(10-ibar)*2,0.0445,'$\geq$1.0',color=cols[iiIhcoeff],rotation=90,fontsize=4,ha='center',va='bottom')
        #axarr[iax].text(ibar*4+2*iiIhcoeff-(10-ibar)*2,0.0445,'{:.3f}'.format(my_y),color=cols[iiIhcoeff],rotation=90,fontsize=4,ha='center',va='bottom')
        print('my_y = '+str(my_y))
      else:
        axarr[iax].text(ibar*4+2*iiIhcoeff-(10-ibar)*2,0.0445,'{:.3f}'.format(my_y),color=cols[iiIhcoeff],rotation=90,fontsize=4,ha='center',va='bottom')
for ibar in [6,7,8]:
  diagvals = diagvals_all[1]
  discont(axarr[1],ibar*4+2,0.035,width=0.5,height=0.0006,col='#FF00FF')
  iiIhcoeff = 2
  iIhcoeff = 2
  my_y = 1e3*mean(diagvals[ibar][iIhcoeff])
  discont(axarr[1],ibar*4+iiIhcoeff,0.035,width=0.5,height=0.0006,col=cols[iiIhcoeff])
  axarr[1].plot([ibar*4+iiIhcoeff,ibar*4-2],[0.0405,0.044],'k-',lw=0.3,color=cols[iiIhcoeff],clip_on=False)
  if my_y  > 1.0:
    axarr[1].text(ibar*4-2,0.0445,'$\geq$1.0',color=cols[iiIhcoeff],rotation=90,fontsize=4,ha='center',va='bottom')
  else:
    axarr[1].text(ibar*4-2,0.0445,'{:.3f}'.format(my_y),color=cols[iiIhcoeff],rotation=90,fontsize=4,ha='center',va='bottom')
for iax in range(2,4):
  axarr[iax].set_xticks([0.5+x for x in range(0,len(dist1s))])
  axarr[iax].set_xticklabels([str(d) for d in dist1s],fontsize=4,rotation=20,ha='right',va='top')
  axarr[iax].set_yticks([0.5+x for x in range(0,len(dist1s))])
  axarr[iax].set_yticklabels([str(d) for d in dist1s],fontsize=4,rotation=20,ha='right',va='top')

axarr[0].set_title('DA modulation of proximal dendrite\n\n\n',fontsize=6)
axarr[1].set_title('ACh modulation of proximal dendrite\n\n\n',fontsize=6)
axarr[2].set_title('DA modulation of proximal dendrite',fontsize=6)
axarr[3].set_title('ACh modulation of proximal dendrite',fontsize=6)
#leg0=axarr[0].legend(fontsize=4,edgecolor='none', facecolor=(1, 1, 1, 0.1),loc='upper left',borderpad=-0.25)
#leg1=axarr[1].legend(fontsize=4,edgecolor='none', facecolor=(1, 1, 1, 0.1),loc='upper left',borderpad=-0.25)
for iax in range(0,2):
  pos = axarr[iax].get_position()
  leg = mytools.mylegend(f,[pos.x0+0.02,pos.y0+0.22,0.22,0.05],['b-','b-','b-'],[titles[iax][i] for i in iIhord],1,3,0.5,0.35,colors=cols,myfontsize=4)
  leg.get_xaxis().set_visible(False)
  leg.get_yaxis().set_visible(False)
  leg.spines['top'].set_visible(False)
  leg.spines['right'].set_visible(False)
  leg.spines['bottom'].set_visible(False)
  leg.spines['left'].set_visible(False)

cols = ['#000000','#0000FF','#FF00FF']

axnew = f.add_axes([0.03,0.01,0.033,0.57])
axnew2 = f.add_axes([0.92,0.01,0.033,0.57])
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

axnew2.text(2.2,1.5,'Threshold larger during (extra) modulation',fontsize=4,color='#FF0000',ha='left',va='bottom',rotation=90,clip_on=False)
axnew2.text(2.2,-0.5,'Threshold smaller during (extra) modulation',fontsize=4,color='#0000FF',ha='left',va='top',rotation=90,clip_on=False)

axnew3 = f.add_axes([0.1,0.002,0.8,0.09])
axnew3.get_xaxis().set_visible(False)
axnew3.get_yaxis().set_visible(False)
for q in ['top','bottom','left','right']:
  axnew3.spines[q].set_visible(False)
axnew3.set_ylim([0,5.6])
axnew3.set_xlim([0,60])

for ipatch in range(0,3):
  ystart = 0.4+1.6*ipatch
  if ipatch == 1:
    polygon = Polygon(array([[0,0+1,0],[ystart,ystart,ystart+1]]).T, True) #lower left triangle Ih modulated proxim
  elif ipatch == 2:
    polygon = Polygon(array([[0,0+1,0+1],[ystart+1,ystart,ystart+1]]).T, True) #upper right triangle control
  else:
    polygon = Polygon(array([[0+x for x in patchxs],[ystart+y for y in patchys]]).T, True) #circle around Ih modulated both
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  p.set_facecolor('#FFFFFF'); p.set_edgecolor('#000000')
  axnew3.add_collection(p)
  
axnew3.text(1.2,0.4+1.6*2+0.5,'$I_h$ non-modulated',fontsize=5,ha='left',va='center')
axnew3.text(1.2,0.4+1.6*1+0.5,'$I_h$ modulated at prox. dend.',fontsize=5,ha='left',va='center')
axnew3.text(1.2,0.4+1.6*0+0.5,'$I_h$ modulated at prox. dend., opposite modulation at distal dend.',fontsize=5,ha='left',va='center')
axarr[2].set_zorder(2)
axarr[3].set_zorder(2)
for ipatch in range(0,2):
  ystart = 0.4+1.6*ipatch
  #if ipatch == 1:
  #  polygon = Polygon(array([[0,0+1,0],[ystart,ystart,ystart+1]]).T, True) #lower left triangle Ih modulated proxim
  #elif ipatch == 2:
  #  polygon = Polygon(array([[0,0+1,0+1],[ystart+1,ystart,ystart+1]]).T, True) #upper right triangle control
  #p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  #p.set_facecolor('#FFFFFF'); p.set_edgecolor('#000000')
  #axnew3.add_collection(p)
  if ipatch == 0:
    axnew3.plot([30+x for x in [0,0+1,0,0]],[ystart,ystart,ystart+1,ystart],'b-')
    axnew3.plot([30+x for x in [0,0+1,0,0]],[ystart,ystart,ystart+1,ystart],'r--',dashes=(1,1))
  else:
    axnew3.plot([30+x for x in [0,0+1,0+1,0]],[ystart+1,ystart,ystart+1,ystart+1],'b-')
    axnew3.plot([30+x for x in [0,0+1,0+1,0]],[ystart+1,ystart,ystart+1,ystart+1],'r--',dashes=(1,1))

axnew3.text(30+1.2,0.4+1.6*1+0.5,'$I_h$ modulated at prox. dend. vs. $I_h$ non-modulated',fontsize=5,ha='left',va='center')
axnew3.text(30+1.2,0.4+1.6*0+0.5,'Combined $I_h$ modulation vs. only proximal modulation',fontsize=5,ha='left',va='center')

for iax in range(0,len(axarr)):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.03, pos.y1 - 0.01, chr(ord('A')+iax), fontsize=9)

f.savefig("fig6.eps")


