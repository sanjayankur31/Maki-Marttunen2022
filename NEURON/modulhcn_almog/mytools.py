import math
import copy
import pylab

def spike_times(time,vrec,V_min_peak=-20,V_max_valley=0):
  valley_reached = 1
  sptime = []  
  for j in range(1,len(time)-1):
    if valley_reached and vrec[j] >= V_min_peak and vrec[j] > vrec[j-1] and vrec[j] >= vrec[j+1]:
      valley_reached = 0
      sptime.append(time[j])
    elif valley_reached==False and vrec[j] <= V_max_valley:
      valley_reached = 1
  return sptime  

#membpotderivs(time,vrec): Given the membrane potentials (vrec) at time points time[0],time[1],...,time[N-1],
#return the derivatives at time points time[1],time[2],...,time[N-3]
def membpotderivs(time,vrec):
  N = len(time)
  tdiff = [x-y for x,y in zip(time[1:N-1],time[0:N-2])]
  vdiff = [x-y for x,y in zip(vrec[1:N-1],vrec[0:N-2])]
  mderiv = [x/y for x,y in zip(vdiff,tdiff)]
  return [0.5*(x+y) for x,y in zip(mderiv[1:N-2],mderiv[0:N-3])]

#limitcyclescaledv(v1,dv1,v2,dv2): Give the coefficient for memb. pot. derivative that one has to use in order to make
#the difference on the derivative axis as significant as the difference on the memb. pot. axis
def limitcyclescaledv(v1,dv1,v2,dv2):
  maxv = max(max(v1),max(v2))
  minv = min(min(v1),min(v2))
  maxdv = max(max(dv1),max(dv2))
  mindv = min(min(dv1),min(dv2))
  return 1.0*(maxv-minv)/(maxdv-mindv)

def limitcyclediff(v1,dv1,v2,dv2,dvcoeff=0.1):
  N1 = len(v1)
  N2 = len(v2)
  dv1 = [dvcoeff*x for x in dv1]
  dv2 = [dvcoeff*x for x in dv2]
  Dmin = N1*[0]
  for i in range(0,N1):
    Dmin[i] = math.sqrt(min([(x-v1[i])**2+(y-dv1[i])**2 for x,y in zip(v2,dv2)]))
  vdiff = [x-y for x,y in zip(v1[1:N1],v1[0:N1-1])] 
  dvdiff = [x-y for x,y in zip(dv1[1:N1],dv1[0:N1-1])] 
  h = [math.sqrt(x**2+y**2) for x,y in zip(vdiff,dvdiff)]
  #use the trapezoid rule for integration:
  Dminmean = [(x+y)/2.0 for x,y in zip(Dmin[1:N1],Dmin[0:N1-1])]
  #print "hsum="+str(sum(h))
  return sum([x*y for x,y in zip(Dminmean,h)])
  

def interpolate(tref,vref,tint): #Assumes that the trefs come sorted!
  vint = len(tint)*[0.0]
  addedOne = False
  #print tref
  #print tint
  #if tref[len(tref)-1] == tint[len(tint)-1]:
  #  tref.append(tref[len(tref)-1]+0.0001)
  #  vref.append(vref[len(tref)-1])
  #  addedOne = True
  if tref[0] > tint[0] or tref[len(tref)-1] < tint[len(tint)-1]:
    print("Extrapolation needed!")
    return len(tint)*[-1]
  indvrecnow = 0  
  for j in range(0,len(tint)):
    while tref[indvrecnow+1] <= tint[j]:
      indvrecnow = indvrecnow + 1
      if indvrecnow == len(tref)-1: # It must be the last index if this happens
        vint[j:len(tint)] = [vref[indvrecnow]]*(len(tint)-j)
        return vint
    vint[j] = vref[indvrecnow] + 1.0*(tint[j]-tref[indvrecnow])/(tref[indvrecnow+1]-tref[indvrecnow])*(vref[indvrecnow+1]-vref[indvrecnow])
  return vint

def interpolate_extrapolate_constant(tref,vref,tint): #Assumes that the trefs come sorted!
  vint = len(tint)*[0.0]
  addedOne = False
  #print tref
  #print tint
  #if tref[len(tref)-1] == tint[len(tint)-1]:
  #  tref.append(tref[len(tref)-1]+0.0001)
  #  vref.append(vref[len(tref)-1])
  #  addedOne = True

  indvrecnow = 0  
  for j in range(0,len(tint)):
    if tint[j] < tref[0]:
      vint[j] = vref[0]
      continue
    if indvrecnow >= len(tref) - 1:
      vint[j] = vref[-1]
    while tref[indvrecnow+1] <= tint[j]:
      indvrecnow = indvrecnow + 1
      if indvrecnow == len(tref)-1: # It must be the last index if this happens
        vint[j:len(tint)] = [vref[indvrecnow]]*(len(tint)-j)
        return vint
    vint[j] = vref[indvrecnow] + 1.0*(tint[j]-tref[indvrecnow])/(tref[indvrecnow+1]-tref[indvrecnow])*(vref[indvrecnow+1]-vref[indvrecnow])
  return vint

#kronecker product of list A and list B
def kron(A,B):
  C = []
  if type(B[0]) is int or type(B[0]) is float:
    for i in range(0,len(A)):
      for j in range(0,len(B)):
        print("asdf")
        print((B[j]))
        C.append(A[i]*B[j])
  elif type(B[0][0]) is int or type(B[0][0]) is float:
    for i in range(0,len(A)):
      for j in range(0,len(B)):
        C.append([x*A[i] for x in B[j]])
  return C
  
def cumprod(A):
  B = len(A)*[0]; B[0]=A[0]
  for j in range(1,len(A)):
    B[j] = B[j-1]*A[j]
  return B

def oscillatorypoissontimeseries(N, minfreq, maxfreq, oscfreq, phase, T, isprinted=0):
  a = 0.5*(minfreq + maxfreq)
  b = 0.5*(maxfreq - minfreq)
  meanT = 1.0/a
  omega = 2*pylab.pi*oscfreq
  sptime = pylab.zeros([int(0.001*2*T*N/meanT)+100,2],dtype='d')
  #sptime = []
  Nplaced = 0
  for iN in range(0,N):
    t = 0
    while t < T:
      e = pylab.log(1-pylab.rand()) 
      Ttrial = meanT #Ttrial in sec, t and T in msec
      # Newton's method:
      for iter in range(0,50):
        Ttrialold = Ttrial
        Ttrial = Ttrial - (e + a*Ttrial - b/omega*(pylab.cos(omega*(0.001*t+Ttrial)+phase)-pylab.cos(omega*0.001*t+phase))) / (a + b*pylab.sin(omega*(0.001*t+Ttrial)+phase))
      myerr = abs(e + a*Ttrial - b/omega*(pylab.cos(omega*(0.001*t+Ttrial)+phase)-pylab.cos(omega*0.001*t+phase)))
      if Ttrial < 0 or myerr > 0.0001:
        Ttrial = meanT
        # If Newton's method didn't give good result, try fixed point method: Find T such that -aT+b/omega(pylab.cos(omega*t+phase)-pylab.cos(phase)) = pylab.log(1-x):
        for iter in range(0,30000):
          Ttrial = (b/omega*(pylab.cos(omega*(0.001*t+Ttrial)+phase)-pylab.cos(omega*0.001*t+phase)) - e)/a
        Ttrialold = Ttrial
        Ttrial = (b/omega*(pylab.cos(omega*(0.001*t+Ttrial)+phase)-pylab.cos(omega*0.001*t+phase)) - e)/a
      if isprinted > 1:
        print(('iN='+str(iN)+', Err(Ttrial) = '+str(abs(Ttrial-Ttrialold))+', Ttrial = '+str(Ttrial)+', t = '+str(t)+', newPhase='+str(omega*(t/1000+Ttrial)+phase)+', myerr='+str(myerr)))
      t = t+Ttrial*1000
      sptime[Nplaced,:] = [iN,t]
      #sptime.append([iN,t])
      Nplaced = Nplaced + 1
    if isprinted > 0:
      print(('spikes for iN='+str(iN)+' complete'))
  sptime = sptime[:Nplaced,:]
  if isprinted > 0:
    print((str(Nplaced)+' spikes in total, '+str(int(0.001*2*T*N/meanT)+100)+' reserved'))
  return sptime

def printlistlen(A):
  #TODO: recursive method might work out but needs some thought...
  #toCheck = A
  #lens = [len(x) for x in A]
  #levels = [0 for x in A]
  #while type(toCheck) is list and len(toCheck) > 0:
  #  while type(toCheck[0]) is list and len(toCheck[0]) > 0:
  #    toCheck.append(toCheck[0][0])
  #    lens.append[len
  #    toCheck[0].pop(0)
  nan = -1
  if type(A) is list:
   B = copy.deepcopy(A)
   listFound0 = False
   for i0 in range(0,len(B)):
    if type(B[i0]) is list:
     listFound0 = True
     listFound1 = False
     for i1 in range(0,len(B[i0])):
      if type(B[i0][i1]) is list:
       listFound1 = True
       listFound2 = False
       for i2 in range(0,len(B[i0][i1])):
        if type(B[i0][i1][i2]) is list:
         listFound2 = True
         listFound3 = False
         for i3 in range(0,len(B[i0][i1][i2])):
          if type(B[i0][i1][i2][i3]) is list:
           listFound3 = True
           listFound4 = False
           for i4 in range(0,len(B[i0][i1][i2][i3])):
            if type(B[i0][i1][i2][i3][i4]) is list:
             listFound4 = True
             listFound5 = False
             for i5 in range(0,len(B[i0][i1][i2][i3][i4])):
              if type(B[i0][i1][i2][i3][i4][i5]) is list:
               listFound5 = True
               listFound6 = False               
               for i6 in range(0,len(B[i0][i1][i2][i3][i4][i5])):
                if type(B[i0][i1][i2][i3][i4][i5][i6]) is list:
                  listFound6 = True
                  B[i0][i1][i2][i3][i4][i5][i6] = len(B[i0][i1][i2][i3][i4][i5][i6])
                else:
                  B[i0][i1][i2][i3][i4][i5][i6] = nan
               if not listFound6:
                 B[i0][i1][i2][i3][i4][i5] = len(B[i0][i1][i2][i3][i4][i5])
              else:
                B[i0][i1][i2][i3][i4][i5] = nan
             if not listFound5:
               B[i0][i1][i2][i3][i4] = len(B[i0][i1][i2][i3][i4])
            else:
              B[i0][i1][i2][i3][i4] = nan
           if not listFound4:
             B[i0][i1][i2][i3] = len(B[i0][i1][i2][i3])
          else:
            B[i0][i1][i2][i3] = nan
         if not listFound3:
           B[i0][i1][i2] = len(B[i0][i1][i2])
        else:
          B[i0][i1][i2] = nan
       if not listFound2:
         B[i0][i1] = len(B[i0][i1])
      else:
        B[i0][i1] = nan
     if not listFound1:
       B[i0] = len(B[i0])
    else:
      B[i0] = nan
   if not listFound0:
     B = len(B)
  else:
   B = nan
  print(B)

def drawarrow(ax,x,y,acoeff=1,prc=0.9,lw=1,lc='#000000'):
  d = [x[1]-x[0], y[1]-y[0]];
  k = pylab.sqrt(d[0]**2 + d[1]**2)
  d = d/k

  dperp = [acoeff*z for z in [-d[1], d[0]]];
  lens = k-k*(1-prc);
  perplen = 0.5*k*(1-prc);

  px = [x[0],x[1],x[0]+lens*d[0]+perplen*dperp[0],x[0]+lens*d[0]-perplen*dperp[0]];
  py = [y[0],y[1],y[0]+lens*d[1]+perplen*dperp[1],y[0]+lens*d[1]-perplen*dperp[1]];

  px = [px[0],px[1],px[2],px[1],px[3]]
  py = [py[0],py[1],py[2],py[1],py[3]]
  #px = reshape([px(:,[1 2 3 2 4]) nan(size(x,1),1)]',size(x,1)*6,1);
  #py = reshape([py(:,[1 2 3 2 4]) nan(size(x,1),1)]',size(x,1)*6,1);
  ax.plot(px,py,'k-',linewidth=lw,color=lc)

def timeseriesmean(times,x):
  return 1.0*sum([(t2-t1)*(x1+x2)/2 for t1,t2,x1,x2 in zip(times[0:-1],times[1:],x[0:-1],x[1:])])/(times[-1]-times[0])

def timeseriessecondmoment(times,x):
  return 1.0*sum([(t2-t1)*(x1**2+x2**2)/2 for t1,t2,x1,x2 in zip(times[0:-1],times[1:],x[0:-1],x[1:])])

def timeseriesstd(times,x,xmean=pylab.nan):
  if pylab.isnan(xmean):
    xmean = timeseriesmean(times,x)
  return pylab.sqrt(1.0*sum([(t2-t1)*((x1-xmean)**2+(x2-xmean)**2)/2 for t1,t2,x1,x2 in zip(times[0:-1],times[1:],x[0:-1],x[1:])])/(times[-1]-times[0]))

def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y-0.5*yoffset,y+0.5*yoffset],'k-',color='#FFFFFF',zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

def firingratecurve(spikes,T=[],dt=1.0,gauss_std=5.0):
  if type(spikes) is list:
    spikes = pylab.array(spikes)
  if type(T) is list and len(T) == 0:
    T = max(spikes[0,:])
  FRs = pylab.zeros([int(T/dt),1])
  FRts = [dt*(i+0.5) for i in range(0,len(FRs))]
  for iFRt in range(0,len(FRs)):
    FRs[iFRt] = sum(pylab.exp(-1/2*(FRts[iFRt]-spikes[0,:])**2/gauss_std**2))
  return [FRs,FRts]

def find(condition):
    "Return the indices where ravel(condition) is true"
    res, = np.nonzero(np.ravel(condition))
    return res

def mylegend(fig,pos,styles,texts,nx=1,dx=2,yplus=0.5,yplustext=0.35,colors=[],dashes=[],linewidths=[],myfontsize=8):
  ny = int(pylab.ceil(1.0*len(styles)/nx))
  print("ny = "+str(ny))
  axnew = fig.add_axes(pos)
  handles = []
  for i in range(0,len(styles)):
    handles.append(axnew.plot([dx*int(i/ny)+0.15,dx*int(i/ny)+0.35],[yplus+ny-1-(i%ny),yplus+ny-1-(i%ny)],styles[i]))
    axnew.text(dx*int(i/ny)+0.5,yplustext+ny-1-(i%ny),texts[i],fontsize=myfontsize)
  axnew.set_xlim([0,dx*nx])
  axnew.set_ylim([0,len(styles)])
  for i in range(0,len(dashes)):
    if len(dashes[i]) > 0:
      handles[i][0].set_dashes(dashes[i])
  for i in range(0,len(colors)):
    if len(colors[i]) > 0:
      handles[i][0].set_color(colors[i])
  for i in range(0,len(linewidths)):
    if type(linewidths[i]) is not list and linewidths[i] > 0:
      handles[i][0].set_linewidth(linewidths[i])
  axnew.get_xaxis().set_visible(False)
  axnew.get_yaxis().set_visible(False)
  axnew.set_xlim([0,dx*nx])
  axnew.set_ylim([0,ny])
  return axnew

import math

def colorsredtolila(N,brightness=1.0,whiteness=0.0):
  if N==1:
    return ['#0000FF']
  C = colorsredtolilaint(N,brightness,whiteness)
  hexlist = []
  for j in range(0,N):
    myhex = '#'
    for k in range(0,3):
      if C[j][k] < 16:
        myhex = myhex + '0' + hex(C[j][k])[2]
      else:
        myhex = myhex + hex(C[j][k])[2:4]
    hexlist.append(myhex)
  return hexlist

def colorsredtolilaint(N,brightness=1.0,whiteness=0.0):

  linchange = [1.0/(N-1)*x for x in range(0,N)]
  linchange_x = [x/0.8 for x in [0.0, 0.36, 0.45, 0.65, 0.7, 0.8]]
  linchange_y = [0.0, 0.29, 0.46, 0.65, 0.75, 0.83]
  for i in range(0,len(linchange)):
    ind = N
    for j in range(0,len(linchange_x)):
      if linchange[i] >= linchange_x[j]:
        ind = j
    if ind < len(linchange_x)-1:
      linchange[i] = linchange_y[ind]+(linchange_y[ind+1]-linchange_y[ind])*(linchange[i]-linchange_x[ind])/(linchange_x[ind+1]-linchange_x[ind])

  C = hsv2rgblist([int(255*x) for x in linchange],[1.0-whiteness]*N,[brightness]*N);
  return C


def hsv2rgblist(h, s, v):
  hs = []
  ss = []
  vs = []
  for i in range(0,len(h)):
    h1, s1, v1 = hsv2rgb(h[i],s[i],v[i])
    hs.append(h1)
    ss.append(s1)
    vs.append(v1)
  return list(zip(hs,ss,vs))

def hsv2rgb(h, s, v):
  h = float(h)
  s = float(s)
  v = float(v)
  h60 = h / 60.0
  h60f = math.floor(h60)
  hi = int(h60f) % 6
  f = h60 - h60f
  p = v * (1 - s)
  q = v * (1 - f * s)
  t = v * (1 - (1 - f) * s)
  r, g, b = 0, 0, 0
  if hi == 0: r, g, b = v, t, p
  elif hi == 1: r, g, b = q, v, p
  elif hi == 2: r, g, b = p, v, t
  elif hi == 3: r, g, b = p, q, v
  elif hi == 4: r, g, b = t, p, v
  elif hi == 5: r, g, b = v, p, q
  r, g, b = int(r * 255), int(g * 255), int(b * 255)
  return r, g, b
    
def rgb2hsv(r, g, b):
  r, g, b = r/255.0, g/255.0, b/255.0
  mx = max(r, g, b)
  mn = min(r, g, b)
  df = mx-mn
  if mx == mn:
    h = 0
  elif mx == r:
    h = (60 * ((g-b)/df) + 360) % 360
  elif mx == g:
    h = (60 * ((b-r)/df) + 120) % 360
  elif mx == b:
    h = (60 * ((r-g)/df) + 240) % 360
  if mx == 0:
    s = 0
  else:
    s = df/mx
  v = mx
  return h, s, v
