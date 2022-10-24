
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

for i in range(0,12):
  axarr[i].set_xlim([-5,75])
  axarr[i].set_ylim([-80,120])
for i in [1,2,5,6]:
  axarr[i].set_ylim([0,1.01])
#axarr[4].set_xlim([-40,328])
for i in [3,7,11]:
  axarr[i].set_visible(False)
  

iamp = 2
#Threshold, normal
#picklelist = [Is,times_all,m_calvas_all,h_calvas_all,Vsomas_all,ecas_all,icalvas_all]
unpicklefile = open('../modulhcn_hay/singlecompartment_cond_recs0.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[0].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'k-',lw=0.5)
axarr[1].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[2][iamp],'k-',lw=0.5)
axarr[2].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[3][iamp],'k-',lw=0.5)
axarr[10].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'k-',lw=0.5)

#Threshold, Ih blocked
#picklelist = [Is,times_all,m_calvas_all,h_calvas_all,Vsomas_all,ecas_all,icalvas_all]
unpicklefile = open('../modulhcn_hay/singlecompartment_cond_recs1.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[4].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'b-',lw=0.5)
axarr[5].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[2][iamp],'b-',lw=0.5)
axarr[6].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[3][iamp],'b-',lw=0.5)
axarr[10].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'b-',lw=0.5)

#Threshold, h_CaLVA replaced by that recorded in absence of Ih
#picklelist = [Is,times_all,m_calvas_all,h_calvas_all,Vsomas_all]
unpicklefile = open('../modulhcn_hay/singlecompartment_cond_artificialCaLVAcondreadeca_recs0_mfrombl1_hfrombl0_dtstim1.0_fixeddt.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[8].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'m-',lw=0.5)
axarr[10].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'m-',lw=0.5)

#Threshold, h_CaLVA replaced by that recorded in absence of Ih
#picklelist = [Is,times_all,m_calvas_all,h_calvas_all,Vsomas_all]
unpicklefile = open('../modulhcn_hay/singlecompartment_cond_artificialCaLVAcond_recs0_mfrombl0_hfrombl1_dtstim1.0_fixeddt.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
axarr[9].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'y-',lw=0.5)
axarr[10].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'y-',lw=0.5)

#(this is not needed)
##Threshold, m_CaLVA and h_CaLVA replaced by those recorded in absence of Ih
##picklelist = [Is,times_all,m_calvas_all,h_calvas_all,Vsomas_all]
#unpicklefile = open('../modulhcn_hay/singlecompartment_cond_artificialCaLVAcond_recs0_mfrombl1_hfrombl1_dtstim1.0_fixeddt.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
#axarr[10].plot([x-800 for x in unpickledlist[1][iamp]],unpickledlist[4][iamp],'g-',lw=0.5)

iCHR = 0
for iax in [0,1,2,4,5,6,8,9,10]:
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.03, pos.y1 - 0.01, chr(ord('A')+iCHR), fontsize=9)
  iCHR = iCHR + 1

f.savefig("figS3.eps")
