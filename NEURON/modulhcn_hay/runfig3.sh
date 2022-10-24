
#Assume the following lines were run from runfig2.sh - the results from these runs are needed for Fig. 3
#for Ihcoeff in 1.0 0.0
#do
#  for dist in `seq 50.0 50.0 1000.0`
#  do
#    echo "python3 strongdendstim_findthresh_absbound.py $Ihcoeff $dist"
#    python3 strongdendstim_findthresh_absbound.py $Ihcoeff $dist
#    echo "python3 strongdendstimcond_findthresh_absbound.py $Ihcoeff $dist"
#    python3 strongdendstimcond_findthresh_absbound.py $Ihcoeff $dist
#  done
#done

for Ihcoeff in 1.0 0.0
do
  for dist in `seq 50.0 50.0 1000.0`
  do
    echo "python3 strongdendstim_findthresh_absbound_nohotLVA.py $Ihcoeff $dist"
    python3 strongdendstim_findthresh_absbound_nohotLVA.py $Ihcoeff $dist
    echo "python3 strongdendstimcond_findthresh_absbound_nohotLVA.py $Ihcoeff $dist"
    python3 strongdendstimcond_findthresh_absbound_nohotLVA.py $Ihcoeff $dist
    echo "python3 strongdendstim_findthresh_additionalblocked_absbound.py $Ihcoeff $dist gCa_LVAstbar_Ca_LVAst"
    python3 strongdendstim_findthresh_additionalblocked_absbound.py $Ihcoeff $dist gCa_LVAstbar_Ca_LVAst
  done
  for dist in 800
  do
    echo "python3 strongdendstimcond_recCaLVA_Ih.py $Ihcoeff $dist"
    python3 strongdendstimcond_recCaLVA_Ih.py $Ihcoeff $dist
  done
done

