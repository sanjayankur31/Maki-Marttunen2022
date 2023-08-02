
#Assume the following lines were run from runfig2.sh - the results from these runs are needed for Fig. 5
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
  for dist in `seq 50 50 1000`
  do
    echo "python3 strongdendstim_findthresh_absbound_withbasalcond.py $Ihcoeff $dist 0.1 50.0 0.0"
    python3 strongdendstim_findthresh_absbound_withbasalcond.py $Ihcoeff $dist 0.1 50.0 0.0
    echo "python3 strongdendstim_findthresh_absbound_withbasalcond.py $Ihcoeff $dist 5.0 50.0 -80.0"
    python3 strongdendstim_findthresh_absbound_withbasalcond.py $Ihcoeff $dist 5.0 50.0 -80.0
  done
done

for myseed in `seq 1 41`
do
  for Ihcoeff in 1.0 0.0
  do
    for dist1 in `seq 200 100 1200`
    do
      for dist2 in `seq $((dist1+100)) 100 1300`
      do
        echo "python3 calcffthreshs_givendists_withbasal.py $Ihcoeff $dist1 $dist2 0.8 0.0 $myseed"
	python3 calcffthreshs_givendists_withbasal.py $Ihcoeff $dist1 $dist2 0.8 0.0 $myseed
        echo "python3 calcffthreshs_givendists_withunijitterbasalgaba.py $Ihcoeff $dist1 $dist2 0.0002 25.0 $myseed"
	python3 calcffthreshs_givendists_withunijitterbasalgaba.py $Ihcoeff $dist1 $dist2 0.0002 25.0 $myseed
      done
    done
  done
done

