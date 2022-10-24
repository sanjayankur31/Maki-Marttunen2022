
for Ihcoeff in 1.0 0.0
do

  for dist in 500.0 800.0
  do
    echo "python3 strongdendstim.py $Ihcoeff $dist"
    python3 strongdendstim.py $Ihcoeff $dist
  done
  for dist in `seq 50.0 50.0 1000.0`
  do
    echo "python3 strongdendstim_findthresh_absbound.py $Ihcoeff $dist"
    python3 strongdendstim_findthresh_absbound.py $Ihcoeff $dist
    echo "python3 strongdendstimcond_findthresh_absbound.py $Ihcoeff $dist"
    python3 strongdendstimcond_findthresh_absbound.py $Ihcoeff $dist
  done
done

    
