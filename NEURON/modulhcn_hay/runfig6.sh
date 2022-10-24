#Assume these have been run already (drawfig6.py uses these data)
#for myseed in `seq 1 41`
#do
#  for dist1 in `seq 200 100 1200`
#  do
#    for dist2 in `seq $((dist1+100)) 100 1300`
#    do
#      echo "python3 calcffthreshs_givendists.py 1.0 $dist1 $dist2 $myseed"
#      python3 calcffthreshs_givendists.py 1.0 $dist1 $dist2 $myseed
#    done
#  done
#done

for myseed in `seq 1 41`
do
  for Ihcoeff in 1.0 0.0
  do
    for dv in 10.0,0.0 -10.0,0.0 10.0,-10.0 -10.0,10.0
    do
	  
      for dist1 in `seq 200 100 1200`
      do
        for dist2 in `seq $((dist1+100)) 100 1300`
        do
          echo "python3 calcffthreshs_givendists_modih2ways.py ${dv},500.0 $dist1 $dist2 $myseed"
          python3 calcffthreshs_givendists_modih2ways.py ${dv},500.0 $dist1 $dist2 $myseed
	done
      done
    done
  done
done




    
