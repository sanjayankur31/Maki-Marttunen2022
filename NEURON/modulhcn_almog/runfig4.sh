echo "python3 drawmorph_multicolor.py"
python3 drawmorph_multicolor.py
for myseed in `seq 1 41`
do
  for Ihcoeff in 1.0 0.0
  do
    echo "python3 calcbasalffthreshs_givendists.py $Ihcoeff 0 356 $myseed"
    python3 calcbasalffthreshs_givendists.py $Ihcoeff 0 356 $myseed
    echo "python3 calcffthreshs_givendists.py $Ihcoeff 0 1300 $myseed"
    python3 calcffthreshs_givendists.py $Ihcoeff 0 1300 $myseed
    for dist1 in `seq 200 100 1200`
    do
      for dist2 in `seq $((dist1+100)) 100 1300`
      do
        echo "python3 calcffthreshs_givendists.py $Ihcoeff $dist1 $dist2 $myseed"
        python3 calcffthreshs_givendists.py $Ihcoeff $dist1 $dist2 $myseed
      done
    done
  done
done



    
