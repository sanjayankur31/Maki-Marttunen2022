
echo "python3 drawmorph_unicolor.py"
python3 drawmorph_unicolor.py
for Ihcoeff in 1.0 0.0 2.0
do
    echo "python3 calcifcurves2.py $Ihcoeff"
    python3 calcifcurves2.py $Ihcoeff
    echo "python3 calcifcurves2_wait_findthresh.py $Ihcoeff"
    python3 calcifcurves2_wait_findthresh.py $Ihcoeff
done

for Ihmod in -5.0 5.0
do
    echo "python3 calcifcurves2_modih.py $Ihmod"
    python3 calcifcurves2_modih.py $Ihmod
    echo "python3 calcifcurves2_wait_findthresh_modih.py $Ihmod"
    python3 calcifcurves2_wait_findthresh_modih.py $Ihmod
done

    
