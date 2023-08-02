
echo "python3 drawmorph_unicolor.py"
python3 drawmorph_unicolor.py
for Ihcoeff in 1.0 0.0 2.0
do
    echo "python3 calcifcurves.py $Ihcoeff"
    python3 calcifcurves.py $Ihcoeff
    echo "python3 calcifcurves_wait_findthresh.py $Ihcoeff"
    python3 calcifcurves_wait_findthresh.py $Ihcoeff
done

for Ihmod in -10.0 10.0
do
    echo "python3 calcifcurves_modih.py $Ihmod"
    python3 calcifcurves_modih.py $Ihmod
    echo "python3 calcifcurves_wait_findthresh_modih.py $Ihmod"
    python3 calcifcurves_wait_findthresh_modih.py $Ihmod
done

    
