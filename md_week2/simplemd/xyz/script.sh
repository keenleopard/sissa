for dt in 0.005 0.01 0.015 0.02 0.025 
do 
#sed s/'$time'/$dt/g in > in_$dt
../src/simplemd.x < in_$dt
done
