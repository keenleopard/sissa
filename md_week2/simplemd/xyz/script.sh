for T in 0.3 2.0 
do 
sed s/'$temp'/$T/g in > in_$T
../src/simplemd.x < in_$T
done


# a compact to write for is 
# for ((temp=0.0;temp<100;temp++))
