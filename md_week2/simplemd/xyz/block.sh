#! /bin/bash
echo "#l  var_Epot  sig_T"

#for leng in 1 2 4 8 16 32 64 128
for leng in $(seq 1 1 20000)
do
    if [ $((50000/leng*leng)) -ne 50000 ] 
    then 
        continue
    else 
        awk -v l=$leng 'BEGIN{
            tot=50000;
            nb=int(tot/l);
        }{
        block[int((NR)/l)] += $4;
        a += $4;
        }
        END{
            mean = a/(tot);
            for (i=2;i<nb;i++){
                block[i] /= l;
                var += (block[i]-mean)**2;
            }
            print nb, var/(nb**2)
        }' < energies_2.0.dat
    fi
done
