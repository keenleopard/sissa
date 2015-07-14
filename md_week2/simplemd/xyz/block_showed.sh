#! /bin/bash

awk '{
a += $4
a2 += $4 **2
data[NR] = $4
}
END{
    var = a2/NR - (a/NR)**2
    avg = a/NR

    for(nb=2;nb<NR;nb++){
#        print nb
        if(int(NR/nb)*nb != NR) continue

        lb=int(NR/nb)
        var_lb = 0
        for (i=0;i<nb;i++){
            avg_block=0
            for (k=i*lb;k<(i+1)*lb;k++){
                avg_block += data[k+1]
            }
            avg_block/=lb
            var_lb += (avg_block-avg)**2
        }
        var_lb /= nb
        print lb, nb, sqrt(var_lb/nb)
    }
}'
