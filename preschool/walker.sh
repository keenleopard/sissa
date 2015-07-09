#! /bin/bash
# Generate 1D random walker by using awk
awk 'BEGIN{
    srand();
    for (i=1; i<=100; i++){
        r = rand()
        if (r>0.5) s++
        else s--
    }
    print "the final position is", s
}'
