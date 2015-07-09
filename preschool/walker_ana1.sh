#! /bin/bash
# find the average and variance of the final positions of walkers 
awk '{D+= ($1*$2); D2 += ($1**2 *$2); M += $2} END{Da=D/M; D2a = D2/M - Da**2; print Da, D2a}'



