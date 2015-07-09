#! /bin/bash
awk '
function f(x){
    y = x - cos(x);
    return y;
}
function abs(x){
    if (x>0) return x;
    else return -x;
}

BEGIN{
    xmin = 0.;
    xmax = 1.;
    xold = 1.;
    for (xnew = 0.5; (abs(xnew-xold)/xold) > 0.00001; xnew = (xmin+xmax)/2.){
        if (f(xmin)*f(xnew)< 0) xmax = xnew;
        else xmin = xnew;
        xold = xnew;
        print xnew
    }
}'
