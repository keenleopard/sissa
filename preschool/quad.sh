#! /bin/bash
# Estimate integral of x^19 over [0,1]

awk '
function f(z) { return z**19}
function quad(a,b){ 
    y=(b-a)* f((a+b)/2.); 
    return y
}

function abs(x){
    if (x>0) return x;
    else return -x;
}

function ad(a, b, tau){
    q = quad(a,(a+b)/2) + quad((a+b)/2,b)
    if (abs(q-quad(a,b))>tau) {
        q = ad(a, (a+b)/2, tau/2) + ad( (a+b)/2, b,tau/2)
    }
    return q
}
    
BEGIN{
    print "the estimate is", ad(0,1,1e-7);
}
'
#    for (q=1./5.; e>tau; i++){
#        m = (a+b)/2.;
#        qnew = ad(a,m) + ad(m,b);
#        e = abs(q-qnew);
#        q = qnew
