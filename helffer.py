from __future__ import division
from functions import classical_lt, surf_nsphere
from scipy.special import gamma
import math
from numpy import *

#sum_i+j+k=n 1
def coeff(i,d):
    prod = 1
    for j in range(1,d):
        prod *= (i+j)/j
    return prod
    if d == 1:
        return 1
    elif d == 2:
        return i+1
    elif d == 3:
        return ((i+1)*(i+2))/2
    elif d == 4:
        return ((i+1)*(i+2)*(i+3))/6
    else:
        prod = 1
        for j in range(1,d):
            prod *= (i+j)/j
        return prod
            
def compute_r(gam,d,h):
    bound = floor((1-d*h)/(2*h)) #last valid i
    i = arange(0,bound+1)
    eig = (1-d*h-(2*i)*h)
    if eig[-1] < 0:
        print "eig[-1]", eig[-1]
        eig[-1] = 0 #if d and h are in an exact ratio, roundoff errors might introduce negative numbers
    s = 0

    return math.fsum(coeff(i,d) * eig ** gam)
    # return sum(coeff(i,d) * eig ** gam)

def integral(gam,d):
    return surf_nsphere(d) / 2 * gamma(gam+d/2+1)*gamma(d/2)/gamma(gam+d+1)

def compute_C(gam,d,h):
    print d,gam,h
    return h**d*compute_r(gam,d,h)/integral(gam,d)/classical_lt(gam,d)

    
# h = .0001
# gammas = linspace(0,1.5,1000);
# d = 2
# Cs = [compute_C(gam, d, h) - 1 for gam in gammas]
# Cs = array(Cs)
# figure()
# plot(gammas,Cs)
# #plot(gammas,-sign(Cs)*log(abs(Cs/1000)))
# plot(gammas,zeros_like(Cs))
# show()


gam = 1.05
# hs = logspace(-6,-7,10);
hs = logspace(-6,-8,10);
hs = [10e-9]
# hs = linspace(1e-7,1e-8, 10)
d = 3
Cs = [compute_C(gam, d, h) - 1 for h in hs]
Cs = array(Cs)
figure()
plot(hs,Cs)
#plot(gammas,-sign(Cs)*log(abs(Cs/1000)))
plot(hs,zeros_like(Cs))
show()
