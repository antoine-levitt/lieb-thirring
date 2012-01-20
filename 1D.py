# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np, matplotlib, scipy as sp, scipy.linalg as la
# from matplotlib.pyplot import *

import functions
reload(functions)
from functions import *
import grid
reload(grid)

Lgamma = []
gammas = np.linspace(1.2,1.25,3)
gammas = np.array([1.2])
#gammas = np.linspace(.84,1,5)
maxiter = 1
tol = 1e-4
sigma = -.005

d = 2
L = 50
Vinit_size = 5
h = 0.1
# N = int(L / h)
# N = 8000
# N = 10
N = 5000
radial = True
init = 'gauss'
# init = 'prev_scaled'
# init = 'harm'

g = grid.Grid(N,L,d,radial)

if init == 'prev':
    Vinit = np.interp(g.x, prev_g.x, V,right=0)
elif init == 'prev_scaled':
    scaling = 3
    Vinit = np.interp(g.x, prev_g.x*scaling, V,right=0)
elif init == 'gauss':
    Vinit = -np.exp(-np.abs((g.x)/Vinit_size)**2)
elif init == 'harm':
    Vinit = (g.x/Vinit_size)**2 - 1
    Vinit[Vinit > 0] = 0
elif init == 'constant':
    Vinit = -np.ones_like(g.x)
else:
    raise Error('init invalide')

V = Vinit / compute_lp(g,Vinit,gammas[0]+d/2)
Vgammas = []
psiss = []
for gamma in gammas:
    if len(gammas) != 1:
        print
        print "********** Calcul avec gamma = ", gamma, "************"
        print
    V,rho,eigs,Vs,rhos,psiss = scf(g,V,gamma,maxiter=maxiter,tol=tol,sigma=sigma)
    sigma = eigs[0]
    Lgamma.append(sum(abs(eigs)**gamma))
    Vgammas.append(V)


Cgamma = np.array(Lgamma) / classical_lt(gammas,d)

# plot(gammas,Cgamma)

# close()
# for rho in rhos:
#     plot(g.x,rho)
# show()
# title ("Rho")

if(len(Cgamma) > 1):
    gammasextr = np.linspace(1,1.5,100)
    Cgammaextr = functions.extrap(gammasextr, gammas,Cgamma)

    figure()
    plot(gammasextr,Cgammaextr)
    show()

# figure()
# for V in Vs:
#     plot(g.x,V)
# show()
# title ("V")

psis = psiss[-1]

# figure()
# plot(g.x,psis)
# show()

print radial_stddev(g,V)

prev_g = g

# figure()
# plot(psiss[-1]);plot(psiss[-8])
# show()

