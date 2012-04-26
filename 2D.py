# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np, matplotlib, scipy as sp, scipy.linalg as la
# from matplotlib.pyplot import *
import functions2D
reload(functions2D)
from functions2D import *
import grid2D
reload(grid2D)
import extrapolation
reload(extrapolation)

Cgamma = []
gammas = np.linspace(1.2,1.25,3)
gammas = np.array([1.2])
maxiter = 200
tol = 1e-10
sigma = -.05

d = 2
L = 60
Vinit_size = 18
N = 75
g = grid2D.Grid2D(N,L,d)

Vinit = -np.exp(-((g.x)/Vinit_size)**2 - (g.y/Vinit_size)**2)
# Vinit += -np.exp(-((g.x+8)/Vinit_size)**2 - (g.y/Vinit_size)**2)
# Vinit *= (1 +(g.x/Vinit_size)**2 + (g.y/Vinit_size)**2)
# Vinit = -np.exp(-np.abs(g.x/Vinit_size) - np.abs(g.y/Vinit_size))
Vinit *= (1+np.cos(4*g.theta))

Vinit = Vinit / compute_lp_2D(g,Vinit,gammas[0]+d/2)
V = Vinit * 1
Vgammas = []
psiss = []
for gamma in gammas:
    if len(gammas) != 1:
        print
        print "********** Calcul avec gamma = ", gamma, "************"
        print
    V,rho,eigs,Vs,rhos,psiss,Cs = scf_2D(g,V,gamma,maxiter=maxiter,tol=tol,sigma=sigma)
    sigma = eigs[0]
    if len(Cs) > 3:
        Rs,rs = extrapolation.richardson(Cs)
        print "Calcul fini, Richardson dit", Rs[-1]
    Cgamma.append(Cs[-1])
    Vgammas.append(V)


# Cgamma = np.array(Cgamma)

# if(len(Cgamma) > 1):
#     gammasextr = np.linspace(1,1.5,100)
#     Cgammaextr = functions.extrap(gammasextr, gammas,Cgamma)

#     figure()
#     plot(gammasextr,Cgammaextr)
#     plot(gammas,Cgamma,'x')
#     show()


figure()
contourf(g.x,g.y,V)
gray()
show()

diffC = Cs[1:] - Cs[:-1]
