# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np, matplotlib, scipy as sp, scipy.linalg as la
# from matplotlib.pyplot import *

import functions
reload(functions)
from functions import *
import grid
reload(grid)
import extrapolation
reload(extrapolation)

Cgamma = []
gammas = np.linspace(1.2,1.25,3)
gammas = np.array([1.167])
#gammas = np.linspace(.84,1,5)
maxiter = 40
tol = 1e-10
sigma = -.1

d = 2
L = 100
Vinit_size = 5
h = 0.1
# N = int(L / h)
# N = 800
# q optimal pour singularité r^a : 5/(1+2*a)
# N = 8208
N = 800
q = 2
init = 'gauss'
# init = 'prev_scaled'
# init = 'harm'

g = grid.Grid(N,L,d,q)
N = g.N

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
elif init == 'linear':
    Vinit = (g.x-L/2)
    Vinit[Vinit > 0] = 0
else:
    raise Error('init invalide')


Vinit = Vinit / compute_lp(g,Vinit,gammas[0]+d/2)
V = Vinit
Vgammas = []
psiss = []
for gamma in gammas:
    if len(gammas) != 1:
        print
        print "********** Calcul avec gamma = ", gamma, "************"
        print
    V,rho,eigs,Vs,rhos,psiss,Cs = scf(g,V,gamma,maxiter=maxiter,tol=tol,sigma=sigma)
    sigma = eigs[0]
    Rs,rs = extrapolation.richardson(Cs)
    Cgamma.append(Cs[-1])
    Vgammas.append(V)

    print "Calcul fini, Richardson dit", Rs[-1]

Cgamma = np.array(Cgamma)

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
    plot(gammas,Cgamma,'x')
    show()

# figure()
# for V in Vs:
#     plot(g.x,V)
# show()
# title ("V")

psis = psiss[-1]
psi = psis[:,0]
phi = functions.psi_to_phi(g,psi)

figure()
plot(g.x,psis,'-x')
show()

prev_g = g

# figure()
# plot(psiss[-1]);plot(psiss[-8])
# show()

