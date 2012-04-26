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
gammas = np.linspace(1,1.03,10)
gammas = np.array([1.0894])
gammas = np.array([1.4])
#gammas = np.linspace(.84,1,5)
gammas = np.linspace(1.1,1.3,8)
# gammas = np.linspace(.5,1.5,100)
gammas = np.array([1.2])
maxiter = 20
tol = 1e-10
sigma = -1

d = 1
# L = 40
Vinit_size = 1
# N = 10000
q = 1
init = 'gauss'
init = 'prev'
# init = 'prev_scaled'
# init = 'harm'
#init = 'prev'

g = grid.Grid(N,L,d,q)
N = g.N

if init == 'prev':
    Vinit = functions.extrap(g.x, prev_g.x, V)
elif init == 'prev_scaled':
    scaling = 3
    Vinit = np.interp(g.x, prev_g.x*scaling, V,right=0)
elif init == 'gauss':
    Vinit = -np.exp(-np.abs((g.x)/Vinit_size)**2) * (1 + (g.x/Vinit_size)**2 + (g.x/Vinit_size)**4)
    Vinit = -np.exp(-np.abs((g.x)/Vinit_size)**2) * (1)
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

# Vinit = -np.exp(-np.abs((g.x - 3)/Vinit_size)**2) * (1)
# Vinit += -np.exp(-np.abs((g.x + 3)/Vinit_size)**2) * (1)


# i = 72
# Vinit = shelf['Vgammas'][i]
# gammas = [shelf['gammas'][i]]
# g = shelf["g"]
# N = g.N
# L = g.L


import Logger
sys.stdout = Logger.Logger('computation.log')
Vinit = Vinit / compute_lp(g,Vinit,gammas[0]+d/2)
V = Vinit
Vgammas = []
psiss = []
try:
    for gamma in gammas:
        if len(gammas) != 1:
            print
            print "********** Calcul avec gamma = ", gamma, "************"
            print
        V /= compute_lp(g,V,gamma+d/2)
        V,rho,eigs,Vs,rhos,psiss,Cs = scf(g,V,gamma,maxiter=maxiter,tol=tol,sigma=sigma)
        sigma = eigs[0]
        Rs,rs = extrapolation.richardson(Cs)
        if len(Rs) > 0:
            print "Calcul fini, Richardson dit", Rs[-1]
            Cgamma.append(Rs[-1])
        else:
            Cgamma.append(Cs[-1])
        Vgammas.append(V)
finally:
    sys.stdout.close()

Cgamma = np.array(Cgamma)

if(len(Cgamma) > 1):
    gammasextr = np.linspace(1,1.5,100)
    Cgammaextr = functions.extrap(gammasextr, gammas,Cgamma)

    figure()
    plot(gammasextr,Cgammaextr)
    plot(gammas,Cgamma,'x')
    title('N = {},L = {}'.format(N,L))
    show()

psis = psiss[-1]
psi = psis[:,0]
phi = functions.psi_to_phi(g,psi)

# figure()
# # plot(g.x,psis,'-x')
# [plot(g.x,V) for V in Vs]
# show()

prev_g = g




# w, v, needmore = get_schrodinger_eigen(g,V,sigma,2,0,"eigsh")
