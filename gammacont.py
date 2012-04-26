# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np, matplotlib, scipy as sp, scipy.linalg as la
# from matplotlib.pyplot import *
import sys, time

import functions
reload(functions)
from functions import *
import grid
reload(grid)
import extrapolation
reload(extrapolation)
import shelve

gammas = np.linspace(1.1,1.13,10)
gammas = np.linspace(.88,.90,10)
gammas = np.linspace(.9,.91,5)
gammas = np.linspace(1.1,1.2,5)
gammas = np.linspace(1.1,.9,20)
gammas = np.linspace(1.1,1.5,20)
gammas = np.linspace(1.1,1.17,30)
gammas = np.linspace(.5, 1.5,100)
# gammas = np.linspace(.8,1,10)
# gammas = np.array([.8])
#gammas = np.linspace(.9,1.1,40)
maxiter = 100
tol = 1e-6
sigma = -1e-2

d = 2
L = 400
Vinit_size = 39
N = 2000
q = 1
init = 'gauss'

g = grid.Grid(N,L,d,q)

if init == 'prev':
    Vinit = np.interp(g.x, prev_g.x, V,right=0)
elif init == 'gauss':
    Vinit = -np.exp(-np.abs(g.x/Vinit_size)**2)
else:
    raise Error('init invalide')

# Vinit = shelve.open("3D/140.shelf")["Vgammas"][-1]

Vinit = Vinit / compute_lp(g,Vinit,gammas[0]+d/2)

V = Vinit
Vgammas = []
Cgamma = []
eigss = []


import Logger
sys.stdout = Logger.Logger('computation.log')
start_time = time.clock()
print "Début du calcul"
print "Grille : d = " + str(d) + "L = " + str(L) + ", N = " + str(N) + ", q = " + str(q)
print "Init : " + str(init) + ", Vinit_size = " + str(Vinit_size)
print "Scf : tol = " + str(tol) + ", sigma = " + str(sigma) + ", maxiter = " + str(maxiter)
print "Gammas de " + str(gammas[0]) + " à " + str(gammas[-1]) + " avec " + str(len(gammas)) + " points."
try:
    for gamma in gammas:
        elapsed = time.clock() - start_time
        estimate = (gammas[-1] - gamma) / (gamma-gammas[0]) * elapsed
        if gamma == gammas[0]:
            estimate = 0
        sec = int(estimate)
        min_, sec = divmod(sec, 60)
        hour, min_ = divmod(min_, 60)
        ETA = '%d:%02d:%02d' % (hour, min_, sec)

        print
        print "********** Calcul avec gamma = ", gamma, ", ", (gamma-gammas[0])/(gammas[-1] - gammas[0])*100, "%, ETA ", ETA,"************"
        print

        V /= compute_lp(g,V,gamma+d/2)
        V,rho,eigs,Vs,rhos,psiss,Cs = scf(g,V,gamma,maxiter=maxiter,tol=tol,sigma=sigma)
        sigma = 2*eigs[0]
        Rs,rs = extrapolation.richardson(Cs)
        if len(Rs) > 0:
            print "Calcul fini, Richardson dit", Rs[-1]
            # Cgamma.append(Rs[-1])
            Cgamma.append(Cs[-1])
            Vgammas.append(V)
        else:
            Cgamma.append(Cs[-1])
            Vgammas.append(V)
        eigss.append(eigs)
finally:
    sys.stdout.close()

bs = np.array([len(eigs) for eigs in eigss])
if any(bs != bs[0]):
    print ""
    print ""
    print "ALERTE, ALERTE, BS ONT CHANGE"
execfile("save_shelf.py")
import os
os.system("mv computation.log " + str(d) + "D/" + str(bs) + ".log")

Cgamma = np.array(Cgamma)

if(len(Cgamma) > 1):
    gammasextr = np.linspace(.7,1.5,100)
    Cgammaextr = functions.extrap(gammasextr, gammas,Cgamma)

    figure()
    plot(gammasextr,Cgammaextr)
    plot(gammas,Cgamma,'x')
    show()

psis = psiss[-1]

prev_g = g


gammasextr = np.linspace(.8,.9,100)
Cgammaextr = functions.extrap(gammasextr, gammas,Cgamma)
figure()
plot(gammasextr,Cgammaextr)
plot(gammas,Cgamma,'x')
show()
