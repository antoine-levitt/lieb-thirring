# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import functions
from functions import classical_lt

def build_rho_2D(g, V, gamma, sigma, numeigs=4):
    '''Construit rho à partir des premières fonctions propres de -Delta+V'''
    bigN = g.N**g.d
    H = g.Delta + functions.tridiag(bigN,0,V.flat,0)
    w,vs = scipy.sparse.linalg.eigsh(H,numeigs,sigma=sigma,which='LM')
    if w[-1] < 0:
        return build_rho_2D(g,V,gamma, sigma, numeigs*2)
    elif w[0] > 0:
        raise Exception('No negative eigenvalue')
    print 'Eigs', w
    vecs = []
    for i in range(0,len(w)):
        if w[i] < 0:
            vecs.append(vs[:,i].reshape((g.N,g.N)))

    rho = np.zeros_like(V)
    for i in range(0,len(vecs)):
        rho = rho + np.abs(w[i])**(gamma-1)*vecs[i]**2
    return rho, w[w < 0], vecs

def build_V_2D(g, rho, gamma):
    '''Construction de V à partir de la densité rho'''
    V = -rho**(1/(gamma+g.d/2-1))
    V /= compute_lp_2D(g, V, gamma+g.d/2)
    return V

def compute_lp_2D(g,V,p):
    """Compute the lp norm of V"""
    V = np.abs(V)
    sumpow = g.delta**g.d * np.sum(V**p)
    return sumpow ** (1/p)


def scf_2D(g,Vinit,gamma,maxiter,tol, sigma):
    Vs = [Vinit]
    rhos = []
    Cs = []
    V = Vinit
    prevC = 0
    eigs = [sigma]
    psiss = []
    for i in range(0,maxiter):
        print 'Itération %d' %(i+1)
        rho,eigs,psis = build_rho_2D(g,V,gamma,sigma=eigs[0])
        V = build_V_2D(g,rho,gamma)
        psiss.append(psis)
        Vs.append(V)
        rhos.append(rho)
        res = scipy.linalg.norm((Vs[-1] - Vs[-2])/(Vs[-1]))

        print "Résidu %s" %res
        C = sum(abs(eigs)**gamma) / classical_lt(gamma,g.d)
        Cs.append(C)
        print "C", C
        print "deltaC", C - prevC
        if C - prevC < 0:
            print "ATTENTION, ENERGIE DECROISSANTE"
        elif np.abs(C - prevC) < tol:
            break
        prevC = sum(abs(eigs)**gamma) / classical_lt(gamma,g.d)
        print ''

    return V,rho,eigs,Vs,rhos,psiss, np.array(Cs)
