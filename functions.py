# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import scipy
import scipy.special
import scipy.linalg
import scipy.sparse
import scipy.integrate
import scipy.sparse.linalg
import scipy.ndimage.measurements
import primme_inter
from matplotlib.pyplot import *

def get_schrodinger_eigen(g, V, sigma, numeigs, l, method):
    '''Retourne les numeigs premières valeurs propres de -Delta + V,
    avec moment angulaire l. Sigma est une approximation du bas du
    spectre.
    '''
    params = g.fem_params(l)
    
    # Matrices FEM
    Delta,Mass = g.fem_mats(l)
    
    # Formation du hamiltonien
    H = Delta
    if params['phi']:
        Vl = (l+(g.d-1)/2)*(l+(g.d-3)/2) * g.x**(-2)
        Vlop = g.fem_mult_op(Vl,l,Mass)
        # Vlop = (l+(g.d-1)/2)*(l+(g.d-3)/2) * g.fem_overxsquare_op()
        H = H + Vlop
    else:
        #shift V
        V0 = extrap(0,g.x,V)

        V = np.roll(V,1)
        V[0] = V[1]
        V[0] = V0

    Vop = g.fem_mult_op(V,l,Mass)
    H = H + Vop
    
    # Diagonalisation
    if method == 'eigsh':
        w,v = scipy.sparse.linalg.eigsh(H,numeigs,sigma=sigma,which='LM',M=Mass)
        needmore = w[-1] < 0
    if method == 'primme':
        w,v = primme_inter.solve_eigen(H,numeigs,Mass)
        needmore = w[-1] < 0
    elif method == 'banded':
        # print "WARNING : ignores mass matrix"
        Hb = tridiag_to_banded(g,H)
        w,v = scipy.linalg.eig_banded(Hb,overwrite_a_band=True,select='v',select_range=(-1000,0))
        needmore = 0
    elif method == 'lobpcg':
        X = np.random.rand(len(V),numeigs)
        w,v = scipy.sparse.linalg.lobpcg(H + tridiag(len(V),0,sigma,0),X,largest=False,verbosityLevel=1,maxiter=800)
        w = w - sigma
        print w
        needmore = w[-1] < 0

        
    # Changement de variable
    for i in range(0,len(w)):
        if params['phi']:
            phi = v[:,i]
            psi = phi_to_psi(g,phi)
        else:
            psi = np.roll(v[:,i],-1)
            psi[-1] = psi[-2]
            # psi = v[:,i]
        # print "Norm by mass matrix", np.dot(psi, Mass * psi)
        # print "Norm by compute_lp", compute_lp(g,psi,2) / np.sqrt(surf_nsphere(2))
        v[:,i] = psi / compute_lp(g,psi,2)

        v[:,i] = v[:,i] * np.sign(v[2,i]) # force positive sign

    return w,v,needmore

def get_all_negative_schrodinger_eigen(g,V,sigma,numeigs, l=0):
    "Retourne la première positive si il n'y a pas de négative"
    w,v,needmore = get_schrodinger_eigen(g,V,sigma,numeigs,l,'eigsh')
    # w,v,needmore = get_schrodinger_eigen(g,V,sigma,numeigs,l,'eigsh')
    if needmore:
        print numeigs, "valeurs propres pas suffisantes, je réessaie avec", numeigs*2
        return get_all_negative_schrodinger_eigen(g,V,w[0],numeigs*2,l)
    if len(w) == 0:
        w = np.ones(1)
        v = np.zeros((g.N,1))
        return w,v
    elif w[0] > 0:
        return w[[0]], v[:,[0]]
        
    # print "Vp pour l=",l, ": ", w[w < 0], ", première positive", w[w>=0][0]
    print "Vp pour l=",l, ": ", w
    return w[w < 0],v[:, w < 0]
    
def build_rho(g,V,gamma, sigmas=[-1], numeigs=4):
    '''Construit rho à partir des premières fonctions propres de -Delta+V'''

    # l = 0
    w,v = get_all_negative_schrodinger_eigen(g,V,sigmas[0],numeigs,0)
    newsigmas = [w[0]]
    shift = w[0]
    if (w > 0).all():
        print "ATTENTION: pas de vp négatives, cutoff", w[0]
        newsigmas = sigmas
    else:
        # l > 0
        if g.radial:
            l = 1
            while True:
                #choix du shift
                if len(sigmas) > l:
                    shift = sigmas[l]
                elif l > 1:
                    shift = wl[0]
                #else, just use shift from l-1

                wl,vl = get_all_negative_schrodinger_eigen(g,V,shift,numeigs,l)
                if (wl > 0).all():
                    break
                else:
                    w = np.hstack((w,np.repeat(wl,int(round(multipl_lb(g.d,l))))))
                    v = np.hstack((v,np.repeat(vl,int(round(multipl_lb(g.d,l))),1)))
                newsigmas = np.hstack((newsigmas,wl[0]))
                l += 1


    # construction de rho
    rho = np.zeros(V.shape)
    for i in range(0,w.size):
        rho = rho + np.abs(w[i])**(gamma-1)*v[:,i]**2
    # rho /= np.sum(np.abs(w)**gamma)**(1-1/gamma)
    
    return rho, w, v, newsigmas

def build_V(g, rho, gamma):
    '''Construction de V à partir de la densité rho'''
    V = -rho**(1/(gamma+g.d/2-1))
    V /= compute_lp(g, V, gamma+g.d/2)
    return V
    

def compute_lp(g,V,p):
    """Compute the lp norm of V"""
    Vs = np.abs(V)

    # if g.radial and (g.d <= 3):
    #     sumpow = 0
    #     h = g.delta
    #     V = Vs[:-1]
    #     dV = Vs[1:] - V
    #     r = g.x[:-1]
    #     special_ind = np.abs(dV) < 1e-10
    #     normal_ind = np.abs(dV) >= 1e-10
    #     integrands = np.zeros_like(V)
    #     if g.d == 1:
    #         integrands[special_ind] = V[special_ind]**p * h
    #     elif g.d == 2:
    #         integrands[special_ind] = V[special_ind]**p * h * (g.x[special_ind] + h/2)
    #     elif g.d == 3:
    #         integrands[special_ind] = V[special_ind]**p * h * (g.x[special_ind]**2 + g.x[special_ind]*h + h**2/3)
            
    #     V = V[normal_ind]
    #     dV = dV[normal_ind]
    #     r = r[normal_ind]
    #     # de maple
    #     if g.d == 1:
    #         integrands[normal_ind] = h * (-V ** (p + 1) + (V + dV) ** p * V + (V + dV) ** p * dV) / dV / (p + 1)
    #     elif g.d == 2:
    #         integrands[normal_ind] = -h * (V ** (p + 1) * p * dV * r - V ** (p + 2) * h + 2 * V ** (p + 1) * r * dV - (V + dV) ** p * V * p * dV * r - (V + dV) ** p * V * p * dV * h - 2 * (V + dV) ** p * V * r * dV + (V + dV) ** p * V ** 2 * h - (V + dV) ** p * dV ** 2 * p * r - (V + dV) ** p * dV ** 2 * p * h - 2 * (V + dV) ** p * dV ** 2 * r - (V + dV) ** p * dV ** 2 * h) / dV ** 2 / (p ** 2 + 3 * p + 2)
    #     elif g.d == 3:
    #         integrands[normal_ind] = h * (2 * (V + dV) ** p * V ** 3 * h ** 2 + 2 * (V + dV) ** p * dV ** 3 * h ** 2 + 6 * (V + dV) ** p * dV ** 3 * r ** 2 + (V + dV) ** p * V * p ** 2 * dV ** 2 * r ** 2 + (V + dV) ** p * V * h ** 2 * p ** 2 * dV ** 2 + (V + dV) ** p * V * h ** 2 * p * dV ** 2 - 2 * (V + dV) ** p * V ** 2 * h ** 2 * dV * p - 6 * (V + dV) ** p * V ** 2 * r * h * dV + 5 * (V + dV) ** p * V * p * dV ** 2 * r ** 2 + 2 * (V + dV) ** p * dV ** 3 * h * p ** 2 * r + 8 * (V + dV) ** p * dV ** 3 * h * p * r + 2 * V ** (p + 2) * h * p * dV * r + 3 * (V + dV) ** p * dV ** 3 * h ** 2 * p + 5 * (V + dV) ** p * dV ** 3 * p * r ** 2 + 6 * (V + dV) ** p * dV ** 3 * h * r + 6 * V ** (p + 2) * r * h * dV - 5 * V ** (p + 1) * p * dV ** 2 * r ** 2 - 6 * V ** (p + 1) * dV ** 2 * r ** 2 - 2 * V ** (p + 3) * h ** 2 - V ** (p + 1) * p ** 2 * dV ** 2 * r ** 2 + 2 * (V + dV) ** p * V * h * dV ** 2 * p ** 2 * r - 2 * (V + dV) ** p * V ** 2 * h * p * dV * r + 6 * (V + dV) ** p * V * h * p * dV ** 2 * r + (V + dV) ** p * dV ** 3 * h ** 2 * p ** 2 + 6 * (V + dV) ** p * V * dV ** 2 * r ** 2 + (V + dV) ** p * dV ** 3 * p ** 2 * r ** 2) / dV ** 3 / (p + 3) / (p + 2) / (p + 1)
            
    #     sumpow = sum(integrands) * surf_nsphere(g.d)
    # elif g.radial:
    #     sumpow = g.delta * np.sum(Vs**p * g.ds)
    # else:
    #     sumpow = g.delta**g.d * np.sum(Vs**p)

    sumpow = np.sum(g.delta* Vs**p * g.ds)
    
    return sumpow ** (1/p)

def surf_nsphere(n):
    '''Surface de la sphère en dimension n'''
    return 2*np.pi**(n/2) / scipy.special.gamma(n/2)

def classical_lt(gamma,d):
    '''Valeur de la constante de LT semi-classique'''
    return 2**(-d) * np.pi**(-d/2) * scipy.special.gamma(gamma+1) / scipy.special.gamma(gamma + 1 + d/2)

def one_d_one_bs(gamma):
    '''Constante optimale en 1D, état 1-bound state'''
    return 2*((gamma-.5)/(gamma+.5))**(gamma-.5)

def multipl_lb(d,l):
    "Multiplicité de la valeur propre l de l'opérateur de Laplace-Beltrami"
    return (2*l+d-2)*scipy.special.gamma(d+l-2)/scipy.special.gamma(d-1)/scipy.special.gamma(l+1)
    
def phi_to_psi(g,phi):
    return phi / np.sqrt(g.ds)

def psi_to_phi(g,psi):
    return psi * np.sqrt(g.ds)
    
def translate(g,V):
    V = V.reshape(g.x.shape)
    com = scipy.ndimage.measurements.center_of_mass(V)
    Vtranslated = V
    for i in range(int(g.d)):
        if int(com[i]-g.N/2) != 0:
            print "Shifted center of mass by", int(com[i]-g.N/2), "in direction", i
        Vtranslated = np.roll(Vtranslated,int(com[i]-g.N/2),i)

    return Vtranslated.flatten()
    

def scf(g,Vinit,gamma,maxiter,tol, sigma):
    Vs = [Vinit]
    rhos = []
    Cs = []
    V = Vinit
    prevC = 0
    # eigs = [-1000/(g.L**2)]
    sigmas = [sigma]
    eigs = [sigma]
    # eigs = [0]
    psiss = []
    for i in range(0,maxiter):
        print 'Itération %d' %(i+1)
        rho,eigs,psis,sigmas = build_rho(g,V,gamma,sigmas=sigmas)
        V = build_V(g,rho,gamma)
        psiss.append(psis)
        # if g.radial == False:
        #     V = translate(g,V)
        Vs.append(V)
        rhos.append(rho)
        res = scipy.linalg.norm((Vs[len(Vs)-1] - Vs[len(Vs)-2])/(Vs[len(Vs)-1]))

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

def plot_radial(g,V):
    plot(g.onedgrid[g.N/2:],V.reshape(g.x.shape)[g.N/2,g.N/2:])

def tridiag_to_banded(g,M):
    ab = np.zeros((2, M.get_shape()[0]))
    ab[0,:] = M[2,3]
    ab[1,:] = M.diagonal()

    return ab

def tridiag(N,a,b,c):
    # replaced by "diags" in scipy 0.12
    diags = scipy.zeros((3,N))
    diags[0,:] = a
    diags[1,:] = b
    diags[2,:] = c
    diags[2,1:] = diags[2,0:-1]

    return scipy.sparse.spdiags(diags, (-1,0,1),N,N,'csr')

def radial_stddev(g,V):
    return np.sqrt(np.sum(g.x**2 * g.ds * np.abs(V)) / np.sum(g.ds * np.abs(V)))

def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y = np.where(x < xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = np.where(x > xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
    return y
