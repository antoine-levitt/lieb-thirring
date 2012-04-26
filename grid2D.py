from __future__ import division
import functions
import numpy as np
import scipy
import scipy.sparse

class Grid2D:
    def __init__(self,N=100,L=10.0,d=1):
        self.N = N
        self.L = L
        self.d = d
        self.onedgrid = np.linspace(-L,L,N)
        self.x,self.y = np.meshgrid(self.onedgrid, self.onedgrid)
        self.r = np.sqrt(self.x**2 + self.y**2)
        self.theta = np.arctan2(self.x,self.y)
        self.delta = self.onedgrid[1] - self.onedgrid[0]

        Delta1D = functions.tridiag(N,-1,2,-1)/self.delta**2
        eye = functions.tridiag(N,0,1,0)
        self.Delta = scipy.sparse.kron(Delta1D,eye,'csr') + scipy.sparse.kron(eye,Delta1D,'csr')
