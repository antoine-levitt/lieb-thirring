from __future__ import division
import functions
import numpy as np
import scipy
import scipy.sparse
class Grid:
    def __init__(self,N=100,L=10.0,d=1,radial=False):
        self.N = N
        self.L = L
        self.d = d
        self.radial = radial
        if radial:
            x = np.linspace(L,0,N,endpoint=False)[::-1]
        else:
            x = np.linspace(-L,L,N)
            x = np.linspace(0,L,N)
        self.onedgrid = x
        self.delta = x[1] - x[0]
        # fd ou fe. d=2, l=0 est hard-code
        self.method = 'fe'

        self.x = x
        if radial:
            self.ds = functions.surf_nsphere(d) * x**(d-1)

    def fem_mats(self,l):
        # returns Delta, Mass
        N = self.N
        
        if self.d == 2 and l == 0:
            h = self.x[1] - self.x[0]
            # the grid has been shifted
            x = self.x - h
            
            Mdiag = h * 2/3 * x
            Moffdiag = h * 1/6 * (x + h/2)
            Mdiag[0] = h**2 / 12

            Ddiag = 2 / h * x #OK
            Doffdiag = -1 / h * (x + h/2) #OK
            Ddiag[0] = 1/2 #OK
            
            Delta = functions.tridiag(N,Doffdiag,Ddiag,Doffdiag)
            Mass = functions.tridiag(N,Moffdiag,Mdiag,Moffdiag)

        else:
            Delta = functions.tridiag(N,-1,2,-1) / self.delta
            if self.method == 'fd':
                Delta = functions.tridiag(N,-1,2,-1) / self.delta**2
                # Mass = functions.tridiag(N,0,1,0)
                Mass = None
            else:
                Delta = functions.tridiag(N,-1,2,-1) / self.delta
                Mass = functions.tridiag(N,1,4,1) * self.delta / 6
            # Mass = functions.tridiag(N,0,1,0) * self.delta
            
        return Delta,Mass

    def fem_mult_op(self,V,l,Mass):
        # Vsqrtop = scipy.sparse.spdiags(np.sqrt(np.abs(V)),0,V.shape[0],V.shape[0],'csr')
        # return np.sign(V[0]) * np.dot(np.dot(Vsqrtop,Mass),Vsqrtop)

        if(self.d == 2 and l == 0):
            # the grid has been shifted
            h = self.delta
            x = self.x - h
            
            diag = V / 2 * x
            diag[1:] += 1/12*V[:-1] * (x[1:] - 2/5*h)
            diag[:-1] += 1/12*V[1:] * (x[:-1] + 2/5*h)
            diag[0] = h*(1/20*V[0] + 1/30 * V[1])
            offdiag = 1/12*V * (x + 2/5*h)
            offdiag[:-1] += 1/12*V[1:] * (x[:-1] + 3/5*h)

            return h*functions.tridiag(self.N,offdiag,diag,offdiag)
        else:
            if self.method == 'fd':
                return functions.tridiag(self.N,0,V,0)
            else:
                diag = V / 2
                diag[1:] += 1/12*V[:-1]
                diag[:-1] += 1/12*V[1:]
                offdiag = 1/12*V
                offdiag[:-1] += 1/12*V[1:]
            
                return self.delta*functions.tridiag(self.N,offdiag,diag,offdiag)
