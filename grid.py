from __future__ import division
import functions
import numpy as np
import scipy
import scipy.sparse
class Grid:
    def __init__(self,N=100,L=10.0,d=1,q=1):
        self.N = N
        self.L = L
        self.d = d
        self.radial = True
        x = np.linspace(L,0,N,endpoint=False)[::-1]
        

        method = 'uniform'
        method = 'geometric'
        if method == 'uniform':
            x = np.linspace(L,0,N,endpoint=False)[::-1]
        elif method == 'geometric':
            x = np.linspace(1,0,N,endpoint=False)[::-1]
            x = x**q
            x = x*L
        else:
            n = int(np.floor(np.log2(N)))
            kappa = 1/4
            x = np.array([kappa*L,L])
            for i in range(n):
                y = (x[1:-1] + x[2:])/2
                x = np.hstack((x,y,kappa*x[0]))
                x.sort()
            self.N = len(x)
            N = self.N
        
        #deltai = xi - xi-1
        self.deltam = np.zeros_like(x)
        self.deltam[0] = x[0]
        self.deltam[1:] = x[1:] - x[:-1]

        self.deltap = np.zeros_like(self.deltam)
        self.deltap[:-1]= self.deltam[1:]
        self.deltap[-1] = self.deltam[-1]
        self.delta = (self.deltap+self.deltam)/2


        # fd ou fe. d=2, l=0 est hard-code
        self.method = 'fe'

        self.x = x
        self.ds = functions.surf_nsphere(d) * x**(d-1)

    def fem_params(self,l):
        # return {'phi': True, 'weight': False, 'fd':False}
        if l == 0 and self.d == 2:
            return {'phi': False, 'weight': True, 'fd':False}
        else:
            return {'phi': True, 'weight': False, 'fd':False}
    
    def fem_mats(self,l):
        # returns Delta, Mass
        N = self.N
        params = self.fem_params(l)
        
        if not params['phi']:
            if self.d != 2:
                raise Exception('Not implemented')

            # the grid has been shifted
            x = np.zeros_like(self.x)
            x[1:] = self.x[:-1]
            deltap = self.deltam
            deltam = np.zeros_like(deltap)
            deltam[1:] = deltap[:-1]
            deltam[0] = np.NaN

            if params['weight']:
                # integrate against r
                # Mdiag = h * 2/3 * x
                # Moffdiag = h * 1/6 * (x + h/2)
                # Mdiag[0] = h**2 / 12

                # Ddiag = 2 / h * x #OK
                # Doffdiag = -1 / h * (x + h/2) #OK
                # Ddiag[0] = 1/2 #OK

                Mdiag = deltam/3*(x-deltam/4) + deltap/3*(x+deltap/4)
                Moffdiag = 1/6*deltap*(x+deltap/2)
                Mdiag[0] = deltap[0]**2/12

                Ddiag = x*(1/deltam+1/deltap)
                Doffdiag = -1/deltap*(x+deltap/2)
                Ddiag[0] = 1/2

                Delta = functions.tridiag(N,Doffdiag,Ddiag,Doffdiag)
                Mass = functions.tridiag(N,Moffdiag,Mdiag,Moffdiag)
            else:
                if l != 0:
                    raise Exception('Not implemented')
                Delta = functions.tridiag(N,-1,2,-1) / self.delta
                Mass = functions.tridiag(N,1,4,1) * self.delta / 6
                #neumann
                Delta[0,0] = 1/self.delta
                Mass[0,0] = 1/3*self.delta

        else:
            if params['fd']:
                Delta = functions.tridiag(N,-1,2,-1) / self.delta**2
                # Mass = functions.tridiag(N,0,1,0)
                Mass = None
                # Delta[0,0] = 1/self.delta**2
            else:
                # Delta = functions.tridiag(N,-1,2,-1) / self.delta
                # Mass = functions.tridiag(N,1,4,1) * self.delta / 6

                # deltamat = functions.tridiag(N,0, np.sqrt(self.delta),0)
                # deltainvmat = functions.tridiag(N,0, 1/np.sqrt(self.delta),0)
                # Delta = np.dot(np.dot(deltainvmat,functions.tridiag(N,-1,2,-1)),deltainvmat)
                # Mass = np.dot(np.dot(deltamat,functions.tridiag(N,1,4,1)),deltamat) / 6

                # if neumann:
                #     Delta[0,0] = 1/self.delta
                #     Mass[0,0] = 1/3*self.delta

                Ddiag = (self.deltap+self.deltam)/self.deltap/self.deltam
                Doffdiag = -1/self.deltap
                Mdiag = 1/3*(self.deltap + self.deltam)
                Moffdiag = 1/6 * self.deltap

                Delta = functions.tridiag(N,Doffdiag,Ddiag,Doffdiag)
                Mass = functions.tridiag(N,Moffdiag,Mdiag,Moffdiag)
            
        return Delta,Mass

    def fem_mult_op(self,V,l,Mass):
        # Vsqrtop = scipy.sparse.spdiags(np.sqrt(np.abs(V)),0,V.shape[0],V.shape[0],'csr')
        # return np.sign(V[1]) * np.dot(np.dot(Vsqrtop,Mass),Vsqrtop)

        params = self.fem_params(l)
        if not params['phi']:
            if not params['weight']:
                raise Exception("Not implemented")
            # the grid has been shifted
            x = np.zeros_like(self.x)
            x[1:] = self.x[:-1]
            deltap = self.deltam
            deltam = np.zeros_like(deltap)
            deltam[1:] = deltap[:-1]
            deltam[0] = np.NaN

            
            h = deltam[5]
            diag = V / 2 * x
            diag[1:] += 1/12*V[:-1] * (x[1:] - 2/5*h)
            diag[:-1] += 1/12*V[1:] * (x[:-1] + 2/5*h)
            diag[0] = h*(1/20*V[0] + 1/30 * V[1])
            offdiag = 1/12*V * (x + 2/5*h)
            offdiag[:-1] += 1/12*V[1:] * (x[:-1] + 3/5*h)
            # debug

            diag = V/4*(deltam*(x-1/5*deltam) + deltap*(x+1/5*deltap))
            diag[1:] += 1/12*V[:-1] * deltam[1:]*(x[1:] - 2/5*deltam[1:])
            diag[:-1] += 1/12*V[1:] * deltap[:-1]*(x[:-1] + 2/5*deltap[:-1])
            diag[0] = deltap[0]**2*(1/20*V[0] + 1/30 * V[1])
            offdiag = 1/12*V * deltap*(x + 2/5*deltap)
            offdiag[:-1] += 1/12*V[1:] * deltap[:-1]*(x[:-1] + 3/5*deltap[:-1])


            return functions.tridiag(self.N,offdiag,diag,offdiag)
        else:
            if params['fd']:
                return functions.tridiag(self.N,0,V,0)
            else:
                diag = V/4*(self.deltap + self.deltam)
                diag[1:] += 1/12*V[:-1]*self.deltam[1:]
                diag[:-1] += 1/12*V[1:]*self.deltap[:-1]
                offdiag = 1/12*V*self.deltap
                offdiag[:-1] += 1/12*V[1:]*self.deltap[:-1]
            
                return functions.tridiag(self.N,offdiag,diag,offdiag)

    def fem_overxsquare_op(self):
        x = self.x

        settings = np.geterrobj()
        np.seterr(all='ignore')
        diag = (2 * np.log((x - self.deltam)) * (x ** 2) - 2 * np.log((x - self.deltam)) * self.deltam * x + (2 * x * self.deltam) - (self.deltam ** 2) - 2 * np.log(x) * (x ** 2) + 2 * np.log(x) * self.deltam * x) / (self.deltam ** 2) / x
        diag[0] = ((2 * x[0] * self.deltam[0]) - (self.deltam[0] ** 2) - 2 * np.log(x[0]) * (x[0] ** 2) + 2 * np.log(x[0]) * self.deltam[0] * x[0]) / (self.deltam[0] ** 2) / x[0]
        np.seterrobj(settings)
        
        diag+= ((2 * x * self.deltap) + (self.deltap ** 2) + 2 * np.log(x) * self.deltap * x + 2 * np.log(x) * (x ** 2) - 2 * np.log((x + self.deltap)) * (x ** 2) - 2 * np.log((x + self.deltap)) * x * self.deltap) / (self.deltap ** 2) / x;

        offdiag =  - ((2 * self.deltap) + np.log(x) * self.deltap + 2 * np.log(x) * x - 2 * np.log((x + self.deltap)) * x - np.log((x + self.deltap)) * self.deltap) / (self.deltap ** 2);
    
        return functions.tridiag(self.N,offdiag,diag,offdiag)
