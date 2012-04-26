from __future__ import division
Ls = linspace(10,50,20)

eigss = []
psis_N = []
g_N = []
C_N = []
oneoverh = 2000
for L in Ls:
    N = L*oneoverh
    N = floor(N/2) * 2
    execfile("1D.py")
    eigss.append(eigs)
    psis_N.append(psis)
    g_N.append(g)
    C_N.append(Cgamma[0])
    init = "prev"


# i = 0

Lfinal = Ls[-1] * 2
L = Lfinal
N = L*oneoverh
N = floor(N/2) * 2
execfile("1D.py")
psi = psis[:,i]

myvecs = []
for j in range(0,len(Ns)):
    # myvecs.append(np.interp(g_N[-1].x,g_N[j].x,psis_N[j][:,i]))
    myvecs.append(np.interp(g.x,g_N[j].x,psis_N[j][:,i]))
    #myvecs.append(functions.extrap(g_N[-1].x,g_N[j].x,phi))

myvecsdiff = [functions.compute_lp(g,myvecs[j+1] - myvecs[j],2) for j in range(0,len(Ns)-1)]
myvecsdiff = [functions.compute_lp(g,v - psi,2) for v in myvecs]
# myvecsgraddiff = []
# for j in range(0,len(Ns)-1):
#     d = myvecs[j+1] - myvecs[j]
#     diff = zeros_like(g.x)
#     diff[:-1] = (d[1:] - d[:-1])/g.deltap[:-1]
#     myvecsgraddiff.append(functions.compute_lp(g,diff,2))
myvecsgraddiff = []
for j in range(0,len(Ns)):
    d = myvecs[j] - psi
    diff = zeros_like(g.x)
    diff[:-1] = (d[1:] - d[:-1])/g.deltap[:-1]
    myvecsgraddiff.append(functions.compute_lp(g,diff,2))

myvecsh1diff = array([la.norm(array([myvecsdiff[j],myvecsgraddiff[j]])) for j in range(0,len(myvecsdiff))])

myvecsdiff = myvecsgraddiff
myvecsdiff = myvecsh1diff
# myvecsdiff = myvecsdiff




pylab.rcParams.update({'text.usetex':False,
                       'axes.labelsize': 18,
                       'text.fontsize': 14,
                       'legend.fontsize': 14,
                       'xtick.labelsize': 14,
                       'ytick.labelsize': 14})



Ls = Ls/2

semilogy(Ls,myvecsdiff,'-x')
p = polyfit(Ls,log(myvecsdiff),1)
print "psi", p[0]


# plot(Ns[:-1],myeigsdiff,'-x')
# p = polyfit(log(Ns[:-1]),log(myeigsdiff),1)
# print p[0]

show()

semilogy(Ls,one_d_one_bs(gamma) - C_N,'--s')
p = polyfit(Ls,log(one_d_one_bs(gamma) - C_N),1)
print "C", p[0]


show()

xlabel("$L$")
ylabel('Error')
legend(("$||\psi_L - \psi||_{H^1}$", "$|E(V_L) - L_{\gamma}|$"))
savefig("plots/conv_L.pdf")

show()
