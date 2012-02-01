#eigsN pour les vp de N points
Ns = np.floor(logspace(log10(4000),log10(8000),20))
#Ns = [2000, 4000,8000]
#Ns = np.floor(logspace(log10(100),log10(200),10))

eigss = []
psis_N = []
g_N = []
for N in Ns:
    execfile("1D.py")
    eigss.append(eigs)
    psis_N.append(psis)
    g_N.append(g)

i = 0

g = g_N[-1]

myeigs = array([eigs[i] for eigs in eigss])
myeigsdiff = abs(myeigs[1:] - myeigs[:-1])

myvecs = []
for j in range(0,len(Ns)):
    myvecs.append(np.interp(g_N[-1].x,g_N[j].x,psis_N[j][:,i]))
    #myvecs.append(functions.extrap(g_N[-1].x,g_N[j].x,phi))

myvecsdiff = [functions.compute_lp(g,myvecs[j+1] - myvecs[j],2) for j in range(0,len(Ns)-1)]
myvecsgraddiff = []
for j in range(0,len(Ns)-1):
    d = myvecs[j+1] - myvecs[j]
    diff = zeros_like(g.x)
    diff[:-1] = (d[1:] - d[:-1])/g.deltap[:-1]
    myvecsgraddiff.append(functions.compute_lp(g,diff,2))

myvecsh1diff = array([la.norm(array([myvecsdiff[j],myvecsgraddiff[j]])) for j in range(0,len(myvecsdiff))])

myvecsdiff = myvecsgraddiff

plot(Ns[:-1],myvecsdiff,'-x')
p = polyfit(log(Ns[:-1]),log(myvecsdiff),1)
print p[0]


plot(Ns[:-1],myeigsdiff,'-x')
p = polyfit(log(Ns[:-1]),log(myeigsdiff),1)
print p[0]

show()


# V0 = extrap(0,g.x,Vinit)
# Vbis = np.roll(Vinit,1)
# Vbis[0] = V0

# phi1 = myvecs[0]
# phi2 = myvecs[1]
# Delta,Mass = g.fem_mats(0)
# Vop = g.fem_mult_op(Vinit,0,Mass)
# H = Delta + Vop
# print dot(phi2,H.dot(phi2)) / dot(phi2,Mass.dot(phi2))
# print myeigs
