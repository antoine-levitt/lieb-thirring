#exploit, 2D avec branche 4
plot(shelf["gammas"], shelf["Cgamma"], 'kx',markersize=6)
xlabel('$\gamma$')
ylabel('$E(V_\gamma)$')
xlim((.8,1.5))
ylim((.95,1.1))
pylab.rcParams.update({'text.usetex':True,
                       'axes.labelsize': 20,
                       'text.fontsize': 16,
                       'legend.fontsize': 16,
                       'xtick.labelsize': 16,
                       'ytick.labelsize': 16})
# pylab.savefig('plots/bifur_R.pdf')
show()


figure()
[plot(shelf["g"].x,V) for V in shelf["Vgammas"][70:74]]
xlabel('$r$')
ylabel('$V(r)$')
xlim((0,80))
legend(['$\gamma = {0:.2f}$'.format(shelf["gammas"][g]) for g in arange(70,74)])
# pylab.savefig('plots/bifur_V.pdf')
show()
