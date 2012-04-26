pylab.rcParams.update({'text.usetex':True,
                       'axes.labelsize': 18,
                       'text.fontsize': 14,
                       'legend.fontsize': 14,
                       'xtick.labelsize': 14,
                       'ytick.labelsize': 14})
figure()
contourf(g.x,g.y,V)
gray()
colorbar()
lim = 35
xlim((-lim,lim))
ylim((-lim,lim))
xlabel('$x$')
ylabel('$y$')
pylab.legend()
pylab.savefig('/home/antoine/Desktop/fig.pdf')
show()
