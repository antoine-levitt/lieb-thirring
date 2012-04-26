import shelve,pickle
bs = len(eigs)
filename = str(d) + "D/" + str(bs) + ".shelf"
#filename='7BS.shelf'
print "Saving into " + filename
my_shelf = shelve.open(filename,'n') # 'n' for new

to_save = ("g", "maxiter", "Vinit", "Vgammas","gammas", "Cgamma","eigs", "eigss")

for key in to_save:
    try:
        my_shelf[key] = globals()[key]
    # except (TypeError,KeyError,pickle.PicklingError):
    except:
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()
