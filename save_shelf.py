import shelve,pickle
filename='7BS.shelf'
my_shelf = shelve.open(filename,'n') # 'n' for new

to_save = ("g", "maxiter", "Vinit", "Vgammas","gammas", "Cgamma","eigs")

for key in to_save:
    try:
        my_shelf[key] = globals()[key]
        print "saved", key
    # except (TypeError,KeyError,pickle.PicklingError):
    except:
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()
