import shelve
def load_shelf():
    my_shelf = shelve.open("shelve.out")
    for key in my_shelf:
        globals()[key]=my_shelf[key]

    print my_shelf['g']
    my_shelf.close()
