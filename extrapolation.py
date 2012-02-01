def aitken(x):
    a = x[:-2]
    b = x[1:-1]
    c = x[2:]
    
    return (c*a-b**2)/(-2*b+c+a)

def richardson(x):
    a = x[:-2]
    b = x[1:-1]
    c = x[2:]

    r = (c-b)/(b-a)
    return b + (c-b)/(1-r), r
