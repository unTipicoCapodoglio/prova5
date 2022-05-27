import numpy  as np

def quad(a,b, c = 0, d = 0, e = 0, f = 0, root = True):
    if (root):
        return (a**2 + b**2 + c**2 + d**2 + e**2 + f**2)**0.5
    else:
        return (a**2 + b**2 + c**2 + d**2 + e**2 + f**2)

def getV(x2,x1,t2,t1,sigma_x2,sigma_x1,sigma_t2,sigma_t1):
    sigma_v = quad(sigma_x2/(t2-t1),
                sigma_x1/(t2-t1),
                (x2-x1)*sigma_t2/(t2-t1)**2,
                (x2-x1)*sigma_t1/(t2-t1)**2)
    v = (x2-x1)/(t2-t1)
    return v, sigma_v
    