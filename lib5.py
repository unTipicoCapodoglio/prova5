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


def getM_a(s, r, r_star, r_f, m_b, rho, B,sigma_s, sigma_r, sigma_r_star, sigma_r_f, sigma_m_b, sigma_rho, sigma_B):
    _pi, _d, _g = sp.symbols('pi, d,g')
    _s, _r, _r_star, _r_f, _m_b, _rho, _B =sp.symbols('s, r, r^*, r_f, m_b, rho, B ')
    _sigma_s, _sigma_r, _sigma_r_star, _sigma_r_f, _sigma_m_b, _sigma_rho, _sigma_B  = sp.symbols('sigma_s, sigma_r, sigma_r^*, sigma_r_f, sigma_m_b, sigma_rho, sigma_B')
    _M_a, _sigma_M_a = sp.symbols('M_a, sigma_M_a')
    _M_a = _s*_pi*_d*_g*_r*(((_r_star+_r)**2)-12*(_r_f**2)) - (2*_m_b*((_r+_rho+_r_f)**2))/(_B*_r)
    _sigma_M_a = quad(sp.diff(_M_a,_s)*_sigma_s,
                      sp.diff(_M_a,_r)*_sigma_r,
                      sp.diff(_M_a,_r_star)*_sigma_r_star,
                      sp.diff(_M_a,_r_f)*_sigma_r_f,
                      sp.diff(_M_a,_m_b)*_sigma_m_b,
                      sp.diff(_M_a,_rho)*_sigma_rho,
                      sp.diff(_M_a,_B)*_sigma_B)

    M_a = _M_a.subs([(_pi,np.pi),(_d,2.70/1000),(_g,9.805),(_s,s),(_r,r),(_r_star,r_star),(_r_f,r_f),(_m_b,m_b),(_rho,rho),(_B,B)])
    sigma_M_a = _sigma_M_a.subs([(_pi,np.pi),(_d,2.70/1000),(_g,9.805),(_s,s),(_r,r),(_r_star,r_star),(_r_f,r_f),(_m_b,m_b),(_rho,rho),(_B,B),
                                 (_sigma_s,sigma_s),(_sigma_r,sigma_r),(_sigma_r_star,sigma_r_star),(_sigma_r_f,sigma_r_f),(_sigma_m_b,sigma_m_b),
                                 (_sigma_rho,sigma_rho),(_sigma_B, sigma_B)])
    #Nonostante i risulati di sopra di V siano espressi in cm^3, il notebook ha comunque lavorato sempre in g/mm^3
    #La conversione g->Kg avviene ora, insieme a far diventare le espresisoni di sympy dei float di Python
    M_a = float(M_a)/1000 #M_a in questo caso è espresso in g*m^2/s^2, ci serve in N*m 
    sigma_M_a = float(sigma_M_a)/1000 # M_a in questo caso è espresso in g*m^2/s^2, ci serve in N*m 
    
    latex1 = sp.latex(_M_a)
    latex2 = sp.latex(sp.simplify(_sigma_M_a))
    return M_a, sigma_M_a, latex1, latex2