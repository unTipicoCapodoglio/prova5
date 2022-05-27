import matplotlib.pyplot as plt
import numpy as np
from math import log10, floor, sqrt
import scipy as sp
from scipy import stats

#########################################################################
# funzione per arrotondare con un certo numero di cifre significative
#########################################################################
def PrintResult(name,mean,sigma,digits,unit):
    mean = round(mean,digits)
    sigma = round(sigma,digits)
    nu = sigma / mean
    result = (name+" = ({0} +/- {1} ) ".format(mean,sigma)+unit+" [{0:.2f}%]".format(nu*100))
    print (result)
    #return ""
    
#########################################################################
# funzioni per fare il fit lineare
#########################################################################
#
def my_mean(x, w):
    return np.sum( x*w ) / np.sum( w )

def my_cov(x, y, w):
    return my_mean(x*y, w) - my_mean(x, w)*my_mean(y, w)

def my_var(x, w):
    return my_cov(x, x, w)

def my_line(x, m=1, c=0):
    return m*x + c

def y_estrapolato(x, m, c, sigma_m, sigma_c, cov_mc):
    y = m*x + c
    uy = np.sqrt(np.power(x, 2)*np.power(sigma_m, 2) +
                   np.power(sigma_c, 2) + 2*x*cov_mc ) 
    return y, uy

def lin_fit(x, y, sd_y, sd_x=0, m_0 = 0., xlabel="x [ux]", ylabel="y [uy]", xm=0., xM=1., ym=0., yM=1.,
            name = "name", title = 'Fit Lineare', verbose=True, plot=True, setrange=False, save=False):
    s = np.arange(float(sd_y.size))
    #pesi
    #
    s = ((sd_x*m_0)**2 + sd_y**2)**0.5
    w_y = np.power(s.astype(float), -2)
    
    #m
    m = my_cov(x, y, w_y) / my_var(x, w_y)
    var_m = 1 / ( my_var(x, w_y) * np.sum(w_y) )
    
    #c
    c = my_mean(y, w_y) - my_mean(x, w_y) * m
    var_c = my_mean(x*x, w_y)  / ( my_var(x, w_y) * np.sum(w_y) )
    
    #cov
    cov_mc = - my_mean(x, w_y) / ( my_var(x, w_y) * np.sum(w_y) ) 
   
    #rho
    rho_mc = cov_mc / ( sqrt(var_m) * sqrt(var_c) )

    if (verbose):
        
        print ('m         = ', m.round(4))
        print ('sigma(m)  = ', np.sqrt(var_m).round(4))
        print ('c         = ', c.round(4))
        print ('sigma(c)  = ', np.sqrt(var_c).round(4))
        print ('cov(m, c) = ', cov_mc.round(4))
        print ('rho(m, c) = ', rho_mc.round(4))
        
    if (plot):
        
        # rappresento i dati
        plt.errorbar(x, y, yerr=sd_y, xerr=sd_x, ls='', marker='.', 
                     color="black",label='dati')

        # costruisco dei punti x su cui valutare la retta del fit              
        xmin = float(np.min(x)) 
        xmax = float(np.max(x))
        xmin_plot = xmin-.2*(xmax-xmin)
        xmax_plot = xmax+.2*(xmax-xmin)
        if (setrange):
            xmin_plot = xm
            xmax_plot = xM  
        x1 = np.linspace(xmin_plot, xmax_plot, 100)
        y1 = my_line(x1, m, c)
        
        # rappresento la retta del fit
        plt.plot(x1, y1, linestyle='-.', color="green", label='fit')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xlim(xmin_plot,xmax_plot)
        if (setrange):
            plt.ylim(ym,yM)
        
        # rappresento le incertezze sulla retta 
        y1_plus_1sigma = y1+y_estrapolato(x1, m, c, np.sqrt(var_m), np.sqrt(var_c), cov_mc)[1]
        y1_minus_1sigma = y1-y_estrapolato(x1, m, c, np.sqrt(var_m), np.sqrt(var_c), cov_mc)[1]         
        plt.plot(x1,y1_plus_1sigma, linestyle='-', color="orange", label=r'fit $\pm 1\sigma$')
        plt.plot(x1,y1_minus_1sigma, linestyle='-', color="orange")
        
        plt.grid()
        
        plt.legend()
        
    if (save):
        plt.savefig(name)
        
    return m, np.sqrt(var_m), c, np.sqrt(var_c), cov_mc, rho_mc

def help_residui():
    print('''x,y,m,c,sigma_y, xlabel = 'x [ux]', ylabel='y [uy]', title = 'title', name = 'residui.png',\n
            normalize = True, save = False''')

def residui(x,y,m,c, sigma_y, sigma_x = 0, xlabel = 'x [ux]', ylabel='y [uy]', title = 'title', name = 'residui.pdf',
            normalize = True, save = False):
    #plt.plot([x.min(), x.max()],[0,0], '-', color = 'orange')  
    if (normalize):
        d = (y - (m*x + c))/sigma_y
        plt.errorbar(x,d, xerr = sigma_x, yerr = 1, marker = '.', ls='')
        if (title == 'title'):
            title = 'Residui normalizzati'                  
    else:
        d = y - (m*x + c)
        plt.errorbar(x,d, xerr = sigma_x, yerr = sigma_y, marker = '.', ls='')       
        if (title == 'title'):
            title = 'Residui'        
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    if (save):
        plt.savefig(name)
    
