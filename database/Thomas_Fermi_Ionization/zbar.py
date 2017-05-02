import numpy as np
import math

#calculates Zbar of a mixture, based on glosli ddcmd calculation

Na = 6.0221415e23   #Avogadro/Avocado number

def zBar(n, Z, T):
    """
    This function calculates the partial ionizations of a mixture of species. 

    Inputs:
    n - list of number densities of the species. Units: 1/cc
    Z - list of charg states when each species is fully ionized (e.g. he has Z=2)
    T - electron temperature. Units: eV

    Output:
    zBar - list of partial ionizations for each species
    
    """
    nspec = len(n)

    Ions = []
    for spec in range(nspec):
        Ions.append(Ion_state(n[spec]/Na,Z[spec],T) )
    
    zBarFunc(Ions)

    zbar = []
    for Ion in Ions:
        zbar.append(Ion.zbar)

    return zbar

                  



class Ion_state(object):
    def __init__(self,n,Z,T):
        self.n = n
        self.T = T
        self.Z = Z
        self.vol = 0.0
        self.nf = 0.0
        self.dnf = 0.0
        self.zbar = 0.0
    
        a1 = 0.003323;
        a2 = 0.9718;
        a3 = 9.26148e-5;
        a4 = 3.10165;
        b0 = -1.7630;
        b1 = 1.43175;
        b2 = 0.31546;
        c1 = -0.366667;
        c2 = 0.983333;

        self.T = self.T + 1e-100; 
        T0 = T/pow(Z,(4./3.)); 
        Tf = T0/(1.+T0);
        self.A = a1*pow(T0,a2)+a3*pow(T0,a4);
        self.B = -math.exp(b0+b1*Tf+b2*pow(Tf,7.));
        self.C = c1*Tf+c2;



def fttfq(I):
    A = I.A
    B = I.B
    C = I.C
    Z = I.Z

    vol = I.vol

    alpha = 14.3139
    beta = 0.6624
    R = (Z*vol)
    Q1 = A*R**-B
    Y1 = R**-C 
    Y2 = Q1**C
    Y = Y1 + Y2; 
    Q = Y**(1.0/C)
    x = alpha*Q**beta
    zbar = Z*x/(1.+x+math.sqrt(1.+2.*x));
    
    dR = Z; 
    dQ1 = -B*Q1/R*dR; 
    dY = C*(-Y1/R*dR + Y2/Q1*dQ1); 
    dQ = Q/(C*Y)*dY; 
    dx = x * beta/Q*dQ; 

    t1 = math.sqrt(1.+2.*x);
    dt1 = 1/t1 * dx; 
    t2 = 1 + x + t1; 
    dt2  = (dx + dt1); 
    zBar =   Z*x/t2; 
    dzBar = zBar/x*dx - zBar/t2*dt2; 
    rho = zbar/vol; 
    drho = (dzBar - rho)/vol; 
    I.zbar = zBar; 
    I.nf = rho; 
    I.dnf = drho; 

    return rho;

    

def zBarFunc(Ions):

    nspec = len(Ions)
    ntotal = 0.0

    for Ion in Ions:
        ntotal += Ion.n

    vol = 1/ntotal
    nf = 0.0
    
    for Ion in Ions:
        Ion.vol = vol
        fttfq(Ion)
        nf = nf + Ion.nf * Ion.n
    
    nf = nf * vol
    
    for loop in range(50):
        f = 0.0
        df = 0.0
        error = 0.0
        for Ion in Ions:
            fttfq(Ion)
            delta = -(Ion.nf - nf)/Ion.dnf
            if (Ion.vol + delta) > 1.0/Ion.n:
                delta = 0.5*(1.0/Ion.n - Ion.vol)
            if (Ion.vol + delta < 0):
                delta = -0.5*Ion.vol
            f = f + Ion.n*Ion.vol
            df = df + Ion.n / Ion.dnf
            error = (Ion.nf / nf - 1.0)**2
        delta_nf = - (f-1.0)/df
        if(delta_nf > 0.2*nf):
            delta_nf = 0.2*nf
        if(delta_nf < -0.2*nf):
            delta_nf = -0.2*nf
        nf = nf + delta_nf
        
        error = error + (f-1.0)*(f-1.0)
        error = math.sqrt(error)
        if(error < 1.0e-12):
            break


