"""Class to find mean ionizations in plasmas, as found in [1].

[1]Stanton, L. G., & Murillo, M. S. (2016). 
Ionic transport in high-energy-density matter. 
Physical Review E, 93(4), 043203. 
https://doi.org/10.1103/PhysRevE.93.043203"""

import numpy as np


class ZbarSolver(object):
    """A simple class that iteratively solves for mean ionization states with
    the ability to be reset.
    
    This class solves for the Thomas-Fermi average atom ionization of a
    plasma system, and has the ability to be reset with new inputs inorder
    to avoid re-initialization when a large number of solves need to be done.
    
    Parameters
    ----------
    ndens : (N,) nd.array
        The number densities (in units of 1/cm^{3}) of each species.
    znuc : (N,) nd.array
        Atomic number of each species.
    Te : float
        Electron temperature in the system, in units of eV.


    Class Attributes
    ----------------
    a : (2,) nd.array
       Fit constants.
    k : (9,) nd.array
       Other fit constants.
    

    Attributes
    ----------
    zbars : (N,) nd.array
       Ionization of each species, in the same order as the inputs.
    A : (N,) nd.array
       Invariant of the mean ionization function with respect to inputs.
    B : (N, nd.array)
       Invariant of the mean ionization function with respect to inputs.
    C : (N, nd.array)
       Invariant of the mean ionization function with respect to inputs.
    s1 : (N, nd.array)
       Storage array. Once the solution has converged, holds volume 
       (in units of cm^{3}/mol) of each species. 
    s2 : (N, nd.array)
       Storage array. Once the solution has converged, it holds the
       mean ionization from the previous step.
    """
    a = np.array([14.3139, 0.6624], dtype=np.double)
    k = np.array([3.323e-3, 0.9718, 9.26148e-5, 3.10165, -1.7630, 
                  1.43175, 0.31546, -0.366667, 0.983333], dtype=np.double)


    def __init__(self, ndens, znuc, Te):
        super(ZbarSolver, self).__init__()
        self.zbars = np.ones(znuc.shape[0], dtype=znuc.dtype)
        self.A = np.ones(znuc.shape[0], dtype=znuc.dtype)
        self.B = np.ones(znuc.shape[0], dtype=znuc.dtype)
        self.C = np.ones(znuc.shape[0], dtype=znuc.dtype)
        self.s1 = np.ones(znuc.shape[0], dtype=znuc.dtype)
        self.s2 = np.ones(znuc.shape[0], dtype=znuc.dtype)
        self.solve_zbar(ndens/6.0221409e23, znuc, Te)

    
    def reset(self, ndens, znuc, Te):
        """Essentially re-initializes the class, but without additional memory
        allocation or initialization overhead.
        
        Parameters
        ----------
        ndens : (N,) nd.array
            The number densities (in units of 1/cm^{3}) of each species.
        znuc : (N,) nd.array
            Atomic number of each species.
        Te : float
            Electron temperature in the new system, in units of eV."""
        self.zbars[:] = 1.0
        self.solve_zbar(ndens/6.0221409e23, znuc, Te)


    def set_constants(self, znuc, Te):
        """Set variables in the zbar function that are independent of 
        the varying ionization states.
        
        Parameters
        ----------
        znuc : (N,) nd.array
            Atomic number of each species.
        Te : float
            Electron temperature in the system.
        """
        self.A = Te*np.power(znuc, -4.0/3.0)
        self.B = self.A/(1.0 + self.A)
        self.C = np.copy(self.B)
        self.A = self.k[0]*np.power(self.A, self.k[1]) + self.k[2]*np.power(self.A, self.k[3])
        self.B = self.k[4] + self.k[5]*self.B + self.k[6]*np.power(self.B, 7.0)
        self.B = -1*np.exp(self.B)
        self.C = self.k[7]*self.C + self.k[8]

    
    def set_zbar(self, znuc):
        """Generates the mean ionization states for the current values of `self.V`.
        
        Parameters
        ----------
        znuc : (N,) nd.array
            Atomic number of each species.
        """
        self.s1[:] = np.power(self.s1, -1*self.C)
        self.s1[:] = self.s1[:] + np.power(self.A[:], self.C[:])\
                    *np.power(self.s1[:], self.B[:])
        self.s1[:] = self.a[0]*np.power(self.s1[:], self.a[1]/self.C[:])
        self.zbars[:] = 1 + self.s1[:] + np.sqrt(1 + 2*self.s1[:])
        self.zbars[:] = znuc*self.s1[:]/self.zbars[:]


    def solve_zbar(self, ndens, znuc, Te):
        """Generate the consistent mean-ionizations for each species
        iteratively.
        
        Parameters
        ----------
        ndens : (N,) nd.array
            Number density of each species, in units of mol/cm^{3}.
        znuc : (N,) nd.array
            Atomic number of each species.
        Te : float
            Electron temperature, in units of eV.
        """
        self.set_constants(znuc, Te)
        self.s2[:] = self.zbars[:]
        self.s1[:] = znuc[:]*self.s1[:]
        self.set_zbar(znuc)
        rhotot = np.sum(ndens*self.zbars)
        self.s1[:] = self.zbars[:]/rhotot
        err = np.max(np.abs(self.zbars - self.s2))
        while err>np.finfo(self.zbars.dtype).eps:
            self.s2[:] = self.zbars[:]
            self.s1[:] = znuc[:]*self.s1[:]
            self.set_zbar(znuc)
            rhotot = np.sum(ndens*self.zbars)
            self.s1[:] = self.zbars[:]/rhotot
            err = np.max(np.abs(self.zbars - self.s2))