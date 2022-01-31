
import numpy as np
import scipy.integrate as integrate

def zbar(Z, AM, rho, T):
    """
    Finite Temperature Thomas Fermi Charge State using
    R.M. More, "Pressure Ionization, Resonances, and the
    Continuity of Bound and Free States", Adv. in Atomic
    Mol. Phys., Vol. 21, p. 332 (Table IV).

    Z = atomic number
    AM = atomic mass
    rho = density (g/cc)
    T = temperature (eV)
    """

    alpha = 14.3139
    beta = 0.6624
    a1 = 0.003323
    a2 = 0.9718
    a3 = 9.26148e-5
    a4 = 3.10165
    b0 = -1.7630
    b1 = 1.43175
    b2 = 0.31546
    c1 = -0.366667
    c2 = 0.983333

    R = rho/(Z*AM)
    T0 = T/Z**(4./3.)
    Tf = T0/(1+T0)
    A = a1*T0**a2+a3*T0**a4
    B = -np.exp(b0+b1*Tf+b2*Tf**7)
    C = c1*Tf+c2
    Q1 = A*R**B
    Q = (R**C+Q1**C)**(1/C)
    x = alpha*Q**beta

    return Z*x/(1 + x + np.sqrt(1 + 2.*x))



def FD_int(eta, order):
    '''Numerically computes the Fermi integral.'''
    return integrate.quad(lambda x: x**order/(1 + np.exp(x - eta)), 0, 100)[0];



def Ichimaru_chem_pot(theta):
    '''Ichimaru's fit to the chemical potential: Ichimaru, Volume II, page 87, (3.147)
    input: theta = T/E_F
    ouput: unitless ideal gas chemical potential
    '''
    A = 0.25954
    B = 0.072
    b = 0.858

    return -1.5*np.log(theta) + np.log(4/(3*np.sqrt(np.pi))) + (A*theta**(-b-1) + B*theta**(-b/2-1/2))/(1+ A*theta**(-b))



def A_alpha_fit(eta):
    '''Fit from the appendix of LM, used to check the numerical results.'''
    a_1 = 3.39
    a_2 = 3.47e-1
    a_3 = 1.29e-1
    b_2 = 5.11e-1
    b_3 = 1.24e-1

    y = np.log(1 + np.exp(eta))

    return (a_1 + a_2*y + a_3*y**2)/(1 + b_2*y + b_3*y**2)



def A_alpha(eta):
    '''Compute A^\alpha directly from numerical Fermi integrals.
    This is called D in the paper to avoid conflict with other uses of A.'''
    return (4./3)*FD_int(eta,2)/((1 + np.exp(-eta))*FD_int(eta,0.5)**2)



def effective_temperature(T_e, eta):
    '''Return the "Thomas-Fermi temperature", which yields
    T_e and (2/3)E_F in the classical and degenerate limits, respectively.
    T_e isin eV and eta is dimensionless.'''
    return 2.0*T_e*FD_int(eta,0.5)/FD_int(eta,-0.5)



def LM_Coulomb_logarithm(n_e, T_e, n_i, T_i, eta, z_bar):
    '''This is the CL from LM, to my best guess of what they did.'''

    # b_max
    # It is not perfectly clear what LM meant by "Fermi temperature", so I used this:
    lambda_e = 1/np.sqrt(4*np.pi*n_e*1.44e-7/effective_temperature(T_e, eta)) # cm
    lambda_i = 1/np.sqrt(4*np.pi*n_i*1.44e-7/T_i) # cm
    a_i = (3/(4*np.pi*n_i))**(1/3) # I use this (a_i) for their R_0, since they didn't specify it.
    b_max = np.max([a_i, 1.0/np.sqrt(1/lambda_e**2 + 1/lambda_i**2)]) # LM rule for minimum CL

    # b_min
    b_min = np.min([z_bar*1.44e-7/(3.0*T_e), 5e-8/np.sqrt(T_e)]) # units are cm
#     b_min = np.min([z_bar*1.44e-7/(3.0*effective_temperature(T_e, eta)), 5e-8/np.sqrt(effective_temperature(T_e, eta))]) # what if LM used T_eff?


    return np.max([0.5*np.log(1 + b_max**2/b_min**2), 2.0])



def MSM_Coulomb_logarithm(n_e, T_e, n_i, T_i, eta, z_bar):
    '''
    This is my modification to the LM CL. Potentially important
    changes are made, such as using the effective electron temperature
    consistently.


    Note that LM state "In our model we require the Coulomb logarithm to have value greater
    than or equal to 2.0." This condition is NOT used here.
    '''

    # b_max
    lambda_e = 1/np.sqrt(4*np.pi*n_e*1.44e-7/effective_temperature(T_e, eta)) # cm
    lambda_i = 1/np.sqrt(4*np.pi*n_i*1.44e-7/T_i) # cm
    a_e = (3/(4*np.pi*n_e))**(1/3)
    a_i = (3/(4*np.pi*n_i))**(1/3)
    b_max = 1/np.sqrt(1/(a_e**2 + lambda_e*2) + 1/(a_i**2 + lambda_i**2))

    # b_min
    b_min = np.sqrt((z_bar*1.44e-7/(3.0*effective_temperature(T_e, eta)))**2 + (5e-8/np.sqrt(effective_temperature(T_e, eta)))**2)

    return 0.5*np.log(1 + b_max**2/b_min**2)



def sigma_MSM(electron_temperature, ion_temperature, ion_density, eta, z_bar):
    '''Lee-More electrical conductivity, using my choices.
    inputs: temperatures in eV, density in 1/cc, others dimensionless
    Requires functions A^\alpha(eta) and TF_zbar.
    Base units are 1/s, or 1/(9*10^9 Ohm-m), but currently
    returns 1/(Ohm-cm).'''

    sigma_0 = ion_density*z_bar*1.44e-7*(2.997925e10)**2/511e3*A_alpha(eta) # units are 1/s^2
    degen_factor =  (1 + np.exp(-eta))*FD_int(eta, 0.5)

    Coulomb_log = MSM_Coulomb_logarithm(z_bar*ion_density, electron_temperature, ion_density, ion_temperature, eta, z_bar)

    tau_temporary = 3*np.sqrt(511e3)*electron_temperature**(3/2)/(2*np.sqrt(2)*2.997925e10*np.pi*z_bar**2*ion_density*(1.44e-7)**2*Coulomb_log) # units of s
    tau = tau_temporary*degen_factor
    sigma = sigma_0*tau # units are 1/s

    unit_factor = 9e11 # convert units to 1/(Ohm-cm)

    return sigma/unit_factor



def sigma_LM(electron_temperature, ion_temperature, ion_density, eta, z_bar):
    '''Lee-More electrical conductivity, to my best guess of what they did.
    inputs: temperatures in eV, density in 1/cc, others dimensionless
    Requires functions A^\alpha(eta) and TF_zbar.
    Base units are 1/s, or 1/(9*10^9 Ohm-m), but currently
    returns 1/(Ohm-cm).'''

    sigma_0 = ion_density*z_bar*1.44e-7*(2.997925e10)**2/511e3*A_alpha(eta) # units are 1/s^2
    degen_factor =  (1 + np.exp(-eta))*FD_int(eta, 0.5)

    Coulomb_log = LM_Coulomb_logarithm(z_bar*ion_density, electron_temperature, ion_density, ion_temperature, eta, z_bar)

    tau_temporary = 3*np.sqrt(511e3)*electron_temperature**(3/2)/(2*np.sqrt(2)*2.997925e10*np.pi*z_bar**2*ion_density*(1.44e-7)**2*Coulomb_log) # units of s
    tau = tau_temporary*degen_factor
    sigma = sigma_0*tau # units are 1/s

    unit_factor = 9e11 # convert units to 1/(Ohm-cm)

    return sigma/unit_factor



def sigma(Z, A, rho, T):
    '''Call this function to get a conductivity prediction. This is a helper
    function that uses the functions above so that an 'external' user does not
    deal with <Z>, chemical potential, and so on.
    Z - no units
    rho - g/cc
    A - no units
    T - eV
    '''

    m_proton = 1.672622e-24
    electron_temperature = T
    ion_temperature = T
    ion_density = rho/(A*m_proton)

    z_bar = zbar(Z, A, rho, T)
    n_e = z_bar*ion_density
    E_F = (0.197326e-4)**2*(3*np.pi**2*n_e)**(2/3)/(2*511e3)
    theta = T/E_F
    eta = Ichimaru_chem_pot(theta)

    return sigma_MSM(electron_temperature, ion_temperature, ion_density, eta, z_bar)
