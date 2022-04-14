import numpy as np
from scipy.integrate import simps
from scipy.optimize import root_scalar, fsolve

from ..constants import speed_of_light, omega_matter, omega_lamda, H0
cc = speed_of_light*1e-3 # Km/s

def distance_from_redshift(z):
    """
    Returns the luminosity distance [Mpc] from the redshift.
    Uses cosmological parameters from the constants.py module.

    :param float z: Redshift.

    :rtype: float
    """
    integrand = lambda x: 1/np.sqrt(omega_matter*(1+x)**3+omega_lamda)
    dz = 1e-5
    X = np.arange(0,z+dz,dz)
    Y = integrand(X)
    d_L = cc/H0*(1+z)*simps(Y,X)
    return d_L

def redshift_from_distance(d_L):
    """
    Returns the redshift from the luminosity distance [Mpc]. 
    Uses cosmological parameters from the constants.py module.

    :param float d_L: Luminosity distance [Mpc].
    
    :rtype: float
    """
    if d_L<=2:
        return d_L*H0/cc
    f = lambda x: d_L - distance_from_redshift(x)
    z = fsolve(f,10).item()
    return z
