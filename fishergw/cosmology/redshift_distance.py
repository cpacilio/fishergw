import numpy as np
from scipy.integrate import simps
from scipy.optimize import root_scalar, fsolve

from ..constants import speed_of_light, omega_matter, omega_lamda, H0
cc = speed_of_light*1e-3 # Km/s

def distance_from_redshift(z):
    """
    Returns the luminosity distance [Mpc] from the redshift.
   
    Uses the Planck 2018 cosmological parameters from Tab. 1 in https://arxiv.org/abs/1807.06209.

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

    Uses the Planck 2018 cosmological parameters from Tab. 1 in https://arxiv.org/abs/1807.06209.

    :param float d_L: Luminosity distance [Mpc]. Must be no less than 1.
    
    :rtype: float
    """
    if d_L<1:
        raise ValueError('d_L must be no less than 1 Mpc!')
    f = lambda x: d_L - distance_from_redshift(x)
    z = fsolve(f,10).item()
    return z
