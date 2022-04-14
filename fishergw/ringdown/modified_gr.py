## contains functions to fit Kerr qnms, amplitudes and phases

import numpy as np
from ..constants import speed_of_light, solar_mass, G, Mpc

cc = speed_of_light
msun = solar_mass
## solar mass in secs
tsun = msun*G/cc**3 ## secs


def qnm_KerrNewman(mass,spin,charge,mode,method='fit2'):
    """
    Returns the frequency and the damping time of a Kerr-Newman black hole
    for a given mass, spin, charge-to-mass ratio and mode.

    :param mass: Mass of the black hole [:math:`M_\odot`].
    :param type: float

    :param spin: Dimensionless spin.
    :param type: float

    :param charge: dimensionless charge-to-mass ratio.
    :param type: float

    :param mode: Angular, azimuthal and overtone number in the form (l,m,n).
    :param type: tuple

:param method: The method used to approximate the Kerr-Newman spectrum. If `fit1`, it uses the fits provided in https://arxiv.org/abs/1307.7315. If `fit2`, it uses the fits provided in the appendix of https://arxiv.org/abs/2109.13961. Deafult: `fit2`.
    """
    l,m,n = mode

    amax = (1-charge**2)**0.5
    y = 1-amax

    methods = ['fit1','fit2','interp']
    if method not in methods:
        raise ValueError('method should be one of '+str(methods))

    if method=='fit1':
        ## use fits as provided in https://arxiv.org/abs/1307.7315
        # add non-rotating contribution
        coeffs_R = {(2,0):(0.3737,0.0525,0.0607,-0.0463,-0.0070),\
                  (2,1):(0.3467,0.0546,0.0709,-0.0292,-0.0433),\
                  (3,0):(0.5994,0.0790,0.1734,-0.2019,0.0700),\
                  (3,1):(0.5826,0.0819,0.1752,-0.1753,0.0304)}

        coeffs_I = {(2,0):(-0.0890,-0.0055,0.0024,0.0214,-0.0084),\
                  (2,1):(-0.2739,-0.0157,0.0099,0.0668,-0.0239),\
                  (3,0):(-0.0927,-0.0043,-0.0013,0.0292,-0.0130),\
                  (3,1):(-0.2813,-0.0123,-0.0050,0.0972,-0.0472)}

        f0, f1, f2, f3, f4 = coeffs_R[(l,n)]
        q0, q1, q2, q3, q4 = coeffs_I[(l,n)]
        omegaR = f0 + f1*y + f2*y**2 + f3*y**3 + f4*y**4
        omegaI = q0 + q1*y + q2*y**2 + q3*y**3 + q4*y**4

        # add rotating contribution
        coeffs_R = {(2,0):(0.0628,0.0676,0.0209,0.0823,-0.0810),\
                  (2,1):(0.0717,0.0764,0.0020,0.1959,-0.2213),\
                  (3,0):(0.0673,0.0693,0.0211,0.0791,-0.0677),\
                  (3,1):(0.0713,0.0736,0.0108,0.1287,-0.1282)}

        coeffs_I = {(2,0):(0.0010,0.0014,0.0091,0.0174,0.0145),\
                  (2,1):(0.0065,0.0070,0.0360,0.0254,0.0905),\
                  (3,0):(0.0006,0.0014,0.0084,0.0058,0.0122),\
                  (3,1):(0.0029,0.0059,0.0194,0.0374,0.0231)}

        f0, f1, f2, f3, f4 = coeffs_R[(l,n)]
        q0, q1, q2, q3, q4 = coeffs_I[(l,n)]
        omegaR += spin*m*(f0 + f1*y + f2*y**2 + f3*y**3 + f4*y**4)
        omegaI += spin*m*(q0 + q1*y + q2*y**2 + q3*y**3 + q4*y**4)

        # f and tau
        f = omegaR/(2*np.pi)/mass/tsun
        tau = -1/omegaI*mass*tsun

    elif method=='fit2':
        ## use fits provided in https://arxiv.org/abs/2109.13961
        BR = {(2,2,0): 0.37367168*np.array([\
                [1,0.537583,-2.990402,1.503421],\
                [-1.899567,-2.128633,6.626680,-2.903790],\
                [1.015454,2.147094,-4.672847,1.891731],\
                [-0.111430,-0.581706,1.021061,-0.414517]]),\
                ##
                (2,2,1): 0.34671099*np.array([\
                [1,-2.918987,2.866252,-0.944554],\
                [-1.850299,7.321955,-8.783456,3.292966],\
                [0.944088,-5.584876,7.675096,-3.039132],\
                [-0.088458,1.198758,-1.973222,0.838109]]),\
                ##
                (3,3,0): 0.59944329*np.array([\
                [1,-0.311963,-1.457057,0.825692],\
                [-1.928277,-0.026433,3.139427,-1.484557],\
                [1.044039,0.545708,-2.188569,0.940019],\
                [-0.112303,-0.226402,0.482482,-0.204299]])}
        ##
        CR = {(2,2,0): np.array([[1,0.548651,-3.141145,1.636377],\
                [-2.238461,-2.291933,7.695570,-3.458474],\
                [1.581677,2.662938,-6.256090,2.494264],\
                [-0.341455,-0.930069,1.688288,-0.612643]]),\
                ##
                (2,2,1): np.array([[1,-2.941138,2.907859,-0.964407],\
                [-2.250169,8.425183,-9.852886,3.660289],\
                [1.611393,-7.869432,9.999751,-3.737205],\
                [-0.359285,2.392321,-3.154979,1.129776]]),\
                ##
                (3,3,0): np.array([[1,-0.299153,-1.591595,0.938987],\
                [-2.265230,0.058508,3.772084,-1.852247],\
                [1.624332,0.533096,-3.007197,1.285026],\
                [-0.357651,-0.300599,0.810387,-0.314715]])}
        ##
        BI = {(2,2,0): 0.08896232*np.array([\
                [1,-2.721789,2.472860,-0.750015],\
                [-2.533958,7.181110,-6.870324,2.214689],\
                [2.102750,-6.317887,6.206452,-1.980749],\
                [-0.568636,1.857404,-1.820547,0.554722]]),\
                ##
                (2,2,1): 0.27391488*np.array([\
                [1,-3.074983,3.182195,-1.105297],\
                [0.366066,4.296285,-9.700146,5.016955],\
                [-3.290350,-0.844265,9.999863,-5.818349],\
                [1.927196,-0.401520,-3.537667,2.077991]]),\
                ##
                (3,3,0): 0.09270305*np.array([\
                [1,-2.813977,2.666759,-0.850618],\
                [-2.163575,6.934304,-7.425335,2.640936],\
                [1.405496,-5.678573,6.621826,-2.345713],
                [-0.241561,1.555843,-1.890365,0.637480]])}
        ##
        CI = {(2,2,0): np.array([[1,-2.732346,2.495049,-0.761581],\
                [-2.498341,7.089542,-6.781334,2.181880],\
                [2.056918,-6.149334,6.010021,-1.909275],\
                [-0.557557,1.786783,-1.734461,0.524997]]),\
                ##
                (2,2,1): np.array([[1,-3.079686,3.191889,-1.110140],\
                [0.388928,4.159242,-9.474149,4.904881],\
                [-3.119527,-0.914668,9.767356,-5.690517],\
                [1.746957,-0.240680,-3.505359,2.049254]]),\
                ##
                (3,3,0): np.array([[1,-2.820763,2.680557,-0.857462],\
                [-2.130446,6.825101,-7.291058,2.583282],\
                [1.394144,-5.533669,6.393699,-2.254239],\
                [-0.261229,1.517744,-1.810579,0.608393]])}
        ##
        qq = np.array([charge**i for i in range(4)])
        ss = np.array([spin**i for i in range(4)])
        #omegaR = np.sum(np.dot(ss,np.dot(BR[mode],qq)))/np.sum(np.dot(ss,np.dot(CR[mode],qq)))
        #omegaI = np.sum(np.dot(ss,np.dot(BI[mode],qq)))/np.sum(np.dot(ss,np.dot(CI[mode],qq)))
        omegaR = np.sum(ss*np.dot(BR[mode],qq))/np.sum(ss*np.dot(CR[mode],qq))
        omegaI = np.sum(ss*np.dot(BI[mode],qq))/np.sum(ss*np.dot(CI[mode],qq))
        # f and tau
        f = omegaR/(2*np.pi)/mass/tsun
        tau = 1/omegaI*mass*tsun

    else:
        raise KeyError('Only the fit methods are implemented!')

    return f, tau

def qnm_EdGB(mass,spin,zeta,mode,method='fit1'):
    """
    Returns the frequency and the damping time of a Kerr-Newman black hole
    for a given mass, spin, charge-to-mass ratio and mode.

    :param mass: Mass of the black hole [:math:`M_\odot`].
    :param type: float

    :param spin: Dimensionless spin.
    :param type: float

    :param charge: dimensionless charge-to-mass ratio.
    :param type: float

    :param mode: Angular, azimuthal and overtone number in the form (l,m,n).
    :param type: tuple

:param method: The method used to approximate the Kerr-Newman spectrum. If `fit1`, it uses the fits provided in https://arxiv.org/abs/1307.7315. If `fit2`, it uses the fits provided in the appendix of https://arxiv.org/abs/2109.13961. Deafult: `fit2`.
    """
    l,m,n = mode

    methods = ['fit1']
    if method not in methods:
        raise ValueError('method should be one of '+str(methods))

    if method=='fit1':
        ## use fits as provided by Lorenzo
        # add non-rotating contribution
        coeffs_R = {(2,0):(0.373671,0,-0.0140594,-0.00754829,\
                0.00152926,0.00642125,-0.000992664),\
                  (3,0):(0.599442,0,-0.0544579,-0.0322951,\
                   0.0733055,-0.00175477,-0.0258144)}

        coeffs_I = {(2,0):(-0.0889625,0,-0.00469285,-0.00627901,\
                -0.00195346,-0.00541852,0.00389027),\
                  (3,0):(-0.0927031,0,-0.00715884,-0.011622,\
                  0.0145533,0.00799485,-0.0189346)}

        fs = coeffs_R[(l,n)]
        qs = coeffs_I[(l,n)]
        omegaR = sum([fs[i]*zeta**i for i in range(len(fs))])
        omegaI = sum([qs[i]*zeta**i for i in range(len(qs))])

        # add O(m*spin) contribution
        coeffs_R = {(2,0):(0.0628978,0,-0.0104803,-0.0107957,\
                -0.00112558,-0.00357684,0.0154108),\
                  (3,0):(0.0673656,0,-0.0214809,-0.0220525,\
                   0.0550205,-0.00731002,-0.0344358)}

        coeffs_I = {(2,0):(0.00100132,0,-0.0000142364,-0.00194089,\
                -0.0151062,0.0245121,-0.054223),\
                  (3,0):(0.000651781,0,0.0012462,-0.00442894,\
                  0.00881815,-0.0366831,0.0363787)}

        fs = coeffs_R[(l,n)]
        qs = coeffs_I[(l,n)]
        omegaR += spin*m*sum([fs[i]*zeta**i for i in range(len(fs))])
        omegaI += spin*m*sum([qs[i]*zeta**i for i in range(len(qs))])

        # add O(spin^2) contribution
        coeffs_R = {(2,0):(0.0359141,0,0.0136826,0.0111001,\
                0.0063752,0.0100457,-0.000194589),\
                  (3,0):(0.0475846,0,0.0300244,0.00172412,\
                   0.0823261,-0.234251,0.223573)}

        coeffs_I = {(2,0):(0.00638133,0,0.00810781,0.0108534,\
                -0.00172558,0.0371894,-0.0265663),\
                  (3,0):(0.00659952,0,0.0106413,0.0907605,\
                  -0.451635,1.03177,-0.864954)}

        fs = coeffs_R[(l,n)]
        qs = coeffs_I[(l,n)]
        omegaR += spin**2*sum([fs[i]*zeta**i for i in range(len(fs))])
        omegaI += spin**2*sum([qs[i]*zeta**i for i in range(len(qs))])

        # add O(m^2*spin^2) contribution
        coeffs_R = {(2,0):(0.00895781,0,-0.0083094,-0.0130798,\
                0.00968677,-0.0783303,0.0996714),\
                  (3,0):(0.00661201,0,-0.009497,-0.00998551,\
                   0.0197499,0.00284495,-0.0411144)}

        coeffs_I = {(2,0):(-0.000312294,0,0.00302198,0.00492451,\
                -0.0241024,0.0709671,-0.110537),\
                  (3,0):(0.0000615071,0,-0.000712262,0.0208969,\
                  -0.13611,0.332366,- 0.305122)}

        fs = coeffs_R[(l,n)]
        qs = coeffs_I[(l,n)]
        omegaR += (m*spin)**2*sum([fs[i]*zeta**i for i in range(len(fs))])
        omegaI += (m*spin)**2*sum([qs[i]*zeta**i for i in range(len(qs))])

        # f and tau
        f = omegaR/(2*np.pi)/mass/tsun
        tau = -1/omegaI*mass*tsun

    return f, tau
