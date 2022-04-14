## contains functions to fit Kerr qnms, amplitudes and phases

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import bisect
from ..constants import speed_of_light, solar_mass, G, Mpc
from .harmonics import sYlm
from os.path import realpath, dirname
#import warnings

full_path = realpath(__file__)
dir_path = dirname(full_path)

cc = speed_of_light
msun = solar_mass
## solar mass in secs
tsun = msun*G/cc**3 ## secs

def spherical_harmonics(angle,mode):
    l, m = mode
    Yp = np.real(sYlm(-2,l,m,angle,0)+(-1)**l*sYlm(-2,l,-m,angle,0))
    Yc = np.real(sYlm(-2,l,m,angle,0)-(-1)**l*sYlm(-2,l,-m,angle,0))
    return Yp, Yc

def ringdown_relative_amplitudes(mass_ratio,spin1,spin2,mode,method=1):
    """
    Returns the ratio of the ringdown amplitudes relative to the amplitude of the (2,2,0) mode.
    
    :param mass_ratio: the mass ratio of the progenitor binary, defined as m1/m2>1.
    :type mass_ratio: float

    :param spin1: the projection of the dimensionless spin of the first black hole component the z axis.
    :type spin1: float

    :param spin2: the projection of the dimensionless spin of the second component along the z axis.
    :type spin2: float

    :param mode: Angular, azimuthal and overtnone number in the form (l,m,n).
    :type mode: tuple

    :param method: if 1, use fits from https://arxiv.org/abs/1207.0399; if 2, use fits from https://arxiv.org/abs/2005.03260 in an updated form; if 3, uses EMOP fits from https://arxiv.org/abs/1710.02156. Default: 1.
    """

    if method == 1:
        ## use https://arxiv.org/abs/1207.0399
        eta = mass_ratio/(mass_ratio+1)**2
        chi_p = (mass_ratio*spin1+spin2)/(1+mass_ratio)
        chi_m = (mass_ratio*spin1-spin2)/(1+mass_ratio)
        delta = (1-4*eta)**0.5
        chi_eff = 0.5*(delta*chi_p+chi_m)
        A = {(2,2,0):1.}
        A[(2,1,0)] = 0.43*(delta-chi_eff)
        A[(3,3,0)] = 0.44*(1-4*eta)**0.45
        A[(4,4,0)] = 5.4*(eta-0.22)**2+0.04
    
    elif method == 2:
        ## use https://arxiv.org/abs/2005.03260
        q = mass_ratio
        A = {(2,2,0):1.}
        A[(2,1,0)] = 0.473846-1.22756/q+1.61047/q**2-0.85676/q**3
        A[(3,3,0)] = 0.439698-0.611581/q+0.199865/q**2-0.0279826/q**3
        A[(2,2,1)] = 0.373974+0.074412/q+0.416288/q**2-0.322963/q**3

    elif method == 3:
        ## use https://arxiv.org/abs/1710.02156
        q = mass_ratio
        eta = mass_ratio/(mass_ratio+1)**2
        delta = (1-4*eta)**0.5
        chi_p = (mass_ratio*spin1+spin2)/(1+mass_ratio)
        chi_m = (mass_ratio*spin1-spin2)/(1+mass_ratio)
        def EMOP():
            modes = ((2,2,0),(3,3,0),(2,1,0),(4,4,0))
            Amode0 = {}
            Amode1 = {}
            ## nonspinning coefficients
            coeffs = {(2,2,0):(0.303,0.571,0),\
                    (3,3,0):(0.157,0.671,0),\
                    (2,1,0):(0.099,0.06,0),\
                    (4,4,0):(0.122,-0.188,-0.964)}
            for mode in modes:
                a,b,c = coeffs[mode]
                Amode0[mode] = a + b*eta + c*eta**2
            ## spinning coefficients
            coeffs = {(2,2,0):(-0.07,0.255,0.189,-0.013,0.084),\
                    (3,3,0):(0.163,-0.187,0.021,0.073,0),\
                    (2,1,0):(-0.067,0,0,0,0),\
                    (4,4,0):(-0.207,0.034,-0.701,1.387,0.122)}
            a,b,c,d,e = coeffs[(2,2,0)]
            Amode1[(2,2,0)] = eta*chi_p*(a + b/q + c*q + d*q**2) + e*delta*chi_m
            a,b,c,d,e = coeffs[(3,3,0)]
            Amode1[(3,3,0)] = eta*chi_m*(a + b/q + c*q) + d*delta*chi_p
            a,b,c,d,e = coeffs[(2,1,0)]
            Amode1[(2,1,0)] = a*chi_m
            a,b,c,d,e = coeffs[(4,4,0)]
            Amode1[(4,4,0)] = eta*chi_p*(a/q + b*q) + delta*eta*chi_m*(c + d/q + e*q)
            ## combine
            out = {}
            for mode in modes:
                if mode[1]%2:
                    out[mode] = eta**2*(Amode0[mode]*delta + Amode1[mode])**2
                else:
                    out[mode] = eta**2*(Amode0[mode] + Amode1[mode])**2
            return out
        emops = EMOP()
        A = {}
        for k in emops.keys():
            A[k] = (emops[k]/emops[(2,2,0)])**0.5
        spin_f = final_spin(mass_ratio,1,spin1,spin2)
        f_22, tau_22 = qnm_Kerr(1,spin_f,(2,2,0),method='interp')
        w_22 = tau_22/(1+4*np.pi**2*f_22**2+tau_22**2)
        for k in emops.keys():
            f_mode, tau_mode = qnm_Kerr(1,spin_f,k,method='interp')
            w_mode = tau_mode/(1+4*np.pi**2*f_mode**2+tau_mode**2)
            A[k] *= (w_mode/w_22)**0.5

    return A[mode]

def ringdown_absolute_amplitudes(mass_f,mass_ratio,spin1,spin2,d_L,mode,method=1):
    eta = mass_ratio/(mass_ratio+1)**2
    A = {(2,2,0):0.864*eta*mass_f*tsun*cc/Mpc/d_L}
    out = ringdown_relative_amplitudes(mass_ratio,spin1,spin2,mode,method)*A[(2,2,0)]
    return out

def qnm_Kerr(mass,spin,mode,method='interp',interp_qnm={}):
    """
    Returns the frequency and the damping time of a Kerr black hole
    for a given mass, spin and mode.

    :param mass: Mass of the black hole [:math:`M_\odot`].
    :type mass: float

    :param spin: Dimensionless spin.
    :type spin: float

    :param mode: Angular, azimuthal and overtone number in the form (l,m,n).
    :type mode: tuple

    :param method: The method used to approximate the Kerr spectrum. If `fit1`, it uses the fits provided in https://arxiv.org/abs/gr-qc/0512160. If `fit2`, it uses the fits provided in the appendix of https://arxiv.org/abs/2109.13961. If `interp`, it interpolates from the numerical tables provided in https://pages.jh.edu/eberti2/ringdown/. Deafault: `interp`.
    """
   
    methods = ['fit1','fit2','interp']

    if method not in methods:
        raise ValueError('method should be one of '+str(methods))

    if method=='fit1':
        ## use https://arxiv.org/abs/gr-qc/0512160
        coeffs = {(2,2,0):(1.5251,-1.1568,0.1292,0.7,1.4187,-0.4990),\
                  (2,2,1):(1.3673,-1.0260,0.1628,0.1,0.5436,-0.4731),\
                  (2,1,0):(0.6,-0.2339,0.4175,-0.3,2.3561,-0.2277),\
                  (3,3,0):(1.8956,-1.3043,0.1818,0.9,2.3430,-0.481),\
                  (4,4,0):(2.3,-1.5056,0.2244,1.1929,3.1191,-0.4825)}
        f1, f2, f3, q1, q2, q3 = coeffs[mode]
        omega = (f1+f2*(1-spin)**f3)/mass/tsun
        f = omega/(2*np.pi)
        Q = q1+q2*(1-spin)**q3
        tau = Q/(np.pi*f)

    elif method=='fit2':
        ## use https://arxiv.org/abs/2109.13961
        BR = {(2,2,0): 0.37367168*np.array([1,-1.899567,1.015454,-0.111430]),\
            (2,2,1): 0.34671099*np.array([1,-1.850299,0.944088,-0.088458]),\
            (3,3,0): 0.59944329*np.array([1,-1.928277,1.044039,-0.112303])}
        ##
        CR = {(2,2,0): np.array([1,-2.238461,1.581677,-0.341455]),\
            (2,2,1): np.array([1,-2.250169,1.611393,-0.359285]),\
            (3,3,0): np.array([1,-2.265230,1.624332,-0.357651])}
        ##
        BI = {(2,2,0): 0.08896232*np.array([1,-2.533958,2.102750,-0.568636]),\
            (2,2,1): 0.27391488*np.array([1,0.366066,-3.290350,1.927196]),\
            (3,3,0): 0.09270305*np.array([1,-2.163575,1.405496,-0.241561])}
        ##
        CI = {(2,2,0): np.array([1,-2.498341,2.056918,-0.557557]),\
            (2,2,1): np.array([1,0.388928,-3.119527,1.746957]),\
            (3,3,0): np.array([1,-2.130446,1.394144,-0.261229])}
        ##
        ss = np.array([spin**i for i in range(4)])
        omegaR = np.sum(BR[mode]*ss)/np.sum(CR[mode]*ss)
        omegaI = np.sum(BI[mode]*ss)/np.sum(CI[mode]*ss)
        # f and tau
        f = omegaR/(2*np.pi)/mass/tsun
        tau = 1/omegaI*mass*tsun

    elif method=='interp':
        ## use numerical tables provided in https://pages.jh.edu/eberti2/ringdown/
        if mode in interp_qnm.keys():
            omega_interp, tau_interp_m1 = interp_qnm[mode]
        else:
            filename = dir_path+'/qnm_data/%s%s%s.dat'%(mode[0],mode[1],mode[2])
            x,y,z,_,_ = np.loadtxt(filename).T
            omega_interp = interp1d(x,y)#,kind='cubic')
            tau_interp_m1 = interp1d(x,-z)#,kind='cubic')
            interp_qnm[mode] = (omega_interp,tau_interp_m1)

        omega = omega_interp(spin)/mass/tsun
        f = omega/(2*np.pi)
        tau = mass*tsun/tau_interp_m1(spin)
    
    return f, tau

def invert_qnm_from_f(f,mode,method='fit2'):
    """
    Returns the masses (in units if solar masses) and dimensionless spins compatible with a given QNM frequency, assuming it is generated from the Kerr spectrum.
    
    :param f: freqency [Hz].
    :type f: float

    :param mode: the QNM numbers (l,m,n).
    :type mode: tuple
    """
    #spins = np.linspace(-0.998,0.998,200)
    spins = np.linspace(0,0.998,100)
    qnms = np.array([qnm_Kerr(1,spin,mode,method)[0] for spin in spins])
    masses = qnms/f
    return masses,spins

def invert_qnm_from_tau(tau,mode,method='fit2'):
    """
    Returns the masses (in units of solar masses) and spins compatible with a given QNM damping time, assuming it is generated from the Kerr spectrum.
    
    :param tau: damping time [s].
    :type tau: float

    :param mode: the QNM numbers (l,m,n).
    :type mode: tuple
    """
    #spins = np.linspace(-0.998,0.998,200)
    spins = np.linspace(0,0.998,100)
    qnms = np.array([qnm_Kerr(1,spin,mode,method)[1] for spin in spins])
    masses = tau/qnms
    return masses,spins


def invert_qnm_from_ff(f1,f2,mode1,mode2,method='fit2',qnm_invert={}):
    """
    Returns the mass (in units of solar masses) and the dimensionless spin corresponding to a pair of frequencies, assuming they are generated from the Kerr spectrum.

    :param f1: The first frequency [Hz].
    :type f1: float

    :param f2: The second frequency [Hz].
    :type f2: float

    :param mode1: The set (l,m,n) of angular and overtone numbers of the first frequency.
    :type mode1: tuple

    :param mode2: The set (l,m,n) of angular oand overtone numbers of the second frequency.
    :type mode2: tuple

    :param method: The method used to approximate the Kerr spectrum. See methods available in the \qnm_Kerr function.

    """
    if mode1 == mode2:
        raise ValueError('mode1 and mode2 cannot be the same')

    def f_ratio(spins):
        f1 = np.array([qnm_Kerr(1,spin,mode1,method)[0] for spin in spins])
        f2 = np.array([qnm_Kerr(1,spin,mode2,method)[0] for spin in spins])
        #f1 = qnm_Kerr(1,spin,mode1,method)[0]
        #f2 = qnm_Kerr(1,spin,mode2,method)[0]
        return f1/f2
    
    if (mode1,mode2) not in qnm_invert.keys():
        #X = np.linspace(-0.998,0.998,200)
        X = np.linspace(0.,0.998,100)
        y = f_ratio(X)
        qnm_invert[(mode1,mode2)] = interp1d(y,X)#,kind='cubic')
  
    try:
        #spin = bisect(lambda x: f_ratio(x)-f1/f2,0.,0.998)
        spin = qnm_invert[(mode1,mode2)](f1/f2)
        spin = spin.item()
        mass = qnm_Kerr(1,spin,mode1,method)[0]/f1
    except:
        spin = 10
        mass = 0
    return mass, spin

def invert_qnm_from_ftau(f1,tau2,mode1,mode2,method='fit2',qnm_invert={}):
    """
    Returns the mass (in units of solar masses) and the dimensionless spin corresponding to a frequency and a damping time, assuming they are generated from the Kerr spectrum.

    :param f1: The frequency [Hz].
    :type f1: float

    :param f2: The damping time [s].
    :type f2: float

    :param mode1: The set (l,m,n) of angular and overtone numbers of the frequency.
    :type mode1: tuple

    :param mode2: The set (l,m,n) of angular oand overtone numbers of the damping time.
    :type mode2: tuple

    :param method: The method used to approximate the Kerr spectrum. See methods available in the \qnm_Kerr function.
    """
    def ftau_product(spins):
        f1 = np.array([qnm_Kerr(1,spin,mode1,method)[0] for spin in spins])
        tau2 = np.array([qnm_Kerr(1,spin,mode2,method)[1] for spin in spins])
        return f1*tau2

    if (mode1,mode2) not in qnm_invert.keys():
        #X = np.linspace(-0.998,0.998,200)
        X = np.linspace(0.,0.998,100)
        y = ftau_product(X)
        qnm_invert[(mode1,mode2)] = interp1d(y,X)#,kind='cubic')

    try:
        spin = qnm_invert[(mode1,mode2)](f1*tau2)
        mass = qnm_Kerr(1,spin,mode1,method)[0]/f1
    except:
        spin = 10
        mass = 0
    return mass, spin

def invert_qnm_from_tautau(tau1,tau2,mode1,mode2,method='fit2',qnm_invert={}):
    """
    Returns the mass (in units of solar masses) and the dimensionless spin corresponding to a pair of damping times, assuming they are generated from the Kerr spectrum.

    :param f1: The first damping time.
    :type f1: float

    :param f2: The second damping time.
    :type f2: float

    :param mode1: The set (l,m,n) of angular and overtone numbers of the first damping time.
    :type mode1: tuple

    :param mode2: The set (l,m,n) of angular oand overtone numbers of the second damping time.
    :type mode2: typle

    :param method: The method used to approximate the Kerr spectrum. See methods available in the \qnm_Kerr function.

    .. note::

    The use of ``fit1`` methods can be badly posed for certain combination of modes, because the ratio tau1/tau2 is not monotonic. This is just an artifact of the ``fit1`` method and it does not happen with other methods.

    """
    if mode1 == mode2:
        raise ValueError('mode1 and mode2 cannot be the same')

    def tautau_ratio(spins):
        tau1 = np.array([qnm_Kerr(1,spin,mode1,method)[1] for spin in spins])
        tau2 = np.array([qnm_Kerr(1,spin,mode2,method)[1] for spin in spins])
        return tau1/tau2

    if (mode1,mode2) not in qnm_invert.keys():
        #X = np.linspace(-0.998,0.998,200)
        X = np.linspace(0.,0.998,100)
        y = tautau_ratio(X)
        qnm_invert[(mode1,mode2)] = interp1d(y,X)#,kind='cubic')

    try:
        #spin = bisect(lambda x: tautau_ratio(x)-tau1/tau2,0.,0.998)
        spin = qnm_invert[(mode1,mode2)](tau1/tau2)
        mass = tau1/qnm_Kerr(1,spin,mode1,method)[1]
    except:
        spin = 10
        mass = 0
    return mass, spin

def r_isco(spin):
    """
    Returns the dimensionless ISCO radius for a given dimensionless spin.

    :param spin: Dimensionless spin.
    :type spin: float
    """
    Z1 = 1+(1-spin**2)**(1/3)*((1+spin)**(1/3)+(1-spin)**(1/3))
    Z2 = (3*spin**2+Z1**2)**0.5
    out = 3+Z2-np.sign(spin)*((3-Z1)*(3+Z1+2*Z2))**0.5
    return out

def E_isco(spin):
    """
    Returns the dimensionless energy for a test particle at the ISCO of a Kerr black hole.

    :param spin: Dimensionless spin.
    :type spin: float
    """
    return (1-2/3/r_isco(spin))**0.5

def L_isco(spin):
    """
    Returns the dimensionless angular momentum for a test particle at the ISCO of a Kerr black hole.

    :param spin: Dimensionless spin.
    :type spin: float
    """
    return 2/3**1.5*(1+2*(3*r_isco(spin)-2)**0.5)


def final_mass(mass1,mass2,spin1,spin2,beta=0.,gamma=0.,method=1):
    """
    Returns the final mass (in units of solar masses) of the Kerr black hole remnant of a quasi-circular binary black hole merger.

    :param mass1: Mass of the primary component [:math:`M_\odot`].
    :type mass1: float

    :param mass2: Mass of the secondary component [:math:`M_\odot`].
    :type mass2: float

    :param spin1: Dimensionless spin of the primary component.
    :type spin1: float

    :param spin2: Dimensionless spin of the secondary component.
    :type spin2: float
    
    :param beta: Angle between spin1 and the z direction. Default: 0.
    :type beta: float
    
    :param gamma: Angle between spin2 and the z direction. Default: 0.
    :type gamma: float

    :param method: If 1, uses the fit in https://arxiv.org/abs/1206.3803. If 2, uses the fit in https://arxiv.org/abs/1508.07250. If 3, uses the fit in https://arxiv.org/abs/0807.2985. Default: 1.
    """
    if mass2>mass1:
        raise ValueError('mass2 must be no more than mass1!')
    q = mass2/mass1
    eta = q/(1+q)**2
    m_tot = mass1+mass2

    if method==1:
        ## use https://arxiv.org/abs/1206.3803
        a_tot = (spin1*np.cos(beta)+spin2*np.cos(gamma)*q**2)/(1+q)**2
        p0, p1 = 0.04827, 0.01707
        E_rad = (1-E_isco(a_tot))*eta+\
               4*eta**2*(4*p0+16*p1*a_tot*(a_tot+1)+\
                E_isco(a_tot)-1)
        m_rad = E_rad*m_tot
        m_fin = m_tot-m_rad

    elif method==2:
        ## use https://arxiv.org/abs/1508.07250
        a_tot = (spin1*np.cos(beta)+spin2*np.cos(beta)*q**2)/(1+q)**2/(1-2*eta)
        E_rad = 0.0559745*eta+0.580951*eta**2-0.960673*eta**3+3.35241*eta**4
        E_rad *= (1+a_tot*(-0.00303023-2.00661*eta+7.70506*eta**2))/(1+\
                a_tot*(-0.67144-1.47569*eta+7.30468*eta**2))
        m_rad = E_rad*m_tot
        m_fin = m_tot-m_rad

    elif method==3:
        ## use https://arxiv.org/abs/0807.2985
        m0 = 0.9515
        m1 = -0.013
        s0 = 0.686
        s1 = 0.15
        E_fin = 1+(m0-1)*4*eta+m1*16*eta**2*(spin1*np.cos(beta)+\
                                             spin2*np.cos(gamma))
        m_fin = m_tot*E_fin

    return m_fin


def final_spin(mass1,mass2,spin1,spin2,beta=0.,gamma=0.,method=2):
    """
    Returns the dimensionless final spin of the Kerr black hole remnant of a quasi-circular binary black hole merger.

    :param mass1: Mass of the primary component [:math:`M_\odot`].
    :type mass1: float

    :param mass2: Mass of the secondary component [:math:`M_\odot`].
    :type mass2: float

    :param spin1: Dimensionless spin of the primary component.
    :type spin1: float

    :param spin2: Dimensionless spin of the secondary component.
    :type spin2: float
    
    :param beta: Angle between spin1 and the z direction. Default: 0.
    :type beta: float
    
    :param gamma: Angle between spin2 and the z direction. Default: 0.
    :type gamma: float

    :param method: If 1, uses the fit in https://arxiv.org/abs/0904.2577. If 2, uses the fit in https://arxiv.org/abs/1605.01938. If 3, uses the fit in https://arxiv.org/abs/1508.07250. If 4, uses the fit in https://arxiv.org/abs/0807.2985. Default: 2.
    """
    if mass2>mass1:
        raise ValueError('mass2 must be no more than mass1!')
    q = mass2/mass1
    eta = q/(1+q)**2

    if method==1:
        ## reproduces https://arxiv.org/abs/0904.2577
        ## deprecated, do not use it
        a_tot = (spin1+spin2*q**2)/(1+q**2)
        kappas = [[2*3**0.5,-2.8904,-0.1229],\
                [-3.5171,0.4537],\
                [2.5763]]
        a_fin = a_tot
        for i in range(len(kappas)):
            for j in range(len(kappas[i])):
                a_fin += kappas[i][j]*eta**(1+i)*a_tot**j

    elif method==2:
        ## reproduces https://arxiv.org/abs/1605.01938
        a_tot = (spin1*np.cos(beta)+spin2*np.cos(gamma)*q**2)/(1+q)**2
        xi = 0.41616
        a_eff = a_tot+xi*eta*(spin1*np.cos(beta)+spin2*np.cos(gamma))
        kappas = [[-3.82,-1.2019,-1.20764],\
                  [3.79245,1.18385,4.90494]]
        a_fin = a_tot + eta*(L_isco(a_eff)-2*a_tot*(E_isco(a_eff)-1))
        for i in range(len(kappas)):
            for j in range(len(kappas[i])):
                a_fin += kappas[i][j]*eta**(2+i)*a_eff**j

    elif method==3:
        ## reproduces https://arxiv.org/abs/1508.07250
        a_tot = (spin1*np.cos(beta)+spin2*np.cos(gamma)*q**2)/(1+q)**2
        kappas = [[2*3**0.5,-0.085,0.102,-1.355,-0.868],\
                  [-4.399,-5.837,-2.097,4.109,2.064],\
                  [9.397],[-13.181]]
        a_fin = a_tot
        for i in range(len(kappas)):
            for j in range(len(kappas[i])):
                a_fin += kappas[i][j]*eta**(1+i)*a_tot**j

    elif method==4:
        ## reproduces https://arxiv.org/abs/0807.2985
        m0 = 0.9515
        m1 = -0.013
        s0 = 0.686
        s1 = 0.15
        c = s1**0.5/(1-s1**0.5)
        g1 = (1+c)**2/(q+c)**2
        g2 = (1+c)**2/(1/q+c)**2
        w = 1.26
        k = 0.008
        a_fin = s0*(4*w*eta+16*(1-w)*eta**2)+s1*(g1*spin1*np.cos(beta)+\
                g2*spin2*np.cos(gamma))+16*k*eta**2*(spin1*np.cos(beta)+\
                                                spin2*np.cos(gamma))**2

    return a_fin

