import numpy as np
import sympy as sp
from sympy import Rational
from copy import deepcopy
from scipy.interpolate import interp1d
from scipy.integrate import simps
from ..cosmology import redshift_from_distance
from ..constants import speed_of_light, solar_mass, G, Mpc
from .utils import spherical_harmonics

class RingdownLM():
    """
    A class to define the ringdown components with the same angular and azimuthal numbers (l,m) and including the fundamental tone (n=0) and possibly overtones (n>0). 
    """
    def __init__(self,modes,freqs,taus,amps,phis,iota=None):
        """
        :param modes: a list specifying the corresponding modes, each as a tuple of angular, azimuthal and overtone numbers (l,m,n).
        :type modes: tuple

        :param freqs: frequencies of the tones in ascending order starting from n=0.
        :type freqs: list

        :param taus: damping times of the tones in ascending order starting from n=0.
        :type taus: list

        :param amps: amplitudes of the tones in ascending order starting from n=0.
        :type taus: list

        :param phis: phases of the tones in ascending order starting from n=0.
        :type taus: list

        :param iota: the inclination angle, which determines the values of the plus and cross spherical harmonics. If ``None``, spherical harmonics are averaged over the inclination angle. Default: ``None``.
        :type iota: float or ``None``

        """
        self.modes = modes
        lm = modes[0][:-1]
        self.lm = lm
        self.dim = len(freqs)
        # compute averages of Yp and Yc over inclination angle
        if iota is not None:
            self.Yp, self.Yc = spherical_harmonics(iota,lm)
        else:
            X = np.linspace(0,np.pi,1000)
            Yp, Yc = np.array([spherical_harmonics(i,lm) for i in X]).T
            self.Yp = (simps(Yp**2*np.sin(X),X)*0.5)**0.5
            self.Yc = (simps(Yc**2*np.sin(X),X)*0.5)**0.5
        # parameters of the component
        self.keys = []
        for i in range(self.dim):
            self.keys += ['f%s'%i,'tau%s'%i,'logA%s'%i,'phi%s'%i]
        # store parameter values
        for i in range(self.dim):
            self.__dict__['f%s'%i] = freqs[i]
            self.__dict__['tau%s'%i] = taus[i]
            self.__dict__['logA%s'%i] = np.log10(amps[i])
            self.__dict__['phi%s'%i] = phis[i]
        # when eval=True, avoid repeating costing evaluations
        self.eval = False

    def __repr__(self):
        out = 'RingdownLM(lm=%s)\n'%str(self.lm)
        for i in range(self.dim):
            out += 'mode%d:\t'%i+str(self.modes[i])+'\n'
            for j in range(4*i,4*i+4):
                k = self.keys[j]
                v = self.__dict__[k]
                out += '\t%s:\t'%k + str(v) + '\n'
        return out

    def __call__(self,f):
        """
        Returns the plus and cross polarizations at a given frequency f, including multiplicative factors from the plus and cross spherical harmonics.
        
        :param f: frequency
        :type f: float

        :rtype: :class:`numpy.array`, shape (2,)
        """
        if not self.eval:
            params = {k:self.__dict__[k] for k in self.keys}
            ## hp
            hp_eval = self.hp().subs(params)
            self.hp_eval = sp.lambdify('f',hp_eval,modules='numpy')
            ## hp
            hc_eval = self.hc().subs(params)
            self.hc_eval = sp.lambdify('f',hc_eval,modules='numpy')
            ## update eval
            self.eval = True
        hp = self.hp_eval(f)*self.Yp
        hc = self.hc_eval(f)*self.Yc
        return hp,hc

    def Nabla(self):
        """
        Returns the gradient of the plus and cross signal polarizations.

        :rtype: list of lambda functions, shape: ``len(self.keys)``
        """
        out = []
        for argument in self.keys:
            x = self.diff_(argument)
            out.append(x)
        return out

    def diff_(self,argument):
        """
        Returns the derivative of the plus and cross signal polarizations w.r.t. argument.

        :param argument: argument of the derivative.
        :type argument: str

        :rtype: tuple of lambda functions, shape (2,)
        """
        keys = deepcopy(self.keys)
        keys.remove(argument)
        args = {k:self.__dict__[k] for k in keys}
        ## hp
        hp_diff = self.hp().subs(args)
        hp_diff = sp.diff(hp_diff,argument)
        hp_diff = hp_diff.subs({argument:self.__dict__[argument]})
        hp_diff = sp.lambdify('f',hp_diff,modules='scipy')
        ## hc
        hc_diff = self.hc().subs(args)
        hc_diff = sp.diff(hc_diff,argument)
        hc_diff = hc_diff.subs({argument:self.__dict__[argument]})
        hc_diff = sp.lambdify('f',hc_diff,modules='scipy')

        return lambda f: (self.Yp*hp_diff(f),self.Yc*hc_diff(f))

    def hp(self):
        """
        Returns the fourier transform of the plus polarization of the signal.
        """
        # fundamental tone
        sv = 'f'
        for k in self.keys[:4]:
            sv += ','+k
        f, f0, tau0, logA0, phi0 = sp.symbols(sv)
        h = (10**logA0)*sp.exp(-sp.I*phi0)*tau0*(sp.exp(2*sp.I*phi0)/(1+4*(f+f0)**2*np.pi**2*tau0**2)+\
                                   1/(1+4*(f-f0)**2*np.pi**2*tau0**2))/np.sqrt(2)
        # overntones
        if self.dim>1:
            for i in range(1,self.dim):
                sv = 'f'
                for k in self.keys[i*4:(i+1)*4]:
                    sv += ','+k
                f, f0, tau0, logA0, phi0 = sp.symbols(sv)
                h += (10**logA0)*sp.exp(-sp.I*phi0)*tau0*(sp.exp(2*sp.I*phi0)/(1+4*(f+\
                        f0)**2*np.pi**2*tau0**2)+1/(1+4*(f-f0)**2*np.pi**2*tau0**2))/np.sqrt(2)
        return h


    def hc(self):
        """
        Returns the fourier transform of the cross polarization of the signal.
        """
        # fundamental tone
        sv = 'f'
        for k in self.keys[:4]:
            sv += ','+k
        f, f0, tau0, logA0, phi0 = sp.symbols(sv)
        h = sp.I*(10**logA0)*sp.exp(-sp.I*phi0)*tau0*(-sp.exp(2*sp.I*phi0)/(1+4*(f+\
                    f0)**2*np.pi**2*tau0**2)+1/(1+4*(f-f0)**2*np.pi**2*tau0**2))/np.sqrt(2)
        # overntones
        if self.dim>1:
            for i in range(1,self.dim):
                sv = 'f'
                for k in self.keys[i*4:(i+1)*4]:
                    sv += ','+k
                f, f0, tau0, logA0, phi0 = sp.symbols(sv)
                h += sp.I*(10**logA0)*sp.exp(-sp.I*phi0)*tau0*(-sp.exp(2*sp.I*phi0)/(1+4*(f+\
                    f0)**2*np.pi**2*tau0**2)+1/(1+4*(f-f0)**2*np.pi**2*tau0**2))/np.sqrt(2)
        return h

class RingdownMultimode():
    """
    A class to define a multimode ringdown signal as a superimposition of damped sinusoids.
    """
    def __init__(self,modes,freqs,taus,amps,phis,iota=None):
        """
        :param freqs: frequencies of each mode.
        :type freqs: list

        :param taus: damping times of each mode.
        :type taus: list

        :param amps: amplitudes of each mode.
        :type taus: list

        :param phis: phases of each mode.
        :type taus: list

        :param iota: inclination angle. If ``None``, polarizations are averaged over the inclination angle. Default: ``None``.
        :type iota: float or None

        :param modes: a list specifying the corresponding modes, each as a set of angular, azimuthal and overtone numbers (l,m,n).
        :type modes: list
        """
        self.components = {}
        ## sort modes by (l,m)
        sorted_modes = sorted(modes, key=lambda tup: (tup[0],tup[1],tup[2]))
        mode = sorted_modes[0]
        tmp = [[mode]]
        l,m = mode[:-1]
        for mode in sorted_modes[1:]:
            if mode[0]==l and mode[1]==m:
                tmp[-1].append(mode)
            else:
                tmp.append([mode])
                l,m = mode[:-1]
        sorted_modes = tmp
        self.sorted_modes = sorted_modes
        ## fill components
        for lm in sorted_modes:
            loc_modes = lm
            indices = [modes.index(lmn) for lmn in loc_modes]
            loc_freqs = [freqs[i] for i in indices]
            loc_taus = [taus[i] for i in indices]
            loc_amps = [amps[i] for i in indices]
            loc_phis = [phis[i] for i in indices]
            loc_mode = lm[0][:-1]
            self.components[loc_mode] = RingdownLM(lm,loc_freqs,loc_taus,loc_amps,loc_phis,iota)
        
        self.dim = len(self.components)

    def __repr__(self):
        out = 'RingdownMultimode\n'
        out += 'sorted modes:\n'
        for lm in self.sorted_modes:
            for lmn in lm:
                out += '\t' + str(lmn)+'\t'
            out += '\n'
        return out

    def __getitem__(self,lm):
        return self.components[lm]

    def __call__(self,f):
        hp, hc = 0., 0.
        for comp in self.components.values():
            loc_hp, loc_hc = comp(f)
            hp += loc_hp
            hc += loc_hc
        return hp, hc
