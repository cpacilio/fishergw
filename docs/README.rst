About
-----
A Python package to compute Fisher matrices for gravitational wave models

Installation
------------
Install from folder::
    
   $ pip install .

Install from pip::

   $ pip install fishergw

Usage of taylorf2
-----------------
    >>> from fishergw.taylorf2 import CompactObject, TaylorF2, Fisher
    >>>
    >>> mass_1, mass_2 = 1.6, 1.4
    >>> luminosity_distance = 100
    >>> spin_1, spin_2 = 0., 0.
    >>> lamda_1, lamda_2 = 200, 350
    >>> obj1 = CompactObject(mass_1,spin_1,Lamda=lamda_1)
    >>> obj2 = CompactObject(mass_2,spin_2,Lamda=lamda_2)
    >>> signal = TaylorF2(obj1,obj2,d_L=luminosity_distance,redshift=False)
    >>>
    >>> keys=['t_c','phi_c','M_c','eta','Lamda_T','chi_s','chi_a']
    >>> logscale_keys = ['M_c','eta']
    >>> fisher = Fisher(signal,detector='etd',\
    >>>         keys=keys,logscale_keys=logscale_keys)
    >>> fmin = 5
    >>> fmax = signal.isco(mode='static')
    >>>
    >>> snr = fisher.snr(fmin,fmax,nbins=1e5)
    >>> priors = {'chi_s':0.05,'chi_a':0.05}
    >>> fisher_matrix = fisher.fisher_matrix(fmin,fmax,nbins=1e5,priors=priors)
    >>> covariance_matrix = fisher.covariance_matrix(fisher_matrix)
    >>> correlation_matrix = fisher.correlation_matrix(fisher_matrix)
    >>> sigmas = fisher.sigma1d(fisher_matrix)
    >>> samples = fisher.sample(covariance_matrix,nsamples=1e5)

Usage of ringdown
-----------------
    >>> from fishergw.ringdown import RingdownMultimode, Fisher
    >>> from fishergw.ringdown.utils import qnm_Kerr, final_mass, final_spin,\
    >>>     ringdown_absolute_amplitudes
    >>> from fishergw.cosmology import redshift_from_distance
    >>>
    >>> mass_1_source, mass_2_source = 36, 30
    >>> spin_1, spin_2 = 0, 0
    >>> mass_f_source = final_mass(mass_1_source,mass_2_source,spin_1,spin_2)
    >>> spin_f = final_spin(mass_1_source,mass_2_source,spin_1,spin_2)
    >>> 
    >>> luminosity_distance = 100
    >>> redshift = redshift_from_distance(luminosity_distance)
    >>> mass_f = mass_f_source*(1+redshift)
    >>>
    >>> modes = [(2,2,0),(3,3,0)]
    >>> freqs, taus, amps = [], [], []
    >>> for mode in modes:
    >>>     f,tau = qnm_Kerr(mass_f,spin_f,mode)
    >>>     freqs.append(f)
    >>>     taus.append(tau)
    >>>     amps.append(ringdown_absolute_amplitude(mass_f,mass_1_source/mass_2_source,\
    >>>         spin_1,spin_2,luminosity_distance)
    >>> phis = [np.random.uniform(0,2*np.pi) for m in modes]
    >>> signal = RingdownMultimode(modes,freqs,taus,amps,phis)
    >>>
    >>> fisher = Fisher(signal,detector='etd')
    >>>
    >>> snr = fisher.snr()
    >>> fisher_matrix = fisher.fisher_matrix(nbins=1e5)
    >>> sigmas = fisher.sigma1d(fisher_matrix)

Usage of cosmology
------------------

    >>> from fishergw.cosmology import redshift_from_distance, distance_from_redshift
    >>> z = redshift_from_distance(100)
    >>> d_L = distance_from_redshift(z)
