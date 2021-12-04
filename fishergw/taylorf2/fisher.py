import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simps
from os.path import realpath, dirname
#import warnings

full_path = realpath(__file__)
dir_path = dirname(full_path)

class Fisher():
    """
    An object to load the power spectral density (PSD) and compute the signal-to-noise ratio (SNR) and the Fisher matrix elements of a gravitational wave signal. The SNR and the Fisher matrix are averaged over orientation and inclination angles.
    
    Attributes:
        **signal** (:class:`TaylorF2`) -- A :class:`TaylorF2` waveform instance.

        **integration_method** (:class:`scipy.integrate`) -- An integration method from the :class:`scipy.integrate` module.

        **psd** (:class:`scipy.interpolate.interp1d`) -- Interpolant of the PSD. An instance of :class:`scipy.interpolate.interp1d`.

        **fmin** *(float)* -- If ``psd_name`` is provided, ``fmin`` is the corresponding minimum frequency. Otherwise, ``fmin = None``.

        **fmax** *(float)* -- If ``psd_name`` is provided, ``fmax`` is the corresponding maximum frequency. Otherwise, ``fmax = None``.

        **Qavg** *(float)* -- Angle-averaging factor :math:`Q`. When :func:`load_psd` is called, the PSD is divided by :math:`Q^2` to ensure that the SNR and the Fisher matrix elements are angle-averaged.
    
        **keys** *(list)* -- Independent variables w.r.t. which the Fisher matrix is evaluated.

        **_detectors_** *(dict)* -- Dictionary mapping built-in detectors to their ``psd_name`` and angle-averaging factor :math:`Q`. Built-in detectors are Advanced Ligo ('aLigo'), Cosmic Explorer ('CE'), Einstein Telescope ('ET') and LISA ('lisa'). The following conventions hold for :math:`Q`:
        
            'aLigo' and 'CE' are mapped to the factor :math:`Q=2/5` for a two-armed 90-degrees detector (see, e.g., Eq. (7.177) in [1]);
        
            'ET' is mapped to :math:`Q=2/5\sqrt{3/2}`, the additional factor :math:`\sqrt{3/2}` coming from the fact that ET is a three-armed 60-degrees detector with two channels (see Eq. (4) in https://arxiv.org/abs/1012.0908);
        
            'lisa' is mapped to :math:`Q=2/\sqrt{5}`. This only accounts for averaging over the inclination angle, because the LISA sensitivity curve is already averaged over orientation and detector channels (see Eq.s (2,8-9) in https://arxiv.org/abs/1803.01944).

References:

    [1] Maggiore, Michele. Gravitational waves: Volume 1: Theory and experiments. Vol. 1. Oxford university press, 2008.
    """

    _detectors_ = {'aLigo':(dir_path+'/../detector/aligo_psd.dat',2/5),\
                 'CE':(dir_path+'/../detector/ce_psd.dat',2/5),\
                 'ET':(dir_path+'/../detector/etd_psd.dat',2/5*np.sqrt(3/2)),\
                 'lisa':(dir_path+'/../detector/lisa_psd.dat',2/5*np.sqrt(5))}
    
    def __init__(self,signal,integration_method=simps,\
            psd_name=None,detector=None,keys=None,log_scale_keys=[]):
        """
        :param signal: A :class:`TaylorF2` waveform instance
        :type signal: :class:`TaylorF2`

        :param integration_method: The integration method to compute the SNR and the Gisher matrix elements.
        :type integration_method: :class:`scipy.integrate` object, default= :class:`scipy.integrate.simps`

        :param psd_name: The path to a text file with the tabulated PSD. If ``None``, psd defaults to 1.0. If ``detector`` is not ``None``, ``psd_name`` is read from the ``_detectors_`` dictionary
        :type psd_name: filepath or None, optional
            
        :param detector: If not ``None``, must be one of ``['aLigo','CE','ET','lisa']``.
        :type detector: str or None, optional

        :param keys: Independent variables w.r.t. which the Fisher matrix is evaluated. If ``None``, defaults to ``self.signal.keys``.
        :type keys: list of str or None
        
        :param log_scale_keys: Subset of ``keys`` to be converted in log scale when computing the Fisher. Defaults to the empty list.
        :type log_scale_keys: list
        """
        self.signal = signal
        if not keys:
            self.keys = self.signal.keys
        else:
            self.keys = keys
        self.log_scale_keys = log_scale_keys
        self.integration_method = simps
        if detector:
            psd_name, Qavg = self._detectors_[detector]
            self.load_psd(psd_name,Qavg)
        elif psd_name:
            self.Qavg = 1.0
            self.load_psd(psd_name)
        else:
            self.Qavg = 1.0
            self.psd = lambda x: 1.0
            self.fmin, self.fmax = None, None
        
    def snr(self,fmin=None,fmax=None,nbins=int(1e5)):
        """
        Returns the SNR of the signal.

        :param fmin: Minimum frequency. If ``None``, defaults to the minimum frequency set by the provided PSD file.
        :type fmin: float or None, optional

        :param fmax: Maximum frequency. If ``None``, defaults to the maximum frequency set by the provided PSD file.
        :type fmax: float or None, optional

        :param nbins: Binning of the integration domain.
        :type nbins: int, default=1e5

        :rtype: float
        """
        if not fmin:
            fmin = self.fmin
        if not fmax:
            fmax = self.fmax
        x = np.linspace(fmin,fmax,int(nbins))
        y = 4*np.abs(self.signal(x))**2/self.psd(x)
        snr = self.integration_method(y,x)
        snr = np.sqrt(snr)
        return snr
    
    def load_psd(self,psd_name,Qavg=1.0):
        """
        Loads the PSD from a text file. Initializes the attributes ``psd``, ``fmin``, ``fmax`` and ``Qavg``.

        :param psd_name: Path to the PSD text file.
        :type psd_name: filepath

        :param Qavg: Value of the angle-averaging factor. See the description in the :class:`Fisher` attributes.
        :type Qavg: float, default=1.0

        :returns: self instance
        :rtype: :class:`Fisher`
        """
        s = np.genfromtxt(psd_name).T
        self.Qavg = Qavg
        self.fmin = s[0].min()
        self.fmax = s[0].max()
        self.psd = interp1d(s[0],s[1]/Qavg**2,kind='cubic')
        return self
    
    def fisher_matrix(self,fmin=None,fmax=None,nbins=int(1e5),priors=None):
        """
        Returns the Fisher matrix. Gaussian priors can be specified based on https://arxiv.org/abs/gr-qc/9502040.
        
        :param fmin: Minimum frequency. If ``None``, defaults to the minimum frequency set by the provided PSD file.
        :type fmax: float or None, optional

        :param fmax: Maximum frequency. If ``None``, defaults to the maximum frequency set by the provided PSD file.
        :type fmax: float or None, optional

        :param nbins: Binning of the integration domain.
        :type nbins: int, default=1e5
        
        :param priors: Dictionary mapping a key with the standard deviation of its Gaussian prior. Defaults to None.
        :type priors: dict or None

        :rtype: :class:`numpy.array`, shape ``(len(keys),len(keys))``
        """
        Nabla = self.signal._evaluate_Nabla_(keys=self.keys,log_scale_keys=self.log_scale_keys)
        dim = len(self.keys)
        if not fmin:
            fmin = self.fmin
        if not fmax:
            fmax = self.fmax
        fm = np.zeros((dim,dim))
        derivatives = list(Nabla.values())
        f = np.linspace(fmin,fmax,int(nbins))
        for i in range(dim):
            for j in range(i,dim):
                y = 4*np.real(derivatives[i](f)*np.conj(derivatives[j](f)))/self.psd(f)
                fm[i,j] = self.integration_method(y,f)
                fm[j,i] = fm[i,j]
        if priors:
            for k,v in priors.items():
                if k in self.keys:
                    i = self.keys.index(k)
                    fm[i,i] += 1/v**2
                else:
                    warning_message = 'Prior specification on %s is invalid and will be ignored.\nValid keys :'%(k)+str(self.keys)+'\n'
                    print(warning_message)
                    #warnings.warn(warning_message)
        return fm
    
    def invert_matrix(self,fm,tol=1e-6):
        """
        Returns the inverse of the \Fisher matrix. A singular value decomposition is applied and directions with singular values smaller than ``tol`` are discarded.

        :param fm: The \Fisher matrix.
        :type fm: :class:`numpy.array`, shape ``(len(keys),len(keys))``

        :param tol: Directions with singular values smaller than ``tol`` are discarded when inverting the matrix. Default: 1e-6.
        :type tol: float

        :rtype: :class:`numpy.array`, shape ``(len(keys),len(keys))``        
        """
        dim = len(fm)
        diag = np.sqrt(fm[range(dim),range(dim)])
        norm = np.outer(diag,diag)
        fm_normalized = fm/norm
        u,s,vh = np.linalg.svd(fm_normalized)
        idx = sum(s<tol)
        if idx:
            s[-idx:] = np.inf
        inverse_fm_normalized = np.dot(vh.T,np.dot(np.diag(1/s),u.T))
        inverse_fm = inverse_fm_normalized/norm
        return inverse_fm
    
    def covariance_matrix(self,fm,svd=True):
        """
        Returns the covariance matrix. If `svd`, the Fisher matrix is inverted with a singular value decomposition.

        :param fm: The \Fisher matrix.
        :type fm: :class:`numpy.array`, shape ``(len(keys),len(keys))``

        :rtype: :class:`numpy.array`, shape ``(len(keys),len(keys))``
        """
        if svd:
            inverse_fm = self.invert_matrix(fm)
        else:
            inverse_fm = np.matrix(fm).I
        return inverse_fm
    
    def correlation_matrix(self,fm,svd=True):
        """
        Returns the correlation matrix. If `svd`, the Fisher matrix is inverted with a singular value decomposition.

        :param fm: The \Fisher matrix.
        :type fm: :class:`numpy.array`, shape ``(len(keys),len(keys))``

        :rtype: :class:`numpy.array`, shape ``(len(keys),len(keys))``
        """
        if svd:
            inverse_fm = self.invert_matrix(fm)
        else:
            inverse_fm = np.matrix(fm).I
        corr = np.zeros_like(inverse_fm)
        dim = len(fm)
        for i in range(dim):
            for j in range(i,dim):
                corr[i,j] = inverse_fm[i,j]/np.sqrt(inverse_fm[i,i]*inverse_fm[j,j])
                corr[j,i] = corr[i,j]
        #sigma = {self.keys[i]:np.sqrt(inverse_fm[i,i]) for i in range(dim)}
        return corr
    
    def sigma1d(self,fm,svd=True):
        """
        Returns the standard deviations of the 1D marginalized posteriors. If `svd`, the Fisher matrix is inverted with a singular value decomposition.


        :param fm: The \Fisher matrix.
        :type fm: :class:`numpy.array`, shape ``(len(keys),len(keys))``

        :returns: Dictionary mapping ``keys`` to the corresponding standard deviations.
        :rtype: dict
        """
        if svd:
            inverse_fm = self.invert_matrix(fm)
        else:
            inverse_fm = np.matrix(fm).I
        dim = len(fm)
        sigma = {self.keys[i]:np.sqrt(inverse_fm[i,i]) for i in range(dim)}
        return sigma

