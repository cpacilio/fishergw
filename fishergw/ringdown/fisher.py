import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simps
#from scipy.linalg import block_diag
from os.path import realpath, dirname
#import warnings

full_path = realpath(__file__)
dir_path = dirname(full_path)

class Fisher():
    """
    An object to load the power spectral density (PSD) and compute the signal-to-noise ratio (SNR) and the Fisher matrix elements of a gravitational wave signal. The SNR and the Fisher matrix are averaged over orientation angles. Recall that average over the inclination angle is already included in the initialization of :class:`RingdownMultimode` .
    
    Attributes:
        **signal** (:class:`RingdownMultimode`) -- A :class:`RingdownMultimode` waveform instance.

        **integration_method** (:class:`scipy.integrate`) -- An integration method from the :class:`scipy.integrate` subpackage.

        **psd** (:class:`scipy.interpolate.interp1d`) -- Interpolant of the PSD. An instance of :class:`scipy.interpolate.interp1d`.

        **fmin** *(float)* -- If ``asd_name`` is provided, ``fmin`` is the corresponding minimum frequency. Otherwise, ``fmin = None``.

        **fmax** *(float)* -- If ``asd_name`` is provided, ``fmax`` is the corresponding maximum frequency. Otherwise, ``fmax = None``.

        **Qavg** *(float)* -- Angle-averaging factor :math:`Q`. When :func:`load_psd` is called, the PSD is divided by :math:`Q^2` to ensure that the SNR and the Fisher matrix elements are angle-averaged.
    
        **_detectors_** *(dict)* -- Dictionary mapping built-in detectors to their ``asd_name`` and angle-averaging factor :math:`Q`. Built-in detectors are Advanced Ligo (``'aligo'``), Cosmic Explorer (``'ce'``), Einstein Telescope (``'etd'``) and LISA (``'lisa'``). The following conventions hold for :math:`Q`:
        
            ``'aLigo'`` and ``'ce'`` are mapped to the factor :math:`Q=1/\sqrt{5}` for a two-armed 90-degrees detector (see, e.g., Eq. (7.177) in [1]);
        
            ``'etd'`` is mapped to :math:`Q=1/\sqrt{5}*\sqrt{3/2}`, the additional factor :math:`\sqrt{3/2}` coming from the fact that ET is a three-armed 60-degrees detector with two channels (see Eq. (4) in https://arxiv.org/abs/1012.0908);
        
            ``'lisa'`` is mapped to :math:`Q=1`, because the LISA sensitivity curve is already averaged over orientation and detector channels (see Eq.s (2,8-9) in https://arxiv.org/abs/1803.01944).

References:

    [1] Maggiore, Michele. Gravitational waves: Volume 1: Theory and experiments. Vol. 1. Oxford university press, 2008.
    """

    _detectors_ = {'aligo':('aligo_asd.txt',1/np.sqrt(5)),\
                 'ce':('ce_asd.txt',1/np.sqrt(5)),\
                 'etd':('etd_asd.txt',1/np.sqrt(3/10)),\
                 'lisa':('lisa_asd.txt',1.)}
    
    def __init__(self,signal,integration_method=simps,\
            asd_name=None,detector=None):
        """
        :param signal: A :class:`RingdownMultimode` waveform instance.
        :type signal: :class:`RingdownMultimode`

        :param integration_method: The integration method to compute the SNR and the Fisher matrix elements. Deafult: :class:`scipy.integrate.simps`.
        :type integration_method: :class:`scipy.integrate` object

        :param asd_name: The path to a text file with the tabulated ASD. If ``None``, psd defaults to 1.0. If ``detector`` is not ``None``, ``asd_name`` is read from the ``_detectors_`` internal dictionary.
        :type asd_name: filepath or ``None``, optional
            
        :param detector: If not ``None``, must be one of ``['aligo','ce','etd,'lisa']``.
        :type detector: str or ``None``, optional

        """
        self.signal = signal
        self.integration_method = simps
        if detector:
            asd_name, Qavg = self._detectors_[detector]
            self.load_psd(asd_name,Qavg)
        elif asd_name:
            self.Qavg = 1.0
            self.load_psd(asd_name)
        else:
            self.Qavg = 1.0
            self.psd = lambda x: 1.0
            self.fmin, self.fmax = None, None
    
    def __repr__(self):
        out = 'Fisher\n'
        out += '\tPSD:\n'
        out += '\t\tname: %s\n'%self.asd_name
        out += '\t\tQavg: '+str(self.Qavg)+'\n'
        try:
            out += '\t\tfmin: %.2f\n'%self.fmin
            out += '\t\tfmax: %.2f\n'%self.fmax
        except:
            out += '\t\tfmin and fmax are None\n'
        return out        
    
    def snr(self,fmin=None,fmax=None,nbins=1e5,logspace=False):
        """
        Returns the SNR of the signal.

        :param fmin: Minimum frequency. If ``None``, defaults to the minimum frequency set by the provided PSD file.
        :type fmin: float or None, optional

        :param fmax: Maximum frequency. If ``None``, defaults to the maximum frequency set by the provided PSD file.
        :type fmax: float or None, optional

        :param nbins: Binning of the integration domain. Default: 1e5.
        :type nbins: float
        
        :param logspace: If ``True``, integration is performed in logspace over frequencies. Default: ``False``.
        :type logspace: bool

        :rtype: float
        """
        if not fmin:
            fmin = self.fmin
        if not fmax:
            fmax = self.fmax
        if not logspace:
            x = np.linspace(fmin,fmax,int(nbins))
        else:
            x = np.logspace(np.log10(fmin),np.log10(fmax),int(nbins))
        out = 0.
        for comp in self.signal.components.values():
            hp, hc = comp(x)
            y = 4*(abs(hp)**2+abs(hc)**2)/self.psd(x)
            out += self.integration_method(y,x)
        out = np.sqrt(out)
        return out

    def snr_BCW(self):
        ## Note: as it is currently implemented,
        ## this approximation is only valid w/o overtones
        ## and it fails when the peak of the lorentian
        ## does not have support in the psd domain
        out = 0.
        for comp in self.signal.components.values():
            Yp, Yc = comp.Yp, comp.Yc
            for i in range(comp.dim):
                freq = comp.__dict__['f%s'%i]
                tau = comp.__dict__['tau%s'%i]
                amp = 10**comp.__dict__['logA%s'%i]
                phi = comp.__dict__['phi%s'%i]
                try:
                    out += 0.5*amp**2*tau/self.psd(freq)*(Yp**2+Yc**2)
                except:
                    out += 1e-20
        out = out**0.5
        return out
    
    def load_psd(self,asd_name,Qavg=1.0,logspace=False):
        """
        Loads the PSD from a text file. Initializes the attributes ``psd``, ``fmin``, ``fmax`` and ``Qavg``.

        :param asd_name: name of the text file containing the ASD table. The file will be first looked into the ``fishergw/detector/`` folder, then into the ``fishergw/detector/unofficial_asd/`` folder. Otherwise, ``asd_name`` must be the full path of the text file.
        :type asd_name: str

        :param Qavg: Value of the angle-averaging factor. See the description in the :class:`Fisher` attributes. Default: 1.
        :type Qavg: float
        
        :param logspace: If ``True``, the ASD is interpolated in logspace along the y-axis. Default: ``False``.
        :type logspace: bool

        """
        try:
            s = np.loadtxt(dir_path+'/../detector/'+asd_name).T
            self.asd_name = 'fishergw/detector/'+asd_name
        except:
            try:
                s = np.loadtxt(dir_path+'/../detector/unofficial_asd/'+asd_name).T
                self.asd_name = 'fishergw/detector/unofficial_asd/'+asd_name
            except:
                s = np.loadtxt(asd_name).T
                self.asd_name = asd_name
        self.Qavg = Qavg
        self.fmin = s[0].min()
        self.fmax = s[0].max()
        if not logspace:
            self.psd = interp1d(s[0],s[1]**2/Qavg**2,kind='linear',\
                                #bounds_error=False,fill_value=1e10,\
                                #fill_value='extrapolate',\
                               )
        else:
            psd = interp1d(s[0],2*np.log10(s[1]),kind='linear',\
                              #bounds_error=False,fill_value=5,\
                              #fill_value='extrapolate',\
                          )
            self.psd = lambda x: 10**psd(x)/Qavg**2
        return None


    def fisher_matrix(self,fmin=None,fmax=None,nbins=1e5,logspace=False):
        """
        Returns the Fisher matrix. Because it is block diagonal w.r.t. (l,m), the output is a dictionary of (l,m) couples and their corresponding blocks.
        
        :param fmin: Minimum frequency. If ``None``, defaults to the minimum frequency set by the provided ASD file.
        :type fmin: float or None, optional

        :param fmax: Maximum frequency. If ``None``, defaults to the maximum frequency set by the provided ASD file.
        :type fmax: float or None, optional

        :param nbins: Binning of the integration domain. Default: 1e5.
        :type nbins: float
        
        :param logspace: If ``True``, integrations are performed in logspace over frequencies. Default: ``False``.
        :type logspace: Bool
        """
        #dim = self.signal.dim
        if not fmin:
            fmin = self.fmin
        if not fmax:
            fmax = self.fmax
        out = {}
        if not logspace:
            f = np.linspace(fmin,fmax,int(nbins))
        else:
            f = np.logspace(np.log10(fmin),np.log10(fmax),int(nbins))
        for comp in self.signal.components.values():
            lm = comp.lm
            ## build block matrix
            loc_Nabla = comp.Nabla()
            loc_dim = len(loc_Nabla)
            loc_fm = np.zeros((loc_dim,loc_dim))
            for i in range(loc_dim):
                for j in range(i,loc_dim):
                    dhp_i, dhc_i = loc_Nabla[i](f)
                    dhp_j, dhc_j = loc_Nabla[j](f)
                    y = 4*np.real(dhp_i*np.conj(dhp_j)+dhc_i*np.conj(dhc_j))/self.psd(f)
                    loc_fm[i,j] = self.integration_method(y,f)
                    loc_fm[j,i] = loc_fm[i,j]
            out[lm] = loc_fm
        return out

    @staticmethod
    def invert_matrix(matrix,svd=False,tol=1e-6):
        """
        Returns the inverse of a matrix.

        :param matrix: The input matrix.
        :type matrix: :class:`numpy.array`

        :param svd: If ``True`` a singular value decomposition is applied and principal directions with singular values smaller than ``tol`` are discarded. Default: ``False``.
        :type svd: bool, optional

        :param tol: Directions with singular values smaller than ``tol`` are discarded when inverting the matrix. Default: 1e-6.
        :type tol: float

        :rtype: :class:`numpy.array` 
        """
        dim = len(matrix)
        if not svd:
            inverse_matrix = np.matrix(matrix).I
        else:
            diag = np.sqrt(matrix[range(dim),range(dim)])
            norm = np.outer(diag,diag)
            matix_normalized = matrix/norm
            u,s,vh = np.linalg.svd(matrix_normalized)
            idx = sum(s/s[0]<tol)
            if idx:
                s[-idx:] = np.inf
            inverse_matrix_normalized = np.dot(vh.T,np.dot(np.diag(1/s),u.T))
            inverse_matrix = inverse_matrix_normalized/norm
        return inverse_matrix

    def covariance_matrix(self,fm,svd=False):
        """
        Returns the covariance matrix for a single (l,m) mode. If ``svd``, the Fisher matrix is inverted with a singular value decomposition.

        :param fm: The \Fisher matrix for a single (l,m) mode.
        """
        return self.invert_matrix(fm,svd)
    
    def correlation_matrix(self,fm,svd=False):
        """
        Returns the correlation matrix for a single (l,m) mode. If ``svd``, the Fisher matrix is inverted with a singular value decomposition.

        :param fm: The \Fisher matrix for a single (l,m) mode.
        """
        inverse_fm = self.invert_matrix(fm,svd)
        corr = np.zeros_like(inverse_fm)
        dim = len(fm)
        for i in range(dim):
            for j in range(i,dim):
                corr[i,j] = inverse_fm[i,j]/np.sqrt(inverse_fm[i,i]*inverse_fm[j,j])
                corr[j,i] = corr[i,j]
        return corr
    
    def sigma1d(self,fm,svd=False):
        """
        Returns the standard deviations of the 1D marginalized posteriors. If ``svd``, the Fisher matrix is inverted with a singular value decomposition.

        :param fm: The \Fisher matrix.

        :returns: Dictionary mapping ``keys`` to the corresponding standard deviations.
        :rtype: dict
        """
        sigma = {}
        for comp in self.signal.components.values():
            loc_fm = fm[comp.lm]
            inverse_loc_fm = self.invert_matrix(loc_fm,svd)
            loc_modes = comp.modes
            for i in range(len(loc_modes)):
                sigma[loc_modes[i]] = [np.sqrt(inverse_loc_fm[j+4*i,j+4*i]) for j in range(4)]
        return sigma

    def sample(self,lm,covariance_matrix,nsamples=1e4,keys=None):
        """
        Returns sampled parameters for a given (l,m) from the corresponding ``covariance matrix``.

        :param covariance_matrix: The covariance matrix for the (l,m) modes.

        :param nsamples: Number of samples to draw. Default: 1e4.
        :type nsamples: float
        """
        comp = self.signal[lm]
        if not keys:
            keys = comp.keys
            cov = covariance_matrix
        else:
            idx = [comp.keys.index(k) for k in keys]
            cov = covariance_matrix[idx][:,idx]
        means = np.array([comp.__dict__[k] for k in keys])
        samples = np.random.multivariate_normal(means,cov,int(nsamples)).T
        #out = {}
        #for i in range(len(keys)):
            #out[keys[i]] = samples[i]
        return samples
