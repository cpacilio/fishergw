import numpy as np
from scipy.interpolate import interp1d
import sys

print('\nThis script highlights the effects of priors and singular value decomposition\n\
on the computation of the Fisher matrix.\n')

## import fishergw objects
sys.path.append('..')
from fishergw.taylorf2 import CompactObject, TaylorF2, Fisher

## define a function to compute sigmas without code repetition
def compute_sigma(signal,svd=False,priors=None):
    ## define fisher matrix
    log_scale_keys = ['M_c','eta']
    fisher = Fisher(signal,keys=keys,log_scale_keys=log_scale_keys,detector='aLigo')
    fmax = signal.isco(mode='static')
    fm = fisher.fisher_matrix(fmax=fmax,nbins=int(1e4),priors=priors)
    ## compute snr
    snr = fisher.snr(fmax=fmax,nbins=int(1e4))
    print('SNR:\t%d'%snr)
    ## compute uncertainties
    sigma = fisher.sigma1d(fm,svd=svd)
    ## scale time in milliseconds and masses in percent
    sigma['t_c'] *= 1000
    sigma['M_c'] *= 100
    sigma['eta'] *= 100
    return sigma

keys = ['t_c','phi_c','M_c','eta','chi_s','chi_a']
string = '\n'
for k in keys:
    string += '%s\t\t'%k
print(string)

## define intrinsic parameters
m1, m2 = 20, 10
chi1, chi2 = 0., 0.
## define binary objects
obj1 = CompactObject(m1,chi1)
obj2 = CompactObject(m2,chi2)
## define signal
d_L = 50.
signal = TaylorF2(obj1,obj2,redshift=True,d_L=d_L)

priors = {'chi_s':0.05,'chi_a':0.05}

## no svd no priors
print('\nsvd: NO, priors: NO')
sigma = compute_sigma(signal)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

## blach hole neutron star case
print('\nsvd: YES, priors: NO')
sigma = compute_sigma(signal,svd=True)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

## binary black hole star case
print('\nsvd: NO, priors: YES')
sigma = compute_sigma(signal,priors=priors)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

## binary black hole star case
print('\nsvd: YES, priors: YES')
sigma = compute_sigma(signal,priors=priors,svd=True)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

