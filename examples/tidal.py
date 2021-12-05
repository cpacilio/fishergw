import numpy as np
from scipy.interpolate import interp1d
import sys

print('\nThis script explores the effects of including the spins and delta_Lamda\n\
in the computation of the Fisher matrix for binary neutron star system.\n\
The effect of low-spin priors is also explored.\n')

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
    #snr = fisher.snr(fmax=fmax,nbins=int(1e4))
    #print('SNR:\t%d'%snr)
    ## compute uncertainties
    sigma = fisher.sigma1d(fm,svd=svd)
    ## scale time in milliseconds and masses in percent
    sigma['t_c'] *= 1000
    sigma['M_c'] *= 100
    sigma['eta'] *= 100
    return sigma

## define intrinsic parameters
m1, m2 = 1.6, 1.4
chi1, chi2 = 0., 0.
l1, l2 = 200, 300
## define binary objects
obj1 = CompactObject(m1,chi1,Lamda=l1)
obj2 = CompactObject(m2,chi2,Lamda=l2)
## define signal
d_L = 50.
signal = TaylorF2(obj1,obj2,redshift=True,d_L=d_L)
print('Lamda_T:\t%.2f'%signal.Lamda_T)

keys = ['t_c','phi_c','M_c','eta','Lamda_T','chi_s','chi_a','delta_Lamda']
string = '\n'
for k in keys:
    string += '%s\t\t'%k
print(string)

print('\ndelta_Lamda: YES, spins: YES')
sigma = compute_sigma(signal)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

print('\ndelta_Lamda: YES, spins: YES, low-spin priors: YES')
priors = {'chi_s':0.05, 'chi_a':0.05}
sigma = compute_sigma(signal,priors=priors)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

print('\ndelta_Lamda: NO, spins: YES')
keys = ['t_c','phi_c','M_c','eta','Lamda_T','chi_s','chi_a']
sigma = compute_sigma(signal)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

print('\ndelta_Lamda: NO, spins: YES, low-spin priors: YES')
keys = ['t_c','phi_c','M_c','eta','Lamda_T','chi_s','chi_a']
priors = {'chi_s':0.05, 'chi_a':0.05}
sigma = compute_sigma(signal,priors=priors)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

print('\ndelta_Lamda: NO, spins: NO')
keys = ['t_c','phi_c','M_c','eta','Lamda_T']
sigma = compute_sigma(signal)
string = ''
for k in keys:
    string += '%.2E\t'%sigma[k]
print(string)

