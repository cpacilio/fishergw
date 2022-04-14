import numpy as np
from scipy.interpolate import interp1d
import sys

print('------------------------------------------------------------------------------------------------.')
print('Reproduces the first block of Table 1 in Arun et al. (2005), https://arxiv.org/abs/gr-qc/0411146.')
print('NOTE: after the 2PN order, uncertainties in phi_c deviate \n\
        because of a different convention in the definition of the 2.5PN term.')
print('------------------------------------------------------------------------------------------------.')

## import fishergw objects
sys.path.append('..')
from fishergw.taylorf2 import CompactObject, TaylorF2, Fisher

## define the sensitivity curve
def S_h(f):
    S0 = 1e-49
    fs = 20
    f0 = 215
    x = f/f0
    out = S0*(x**(-4.14)-5/x**2+111*(1-x**2+0.5*x**4)/(1+0.5*x**2))
    return out

## define a function to compute sigmas without code repetition
def compute_sigma(signal):
    ## define fisher matrix
    logscale_keys = ['M_c','eta']
    fisher = Fisher(signal,keys=keys,logscale_keys=logscale_keys)
    fisher.psd = S_h
    fmin = 20
    fmax = signal.isco(mode='static')
    fm = fisher.fisher_matrix(fmin,fmax,nbins=1e4)
    ## compute snr
    snr = fisher.snr(fmin,fmax,nbins=1e4)
    ## renormalize to snr=10
    rho0 = 10
    fm *= rho0**2/snr**2
    ## compute uncertainties
    sigma = fisher.sigma1d(fm,svd=False)
    ## scale time in milliseconds and masses in percent
    sigma['t_c'] *= 1000
    sigma['M_c'] *= 100
    sigma['eta'] *= 100
    return sigma

keys = ['t_c','phi_c','M_c','eta']
string = '\nPN order\t'
for k in keys:
    string += '%s\t\t'%k
print(string)

## binary neutron star case
print('\nNS-NS\n')

## define intrinsic parameters
m1, m2 = 1.4, 1.4
chi1, chi2 = 0., 0.
## define binary objects
obj1 = CompactObject(m1,chi1)
obj2 = CompactObject(m2,chi2)
## define signal
for PN in [1,1.5,2,2.5,3,3.5]:
    signal = TaylorF2(obj1,obj2,redshift=False,PN_phase=PN)
    sigma = compute_sigma(signal)

    string = '%s\t\t'%PN
    for k in keys:
        string += '%.2E\t'%sigma[k]
    print(string)

## black hole neutron star case
print('\nNS-BH\n')

## define intrinsic parameters
m1, m2 = 10, 1.4
chi1, chi2 = 0., 0.
## define binary objects
obj1 = CompactObject(m1,chi1)
obj2 = CompactObject(m2,chi2)
## define signal
for PN in [1,1.5,2,2.5,3,3.5]:
    signal = TaylorF2(obj1,obj2,redshift=False,PN_phase=PN)
    sigma = compute_sigma(signal)

    string = '%s\t\t'%PN
    for k in keys:
        string += '%.2E\t'%sigma[k]
    print(string)

## binary black hole case
print('\nBH-BH\n')

## define intrinsic parameters
m1, m2 = 10, 10
chi1, chi2 = 0., 0.
## define binary objects
obj1 = CompactObject(m1,chi1)
obj2 = CompactObject(m2,chi2)
## define signal
for PN in [1,1.5,2,2.5,3,3.5]:
    signal = TaylorF2(obj1,obj2,redshift=False,PN_phase=PN)
    sigma = compute_sigma(signal)

    string = '%s\t\t'%PN
    for k in keys:
        string += '%.2E\t'%sigma[k]
    print(string)

