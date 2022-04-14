import sys
import numpy as np
## import fishergw objects
sys.path.append('..')
from fishergw.ringdown import RingdownMultimode, Fisher
from fishergw.ringdown.utils import final_mass, final_spin, qnm_Kerr, ringdown_absolute_amplitudes
from fishergw.cosmology import redshift_from_distance

print('------------------------------------------------------------------')
print('Compute the 1-D uncertainties on the parameters of a BH ringdown')
print('using Fisher matrix formalism on a template of damped sinusoids.\n')
print('------------------------------------------------------------------')
print('')

mass_1 = 36
mass_2 = 30
mass_ratio = mass_1/mass_2
spin_1 = 0
spin_2 = 0
luminosity_distance = 100

mass_f_source = final_mass(mass_1,mass_2,spin_1,spin_2)
redshift = redshift_from_distance(luminosity_distance)
mass_f = mass_f_source*(1+redshift)
spin_f = final_spin(mass_1,mass_2,spin_1,spin_2)

print('final mass: %.2f'%mass_f)
print('final spin: %.2f'%spin_f)
print('')

## define signal params
modes = [(2,2,0),(3,3,0)]
# uncomment to use overtone
#modes = [(2,2,0),(2,2,1)]
freqs = []
taus = []
amps = []
phis = [0 for m in modes]
for mode in modes:
    f, tau = qnm_Kerr(mass_f,spin_f,mode)
    freqs.append(f)
    taus.append(tau)
    amp = ringdown_absolute_amplitudes(mass_f,mass_ratio,spin_1,spin_2,\
            luminosity_distance,mode,method=2)
    amps.append(amp)

print('-------------------------------------------------------------------------')
print('Ringdown parameters')
print('mode\t\tfrequency\tdamping time\tA/A_220\t\tphase')
print('-------------------------------------------------------------------------')
for mode in modes:
    i = modes.index(mode)
    print('%s\t%.2f\t\t%.2E\t%.2E\t%.2f'\
            %(mode,freqs[i],taus[i],amps[i]/amps[0],phis[i]))
print('')

## define signal
signal = RingdownMultimode(modes,freqs,taus,amps,phis)

## define Fisher
fisher = Fisher(signal,detector='aligo')

## compute SNR
snr = fisher.snr(nbins=1e5)
print('SNR: %.2f'%snr)
print('')

## compute fisher matrix
fisher_matrix = fisher.fisher_matrix(nbins=1e5)

## compute 1-D uncertainties
sigmas = fisher.sigma1d(fisher_matrix)

print('-------------------------------------------------------------------------')
print('1-D uncertainties:')
print('mode\t\tfrequency\tdamping time\tlog10(A)\tphase')
print('-------------------------------------------------------------------------')
for k,v in sigmas.items():
    print('%s\t%.2E\t%.2E\t%.2E\t%.2E'%(k,*v))
print('')

## compute measurabilities
## defined as relative percentage uncertainties
## over frequencies and damping times
print('--------------------------------------------')
print('Measurabilities [%]:')
print('mode\t\tfrequency\tdamping time')
print('--------------------------------------------')
for k,v in sigmas.items():
    i = modes.index(k)
    a = v[0]/freqs[i]*100
    b = v[1]/taus[i]*100
    print('%s\t%.2f\t\t%.2f'%(k,a,b))
print('')

