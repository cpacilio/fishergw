import sys
## import fishergw objects
sys.path.append('..')
from fishergw.taylorf2 import CompactObject, TaylorF2, Fisher

print('--------------------------------------------------------')
print('Compute the 1-D uncertainties on the parameters of a BNS')
print('using Fisher matrix formalism on a TaylorF2 template.\n')
print('--------------------------------------------------------')

## define binary components
params_1 = {}
params_2 = {}
params_1['mass'] = 1.6 # M_sun
params_2['mass'] = 1.4
params_1['spin'] = 0.
params_2['spin'] = 0.
params_1['Lamda'] = 300
params_2['Lamda']= 200
params_1['radius'] = 12e3
params_2['radius'] = 12e3 # meters
luminosity_distance = 100 # Mpc

for k,v in params_1.items():
    print('%s: %.2f'%(k,v))
print('')
for k,v in params_2.items():
    print('%s: %.2f'%(k,v))
print('')

obj1 = CompactObject(**params_1)
obj2 = CompactObject(**params_2)

## define signal
## set redshift=True to redshift the component masses
signal = TaylorF2(obj1,obj2,d_L=luminosity_distance,redshift=True)
print('luminosity distance: %.2f'%luminosity_distance)
print('redshift: %.2f'%signal.redshift)
print('Lamda tilde: %.2f'%signal.Lamda_T)
print('')

## define Fisher
keys = ['t_c','phi_c','M_c','eta','Lamda_T','chi_s','chi_a']
logscale_keys = ['M_c','eta'] # params differentiated in logscale
detector = 'aligo'
fisher = Fisher(signal,keys=keys,logscale_keys=logscale_keys,detector=detector)
## the above line is a shortcut for
# fisher = Fisher(signal,keys=keys,logscale_keys=logscale_keys)
# fisher.load_psd(asd_name='fishergw/detector/aligo_asd.txt',Qavg=2/5)

## define isco frequency
# mode=`static` uses the Schwarzschild isco
# mode=`contact` uses contact frequency
fmax = signal.isco(mode='static')

## define minimum frequency
fmin = 10

## compute snr
snr = fisher.snr(fmin,fmax,nbins=1e5)
print('SNR: %.2f\n'%snr)

## compute fisher matrix
fisher_matrix = fisher.fisher_matrix(fmin,fmax,nbins=1e5)

## compute 1-D uncertainties
sigmas = fisher.sigma1d(fisher_matrix)
print('1-D uncertainties:')
for k,v in sigmas.items():
    print('\t%s: %.2E'%(k,v))
print('')

## impose priors over the spins
priors = {'chi_s':0.05,'chi_a':0.05}
fisher_matrix = fisher.fisher_matrix(fmin,fmax,nbins=1e5,priors=priors)
sigmas = fisher.sigma1d(fisher_matrix)
print('1-D uncertainties with restrictive priors over the spins:')
for k,v in sigmas.items():
    print('\t%s: %.2E'%(k,v))


