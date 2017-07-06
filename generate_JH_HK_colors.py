import pynrc

from pylab import *;ion()
import pysynphot as S

from scipy.interpolate import CubicSpline

def get_magnitudes(Teff, FeH=0.0, logg=4.5, Vmag=10):
    # Load the Filters
    bp_j = S.ObsBandpass('j')
    bp_h = S.ObsBandpass('h')
    bp_k = S.ObsBandpass('k')
    bp_v = S.ObsBandpass('johnson,v')
    
    # Stellar spectrum normalized to V=10 mags (default Phoenix models)
    sp = S.Icat(Teff, FeH, logg)#pynrc.stellar_spectrum(stellarType, Vmag, 'vegamag', bp_v)
    
    sp_norm = sp.renorm(Vmag, 'vegamag', bp_v)
    sp_norm.name = sp.name
    sp = sp_norm
    
    # Observe in J, H, and K
    obs_j = S.Observation(sp, bp_j, binset=sp.wave)
    obs_h = S.Observation(sp, bp_h, binset=sp.wave)
    obs_k = S.Observation(sp, bp_k, binset=sp.wave)
    
    # Magnitudes in each filter
    mag_j = obs_j.effstim('vegamag')
    mag_h = obs_h.effstim('vegamag')
    mag_k = obs_k.effstim('vegamag')
    
    return mag_j, mag_h, mag_k

# Loop over the 30 models
nModels   = 30
Jmags     = np.zeros(nModels)
Hmags     = np.zeros(nModels)
Kmags     = np.zeros(nModels)

teff_list = [x for x in range(2800,5500+100,100)] + [5800] + [6000]

# check if I did that write
assert(len(teff_list) == nModels)

for kt, teff in enumerate(teff_list):
    modelTeff = S.Icat('phoenix', teff, 0.0, 4.5) # load the model
    Jmags[kt], Hmags[kt], Kmags[kt] = get_magnitudes(stellarType = 'G2V')

jhmod = Jmags - Hmags
hkmod = Hmags - Kmags
