import pynrc

from pylab import *
import pysynphot as S

from scipy.interpolate import CubicSpline

def get_magnitudes(Teff, FeH=0.0, logg=4.5, Vmag=10, magScale='vegamag', modelgrid='ck04models'):
    # Load the Filters
    bp_j = S.ObsBandpass('j')
    bp_h = S.ObsBandpass('h')
    bp_k = S.ObsBandpass('k')
    bp_v = S.ObsBandpass('johnson,v')
    
    # Stellar spectrum normalized to V=10 mags (default Castelli & Kurucz 2004 models)
    sp = S.Icat(modelgrid, Teff, FeH, logg)#pynrc.stellar_spectrum(stellarType, Vmag, 'vegamag', bp_v)
    
    sp_norm = sp.renorm(Vmag, magScale, bp_v)
    sp_norm.name = sp.name
    sp = sp_norm
    
    # Observe in J, H, and K
    obs_j = S.Observation(sp, bp_j, binset=sp.wave)
    obs_h = S.Observation(sp, bp_h, binset=sp.wave)
    obs_k = S.Observation(sp, bp_k, binset=sp.wave)
    
    # Magnitudes in each filter
    mag_j = obs_j.effstim(magScale)
    mag_h = obs_h.effstim(magScale)
    mag_k = obs_k.effstim(magScale)
    
    return mag_j, mag_h, mag_k

# Loop over the 30 models
nModels   = 30
Jmags     = np.zeros(nModels)
Hmags     = np.zeros(nModels)
Kmags     = np.zeros(nModels)

Teff_list = [x for x in range(3500,5500+100,100)] + [5800] + [6000]

# check if I did that write
# assert(len(Teff_list) == nModels)

for kt, Teff in enumerate(Teff_list):
    Jmags[kt], Hmags[kt], Kmags[kt] = get_magnitudes(Teff, \
                                                     FeH=0.0, \
                                                     logg=4.5, Vmag=10, magScale='vegamag',modelgrid='ck04models')

jhmod = Jmags - Hmags
hkmod = Hmags - Kmags

### ----

import os 

from numpy import float32
from astroquery.irsa import Irsa
import astropy.coordinates as coords
import astropy.units as u

def deg2HMS(ra='', dec='', round=False):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    if str(dec)[0] == '-':
      ds, dec = '-', abs(dec)
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    if round:
      decS = int((abs((dec-deg)*60)-decM)*60)
    else:
      decS = (abs((dec-deg)*60)-decM)*60
    DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)
  
  if ra:
    if str(ra)[0] == '-':
      rs, ra = '-', abs(ra)
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    if round:
      raS = int(((((ra/15)-raH)*60)-raM)*60)
    else:
      raS = ((((ra/15)-raH)*60)-raM)*60
    RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC

def hmsdms2deg(ra,dec):
    hours2degs = 15.0
    mins2hours = 1/60.
    secs2hours = 1/3600.
    
    raOut   = float32(ra.split(':'))
    decOut  = float32(dec.split(':'))
    
    raOut   = raOut[0] + raOut[1]*mins2hours + raOut[2]*secs2hours
    raOut   = raOut*hours2degs
    decOut  = decOut[0] + decOut[1]*mins2hours + decOut[2]*secs2hours
    
    return raOut, decOut

# RA and Dec for HAT-P-11
ra     = '19:50:50.25'
dec    = '+48:04:51.1'
binComp=None

deg2rad     = np.pi/180
rad2deg     = 180/np.pi
deg2arcsec  = 3600
# binComp=[sourceDecRA,sourceDecDEC,J,H,K]

# stars in large field around target
targetCoords  = coords.SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg), frame='icrs')
fieldSources  = Irsa.query_region(targetCoords, catalog="fp_psc", spatial='Box', width=4 * u.arcmin)

sourceRA    = fieldSources['ra'].data.data  # in degrees
sourceDec   = fieldSources['dec'].data.data # in degrees
sourceJmag  = fieldSources['j_m'].data.data
sourceHmag  = fieldSources['h_m'].data.data
sourceKmag  = fieldSources['k_m'].data.data

# target coords
sky_distance= sqrt((sourceRA-targetCoords.ra.value)**2. + (sourceDec-targetCoords.dec.value)**2.)
targetIndex = np.argmin(sky_distance)

sourceJ_H   = (sourceJmag - sourceHmag)
sourceH_K   = (sourceHmag - sourceKmag)

ra,dec      = deg2HMS(sourceRA[targetIndex], sourceDec[targetIndex], round=True)

rcParams['figure.dpi'] = 300

# scatter(sourceJ_H, sourceH_K, s=50, c='k', alpha=0.5, lw = 0)
# sc = scatter(jhmod, hkmod,s=50,c=Jmags,cmap=cm.plasma)
scatter(sourceH_K, sourceJ_H, s=50, c='k', alpha=0.5, lw = 0)

sc = scatter(hkmod, jhmod, s=50, c=Jmags, cmap=cm.plasma)
xlabel('J-H')
ylabel('H-K')
title('Color-Color Synthetic vs Catalog')
cb = colorbar(sc)
cb.set_label('J-Mag')
savefig('JH_HK_color_color_plots_sunlike_stars.png')
# plt.show()
