# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 15:25:12 2024

@author: Özgür Özcan
"""

from prysm.coordinates import make_xy_grid, cart_to_polar
from prysm.geometry import regular_polygon

from matplotlib import pyplot as plt

efl = 50   
fno = 5

x, y = make_xy_grid(256, diameter=efl/fno)
dx = x[0,1]-x[0,0]
r, t = cart_to_polar(x, y)
radius = efl/fno/2
rho = r / radius
n_sides = 14

aperture = regular_polygon(n_sides, radius, x, y)

plt.imshow(aperture, origin='lower')

from prysm.polynomials import zernike_nm
from prysm.propagation import Wavefront
wvl = 0.555 # mid visible band, um

wfe_nm_rms = wvl/14*1e3 # nm, 3/4 of a wave, 1e3 = um to nm
mode = zernike_nm(4, 0, rho, t)
opd = mode * wfe_nm_rms
pup = Wavefront.from_amp_and_phase(aperture, opd, wvl, dx)
coherent_psf = pup.focus(efl, Q=2)

from prysm.otf import mtf_from_psf, diffraction_limited_mtf
psf = coherent_psf.intensity
mtf = mtf_from_psf(psf)


fx, _ = mtf.slices().x
fig, ax = mtf.slices().plot(['x', 'y', 'azavg'], xlim=(0,400))
difflim = diffraction_limited_mtf(fno, wvl, fx)

ax.plot(fx, difflim, ls=':', c='k', alpha=0.75, zorder=1)
ax.set(xlabel='Spatial frequency, cy/mm', ylabel='MTF')

wfe_nm_rms = wvl/14*1e3
mode = zernike_nm(3, 1, rho, t) # only this line changed
opd = mode * wfe_nm_rms
pup = Wavefront.from_amp_and_phase(aperture, opd, wvl, dx)
coherent_psf = pup.focus(efl, Q=2)
psf = coherent_psf.intensity
mtf = mtf_from_psf(psf, psf.dx)

fig, ax = mtf.slices().plot(['x', 'y', 'azavg'], xlim=(0,400))

ax.plot(fx, difflim, ls=':', c='k', alpha=0.75, zorder=1)
ax.set(xlabel='Spatial frequency, cy/mm', ylabel='MTF')
