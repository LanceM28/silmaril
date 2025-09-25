import yt
import numpy as np
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from silmaril.lens import Lens
from astropy.io import fits
from silmaril.utilities import Grid
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

#TODO:
#Add correct particle fields
#understand what wcs is 
#what is the actual output

# Load dataset
f1 = "/Users/lamoreau/python/ASpec/SimulationFiles/output_00273"
ds = yt.load(f1)

# Example projection (along z-axis, density)
proj = yt.ProjectionPlot(ds, "z", ("gas", "density"))

# Convert to a Fixed Resolution Buffer (FRB)
# width = physical size of the projection, resolution = output pixels
width = 50 * u.kpc
res   = (512, 512)
frb = proj.data_source.to_frb(width, res)

# Pull out multiple fields as numpy arrays
density     = np.array(frb[("gas", "density")])
temperature = np.array(frb[("gas", "temperature")])
pressure    = np.array(frb[("gas", "pressure")])

print("Density FRB stats:")
print("Min:", np.nanmin(density))
print("Max:", np.nanmax(density))
print("Number of zeros:", np.sum(density==0))
print("Number of NaNs:", np.sum(np.isnan(density)))

# Stack into layers: shape (n_layers, ny, nx)
layers = np.stack([density, temperature, pressure], axis=0)

#imports for part 2 formerly here 

# Cosmology for angular scale conversion
cosmo = FlatLambdaCDM(H0=70, Om0=0.3) #DUCK

# Suppose your source is at z = 2
z_source = 10.0
angular_diam_dist = cosmo.angular_diameter_distance(z_source)  # in Mpc

# FRB dimensions
ny, nx = res
width_kpc = width.to("kpc").value

# Pixel scale in radians
pixel_size_kpc = width_kpc / nx
pixel_scale_rad = ((pixel_size_kpc * u.kpc / angular_diam_dist)
                   .to(u.rad, equivalencies=u.dimensionless_angles()).value)
pixel_scale_arcsec = np.degrees(pixel_scale_rad) * 3600

# Build WCS #DUCK
wcs = WCS(naxis=2)
wcs.wcs.crpix = [nx/2, ny/2]          # reference pixel = center
wcs.wcs.cdelt = [-pixel_scale_arcsec, pixel_scale_arcsec]  # arcsec/pixel
wcs.wcs.crval = [0.0, 0.0]            # reference coord (arbitrary center)
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

#former silmaril loading

# Example deflection maps shipped in silmaril
x_defl = fits.getdata("docs/source/notebooks/hlsp_relics_model_model_whl0137-08_glafic_v1_x-arcsec-deflect.fits")
y_defl = fits.getdata("docs/source/notebooks/hlsp_relics_model_model_whl0137-08_glafic_v1_y-arcsec-deflect.fits")
lens_wcs = WCS(fits.getheader("docs/source/notebooks/hlsp_relics_model_model_whl0137-08_glafic_v1_x-arcsec-deflect.fits"))

# Create lens object
lens = Lens(x_defl, y_defl, wcs=lens_wcs, redshift=0.3) #TODO: Fix what redshift this is, is this in the .FITS file


#grid stuff

#the original GPT build, 

# Build a Grid with the FRB dimensions and WCS
#lens = Lens(nx, ny, wcs=wcs, redshift = )
"""
# Trace through the lens
x_src, y_src = lens.trace_grid(grid, source_redshift=z_source)

# Flatten for interpolation
coords = np.array([y_src.ravel(), x_src.ravel()])

lensed_layers = []
for layer in layers:
    lensed = map_coordinates(layer, coords, order=1, mode="constant", cval=0.0)
    lensed = lensed.reshape(ny, nx)
    lensed_layers.append(lensed)

lensed_layers = np.stack(lensed_layers, axis=0)  # (n_layers, ny, nx)
"""

#new one with grid doing things as designed
# Example: grid centered at RA=150 deg, Dec=2 deg
center = SkyCoord(ra=150*u.deg, dec=2*u.deg, frame="icrs")

num_pix = nx   # number of pixels per side
scale = 0.2    # arcsec per pixel, replace with your actual pixel scale

grid = Grid(center, num_pix, scale)

# Trace through the lens
x_src, y_src = lens.trace_grid(grid, source_redshift=z_source)

# Flatten for interpolation
coords = np.array([y_src.ravel(), x_src.ravel()])

lensed_layers = []
for layer in layers:
    lensed = map_coordinates(layer, coords, order=1, mode="constant", cval=0.0)
    lensed = lensed.reshape(ny, nx)
    lensed_layers.append(lensed)

lensed_layers = np.stack(lensed_layers, axis=0)  # (n_layers, ny, nx)

from astropy.visualization.wcsaxes import WCSAxes


layer_names = ["Density", "Temperature", "Pressure"]

fig = plt.figure(figsize=(15,5))

for i, layer in enumerate(lensed_layers):
    print(np.max(layer))
    print(np.min(layer))
    ax = fig.add_subplot(1, len(lensed_layers), i+1, projection=grid.wcs)
    im = ax.imshow(layer, origin="lower", cmap="viridis")
    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")
    ax.set_title(layer_names[i])
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()
