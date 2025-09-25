import yt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from silmaril.lens import Lens  # or whichever API you use
from silmaril.utilities import Grid          # your pasted Grid class

# --- Step 1: Load dataset and make a projection ---
ds = yt.load("your_dataset")   # replace with your file
width = (100, "kpc")           # example width

proj = yt.ProjectionPlot(ds, "z", ("gas", "density"), width=width)
frb = proj.to_frb(width=(100, "kpc"), resolution=(256, 256))

# --- Step 2: Extract the FRB as a numpy array ---
image_data = np.array(frb[("gas", "density")])

# --- Step 3: Compute pixel scale ---
# Physical width in kpc
width_kpc = width[0] * ds.quan(1, width[1]).to("kpc").value
pixel_size_kpc = width_kpc / image_data.shape[0]

# Convert to angular scale given angular diameter distance
angular_diam_dist = (1000 * u.Mpc)  # replace with real z/distance
pixel_scale_rad = (
    (pixel_size_kpc * u.kpc / angular_diam_dist)
    .to(u.rad, equivalencies=u.dimensionless_angles())
    .value
)
pixel_scale_arcsec = pixel_scale_rad * (180/np.pi) * 3600  # rad → arcsec

# --- Step 4: Build Grid centered on some RA/Dec ---
center = SkyCoord(ra=150*u.deg, dec=2*u.deg, frame="icrs")  # pick center
grid = Grid(center, image_data.shape[0], pixel_scale_arcsec)

# --- Step 5: Apply lensing transformation ---
# Silmaril expects: (image array, grid of coordinates, lens model params)
# Example API: lens_image(source_image, grid, lens)
lensed_image = Lens(image_data, grid.as_2d_array(), lens_model)

# Now you can save or plot lensed_image