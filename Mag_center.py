# Author: Jennifer Vezilj
# Editor: Alex Frank
# This is a cleaned copy of clio_sub.py

from astropy.io import fits
import numpy as np
import sys

inpath = str(sys.argv[1])
file_name = str(sys.argv[2])
outpath = str(sys.argv[3])
fwhm = float(sys.argv[4])
r = float(sys.argv[5])
# threshold = float(sys.argv[5])
# r = float(sys.argv[6])

hdulist = fits.open(inpath + file_name)
image = hdulist[0].data

#calc threshold
sum = 0.0
n = 0
#first corner
for i in range(99):
    for j in range(99):
        sum += image[i][j]
        n += 1
#second corner
for i in range(99):
    for j in range(99):
        sum += image[len(image)-i][j]
        n += 1
#thrid corner
for i in range(99):
    for j in range(99):
        sum += image[i][len(image)-j]
        n += 1
#fourth corner
for i in range(99):
    for j in range(99):
        sum += image[len(image)-i][len(image)-j]
        n += 1

mu = sum/n
sum = 0

#first corner
for i in range(99):
    for j in range(99):
        sum += (image[i][j] - mu)**2
#second corner
for i in range(99):
    for j in range(99):
        sum += (image[len(image)-i][j] - mu)**2
#thrid corner
for i in range(99):
    for j in range(99):
        sum += (image[i][len(image)-j] - mu)**2
#fourth corner
for i in range(99):
    for j in range(99):
        sum += (image[len(image)-i][len(image)-j] - mu)**2

sigma2 = sum/(n-1)
sigma = sigma2**(0.5)
threshold = 5*sigma

from photutils import daofind
from astropy.stats import mad_std
bkg_sigma = mad_std(image)
sources = daofind(image, fwhm, threshold*bkg_sigma)

sources_2 = np.array(sources["id", "xcentroid", "ycentroid", "sharpness", "roundness1", "roundness2", "npix", "sky", "peak", "flux", "mag"])
print_line= (file_name+","+str(sources_2))
file= open(outpath, "a")
file.write(print_line)
file.close()

from photutils import aperture_photometry, CircularAperture
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r)
phot_table = aperture_photometry(image, apertures)
phot_table_2 = np.array(phot_table["aperture_sum", "xcenter", "ycenter"])
print_line= (","+str(phot_table_2)+"\n")
file= open(outpath, "a")
file.write(print_line)
file.close()

import matplotlib.pylab as plt
image[image<=0]=0.0001
plt.imshow(image, cmap='gray', origin='lower')
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()