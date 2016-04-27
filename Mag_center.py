# Author: Jennifer Vezilj
# Editor: Alex Frank
# This is a cleaned copy of clio_sub.py

from astropy.io import fits
import numpy as np
import sys
import os
from photutils import daofind
from astropy.stats import mad_std
from photutils import aperture_photometry, CircularAperture
import matplotlib.pylab as plt

def calc_threshold(image, corner_width, corner_height):
    sum = 0.0
    n = 0
    #first corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += image[i][j]
            n += 1
    #second corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += image[len(image)-i-1][j]
            n += 1
    #thrid corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += image[i][len(image)-j-1]
            n += 1
    #fourth corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += image[len(image)-i-1][len(image)-j-1]
            n += 1

    mu = sum/n
    sum = 0

    #first corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += (image[i][j] - mu)**2
    #second corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += (image[len(image)-i-1][j] - mu)**2
    #thrid corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += (image[i][len(image)-j-1] - mu)**2
    #fourth corner
    for i in range(corner_width):
        for j in range(corner_height):
            sum += (image[len(image)-i-1][len(image)-j-1] - mu)**2

    sigma2 = sum/(n-1)
    sigma = sigma2**(0.5)
    threshold = 5*sigma
    return threshold


def process_file(inpath, file_name, corner_width, corner_height, fwhm, r, outpath):
    hdulist = fits.open(inpath + file_name)
    image = hdulist[0].data

    threshold = calc_threshold(image, corner_width, corner_height)

    bkg_sigma = mad_std(image)
    sources = daofind(image, threshold*bkg_sigma, fwhm)

    sources_2 = np.array(sources["id", "xcentroid", "ycentroid", "sharpness", "roundness1", "roundness2", "npix", "sky", "peak", "flux", "mag"])
    print_line= (file_name+","+str(sources_2))
    file= open(outpath, "a")
    file.write("\n")
    file.write(print_line)
    file.close()

    positions = (sources['xcentroid'], sources['ycentroid'])
    print positions
    apertures = CircularAperture(positions, r)
    phot_table = aperture_photometry(image, apertures)
    phot_table_2 = np.array(phot_table["aperture_sum", "xcenter", "ycenter"])
    print_line= (","+str(phot_table_2)+"\n")
    file= open(outpath, "a")
    file.write(print_line)
    file.close()

    image[image<=0]=0.0001
    plt.imshow(image, cmap='gray', origin='lower')
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()


def main():
    inpath = str(sys.argv[1])
#    file_name = str(sys.argv[2])
    outpath = str(sys.argv[2])
    fwhm = float(sys.argv[3])
    r = float(sys.argv[4])
    corner_width = int(sys.argv[5])
    corner_height = int(sys.argv[6])

    if os.path.isdir(inpath):
        file_names = os.listdir(inpath)
        for file_name in file_names:
            process_file(inpath, file_name, corner_width, corner_height, fwhm, r, outpath)

    else:
        process_file(inpath, "", corner_width, corner_height, fwhm, r, outpath)

if __name__ == "__main__":
	main()
