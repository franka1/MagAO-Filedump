# Author: Alex Frank
# Contributor: Jennifer Vezilj

from astropy.io import fits
import numpy as np
import sys
import os
from scipy import signal
from photutils import daofind
#from astropy.stats import mad_std
from photutils import aperture_photometry, CircularAperture
import matplotlib.pylab as plt
import argparse

def calc_sigma(image, corner_width, corner_height):
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
    return sigma


def process_file(inpath, file_name, t_constant, sigma, fwhm, r, kernel_size, outpath, plot):
    print "Processing " + file_name
    hdulist = fits.open(inpath + file_name)
    image = hdulist[0].data

    if isinstance(sigma, list):
        threshold = calc_sigma(image, sigma[0], sigma[1]) * t_constant
    else:
        threshold = t_constant*sigma

    median_out = signal.medfilt(image,kernel_size)
    median_sub = np.subtract(image,median_out)
    sources = daofind(median_sub, threshold, fwhm)

    sources_2 = np.array(sources["id", "xcentroid", "ycentroid", "sharpness", "roundness1", "roundness2", "npix", "sky", "peak", "flux", "mag"])
    print_line= (file_name+","+str(sources_2))
    base_name = os.path.splitext(file_name)[0]
    file = open(outpath + base_name + ".out", "a")
    file.write(print_line)
    file.close()

    positions = (sources['xcentroid'], sources['ycentroid'])
#    print positions
    apertures = CircularAperture(positions, r)
    phot_table = aperture_photometry(median_sub, apertures)
    phot_table_2 = np.array(phot_table["aperture_sum", "xcenter", "ycenter"])
    print_line= (","+str(phot_table_2)+"\n")
    file = open(outpath + base_name + ".out", "a")
    file.write(print_line)
    file.write("\n")
    file.close()

    hdulist[0].data = median_sub
    file = open(outpath + base_name + ".fits", "w")
    hdulist.writeto(file)
    file.close()

    if plot:
        median_sub[median_sub<=0]=0.0001
        plt.imshow(median_sub, cmap='gray', origin='lower')
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Find stars in an image.')
    parser.add_argument('-i','--input', type=str, required=True,
                        help='The path to input file or directory. Trailing / is require for directories.')
    parser.add_argument('-o','--output', type=str, required=True, help='The path to output directory.')
    parser.add_argument('-f','--fwhm', type=float, default=4.7, help='Specifies fwhm. (default: 4.7)')
    parser.add_argument('-r','--radius', type=int, default=12, help='Specifies radius of stars in image. (default:12)')
    parser.add_argument('-t','--threshold', type=int, default=1,
                        help='Specifies the threshold constant to be used in concert with sigma. (default: 1)')
    parser.add_argument('-s','--calc_sigma', nargs=2, metavar=('CORNER_WIDTH','CORNER_HEIGHT'),
                        help='Calculate sigma using specified corner width and height.')
    parser.add_argument('-S','--sigma', type=float, help='Use the supplied sigma value.')
    parser.add_argument('-k','--kernel', type=int, default=17, help='Specifies kernel size for median filter. (default: 17)')
    parser.add_argument('-p','--plot', action='store_true', help='Show plotted stars after for each image')

    args = parser.parse_args()

    if args.sigma:
        sigma = args.sigma
    else:
        sigma = [int(args.calc_sigma[0]),int(args.calc_sigma[1])]

    if os.path.isdir(args.input):
        file_names = os.listdir(args.input)
        for file_name in file_names:
            process_file(args.input, file_name, args.threshold, sigma, args.fwhm, args.radius, args.kernel, args.output, args.plot)
    else:
        process_file(args.input, "", args.threshold, sigma, args.fwhm, args.radius, args.kernel, args.output, args.plot)

if __name__ == "__main__":
	main()
