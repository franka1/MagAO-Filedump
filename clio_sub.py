# Author: Jennifer Vezilj
# Editor: Alex Frank

from astropy.io import fits
import numpy as np
import sys

# inputs: inpath, filename, outpath, fwhm, threshold, r


#read1 = []

#star=[9, 10, 15, 16]

#for i in range(0, 4):
#    if star[i] < 10:
#        path = "/Users/jennifervezilj/mypy/Python/data/clio_20141201_02/sub/"
#        file_name = "Test33_4ND0000%d.fit" %star[i]
#        read1.append(file_name)
#    else:
#        path = "/Users/jennifervezilj/mypy/Python/data/clio_20141201_02/sub/"
#        file_name = "Test33_4ND000%d.fit" %star[i]
#        read1.append(file_name)

inpath = str(sys.argv[1])
filename = str(sys.argv[2])
outpath = str(sys.argv[3])
fwhm = float(sys.argv[4])
threshold = float(sys.argv[5])
r = float(sys.argv[6])

#len(read1)

#for i in range(0, 4):
#    file_name = read1[i]

    hdulist = fits.open(path+file_name)
    image = hdulist[0].data
    #image = image.astype(float) - np.median(image)
    from photutils import daofind
    from astropy.stats import mad_std
    bkg_sigma = mad_std(image)
#    sources = daofind(image, fwhm=10., threshold=20.*bkg_sigma)
    sources = daofind(image, fwhm, threshold*bkg_sigma)
    #print_line= (file_name+","+str(sources_2)+"\n")
    sources_2 = np.array(sources["id", "xcentroid", "ycentroid", "sharpness", "roundness1", "roundness2", "npix", "sky", "peak", "flux", "mag"])
    print_line= (file_name+","+str(sources_2))
#    file= open("/Users/jennifervezilj/mypy/Python/data/clio_20141201_02/code/test33_4nd_phot.txt", "a")
    file= open(outpath, "a")
    file.write(print_line)
    file.close()

    from photutils import aperture_photometry, CircularAperture
    positions = (sources['xcentroid'], sources['ycentroid'])
    apertures = CircularAperture(positions, r)
    phot_table = aperture_photometry(image, apertures)
    phot_table_2 = np.array(phot_table["aperture_sum", "xcenter", "ycenter"])
    print_line= (","+str(phot_table_2)+"\n")
#    file= open("/Users/jennifervezilj/mypy/Python/data/clio_20141201_02/code/test33_4nd_phot.txt", "a")
    file= open(outpath, "a")
    file.write(print_line)
    file.close()

    import matplotlib.pylab as plt
    im2 = image
    im2[im2<=0]=0.0001
    plt.imshow(im2, cmap='gray', origin='lower')
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()
