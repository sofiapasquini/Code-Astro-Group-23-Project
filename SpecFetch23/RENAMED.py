import Extractor
import matplotlib.pyplot as plt
import numpy as np

def plot_spec(dict, radec, z): # takes as input the dictionary holding the data, the radec, and the redshift

    for i in range(len(dict['wavelength'])):
        #extract data
        wavelength = dict['wavelength'][i]
        sigma = dict['1sigma'][i]
        flux = dict['flux'][i]

        # instantiate a figure object
        fig=plt.figure()
        plt.title(str(radec)+str('; ')+'z={}'.format(z))
        plt.xlabel("Rest-frame Wavelength [$\AA$]")
        plt.ylabel("Flux [$10^{-17}$ erg$^{-1}$s$^{-1}$cm$^{-2}$$\AA^{-1}$]")
        plt.plot(wavelength, flux) # plot the actual data
        # now create upper and lower bounds on the uncertainty regions
        sigmaUpper=np.add(flux,sigma)
        sigmaLower=np.subtract(flux,sigma)
        plt.fill_between(wavelength, sigmaLower, sigmaUpper, color='grey', alpha=0.5)

        plt.show()

radec='22h38m12.39s +21d32m03.4s'
z=Extractor.redshift(radec)
data=Extractor.extractor(radec)
spec_list=Extractor.downloader(data)
dic = Extractor.transform_data(spec_list,z)
# print(dic)

plot_spec(dic, radec, z)
