import Extractor
import matplotlib.pyplot as plt
import numpy as np

##TEST  
z=0  
radec='00h53m13.81s +13d09m55.0s'
#data=extractor(radec)
#spec_list= downloader(data)
#dic = transform_data(spec_list,z)

data = Extractor.extractor(radec)
spec_list = Extractor.downloader(data)
dic = Extractor.transform_data(spec_list,z)

#print(dic)

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

plot_spec(dic, radec, z)