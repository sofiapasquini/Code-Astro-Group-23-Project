#define functions that will extract the data from SDSS based on an input RA/DEC

from astroquery.sdss import SDSS
from astropy import coordinates as coords
import pandas as pd 
from astroquery.ned import Ned 
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Box1DKernel
import numpy as np
from astropy import units as u


def ra_dec_format(val):
    """ Ra/Dec string formatting

    Converts the input string format of a right ascension/ declination coordinate
    to one recognizable by astroquery

    Args:
      val (str): string; an ra/dec expression formatted as "005313.81 +130955.0".

    Returns:
      string: the ra/dec coordinates re-formatted as "00h53m13.81s +13d09m55.0s"
    """
    #ra
    hour = val[0:2]
    min_ = val[2:4]
    sec = val[4:9]
    ra = hour+'h'+min_+'m'+sec+'s'
    #dec
    deg = val[9:13]
    min_d = val[13:15]
    sec_d = val[15:]
    dec = deg+'d'+min_d+'m'+sec_d+'s'
    return ra+" "+dec

def extractor(position):
  """Spectral data extraction from SDSS.

  Extract the information from the SDSS database (query by a region surrounding the input ra/dec) and return
  a pandas dataframe holding the query results. Please ensure that the 'position'
  input is formatted as '005313.81 +130955.0

  Args:
    position (str) : string; an ra/dec expression formatted as "005313.81 +130955.0".

  Returns:
    pandas DataFrame: the results returned by astroquery.SDSS.query_region().
  """

  # convert the input position argument to the format recognized by astroquery.SDSS
#   position=ra_dec_format(position)

  # query the region and get the data
  position = ra_dec_format(position)
  pos = coords.SkyCoord(position, frame='icrs')
  data = SDSS.query_region(pos, spectro=True)
  return data.to_pandas()


def downloader(data):
  """ Download database query results

  Uses extracted information from extractor() in order to download spectra availible for 
  the queried ra/dec.

  Args:
    data (pd.DataFrame): the output results from extractor(); a dataframe holding the results returned
    by astroquery.SDSS.query_region(). Returns the data portions of the original output fits files in a list format.

  Returns:
    list: a list of fits files holding the data sections of the resuls files from the SDSS database.
  """
  #create a empty list
  spec_list=[]

  # iteration over the pandas
  for i in range(len(data)):
    results = SDSS.query_specobj(plate   = data['plate'][i],
                                 mjd     = data['mjd'][i],
                                 fiberID = data['fiberID'][i])
    
    # try if it can download the data (SDSS)
    try:
      spec    = SDSS.get_spectra(matches=results)[0]
      spec_list.append(spec)

    # if it cant download, is because is from (BOSS)
    except:
      results.remove_column("instrument")
      results.add_column(name="instrument", col="eboss") # replace the instrument column
      spec    = SDSS.get_spectra(matches=results)[0]
      spec_list.append(spec)

  return spec_list 



# test=downloader(result)
# print(test)

# define a function which grabs the object's redshift from the Ned database (better calibration)- needed for plotting in the object's rest-frame
def redshift(position):
  """ Query Ned for object's redshift

  Query the Ned database using astroquery.Ned by a right ascension/ declination coordinate for the corresponding redshift value.

  Args:
    position (str): string; an ra/dec expression formatted as "005313.81 +130955.0".

  Returns:
    float: the redshift of the object located at the input ra/dec.


  """

    # make sure to format the input position argument such that it is recognizable by astroquery.Ned
    # position=ra_dec_format(position)
    position = ra_dec_format(position)
    pos=coords.SkyCoord(position, frame='icrs') # create a position object
    ned_results=Ned.query_region(pos,equinox="J2000", radius=2*u.arcsecond) # query the database
    z=float(ned_results[0][6]) # grab the redshift value from the query results
    return z

# define a function that transforms an objects wavelength array into the object's rest-frame
def redshift_correct(z, wavelengths): 
    """ Redshift correction

    Correct a given numpy array of wavelengths for redshift effects using an accepted redshift value.

    Args:
      z (float): a redshift value

      wavelengths (array): a numpy array containing wavelenth values

    Returns:
      array: a numpy array of wavelength values of the same size as input wavelengths array which has been 
      corrected to the emitted wavelength values for the rest-frame of the object in quetion.
    """



    wavelengths_corrected = wavelengths/(z+1)
    return wavelengths_corrected

# define a function that transforms the results of downloader() into an array of data which will be plotted
def transform_data(spec_list, z): 
    """Organize output of downloader()

    Sort the information from the output list of spectra from downloader() (results of a 
    query by astroquery.SDSS.query_region) into numpy arrays containing flux, wavelength, and one-sigma (one-sigma uncertainty associated with each flux value) values
    in a nested dictionary (each of the three data arrays per epoch availible).

    Args:
      spec_list (list): a list of fits files holding the data sections of the resuls files from the SDSS database; output from
      downloader().

      z (float): a redshift value

    Returns:
      dictionary : a (nested) dictonary containing flux, wavelength, and one-sigma arrays for each availible epoch. 
      The numerical keys representing number of availible epochs (1,2,3 ... etc), internal dictionaries
      for each contain numpy arrays for flux, wavelength, and one-sigma values.


    """
    # iterate over each file and grab the important data
    #fluxes={} # containers for each of the data arrays to be plotted ( will be lists of lists/arrays)
    #wavelengths={}
    #inverse_variances={} # <- dictionaries!

    dict={}

    for spec in spec_list:
        
        flux_array=[]
        wavelength_array=[]
        sigma_array=[]

        data=spec[1].data # this is the data part of the file
        #print(data.shape[0])
        #print(data)

        # store the appropriate columns in the designated containers- each row is a single spectrum?
        # SOFIA- try a nested dictionary?!?! 
        for j in range(data.shape[0]):
            #print(data[j][0])

            #smoothedFlux=convolve(data[0],Box1DKernel(9)) # smooth the fluxes using a boxcar
            #print(smoothedFlux)
            flux_data = data[j][0]
            flux_array.append(flux_data)
            
            wavelengths_uncorrected=10**data[j][1] # the wavelengths (transformed from the log scale)
            #print(wavelengths_uncorrected)
            wavelengths_corrected=redshift_correct(z, wavelengths_uncorrected) # save the wavelengths after they have been scaled to the rest-frame
            #print(wavelengths_corrected)
            wavelength_array.append(wavelengths_corrected)

            inverse_variance=data[j][2] # the inverse variance of the flux
            one_over_sigma=inverse_variance**0.5
            sigma=1/one_over_sigma # the one-sigma  uncertainty associated with the flux array
            sigma_array.append(sigma)
        
    smoothedFlux = convolve(flux_array,Box1DKernel(9))
    if 'flux' in dict:
        dict['flux'].append([smoothedFlux])
    else:
        dict['flux'] = [smoothedFlux]
        
    if 'wavelength' in dict:
        dict['wavelength'].append([wavelength_array])
    else:
        dict['wavelength'] = [wavelength_array]
        
    if '1sigma' in dict:
        dict['1sigma'].append([sigma_array])
    else:
        dict['1sigma'] = [sigma_array]

    # now return the nested dictionary with three keys:(flux, wavelength and sigma)
    # each key should have data.shape[0] number of arrays with all fluxes, wavelength and sigmas for every spec in spec_list
    return dict


def plot_spec(dict, radec, z): 
  """Plot availible spectra from the SDSS database

  Plots the spectra availible for a given object in the SDSS database queried by right ascension/ declination coordinate.
  Fluxes [10^-17 ergs/s/cm^2/Angstrom] are plotted with corresponding one-sigma uncertainty regsions as functions of Wavelength
  [Angstroms] in the object rest-frame.

  Args:
    dict (dictionary): a (nested) dictonary containing flux, wavelength, and one-sigma arrays for each availible epoch. 
      The numerical keys representing number of availible epochs (1,2,3 ... etc), internal dictionaries
      for each contain numpy arrays for flux, wavelength, and one-sigma values. The output of transform_data().

    radec (int/str): an ra/dec (of the object queried), used for the plot title

    z (int): the redshift (of the object queried), used for the plot title

  Returns:

    fig: A plot of the object's flux as a function of rest-frame wavelength. Associated one-sigma
    uncertainty regions are highlighted in grey above/below each spectrum, all availible epochs (from the SDSS database) are
    represented.
  """



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



#TEST
radec='223812.39 +213203.4'
z=redshift(radec)
data=extractor(radec)
spec_list=downloader(data)
dic = transform_data(spec_list,z)
plot_spec(dic, radec, z)