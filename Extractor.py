#define functions that will extract the data from SDSS based on an input RA/DEC

from astroquery.sdss import SDSS
from astropy import coordinates as coords
import pandas as pd 
from astroquery.ned import Ned 

# this function extracts the information from the database 
def extractor(position): # the format of the positon must be as follows: '0h8m05.63s +14d50m23.3s'
	pos=coords.SkyCoord(position, frame='icrs') # create a position object
	xid = SDSS.query_region(pos, spectro=True) # query the database
	xid.to_pandas() # convert to a pandas data frame
	return xid

# testing it out
radec='1h11m17.36s +14d26m53.6s'
result=extractor(radec)
print(result)

#this function uses extracted information in order to dwonaload spectra
def downloader(xid):
    # start by splitting up instruments (SDSS vs eBOSS)
    boss_dataframe = xid[xid['instrument'] == b'BOSS']
    sdss_dataframe = xid[xid['instrument'] == b'SDSS']


    # make sure to fix the instrument column in the eBOSS rows
    boss_dataframe.remove_column("instrument")
    boss_dataframe.add_column(name="instrument", col="eboss")

    # change the tables to pandas dataframes
    boss_dataframe=boss_dataframe.to_pandas()
    sdss_dataframe=sdss_dataframe.to_pandas()

    # combine the two data frames now that they have been edited
    spec_frame=pd.concat([boss_dataframe, sdss_dataframe], ignore_index=False)

    # now we want to get the spectra for each entry in the spec_frame
    
    # create an empty list to hold all of the results files
    spec_list=[]

    # print(spec_frame)


    # query the database
    for i in range(len(spec_frame)):
        plateID=spec_frame['plate'][i] # get the plate number
        date=spec_frame['mjd'][i] # get the date (mjd)
        fiber=spec_frame['fiberID'][i] # get the fiber ID
        print(plateID)
        print(date)
        print(fiber)

        results=SDSS.query_specobj(plate=plateID, mjd=date, fiberID=fiber) # querying for matching spectra
        spec=SDSS.get_spectra(matches=results) # now actually download/get the spectra

        # add each query result to the list holding the results files
        spec_list.append(spec)

    return spec_list



# test=downloader(result)
# print(test)

# define a function which grabs the object's redshift from the Ned database (better calibration)- needed for plotting in the object's rest-frame
def redshift(position):
    pos=coords.SkyCoord(position, frame='icrs') # create a position object
    ned_results=Ned.query_region(pos,equinox="J2000", radius=2*u.arcsecond) # query the database
    z=ned_results[0][6] # grab the redshift value from the query results
    return z

# define a function that transforms an objects wavelength array into the object's rest-frame
def redshift_correct(z, wavelengths): # takes as input the redshift and the array of wavelengths
    wavelengths_corrected = wavelengths/(z+1)
    return wavelengths_corrected

# define a function that transforms the results of downloader() into an array of data which will be plotted
def transform_data(spec_list, z): # takes as input a list of (I think?) fits files results and the redshift of the object
    
    # iterate over each file and grab the important data
    fluxes={} # containers for each of the data arrays to be plotted ( will be lists of lists/arrays)
    wavelengths={}
    inverse_variances={} # <- dictionaries!

    dict={}

    for spec in spec_list:

        data=spec[0][1].data # this is the data part of the file

        # store the appropriate columns in the designated containers- each row is a single spectrum?
        # SOFIA- try a nested dictionary?!?! 
        for j in range(data.shape[0]):

            dict[j]['flux']=data[j][0] # the fluxes

            wavelengths_uncorrected=10**data[j][1] # the wavelengths (transformed from the log scale)
            wavelengths_corrected=redshift_correct(z, wavelengths_uncorrected) # save the wavelengths after they have been scaled to the rest-frame
            dict[j]['wavelength']=wavelengths_corrected

            inverse_variance=data[j][2] # the inverse variance of the flux
            one_over_sigma=inverse_variance**0.5
            sigma=1/one_over_sigma # the one-sigma  uncertainty associated with the flux array
            dict[j]['1sigma']=sigma

    # now return the nested dictionary, each key should have three arrays (flux, wavelength, and sigma)
    return dict

# now define a function which will plot the actual spectra given a spec dictionary
def plot_spec(dict, radec, z): # takes as input the dictionary holding the data, the radec, and the redshift
    # instantiate a figure object
    fig=plt.figure()
    plt.title(str(radec)+str('; ')+str("z"))
    plt.xlabel("Rest-frame Wavelength [$\AA$]")
    plt.ylabel("Flux [$10^{-17}$ erg$^{-1}$s$^{-1}$cm$^{-2}$$\AA^{-1}$]")
    for epoch in range(len(dict)):
        plt.plot(dict[epoch]['wavelength'], dict[epoch]['flux']) # plot the actual data
        # now create upper and lower bounds on the uncertainty regions
        sigmaUpper=np.add(dict[epoch]['flux'],dict[epoch]['1sigma'])
        sigmaLower=np.subtract(dict[epoch]['flux'],dict[epoch]['1sigma'])
        plt.fill_between(dict[epoch]['wavelength'],sigmaLower, sigmaUpper, color='grey', alpha='0.5')

    plt.show()