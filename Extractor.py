#define functions that will extract the data from SDSS based on an input RA/DEC

from astroquery.sdss import SDSS
from astropy import coordinates as coords
import pandas as pd 

# this function extracts the information from the database 
def extractor(position): # the format of the positon must be as follows: '0h8m05.63s +14d50m23.3s'
	pos=coords.SkyCoord(position, frame='icrs') # create a position object
	xid = SDSS.query_region(pos, spectro=True) # query the database
	xid.to_pandas() # convert to a pandas data frame
	return xid

# testing it out
radec='1h11m17.36s +14d26m53.6s'
result=extractor(radec)
# print(result)

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

        sdss_results=SDSS.query_specobj(plate=plateID, mjd=date, fiberID=fiber) # querying for matching spectra
        spec=SDSS.get_spectra(matches=results) # now actually download/get the spectra

        # add each query result to the list holding the results files
        spec_list.append(spec)

    return spec_list



test=downloader(result)
print(test)
