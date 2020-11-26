# ** THE CODE IN THIS REPOSITORY IS CURRENTLY BEING REFORMATTED AND IS NOT YET COMPLETE **

# Goldblatt_etal_2020


This repository contains the code used to analyze model output and make the figures for Goldblatt et al 2020: 'Clouds Stabilize Earth's Long Term Climate', doi:_________

Model output from this experiment is available and can be downloaded from doi:_________

Code is written in python and makes use of the Climate Data Operators (CDO) version 1.7.0 for processing the model output. For most scripts, CDO is used to extract the variable(s) from the model output, process it by time averaging, selecting a pressure level, or calculating a new variable, and then saves this processed output to a new NetCDF file, which is then plotted using Matplotlib.  
These processing steps can be time intensive, taking up to ~20 minutes before the figure will be plotted. Intervening steps will be written to the terminal to track progress. There are probably lots of ways this could be better optimized.

Since processing the model output is the more time-intensive step, saving the processed variables to their own files means that the processing only needs to be done once, and then subsequent changes to the figures are *much* faster. 

Also included is an environment.yml file - this describes the python libraries required to run the scripts and it is recommended to make a new conda
environment using this yaml file before you run the scripts. 

Importantly, you will also need to download and install the Climate Data Operators (https://code.mpimet.mpg.de/projects/cdo/wiki#Installation-and-Supported-Platforms),
NetCDF (https://www.unidata.ucar.edu/downloads/netcdf/), and HDF5 (https://www.hdfgroup.org/downloads/hdf5/source-code/) manually.
These three installations can be a bit troublesome but if you install HDF5 first, 
then NetCDF, then finally configure CDO before installing, it should go smoothly. 
** Please consult the support pages for these software packages if you have trouble.

Note that to read the large NetCDF2 files ***version 1.7.0*** of CDO must be used. If you have other versions of CDO installed, you may need to change the commands in processing_functions.py to the location of version 1.7.0 instead, ie:
replace instances of r"//usr//bin//cdo " with "cdo" if that is the only installed version, or the path to version 1.7.0 if different.

Also, the processed model output will be saved to new netCDF files that are then used to make the figures. The variable outfileloc in each figure script points to the location where these processed files will be saved, and the directory must already exist on your system, so please create a new directory (we recommend outside of this repository), or change to point to an existing directory.

General correspondence about the paper and methods should be directed to corresponding author Colin Goldblatt (czg@uvic.ca).  
Questions about the code in this repository can be sent to Victoria McDonald (vmcd@atmos.washington.edu), or feel free to open an issue in the GitHub repository.