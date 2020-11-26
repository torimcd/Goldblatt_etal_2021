This repository contains the code used to make the figures for Goldblatt et al 2020: 'Clouds Stabilize Earth's Long Term Climate', doi:________

Model output from this experiment is available and can be downloaded from doi:_________

Code is written in python and makes use of the Climate Data Operators (CDO) for processing the model output. For most scripts, CDO is used to extract
the variable(s) from the model output, process it by time averaging, selecting a pressure level, or calculating a new variable, and then saves this processed
output to a new NetCDF file, which is then plotted using Matplotlib. Thus, before running each figure plotting script, the processing scripts need to be run first.
Processing the model output is the more time-intensive step, so this way the processing only needs to be done once, and then subsequent changes to the figures 
are much faster. Please note that to support the large file sizes the model output is in NetCDF2 format, so CDO must be built with NetCDF2 support.

Also included is an environment.yml file - this describes the python libraries required to run the scripts and it is recommended to make a new conda
environment using this yaml file before you run the scripts. 

Importantly, you will also need to download and install the Climate Data Operators (https://code.mpimet.mpg.de/projects/cdo/wiki#Installation-and-Supported-Platforms),
NetCDF (https://www.unidata.ucar.edu/downloads/netcdf/), and HDF5 (https://www.hdfgroup.org/downloads/hdf5/source-code/) manually.
These three installations can be a bit troublesome but if you install HDF5 first, 
then NetCDF, then finally configure CDO with the --with-netcdf=<NetCDF2_root_directory> flag before installing, it should go smoothly. 
** Please consult the support pages for these software packages if you have trouble.


General correspondence about the paper and methods should be directed to corresponding author Colin Goldblatt (czg@uvic.ca). 
Questions about the code in this repository can be sent to Victoria McDonald (vmcd@atmos.washington.edu), or feel free to open an issue in the GitHub repository.