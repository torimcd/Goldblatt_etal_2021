% -----Updated by Victoria to use CAM output------
% Colin learns to read and process reanalysis data

% I put some data in my sandbox 
% this is one month-average of air temperatures

% read in Colin's test data
%T = ncread('~/projects/wetbulb/MERRA2_100.instM_3d_ana_Np.198001.nc4.nc','T');
%lat = ncread('~/projects/wetbulb/MERRA2_100.instM_3d_ana_Np.198001.nc4.nc','lat');
%lon = ncread('~/projects/wetbulb/MERRA2_100.instM_3d_ana_Np.198001.nc4.nc','lon');
%lev = ncread('~/projects/wetbulb/MERRA2_100.instM_3d_ana_Np.198001.nc4.nc','lev');
%time = ncread('~/projects/wetbulb/MERRA2_100.instM_3d_ana_Np.198001.nc4.nc','time');


% read in CAM output 1.0
T = ncread('~/projects/modelOutput/gcmexpt/CAM5/1.0/wbpt_pressure_.nc','T');
lat = ncread('~/projects/modelOutput/gcmexpt/CAM5/1.0/wbpt_pressure_.nc','lat');
lon = ncread('~/projects/modelOutput/gcmexpt/CAM5/1.0/wbpt_pressure_.nc','lon');
lev = ncread('~/projects/modelOutput/gcmexpt/CAM5/1.0/wbpt_pressure_.nc','lev');
time = ncread('~/projects/modelOutput/gcmexpt/CAM5/1.0/wbpt_pressure_.nc','time');


%% make zonal mean 1.0
Tmean = nanmean(T,1);
Tmean = squeeze(Tmean);


% %% plot temperature
% figure(10)
% clf
% contourf(lat,lev,Tmean')
% shading flat
% axis ij
% xlabel('latitude')
% ylabel('pressure')
% title('temperature')


%% calculate and plot potential temperature

% constants
R= 287.1;
cp = 1004;
p0 = 1e5;

% unit conversion
p = lev*100;
pp = ones(size(lat))*p';

sigma = Tmean.*(p0./pp).^(R/cp);

% figure(11)
% clf
% contourf(lat,lev,sigma',240:10:500)
% shading flat
% axis ij
% xlabel('latitude')
% ylabel('pressure')
% title('potential temperature')
% %caxis([240 500])
% colorbar


%% find and plot wet bulb potential temperature

% use some prevuiously wirtten code
%addpath('~/Documents/work/projects/runawaygreenhouse/tools/')
%addpath('~/Documents/work/projects/radcode/atmosutils/thermoprop/h2o/')

% inititalize
sigmaw = 999*ones(size(sigma));

% function we will integrate is
% dTdp = pseudoadiabatig(p,T,Rd,Md,cpd,condensablegas,condensedphase)
Rd = R;
cpd = cp;
Md = 0.02897;

options = odeset('reltol',1e-6');

parfor ii = 1:numel(sigmaw)
    pspan = [pp(ii) 1e5];
    if pspan(1) == pspan(2)
        sigmaw(ii) = Tmean(ii);
    else
        To = Tmean(ii);
        [~,Tw] = ode45(@(p,T) pseudoadiabatig(p,T,Rd,Md,cpd,'h2o','l'), pspan,To, options);
        sigmaw(ii) = Tw(end);
    end
end
%% plot
figure(12)
clf

contourf(lat,lev,sigmaw',240:5:320)
shading flat
axis ij
xlabel('latitude')
ylabel('pressure')
title('wet bulb potential temperature')
%caxis([240 500])
colorbar

%display(sigmaw)

%% write potential temp as nc file

%Open the file
ncid = netcdf.create(['./wbpt_1.0.nc'],'NC_WRITE');
 
%Define the dimensions
%dimidt = netcdf.defDim(ncid,'time',size(time));
%dimidp = netcdf.defDim(ncid,'sigmaw',size(sigmaw));
dimidlat = netcdf.defDim(ncid,'latitude',96);
dimidlon = netcdf.defDim(ncid,'lev',30);
 
%Define IDs for the dimension variables (pressure,time,latitude,...)
%time_ID=netcdf.defVar(ncid,'time','double',[dimidt]);
%pressure_ID=netcdf.defVar(ncid,'pressure','double',[dimidp]);
latitude_ID=netcdf.defVar(ncid,'latitude','double',[dimidlat]);
level_ID=netcdf.defVar(ncid,'lev','double',[dimidlon]);
 
%Define the main variable ()
sigmaw_ID = netcdf.defVar(ncid,'sigmaw','double',[dimidlat dimidlon]);
 
%We are done defining the NetCdf
netcdf.endDef(ncid);
 
%Then store the dimension variables in
%netcdf.putVar(ncid,time_ID,mytimearray);
%netcdf.putVar(ncid,pressure_ID,mypressurearray);
netcdf.putVar(ncid,latitude_ID,lat);
netcdf.putVar(ncid,level_ID,lev);
 
%Then store my main variable
netcdf.putVar(ncid,sigmaw_ID,sigmaw);
 
%We're done, close the netcdf
netcdf.close(ncid)
