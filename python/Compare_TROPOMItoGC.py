#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import numpy as np
import xarray as xr
import re
import pickle
import os
import pandas as pd
import datetime
import scipy.integrate as integrate
from scipy.integrate import quad

#----- define function -------
def save_obj(obj, name):
    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open( name, 'rb') as f:
        return pickle.load(f)

def read_tropomi(filename):
    met={}
    #--- read methane, qa_value, longitude, latitude ---
    data=xr.open_dataset(filename,group="PRODUCT")

    # Methane concentration (ppb)
    met['methane']=data['methane_mixing_ratio_bias_corrected'].values[0,:,:]

    # qa value
    met['qa_value']=data['qa_value'].values[0,:,:]

    # Lat and lon
    met['longitude']=data['longitude'].values[0,:,:]
    met['latitude']=data['latitude'].values[0,:,:]

    # Precision (ppb)
    met['precision']=data['methane_mixing_ratio_precision'].values[0,:,:]

    # Times
    referencetime=data['time'].values
    delta_time=data['delta_time'][0].values
    strdate=[]
    if delta_time.dtype=='<m8[ns]':
        strdate=referencetime+delta_time
    elif delta_time.dtype=='<M8[ns]':
        strdate=delta_time
    else:
        print(delta_time.dtype)
        pass
    timeshift=np.array(met['longitude']/15*60,dtype=int)#convert to minutes
    met['utctime']=strdate

    localtimes=np.zeros(shape=timeshift.shape,dtype='datetime64[ns]')
    for kk in range(timeshift.shape[0]):
        item=pd.to_timedelta(timeshift[kk,:], unit='m')
        localtimes[kk,:]=strdate[kk]+item.values
    data.close()
    met['localtime']=localtimes

    # Albedo, AOD, TROPOMI averaging kernel
    data=support = xr.open_dataset(filename, group="PRODUCT/SUPPORT_DATA/DETAILED_RESULTS")
    met['albedo']=data['surface_albedo_SWIR'].values[0,:,:]
    met['aerosol_optical_depth']=data['aerosol_optical_thickness_SWIR'].values[0,:,:]
    met['column_AK']=data['column_averaging_kernel'].values[0,:,:,::-1]
    data.close()

    data=xr.open_dataset(filename,group="PRODUCT/SUPPORT_DATA/INPUT_DATA")

    # Methane a priori profile (mol m-2)
    met['methane_profile_apriori']=data['methane_profile_apriori'].values[0,:,:,::-1]

    # Pressure information (Pa)
    pressure_interval=data['pressure_interval'].values[0,:,:]
    surface_pressure=data['surface_pressure'].values[0,:,:]

    # Dry air subcolumns (mol m-2)
    met['dry_air_subcolumns']=data['dry_air_subcolumns'].values[0,:,:,::-1]
    data.close()

    # Get data shape
    N1=met['methane'].shape[0]
    N2=met['methane'].shape[1]

    # Create a pressures array corresponding to 13 vertical levels?
    pressures=np.zeros([N1,N2,13], dtype=np.float)
    pressures.fill(np.nan)
    pressures_mid=np.zeros([N1,N2,12], dtype=np.float)
    pressures_mid.fill(np.nan)
    for i in range(12+1):
        pressures[:,:,i]=surface_pressure-i*pressure_interval
        if i < 12:
            pressures_mid[:,:,i]=surface_pressure-(i+0.5)*pressure_interval
    # HN 2020/07/07 : converted from for loop to numpy


    met['pressures']=pressures
    met['pressures_mid']=pressures_mid

    return met

def get_diagnostic(diag_name, date):
    return os.path.join(GC_datadir,
                        'GEOSChem.' + diag_name + '.' + date + '.0000z.nc4')

def read_GC(date):
    month=int(date[4:6])
    short_date = date[:8]
    hour = int(date[-2:])

    #--- read base GC -----
##
#These L,I,J dims have I and J flipped for lon&lat to lat&lon.
#However, if I correct this, then TOP isn't able to concatenate corectly b/c it has the wrong dimensions
##

    # Start by downloading methane data (ppb), latitude, and longitude
    species_filename=get_diagnostic('SpeciesConc', short_date)
    # species_filename=GC_datadir+'/GEOSChem.SpeciesConc.'+date[0:8] +'.0000z.nc4'
    species_data=xr.open_dataset(species_filename)
    #subset for the relevant hours
    species_data=species_data.where(species_data.time.dt.hour == hour,
                                    drop=True).squeeze()
    CH4=species_data['SpeciesConc_CH4'].values*1e9#[0,:,:,:]*1e9 #mrew
    CH4=np.einsum('lij->jil',CH4)
    LON=species_data['lon'].values
    LAT=species_data['lat'].values
    LEV=species_data.lev.values[:]

    species_data.close()

    # Now get the other varaibles
    #making date only consider month and day, not time
    met_filename=get_diagnostic('StateMet', short_date)
    # met_filename=GC_datadir+'/GEOSChem.StateMet.'+date[0:8] +'.0000z.nc4'
    #print('GC file', filename)
    met_data=xr.open_dataset(met_filename)
    met_data=met_data.where(met_data.time.dt.hour == hour,
                            drop=True).squeeze()

    # PBL height (m)
    TROPP=met_data['Met_PBLH'].values#[0,:,:]
    TROPP=np.einsum('ij->ji',TROPP)

    # Air density (kg m-3)
    AIRDEN=met_data['Met_AIRDEN'].values#[0,:,:,:]
    AIRDEN=np.einsum('lij->jil',AIRDEN)

    # Box height wrt dry air (m)
    BXHGHT=met_data['Met_BXHEIGHT'].values#[0,:,:,:]
    BXHGHT=np.einsum('lij->jil',BXHGHT)

    # Dry air column (kg m-2)
    DRYAIR = np.multiply(AIRDEN, BXHGHT) #* 100

    # Albedo
    ALBD=met_data['Met_ALBD'].values#[0,:,:]
    ALBD=np.einsum('ij->ji',ALBD)

    # Dry air mass
    DRYAIRMASS = met_data['Met_AD'].values
    DRYAIRMASS = np.einsum('lij->jil', DRYAIRMASS)

    # Total column
    TOTCOL = (CH4*DRYAIRMASS).sum(axis=2)/DRYAIRMASS.sum(axis=2)

    met_data.close()

    # Pressure information (hPa)
    # pedge_filename=GC_datadir+'/GEOSChem.LevelEdgeDiags.'+date[0:8] +'.0000z.nc4'
    pedge_filename=get_diagnostic('LevelEdgeDiags', short_date)
    pedge_data=xr.open_dataset(pedge_filename)
    pedge_data=pedge_data.where(pedge_data.time.dt.hour == hour,
                               drop=True).squeeze()

    PEDGE=pedge_data['Met_PEDGE'].values#[0,:,:,:]*100
    PEDGE=np.einsum('lij->jil',PEDGE)
    # LON=pedge_data['lon'].values
    # LAT=pedge_data['lat'].values

    pedge_data.close()

    #what is this
    TOP=np.ones([len(LON),len(LAT)],dtype=float)*0.01
    PEDGE=np.dstack((PEDGE,TOP))

    #CH4_adjusted=CH4.copy()

    # Latitudinal stratosphere correction?
    # for i in range(len(LON)):
    #     for j in range(len(LAT)):
    #         l=int(TROPP[i,j])
    #         CH4_adjusted[i,j,l:]=CH4[i,j,l:]*lat_ratio[j,month-1]

    met={}
    met['lon']=LON
    met['lat']=LAT
    met['PEDGE']=PEDGE
    met['CH4']=CH4
    met['LEV']=LEV

    #met['CH4_adjusted']=CH4_adjusted

    met['TROPP']=TROPP
    met['DRYAIR']=DRYAIR
    met['ALBD']=ALBD
    met['GCCOL']=TOTCOL

    #--- read sensitivity ---

    #filename=Sensi_datadir+'/'+date+'0000.nc'
    #print,'sensfile',filename
    #data=xr.open_dataset(filename)
    #Sensi=data['Sensi'].values
    #Sensi=np.einsum('klji->ijlk',Sensi)
    #data.close()
    # Latitudinal stratosphere correction
    # for i in range(len(LON)):
    #     for j in range(len(LAT)):
    #         l=int(TROPP[i,j])
    #         Sensi[i,j,l:,:]=Sensi[i,j,l:,:]*lat_ratio[j,month-1]
    #met['Sensi']=Sensi

    return met

# def read_pert(date,ipert):
#     month=int(date[4:6])

#     #--- read base GC -----
#     pertdir = '/n/holyscratch01/jacob_lab/zhenqu/inv2019/PertA/'

#     filename=pertdir+'PertA_'+str(ipert).zfill(4)+ '/nc/nc_ts'+date +'0000.nc'
#     #print('pert file', filename)
#     data=xr.open_dataset(filename)
#     LON=data['LON'].values
#     LAT=data['LAT'].values
#     CH4=data['IJ-AVG-S__CH4'].values
#     CH4=np.einsum('lij->jil',CH4)

#     AIRDEN=data['TIME-SER__AIRDEN'].values
#     AIRDEN=np.einsum('lij->jil',AIRDEN)
#     BXHGHT=data['BXHGHT-S__BXHEIGHT'].values
#     BXHGHT=np.einsum('lij->jil',BXHGHT)

#     DRYAIR = np.multiply(AIRDEN, BXHGHT) * 100

#     data.close()

#     filename=GC_datadir+'/nc_ts'+date +'0000.nc'
#     #print('GC file', filename)
#     data=xr.open_dataset(filename)
#     TROPP=data['PBLDEPTH__PBL-L'].values[0,:,:]
#     TROPP=np.einsum('ij->ji',TROPP)
#     data.close()

#     TOP=np.ones([len(LON),len(LAT)],dtype=float);TOP.fill(0.01)

#     CH4_adjusted=CH4.copy()
#     # Latitudinal stratospheric correction?
#     # for i in range(len(LON)):
#     #     for j in range(len(LAT)):
#     #         l=int(TROPP[i,j])
#     #         CH4_adjusted[i,j,l:]=CH4[i,j,l:]*lat_ratio[j,month-1]

#     met={}
#     met['lon']=LON
#     met['lat']=LAT
#     met['CH4']=CH4
#     met['CH4_adjusted']=CH4_adjusted
#     met['DRYAIR']=DRYAIR

#     return met


def read_all_GC(all_strdate):
    met={}

    # reduce all_strdate to daily list (no time, nt)
    # all_strdate_nt = [d[:8] for d in all_strdate]
    # all_strdate_nt = list(set(all_strdate_nt))
    for strdate in all_strdate:
        met[strdate]= read_GC(strdate)
    return met

# def read_allpert_GC(all_strdate):
#     met={}
#     for strdate in all_strdate:
#         for ipert in range(1009):
#         #for ipert in range(99):
#             met[(strdate,ipert)]= read_pert(strdate,ipert+1)
#     return met

# quzhen 2020/2/13
def get_intmap(Sat_p, GC_p):

    intmap = np.zeros((len(GC_p)-1, len(Sat_p)))
    count = 0
    for ltm in range(len(Sat_p)):
        if Sat_p[ltm] > GC_p[0]:
            count+=1
            intmap[:,ltm] = 1e0

            if ltm > 0:
                for ll in range(ltm):
                  intmap[1:,ll] = 0e0
                  intmap[0,ll] = 1.0/count
                intmap[0,ltm] = 1.0/count
        elif Sat_p[ltm]<GC_p[len(GC_p)-1]:
            intmap[len(GC_p)-2,ltm] = 1e0
        else:
            for lgc in range(len(GC_p)-1):
              low = GC_p[lgc+1]
              hi = GC_p[lgc]
              if (Sat_p[ltm]<=hi and Sat_p[ltm]>low):
                diff = hi - low

                if (ltm == 0):
                  intmap[0:lgc+1,ltm] = 1e0

                if (ltm > 0):
                  intmap[lgc,ltm-1] = (hi-Sat_p[ltm])/diff
                  if(lgc < len(GC_p)-1):
                    intmap[lgc+1:len(GC_p)-1,ltm-1] = 0e0
                  intmap[lgc,ltm] = (Sat_p[ltm]-low)/diff

                if (lgc < len(GC_p)-1):
                  intmap[lgc+1:len(GC_p),ltm] = 1e0
    return intmap


def newmap(intmap,lgos, GC_p, Sat_p,gc_ch4_native,dryair):
    gc_ch4 = np.zeros(lgos-1)
    gc_weight = np.zeros(lgos-1)
    count = 0e0
    for ll in range(lgos-1):
      temp_gc = 0e0
      temp_dry = 0e0
      for l in range(len(GC_p)-2): #changed from len()-1 to len()-2 mrew
        temp_gc += abs(intmap[l,ll]) * gc_ch4_native[l] * dryair[l]
        temp_dry += abs(intmap[l,ll]) * dryair[l]
        count += abs(intmap[l,ll]) * dryair[l]
      gc_ch4[ll] = temp_gc / temp_dry
      gc_weight[ll] = temp_dry / np.sum(dryair)
      #gc_weight[ll] = temp_dry / np.sum(dryair*intmap[:,ll])


    met={}
    met['GC_CH4']=gc_ch4
    met['GC_WEIGHT']=gc_weight

    return met

def newmap2(intmap,lgos, GC_p, Sat_p,gc_sens,dryair):
    gc_ch4 = np.zeros((lgos-1,1009))
    count = 0e0
    for ll in range(lgos-1):
      temp_gc = 0e0
      temp_dry = 0e0
      for l in range(len(GC_p)-1):
        temp_gc += abs(intmap[l,ll]) * gc_sens[l,:] * dryair[l]
        temp_dry += abs(intmap[l,ll]) * dryair[l]
        count += abs(intmap[l,ll]) * dryair[l]
      gc_ch4[ll,:] = temp_gc / temp_dry
    met={}
    met['Sens']=gc_ch4

    return met

def remap2(Sensi, data_type, Com_p, location, first_2):
    MM=Sensi.shape[1]
    conc=np.zeros((len(Com_p)-1,MM))
    conc.fill(np.nan)
    k=0
    for i in range(first_2,len(Com_p)-1):
        conc[i,:]=Sensi[k,:]
        if data_type[i+1]==2:
            k=k+1
    if first_2>0:
        conc[:first_2,:]=conc[first_2,:]

    Sat_CH4=np.zeros((12,MM));Sat_CH4.fill(np.nan)

    delta_p=Com_p[:-1]-Com_p[1:]
    delta_ps=np.transpose(np.tile(delta_p,(MM,1)))
    for i in range(len(location)-1):
        start=location[i]
        end=location[i+1]
        fenzi=np.sum(conc[start:end,:]*delta_ps[start:end,:],0)
        fenmu=np.sum(delta_p[start:end])
        Sat_CH4[i,:]=fenzi/fenmu

    return Sat_CH4

#==============================================================================
#===========================Define functions ==================================
#==============================================================================
#Sat_datadir="/n/seasasfs02/hnesser/TROPOMI/downloads_201910/"
Sat_datadir="/n/holyscratch01/jacob_lab/mwinter/newTROPOMI/"
GC_datadir="/n/holyscratch01/jacob_lab/mwinter/Nested_NA/run_dirs/Hannah_NA_0000/OutputDir_2018/"
outputdir="/net/seasasfs02/srv/export/seasasfs02/share_root/mwinter/TROPOMI_processed/data_2018_ch4col/"
biasdir="/net/seasasfs02/srv/export/seasasfs02/share_root/mwinter/TROPOMI_processed/bias/"
#Sensi_datadir="/n/holyscratch01/jacob_lab/zhenqu/aggregate/data/"

# #==== read lat_ratio ===
# df=pd.read_csv("./lat_ratio.csv",index_col=0)
# lat_mid=df.index
# lat_ratio=df.values

#==== read Satellite ===

# List all raw netcdf TROPOMI files
allfiles=glob.glob(Sat_datadir+'*.nc')

# Create empty list
Sat_files=[]

# Iterate through the raw TROPOMI data
for index in range(len(allfiles)):
    filename=allfiles[index]

    # Get the date (YYYY, MM, and DD) of the raw TROPOMI file
    shortname=re.split('\/', filename)[-1]
    shortname=re.split('\.', shortname)[0]
    strdate=re.split('\.|_+|T',shortname)[4]
    YYYY=int(strdate[:4])
    MM=int(strdate[4:6])
    DD=int(strdate[6:8])

    # Skip observations not in range
   # number = [i for i in range(2,9)]
    if not ((YYYY==2018 and MM==5)):
        continue

    # Add the file to the list of Sat_files
    Sat_files.append(filename)

# Sort by date and print the number of files
Sat_files.sort()
print("Number of files",len(Sat_files))

# Create an array that corresponds to the state vector
b = np.zeros((46,72))
bcount = np.zeros((46,72))

# Iterate throught the Sat_files we created
#for index in range((1-1)*400,1*400):
#for index in range(0,1):
for index in range(1):#range(0,len(Sat_files)):
    print('========================')

    # Again, get the date of the sat file in question
    # (could potentially be improved by using a dictionary
    # instead of a list)
    filename=Sat_files[index]
    temp=re.split('\/', filename)[-1]
    print(temp)
    date=re.split('\.',temp)[0]

    # If already processed, skip the rest of the processing
    # within this loop
    # if os.path.isfile(outputdir+date+'_GCtoTROPOMI.pkl'):
    #     continue

    #--- read TROPOMI ---
    # This function creates a dictionary that saves out the
    # bias corrected methane mixing ratio, qa, lon, lat, precision,
    # utctime, localtime, column averaging kernel, the a priori
    # methane profile, the dry air sub columns, pressures, and
    # pressure midpoints
    TROPOMI=read_tropomi(filename)#read TROPOMI data

    # Get the indices of the data with sufficiently high quality
    # data (really this should just be a NN x 1 array...)
    #sat_ind=np.where((TROPOMI['qa_value']>=0.5) & (TROPOMI['utctime']>=GC_startdate) & (TROPOMI['utctime']<=GC_enddate))
    sat_ind=np.where((TROPOMI['qa_value']>=0.5))

    # observation dimension (number of good observations in that single
    # observation file)
    NN=len(sat_ind[0])

    # state vector dimension (we will need to change this)
    # MM=1009

    # create an empty matrix for the Jacobian
    # temp_KK=np.zeros([NN,MM],dtype=np.float32)#Store the K

    # create an empty matrix to store TROPOMI CH4, GC CH4,
    # lon, lat, II, and JJ (GC indices)
    temp_obs_GC=np.zeros([NN,10],dtype=np.float32)#TROPOMI-CH4, GC-CH4, longitude,latitude, II, JJ

    #================================
    #--- now compute sensitivity ---
    #================================
    #--- generate all strdate----

    # Create a list to contain all the dates and hours
    # in YYYYMMDD.HH format as strings
    all_strdate=[]

    # Iterate through observation dimension
    for iNN in range(NN):

        # Get indices of the particular good satellite observation
        # in question (i.e. qa > 0.5)
        iSat=sat_ind[0][iNN] # row
        jSat=sat_ind[1][iNN] # col

        # Get the utc time and convert it to YYYYMMDD.HH
        # format. Add this date and time to all_strdate
        utctime=TROPOMI['utctime'][iSat]
        utctime=pd.to_datetime(str(utctime))
        strdate=utctime.strftime("%Y%m%d.%H")
        all_strdate.append(strdate)

    # Remove duplicates
    # This is now a list of the day/hour at which every observation
    # in a single TROPOMI file occurs at
    all_strdate=list(set(all_strdate))

    # Then, read in the GC data for these dates. This works by
    # reading the lon, lat, pressure edge, xch4, xch4_adjusted
    # (which I believe is the stratospheric corrected data), TROPP
    # (which is the planetary boundary layer info), and dry air.
    # It also reads sensitivity data
    all_date_GC=read_all_GC(all_strdate)
    #all_pert_GC=read_allpert_GC(all_strdate)

    #print('all_strdate', all_strdate)
    #temp = all_pert_GC[('20190101.02',0)]
    #print(temp)
    #quit()
    # pert = np.zeros([NN,1009],dtype=np.float)

    # Iterate through observations (individual pixels)
    for iNN in range(NN):
        # Get indices of good TROPOMI observations
        iSat=sat_ind[0][iNN]
        jSat=sat_ind[1][iNN]

        # Get the pressures, dry air subcolumns, prior methane
        # profile, pressure mid points, and averaging kernel
        # of the true TROPOMI retrieval
        Sat_p=TROPOMI['pressures'][iSat,jSat,:]
        dry_air_subcolumns=TROPOMI['dry_air_subcolumns'][iSat,jSat,:]#mol m-2
        priori=TROPOMI['methane_profile_apriori'][iSat,jSat,:]
        # quzhen 2020/2/5
        Sat_pmid=TROPOMI['pressures_mid'][iSat,jSat,:]
        AK=TROPOMI['column_AK'][iSat,jSat,:]

        # I don't know what this does
        timeshift=int(TROPOMI['longitude'][iSat,jSat]/15*60)
        utctime=TROPOMI['utctime'][iSat]
        utctime=pd.to_datetime(str(utctime))
        hour = str(utctime.hour)
        strdate=utctime.strftime("%Y%m%d.%H")

        # Get the GC data associated with the data we found above
        # (lon, lat, pedge, ch4, tropp (pbl depth?), dryair)
        GC=all_date_GC[strdate]

        # Find the grid box indices corresponding to the good TROPOMI
        # observations
        # Individual TROPOMI obs: TROPOMI['lat/lon'][iSat, jSat] (scalar)
        # Enitre GC grid: GC['lon'] (vector)
        iGC=np.abs(GC['lon']-TROPOMI['longitude'][iSat,jSat]).argmin()
        jGC=np.abs(GC['lat']-TROPOMI['latitude'][iSat,jSat]).argmin()

        # Then find the pressure and dry air (kg/m2?) for the column associated
        # with that grid box (i.e. i=iGC, j=jGC, lev=ALL)
        GC_p=GC['PEDGE'][iGC,jGC,:]
        dryair = GC['DRYAIR'][iGC,jGC,:]

	#mrew
        GC_albd=GC['ALBD'][iGC,jGC]

        # Column calculation
        #GC_CH4=GC['CH4_adjusted'][iGC,jGC,:]
        # GC_CH4=GC['CH4'][iGC,jGC,:]
        # GC_CH4_col=integrate.simps(GC['CH4'][iGC,jGC,:]*((GC['LEV']*10**6*6.023*10**23)/(8.31)), x=None, dx=1)/integrate.simps((GC['LEV']*10**6*6.023*10**23)/(8.31), x=None, dx=1)
        GC_CH4 = GC['CH4'][iGC, jGC, :]
        GC_COL = GC['GCCOL'][iGC, jGC]

        # quzhen 2020/2/13
        # I suspect that she is creating an interpolated pressure grid
        intmap = get_intmap(Sat_p, GC_p)
        temp = newmap(intmap, len(Sat_p), GC_p, Sat_p, GC_CH4, dryair)
        Sat_CH4 = temp['GC_CH4']
        GC_WEIGHT = temp['GC_WEIGHT']

        # quzhen 2020/2/5
        temp_gc = 0e0
        temp_gcpri = 0e0

        # Applies the averaging kernel to the geos-chem observations
        # This is the TROPOMI operator
        # For reference: see GEOS-Chem GOSAT Obs module
        # ! Apply GOSAT observation operator
        # !
        # !   Xch4_m = Xch4_a + SUM_j( h_j * a_j * (x_m - x_a) )
        # !
        # !   Xch4_m  - model XCH4
        # !   Xch4_a  - apriori XCH4 = h^T * x_a
        # !   h       - pressure weighting function
        # !   a       - column averaging kernel
        # !   x_m     - model CH4 [v/v]
        # !   x_a     - apriori CH4 [v/v]
        # !
        # !   The pressure weighting function is defined in Connor et al. 2008
        # !     and the OCO-2 ATBD
        # !--------------------------------------------------------------
        for ll in range(12):
            temp_gc += GC_WEIGHT[ll]*(priori[ll]/dry_air_subcolumns[ll]*1e9+AK[ll]*(Sat_CH4[ll] - priori[ll]/dry_air_subcolumns[ll]*1e9))
            temp_gcpri += GC_WEIGHT[ll]*(priori[ll]/dry_air_subcolumns[ll]*1e9+(Sat_CH4[ll] - priori[ll]/dry_air_subcolumns[ll]*1e9))

        # I suspect(though I don't know) that GC_base_posteri is the
        # GC observations with averaging kernel applied
        GC_base_posteri = temp_gc
        GC_base_pri = temp_gcpri

        #Sensi=GC['Sensi'][iGC,jGC,:,:]
        #temp = newmap2(intmap, len(Sat_p), GC_p, Sat_p, Sensi, dryair)
        #Sens = temp['Sens']
        #print(Sens.shape)
        #temp_gcsens = np.zeros(1009)
        #for ll in range(12):
        #    temp_gcsens[:] += GC_WEIGHT[ll]*AK[ll]*Sens[ll,:]

            # perturbation = temp_gcsens-temp_gc
            # pert[iGC,jGC,isens] += (temp_gcsens-temp_gc)/0.5 # for grid aggregate
        #pert[iNN,:] = temp_gcsens/0.5 # for observation individual



        #print('GC_pos', GC_base_posteri)
        #print('GC_pri', GC_base_pri)

        temp_obs_GC[iNN,0]=TROPOMI['methane'][iSat,jSat]
        temp_obs_GC[iNN,1]=GC_base_posteri
        temp_obs_GC[iNN,2]=TROPOMI['longitude'][iSat,jSat]
        temp_obs_GC[iNN,3]=TROPOMI['latitude'][iSat,jSat]
        temp_obs_GC[iNN,4]=iGC
        temp_obs_GC[iNN,5]=jGC
        temp_obs_GC[iNN,6]=TROPOMI['precision'][iSat,jSat]
        temp_obs_GC[iNN,7]=TROPOMI['albedo'][iSat,jSat]
        temp_obs_GC[iNN,8]=TROPOMI['aerosol_optical_depth'][iSat,jSat]
	#mrew - not sure
        # temp_obs_GC[iNN,9]=GC_albd
        temp_obs_GC[iNN,9]=GC_COL

        # Total error (unclear why not abs) in each grid cell
        #b[jGC, iGC] += GC_base_posteri - TROPOMI['methane'][iSat,jSat]

        # Number of observations in each grid cell
        #bcount[jGC, iGC] += 1

    result={}
    result['obs_GC']=temp_obs_GC
    # result['KK']=pert


    save_obj(result,outputdir+date+'_GCtoTROPOMI.pkl')
#b[bcount>0] = b[bcount>0]/bcount[bcount>0]
#save_obj(b,biasdir+'1.pkl')
print('CODE FINISHED')
