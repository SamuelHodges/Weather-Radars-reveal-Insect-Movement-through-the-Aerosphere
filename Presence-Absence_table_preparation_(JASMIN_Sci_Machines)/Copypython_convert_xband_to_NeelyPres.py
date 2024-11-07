# Neely Filter WSR to insects

# Import packages
import os
import sys
import datetime
import wradlib
import pyart
import numpy as np
import netCDF4
from numpy import emath as e
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import pylab
from matplotlib import colors as c
import cartopy.crs as ccrs
import warnings
import glob
import cmcrameri.cm as cm

# sys.argvs #
# 1 - full filename, can set "all" to loop through directory
# 2 - working directory

#change to working directory and check arguments are set
import sys
import os

if (len(sys.argv)-1)<1 :
    print("Error, filename must be defined")
    sys.exit()

if (len(sys.argv)-1)==2 :
    os.chdir(sys.argv[2])

badfiles=[]

filename=sys.argv[1]
if filename=="all":
    filename=glob.glob("*SUR_v1.nc")
elif len(filename)==8:
    filename=glob.glob("*"+filename+"*SUR_v1.nc")
else:
    filename=[filename]


if len(filename)==0:
    print("Error, no files found with user-given tag, exiting")
    exit()

for i in filename:
    # report place
    if len(glob.glob(i[0:-9]+"NeelyPresences.nc")) > 0 :
        print("Presence-Absence variant of "+i+" already exists, skipping")
        continue
    
    print(i)
    
    # get data and set up workspace
    
    
    os.environ['WRADLIB_DATA']=os.getcwd()
    file=wradlib.util.get_wradlib_data_file(i)
    try:
        radar=wradlib.io.read_generic_netcdf(file)	#wradlib for general processing and plots
        pyart_radar=pyart.io.read(i)	#pyart for specialised operations eg. despeckling, dealiasing
    except OSError:
        print("OSError exception raised, check "+i+" for data corruption or incorrect structure")
        badfiles.append(i)
    #dict.keys(radar['variables'])
    
    Range=np.ma.getdata(radar['variables']['range']['data'])
    Azimuth=np.ma.getdata(radar['variables']['azimuth']['data'])
    Elevation=np.ma.getdata(radar['variables']['elevation']['data'])
    sitecoords=(float(np.ma.getdata(radar['variables']['longitude']['data'])), float(np.ma.getdata(radar['variables']['latitude']['data'])), float(np.ma.getdata(radar['variables']['altitude']['data'])))

    dBZ=np.ma.getdata(radar['variables']['dBZ']['data'])
    Time=np.ma.getdata(radar['variables']['time']['data'])
    ZDR=np.ma.getdata(radar['variables']['ZDR']['data'])
    RhoHV=np.ma.getdata(radar['variables']['RhoHV']['data'])
    V=np.ma.getdata(radar['variables']['V']['data'])
    SNR=np.ma.getdata(radar['variables']['SNR']['data'])

    # calculate lon lat and alt matrices
    guide=wradlib.georef.spherical_to_proj(r=Range, phi=Azimuth, theta=Elevation, sitecoords=sitecoords)
    
    Longitude=guide[:,:,0].transpose()
    Latitude=guide[:,:,1].transpose()
    Altitude=guide[:,:,2].transpose()
    
    dBZ=dBZ.transpose()
    ZDR=ZDR.transpose()
    RhoHV=RhoHV.transpose()
    V=V.transpose()
    SNR=SNR.transpose()

    # get mask fill value
    dimfill=radar['variables']['azimuth']['data'].get_fill_value()
    maskfill=radar['variables']['ZDR']['data'].get_fill_value()
    
    #set nans
    dBZnanindex=np.argwhere(np.ma.getmask(radar['variables']['dBZ']['data']).transpose())
    ZDRnanindex=np.argwhere(np.ma.getmask(radar['variables']['ZDR']['data']).transpose())
    RhoHVnanindex=np.argwhere(np.ma.getmask(radar['variables']['RhoHV']['data']).transpose())
    Vnanindex=np.argwhere(np.ma.getmask(radar['variables']['V']['data']).transpose())
    SNRnanindex=np.argwhere(np.ma.getmask(radar['variables']['SNR']['data']).transpose())

    dBZ[dBZnanindex[0:dBZnanindex.shape[0],0], dBZnanindex[0:dBZnanindex.shape[0],1]]=np.nan
    ZDR[ZDRnanindex[0:ZDRnanindex.shape[0],0], ZDRnanindex[0:ZDRnanindex.shape[0],1]]=np.nan
    RhoHV[RhoHVnanindex[0:RhoHVnanindex.shape[0],0], RhoHVnanindex[0:RhoHVnanindex.shape[0],1]]=np.nan
    V[Vnanindex[0:Vnanindex.shape[0],0], Vnanindex[0:Vnanindex.shape[0],1]]=np.nan
    SNR[SNRnanindex[0:SNRnanindex.shape[0],0], SNRnanindex[0:SNRnanindex.shape[0],1]]=np.nan

    # Remove data from elevations below 2 degrees to exclude ground clutter
    elenan=np.argwhere(Elevation<2)

    dBZ[:,elenan[:,0]]=np.nan
    ZDR[:,elenan[:,0]]=np.nan
    RhoHV[:,elenan[:,0]]=np.nan
    V[:,elenan[:,0]]=np.nan
    SNR[:,elenan[:,0]]=np.nan
    
    # Kilambi et al. 2018 Depolarisation Ratio Filter
    lZDR=pow(10, ZDR/10)
    dr = 10 * e.logn(10, ( lZDR + 1 - 2 * e.power(lZDR, 1/2) * RhoHV ) / ( lZDR + 1 + 2 * e.power(lZDR, 1/2) * RhoHV ))

    dr_isnan = np.isnan(dr)
    dbz_isnan = np.isnan(dBZ)
    print(np.unique(dr_isnan, return_counts=True))

    Presence = np.logical_and((dr > -12),  (dBZ < 35), (SNR > 50)).astype(float)
    
    DR=dr
    
    # set absences to correct types (True Absence, Weather)
    Weather = np.argwhere(np.logical_and(Presence==0.0, np.isnan(dBZ)==False))

    Presence[dr_isnan] = 0
    Presence[dbz_isnan] = 0

    #Presence[NegligableSignals[:,0],NegligableSignals[:,1]] = 0
    #Presence[:,(MaxEffectiveRangeIndice+1):] = -2
    Presence[Weather[:,0],Weather[:,1]] = -1
    
    ## Despeckle true presence
    pyart_Presence=np.copy(Presence)
    pyart_FalsePresence=np.argwhere(pyart_Presence>1.5)
    pyart_Presence[pyart_FalsePresence[:,0], pyart_FalsePresence[:,1]]=-5 #temporarily set false presences to -ve so as not to interfere with the despeckling filter


    pyart_radar.add_field('Kilambi_Presence', {'data':Presence.transpose()}) #this is the original Kilambi filtered data
    PresenceDespeckle=pyart.correct.despeckle_field(pyart_radar, 'Kilambi_Presence', threshold=0.5, size=30)
    pyart_radar.add_field('Kilambi_Presence_Despeckle', {'data':PresenceDespeckle.gate_included}) #this is just the despeckle filter

    SpeckleMask=PresenceDespeckle.gate_included
    Presence=Presence*SpeckleMask.transpose() # Set speckle as 'No Significant Signal'
    pyart_radar.add_field('Kilambi_Presence_NoSpeckle', {'data':Presence.transpose()}) #this is with the despeckle filter on non-meteorological presence applied

    
    ## Remove Known Radar shadows, their edge effects, and the ZDR ground clutter cone (manually remove as part of Effective Range to avoid false positives IDs from Kilambi)

    ZDR_Cone=np.argwhere(np.logical_and( Azimuth > 175, Azimuth < 185 ))
    Radar_shadows1=np.argwhere(np.logical_and( Azimuth > 54, Azimuth < 65 ))
    Radar_shadows2=np.argwhere(np.logical_and( Azimuth > 82,  Azimuth < 91 ))
    Radar_shadows3=np.argwhere(np.logical_and( Azimuth > 128, Azimuth < 155 ))

    Mask=np.concatenate((ZDR_Cone,Radar_shadows1,Radar_shadows2,Radar_shadows3), axis=0)

    Presence[:,Mask[:,0]] = -2
    
    # write file
    Presence_rawpres=np.copy(Presence)
    
    writenetcdf=netCDF4.Dataset(i[0:-9]+"NeelyPresences.nc", "w", format="NETCDF4")

    # dims
    range=writenetcdf.createDimension("range", int(Range.shape[0]))
    azimuth=writenetcdf.createDimension("azimuth", int(Azimuth.shape[0]))
    elevation=writenetcdf.createDimension("elevation", int(Elevation.shape[0]))
    time=writenetcdf.createDimension("time", int(Time.shape[0]))
    sitecoordinatedim=writenetcdf.createDimension("sitecoordinates", int(3))

    ranges=writenetcdf.createVariable("range", "f4",("range",), fill_value=dimfill)
    azimuths=writenetcdf.createVariable("azimuth", "f4",("azimuth",), fill_value=dimfill)
    elevations=writenetcdf.createVariable("elevation", "f4",("elevation",), fill_value=dimfill)
    times=writenetcdf.createVariable("time", "f4",("time",), fill_value=dimfill)
    sitecoordinates=writenetcdf.createVariable("sitecoordinates", "f4",("sitecoordinates"))

    ranges[:]=Range
    azimuths[:]=Azimuth
    elevations[:]=Elevation
    times[:]=Time
    
    print("dimensions written")
    
    sitecoordinates[:]=sitecoords

    # main vars
    dBZs=writenetcdf.createVariable("dBZ", "f4",("range","time"), fill_value=maskfill)
    ZDRs=writenetcdf.createVariable("ZDR", "f4",("range","time"), fill_value=maskfill)
    RhoHVs=writenetcdf.createVariable("RhoHV", "f4",("range","time"), fill_value=maskfill)
    DRs=writenetcdf.createVariable("DR", "f4",("range","time"), fill_value=maskfill)
    Longitudes=writenetcdf.createVariable("Longitude", "f4",("range","time"), fill_value=maskfill)
    Latitudes=writenetcdf.createVariable("Latitude", "f4",("range","time"), fill_value=maskfill)
    Altitudes=writenetcdf.createVariable("Altitude", "f4",("range","time"), fill_value=maskfill)
    
    Presences=writenetcdf.createVariable("Presence", "f4",("range", "time"), fill_value=-2)
    #MaxEffectiveRanges=writenetcdf.createVariable("MaxEffectiveRange", "i4")
    #MinimumDetectableSignals=writenetcdf.createVariable("MinimumDetectableSignals", "f4")
    
    dBZs[:]=dBZ
    ZDRs[:]=ZDR
    RhoHVs[:]=RhoHV
    DRs[:]=dr
    Longitudes[:]=Longitude
    Latitudes[:]=Latitude
    Altitudes[:]=Altitude
    
    Presences[:]=Presence_rawpres
    #MaxEffectiveRanges=MaxEffectiveRange
    #MinimumDetectableSignals=MinimumDetectableSignal

    print("Key variables written")

    Vs=writenetcdf.createVariable("radial_velocity_h", "f4",("range","time"), fill_value=maskfill)
    Vs[:]=V

    print("Optional variables written")

    # close new netcdf
    writenetcdf.close()

    # cleanup
    del writenetcdf
    
    # plot Neely PPIs
    
    rangeindicemax=np.count_nonzero(np.less_equal(Range, 40000))
    RangesLim=Range[0:rangeindicemax]/1000

    ScanElevation1=np.unique(Elevation)[3]
    ScanElevation2=np.unique(Elevation)[5]

    strElevation1=str(ScanElevation1)
    strElevation2=str(ScanElevation2)


    elevationindices1=np.argwhere(np.equal(Elevation, ScanElevation1))
    elevationindices2=np.argwhere(np.equal(Elevation, ScanElevation2))

    Azimuths1=Azimuth[elevationindices1.flatten()]

    dBZLim1=dBZ[0:rangeindicemax, elevationindices1.flatten()]

    ZDRLim1=ZDR[0:rangeindicemax, elevationindices1.flatten()]
    RhoHVLim1=RhoHV[0:rangeindicemax, elevationindices1.flatten()]
    DRLim1=DR[0:rangeindicemax, elevationindices1.flatten()]

    VLim1=V[0:rangeindicemax, elevationindices1.flatten()]
    SNRLim1=SNR[0:rangeindicemax, elevationindices1.flatten()]

    Presence_Lim1=Presence_rawpres[0:rangeindicemax, elevationindices1.flatten()]

    Azimuths2=Azimuth[elevationindices2.flatten()]

    dBZLim2=dBZ[0:rangeindicemax, elevationindices2.flatten()]

    ZDRLim2=ZDR[0:rangeindicemax, elevationindices2.flatten()]
    RhoHVLim2=RhoHV[0:rangeindicemax, elevationindices2.flatten()]
    DRLim2=DR[0:rangeindicemax, elevationindices2.flatten()]

    VLim2=V[0:rangeindicemax, elevationindices2.flatten()]
    SNRLim2=SNR[0:rangeindicemax, elevationindices2.flatten()]

    Presence_Lim2=Presence_rawpres[0:rangeindicemax, elevationindices2.flatten()]

    del radar
    del pyart_radar
    
    # Plot figure settings (raw Radar variables)

    fig1 = plt.figure(figsize=[48, 14])

    title=i[38:53]
    fig1_title="Raw Radar Variables - "+title[0:4]+"/"+title[4:6]+"/"+title[6:8]+" "+title[9:11]+":"+title[11:13]+":"+title[13:15]
    fig1.suptitle(fig1_title, size=28, weight=5)

    cmap=c.ListedColormap(np.flip(['limegreen', 'black', 'mediumblue', 'whitesmoke']))
    RhoHVnorm=c.TwoSlopeNorm(vcenter=0.8, vmin=0, vmax=1)

    ax1, pm = wradlib.vis.plot_ppi(data=dBZLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig1, ax=pylab.subplot(2,4,1), cmap=plt.get_cmap('pyart_HomeyerRainbow', 20), vmax=20, vmin=-20)
    #cb1=plt.colorbar(pm, ax=ax1, pad=0.12, shrink=0.8, label="decibels")
    ax7, pm = wradlib.vis.plot_ppi(data=dBZLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig1, ax=pylab.subplot(2,4,2), cmap=plt.get_cmap('pyart_HomeyerRainbow', 20), vmax=20, vmin=-20)
    cb7=plt.colorbar(pm, ax=ax7, pad=0.12, shrink=0.8, label="decibels (dB)")

    ax2, pm = wradlib.vis.plot_ppi(data=ZDRLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig1, ax=pylab.subplot(2,4,3), cmap=plt.get_cmap('viridis', 12), vmax=10, vmin=-2)
    #cb2=plt.colorbar(pm, ax=ax2, pad=0.12, shrink=0.8, label="decibels (dB)")
    ax8, pm = wradlib.vis.plot_ppi(data=ZDRLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig1, ax=pylab.subplot(2,4,4), cmap=plt.get_cmap('viridis', 12), vmax=10, vmin=-2)
    cb8=plt.colorbar(pm, ax=ax8, pad=0.12, shrink=0.8, label="decibels (dB)")

    ax3, pm = wradlib.vis.plot_ppi(data=RhoHVLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig1, ax=pylab.subplot(2,4,5), norm=RhoHVnorm, cmap=plt.get_cmap(cm.managua, 10), vmax=1, vmin=0) # norm=colors.Normalise(vmin=0, vmin=1)
    #cb3=plt.colorbar(pm, ax=ax3, pad=0.12, shrink=0.8, label="unitless")
    ax9, pm = wradlib.vis.plot_ppi(data=RhoHVLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig1, ax=pylab.subplot(2,4,6), norm=RhoHVnorm, cmap=plt.get_cmap(cm.managua, 10), vmax=1, vmin=0) # norm=colors.LogNorm(vmin = 1e-10, vmax=1, clip=True) , replace vmin and vmax with this line # see also norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.8, vmax=1), cmap=cm.managua)
    cb9=plt.colorbar(pm, ax=ax9, pad=0.12, shrink=0.8, label="unitless")

    ax4, pm = wradlib.vis.plot_ppi(data=Presence_Lim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig1, ax=pylab.subplot(2,4,7), cmap=plt.get_cmap(cmap, 4), vmax=1.5, vmin=-2.5)
    #cb1=plt.colorbar(pm, ax=ax1, pad=0.12, shrink=0.8, ticks=np.arange(-3,7))
    ax10, pm = wradlib.vis.plot_ppi(data=Presence_Lim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation1, fig=fig1, ax=pylab.subplot(2,4,8), cmap=plt.get_cmap(cmap, 4), vmax=1.5, vmin=-2.5)
    cb10=plt.colorbar(pm, ax=ax10, pad=0.12, shrink=0.8, ticks=np.arange(-2,2))


    #text1=ax1.set_title("2° Elevation\nHorizontal Reflectivity Factor (dBZ)", fontsize=18)
    text1=ax1.set_title("2° Elevation", fontsize=18)
    text2=ax2.set_title("2° Elevation", fontsize=18)
    #text3=ax3.set_title("2° Elevation", fontsize=18)
    #text4=ax4.set_title("2° Elevation", fontsize=18)

    #text7=ax7.set_title("4.5° Elevation\nHorizontal Reflectivity Factor (dBZ)", fontsize=18)
    text7=ax7.set_title("4.5° Elevation", fontsize=18)
    text8=ax8.set_title("4.5° Elevation", fontsize=18)
    #text9=ax9.set_title("Correlation Coefficient (ρHV)", fontsize=18)

    ax1.set_ylabel('Polar Range Northward (km)')
    ax1.set_xlabel('Polar Range Eastward (km)')
    ax2.set_ylabel('Polar Range Northward (km)')
    ax2.set_xlabel('Polar Range Eastward (km)')
    ax3.set_ylabel('Polar Range Northward (km)')
    ax3.set_xlabel('Polar Range Eastward (km)')
    ax7.set_ylabel('Polar Range Northward (km)')
    ax7.set_xlabel('Polar Range Eastward (km)')
    ax8.set_ylabel('Polar Range Northward (km)')
    ax8.set_xlabel('Polar Range Eastward (km)')
    ax9.set_ylabel('Polar Range Northward (km)')
    ax9.set_xlabel('Polar Range Eastward (km)')

    #ax5.set_facecolor('black')
    #ax11.set_facecolor('black')


    cb7.ax.tick_params(labelsize=14)
    #cb7.set_label(label="decibels (dB)", size=18)
    cb7.set_label(label="Horizontal\nReflectivity Factor (dB)", size=18)
    cb8.ax.tick_params(labelsize=14)
    #cb8.set_label(label="decibels (dB)", size=18)
    cb8.set_label(label="Differential Reflectivity (dB)", size=18)
    cb9.ax.tick_params(labelsize=14)
    #cb9.set_label(label="unitless", size=18)
    cb9.set_label(label="Correlation Coefficient\n(unitless)", size=18)
    cb10.ax.set_yticklabels(np.flip(["Non-Meteorological","No Significant Signal", "Meteorological", "Radar Shadow"]))

    #text=fig1.suptitle("Overview of Raw Radar Variables on"+"\n"+"2017-07-05 at 15:46:57, "+strElevation1+"° & "+strElevation2+
    #"° Elevation PPIs", fontsize=24)

    fig1.tight_layout(h_pad=0.0,w_pad=0.0)
    fig1.subplots_adjust(wspace=0.0, hspace=0.25, left=0.15, right=0.85, top=0.92)

    fig1.savefig(fname=i[38:53]+"_RawAndKilambi.png",format="png")
    
    del guide
    del Range
    del range
    del ranges
    del Azimuth
    del azimuth
    del azimuths
    del Elevation
    del elevation
    del elevations
    del time
    del Time
    del times
    del dBZ
    del dBZs
    del ZDR
    del ZDRs
    del RhoHV
    del RhoHVs   
    del V
    del Vs
   
    del Longitude
    del Longitudes
    del Latitude
    del Latitudes
    del Altitude
    del Altitudes
    
    del dBZnanindex
    del ZDRnanindex
    del RhoHVnanindex
    del lZDR
    del dr
    del dr_isnan
    del dbz_isnan
    del Weather
    del Presence
    del Presences