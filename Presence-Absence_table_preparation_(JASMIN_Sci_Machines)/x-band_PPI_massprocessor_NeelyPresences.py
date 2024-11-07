# X-band PPI mass processor
# args 1 - target directory
#      2 - target file (if parallelised variant)

# set up

import sys
import os
import datetime
import wradlib
import pyart
import numpy as np
from numpy import emath as e
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import pylab
from matplotlib import colors as c
import cartopy.crs as ccrs
import warnings
import glob
import cmcrameri.cm as cm

warnings.simplefilter("ignore")

os.environ['WRADLIB_DATA']=os.getcwd()
matplotlib.use('Agg')

# Look for cases where all three filtering stages are present
os.chdir(sys.argv[1])
#os.chdir("/gws/nopw/j04/ncas_radar_vol1/eesjir/transfer2_data_prep/data/x-band")

#files_xyz=glob.glob("*XYZ_raw.nc")
#files_rawpres=glob.glob("*XYZ_rawpresences.nc")
#files_mdspres=glob.glob("*XYZ_MDSpresences.nc")
#files_bipres=glob.glob("*XYZ_MDS_BirdInsect_presences.nc")

#loadindex=np.empty(len(files_xyz))
#entry=0
#for f in files_xyz:
#    if len(glob.glob("*"+f[38:53]+"*.nc"))<5:
#        loadindex[entry]=0
#    else:
#        loadindex[entry]=1
#    entry+=1

#loadindex=np.loadtxt(fname="PPI_loadindex.csv", delimiter=",")

# Begin loop

files_xyz=glob.glob(sys.argv[2]) #for a for loop variant (not parallelised), uncomment the above and comment out this line. 
files_dates=[f[38:53] for f in files_xyz]

#entry=0
for f in files_dates:
    #if loadindex[entry]==0:
    #    continue
    #entry+=1
    
    # Load files and key variables
    files=glob.glob("*"+f+"*.nc")
    #filesraw=glob.glob("*"+f+"*.nc")

    radar_rawpres=wradlib.io.read_generic_netcdf(files[0])
    radar_xyz=wradlib.io.read_generic_netcdf(files[0])

    Range=np.ma.getdata(radar_rawpres['variables']['range']['data'])
    Azimuth=np.ma.getdata(radar_rawpres['variables']['azimuth']['data'])
    Elevation=np.ma.getdata(radar_rawpres['variables']['elevation']['data'])
    dBZ=np.ma.getdata(radar_rawpres['variables']['dBZ']['data'])
    dBZv=np.ma.getdata(radar_rawpres['variables']['dBZ']['data'])

    sitecoords=tuple(np.ma.getdata(radar_rawpres['variables']['sitecoordinates']['data']))

    ZDR=np.ma.getdata(radar_rawpres['variables']['ZDR']['data'])
    RhoHV=np.ma.getdata(radar_rawpres['variables']['RhoHV']['data'])
    DR=np.ma.getdata(radar_rawpres['variables']['DR']['data'])

    V=np.ma.getdata(radar_rawpres['variables']['radial_velocity_h']['data'])
    #Airspeed=np.ma.getdata(radar_BirdInsect['variables']['Airspeed']['data'])

    Presence_rawpres=np.ma.getdata(radar_rawpres['variables']['Presence']['data'])
 
    #Set nans for raw observations
    dBZnanindex=np.argwhere(np.less_equal(np.ma.getdata(radar_rawpres['variables']['dBZ']['data']), -100000))
    ZDRnanindex=np.argwhere(np.less_equal(np.ma.getdata(radar_rawpres['variables']['ZDR']['data']), -100000))
    RhoHVnanindex=np.argwhere(np.isnan(np.ma.getdata(radar_rawpres['variables']['RhoHV']['data'])))
    DRnanindex=np.argwhere(np.isnan(np.ma.getdata(radar_rawpres['variables']['DR']['data'])))
    Vnanindex=np.argwhere(np.less_equal(np.ma.getdata(radar_rawpres['variables']['radial_velocity_h']['data']), -900))
    #Airspeednanindex=np.argwhere(np.greater_equal(np.ma.getdata(radar_BirdInsect['variables']['Airspeed']['data']), 100000))

    dBZ[dBZnanindex[0:dBZnanindex.shape[0],0], dBZnanindex[0:dBZnanindex.shape[0],1]]=np.nan
    ZDR[ZDRnanindex[0:ZDRnanindex.shape[0],0], ZDRnanindex[0:ZDRnanindex.shape[0],1]]=np.nan
    RhoHV[RhoHVnanindex[0:RhoHVnanindex.shape[0],0], RhoHVnanindex[0:RhoHVnanindex.shape[0],1]]=np.nan
    DR[ZDRnanindex[0:ZDRnanindex.shape[0],0], ZDRnanindex[0:ZDRnanindex.shape[0],1]]=np.nan
    V[Vnanindex[0:Vnanindex.shape[0],0], Vnanindex[0:Vnanindex.shape[0],1]]=np.nan
    #Airspeed[Airspeednanindex[0:Airspeednanindex.shape[0],0], Airspeednanindex[0:Airspeednanindex.shape[0],1]]=np.nan
    
    
    ### Temporary additional raw presence processing, remove once filters are fully updated ###
    
    ## Remove data from elevations below 2 degrees to exclude ground clutter
    elenan=np.argwhere(Elevation<2)
    
    dBZ[:,elenan[:,0]]=np.nan
    dBZv[:,elenan[:,0]]=np.nan
    ZDR[:,elenan[:,0]]=np.nan
    RhoHV[:,elenan[:,0]]=np.nan
    V[:,elenan[:,0]]=np.nan
    
    ## Remove Known Radar shadows, their edge effects, and the ZDR ground clutter cone (manually remove as part of Effective Range to avoid false positives IDs from Kilambi)
    
    ZDR_Cone=np.argwhere(np.logical_and( Azimuth > 175, Azimuth < 185 ))
    Radar_shadows1=np.argwhere(np.logical_and( Azimuth > 54, Azimuth < 65 ))
    Radar_shadows2=np.argwhere(np.logical_and( Azimuth > 82,  Azimuth < 91 ))
    Radar_shadows3=np.argwhere(np.logical_and( Azimuth > 128, Azimuth < 155 ))

    Mask=np.concatenate((ZDR_Cone,Radar_shadows1,Radar_shadows2,Radar_shadows3), axis=0)

    Presence_rawpres[:,Mask[:,0]] = -2
    Presence_rawpres[:,elenan[:,0]]= -2

    Presence_rawpres[np.logical_and((Presence_rawpres==-1),(dBZ!=np.nan),(dBZv==np.nan))]=-3 # set areas where vertical velocity has attenuated to false absence
    
    ###     ###
    
    
    rangeindicemax=np.count_nonzero(np.less_equal(Range, 80000))
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
    #AirspeedLim1=Airspeed[0:rangeindicemax, elevationindices1.flatten()]

    Presence_Lim1=Presence_rawpres[0:rangeindicemax, elevationindices1.flatten()]

    Azimuths2=Azimuth[elevationindices2.flatten()]

    dBZLim2=dBZ[0:rangeindicemax, elevationindices2.flatten()]

    ZDRLim2=ZDR[0:rangeindicemax, elevationindices2.flatten()]
    RhoHVLim2=RhoHV[0:rangeindicemax, elevationindices2.flatten()]
    DRLim2=DR[0:rangeindicemax, elevationindices2.flatten()]

    VLim2=V[0:rangeindicemax, elevationindices2.flatten()]
    #AirspeedLim2=Airspeed[0:rangeindicemax, elevationindices2.flatten()]

    Presence_Lim2=Presence_rawpres[0:rangeindicemax, elevationindices2.flatten()]
     
    
    # Figure 4 (widescale plot of raw and Kilambi-only filtering)
    
    fig4 = plt.figure(figsize=[26, 20])
    
    title=files[0][38:53]
    fig4_title="Raw Radar Variables and Kilambi Filter - "+title[0:4]+"/"+title[4:6]+"/"+title[6:8]+" "+title[9:11]+":"+title[11:13]+":"+title[13:15]
    fig4.suptitle(fig4_title, size=28, weight=5)
    
    Vcmap=c.LinearSegmentedColormap.from_list(name='Vcmap', colors=np.flip(['red', 'black', 'blue']))

    cmap=c.ListedColormap(np.flip(['limegreen', 'black', 'mediumblue', 'whitesmoke','firebrick']))
    RhoHVnorm=c.TwoSlopeNorm(vcenter=0.8, vmin=0, vmax=1)

    ax1, pm = wradlib.vis.plot_ppi(data=dBZLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,1), cmap=plt.get_cmap('pyart_HomeyerRainbow', 20), vmax=20, vmin=-20)
    #cb1=plt.colorbar(pm, ax=ax1, pad=0.12, shrink=0.8, label="decibels")
    ax7, pm = wradlib.vis.plot_ppi(data=dBZLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig4, ax=pylab.subplot(2,4,2), cmap=plt.get_cmap('pyart_HomeyerRainbow', 20), vmax=20, vmin=-20)
    cb7=plt.colorbar(pm, ax=ax7, pad=0.12, shrink=0.8, label="decibels (dB)")

    ax2, pm = wradlib.vis.plot_ppi(data=ZDRLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,3), cmap=plt.get_cmap('viridis', 12), vmax=10, vmin=-2)
    #cb2=plt.colorbar(pm, ax=ax2, pad=0.12, shrink=0.8, label="decibels (dB)")
    ax8, pm = wradlib.vis.plot_ppi(data=ZDRLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig4, ax=pylab.subplot(2,4,4), cmap=plt.get_cmap('viridis', 12), vmax=10, vmin=-2)
    cb8=plt.colorbar(pm, ax=ax8, pad=0.12, shrink=0.8, label="decibels (dB)")

    ax3, pm = wradlib.vis.plot_ppi(data=RhoHVLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,5), cmap=plt.get_cmap('pyart_LangRainbow12', 10), vmax=1, vmin=0) # norm=colors.Normalise(vmin=0, vmin=1)
    #cb3=plt.colorbar(pm, ax=ax3, pad=0.12, shrink=0.8, label="unitless")
    ax9, pm = wradlib.vis.plot_ppi(data=RhoHVLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig4, ax=pylab.subplot(2,4,6), cmap=plt.get_cmap('pyart_LangRainbow12', 10), vmax=1, vmin=0) # norm=colors.LogNorm(vmin = 1e-10, vmax=1, clip=True) , replace vmin and vmax with this line # see also norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.8, vmax=1), cmap=cm.managua)
    cb9=plt.colorbar(pm, ax=ax9, pad=0.12, shrink=0.8, label="unitless")

    ax4, pm = wradlib.vis.plot_ppi(data=Presence_Lim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,7), cmap=plt.get_cmap(cmap, 4), vmax=1.5, vmin=-3.5)
    #cb1=plt.colorbar(pm, ax=ax1, pad=0.12, shrink=0.8, ticks=np.arange(-3,7))
    ax10, pm = wradlib.vis.plot_ppi(data=Presence_Lim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,8), cmap=plt.get_cmap(cmap, 4), vmax=1.5, vmin=-3.5)
    cb10=plt.colorbar(pm, ax=ax10, pad=0.12, shrink=0.8, ticks=np.arange(-3,2))


    #text1=ax1.set_title("2° Elevation\nHorizontal Reflectivity Factor (dBZ)", fontsize=18)
    text1=ax1.set_title("2° Elevation", fontsize=18)
    text2=ax2.set_title("2° Elevation", fontsize=18)

    #text7=ax7.set_title("4.5° Elevation\nHorizontal Reflectivity Factor (dBZ)", fontsize=18)
    text7=ax7.set_title("4.5° Elevation", fontsize=18)
    text8=ax8.set_title("4.5° Elevation", fontsize=18)

    ax1.set_ylabel('Polar Range Northward (km)')
    ax1.set_xlabel('Polar Range Eastward (km)')
    ax2.set_ylabel('Polar Range Northward (km)')
    ax2.set_xlabel('Polar Range Eastward (km)')
    ax3.set_ylabel('Polar Range Northward (km)')
    ax3.set_xlabel('Polar Range Eastward (km)')
    ax4.set_ylabel('Polar Range Northward (km)')
    ax4.set_xlabel('Polar Range Eastward (km)')
    ax7.set_ylabel('Polar Range Northward (km)')
    ax7.set_xlabel('Polar Range Eastward (km)')
    ax8.set_ylabel('Polar Range Northward (km)')
    ax8.set_xlabel('Polar Range Eastward (km)')
    ax9.set_ylabel('Polar Range Northward (km)')
    ax9.set_xlabel('Polar Range Eastward (km)')
    ax10.set_ylabel('Polar Range Northward (km)')
    ax10.set_xlabel('Polar Range Eastward (km)')

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

    # set axis labels (1 - original ticks, 2 - reduced ticks (no Nussbaumer filtering), 3 - tick labels with counts of radar sectors affiliated
    # 4 - tick labels with percents of radar sectors affiliated)
    
    cb10.ax.set_yticklabels(np.flip(["Non-Meteorological","No Significant Signal", "Meteorological", "Radar Shadow", "Indeterminate Scatter"]))
    #cb10.ax.set_yticklabels(np.flip(["","Non-Meteorological","No Significant Signal", "Meteorological", ""]))
    #cb10.ax.set_yticklabels(np.flip(["","Non-Meteorological"+"\nTotals (left) = "+str(nNonMeta)+"\nTotals (right) = "+str(nNonMetb),"No Significant Signal"+"\nTotals (left) = "+str(nNoSiga)+"\nTotals (right) = "+str(nNoSigb), "Meteorological"+"\nTotals (left) = "+str(nMeta)+"\nTotals (right) = "+str(nMetb), ""]))
    #cb10.ax.set_yticklabels(np.flip(["","Non-Meteorological"+"\nTotals (left) = "+str(np.around(nNonMeta/np.size(Presence_rawpresLim1)*100,1))+"%\nTotals (right) = "+str(np.around(nNonMetb/np.size(Presence_rawpresLim2)*100,1))+"%","No Significant Signal"+"\nTotals (left) = "+str(np.around(nNoSiga/np.size(Presence_rawpresLim1)*100,1))+"%\nTotals (right) = "+str(np.around(nNoSigb/np.size(Presence_rawpresLim1)*100,1))+"%", "Meteorological"+"\nTotals (left) = "+str(np.around(nMeta/np.size(Presence_rawpresLim1)*100,1))+"%\nTotals (right) = "+str(np.around(nMetb/np.size(Presence_rawpresLim1)*100,1))+"%", ""]))

    #text=fig4.suptitle("Overview of Raw Radar Variables on"+"\n"+"2017-07-05 at 15:46:57, "+strElevation1+"° & "+strElevation2+
    #"° Elevation PPIs", fontsize=24)

    fig4 = plt.figure(figsize=[48, 14])

    title=f
    fig4_title="Raw Radar Variables and Kilambi Classification - "+title[0:4]+"/"+title[4:6]+"/"+title[6:8]+" "+title[9:11]+":"+title[11:13]+":"+title[13:15]
    fig4.suptitle(fig4_title, size=28, weight=5)

    cmap=c.ListedColormap(np.flip(['limegreen', 'black', 'mediumblue', 'whitesmoke','firebrick']))
    RhoHVnorm=c.TwoSlopeNorm(vcenter=0.8, vmin=0, vmax=1)

    ax1, pm = wradlib.vis.plot_ppi(data=dBZLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,1), cmap=plt.get_cmap('pyart_HomeyerRainbow', 20), vmax=20, vmin=-20)
    #cb1=plt.colorbar(pm, ax=ax1, pad=0.12, shrink=0.8, label="decibels")
    ax7, pm = wradlib.vis.plot_ppi(data=dBZLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig4, ax=pylab.subplot(2,4,2), cmap=plt.get_cmap('pyart_HomeyerRainbow', 20), vmax=20, vmin=-20)
    cb7=plt.colorbar(pm, ax=ax7, pad=0.12, shrink=0.8, label="decibels (dB)")

    ax2, pm = wradlib.vis.plot_ppi(data=ZDRLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,3), cmap=plt.get_cmap('viridis', 12), vmax=10, vmin=-2)
    #cb2=plt.colorbar(pm, ax=ax2, pad=0.12, shrink=0.8, label="decibels (dB)")
    ax8, pm = wradlib.vis.plot_ppi(data=ZDRLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig4, ax=pylab.subplot(2,4,4), cmap=plt.get_cmap('viridis', 12), vmax=10, vmin=-2)
    cb8=plt.colorbar(pm, ax=ax8, pad=0.12, shrink=0.8, label="decibels (dB)")

    ax3, pm = wradlib.vis.plot_ppi(data=RhoHVLim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,5), cmap=plt.get_cmap('pyart_LangRainbow12', 10), vmax=1, vmin=0) # norm=colors.Normalise(vmin=0, vmin=1)
    #cb3=plt.colorbar(pm, ax=ax3, pad=0.12, shrink=0.8, label="unitless")
    ax9, pm = wradlib.vis.plot_ppi(data=RhoHVLim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation2, fig=fig4, ax=pylab.subplot(2,4,6), cmap=plt.get_cmap('pyart_LangRainbow12', 10), vmax=1, vmin=0) # norm=colors.LogNorm(vmin = 1e-10, vmax=1, clip=True) , replace vmin and vmax with this line # see also norm = colors.TwoSlopeNorm(vmin=0, vcenter=0.8, vmax=1), cmap=cm.managua)
    cb9=plt.colorbar(pm, ax=ax9, pad=0.12, shrink=0.8, label="unitless")

    ax4, pm = wradlib.vis.plot_ppi(data=Presence_Lim1.transpose(), r=RangesLim, az=Azimuths1, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,7), cmap=plt.get_cmap(cmap, 4), vmax=1.5, vmin=-3.5)
    #cb1=plt.colorbar(pm, ax=ax1, pad=0.12, shrink=0.8, ticks=np.arange(-3,7))
    ax10, pm = wradlib.vis.plot_ppi(data=Presence_Lim2.transpose(), r=RangesLim, az=Azimuths2, elev=ScanElevation1, fig=fig4, ax=pylab.subplot(2,4,8), cmap=plt.get_cmap(cmap, 4), vmax=1.5, vmin=-3.5)
    cb10=plt.colorbar(pm, ax=ax10, pad=0.12, shrink=0.8, ticks=np.arange(-3,2))


    #text1=ax1.set_title("2° Elevation\nHorizontal Reflectivity Factor (dBZ)", fontsize=18)
    text1=ax1.set_title(strElevation1+"° Elevation", fontsize=18)
    text2=ax2.set_title(strElevation1+"° Elevation", fontsize=18)
    #text3=ax3.set_title("2° Elevation", fontsize=18)
    #text4=ax4.set_title("2° Elevation", fontsize=18)

    #text7=ax7.set_title("4.5° Elevation\nHorizontal Reflectivity Factor (dBZ)", fontsize=18)
    text7=ax7.set_title(strElevation2+"° Elevation", fontsize=18)
    text8=ax8.set_title(strElevation2+"° Elevation", fontsize=18)
    #text9=ax9.set_title("Correlation Coefficient (ρHV)", fontsize=18)

    ax1.set_ylabel('Polar Range Northward (km)')
    ax1.set_xlabel('Polar Range Eastward (km)')
    ax2.set_ylabel('Polar Range Northward (km)')
    ax2.set_xlabel('Polar Range Eastward (km)')
    ax3.set_ylabel('Polar Range Northward (km)')
    ax3.set_xlabel('Polar Range Eastward (km)')
    ax4.set_ylabel('Polar Range Northward (km)')
    ax4.set_xlabel('Polar Range Eastward (km)')
    ax7.set_ylabel('Polar Range Northward (km)')
    ax7.set_xlabel('Polar Range Eastward (km)')
    ax8.set_ylabel('Polar Range Northward (km)')
    ax8.set_xlabel('Polar Range Eastward (km)')
    ax9.set_ylabel('Polar Range Northward (km)')
    ax9.set_xlabel('Polar Range Eastward (km)')
    ax10.set_ylabel('Polar Range Northward (km)')
    ax10.set_xlabel('Polar Range Eastward (km)')

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
    #cb9.ax.set_yticklabels(["0.00","0.20","0.40","0.60","0.80","0.85","0.90","0.95","1.00"])
    cb10.ax.set_yticklabels(np.flip(["Non-Meteorological","No Significant Signal", "Meteorological", "Radar Shadow", "Indeterminate Scatter"]))

    #labels=[str(item) for item in np.arange(0,1.1,0.1)]
    #pos=[item[0:3:1] for item in labels]
    #cb9.ax.set_yticks(np.arange(-0.1,1.0,0.1))
    #cb9.ax.set_yticklabels(pos)

    #text=fig4.suptitle("Overview of Raw Radar Variables on"+"\n"+"2017-07-05 at 15:46:57, "+strElevation1+"° & "+strElevation2+
    #"° Elevation PPIs", fontsize=24)

    fig4.tight_layout(h_pad=0.0,w_pad=0.0)
    fig4.subplots_adjust(wspace=0.0, hspace=0.25, left=0.15, right=0.85, top=0.92)

    print(fig4)
    # Print and label radar figures
    
    fig4.savefig(fname=f+"_RawAndKilambi.png", format="png")
    
    print(f+" plots completed at time:",datetime.datetime.now())

