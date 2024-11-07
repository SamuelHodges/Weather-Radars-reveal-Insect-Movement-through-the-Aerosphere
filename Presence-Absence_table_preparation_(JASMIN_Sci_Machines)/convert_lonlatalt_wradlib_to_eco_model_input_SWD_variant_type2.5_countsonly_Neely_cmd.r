# args
# 1 - date in yyyymmdd or 'all', to iterate through directory
# 2 - start hour in 24hr clock
# 3 - end hour in 24hr clock
# no target directory argument, run at top level for file structure (biomod2_transfer_results)

Args=commandArgs(trailingOnly=TRUE)

#set date & working directory

if (Args[1]=="all") {
	date=c()
	#check for which dates x-band and ecmwf files both exist
	xband_check=names(table(substr(list.files(path="./data/x-band/NeelyFilterTest/NeelyFinal", pattern="NeelyPresences.nc"), start=39, stop=46)))
	ecmwf_check=gsub("-", "", names(table(substr(list.files(path="./data/ecmwf-operational-archive", pattern=".nc"), start=11, stop=20))))
	for (xband_day in xband_check) {
		day=ecmwf_check[which(ecmwf_check==xband_day)]
		date=c(date, day) }
	if (length(which(is.na(date)==TRUE))>0) { date=as.numeric(date[-which(is.na(date)==TRUE)]) }
} else { date=as.numeric(Args[1]) }

#setwd("/home/users/eesjir/eesjir_group_workspace/biomod2_transfer_results")

#load packages
library(ncdf4)
library(sp)

#loop to load relevant files

hour=seq(as.numeric(Args[2]),as.numeric(Args[3]),1)

for (day in date) {
	
	daycode=paste("OP-sam_fc-", paste(substr(day, start=1, stop=4), substr(day, start=5, stop=6), substr(day, start=7, stop=8), sep="-"), sep="")
	
	if (length(list.files(pattern=paste("Bioscatter_biomod", daycode, "table.csv", sep="_")))>0) { print(paste("Output for", daycode, "already exists, skipping")) ; next
	} else { print(day) }

	ecmwf_file=nc_open(list.files(pattern=daycode, recursive=TRUE))

	#expand ecmwf
	lon=ncvar_get(ecmwf_file, "lon")
	lat=ncvar_get(ecmwf_file, "lat")
	plev=ncvar_get(ecmwf_file, "plev")
	time=ncvar_get(ecmwf_file, "time")
	#temp=ncvar_get(ecmwf_file, "t")-273.15
	#rh=ncvar_get(ecmwf_file, "r")
	#u_wind=ncvar_get(ecmwf_file, "u")
	#v_wind=ncvar_get(ecmwf_file, "v")
	#windspeed=sqrt(u_wind^2+v_wind^2)
	
	#convert pressure levels to altitude
	alt=plev
	alt[which(plev==100000)]=0
	alt[which(plev==90000)]=1000
	alt[which(plev==80000)]=2000
	alt[which(plev==70000)]=3000
	
	#correct time
	#time=seq(0, length(time), 1)

	date_files=list.files(pattern=as.character(day), recursive=TRUE)
	xbands=date_files[grep(x=date_files, pattern="NeelyPresences.nc")]

	dir_offset_radar=lapply(gregexpr(xbands, pattern="/"), max) # calculate folder offset for strings
	dir_offset_radar=dir_offset_radar[[1]]
	
	#initialise bioscatter table
	bioscatter_full_csv=data.frame(Species=character(), Presence.Absence=logical(0), PresenceCount=numeric(0), TrueAbsenceCount=numeric(0), Hour=character(),
	env_Longitude=numeric(0), env_Latitude=numeric(0), env_Altitude=numeric(0), Temperature=numeric(0), Relative_Humidity=numeric(0), Windspeed=numeric(0))

	hibioscatter_full_csv=data.frame(Longitude=numeric(0), Latitude=numeric(0), Altitude=numeric(0), Presence.Absence=logical(0)) # spare dataframe for values above altitude limit

	#loadflag=0

	# begin hour loop
	for (i in hour) {
		#get files
		#if (i>11 & loadflag==0) { ecmwf_file=nc_open(list.files(pattern=paste("ECMWF_UK", day, "1200-2300", sep="_"), recursive=TRUE)) ; loadflag=1}

		if (nchar(as.character(i))==1) { select=paste("0",i,sep="") } else { select=i }
		#if (which(i==hour)>length(xbands)) { print("last available hour indexed") ; break }
		if (length(which(substr(xbands, start=48+dir_offset_radar, stop=49+dir_offset_radar)==select))==0) { print(paste("No data for hour", select, "on date", day, sep=" ")) ; next }
		
		print(paste("Hour", i))
		
		xband_file=nc_open(xbands[which(substr(xbands, start=48+dir_offset_radar, stop=49+dir_offset_radar)==select)])


		#if (loadflag==1) {
		#	#reload ecmwf from new file
		#	lon=ncvar_get(ecmwf_file, "lon")
		#	lat=ncvar_get(ecmwf_file, "lat")
		#	plev=ncvar_get(ecmwf_file, "plev")
		#	time=ncvar_get(ecmwf_file, "time")
		#	temp=ncvar_get(ecmwf_file, "var130")-273.15
		#	rh=ncvar_get(ecmwf_file, "var157")
		#	u_wind=ncvar_get(ecmwf_file, "var131")
		#	v_wind=ncvar_get(ecmwf_file, "var132")
		#	windspeed=sqrt(u_wind^2+v_wind^2)
		#	
		#	#convert pressure levels to altitude
		#	alt=plev
		#	alt[which(plev==100000)]=0
		#	alt[which(plev==90000)]=1000
		#	alt[which(plev==80000)]=2000
		#	alt[which(plev==70000)]=3000
		#	
		#	#correct time
		#	time=seq(0, length(time), 1)
		#	
		#	loadflag==2 }
		
		
		#expand x-band data
		dbz=as.vector(ncvar_get(xband_file, "dBZ"))
		#zdr=as.vector(ncvar_get(xband_file, "ZDR"))
		#rhohv=as.vector(ncvar_get(xband_file, "RhoHV"))

		#lzdr=10^(zdr/10)
		#dr = 10 * log10( ( lzdr + 1 - 2 * sqrt(lzdr) * rhohv ) / ( lzdr + 1 + 2 * sqrt(lzdr) * rhohv ) )

		#Presence.Absence=dr > -12 & dbz < 35
		Presence.Absence=as.vector(ncvar_get(xband_file, "Presence"))

		radar_Longitude=as.vector(ncvar_get(xband_file, "Longitude"))
		radar_Latitude=as.vector(ncvar_get(xband_file, "Latitude"))
		radar_Altitude=as.vector(ncvar_get(xband_file, "Altitude"))
		
		#remove out-of-range and weather absences
		
		radar_Longitude=radar_Longitude[-which(is.na(Presence.Absence)==TRUE | Presence.Absence==-1)]
		radar_Latitude=radar_Latitude[-which(is.na(Presence.Absence)==TRUE | Presence.Absence==-1)]
		radar_Altitude=radar_Altitude[-which(is.na(Presence.Absence)==TRUE | Presence.Absence==-1)]

		Presence.Absence=Presence.Absence[-which(is.na(Presence.Absence)==TRUE | Presence.Absence==-1)]
		
		#set aside values above env_altitude limit

		high_scatter=which(radar_Altitude>3500)
		
		hibioscatter_csv=data.frame(Longitude=radar_Longitude[high_scatter], Latitude=radar_Latitude[high_scatter], Altitude=radar_Altitude[high_scatter], Presence.Absence=Presence.Absence[high_scatter])
		hibioscatter_full_csv=rbind(hibioscatter_full_csv, hibioscatter_csv)
		
		if (length(high_scatter) != 0) {
			radar_Longitude=radar_Longitude[-high_scatter]
			radar_Latitude=radar_Latitude[-high_scatter]
			radar_Altitude=radar_Altitude[-high_scatter]

			Presence.Absence=Presence.Absence[-high_scatter]
		}
		
		#subset env data to Chilbolton area (radar data area) at hour
		#if (i<12) { hourindex=which(time==i) } else { hourindex=which(time==i-12) }
		hourindex=which(time==i)
		
		chillon=which(lon > min(radar_Longitude) & lon < max(radar_Longitude))	# note that this defines the IFS size by absence extent, see above lines 132-134
		chillat=which(lat > min(radar_Latitude) & lat < max(radar_Latitude))
		
		#temp1=temp[chillon, chillat, , hourindex]
		#rh1=rh[chillon, chillat, , hourindex]
		#windspeed1=windspeed[chillon, chillat, , hourindex]


		#build initial presence/absence csv
		lonlatalt=expand.grid(lon[chillon], lat[chillat], alt)

		#bioscatter_csv=data.frame(Species="Bioscatter", Presence.Absence=numeric(1), PresenceCount=numeric(1), TrueAbsenceCount=numeric(1), Hour=select,
		#env_Longitude=as.numeric(lonlatalt[,1]), env_Latitude=as.numeric(lonlatalt[,2]), env_Altitude=as.numeric(lonlatalt[,3]), 
		#Temperature=as.numeric(temp1), Relative_Humidity=as.numeric(rh1), Windspeed=as.numeric(windspeed1))
		
		bioscatter_csv=data.frame(Species="Bioscatter", Presence.Absence=numeric(1), PresenceCount=numeric(1), TrueAbsenceCount=numeric(1), Hour=select,
		env_Longitude=as.numeric(lonlatalt[,1]), env_Latitude=as.numeric(lonlatalt[,2]), env_Altitude=as.numeric(lonlatalt[,3]))

		namelist=names(bioscatter_csv)

		for (n in 1:length(ecmwf_file$var)){

			varname=ecmwf_file$var[[n]]$name
			varvalues=ncvar_get(ecmwf_file, varname)
			if ( length(dim(varvalues))==4 ) { varsubset=varvalues[chillon, chillat, , hourindex] ; varcolumn=as.numeric(varsubset)}
			
			if ( length(dim(varvalues))==3 ) { # if varvalues does not contain an altitude dimension (surface), expand it to include this dimension
				varsubset=varvalues[chillon, chillat, hourindex] 
				varcolumn=as.numeric(varsubset)
				varsubset=rep(varsubset, length(alt))}

			bioscatter_csv=cbind(bioscatter_csv, varcolumn)
			namelist=c(namelist, varname)

			next
			}
		
		names(bioscatter_csv)=namelist
		
		#assign data from x-band radar voxels to ecmwf raster layers
		
		# get longitude, latitude and altitude gates
		londiff=abs(mean(lon[2:length(lon)]-lon[1:(length(lon)-1)])/2)
		latdiff=abs(mean(lat[2:length(lat)]-lat[1:(length(lat)-1)])/2)
		altdiff=abs(mean(alt[2:length(alt)]-alt[1:(length(alt)-1)])/2)

		# create a 'holding pen' with a fill value, to hold selected data to avoid having to subset iteratively. nrow is set by radar data array
		#Presence.Absence.pen=rep(666, length(dbz))
		
		for (j in 1:nrow(bioscatter_csv)){
			
			envlon=bioscatter_csv$env_Longitude[j]
			envlat=bioscatter_csv$env_Latitude[j]
			envalt=bioscatter_csv$env_Altitude[j]
			
			Presence.Absence.bag=Presence.Absence[which( radar_Longitude < (envlon+londiff) & radar_Longitude > (envlon-londiff) & 
									radar_Latitude < (envlat+latdiff) & radar_Latitude > (envlat-latdiff) & 
									radar_Altitude < (envalt+altdiff) & radar_Altitude > (envalt-altdiff) )]
			
			bioscatter_csv$PresenceCount[j]=length(which(Presence.Absence.bag==1))
			bioscatter_csv$TrueAbsenceCount[j]=length(which(Presence.Absence.bag==0))
			if( length(which(Presence.Absence.bag==1)) > 0 ) { bioscatter_csv$Presence.Absence[j] = 1 }
			if( ( length(which(Presence.Absence.bag==1)) + length(which(Presence.Absence.bag==0)) ) < 1 ) { bioscatter_csv$Presence.Absence[j] = NA }
			
			if (j/500==as.integer(j/500)) {print(paste(j,"/",nrow(bioscatter_csv), " ", Sys.time(), sep=""))}
			
			next }

		bioscatter_full_csv=rbind(bioscatter_full_csv, bioscatter_csv)
		
	}

	
	#write both files
	write.csv(bioscatter_full_csv, paste("Bioscatter_biomod", day, "binned.csv", sep="_"))
	#write.csv(hibioscatter_full_csv, paste("Bioscatter_biomod", day, "high_altitude.csv", sep="_"))
}
