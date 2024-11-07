# projection script, using median accuracy model. Also saves a simple plot.
#
# Run from the r4 environment on SAM
#
# Args:
# 1 - Model directory (full path)
# 2 - Model name
# 3 - Atmospheric netcdf4 data (ECMWF)
# 4 - (Optional) Altitude levels to use, given as a vector index, separated by ','
# 5 - (Optional) Hours (index) to use in projection, separated by ','
# 6 - (Optional) algorithm argument to pass to get_built_models
# for optional arguments, set the position to 'NA'
# for particular days, make a sub directory containing the desired days. These should all have long,lat,alt and time dimensions!

library(raster)
library(ncdf4)
library(biomod2)

# Model directory

Args=commandArgs(trailingOnly=TRUE)

setwd(Args[1])
setwd("..")
toplevel=getwd() # save level code is run from to keep setwd() calls consistent (this is the level biomod calls must be run from)

# load in model and selection of Models
biomod_projection_dir=Args[1]
model_name=Args[2]

setwd(biomod_projection_dir)
models = get(load(model_name))
setwd(toplevel)

# log location of CEH pre-gridded data (currently directory is hardcoded, prepare these tables using the pairing procedure in add_CEH_land_cover_percents)
CEHdir="/data/transfer2_biomod/CEH-land-cover"
CEHtabs=list.files(CEHdir, pattern=".asc")

# determine the run with the median evaluation score and prepare it for biomod2 loading
if(is.na(Args[6])==TRUE) { evaluations = get_evaluations(models) } else { evaluations = get_evaluations(models, algo=Args[6]) }
TSS_evaluations = subset(evaluations, metric.eval=="TSS")

median_dists = abs(TSS_evaluations$validation - median(TSS_evaluations$validation)) 

median_model_name = evaluations$full.name[order(median_dists)][1]
median_model = get_built_models(models, full.name = median_model_name)

# loop the projection, to avoid building large dataframes in memory.

# prepare selection variables
atmo_dir=Args[3]
alts=as.integer(unlist(strsplit(Args[4], ",")))
hours=as.integer(unlist(strsplit(Args[5], ",")))

setwd(atmo_dir)
days=list.files(pattern=".nc")

for (d in days) {
	setwd(atmo_dir)
	atmo_nc=nc_open(d) # get atmo data then revert to toplevel of Biomod repo
	setwd(toplevel)
	#setwd(biomod_projection_dir)
	#setwd("..")
	
	for (h in hours) {
	
		for (a in alts) {
		
		atmo_frame=expand.grid(Longitude=ncvar_get(atmo_nc, "lon"), Latitude=ncvar_get(atmo_nc, "lat"), env_Altitude=ncvar_get(atmo_nc, "plev")[a], Hour=ncvar_get(atmo_nc, "time")[h])
		#Altitude=ncvar_get(atmo_nc, "plev")[a], Hours=ncvar_get(atmo_nc, "time")[h]) #set/reset atmo_frame to take new data
		
		#convert altitude layers from pressure levels to altitude, using MetPy's conversion under standard atmosphere as a standard.
		if(ncvar_get(atmo_nc, "plev")[a]==100000) { atmo_frame$env_Altitude=0 }
		if(ncvar_get(atmo_nc, "plev")[a]==90000) { atmo_frame$env_Altitude=1000 }
		if(ncvar_get(atmo_nc, "plev")[a]==80000) { atmo_frame$env_Altitude=2000 }
		if(ncvar_get(atmo_nc, "plev")[a]==70000) { atmo_frame$env_Altitude=3000 }

		
			for (v in names(atmo_nc$var)) {
			variable=ncvar_get(atmo_nc, v)
			
			if (length(dim(variable))==4) {  atmo_frame=cbind(atmo_frame, as.vector(variable[,,a,h]) ) }
			if (length(dim(variable))==3) {  atmo_frame=cbind(atmo_frame, as.vector(variable[,,h]) ) }
			}
		
		names(atmo_frame)=c("Longitude", "Latitude", "env_Altitude", "Hour", names(atmo_nc$var))
		
		#unpack and bind CEH data to projection frame
		lon=ncvar_get(atmo_nc, "lon")
		lat=ncvar_get(atmo_nc, "lat")
		londiff=abs(mean(lon[2:length(lon)]-lon[1:(length(lon)-1)])/2)
		latdiff=abs(mean(lat[2:length(lat)]-lat[1:(length(lat)-1)])/2)
		originalncol=ncol(atmo_frame)

		#make empty space in dataframe to fill with CEH data
		dummy=matrix(c(0), ncol=length(CEHtabs), nrow=nrow(atmo_frame))
		atmo_frame=cbind(atmo_frame, as.data.frame(dummy))
		names(atmo_frame)[(originalncol+1):(originalncol+length(CEHtabs))]=substr(CEHtabs, start=23, stop=nchar(CEHtabs)-4)

		setwd(CEHdir)
		for (C in 1:length(CEHtabs)){
			
			print(CEHtabs[C])
			
			CEH=raster(CEHtabs[C])
			CEH=as.data.frame(CEH, xy=TRUE)
			
			CEH_varname=names(CEH)[3]
			CEH_varname=substr(CEH_varname, start=23, stop=nchar(CEH_varname))
					
			for (j in 1:nrow(atmo_frame)){
				
				envlon=atmo_frame$Longitude[j]
				envlat=atmo_frame$Latitude[j]
				
				CEH.bag=CEH[,3][which( CEH$x < (envlon+londiff) & CEH$x > (envlon-londiff) & 
										CEH$y < (envlat+latdiff) & CEH$y > (envlat-latdiff))]
				
				atmo_frame[,(originalncol+C)][j]=mean(CEH.bag)
				
				if (j/500==as.integer(j/500)) {print(paste(j,"/",nrow(atmo_frame), " ", Sys.time(), sep=""))}
				
				next }
			
			}
			
		#rename band variables
		names(atmo_frame)[which(names(atmo_frame)=="Band1")]="Perc_Broadleaf_Wood"
		names(atmo_frame)[which(names(atmo_frame)=="Band10")]="Perc_Urban.Suburban"
		names(atmo_frame)[which(names(atmo_frame)=="Band2")]="Perc_Coniferous_Wood"
		names(atmo_frame)[which(names(atmo_frame)=="Band3")]="Perc_Arable"
		names(atmo_frame)[which(names(atmo_frame)=="Band5")]="Perc_Grassland"

		setwd(toplevel)
		
		#set projection directory name, based on model name, the get_built_models call and the date, time and altitude (yyyymmdd_hh_aaaa)
		mname=unlist(strsplit(model_name, "\\."))
		mname=mname[length(mname)-2]
		
		if(is.na(Args[6])==FALSE) { mname=paste(mname, Args[6], sep="_") }
		
		dname=unlist(strsplit(d, "")) # split date file string into single characters
		dname=dname[which(is.na(as.integer(dname))==FALSE)] # find all numbers inside the string
		dname=paste(dname, collapse="")
		
		projection_dir=paste(mname, dname, paste0(atmo_frame$Hour[1], "hrs"), paste0(atmo_frame$env_Altitude[1], "metres"), sep="_")
		
		#check for variable names beginning with an integer, and add an 'X' to the front, in keeping with Biomod2
		numbernames = which(is.na(as.integer(substr(names(atmo_frame), start=1, stop=1)))==FALSE)
		names(atmo_frame)[ numbernames ] = paste0( "X", names(atmo_frame)[ numbernames ] ) 
		
		projection = BIOMOD_Projection(bm.mod=models, proj.name=projection_dir, new.env=atmo_frame[,c(-1,-2)],
		new.env.xy=atmo_frame[,c(1,2)], models.chosen=median_model, build.clamping.mask=FALSE)
		
		# write simple plot image into projection directory
		setwd( paste0(biomod_projection_dir, "/proj_", projection_dir ) )
		png(width=960, height=960)
		plot(projection)
		dev.off()
		
		# write file of plot dimensions for later panelling
		panels_and_axes=list("Longitude"=ncvar_get(atmo_nc, "lon"), "Latitude"=ncvar_get(atmo_nc, "lat"),
		"Altitude"=unique(atmo_frame$env_Altitude), "Time"=unique(atmo_frame$Hour))
		
		saveRDS(panels_and_axes, "panels_and_axes.rds")
		
		# write data frame of predictions for loading into sf environment, during plotting
		preds=get_predictions(projection)
		write.csv(preds, "predictions.csv")
		
		setwd(toplevel)
		
		}
	}
}


