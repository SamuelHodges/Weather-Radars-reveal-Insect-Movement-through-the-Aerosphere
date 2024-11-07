# plot biomod projection panels, per date, hour and altitude
#
# Run from the sf2 environment on SAM
#

#Args list:
#
# 1 - Directory with 'proj_*' directories
# 2 - Model name (based on .out file)
# 3 - (Optional) Algorithm (where multiple were used)

library(geodata)
library(sf)
library(ggplot2)

Args=commandArgs(trailingOnly=TRUE)

setwd(Args[1])
setwd("..")
toplevel=getwd() # set top level to the directory above the folder containing projections
setwd(Args[1])

allprojections=list.dirs()[grep(list.dirs(), pattern="proj_")]

splits=strsplit(allprojections, "_")

#get requested model-algorithm's dates
models=which( do.call(c, lapply(splits, FUN="[", 2)) == Args[2] )

if (nchar(Args[3])<8) {	#determine whether Args[3] is set as a date or an algorithm name
model_algos=which( do.call(c, lapply(splits, FUN="[", 2)) == Args[2] & do.call(c, lapply(splits, FUN="[", 3)) == Args[3] )

dates=unique(do.call(c, lapply(splits[model_algos], FUN="[", 4)))
#hours=unique(do.call(c, lapply(splits, FUN="[", 3)))
#alts=unique(do.call(c, lapply(splits, FUN="[", 4)))
} else {
dates=unique(do.call(c, lapply(splits[models], FUN="[", 3)))

}

# get world map data in a spatvector file, then convert to simple features
countries = geodata::world(path = ".")
map = st_as_sf(countries)

for (d in dates) {
	
	print(d)
	# ready base projections
	if (nchar(Args[3])<8) {	#determine whether Args[3] is set as a date or an algorithm name
	projections=list.files(pattern=paste("proj",Args[2],Args[3],d, sep="_"))
	} else {
	projections=list.files(pattern=paste("proj",Args[2],d, sep="_")) }
	
	# load projections and label with date, altitude and time 
	setwd(projections[1])
	panels_and_axes=readRDS("panels_and_axes.rds")
	lonlat=expand.grid(Longitude=panels_and_axes$Longitude, Latitude=panels_and_axes$Latitude)
	#lonlat=lonlat[order(lonlat$Longitude, lonlat$Latitude),]
	
	preds=read.csv("predictions.csv")
	#preds=preds[order(preds$points, decreasing=TRUE),]
	
	pred_base = cbind( Date=d, Time=panels_and_axes$Time, Altitude=panels_and_axes$Altitude,
	Longitude=lonlat$Longitude, Latitude=lonlat$Latitude, preds )
	
	setwd(Args[1])

	for (p in projections[-1]) {

		setwd(p)
		panels_and_axes=readRDS("panels_and_axes.rds")
		lonlat=expand.grid(Longitude=panels_and_axes$Longitude, Latitude=panels_and_axes$Latitude)
		#lonlat=lonlat[order(lonlat$Longitude, lonlat$Latitude),]
		
		preds=read.csv("predictions.csv")
		#preds=preds[order(preds$points, decreasing=TRUE),]

		pred_table = cbind( Date=d, Time=panels_and_axes$Time, Altitude=panels_and_axes$Altitude,
		Longitude=lonlat$Longitude, Latitude=lonlat$Latitude, preds )
		
		setwd(Args[1])
		
		pred_base=rbind(pred_base, pred_table)

	}
	
	#transform to probability scale and correct labels
	pred_base$pred=pred_base$pred/1000
	pred_base$Altitude=paste0(pred_base$Altitude, "m")
	
	pred_base$Time=paste0(pred_base$Time, ":00")
	pred_base$Time[which(nchar(pred_base$Time)==4)]=paste0("0",pred_base$Time[which(nchar(pred_base$Time)==4)])
	pred_base$Time=paste0(pred_base$Time, " UTC+0")
	
	#plot mapgrids faceted by Time and Altitude
	plot=ggplot()+geom_raster(aes(Longitude, Latitude, fill=pred), pred_base)+facet_grid(Altitude~Time)+
	scale_fill_continuous(limits=c(0,1), type="viridis", name="Probability\nof Presence")+geom_sf(colour="white", fill=NA, data=map)+
	xlim(c(min(pred_base$Longitude),max(pred_base$Longitude)))+ylim(c(min(pred_base$Latitude),max(pred_base$Latitude)))
	
	if (nchar(Args[3])<8) {	#determine whether Args[3] is set as a date or an algorithm name
	ggsave(paste(Args[2], Args[3], d, "Projection_Map.png", sep="_"), plot, width=11, height=8, scale=0.7) } else {
	ggsave(paste(Args[2], d, "Projection_Map.png", sep="_"), plot) }
	
	setwd(Args[1])
}