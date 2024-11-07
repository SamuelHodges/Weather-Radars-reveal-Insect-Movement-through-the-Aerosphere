# Create full ENM data tables

# Place script in top level of biomod repo folders
# 1 - Sub-Directory to target

Args=commandArgs(trailingOnly=TRUE)

library(biomod2)
library(stringr)

toplevel=paste(getwd(), Args[1], sep="/")
setwd(toplevel)


# list stored model data

files=list.dirs(getwd())
slashcount=str_count(files, "/")
files=files[which(slashcount==4)]
#files=files[grep(files, pattern="models")]

models=lapply(files, list.files)


# prepare path and .out information
modelsfullpath=character(0)
#for (i in 1:length(models)){
for (i in 1:length(files)) {
	modelsfullpathi=rep(files[i], length(models[[i]])-1)
	#modelsfullpathi=paste(files[i], models[[i]], sep="/")
	modelsfullpath=c(modelsfullpath, modelsfullpathi)
	models[[i]]=models[[i]][-length( models[[i]] )]
	}

models=do.call(c, models)

# loop to get evaluation metric tables

#modelsfullpath=modelsfullpath[-grep(models, pattern="RF")] # temporary subset while RF models complete
#models=models[-grep(models, pattern="RF")]

# get first model as base

setwd(modelsfullpath[1])
m1=get(load(models[1]))
setwd(toplevel)

evalmean_base=bm_PlotEvalBoxplot(m1, group.by=c("algo", "run"), do.plot=FALSE)
varimp_base=bm_PlotVarImpBoxplot(m1, group.by=c('expl.var', 'algo', "run"), do.plot=FALSE)

evalmean_base=evalmean_base$tab
varimp_base=varimp_base$tab

evalmean_base=cbind(savename=modelsfullpath[1], model=models[1], evalmean_base)
varimp_base=cbind(savename=modelsfullpath[1], model=models[1], varimp_base)


for (m in 2:length(modelsfullpath)) {
	print(m)

	# get model described by iteration
	setwd(modelsfullpath[m])
	modeldetailsm=get(load(models[m]))
	setwd(toplevel)
	
	evalmean_m=bm_PlotEvalBoxplot(modeldetailsm, group.by=c("algo","run"), do.plot=FALSE)
	varimp_m=bm_PlotVarImpBoxplot(modeldetailsm, group.by=c('expl.var', 'algo', "run"), do.plot=FALSE)

	evalmean_m=evalmean_m$tab
	varimp_m=varimp_m$tab

	evalmean_base=rbind(evalmean_base, cbind(savename=modelsfullpath[m], model=models[m], evalmean_m))
	varimp_base=rbind(varimp_base, cbind(savename=modelsfullpath[m], model=models[m], varimp_m))

	}
	
# write full table
write.csv(evalmean_base, "collected_model_evaluations.csv")
write.csv(varimp_base, "collected_variable_importances.csv")

# Response Curves (requires evalmeans to be calculated first, currently built around ROC and TSS only)

setwd(toplevel)

evalmean=evalmean_base

dirrunlist=character(0)
modelrunlist=character(0)
runlist=character(0)

for (s in levels(factor(evalmean$savename))) {
	for (a in levels(factor(evalmean$algo))) {
		print(paste(s,a))
		
		evalmeansub=subset(evalmean, savename==s & algo==a)	#get each combination of dataset and model algorithm
		
		evalmeanROCs=subset(evalmeansub, select=validation, metric.eval=="ROC")
		evalmeanTSSs=subset(evalmeansub, select=validation, metric.eval=="TSS")
		metric.evals=cbind(1:(nrow(evalmeansub)/2), evalmeanTSSs, evalmeanROCs)	# build a data frame describing combinations of TSS and ROC with a row index
		names(metric.evals)=c("Index", "TSS", "ROC")

		attach(metric.evals)
		metric.evals=metric.evals[order(TSS, ROC, decreasing=TRUE),]	#get the greatest value combination, prioritising TSS
		detach(metric.evals)

		dirrunlist=c(dirrunlist, evalmeansub$savename[((metric.evals[1,][,1]-1)*2)+1])	# append the directory locations for modelrunlist
		modelrunlist=c(modelrunlist, evalmeansub$model[((metric.evals[1,][,1]-1)*2)+1])	# append the models to call
		runlist=c(runlist, evalmeansub$run[((metric.evals[1,][,1]-1)*2)+1])	# append the list of model runs to be used
		}
	}
	
#runlist=as.numeric(substr(runlist, start=4, stop=nchar(runlist)))	# get run with best evaluations by TSS and ROC
runlist=paste(runlist, levels(factor(evalmean$algo)), sep="_")

# model first response curves as a base

setwd(dirrunlist[1])
r1=get(load(modelrunlist[1]))
rmodel1=get_built_models(r1)
rrun1=runlist[1]
rmodel1=rmodel1[grep(rmodel1, pattern=rrun1)]
setwd(toplevel)

response_base=bm_PlotResponseCurves(r1, models.chosen=rmodel1, fixed.var="median", do.plot=FALSE)
response_base=response_base$tab
response_base=cbind(savename=dirrunlist[1], model=runlist[1], response_base)

# model response curves and get the run with the best evaluation, then append into a table

for (r in 2:length(runlist)) {

	setwd(dirrunlist[r])
	rn=get(load(modelrunlist[r]))
	rmodeln=get_built_models(rn)
	rrunn=runlist[r]
	rmodeln=rmodeln[grep(rmodeln, pattern=rrunn)]
	setwd(toplevel)

	response_r=bm_PlotResponseCurves(rn, rmodeln, fixed.var="median", do.plot=FALSE)
	response_r=response_r$tab
	response_r=cbind(savename=dirrunlist[r], model=runlist[r], response_r)

	response_base=rbind(response_base, response_r)

	}

write.csv(response_base, "collected_best_fit_response_curves.csv")
