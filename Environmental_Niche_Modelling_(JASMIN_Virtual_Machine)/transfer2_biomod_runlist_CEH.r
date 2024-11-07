#biomod2 run script

#load in data
setwd("/data/transfer2_biomod/BirdInsect_AerialTerrestrial")
Bioscatter_2017_0_001=read.csv(list.files(pattern="0,001_retained.csv"))[,-1:-2]
Bioscatter_2017_0_005=read.csv(list.files(pattern="0,005_retained.csv"))[,-1:-2]
Bioscatter_2017_0_01=read.csv(list.files(pattern="0,01_retained.csv"))[,-1:-2]

#biomod2 code

library(biomod2)
Options=BIOMOD_ModelingOptions(GAM=list(algo="GAM_gam"))


#Full dataset # Random Forest Runs are carried out seperately owing to intensity of computation # # GAMs were added later owing to a programming fault with the R-base gam model #

##Retention 0.001

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001_t_only", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001_a_only", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


## Retention 0.005

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005_t_only", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005_a_only", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


## Retention 0.01

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01_t_only", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01_a_only", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("CTA-GLM"), models=c("CTA","GLM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


# Random Forest Runs

##Retention 0.001

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001_t_only", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001_a_only", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


## Retention 0.005

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005_t_only", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005_a_only", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


## Retention 0.01

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01_t_only", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01_a_only", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("RF_only"), models=c("RF"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


# GAM Runs

## Retention 0.001

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001_t_only", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001_a_only", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_001", resp.var=Bioscatter_2017_0_001[,2], expl.var=Bioscatter_2017_0_001[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


## Retention 0.005

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005_t_only", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005_a_only", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_005", resp.var=Bioscatter_2017_0_005[,2], expl.var=Bioscatter_2017_0_005[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)


## Retention 0.01

###Terrestrial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01_t_only", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=8, do.full.models=FALSE, seed.val=342)

###Aerial_only

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01_a_only", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,8,9:16)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=10, do.full.models=FALSE, seed.val=342)

###Terrestrial-Aerial

bioprep=BIOMOD_FormatingData(resp.name="Bioscatter_2017_0_01", resp.var=Bioscatter_2017_0_01[,2], expl.var=Bioscatter_2017_0_01[,c(5,8,9:16,17:28)])

BioModel=BIOMOD_Modeling(bm.format=bioprep, modeling.id=c("GAM_only"), models=c("GAM"), bm.options=Options, nb.rep=100, data.split.perc=80, metric.eval=c("TSS", "ROC"), var.import=17, do.full.models=FALSE, seed.val=342)
