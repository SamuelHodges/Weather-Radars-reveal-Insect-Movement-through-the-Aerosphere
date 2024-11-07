
#list structure:
# variable
#	iteration
#		row index / values

#Args
#
# 1 - table date, yyyymmdd, yyyymm, or yyyy
# 2 - fraction of table entries to retain
# 3 - number of random subset iterations, positive integer
# 4 - number of variables to consider, index starts after the env_Altitude variable, integer value of at least 1
# 5 - directory to find tables

# note that if a variable contains a large number of 0s (Q1, Q2 or Q3 are 0), this method will fail!

Args=commandArgs(trailingOnly=TRUE)
setwd(as.character(Args[5]))

# load data
if (length(list.files(pattern=paste("Bioscatter_biomod_", Args[1], sep="")))<2) {
	original_table=read.csv(list.files(pattern=paste("Bioscatter_biomod_", Args[1], sep="")))[,-1]
} else {
	original_table=list.files(pattern=paste("Bioscatter_biomod_", Args[1], sep=""))
	original_table=lapply(original_table[-grep(original_table, pattern="retained")], read.csv) # grep out exisiting subsampled tables
	original_table=do.call(rbind, original_table) 
	original_table=original_table[,-1] }

retention=round( as.numeric(Args[2]) * nrow(original_table) )
iterations=as.numeric(Args[3])

# create randomly subsetted tables, stored as a list of iterations

tabindex=1:nrow(original_table)
subtabs=rep( list(tabindex), iterations ) # set hierarchical level for no. of random subsamples, initialise index

# randomly subsample the index to the retention % , while preserving proportion of presences to absences
subtabs=lapply(subtabs, function(x) { 

x1=which(original_table[,2]==TRUE);
x0=which(original_table[,2]==FALSE);

retention1=round( as.numeric(Args[2]) * length(x1) )
retention0=round( as.numeric(Args[2]) * length(x0) )

c(sample(x1, retention1), sample(x0, retention0)) } )

subtabindices=subtabs # backup subtabs for later use
#subtabs=rep( list(subtabs), Args[4]) # set hierarchical level for variable

# create corresponding list hierarchy for storing summary statistics
subsums=rep( list(0), iterations )
subsums=rep( list(subsums), Args[4])

# redefine random index to summary statistics of each variable, per iteration
for (v in 1:(Args[4]) ) {
	for (i in 1:iterations) {
		subsums[[v]][[i]] = summary(original_table[,v+8][ subtabs[[i]] ])
		subsums[[v]][[i]] = c(subsums[[v]][[i]], var(original_table[,v+8][ subtabs[[i]] ]))
		names(subsums[[v]][[i]])[7]="Variance"
		next }
	next }

subsums_original=subsums # backup subsums

# take ratios compared to original table, then take the difference between these ratios and 1 squared (see sum of squares approximation)
for (v in 1:(Args[4])) {
	original_v_summary=c(summary(original_table[,v+8]), var(original_table[,v+8]))
	subsums[[v]]=lapply(subsums[[v]], function(x,y){ ((x/y)-1)^2 }, original_v_summary)
	subsums[[v]]=lapply(subsums[[v]], function(x){sum(x[c(2:5, 7)])})		#create representivity scores
	next }

# sum scores across variables and get lowest total score (this is the best one)
scores=do.call(rbind, subsums)
scores=matrix(as.numeric(scores), nrow(scores), ncol(scores)) # correct out of list format
scores=apply(scores, 2, sum) # get total score

table_to_keep = original_table[ subtabindices[[which(scores==min(scores))]], ] # this iteration is the table to keep

# print and compare summaries using stored subsums_original
print(paste(Args[1], length(subtabs[[which(scores==min(scores))]]), "/", nrow(original_table), "entries retained" ))
print("")
for (v in 1:(Args[4]) ) {
	#print
	original_v_summary=c(summary(original_table[,v+8]), var(original_table[,v+8]))
	names(original_v_summary)[7]="Variance"
	print(paste("Random Sample", names(original_table)[v+8]) )
	print(subsums_original[[v]][[which(scores==min(scores))]])
	print(paste("Original", names(original_table)[v+8]))
	print(original_v_summary) 
	print("") }
	print("")

filename_retention=gsub("\\.",",",Args[2]) # replace the . with a , for easier storage
write.csv(table_to_keep, paste("Bioscatter_biomod", Args[1], "binned", filename_retention, "retained.csv", sep="_"))
