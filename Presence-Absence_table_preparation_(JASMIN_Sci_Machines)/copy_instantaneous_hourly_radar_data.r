# fileselect - target folder to copy from

Args=commandArgs(trailingOnly=TRUE)

fileselect=Args[1]

hours=seq(0, 23, 1)
for (folder in fileselect) {

#setwd("20170717")
setwd(as.character(folder))
allfiles=list.files()

filehours=substr(allfiles, start=nchar(allfiles)-15, stop=nchar(allfiles)-14)
filetime=as.numeric(substr(allfiles, start=nchar(allfiles)-15, stop=nchar(allfiles)-10))

for (i in hours) {
possiblefiles=grep(x=filehours, pattern=as.character(i))
targetfile=allfiles[which(filetime==min(filetime[possiblefiles]))]

if (length(list.files(path="/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/data/x-band", pattern=targetfile))>0){
print(paste(targetfile, "is already copied to data directory, skipping")) } else {

print(paste("Copying ", targetfile, "...", sep=""))
file.copy(from=paste(getwd(), targetfile, sep="/"), to=paste("/home/users/eesjir/eesjir_group_workspace/transfer2_data_prep/data/x-band", targetfile, sep="/")) }
}
setwd("..")

}
print("Done!")