#=================================================#
# Generate Monte Carlo samples for cross-sections #
#=================================================#

setwd("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcsTmp")

# Monte Carlo parameters
nMC = 500
fMC = 0.2 # 20% uncertainty for all cross-sections

# Target directory for chemistry samples
targetMCDir = '../ChemDBPublic/PhotoProcs/'
sourceDir   = './Cross-sections/'

fileList = list.files(path = sourceDir)

# # Clean target dir
# command = paste0('rm -rf ',targetMCDir,'*')
# dummy   = system(command)

for (iMC in 0:nMC) {
  prefix=paste0(sprintf('%04i',iMC),'_')
  for(file in fileList) {
    rnd = rlnorm(1,meanlog = 0, sdlog = 0.2)
    if(substr(file,1,2)== 'se') {
      xS = read.table(paste0(sourceDir,file),
                      header=FALSE,fill=TRUE)
      xS[,2] = xS[,2] * rnd
      write.table(xS[,1:2], sep=' ', row.names=FALSE, col.names=FALSE,
                file=paste0(targetMCDir,prefix,file))
    } else {
      file.copy(from = paste0(sourceDir,file),
                to = paste0(targetMCDir,prefix,file))
    }
  }
}

