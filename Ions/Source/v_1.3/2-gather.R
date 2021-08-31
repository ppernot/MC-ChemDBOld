# 1/ Gather samples
# 2/ Generate html index for reactions summaries
# 3/ Generate PDF doc of database

version   = '1.3'
rootDir   = '/home/pernot/Bureau/Titan-APSIS/MC-ChemDB/'
sourceDir = paste0(rootDir,'Ions/Source/v_',version,'/')
tmpDir    = paste0(rootDir,'Ions/Tmp/v_',version,'/')
publicDir = paste0(rootDir,'Ions/Public/v_',version,'/')

setwd(sourceDir)

# Load libraries #####
libs =c('stringr','xtable')
for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

source('./massCalc.R')

# Parameters
sampleSize    = 500 # Number max of random samples to gather
randomSamples = TRUE # Gather random databases (vs. nominal only)

setwd(paste0(tmpDir,'Reactions/')) # Contains generated samples
listReacs = list.dirs(full.names=FALSE, recursive=FALSE)
listReacs = gsub("./","",listReacs)

# Reorder by increasing mass of (1) ion & (2) other reactants
tabSpecies  = t(data.frame(sapply(listReacs,
                                  function(x)
                                    strsplit(x, ' + ', fixed = TRUE))))
mass=matrix(0,ncol=ncol(tabSpecies),nrow=nrow(tabSpecies))
for (i in 1:ncol(tabSpecies))
  mass[, i] = as.numeric(unlist(sapply(
    tabSpecies[, i],
    FUN = function(x)
      getMassList(x, excludeList = dummySpecies)
  )))
tabSpecies = tabSpecies[order(mass[,1],mass[,2],na.last=FALSE),]
listReacs = apply(tabSpecies,1,function(x) paste0(x,collapse=' + '))

# Clean target files to be appended to
indexFile=paste0(sourceDir,'index.html')
file.remove(indexFile)

spIndexFile=paste0(sourceDir,'spIndex.html')
file.remove(spIndexFile)

dbDir = paste0(publicDir, 'Databases')
if (!dir.exists(dbDir))
  dir.create(dbDir)
fileList = list.files(path = dbDir, full.name = TRUE)
file.remove(fileList)

dataTableFile=paste0(publicDir,'dataTable.html')
sink(file=dataTableFile,append=FALSE)
cat('<TABLE BORDER=0>')
sink(file=NULL)


allSpecies = allBibKeys = c()

for (reac in listReacs) {
  sink(file=dataTableFile,append=TRUE)
  cat('<TR><TD COLSPAN=5><HR size=1></TD></TR>\n')
  sink(file=NULL)
  file.append(
    file1=dataTableFile,
    file2=paste0(tmpDir,'Reactions/',reac,'/dataTable.html')) 
  
  # Generate Html index to summary files
  cat(paste0(reac,'\n'))
  sink(file=indexFile,append=TRUE)
  cat(paste0('<BR><A HREF="./Data/',reac,'/summary.html">',reac,'</A>\n')) 
  sink(file=NULL)
   
  # Generate collated Monte Carlo samples
  maxNum=ifelse(randomSamples,sampleSize,0)
  for (i in 0:maxNum) {
    file.append(
      file1=paste0(publicDir,'Databases/run_',sprintf('%04i',i),'.csv'),
      file2=paste0(tmpDir,'Reactions/',reac,'/Samples/run_',sprintf('%04i',i),'.csv')) 
  }
  
  # Collate full species list
  species = read.csv(file=paste0(tmpDir,'Reactions/',reac,'/species.txt'),
                     sep=' ',header=FALSE,stringsAsFactors = FALSE)
  allSpecies=c(allSpecies,unlist(species))                   

  # Collate full biblio
  file=paste0(tmpDir,'Reactions/',reac,'/bibKeys.txt')
  if( file.info(file)$size != 0) {
    bibKeys = read.csv(file,sep=' ',header=FALSE,stringsAsFactors = FALSE)
    allBibKeys=c(allBibKeys,unlist(bibKeys))                   
  }
  
}
allSpecies = unique(allSpecies)
masses = sapply(allSpecies, getMassList)

sink(file=dataTableFile,append=TRUE)
cat('</TABLE>')
sink(file=NULL)

# Generate prod-loss file 
sink(file=spIndexFile,append=FALSE)

cat('<H2>Neutrals</H2>')
selIons=grepl('\\+$',allSpecies)
spec = allSpecies[!selIons]
mass  = masses[!selIons]
mord=order(mass)
for (sp in spec[mord]) {
  specDir=paste0(tmpDir,'Species/',sp)
  
  cat(paste0('<BR><B>',sp,'</B> ')) 
  pFile=paste0(specDir,'/prod.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Productions</A>'))
  
  pFile=paste0(specDir,'/loss.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Losses</A>'))
}

cat('<H2>Ions</H2>')
spec = allSpecies[selIons]
mass  = masses[selIons]
mord=order(mass)
for (sp in spec[mord]) {
  specDir=paste0(tmpDir,'Species/',sp)

  cat(paste0('<BR><B>',sp,'</B> ')) 
  pFile=paste0(specDir,'/prod.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Productions</A>'))

  pFile=paste0(specDir,'/loss.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Losses</A>'))
}
sink(file=NULL)

targetHtml = paste0(sourceDir,'speciesList.html')
sink(file=targetHtml, append=FALSE)
cat('<H1>Species List</H1>')

# Dummies 
selAux = is.na(masses) | masses < 1
cat('<H2>Auxiliary species</H2>')
cat(paste0(allSpecies[selAux],collapse='<BR>'))

# Neutrals
trueSpecies=allSpecies[!selAux]
trueMasses=masses[!selAux]
selIons=grepl('\\+$',trueSpecies)
species = trueSpecies[!selIons]
spMass  = trueMasses[!selIons]
mord=order(spMass)
listSp = c()
for (i in seq_along(mord)) 
  listSp[i] = paste0('<font color="blue">',species[mord[i]],
                     '</font> (',signif(spMass[mord[i]],4),')')
cat('<H2>Neutrals</H2>')
cat(paste0(listSp,collapse='<BR>'))

# Cations
species = trueSpecies[selIons]
spMass  = trueMasses[selIons]
mord=order(spMass)
listSp = c()
for (i in seq_along(mord)) 
  listSp[i] = paste0('<font color="blue">',species[mord[i]],
                     '</font> (',signif(spMass[mord[i]],4),')')
cat('<H2>Cations</H2>')
cat(paste0(listSp,collapse='<BR>'))

sink(file = NULL)

setwd(sourceDir)
listHtml=paste0('"Data/',listReacs,'/summary.html"')
listHtml=paste(listHtml,collapse=' ')
listHtml=paste0('Doc/ReleaseNotes.html ','speciesList.html ',listHtml)
command=paste0('htmldoc --book --toclevels 1 --size A4 ',
               '--compression=5 --fontsize 10 --linkcolor purple ',
               '-f summary.pdf ',listHtml, collapse=' ')
system(command,intern=FALSE)
file.rename(from='summary.pdf',to=paste0(publicDir,'summary.pdf'))

# Full biblio
targetHtml = paste0(tmpDir,'bibliography.html')
sink(file=targetHtml, append=FALSE)
allBibKeys=sort(unique(allBibKeys))
printBib(allBibKeys,bib)
sink(file = NULL)
publicPDF = paste0(publicDir,'bibliography.pdf')
command=paste0('pandoc -V geometry:margin=2cm ',targetHtml,' -o ',publicPDF)
system(command,intern=FALSE)



