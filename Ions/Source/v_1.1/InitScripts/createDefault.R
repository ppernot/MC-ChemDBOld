listIons=c("CH3CN+", "C2H5N+", "C2HN2+", "C3H5N+", "C4HN+", "C4H2N+", "C4H4N+", 
"C4H5N+", "C4H6N+", "C4N2+", "C4HN2+", "HC5N+", "C5H2N+", "C5H3N+", 
"C5H4N+", "C5H5N+", "C5H6N+", "C5HN2+", "C6H2N+", "C6H7N+", "C6H8N+", 
"CxHyNz+", "C7H3N+", "C6HN2+", "C6H2N2+", "C7H5N+", "C3HN3+", 
"C2H5CNH+")

version   = '1.0'
sourceDir = paste0('~/Bureau/ChemDB/ChemDB_Source/v_',version,'/')
dataDir   = paste0(sourceDir,'DRData/')

setwd(sourceDir)
X = as.matrix(read.csv(file='Doc/defaultTemplate.csv',
                       header=FALSE, sep='\t', fill=TRUE,
                       stringsAsFactors = FALSE,
                       na.strings="")
)
X = str_trim(X) # remove unwanted spaces

for (ion in listIons) {
  X[1,1] = paste0(ion,' + E')
  targetDir = paste0(dataDir,ion)
  dir.create(targetDir)
  targetFile = paste0(dataDir,ion,'/data.csv')
  write.table(X,file=targetFile, sep='\t', na="",
            row.names = FALSE,
            col.names = FALSE)
}
