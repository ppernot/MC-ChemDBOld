# Generate fixed resolution absorption cross-sections
# from SWRI data

setwd("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcs/v1.1")
source_dir = 'Data/Plessis/'
target_dir = paste0('Generated/Plessis/')

source('R/functions.R')

# Get species list in photoprocess scheme
species = 'CH4'

for (reso in c(1, 0.1)) {
  for (sp in species) {
    path    = paste0(source_dir)
    pattern = paste0('qy', sp)
    files = list.files(path, pattern)
    for (file in files) {
      S = read.table(file = paste0(path,file),
                     header = TRUE,
                     check.names = FALSE)
      wl  = S[, 1]
      xs  = S[, 2]
      xsl = downSample(wl, xs, reso = reso)
      wl1  = xsl$wl
      xs1  = xsl$xs
      write.table(
        cbind(wl1, xs1),
        file = paste0(target_dir, reso, 'nm/',file),
        sep = ' ',
        col.names = FALSE,
        row.names = FALSE
      )
    }
  }
}
