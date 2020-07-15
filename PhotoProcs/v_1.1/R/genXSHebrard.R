# Generate fixed resolution absorption cross-sections
# from SWRI data

library(here)
setwd(here::here())

resolutions = c(1, 0.1)

source_dir = './Data/Hebrard/'
target_dir = paste0('./Generated/Hebrard/')

source('./R/functions.R')

processFileHebrard <- function(file_from,file_to,reso) {
  S = read.table(file = file_from,
                 header = FALSE,
                 check.names = FALSE,
                 fill = TRUE)
  wl  = S[, 1]
  xs  = S[, 2]
  xsl = downSample(wl, xs, reso = reso)
  wl1  = xsl$wl
  xs1  = xsl$xs
  write.table(
    cbind(wl1, xs1),
    file = file_to,
    sep = ' ',
    col.names = FALSE,
    row.names = FALSE
  )
}

# Get species list in photoprocess scheme
l = getSpecies('PhotoSchemeGen.dat')
species = l$species
cat('Species List:', species, '\n')


for (sp in species) {
  print(sp)

  # Cross-section
  file_from = paste0(source_dir, 'se', sp, '.dat')
  if (!file.exists(file_from)) {
    cat('No data for:', sp, ', skipping...\n')
    next
  }
  for (reso in resolutions) {
    file_to = paste0(target_dir, reso, 'nm/se', sp, '.dat')
    processFileHebrard(file_from,file_to,reso)
  }

  # Quantum yields
  pattern = paste0('qy', sp)
  files = list.files(source_dir, pattern)
  for (file in files) {
    file_from = paste0(source_dir, file)
    for (reso in resolutions) {
      file_to   = paste0(target_dir, reso, 'nm/', file)
      processFileHebrard(file_from,file_to,reso)
    }
  }
}
