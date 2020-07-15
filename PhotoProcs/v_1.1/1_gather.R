# Copy files from generated to cross-sections
# according to preference pattern
library(here)
# setwd(here::here())

source('R/functions.R')
source_dir = './Generated/'
target_dir = './Cross-sections/'

# Define scheme
leiden  = TRUE # If TRUE, take prefrentially Leiden cross-sections
                # Except for C2H4, which is restricted in wavelength range
plessis = TRUE # Prefer Plessis's branching ratios for CH4

# Get species list in photolysis scheme
l = getSpecies('PhotoSchemeGen.dat')
species = l$species
cat('Species List:', species, '\n')

# Process files
prov = c()
for (reso in c(1, 0.1)) {
  ifile = 0
  for (sp in species) {
    cat(sp,'\n')
    # Cross-sections

    ## SWRI > Hebrard
    file = paste0(source_dir, 'Hebrard/', reso, 'nm/se', sp, '.dat')
    file1 = paste0(source_dir, 'SWRI/', reso, 'nm/se', sp, '.dat')
    if (file.exists(file1))
      file = file1

    ## Leiden > (SWRI,Hebrard)
    if (leiden & sp != 'C2H4') {
      file1 = paste0(source_dir, 'Leiden/', reso, 'nm/se', sp, '.dat')
      if (file.exists(file1))
        file = file1
    }
    if (!file.exists(file))
      stop(paste0('*** Pb with se',sp))

    ## Copy
    file_to = paste0(target_dir, reso, 'nm/se', sp, '.dat')
    ifile = ifile + 1
    prov[ifile] = file
    file.copy(from = file,
              to = file_to,
              overwrite = TRUE)

    # ## Get wavelength range to check compatibility with qy files
    # wavlXS = read.table(file,header=FALSE,fill=TRUE)[,1]

    # Quantum yields
    path    = paste0(source_dir, 'SWRI/', reso, 'nm')
    pattern = paste0('qy', sp)
    files = list.files(path, pattern)

    ## Plessis > SWRI (for CH4)
    if(plessis) {
      path1   = paste0(source_dir, 'Plessis/', reso, 'nm')
      files1 = list.files(path1, pattern)
      if (length(files1) != 0) {
        path = path1
        files = files1
      }
    }

    ## copy
    for (file in files) {
      file_from = paste0(path,'/',file)
      ifile = ifile + 1

      # ## Get wavelength range to check compatibility with se files
      # wavlqy = read.table(file_from,header=FALSE,fill=TRUE)[,1]
      # if (min(wavlqy) > min(wavlXS) )
      #   cat('Min qy wavl too large\n')
      # if (max(wavlqy) < max(wavlXS) )
      #   cat('Max qy wavl too small\n')

      prov[ifile] = file_from
      file_to   = paste0(target_dir, reso, 'nm/', file)
      file.copy(from = file_from,
                to = file_to,
                overwrite = TRUE)
    }
  }

  # Record provenance track
  sink(file=paste0(target_dir,reso, 'nm/0_Provenance.txt'))
  for(i in 1:length(prov))
    cat(' -',prov[i],'\n')
  sink()
}

