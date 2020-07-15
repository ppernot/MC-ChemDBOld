# Generate fixed resolution absorption cross-sections
# from SWRI data

library(here)
setwd(here::here())

source_dir = './Data/SWRI/'
target_dir = paste0('./Generated/SWRI/')

source('./R/functions.R')

# Get species list in photoprocess scheme
l = getSpecies('PhotoSchemeGen.dat')
species = l$species
cat('Species List:', species, '\n')

for (reso in c(1, 0.1)) {
  for (sp in species) {
    file = paste0(source_dir, sp, '.dat')
    if (!file.exists(file)) {
      cat('No data for:', sp, ', skipping...\n')
      next
    }
    S = read.table(file = file,
                   header = TRUE,
                   check.names = FALSE)
    # Convert A to nm
    wl  = S[, 1] / 10

    # Absorption cross-section
    xs  = S[, 2]

    xsl = downSample(wl, xs, reso = reso)
    wl1  = xsl$wl
    xs1  = xsl$xs

    # Remove tailing zeroes
    first = which(xs1 != 0)[1]
    last  = length(xs1)-which(rev(xs1) !=0)[1] + 1
    wl1   = wl1[first:last]
    xs1   = xs1[first:last]

    write.table(
      cbind(wl1, xs1),
      file = paste0(target_dir, reso, 'nm/se', sp, '.dat'),
      sep = ' ',
      col.names = FALSE,
      row.names = FALSE
    )

    # Quantum yields
    np = ncol(S) - 2
    if (np == 0) {
      # Single channel: compute qy from total xs => qy=1
      np = 1
      i0 = 1
    } else {
      # Normal case: several channels
      i0 = 2
    }
    for (i in 1:np) {
      qy = downSample(wl, S[, i+i0]/xs, reso = reso)$xs
      write.table(
        cbind(wl1, qy[first:last]),
        file = paste0(target_dir, reso, 'nm/qy', sp, '_', i, '.dat'),
        sep = ' ',
        col.names = FALSE,
        row.names = FALSE
      )
    }
  }
}
