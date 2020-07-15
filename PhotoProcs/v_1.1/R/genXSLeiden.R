# Generate fixed resolution absorption cross-sections
# from Leiden hdf5 data

test = FALSE

setwd("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcs/v1.1")

libs = c('hdf5r','repmis')
repmis::LoadandCite(libs,file=NULL)

source_dir ='Data/Leiden/'
target_dir =paste0('Generated/Leiden/')

source('R/functions.R')

# Get species list in photoprocess scheme
l = getSpecies('PhotoScheme.dat')
species = l$species
cat('Species List:',species,'\n')

# Generate fixed-res data ####
for(reso in c(1,0.1)) {
  for( sp in species) {
    xsl  = getXShdf5( sp, source_dir = source_dir)
    if(is.null(xsl)) {
      cat('No data for:',sp,', skipping...\n')
      next
    }
    uF   = xsl$uncF
    wl   = xsl$wavelength
    xs   = xsl$photoabsorption
    xsl  = downSample(wl,xs,reso=reso)
    wl   = xsl$wl
    xs   = xsl$xs
    cat('Creating ',paste0(target_dir,reso,'nm/se',sp,'.dat\n'))
    cat('Creating ',paste0(target_dir,reso,'nm/uF',sp,'.dat\n'))
    if(!test) {
      write.table(
        cbind(wl,xs),
        file = paste0(target_dir,reso,'nm/se',sp,'.dat'),
        sep = ' ',
        col.names = FALSE,
        row.names = FALSE)
      write.table(
        uF,
        file = paste0(target_dir,reso,'nm/uF',sp,'.dat'),
        sep = ' ',
        col.names = FALSE,
        row.names = FALSE)
    }
  }
}
