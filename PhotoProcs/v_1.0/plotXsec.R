setwd("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcsTmp")
options(warn=2)

leiden = FALSE
reso = 1

libs = c('hdf5r','repmis','RandomFields','inlmisc')
repmis::LoadandCite(libs,file='packages.bib')
RFoptions(spConform=FALSE)
source_dir ='./Leiden/'

cols     = rev(inlmisc::GetColors(11))[1:10]
cols_tr  = rev(inlmisc::GetColors(11, alpha = 0.1))[1:10]
cols_tr2 = rev(inlmisc::GetColors(11, alpha = 0.5))[1:10]

# Misc functions ####
getXShdf5 = function (species,source_dir='./Leiden/') {

  # Misc. info
  info = read.csv(
    file = paste0(source_dir,'cross_section_properties.csv'),
    header=TRUE,
    skip=18,
    stringsAsFactors = FALSE)

  # Get out if species not available in dataset
  if(!(species %in% trimws(info$species)))
    return(NULL)

  # Species OK: proceed
  rownames(info) = trimws(info$species)
  uncCode =  trimws(info[species,'cross_sec_unc'])
  uF = switch(
    uncCode,
    'A+' = 1.2,
    'B'  = 1.3,
    'C'  = 2,
    'D'  = 10
  )

  # Cross sections
  file = paste0(source_dir,'cross_sections/',species,'/',species,'.hdf5')
  fh5 = hdf5r::H5File$new(file,mode='r')
  xs = list()
  for (key in c('photoabsorption',
                'photodissociation',
                'photoionisation',
                'wavelength')       ) {
    xs[[key]] =fh5[[key]]$read()
  }
  # Add XS uncertainty factor to list
  xs[['uncF']] = uF

  return(xs)
}
downSample = function(wl, xs, reso=1) {

  # Define new grid
  lims = round(range(wl))
  wl1 = seq(min(lims), max(lims), by = reso)
  nl1 = length(wl1)

  # Compute increments on original grid
  dwl = diff(wl)
  ## Remove intervals of null width
  sel = dwl != 0
  dwl = dwl[sel]
  sel = c(sel,FALSE)
  pI = xs[sel] * dwl
  wl = wl[sel]

  # Interpolate from irregular to regular grid
  p1 = rep(NA, length(wl1) - 1)
  i = 0
  for (il in wl1[1:(nl1 - 1)]) {
    i = i + 1
    # Which points within interval
    sel = wl >= il - 0.5 * reso &
      wl <  il + 0.5 * reso
    if (sum(sel) == 0) {
      # No point in interval: increase upper limit
      sel1 = which(wl > il + 0.5 * reso)[1]
      if (is.na(sel1)) {
        # No point within full range: contribution is null
        p1[i] = 0
      } else {
        # Found some points
        if (i == 1) {
          # If first point on grid assign first value
          p1[i] = pI[sel1] / dwl[sel1]
        } else {
          # If not first point, linear interpolation from previous value
          x0 = il - reso
          v0 = p1[i - 1]
          x1 = wl[sel1]
          v1 = pI[sel1] / dwl[sel1]
          p1[i] = v0 + (v1 - v0) / (x1 - x0) * reso
        }
      }
    } else {
      # At least one point in gegular interval: sum contributions
      p1[i] = sum(pI[sel]) / sum(dwl[sel])
    }
  }
  # Remove last point
  wl1 = wl1[1:(length(wl1) - 1)]

  return(list(wl = wl1, xs = p1))
}

# Get species list ####
# from photochemistry scheme
nbReac=0
reactants = products = params = type = orig = locnum = list()

filename= 'PhotoScheme.dat'
scheme  = read.fwf(file=filename, widths= rep(11,12))
scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
for (i in 1:nrow(scheme)) {
  nbReac = nbReac + 1
  terms=scheme[i,1:2]
  reactants[[nbReac]] = terms[!is.na(terms) & terms!="" & terms!="HV"]
  terms=scheme[i,3:6]
  products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
  terms=scheme[i,7:12]
  params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
  type[[nbReac]]      = 'photo'
  locnum[[nbReac]]    = i
  orig[[nbReac]]      = filename
}
species=levels(as.factor(unlist(reactants)))
nbSpecies=length(species)
print("Species List:")
print(species)

# Build loss matrix ####
L = matrix(0,ncol=nbSpecies,nrow=nbReac)
for (m in 1:nbReac) {
  reac = unlist(reactants[m])
  for (n in 1:nbSpecies) {
    search=species[n]
    L[m,n] = length(which( search == reac )) # Loss
  }
}

# Generate XS ####
for (i in seq_along(species)) {
  sp=species[i]; print(sp)
  if(leiden) {
    # Check if xs available
    xs  = getXShdf5(sp)
    if(!is.null(xs)) {
      wl   = xs$wavelength
      xsa  = xs$photoabsorption
      xs1  = downSample(wl,xsa,reso=1)
      wavl = round(xs1$wl)
      xsec = matrix(0,ncol=1,nrow=max(wavl))
      xsec[wavl] = xs1$xs
    }
  } else {
    x = read.table(paste0('Cross-sections/se',sp,'.dat'),
                   col.names=1:10,
                   header=FALSE, fill=TRUE)
    wavl = x[,1]
    xsec = matrix(0,ncol=1,nrow=max(wavl))
    xsec[wavl] = x[,2]

  }
  file=paste0('figXsec_',sp, ifelse(leiden,'_Leiden',''), '.png')
  png(file,width=1200,height=1200)
  par(
    pty = 'm',
    mar = c(3.5, 3.5, 1.6, .2),
    mgp = c(2, .75, 0),
    tcl = -0.5,
    lwd = 2,
    cex = 2,
    yaxs = 'i'
  )
  plot(
    xsec,
    log = '',
    type = 'l',
    col = 'gold',
    lwd = 3,
    xlim = c(50, 250),
    # ylim = c(0, 1.2*max(xsec)),
    ylim = c(0, 2e-16),
    main = sp,
    xlab = 'Wavelength [nm]',
    ylab = expression(paste('Cross-section [', cm ^ 2, ']'))
  )
  grid()
  # Partial cross-sections
  reacs = which(L[, i] != 0)
  leg = c()
  for (j in 1:length(reacs)) {
    jreac   = reacs[j]
    channel = params[[jreac]][1]
    if (channel != 0) {
      file = paste0('Cross-sections/qy', sp, '_', channel, '.dat')
      x = read.table(file)
      w  = x[, 1]
      qy = matrix(1e-25, ncol = 1, nrow = max(c(w, wavl)))
      qy[w] = x[, 2] * xsec[w]
      if(max(w) < max(wavl))
        for(i in (max(w)+1):max(wavl))
          qy[i] = x[nrow(x),2] * xsec[i]
      lines(qy,
            col = cols[j],
            lwd = 2,
            lty = 2)
      leg[j] = paste(file, ':',
                     paste(unlist(products[jreac]),
                           collapse = ' + ')
      )
    }
  }
  if (length(leg) != 0)
    legend(
      'topright',
      leg,
      col = cols[1:length(leg)],
      lwd = 2,
      lty = 2
    )
  box()
  dev.off()
}

