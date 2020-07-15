setwd("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcs/v1.1")
options(warn=2)

leiden = FALSE
reso = 1
source_dir ='./Data/Leiden/'
target_dir ='./Generated/Leiden'


libs = c('hdf5r','repmis','RandomFields','inlmisc')
repmis::LoadandCite(libs,file='packages.bib')
RFoptions(spConform=FALSE)

cols     = rev(inlmisc::GetColors(11))[1:10]
cols_tr  = rev(inlmisc::GetColors(11, alpha = 0.1))[1:10]
cols_tr2 = rev(inlmisc::GetColors(11, alpha = 0.5))[1:10]

source('R/functions.R')

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

