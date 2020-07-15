setwd("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcsTmp")

libs = c('hdf5r','repmis','RandomFields','inlmisc')
repmis::LoadandCite(libs,file='packages.bib')
RFoptions(spConform=FALSE)
source_dir ='./Leiden/'

cols     = rev(inlmisc::GetColors(8))[1:7]
cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.1))[1:7]
cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.5))[1:7]

# Misc functions ####
getXShdf5 = function (species,source_dir='./Leiden/') {
  # Misc. info
  info = read.csv(
    file = paste0(source_dir,'cross_section_properties.csv'),
    header=TRUE,
    skip=18,
    stringsAsFactors = FALSE)
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

nMC = 500
wlScaleFac = 0.2
reso = 1

# Get data ####
species = 'CH4'
xs  = getXShdf5(species)

# Downsample to lower resolution grid ####
wl   = xs$wavelength
xsa  = xs$photoabsorption
xs1  = downSample(wl,xsa,reso=reso)
wl1  = xs1$wl
xsd  = xs$photodissociation
xsd1 = downSample(wl,xsd,reso=reso)
xsi  = xs$photoionisation
xsi1 = downSample(wl,xsi,reso=reso)

# Simulate standard GP to add noise
rho = wlScaleFac*diff(range(wl1))
cond = RandomFields::RFsimulate(
  model = RandomFields::RMgauss(
    var = 1,
    scale = rho,
    Aniso = NULL,
    proj = NULL
  ),
  x     = wl1,
  n     = nMC,
)

# Noisy cross-sections
luF = log(xs$uncF)
lxs = log(xs1$xs)
br  = xsd1$xs / xs1$xs # photodissoc branching ratio
xsu = xsdu = xsiu = matrix(NA, nrow = length(wl1), ncol = nMC)
for (i in 1:nMC) {
  xsu[, i]  = exp(lxs + luF * cond[, i])
  xsdu[, i] = br * xsu[, i]
  xsiu[, i] = (1 - br) * xsu[, i]
}


# Plots ####
par(
  pty = 'm',
  mar = c(3.5, 3.5, 1.6, .2),
  mgp = c(2, .75, 0),
  tcl = -0.5,
  lwd = 1.5,
  cex = 1.2,
  yaxs = 'i'
)
plot(
  wl, xsa,
  type ='n',
  xlab = 'Wavelength [nm]',
  ylab = expression(paste('Cross-section [',cm^2,']')),
  ylim = c(0, xs$uncF^2*max(xsa)),
  main = paste0(species,' / uF = ',xs$uncF)
)
grid()
lines(xs1$wl,xs1$xs,col=2)
for (i in 1:nMC)
  lines(xs1$wl, xsdu[, i], col = cols_tr[5])
points(wl, xs$photodissociation,
  col=1, pch=16, cex=0.5)
for (i in 1:nMC)
  lines(xs1$wl, xsiu[, i], col = cols_tr[2])
points(wl, xs$photoionisation,
  col=1, pch=16, cex=0.5)
for (i in 1:nMC)
  lines(xs1$wl, xsdu[, i] + xsiu[, i], col = cols_tr[3])
points(wl, xs$photoabsorption,
       col=1, pch=16, cex=0.5)
box()
