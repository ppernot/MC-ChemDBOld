library(here)
setwd(here::here())

reso = 1
nMC  = 200
test = TRUE

source_dir =paste0('../../ChemDBPublic/PhotoProcs/',reso,'nm/')
if(test)
  source_dir =paste0('./Test/',reso,'nm/')

libs = c('repmis','RandomFields','inlmisc')
repmis::LoadandCite(libs,file='packages.bib')
RFoptions(spConform=FALSE)

cols     = rev(inlmisc::GetColors(11))[1:10]
cols_tr  = rev(inlmisc::GetColors(11, alpha = 0.1))[1:10]
cols_tr2 = rev(inlmisc::GetColors(11, alpha = 0.5))[1:10]

source('R/functions.R')

# Get species list and network ####
spL = getSpecies('PhotoSchemeGen.dat')
species  = spL$species
L        = spL$L
params   = spL$params
products = spL$products

# if(test)
#   species = 'N2' # species[1]

# XS ####
for (sp in species) {
  print(sp)
  iMC = 0
  prefix=paste0(sprintf('%04i',iMC),'_')
  x = read.table(paste0(source_dir,prefix,'se',sp,'.dat'),
                 col.names=1:10,
                 header=FALSE, fill=TRUE)
  wavl = x[,1]
  X = matrix(NA,nrow=nrow(x),ncol=nMC+1)
  X[,1] = x[,2]
  for (iMC in 1:nMC) {
    prefix=paste0(sprintf('%04i',iMC),'_')
    x = read.table(paste0(source_dir,prefix,'se',sp,'.dat'),
                   col.names=1:10,
                   header=FALSE, fill=TRUE)
    X[,iMC+1] = x[,2]
  }
  file=paste0('Figs/figXsecMC_',sp,'.png')
  png(file,width=2000,height=1200)
  par(
    mfrow=c(1,2),
    pty = 'm',
    mar = c(3.5, 3.5, 1.6, .2),
    mgp = c(2, .75, 0),
    tcl = -0.5,
    lwd = 2,
    cex = 3.5,
    xaxs = 'i',
    yaxs = 'i'
  )
  matplot(
    wavl, X,
    log = '',
    type = 'l',
    col = cols_tr[8],
    lty = 1,
    lwd = 3,
    xlim = c(50, max(wavl)),
    ylim = c(0, 1.1*quantile(X[,1],0.99)),
    main = sp,
    xlab = 'Wavelength [nm]',
    ylab = expression(paste('Cross-section [', cm ^ 2, ']'))
  )
  grid()
  lines(wavl,X[,1],col=2)
  box()

  # Partial cross-sections ####
  reacs = which(L[, sp] != 0)
  leg = c()
  iplot = TRUE
  nc = length(reacs)
  ch_cols     = rev(inlmisc::GetColors(nc+1))[1:nc]
  ch_cols_tr  = rev(inlmisc::GetColors(nc+1, alpha = 0.02))[1:nc]
  for (j in 1:length(reacs)) {
    jreac   = reacs[j]
    channel = as.numeric(params[[jreac]][1])
    if (channel != 0) {
      iMC = 0
      prefix=paste0(sprintf('%04i',iMC),'_')
      x = read.table(paste0(source_dir,prefix,'qy',sp,'_', channel,'.dat'))
      wavq  = x[, 1]
      qy = matrix(NA, nrow = nrow(x), ncol = nMC+1)
      qy[,1] = x[, 2]
      for (iMC in 1:nMC) {
        prefix=paste0(sprintf('%04i',iMC),'_')
        x = read.table(paste0(source_dir,prefix,'qy',sp,'_', channel,'.dat'))
        qy[,iMC+1] = x[,2]
      }

      if(iplot){
        iplot = FALSE
        matplot(
          wavq, qy,
          log = '',
          type = 'l',
          col = cols_tr[channel],
          lty = 1,
          lwd = 3,
          xlim = c(50, max(wavl)),
          ylim = c(-0.01, 1.19),
          main = sp,
          xlab = 'Wavelength [nm]',
          ylab = expression(paste('Branching ratios'))
        )
        grid(); abline(h=1)
      } else {
        matlines(
          wavq, qy,
          log = '',
          col = ch_cols_tr[channel],
          lty = 1,
          lwd = 3
        )
      }
      lines(wavq,qy[,1],col=ch_cols[channel],lty=2,lwd=3)

      leg[j] = paste('qy',sp,'_', channel, ':',
                     paste(unlist(products[jreac]),
                           collapse = ' + ')
      )
    }
  }
  if (length(leg) != 0)
    legend(
      'top', ncol=2, bty='n',
      leg,
      cex = 0.5,
      col = ch_cols[1:length(leg)],
      lwd = 10,
      lty = 1
    )
  box()
  dev.off()
  hist(apply(qy,1,sd),nclass=33,main=sp)
}

