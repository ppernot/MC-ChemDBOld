#=================================================#
# Generate Monte Carlo samples for cross-sections #
#=================================================#
library(DiceKriging)
library(here)

# setwd(here::here()) # Define working directory

test = TRUE
GP_Fit = FALSE
eps = 5e-3 # Threshold for zero in compositions

resolutions = c(1, 0.1)[1:2]

# Monte Carlo parameters
nMC = 500

# Default uncertainties
uF0    = 1.2  # 20% default uncertainty factor for cross-sections
ruBRN  = 0.2  # relative uncertainty for Neutral branching ratios
ruBRI  = 0.2  # relative uncertainty for Ionic branching ratios
ruBRNI = 0.03 # relative uncertainty on Ionic vs. Neutral channels

# Target directory for chemistry samples
sourceDir   = './Cross-sections/'
uncFDir     = './Generated/Leiden/'
targetMCDir = '../../../Reactor_Runs/ChemDBPublic/PhotoProcs_v_1.1/'
# if(test)
#   targetMCDir = './Test/'

cols     = rev(inlmisc::GetColors(11))[1:10]
cols_tr  = rev(inlmisc::GetColors(11, alpha = 0.1))[1:10]
cols_tr2 = rev(inlmisc::GetColors(11, alpha = 0.5))[1:10]
source('./R/functions.R')

# Get species list in photoprocess scheme
l = getSpecies('PhotoSchemeGen.dat')
species = l$species
L        = l$L
params   = l$params
products = l$products


if(test)
  species = 'CO2'

cat('Species List:', species, '\n')

for (reso in resolutions) {
  sourceDir1 = paste0(sourceDir,reso,'nm/')
  targetMCDir1 = paste0(targetMCDir,reso,'nm/')
  for (sp in species) {

    # se ####
    file = paste0(sourceDir1,'se', sp, '.dat')
    if (!file.exists(file))
      stop('No data for:', sp, '\n')
    x = read.table(file,header=FALSE,fill=TRUE)
    wavl = x[,1]
    se   = x[,2]

    ## Define common wavl range for se and qy
    fileList = list.files(
      path = sourceDir1,
      pattern = paste0('qy',sp,'_')
    )
    x = read.table(paste0(sourceDir1,fileList[1]))
    wavlqy = x[,1]

    selXS = wavl %in% wavlqy
    selBR = wavlqy %in% wavl

    ## Adjust XS tables
    wavl = wavl[selXS]
    se   = se[selXS]

    cat('Generating se samples for ',sp,' at ', reso,'nm\n')
    # se uncertainty factor from Leiden data or uF0
    uFile = paste0(uncFDir,reso,'nm/uF',sp,'.dat')
    if(file.exists(uFile)) {
      uF = unlist(read.table(uFile,header=FALSE,fill=TRUE))
    } else {
      uF = uF0
    }

    for (iMC in 0:nMC) {
      prefix=paste0(sprintf('%04i',iMC),'_')
      if(iMC == 0) {
        rnd = 1 # Nominal value
      } else {
        rnd = rlnorm(1, meanlog = 0, sdlog = log(uF))
      }
      seMC = se * rnd
      write.table(
        cbind(wavl, seMC),
        sep = ' ', 
        row.names = FALSE, 
        col.names = FALSE,
        file = gzfile(
          paste0(targetMCDir1, prefix, 'se', sp, '.dat.gz')
        )
      )
    }

    # qy ####
    cat('Generating qy samples for ',sp,' at ', reso,'nm\n')

    ## Identify ionization channels
    reacs = which(L[, sp] != 0)
    ionic = c()
    for (j in 1:length(reacs)) {
      ch = as.numeric(params[[ reacs[j] ]][1])
      if (ch != 0)
        ionic[ch] = 'E' %in% products[[ reacs[j] ]]
    }

    ## Get BR data
    
    #### Build file list in numeric order (fileList is alphabetic) 
    nch = length(fileList)
    qyFiles = paste0('qy',sp,'_',1:nch,'.dat') 
    
    x = read.table(paste0(sourceDir1,qyFiles[1]))
    qy = matrix(NA,nrow=sum(selBR),ncol=nch)
    wavl   = x[selBR,1]
    qy[,1] = x[selBR,2]
    for(i in 2:nch) {
      x = read.table(paste0(sourceDir1,qyFiles[i]))
      qy[,i] = x[selBR,2]
    }
    nw = length(wavl)
    nc = ncol(qy)

    ## nominal run
    iMC = 0
    prefix = paste0(sprintf('%04i',iMC),'_')
    # save
    for(i in 1:nc) {
      fileOut = paste0(targetMCDir1, prefix, 'qy', sp, '_', i, '.dat')
      write.table(
        cbind(wavl, qy[, i]),
        sep = ' ',
        row.names = FALSE,
        col.names = FALSE,
        file = gzfile(paste0(fileOut, '.gz'))
      )
    }


    if(nc == 1) {
      # A single channel: no uncertainty in BR
      qySample = array(
        data = 1,
        dim = c(nMC, nw, nc)
      )

    } else {
      
      qy = qy / rowSums(qy)

      if(GP_Fit) {
        # Experimental...

        qySample = array(
          data=NA,
          dim=c(nMC,nw,nc)
        )
        ## Apply threshold
        mask = qy < eps
        qy[mask] = eps

        # Transform to logit-space
        qy1  = log(qy/(1-qy))
        noise = 0.2 * sqrt(qy*(1-qy)) # Distrib of uncert peaks at 0.06
        unc1 = noise/(qy*(1-qy))

        par(
          mfrow=c(1,1),
          pty = 'm',
          mar = c(3.5, 3.5, 1.6, .2),
          mgp = c(2, .75, 0),
          tcl = -0.5,
          lwd = 2,
          cex = 2,
          xaxs = 'i',
          yaxs = 'i'
        )
        matplot(
          wavl, qy1,
          log = '',
          type = 'l',
          col = cols,
          lty = 1,
          lwd = 3,
          xlim = c(50, max(wavl)),
          # ylim = c(0, 1.1*quantile(X[,1],0.99)),
          main = sp,
          xlab = 'Wavelength [nm]',
          ylab = 'log-BR'
        )

        ##  GP
        gp=list()
        kernel = c('gauss','matern5_2','matern3_2','exp')[2]
        coef.var = 20.0
        coef.cov = 10 # nm
        sel = seq(1, nw, length.out = 50)
        for (i in 1:nc) {
          gp[[i]] =km(design   = data.frame(x=wavl[sel]),
                      response = data.frame(y=qy1[sel,i]),
                      covtype  = kernel,
                      coef.trend = mean(qy1[,i]),
                      coef.cov = coef.cov,
                      coef.var = coef.var,
                      noise.var= unc1[sel,i]^2
          )
        }
        matpoints(
          wavl[sel],
          qy1[sel,],
          pch=16,
          col=cols
        )
        for (i in 1:nc)
          segments(
            wavl[sel],qy1[sel,i]+2*unc1[sel,i],
            wavl[sel],qy1[sel,i]-2*unc1[sel,i],
            lty =1, lwd=3,
            col=cols[i]
          )

        for (iMC in 1:nMC) {
          pred=matrix(nrow=nw,ncol=nc)
          for(i in 1:nc) {
            p = simulate(
              gp[[i]],
              newdata=data.frame(x=wavl),
              cond=TRUE
            )
            pred[,i] = p
          }
          matlines(
            wavl,
            pred,
            col = cols_tr,
            lty = 1,
            lwd = 3)

          ## back-transform
          brPred = 1/(1+exp(-pred))
          brPred[mask] = 0 # Restore true 0s
          qySample[iMC,,] = brPred/rowSums(brPred)
        }

      } else {
        # Generate ordered Diri-based samples at each wavelength

        if (sum(ionic) * sum(!ionic) == 0) {
          # Diri sampling
          qySample = diriSample(qy,
                                ru = ifelse(
                                  sum(ionic)==0,
                                  ruBRN,
                                  ruBRI),
                                nMC,
                                eps)
        } else {
          # Nested sampling
          qySample = hierSample(qy, ionic,
                                ru = c(ruBRNI, ruBRN, ruBRI),
                                nMC,
                                eps)
        }
      }
    }

    # Save sets in random order
    shuffle = sample.int(n=nMC,size=nMC)
    for(iMC in 1:nMC) {
      for(i in 1:nc) {
        # Randomize global sample through file prefix
        prefix = paste0(sprintf('%04i',shuffle[iMC]),'_')
        fileOut = paste0(targetMCDir1,prefix,'qy',sp,'_',i,'.dat')
        write.table(
          cbind(wavl,qySample[iMC,,i]),
          sep=' ', row.names=FALSE, col.names=FALSE,
          file = gzfile(paste0(fileOut, '.gz'))
        )
      }
    }
  }
}

