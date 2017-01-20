version   = '1.0'
rootDir   = '/home/pernot/Bureau/Titan-APSIS/ChemDB/'
subDB     = 'Photo/'
sourceDir = paste0(rootDir,subDB,'Source/v_',version,'/')
targetDir = paste0(rootDir,subDB,'Tmp/v_',version,'/')
publicDir = paste0(rootDir,subDB,'Public/v_',version,'/')
docDir    = paste0(sourceDir,'Doc/')

# Load libraries #####
library(stringr)
library(xtable)
library(bibtex)
library(RColorBrewer)

# Load functions #####
setwd(sourceDir)
source('Scripts/myApe.R')         # Tree plots
source('Scripts/stoechiometry.R') # Stoechiometry

# Additional functions
col2tr = function(col,alpha=20){  
  rgb(t(col2rgb(col)),alpha=alpha,maxColorValue=255)
} 
capwords = function(s, strict = TRUE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s = substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
sanitize = function(str) {
  gsub('([#$%&~_\\^\\\\{}\\s\\(\\)])', '\\\\\\1', str, perl = TRUE)
}

# Construction of pdf strings
getDistString = function (X,kwd) {
  topLeft = which(X==kwd,arr.ind=TRUE)
  if(length(topLeft)==0) {
    dist = 'Delta' # default
    pars = '0'
  } else {  
    dist = X[topLeft[1],topLeft[2]+1]
    pars = X[topLeft[1],(topLeft[2]+2):ncol(X)]
    pars = pars[!is.na(pars)]
  }
  string=paste0(dist,'(',paste0(pars,collapse=','),')')
}
getParams = function (X,kwd) {
  topLeft = which(X==kwd,arr.ind=TRUE)
  if(length(topLeft)==0) {
    pars = NA
  } else {  
    pars= X[topLeft[1],(topLeft[2]+1):ncol(X)]
    if(any(!is.na(pars))) {
      pars = pars[!is.na(pars)]
    } else {
      pars = NA
    }
  }
}
sampleDistString = function (distString,sampleSize) {
  spl = unlist(strsplit(distString,'(',fixed=TRUE))
  if(spl[1] =='Delta') {
    par1 = strsplit(spl[2],',',fixed=TRUE)[1]
    par1 = sub(')','',par1)
    sample = matrix(as.numeric(par1),ncol=1,nrow=sampleSize)
  } else {
    sample = nds(sampleSize,distString)    
  }
  return(sample)
}
nds = function(ns,dist) {
  command=paste("echo ",ns," '",dist,"' | ../Rnested.x") 
  # quotes around dist avoid shell interpretation
  tc=textConnection(system(command,intern=T))
  liste=scan(tc,quiet=TRUE)
  close(tc)
  nleaves=liste[1]
  nlast=nleaves*ns+1
  nds=matrix(liste[2:nlast],ncol=nleaves,byrow=T)
  return(nds)
}
oneDist = function(ic,d,tags,tagged=FALSE) {
  dSub=d[[ic]]
  string=paste0(dSub$dist,'(')
  sSub=c()
  for (ip in 1:length(dSub$elem)) {
    sSub[ip]=paste0(dSub$mu[ip],
                    ifelse(dSub$link[ip]==0,
                           ifelse(tagged,
                                  paste0(":'",tags[dSub$elem[ip]],"'"),
                                  ''),
                           paste0('*LINK/',dSub$link[ip],'/')
                    ),
                    collapse='')
  }
  sig=dSub$sig[which(!is.na(dSub$sig))]
  return(paste0(string,
                paste0(sSub,collapse=','),
                ifelse(length(sig)==0,
                       '',
                       paste0(';',paste0(sig,collapse=','))
                ),                
                ')')
  )
}
oneNewick = function(ic,d,tags,tagged=TRUE) {
  dSub=d[[ic]]
  sSub=c()
  for (ip in 1:length(dSub$elem)) {
    sSub[ip]=paste0(
                    ifelse(dSub$link[ip]==0,
                           tags[dSub$elem[ip]],
                           paste0('LINK/',dSub$link[ip],'/')
                    ),
                    collapse='')
  }
  return(
    paste0("(",paste(sSub,collapse=","),")")
  )
}
oneEdgeTag = function(ic,d,tags) {
  dSub=d[[ic]]
#   sig=dSub$sig[which(!is.na(dSub$sig))]
  mu = as.numeric(dSub$mu)
  sig= as.numeric(dSub$sig)   
  sSub=c()
  for (ip in 1:length(dSub$elem)) {
    label = switch(dSub$dist,
                   Dirg = {paste0(signif(mu[ip],2),' +/- ',
                                  signif(sig[ip],2))},
                   Mlgn = {paste0(signif(mu[ip],2),' */: ',
                                  signif(sig[ip],2))},
                   Diut = {paste0('[',signif(mu[ip],2),'; ',
                                  signif(sig[ip],2),']')},
                   Diri = {paste0(signif(mu[ip]/sum(mu),2))},
                   Dior = {paste0(signif((sum(mu)-mu[ip])/sum(mu),2))}
    )
    if(dSub$link[ip]==0) {
      sSub[ip]=label  
    } else {
      sSub[ip]=paste0(label,',LINK/',dSub$link[ip],'/')
    }
  }
  return(paste0(sSub,collapse=","))
}
oneNodeTag = function(ic,d,tags) {
  dSub=d[[ic]]
  #   sig=dSub$sig[which(!is.na(dSub$sig))]
  sSub=c()
  for (ip in 1:length(dSub$elem)) {
    if(dSub$link[ip]==0) {
      sSub[ip]=dSub$dist[ip]  
    } else {
      sSub[ip]=paste0(dSub$dist[ip],',LINK/',dSub$link[ip],'/')
    }
  }
  return(paste0(sSub,collapse=","))
}
getSpecies = function (chain, maxProds) {
  species = unlist(strsplit(chain,' + ',fixed=TRUE))
  if(length(species) > maxProds) stop (paste0('Nb of products exceds maxProd=',maxProds))
  species = as.vector(sapply(species,str_trim))
  return(species)
}
writeSample = function(i,dir,reac,tags,drawPars,drawBR,type,maxReacts,maxProds) {
# Generate csv file of a database draw
  dbOut=data.frame()
  reacts = getSpecies(reac,maxReacts)
  signBeta=1
  if('E' %in% reactants) signBeta=-1
  for (ip in 1:length(tags)) {
    prods = getSpecies(tags[ip],maxProds)
    db1=data.frame(R1=reacts[1],R2=reacts[2],R3=reacts[3],
                   P1=prods[1],P2=prods[2],P3=prods[3],P4=prods[4],
                   a=drawPars['ALPHA']*drawBR[ip],
                   b=drawPars['BETA']*signBeta,
                   c=drawPars['GAMMA'],
                   f=1.0,g=0.0,
                   type=type) 
    dbOut=rbind(dbOut,db1)
  }
  dbOutm <- within(dbOut, {
    b <- sprintf("%6.3f", b)
    a  <- sprintf("%6.3e", a)
  })
  write.table(dbOutm,
              file=paste0(dir,'run_',sprintf('%04i',i),'.csv'),
              quote=TRUE, sep=';', na=" ",
              row.names=FALSE,col.names=FALSE)  
}  
writeSampleWL = function(i,dir,reac,tags,wl,xs,br,type,maxReacts,maxProds) {
  # Generate reacs file of a database random draw
  dbOut=data.frame()
  reacts = getSpecies(reac,maxReacts)
  signBeta=1
  if('E' %in% reactants) signBeta=-1
  for (ip in 1:length(tags)) {
    prods = getSpecies(tags[ip],maxProds)
    db1=data.frame(R1=reacts[1],R2=reacts[2],R3=reacts[3],
                   P1=prods[1],P2=prods[2],P3=prods[3],P4=prods[4],
                   # a=drawPars['ALPHA']*drawBR[ip],
                   # b=drawPars['BETA']*signBeta,
                   # c=drawPars['GAMMA'],
                   # f=1.0,g=0.0,
                   type=type) 
    dbOut=rbind(dbOut,db1)
  }
  # dbOutm <- within(dbOut, {
  #   b <- sprintf("%6.3f", b)
  #   a  <- sprintf("%6.3e", a)
  # })
  write.table(dbOut,
              file=paste0(dir,'run_',sprintf('%04i',i),'.csv'),
              quote=TRUE, sep=';', na=" ",
              row.names=FALSE,col.names=FALSE)  

  # Cross-section
  locTag = gsub(' ','_',reac)
  write.table(data.frame(wl=wl,xs=xs),
              file=paste0(dir,sprintf('%04i',i),'_se_',locTag,'.dat'),
              quote=TRUE, sep=' ',
              row.names=FALSE,col.names=FALSE) 

  # Branching ratios
  for (ic in 1:ncol(br)) {
    locTag = gsub(' ','_',tags[ic])
    write.table(data.frame(wl=wl,br=br[,ic]),
                file=paste0(dir,sprintf('%04i',i),'_qy_',locTag,'.dat'),
                quote=TRUE, sep=' ',
                row.names=FALSE,col.names=FALSE) 
  }
  
}  
getXSData = function(reac) {
  locTag = gsub(' ','_',reac)
  tmp = try(read.table(file=paste0(reac,'/se_',locTag,'.dat'), 
                       header=FALSE,stringsAsFactors = FALSE,
                       colClasses='numeric'),
            silent = FALSE
  )
  if (class(tmp) == "try-error")
      stop()

  return(list(wl=tmp[,1],xs=tmp[,2]))
}
getBRData = function(reac,tags,thresh=0.0001) {
  nb = length(tags)
  br = list(); wl = c()
  for(ib in 1:nb) {
    locTag = gsub(' ','_',tags[ib])
    tmp = try(read.table(file=paste0(reac,'/qy_',locTag,'.dat'), 
                         header=FALSE,stringsAsFactors = FALSE,
                         colClasses='numeric'),
              silent = FALSE
    )
    if (class(tmp) == "try-error")
      stop()
    else
      br[[ib]] = tmp
  }
  wl = tmp[,1]
  nl = length(wl)
  
  BR = matrix(NA,ncol=nb,nrow=nl)
  for(ib in 1:nb)
    BR[,ib] = br[[ib]][,2]
 
  # Threshold and normalize
  BR[BR<thresh] = thresh
  BR[(1-BR)<thresh] = 1-thresh
  BR=BR/rowSums(BR)
  
  rm(br,tmp)
  return(list(wl=wl,br=BR))
}
# Bibliography
printBib = function(keys,bibList) {
  if(length(keys) != 0) {
    cat('<H2>References</H2><DL>\n')
    for (key in keys) {
      cat(paste0('<DT>[',key,']</DT><DD>'))
      print(bibList[key],style="html")
      cat('</DD>')
    }
    cat('</DL>')
  }
}
printBibKeys = function(keys) {
  if(any(!is.na(keys))) {
    refs=paste0(sort(keys),collapse=',')
    paste0('[',refs,']')
  } else {
    NA
  }
}
printRQ = function(comments) {
  if(!is.na(comments)) {
    cat('<H2>Comments</H2>\n')
    cat(paste0('<font color="red">',comments,'</font>\n'))
  }
}


# Manage bibliography file to avoid repeated processing #####
setwd(docDir)
bibFile='refsDR.bib'
if(!file.exists('bib.Rdata')) {
  cat('*** Processing .bib file\n')
  bib = read.bib(file=bibFile)
  save(bib, file='bib.Rdata')
} else {
  sourceTime = file.info(bibFile)["mtime"]
  bibTime    = file.info('bib.Rdata')["mtime"]
  if(sourceTime > bibTime) {
    cat('*** Processing .bib file\n')
    bib = read.bib(file=bibFile)
    save(bib, file='bib.Rdata')    
  } else {
    cat('*** Loading  processed .bib file\n')
    load('bib.Rdata')
  }
}

# Parameters #####
maxReacts      = 3     # Max number of reactants slots in generated dBases
maxProds       = 4     # Max number of product slots in generated dBases
sampleSize     = 100   # Number of random samples to generate
writeSamples   = TRUE  # Output samples to disk (slow) or not (nominal sample is written)
tagged         = FALSE # Print tagged strings for debug
checkFragments = FALSE  # Partial run of script to check mass of fragments (no sampling)

dataDir = paste0(sourceDir,'Data/')
setwd(dataDir)
listReacs = list.dirs(full.names=FALSE, recursive=FALSE)
listReacs = gsub("./","",listReacs); cleanTmp=TRUE
# listReacs=c('N2 + HV','H2 + HV'); cleanTmp = TRUE

if (cleanTmp) {
  # Clean Tmp
  command=paste0('rm -rf ',targetDir,'Reactions/*') 
  system(command,intern=FALSE)
  command=paste0('rm -rf ',targetDir,'Species/*') 
  system(command,intern=FALSE)  
}

# Main loop over reactions #####
for (reac in listReacs) {
  allBibKeys = allSpecies = c()
  
  cat(paste0(reac,'\n'))
  
  reactants = getSpecies(reac,maxReacts)
  allSpecies = c(allSpecies,reactants)
  massReactants = getMassList(reactants)
  
  # Get data for this reaction #
  X = as.matrix(read.csv(paste0(reac,'/data.csv'),
                         header=FALSE, sep='\t', fill=TRUE,
                         stringsAsFactors = FALSE,
                         na.strings="")
  )
  X[] = str_trim(X) # remove unwanted spaces
  reacName = X[1,1]
  if(reacName != reac) stop(paste0('Pb reac identity:',reacName))
  
  # Locate Rate Info in dataFrame by keywords #
  
  # Reaction rate expression
  topLeft = which(X=='TYPE',arr.ind=TRUE)
  if(length(topLeft)==0) {
    reacType = 'kooij' # default
    if('E' %in% reactants) reacType ='dr'
    if('HV' %in% reactants) reacType ='photo'
  } else {  
    reacType = X[topLeft[1],topLeft[2]+1]
    if(! reacType %in% c('dr','kooij','ionpol1','ionpol2')) 
      stop(paste0('Improper rate type:',reacType))
  }
  
  # Rate parameters 
  rateParKwdList = c('ALPHA','BETA','GAMMA') 
  sampleRateParams = matrix(NA,nrow=sampleSize,ncol=length(rateParKwdList))
  colnames(sampleRateParams) = rateParKwdList
  rateParDistStrings = rep(NA,length(rateParKwdList))
  names(rateParDistStrings) = rateParKwdList
  for (kwd in rateParKwdList) {
    stringDist = getDistString(X,kwd)
    sampleDist = sampleDistString(stringDist,sampleSize)
    sampleRateParams[1:sampleSize,kwd]=sampleDist
    rateParDistStrings[kwd]=stringDist
  }
  
  # Locate BR Info #
  topLeft =  which(X=='BR',arr.ind=TRUE)
  X0=X[topLeft[1]:nrow(X),topLeft[2]:ncol(X)]
  tags=X0[2:nrow(X0),1] # List of branches
  
  # Check mass compatibility
  allProds = c()
  for (ip in 1:length(tags)) {
    prods=getSpecies(tags[ip],maxProds)
    allProds = c(allProds,prods)
    allSpecies = c(allSpecies,prods)
    massFrags = getMassList(prods)
    if(!is.na(massFrags)) {
      if(abs(massFrags-massReactants) > 0.01) {
        setwd("..")
        stop(paste0('Pb mass of fragments: ',tags[ip]))
      }      
    }
  }
  allProds = unique(allProds)
  
  # Build Prod-loss files
  for (sp in reactants) {
    specDir=paste0(targetDir,'Species/',sp)
    if(!file.exists(specDir)) dir.create(specDir)
    sink(file=paste0(specDir,'/loss.html'),append=TRUE)
    cat(paste0('<BR><A HREF="',dataDir,reac,'/summary.html">',reac,'</A>\n')) 
    sink(file=NULL)
  }
  for (sp in allProds) {
    specDir=paste0(targetDir,'Species/',sp)
    if(!file.exists(specDir)) dir.create(specDir)
    sink(file=paste0(specDir,'/prod.html'),append=TRUE)
    cat(paste0('<BR><A HREF="',dataDir,reac,'/summary.html">',reac,'</A>\n')) 
    sink(file=NULL)
  }
  
  if (checkFragments) next
  
  # Generate random XS ####
  XS = getXSData(reac)
  wl = XS$wl # Wavelength range
  xs = XS$xs # wl dependent cross-section
  sampleXS = matrix(NA,nrow=sampleSize,ncol=length(wl))
  for(iMC in 1:sampleSize) {
    alp = sampleRateParams[iMC,'ALPHA']
    sampleXS[iMC,] = xs * alp
  }
  # matplot(wl,t(sampleXS),type='l',lty=1,lwd=2,col=col_tr[1])
  # lines(wl,xs)
  
  if(length(tags) >=2) {
    
    # Get BR data
    BR = getBRData(reac,tags)
    if(!identical(BR$wl,wl)) 
      stop('Inconsistent wl range between XS and BR')
    br = BR$br # wl dependent BR
    
    # Build BR tree #####
    dist=X0[1,2:ncol(X0)] # List of distributions in tree
    dist=dist[!is.na(dist)]
    XT=X0[-1,-1]          # Matrix of parameters
    nc=1
    nl=nrow(XT)
    X1=matrix(XT[1:nl,1:nc],nrow=nl,ncol=nc)

    for (ic in 2:ncol(XT)) 
      if(sum(is.na(XT[,ic]))!=nl) {
        nc=nc+1
        X1=cbind(X1,XT[,ic])
      }
    
    d=list()
    dSub=list()
    for(ic in 1:nc){
      dSub$elem = which(X1[,ic] !='')
      dSub$dist = dist[ic]
      pars = as.vector(X1[dSub$elem,ic])
      dSub$mu = dSub$sig = c()
      for (ip in 1:length(pars)) {
        if( grepl("/",pars[ip]) ) {      
          loc  = gregexpr("(?<mu>.*)/(?<sig>.*)",pars,perl=TRUE)      
          start = attr(loc[[ip]],'capture.start')[1]
          stop  = start + attr(loc[[ip]],'capture.length')[1] -1
          dSub$mu[ip] = as.numeric(substr(pars[ip],start,stop))
          start = attr(loc[[ip]],'capture.start')[2]
          stop  = start + attr(loc[[ip]],'capture.length')[2] -1
          dSub$sig[ip] = as.numeric(substr(pars[ip],start,stop))  
        } else {
          dSub$mu[ip]  = pars[ip]
          dSub$sig[ip] = NA
        }    
      }
      dSub$link = rep(0,length(dSub$elem))
      if(ic!=1) {
        for (iel in 1:length(dSub$elem)) {
          ip = dSub$elem[iel]
          for (ic1 in 1:(ic-1)) {
            links = ip %in% d[[ic1]]$elem
            if(sum(links)!=0) dSub$link[iel] = ic1  
          }
        }
      }
      d[[ic]] = dSub
    }
    
    # Scheme for random errors #####    
#     dRef = d
#     matplot(wl,sampleBR,type='n',pch=19, col=col_tr,ylim=c(0,1))
#     for(iMC in 1:10) {
#   sampleBR = matrix(0,nrow=length(wl),ncol=ncol(br))    
#   # Loop over wavelengths #
#   for (iwl in 1:length(wl)) {
#     # Update d$mu with wl-dept values #
#     for (ic in 1:nc) {
#       elem = d[[ic]]$elem
#       if(sum(d[[ic]]$link)== 0){
#         d[[ic]]$sum = sum(br[iwl,elem])
#         d[[ic]]$mu  = br[iwl,elem] / d[[ic]]$sum
#         d[[ic]]$sig = dRef[[ic]]$sig * d[[ic]]$mu # Interpret as %
#         
#       } else {
#         for (ie in 1:length(elem)) {
#           if(d[[ic]]$link[ie] == 0) 
#             d[[ic]]$mu[ie] = br[iwl,ie]
#           else
#             d[[ic]]$mu[ie] = d[[ d[[ic]]$link[ie] ]]$sum
#         }
#         d[[ic]]$sum = sum(d[[ic]]$mu)
#         d[[ic]]$mu  = d[[ic]]$mu / d[[ic]]$sum
#         d[[ic]]$sig = dRef[[ic]]$sig * d[[ic]]$mu # Interpret as %
#       }
#     }
#     
#     # Build probabilistic tree string for sampler #
#     stringBR = oneDist(nc,d,tags,tagged)
#     while(grepl("LINK/",stringBR) ) {
#       poc = regmatches(stringBR,gregexpr('LINK/[0-9]+/',stringBR))
#       po  = sapply(poc[[1]],
#                    function(x) as.numeric(
#                      sub('LINK','',gsub('/','',x)) ))
#       for (ip in 1:length(po)) {
#         str = oneDist(po[ip],d,tags,tagged)
#         stringBR = sub(poc[[1]][ip],str,stringBR)
#       } 
#     }   
#     
#     # Generate BR sample #
#     # sampleBR = nds(sampleSize,stringBR)  
#     sampleBR[iwl,] = nds(1,stringBR)  
#   }
#   matplot(wl,sampleBR,type='p',pch=19, col=col_tr,add=TRUE)
# }
#     matplot(wl,br,type='l', col=cols,add=TRUE,lty=1,lwd=4)

    
    # Scheme for systematic errors #####   
    # Build probabilistic tree string for sampler #
    stringBR = oneDist(nc,d,tags,tagged)
    while(grepl("LINK/",stringBR) ) {
      poc = regmatches(stringBR,gregexpr('LINK/[0-9]+/',stringBR))
      po  = sapply(poc[[1]],
                   function(x) as.numeric(
                     sub('LINK','',gsub('/','',x)) ))
      for (ip in 1:length(po)) {
        str = oneDist(po[ip],d,tags,tagged)
        stringBR = sub(poc[[1]][ip],str,stringBR)
      } 
    }
    stringBR1 = stringBR

    # Build probabilistic tree string for sampler (no uncertainty) #
    d0 = d
    for (ic in 1:nc) 
      d0[[ic]]$sig = d0[[ic]]$sig/d0[[ic]]$sig *0.00001
    stringBR = oneDist(nc,d0,tags,tagged)
    while(grepl("LINK/",stringBR) ) {
      poc = regmatches(stringBR,gregexpr('LINK/[0-9]+/',stringBR))
      po  = sapply(poc[[1]],
                   function(x) as.numeric(
                     sub('LINK','',gsub('/','',x)) ))
      for (ip in 1:length(po)) {
        str = oneDist(po[ip],d,tags,tagged)
        stringBR = sub(poc[[1]][ip],str,stringBR)
      } 
    }   
    stringBR0 = stringBR
    rm(d0)
    stringBR = stringBR1 # For further treatment 
    
    # Power-perturbate reference BR #
    sampleBR = array(NA,dim=c(sampleSize,length(wl),length(tags)))
      # matplot(wl,br,type='n',ylim=c(0,1))
    sampleBR1 = nds(sampleSize,stringBR1)  
    sampleBR0 = nds(1,stringBR0)  
    for(iMC in 1:sampleSize) {
      pow = exp((sampleBR1[iMC,]-sampleBR0)/sampleBR0)
      br1 = t(apply(br,1,function(x) x^pow/sum(x^pow)))
      sampleBR[iMC,,] = br1
        # matplot(wl,br1,type='l',lty=1, lwd=3,col=col_tr,add=TRUE)
    }
      # matplot(wl,br,type='l', col=1,add=TRUE,lty=1,lwd=1)
      # legend('right',legend = tags, col=cols, lty=1,lwd=3)
    
    
    # Build Newick string for tree plotting #####
    newickBR = oneNewick(nc,d,tags)
    while(grepl("LINK/",newickBR) ) {
      poc=regmatches(newickBR,gregexpr('LINK/[0-9]+/',newickBR))
      po=sapply(poc[[1]],
                function(x) as.numeric(
                  sub('LINK','',gsub('/','',x)) ))
      for (ip in 1:length(po)) {
        str = oneNewick(po[ip],d,tags)
        newickBR = sub(poc[[1]][ip],str,newickBR)
      } 
    }  
    newickBR = paste0(newickBR,";")
    mytree <- read.tree(text=newickBR)
    
    # Build edge tags for tree annotation #
    edgeTags=oneEdgeTag(nc,d,tags)
    while(grepl("LINK/",edgeTags) ) {
      poc=regmatches(edgeTags,gregexpr('LINK/[0-9]+/',edgeTags))
      po=sapply(poc[[1]],
                function(x) as.numeric(
                  sub('LINK','',gsub('/','',x)) ))
      for (ip in 1:length(po)) {
        str = oneEdgeTag(po[ip],d,tags)
        edgeTags = sub(poc[[1]][ip],str,edgeTags)
      } 
    }
    edgeTags=unlist(strsplit(edgeTags,','))
    
    # Build node tags for tree annotation #
    nodeTags=oneNodeTag(nc,d,tags)
    while(grepl("LINK/",nodeTags) ) {
      poc=regmatches(nodeTags,gregexpr('LINK/[0-9]+/',nodeTags))
      po=sapply(poc[[1]],
                function(x) as.numeric(
                  sub('LINK','',gsub('/','',x)) ))
      for (ip in 1:length(po)) {
        str = oneNodeTag(po[ip],d,tags)
        nodeTags = sub(poc[[1]][ip],str,nodeTags)
      } 
    }
    nodeTags=unlist(strsplit(nodeTags,','))
    nodeTags=nodeTags[nodeTags!='NA']
    
    
    
  } else {
    
    # Single pathway with BR=1
    sampleBR = array(1,dim=c(sampleSize,length(wl),1))
    
  }
  
  # Generate output for kinetics code #####
  
  # Nominal/mean/median values from samples
  meanPars = rep(NA,ncol(sampleRateParams))
  names(meanPars)=colnames(sampleRateParams)
  sigPars = rep(NA,ncol(sampleRateParams))
  names(sigPars)=colnames(sampleRateParams)
  for (kwd in rateParKwdList) {
    sample = sampleRateParams[,kwd]
    meanPars[kwd] = exp(mean(log(sample)))
    sigPars[kwd]  = exp(sd(log(sample)))
  }
  rm(sample)
    
  # meanBR    = colMeans(sampleBR)
  # meanBR    = meanBR/sum(meanBR) 
  # sigBR     = apply(sampleBR,2,sd)
  
  #Generate kinetic databases
  samplesDir=paste0(targetDir,'Reactions/',reac)
  if(!file.exists(samplesDir)) dir.create(samplesDir)
  samplesDir=paste0(targetDir,'Reactions/',reac,'/Samples/')
  if(!file.exists(samplesDir)) dir.create(samplesDir)
  
  # Nominal run
  writeSampleWL(0, samplesDir, reac, tags, 
                wl, xs, br, 
                reacType, maxReacts, maxProds)     
  # Random runs
  if(writeSamples)
    for (i in 1:sampleSize) 
      writeSampleWL(i, samplesDir, reac, tags, 
                    wl, sampleXS[i,], sampleBR[i,,],
                    reacType, maxReacts, maxProds)

  
  
  # Plots #####
  cols = brewer.pal(8,'Dark2')
  col_tr=c()
  for (i in 1:length(cols))
    col_tr[i] = col2tr(cols[i],20)
  trBlue=col2tr('blue')
  np = min(sampleSize,100) # nb of plotted samples
  
  # Branching ratios
  nt=length(tags)
  if (nt >=2) {
    png(file=paste0(reac,'/figTreeBR.png'),
        width=max(800,nt*30),height=max(300,nt*30))
    tagStat=tags
    for (ip in 1:nt) {
      # tagStat[ip]=paste0(tags[ip],' (',
      #                    signif(meanBR[ip],2),' +/- ',
      #                    signif(sigBR[ip],1),')')   
      tagStat[ip] = tags[ip]
    }
    par(mar=c(1,1,1,1),cex=1)
    mytree$tip.label=tagStat
    mytree$edge.length=rep(1,dim(mytree$edge)[1])
    plot(mytree,type='clado',y.lim=c(length(tags),1),
         show.tip.label=TRUE,tip.color='blue',
         use.edge.length=TRUE,
         root.edge=TRUE, edge.width=2,
         edge.color='orange',
         font=1, main='')
    mynodelabels(nodeTags,bg='gold')
    myedgelabels(edgeTags)
    dev.off()

    png(file=paste0(reac,'/figBR.png'),width=600,height=400)
    par(mar=c(3,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5)
    matplot(wl,br,type='n',ylim=c(0,1),
            xlab='Wavelength / nm',
            ylab = 'Branching ratios', main=reac)
    grid(col='darkgray'); box()
    for(iMC in 1:np)
      matplot(wl,sampleBR[iMC,,],type='l',lty=1, lwd=3,col=col_tr,add=TRUE)
    matplot(wl,br,type='l', col=1,add=TRUE,lty=1,lwd=1)
    legend('topleft',legend = tags, col=cols, lty=1, lwd=3, bty='o')
    dev.off()

  }
  
  png(file=paste0(reac,'/figXS.png'),width=600,height=400)
  par(mar=c(3,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5)
  matplot(wl,xs,type='n', log='y',
          xlab='Wavelength / nm',
          ylab='Cross-Section', main=reac)
  grid(col='darkgray'); box()
  for(iMC in 1:np)
    matplot(wl,sampleXS[iMC,],type='l',lty=1, lwd=3,col=col_tr,add=TRUE)
  matplot(wl,xs,type='l', col=1,add=TRUE,lty=1,lwd=1)
  dev.off()
  
  # Generate summary document #####
  # Get Comments
  comments = getParams(X,'RQ')
  if(!is.na(comments)) 
    comments = paste0(comments,collapse=';')
  
  # Get bibliography
  bibKwd = paste0('REF_',c(rateParKwdList,'BR'))
  bibKeys=matrix(NA,nrow=ncol(X),ncol=length(bibKwd))
  colnames(bibKeys)=bibKwd
  for (kwd in bibKwd) {
    refs = getParams(X,kwd)
    if(sum(!is.na(refs))!=0) bibKeys[1:length(refs),kwd]=refs
  }
  
  sink(file=paste0(reac,'/summary.html'),append=FALSE)
  cat(paste0('<H1>',reacName,'</H1>\n'))
  
  printRQ(comments)
 
  cat('<H2>Global rate constant</H2>\n')
  
  cat(paste0('<H3>Cross Section: ',reacType,'</H3>\n'))

  parTable=matrix(NA,ncol=length(rateParKwdList),nrow=4)
  colnames(parTable)=c(capwords(rateParKwdList))
  rownames(parTable)=c('Dist.','Median','F','Refs')
  icol=0
  for (kwd in rateParKwdList) {
    icol=icol+1
    parTable[1,icol]=rateParDistStrings[kwd]
    if(is.finite(sigPars[kwd])){
      parTable[2,icol]=signif(meanPars[kwd],2)
      parTable[3,icol]=signif(sigPars[kwd] ,2)
    }
    refKey=paste0('REF_',kwd)
    parTable[4,icol]=printBibKeys(bibKeys[,refKey])
  }
   print(xtable(parTable), type='html', 
         html.table.attributes = "border = '1', align = 'center', width = 600",
         include.rownames = TRUE, include.colnames = TRUE)
  
   cat('<table border=1 align="center">\n')
   cat('<tr><th><I>Cross section sample</I></th></tr>')
   cat('<tr><td><img src="./figXS.png" width=600></td></tr>')
   cat('</table>\n')
  
  if( length(tags) >= 2 ) cat('<!-- PAGE BREAK -->\n')
  cat('<H2>Branching Ratios</H2>\n')
  refs=printBibKeys(bibKeys[,'REF_BR'])
  if(!is.na(refs)) cat(paste0('<B>Refs: ',refs,'</B>\n'))
  
  if( length(tags) < 2 ) {
    cat('<H3>Data</H3>\n')
    print(xtable(X0),type="html",
          html.table.attributes = "border = '1', align = 'center', width = 600",
          include.rownames = FALSE, include.colnames = FALSE)   
  } else {
    cat('<P><table border=1 align="center">\n')
    cat('<tr><th><I>Probabilistic tree of branching ratios</I></th></tr>')
    cat('<tr><td><img src="./figTreeBR.png" width=600></td></tr>')
    cat('</table>\n')
    
    if( length(tags) >= 10 )cat('<!-- PAGE BREAK -->\n')
    cat('<table border=1 align="center">\n')
    cat('<tr><th><I>Branching ratios sample</I></th></tr>')
    cat('<tr><td><img src="./figBR.png" width=600></td></tr>')
    cat('</table>\n')
    
  }

  thisBibKeys=sort(unique(as.vector(bibKeys)))
  if(length(thisBibKeys)>0) {
    printBib(thisBibKeys,bib)
    allBibKeys = c(allBibKeys,thisBibKeys) 
  }

  sink(file = NULL)

#   # Generate PDF
#   setwd(paste0(dataDir,reac))
#   targetPDF = paste0(targetDir,reac,'/summary.pdf')
#   command='pandoc -V geometry:margin=2cm summary.html -o summary.pdf' 
#   system(command,intern=FALSE)
#   file.rename(from='summary.pdf',to=targetPDF)
#   setwd(dataDir)

  sink(file=paste0(targetDir,'Reactions/',reac,'/species.txt'),append=FALSE)
  cat(unique(allSpecies))
  sink(file = NULL)

  sink(file=paste0(targetDir,'Reactions/',reac,'/bibKeys.txt'),append=FALSE)
  cat(unique(allBibKeys))
  sink(file = NULL)

#   sink(file=paste0(targetDir,'Reactions/',reac,'/dataTable.html'),append=FALSE)
#   cat('<TR>\n')
#   cat(paste0('<TD>',reac,'</TD>\n'))
#   cat(paste0('<TD>--></TD>\n'))
#   cat(paste0('<TD>',tags[1],'</TD>\n'))
#   if(length(tags) <2) {
#     mbr = 1
#     sbr = 0
#   } else{
#     mbr = meanBR[1]
#     sbr = sigBR[1]    
#   }
#   cat(paste0('<TD>',paste0(signif(mbr,2),' +/- ',
#                            signif(sbr,1))    ,'</TD>\n'))
#   cat(paste0('<TD>',paste0(signif(meanPars['ALPHA'],2),' */ ',
#                            signif(sigPars['ALPHA'],2))    ,'</TD>\n'))
#   cat('</TR>\n')
#   if(length(tags) >=2) {
#     for(i in 2:length(tags)) {
# #       cat('<TR><TD COLSPAN=5>HR</TD></TR>\n')
#       cat('<TR>\n')
#       cat(paste0('<TD> </TD>\n'))
#       cat(paste0('<TD> </TD>\n'))
#       cat(paste0('<TD>',tags[i],'</TD>\n'))
#       cat(paste0('<TD>',paste0(signif(meanBR[i],2),' +/- ',
#                                signif(sigBR[i],1))    ,'</TD>\n'))
#       cat(paste0('<TD> </TD>\n'))
#       cat('</TR>\n')      
#     }
#   }
#       
#   sink(file = NULL)

}
