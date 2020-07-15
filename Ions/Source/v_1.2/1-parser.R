version   = '1.2'
rootDir   = '/home/pernot/Bureau/Titan-APSIS/MC-ChemDB/'
sourceDir = paste0(rootDir,'Ions/Source/v_',version,'/')
tmpDir    = paste0(rootDir,'Ions/Tmp/v_',version,'/')
publicDir = paste0(rootDir,'Ions/Public/v_',version,'/')
docDir    = paste0(sourceDir,'Doc/')
targetDir = tmpDir

setwd(sourceDir)

# Load libraries #####
libs =c('stringr','xtable','bibtex','ape','CHNOSZ')
for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

# Stoechiometry functions
source('./massCalc.R')

# Misc functions #####
# Modif of ape function to enlarge "rect" frames in labels
mynodelabels=function (text, node, adj = c(0.5, 0.5), frame = "rect", pch = NULL, 
          thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
          bg = "lightblue", horiz = FALSE, width = NULL, height = NULL, 
          ...) 
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (missing(node)) 
    node <- (lastPP$Ntip + 1):length(lastPP$xx)
  XX <- lastPP$xx[node]
  YY <- lastPP$yy[node]
  myBOTHlabels(text, node, XX, YY, adj, frame, pch, thermo, pie, 
             piecol, col, bg, horiz, width, height, ...)
}
myedgelabels=function (text, edge, adj = c(0.5, 0.5), frame = "rect", pch = NULL, 
          thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
          bg = "lightgreen", horiz = FALSE, width = NULL, height = NULL, 
          date = NULL, ...) 
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (missing(edge)) {
    sel <- 1:dim(lastPP$edge)[1]
    subedge <- lastPP$edge
  }
  else {
    sel <- edge
    subedge <- lastPP$edge[sel, , drop = FALSE]
  }
  if (lastPP$type == "phylogram") {
    if (lastPP$direction %in% c("rightwards", "leftwards")) {
      XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[, 
                                                         2]])/2
      YY <- lastPP$yy[subedge[, 2]]
    }
    else {
      XX <- lastPP$xx[subedge[, 2]]
      YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[, 
                                                         2]])/2
    }
  }
  else {
    XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[, 
                                                       2]])/2
    YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[, 
                                                       2]])/2
  }
  if (!is.null(date)) 
    XX[] <- max(lastPP$xx) - date
  myBOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo, pie, 
             piecol, col, bg, horiz, width, height, ...)
}
myBOTHlabels = function (text, sel, XX, YY, adj, frame, pch, thermo, pie, piecol, 
          col, bg, horiz, width, height, ...) 
{
  if (missing(text)) 
    text <- NULL
  if (length(adj) == 1) 
    adj <- c(adj, 0.5)
  if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie)) 
    text <- as.character(sel)
  frame <- match.arg(frame, c("rect", "circle", "none"))
  args <- list(...)
  CEX <- if ("cex" %in% names(args)) 
    args$cex
  else par("cex")
  if (frame != "none" && !is.null(text)) {
    if (frame == "rect") {
      width <- strwidth(text, units = "inches", cex = CEX)
      height <- strheight(text, units = "inches", cex = CEX)
      if ("srt" %in% names(args)) {
        args$srt <- args$srt%%360
        if (args$srt == 90 || args$srt == 270) {
          tmp <- width
          width <- height
          height <- tmp
        }
        else if (args$srt != 0) 
          warning("only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n")
      }
      width <- xinch(width)
      height <- yinch(height)
      xl <- XX - width * adj[1] - xinch(0.04)
      xr <- xl + width + xinch(0.08)
      yb <- YY - height * adj[2] - yinch(0.04)
      yt <- yb + height + yinch(0.08)
      rect(xl, yb, xr, yt, col = bg)
    }
    if (frame == "circle") {
      radii <- 0.8 * apply(cbind(strheight(text, units = "inches", 
                                           cex = CEX), strwidth(text, units = "inches", 
                                                                cex = CEX)), 1, max)
      symbols(XX, YY, circles = radii, inches = max(radii), 
              add = TRUE, bg = bg)
    }
  }
  if (!is.null(thermo)) {
    parusr <- par("usr")
    if (is.null(width)) {
      width <- CEX * (parusr[2] - parusr[1])
      width <- if (horiz) 
        width/15
      else width/40
    }
    if (is.null(height)) {
      height <- CEX * (parusr[4] - parusr[3])
      height <- if (horiz) 
        height/40
      else height/15
    }
    if (is.vector(thermo)) 
      thermo <- cbind(thermo, 1 - thermo)
    thermo <- if (horiz) 
      width * thermo
    else height * thermo
    if (is.null(piecol)) 
      piecol <- rainbow(ncol(thermo))
    xl <- XX - width/2 + adj[1] - 0.5
    xr <- xl + width
    yb <- YY - height/2 + adj[2] - 0.5
    yt <- yb + height
    if (horiz) {
      rect(xl, yb, xl + thermo[, 1], yt, border = NA, col = piecol[1])
      for (i in 2:ncol(thermo)) rect(xl + rowSums(thermo[, 
                                                         1:(i - 1), drop = FALSE]), yb, xl + rowSums(thermo[, 
                                                                                                            1:i]), yt, border = NA, col = piecol[i])
    }
    else {
      rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
      for (i in 2:ncol(thermo)) rect(xl, yb + rowSums(thermo[, 
                                                             1:(i - 1), drop = FALSE]), xr, yb + rowSums(thermo[, 
                                                                                                                1:i]), border = NA, col = piecol[i])
    }
    s <- apply(thermo, 1, function(xx) any(is.na(xx)))
    xl[s] <- xr[s] <- NA
    rect(xl, yb, xr, yt, border = "black")
    if (!horiz) {
      segments(xl, YY, xl - width/5, YY)
      segments(xr, YY, xr + width/5, YY)
    }
  }
  if (!is.null(pie)) {
    if (is.vector(pie)) 
      pie <- cbind(pie, 1 - pie)
    xrad <- CEX * diff(par("usr")[1:2])/50
    xrad <- rep(xrad, length(sel))
    XX <- XX + adj[1] - 0.5
    YY <- YY + adj[2] - 0.5
    for (i in seq_along(sel)) {
      if (any(is.na(pie[i, ]))) 
        next
      floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], 
                       col = piecol)
    }
  }
  if (!is.null(text)) 
    text(XX, YY, text, adj = adj, col = col, ...)
  if (!is.null(pch)) 
    points(XX + adj[1] - 0.5, YY + adj[2] - 0.5, pch = pch, 
           col = col, bg = bg, ...)
}

# Additional functions
capwords = function(s, strict = TRUE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s = substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
sanitize = function(str) {
  gsub('([#$%&~_\\^\\\\{}\\s\\(\\)])', '\\\\\\1', str, perl = TRUE)
}
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
col2tr =function(x,alpha=80){
  rgb(unlist(t(col2rgb(x))),alpha=alpha,maxColorValue = 255)
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
sampleSize     = 500   # Number of random samples to generate
writeSamples   = TRUE  # Output samples to disk (slow) or not (nominal sample)
tagged         = FALSE  # Print tagged strings for debug
checkFragments = FALSE  # Partial run of script to check mass of fragments (no sampling)
cleanTmp       = TRUE  # Clean all tmp files (full update)


dataDir = paste0(sourceDir,'Data/'); setwd(dataDir)
listReacs = list.dirs(full.names=FALSE, recursive=FALSE)
listReacs = gsub("./","",listReacs)

# # Select by date (same day)
# finf <- file.info(listReacs, extra_cols = FALSE)
# selByDate = which(difftime(Sys.time(), finf[,"mtime"], units = "days") <= 1 )
# sel = selByDate
# 
# # # Select files where data.csv is newer than summary.html
# # finf1 <- file.info(paste0(listReacs,'/data.csv'), extra_cols = FALSE)
# # finf2 <- file.info(paste0(listReacs,'/summary.html'), extra_cols = FALSE)
# # selByMod = which(difftime(finf1[,"mtime"],finf2[,"mtime"]) > 0)
# # sel = unique(c(selByDate,selByMod))
# 
# listReacs = listReacs[sel]
# cleanTmp=FALSE

#listReacs=c('H3+ + E') ; cleanTmp=FALSE #################

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
  massReactants = getMassList(reactants, excludeList = dummySpecies)
  
  # Get data for this reaction #
  fName = paste0(reac,'/data.csv')
  X = as.matrix(read.csv(fName,
                         header=FALSE, sep='\t', fill=TRUE,
                         stringsAsFactors = FALSE,
                         na.strings="")
  )
  X = trimws(X) # remove unwanted spaces
  reacName = X[1,1]
  if(reacName != reac) stop(paste0('Pb reac identity:',reacName))
  
  lastMod = file.info(paste0(dataDir,fName))$mtime
  
  # Locate Rate Info in dataFrame by keywords #
  
  # Reaction rate expression
  topLeft = which(X=='TYPE',arr.ind=TRUE)
  if(length(topLeft)==0) {
    reacType = 'kooij' # default
    if('E' %in% reactants) reacType ='dr'
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
    massFrags = getMassList(prods, excludeList = dummySpecies)
    if(!is.na(massReactants) & !is.na(massFrags)) {
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
  
  if(length(tags) >=2) {
    # Build tree #####
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
          dSub$mu[ip] = substr(pars[ip],start,stop)
          start = attr(loc[[ip]],'capture.start')[2]
          stop  = start + attr(loc[[ip]],'capture.length')[2] -1
          dSub$sig[ip] = substr(pars[ip],start,stop)  
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
    
    # Build probabilistic tree string for sampler #
    stringBR = oneDist(nc, d, tags, tagged)
    while (grepl("LINK/", stringBR)) {
      poc = regmatches(stringBR, gregexpr('LINK/[0-9]+/', stringBR))
      po = sapply(poc[[1]],
                  function(x)
                    as.numeric(sub('LINK', '', gsub('/', '', x))))
      for (ip in 1:length(po)) {
        str = oneDist(po[ip], d, tags, tagged)
        stringBR = sub(poc[[1]][ip], str, stringBR)
      }
    }    
    # Generate BR sample #
    sampleBR = nds(sampleSize,stringBR)  

    # Build Newick string for tree plotting #####
    newickBR = oneNewick(nc, d, tags)
    while (grepl("LINK/", newickBR)) {
      poc = regmatches(newickBR, gregexpr('LINK/[0-9]+/', newickBR))
      po = sapply(poc[[1]],
                  function(x)
                    as.numeric(sub('LINK', '', gsub('/', '', x))))
      for (ip in 1:length(po)) {
        str = oneNewick(po[ip], d, tags)
        newickBR = sub(poc[[1]][ip], str, newickBR)
      }
    }
    newickBR = paste0(newickBR, ";")
    mytree <- read.tree(text = newickBR)
    
    # Build edge tags for tree annotation #
    edgeTags = oneEdgeTag(nc, d, tags)
    while (grepl("LINK/", edgeTags)) {
      poc = regmatches(edgeTags, gregexpr('LINK/[0-9]+/', edgeTags))
      po = sapply(poc[[1]],
                  function(x)
                    as.numeric(sub('LINK', '', gsub('/', '', x))))
      for (ip in 1:length(po)) {
        str = oneEdgeTag(po[ip], d, tags)
        edgeTags = sub(poc[[1]][ip], str, edgeTags)
      }
    }
    edgeTags = unlist(strsplit(edgeTags, ','))
    
    # Build node tags for tree annotation #
    nodeTags = oneNodeTag(nc, d, tags)
    while (grepl("LINK/", nodeTags)) {
      poc = regmatches(nodeTags, gregexpr('LINK/[0-9]+/', nodeTags))
      po = sapply(poc[[1]],
                  function(x)
                    as.numeric(sub('LINK', '', gsub('/', '', x))))
      for (ip in 1:length(po)) {
        str = oneNodeTag(po[ip], d, tags)
        nodeTags = sub(poc[[1]][ip], str, nodeTags)
      }
    }
    nodeTags = unlist(strsplit(nodeTags, ','))
    nodeTags = nodeTags[nodeTags != 'NA']
    
  } else {
    # Single pathway with BR=1
    sampleBR = matrix(1,ncol=1,nrow=sampleSize)
    
  }
  
  # Generate output for kinetics code #
  
  # Nominal/mean/median values from samples
  meanPars = rep(NA, ncol(sampleRateParams))
  names(meanPars) = colnames(sampleRateParams)
  sigPars = rep(NA, ncol(sampleRateParams))
  names(sigPars) = colnames(sampleRateParams)
  for (kwd in rateParKwdList) {
    sample = sampleRateParams[, kwd]
    if(substr(rateParDistStrings[kwd],1,5)!='Delta') {
      meanPars[kwd] = exp(mean(log(sample)))
      sigPars[kwd]  = exp(sd(log(sample)))
    } else {
      meanPars[kwd] = mean(sample)
      sigPars[kwd]  = 1
    }
  }
  rm(sample)
  
  meanBR    = colMeans(sampleBR)
  meanBR    = meanBR / sum(meanBR)
  sigBR     = apply(sampleBR, 2, sd)
  
  #Generate kinetic databases
  samplesDir=paste0(targetDir,'Reactions/',reac)
  if(!file.exists(samplesDir)) dir.create(samplesDir)
  samplesDir=paste0(targetDir,'Reactions/',reac,'/Samples/')
  if(!file.exists(samplesDir)) dir.create(samplesDir)
  
  # Nominal run
  writeSample(0, samplesDir, reac, tags, meanPars, meanBR, 
              reacType, maxReacts, maxProds)     
  # Random runs
  if(writeSamples)
    for (i in 1:sampleSize) writeSample(i, samplesDir, reac, tags, 
                                        sampleRateParams[i,], sampleBR[i,], 
                                        reacType, maxReacts, maxProds)     
  
  # Plots #####
  trBlue=col2tr('blue')
  np = min(sampleSize,100) # nb of plotted samples
  # Branching ratios
  nt=length(tags)
  if (nt >=2) {
    png(file=paste0(reac,'/figTreeBR.png'),width=max(800,nt*30),height=max(300,nt*30))
    tagStat=tags
    for (ip in 1:nt) {
      tagStat[ip]=paste0(tags[ip],' (',
                         signif(meanBR[ip],2),' +/- ',
                         signif(sigBR[ip],1),')')   
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
         
    png(file=paste0(reac,'/figBR.png'),width=600,height=min(max(200,nt*20),800))
    par(mar=c(1,10,4,1),cex.main=1)
    matplot(t(sampleBR)[nt:1,1:np],1:nt,
            main='',
            type='l',lty=1,col=trBlue,lwd=2,
            yaxt='n',ylab='',yaxs='i',
            xlim=c(0,1),xlab='',xaxt='n',xaxs='i')
    lines(meanBR[nt:1],1:nt,lwd=3,col='red')
    # lines(1:length(tags),cumsum(vMean),lwd=3,col='green')
    text(x=seq(0,1,by=0.1),y=nt,labels=seq(0,1,by=0.1),
         xpd=TRUE,cex=1,pos=3)
    text(x=-0.02,y=nt:1,labels=tags,srt=0,adj=1,xpd=TRUE,cex=1)
    grid(col='darkgray')
    abline(h=1:length(tags),col='darkgray',lwd=3,lty=3)    
    #     text(meanBR,nt:1,labels=tagStat,srt=0,adj=1,pos=4,xpd=TRUE,cex=1)
    dev.off()
  }
  
  # Rate
  png(file=paste0(reac,'/figRate.png'),width=800,height=600)
  
  split.screen(c(2, 1))
  split.screen(c(1,3),1)
  iscreen=2
  for (kwd in rateParKwdList) {
    iscreen = iscreen+1; screen(iscreen)
    if(!is.finite(sigPars[kwd]) | sigPars[kwd] == 1 ) next
    par(mar=c(4,5,1,1)) 
    hist(sampleRateParams[,kwd],
         xlab=paste0(kwd,'~',rateParDistStrings[kwd]),
         col=trBlue,main='')
    abline(v=meanPars[kwd],col='red',lwd=2)
    abline(v=range(sampleRateParams[,kwd]),col='red',lty=2,lwd=1)
  }
  
  screen(2)
  par(mar=c(4,5,1,1))
  temp=seq(150,1000,by=10)
  nt=length(temp)
  if(reacType =='kooij') {
    if('E' %in% reactants) { 
      rateFun = function(t,pars) 
        pars['ALPHA']*(300/t)^pars['BETA']
    } else {
      rateFun = function(t,pars) 
        pars['ALPHA']*(t/300)^pars['BETA']*exp(-pars['GAMMA']/t)          
    }
  } else { 
    if(reacType =='ionpol1') {
      # KIDA::ionpol1 type 
      # BEWARE: notation change because here branching ratios are stored in BR
      rateFun = function(t,pars) 
        pars['ALPHA']*(0.62 + 
                       0.4767*pars['BETA']*(300/t)^0.5
                      )          
    } else {
      # KIDA::ionpol2 type 
      # BEWARE: notation change because here branching ratios are stored in BR
      rateFun = function(t,pars) 
        pars['ALPHA']*(1 + 
                       0.0967*pars['BETA']*(300/t)^0.5 +
                       pars['BETA']^2*300/(10.526*t)
                      )          
    }
  }
  
  Y=matrix(NA,nrow=np,ncol=nt)
  for (ip in 1:np) 
    Y[ip,1:nt]=rateFun(temp,sampleRateParams[ip,])
  matplot(temp,t(Y),type='l',lty=1,col=trBlue,lwd=2,log='y',
          xlab='T / K', ylab='rate ct. / cm^3.s^-1',
          ylim=c(min(Y)/1.5,max(Y)*1.5),
          main='' )
  lines(temp,rateFun(temp,meanPars),col='red',lwd=2)
  grid(col='darkgray')
  close.screen(all = TRUE)
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
  
  cat(paste0('<B>Last mod. : </B>',lastMod,'\n'))
  
  printRQ(comments)
 
  cat('<H2>Global rate constant</H2>\n')
  
  cat(paste0('<H3>Rate Formula: ',reacType,'</H3>\n'))

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
  
  cat('<P><table border=1 width=600 align="center">\n')
  cat('<tr><th><I>(Top) Histograms of the samples rate parameters;
      (Bottom) Samples of the T[Te]-dependent rate curves</I></th></tr>')
  cat('<tr><td><img src="./figRate.png" width=600></td></tr>')
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
    
#     cat('<H3>Sample stats</H3>\n')
#     tabBR=data.frame(Fragments=NULL,Stats=NULL,CumSum=NULL)
#     cSum=cumsum(meanBR)
#     tagStat=tags
#     for (ip in 1:length(tags)) {
#       tagStat[ip]=paste0(signif(meanBR[ip],2),' +/- ',
#                          signif(sigBR[ip],1))   
#       tabBR=rbind(tabBR,data.frame(Fragments=tags[ip],Stats=tagStat[ip],
#                                    CumSum=signif(cSum[ip],2)))
#     }
#     print(xtable(tabBR), type='html', 
#           html.table.attributes = "border = '1', align = 'center', width = 600",
#           include.rownames = FALSE, include.colnames = TRUE)

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

  sink(file=paste0(targetDir,'Reactions/',reac,'/dataTable.html'),append=FALSE)
  cat('<TR>\n')
  cat(paste0('<TD>',reac,'</TD>\n'))
  cat(paste0('<TD>--></TD>\n'))
  cat(paste0('<TD>',tags[1],'</TD>\n'))
  if(length(tags) <2) {
    mbr = 1
    sbr = 0
  } else{
    mbr = meanBR[1]
    sbr = sigBR[1]    
  }
  cat(paste0('<TD>',paste0(signif(mbr,2),' +/- ',
                           signif(sbr,1))    ,'</TD>\n'))
  cat(paste0('<TD>',paste0(signif(meanPars['ALPHA'],2),' */ ',
                           signif(sigPars['ALPHA'],2))    ,'</TD>\n'))
  cat('</TR>\n')
  if(length(tags) >=2) {
    for(i in 2:length(tags)) {
#       cat('<TR><TD COLSPAN=5>HR</TD></TR>\n')
      cat('<TR>\n')
      cat(paste0('<TD> </TD>\n'))
      cat(paste0('<TD> </TD>\n'))
      cat(paste0('<TD>',tags[i],'</TD>\n'))
      cat(paste0('<TD>',paste0(signif(meanBR[i],2),' +/- ',
                               signif(sigBR[i],1))    ,'</TD>\n'))
      cat(paste0('<TD> </TD>\n'))
      cat('</TR>\n')      
    }
  }
      
  sink(file = NULL)

}
