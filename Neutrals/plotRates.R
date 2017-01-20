rm(list = ls())

version   = '1.0'
sourceDir = paste0('/home/pernot/Bureau/ChemDB/Neutrals/Source/v_',version,'/')
publicDir = paste0('/home/pernot/Bureau/ChemDB/Neutrals/Public/v_',version,'/')
samplesDir= paste0(publicDir,'Databases/')
figsDir   = paste0(publicDir,'Figs/')
tmpDir    = paste0('/home/pernot/Bureau/ChemDB/Neutrals/Tmp/')

# T, P
T0 = 300
tempRange = seq(50,350,by=5)
pressure  = 1e18 # molec/cm^3
  
# Clean target dirs
command = paste0('rm -rf ',tmpDir,'*')
dummy = system(command)

# List sample files
samplesList=list.files(path=paste0(publicDir,'Databases'),pattern='run_',full.name=TRUE)

# Generate individual curves
irun=-1
for (file in samplesList) {
  irun=irun+1

  nbReac=0  
  reactants = products = params = type = tag = list()
  scheme  = read.csv(file=file,header=FALSE,sep=';')
  scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
  for (i in 1:nrow(scheme)) {
    nbReac = nbReac + 1
    terms=scheme[i,1:3]
    reactants[[nbReac]] = terms[!is.na(terms) & terms!=""]
    terms=scheme[i,4:7]
    products[[nbReac]]  = terms[!is.na(terms) & terms!=""]
    terms=scheme[i,8:18]
    params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
    tag[[nbReac]]       = paste0(nbReac,': ',
                                 paste0(reactants[[nbReac]],collapse='+'),
                                 '->',
                                 paste0(products[[nbReac]],collapse='+'))
  }
  
  alerts = c()
  for (i in 1:nbReac) {
    tdir=paste0(tmpDir,tag[[i]])
    if(!file.exists(tdir)) dir.create(tdir)
    type=params[[i]][length(params[[i]])]
    
    if(type == 'kooij') {
      cf = as.numeric(params[[i]][1:5])
      k = cf[1] * (tempRange/T0)**cf[2] * exp(-cf[3]/tempRange) *
        cf[4] * exp( cf[5]*abs(1/tempRange-1/T0) )  
      if(sum(k<=0) != 0) alerts=c(alerts,paste0('Null RC: ',tag[[i]],'\n'))
      
    } else { # 3-body
      cf = as.numeric(params[[i]][1:5])
      k0 = cf[1] * (tempRange/T0)**cf[2] * exp(-cf[3]/tempRange) *
        cf[4] * exp( cf[5]*abs(1/tempRange-1/T0) )
      cf = as.numeric(params[[i]][6:10])
      kInf = cf[1] * (tempRange/T0)**cf[2] * exp(-cf[3]/tempRange) *
        cf[4] * exp( cf[5]*abs(1/tempRange-1/T0) )
      fc = 0.64
      Pr = k0*pressure/kInf
      cExp = -0.4 -0.67*log10(fc)
      NExp = 0.75 -1.27*log10(fc)
      dExp = 0.14
      fExp = 1+((log10(Pr)+cExp)/(NExp-dExp*(log10(Pr)+cExp)))^2
      broadF = fc^(1/fExp)
      k = kInf*(Pr/(1+Pr))*broadF    
      if(sum(k<=0) != 0) alerts=c(alerts,paste0('Null RC: ',tag[[i]],'\n'))
      
    }
    write.table(data.frame(T=tempRange,kval=k),
                file=paste0(tdir,'/curve_',sprintf('%04i',irun),'.csv'),
                row.names=FALSE,col.names=FALSE)
  }
}
if(length(alerts != 0)) print(unique(alerts))

# Gather info from source database for legending plots

# Kinetic Parser ############################################################
nbReac0=0  
reactants = products = params = type = orig = legText = list()
filename='Titan - Réactions bimoléculaires.csv'
scheme  = read.csv(file=paste0(sourceDir,filename),header=FALSE,
                   stringsAsFactors = FALSE)
comments = scheme[,ncol(scheme)]
scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
for (i in 1:nrow(scheme)) {
  nbReac0 = nbReac0 + 1
  terms=scheme[i,1:3]
  reactants[[nbReac0]] = terms[!is.na(terms) & terms!=""]
  terms=scheme[i,4:8]
  products[[nbReac0]]  = terms[!is.na(terms) & terms!=""]
  terms=scheme[i,9:13]
  params[[nbReac0]]    = terms[!is.na(terms) & terms!=""]
  type[[nbReac0]]      = 'kooij'
  orig[[nbReac0]]      = filename
}

filename = 'Titan - Réactions trimoléculaires.csv'
scheme  = read.csv(file=paste0(sourceDir,filename),header=FALSE,
                   stringsAsFactors = FALSE)
comments = c(comments,scheme[,ncol(scheme)])
scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
for (i in 1:nrow(scheme)) {
  nbReac0 = nbReac0 + 1
  terms=scheme[i,1:3]
  reactants[[nbReac0]] = terms[!is.na(terms) & terms!=""]
  terms=scheme[i,4:8]
  products[[nbReac0]]  = terms[!is.na(terms) & terms!=""]
  terms=scheme[i,9:18]
  params[[nbReac0]]    = terms[!is.na(terms) & terms!=""]
  type[[nbReac0]]      = '3-body'
  orig[[nbReac0]]      = filename
}

if(nbReac0 != nbReac) {cat('Pb. of database version\n'); stop()}

# Generate plots from curves
fileList = list.files(path=figsDir,full.name=TRUE)
dummy = file.remove(fileList)

col2tr =function(x,alpha=80){
  rgb(unlist(t(col2rgb(x))),alpha=alpha,maxColorValue = 255)
}
trBlue=col2tr('blue',60)
for (ireac in 1:nbReac) {
  reac=tag[[ireac]]
  
  legText = paste0(tag,'\n',
                   'Rate law: ',type[[ireac]],'\n')
  if(type[[ireac]] == 'kooij')
    legText=paste0(legText,
                   'Parameters: ',paste0(params[[ireac]][1:5],collapse=' / '),'\n')
  else
    legText=paste0(legText,
                   'Parameters k0  : ',paste0(params[[ireac]][1:5],collapse=' / '),'\n',
                   'Parameters kInf: ',paste0(params[[ireac]][6:10],collapse=' / '),'\n')
  com = comments[ireac]
  if(!is.na(com)) {
    splCom = unlist(strsplit(com,' '))
    while (length(splCom)!=0) {
      sel = cumsum(nchar(splCom)) <= 65
      line= paste0(splCom[sel],collapse=' ')
      legText = paste0(legText,'! ',line,'\n')
      splCom = splCom[!sel]
    }
  }
  
  tdir=paste0(tmpDir,reac)
  fileList=list.files(path=tdir,full.name=TRUE)
  
  ktab=matrix(0,nrow=length(tempRange),ncol=length(fileList))
  for(i in seq_along(fileList)) 
    ktab[,i]=read.csv(fileList[i],header=FALSE,sep=' ')[,2]

  if(diff(range(ktab))!=0) {
    png(file=paste0(figsDir,reac,'.png'),width=1000,height=1000)
    par(mar=c(4,5,20,1),cex.lab=2,cex.axis=2,cex.main=3)
    matplot(tempRange,ktab,type='l',lty=1,col=trBlue,lwd=3,log='y',
            xlab='T / K', ylab='rate ct. / cm^3.s^-1',
            ylim=c(1e-28,1e-7),
            #           ylim=c(min(ktab)/1.5,max(ktab)*1.5),
            main='' )
    lines(tempRange,ktab[,1],col='red',lwd=3)
    grid(col='darkgray')
    mtext(legText[[ireac]],side=3,cex=2,
          adj=0,line=18,padj=1,
          col='darkgreen')
    box(lwd=4)
    dev.off()      
  }

}
