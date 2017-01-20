setwd("~/Bureau/Reactor/General/Photo")
options(warn=2)

nbReac=0  
reactants = products = params = type = list()
filename='../PhotoScheme.dat'
scheme  = read.fwf(file=filename, widths= rep(11,12))
scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
for (i in 1:nrow(scheme)) {
  nbReac = nbReac + 1
  terms=scheme[i,1:2]
  reactants[[nbReac]] = terms[!is.na(terms) & terms!="" & terms!="HV"]
  terms=scheme[i,3:6]
  products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
  params[[nbReac]]    = scheme[i,7]
  type[[nbReac]]      = 'photo'
}

species=levels(as.factor(unlist(reactants)))
nbSpecies=length(species)
print("Species List:")
print(species)

L = matrix(0,ncol=nbSpecies,nrow=nbReac)
for (m in 1:nbReac) {
  reac = unlist(reactants[m])
  for (n in 1:nbSpecies) {
    search=species[n]
    L[m,n] = length(which( search == reac )) # Loss
  }
}

for (i in seq_along(species)) {
  sp=species[i]; print(sp)
  png(file=paste0('figXsec_',sp,'.png'),width=800,height=800)
  x = read.table(paste0('se',sp,'.dat'),
                 col.names=1:10,
                 header=FALSE, fill=TRUE)
  wavl = x[,1]
  xsec = matrix(1e-25,ncol=1,nrow=max(wavl))
  xsec[wavl] = x[,2]
  plot(xsec,log='y',type='l',col='pink',lwd=3,
       xlim=c(50,250),ylim=c(1e-20,1e-15), main=sp,
       xlab="Wavelength / nm", ylab="Cross-section / cm^2")

  reacs = which(L[,i]!=0) 
  leg = c()
  for (j in seq_along(reacs)) {
    jreac=reacs[j]
    channel = params[[jreac]][1]
    if(channel != 0) {
      file = paste0('qy',sp,'_',channel,'.dat')
      x = read.table(file)
      w  = x[,1]
      qy=matrix(1e-25,ncol=1,nrow=max(c(w,wavl)) )  
      qy[w] = x[,2]*xsec[w]
      lines(qy,col=j,lwd=2,lty=2)        
      leg[j] = paste(file,':',paste(unlist(products[jreac]),collapse=' + '))
    }
  }
  grid(col='gray30')
  if(length(leg)!=0) legend('topright',leg,col=1:length(leg),lwd=2,lty=2)
  dev.off()
}

  