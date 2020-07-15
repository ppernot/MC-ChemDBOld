library('CHNOSZ')
# data('thermo') # Obsolete

calcAtoms = function(formula) {
  atoms=i2A(formula)
  nC = tryCatch(atoms[1,'C' ],error= function(x) 0)
  nH = tryCatch(atoms[1,'H' ],error= function(x) 0)
  nN = tryCatch(atoms[1,'N' ],error= function(x) 0)
  nO = tryCatch(atoms[1,'O' ],error= function(x) 0)
  nAr= tryCatch(atoms[1,'Ar'],error= function(x) 0)
  list(C=nC,H=nH,N=nN,O=nO,Ar=nAr)
}
get.atoms = function (sp) {
  sp1=sub('^l','',sub('^c','',sub('^t','',sp)))
  sp1=sub('^X','',sp1)
  sp1=sub('X$','',sp1)
  sp1=sub('^1','',sp1)
  sp1=sub('^2','',sp1)
  sp1=sub('^3','',sp1)
  sp1=sub('^4','',sp1)
  sp1=sub('HV','',sp1)
  sp1=sub('E','',sp1)
  sp1=sub('N4S','N',sp1)
  sp1=sub('N2D','N',sp1)
  sp1=sub('N2P','N',sp1)
  sp1=sub('N3P','N',sp1)
  sp1=sub('N1D','N',sp1)
  sp1=sub('C1D','C',sp1)
  sp1=sub('C3P','C',sp1)
  sp1=sub('C1S','C',sp1)
  sp1=sub('O3P','O',sp1)
  sp1=sub('\\(.*\\)','',sp1)
  tryCatch(calcAtoms(sp1),
           error= function(x) list(C=NA,H=NA,N=NA,O=NA,Ar=NA))
}
massCalc = function(sto) {
  sto$C  * 12.0107 + 
    sto$H  *  1.00794 + 
    sto$N  * 14.00674 + 
    sto$O  * 15.994915 + 
    sto$Ar * 39.962383
}
# getMass = function(x) massCalc(get.atoms(x))
getMassList = function (species,
                        excludeList=c('Products','CxHyNz+','SOOTC+','SOOTN+',
                                      'CxHy','CxHy+','dummy')) {
  if( any(species %in% excludeList) ) {
    mass = NA
  } else {
    mass=sum(sapply(species,
                    function(x) massCalc(get.atoms(x))),
             na.rm=TRUE)
  }
  return(mass)
}
