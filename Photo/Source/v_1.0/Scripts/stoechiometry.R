library(CHNOSZ)

calcAtoms = function(formula) {
  atoms=i2A(formula)
  nC = tryCatch(atoms[1,'C'],error= function(x) 0)
  nH = tryCatch(atoms[1,'H'],error= function(x) 0)
  nN = tryCatch(atoms[1,'N'],error= function(x) 0)
  nO = tryCatch(atoms[1,'O'],error= function(x) 0)
  c(nC,nH,nN,nO)
}
get.atoms = function (sp) {
  sp1=sub('^l','',sub('^c','',sub('^t','',sp)))
  sp1=sub('^X','',sp1)
  sp1=sub('^1','',sp1)
  sp1=sub('^3','',sp1)
  sp1=sub('HV','',sp1)
  sp1=sub('E','',sp1)
  sp1=sub('H2S','H',sp1)
  sp1=sub('N4S','N',sp1)
  sp1=sub('N2D','N',sp1)
  sp1=sub('N2P','N',sp1)
  sp1=sub('N3P','N',sp1)
  sp1=sub('N1D','N',sp1)
  sp1=sub('C1D','C',sp1)
  sp1=sub('C3P','C',sp1)
  sp1=sub('C1S','C',sp1)
  sp1=sub('O3P','O',sp1)
  sp1=sub('O1D','O',sp1)
  sp1=sub('\\(.*\\)','',sp1)
  tryCatch(calcAtoms(sp1),error= function(x) NA)
}
massCxHyNzOw = function(sto) {
  sto[1]*12.0107 + sto[2]*1.00794 + 
    sto[3]*14.00674 + sto[4]*15.994915
}
getMassList = function (species,
                        excludeList=c('Products','CxHyNz+','SOOTC+','SOOTN+',
                                      'CxHy','CxHy+','dummy')) {
  if( any(species %in% excludeList) ) {
    mass = NA
  } else {
    mass=sum(sapply(species,
                    function(x) massCxHyNzOw(get.atoms(x))),
             na.rm=TRUE)
  }
  return(mass)
}
