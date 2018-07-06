# Read datafiles in Source and generate Monte Carlo
# realizations in Publis/Databases, using
# transformation of the lognormal rate uncertainty params
# to account transparently for temperature effects:
# f=f^a; g=g*a, with a~norm(0,1)

version   = '1.0'
rootDir   = '/home/pernot/Bureau/Titan-APSIS/MC-ChemDB/'
sourceDir = paste0(rootDir,'Neutrals/Source/v_',version,'/')
publicDir = paste0(rootDir,'Neutrals/Public/v_',version,'/')
samplesDir= paste0(publicDir,'Databases/')
docDir    = paste0(sourceDir,'Doc/')

# Parameters #####
maxReacts      = 3    # Max number of reactants slots in generated dBases
maxProds       = 4    # Max number of product slots in generated dBases
sampleSize     = 500  # Number of random samples to generate

# Function to get stoichiometry and mass #####
library('CHNOSZ')
data('thermo')

elements = c('H','C','N','O','S','Ar')
massElem = mass(elements)

filterFormula <- function (sp) {
  sp1=sub('^l-','',sub('^c-','',sub('^t-','',sp )))
  sp1=sub('^l' ,'',sub('^c' ,'',sub('^t' ,'',sub('^i' ,'',sp1))))
  sp1=sub('\\(.*\\)','',sp1)
  sp1=sub('^X','',sp1)
  sp1=sub('^1','',sp1)
  sp1=sub('^3','',sp1)
  sp1=sub('C3P','C',sp1)
  sp1=sub('N4S','N',sp1)
  sp1=sub('N2D','N',sp1)
  sp1=sub('N1D','N',sp1)
  sp1=sub('N2P','N',sp1)
  sp1=sub('N3P','N',sp1)
  sp1=sub('O1D','O',sp1)
  sp1=sub('O3P','O',sp1)
  sp1=sub('S1D','S',sp1)
  sp1=sub('S3P','S',sp1)
  sp1=sub('\\+','',sp1)
  return(sp1)
}

calcAtoms = function(formula) {
  atoms=i2A(formula)  
  compo=matrix(0,nrow=1,ncol=length(elements))
  colnames(compo)=elements
  compo[1,colnames(atoms)]=atoms
  return(compo[1,])
}
get.atoms <- function (sp) {
  sp1=filterFormula(sp)
  tryCatch(calcAtoms(sp1),error= function(x) rep(NA,length(elements)))
}
massFormula = function(sto) sum(sto*massElem)
getMassList = function (species,
                        excludeList=c('Products','CxHyNz+','SOOTC+','SOOTN+',
                                      'CxHy','CxHy+','dummy')) {
  if( any(species %in% excludeList) ) {
    mass = NA
  } else {
    mass=sum(sapply(species,
                    function(x) massFormula(get.atoms(x))),
             na.rm=TRUE)
  }
  return(mass)
}
 
# Kinetic Parser ############################################################
nbReac=0  
reactants = products = params = type = orig = list()

filename='Titan - Réactions bimoléculaires.csv'
scheme  = read.csv(file=paste0(sourceDir,filename),header=FALSE)
# filename='reactions_1.dat'
# scheme  = read.fwf(file=filename, widths= rep(11,12))
scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
for (i in 1:nrow(scheme)) {
  nbReac = nbReac + 1
  terms=scheme[i,1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms!="" & terms!="HV"]
  terms=scheme[i,4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
  terms=scheme[i,9:13]
  params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
  params[[nbReac]][6] = 'kooij'
  type[[nbReac]]      = 'kooij'
  orig[[nbReac]]      = filename
}

filename = 'Titan - Réactions trimoléculaires.csv'
scheme  = read.csv(file=paste0(sourceDir,filename),header=FALSE)
# filename='reactions_2.dat'
# scheme  = read.fwf(file=filename, widths= rep(11,16))
scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
for (i in 1:nrow(scheme)) {
  nbReac = nbReac + 1
  terms=scheme[i,1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms!="" & terms!="HV"]
  terms=scheme[i,4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
  terms=scheme[i,9:18]
  params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
  params[[nbReac]][11] = '3-body'
  type[[nbReac]]      = '3-body'
  orig[[nbReac]]      = filename
}

# filename = 'Titan - Réactions bimol_trimol_association.csv'
# scheme  = read.csv(file=paste0(sourceDir,filename),header=FALSE)
# scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
# for (i in 1:nrow(scheme)) {
#   nbReac = nbReac + 1
#   terms=scheme[i,1:3]
#   reactants[[nbReac]] = terms[!is.na(terms) & terms!="" & terms!="HV"]
#   terms=scheme[i,4:8]
#   products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
#   terms=scheme[i,9:23]
#   params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
#   params[[nbReac]][16] = 'assoc'
#   type[[nbReac]]      = 'assoc'
#   orig[[nbReac]]      = filename
# }
  
species=levels(as.factor(unlist(c(reactants,products))))
nbSpecies=length(species)
print("Species List:")
print(species)

# Stoechiometry ######
compo = t(apply(as.matrix(species,ncol=1),1,get.atoms))
colnames(compo)=elements
mass  = apply(compo,1,massFormula)
names(mass) = species
dummyMass   = round(max(mass,na.rm=TRUE)+2)
dummySpecies = c('CxHy','CxHy+','CxHyNz','CxHyNz+',
                 'Products','HV','E','SOOTC','SOOTC+',
                 'C4H2X','HC3NX','HC5NX','SOOTS','SOOTN')
mass[dummySpecies] = dummyMass

# Check mass compatibility #####
for (m in 1:nbReac) {
  reac = reactants[m]
  prod = products[m]
  massReacs = getMassList(reac,excludeList = dummySpecies)
  massFrags = getMassList(prod,excludeList = dummySpecies)
  if(!is.na(massReacs) & !is.na(massFrags)) {
    if(abs(massFrags-massReacs) > 0.01) {
      setwd("..")
      stop(paste0('Pb mass of fragments vs. reactants: ',m))
    }      
  }
}

# Generate random samples #####

for (i in 0:sampleSize) {
  if(i%%10 == 0) cat(i,' over',sampleSize,'\n')
  
  if (i==0) # Nominal database
    rnd = rep(0,2*nbReac)
  else
    rnd = rnorm(2*nbReac,0,1)
  
  dbOut=data.frame(R1='x',R2='x',R3='x',
                   P1='x',P2='x',P3='x',P4='x',
                   a =NA,b =NA,c =NA,f =NA,g =NA,
                   a1=NA,b1=NA,c1=NA,f1=NA,g1=NA,
#                    a2=NA,b2=NA,c2=NA,f2=NA,g2=NA,
                   ttype='x',stringsAsFactors=FALSE)  

  for (m in 1:nbReac) {
    reac = unlist(reactants[m])
    prod = unlist(products[m])
    typ  = type[m]
    pars = unlist(params[m])
    
    if (typ == 'kooij') 
      db1=data.frame(R1=reac[1],R2=reac[2],R3=reac[3],
                     P1=prod[1],P2=prod[2],P3=prod[3],P4=prod[4],
                     a=as.numeric(pars[1]),
                     b=as.numeric(pars[2]),
                     c=as.numeric(pars[3]),
                     f=as.numeric(pars[4])^rnd[m],
                     g=as.numeric(pars[5])*rnd[m],
                     a1=NA,b1=NA,c1=NA,f1=NA,g1=NA,
#                      a2=NA,b2=NA,c2=NA,f2=NA,g2=NA,
                     ttype=typ,stringsAsFactors=FALSE)    
    else      
      if (typ == '3-body') 
        db1=data.frame(R1=reac[1],R2=reac[2],R3=reac[3],
                       P1=prod[1],P2=prod[2],P3=prod[3],P4=prod[4],
                       a=as.numeric(pars[1]),
                       b=as.numeric(pars[2]),
                       c=as.numeric(pars[3]),
                       f=as.numeric(pars[4])^rnd[m],
                       g=as.numeric(pars[5])*rnd[m],
                       a1=as.numeric(pars[6]),
                       b1=as.numeric(pars[7]),
                       c1=as.numeric(pars[8]),
                       f1=as.numeric(pars[9] )^rnd[nbReac+m],
                       g1=as.numeric(pars[10])*rnd[nbReac+m],
#                        a2=NA,b2=NA,c2=NA,f2=NA,g2=NA,
                       ttype=typ,stringsAsFactors=FALSE)        
#     else
#       db1=data.frame(R1=reac[1],R2=reac[2],R3=reac[3],
#                      P1=prod[1],P2=prod[2],P3=prod[3],P4=prod[4],
#                      a=as.numeric(pars[1]),
#                      b=as.numeric(pars[2]),
#                      c=as.numeric(pars[3]),
#                      f=as.numeric(pars[4])^rnd[m],
#                      g=as.numeric(pars[5])*rnd[m],
#                      a1=as.numeric(pars[6]),
#                      b1=as.numeric(pars[7]),
#                      c1=as.numeric(pars[8]),
#                      f1=as.numeric(pars[9] )^rnd[nbReac+m],
#                      g1=as.numeric(pars[10])*rnd[nbReac+m],
#                      a2=as.numeric(pars[6]),
#                      b2=as.numeric(pars[7]),
#                      c2=as.numeric(pars[8]),
#                      f2=as.numeric(pars[9] )^rnd[2*nbReac+m],
#                      g2=as.numeric(pars[10])*rnd[2*nbReac+m],
#                      ttype=typ,stringsAsFactors=FALSE)        
    
    dbOut[m,]=db1
  }
#   dbOutm <- within(dbOut, {
#     b <- sprintf("%6.3f", b)
#     a <- sprintf("%6.3e", a)
#     f <- sprintf("%6.3e", f)
#     g <- sprintf("%6.3f", g)
#   })
  write.table(dbOut,
              file=paste0(samplesDir,'run_',sprintf('%04i',i),'.csv'),
              quote=TRUE, sep=';', na=" ",
              row.names=FALSE,col.names=FALSE)  
}
