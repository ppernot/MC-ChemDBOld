# Read datafiles in Source and generate Monte Carlo
# realizations in Public/Databases, using
# transformation of the lognormal rate uncertainty params
# to account transparently for temperature effects:
# f=f^a; g=g*a, with a~norm(0,1)

options(stringsAsFactors = FALSE)

version   = '1.2'
rootDir   = '/home/pernot/Bureau/Titan-APSIS/MC-ChemDB/'
sourceDir = paste0(rootDir,'Neutrals/Source/v_',version,'/')
publicDir = paste0(rootDir,'Neutrals/Public/v_',version,'/')
samplesDir= paste0(publicDir,'Databases/')
docDir    = paste0(sourceDir,'Doc/')

# Parameters #####
maxReacts      = 3    # Max number of reactants slots in generated dBases
maxProds       = 4    # Max number of product slots in generated dBases
sampleSize     = 500  # Number of random samples to generate

# Stoechiometry functions
source('./massCalc.R')

# Kinetic Scheme
source('./kinParse.R')
  
# species=levels(as.factor(unlist(c(reactants,products))))
species = sort(unique(unlist(c(reactants, products))))
nbSpecies = length(species)
print("Species List:")
print(species)

# Stoechiometry ######
compo = t(apply(as.matrix(species,ncol=1),1,get.atoms))
colnames(compo)=elements
mass  = apply(compo,1,massFormula)
names(mass) = species
dummyMass   = round(max(mass,na.rm=TRUE)+2)
mass[dummySpecies] = dummyMass

# Check mass compatibility #####
massProds = reacTags = c()
for (m in 1:nbReac) {
  reac = unlist(reactants[m])
  prod = unlist(products[m])
  massReacs = getMassList(reac,excludeList = dummySpecies)
  massFrags = getMassList(prod,excludeList = dummySpecies)
  reacTags[m] = paste(
    paste(reac,collapse = ' + '),
    '-->',
    paste(prod,collapse = ' + ')
  )
  if(!is.na(massReacs) & !is.na(massFrags)) {
    if(abs(massFrags-massReacs) > 0.01) {
      # setwd("..")
      stop(
        paste0(
          'Pb mass of fragments vs. reactants: ',m,' : \n',
          reacTags[m],'\n',
          massReacs, ' /= ', massFrags
        )
      )
    }      
  }
  massProds[m] = round(massReacs,4)
}
names(massProds) = reacTags

# Reactions per mass to identify unsuitable dummy species
sink(file = 'reacsMassList.txt')
mu = sort(unique(massProds))
for (m in mu) {
  cat(
    m,' : ',
    paste(names(massProds)[massProds==m],
          collapse = '\n\t\t\t\t\t\t '),
    '\n\n'
  )
}
sink()


# Species per mass
sink(file = 'speciesList.txt')
mu = sort(unique(mass))
for (m in mu) {
  cat(m,' : ',names(mass)[mass==m],'\n')
}
sink()

stop()

# Generate random samples #####

for (i in 0:sampleSize) {
  if(i%%10 == 0) cat(i,' over',sampleSize,'\n')
  
  if (i==0) # Nominal database
    rnd = rep(0,3*nbReac)
  else
    rnd = rnorm(3*nbReac,0,1)
  
  dbOut = data.frame(
    R1 = 'x', R2 = 'x', R3 = 'x',
    P1 = 'x', P2 = 'x', P3 = 'x', P4 = 'x',
    a = NA,
    b = NA,
    c = NA,
    f = NA,
    g = NA,
    a1 = NA,
    b1 = NA,
    c1 = NA,
    f1 = NA,
    g1 = NA,
    a2 = NA,
    b2 = NA,
    c2 = NA,
    f2 = NA,
    g2 = NA,
    fc = NA,
    ttype = 'x'
  )
  
  for (m in 1:nbReac) {
    reac = unlist(reactants[m])
    prod = unlist(products[m])
    typ  = type[m]
    pars = unlist(params[m])
    
    if (typ == 'kooij') {
      db1 = data.frame(
        R1 = reac[1], R2 = reac[2], R3 = reac[3],
        P1 = prod[1], P2 = prod[2], P3 = prod[3], P4 = prod[4],
        a = as.numeric(pars[1]),
        b = as.numeric(pars[2]),
        c = as.numeric(pars[3]),
        f = as.numeric(pars[4]) ^ rnd[m],
        g = as.numeric(pars[5]) * rnd[m],
        a1 = NA,
        b1 = NA,
        c1 = NA,
        f1 = NA,
        g1 = NA,
        a2 = NA,
        b2 = NA,
        c2 = NA,
        f2 = NA,
        g2 = NA,
        fc = NA,
        ttype = typ
      )
    } else {
      if (typ == 'assocMD' | 
          typ == 'assocVV'  ) {
        db1 = data.frame(
          R1 = reac[1], R2 = reac[2], R3 = reac[3],
          P1 = prod[1], P2 = prod[2], P3 = prod[3], P4 = prod[4],
          a = as.numeric(pars[1]),
          b = as.numeric(pars[2]),
          c = as.numeric(pars[3]),
          f = as.numeric(pars[4]) ^ rnd[m],
          g = as.numeric(pars[5]) * rnd[m],
          a1 = as.numeric(pars[6]),
          b1 = as.numeric(pars[7]),
          c1 = as.numeric(pars[8]),
          f1 = as.numeric(pars[9]) ^ rnd[nbReac + m],
          g1 = as.numeric(pars[10]) * rnd[nbReac + m],
          a2 = as.numeric(pars[11]),
          b2 = as.numeric(pars[12]),
          c2 = as.numeric(pars[13]),
          f2 = as.numeric(pars[14]) ^ rnd[2 * nbReac + m],
          g2 = as.numeric(pars[15]) * rnd[2 * nbReac + m],
          fc = as.numeric(pars[16]),
          ttype = typ
        )
      }
    }      
    dbOut[m,]=db1
  }
  write.table(dbOut,
              file=paste0(samplesDir,'run_',sprintf('%04i',i),'.csv'),
              quote=TRUE, sep=';', na=" ",
              row.names=FALSE,col.names=FALSE)  
}
