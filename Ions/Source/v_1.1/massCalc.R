library('CHNOSZ')

elements = c('H', 'C', 'N', 'O', 'S', 'Ar')
massElem = mass(elements)

dummySpecies = c(
  # 'HV', # Treated in calcAtoms
  # 'E',  # Treated in calcAtoms
  'CxHy',
  'CxHy+',
  'CxHyNz',
  'CxHyNz+',
  'Products',
  'SOOTC',
  'SOOTC+',
  'C4H2X',
  'HC3NX',
  'HC5NX',
  'SOOTS',
  'SOOTN'
)
filterFormula <- function (sp) {
  sp1 = sub('^l-', '', sub('^c-', '', sub('^t-', '', sp)))
  sp1 = sub('^l' , '', sub('^c' , '', sub('^t' , '', sub('^i' , '', sp1))))
  sp1 = sub('\\(.*\\)', '', sp1)
  sp1 = sub('^X', '', sp1)
  sp1 = sub('X$', '', sp1)
  sp1 = sub('^1', '', sp1)
  sp1 = sub('^2', '', sp1)
  sp1 = sub('^3', '', sp1)
  sp1 = sub('^4', '', sp1)
  sp1 = sub('C1D', 'C', sp1)
  sp1 = sub('C1S', 'C', sp1)
  sp1 = sub('C3P', 'C', sp1)
  sp1 = sub('N4S', 'N', sp1)
  sp1 = sub('N2D', 'N', sp1)
  sp1 = sub('N1D', 'N', sp1)
  sp1 = sub('N2P', 'N', sp1)
  sp1 = sub('N3P', 'N', sp1)
  sp1 = sub('O1D', 'O', sp1)
  sp1 = sub('O3P', 'O', sp1)
  sp1 = sub('S1D', 'S', sp1)
  sp1 = sub('S3P', 'S', sp1)
  sp1 = sub('\\+', '', sp1)
  return(sp1)
}
calcAtoms = function(formula) {
  compo = matrix(0, nrow = 1, ncol = length(elements))
  colnames(compo) = elements
  if (!(formula %in% c("E","HV"))) { # Mass = 0
    atoms = CHNOSZ::i2A(formula)
    compo[1, colnames(atoms)] = atoms
  }
  return(compo[1, ])
}
get.atoms <- function (sp) {
  sp1 = filterFormula(sp)
  tryCatch(
    calcAtoms(sp1),
    error = function(x)
      rep(NA, length(elements))
  )
}
massFormula = function(sto) {
  sum(sto * massElem)
}
getMassList = function (species,excludeList = 'Products') {
  if (any(species %in% excludeList)) {
    mass = NA
  } else {
    mass = sum(sapply(species,
                      function(x)
                        massFormula(get.atoms(x))),
               na.rm = TRUE)
  }
  return(mass)
}
