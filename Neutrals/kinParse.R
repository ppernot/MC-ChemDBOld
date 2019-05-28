# Kinetic Parser ############################################################

nbReac = 0
reactants = products = params = type = orig =  list()
filename = 'Titan - Réactions bimoléculaires.csv'
scheme  = read.csv(
  file = paste0(sourceDir, filename),
  header = FALSE,
  stringsAsFactors = FALSE
)
comments = scheme[, ncol(scheme)]
scheme  = t(apply(scheme, 1, function(x)
  gsub(" ", "", x)))
for (i in 1:nrow(scheme)) {
  nbReac = nbReac + 1
  terms = scheme[i, 1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 9:13]
  params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
  params[[nbReac]][6] = 'kooij'
  type[[nbReac]]      = 'kooij'
  orig[[nbReac]]      = filename
}

filename = 'Titan - Réactions bimoléculaires_supp.csv'
scheme  = read.csv(
  file = paste0(sourceDir, filename),
  header = FALSE,
  stringsAsFactors = FALSE
)
comments = scheme[, ncol(scheme)]
scheme  = t(apply(scheme, 1, function(x)
  gsub(" ", "", x)))
for (i in 1:nrow(scheme)) {
  nbReac = nbReac + 1
  terms = scheme[i, 1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 9:13]
  params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
  params[[nbReac]][6] = 'kooij'
  type[[nbReac]]      = 'kooij'
  orig[[nbReac]]      = filename
}

filename = 'Titan - Réactions trimoléculaires.csv'
scheme  = read.csv(
  file = paste0(sourceDir, filename),
  header = FALSE,
  stringsAsFactors = FALSE
)
comments = c(comments, scheme[, ncol(scheme)])
scheme  = t(apply(scheme, 1, function(x)
  gsub(" ", "", x)))
for (i in 1:nrow(scheme)) {
  nbReac = nbReac + 1
  terms = scheme[i, 1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 9:24]
  params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
  params[[nbReac]][17]= '3-body'
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
