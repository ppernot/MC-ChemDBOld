lines = readLines(con="list_reac_ion.dat")
scheme = lapply(lines,
                 FUN=function(x) scan(text=x, what=character(), 
                                      strip.white=TRUE, quiet=TRUE))
nbReac=length(scheme)

reactants = products = list()
for (i in 1:nbReac) {
  terms=scheme[[i]][1:2]
  reactants[[i]] = terms[!is.na(terms) & terms!=""]
  terms=scheme[[i]][3:length(scheme[[i]])]
  products [[i]] = terms[!is.na(terms) & terms!=""]
}

tabR=unlist(lapply(reactants,function(x) paste(x,collapse=" + ")))
tabP=unlist(lapply(products,function(x) paste(x,collapse=" + ")))
tab= cbind(tabR,tabP)
write.table(tab,file='list_reac_ion_collapsed.csv',
            sep=';',quote=FALSE,col.names=FALSE,row.names=FALSE)
