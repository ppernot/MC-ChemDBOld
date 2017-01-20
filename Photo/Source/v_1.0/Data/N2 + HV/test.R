library(RColorBrewer)
col2tr = function(col,alpha) 
  rgb(t(col2rgb(col)),alpha=alpha,maxColorValue=255)
cols = brewer.pal(8,'Dark2')
col_tr=c()
for (i in 1:length(cols))
  col_tr[i] = col2tr(cols[i],10)

setwd("~/Bureau/Titan-APSIS/ChemDB/Photo/Source/v_1.0/Data/N2 + HV")

# Gather BR from data files
files = list.files(pattern = 'qy_')
br = list()
wl = c()
for(file in files) {
  tmp = read.table(file, 
                   header=FALSE,
                   stringsAsFactors = FALSE)
  br[[file]] = tmp
  wl = c(wl,tmp[,1])
}
wl = sort(unique(wl))
nl = length(wl)

BR = matrix(NA,ncol=length(files),nrow=nl)
i=0
for(file in files) {
  i = i+1
  BR[,i] = br[[file]][,2]
}

# ################################################################
# #### PB !!!!!  Request BR with hierarchical data...
# pat = "Dirg(BR1,(1-BR1)*Dirg(BR2,(1-BR2);0.2,0.2);0.03,0.03)"
# 
# # Generate stats per wl
# for (i in 1:nl) {
#   # 1 - Plug BR in nested pattern
#   
#   # 2 - Generate random values (nds)
#   
#   # 4 - Clr transform
#   
#   # 3 - Statistical summaries per channel
#   
# }
# 
# # Generate random samples
# for (iSample in 1:nSample)
#   for (i in 1:nBR) {
#     # 1- Generate random sample for each channel
#     # GP or mvnorm ?
#     
#     for (i in 1:nl)
#       # 2- ClrInv
#   }

##################################################
# Basic approach, without nesting....

# Threshold, normalize and log-ratio
thresh = 0.001
BR[BR<thresh] = thresh
BR[(1-BR)<thresh] = 1-thresh
BR=BR/rowSums(BR)
BR1 = compositions::clr(BR)

# Generate uncertainties
uq = c(0.03,0.2,0.2) # Not realistic without nesting....
taul = 10 # Correlation length

CM = matrix(NA,ncol=nl,nrow=nl)
for (i in 1:nl)
  for (j in 1:nl)
    CM[i,j] = exp(-abs(wl[i]-wl[j])/taul)

matplot(wl,BR1,col=cols,pch=19)
BRN = BR1
for (iMC in 1:200) {
  for (i in 1:ncol(BR1)) {
    BRN[,i] = mvtnorm::rmvnorm(1,BR1[,i],CM*uq[i]^2)
  }
  BRP=compositions::clrInv(BRN)
  matplot(wl,BRN,add=TRUE,type='l',col=col_tr)
}

