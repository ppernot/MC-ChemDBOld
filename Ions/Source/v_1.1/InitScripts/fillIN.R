version   = '1.0'
sourceDir = paste0('~/Bureau/ChemDB/ChemDB_Source/v_',version,'/')
dataDir   = paste0(sourceDir,'INData/')

library(stringr)
library(CHNOSZ)
data('thermo')

setwd(sourceDir)
X = as.matrix(read.csv(file='Doc/dbReacs.csv',
                       header=FALSE, sep=';', fill=TRUE,
                       stringsAsFactors = FALSE,
                       na.strings="")
)
X = str_trim(X) # remove unwanted spaces
codes = apply(X[,1:2], 1,function(x) str_split_fixed(x[1],'\\.',4)[2])
X[,1] = codes
X=as.data.frame(X,stringsAsFactors = FALSE)
colnames(X)=c("id","reac","dummy","rate","br1","br2")

# Convert data in matrix shape
for (f in unique(codes)) {
  print(f)
  X0 = X[X$id==f,]
  nBR = nrow(X0)
  nLines = nBR + 7
  nCols  = 10
  Y = matrix("",ncol=nCols,nrow=nLines)
  reac=X0$reac[1]
  reactants=str_trim(str_split_fixed(str_split_fixed(reac,' -> ',2)[1],' ',2))
  Y[1,1]=paste(reactants,collapse=' + ')
  
  rate=X0$rate[1]
  vmean=str_split_fixed(str_split_fixed(rate,'\\(',2)[2],',',2)[1]
  vF=sub('\\)','',str_split_fixed(str_split_fixed(rate,'\\(',2)[2],',',2)[2])
  Y[3,1]='ALPHA'  
  Y[3,2]='Logn'  
  Y[3,3]=vmean  
  Y[3,4]=vF  
  
  Y[5,1]='BETA'
  Y[5,2]='Delta'
  Y[5,3]='0'
  
  Y[7,1]='BR'
  for (i in 1:nBR) {
    reac=X0$reac[i]
    products=str_trim(unlist(str_split(str_split_fixed(reac,' -> ',2)[2],' ')))
    Y[7+i,1]=paste0(products,collapse=' + ')
    if (nBR==1) break
    Y[7,2]='Dirg'
    br = X0$br1[i]
    b1 = str_split_fixed(str_split_fixed(br,'\\(',2)[2],',',2)[1]
    b2 = sub('\\)','',str_split_fixed(str_split_fixed(br,'\\(',2)[2],',',2)[2])
    Y[7+i,2]=paste(b1,b2,sep='/')
  }    
  
  targetDir = paste0(dataDir,Y[1,1])
  dir.create(targetDir)
  targetFile = paste0(targetDir,'/data.csv')
  write.table(Y,file=targetFile, sep='\t', na="",quote=FALSE,
            row.names = FALSE,
            col.names = FALSE)
}
