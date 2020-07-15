species = 'H2O'

# setwd(paste0("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcsTmp/",species))

S = read.table(file = 'extract_swri.dat', header=TRUE)

# Convert A to nm
lambda = S[,1]/10

matplot(lambda, log10(S[,2:ncol(S)]), type = 'l',
        xlim= c(50,150), ylim = c(-20,-15))

# Interpolate on 1 nm resolution, [5 - lmax] grid
lmax = round(max(lambda))
# lmax = 250 # !!!!!!!!!!!!!!!!!!!!!!!!!!
lg = seq(5,lmax,by=1)

# Total XS
#
xs = pmax(min(S[,2]),spline(lambda,S[,2],xout = lg)$y) # Splines. BOF

#
# xs = c()
# for(i in 1:length(lg)) {
#   np = which.min(abs(lambda - lg[i])) # Select nearest point. BOF
#   xs[i] = S[np,2]
# }

## Average on 1 nm interval
average_1nm <- function(lambda, xin, lg) {
  xs = c()
  for(i in 1:length(lg)) {
    sel = which( (lambda >= lg[i]-0.5) & (lambda < lg[i]+0.5))
    if(sum(sel) == 0)
      if(i==1) {
        xs[i] = 0
      } else {
        # Lin. Average between nearest points
        i1 = which(lambda >= lg[i])[1]-1
        x1 = xin[i1]
        i2 = i1 + 1
        x2 = xin[i2]
        xs[i] = x1 + (x2-x1)/(lambda[i2]-lambda[i1])*(lg[i]-lambda[i1])
      }
    else
      xs[i] = mean(xin[sel])
  }
  return(xs)
}
xs = average_1nm(lambda, S[,2], lg)

plot(lambda, S[,2], type = 'l', pch=1, xlim=range(lg))
points(lg,xs,pch=16,cex=0.5,col=2)

write.table(
  cbind(lg,signif(xs,5)), sep=' ',
  row.names=FALSE, col.names=FALSE,
  file=paste0('se',species,'.dat')
)


# Quantum yields
np = ncol(S) - 2
qy = matrix(NA, nrow = length(lg), ncol = np)
for (i in 1:np) {
  xs = S[, i+2] / S[, 2]
  qy[,i] = average_1nm(lambda, xs, lg)
  # qy[,i] =pmin(1,pmax(0,spline(lambda,xs,xout = lg)$y))
}
qy = qy /rowSums(qy) #!!!!!!!!!!!!!!
qy[qy < 1e-25] = 0

matplot(lambda, S[,3:(np+2)]/S[,2], type = 'l', xlim=range(lg))
matpoints(lg, qy, pch = 16)
lines(lg,rowSums(qy))

for (i in 1:np) {
  write.table(
    cbind(lg,signif(qy[,i],5)), sep=' ',
    row.names=FALSE, col.names=FALSE,
    file=paste0('qy',species,'_',i,'.dat')
  )
}


