setwd("~/Bureau/Titan-APSIS/Reactor_Runs/PhotoProcsTmp/CO")

S = read.table(file = 'extract_swri.dat', header=TRUE)

# Convert A to nm
lambda = S[,1]/10

matplot(lambda, log10(S[,2:ncol(S)]), type = 'l',
        xlim= c(50,150), ylim = c(-20,-15))

# Interpolate on 1 nm resolution, [5 - round(max(lambda))] grid
lg = seq(5,round(max(lambda)),by=1)

# Total XS
#
# xs = pmax(0,spline(lambda,S[,2],xout = lg)$y) # Splines. BOF
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
      if(i==1)
        xs[i] = 0
      else
        xs[i] = xs[i-1]
    else
      xs[i] = mean(xin[sel])
  }
  return(xs)
}
xs = average_1nm(lambda, S[,2], lg)

plot(lambda, S[,2], type = 'l')
points(lg,xs,pch=16,col=2)

write.table(cbind(lg,signif(xs,5)), sep=' ',
            row.names=FALSE, col.names=FALSE,
            file='seCO.dat')


# Quantum yields
np = ncol(S) - 2
qy = matrix(NA, nrow = length(lg), ncol = np)
for (i in 1:np) {
  xs = S[, i+2] / S[, 2]
  qy[,i] = average_1nm(lambda, xs, lg)

}

matplot(lambda, S[,3:(np+2)]/S[,2], type = 'l')
matpoints(lg, qy, pch = 16)
lines(lg,rowSums(qy))

for (i in 1:np) {
  write.table(cbind(lg,signif(qy[,i],5)), sep=' ',
              row.names=FALSE, col.names=FALSE,
              file=paste0('qyCO_',i,'.dat'))
}


