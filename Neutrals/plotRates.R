rm(list = ls())

version   = '1.1'
rootDir   = '/home/pernot/Bureau/Titan-APSIS/MC-ChemDB/'
sourceDir = paste0(rootDir, 'Neutrals/Source/v_', version, '/')
publicDir = paste0(rootDir, 'Neutrals/Public/v_', version, '/')
tmpDir    = paste0(rootDir, 'Neutrals/Tmp/')
samplesDir = paste0(publicDir, 'Databases/')
docDir    = paste0(sourceDir, 'Doc/')
figsDir   = paste0(publicDir, 'Figs/')


# Clean target dirs
command = paste0('rm -rf ', tmpDir, '*')
dummy = system(command)

# List sample files
samplesList = list.files(
  path = paste0(publicDir, 'Databases'),
  pattern = 'run_',
  full.name = TRUE
)[1:51]

kooij = function(p)
  p[1] * (tempRange / T0)^p[2] * exp(-p[3]/tempRange) *
  p[4] * exp(p[5] * abs(1/tempRange -1/ T0))

# T varies, fixed density (M)
T0 = 300
tempRange = seq(50, 350, by = 5)
M  = 1e18 # molec/cm^3

# Generate individual curves
irun = -1
for (file in samplesList) {
  irun = irun + 1
  
  nbReac = 0
  reactants = products = params = type = tag = list()
  scheme  = read.csv(file = file,
                     header = FALSE,
                     sep = ';')
  scheme  = t(apply(scheme, 1, function(x)
    gsub(" ", "", x)))
  for (i in 1:nrow(scheme)) {
    nbReac = nbReac + 1
    terms = scheme[i, 1:3]
    reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
    terms = scheme[i, 4:7]
    products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
    terms = scheme[i, 8:24]
    params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
    tag[[nbReac]]       = paste0(
      nbReac,
      ': ',
      paste0(reactants[[nbReac]], collapse =
               '+'),
      '->',
      paste0(products[[nbReac]], collapse = '+')
    )
  }
  
  alerts = c()
  for (i in 1:nbReac) {
    tdir = paste0(tmpDir, tag[[i]])
    if (!file.exists(tdir))
      dir.create(tdir)
    l = length(params[[i]])
    pars = as.numeric(params[[i]][1:(l-1)])
    type = params[[i]][l]
    
    if (type == 'kooij') {
      k = kooij(pars[1:5])
      if (sum(k <= 0) != 0)
        alerts = c(alerts, paste0('Null RC: ', tag[[i]], '\n'))
      
    } else {
      # 3-body
      k0   = kooij(pars[1:5])
      kInf = kooij(pars[6:10])
      kr   = kooij(pars[11:15])
      Fc   = pars[16] 
      
      # kAss from Dobrijevic2016
      Ni = 1
      F1 = exp( log(Fc) / (1 + (log(k0 * M + kInf) / Ni) ^ 2) )
      k  = kInf*(k0*M*F1 + kr) / (k0*M + kInf)
      
      # # Troe formula
      # Pr   = k0 * pressure / kInf
      # cExp = -0.4 - 0.67 * log10(fc)
      # NExp = 0.75 - 1.27 * log10(fc)
      # dExp = 0.14
      # fExp = 1 + ((log10(Pr) + cExp) / (NExp - dExp * (log10(Pr) + cExp)))^2
      # broadF = fc ^ (1 / fExp)
      # k = kInf * (Pr / (1 + Pr)) * broadF

      
      if (sum(k <= 0) != 0)
        alerts = c(alerts, paste0('Null RC: ', tag[[i]], '\n'))
      
    }
    write.table(
      data.frame(T = tempRange, kval = k),
      file = paste0(tdir, '/curve_T_', sprintf('%04i', irun), '.csv'),
      row.names = FALSE,
      col.names = FALSE
    )
  }
}
if (length(alerts != 0))
  print(unique(alerts))

# fixed T, M varies
T0 = 300
tempRange = 150
M  = 10^seq(8,20,by=1) # molec/cm^3

# Generate individual curves
irun = -1
for (file in samplesList) {
  irun = irun + 1
  
  nbReac = 0
  reactants = products = params = type = tag = list()
  scheme  = read.csv(file = file,
                     header = FALSE,
                     sep = ';')
  scheme  = t(apply(scheme, 1, function(x)
    gsub(" ", "", x)))
  for (i in 1:nrow(scheme)) {
    nbReac = nbReac + 1
    terms = scheme[i, 1:3]
    reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
    terms = scheme[i, 4:7]
    products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
    terms = scheme[i, 8:24]
    params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
    tag[[nbReac]]       = paste0(
      nbReac,
      ': ',
      paste0(reactants[[nbReac]], collapse =
               '+'),
      '->',
      paste0(products[[nbReac]], collapse = '+')
    )
  }
  
  alerts = c()
  for (i in 1:nbReac) {
    tdir = paste0(tmpDir, tag[[i]])
    if (!file.exists(tdir))
      dir.create(tdir)
    l = length(params[[i]])
    pars = as.numeric(params[[i]][1:(l-1)])
    type = params[[i]][l]
    
    if (type == 'kooij') {
      k = kooij(pars[1:5])
      if (sum(k <= 0) != 0)
        alerts = c(alerts, paste0('Null RC: ', tag[[i]], '\n'))
      
    } else {
      # 3-body
      k0   = kooij(pars[1:5])
      kInf = kooij(pars[6:10])
      kr   = kooij(pars[11:15])
      Fc   = pars[16] 
      
      # kAss from Dobrijevic2016
      Ni = 1
      F1 = exp( log(Fc) / (1 + (log(k0 * M + kInf) / Ni) ^ 2) )
      k  = kInf*(k0*M*F1 + kr) / (k0*M + kInf)
      
      # # Troe formula
      # Pr   = k0 * M / kInf
      # cExp = -0.4 - 0.67 * log10(fc)
      # NExp = 0.75 - 1.27 * log10(fc)
      # dExp = 0.14
      # fExp = 1 + ((log10(Pr) + cExp) / (NExp - dExp * (log10(Pr) + cExp)))^2
      # broadF = fc ^ (1 / fExp)
      # k = kInf * (Pr / (1 + Pr)) * broadF
      
      if (sum(k <= 0) != 0)
        alerts = c(alerts, paste0('Null RC: ', tag[[i]], '\n'))
      
    }
    write.table(
      data.frame(P = M, kval = k),
      file = paste0(tdir, '/curve_P_', sprintf('%04i', irun), '.csv'),
      row.names = FALSE,
      col.names = FALSE
    )
  }
}
if (length(alerts != 0))
  print(unique(alerts))


# Gather info from source database for legending plots
nbReac_old = nbReac
source('./kinParse.R')
if (nbReac_old != nbReac) {
  cat('Pb. of database version\n')
  stop()
}

# Generate plots from curves
fileList = list.files(path = figsDir, full.name = TRUE)
dummy = file.remove(fileList)

col2tr = function(x, alpha = 80) {
  rgb(unlist(t(col2rgb(x))), alpha = alpha, maxColorValue = 255)
}
trBlue = col2tr('blue', 60)
for (ireac in 1:nbReac) {
  reac = tag[[ireac]]
  
  legText = paste0(reac, '\n',
                   'Rate law: ', type[[ireac]])
  if (type[[ireac]] == 'kooij') {
    legText = paste0(legText,
                     '\n',
                     'Parameters: ', 
                     paste0(params[[ireac]][1:5],
                            collapse = ' / '),
                     '\n')
  } else {
    legText = paste0(
      legText,'; Fc : ', params[[ireac]][16],
      '\n',
      'Parameters k0  : ',
      paste0(params[[ireac]][1:5], collapse = ' / '),
      '\n',
      'Parameters kInf: ',
      paste0(params[[ireac]][6:10], collapse = ' / '),
      '\n',
      'Parameters kr: ',
      paste0(params[[ireac]][11:15], collapse = ' / '),
      '\n'
    )
  }
  
  com = comments[ireac]
  if (!is.na(com)) {
    splCom = unlist(strsplit(com, ' '))
    while (length(splCom) != 0) {
      sel = cumsum(nchar(splCom)) <= 65
      line = paste0(splCom[sel], collapse = ' ')
      legText = paste0(legText, '! ', line, '\n')
      splCom = splCom[!sel]
    }
  }
  
  
  png(
    file = paste0(figsDir, reac, '.png'),
    width = 2000,
    height = 1000
  )
  par(
    mfrow = c(1, 2),
    mar = c(4, 6, 20, 1),
    cex.lab = 2,
    cex.axis = 2,
    cex.main = 3
  )
  
  tdir = paste0(tmpDir, reac)
  
  # T-dep
  fileList = list.files(path = tdir,
                        pattern = 'curve_T',
                        full.names = TRUE)
  dat = read.csv(fileList[1], header = FALSE, sep = ' ')
  tempRange = dat[, 1]
  ktab = matrix(0, nrow = length(tempRange), ncol = length(fileList))
  for (i in seq_along(fileList))
    ktab[, i] = read.csv(fileList[i], header = FALSE, sep = ' ')[, 2]
  # if (diff(range(ktab)) != 0) {
  matplot(
    tempRange,
    ktab,
    type = 'l',
    lty = 1,
    col = trBlue,
    lwd = 3,
    log = 'y',
    xlab = 'T / K',
    ylab = 'rate ct. / cm^3.s^-1',
    ylim = c(1e-18, 1e-7),
    main = ''
  )
  lines(tempRange, ktab[, 1], col = 'red', lwd = 3)
  grid(col = 'darkgray')
  legend('topright',title = 'M = 1e18 cm^-3', legend = NA, bty='n',cex=2)
  mtext(
    legText,
    side = 3,
    cex = 2,
    adj = 0,
    line = 18,
    padj = 1,
    col = 'darkgreen'
  )
  box(lwd = 4)
  
  # P-dep
  fileList = list.files(path = tdir,
                        pattern = 'curve_P',
                        full.names = TRUE)
  dat = read.csv(fileList[1], header = FALSE, sep = ' ')
  tempRange = dat[, 1]
  ktab = matrix(0, nrow = length(tempRange), ncol = length(fileList))
  for (i in seq_along(fileList))
    ktab[, i] = read.csv(fileList[i], header = FALSE, sep = ' ')[, 2]
  # if (diff(range(ktab)) != 0) {
  matplot(
    tempRange,
    ktab,
    type = 'l',
    lty = 1,
    col = trBlue,
    lwd = 3,
    log = 'xy',
    xlab = 'M / cm^-3',
    ylab = 'rate ct. / cm^3.s^-1',
    ylim = c(1e-18, 1e-7),
    main = ''
  )
  lines(tempRange, ktab[, 1], col = 'red', lwd = 3)
  grid(col = 'darkgray')
  legend('topright',title = 'T = 150 K', legend = NA, bty='n',cex=2)
  box(lwd = 4)
  
  dev.off()
  
}
