library(devtools)
# install_github('ramnathv/htmlwidgets')
# devtools::install_github("bokeh/rbokeh")

library(htmlwidgets)
library(rbokeh)
library(RColorBrewer)

files = c('HCN_Titan.dat','HCN_meudon.csv',
          'bessy_new_195K.txt','bessy_new_295K.txt',
          'huebner_1992.txt')

p = figure(width = 800, height = 600)
cols=brewer.pal(8,'Dark2')
i=0
for(file in files) {
  X = read.table(file,header=FALSE)
  names(X)=c('Wavelength','CrossSection')
  i=i+1
  p = ly_lines(p,Wavelength,CrossSection,data=X,
               color=cols[i],legend=file)
  p = ly_points(p,Wavelength,CrossSection,data=X,
                color=cols[i],size=5)   
}
p