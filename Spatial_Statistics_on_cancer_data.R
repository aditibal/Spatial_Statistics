library(spdep)
library(maps)
library(maptools)
library(classInt)
library(RColorBrewer)
library(rgdal)
library(spatialreg)
library(stats)
library(ggplot2)
library(ggpubr)



usa.state=map(database="state",fill=TRUE,plot=FALSE)
state.ID<-sapply(strsplit(usa.state$names,":"),function(x) x[1])
usa.poly=map2SpatialPolygons(usa.state,IDs=state.ID)
usa.nb=poly2nb(usa.poly)
usa.listw_binary=nb2listw(usa.nb,style="B",zero.policy = TRUE)
usa.listw_rownorm=nb2listw(usa.nb,style="W")


filename=paste0("project_data_childhood_cancer.csv")
z<-read.csv(filename, header = FALSE, sep = ",")
z<-as.matrix(z)
z<-z[2:50,]
count<-as.numeric(z[,5])
age<-as.numeric(z[,2])
sex<-as.numeric(z[,3])
race<-as.numeric(z[,4])

##choropleth map
brks=c(0,25,50,75,100,125,150,175,200)
color.pallete=rev(brewer.pal(8,"RdBu"))
class.fitted = classIntervals(var=count,n=8,style="fixed", fixedBreaks=brks,dataPrecision=8)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,25,50,75,100,125,150,175,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Counts")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

# Moran's I Geary's C 
moran.test(count,listw=usa.listw_rownorm)
geary.test(count,listw=usa.listw_rownorm)


############ SAR
outBin = spautolm(count~age+sex+race,family="SAR", listw =usa.listw_binary, zero.policy = TRUE)
summary(outBin)

outRowN = spautolm(count~age+sex+race,family="SAR", listw =usa.listw_rownorm)
summary(outRowN)

#plot SAR predictions for outBin
fittedBin = fitted(outBin)
brks=c(-25,0,25,50,75,100)
color.pallete=rev(brewer.pal(5,"RdBu"))
class.fitted = classIntervals(var=fittedBin,n=5,style="fixed", fixedBreaks=brks,dataPrecision=5)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,0,25,50,75,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Fitted Values")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

resBin<-residuals.spautolm(outBin)
brks=c(-56,-10,0,10,150)
color.pallete=rev(brewer.pal(4,"RdBu"))
class.fitted = classIntervals(var=resBin,n=4,style="fixed", fixedBreaks=brks,dataPrecision=4)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,-10,0,10,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Residuals")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

#plot SAR predictions for outRowN
fittedRow = fitted(outRowN)
brks=c(-25,0,25,50,75,102)
color.pallete=rev(brewer.pal(5,"RdBu"))
class.fitted = classIntervals(var=fittedRow,n=5,style="fixed", fixedBreaks=brks,dataPrecision=5)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,0,25,50,75,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Fitted Values")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

resRow<-residuals.spautolm(outRowN)
brks=c(-55,-10,0,10,140)
color.pallete=rev(brewer.pal(4,"RdBu"))
class.fitted = classIntervals(var=resRow,n=4,style="fixed", fixedBreaks=brks,dataPrecision=4)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,-10,0,10,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Residuals")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

############ CAR
outBin = spautolm(count~age+sex+race,family="CAR", listw =usa.listw_binary, zero.policy = TRUE)
summary(outBin)

outRowN = spautolm(count~age+sex+race,family="CAR", listw =usa.listw_rownorm)
summary(outRowN)

#plot CAR predictions for outBin
fittedBin = fitted(outBin)
brks=c(-25,0,25,50,75,100)
color.pallete=rev(brewer.pal(5,"RdBu"))
class.fitted = classIntervals(var=fittedBin,n=5,style="fixed", fixedBreaks=brks,dataPrecision=5)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,0,25,50,75,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Fitted Values")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

resBin<-residuals.spautolm(outBin)
brks=c(-56,-10,0,10,150)
color.pallete=rev(brewer.pal(4,"RdBu"))
class.fitted = classIntervals(var=resBin,n=4,style="fixed", fixedBreaks=brks,dataPrecision=4)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,-10,0,10,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Residuals")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

#plot CAR predictions for outRowN
fittedRow = fitted(outRowN)
brks=c(-26,0,25,50,75,113)
color.pallete=rev(brewer.pal(5,"RdBu"))
class.fitted = classIntervals(var=fittedRow,n=5,style="fixed", fixedBreaks=brks,dataPrecision=5)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,0,25,50,75,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Fitted Values")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))

resRow<-residuals.spautolm(outRowN)
brks=c(-55,-10,0,10,130)
color.pallete=rev(brewer.pal(4,"RdBu"))
class.fitted = classIntervals(var=resRow,n=4,style="fixed", fixedBreaks=brks,dataPrecision=4)
color.code.fitted = findColours(class.fitted, color.pallete)
brks = c(-Inf,-10,0,10,Inf)
plot(usa.poly,col = color.code.fitted)
title(main="Residuals")
legend("bottomleft", fill=color.pallete, legend=leglabs(brks))
