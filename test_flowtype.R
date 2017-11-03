library(flowType)
library(flowViz)
library(flowCore)
library(flowMeans)
setwd("C:/Users/xyshe/OneDrive/MyProject/c2g")
file.name <- "testdata/ctr.fcs"
d <- read.FCS(file.name, transformation=FALSE)
summary(d)
PropMarkers <- c(3,4,6,8,9,11,12,13,22,24,25,27);
tl <- transformList( featureNames(d)[PropMarkers], arcsinhTransform(b=0.2))
td <- transform(d, tl)
MFIMarkers <- PropMarkers
MarkerNames <- as.character(featureNames(td)[PropMarkers])
Res <- flowType(td, PropMarkers, MFIMarkers, 'flowMeans', MarkerNames)
plot(Res, td)

threshold = unlist(Res@Thresholds)
marker = Res@MarkerNames

asinht(a=5)


