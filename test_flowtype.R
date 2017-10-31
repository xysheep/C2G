library(flowType)
data(DLBCLExample)
PropMarkers <- 3:5
MFIMarkers <- PropMarkers
MarkerNames <- c('FS', 'SS','CD3','CD5','CD19')
Res <- flowType(DLBCLExample, PropMarkers, MFIMarkers, 'kmeans', MarkerNames)







