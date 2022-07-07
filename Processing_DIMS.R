library("MALDIquant")
library("MALDIquantForeign")

# Set directiry with profile mass-spectra with multiple scans
setwd("/your/directory")
file.prof <- list.file()

# Set path for collectin consensus spectra
pth2 <- "/your/directory/consensus/spectra"

# Make consensus spectra with parallel
myCluster <- makeCluster(4, type="FORK")
registerDoParallel(myCluster)
foreach(i=1:length(file.prof)) %dopar% {
  pf <- readMSData(file.prof[i], pdata = NULL, mode = "onDisk")
  spec <- MSpectra(spectra(pf))
  comS <- combineSpectra(spec, timeDomain = F, mzd = 0.001, ppm = 1, weighted = TRUE)
  writeMSData(object = as(comS, "MSnExp"), file = paste0(pth2, basename(file.prof[i])))
}
stopCluster(myCluster)

# Import your consensus spectra with parallel
mzml <- importMzMl("/your/directory/consensus/spectra", centroided=FALSE, mc.cores=4,
                     verbose=FALSE)

# Set th method for transform intensity
t_spec <- transformIntensity(mzml, method="sqrt")

# Smoothing intensities with Savitzky-Golay method
smooth_spec <- smoothIntensity(t_spec, method="SavitzkyGolay", halfWindowSize=4)

# Set the method for baseline correction and number of iterations
bsln_off_spec <- removeBaseline(smooth_spec, method="SNIP", iterations=80)

# Make the spectra alignment
my_peaks <- alignSpectra(bsln_off_spec, halfWindowSize=2, noiseMethod="MAD", SNR=1,
                         tolerance=0.001, warpingMethod="lowess")

# Detect peaks with noise filtering method "median absolute deviation"
my_peaks <- detectPeaks(my_peaks, method="MAD", halfWindowSize=4, SNR=1)

# Make peak bins (all aligned peaks will have the same m/z values)
my_peaks <- binPeaks(my_peaks, method="strict", tolerance=0.001)


setwd("/home/rsuser/")
dir.create("/your/results/")
setwd("/your/results/")

#Export your results with MALDIquantForeign in csv (exportCsv) or txt (exportTxt) formats