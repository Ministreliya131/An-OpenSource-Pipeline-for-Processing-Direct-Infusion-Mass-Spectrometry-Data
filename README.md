<h1 align="center">An-OpenSource-Pipeline-for-Processing-Direct-Infusion-Mass-Spectrometry-Data</h1>

<h2 align="center">
<p align="right"><i>Every year is getting shorter, never seem to find the time,</i></p>
<p align="right"><i>Plans that either come to naught or half a page of scribbled lines.</i></p>
</h2>
<h3 align="right"><i>Pink Floyd</i></h3>

<img src="./h1.jpg" width="100%">

## Description

The entire pipeline was built on the publicly available R language packages most commonly used to process, analyze and visualize mass spectrometry data (MSnbase, MALDIquant and MetaboAnalyst). 
Bioinformaticians can now reanalyze the direct infusion mass spectra (DIMS) without access to commercial software.

## Input format

Before processing your data should be converted into mzML or mzXML with MSconvert (a part of ProteoWizard, can be downloaded here: https://proteowizard.sourceforge.io/download.html )

<img src="./p1.png" width="80%">

## Make the consensus mass spectum

If your .mzML files contain > 1 spectra, then a consensus spectrum should be obtained. Use the MSnbase package in your RStudio IDE.
<p>For example:</p>

```R
library("MSnbase")
library("foreach")
library("parallel")
library("doParallel")
library("iterators")

file.prof <- "/path/to/your/files/"

# make consensus spectra with parallel
myCluster <- makeCluster(4, type="FORK")
registerDoParallel(myCluster)
foreach(i=1:length(file.prof)) %dopar% {
  pf <- readMSData(file.prof[i], pdata = NULL, mode = "onDisk")
  spec <- MSpectra(spectra(pf))
  cons.spec <- combineSpectra(spec, timeDomain = F, mzd = 0.001, ppm = 1, weighted = TRUE)
  writeMSData(object = as(cons.spec, "MSnExp"), file = paste0("/your/path/to/consensus/spectra/", basename(file.prof[i])))
}
stopCluster(myCluster)
```

<p>MSnbase documentation: https://www.bioconductor.org/packages/devel/bioc/manuals/MSnbase/man/MSnbase.pdf</p>

## Processing consensus spectra with MALDIquant

MALDIquant was originally developed for processing MALDI mass spectra. There is some superficial similarity between the MALDI spectrum and DIMS as shown in the figure below. Thus, the MALDIquant can be used to process direct input mass spectra with a standard pipeline for MALDI data.


