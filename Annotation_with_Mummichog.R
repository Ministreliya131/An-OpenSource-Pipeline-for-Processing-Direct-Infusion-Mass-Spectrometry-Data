library("MetaboAnalystR")

# set working directory with mass-lists in csv format
pth <- "/your/directory/with/mass_lists/"
setwd(pth)
mumFles <- list.files()

# Initialize internal object for Mummichog 
mSet<-InitDataObjects("mass_all", "mummichog", FALSE)
# Set the internal structure of mass-list. Set "rmp" if retention time, m/z value and intensity are presented.
mSet<-SetPeakFormat(mSet, "rmp")
# Set the mode of obtainin mass-spectrum, mass tolerance and resolution
mSet<-UpdateInstrumentParameters(mSet, 10.0, "positive", "yes", 0.02);
# Set the adduct list
add.vec <- c("M [1+]","M+H [1+]","M+2H [2+]","M+3H [3+]","M+Na [1+]","M+H+Na [2+]","M+K [1+]")

dir.create("/your/results")
setwd("/your/results")

# Get putative annotations for m/z values
for (i in 1:length(mumFles)){
  mSet <- Read.PeakListData(mSet, paste0(pth,mumFles[i]))
  mSet<-SanityCheckMummichogData(mSet)
  mSet<-Setup.AdductData(mSet, add.vec)
  mSet<-PerformAdductMapping(mSet, "positive")
  mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
  mSet<-SetMummichogPval(mSet, 0.99)
  mSet<-PerformPSEA(mSet, "hsa_kegg", "current", 3 , 100)
  
  file.rename("mummichog_matched_compound_all.csv", 
              paste0(smplNms[i], ".csv"))
  
  
  file.remove("mum_raw.qs", "mum_res.qs",
              "mummichog_pathway_enrichment.csv", "mummichog_query.json",
              "hsa_kegg.qs", "pos_adduct.qs")
}

