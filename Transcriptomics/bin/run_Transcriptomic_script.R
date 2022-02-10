########################################
## Run script from the Transcriptomic cohort
########################################

unlink( "../result/CorrPlot/" , recursive = TRUE )
dir.create( "../result/CorrPlot/"  )


############################
## Soluble CD73 analysis
############################

source( "CorrPlot_ICGC.R" )
source( "CorrPlot_TCGA.R" )


############################
## TMA CD73, CD39, CD8 analysis
############################

source( "Meta-analysis_CD39_CD73.R" )
