########################################
## Run script from the crCHUM cohort
########################################

unlink( "../result/Soluble/" , recursive = TRUE )
unlink( "../result/TMA/" , recursive = TRUE )

dir.create( "../result/Soluble/" )
dir.create( "../result/TMA/" )

############################
## Soluble CD73 analysis
############################

source( "Soluble_KMplot_CD73.R" )
source( "Boxplot_CD73_soluble.R" )


############################
## TMA CD73, CD39, CD8 analysis
############################

source( "KMplot_1_variable_analysis.R" )
source( "KMplot_2_variables_median.R" )
