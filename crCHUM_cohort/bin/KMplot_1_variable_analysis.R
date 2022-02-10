
library(Cairo)
library(forestplot)
library(DescTools)
library(sva)
library(genefu)


############################################################################
############################################################################
## Creating the directories

unlink( "../result/TMA/KMplot_1_variable/" , recursive = TRUE )

dir.create( "../result/TMA/KMplot_1_variable" )

############################################################################
############################################################################
## Reading the file

data = read.csv( file = "../data/TMA_CD73_CD8_CD39.txt" , stringsAsFactors = FALSE , sep = "\t" , header = TRUE , dec = ',' )

id = c( "MFI_CD73_epithelium_tumeur" ,
		"MFI_CD39_Stroma_tumor" , 
		"CD8_cm2_total_tumor" , 
		"Overall_survival_time" , "Death_Status" )

data = data[ , id ]
colnames( data ) = c( "CD73_E_T" , "CD39_S_T" , "CD8_cm_tot" , "os.time" , "os" )

############################################################################
############################################################################
## Censuring survival data to 5-years

data[ , "os.time" ] = as.numeric( as.character( data[ , "os.time" ] ) ) / 12

for( i in 1:nrow( data ) ){

	if( !is.na( as.numeric( as.character( data[ i , "os.time" ] ) ) ) && as.numeric( as.character( data[ i , "os.time" ] ) ) > 5 ){
		data[ i , "os.time" ] = 5
		data[ i , "os" ] = ifelse( is.na( data[ i , "os" ] ) , NA , 0 )
	}
}


############################################################################
############################################################################
## Use of the top-tertile cut-point for CD73 & CD39 variables

id = c( "CD73_E_T" , "CD39_S_T" )

for(i in 1: length( id ) ){
	data[ , id[ i ] ] = ifelse( round( data[ , id[ i ] ] ) > round( quantile( data[ , id[ i ] ] , probs = .66 , na.rm = TRUE ) ) , 1 , 
												ifelse( is.na( data[ , id[ i ] ] ) , NA , 0 ) )
}

############################################################################
############################################################################
## Use of the optimal cut-point for CD8 variable

data$CD8_cm_tot = as.numeric( as.character( data$CD8_cm_tot ) )

library(survminer)

os.cut <- surv_cutpoint(
	data,
	time = "os.time",
	event = "os",
	variables = "CD8_cm_tot"
)
summary(os.cut)

data$CD8_cm_tot = ifelse( round( data$CD8_cm_tot ) > summary(os.cut)[1,1] , 1 , 
 												ifelse( round( data$CD8_cm_tot ) <= summary(os.cut)[1,1] , 0 , NA ) )


############################################################################
############################################################################
## Generating the figures

id = c( "CD73_E_T" , "CD39_S_T" , "CD8_cm_tot" )

for( i in 1 : length( id ) ){

	pdf( paste( "../result/TMA/KMplot_1_variable/KMplot_" , id[ i ] , "_OS_toptertile.pdf" , sep = "" ) , 
				height = 5 , width = 7 , bg = "transparent" , onefile = FALSE )
		# par( mar = c( 5 , 12 , 1 , .5 ) )
		km.coxph.plot( formula.s = formula( paste( "Surv( os.time , os ) ~" , id[i] ) ) , data.s = data , x.label = "Time (Years)" , y.label="Overall Survival" , 
			main.title = id[ i ] , sub.title = "" , 
			leg.text = c( paste( id[ i ] , "Low" ) , paste( id[ i ] , "High" ) ) , 
			leg.pos = "topright" , .col = c( "#1976d2" , "#e53935" ) , show.n.risk = TRUE , n.risk.step = 1 , n.risk.cex = 0.85 , ylim = c( 0 , 1 ) , 
			leg.inset = 0 , .lwd = 3.5 , verbose = FALSE )
	dev.off()

}


