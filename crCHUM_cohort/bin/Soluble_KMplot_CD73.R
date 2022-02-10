
library(Cairo)
library(forestplot)
library(sva)
library(genefu)
 
data = read.csv( file = "../data/CD73_soluble.txt" , stringsAsFactors = FALSE , sep = "\t" , header = TRUE , dec = ',' )
rownames(data) = paste( "P" , data[ , "Nbr_bank" ] , sep = "" )


id = c( "CD73_level_homemade_ELISA" , 
		"OS_months" , "OS_censure" )

data = data[ , id ]
colnames( data ) = c( "sCD73" ,
		"os.time" , "os" )

data[ , "os.time" ] = as.numeric( as.character( data[ , "os.time" ] ) ) / 12

for( i in 1:nrow( data ) ){

	if( !is.na( as.numeric( as.character( data[ i , "os.time" ] ) ) ) && as.numeric( as.character( data[ i , "os.time" ] ) ) > 5 ){
		data[ i , "os.time" ] = 5
		data[i , "os" ] = ifelse( is.na( data[ i , "os" ] ) , NA , 0 )
	}

}

data[ , "sCD73" ] = ifelse( data[ , "sCD73" ] >= quantile( data[ , "sCD73" ] , probs = .75 , na.rm = TRUE ) , 1 , 0 )

	
pdf( "../result/Soluble/KMplot_sCD73_OS_topquartile.pdf" , 
			height = 5 , width = 7 , bg = "transparent" , onefile = FALSE )
	km.coxph.plot( formula.s = Surv( os.time , os ) ~ sCD73 , data.s = data , x.label = "Time (Years)" , y.label="Overall Survival" , 
		main.title = "sCD73" , sub.title = "" , 
		leg.text = c(  "sCD73 Low" , "sCD73 High" ) , 
		leg.pos = "topright" , .col = c( "#1976d2" , "#e53935" ) , show.n.risk = TRUE , n.risk.step = 1 , n.risk.cex = 0.85 , ylim = c( 0 , 1 ) , 
		leg.inset = 0 , .lwd = 3.5 , verbose = FALSE )
dev.off()



