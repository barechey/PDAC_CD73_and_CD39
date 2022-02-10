
library(beeswarm)
library(Cairo)

tumor = as.numeric( as.character( read.csv( file = "../data/CD73_soluble.txt" , stringsAsFactors = FALSE , sep = "\t" , header = TRUE , dec = ',' )[ , "CD73_level_homemade_ELISA" ] ) )

hb = read.csv( file = "../data/CD73_soluble_H-B.txt" , stringsAsFactors = FALSE , sep = "\t" , header = TRUE , dec = ',' )
hb$rep1 = as.numeric( as.character( hb$rep1 ) )
hb$rep2 = as.numeric( as.character( hb$rep2 ) )

healthy =  apply( hb[ hb$patient %in% "H" , c( "rep1" , "rep2" ) ] , 1 , mean )
benign =  apply( hb[ hb$patient %in% "B" , c( "rep1" , "rep2" ) ] , 1 , mean )

############################################################################
############################################################################

data = as.data.frame( rbind( cbind( "BT" , tumor ) ,
			cbind( "BT" , benign ) ,
			cbind( "H" , healthy )
			) )
colnames(data) = c( "sample" , "CD73" )
data$sample = as.character( data$sample )
data$sample = factor( data$sample , levels = c( "H" , "BT" ) )
data$CD73 = as.numeric( as.character( data$CD73 ) )

CairoPDF( "../result/Soluble/Boxplot_sCD73.pdf" , height=4, width=3.5, bg="transparent")
	
	################################################
	################################################
	## Plot

	xLabels <- sort( unique( data$sample ) )
    yLabels <- seq( round( min( c( data$CD73 , .4 ) , na.rm=TRUE ) , 1 ) ,  round( max( data$CD73 , na.rm=TRUE ) , 1 ) , by=5 ) 
    boxplot( CD73 ~ sample , data= data , ylab= "sCD73 level" , xlab="" , main="" , 
    	col= "white" , 	    	 
    	boxlty = 1 , outline= FALSE , axes= FALSE , ylim= c( min( c( data$CD73 , .4 ) ) ,  max( data$CD73 )))

	beeswarm( CD73 ~ sample , data= data , 
	        pch = 19, cex = .5 , 
	        corral="wrap", 
	        col=c( adjustcolor("#039be5", alpha.f = 0.9) , 
	           	 adjustcolor("#e53935", alpha.f = 0.9) ), 
	        add = TRUE)

    axis(side = 2, at=yLabels, labels= as.character( yLabels ), las= 2,
             cex.axis=1,tick=1 , col="black")
    axis(side = 1, at=seq( 1 , length(xLabels) , 1 ) , labels= xLabels, las= 1,
             cex.axis=1,tick=1 , col="black")

	wil = wilcox.test( CD73 ~ sample , data = data ) 
	mtext( paste( "Wilcoxon Rank Sum Test " , 
			ifelse( round( wil$p.value , 3 ) < 0.001 ,  paste( "P≤0.001" ) ,
				 paste( "P=" , round( wil$p.value , 3 ) )) , sep="" ) , 
			col="#757575" , 
			cex=.8 
		)

dev.off() 
