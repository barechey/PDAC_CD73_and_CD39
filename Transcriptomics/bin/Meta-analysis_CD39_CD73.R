library(survcomp)
library(genefu)
library(Cairo)
library(meta)
library(metafor)
library(grid)

####################################################################################
####################################################################################

unlink( "../result/KMplot/" , recursive = TRUE )
unlink( "../result/Forestplot/" , recursive = TRUE )

dir.create( "../result/KMplot/" )
dir.create( "../result/Forestplot/" )


####################################################################################
####################################################################################

Get_HR_median = function( data , geneID ){
	d = data[ !is.na( data$os ) , ]
	d[ , geneID ] = ifelse( d[ , geneID ] <= median( d[ , geneID ] , na.rm=TRUE) , 0 , 
					ifelse( d[ , geneID ] > median( d[ , geneID ] , na.rm=TRUE) , 1 , NA ) )
	cox = coxph( formula= formula( paste( "Surv( t.os , os ) ~" , geneID ) ) , data = d )
	mean <- round( coef( cox )[1] , 2 )
	low <- round( confint( cox , level=.95 )[ 1 , 1 ] , 2 )
	up <- round( confint( cox , level=.95 )[ 1 , 2 ] , 2 )
	pval <- summary(cox)$coef[ 1 , 5 ]

	c( summary(cox)$n , mean , round( summary(cox)$coefficients[3] , 2 ) , low , up , pval )
}

Generate_Forestplot = function( data , label ){
	data = as.data.frame( data )
	data$study = as.character( data$study )
	data$N = as.numeric(as.character( data$N ))
	data$HR = as.numeric(as.character( data$HR ))
	data$SE = as.numeric(as.character( data$SE ))
	data$Pval = as.numeric(as.character( data$Pval ))

	data = data[ order( data$HR ) , ]

	# Get xlim
	m = c( ( max( abs(data$HR) , na.rm=TRUE) * -1 ) - 1 , ( max( abs(data$HR) , na.rm=TRUE) ) + 1 )

	d = data
	d$study = paste( d$study , ", n = " , data$N , sep= "" ) 

	meta <- metagen( TE = HR,
	                  seTE = SE,
	                  data = d,
	                  studlab = paste(study) ,
	                  fixed = FALSE ,
	                  random = TRUE ,
	                  control = list( maxiter = 10000 , stepadj=0.5 )
	                 ) 

	CairoPDF( label , height= 4, width= 8 , bg="transparent" )
		forest( meta , 
            leftcols = c("studlab", "effect.ci" , "Pval" ),
			leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
   			xlab = "logHR estimate",
			digits.se = 2 ,
   			colgap.forest=unit(10, "mm") ,
	      	plotwidth = unit( 30 , "mm") , 
	       	pooled.totals = TRUE,
	       	smlab = " ",
	       	comb.random =TRUE,
	       	comb.fixed = FALSE,
	       	text.fixed.w = FALSE,
		    layout = "JAMA",
		    print.I2.ci = TRUE,
		    print.Q = FALSE,
		    print.pval.Q = TRUE,
		    print.I2 = TRUE,
		    print.tau2 = FALSE,
		    resid.hetstat = FALSE,
	       	test.overall.random = TRUE,
	       	test.overall.fixed = FALSE,
	       	xlim = m ,  
	       	col.square= "black" ,  
	       	col.study= "black" ,  
	       	col.square.lines = "black" ,
	       	col.diamond.random  = "#1565c0"  ,
	       	col.diamond.lines.random  ="#1565c0" ,
	       	col.by = "#1565c0",
		    addrow.subgroups=TRUE 
	    )
	dev.off()
}

Get_Gene_Kaplan_OS = function( data , geneID1 , geneID2 , study ){

	data = data[ !is.na( data[ ,"t.os" ] ) , ]
	data[ , geneID2 ] = ifelse( data[ , geneID2 ] <= median( data[ , geneID2 ] , na.rm=TRUE) , 0 , 
			ifelse( data[ , geneID2 ] > median( data[ , geneID2 ] , na.rm=TRUE) , 1 , NA ) )

	##########################################
	##########################################
	CairoPDF( paste( "../result/KMplot/Kaplan_Combi_" , geneID1 , "_in_" , geneID2 , "_median_" , study , ".pdf" , sep = "" ) , height = 4 , width = 12 , bg = "transparent")
	  
	  layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))

		##########################################
		##########################################
		par(mar = c(5, 2, 5, 2))
		d = data
		d[ , geneID1 ] = ifelse( d[ , geneID1 ] <= median( d[ , geneID1 ] , na.rm=TRUE) , 0 , 
				ifelse( d[ , geneID1 ] > median( d[ , geneID1 ] , na.rm=TRUE) , 1 , NA ) )

	   	km.coxph.plot( formula.s = formula( paste( "Surv(t.os, os) ~" , geneID1 ) ) , data.s = d , x.label = "Time (Years)" , y.label = "Overall Survival\nProbability" , main.title = paste( geneID1 , "overall" ), sub.title="",
		    leg.text = c( "Low" , "High" ) , leg.pos = "topright" , .col = c( "#3f51b5" , "#e53935" ) , show.n.risk = TRUE , 
		    n.risk.step = 1 , n.risk.cex = 0.85 , ylim = c( 0 , 1 ) , leg.inset = 0 , .lwd = 3 , verbose = FALSE )  


	  	##########################################
		##########################################
		par(mar = c(5, 2, 5, 2))
		
		d = data[ data[ , geneID2 ] %in% 0 , ]
		d[ , geneID1 ] = ifelse( d[ , geneID1 ] <= median( d[ , geneID1 ] , na.rm=TRUE) , 0 , 
				ifelse( d[ , geneID1 ] > median( d[ , geneID1 ] , na.rm=TRUE) , 1 , NA ) )

		km.coxph.plot( formula.s = formula( paste( "Surv(t.os, os) ~" , geneID1 ) ) , data.s = d , x.label = "Time (Years)" , y.label = "Overall Survival\nProbability" , main.title = paste( geneID1 , " in ", geneID2 , " Low", sep = "" ) , sub.title="",
		    leg.text = c( "Low" , "High" ) , leg.pos = "topright" , .col = c( "#3f51b5" , "#e53935" ) , show.n.risk = TRUE , 
		    n.risk.step = 1 , n.risk.cex = 0.85 , ylim = c( 0 , 1 ) , leg.inset = 0 , .lwd = 3 , verbose = FALSE )  

	  	##########################################
	  	##########################################
	  	
	  	par(mar = c(5, 2, 5, 2))
		
		d = data[ data[ , geneID2 ] %in% 1 , ]
		d[ , geneID1 ] = ifelse( d[ , geneID1 ] <= median( d[ , geneID1 ] , na.rm=TRUE) , 0 , 
				ifelse( d[ , geneID1 ] > median( d[ , geneID1 ] , na.rm=TRUE) , 1 , NA ) )

	  	km.coxph.plot( formula.s = formula( paste( "Surv(t.os, os) ~" , geneID1 ) ) , data.s = d , x.label = "Time (Years)" , y.label = "Overall Survival\nProbability" , main.title = paste( geneID1 , " in ", geneID2 , " High" , sep = "" ) , sub.title="",
		    leg.text = c( "Low" , "High" ) , leg.pos = "topright" , .col = c( "#3f51b5" , "#e53935" ) , show.n.risk = TRUE , 
		    n.risk.step = 1 , n.risk.cex = 0.85 , ylim = c( 0 , 1 ) , leg.inset = 0 , .lwd = 3 , verbose = FALSE )  

	dev.off()

}
################################################
################################################

load( "../data/MetaGxPancreas-ICGC.RData" )
rownames(clinical) = clinical$sample

patient = intersect( colnames( data ) , rownames( clinical ) )

mat = as.data.frame( cbind( t( data[ c( "ENTPD1" , "NT5E" ) , patient ] ) , clinical[ patient , ] ) )
colnames(mat) = c( "ENTPD1" , "NT5E" , colnames( clinical ) )

mat[ , "ENTPD1" ] = as.numeric( as.character( mat[ , "ENTPD1" ] ) )
mat[ , "NT5E" ] = as.numeric( as.character( mat[ , "NT5E" ] ) )
mat[ , "dataset" ] = as.character( mat[ , "dataset" ] )
mat[ , "os" ] = as.numeric( as.character( mat[ , "os" ] ) )
mat[ , "t.os" ] = as.numeric( as.character( mat[ , "t.os" ] ) )

####################################################################################################################################################
####################################################################################################################################################

dataset = sort( unique( mat$dataset ) )

CD39 = CD73 = NULL
for( i in 1:length( dataset ) ){
	m = mat[ mat$dataset %in% dataset[i] , ]

	CD39 = rbind( CD39 , c( dataset[i] , Get_HR_median( data = m , geneID = "ENTPD1" ) ) )

	CD73 = rbind( CD73 , c( dataset[i] , Get_HR_median( data = m , geneID = "NT5E" ) ) )

	if( dataset[i] %in% c( "ICGCMICRO_SumExp" , "TCGA_SumExp" ) ){
		
		Get_Gene_Kaplan_OS( data = m , geneID1 = "NT5E" , geneID2 = "ENTPD1" , study = dataset[i] )

	}

}

colnames( CD39 ) = colnames( CD73 ) = c( "study" , "N" , "HR" , "SE" , "95low" , "95high" , "Pval")

##########################################
##########################################

Generate_Forestplot( data = CD39 , label = "../result/Forestplot/Meta-analysis_CD39.pdf" )

Generate_Forestplot( data = CD73 , label = "../result/Forestplot/Meta-analysis_CD73.pdf" )


