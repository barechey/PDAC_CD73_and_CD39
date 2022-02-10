rm(list=ls(all=TRUE))
library(genefu)
library(MetaGxPancreas)

pancreasData <- loadPancreasDatasets()[1]

esets = pancreasData$SummarizedExperiments

datasets = names(esets)

data = geneID = clinical = NULL

for(i in 1:length(datasets)){
  expr = assays(esets[[i]])[[1]]
  colnames( expr ) = paste( datasets[i] , colnames( expr ) , sep = "_" )
  #rownames(expr) = featureData(esets[[i]])$gene

  clinic = as.data.frame(colData(esets[[i]]))
  rownames( clinic ) = paste( datasets[i] , rownames( clinic ) , sep = "_" )

  expr = expr[,ifelse( clinic$sample_type%in%"tumour", TRUE, FALSE)]
      
  if(ncol(expr)>=40){
    #############################################################################
    #############################################################################

    expr.sd <- apply(expr,1,sd, na.rm=T)

    t_uniq <- expr[!(rownames(expr)%in%rownames(expr[duplicated(rownames(expr)),])),]
    t_dup <- expr[(rownames(expr)%in%rownames(expr[duplicated(rownames(expr)),])),]
    
    if(nrow(t_dup)>0){
      t_dup <- t_dup[order(rownames(t_dup)),]
      id <- unique(rownames(t_dup))

      t.dup.rm <- NULL
      names <- NULL
      for(j in 1:length(id)){
        tmp <- t_dup[which(rownames(t_dup)%in%id[j]),]
        tmp.sd <- apply(tmp,1,function(x){sd(as.numeric(as.character(x)),na.rm=T)})
        tmp <- tmp[which(tmp.sd%in%max(tmp.sd,na.rm=T)),]

        if( is.null(dim(tmp)) ){
          t.dup.rm <- rbind(t.dup.rm,tmp) 
          names <- c(names,names(tmp.sd)[1])
        }   
      }
      expr <- rbind(t_uniq,t.dup.rm)
      rownames(expr) <- c(rownames(t_uniq),names)
    }


    #############################################################################
    #############################################################################
    clin = as.data.frame( cbind(
                  as.numeric(as.character( clinic$days_to_death)), 
                  ifelse(as.character( clinic$vital_status)%in%"living",0,
                  ifelse(as.character( clinic$vital_status)%in%"deceased",1,NA))))
    colnames(clin) = c("t.os","os")

    if(!is.null(geneID)){
      geneID = intersect(rownames(expr),rownames(data))
    } else { 
      geneID = rownames(expr)
    }


    
    if(sum(clin[,"t.os"],na.rm=T)>0){

      clin = cbind( clinic$sample_name , datasets[i] , clin )
      colnames( clin ) = c( "sample" , "dataset" , "t.os" , "os" )
      clin[ , "sample" ] = as.character( clin[ , "sample" ] )
      clin[ , "dataset" ] = as.character( clin[ , "dataset" ] )
      clin[ , "t.os" ] = as.numeric( as.character( clin[ , "t.os" ] ) ) / 365
      clin[ , "os" ] = as.numeric( as.character( clin[ , "os"] ) ) 

      clin = clin[ clinic$sample_type %in% "tumour" , ]
      clin = clin[ !is.na( clin[ , "t.os" ] ) , ]

      if( nrow( clin ) >= 40 ){ 

        if( is.null( data ) ){
          data = expr[ geneID , ]
        } else{
          data = rbind(
                    cbind(
                      data[ geneID , ] , expr[ geneID , ] 
                    ) ,
                    cbind( 
                      data[ !rownames( data ) %in% geneID , ] , 
                      matrix( NA , nrow = nrow( data[ !rownames( data ) %in% geneID , ] ) , ncol = ncol( expr[ !rownames( expr ) %in% geneID , ] ) )
                    ) ,
                    cbind(
                      matrix( NA , nrow = nrow( expr[ !rownames( expr ) %in% geneID , ] ) , ncol = ncol( data[ !rownames( data )%in% geneID , ] ) ) ,
                      expr[!rownames(expr)%in%geneID,]
                    )
                )
        } 

        clinical = rbind( clinical , clin )
      }
    }
  }
}

clinical[ , "sample" ] = paste( clinical[ , "dataset" ] , clinical[ , "sample" ] , sep = "_" ) 
save( data , clinical , file = "../data/MetaGxPancreas.RData" )



##########################################
##########################################

load( "../data/MetaGxPancreas.RData" )
library(genefu)
data_meta = data
clinical_meta = clinical

##########################################
##########################################
load( "../data/icgc_PAAD_ca_expression.RData" )

data_ca = normCount

clin_ca = as.data.frame( cbind( rownames( clin ) , "ICGC_CA" , clin[ colnames( data_ca ) , c( "donor_survival_time" , "donor_vital_status" ) ] ) )
colnames( clin_ca ) <- c( "sample","dataset","t.os","os")

clin_ca[ , "os" ] <- ifelse( as.character( clin_ca[ , "os" ] ) %in% "deceased" , 1 , 
              ifelse( as.character( clin_ca[ , "os" ] ) == "alive" , 0 , NA ) )

clin_ca[ , "sample" ] <- as.character( clin_ca[ , "sample" ] )
clin_ca[ , "dataset" ] <- as.character( clin_ca[ , "dataset" ] )
clin_ca[ , "t.os" ] <- as.numeric( as.character( clin_ca[ , "t.os" ] ) ) / 365

colnames( data_ca ) = paste( "ICGC_CA" , colnames( data_ca ) , sep = "_" )
rownames( clin_ca ) = paste( "ICGC_CA" , rownames( clin_ca ) , sep = "_" )
clin_ca$sample =  rownames( clin_ca ) 

######################################################################################################################################################
######################################################################################################################################################

geneID = intersect( rownames( data_meta ) , rownames( data_ca ) )

data = cbind( data_meta[ geneID , ] , data_ca[ geneID , ] )
clinical = rbind( clinical_meta , clin_ca )

######################################################################################################################################################
######################################################################################################################################################

for(i in 1:nrow( clinical ) ){
  if( !is.na( as.numeric( as.character( clinical[ i , "t.os" ] ) ) ) && as.numeric( as.character( clinical[ i , "t.os" ] ) ) > 5 ){
    clinical[ i , "t.os" ] = 5
    clinical[ i , "os" ] = 0
  }
}

rownames( clinical ) = clinical$sample

######################################################################################################################################################
######################################################################################################################################################
save(data, clinical, file="../data/MetaGxPancreas-ICGC.RData")

