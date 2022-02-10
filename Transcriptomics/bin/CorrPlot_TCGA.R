rm(list=ls(all=TRUE))
load("../data/tcga_PAAD_expression.RData")

library(RColorBrewer)
library(Cairo)
library(corrplot)

######################################################################################
######################################################################################

genes = c("PDCD1", "LAG3" , "CTLA4" , "TIGIT")

## CYT score signature from
## https://clincancerres.aacrjournals.org/content/clincanres/early/2018/09/14/1078-0432.CCR-18-0599.full.pdf
cyt = c("CD8A", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1")
cyt_signature = as.numeric( apply( data[ cyt , ] , 2 , function(x){ mean( x , na.rm=T ) } ) )

## CD8 signature from 
## Jiang et al., Nature Med 2018 : https://www.nature.com/articles/s41591-018-0136-1
cd8 = c("CD8A","CD8B", "GZMA", "GZMB","PRF1")
cd8_signature = as.numeric( apply( data[ cd8 , ] , 2 , function(x){ mean( x , na.rm=T ) } ) )

## Treg signature 
treg = c("FOXP3" , "CCR8")
treg_signature = as.numeric( apply( data[ treg , ] , 2 , function(x){ mean( x , na.rm=T ) } ) )

######################################################################################
######################################################################################
## Corrplot with Genes Signature

dat = rbind( data[ c( "ENTPD1" , "NT5E" ) , ] , cyt_signature , cd8_signature , treg_signature , data[ genes , ] )
rownames(dat) = c(  "ENTPD1" , "NT5E"  , "CYT_signature", "CD8_signature" , "Treg_signature" , genes)


cor.test.estimate = function(mat){
    mat <- as.matrix(mat)
    n <- ncol(mat)
    estimate <- matrix(NA, n, n)
    diag(estimate) <- 1
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            m = mat[ !is.na( mat[ , i ] )& !is.na( mat[ , j ] ) , c( i , j ) ]
            if( nrow( m ) ){
                tmp <- cor.test(m[, 1], m[, 2] , method = 's' , na.rm=TRUE )
                estimate[i, j] <- estimate[j, i] <- tmp$estimate
            }
        }
    }
  colnames(estimate) <- rownames(estimate) <- colnames(mat)
  estimate
}


cor.test.pval = function(mat){
    mat <- as.matrix(mat)
    n <- ncol(mat)
    pval <- matrix(NA, n, n)
    diag(pval) <- 1
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            m = mat[ !is.na( mat[ , i ] )& !is.na( mat[ , j ] ) , c( i , j ) ]
            if( nrow( m ) ){
                tmp <- cor.test(m[, 1], m[, 2] , method = 's' , na.rm=TRUE )
                pval[i, j] <- pval[j, i] <- tmp$p.value
            }
        }
    }
  colnames(pval) <- rownames(pval) <- colnames(mat)
  pval
}
# Matrice de p-value de la corrélation
p.mat <- cor.test.pval( mat = t( dat ) ) 
M = cor.test.estimate( mat = t( dat ) ) 


CairoPDF("../result/CorrPlot/CorrPlot_TCGA_GeneSignature.pdf",height=6,width=6,bg="transparent")
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    corrplot( M, method = "color" , col = col(200) ,  
         type = "upper" ,  
         addCoef.col = "black", 
         tl.col = "black", tl.srt = 45,
         diag = FALSE 
    )
dev.off()

