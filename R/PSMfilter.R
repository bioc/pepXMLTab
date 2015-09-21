##' The PSMfilter() function filter the peptide identification based on user 
##' chosen paramter.
##' 
##' Filter the peptide identification based on FDR, hit rank, or peptide length.
##' @title Filter the peptide indentification.
##' @param PSMtab a data frame contain peptide identification from a pepXML file
##' @param pepFDR filter the peptides based on this chosen FDR, 
##' default is 0.01. 
##' @param scorecolumn which column is chosen to calculate FDR 
##' @param hitrank an integer indicates how many peptides to retain for a 
##' spectrum. A spectrum can match to multiple peptides. Default is 1.
##' @param minpeplen an integer of minimum peptide length 
##' @param decoyprefix a character indicates decoy sequence in the 'protein' 
##' column. Usually is 'rev_' or 'DECOY_'.       
##' @param ... additional arguments
##' @return a data frame object, contain PSMs (peptide spectrum match) passed 
##' the filters.
##' @author Xiaojing Wang
##' @examples
##' ##MyriMatch example
##' pepxml <- system.file("extdata/pepxml", "Myrimatch.pepXML", 
##'             package="pepXMLTab")
##' tttt <- pepXML2tab(pepxml)
##' passed <- PSMfilter(tttt, pepFDR=0.01, scorecolumn='mvh', hitrank=1, 
##'         minpeplen=6, decoyprefix='rev_')
##' 

PSMfilter <- function(PSMtab, pepFDR=0.01, scorecolumn='mvh', hitrank=1, 
        minpeplen=6, decoyprefix='rev_', ...)
{
    ##hitrank
    if(dim(PSMtab)[1] == 0){
        message("No input data")
    }else{
        #PSMtab <- subset(PSMtab, hit_rank <= hitrank)
        PSMtab <- PSMtab[PSMtab$hit_rank <= hitrank, ]
        
        ##peptide length
        peptide <- PSMtab[, 'peptide']
        idx <- which(nchar(peptide) < minpeplen)
        if(length(idx) > 0) PSMtab <- PSMtab[-idx, ]
        
        ##peptide hit to reverse protein
        ## too slow
        #peptide <- PSMtab[, 'peptide']
        #decoy <- AAString(paste(reverse(proteinseq[, 'peptide']), 
        #    collapse=';'))
        #system.time(idx <- lapply(peptide[1:200], function(x) 
        #                countPattern(x, decoy, fixed=TRUE)))
        
        NTT <- rep(0, dim(PSMtab)[1])
        cut1 <- grep('[KR-]', PSMtab[, 'peptide_prev_aa'], fixed=FALSE)
        NTT[cut1] <- 1
        cut2 <- union(grep('[KR]', 
                    substring(PSMtab[, 'peptide'], nchar(PSMtab[, 'peptide'])), 
                    fixed=FALSE),
                grep('-', PSMtab[, 'peptide_next_aa'], fixed=TRUE))
        NTT[cut2] <- 1
        NTT[intersect(cut1, cut2)] <- 2
        PSMtab <- cbind(PSMtab, NTT)
        
        PSMtabssubs <- split(PSMtab, paste(PSMtab$assumed_charge, 
                                PSMtab$NTT))
        
        PSMsubpass <- lapply(PSMtabssubs, function(y){
            protein <- y[, 'protein']
            index <- grep(decoyprefix, protein, fixed=TRUE)

            if(dim(y)[1] == 1){
                if(length(index) == 0) res <- y
            }else{
                proORrev <- rep(0, dim(y)[1])
                index <- grep(decoyprefix, protein, fixed=TRUE)
                proORrev[index] <- 1
                
                #y <- cbind(y, proORrev)
                #y <- y[order(as.numeric(y[, scorecolumn]), decreasing=TRUE),]
                #calculate FDR based on score
                score <- y[, scorecolumn]
                score <- as.numeric(score)
            
                tmp <- cbind(score, proORrev)
                tmp <- tmp[order(score, decreasing=TRUE),]
                FP <- cumsum(tmp[, 2])
                tmp <- cbind(tmp, FP)
                calFDR <- unlist(lapply(1:dim(y)[1], function(x) 
                        (2*tmp[x, 'FP'])/x))
                #y <- cbind(y, calFDR)
                #res <- subset(y, calFDR <= pepFDR)
                if(TRUE %in% (calFDR <= pepFDR)){
                    cutoff <- tmp[max(which(calFDR <= pepFDR)), 'score']
                    #print(cutoff)
                    res <- y[which(as.numeric(y[, scorecolumn]) >= cutoff),]
                    res
                }
            }
        })
        
        PSMpass <- do.call(rbind.data.frame, PSMsubpass)
        rownames(PSMpass) <- NULL
        PSMpass
    }
}
