##' The pepXML2tab() function generates a data frame from a pepXML file.
##' 
##' Read peptide identification from pepXML file into an data frame object.
##' @title Generate a data frame objects from a pepXML file.
##' @param pepxml a character contains the path and name of a pepXML file 
##' @param ... additional arguments
##' @return a data frame object, each row represent a PSM 
##'     (peptide spectrum match) from the pepXML file
##' @author Xiaojing Wang
##' @examples
##' ## MyriMatch example
##' pepxml <- system.file("extdata/pepxml", "Myrimatch.pepXML", 
##'     package="pepXMLTab")
##' tttt <- pepXML2tab(pepxml)
##' 
##' ## Mascot example
##' pepxml <- system.file("extdata/pepxml", "Mascot.pepXML", 
##'     package="pepXMLTab")
##' tttt <- pepXML2tab(pepxml)
##' 
##' ## SEQUEST example
##' pepxml <- system.file("extdata/pepxml", "SEQUEST.pepXML", 
##'     package="pepXMLTab")
##' tttt <- pepXML2tab(pepxml)
##' 
##' ## XTandem example
##' pepxml <- system.file("extdata/pepxml", "XTandem.pepXML", 
##'     package="pepXMLTab")
##' tttt <- pepXML2tab(pepxml)
##' 


pepXML2tab <- function(pepxml, ...) 
{
    options(stringsAsFactors=FALSE)
    #ifelse(verbose,progress <- "text", progress <- "none")
    doc <- xmlRoot(xmlTreeParse(pepxml))
    msms_run_summary <- doc[['msms_run_summary']]

    sampleEnzyme <- xmlGetAttr(msms_run_summary[["sample_enzyme"]], "name")
    searchEngine <- xmlGetAttr(msms_run_summary[["search_summary"]], 
                    "search_engine")
    searchDB <- msms_run_summary[["search_summary"]][['search_database']]
    
    n <- which(names(msms_run_summary) == "spectrum_query")
    spqrieslist <- msms_run_summary[n]


    spectrumQueries<- lapply(spqrieslist, function(spqu) {
    #spqrieslist[[6]]
        attrs <- xmlAttrs(spqu)
        if(xmlSize(spqu) > 0){
            hit <- xmlApply(spqu[[1]], xmlAttrs)
            
            #xmlApply(spqu[[1]], xmlAttrs)
            
            hitTab <- do.call(rbind, hit)
            
            #totalnum <- xmlApply(spqu, xmlSize)[[1]]
            tmp <- xmlApply(spqu[[1]], function(x) {
                #scores <- xmlApply(x, xmlAttrs)
                scores <- xmlSApply(x, xmlGetAttr, "value")
                scorenames <- xmlSApply(x, xmlGetAttr, "name")
                
                idx <- which(names(scores)=="search_score")
                scoreslist <- unlist(scores[idx])
                names(scoreslist) <- unlist(scorenames[idx])
                
                idxx <- which(names(scores)=="modification_info")
                if(length(idxx) > 0){
                    modification <- as.data.frame(
                            xmlSApply(x[["modification_info"]],
                            xmlAttrs))
                   
                }else{
                    modification <- NA
                }
                res <- c(scoreslist, modification=paste(unlist(modification),
                        collapse=';'))
                res
            })
            
            
            PSM <- cbind(hitTab, do.call(rbind, tmp))
            if(!is.null(PSM)){
                PSM <- t(apply(PSM, 1, function(x) c(attrs, x)))
            }
        }
    })
    
    PSMtab <- do.call(rbind.data.frame, spectrumQueries)
    rownames(PSMtab) <- NULL
    PSMtab
}

