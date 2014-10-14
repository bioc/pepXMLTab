test_pepXML2tab <- function() {
    pepxml <- system.file("unitTests", "tmp.pepXML", 
        package="pepXMLTab")
    checkEquals(dim(pepXML2tab(pepxml)), c(0, 0))
    pepxml <- system.file("extdata/pepxml", "Myrimatch.pepXML", 
        package="pepXMLTab")
    checkTrue(!is.na(pepXML2tab(pepxml)[1,1]))
    checkEqualsNumeric(dim(pepXML2tab(pepxml))[2], 25)
}


