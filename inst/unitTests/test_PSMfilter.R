test_PSMfilter <- function() 
{
    pepxml <- system.file("unitTests", "tmp.pepXML", 
        package="pepXMLTab")
    checkEquals(dim(pepXML2tab(pepxml)), c(0, 0))
    checkEquals(PSMfilter(pepXML2tab(pepxml)), message('No input data'))
}

