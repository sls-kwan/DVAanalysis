require(RSelenium)
require(RUnit)
require(audio)
##Find the RSelenium webdriver
##Run first time
RSelenium::checkForServer()

#You need to load the selenium driver if you want to load to server
#Driver can be found @ https://sites.google.com/a/chromium.org/chromedriver/downloads

##Windows Version
#Must use the path to chrome driver
# RSelenium::startServer(args = c("-Dwebdriver.chrome.driver=C:/Users/SLS/Downloads/chromedriver.exe")
#                        , log = FALSE, invisible = F)
##Linux/Mac Version, more reliable to run "java -jar selenium-server-standalone-2.48.2.jar"
RSelenium::startServer(args = c("-Dwebdriver.chrome.driver=/home/shannon/R/x86_64-pc-linux-gnu-library/3.2/RSelenium/bin/chromedriver")
                       , log = FALSE, invisible = FALSE)
remDr <- remoteDriver(remoteServerAddr = "localhost" 
                      , port = 4444
                      , browserName = "chrome")


##Function allows user to do GO analysis using the DAVID website
##inputFile is the inputfile that contains the list of genes that need to have analysis on them
##Make sure the file has content
##geneType is the gene identifier
##Species is the species of interest(M for mouse and everthing else will be use human background )
AutoSubmit("bothDownGSE51372", "OFFICIAL_GENE_SYMBOL", "M")
# AutoSubmit("bothUpE-MTAB-2805", "ENSEMBL_GENE_ID", "M")
dataLookat <- read.csv("~/IT/DVA_2.0/EMTAB2805", header=FALSE)
AutoSubmit("bothDownGSE67835", "OFFICIAL_GENE_SYMBOL", "H")
AutoSubmit("UpdvE-MTAB-2805", "ENSEMBL_GENE_ID", "H")
apply(dataLookat, 1, function(x){
  x <- as.character(x)
  AutoSubmit(x[1], x[2], x[3])
})
AutoSubmit <- function(inputFile, geneType, species){
  #Opens browser window and navigates to DAVID site
  remDr$open()
  remDr$navigate("https://david.ncifcrf.gov/tools.jsp")
  #Selects upload tab
  UploadHome <- remDr$findElement(using='xpath', '//*[@id="tablist"]/li[1]/a')$clickElement()
  #Upload inputFile into and selects the gene name based on geneType
  textBox <- remDr$findElement(using ='name', value='fileBrowser')
  filesIn <- as.list(paste("/home/shannon/IT/DVA_2.0/Results/", inputFile,".txt", sep=""))
  textBox$sendKeysToElement(filesIn)
  GeneSymbol <- remDr$findElement(using ='xpath', paste("//*/option[@value='", geneType, "']" ,sep=""))$clickElement()
  GeneList <- remDr$findElement(using='name', value='rbUploadType')$clickElement()
  SubmitList <- remDr$findElement(using='xpath', value='//*[@id="divUpload"]/table/tbody/tr[16]/td/input')$clickElement()
  ##Test to see if the valid genes have been selected
  testValid <- try(unlist(remDr$findElement(using='xpath', '//*[@id="conversion"]/table/tbody/tr[4]/td/b/input')))
  if(!("try-error" %in% class(testValid))){
      validGenes <- unlist(remDr$findElement(using='xpath', '//*[@id="conversion"]/table/tbody/tr[4]/td/b/input'))
      validGenes$clickElement()
  } else {
    testExist <- try(remDr$findElement(using='xpath', '//*[@id="conversion"]/table/tbody/tr[4]/td/span/input'))
    if(!("try-error" %in% class(testExist))){
        print("There is no gene IDS in file")
        return()
    }
  }
  alertExist <- try(remDr$getAlertText())
  if(!("try-error" %in% class(alertExist))){
    remDr$acceptAlert()
  }
  #Select spcies and submit

  #Select the most recent entry of the background
  BackgroundTab <- remDr$findElement(using ='xpath', "//*[@id='tablist']/li[3]/a")$clickElement()
  ##This selects the mouse background
  if(species =="M"){
    BackgroundSelect <- remDr$findElement(using='xpath', '(//*[@id="affyBGs"])[5]/font/font/input[8]')$clickElement()
  } else {
  ##This selects the human background
    BackgroundSelect <- remDr$findElement(using='xpath', '(//*[@id="affyBGs"])[5]/font/font/input[1]')$clickElement()
  }
  #Get the functional chart
  FunctionChart <- remDr$findElement(using='xpath', '/html/body/table[2]/tbody/tr/td[3]/table/tbody/tr/td/table/tbody/tr/td/table/tbody/tr/td/table/tbody/tr[5]/td[2]/a/ul/li[2]/font/small/big')$clickElement()
  #Use mouse click instead of clickElement() as it can select only that element, not only not that element
  SelectSpecies <- remDr$findElement(using='xpath', '//*[@id="divManager"]/table[1]/tbody/tr[3]/td/p[1]/font/select/option[@value= "0"]')
  remDr$mouseMoveToLocation(webElement = SelectSpecies)
  remDr$click(0)
  submitSpecies <- remDr$findElement(using='name', value='B1')$clickElement()
  clearAll <- remDr$findElement(using='xpath', '/html/body/table[2]/tbody/tr/td[3]/table/tbody/tr/td/table/tbody/tr/td/table/tbody/tr/td/form/table[2]/tbody/tr[2]/td[3]/input')$clickElement()
  #Clear all selection and select GO default terms and all pathway terms
  GeneOntology <- try(remDr$findElement(using='partial link text', value='Gene_Ontology'))
  #Selects default GO terms
  if(!("try-error" %in% class(GeneOntology))){
    GeneOntology <- remDr$findElement(using='partial link text', value='Gene_Ontology')
    GeneOntology$clickElement()
    GOTerms <- remDr$findElements(using='id', value='Gene_Ontology')
    GOindexes <- c("25","32","39")
    lapply(GOTerms, function(x){
      value <- x$getElementAttribute("value")
      presence <- any(grepl(unlist(value),GOindexes, fixed=T))
      if(presence){x$clickElement()}
    })
  }
  #Selects all GO pathway terms
  PathwaysOpen <- try(remDr$findElement(using='id', 'Pathwaystd')$clickElement())
  if(!("try-error" %in% class(PathwaysOpen))){
    Pathways <-remDr$findElements(using='id', value='Pathways')
    lapply(Pathways, function(x){
      x$clickElement()
    })
  }
  #Get chart
  GetChart <- remDr$findElement(using='xpath', '/html/body/table[2]/tbody/tr/td[3]/table/tbody/tr/td/table/tbody/tr/td/table/tbody/tr/td/form/table[3]/tbody/tr[3]/td/button')$clickElement()
  #Need to switch to new window
  allWindows <- unlist(remDr$getWindowHandles())
  switchWindow <- remDr$switchToWindow(allWindows[[2]])
  #Waits for the page to load by checking url and progress bar status
  x <- remDr$getCurrentUrl()[[1]]
  while(x == ""){x <- remDr$getCurrentUrl()[[1]]}
  progress <- remDr$findElement(using="id", value='progress')
  percentProgress <- progress$getElementText()[[1]]
  while(percentProgress != ""){
    progress <- remDr$findElement(using="id", value='progress')
    percentProgress <- progress$getElementText()[[1]]
  }
  RunOptions <- remDr$findElement(using='xpath', '//*[@id="summaryTree"]/li/a')$clickElement()
  #Changes Gene Count to 10
  CountOptions <- remDr$findElement(using='name', value='count')
  CountOptions$sendKeysToElement(list('\ue009', 'a', '\ue003'))
  CountOptions$sendKeysToElement(list('10'))
  #Get the gene list and saves it
  GetList <- remDr$findElement(using='xpath', '/html/body/table/tbody/tr[2]/td/table[2]/tbody/tr/td/input[1]')$clickElement()
  SaveList <- try(remDr$findElement(using='xpath', '/html/body/table/tbody/tr[2]/td/table[3]/tbody/tr/td[2]/font/a'))
  if("try-error" %in% class(SaveList)){
    print("No results");remDr$close();return();
  }
  SaveList <- remDr$findElement(using='xpath', '/html/body/table/tbody/tr[2]/td/table[3]/tbody/tr/td[2]/font/a')
  SaveList$clickElement()
  allWindows <- unlist(remDr$getWindowHandles())
  remDr$switchToWindow(allWindows[3])
  x <- remDr$getCurrentUrl()[[1]]
  while(x == ""){x <- remDr$getCurrentUrl()[[1]]}
  statusCheck <- function(){
   return(remDr$executeScript('return document.readyState')[[1]])
  }
  statuscheck <- statusCheck()
  while (statuscheck != "complete"){
    statuscheck <- statusCheck();
  }
  output <- remDr$getPageSource()[[1]]
  name <- paste("GO", inputFile,".tsv", sep="")
  source("htmlToText.R")
  output <- htmlToText(output)
  fileName <- paste("~/IT/DVA_2.0/GO/",name, sep='')
  fileConn <- file.create(fileName)
  write(output, file=fileName)
  remDr$deleteAllCookies()
  remDr$close()
}
