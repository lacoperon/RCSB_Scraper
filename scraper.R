library(XML)
library(readr)
library(httr)
library(dplyr)

# Opens pre-saved version of the search results
searchResults <- read_file("data/searchResults.html")
d <- htmlParse(searchResults)

# Gets the number of structures returned
node_num_struct <- getNodeSet(doc=d, path="//ul/li[@class='active']//a")[[1]]
num_structure <- trimws(xmlToList(node_num_struct)$text)

# Scrapes all search result HTML nodes (and their PDB accession numbers)
node_set_path <- "//ul[@id='SearchResultsDetails-MainContent']/li"
result_nodes <- getNodeSet(doc=d, path=node_set_path)
# Returns the class element of each search result item corresponding to PDB code
cand_structs  <- data.frame(sapply(result_nodes, 
                                   function(node) {
                                     substr(xmlAttrs(node)[[1]], 11, 1000)
                                   }))
colnames(cand_structs) <- c("PDB_Number")

# Scrapes brief description of structure
desc_path <- "//ul[@id='SearchResultsDetails-MainContent']/li/div[2]/h4/a"
desc_nodes <- getNodeSet(doc=d, path=desc_path)
cand_structs$Description  <- sapply(desc_nodes, 
                                    xmlValue)

# Scrapes all search result links (to the PDB structures themselves)
file_link_path <- "//ul[@id='SearchResultsDetails-MainContent']/li/div[2]/div[1]/div/a[1]"
cand_structs$PDBFileLink   <- xpathSApply(doc=d, path=file_link_path, fun=xmlGetAttr, name="href")


# Scrapes Citation of Paper
cit_path <- "//ul[@id='SearchResultsDetails-MainContent']/li/div[2]/p[2]"
cit_nodes <- getNodeSet(doc=d, path=cit_path)
cand_structs$PaperDate  <- sapply(cit_nodes, xmlValue)

# Scrapes Link to Individual Structure Page
link_path <- "//ul[@id='SearchResultsDetails-MainContent']/li/div[2]/h3/a"
cand_structs$PageLink <- xpathSApply(doc=d, path=link_path, fun=xmlGetAttr, name="href")

# This function returns FASTA link associated with a particular rcsb page
getFASTALink <- function(url) {
  # Gets HTML associated with particular site
  print(url)
  s <- GET(url)
  w <- content(s, as='text') # converts s to plaintext of HTML
  d <- htmlParse(file=content(s, as="text", asText=T))
  
  #XPath associated with FASTA download from page
  fl_path <- '//*[@id="DownloadFilesButton"]/ul/li[1]/a'
  fasta_link <- xpathSApply(doc=d, path=fl_path, fun=xmlGetAttr, name="href")[1]
  if (nchar(fasta_link) == 0) {
    # Tries again (one more time), if necessary
    fasta_link <- xpathSApply(doc=d, path=fl_path, fun=xmlGetAttr, name="href")[1]
  }
  fasta_link <- paste("http://www.rcsb.org", fasta_link, sep="")
  return(fasta_link)
}

# Uses the above function to get the FASTA links for all candidate structures
cand_structs$FastaLink <- sapply(cand_structs$PageLink, getFASTALink)









