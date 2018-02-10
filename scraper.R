library(XML)
library(readr)
library(httr)
library(dplyr)
library(purrr)

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

getLinkFromPage <- function(url, path, prepend) {
  # Gets HTML associated with particular site
  print(url)
  s <- GET(url)
  w <- content(s, as='text') # converts s to plaintext of HTML
  d <- htmlParse(file=content(s, as="text", asText=T))
  
  # Gets link from page
  link <- xpathSApply(doc=d, path=path, fun=xmlGetAttr, name="href")[1]
  # Sees that it's valid (ie is not empty)
  if (nchar(link) == 0) {
    # Tries again (one more time), if necessary
    s <- GET(url)
    w <- content(s, as='text') # converts s to plaintext of HTML
    d <- htmlParse(file=content(s, as="text", asText=T))
    link <- xpathSApply(doc=d, path=path, fun=xmlGetAttr, name="href")[1]
  }
  link <- paste(prepend, link, sep="")
  return(link)
}

fl_path <- '//*[@id="DownloadFilesButton"]/ul/li[1]/a'
fl_prepend <- "http://www.rcsb.org"
getFASTALink <- partial(getLinkFromPage, path=fl_path, prepend=fl_prepend)
# Uses the above function to get the FASTA links for all candidate structures
cand_structs$FastaLink <- sapply(cand_structs$PageLink, getFASTALink)


# TODO: Method, Resolution

getAttrsFromPage <- function(url, path) {
  print(url)
  s <- GET(url)
  w <- content(s, as='text') # converts s to plaintext of HTML
  d <- htmlParse(file=content(s, as="text", asText=T))
  nodes <- getNodeSet(doc=d, path=path)
  values <- sapply(nodes, xmlValue)
  return(values)
}

getAttrFromPage <- function(url, path) {
  values <- getAttrsFromPage(url, path)
  if (is.null(values[1])) {
    values <- getAttrsFromPage(url, path)
  }
  return(values[1])
}

# Gets Authors associated with particular structure
author_path <- '//*[@id="header_deposition-authors"]'
getAuthorNames <- partial(getAttrFromPage, path=author_path)
cand_structs$Authors <- substr(sapply(cand_structs$PageLink, getAuthorNames), 27, 200)

# Gets Method used to derive structure
method_path <- '//*[@id="exp_header_0_method"]'
getMethod <- partial(getAttrFromPage, path=method_path)
cand_structs$Method <- substr(sapply(cand_structs$PageLink, getMethod), 13, 100)

# Gets Resolution of Structure
diff_path <- '//*[@id="exp_header_0_diffraction_resolution"]'
em_path <- '//*[@id="exp_header_0_em_resolution"]'
getDiff <- partial(getAttrFromPage, path=diff_path)
getEM   <- partial(getAttrFromPage, path=em_path)
# Gets resolution (the tagged Angstrom resolution varies by diffraction or Em)
getRes <- function(url) {
  diff_path <- '//*[@id="exp_header_0_diffraction_resolution"]'
  res <- getDiff(url)
  if (is.null(res[[1]])) {
    em_path <- '//*[@id="exp_header_0_em_resolution"]'
    res <- getEM(url)
  }
  return(res)
}

# Parses Resolution Text to numeric
res <- sapply(cand_structs$PageLink, getRes)
res_parsed <- sub("Resolution:", "", res)
res_parsed <- gsub("&nbsp", "", res_parsed)
res_parsed <- gsub("Å", "", res_parsed)
res_parsed <- as.numeric(res_parsed)
cand_structs$Resolution <- res_parsed

write_csv(cand_structs, "./data/candStructs.csv")

# Path associated with 
DNA_path <- '//*[@id="MacromoleculeTableDNA"]/div/table/tbody/tr[3]'
# Function to get the constituent Macromolecules in the structure
# (In this case, I'll only use it for the RNA Parts, but it can also
# be applied to Protein as well)
getConstituentEntities <- function(url, path) {
  print(url)
  s <- GET(url)
  w <- content(s, as='text') # converts s to plaintext of HTML
  d <- htmlParse(file=content(s, as="text", asText=T))
  nodes <- getNodeSet(doc=d, path=path)
  nodeChildren <- sapply(nodes, xmlChildren)
  data_list <- sapply(nodeChildren, xmlValue)
  if (length(data_list) == 0) {
    return(NULL)
  }
  EntityName <- data_list[seq(1, length(data_list), 4)]
  entities <- data.frame(EntityName)
  entities$ChainNames <- data_list[seq(2, length(data_list), 4)]
  entities$Length <- data_list[seq(3, length(data_list), 4)]
  entities$Organism <- data_list[seq(4, length(data_list), 4)]
  return(entities)
}

# Path for DNA details
getRNA <- partial(getConstituentEntities, path=DNA_path)
cand_structs$NucleotideEntities <- sapply(cand_structs$PageLink, getRNA)

# Some of the above don't actually contain Nucleotides with protein
# (And, based on a cursory look at the dataset, don't contain Ribosome structs)
cand_structs$hasNoEntities <- sapply(cand_structs$NucleotideEntities, is.null)
cand_structs2 <- cand_structs2 %>% 
  filter(hasNucleotideEntities != TRUE) %>% 
  select(-hasNucleotideEntities)

# Maybe we want to look at the structures that don't contain nucleotideEntities
cand_structsFail <- cand_structs %>%
  filter(hasNucleotideEntities == T) %>%
  select(-hasNucleotideEntities)

