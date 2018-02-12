library(XML)
library(readr)
library(httr)
library(dplyr)
library(purrr)
library(stringr)

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

# Path associated with 
DNA_path <- '//*[@id="MacromoleculeTableDNA"]/div/table/tbody/tr[3]'
# Function to get the constituent Macromolecules in the structure
# (In this case, I'll only use it for the RNA Parts, but it can also
# be applied to Protein as well)
getConstituentEntities <- function(url, path, mode="Nucleo") {
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
  if (mode == "Nucleo") {
    EntityName <- data_list[seq(1, length(data_list), 4)]
    entities <- data.frame(EntityName)
    entities$ChainNames <- data_list[seq(2, length(data_list), 4)]
    entities$Length <- data_list[seq(3, length(data_list), 4)]
    entities$Organism <- data_list[seq(4, length(data_list), 4)]
  }
  
  if (mode == "Protein") {
    EntityName <- data_list[seq(1, length(data_list), 5)]
    entities <- data.frame(EntityName)
    entities$ChainNames <- data_list[seq(2, length(data_list), 5)]
    entities$Length <- data_list[seq(3, length(data_list), 5)]
    entities$Organism <- data_list[seq(4, length(data_list), 5)]
    entities$Details <- data_list[seq(5, length(data_list), 5)]
    entities$EntityName <- as.character(entities$EntityName)
  }
  
  return(entities)
}

# Path for RNA details
getRNA <- partial(getConstituentEntities, path=DNA_path)
NucleotideEntities <- sapply(cand_structs$PageLink, getRNA)

# Pointer to original index
cand_structs$OriginalIndex <- seq(1, dim(cand_structs)[1])


# Some of the above don't actually contain Nucleotides with protein
# (And, based on a cursory look at the dataset, don't contain Ribosome structs)
cand_structs$noNucleotideEntities <- sapply(NucleotideEntities, is.null)

cand_structs2 <- cand_structs %>% 
  filter(noNucleotideEntities != T) %>% 
  select(-noNucleotideEntities)

# Maybe we want to look at the structures that don't contain nucleotideEntities
# Went through by hand -- they actually don't have nucleotide elements,
# and often are irrelevant
cand_structsFailNucleotideEntities <- cand_structs %>%
  filter(noNucleotideEntities == T) %>%
  select(-noNucleotideEntities)

# Function that returns whether or not a particular structure has mRNA entity
ismRNA <- function(nucleotideDF) {
  containsmRNA <- grepl("[mM]RNA", nucleotideDF$EntityName)
  return(sum(containsmRNA) > 0)
}

# Vector of whether or not the structure contains mRNA
cand_structs$hasMRNA <- sapply(NucleotideEntities, ismRNA)

# cand_structs2 contains all structures which have mRNA
cand_structs2 <- cand_structs %>% 
  filter(hasMRNA == TRUE, noNucleotideEntities != TRUE) %>% 
  select(-noNucleotideEntities, hasMRNA)

# This contains all structures sthat don't contain mRNA. 
# Went through by hand -- they actually don't have mRNA
cand_structsNoMRNA <- cand_structs %>%
  filter(hasMRNA == FALSE, noNucleotideEntities != TRUE) %>%
  select(-noNucleotideEntities, hasMRNA)


# Path for Protein details
prot_path <- '//*[@id="MacromoleculeTable"]/div/table/tbody/tr[3]'
getProtein <- partial(getConstituentEntities, path=prot_path, mode="Protein")
ProteinEntities <- sapply(cand_structs$PageLink, getProtein)
ProteinEntities_transp <- t(ProteinEntities)

prot_seq <- seq(1, dim(ProteinEntities_transp)[1])
containsS12 <- function(index) {
  temp <- ProteinEntities_transp[index,]
  contS12 <- (sum(grepl("[sS]12", temp$EntityName)) > 0)
  return(contS12)
}
cand_structs$containsS12 <- sapply(prot_seq, containsS12)

# cand_structs2 contains all structures which have mRNA, and S12
cand_structs2 <- cand_structs %>% 
  filter(hasMRNA == TRUE, 
         noNucleotideEntities != TRUE,
         containsS12 == TRUE) %>% 
  select(-noNucleotideEntities, -containsS12, -hasMRNA)

# These don't seem to be structures we're interested in, either
cand_structsNoS12 <- cand_structs %>% 
  filter(containsS12 != TRUE) %>% 
  select(-containsS12)

# Now, to try to isolate the mRNA information

# This function returns the chain names associated with the mRNA in its
# FASTA file
get_mRNA_Chains <- function(index, filtered_ne, urls) {
  
  url <- urls[index]
  NucleoEntity <- data.frame(filtered_ne[index])
  print(url)
  s <- GET(url)
  w <- content(s, as='text') # converts FASTA response to plaintext
  fasta_list <- strsplit(w, ">")[[1]]
  df_fasta <- data.frame(fasta_list)
  df_fasta <- filter(df_fasta, fasta_list != "")
  df_fasta$chainName <- sub(".{4}:", "", df_fasta$fasta_list, perl=T)
  df_fasta$chainName <- trimws(sub("\\|.*+([[:space:]].*)*", 
                                   "", 
                                   df_fasta$chainName, perl=T))
  df_fasta$sequence  <- sub("(.*SEQUENCE)[[:space:]]", "", df_fasta$fasta_list)
  
  
  colnames(NucleoEntity) <- c("EntityName","ChainNames", "Length", "Organism")
  NucleoEntity$ismRNA <- sapply(NucleoEntity$EntityName, 
                                partial(grepl, pattern="[mM]RNA"))
  mRNAChainNames <- filter(NucleoEntity, ismRNA)$ChainNames
  return(mRNAChainNames)
}

indices <- seq(1, dim(cand_structs2)[1])
filtered_ne <- NucleotideEntities[cand_structs2$OriginalIndex]
cand_structs2$mRNA_ChainNames <- sapply(indices, partial(get_mRNA_Chains, 
                                          urls=cand_structs2$FastaLink,
                                          filtered_ne = filtered_ne))
cand_structs2$mRNA_ChainNames <- sapply(cand_structsw$mRNA_ChainNames, 
                                        function(elt) return(elt[[1]]))


# Gets number of mRNA chains
cand_structs2 <- cand_structs2 %>%
  mutate(num_mRNA_Chains=str_count(mRNA_ChainNames, ",") + 1)

# Gets max number of chains
max(cand_structs2$numChains)
# Gets min number of chains (should be 1)
min(cand_structs2$numChains)

get_mRNA_Seqs <- function(index, filtered_nu=filtered_ne, urls=cand_structs2$FastaLink) {
  url <- urls[index]
  NucleoEntity <- data.frame(filtered_nu[index])
  s <- GET(url)
  print(url)
  w <- content(s, as='text') # converts FASTA response to plaintext
  fasta_list <- strsplit(w, ">")[[1]]
  df_fasta <- data.frame(fasta_list)
  df_fasta <- filter(df_fasta, fasta_list != "")
  df_fasta$chainName <- sub(".{4}:", "", df_fasta$fasta_list, perl=T)
  df_fasta$chainName <- trimws(sub("\\|.*+([[:space:]].*)*", 
                                   "", 
                                   df_fasta$chainName, perl=T))
  df_fasta$sequence  <- sub("(.*SEQUENCE)[[:space:]]", "", df_fasta$fasta_list)
  
  colnames(NucleoEntity) <- c("EntityName","ChainNames", "Length", "Organism")
  NucleoEntity$ismRNA <- sapply(NucleoEntity$EntityName, partial(grepl, pattern="[mM]RNA"))
  mRNAChainNames <- filter(NucleoEntity, ismRNA)$ChainNames
  chainNames <- strsplit(mRNAChainNames, ",")[[1]]
  mrna_fasta <- filter(df_fasta, chainName %in% chainNames)
  mrna_fasta$sequence <- trimws(mrna_fasta$sequence)
  return(mrna_fasta$sequence)
}

cand_structs2$mRNA_Sequences <- sapply(indices, get_mRNA_Seqs)
cand_structs2$mRNA_Length <- sapply(cand_structs2$mRNA_Sequences, 
                                    function(seqs) (sum(nchar(seqs[[1]]))))

# Writing out logic
write_csv(cand_structs,  path="./data/candStructs.csv")
cand_structs2w <- cand_structs2

cand_structs2w$mRNA_Sequences <- sapply(cand_structs2w$mRNA_Sequences, flattenSeqList)
cand_structs2w$mRNA_ChainNames <- sapply(cand_structs2w$mRNA_ChainNames, function(elt) return(elt[[1]]))
write_csv(cand_structs2w,  path="./data/candStructs_filtered.csv")




