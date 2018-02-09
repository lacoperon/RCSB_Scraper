library(XML)
library(readr)

# Opens pre-saved version of the search results
searchResults <- read_file("data/searchResults.html")
d <- htmlParse(searchResults)

# Gets the number of structures returned
node_num_struct <- getNodeSet(doc=d, path="//ul/li[@class='active']//a")[[1]]
num_structure <- trimws(xmlToList(node_num_struct)$text)

# Scrapes all search result HTML nodes (and their PDB accession numbers)
node_set_path <- "//ul[@id='SearchResultsDetails-MainContent']/li"
result_nodes <- getNodeSet(doc=d, path=node_set_path)
cand_structs  <- data.frame(sapply(result_nodes, function(node) {substr(xmlAttrs(node)[[1]], 11, 1000)}))
colnames(cand_structs) <- c("PDB_Number")

# Scrapes all search result links (to the PDB structures themselves)
file_link_path <- "//ul[@id='SearchResultsDetails-MainContent']/li/div[2]/div[1]/div/a[1]"
cand_structs$Link   <- xpathSApply(doc=d, path=file_link_path, fun=xmlGetAttr, name="href")
