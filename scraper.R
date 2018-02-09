library(httr)
library(XML)
library(readr)
library(purrr)

url <- "http://www.rcsb.org/pdb/results/results.do?tabtoshow=Current&qrid=58258A08"
s <- GET(url)
w <- content(s, as='text') # converts s to plaintext of HTML
d <- htmlParse(file=content(s, as="text", asText=T))

# Gets the number of structures returned
node_num_struct <- getNodeSet(doc=d, path="//ul/li[@class='active']//a")[[1]]
num_structure <- trimws(xmlToList(node_num_struct)$text)

# Scrapes all search result HTML nodes
result_nodes <- getNodeSet(doc=d, path="//ul[@class='list-unstyled SearchResultsBlock']")
