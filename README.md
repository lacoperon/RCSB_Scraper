# Web Scraper for Manipulating the 238 Candidate Structures



(Within `scraper.R`)

Firstly, I downloaded the HTML associated with this query: http://www.rcsb.org/pdb/results/results.do?tabtoshow=Current&qrid=58258A08&startAt=0&resultsperpage=238, to get all candidate PDB structures of our desired search. The search itself was dynamically generated using Javascript, so downloading it was necessary.

Then, I loaded that HTML into R, and extracted some relevant information embedded in the HTML, as well as the links to the individual structure sites returned by the search.

Finally, I did some web scraping on the individual structure sites, to get more information.

## Discussion of Measures Collected

I collected information on the PDB Number of the structure, a basic description of the structure, a link to the PDB structure file, the year + journal the paper associated with the structure was published in, the link to the page of the structure candidate, the FASTA link of the candidate structure, the authors of the paper associated with the structures, the method of obtaining the structure, and the resolution of the structures (in Angstroms).

## Graphing These Measures

Code within `graphing.R`, figures within `./figures`.