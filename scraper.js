// Uses phantomjs to scrape the search result page
var page = require('webpage').create();
var num_pages = 0;
page.open('http://www.rcsb.org/pdb/results/results.do?tabtoshow=Current&qrid=58258A08', function(status) {
  num_pages += 1;
  console.log("Status: " + status +"\nNumScraped: " + num_pages);
  if(status === "success") {
    // console.log(page.content);
  }
  phantom.exit();
});
