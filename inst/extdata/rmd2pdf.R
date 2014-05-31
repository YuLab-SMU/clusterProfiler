require(knitr)
require(markdown)

knit('for_unsupported.rmd', 
	 'for_unsupported.md')

markdownToHTML('for_unsupported.md',
			   'for_unsupported.html', 
			   options=c('use_xhml'))

system("wkhtmltopdf for_unsupported.html clusterProfiler_for_unsupported_organisms.pdf")
