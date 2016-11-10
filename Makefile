PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd readme check clean

alldocs: rd readme mkdocs

rd:
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd")'

build:
	cd ..;\
	R CMD build $(PKGSRC)

build2:
	cd ..;\
	R CMD build --no-build-vignettes $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

site: mkdocs

mkdocs: mdfiles
	cd mkdocs;\
	mkdocs build;\
	cd ../docs;\
	rm -rf fonts;\
	rm -rf css/font-awesome*;\
	Rscript -e 'library(ypages); add_biobabble("index.html")'


mdfiles:
	cd mkdocs;\
	Rscript -e 'library(ypages); gendoc("src/index.md", "blue", "docs/index.md")';\
	Rscript -e 'library(ypages); gendoc("src/documentation.md", "blue", "docs/documentation.md")';\
	Rscript -e 'library(ypages); gendoc("src/featuredArticles.md", "blue", "docs/featuredArticles.md")';\
	cd docs;\
	ln -f -s ../mysoftware/* ./

svnignore:
	svn propset svn:ignore -F .svnignore .

