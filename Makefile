PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: alldocs check clean

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
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

mkdocs: featuredArticles.md index.md documentation.md
	cd mkdocs;\
	mkdocs build

index.md:
	cd mkdocs;\
	Rscript -e 'library(ypages); gendoc("private/index.md", "blue", "docs/index.md")'

documentation.md:
	cd mkdocs;\
	Rscript -e 'library(ypages); gendoc("private/documentation.md", "blue", "docs/documentation.md")'

featuredArticles.md:
	cd mkdocs;\
	Rscript -e 'library(ypages); gendoc("private/featuredArticles.md", "blue", "docs/featuredArticles.md")'

svnignore:
	svn propset svn:ignore -F .svnignore .

