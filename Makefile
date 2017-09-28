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

site:
	cd site_src;\
	ln -s ../../software/themes themes;\
	Rscript -e 'blogdown::build_site()';\
	rm themes;\
	cd ..

preview:
	cd site_src;\
	ln -s ../../software/themes themes;\
	Rscript -e 'blogdown::serve_site()';\
	rm themes;\
	cd ..

update:
	git fetch --all;\
	git checkout master;\
	git merge upstream/master;\
	git merge origin/master


push: update
	git push upstream master;\
	git push origin master


