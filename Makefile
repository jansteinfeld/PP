# Makefile for generating R packages.
# 2011 Andrew Redd
#
# Assumes Makefile is in a folder where package contents are in a subfolder PP.
# Roxygen uses the roxygen2 package, and will run automatically on check and all.

## run it like:
# > cd .. && make 

PKG_VERSION=$(shell grep -i ^version PP/DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package PP/DESCRIPTION | cut -d : -d \  -f 2)

R_FILES := $(wildcard PP/R/*.R)
SRC_FILES := $(wildcard PP/src/*) $(addprefix PP/src/, $(COPY_SRC))
PKG_FILES := PP/DESCRIPTION PP/NAMESPACE $(R_FILES) $(SRC_FILES)

.PHONY: tarball install check clean build

tarball: $(PKG_NAME)_$(PKG_VERSION).tar.gz
$(PKG_NAME)_$(PKG_VERSION).tar.gz: $(PKG_FILES)
	R CMD build PP


check: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD check --as-cran $(PKG_NAME)_$(PKG_VERSION).tar.gz

build: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL --build $(PKG_NAME)_$(PKG_VERSION).tar.gz

install: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz

PP/NAMESPACE: $(R_FILES)
	Rscript -e "library(roxygen2);roxygenize('PP')"

clean:
	-rm -f $(PKG_NAME)_*.tar.gz
	-rm -r -f $(PKG_NAME).Rcheck
	-rm -r -f PP/man/*
	-rm -r -f PP/NAMESPACE
.PHONY: list
list:
	@echo "R files:"
	@echo $(R_FILES)
	@echo "Source files:"
	@echo $(SRC_FILES)