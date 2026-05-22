# Makefile for building MaroonX DRAGONS unified documentation

.PHONY: help clean html unified user programmer tutorial all

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make unified HTML documentation (recommended)"
	@echo "  unified    to make unified HTML documentation (same as html)"
	@echo "  user       to make standalone HTML files for user manual only"
	@echo "  programmer to make standalone HTML files for programmer manual only"
	@echo "  tutorial   to make standalone HTML files for tutorial only"
	@echo "  clean      to clean all build directories"

clean:
	rm -rf _build/*
	rm -rf usermanuals/MAROONXDR_UserManual/build/*
	rm -rf progmanuals/MAROONXDR_ProgManual/build/*
	rm -rf tutorials/MAROONXDR_Tutorial/build/*

# Build unified documentation (all three manuals in one)
unified:
	sphinx-build -b html . _build/html

# Default target builds unified docs
html: unified

# Individual manual builds (legacy support)
user:
	$(MAKE) -C usermanuals/MAROONXDR_UserManual html

programmer:
	$(MAKE) -C progmanuals/MAROONXDR_ProgManual html

tutorial:
	$(MAKE) -C tutorials/MAROONXDR_Tutorial html

all: html