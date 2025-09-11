# Makefile for building both documentation sets

.PHONY: help clean html user programmer all

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files for both manuals"
	@echo "  user       to make standalone HTML files for user manual only"
	@echo "  programmer to make standalone HTML files for programmer manual only"
	@echo "  clean      to clean all build directories"

clean:
	rm -rf usermanuals/MAROONXDR_UserManual/build/*
	rm -rf progmanuals/MAROONXDR_ProgManual/build/*

user:
	$(MAKE) -C usermanuals/MAROONXDR_UserManual html

programmer:
	$(MAKE) -C progmanuals/MAROONXDR_ProgManual html

tutorial:
	$(MAKE) -C tutorials/MAROONXDR_Tutorial html


html: user programmer tutorial

all: html