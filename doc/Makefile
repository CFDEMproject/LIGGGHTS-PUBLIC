# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS     ?=
SPHINXBUILD    ?= sphinx-build
PAPER          ?=
GIT            ?= git
TXT2RST        ?= txt2rst
BUILDDIR       ?= _build
TMPDIR         ?= _tmp
SPHINXTHEMEDIR ?= $(PWD)/$(TMPDIR)/sphinx_rtd_theme
# Downloads for external utilities
TXT2RST_URL     = https://github.com/CFDEMproject/lammps-doc-utils.git
SPHINXTHEME_URL = https://github.com/CFDEMproject/sphinx_rtd_theme.git

# User-friendly check for sphinx-build
ifneq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 0)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don`t have Sphinx installed, grab it from http://sphinx-doc.org/. If you use Ubuntu or Debian you might be able to install it using 'sudo apt-get install python-sphinx')
endif

# Check for installed git
GITINSTALLED = 1
ifneq ($(shell which $(GIT) >/dev/null 2>&1; echo $$?), 0)
GITINSTALLED = 0
endif

UTILITY_CHECK = 1
ifeq ($(MAKECMDGOALS), )
UTILITY_CHECK = 0
else
ifneq (,$(filter $(MAKECMDGOALS), help clean clean-all))
UTILITY_CHECK = 0
endif
endif

ifeq ($(UTILITY_CHECK), 1)
# Find txt2rst or install it if possible
ifneq ($(shell which $(TXT2RST) >/dev/null 2>&1; echo $$?), 0)
# Check if it exists in default temporary directory
ifneq ("$(wildcard $(TMPDIR)/lammps-doc-utils/bin/txt2rst)", "")
TXT2RST=$(TMPDIR)/lammps-doc-utils/bin/txt2rst
$(warning Found txt2rst in '$(TMPDIR)/lammps-doc-utils/bin/txt2rst' and using this version. Use 'make clean-all' to remove it)
else
# If git is not available we can't do anything
ifeq ($(GITINSTALLED), 0)
$(error Could not find '$(TXT2RST)' and '$(GIT)' is not present either. Please install either the txt2rst utility from $(TXT2RST_URL), change the TXT2RST environment variable or add the directory containing your txt2rst to the PATH variable. Alternatively, you can install '$(GIT)' and run this makefile again to install txt2rst automatically. If '$(GIT)' is already installed ensure that it points to the correct path or add the directory containing your git executable to the PATH variable.)
# else we can download it
else
# create temporary folder
ifneq ($(shell mkdir -p $(TMPDIR) >/dev/null 2>&1; echo $$?), 0)
$(error Could not create $(TMPDIR))
endif
# download git repository
ifneq ($(shell git clone $(TXT2RST_URL) $(TMPDIR)/lammps-doc-utils >/dev/null 2>&1; echo $$?), 0)
$(error Could not download lammps-doc-utils that provides txt2rst. Please install the utility manually)
endif
# check if file is available now
TXT2RST=$(TMPDIR)/lammps-doc-utils/bin/txt2rst
ifneq ($(shell which $(TXT2RST) >/dev/null 2>&1; echo $$?), 0)
$(error Unknown error in txt2rst installation. Please install the utility manually)
endif
$(warning Successfully downloaded txt2rst and installed it in the $(TMPDIR) directory)
endif
# end of download (git available)
endif
# end of not in default dir
endif
# end of not in specified dir

# Check if the sphinx theme is available
ifeq ("$(wildcard $(SPHINXTHEMEDIR)/sphinx_rtd_theme/__init__.py)", "")
# Set to default directory
SPHINXTHEMEDIR=$(PWD)/$(TMPDIR)/sphinx_rtd_theme
ifneq ("$(wildcard $(SPHINXTHEMEDIR)/sphinx_rtd_theme/__init__.py)", "")
$(warning Found sphinx_rtd_theme in default $(TMPDIR) directory and using this version. Use 'make clean-all' to remove it)
else
# If git is not available we can't do anything
ifeq ($(GITINSTALLED), 0)
$(error Could not find sphinx_rtd_theme and '$(GIT)' is not present either. Please install either sphinx_rtd_theme from $(SPHINXTHEME_URL) or change the SPHINXRTDTHEMEDIR environment variable. Alternatively, you can install '$(GIT)' and run this makefile again to install sphinx_rtd_theme automatically. If '$(GIT)' is already installed ensure that it points to the correct path or add the directory containing your git executable to the PATH variable.)
# else we can download it
else
# create temporary folder
ifneq ($(shell mkdir -p $(TMPDIR) >/dev/null 2>&1; echo $$?), 0)
$(error Could not create $(TMPDIR))
endif
# download git repository
ifneq ($(shell git clone $(SPHINXTHEME_URL) $(TMPDIR)/sphinx_rtd_theme >/dev/null 2>&1; echo $$?), 0)
$(error Could not download sphinx_rtd_theme. Please install the theme manually from $(SPHINXTHEME_URL))
endif
# check if file is available now
#ifeq ("$(wildcard $(SPHINXTHEMEDIR)/sphinx_rtd_theme/__init__.py)", "")
#$(error Unknown error in the sphinx_rtd_theme installation. Please install the theme manually from $(SPHINXTHEME_URL))
#endif
$(warning Successfully downloaded the sphinx_rtd_theme and installed it in the $(TMPDIR) directory)
endif
# end of download (git available)
endif
# end of not in default dir
endif
# end of not in specified dir
endif
# end of UTILITY_CHECK

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean clean-all html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext generate-rst

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  clean      to remove all generated files"
	@echo "  clean-all  to remove all generated and temporary files"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  latexpdfja to make LaTeX files and run them through platex/dvipdfmx"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  xml        to make Docutils-native XML files"
	@echo "  pseudoxml  to make pseudoxml-XML files for display purposes"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"

clean:
	rm -rf $(BUILDDIR)/* *.rst

clean-all: clean
	rm -rf $(TMPDIR)/lammps-doc-utils $(TMPDIR)/sphinx_rtd_theme _themes/sphinx_rtd_theme

generate-rst:
	@mkdir -p _themes
	@rm -rf _themes/sphinx_rtd_theme
	@ln -s $(SPHINXTHEMEDIR)/sphinx_rtd_theme _themes/
	$(TXT2RST) *txt

html: generate-rst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml: generate-rst
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml: generate-rst
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle: generate-rst
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json: generate-rst
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp: generate-rst
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp: generate-rst
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(BUILDDIR)/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/LIGGGHTS.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/LIGGGHTS.qhc"

devhelp: generate-rst
	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(BUILDDIR)/devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/LIGGGHTS"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/LIGGGHTS"
	@echo "# devhelp"

epub: generate-rst
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

latex: generate-rst
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf: generate-rst
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

latexpdfja: generate-rst
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through platex and dvipdfmx..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf-ja
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

text: generate-rst
	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man: generate-rst
	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(BUILDDIR)/man
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo: generate-rst
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info: generate-rst
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext: generate-rst
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes: generate-rst
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck: generate-rst
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest: generate-rst
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

xml: generate-rst
	$(SPHINXBUILD) -b xml $(ALLSPHINXOPTS) $(BUILDDIR)/xml
	@echo
	@echo "Build finished. The XML files are in $(BUILDDIR)/xml."

pseudoxml: generate-rst
	$(SPHINXBUILD) -b pseudoxml $(ALLSPHINXOPTS) $(BUILDDIR)/pseudoxml
	@echo
	@echo "Build finished. The pseudo-XML files are in $(BUILDDIR)/pseudoxml."
