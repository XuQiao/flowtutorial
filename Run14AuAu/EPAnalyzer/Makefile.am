AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libEPAnalyzer.la

noinst_HEADERS = \
  EPAnalyzer.h \
  EPAnalyzerLinkDef.h

libEPAnalyzer_la_SOURCES = \
  EPAnalyzer.C \
  EPAnalyzer_Dict.C 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4allfuncs \
  -lreactionplane \
  -lTOAD

BUILT_SOURCES = \
  testexternals.C

testexternals_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libEPAnalyzer.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


# Rule for generating CINT dictionaries from class headers.
EPAnalyzer_Dict.C: \
  EPAnalyzer.h \
  EPAnalyzerLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^
