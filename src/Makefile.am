AM_LDFLAGS = -no-undefined
AM_CPPFLAGS = -Wall -Wextra -I$(top_builddir)/libs/CppFasta/src/ -I$(top_builddir)/libs/kyotocabinet/

bin_PROGRAMS = matataki matataki-builddb matataki-dbsummary refseqextract
matataki_SOURCES = quickexpression.cpp
matataki_LDADD = libmatataki.la

matataki_builddb_SOURCES = quickexpression-builddb.cpp
matataki_builddb_LDADD = libmatataki.la

matataki_dbsummary_SOURCES = quickexpression-dbsummary.cpp
matataki_dbsummary_LDADD = libmatataki.la 

refseqextract_SOURCES = refseqextract.cpp
refseqextract_LDADD = libmatataki.la

noinst_LTLIBRARIES = libmatataki.la
libmatataki_la_SOURCES = dbbuilder.cpp gene2refseq.cpp tabtable.cpp quickcommon.cpp mapping.cpp dbcommon.cpp dbsummary.cpp dbcoverage.cpp genomehash.cpp MurmurHash3.cpp spawn.c pipedopen.c
libmatataki_la_CPPFLAGS = $(AM_CPPFLAGS)
libmatataki_la_LIBADD = $(top_builddir)/libs/CppFasta/src/libcppfasta.la $(top_builddir)/libs/kyotocabinet/libkyotocabinet.la

noinst_HEADERS = MurmurHash3.h config.h dbbuilder.hpp dbcommon.hpp dbcoverage.hpp dbsummary.hpp gene2refseq.hpp genomehash.hpp mapping.hpp pipedopen.h quickcommon.hpp spawn.h tabtable.hpp
