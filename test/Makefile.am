AUTOMAKE_OPTIONS = subdir-objects
SUBDIRS = gtest data
AM_CPPFLAGS = -Wall -Wextra

noinst_PROGRAMS = test_builder test_tabtable test_common test_mapping test_coverage test_genomehash test_pipedopen
TESTS = $(noinst_PROGRAMS)

LDADD = -lpthread $(top_builddir)/src/libmatataki.la $(top_builddir)/test/gtest/libgtest.la 

test_builder_SOURCES = test_builder.cpp
test_builder_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)/libs/CppFasta/src $(AM_CPPFLAGS)

test_common_SOURCES = test_common.cpp

test_tabtable_SOURCES = test_tabtable.cpp
test_tabtable_CPPFLAGS = -I$(top_builddir)/libs/CppFasta/src $(AM_CPPFLAGS)

test_mapping_SOURCES = test_mapping.cpp
test_mapping_CPPFLAGS =  $(AM_CPPFLAGS) -I$(top_builddir)/libs/CppFasta/src $(AM_CPPFLAGS)

test_coverage_SOURCES = test_coverage.cpp
test_coverage_CPPFLAGS =  $(AM_CPPFLAGS) -I$(top_builddir)/libs/CppFasta/src $(AM_CPPFLAGS)

test_genomehash_SOURCES = test_genomehash.cpp

test_pipedopen_SOURCES = test_pipedopen.cpp
