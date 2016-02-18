AM_LDFLAGS = -no-undefined
AM_CPPFLAGS = -Wall -Wextra

noinst_LTLIBRARIES = libkyotocabinet.la
libkyotocabinet_la_SOURCES = kcutil.cc kcthread.cc kcfile.cc kccompress.cc kccompare.cc kcmap.cc kcregex.cc kcdb.cc kcplantdb.cc kcprotodb.cc kcstashdb.cc kccachedb.cc kchashdb.cc kcdirdb.cc kctextdb.cc kcpolydb.cc kcdbext.cc kclangc.cc
libkyotocabinet_la_CFLAGS = -ansi -pedantic -fsigned-char
libkyotocabinet_la_CXXFLAGS = -fsigned-char
libkyotocabinet_la_CPPFLAGS = $(AM_CPPFLAGS) -DNDEBUG -D_GNU_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_REENTRANT -D__EXTENSIONS__ -D_MYZLIB -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-conversion-null

noinst_HEADERS = cmdcommon.h kccachedb.h kccommon.h kccompare.h kccompress.h kcdb.h kcdbext.h kcdirdb.h kcfile.h kchashdb.h kclangc.h kcmap.h kcplantdb.h kcpolydb.h kcprotodb.h kcregex.h kcstashdb.h kctextdb.h kcthread.h kcutil.h myconf.h

EXTRA_DIST = COPYING README ChangeLog FOSSEXCEPTION LINKEXCEPTION
