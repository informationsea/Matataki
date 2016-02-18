SUBDIRS = libs src test
EXTRA_DIST = README.md build-binary-release.sh autogen.sh

DIST_SUFFIX ?= unknown

dist-bin: all check
	$(am__remove_distdir)
	test -d "$(distdir)" || mkdir "$(distdir)"
	cp ./src/matataki "$(distdir)"
	cp ./src/matataki-builddb "$(distdir)"
	cp ./src/matataki-dbsummary "$(distdir)"
	cp ./src/refseqextract "$(distdir)"
	cp README.md "$(distdir)"
	cp COPYING "$(distdir)"
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir)-bin-$(DIST_SUFFIX).tar.gz
	$(am__post_remove_distdir)

