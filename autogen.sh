#!/bin/sh

run()
{
    $@
    if test $? -ne 0; then
        echo "Failed $@"
        exit 1
    fi
}

case $OSTYPE in
    darwin* ) LIBTOOLIZE=glibtoolize;;
    *) LIBTOOLIZE=libtoolize;;
esac

mkdir -p m4

run aclocal ${ACLOCAL_ARGS}
run ${LIBTOOLIZE} --copy --force
run autoheader
run automake --add-missing --foreign --copy
run autoconf

(cd libs/CppFasta; ./autogen.sh)
