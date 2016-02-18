#!/bin/sh

function run() {
    echo "$@"
    "$@" || exit 1
}

if [ -f Makefile ];then
    make distclean
fi

case `uname` in
    "Darwin" ) export DIST_SUFFIX="OSX" ;;
    "Linux" )
        if [ -f /etc/redhat-release ];then
            RHEL_VERSION=`sed -e 's/^[^[:digit:]]*\([[:digit:]]\).*$/\1/g' /etc/redhat-release`
            DIST_SUFFIX="EL${RHEL_VERSION}"
        else
            DIST_SUFFIX="linux"
        fi
        export CPPFLAGS="-static"
        ;;
    * )
        DIST_SUFFIX="unknown";;
esac

run ./autogen.sh
run ./configure --disable-shared
run make dist-bin

