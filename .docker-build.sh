#! /bin/bash
set -e
yum install -y gmp-devel gcc gcc-c++ || ( apt-get update && apt-get -y install g++ libgmp-dev make)
mkdir -p /build
cd /build
/src/configure --srcdir=/src  || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
make -j2 && CAT_CHECKDIR_ON_ERROR=t make check
