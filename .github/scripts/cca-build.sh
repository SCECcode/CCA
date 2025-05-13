#!/bin/bash

tmp=`uname -s`

if [ $tmp == 'Darwin' ]; then
##for macOS, make sure have automake/aclocal
  brew install automake
  brew reinstall gcc
fi

aclocal -I m4
autoreconf -i -f
automake --add-missing --force-missing
./configure --prefix=$UCVM_INSTALL_PATH/model/cca --with-etree-libdir=${UCVM_INSTALL_PATH}/lib/euclid3/lib --with-etree-incdir=${UCVM_INSTALL_PATH}/lib/euclid3/include --with-proj-libdir=${UCVM_INSTALL_PATH}/lib/proj/lib --with-proj-incdir=${UCVM_INSTALL_PATH}/lib/proj/include --with-tiff-libdir=${UCVM_INSTALL_PATH}/lib/tiff/lib --with-tiff-incdir=${UCVM_INSTALL_PATH}/lib/tiff/include --with-sqlite-libdir=${UCVM_INSTALL_PATH}/lib/sqlite/lib --with-sqlite-incdir=${UCVM_INSTALL_PATH}/lib/sqlite/include
make
make install

