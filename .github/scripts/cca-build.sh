#!/bin/bash

tmp=`uname -s`

if [ $tmp == 'Darwin' ]; then
##for macOS, make sure have automake/aclocal
  brew install automake
  brew reinstall gcc
fi

aclocal -I m4
automake --add-missing --force-missing
autoconf
./configure --prefix=$UCVM_INSTALL_PATH/model/cca
make
make install

