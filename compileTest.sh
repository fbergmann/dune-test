#!/bin/bash
DIRECTORY=$(cd `dirname $0` && pwd)
cd $DIRECTORY

ALUGRID_VERSION=-DALUGRID_VERSION_CMD=$DIRECTORY/install/bin/alugridversion
if [ -f "$DIRECTORY/install/bin/alugridversion.exe" ];
then
  ALUGRID_VERSION=-DALUGRID_VERSION_CMD=$DIRECTORY/install/bin/alugridversion.exe
fi

CMAKE_OPTS="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install "

mkdir build-test
cd build-test
cmake $CMAKE_OPTS $ALUGRID_VERSION ../dune-test
make -i
cd ..
