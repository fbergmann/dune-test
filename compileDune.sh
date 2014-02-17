#!/bin/bash
DIRECTORY=$(cd `dirname $0` && pwd)
cd $DIRECTORY

ALUGRID_VERSION=-DALUGRID_VERSION_CMD=$DIRECTORY/install/bin/alugridversion
if [ -f "$DIRECTORY/install/bin/alugridversion.exe" ];
then
  ALUGRID_VERSION=-DALUGRID_VERSION_CMD=$DIRECTORY/install/bin/alugridversion.exe
fi


CMAKE_OPTS="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install "

# build alberta if it does not exist
if [ ! -d "$DIRECTORY/alberta-2.0.1" ]; then
  if [ ! -d "$DIRECTORY/alberta-2.0.1.tar.gz" ]; then
    curl http://dl.dropboxusercontent.com/u/7222043/alberta-2.0.1.tar.gz -o alberta-2.0.1.tar.gz
  fi 
  tar zxf alberta-2.0.1.tar.gz
  chmod -R 755 alberta-2.0.1
  cd alberta-2.0.1
  ./configure --prefix=$DIRECTORY/install
  make 
  make install
  cd ..  
fi


# build ALUgrid if it does not exist
if [ ! -d "$DIRECTORY/ALUGrid-1.50" ]; then
  if [ ! -d "$DIRECTORY/ALUGrid-1.50.tar.gz " ]; then
  curl http://dl.dropboxusercontent.com/u/7222043/ALUGrid-1.50.tar.gz -o ALUGrid-1.50.tar.gz 
  fi 
  tar zxf ALUGrid-1.50.tar.gz
  chmod -R 755 ALUGrid-1.50
  cd ALUGrid-1.50
  ./configure --prefix=$DIRECTORY/install
  make 
  make install
  cd ..
fi


# build SuperLU_4.3 if it does not exist
if [ ! -d "$DIRECTORY/SuperLU_4.3" ]; then
  if [ ! -d "$DIRECTORY/superlu_4.3.tar.gz " ]; then
  curl http://dl.dropboxusercontent.com/u/7222043/superlu_4.3.tar.gz -o superlu_4.3.tar.gz
  fi 
  tar zxf superlu_4.3.tar.gz
  chmod -R 755 SuperLU_4.3
  mkdir build_superlu
  cd build_superlu
  cmake $CMAKE_OPTS ../SuperLU_4.3
  make 
  make install
  cd ..
fi



mkdir build-common
cd build-common
cmake $CMAKE_OPTS ../dune-common
make
make install
cd ..

mkdir build-geometry
cd build-geometry
cmake $CMAKE_OPTS ../dune-geometry
make
make install
cd ..


mkdir build-istl
cd build-istl
cmake $CMAKE_OPTS ../dune-istl
make
make install
cd ..

mkdir build-grid
cd build-grid
cmake $CMAKE_OPTS $ALUGRID_VERSION ../dune-grid
make
make install
cd ..

mkdir build-localfunctions
cd build-localfunctions
cmake $CMAKE_OPTS ../dune-localfunctions
make
make install
cd ..

mkdir build-typetree
cd build-typetree
cmake $CMAKE_OPTS ../dune-typetree
make
make install
cd ..

mkdir build-grid-howto
cd build-grid-howto
cmake $CMAKE_OPTS ../dune-grid-howto
make
make install
cd ..


mkdir build-pdelab
cd build-pdelab
cmake $CMAKE_OPTS $ALUGRID_VERSION ../dune-pdelab
make
make install
cd ..

mkdir build-pdelab-howto
cd build-pdelab-howto
cmake $CMAKE_OPTS $ALUGRID_VERSION ../dune-pdelab-howto
make -i
make install
cd ..
