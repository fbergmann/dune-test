DIRECTORY=$(cd `dirname $0` && pwd)
cd ALUGrid-1.50
./configure --prefix=$DIRECTORY/install
make
make install
cd ..
