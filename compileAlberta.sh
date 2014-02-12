DIRECTORY=$(cd `dirname $0` && pwd)
cd $DIRECTORY
cd alberta-2.0.1
./configure --prefix=$DIRECTORY/install
make
make install
cd ..
