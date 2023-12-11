#!/bin/bash
INSTALL_DIR=/usr/local

git clone https://github.com/HDFGroup/hdf5.git
cd hdf5
git checkout hdf5_1_10_11
./configure --prefix=$INSTALL_DIR
make -j 2
make check
sudo make install
sudo make check-install
cd ..

wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.1.tar.gz
tar -xzvf v3.10.1.tar.gz
cd lapack-3.10.1
mkdir build
cd build
sudo cmake -DCMAKE_INSTALL_LIBDIR=$INSTALL_DIR/lib/ ..
sudo cmake --build . -j --target install
cd ../..
