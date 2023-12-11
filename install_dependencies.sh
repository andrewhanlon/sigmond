#!/bin/bash
INSTALL_DIR=/usr/local

git clone https://github.com/HDFGroup/hdf5.git
git checkout hdf5_1_10_11
cd hdf5
./configure --prefix=$INSTALL_DIR
make -j 2
make check
sudo make install
sudo make check-install
cd ..

wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.tar.gz
tar -xzvf v3.11.tar.gz
cd lapack-3.11
mkdir build
cd build
sudo cmake -DCMAKE_INSTALL_LIBDIR=($INSTALL_DIR)/lib/ ..
sudo cmake --build . -j --target install
cd ../..
