curl -O https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar -xf gmp-6.3.0.tar.xz
mv gmp-6.3.0/* .
./configure enable-cxx
make -j16
export LD_LIBRARY_PATH="./.libs:$LD_LIBRARY_PATH"
g++ -o exe -I./ -L./.libs main.cpp -O3 -march=native -lgmpxx -lgmp
./exe 10000000
