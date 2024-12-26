export LD_LIBRARY_PATH="./.libs:$LD_LIBRARY_PATH"
g++ -o exe -I./ -L./.libs main.cpp -O3 -march=native -lgmpxx -lgmp -lpthread
./exe $1 $2
