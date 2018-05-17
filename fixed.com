ulimit -s unlimited
gcc fixed.c -o fixed -C -lgsl -L/usr/lib -lgslcblas  -lm
#LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
./fixed -i 50 -c 5.0
