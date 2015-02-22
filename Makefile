all:
	clang++ -O2 -std=c++11 -fPIC -shared -L/usr/lib/x86_64-linux-gnu/ -lboost_python-py27 -I /usr/include/python2.7 python/bindings.cpp -o suffix_array.so
