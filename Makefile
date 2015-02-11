all:
	clang++-3.5 -g -std=c++11 test/test.cpp -o test_bin
	clang++-3.5 -fPIC -shared -L/usr/lib/x86_64-linux-gnu/ -lboost_python-py27 -I /usr/include/python2.7 python/bindings.cpp -o python_bin.so


