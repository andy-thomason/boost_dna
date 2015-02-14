all:
	clang++ -std=c++11 -fPIC -shared -L/usr/lib/x86_64-linux-gnu/ -lboost_python-py27 -I /usr/include/python2.7 python/bindings.cpp -o suffix_array.so

win:
	"C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\BIN\x86_amd64\cl.exe" /MD /LD -I c:/Python27/include -I C:/tmp/boost/boost_1_57_0 python/bindings.cpp /Fe:python_bin.dll

