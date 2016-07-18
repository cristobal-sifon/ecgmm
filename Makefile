py: _ecgmm.so
	mv bin/_ecgmm.so py/ 
	mv src/ecgmm.py* py/

_ecgmm.so:ecgmm_wrap.o
	g++ -shared bin/ecgmm_wrap.o -o bin/_ecgmm.so -Ipython

ecgmm_wrap.o:ecgmm_wrap.cxx
	g++ -c -fPIC src/ecgmm_wrap.cxx -I/usr/include/python2.7 -o bin/ecgmm_wrap.o

ecgmm_wrap.cxx:
    mkdir -p bin
	swig -c++ -python src/ecgmm.i 

cpp: ecGMMexample.o
	g++ bin/ecGMMexample.o -Wall -lgsl -lgslcblas -lm -o bin/ecGMMexample

ecGMMexample.o:
    mkdir -p bin
	g++ -c src/ecGMMexample.cpp -o bin/ecGMMexample.o

clean:
	rm -rf bin/ecGMMexample bin/ecGMMexample.o bin/ecgmm_wrap src/ecgmm.py src/ecgmm_wrap.cxx py/ecgmm.py py/ecgmm.pyc py/_ecgmm.so bin/ecgmm_wrap.o 