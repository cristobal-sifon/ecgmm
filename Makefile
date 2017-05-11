pyver=$(shell python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)")
py_include_path=$(shell python -c "import sys; from distutils import sysconfig; sys.stdout.write(sysconfig.get_python_inc())")

$(info Python version: $(pyver))
$(info Python include directory: $(py_include_path))

py: _ecgmm.so
	mv bin/_ecgmm.so py/ 
	mv src/ecgmm.py* py/

_ecgmm.so:ecgmm_wrap.o
	g++ -shared bin/ecgmm_wrap.o -o bin/_ecgmm.so -Ipython

ecgmm_wrap.o:ecgmm_wrap.cxx
	g++ -c -fPIC src/ecgmm_wrap.cxx -I$(py_include_path) -o bin/ecgmm_wrap.o

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
