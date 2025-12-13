GCC = g++
INC = C:\msys64\SuiteSparse\build\install\include\suitesparse
LIB = C:\msys64\SuiteSparse\build\install\lib
LPK = C:\msys64\SuiteSparse\lapack_windows\x64

FLAGS = -l cholmod -l suitesparseconfig -l umfpack

HEX = --elemType hex8
REM = --remesh

# Extra Arguments
#-l amd -l colamd -l camd -l ccolamd -l suitesparseconfig -l blas -l lapack


main: main.cpp libs/HexElem.h
	@echo "Changes found in sourcecode, compiling..."
	rm -f $(wildcard *.exe)
	${GCC} main.cpp -o main.exe -O2 -I$(INC) -L$(LIB) -L$(LPK) $(FLAGS) -Wall

clean:
	rm -f $(wildcard *.txt)

run: main clean
	.\main.exe $(HEX)

run_mesh: main clean
	.\main.exe $(HEX) $(REM)

stress: main clean
	.\main.exe $(HEX) $(REM) --stress

transient: main clean
	.\main.exe $(HEX) $(REM) --transient $(T_STOP) $(DT)