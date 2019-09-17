# Bellow the code is being linked to rtklib.h library
# In case need to use other source files, they need to be
# compiled and linked as well, in this case by mentioning
# it in the line of compilation.
# See rtklib manual, page 89 for details.

# Example:
# Suppose that your program has two source files: file1.cpp,
# file2.cpp. You could compile all of them in a single command:

# > g++ -o myprog.exe file1.cpp file2.cpp (C++ in this case)

# However, we usually compile each of the source files separately
# into object file, and link them together in the later stage.
# In this case, changes in one file does not require re-compilation
# of the other files.

SRC     = ../lib/RTKLIB/src
SRC1     = ../lib/gnssins
LIB	= ../lib

all:	satinsmap

satinsmap:	satinsmap.c
	gcc -Wall -g -w -o satinsmap satinsmap.c plots.c $(SRC1)/*.c $(SRC)/*.c $(SRC)/rcv/*.c -I$(SRC) -DENAGLO -DLAPACK -llapack -lblas -lm -lpthread
