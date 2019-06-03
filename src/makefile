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
	gcc -Wall -g -w -o satinsmap satinsmap.c mapmatch.c plots.c ins.c $(SRC1)/INS_GNSS.c $(SRC1)/ins_gnss_index.c $(SRC)/stream.c $(SRC)/rtkcmn.c $(SRC)/rinex.c $(SRC)/streamsvr.c $(SRC)/solution.c $(SRC)/options.c $(SRC)/geoid.c $(SRC)/postpos.c $(SRC)/rtkpos.c $(SRC)/ppp.c $(SRC)/ephemeris.c $(SRC)/sbas.c $(SRC)/qzslex.c $(SRC)/preceph.c $(SRC)/pntpos.c $(SRC)/ionex.c $(SRC)/ppp_ar.c $(SRC)/lambda.c $(SRC)/rtcm.c $(SRC)/rtcm2.c $(SRC)/rtcm3.c $(SRC)/rtcm3e.c $(SRC)/rcvraw.c $(SRC)/rcv/ublox.c $(SRC)/rcv/ss2.c $(SRC)/rcv/skytraq.c $(SRC)/rcv/rt17.c $(SRC)/rcv/rcvlex.c $(SRC)/rcv/septentrio.c $(SRC)/rcv/nvs.c $(SRC)/rcv/novatel.c $(SRC)/rcv/javad.c $(SRC)/rcv/gw10.c $(SRC)/rcv/crescent.c $(SRC)/rcv/binex.c -I$(SRC) -DENAGLO -DLAPACK -DTRACE -llapack -lblas -lm -lpthread
