all:	INS_GNSS
	gcc -Wall -g -w -o INS_GNSS INS_GNSS.c -lm -lpthread
