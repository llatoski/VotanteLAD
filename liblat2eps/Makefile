CC = gcc
LD = gcc
AR = ar

CFLAGS = -Wall -O3 -I.
LFLAGS = -L. 
LIBS   = -llat2eps -lm

all: liblat2eps.a lat2eps test test2

clean:
	rm -f *.o

cleaner:
	rm -f liblat2eps.a lat2eps test test2 *.o

liblat2eps.a:	lat2eps_lib.o
	$(AR) rcs liblat2eps.a lat2eps_lib.o

lat2eps:	lat2eps_cmd.o liblat2eps.a
	$(LD) $(LFLAGS) lat2eps_cmd.o -o lat2eps $(LIBS)

test:	test.o liblat2eps.a
	$(LD) $(LFLAGS) test.o -o test $(LIBS)

test2:	test2.o
	$(LD) test2.o -o test2

lat2eps_lib.o:	lat2eps_lib.c lat2eps.h
	$(CC) $(CFLAGS) -c lat2eps_lib.c -o lat2eps_lib.o

lat2eps_cmd.o:	lat2eps_cmd.c lat2eps.h
	$(CC) $(CFLAGS) -c lat2eps_cmd.c -o lat2eps_cmd.o

test.o:	test.c lat2eps.h
	$(CC) $(CFLAGS) -c test.c -o test.o

test2.o:	test2.c
	$(CC) $(CFLAGS) -c test2.c -o test2.o

