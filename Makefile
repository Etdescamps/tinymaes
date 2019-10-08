CC = gcc
AR = ar rcs
#CFLAGS = -O2 -funroll-loops -Wall
CFLAGS = -ggdb -fsanitize=address -Wall
LDFLAGS = -lm
RM = rm -f
PREFIX = /usr/local
INSTALL = install -D -m644 -t
INSTALLEXEC = install -D -m755 -t

all: libtinymaes.so libtinymaes.a

install: libtinymaes.so libtinymaes.a
	$(INSTALLEXEC) $(PREFIX)/lib libtinymaes.so
	$(INSTALL) $(PREFIX)/lib libtinymaes.a
	$(INSTALL) $(PREFIX)/include/tinymaes/ mt64.h tinymaes.h heapsort.h

libtinymaes.so: mt64.c tinymaes.c heapsort.c
	$(CC) $(CFLAGS) $(LDFLAGS) -fPIC -shared mt64.c tinymaes.c heapsort.c -o libtinymaes.so

libtinymaes.a: mt64.o tinymaes.o heapsort.o
	$(AR) libtinymaes.a mt64.o tinymaes.o heapsort.o

clean:
	$(RM) mt64.o tinymaes.o heapsort.o libtinymaes.so libtinymaes.a

.c.o:
	$(CC) $(CFLAGS) -c $<



