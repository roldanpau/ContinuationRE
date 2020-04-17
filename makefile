SHELL = /bin/sh
prefix = /home/pau
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -g
LDFLAGS = -g
LDLIBS = -lm -lgsl -lgslcblas

all : Lagr

install : Lagr
#	ar rv $(libdir)/libds.a ctbp.o
#	cp ctbp.h $(includedir)
	cp ctbp $(bindir)

Lagr : Lagr.o

clean : 
	rm Lagr Lagr.o
