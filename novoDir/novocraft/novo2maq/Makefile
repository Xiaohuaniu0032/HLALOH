CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 -m64 # comment out `-m64' for 32-bit compilation
CXXFLAGS=	$(CFLAGS)
#DFLAGS=		-D_FASTMAP	
DFLAGS=		-D_FASTMAP -DMAQ_LONGREADS
OBJS=		const.o \
			maqmap.o \
			novo2maq.o 
PROG=		novo2maq
VERSION=	MAQV0.6.8/NovoV2
LIBS=		-lz -lm

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $< -o $@

all:$(PROG) $(MANPAGE)

novo2maq:$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

main.o:main.c main.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DPACKAGE_VERSION=\"$(VERSION)\" main.c -o main.o

const.o:const.h

clean:
		rm -f *.o a.out *~ *.a $(PROG) $(MANPAGE)
