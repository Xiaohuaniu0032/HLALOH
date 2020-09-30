
#ifndef OSTYPE
OSTYPE=$(shell uname -s)
#endif
OS=$(shell cat /etc/issue | cut -f 1 "--delimiter= " | head -1)

ifdef CRAY_MPICH2
	INCLUDE=-I ${MPICH_DIR}/include
	CFLAGS=-m64 -c -O3 -DCRAY_MPICH2
	MPICC=g++
	LDLIBSOPTIONS=./bamtools/bamtools.a -pthread -lcrypto -lz -lbz2 -ldl -Wl,-rpath,/opt/gcc/4.8.1/snos/lib64 -l mpich -L ${MPICH_DIR}/lib
	MAP=-Wl,-Map=
else
	MPICC=mpicxx
	ifeq ($(OSTYPE),Darwin)
		MACHTYPE=$(shell echo $MACHTYPE)
		ifeq ($(MACHTYPE), x86_64-apple-darwin10.0)
			CFLAGS=-arch x86_64 -m64 -c -O3 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5 -I. -I/usr/local/include
			LNFLAGS=-m64 -mmacosx-version-min=10.5  -isysroot /Developer/SDKs/MacOSX10.5.sdk 
		else
		    CFLAGS=-arch x86_64 -m64 -c -O3 -I. -I/usr/local/include
		    LNFLAGS=-m64
		endif
		MAP=-Wl,-map 
		LDLIBSOPTIONS= ./bamtools/bamtools.a -lhdf5 -lhdf5_cpp -lpthread -lcrypto -lz -lbz2
	else
		ifeq ($(OSTYPE),SunOS)
			CFLAGS=-m64 -c -O3 -I/usr/sfw/include
			LNFLAGS=-m64
			MAP=-Wl,-m >
			LDLIBSOPTIONS= -lhdf5 -lhdf5_cpp -lpthread -lrt -lz -lbz2 -L /usr/sfw/lib/amd64 -R /usr/sfw/lib/amd64 -lcrypto -ldl
		else
			LNFLAGS=-m64 -Wl,--eh-frame-hdr
			MAP=-Wl,-Map=
			LDLIBSOPTIONS= ./tcmalloc/libtcmalloc.a ./unwind/libunwind.a -lhdf5 -lhdf5-cpp -pthread -lcrypto -lz -lbz2 -ldl
			ifeq ($(OS),Scientific)
				CFLAGS=-m64 -c -O3 -I./unwind -I/usr/lib64 
			else
				CFLAGS=-m64 -c -O3 -I./unwind
			endif
		endif
	endif
endif


# Build Targets
all: novoalignMPI novoalignCSMPI

novoalignMPI:	novoalignMPI.a mpidriver.o bamtools/bamtools.a trace.o
	echo >novoalignMPI.map # for Mac ld as map file must exist!!!
	${MPICC}  ${LNFLAGS} -o novoalignMPI mpidriver.o trace.o novoalignMPI.a bamtools/bamtools.a ${LDLIBSOPTIONS} ${MAP}novoalignMPI.map 

%.o: %.cc 
	${MPICC}  ${CFLAGS} ${INCLUDE}  -o $@ $?

novoalignCSMPI: novoalignCSMPI.a mpidriver.o trace.o
	echo >novoalignCSMPI.map # for Mac ld as map file must exist!!!
	${MPICC}  ${LNFLAGS} -o novoalignCSMPI  mpidriver.o trace.o novoalignCSMPI.a ${LDLIBSOPTIONS} ${MAP}novoalignCSMPI.map 

# Clean Targets
clean:
	rm novoalignMPI novoalignCSMPI mpidriver.o trace.o
