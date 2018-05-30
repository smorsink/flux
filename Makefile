# Source (C) Coire Cadeau 2007, all rights reserved.

# Permission is granted for private use only, and not
# distribution, either verbatim or of derivative works,
# in whole or in part.

# The code is not thoroughly tested or guaranteed for
# any particular use.

CC=g++
CCFLAGS=-Wall -pedantic -O3
LDFLAGS=-lm

# Location of matpack.a static library
# MATPACK=/Users/coire/build/matpack/matpack.a -L/usr/X11R6/lib -lXpm -lX11 -lm
# MPINCLUDE=-DXPM_INCLUDE="<X11/xpm.h>" -I /Users/coire/build/matpack/include

NAMES=flux

OBJ=PolyOblModelBase.o  PolyOblModelCFLQS.o PolyOblModelNHQS.o Units.o OblDeflectionTOA.o \
	Chi.o SphericalOblModel.o matpack.o

APPOBJ=Flux.o

all: $(NAMES)

flux: Solid.o $(OBJ)
	$(CC) $(CCFLAGS) Solid.o $(OBJ) $(LDFLAGS) -o flux

Solid.o: \
	Solid.cpp \
	OblDeflectionTOA.h \
	Chi.h \
	Struct.h \
	PolyOblModelNHQS.h \
	PolyOblModelCFLQS.h \
	SphericalOblModel.h \
	OblModelBase.h \
	Units.h \
	Makefile
	$(CC) $(CCFLAGS) -c Solid.cpp

PolyOblModelBase.o: \
	PolyOblModelBase.h \
	PolyOblModelBase.cpp \
	OblModelBase.h
	$(CC) $(CCFLAGS) -c PolyOblModelBase.cpp

PolyOblModelCFLQS.o: \
	PolyOblModelCFLQS.h \
	PolyOblModelCFLQS.cpp \
	PolyOblModelBase.h
	$(CC) $(CCFLAGS) -c PolyOblModelCFLQS.cpp

PolyOblModelNHQS.o: \
	PolyOblModelNHQS.h \
	PolyOblModelNHQS.cpp \
	PolyOblModelBase.h
	$(CC) $(CCFLAGS) -c PolyOblModelNHQS.cpp


SphericalOblModel.o: \
	SphericalOblModel.h \
	SphericalOblModel.cpp \
	OblModelBase.h
	$(CC) $(CCFLAGS) -c SphericalOblModel.cpp

OblDeflectionTOA.o: \
	OblDeflectionTOA.h \
	OblDeflectionTOA.cpp \
	OblModelBase.h \
	Units.h \
	matpack.h
#	$(CC) $(CCFLAGS) $(MPINCLUDE) -c OblDeflectionTOA.cpp
	$(CC) $(CCFLAGS) -c OblDeflectionTOA.cpp


Chi.o: \
	Chi.h \
	OblDeflectionTOA.h \
	Chi.cpp \
	OblModelBase.h \
	Units.h \
	matpack.h
#	$(CC) $(CCFLAGS) $(MPINCLUDE) -c OblDeflectionTOA.cpp
	$(CC) $(CCFLAGS) -c Chi.cpp


Units.o: \
	Units.h \
	Units.cpp
	$(CC) $(CCFLAGS) -c Units.cpp

matpack.o: \
	matpack.h \
	matpack.cpp \
	Exception.h
	$(CC) $(CCFLAGS) -c matpack.cpp

clean:
	rm -f core *~ $(OBJ) $(APPOBJ)

veryclean:
	rm -f core *~ $(OBJ) $(APPOBJ) $(NAMES)
