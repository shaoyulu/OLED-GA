# Use this Makefile with make

# Executable name
CMD = OLEDsearch.exe

# -------- description of DFLAGS ---------------


# -------- Define environmental variable C_COMPILER -----------
# Make sure it is defined
#          ifeq ($(strip$(FORTRAN_COMPILER)),)
# Otherwise you can define it here also by uncommenting next line
FC = g++ -std=c++11 -I$(MKLROOT)/include
CFLAGS = -c -I /export/zimmerman/shaoyulu/local/include/openbabel-2.0
LDFLAGS=-rdynamic $(HOME)/local/lib/libopenbabel.so.4.0.2 -lm -ldl -Wl,-rpath,$(HOME)/local/lib -Wl,--enable-new-dtags -Wl,--fatal-warnings -Wl,--no-undefined -lc -Wl,--enable-new-dtags -Wl,--fatal-warnings -Wl,--no-undefined

DFLAGS =  -DGEOMETRIC #-Define the cpp flags to be used
#DFLAGS = -DGEOMETRIC #-Define the cpp flags to be used
OFLAGS =  # optimization
F95ROOT = $(MKLROOT)
#Intel Linkers
LINKERFLAGS =  -L$(MKLROOT)/lib/intel64 $(F95ROOT)/lib/intel64/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#Intel parallel openmp (only w/icpc compiler)
#LINKERFLAGS =  -L$(MKLROOT)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
# MAC OS linkers
#LINKERFLAGS = -lm -framework Accelerate



#
# Implicit rules for handling src files
#  ( uses static pattern rules, see info for make )
.c.o:
	$(FC) -c -g $(DFLAGS) -Wimplicit $<
.cpp.o:
	$(FC) -g -O2 $(LDFLAGS) $(DFLAGS) $(CFLAGS) $<

OBJECTS = main.o pTable.o mol2xyz.o OBabel.o nbo.o utils.o stringtools.o read_inp.o fragment.o mopac.o gcode.o  simple_genetic.o

$(CMD) : $(OBJECTS)
	$(FC) $(LDFLAGS) $(DEBUG_FLAGS) $(OFLAGS) $(LINKERFLAGS) $(OBJECTS)  -o $(CMD)

clean:
	/bin/rm -f *.o *.i *.mod *.exe a.out make.log

cleano:
	rm -f *.o *.i

depend :
	g++ -MM *.cpp *.c >> Makefile 

# DO NOT DELETE created with g++ -MM *.cpp *.c
main.o: main.cpp
mol2xyz.o: mol2xyz.cpp mol2xyz.h
read_inp.o: read_inp.cpp read_inp.h
fragment.o: fragment.cpp fragment.h
#DFT.o: DFT.cpp DFT.h 
#fitness_f.o: fitness.cpp fitness.h read_inp.h xyzlist.h CatsGen.h DFT.h 
#CatsGen.o: CatsGen.cpp CatsGen.h icoord.h utils.h 
#linkedlist.o: linkedlist.h
#node.o: node.h
#xyzlist.o: xyzlist.cpp xyzlist.h CatsGen.h
#genetic.o: genetic.cpp genetic.h xyzlist.h CatsGen.h DFT.h utils.h mopac.h fitness.h
#dft.o: dft.cpp dft.h constants.h
#enhanced_init_mut.o: enhanced_init_mut.cpp
simple_genetic.o: simple_genetic.cpp simple_genetic.h
#genetic.o: genetic.cpp genetic.h
OBabel.o: OBabel.cpp OBabel.h
#icoord.o: icoord.cpp icoord.h 
gcode.o: gcode.cpp gcode.h 
#iso.o: iso.cpp icoord.h genetic.h
#mm_grad.o: mm_grad.cpp icoord.h
mopac.o: mopac.cpp mopac.h
nbo.o: nbo.cpp nbo.h
#mem.o: mem.cpp icoord.h
#opt.o: opt.cpp icoord.h
#optic.o: optic.cpp icoord.h
pTable.o: pTable.cpp pTable.h
#print.o: print.cpp icoord.h
stringtools.o: stringtools.cpp stringtools.h
utils.o: utils.cpp utils.h 

