SHELL = /bin/sh
ARCH = macx
TARGET1= main
TARGETS= $(TARGET1)
SRCS1= SUB_main.f  \
dlsode.f  \
rg.f  \
balanc.f  \
balbak.f  \
cdiv.f  \
d1mach.f  \
elmhes.f  \
eltran.f  \
hqr.f  \
hqr2.f  \
SUB_eigen.f  \
SUB_sinve.f\
SUB_smult.f\
SUB_unitary.f\
BLM_ICs.f \
BLM_mol2conc.f \
BLM_chparam.f \
SUB_BLM.f \
SUB_sort.f \
CSP_kern.f \
SUB_diagnostics.f
CSRCS1= 
OBJS1= $(SRCS1:.f=.o) $(CSRCS1:.c=.o)
DEPSRC = $(SRCS1) 
DEPCSRC = $(CSRCS1)  
LIBDIR =
LIBS1=
DEFS = 
#FFLAGS =   -O -fno-automatic -Wunused -Wsurprising $(DEFS) 
FFLAGS =  -O3 $(DEFS) 
LDFLAGS = 
INCL = 
F77 = gfortran	
FC = $(F77)	
CC = cc	
CFLAGS =   $(DEFS) -O	
CPP = gcc	
MKDEP = 	

all: $(TARGETS)	

$(TARGET1): $(OBJS1)	
	$(F77) $(LDFLAGS) $(FFLAGS) $(INCL) -o $(TARGET1) $(OBJS1) $(LIBDIR) $(LIBS1)

clean:	
	rm -f $(OBJS1)  Makefile.bak

depend:	
	@(sed '/^# DO NOT DELETE THIS LINE/q' Makefile && \
	 $(MKDEP) CPP=$(CPP) $(DEPCSRC) $(DEPSRC) \
	) > Makefile.new
	@mv Makefile Makefile.bak
	@mv Makefile.new Makefile

.f.o:	
	$(F77) $(FFLAGS) -c $*.f

#	Dependencies of F77 object files
SUB_main.o:		paramet.i

# DO NOT DELETE THIS LINE	


