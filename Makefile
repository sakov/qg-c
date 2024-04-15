PROGRAMS =\
qg

CFLAGS = -g -O2
FFLAGS = -g -frecursive -O2
VERSION := $(shell grep QG_VERSION version.c | tr "\"" " " | awk '{print $$5}')

LIB_NC = -lnetcdf -lhdf5 -lhdf5_hl
# On Ubuntu the above should be set as follows
# LIB_NC = -lnetcdf -lhdf5_serial -lhdf5_serial_hl
LIBS = $(LIB_NC) -lgfortran -lm

SRC =\
ncw.c\
version.c\
utils.c\
qgprm.c\
model.c\
calcs.c\
qgstep.c\
qg.c

HDR=\
ncw.h\
utils.h\
qg.h\
qgprm.h\
model.h\
calcs.h\
qgstep.h\
helmholtz.h

HELMHOLTZ = helmholtz.o

default: $(PROGRAMS)

helmholtz.o: helmholtz.f90
	gfortran $(FFLAGS) -c helmholtz.f90

qg: $(SRC) $(HDR) $(HELMHOLTZ) Makefile
	gcc $(CFLAGS) -pedantic -Wall -std=c99 -D_GNU_SOURCE -o qg $(SRC) $(HELMHOLTZ) $(LIBS)

clean:
	rm -f $(PROGRAMS) *.o *.a *.mod *.aux *.log

tar:
	make clean; cd ..; tar --exclude=redundant --exclude='*.nc*' -czvf qg-c-v$(VERSION).tar.gz qg-c; echo "  ->../qg-c-v$(VERSION).tar.gz"

pdf:
	pdflatex qg.tex; rm -f qg.aux qg.log

indent:
	indent -T FILE -T size_t -T nc_type -T qgprm -T model *.[ch]; rm -f *.[ch]~
