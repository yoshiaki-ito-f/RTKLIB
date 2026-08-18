# makefile for snrmp
BINDIR  = /usr/local/bin
SRC     = ../../../../src

OPTS    = -DTRACE -DENAGLO -DENAQZS -DENAGAL -DENACMP -DENAIRN -DNFREQ=4 -DNEXOBS=3

CFLAGS  = -std=c99 -Wall -O3 -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTS)
LDLIBS  = -lm

OBJS    = snrmp.o rtkcmn.o trace.o rinex.o pntpos.o ephemeris.o preceph.o \
          ionex.o sbas.o options.o tides.o sofa.o

all        : snrmp
snrmp      : $(OBJS)

snrmp.o    : ../snrmp.c
	$(CC) -c $(CFLAGS) ../snrmp.c
rtkcmn.o   : $(SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(SRC)/rtkcmn.c
trace.o    : $(SRC)/trace.c
	$(CC) -c $(CFLAGS) $(SRC)/trace.c
rinex.o    : $(SRC)/rinex.c
	$(CC) -c $(CFLAGS) $(SRC)/rinex.c
pntpos.o   : $(SRC)/pntpos.c
	$(CC) -c $(CFLAGS) $(SRC)/pntpos.c
ephemeris.o: $(SRC)/ephemeris.c
	$(CC) -c $(CFLAGS) $(SRC)/ephemeris.c
preceph.o  : $(SRC)/preceph.c
	$(CC) -c $(CFLAGS) $(SRC)/preceph.c
ionex.o    : $(SRC)/ionex.c
	$(CC) -c $(CFLAGS) $(SRC)/ionex.c
sbas.o     : $(SRC)/sbas.c
	$(CC) -c $(CFLAGS) $(SRC)/sbas.c
options.o  : $(SRC)/options.c
	$(CC) -c $(CFLAGS) $(SRC)/options.c
tides.o    : $(SRC)/tides.c
	$(CC) -c $(CFLAGS) $(SRC)/tides.c
sofa.o     : $(SRC)/sofa.c
	$(CC) -c $(CFLAGS) $(SRC)/sofa.c

snrmp.o    : $(SRC)/rtklib.h
rtkcmn.o   : $(SRC)/rtklib.h
trace.o    : $(SRC)/rtklib.h
rinex.o    : $(SRC)/rtklib.h
pntpos.o   : $(SRC)/rtklib.h
ephemeris.o: $(SRC)/rtklib.h
preceph.o  : $(SRC)/rtklib.h
ionex.o    : $(SRC)/rtklib.h
sbas.o     : $(SRC)/rtklib.h
options.o  : $(SRC)/rtklib.h
tides.o    : $(SRC)/rtklib.h
sofa.o     : $(SRC)/rtklib.h

test :
	./snrmp ../../../../test/data/rinex/07590920.05o \
	        ../../../../test/data/rinex/07590920.05n test.txt

clean :
	rm -f snrmp snrmp.exe *.o *.txt *.trace

install :
	cp snrmp $(BINDIR)
