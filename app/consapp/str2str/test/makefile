# makefile for obscorr unit tests

SRC    = ../../../../src

OPTION = -DENAGLO -DENAGAL -DENAQZS -DENACMP -DENAIRN -DTRACE -DNFREQ=4 -DNEXOBS=3 -DSVR_REUSEADDR
CFLAGS = -std=c99 -Wall -O0 -g -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTION)
LDLIBS = -lm -lpthread

OBJS   = test_obscorr.o obscorr.o rtkcmn.o trace.o ephemeris.o rinex.o \
         preceph.o sbas.o geoid.o sofa.o solution.o stream.o streamsvr.o \
         rcvraw.o rtcm.o rtcm2.o rtcm3.o rtcm3e.o \
         novatel.o ublox.o crescent.o skytraq.o javad.o nvs.o binex.o \
         rt17.o septentrio.o swiftnav.o unicore.o

# test_refsatmode.c #includes rtkpos.c directly (ddres() has internal
# linkage), so rtkpos.o must NOT be linked in separately.
REFSATMODE_OBJS = test_refsatmode.o options.o rtkcmn.o trace.o ephemeris.o \
         rinex.o preceph.o sbas.o geoid.o sofa.o solution.o stream.o \
         streamsvr.o rcvraw.o rtcm.o rtcm2.o rtcm3.o rtcm3e.o \
         novatel.o ublox.o crescent.o skytraq.o javad.o nvs.o binex.o \
         rt17.o septentrio.o swiftnav.o unicore.o \
         pntpos.o ppp.o ppp_ar.o lambda.o tides.o obscorr.o ionex.o

all: test_obscorr test_refsatmode

test_obscorr: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDLIBS)

test_refsatmode: $(REFSATMODE_OBJS)
	$(CC) -o $@ $(REFSATMODE_OBJS) $(LDLIBS)

test_obscorr.o: test_obscorr.c $(SRC)/rtklib.h
	$(CC) -c $(CFLAGS) test_obscorr.c

test_refsatmode.o: test_refsatmode.c $(SRC)/rtkpos.c $(SRC)/rtklib.h
	$(CC) -c $(CFLAGS) test_refsatmode.c

%.o: $(SRC)/%.c $(SRC)/rtklib.h
	$(CC) -c $(CFLAGS) $<

%.o: $(SRC)/rcv/%.c $(SRC)/rtklib.h
	$(CC) -c $(CFLAGS) $<

test: test_obscorr test_refsatmode
	./test_obscorr
	./test_refsatmode

clean:
	rm -f test_obscorr test_refsatmode *.o test_obscorr_*.csv test_obscorr_*.bin test_refsatmode.conf
