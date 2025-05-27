# Autoconf/automake file

bin_PROGRAMS = test_cca

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${ETREE_INCLUDES} ${PROJ_INCLUDES} -I../src
AM_LDFLAGS = ${LDFLAGS} ${ETREE_LDFLAGS} ${PROJ_LDFLAGS} -L../src -lcca -lm

objects = test_api.o
TARGETS = $(bin_PROGRAMS)

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/tests
	cp test_cca ${prefix}/tests

test_cca$(EXEEXT): $(objects)
	$(CC) -o $@ $^ $(AM_LDFLAGS)

$(objects): %.o: %.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

run_unit : test_cca
	./run_unit

clean :
	rm -rf *~ *.o test_cca

