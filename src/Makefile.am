# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libcca.a libcca.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcca.so ${prefix}/lib
	cp libcca.a ${prefix}/lib
	cp cca.h ${prefix}/include

libcca.a: cca_static.o utm_geo.o
	$(AR) rcs $@ $^

libcca.so: cca.o utm_geo.o
	$(FC) -shared $(AM_FCFLAGS) -o libcca.so $^ $(AM_LDFLAGS)

cca.o: cca.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
cca_static.o: cca.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
utm_geo.o: utm_geo.f90
	$(FC) -fPIC -o $@ -c $^ $(AM_FCFLAGS)