# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${ETREE_INCLUDES} ${PROJ_INCLUDES}
AM_LDFLAGS = ${LDFLAGS} ${ETREE_LDFLAGS} ${PROJ_LDFLAGS}


TARGETS = libcca.a libcca.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcca.so ${prefix}/lib
	cp libcca.a ${prefix}/lib
	cp cca.h ${prefix}/include

libcca.a: cca_static.o 
	$(AR) rcs $@ $^

libcca.so: cca.o 
	$(CC) -shared $(AM_CCFLAGS) -o libcca.so $^ $(AM_LDFLAGS)

cca.o: cca.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
cca_static.o: cca.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
clean:
	rm -rf $(TARGETS)
	rm -rf cca.o 

