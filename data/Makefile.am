# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
if WITH_CCA_LARGEDATA_DIR
	echo "USING WITH_CCA_LARGEDATA_DIR in cca"
else
	if [ ! -d "${prefix}/data/i06" ]; then \
		./make_data_files.py; \
	fi
endif

install:
if WITH_CCA_LARGEDATA_DIR
	if [ ! -d "${prefix}/data/i06" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		$(LN_S) -f $(CCA_LARGEDATA_DIR) ${prefix}/data/i06; \
		echo "USING WITH_CCA_LARGEDATA_DIR in cca"; \
	fi
	cp config ${prefix}/data
else
	if [ ! -d "${prefix}/data/i06" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		mv i06 ${prefix}/data; \
		$(LN_S) -f ${prefix}/data/i06 .; \
	fi
	cp config ${prefix}/data
endif

clean:
	echo "clean"

