# GNU Automake config

SQLITE_URL=https://www.sqlite.org/2022/sqlite-autoconf-3370200.tar.gz
SQLITE_VER=3370200

TIFF_URL=http://download.osgeo.org/libtiff
TIFF_VER=4.3.0

PROJ_URL=https://download.osgeo.org/proj
PROJ_VER=8.2.1

EUCLID_URL=https://github.com/baagaard-usgs/euclid-etree/releases/download/v$(EUCLID_VER)
EUCLID_VER=3.1.4


all: installed_etree installed_proj

# ------------------------------------------------------------------------------
# ETREE/EUCLID
# ------------------------------------------------------------------------------
etree:
if INSTALL_ETREE
EUCLID_TARBALL=euclid-$(EUCLID_VER).tar.gz
	if [ ! -f "$(EUCLID_TARBALL)" ]; then $(CURL) -L -O $(EUCLID_URL)/$(EUCLID_TARBALL); fi
	$(TAR) -xf $(EUCLID_TARBALL)
	$(MKDIR_P) euclid-build
	cd euclid-build && \
		../euclid-$(EUCLID_VER)/configure \
			--prefix=$(prefix) \
			CC=$(CC) CPPFLAGS="$(CPPFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"; \
		make -C libsrc && make install -C libsrc && \
		touch $(abs_builddir)/installed_$@
else
	@echo "Using existing $@ installation."
	touch $(abs_builddir)/installed_$@
endif

installed_etree:
	$(MAKE) $(AM_MAKEFLAGS) etree


# ------------------------------------------------------------------------------
# SQLITE
# ------------------------------------------------------------------------------
sqlite:
if INSTALL_SQLITE
SQLITE_TARBALL=sqlite-autoconf-$(SQLITE_VER).tar.gz
	if [ ! -f "$(SQLITE_TARBALL)" ]; then $(CURL) -O $(SQLITE_URL)/$(SQLITE_TARBALL); fi
	$(TAR) -xf $(SQLITE_TARBALL)
	$(MKDIR_P) sqlite-build
	cd sqlite-build && \
		../sqlite-autoconf-$(SQLITE_VER)/configure \
			--prefix=$(prefix) --enable-shared \
			CC=$(CC) CPPFLAGS="$(CPPFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"; \
		make && \
		make install && \
		touch $(abs_builddir)/installed_sqlite

else
	@echo "Using existing $@ installation."
	touch $(abs_builddir)/installed_$@
endif

installed_sqlite:
	$(MAKE) $(AM_MAKEFLAGS) sqlite


# ------------------------------------------------------------------------------
# TIFF
# ------------------------------------------------------------------------------
tiff:
if INSTALL_TIFF
TIFF_TARBALL=tiff-$(TIFF_VER).tar.gz
	if [ ! -f "$(TIFF_TARBALL)" ]; then $(CURL) -O $(TIFF_URL)/$(TIFF_TARBALL); fi
	$(TAR) -xf $(TIFF_TARBALL)
	$(MKDIR_P) tiff-build
	cd tiff-build && \
		../tiff-$(TIFF_VER)/configure \
			--prefix=$(prefix) ; \
			CC=$(CC) CPPFLAGS="$(CPPFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"; \
		make -j $(make_threads) && \
		make install && \
		touch $(abs_builddir)/installed_$@

else
	@echo "Using existing $@ installation."
	touch $(abs_builddir)/installed_$@
endif

installed_tiff:
	$(MAKE) $(AM_MAKEFLAGS) tiff


# ------------------------------------------------------------------------------
# PROJ
# ------------------------------------------------------------------------------
proj: installed_sqlite installed_tiff
if INSTALL_PROJ
PROJ_TARBALL=proj-$(PROJ_VER).tar.gz
	if [ ! -f "$(PROJ_TARBALL)" ]; then $(CURL) -O $(PROJ_URL)/$(PROJ_TARBALL); fi
	$(TAR) -xf $(PROJ_TARBALL)
	$(MKDIR_P) proj-build
	cd proj-build && \
		../proj-$(PROJ_VER)/configure \
			--prefix=$(prefix) \
			SQLITE3_CFLAGS=$(SQLITE3_INCLUDES) SQLITE3_LIBS="$(SQLITE3_LDFLAGS) -lsqlite3" \
			TIFF_CFLAGS=$(TIFF_INCLUDES) TIFF_LIBS="$(TIFF_LDFLAGS) -ltiff" \
			CC=$(CC) CPPFLAGS="$(CPPFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"; \
		make && make install && \
		projsync --system-directory --all && \
		touch $(abs_builddir)/installed_$@
else
	@echo "Using existing $@ installation."
	touch $(abs_builddir)/installed_$@
endif

installed_proj:
	$(MAKE) $(AM_MAKEFLAGS) proj


# End of file
