
help:
	@echo "crust: "
	@echo "clean: "
	@echo "check: "
	@echo "doc: "
	@echo "sync-doc: "
	@echo "test-sync: "

# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------

# LIBS is the list of libraries
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags

# Default settings
LCXX = $(CXX)
LIBS = -L/usr/local/lib -lo2scl -lhdf5 -lgsl -lreadline $(LDFLAGS) 
LCFLAGS = -O3 -std=c++11 -DNO_MPI -fopenmp $(CFLAGS) \
	-DO2SCL_NO_BOOST_MULTIPRECISION

# ----------------------------------------------------------------
# UTK-specific settings
# ----------------------------------------------------------------

ifdef UTKNA_MAKEFILE

include $(UTKNA_MAKEFILE)

# UTK configuration

LIBS = $(UTKNA_O2SCL_LIBS) $(UTKNA_PYTHON_LDFLAGS)

LCXX = $(UTKNA_CXX) 
LCFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) -DNO_MPI \
        $(UTKNA_OPENMP_FLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION

endif

#----------------------------------------------------------------------

OBJS = crust.o main.o ldrop_crust.o pyc_rates.o eigen.o \
		matter.o rxns.o dist_thermo.o crust_fit.o multi_zone.o \
		nm_thermo.o sna_thermo.o

crust.o: crust.cpp crust.h
	$(LCXX) $(LCFLAGS) -o crust.o -c crust.cpp

crust_fit.o: crust_fit.cpp crust_fit.h
	$(LCXX) $(LCFLAGS) -o crust_fit.o -c crust_fit.cpp

dist_thermo.o: dist_thermo.cpp dist_thermo.h
	$(LCXX) $(LCFLAGS) -o dist_thermo.o -c dist_thermo.cpp

eigen.o: eigen.cpp eigen.h
	$(LCXX) $(LCFLAGS) -o eigen.o -c eigen.cpp

ldrop_crust.o: ldrop_crust.cpp ldrop_crust.h
	$(LCXX) $(LCFLAGS) -o ldrop_crust.o -c ldrop_crust.cpp

main.o: main.cpp 
	$(LCXX) $(LCFLAGS) -o main.o -c main.cpp

matter.o: matter.cpp matter.h
	$(LCXX) $(LCFLAGS) -o matter.o -c matter.cpp

multi_zone.o: multi_zone.cpp multi_zone.h
	$(LCXX) $(LCFLAGS) -o multi_zone.o -c multi_zone.cpp

nm_thermo.o: nm_thermo.cpp nm_thermo.h
	$(LCXX) $(LCFLAGS) -o nm_thermo.o -c nm_thermo.cpp

pyc_rates.o: pyc_rates.cpp pyc_rates.h
	$(LCXX) $(LCFLAGS) -o pyc_rates.o -c pyc_rates.cpp

rxns.o: rxns.cpp rxns.h
	$(LCXX) $(LCFLAGS) -o rxns.o -c rxns.cpp

sna_thermo.o: sna_thermo.cpp sna_thermo.h
	$(LCXX) $(LCFLAGS) -o sna_thermo.o -c sna_thermo.cpp

#----------------------------------------------------------------------

crust: $(OBJS)
	$(LCXX) $(LCFLAGS) -o crust $(OBJS) $(LIBS)

#----------------------------------------------------------------------

check: crust
	./crust -check 1 > check_1.scr
	tail -n 2 check_1.scr
	./crust -check 2 > check_2.scr
	tail -n 2 check_2.scr
	./crust -check 3 > check_3.scr
	tail -n 2 check_3.scr
	./crust -check 4 > check_4.scr
	tail -n 2 check_4.scr
	./crust -check 5 > check_5.scr
	tail -n 2 check_5.scr
	./crust -check 8 > check_8.scr
	tail -n 2 check_8.scr
	./crust -check 9 > check_9.scr
	tail -n 2 check_9.scr
	./crust -check 11 > check_11.scr
	tail -n 2 check_11.scr
	./crust -check 12 > check_12.scr
	tail -n 2 check_12.scr
	./crust -check 13 > check_13.scr
	tail -n 2 check_13.scr

#----------------------------------------------------------------------

doc: empty
# Copy most recent tag files
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .
# Get most recent commit hash
	git rev-parse HEAD | awk \
		'{print "`" $$1 " <http://github.com/awsteiner/crust/tree/" $$1 ">`_"}' \
		 > sphinx/commit.rst
# Parse bibliography
	cd sphinx/static; cat bib_header.txt > ../bib.rst
	cd sphinx/static; btmanip -parse crust.bib -rst ../bib_temp.rst
	cd sphinx; cat bib_temp.rst >> bib.rst; rm -f bib_temp.rst
# Run Doxygen
	cd doc; doxygen doxyfile
# Run sphinx
	cd sphinx; make html

sync-doc:
	rsync -Cavzu sphinx/build/html/* $(STATIC_DOC_DIR)/crust

test-sync:
	rsync -Cavzun sphinx/build/html/* $(STATIC_DOC_DIR)/crust


empty:

#----------------------------------------------------------------------

clean:
	rm -f *.o crust

