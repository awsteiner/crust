
ifeq ($(USER),awsteiner)

else

GSL_INC = .
O2SCL_INC = .
EIGEN_INC = .
HDF5_INC = .
GSL_LIB = .
O2SCL_LIB = .
HDF5_LIB = .

endif

FLAGS = -std=c++0x -ggdb -Wreturn-type -Wparentheses -Wall \
	-Wno-unused -Wno-array-bounds -O3 -DO2SCL_READLINE \
	-I$(O2SCL_INC) -I$(GSL_INC) -I$(HDF5_INC) -I$(BOOST_INC) \
	-I$(EIGEN_INC) -DO2SCL_HDF_SVAR \
	-Wno-deprecated-declarations 

LIB = -L$(O2SCL_LIB) -L$(GSL_LIB) -L$(HDF5_LIB) -lo2scl_eos \
	-lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 \
	-lgsl -lgslcblas -lreadline -lncurses

OBJS = crust.o main.o ldrop_crust.o pyc_rates.o eigen.o \
		matter.o rxns.o dist_thermo.o crust_fit.o multi_zone.o \
		nm_thermo.o sna_thermo.o 

#----------------------------------------------------------------------

crust.o: crust.cpp crust.h
	$(CXX) $(FLAGS) -o crust.o -c crust.cpp

crust_fit.o: crust_fit.cpp crust_fit.h
	$(CXX) $(FLAGS) -o crust_fit.o -c crust_fit.cpp

dist_thermo.o: dist_thermo.cpp dist_thermo.h
	$(CXX) $(FLAGS) -o dist_thermo.o -c dist_thermo.cpp

eigen.o: eigen.cpp eigen.h
	$(CXX) $(FLAGS) -o eigen.o -c eigen.cpp

ldrop_crust.o: ldrop_crust.cpp ldrop_crust.h
	$(CXX) $(FLAGS) -o ldrop_crust.o -c ldrop_crust.cpp

main.o: main.cpp 
	$(CXX) $(FLAGS) -o main.o -c main.cpp

matter.o: matter.cpp matter.h
	$(CXX) $(FLAGS) -o matter.o -c matter.cpp

multi_zone.o: multi_zone.cpp multi_zone.h
	$(CXX) $(FLAGS) -o multi_zone.o -c multi_zone.cpp

nm_thermo.o: nm_thermo.cpp nm_thermo.h
	$(CXX) $(FLAGS) -o nm_thermo.o -c nm_thermo.cpp

pyc_rates.o: pyc_rates.cpp pyc_rates.h
	$(CXX) $(FLAGS) -o pyc_rates.o -c pyc_rates.cpp

rxns.o: rxns.cpp rxns.h
	$(CXX) $(FLAGS) -o rxns.o -c rxns.cpp

sna_thermo.o: sna_thermo.cpp sna_thermo.h
	$(CXX) $(FLAGS) -o sna_thermo.o -c sna_thermo.cpp

#----------------------------------------------------------------------

crust: $(OBJS)
	$(CXX) $(FLAGS) -o crust $(OBJS) $(LIB)

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
		'{print "`" $$1 " <http://github.com/awsteiner/bamr/tree/" $$1 ">`_"}' \
		 > sphinx/commit.rst
# Parse bibliography
	cd sphinx/static; cat bib_header.txt > ../bib.rst
	cd sphinx/static; btmanip -parse crust.bib -rst ../bib_temp.rst
	cd sphinx; cat bib_temp.rst >> bib.rst; rm -f bib_temp.rst
# Run Doxygen
	cd doc; doxygen doxyfile
# Run sphinx
	cd sphinx; make html
# Copy to web
	cp -r sphinx/build/html/* $(HOME)/wcs/int4/web/utk/crust

empty:

#----------------------------------------------------------------------

clean:
	rm -f *.o crust

