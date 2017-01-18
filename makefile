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

fit_test:
	crust -model SLy4 -mf

#----------------------------------------------------------------------

feq_all:
	crust -model SLy4 -mf data/SLy4_moller.fit \
		-rf data/SLy4_moller.fit -feq data/SLy4_moller \
		> data/SLy4_moller.scr 2> data/SLy4_moller.err &
	crust -set fit_moller 0 -model SLy4 -mf data/SLy4_ame.fit \
		-rf data/SLy4_ame.fit -feq data/SLy4_ame \
		> data/SLy4_ame.scr 2> data/SLy4_ame.err &

acc:
	crust -set dist_type ashes \
		-model SLy4 -rf data/SLy4_moller.fit -acc data/SLy4_moller
#		> data/SLy4_moller_acc.scr 2> data/SLy4_moller_acc.err &

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
	crust -model SLy4 -rf data/SLy4_moller.fit -check 1 > check_1.scr
	tail -n 2 check_1.scr
	crust -model SLy4 -rf data/SLy4_moller.fit -check 3 > check_3.scr
	tail -n 2 check_3.scr

check3: crust
	crust -model SLy4 -rf data/SLy4_moller.fit -check 3

#----------------------------------------------------------------------

doc: empty
# Copy most recent tag files
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .
# Get most recent commit hash
	git rev-parse HEAD | awk \
		'{print "<a href=\"http://github.com/awsteiner/crust/tree/" $$1 "\">" $$1 "</a> ."}' \
		 > doc/rev.txt
# Parse bibliography
	cd doc; cat refs_head.txt > refs.dox
	cd doc; btmanip -p crust.bib -dox /tmp/btmanip
	cd doc; cat /tmp/btmanip >> refs.dox
	cd doc; cat refs_foot.txt >> refs.dox
# Run Doxygen
	cd doc; doxygen doxyfile

empty:

#----------------------------------------------------------------------

clean:
	rm -f *.o crust

