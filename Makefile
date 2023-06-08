include ./make.inc

SOURCES=Definition.f90 \
	Main.f90 \
	Buildhydro.f90 \
	Eostable.f90 \
	Fermimo.f90 \
	Getrho.f90 \
	Getrho_eosptor.f90 \
	Getrho_eosrtop.f90 \
	Getepsilon.f90 \
	Ini_der.f90 \
	Tabulareos.f90 

MODULES=Definition.mod
OBJECTS=$(SOURCES:.f90=.o )

TOV: Definition.o $(OBJECTS)  
	$(F90) $(LDFLAGS) -o ./TOV $(OBJECTS) 

$(OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) -c $< -o $@

TOV.o: Definition.o

clean:
	rm -rf Definition
	rm -rf *.o 	
	rm -rf *.mod
cleanfile:
	rm -rf ./outfile/*.dat
	rm -rf ./outfile/*.png
	rm -rf TOV
	rm -rf tmp.txt
