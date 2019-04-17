
OBJECTS = routines.o math.o constants.o  main.o 
        

MODULES = math.mod constants.mod
          

.PHONY: clean

#output.txt: main.exe
#	./main.exe > output.txt

LAPACK =  -llapack -lblas
FLAGS  = 
#FLAGS  = -g -Wall -Wextra -Warray-temporaries -Wconversion  -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow # -finit-real=nan 
main.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS)  $(LAPACK) $(FLAGS) -o main.exe

%.o: %.f90
	gfortran -c $(LAPACK) $(FLAGS)$<

%.mod: %.f90
	gfortran -c $(LAPACK) $(FLAGS) $<

clean:
	rm -f $(OBJECTS) $(MODULES) main.exe
