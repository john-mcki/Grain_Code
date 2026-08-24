FC = mpif90 -O3
MODFILES = tensor_mod.f90 make_props.f90
OBJS = $(MODFILES:.f90=.o)
TARGET = make_props

$(TARGET): $(OBJS)
	$(FC) -o $@ $^

%.o: %.f90
	$(FC) -c -o $@ $<


clean:
	rm -f *.o *.mod $(EXE)
