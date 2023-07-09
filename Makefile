include   Mkinclude
##################################################
TARGET=calc_rivseq_link

all: $(TARGET)

clean:
	$(RM) -rf *.o *.s core *~ *trace *.mod *.dSYN $(TARGET)

.SUFFIXES : .F90
.F90:
	$(FC) $(FFLAGS) $(LFLAG) $(INC) $^ -o $@

