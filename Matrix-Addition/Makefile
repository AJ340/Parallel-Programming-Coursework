

MPICC      = mpicc

SRCS       = utilities.c  matrix_io.c main.c
OBJS       = $(patsubst %.c,%.o,$(SRCS))
EXECS      = addmatrices
CFLAGS     = -Wall 


%.o : %.c
	$(MPICC) -c $< -o $@ $(CFLAGS)

all:
	make $(EXECS)


addmatrices: $(OBJS)
	$(MPICC) $(C++FLAG) -o $@ $(OBJS)


clean:
	-rm -f *.o


cleanall:
	-rm -f $(OBJS) $(EXECS) 

