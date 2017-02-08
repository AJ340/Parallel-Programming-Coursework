
MPICC      = mpicc

SRCS       = utilities.c  matrix_io.c main.c
OBJS       = $(patsubst %.c,%.o,$(SRCS))
EXECS      = sunShine
CFLAGS     = -Wall 


%.o : %.c
	$(MPICC) -c $< -o $@ $(CFLAGS)

all:
	make $(EXECS)


sunShine: $(OBJS)
	$(MPICC) $(C++FLAG) -o $@ $(OBJS) -lm


clean:
	-rm -f *.o


cleanall:
	-rm -f $(OBJS) $(EXECS) 

