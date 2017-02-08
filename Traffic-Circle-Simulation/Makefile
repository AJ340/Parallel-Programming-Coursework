
MPICC      = mpicc

SRCS       = main.c
OBJS       = $(patsubst %.c,%.o,$(SRCS))
EXECS      = traffic_circle 

CFLAGS     = -Wall -lm


%.o : %.c
	$(MPICC) -c $< -o $@ $(CFLAGS) 

all:
	make $(EXECS)

traffic_circle: $(OBJS)
	$(MPICC) -o $@ $(OBJS) $(CFLAGS) 


clean:
	-rm -f *.o


cleanall:
	-rm -f $(OBJS) $(EXECS) 

