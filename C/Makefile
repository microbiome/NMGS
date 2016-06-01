CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = NMGS
LIBS   = -lm -lgsl -lgslcblas
OBJS   = NMGS.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : NMGS.c NMGS.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o NMGS
