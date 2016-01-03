# Include files
SOURCES=initLB.c visualLB.c  boundary.c streaming.c computeCellValues.c main.c helper.c communication.c

# Compiler
# --------
CC=mpicc

CFLAGS=-pedantic -Wall -std=gnu99 -fopenmp -Ofast -msse3 -funroll-loops -g $(PAPI_INC)



# Linker flags
# ------------
LDFLAGS= -fopenmp $(PAPI_LIB)
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
