# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=mpicc

CFLAGS=-pedantic -Wall -Werror -std=gnu99  -Ofast -msse3 -mavx -funroll-loops



# Linker flags
# ------------
LDFLAGS=

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
