CC = gcc
CFLAGS =-c -Wall
LDFLAGS = -lm -lgsl -lgslcblas
LIBS = -L/opt/local/lib/
INCLUDE = -I/opt/local/include
SOURCES = findEta.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = findEta

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LIBS) $(INCLUDE) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o $(EXECUTABLE)