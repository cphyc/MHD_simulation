CC=g++ -std=c++11 -g
CFLAGS=-c -Wall -ffast-math -fomit-frame-pointer -O2 -g
LFLAGS=-O2 -g
SOURCES=misc.cpp constant.cpp vector.cpp simulation.cpp solve.cpp 
OBJECTS=$(SOURCES:.cpp=.o) 
EXECUTABLE=simul

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS) *.o

rebuild: clean all
