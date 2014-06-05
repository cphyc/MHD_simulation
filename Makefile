CC=g++ -std=c++11 -lm -lfftw3 -O2
CFLAGS=-c -Wall
LFLAGS=
SOURCES=misc.cpp constant.cpp vector.cpp simulation.cpp solve.cpp
OBJECTS=$(SOURCES:.cpp=.o)
CLEAN= $(OBJECTS)
EXECUTABLE=simul

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(CLEAN)

rebuild: clean all
