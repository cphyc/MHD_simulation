CC=g++ -std=c++11 -g
CFLAGS=-c -Wall
LFLAGS=
SOURCES=misc.cpp constant.cpp classes.cpp solve.cpp
OBJECTS=$(SOURCES:.cpp=.o) 
EXECUTABLE=simul

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)

rebuild: clean all
