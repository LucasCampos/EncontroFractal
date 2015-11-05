CC          = g++
MAKE        = make
SHELL       = /bin/sh

CFLAGS  = -O3 -std=c++11 -march=native
FOLDERS = 
HEADERS = ../inc
LIBS = 
EXECS = NPH
SOURCES = main.cpp
OBJECTS=$(SOURCES:.cpp=.o)

$(EXECS): $(OBJECTS)
	$(CC) $^ -o $(EXECS) $(LIBS) $(CFLAGS) 

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS) $(FOLDERS) -I $(HEADERS)

.PHONY:clean
clean:
	/bin/rm -f $(OBJECTS) $(EXECS) *.gif *.dat

.PHONY:run
run: $(EXECS)
	./$(EXECS)
