EXEC = sequence-translator

CC = mpiCC --std=c++17 -lstdc++fs
CFLAGS = -Wall -Wextra -O3 -march=native
BOOST_LIBS = -lboost_program_options -lboost_filesystem -lboost_system
LDFLAGS = $(BOOST_LIBS) 

SRC = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, $(SRC))
OBJS = timer.o

.PHONY : clean

all: $(EXEC)

timer.o: timer.c
	$(CC) -O3 -o $@ -c $<

$(EXEC): $(OBJECTS) $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXEC) *.o
