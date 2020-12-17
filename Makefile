EXEC = sequence-translator sequence-translatorOPT

CC = mpiCC --std=c++17 -lstdc++fs
CFLAGS = -Wall -Wextra -O3 -march=native
BOOST_LIBS = -lboost_program_options -lboost_filesystem -lboost_system
LDFLAGS = $(BOOST_LIBS) 

SRC = $(wildcard *.cc)
OBJS = timer.o AminoAcid.o FastaParser.o FastaRecord.o FastaWriter.o GeneticCode.o NucleicAcid.o Parameters.o SequenceTranslator.o
H_FILES = timer.h AminoAcid.hh FastaParser.hh FastaRecord.hh FastaWriter.hh GeneticCode.hh NucleicAcid.hh Parameters.hh SequenceTranslator.hh

.PHONY : clean

all: $(EXEC)

timer.o: timer.c
	$(CC) -O3 -o $@ -c $<

sequence-translator: main.cc $(OBJS) $(H_FILE)
	g++ --std=c++17 -lstdc++fs -O3 -o $@ main.cc $(OBJS) $(LDFLAGS)

sequence-translatorOPT: opt.cc $(OBJECTS) $(OBJS)
	$(CC) -O3 -o $@ opt.cc $(OBJS) $(LDFLAGS)

clean:
	rm -f $(EXEC) *.o
