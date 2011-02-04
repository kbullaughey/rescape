CC = g++
HEADERS = landscape.h error_handling.h command_line.h popsim.h sim_rand.h \
	scapetypes.h stats.h array.h sequence.h configuration.h chain.h
OBJS = landscape.o error_handling.o sim_rand.o command_line.o popsim.o \
	stats.o segal.o array.o sequence.o main.o configuration.o chain.o \
	scapetypes.o
PLATFORM := $(shell uname -s)
CFLAGS = -Wall -O3 -fopenmp -fno-math-errno

ifeq ($(profile),yes)
CFLAGS += -g -pg 
else
CFLAGS += -O3
endif

ifeq ($(BENCHMARK),yes)
CFLAGS += -DBENCHMARK
endif

ifeq ($(PLATFORM),Linux)
LIBS = -lgsl -lm -lgslcblas
CFLAGS += -Wno-deprecated
else ifeq ($(PLATFORM),Darwin)
INC = -I/sw/include
LIBS = -L/sw/lib -lgsl -lm
else
	$(error Error: Unsupported platform)
endif

all: regevscape
	@echo benchmark=$(BENCHMARK)

regevscape: $(OBJS) $(HEADERS)
	$(CC) $(CFLAGS) $(INC) $(LIBS) -o $@ $(OBJS)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INC) -c $<

clean:
	-rm *.o regevscape

install:
	cp regevscape ../bin
