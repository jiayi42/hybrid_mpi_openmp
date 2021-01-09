CC = gcc
CXX = g++
LDLIBS = -lpng -lpthread
CFLAGS = -lm -O3
pp2a: CFLAGS += -pthread 
pp2b: CXX = mpicxx
pp2b: CFLAGS += -fopenmp

CXXFLAGS = $(CFLAGS)
TARGETS = pp2a pp2b 

.PHONY: all
all: $(TARGETS)

.PHONY: clean
clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
