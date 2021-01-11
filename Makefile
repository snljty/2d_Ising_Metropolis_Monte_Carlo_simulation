# Makefile for 2d_Ising_Metropolis_Monte_Carlo_simulation.c

CC = gcc
CLINKER = gcc
CC_FLAGS = -O3
CLINKER_FLAGS = -O3

.PHONY: all

all: 2d_Ising_Metropolis_Monte_Carlo_simulation.exe

2d_Ising_Metropolis_Monte_Carlo_simulation.exe: 2d_Ising_Metropolis_Monte_Carlo_simulation.o
	$(CLINKER) $(CLINKER_FLAGS) -o $@ $< -l m

2d_Ising_Metropolis_Monte_Carlo_simulation.o: 2d_Ising_Metropolis_Monte_Carlo_simulation.c
	$(CC) $(CC_FLAGS) -o $@ -c $<

.PHONY: clean

clean:
	-del /q 2d_Ising_Metropolis_Monte_Carlo_simulation.exe 2d_Ising_Metropolis_Monte_Carlo_simulation.o 1>NUL 2>NUL
# -rm -f 2d_Ising_Metropolis_Monte_Carlo_simulation.exe 2d_Ising_Metropolis_Monte_Carlo_simulation.o 2>/dev/null

