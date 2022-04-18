# Makefile for 2d_Ising_Metropolis_Monte_Carlo_simulation.c

SHELL = cmd
CC = gcc
CLINKER = gcc
CC_FLAGS = -O2
CLINKER_FLAGS = -O2 -static -s

.PHONY: all

all: 2d_Ising_Metropolis_Monte_Carlo_simulation

.PHONY: 2d_Ising_Metropolis_Monte_Carlo_simulation

2d_Ising_Metropolis_Monte_Carlo_simulation: 2d_Ising_Metropolis_Monte_Carlo_simulation.exe

2d_Ising_Metropolis_Monte_Carlo_simulation.exe: 2d_Ising_Metropolis_Monte_Carlo_simulation.obj
	$(CLINKER) $(CLINKER_FLAGS) -o $@ $< -l m

2d_Ising_Metropolis_Monte_Carlo_simulation.obj: 2d_Ising_Metropolis_Monte_Carlo_simulation.c
	$(CC) $(CC_FLAGS) -o $@ -c $<

.PHONY: clean

clean:
	-del /q 2d_Ising_Metropolis_Monte_Carlo_simulation.exe 2d_Ising_Metropolis_Monte_Carlo_simulation.obj 2> NUL

