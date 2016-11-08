CC = g++ -c
LD = g++
OPT = -lm -g -pg -O3 --fast-math
C_FLAGS = -g -Wall -W $(OPT)
O_FLAGS = $(C_FLAGS)
NAME = a.out

all: $(NAME)

$(NAME): main.o funcs.o
	$(LD) $(O_FLAGS) $^ -o $@
   
main.o: main.cpp funcs.h
	$(CC) $(C_FLAGS) $< -o $@
   
funcs.o: funcs.cpp 	
  

clean:
	rm -rf main.o funcs.o$(NAME)
