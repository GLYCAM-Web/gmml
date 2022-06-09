# Makefile for GlycoProtein Builder
# Created by Davis Templeton 2018-07-23
# Last Modified on 2021-10-06 by Oliver Grant

CC = g++
CFLAGS = -std=c++17 -I ${GEMSHOME}/gmml/
RFLAGS = -Wl,-rpath,${GEMSHOME}/gmml/lib/
LFLAGS = -pthread -I ${GEMSHOME}/gmml/includes/ -L ${GEMSHOME}/gmml/lib/ -lgmml
COMPILE = $(CC) $(CFLAGS) -c
LINK = $(LFLAGS)
RUNTIME = $(RFLAGS)
DEBUG = -g 2
RM = rm -f

SRC = ${GEMSHOME}/gmml/src/InternalPrograms/CarbohydrateBuilder/
INC = ${GEMSHOME}/gmml/includes/InternalPrograms/CarbohydrateBuilder/
BUILD = ./build
BIN = ./bin

all:$(BIN)/carbBuilder \
	$(BUILD)/main.o \

$(BIN)/carbBuilder: \
		$(BUILD)/main.o 
	$(CC) \
		$(BUILD)/main.o \
	$(LINK) $(RUNTIME) -o $(BIN)/carbBuilder


$(BUILD)/main.o: ./main.cpp \
		$(INC)/carbohydrateBuilder.hpp 
	$(COMPILE) ./main.cpp -o $(BUILD)/main.o

clean:
	$(RM) $(BUILD)/*.o
	$(RM) $(BIN)/carbBuilder
