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

SRC = ${GEMSHOME}/gmml/src/InternalPrograms/
SRCGP = ${GEMSHOME}/gmml/src/InternalPrograms/GlycoproteinBuilder/
INC = ${GEMSHOME}/gmml/includes/InternalPrograms/
INCGP = ${GEMSHOME}/gmml/includes/InternalPrograms/GlycoproteinBuilder/
BUILD = ./build
BIN = ./bin

all:$(BIN)/gpBuilder \
	$(BUILD)/beadResidues.o \
	$(BUILD)/glycoproteinBuilder.o \
	$(BUILD)/glycosylationSite.o \
    $(BUILD)/gpInputStructs.o \
	$(BUILD)/io.o \
	$(BUILD)/main.o \
	$(BUILD)/functionsForGMML.o \

$(BIN)/gpBuilder: \
		$(BUILD)/beadResidues.o \
		$(BUILD)/glycoproteinBuilder.o \
		$(BUILD)/glycosylationSite.o \
        $(BUILD)/gpInputStructs.o \
		$(BUILD)/io.o \
		$(BUILD)/functionsForGMML.o \
		$(BUILD)/main.o 
	$(CC) \
		$(BUILD)/beadResidues.o \
		$(BUILD)/glycoproteinBuilder.o \
		$(BUILD)/glycosylationSite.o \
        $(BUILD)/gpInputStructs.o \
		$(BUILD)/io.o \
		$(BUILD)/functionsForGMML.o \
		$(BUILD)/main.o \
	$(LINK) $(RUNTIME) -o $(BIN)/gpBuilder

$(BUILD)/beadResidues.o: $(SRC)/beadResidues.cpp \
		$(INC)/beadResidues.hpp \
		$(INCGP)/glycosylationSite.hpp 
	$(COMPILE) $(SRC)/beadResidues.cpp -o $(BUILD)/beadResidues.o

$(BUILD)/glycoproteinBuilder.o: $(SRCGP)/glycoproteinBuilder.cpp \
		$(INCGP)/glycoproteinBuilder.hpp \
		$(INC)/io.hpp \
        $(INCGP)/gpInputStructs.hpp \
		$(INCGP)/glycosylationSite.hpp \
		$(INC)/beadResidues.hpp \
		$(INC)/functionsForGMML.hpp
	$(COMPILE) $(SRCGP)/glycoproteinBuilder.cpp -o $(BUILD)/glycoproteinBuilder.o

$(BUILD)/glycosylationSite.o: $(SRCGP)/glycosylationSite.cpp \
		$(INCGP)/glycosylationSite.hpp \
		$(INC)/functionsForGMML.hpp 
	$(COMPILE) $(SRCGP)/glycosylationSite.cpp -o $(BUILD)/glycosylationSite.o

$(BUILD)/io.o: $(SRC)/io.cpp \
		$(INC)/io.hpp
	$(COMPILE) $(SRC)/io.cpp -o $(BUILD)/io.o

$(BUILD)/gpInputStructs.o: $(SRCGP)/gpInputStructs.cpp \
		$(INCGP)/gpInputStructs.hpp
	$(COMPILE) $(SRCGP)/gpInputStructs.cpp -o $(BUILD)/gpInputStructs.o

$(BUILD)/functionsForGMML.o: $(SRC)/functionsForGMML.cpp \
		$(INC)/functionsForGMML.hpp
	$(COMPILE) $(SRC)/functionsForGMML.cpp -o $(BUILD)/functionsForGMML.o
		
$(BUILD)/main.o: ./main.cpp \
		$(INC)/io.hpp \
		$(INCGP)/glycoproteinBuilder.hpp \
        $(INCGP)/gpInputStructs.hpp
	$(COMPILE) ./main.cpp -o $(BUILD)/main.o

clean:
	$(RM) $(BUILD)/*.o
	$(RM) $(BIN)/gpBuilder
