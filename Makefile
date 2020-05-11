GSL_INCLUDE = /usr/local/include
GSL_LIBRARY= /usr/local/lib
CFLAGS = -DNDEBUG -O -I$(GSL_INCLUDE) -I. -std=c++11 -Wall -msse2
LFLAGS = -lm -lgsl -lgslcblas -L$(GSL_LIBRARY)
GPP = g++

OBJECT_HIV = examples/HIV/HIVlangil.o
OBJECT_CYCLE = examples/CellCycleVariability/cyclo_2states_pop.o
OBJECT_CYCLE_POP = examples/CellCycleVariability/cyclo_2states.o

SRC_HIV = examples/HIV/HIVlangil.cpp
SRC_CYCLE = examples/CellCycleVariability/cyclo_2states_pop.cpp
SRC_CYCLE_POP = examples/CellCycleVariability/cyclo_2states.cpp


#############
# The default target only creates the main object that contains the 
# definition of the main routines
langil.o: langil.cpp langil.h
	@echo "Compiling langil..."
	$(GPP) $(CFLAGS) -c langil.cpp


#########################
#The target all, includes also the codes for the different examples
all: HIVlangil cyclo_2states cyclo_2states_pop

########
HIVlangil: $(OBJECT_HIV) langil.o
	$(GPP) $(CFLAGS) -o HIVlangil $(OBJECT_HIV) langil.o $(LFLAGS)
	mv HIVlangil examples/HIV/

$(OBJECT_HIV): $(SRC_HIV)
	@echo "Compiling HIV example..."
	$(GPP) $(CFLAGS) -c -o $(OBJECT_HIV) $(SRC_HIV)

########
cyclo_2states: $(OBJECT_CYCLE) langil.o
	$(GPP) $(CFLAGS) $(LFLAGS) -o cyclo_2states $(OBJECT_CYCLE) langil.o
	mv cyclo_2states examples/CellCycleVariability/

$(OBJECT_CYCLE): $(SRC_CYCLE)
	@echo "Compiling cyclo_2states example..."
	$(GPP) $(CFLAGS) -c -o $(OBJECT_CYCLE) $(SRC_CYCLE)

#############
cyclo_2states_pop: $(OBJECT_CYCLE_POP) langil.o
	$(GPP) $(CFLAGS) $(LFLAGS) -o cyclo_2states_pop $(OBJECT_CYCLE_POP) langil.o
	mv cyclo_2states_pop examples/CellCycleVariability/

$(OBJECT_CYCLE_POP): $(SRC_CYCLE_POP)
	@echo "Compiling cyclo_2states_pop example..."
	$(GPP) $(CFLAGS) -c -o $(OBJECT_CYCLE_POP) $(SRC_CYCLE_POP)


#############

.PHONY: clean

clean:
	$(RM) langil.o
	$(RM) $(OBJECT_HIV)
	$(RM) $(OBJECT_CYCLE)
	$(RM) $(OBJECT_CYCLE_POP)

