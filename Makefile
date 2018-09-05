#BASE=$(shell pwd)
BASE=/Users/alifarnudi/Documents/Dr\ Ejtehadi/Programming/Latest\ Tiam/Ali/Membrane_OBJ
CC=g++

#CFLAGS= -std=c++11 -O2 

INCLUDE= -I$(BASE)/include

SRC=$(BASE)/source

HDRS_BASE=$(BASE)/include

EXEFILE = $(BASE)/bin/Cell

HDRS_FILES=General_constants.h\
Membrane.h

#18
HDRSPP_FILES=ECM.hpp\
General_functions.hpp\
write_functions.hpp\
interaction.hpp


MAIN_SRC=main.cpp \
General_functions.cpp\
Average_Membrane_Node_Distance.cpp\
check.cpp\
ConstantSurfaceForceLocalTriangles.cpp\
Elastic_Force_Calculator.cpp\
Membrane_bend_force.cpp\
Membrane_MD_Evolution.cpp\
Membrane_node_neighbour_list.cpp\
Membrane_node_pair_identifier.cpp\
Membrane_normal_direction_identifier.cpp\
Membrane_num_of_Node_Pair_Counter_2.cpp\
Membrane_num_of_Node_Pair_Counter.cpp\
Membrane_read_functions.cpp\
Membrane_spring_force.cpp\
Membrane_Triangle_Pair_and_Edges_Identifier.cpp\
Membrane_triangle_pair_counter.cpp\
Membrane_triangle_pair_identifier.cpp\
write_functions.cpp\
ECM_Node_Pair_Identifier.cpp\
ECM_Normal_direction_Identifier.cpp\
ECM_read_functions.cpp\
ECM.cpp\
Interaction_functions.cpp\
interaction_neighbour_list_updaters.cpp\
interaction.cpp\
neighbour_pool_constructor.cpp\
triangle_interaction_neighbour_list_updaters.cpp\

HDRS=$(HDRS_FILES:%.h=$(HDRS_BASE)/%.h)
HDRSPP=$(HDRSPP_FILES:%.hpp=$(HDRS_BASE)/%.hpp)

OBJS=$(MAIN_SRC:.cpp=.o)

MAIN_CMP=$(MAIN_SRC:%.cpp=$(SRC)/%.cpp)

all:main
main:$(OBJS)
	@echo "Linking ...."
	$(CC) -o $(EXEFILE) $(OBJS) 
$(OBJS): $(HDRS HDRSPP)
	$(CC) $(INCLUDE) -c $(MAIN_CMP)

clean:
	@rm -f $(OBJS)