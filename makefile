CC:= g++
INCLUDE:= -I./include -I./include/gmp
LIBS:= -lstdc++ -lgmpxx -lgmp -L./lib
CXXFLAGS:= -std=c++11 -g
DIR_SRC:= ./src
DIR_OBJ:= ./obj
DIR_OBJ_WIN:= .\obj
TARGET:= main.exe
OBJECTS := main.o mpmat.o mpmat_utils.o euler.o

OBJECTS := $(addprefix $(DIR_OBJ)/,$(OBJECTS))

all: $(TARGET)

$(shell mkdir obj)

$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(LIBS)

$(DIR_OBJ)/%.o: $(DIR_SRC)/%.cpp   
	$(CC) -c $(CXXFLAGS) $(INCLUDE) $< -o $@
  
.PHONY : clean
clean:   
	-del $(DIR_OBJ_WIN)\*.o
	-del $(TARGET) 
	rmdir obj
