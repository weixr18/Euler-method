CC:= g++
INCLUDE:= -I./include -IE:/OPEN/libs/gmp-6.2.1/gmp/include -IE:/OPEN/MinGW/mingw64/mingw-w64-gcc-10.2/include
LIBS:= -lstdc++ -lgmpxx -lgmp -LE:/OPEN/libs/gmp-6.2.1/gmp/lib
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
