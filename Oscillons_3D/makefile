EXEC = Oscillon3D
#
CC = g++
#
C++FLAGS = -Wno-deprecated -fopenmp
#
LIBS = -lm -fopenmp
#
#
C++FILES = Oscillon3D.o

$(EXEC): $(C++FILES)
	@echo "    "
	@echo "Linking..."
	@echo "----------"
	$(CC) -o $(EXEC) $(C++FILES) $(LIBS)

.C.o:
	@echo "   "
	@echo "Compiling Source File --" $<
	@echo "------------------------"
	$(CC) $(C++FLAGS) -c $<

.cpp.o:
	@echo "   "
	@echo "Compiling Numerical Recipe File --" $<
	@echo "----------------------------------"
	$(CC) $(C++FLAGS) -c $<

# dependencies:
Oscillon3D.o: Read_Input.h Oscillon.h dumper.h gridfunction.h nr3.h
