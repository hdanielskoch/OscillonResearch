EXEC = main
#
CC = g++
#
C++FLAGS = -Wno-deprecated -fopenmp
#
LIBS = -lm -fopenmp
#
#
C++FILES = main.o Wave_Pred.o

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
main.o: Wave.h
Wave_Pred.o: deviates.h gamma.h ran.h Wave.h nr3.h
