EXEC = g++
FLAGS =

qubit.o: qubit.cpp qubit.h
	$(EXEC) $(FLAGS) -c qubit.cpp

q_circuit.o: q_circuit.cpp q_circuit.h
	$(EXEC) $(FLAGS) -c q_circuit.cpp

test.o: test.cpp q_circuit.h
	$(EXEC) $(FLAGS) -c test.cpp

test: qubit.o q_circuit.o test.o  
	$(EXEC) $(FLAGS) -o test qubit.o q_circuit.o test.o 

clean:
	rm test *.o *.gch

