CXXFLAGS=-std=c++11 -g $(shell root-config --cflags) -fopenmp
LIBS=$(shell root-config --libs) -lMathMore -fopenmp

run : analyze_light
			@echo "Finished Compiling..."
			@echo "To run: ./analyze_light"

analyze_light : analyze_light.o data_output.o semi_analytic_hits.o time_parameterisation.o utility_functions.o

	g++ -o $@ $^ ${LIBS}

%.o : %.cc
	g++ ${CXXFLAGS} -o $@ -c $^
