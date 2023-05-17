SOURCES=Block.cpp ChebyShev.cpp EstimatePDFMain.cpp InputData.cpp InputParameters.cpp JointProbability.cpp MinimizeScore.cpp OutputControl.cpp Score.cpp ScoreQZ.cpp StitchPDF.cpp Variable.cpp WriteResults.cpp
OBJ := $(SOURCES:%.cpp=%.o)
CFLAGS=-fopenmp -ltbb -g 
LDFLAGS=-fopenmp -g

%.o: %.cpp
	g++ -c $(CFLAGS) $< -o $@ 

all: $(OBJ)
	g++ $(LDFLAGS) $(OBJ) -o estimate -ltbb

clean:
	rm $(OBJ)
	rm estimate