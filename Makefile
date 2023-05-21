SOURCES=Block.cpp ChebyShev.cpp EstimatePDFMain.cpp InputData.cpp InputParameters.cpp JointProbability.cpp MinimizeScore.cpp OutputControl.cpp Score.cpp ScoreQZ.cpp StitchPDF.cpp Variable.cpp WriteResults.cpp
OBJ := $(SOURCES:%.cpp=%.o)
CFLAGS=-fopenmp -g -DPDF_BINARY_1D_READ=1 -Wall
LDFLAGS=-fopenmp 

%.o: %.cpp
	mpic++ -c $(CFLAGS) $< -o $@ 

all: $(OBJ)
	mpic++ $(LDFLAGS) $(OBJ) -o estimate 

clean:
	rm $(OBJ)
	rm estimate