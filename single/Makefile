CXX = g++
CXXFLAGS = -O3 -Wall
LIBS = -lm -lgmp -lflint -lntl -lpthread
TARGET = lfr

$(TARGET):
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(TARGET).cpp $(LIBS)

clean:
	rm -rf $(TARGET)
