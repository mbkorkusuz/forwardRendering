SRC_DIR = lib/
CXX = g++
CXXFLAGS = -O3 -std=c++11
TARGET = rasterizer

.PHONY: all rasterizer clean

all: $(TARGET)

$(TARGET): $(wildcard $(SRC_DIR)*.cpp)
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	rm -f $(TARGET)
