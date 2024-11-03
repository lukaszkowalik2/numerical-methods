TARGET = build/NUM3

CXX = g++
CXXFLAGS = -I ./lib/eigen-3.4.0 -O2 -std=c++11

BUILD_DIR = build
SRC = NUM3.cpp

all: $(TARGET)

$(TARGET): $(SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean run
