TARGET = NUM2

BUILD_DIR = build
SRC = NUM2.cpp

CXX = g++
CXXFLAGS = -I ./lib/eigen-3.4.0 -O2 -std=c++11

$(BUILD_DIR)/$(TARGET): $(SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/$(TARGET) $(SRC)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

all: $(BUILD_DIR)/$(TARGET)

run: $(BUILD_DIR)/$(TARGET)
	./$(BUILD_DIR)/$(TARGET)
