CXX = g++
CXXFLAGS = -Wall -std=c++17 -I/usr/include/eigen3
TARGET = NUM7
SRC = NUM7.cpp
OBJ = $(SRC:.cpp=.o)
BUILD_DIR = build

all: $(BUILD_DIR)/$(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/$(TARGET): $(SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(SRC) -o $@

clean:
	rm -rf $(BUILD_DIR)

run: $(BUILD_DIR)/$(TARGET)
	./$(BUILD_DIR)/$(TARGET)

zip:
	zip -r $(TARGET).zip * -x "build/*" ".vscode/*"

.PHONY: all clean run zip