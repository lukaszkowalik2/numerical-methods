CXX = g++
CXXFLAGS = -I /usr/include/eigen3 -O2 -std=c++11
TARGET = NUM5
SRC = NUM5.cpp
OBJ = $(SRC:.cpp=.o)
BUILD_DIR = build

all: $(BUILD_DIR)/$(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/$(TARGET): $(SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(BUILD_DIR)

run: $(BUILD_DIR)/$(TARGET)
	./$(BUILD_DIR)/$(TARGET)

zip:
	zip -r $(TARGET).zip * -x "build/*" ".vscode/*"

.PHONY: all clean run zip
