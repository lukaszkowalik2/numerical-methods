CXX = g++
CXXFLAGS = -Wall -O2 -I./lib/eigen-3.4.0
TARGET = NUM4
SRC = NUM4.cpp
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
	zip -r projekt.zip * -x "build/*" ".vscode/*"

.PHONY: all clean run zip
