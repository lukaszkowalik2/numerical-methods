CXX = g++
CXXFLAGS = -std=c++11 -Wall

SRC1 = NUM1a.cpp
SRC2 = NUM1b.cpp

BUILD_DIR = build
EXEC1 = $(BUILD_DIR)/NUM1a
EXEC2 = $(BUILD_DIR)/NUM1b

all: $(EXEC1) $(EXEC2)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(EXEC1): $(SRC1) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(EXEC1) $(SRC1)

$(EXEC2): $(SRC2) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(EXEC2) $(SRC2)

run1: $(EXEC1)
	./$(EXEC1)

run2: $(EXEC2)
	./$(EXEC2)

clean:
	rm -rf $(BUILD_DIR)
