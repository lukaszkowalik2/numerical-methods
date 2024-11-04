TARGET = build/NUM3
TARGET_TIMED = build/NUM3_timed
CXX = g++
CXXFLAGS = -I ./lib/eigen-3.4.0 -O2 -std=c++11

BUILD_DIR = build
SRC = NUM3.cpp
SRC_TIMED = NUM3_timed.cpp
HEADER = algorithms.h
DATA_FILE = times.csv
PYTHON_SCRIPT = generate_chart.py

all: $(TARGET) $(TARGET_TIMED)

$(TARGET): $(SRC) $(HEADER) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

$(TARGET_TIMED): $(SRC_TIMED) $(HEADER) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(TARGET_TIMED) $(SRC_TIMED)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

run: $(TARGET)
	./$(TARGET)

generate_data: $(TARGET_TIMED)
	./$(TARGET_TIMED)

chart: generate_data $(PYTHON_SCRIPT)
	python3 $(PYTHON_SCRIPT)

clean:
	rm -rf $(BUILD_DIR) $(DATA_FILE)

.PHONY: all clean run generate_data chart
