# Makefile for test_trace (your executable)

CXX = g++
CXXFLAGS = -I./include
LDFLAGS = -L./lib -lmcl -lmclbn256
TARGET = test_trace

SRCS = main.cpp tracing.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(OBJS) $(TARGET)