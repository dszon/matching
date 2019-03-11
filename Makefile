CC = g++

# compiler flags:
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -Wall -Wno-c++11-extensions -std=c++11

# the build target executable:
TARGET = matching_0.3

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CC) $(CFLAGS) -o matching $(TARGET).cc

clean:
	$(RM) $(TARGET)
