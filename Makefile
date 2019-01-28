CC = g++

# compiler flags:
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -Wall -Wno-c++11-extensions

# the build target executable:
TARGET = matching_0.2

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CC) $(CFLAGS) -o matching $(TARGET).cc

clean:
	$(RM) $(TARGET)
