CC = g++

# compiler flags:
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -Wall -Wno-c++11-extensions

# the build target executable:
TARGET = matching

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cc

clean:
	$(RM) $(TARGET)
