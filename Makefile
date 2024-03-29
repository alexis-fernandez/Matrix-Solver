CC = gcc
CFLAGS = -Wall -Wextra -std=c99
SOURCES = functions.c main.c
EXECUTABLE = main
LIBRARIES = -lm

all: $(EXECUTABLE)

$(EXECUTABLE):$(SOURCES)
	$(CC) $(SOURCES) -o $(EXECUTABLE) $(LIBRARIES)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE)
