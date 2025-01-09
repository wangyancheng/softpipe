SYSCONF_LINK = gcc
CFLAGS     = -g
LDFLAGS      =
LIBS         = -lm

DESTDIR = ./
TARGET  = main

OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))

all: $(DESTDIR)$(TARGET)

$(DESTDIR)$(TARGET): $(OBJECTS)
	$(SYSCONF_LINK)  -W -Wall  $(LDFLAGS) -o $(DESTDIR)$(TARGET) $(OBJECTS) $(LIBS)

$(OBJECTS): %.o: %.c
	$(SYSCONF_LINK)  -W -Wall -c $(CFLAGS) $< -o $@

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)
	-rm -f *.tga
	-rm -f *.bmp
	-rm -f *.tile
