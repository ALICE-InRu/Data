CC=g++

INCLUDES= -I.

LIBRARIES= -lm

#CCFLAGS= -g

CCFLAGS= -O3 -DNDEBUG

all: fspgen

fspgen: fspgen.cc
	$(CC) -o fspgen $(INCLUDES) $(CCFLAGS) $(LIBRARIES) fspgen.cc

cleanemacs:
	rm -f *~ \#*

clean: cleanemacs
	rm -rf *.o fspgen core
