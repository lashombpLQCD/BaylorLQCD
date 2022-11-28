CC=gcc
CFLAGS= -I -Wall -g
DEPS = convert_from_milc_spm.c latparams.h
OBJFILES = convert_from_milc_spm.o
TARGET = sub_convert_from_milc_spm 

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) 
$(OBJFILES): $(DEPS)
	$(CC) $(CFLAGS) -c convert_from_milc_spm.c 
clean:
	rm -f $(OBJFILES) $(TARGET) *~
