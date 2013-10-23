CC = g++
CFLAGS = 
LDLIBS =

build: C

t_consc.o: t_consc.h t_consc.cpp helpers.h

t_partition.o: t_partition.h t_partition.cpp helpers.h

t_subset.o: t_subset.h t_subset.cpp

PartitionEnumerator.o: PartitionEnumerator.h PartitionEnumerator.cpp

main.o: main.cpp

C: t_consc.o t_partition.o main.o PartitionEnumerator.o t_subset.o
	$(CC) -o C t_consc.o t_partition.o main.o

clean:
	rm -rf *~ *.o c C

