CC = gcc
CFLAGS = -O2 -funroll-loops -Wall
#CFLAGS = -ggdb -fsanitize=address -Wall
RM = rm -f
ALL_TESTS = test_sort test_maes

all: $(ALL_TESTS)

test_sort: test_sort.c
	$(CC) $(CFLAGS) test_sort.c ../libtinymaes.a -lm -o test_sort -I..

test_maes: test_maes.c
	$(CC) $(CFLAGS) test_maes.c ../libtinymaes.a -lm -o test_maes -I..

clean:
	$(RM) $(ALL_TESTS)

