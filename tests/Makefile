CC = gcc
CFLAGS = -ggdb -fsanitize=address
RM = rm -f
ALL_TESTS = test_sort test_maes

all: $(ALL_TESTS)

test_sort: test_sort.c
	$(CC) $(CFLAGS) test_sort.c ../libtinymaes.a -lm -o test_sort -I..

test_maes: test_maes.c
	$(CC) $(CFLAGS) test_maes.c ../libtinymaes.a -lm -o test_maes -I..

clean:
	$(RM) $(ALL_TESTS)

