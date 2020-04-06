fmm : main.c
	gcc -O4 -o fmm main.c -lm  -I .
.PHONY: run
run : fmm
	./fmm
	rm -f *.o
