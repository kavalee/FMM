tester : fmm.c util.c testing.c
	gcc -O4 -o tester testing.c fmm.c util.c -lm  -I .
.PHONY: run
test : tester
	./tester
	rm -f *.o
