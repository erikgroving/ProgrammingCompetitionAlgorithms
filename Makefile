SOURCE=main.c contest.c

all: contest

contest: $(SOURCE) contest.h
	gcc $(SOURCE) -g -Wall -lm -o contest

clean:
	rm contest
