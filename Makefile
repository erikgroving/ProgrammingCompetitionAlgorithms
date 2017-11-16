SOURCE=main.c contest.c

all: contest

contest: $(SOURCE) contest.h
	gcc $(SOURCE) -Wall -lm -o contest

clean:
	rm contest
