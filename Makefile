SOURCE=main.c contest.c

all: contest

contest: $(SOURCE)
	gcc $(SOURCE) -Wall -lm -o contest

clean:
	rm contest
