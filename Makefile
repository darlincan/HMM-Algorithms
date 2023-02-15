hmm:main_hmm.c
	gcc hmm.c main.c -o hmm

run_hmm:hmm
	./hmm A.txt B.txt pi.txt

clean:
	rm hmm