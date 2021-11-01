CC = g++
CFLAGS = -std=c++11

#test: main.o hung.o
#	$(CC) -o test main.o hung.o

#hung.o: Hungarian.cpp Hungarian.h
#	$(CC) -c Hungarian.cpp -o hung.o
	
#main.o: testMain.cpp Hungarian.h
#	$(CC) $(CFLAGS) -c testMain.cpp -o main.o

#clean:
#	-rm main.o hung.o

test: main.o hung.o
	$(CC) -o test main.o hung.o

hung.o: hungarian.cpp hungarian.h
	$(CC) -c hungarian.cpp -o hung.o -I/usr/include/eigen3/
	
main.o: test.cpp hungarian.h
	$(CC) $(CFLAGS) -c test.cpp -o main.o -I/usr/include/eigen3/

clean:
	-rm main.o hung.o
