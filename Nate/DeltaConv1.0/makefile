CC=g++
CFLAGS=-O1 -larmadillo

all: clean
	#$(CC) $(CFLAGS) -L/storage/group/mnh5174_collab/lab_resources/usr/lib -I/storage/group/mnh5174_collab/lab_resources/usr/include -o deltacon DeltaCon.cpp
	$(CC) $(CFLAGS) -o deltacon_mod DeltaCon.cpp
demo: all
	./deltacon_mod filelist 1 naive 0.1
clean:
	rm -f deltacon_mod
