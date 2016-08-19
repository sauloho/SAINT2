all:
	cd ..; make

tt: transform.cpp transform.h point.o stream_printf.o 
	gcc transform.cpp point.o stream_printf.o -DTEST -o tt -lstdc++ -lm -pedantic -Wall

g: geom.cpp point.h matrix.h
	gcc geom.cpp point.o matrix.o -DTEST -o g -lstdc++ -lm -pedantic -Wall

d: distribution.cpp distribution.h
	gcc -DSELF_TEST distribution.cpp random.cpp -o d -Wall -pedantic -lstdc++

t: t.cpp point.h geom.h
	gcc -o t t.cpp point.o geom.o -Wall -pedantic -lstdc++ -lm

