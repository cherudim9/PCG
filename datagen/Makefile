all: datagen

debug: datagen.cpp
	g++ -g -Wall -std=c++0x datagen.cpp -o datagen

datagen: datagen.cpp
	g++ -O3 -Wall -std=c++0x datagen.cpp -o datagen

clean:
	rm -f datagen

ibmpg1: datagen
	./datagen ibmpg1.spice ibmpg1.solution ibmpg1.data

ibmpg2: datagen
	./datagen ibmpg2.spice ibmpg2.solution ibmpg2.data

ibmpg3: datagen
	./datagen ibmpg3.spice ibmpg3.solution ibmpg3.data

ibmpg4: datagen
	./datagen ibmpg4.spice ibmpg4.solution ibmpg4.data

ibmpg5: datagen
	./datagen ibmpg5.spice ibmpg5.solution ibmpg5.data

ibmpg6: datagen
	./datagen ibmpg6.spice ibmpg6.solution ibmpg6.data

sample1: datagen
	./datagen sample1.spice sample1.solution sample1.data