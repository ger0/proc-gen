LIBS=-lGL -lglfw -lGLEW
HEADERS=*.hpp *.h
FILES=*.cpp

main_file: $(FILES) $(HEADERS)
	g++ -o worms $(FILES)  $(LIBS) -I.