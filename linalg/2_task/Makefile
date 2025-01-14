# Имя исполняемого файла
all: prog

# Правило для сборки исполняемого файла
prog: main2.o solve2.o
	g++  -o prog main2.o solve2.o

# Правило для создания объектного файла main.o
main2.o: main2.cpp solve2.hpp
	g++  -O3 -c main2.cpp

# Правило для создания объектного файла solve.o
solve2.o: solve2.cpp solve2.hpp
	g++  -O3 -c solve2.cpp

# Правило для очистки сгенерированных файлов
clean:
	rm -f *.o prog