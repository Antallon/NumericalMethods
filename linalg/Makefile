# Имя исполняемого файла
all: prog

# Правило для сборки исполняемого файла
prog: main.o solve.o
	g++  -o prog main.o solve.o

# Правило для создания объектного файла main.o
main.o: main.cpp solve.hpp
	g++  -O3 -c main.cpp

# Правило для создания объектного файла solve.o
solve.o: solve.cpp solve.hpp
	g++  -O3 -c solve.cpp

# Правило для очистки сгенерированных файлов
clean:
	rm -f *.o prog