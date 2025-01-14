# Имя исполняемого файла
all: prog

# Правило для сборки исполняемого файла
prog: main3.o solve3.o
	g++  -o prog main3.o solve3.o

# Правило для создания объектного файла main.o
main3.o: main3.cpp solve3.hpp
	g++  -O3 -c main3.cpp

# Правило для создания объектного файла solve.o
solve3.o: solve3.cpp solve3.hpp
	g++  -O3 -c solve3.cpp

# Правило для очистки сгенерированных файлов
clean:
	rm -f *.o prog