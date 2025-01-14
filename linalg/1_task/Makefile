# Имя исполняемого файла
all: prog

# Правило для сборки исполняемого файла
prog: main1.o solve1.o
	g++  -o prog main1.o solve1.o

# Правило для создания объектного файла main.o
main1.o: main1.cpp solve1.hpp
	g++  -O3 -c main1.cpp

# Правило для создания объектного файла solve.o
solve1.o: solve1.cpp solve1.hpp
	g++  -O3 -c solve1.cpp

# Правило для очистки сгенерированных файлов
clean:
	rm -f *.o prog