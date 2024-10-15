# Имя исполняемого файла
all: a.exe

# Правило для сборки исполняемого файла
a.exe: main.o solve.o law.o
	g++  -o a.exe main.o solve.o law.o

# Правило для создания объектного файла main.o
main.o: main.cpp solve.hpp law.hpp
	g++  -c main.cpp

# Правило для создания объектного файла solve.o
solve.o: solve.cpp solve.hpp
	g++  -c solve.cpp

law.o: law.cpp law.hpp
	g++ -c law.cpp

# Правило для очистки сгенерированных файлов
clean:
	rm -f *.o a.exe
