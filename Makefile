export OMP_NUM_THREADS:=8

build:
	@clang++ -o bin/${name} -fopenmp ${path}
lab1:
	@make name=lab1 path=src/lab1/main.cpp build
	@./bin/lab1