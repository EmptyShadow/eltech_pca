export OMP_NUM_THREADS:=1

build:
	@clang++ -o bin/pca -v -fopenmp src/main.cpp
run:
	@make build
	@./bin/pca