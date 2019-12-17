#include <stdio.h>  /* printf, scanf, puts, NULL */
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <omp.h>
#include <math.h>

/**
 * Инициализация и считывания данных из stdin
 * @param v - сгенерированный вектор
 * @param min - минимальное значение элемента из stdin
 * @param max - максимальное значение элемента из stdin
 * @param size - размерность вектора из stdin
 */
void init(double *&v, int &min, int &max, int &size);

/**
 * Процедура для инициализирования вектора результата
 * @param r - вектор результата
 * @param size - размерность вектора
 */
void init_result_vector(double *&r, int size);

/**
 * Функция для получения рандомного вектора заданной длины
 * @param min - минимальное значение элемента
 * @param max - максимальное значение элемента
 * @param size - длина вектора
 * @return сгенерированный вектор
 */
double *get_rand_vector(int min, int max, int size);

/**
 * Процедура для печати вектора
 * @param v - вектор
 * @param size - размер вектора
 */
void print_vector(const double *v, int size);

int main() {
    double *v; // исходный вектор
    int min;   // минимальное значение элемента
    int max;   // максимальное значение элемента
    int size;  // размерность матрицы и векторов
    clock_t t; // переменная времени для подсчета времени вычисления

    // инициализация начальных данных
    init(v, min, max, size);
    print_vector(v, size);

    return 0;
}

void init(double *&v, int &min, int &max, int &size) {
    // получение размерности вектора
    do {
        printf("\nВведите размер вектора: ");
        scanf("%d", &size);
        if (size <= 0) {
            printf("Размер должен быть больше 0!\n");
        }
    } while (size <= 0);

    // ввод минимума
    printf("Введите минимальное значение элемента: ");
    scanf("%d", &min);

    // ввод максимума
    do {
        printf("Введите максимальное значение элемента: ");
        scanf("%d", &max);
        if (max <= min) {
            printf("Максимум должен быть больше %d!\n", min);
        }
    } while (max <= min);

    // инициализация вектора рандомными данными
    v = get_rand_vector(min, max, size);
}

void init_result_vector(double *&r, int size) {
    r = new double[size]();
}

double *get_rand_vector(int min, int max, int size) {
    double *v = new double[size];

    for (int i = 0; i < size; i++) {
        v[i] = rand() % (max - min + 1) + min;
    }

    return v;
}

void print_vector(const double *v, int size) {
    printf("\nvector = |");
    for (int i = 0; i < size - 1; i++) {
        printf("%.2f, ", v[i]);
    }
    printf("%.2f|\n", v[size - 1]);
}