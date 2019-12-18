#include <cstdio>  /* printf, scanf, puts, NULL */
#include <cstdlib> /* srand, rand */
#include <ctime>   /* time */
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <algorithm>    // std::copy

#define N 100000
#define MIN -32
#define MAX 99

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
 * Функция копирования вектора
 * @param v - исходный вектор
 * @param size - размер вектора
 * @return - копия
 */
double *copy_vector(const double *v, int size);

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

/**
 * Процедура для отображения времени
 * @param t - время
 */
void print_clock(double t);

/**
 * Процедура проверки правильности сортировки
 * @param v
 * @param size
 * @return
 */
void check_sort(const double *v, int size);

/**
 * Последовательная сортировка метом Шелла
 * @param v - исходный вектор
 * @param size - размер вектора
 * @return отсортированный вектор
 */
double *serial_sort_shell(const double *v, int size);

/**
 * Параллельная сортировка метом Шелла (openmp)
 * @param v - исходный вектор
 * @param size - размер вектора
 * @return отсортированный вектор
 */
double *openmp_sort_shell(const double *v, int size);

int main(int argc, char **argv) {
    double *v; // исходный вектор
    double *r; // результат
    int min = MIN;   // минимальное значение элемента
    int max = MAX;   // максимальное значение элемента
    int size = N;  // размерность матрицы и векторов
    double t; // переменная времени для подсчета времени вычисления

    // инициализация начальных данных
//    init(v, min, max, size);
    v = get_rand_vector(min, max, size);
//    omp_set_num_threads(4);
//    print_vector(v, size);

    // последовательно
    t = omp_get_wtime();
    r = serial_sort_shell(v, size);
    t = omp_get_wtime() - t;
    printf("\nПоследовательная сортировка методом Шелла");
    print_clock(t);
//    print_vector(r, size);
    check_sort(r, size);

    // параллельно openmp
    t = omp_get_wtime();
    r = openmp_sort_shell(v, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельная сортировка методом Шелла (openmp)");
    print_clock(t);
//    print_vector(r, size);
    check_sort(r, size);

    printf("\n");

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

double *copy_vector(const double *v, int size) {
    double *a;
    init_result_vector(a, size);
    std::copy(v, v + size, a);
    return a;
}

double *get_rand_vector(int min, int max, int size) {
    auto *v = new double[size];

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

void print_clock(double t) {
    printf("\nВремя исполнения: %.16g", t);
}

void check_sort(const double *v, int size) {
    for (int i = 1; i < size; i++) {
        if (v[i - 1] > v[i]) {
            printf("\nСортировка плохая\n");
            return;
        };
    }

    printf("\nСортировка верная");
}

double *serial_sort_shell(const double *v, int size) {
    double *a = copy_vector(v, size);

    //инициализируем шаг
    int step = size / 2;

    //пока шаг не 0
    while (step > 0) {
        for (int i = 0; i < (size - step); i++) {
            int j = i;

            // будем идти начиная с i-го элемента
            // пока не пришли к началу массива
            // и пока рассматриваемый элемент больше
            // чем элемент находящийся на расстоянии шага
            while (j >= 0 && a[j] > a[j + step]) {
                // меняем их местами
                double temp = a[j];
                a[j] = a[j + step];
                a[j + step] = temp;
                j--;
            }
        }

        //уменьшаем шаг
        step = step / 2;
    }

    return a;
}

double *openmp_sort_shell(const double *v, int size) {
    double *a = copy_vector(v, size);

#pragma omp parallel firstprivate(size)
    {
        //инициализируем шаг
        int step = size / 2;

        //пока шаг не 0
        while (step > 0) {
#pragma omp parallel for
            for (int i = 0; i < (size - step); i++)
#pragma omp critical
            {
                int j = i;

                // будем идти начиная с i-го элемента
                // пока не пришли к началу массива
                // и пока рассматриваемый элемент больше
                // чем элемент находящийся на расстоянии шага
                while (j >= 0 && a[j] > a[j + step]) {
                    // меняем их местами
                    double temp = a[j];
                    a[j] = a[j + step];
                    a[j + step] = temp;
                    j--;
                }
            }

            //уменьшаем шаг
            step = step / 2;
        }
    }

    return a;
}
