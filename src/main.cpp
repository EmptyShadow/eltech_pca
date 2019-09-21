#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <omp.h>
#include <math.h>

// Функция для инициализации переменных и считывания данных из stdin
// m - сюда будет помещена матрица
// v - сюда будет помещен вектор
// min - минимальное значение элемента из stdin
// max - максимальное значение элемента из stdin
// s - размерность m и v из stdin
void init(double* &m, double* &v, int &min, int &max, int &size);

// Функция для инициализирования вектора результата
// r - вектор результата
// size - размерность вектора
void init_result(double* &r, int size);

// Функция для получения рандомного вектора заданной длины
// min - минимальное значение элемента
// max - максимальное значение элемента
// size - длина вектора
double* get_rand_vector(int min, int max, int size);

// Функция для печати вектора
// v - вектор
// size - размер вектора
void print_vector(const double* v, int size);

// Функция для печати матрицы
// m - матрица
// size - размер матрицы
void print_matrix(const double* m, int size);

// Функция для отображения времени
// t - время
void print_clock(clock_t t);

// Функция для перемножения матрицы на вектор последовательно
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void serial_matrix_multiply_vector(double* m, double* v, double* r, int size);

// Функция для перемножения матрицы на вектор параллельно по строкам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_rows_multiply_vector(double* m, double* v, double* r, int size);

// Функция для перемножения матрицы на вектор параллельно по столбцам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_cols_multiply_vector(double* m, double* v, double* r, int size);

// Функция для перемножения матрицы на вектор параллельно по блокам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_blocks1_multiply_vector(double* m, double* v, double* r, int size);

// Функция для перемножения матрицы на вектор параллельно по блокам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_blocks2_multiply_vector(double* m, double* v, double* r, int size);

// Функция разности двух векторов
double* v_diff_v(double* v1, double* v2, int size);

// Функция нормы вектора
double v_norma(double* v, int size);

int main() {
    double* m; // матрица
    double* v; // вектор
    double *r1, *r2, *r3, *r4, *r5;
    int min;   // минимальное значение элемента
    int max;   // максимальное значение элемента
    int size;  // размерность матрицы и векторов
    clock_t t; // переменная времени для подсчета времени вычисления

    // инициализация начальных данных
    init(m, v, min, max, size);

    /** Последовательно **/
    init_result(r1, size);
    t = clock();
    serial_matrix_multiply_vector(m, v, r1, size);
    t = clock() - t;
    printf("\nПоследовательное умножение строк матрицы на вектор");
    print_clock(t);
    /********************/

    /** Параллельно **/
    init_result(r2, size);
    t = clock();
    parallel_matrix_rows_multiply_vector(m, v, r2, size);
    t = clock() - t;
    printf("\nПараллельное умножение строк матрицы на вектор");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r1, r2, size), size));

    init_result(r3, size);
    t = clock();
    parallel_matrix_cols_multiply_vector(m, v, r3, size);
    t = clock() - t;
    printf("\nПараллельное умножение столбцов матрицы на вектор");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r1, r3, size), size));

    init_result(r4, size);
    t = clock();
    parallel_matrix_blocks1_multiply_vector(m, v, r4, size);
    t = clock() - t;
    printf("\nПараллельное умножение блоков матрицы на вектор без вложенного расспараллеливания");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r1, r4, size), size));

    init_result(r5, size);
    t = clock();
    parallel_matrix_blocks2_multiply_vector(m, v, r5, size);
    t = clock() - t;
    printf("\nПараллельное умножение блоков матрицы на вектор с вложенным расспараллеливанием в блоках");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r1, r5, size), size));
    /********************/

    return 0;
}

void init(double* &m, double* &v, int &min, int &max, int &size) {
    // получение размерности матрицы, вектора и результирующего вектора
    do {
        printf("\nВведите размер матрицы и вектора: ");
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

    // инициализация матрицы и вектора рандомными данными
    m = get_rand_vector(min, max, size*size);
    v = get_rand_vector(min, max, size);
}

void init_result(double* &r, int size) {
    r = new double [size]();
}

double* get_rand_vector(int min, int max, int size) {
    double* v = new double[size];
    srand(time(NULL)/2);

    for (int i = 0; i < size; i++) {
        v[i] = rand() % (max - min + 1) + min;
    }

    return v;
}

void print_vector(const double* v, int size) {
    printf("\nvector = |");
    for (int i = 0; i < size - 1; i++) {
        printf("%.2f, ", v[i]);
    }
    printf("%.2f|\n", v[size - 1]);
}

void print_matrix(const double* m, int size) {
    int i = 0, j = 0;
    printf("\nmatrix = ");
    for (; i < size; i++) {
        int index = i * size;
        printf("\n|");
        for (j = 0 ; j < size - 1; j++) {
            printf("%.2f, ", m[index + j]);
        }
        printf("%.2f|\n", m[index + j]);
    }
}

void print_clock(clock_t t) {
    printf("\nмикросекунды: %ld, секунды: %f\n", t, ((float)t)/CLOCKS_PER_SEC);
}

void serial_matrix_multiply_vector(double* m, double* v, double* r, int size) {
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            r[i] += m[i * size + j] * v[j];
        }
    }
}

void parallel_matrix_rows_multiply_vector(double* m, double* v, double* r, int size) {
    int i, j;
    // parallel for - расспараллеливает первый цикл 
    // private (j) - говорит, что j должна копироваться для каждого потока
    #pragma omp parallel for private (j)
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            r[i] += m[i * size + j] * v[j];
        }
    }
}

void parallel_matrix_cols_multiply_vector(double* m, double* v, double* r, int size) {
    int i, j;
    omp_lock_t lock;

    // инициализация замка для потоков
    omp_init_lock(&lock);
    for (i = 0; i < size; i++) {
        // распараллеливание
        #pragma omp parallel for
        for (j = 0; j < size; j++) {
            // блокировка текущим потоком доступа другим потоком
            omp_set_lock(&lock);
            r[i] += m[i * size + j]*v[j];
            // разблокировка блокировки текущим потоком для других потоков
            omp_unset_lock(&lock);
        }
    }
    omp_destroy_lock(&lock);
}

void parallel_matrix_blocks1_multiply_vector(double* m, double* v, double* r, int size) {
    int threadID;
    int gridThreadsNum = 4;
    int gridSize = int(sqrt(double(gridThreadsNum)));
    // Предполагается, что размер матрицы кратен
    // размеру сетки потоков
    int blockSize = size / gridSize; // размер одного блока
    omp_set_num_threads(gridThreadsNum);

    #pragma omp parallel private (threadID)
    {
        threadID = omp_get_thread_num(); // получение номера потока

        // создание массива куда будет помещен результат по потоку
        double * pThreadResult;
        init_result(pThreadResult, size);

        // стартовая координата элемента в блоке текущего потока
        int i_start = (int(threadID / gridSize)) * blockSize;
        int j_start = (threadID % gridSize) * blockSize;

        double iterResult;
        // умножение строк блока на нужные элементы вектора
        for (int i = 0; i < blockSize; i++) {
            iterResult = 0;
            for (int j = 0; j < blockSize; j++) {
                iterResult += m[(i + i_start) * size + (j + j_start)] * v[j + j_start];
            }
            pThreadResult[i + i_start] = iterResult;
        }

        #pragma omp critical
        for (int i = 0; i < size; i++) {
            r[i] += pThreadResult[i];
        } // pragma omp parallel
        delete [] pThreadResult;
    }
}

void parallel_matrix_blocks2_multiply_vector(double* m, double* v, double* r, int size) {
    int nestedThreadsNum = 2; // количество потоков для умножения строки на вектор
    omp_set_num_threads(nestedThreadsNum);
    omp_set_nested(true);
    // расспараллелливания строк
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        double threadResult = 0.0;
        // параллельная сумма умножения строки матрицы на вектор при этом применяется
        // nestedThreadsNum количество вложенных потоков
        #pragma omp parallel for reduction(+:threadResult)
        for (int j = 0; j < size; j++) {
            threadResult += m[i * size + j] * v[j];
        }
        r[i] = threadResult;
    }
}

double* v_diff_v(double* v1, double* v2, int size) {
    double* r;
    init_result(r, size);

    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        r[i] = v1[i] - v2[i];
    }

    return r;
}

double v_norma(double* v, int size) {
    double summ = 0.0;
    #pragma omp parallel for reduction(+:summ)
    for (int i = 0; i < size; i++) {
        summ += v[i] * v[i];
    }

    return sqrt(summ);
}