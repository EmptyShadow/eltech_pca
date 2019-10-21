#include <stdio.h>  /* printf, scanf, puts, NULL */
#include <stdlib.h> /* srand, rand */
#include <omp.h>
#include <math.h>

// Функция для инициализации переменных и считывания данных из stdin
// m1 - сюда будет помещена матрица 1
// m2 - сюда будет помещена матрица 2
// v - сюда будет помещен вектор
// min - минимальное значение элемента из stdin
// max - максимальное значение элемента из stdin
// s - размерность m1, m2 и v из stdin
void init(double *&m1, double *&m2, double *&v, int &min, int &max, int &size);

// Функция для инициализирования вектора результата
// r - вектор результата
// size - размерность вектора
void init_result_vector(double *&r, int size);

// Функция для получения рандомного вектора заданной длины
// min - минимальное значение элемента
// max - максимальное значение элемента
// size - длина вектора
double *get_rand_vector(int min, int max, int size);

// Функция для печати вектора
// v - вектор
// size - размер вектора
void print_vector(const double *v, int size);

// Функция для печати матрицы
// m - матрица
// size - размер матрицы
void print_matrix(const double *m, int size);

// Функция для отображения времени
// t - время
void print_clock(double t);

// Функция для перемножения матрицы на вектор последовательно
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void serial_matrix_multiply_vector(double *m, double *v, double *r, int size);

// Функция для перемножения матрицы на матрицу последовательно
// m1 - матрица 1
// m2 - матрица 2
// r - матрица результата
// size - размерность матриц
void serial_matrix_multiply_matrix(double *m1, double *m2, double *r, int size);

// Функция для перемножения матрицы на вектор параллельно по строкам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_rows_multiply_vector(double *m, double *v, double *r, int size);

// Функция для перемножения матрицы на матрицу параллельно по строкам матрицы
// m1 - матрица 1
// m2 - матрица 2
// r - вектор результата
// size - размерность матриц
void parallel_matrix_rows_multiply_matrix(double *m1, double *m2, double *r, int size);

// Функция для перемножения матрицы на вектор параллельно по столбцам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_cols_multiply_vector(double *m, double *v, double *r, int size);

// Функция для перемножения матрицы на матрицу параллельно строки со столбцами
// m1 - матрица 1
// m2 - матрица 2
// r - вектор результата
// size - размерность матриц
void parallel_matrix_rows_multiply_matrix_cols(double *m1, double *m2, double *r, int size);

// Функция для перемножения матрицы на вектор параллельно по блокам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_blocks1_multiply_vector(double *m, double *v, double *r, int size);

// Функция для перемножения матрицы на вектор параллельно по блокам матрицы
// m - матрица
// v - вектор
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_blocks2_multiply_vector(double *m, double *v, double *r, int size);

// Функция для перемножения матрицы на матрицу параллельно по блокам матрицы
// m1 - матрица 1
// m2 - матрица 2
// r - вектор результата
// size - размерность матрицы и векторов
void parallel_matrix_blocks_multiply_matrix(double *m1, double *m2, double *r, int size);

// Функция разности двух векторов
double *v_diff_v(double *v1, double *v2, int size);

// Функция нормы вектора
double v_norma(double *v, int size);

int main()
{
    double *m1, *m2; // матрицы
    double *v;       // вектор
    double *r11, *r21, *r31, *r41, *r51;
    double *r12, *r22, *r32, *r42;
    int min;  // минимальное значение элемента
    int max;  // максимальное значение элемента
    int size; // размерность матрицы и векторов
    double t; // переменная времени для подсчета времени вычисления

    // инициализация начальных данных
    init(m1, m2, v, min, max, size);
    // print_vector(v, size);
    // print_matrix(m1, size);
    // print_matrix(m2, size);

    /** Последовательно **/
    init_result_vector(r11, size);
    t = omp_get_wtime();
    serial_matrix_multiply_vector(m1, v, r11, size);
    t = omp_get_wtime() - t;
    printf("\nПоследовательное умножение строк матрицы на вектор");
    print_clock(t);
    // print_vector(r11, size);

    init_result_vector(r12, size * size);
    t = omp_get_wtime();
    serial_matrix_multiply_matrix(m1, m2, r12, size);
    t = omp_get_wtime() - t;
    printf("\nПоследовательное умножение матрицы на матрицу");
    print_clock(t);
    // print_matrix(r12, size);
    /********************/

    printf("Максимальное количество потоков: %d\n", omp_get_max_threads());
    omp_set_num_threads(omp_get_max_threads());

    /** Параллельно **/
    init_result_vector(r21, size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_vector(m1, v, r21, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на вектор");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r21, size), size));

    init_result_vector(r31, size);
    t = omp_get_wtime();
    parallel_matrix_cols_multiply_vector(m1, v, r31, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение столбцов матрицы на вектор");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r31, size), size));

    init_result_vector(r41, size);
    t = omp_get_wtime();
    parallel_matrix_blocks1_multiply_vector(m1, v, r41, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение блоков матрицы на вектор без вложенного расспараллеливания");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r41, size), size));

    init_result_vector(r51, size);
    t = omp_get_wtime();
    parallel_matrix_blocks2_multiply_vector(m1, v, r51, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение блоков матрицы на вектор с вложенным расспараллеливанием в блоках");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r51, size), size));

    init_result_vector(r22, size * size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_matrix(m1, m2, r22, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на матрицу");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r12, r22, size * size), size * size));

    init_result_vector(r32, size * size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_matrix_cols(m1, m2, r32, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на матрицу ленточное разделение данных строк и параллельно столбцы");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r12, r32, size * size), size * size));

    init_result_vector(r42, size * size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_matrix_cols(m1, m2, r42, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на матрицу по блокам");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r12, r42, size * size), size * size));
    /********************/

    printf("Максимальное количество потоков: %d\n", 1);
    omp_set_num_threads(1);

    /** Параллельно **/
    init_result_vector(r21, size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_vector(m1, v, r21, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на вектор");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r21, size), size));

    init_result_vector(r31, size);
    t = omp_get_wtime();
    parallel_matrix_cols_multiply_vector(m1, v, r31, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение столбцов матрицы на вектор");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r31, size), size));

    init_result_vector(r41, size);
    t = omp_get_wtime();
    parallel_matrix_blocks1_multiply_vector(m1, v, r41, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение блоков матрицы на вектор без вложенного расспараллеливания");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r41, size), size));

    init_result_vector(r51, size);
    t = omp_get_wtime();
    parallel_matrix_blocks2_multiply_vector(m1, v, r51, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение блоков матрицы на вектор с вложенным расспараллеливанием в блоках");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r11, r51, size), size));

    init_result_vector(r22, size * size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_matrix(m1, m2, r22, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на матрицу");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r12, r22, size * size), size * size));

    init_result_vector(r32, size * size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_matrix_cols(m1, m2, r32, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на матрицу ленточное разделение данных строк и параллельно столбцы");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r12, r32, size * size), size * size));

    init_result_vector(r42, size * size);
    t = omp_get_wtime();
    parallel_matrix_rows_multiply_matrix_cols(m1, m2, r42, size);
    t = omp_get_wtime() - t;
    printf("\nПараллельное умножение строк матрицы на матрицу по блокам");
    print_clock(t);
    printf("Ошибка: %f\n", v_norma(v_diff_v(r12, r42, size * size), size * size));
    /********************/

    return 0;
}

void init(double *&m1, double *&m2, double *&v, int &min, int &max, int &size)
{
    // получение размерности матрицы, вектора и результирующего вектора
    do
    {
        printf("\nВведите размер матрицы и вектора: ");
        scanf("%d", &size);
        if (size <= 0)
        {
            printf("Размер должен быть больше 0!\n");
        }
    } while (size <= 0);

    // ввод минимума
    printf("Введите минимальное значение элемента: ");
    scanf("%d", &min);

    // ввод максимума
    do
    {
        printf("Введите максимальное значение элемента: ");
        scanf("%d", &max);
        if (max <= min)
        {
            printf("Максимум должен быть больше %d!\n", min);
        }
    } while (max <= min);

    // инициализация матрицы и вектора рандомными данными
    m1 = get_rand_vector(min, max, size * size);
    m2 = get_rand_vector(min, max, size * size);
    v = get_rand_vector(min, max, size);
}

void init_result_vector(double *&r, int size)
{
    r = new double[size]();
}

double *get_rand_vector(int min, int max, int size)
{
    double *v = new double[size];
    // srand(time(timer_t)/2);

    for (int i = 0; i < size; i++)
    {
        v[i] = rand() % (max - min + 1) + min;
    }

    return v;
}

void print_vector(const double *v, int size)
{
    printf("\nvector = |");
    for (int i = 0; i < size - 1; i++)
    {
        printf("%.2f, ", v[i]);
    }
    printf("%.2f|\n", v[size - 1]);
}

void print_matrix(const double *m, int size)
{
    int i = 0, j = 0;
    printf("\nmatrix = ");
    for (; i < size; i++)
    {
        int index = i * size;
        printf("\n|");
        for (j = 0; j < size - 1; j++)
        {
            printf("%.2f, ", m[index + j]);
        }
        printf("%.2f|\n", m[index + j]);
    }
}

void print_clock(double t)
{
    printf("\nсекунды: %f\n", t);
}

void serial_matrix_multiply_vector(double *m, double *v, double *r, int size)
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            r[i] += m[i * size + j] * v[j];
        }
    }
}

void serial_matrix_multiply_matrix(double *m1, double *m2, double *r, int size)
{
    // индекс строк m1, индекс столбцов m2, индекс элемента в строке m1 и столбце m2
    int i, j, k;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            for (k = 0; k < size; k++)
            {
                r[i * size + j] += m1[i * size + k] * m2[k * size + j];
            }
        }
    }
}

void parallel_matrix_rows_multiply_vector(double *m, double *v, double *r, int size)
{
    int i, j;
// parallel for - расспараллеливает первый цикл
// private (j) - говорит, что j должна копироваться для каждого потока
#pragma omp parallel for shared(m, v, r, size) private(i, j)
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            r[i] += m[i * size + j] * v[j];
        }
    }
}

void parallel_matrix_rows_multiply_matrix(double *m1, double *m2, double *r, int size)
{
    // индекс строк m1, индекс столбцов m2, индекс элемента в строке m1 и столбце m2
    int i, j, k;
#pragma omp parallel for shared(m1, m2, r, size) private(i, j, k)
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            for (k = 0; k < size; k++)
            {
                r[i * size + j] += m1[i * size + k] * m2[k * size + j];
            }
        }
    }
}

void parallel_matrix_cols_multiply_vector(double *m, double *v, double *r, int size)
{
    int i, j;
    omp_lock_t lock;

    // инициализация замка для потоков
    omp_init_lock(&lock);
    for (i = 0; i < size; i++)
    {
// распараллеливание
#pragma omp parallel for
        for (j = 0; j < size; j++)
        {
            // блокировка текущим потоком доступа другим потоком
            omp_set_lock(&lock);
            r[i] += m[i * size + j] * v[j];
            // разблокировка блокировки текущим потоком для других потоков
            omp_unset_lock(&lock);
        }
    }
    omp_destroy_lock(&lock);
}

void parallel_matrix_rows_multiply_matrix_cols(double *m1, double *m2, double *r, int size)
{
    // индекс строк m1, индекс столбцов m2, индекс элемента в строке m1 и столбце m2
    int i, j, k;
    int nestedThreadsNum = 2; // количество потоков для умножения строки на столбцы
    omp_set_nested(true);
    omp_set_num_threads(nestedThreadsNum);
#pragma omp parallel for private(j, k)
    for (i = 0; i < size; i++)
    {
#pragma omp parallel for private(k)
        for (j = 0; j < size; j++)
        {
            for (k = 0; k < size; k++)
            {
                r[i * size + j] += m1[i * size + k] * m2[k * size + j];
            }
        }
    }
}

void parallel_matrix_blocks1_multiply_vector(double *m, double *v, double *r, int size)
{
    int threadID;
    int gridThreadsNum = 4;
    int gridSize = int(sqrt(double(gridThreadsNum)));
    // Предполагается, что размер матрицы кратен
    // размеру сетки потоков
    int blockSize = size / gridSize; // размер одного блока
    omp_set_num_threads(gridThreadsNum);

#pragma omp parallel private(threadID)
    {
        threadID = omp_get_thread_num(); // получение номера потока

        // создание массива куда будет помещен результат по потоку
        double *pThreadResult;
        init_result_vector(pThreadResult, size);

        // стартовая координата элемента в блоке текущего потока
        int i_start = (int(threadID / gridSize)) * blockSize;
        int j_start = (threadID % gridSize) * blockSize;

        double iterResult;
        // умножение строк блока на нужные элементы вектора
        for (int i = 0; i < blockSize; i++)
        {
            iterResult = 0;
            for (int j = 0; j < blockSize; j++)
            {
                iterResult += m[(i + i_start) * size + (j + j_start)] * v[j + j_start];
            }
            pThreadResult[i + i_start] = iterResult;
        }

#pragma omp critical
        for (int i = 0; i < size; i++)
        {
            r[i] += pThreadResult[i];
        } // pragma omp parallel
        delete[] pThreadResult;
    }
}

void parallel_matrix_blocks2_multiply_vector(double *m, double *v, double *r, int size)
{
    int nestedThreadsNum = 2; // количество потоков для умножения строки на вектор
    omp_set_num_threads(nestedThreadsNum);
    omp_set_nested(true);
// расспараллелливания строк
#pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        double threadResult = 0.0;
// параллельная сумма умножения строки матрицы на вектор при этом применяется
// nestedThreadsNum количество вложенных потоков
#pragma omp parallel for reduction(+ \
                                   : threadResult)
        for (int j = 0; j < size; j++)
        {
            threadResult += m[i * size + j] * v[j];
        }
        r[i] = threadResult;
    }
}

void parallel_matrix_blocks_multiply_matrix(double *m1, double *m2, double *r, int size)
{
    int gridThreadsNum = 4;
    int gridSize = int(sqrt(double(gridThreadsNum)));
    // Предполагается, что размер матрицы кратен
    // размеру сетки потоков
    int blockSize = size / gridSize; // размер одного блока
    omp_set_num_threads(gridThreadsNum);

#pragma omp parallel
    {
        int threadID = omp_get_thread_num();
        int rowIndex = threadID / gridSize;
        int colIndex = threadID % gridSize;
        for (int iter = 0; iter < gridSize; iter++)
        {
            for (int i = rowIndex * blockSize; i < (rowIndex + 1) * blockSize; i++)
            {
                for (int j = colIndex * blockSize; j < (colIndex + 1) * blockSize; j++)
                {
                    for (int k = iter * blockSize; k < (iter + 1) * blockSize; k++)
                        r[i * size + j] += m1[i * size + k] * m2[k * size + j];
                }
            }
        }
    } // pragma omp parallel
}

double *v_diff_v(double *v1, double *v2, int size)
{
    double *r;
    init_result_vector(r, size);

#pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        r[i] = v1[i] - v2[i];
    }

    return r;
}

double v_norma(double *v, int size)
{
    double summ = 0.0;
#pragma omp parallel for reduction(+ \
                                   : summ)
    for (int i = 0; i < size; i++)
    {
        summ += v[i] * v[i];
    }

    return sqrt(summ);
}