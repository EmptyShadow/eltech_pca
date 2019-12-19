#include <omp.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

#define N_MAX 1000

/**
 * Инициализация данных
 * @param size - размер квадратной матрицы
 * @param matrix - система уравнений - симметрическая положительно определенная матрица A
 * @param f - вектор f
 * @param eps - погрешность
 */
void init(int &size, double *&matrix, double *&f, double &eps);

/**
 * Функция последовательного метода Gaus-seidel
 * @param size - размер
 * @param matrix - матрица A
 * @param f - вектор f
 * @param om - параметр ω
 * @return
 */
double *serial_gaus_seidel(const int size, const double *matrix, const double *f, double eps);

/**
 * Функция параллельного метода Gaus-seidel
 * @param size - размер
 * @param matrix - матрица A
 * @param f - вектор f
 * @param om - параметр ω
 * @return
 */
double *openmp_gaus_seidel(const int size, const double *matrix, const double *f, double eps);

/**
 * Процедура для печати вектора
 * @param v - вектор
 * @param size - размер вектора
 */
void print_vector(const double *v, int size);

int main(int argc, char **argv) {
    int size;
    double eps, wtime;
    double *matrix, *f;

    init(size, matrix, f, eps);

    wtime = omp_get_wtime();
    auto x = serial_gaus_seidel(size, matrix, f, eps);
    print_vector(x, size);
    printf("%lf", omp_get_wtime() - wtime);

    omp_set_num_threads(2);
    wtime = omp_get_wtime();
    auto x2 = openmp_gaus_seidel(size, matrix, f, eps);
    print_vector(x2, size);
    printf("%lf", omp_get_wtime() - wtime);

    omp_set_num_threads(4);
    wtime = omp_get_wtime();
    auto x3 = openmp_gaus_seidel(size, matrix, f, eps);
    print_vector(x3, size);
    printf("%lf", omp_get_wtime() - wtime);

    omp_set_num_threads(8);
    wtime = omp_get_wtime();
    auto x4 = openmp_gaus_seidel(size, matrix, f, eps);
    print_vector(x4, size);
    printf("%lf", omp_get_wtime() - wtime);

    delete x;
    delete x2;
    delete x3;
    delete x4;
    delete matrix;
    delete f;

    return 0;
}

void init(int &size, double *&matrix, double *&f, double &eps) {
    do {
        printf("\nВведите размер матрицы: ");
        scanf("%d", &size);
        if (size <= 0) {
            printf("Размер должен быть больше 0!\n");
        }
    } while (size <= 0);

    printf("\nAx = f\n");

    printf("\nВведите коэффициенты матрицы A по строчно. Матрица должна быть симметрической положительно определенной.\n");
    matrix = new double[size * size]();
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            scanf("%lf", &matrix[size * row + col]);
        }
    }

    printf("\nВведите вектор f: ");
    f = new double[size]();
    for (int i = 0; i < size; i++) {
        scanf("%lf", &f[i]);
    }

    printf("\nЗадайте погрешность Eps: ");
    scanf("%lf", &eps);
}

void print_vector(const double *v, int size) {
    printf("\nvector = |");
    for (int i = 0; i < size - 1; i++) {
        printf("%lf, ", v[i]);
    }
    printf("%lf|\n", v[size - 1]);
}

double *serial_gaus_seidel(const int size, const double *matrix, const double *f, double eps) {
    int step = 0;
    int i, j;
    double *x, *xPrev, y, xLast, sum, dx;

    x = new double[size];
    xPrev = new double[size];
    std::fill_n(x, size, 0);
    std::fill_n(xPrev, size, 0);

    do {
        sum = 0.0;

        for (i = 0; i < size; i++) {
            xPrev[i] = x[i];
        }

        for (i = 0; i < size; i++) {
            xLast = x[i];
            y = 0.0;

            for (j = 0; j < i; j++) {
                y += matrix[size * i + j] * x[j];
            }

            for (j = i + 1; j < size; j++) {
                y += matrix[size * i + j] * xPrev[j];
            }

            x[i] = (f[i] - y) / matrix[size * i + i];

            dx = x[i] - xLast;
            sum += dx * dx;
        }

        step++;
    } while (std::sqrt(sum) > eps && step < N_MAX);

    delete[] xPrev;

    return x;
}

double *openmp_gaus_seidel(const int size, const double *matrix, const double *f, double eps) {
    int step = 0;
    int i, j;
    double *x, *xPrev, y, sum, dx;

    x = new double[size];
    xPrev = new double[size];
    std::fill_n(x, size, 0);
    std::fill_n(xPrev, size, 0);

    do {
        sum = 0.0;

#pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            xPrev[i] = x[i];
        }

        for (i = 0; i < size; i++) {
            y = 0.0;

#pragma omp parallel for reduction(+:y) private(j)
            for (j = 0; j < i; j++) {
                y += matrix[size * i + j] * x[j];
            }

#pragma omp parallel for reduction(+:y) private(j)
            for (j = i + 1; j < size; j++) {
                y += matrix[size * i + j] * xPrev[j];
            }

            x[i] = (f[i] - y) / matrix[size * i + i];

            dx = x[i] - xPrev[i];
            sum += dx * dx;
        }

        step++;
    } while (std::sqrt(sum) > eps && step < N_MAX);

    delete[] xPrev;

    return x;
}
