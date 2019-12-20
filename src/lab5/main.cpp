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
 * @param lmin - нижняя оценка для нижнего собственного значения
 * @param lmax - верхняя оценка для верхнего собственного значения
 */
void init(int &size, double *&matrix, double *&f, double &eps, double &lmin, double &lmax);

/**
 * Функция последовательного метода Chebyshev
 * @param size - размер
 * @param matrix - матрица A
 * @param f - вектор f
 * @param om - параметр ω
 * @return
 */
double *serial_chebyshev(const int size, const double *matrix, const double *f, double eps, double &lmin, double &lmax);

/**
 * Функция параллельного метода Chebyshev
 * @param size - размер
 * @param matrix - матрица A
 * @param f - вектор f
 * @param om - параметр ω
 * @return
 */
double *openmp_chebyshev(const int size, const double *matrix, const double *f, double eps, double &lmin, double &lmax);

/**
 * Процедура для печати вектора
 * @param v - вектор
 * @param size - размер вектора
 */
void print_vector(const double *v, int size);

int main(int argc, char **argv) {
    int size;
    double eps, wtime, lmin, lmax;
    double *matrix, *f;

    init(size, matrix, f, eps, lmin, lmax);

    wtime = omp_get_wtime();
    auto x = serial_chebyshev(size, matrix, f, eps, lmin, lmax);
    print_vector(x, size);
    printf("%lf", omp_get_wtime() - wtime);

    omp_set_num_threads(2);
    wtime = omp_get_wtime();
    auto x2 = openmp_chebyshev(size, matrix, f, eps, lmin, lmax);
    print_vector(x2, size);
    printf("%lf", omp_get_wtime() - wtime);

    omp_set_num_threads(4);
    wtime = omp_get_wtime();
    auto x3 = openmp_chebyshev(size, matrix, f, eps, lmin, lmax);
    print_vector(x3, size);
    printf("%lf", omp_get_wtime() - wtime);

    omp_set_num_threads(8);
    wtime = omp_get_wtime();
    auto x4 = openmp_chebyshev(size, matrix, f, eps, lmin, lmax);
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

void init(int &size, double *&matrix, double *&f, double &eps, double &lmin, double &lmax) {
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

    printf("\nВведите нижнию оценку для нижнего собственного значения: ");
    scanf("%lf", &lmin);

    printf("\nВведите верхнию оценку для верхнего собственного значения: ");
    scanf("%lf", &lmax);

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

double *serial_chebyshev(
        const int size,
        const double *matrix,
        const double *f,
        double eps,
        double &lmin,
        double &lmax) {

    double d, c, err, sum = 0.0, alpha = 0.0, beta = 0.0;
    int i, j, step = 0;
    double *x, *r, *z, *p;

    d = (lmax + lmin) / 2;
    c = (lmax - lmin) / 2;

    x = new double[size];
    std::fill_n(x, size, 0);

    // r=f-matrix*x;
    r = new double[size];
    for (i = 0; i < size; i++) {
        r[i] = f[i];
    }

    z = new double[size];
    p = new double[size];

    do {
        step++;

        // z = linsolve(E,r);
        for (i = 0; i < size; i++) {
            z[i] = r[i];
        }

        if (step == 1) {
            for (i = 0; i < size; i++) {
                p[i] = z[i];
            }
            alpha = 1.0 / d;
        } else if (step == 2) {
            beta = 0.5 * std::pow(c * alpha, 2);
            alpha = 1.0 / (d - beta / alpha);
            for (i = 0; i < size; i++) {
                p[i] = z[i] + beta * p[i];
            }
        } else {
            beta = std::pow(0.5 * c * alpha, 2);
            alpha = 1.0 / (d - beta / alpha);
            for (i = 0; i < size; i++) {
                p[i] = z[i] + beta * p[i];
            }
        }

        // x=x+alpha*p;
        for (i = 0; i < size; i++) {
            x[i] = x[i] + alpha * p[i];
        }

        // r=f-matrix*x;
        for (i = 0; i < size; i++) {
            sum = 0.0;
            for (j = 0; j < size; j++) {
                sum += matrix[size * i + j] * x[j];
            }

            r[i] = f[i] - sum;
        }

        // norm(r)
        err = 0.0;
        for (i = 0; i < size; i++) {
            err += z[i] * z[i];
        }
    } while (std::sqrt(err) > eps && step < N_MAX);

    delete[] r;
    delete[] z;
    delete[] p;

    return x;
}

double *openmp_chebyshev(
        const int size,
        const double *matrix,
        const double *f,
        double eps,
        double &lmin,
        double &lmax) {

    double d, c, err, sum = 0.0, alpha = 0.0, beta = 0.0;
    int i, j, step = 0;
    double *x, *r, *z, *p;

    d = (lmax + lmin) / 2;
    c = (lmax - lmin) / 2;

    x = new double[size];
    std::fill_n(x, size, 0);

    // r=f-matrix*x;
    r = new double[size];
    for (i = 0; i < size; i++) {
        r[i] = f[i];
    }

    z = new double[size];
    p = new double[size];

    do {
        step++;

        // z = linsolve(E,r);
#pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            z[i] = r[i];
        }

        if (step == 1) {
#pragma omp parallel for private(i)
            for (i = 0; i < size; i++) {
                p[i] = z[i];
            }
            alpha = 1.0 / d;
        } else if (step == 2) {
            beta = 0.5 * std::pow(c * alpha, 2);
            alpha = 1.0 / (d - beta / alpha);
#pragma omp parallel for private(i)
            for (i = 0; i < size; i++) {
                p[i] = z[i] + beta * p[i];
            }
        } else {
            beta = std::pow(0.5 * c * alpha, 2);
            alpha = 1.0 / (d - beta / alpha);
#pragma omp parallel for private(i)
            for (i = 0; i < size; i++) {
                p[i] = z[i] + beta * p[i];
            }
        }

        // x=x+alpha*p;
#pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            x[i] = x[i] + alpha * p[i];
        }

        // r=f-matrix*x;
#pragma omp parallel for private(i, j)
        for (i = 0; i < size; i++) {
            sum = 0.0;
#pragma omp parallel for private(j) reduction(+:sum)
            for (j = 0; j < size; j++) {
                sum += matrix[size * i + j] * x[j];
            }

            r[i] = f[i] - sum;
        }

        // norm(r)
        err = 0.0;
#pragma omp parallel for private(i) reduction(+:err)
        for (i = 0; i < size; i++) {
            err += z[i] * z[i];
        }
    } while (std::sqrt(err) > eps && step < N_MAX);

    delete[] r;
    delete[] z;
    delete[] p;

    return x;
}
