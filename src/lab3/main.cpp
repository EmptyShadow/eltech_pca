#include <omp.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

#define N_MAX 1000

/**
 * Инициализация данных
 * @param size - размер квадратной матрицы
 * @param om - параметр 0 < ω <= 2
 * @param matrix - система уравнений - симметрическая положительно определенная матрица A
 * @param f - вектор f
 * @param eps - погрешность
 */
void init(int &size, double &om, double *&matrix, double *&f, double &eps);

/**
 * Процедура разности
 * @param dr - результат
 * @param v1 - первый вектор
 * @param v2 - второй вектор
 * @param size - размер векторов
 * @return - результат
 */
void diff(double *dr, const double *v1, const double *v2, int size);

/**
 * Преобразование входных данных
 * @param a - А
 * @param f - F
 * @param size - размерность
 * @param om - параметр ω
 * @param c - -A[i,j] * om / A[i,i]
 * @param d - F[i] * om / A[i,i]
 */
void transform_input_data(const double *a, const double *f, const int size, const double om, double *c, double *d);

/**
 * Функция последовательного метода верхней релаксации
 * @param size - размер
 * @param matrix - матрица A
 * @param f - вектор f
 * @param om - параметр ω
 * @return
 */
double *serial_sor(const int size, const double *matrix, const double *f, const double om, double eps);

/**
 * Функция параллельного метода верхней релаксации
 * @param size - размер
 * @param matrix - матрица A
 * @param f - вектор f
 * @param om - параметр ω
 * @return
 */
double *openmp_sor(const int size, const double *matrix, const double *f, const double om, double eps);

/**
 * Норма вектора
 * @param v
 * @param size
 * @return
 */
double norma(const double *v, int size);

/**
 * Процедура для печати вектора
 * @param v - вектор
 * @param size - размер вектора
 */
void print_vector(const double *v, int size);

int main(int argc, char **argv) {
    int size;
    double om, eps, wtime;
    double *matrix, *f;

    init(size, om, matrix, f, eps);

    wtime = omp_get_wtime();
    auto x = serial_sor(size, matrix, f, om, eps);
    print_vector(x, size);
    printf("%lf", omp_get_wtime() - wtime);

    wtime = omp_get_wtime();
    auto x2 = openmp_sor(size, matrix, f, om, eps);
    print_vector(x2, size);
    printf("%lf", omp_get_wtime() - wtime);

    delete x;
    delete x2;
    delete matrix;
    delete f;

    return 0;
}

void init(int &size, double &om, double *&matrix, double *&f, double &eps) {
    do {
        printf("\nВведите размер матрицы: ");
        scanf("%d", &size);
        if (size <= 0) {
            printf("Размер должен быть больше 0!\n");
        }
    } while (size <= 0);

    do {
        printf("\nВведите параметр 0 < ω <= 2: ");
        scanf("%lf", &om);
        if (om < 0.0 || om > 2.0) {
            printf("0 < ω <= 2!\n");
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

void diff(double *dr, const double *v1, const double *v2, int size) {
    for (int i = 0; i < size; i++) {
        dr[i] = v1[i] - v2[i];
    }
}

double norma(const double *v, int size) {
    double summ = 0.0;

    for (int i = 0; i < size; i++) {
        summ += v[i] * v[i];
    }

    return std::sqrt(summ);
}

void print_vector(const double *v, int size) {
    printf("\nvector = |");
    for (int i = 0; i < size - 1; i++) {
        printf("%lf, ", v[i]);
    }
    printf("%lf|\n", v[size - 1]);
}

void transform_input_data(const double *a, const double *f, const int size, const double om, double *c, double *d) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            c[size * i + j] = -a[size * i + j] / a[size * i + i];
        }

        d[i] = f[i] / a[size * i + i];
    }
}

double *serial_sor(const int size, const double *matrix, const double *f, const double om, double eps) {
    int step = 0;
    int i, j;
    double s1, s2;
    double *c, *d;
    double *x_prev, *x_next, *diffA, *tmp;

    x_prev = new double[size];
    x_next = new double[size];
    diffA = new double[size];
    c = new double[size * size];
    d = new double[size];

    std::fill_n(x_prev, size, 0);
    std::fill_n(x_next, size, 0);
    transform_input_data(matrix, f, size, om, c, d);

    do {
        step++;

        tmp = x_next;
        x_next = x_prev;
        x_prev = tmp;

        for (i = 0; i < size; i++) {
            s1 = 0;
            s2 = 0;

            for (j = 0; j < i; j++) {
                s1 = s1 + c[size * i + j] * om * x_next[j];
            }

            for (j = i + 1; j < size; j++) {
                s2 = s2 + c[size * i + j] * om * x_prev[j];
            }

            x_next[i] = s1 + s2 + d[i] * om - x_prev[i] * (om - 1);
        }

        diff(diffA, x_prev, x_next, size);
    } while (norma(diffA, size) > eps && step < N_MAX);


    delete x_prev;
    delete diffA;
    delete c;
    delete d;

    return x_next;
}

double *openmp_sor(const int size, const double *matrix, const double *f, const double om, double eps) {
    int step = 0;
    int i, j;
    double s1, s2;
    double *c, *d;
    double *x_prev, *x_next, *diffA, *tmp;

    x_prev = new double[size];
    x_next = new double[size];
    diffA = new double[size];
    c = new double[size * size];
    d = new double[size];

    std::fill_n(x_prev, size, 0);
    std::fill_n(x_next, size, 0);
    transform_input_data(matrix, f, size, om, c, d);

    do {
        step++;

        tmp = x_next;
        x_next = x_prev;
        x_prev = tmp;

        {
#pragma omp parallel for private(i, j, s1, s2) shared(size, om, x_prev, x_next, c, d) schedule(dynamic)
            for (i = 0; i < size; i++) {
                s1 = 0.0;
                s2 = 0.0;

                for (j = 0; j < i; j++) {
                    double q;
#pragma omp atomic read
                    q = x_next[j];
                    s1 = s1 + c[size * i + j] * om * q;
                }

                for (j = i + 1; j < size; j++) {
                    double q;
#pragma omp atomic read
                    q = x_prev[j];
                    s2 = s2 + c[size * i + j] * om * q;
                }

#pragma omp atomic write
                x_next[i] = s1 + s2 + d[i] * om - x_prev[i] * (om - 1);
            }
        }

        diff(diffA, x_prev, x_next, size);
    } while (norma(diffA, size) > eps && step < N_MAX);

    delete x_prev;
    delete diffA;
    delete c;
    delete d;

    return x_next;
}
