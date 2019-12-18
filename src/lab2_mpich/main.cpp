#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <cstdlib>
#include <algorithm>

#define N 100000
#define MIN -32
#define MAX 99

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
 * Функция слияния двух массивов
 * @param v1 - первый
 * @param n1 - размер первого
 * @param v2 - второй
 * @param n2 - размер второго
 * @return v1.v2
 */
double *merge(double *v1, int n1, double *v2, int n2);

/**
 * Процедура сортировки вектора методом Шелла
 * @param a - исходный вектор, который будет отсортирован
 * @param size - размер массива
 */
void shellsort(double *a, int size);

int main(int argc, char **argv) {
    double *src; // исходный вектор
    double *chunk; // частички вектора, коотрые обрабатываются параллельно
    double *other; // требуется для слияния всех частичек
    int m, min = MIN, max = MAX, sizeSrc = N;
    int processId, countProcesses; // идентификатор текущего процесса, количество процессов в кластере
    int sizeChunk; // размер частички вектора
    int i;
    int step;
    MPI_Status status; // mpich атрибуты сообщения
    int mpiError; // mpich ошибка
    double wtime;

    /**
     * MPI_COMM_WORLD - зарезервированный идентификатор группы, состоящей их всех процессов приложения
     */

    // инициализация mpich библиотеки
    if ((mpiError = MPI_Init(&argc, &argv)) != MPI_SUCCESS) {
        printf("MPI_Init: error = %d", mpiError);
        return mpiError;
    }

    // определение номера процесса в группе
    if ((mpiError = MPI_Comm_rank(MPI_COMM_WORLD, &processId)) != MPI_SUCCESS) {
        printf("MPI_Comm_rank: error = %d", mpiError);
        return mpiError;
    }

    // определение числа процессов в группе
    if ((mpiError = MPI_Comm_size(MPI_COMM_WORLD, &countProcesses)) != MPI_SUCCESS) {
        printf("MPI_Comm_size: error = %d", mpiError);
        return mpiError;
    }

    if (processId == 0) {
        // задание для root процесса
        printf("id: %d, количество процессов: %d", processId, countProcesses);

        int r;
        sizeChunk = sizeSrc / countProcesses;
        r = sizeSrc % countProcesses;

        // инициализация вектора рандомными данными
        src = (double *) malloc((sizeSrc + countProcesses - r) * sizeof(double));
        for (i = 0; i < sizeSrc; i++) {
            src[i] = rand() % (max - min + 1) + min;
        }

        // обнуление хвоста если нужно
        if (r != 0) {
            for (i = sizeSrc; i < sizeSrc + countProcesses - r; i++) {
                src[i] = 0;
            }

            sizeChunk++;
        }

//        print_vector(src, sizeSrc + countProcesses - r);
        wtime = MPI_Wtime();

        // Отправка широковещательных сообщений с размером вектора, который придет следующим сообщением
        MPI_Bcast(&sizeChunk, 1, MPI_INT, 0, MPI_COMM_WORLD);

        init_result_vector(chunk, sizeChunk);

        // Распределения частичек данных по всем процессам группы включая root
        MPI_Scatter(src, sizeChunk, MPI_DOUBLE, chunk, sizeChunk, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // сортировка частички, которая досталась root
        shellsort(chunk, sizeChunk);
    } else {
        // задание для workers

        MPI_Bcast(&sizeChunk, 1, MPI_INT, 0, MPI_COMM_WORLD);
        init_result_vector(chunk, sizeChunk);
        MPI_Scatter(src, sizeChunk, MPI_DOUBLE, chunk, sizeChunk, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        shellsort(chunk, sizeChunk);
    }

    step = 1;
    while (step < countProcesses) {
        if (processId % (2 * step) == 0) {
            if (processId + step < countProcesses) {
                MPI_Recv(&m, 1, MPI_INT, processId + step, 0, MPI_COMM_WORLD, &status);
                init_result_vector(other, m);
                MPI_Recv(other, m, MPI_DOUBLE, processId + step, 0, MPI_COMM_WORLD, &status);
                chunk = merge(chunk, sizeChunk, other, m);
                sizeChunk = sizeChunk + m;
            }
        } else {
            int near = processId - step;
            MPI_Send(&sizeChunk, 1, MPI_INT, near, 0, MPI_COMM_WORLD);
            MPI_Send(chunk, sizeChunk, MPI_DOUBLE, near, 0, MPI_COMM_WORLD);
            break;
        }
        step = step * 2;
    }

    if (processId == 0) {
//        print_vector(chunk, sizeSrc);
        wtime = MPI_Wtime() - wtime;
        print_clock(wtime);
        check_sort(chunk, sizeSrc);
    }
    MPI_Finalize();
}

void init_result_vector(double *&r, int size) {
    r = (double *) malloc(size * sizeof(double));
}

double *copy_vector(const double *v, int size) {
    double *a;
    init_result_vector(a, size);
    std::copy(v, v + size, a);
    return a;
}

double *get_rand_vector(int min, int max, int size) {
    auto *v = (double *) malloc(size * sizeof(double));

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
    printf("\nВремя исполнения: %.16g\n", t);
}

void check_sort(const double *v, int size) {
    for (int i = 1; i < size; i++) {
        if (v[i - 1] > v[i]) {
            printf("\nСортировка плохая\n");
            return;
        };
    }

    printf("\nСортировка верная\n");
}

double *merge(double *v1, int n1, double *v2, int n2) {
    int i, j, k;
    double *result;

    init_result_vector(result, n1 + n2);

    i = 0;
    j = 0;
    k = 0;
    while (i < n1 && j < n2) {
        if (v1[i] < v2[j]) {
            result[k] = v1[i];
            i++;
            k++;
        } else {
            result[k] = v2[j];
            j++;
            k++;
        }
    }

    if (i == n1) {
        while (j < n2) {
            result[k] = v2[j];
            j++;
            k++;
        }
    } else {
        while (i < n1) {
            result[k] = v1[i];
            i++;
            k++;
        }
    }

    return result;
}

void shellsort(double *a, int size) {
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
}