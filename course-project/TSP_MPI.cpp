#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N_COLONY 100
#define CITY     442

int     xColony = 100;
int     xCity = CITY;
double  probab1 = 0.02;
long    maxGen = 200000;

int     colony[N_COLONY * 2][CITY];
double  city_dis[CITY][CITY];
double  dis_p[N_COLONY * 2];
int     temp[CITY];
double  sumbest;
long    GenNum;
double  timeStart;

int     position(int *tmp, int C);
void    invert(int pos_start, int pos_end);
void    printBest(long GenNum, int rank, double elapsed);
void    select1();

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    timeStart = MPI_Wtime();
    srand((unsigned)time(NULL) + rank);

    /* ---------- rank 0 reads file and broadcasts ---------- */
    double cityXY[CITY][2];
    if (rank == 0)
    {
        FILE *fp;
        if ((fp = fopen("pcb442.tsp", "r")) == NULL)
        {
            printf("Rank 0: cannot open pcb442.tsp\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(fp, "%d", &xCity);
        for (int i = 0; i < xCity; i++)
            fscanf(fp, "%*d%lf%lf", &cityXY[i][0], &cityXY[i][1]);
        fclose(fp);
        printf("init: read %d cities\n", xCity);
    }

    MPI_Bcast(&xCity, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cityXY, CITY * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* ---------- all compute distance matrix ---------- */
    for (int i = 0; i < xCity; i++)
        for (int j = 0; j < xCity; j++)
        {
            if (j > i)
            {
                double d = (cityXY[i][0] - cityXY[j][0]) * (cityXY[i][0] - cityXY[j][0])
                         + (cityXY[i][1] - cityXY[j][1]) * (cityXY[i][1] - cityXY[j][1]);
                city_dis[i][j] = (int)(sqrt(d) + 0.5);
            }
            else if (j == i)
                city_dis[i][j] = 0;
            else
                city_dis[i][j] = city_dis[j][i];
        }

    /* ---------- rank 0 initialises colony ---------- */
    if (rank == 0)
    {
        int array[CITY];
        for (int i = 0; i < xCity; i++)
            array[i] = i;

        for (int i = 0; i < xColony; i++)
        {
            int mod = xCity;
            for (int j = 0; j < xCity; j++)
            {
                int sign = rand() % mod;
                colony[i][j] = array[sign];
                int t = array[mod - 1];
                array[mod - 1] = array[sign];
                array[sign] = t;
                mod--;
                if (mod == 1) { j++; colony[i][j] = array[0]; break; }
            }
        }

        for (int i = 0; i < xColony; i++)
        {
            dis_p[i] = 0;
            for (int j = 0; j < xCity - 1; j++)
                dis_p[i] += city_dis[colony[i][j]][colony[i][j + 1]];
            dis_p[i] += city_dis[colony[i][0]][colony[i][xCity - 1]];
        }

        sumbest = dis_p[0];
        for (int i = 1; i < xColony; i++)
            if (dis_p[i] < sumbest) sumbest = dis_p[i];

        printf("init success!!!  initial best = %f\n", sumbest);
    }

    /* ---------- divide tours among processes ---------- */
    int chunk = xColony / size;
    int start = rank * chunk;
    int end = (rank == size - 1) ? xColony : start + chunk;
    int local_count = end - start;

    /* prepare Gatherv arrays on root */
    int *recv_counts = NULL, *displs = NULL;
    int *recv_counts_d = NULL, *displs_d = NULL;
    if (rank == 0)
    {
        recv_counts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
        recv_counts_d = (int *)malloc(size * sizeof(int));
        displs_d = (int *)malloc(size * sizeof(int));
        for (int i = 0; i < size; i++)
        {
            int c = (i == size - 1) ? (xColony - i * chunk) : chunk;
            recv_counts[i] = c * CITY;
            displs[i] = i * chunk * CITY;
            recv_counts_d[i] = c;
            displs_d[i] = i * chunk;
        }
    }

    /* ---------- generation loop ---------- */
    for (GenNum = 0; GenNum < maxGen; GenNum++)
    {
        MPI_Bcast(colony, N_COLONY * CITY, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(dis_p, N_COLONY, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /* each process mutates its chunk */
        for (int i = start; i < end; i++)
        {
            for (int j = 0; j < xCity; j++)
                temp[j] = colony[i][j];

            double disChange = 0;
            int pos_flag = 0;
            int pos_C = rand() % xCity;

            for (;;)
            {
                int pos_C1, C1, j, k1, k2, l1, l2;

                if ((rand() / 32768.0) < probab1)
                {
                    do pos_C1 = rand() % xCity;
                    while (pos_C1 == pos_C);
                    C1 = colony[i][pos_C1];
                }
                else
                {
                    do j = rand() % xColony;
                    while (j == i);
                    int k = position(colony[j], temp[pos_C]);
                    C1 = colony[j][(k + 1) % xCity];
                    pos_C1 = position(temp, C1);
                }

                if ((pos_C + 1) % xCity == pos_C1 || (pos_C - 1 + xCity) % xCity == pos_C1)
                    break;

                k1 = temp[pos_C];
                k2 = temp[(pos_C + 1) % xCity];
                l1 = temp[pos_C1];
                l2 = temp[(pos_C1 + 1) % xCity];
                disChange += city_dis[k1][l1] + city_dis[k2][l2]
                           - city_dis[k1][k2] - city_dis[l1][l2];

                invert(pos_C, pos_C1);
                pos_flag++;
                if (pos_flag > xCity - 1)
                    break;
                pos_C++;
                if (pos_C >= xCity) pos_C = 0;
            }

            dis_p[N_COLONY + i] = dis_p[i] + disChange;
            for (int j = 0; j < xCity; j++)
                colony[N_COLONY + i][j] = temp[j];
        }

        /* gather children to root */
        MPI_Gatherv(&colony[N_COLONY + start][0], local_count * CITY, MPI_INT,
                    &colony[N_COLONY][0], recv_counts, displs, MPI_INT,
                    0, MPI_COMM_WORLD);
        MPI_Gatherv(&dis_p[N_COLONY + start], local_count, MPI_DOUBLE,
                    &dis_p[N_COLONY], recv_counts_d, displs_d, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        /* root: select + print */
        if (rank == 0)
        {
            select1();
            sumbest = dis_p[0];
            for (int j = 1; j < xColony; j++)
                if (dis_p[j] < sumbest) sumbest = dis_p[j];

            printf("%ld:%f\n", GenNum + 1, sumbest);

            if ((GenNum + 1) % 2000 == 0 && GenNum + 1 < maxGen)
            {
                double elapsed = MPI_Wtime() - timeStart;
                printBest(GenNum + 1, rank, elapsed);
            }
        }
    }

    /* final output */
    if (rank == 0)
    {
        double elapsed = MPI_Wtime() - timeStart;
        printf("Finial solution:");
        printBest(maxGen, rank, elapsed);
        printf("  distance = %f, time = %.2f s\n", sumbest, elapsed);
    }

    free(recv_counts); free(displs);
    free(recv_counts_d); free(displs_d);
    MPI_Finalize();
    return 0;
}

/* ---------- helper functions (identical to serial) ---------- */

void select1()
{
    for (int j = 0; j < N_COLONY; j++)
        if (dis_p[N_COLONY + j] < dis_p[j])
        {
            dis_p[j] = dis_p[N_COLONY + j];
            for (int k = 0; k < CITY; k++)
                colony[j][k] = colony[N_COLONY + j][k];
        }
}

void invert(int pos_start, int pos_end)
{
    int j, k, t;
    if (pos_start < pos_end)
    {
        j = pos_start + 1;
        k = pos_end;
        for (; j <= k; j++, k--)
        { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
    }
    else
    {
        if (xCity - 1 - pos_start <= pos_end + 1)
        {
            j = pos_end; k = pos_start + 1;
            for (; k < xCity; j--, k++)
            { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
            k = 0;
            for (; k <= j; k++, j--)
            { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
        }
        else
        {
            j = pos_end; k = pos_start + 1;
            for (; j >= 0; j--, k++)
            { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
            j = xCity - 1;
            for (; k <= j; k++, j--)
            { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
        }
    }
}

int position(int *tmp, int C)
{
    for (int j = 0; j < xCity; j++)
        if (tmp[j] == C) return j;
    return -1;
}

void printBest(long gen, int rank, double elapsed)
{
    FILE *fp = fopen("tsp0.txt", "a");
    if (!fp) return;
    fprintf(fp, "%ld\t%.2f\t%d\n", gen, elapsed, (int)sumbest);
    fclose(fp);
}
