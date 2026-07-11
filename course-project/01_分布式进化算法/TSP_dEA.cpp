// Algorithm 1: dEA (Distributed Evolutionary Algorithm)
// Master-worker: rank 0 maintains population, workers evolve chunks
// Compile: cl /EHsc /I"%MSMPI_INC%" TSP_dEA.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
// Run: mpiexec -n 4 TSP_dEA.exe [runs] [maxGen]

#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define CITY     442
#define N_COLONY 100

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

#include "../TSP_base.h"

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int numRuns = 30;
    if (argc > 1) numRuns = atoi(argv[1]);
    if (argc > 2) maxGen = atol(argv[2]);

    timeStart = MPI_Wtime();
    srand((unsigned)time(NULL) + rank);

    // ---- read city data on rank 0 ----
    double cityXY[CITY][2];
    if (rank == 0) {
        FILE *fp;
        if ((fp = fopen("pcb442.tsp", "r")) == NULL) {
            printf("Rank 0: cannot open pcb442.tsp\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(fp, "%d", &xCity);
        for (int i = 0; i < xCity; i++)
            fscanf(fp, "%*d%lf%lf", &cityXY[i][0], &cityXY[i][1]);
        fclose(fp);
    }

    MPI_Bcast(&xCity, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cityXY, CITY * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // ---- distance matrix ----
    for (int i = 0; i < xCity; i++)
        for (int j = 0; j < xCity; j++) {
            if (j > i) {
                double d = (cityXY[i][0] - cityXY[j][0]) * (cityXY[i][0] - cityXY[j][0])
                         + (cityXY[i][1] - cityXY[j][1]) * (cityXY[i][1] - cityXY[j][1]);
                city_dis[i][j] = (int)(sqrt(d) + 0.5);
            } else if (j == i) city_dis[i][j] = 0;
            else city_dis[i][j] = city_dis[j][i];
        }

    // ---- work division ----
    int chunk = xColony / size;
    int start = rank * chunk;
    int end = (rank == size - 1) ? xColony : start + chunk;
    int local_count = end - start;

    // ---- Gatherv arrays ----
    int *recv_counts = NULL, *displs = NULL;
    int *recv_counts_d = NULL, *displs_d = NULL;
    if (rank == 0) {
        recv_counts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
        recv_counts_d = (int *)malloc(size * sizeof(int));
        displs_d = (int *)malloc(size * sizeof(int));
        for (int i = 0; i < size; i++) {
            int c = (i == size - 1) ? (xColony - i * chunk) : chunk;
            recv_counts[i] = c * CITY;
            displs[i] = i * chunk * CITY;
            recv_counts_d[i] = c;
            displs_d[i] = i * chunk;
        }
    }

    // ---- results & convergence files ----
    FILE *resFile = NULL;
    FILE *convFile = NULL;
    if (rank == 0) {
        resFile = fopen("results_dEA.txt", "w");
        if (resFile) {
            fprintf(resFile, "# dEA  pcb442  pop=%d  maxGen=%ld  procs=%d\n", xColony, maxGen, size);
            fprintf(resFile, "# Run\\tBest\\tTime(s)\\n");
        }
        convFile = fopen("convergence_dEA.txt", "w");
        if (convFile)
            fprintf(convFile, "# Gen\\tTime\\tBest\\n");
    }

    // ---- multi-run ----
    for (int run = 0; run < numRuns; run++) {
        srand((unsigned)time(NULL) + rank + run * 137);

        // ---- init on master ----
        if (rank == 0) {
            init_population(colony, dis_p, xColony);
            sumbest = best_distance(dis_p, xColony);
            if (run == 0) printf("init success! initial best = %.0f\n", sumbest);
        }

        double runStart = MPI_Wtime();
        long prevCvGen = 0;
        double prevCvBest = 0;

        for (GenNum = 0; GenNum < maxGen; GenNum++) {
            MPI_Bcast(colony, N_COLONY * CITY, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(dis_p, N_COLONY, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            // ---- inver-over on local range ----
            for (int i = start; i < end; i++) {
                for (int j = 0; j < xCity; j++) temp[j] = colony[i][j];
                double disChange = 0;
                int pos_flag = 0;
                int pos_C = rand() % xCity;

                for (;;) {
                    int pos_C1, C1, j, k1, k2, l1, l2;
                    if ((rand() / 32768.0) < probab1) {
                        do pos_C1 = rand() % xCity;
                        while (pos_C1 == pos_C);
                        C1 = colony[i][pos_C1];
                    } else {
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
                    if (pos_flag > xCity - 1) break;
                    pos_C++;
                    if (pos_C >= xCity) pos_C = 0;
                }

                dis_p[N_COLONY + i] = dis_p[i] + disChange;
                for (int j = 0; j < xCity; j++)
                    colony[N_COLONY + i][j] = temp[j];
            }

            // ---- gather children to master ----
            if (rank == 0) {
                MPI_Gatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                            &colony[N_COLONY][0], recv_counts, displs, MPI_INT,
                            0, MPI_COMM_WORLD);
                MPI_Gatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                            &dis_p[N_COLONY], recv_counts_d, displs_d, MPI_DOUBLE,
                            0, MPI_COMM_WORLD);
            } else {
                MPI_Gatherv(&colony[N_COLONY + start][0], local_count * CITY, MPI_INT,
                            NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gatherv(&dis_p[N_COLONY + start], local_count, MPI_DOUBLE,
                            NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }

            // ---- master selects ----
            if (rank == 0) {
                select_elitist(colony, dis_p, xColony);
                sumbest = best_distance(dis_p, xColony);

                // critical velocity check
                if (GenNum % 2000 == 0) {
                    double currBest = sumbest;
                    double elapsed = MPI_Wtime() - timeStart;
                    if (prevCvGen == 0) {
                        prevCvGen = GenNum;
                        prevCvBest = currBest;
                    } else {
                        double velocity = (prevCvBest - currBest) / (elapsed > 0 ? elapsed : 1);
                        if (velocity < 5000 && (rand() / 32768.0) < 0.05)
                            mapping_operator(colony, dis_p, xColony);
                        prevCvGen = GenNum;
                        prevCvBest = currBest;
                    }
                }

                // convergence logging
                if (GenNum % 2000 == 0 && convFile) {
                    double elapsed = MPI_Wtime() - timeStart;
                    log_convergence(convFile, GenNum, elapsed, sumbest);
                }
            }
        }

        if (rank == 0) {
            double runTime = MPI_Wtime() - runStart;
            sumbest = best_distance(dis_p, xColony);
            printf("Run %d: best=%.0f  time=%.2f s\n", run + 1, sumbest, runTime);
            if (resFile)
                fprintf(resFile, "%d\t%.0f\t%.2f\n", run + 1, sumbest, runTime);
        }
    }

    if (rank == 0) {
        double totalTime = MPI_Wtime() - timeStart;
        printf("Total: %.2f s (%d runs)\n", totalTime, numRuns);
        if (resFile) fclose(resFile);
        if (convFile) fclose(convFile);
    }

    free(recv_counts); free(displs);
    free(recv_counts_d); free(displs_d);
    MPI_Finalize();
    return 0;
}
