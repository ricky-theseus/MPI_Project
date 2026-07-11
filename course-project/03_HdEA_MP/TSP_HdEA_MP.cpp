// Algorithm 3: HdEA-MP (Moving Colony HdEA)
// 4 islands (50 each: 25 core + 25 migrant), ring topology moving colony
// Compile: cl /EHsc /I"%MSMPI_INC%" TSP_HdEA_MP.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
// Run: mpiexec -n 4 TSP_HdEA_MP.exe [runs] [maxGen]

#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define CITY     442
#define POP      50
#define CORE     25
#define MIG      25

int     xCity = CITY;
double  probab1 = 0.02;
long    maxGen = 200000;
int     xColony = POP;

double  city_dis[CITY][CITY];
int     temp[CITY];

#include "../TSP_base.h"

int main(int argc, char **argv)
{
    int rk, np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int runs = 30;
    if (argc > 1) runs = atoi(argv[1]);
    if (argc > 2) maxGen = atol(argv[2]);
    srand((unsigned)time(NULL) + rk);

    // ---- read city data ----
    double xy[CITY][2];
    if (rk == 0) {
        FILE *f = fopen("pcb442.tsp", "r");
        if (!f) { printf("cannot open pcb442.tsp\\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        fscanf(f, "%d", &xCity);
        for (int i = 0; i < xCity; i++) fscanf(f, "%*d%lf%lf", &xy[i][0], &xy[i][1]);
        fclose(f);
    }
    MPI_Bcast(&xCity, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(xy, CITY * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // ---- distance matrix ----
    for (int i = 0; i < xCity; i++)
        for (int j = 0; j < xCity; j++) {
            if (j > i) {
                double d = (xy[i][0] - xy[j][0]) * (xy[i][0] - xy[j][0])
                         + (xy[i][1] - xy[j][1]) * (xy[i][1] - xy[j][1]);
                city_dis[i][j] = (int)(sqrt(d) + 0.5);
            } else if (j == i) city_dis[i][j] = 0;
            else city_dis[i][j] = city_dis[j][i];
        }

    int colony[POP * 2][CITY];
    double dis[POP * 2];

    // ---- files ----
    FILE *resFile = NULL;
    FILE *convFile = NULL;
    if (rk == 0) {
        resFile = fopen("results_HdEA_MP.txt", "w");
        if (resFile) fprintf(resFile, "# HdEA-MP  pcb442  islands=%d  core=%d  mig=%d  maxGen=%ld\\n", np, CORE, MIG, maxGen);
        if (resFile) fprintf(resFile, "# Run\\tBest\\tTime(s)\\n");
        convFile = fopen("convergence_HdEA_MP.txt", "w");
        if (convFile) fprintf(convFile, "# Gen\\tTime\\tBest\\n");
        printf("HdEA-MP: %d procs, %d runs, %ld gens, POP=%d (core=%d + mig=%d)\\n",
               np, runs, maxGen, POP, CORE, MIG);
    }

    MPI_Datatype MPI_TOUR;
    MPI_Type_contiguous(CITY, MPI_INT, &MPI_TOUR);
    MPI_Type_commit(&MPI_TOUR);

    double t0 = MPI_Wtime();
    long MIG_INTERVAL = 20;
    long MOVING_INTERVAL = MIG_INTERVAL * 20;

    for (int run = 0; run < runs; run++) {
        srand((unsigned)time(NULL) + rk + run * 137);
        double run0 = MPI_Wtime();

        init_population(colony, dis, POP);

        long prevCvGen = 0;
        double prevCvBest = 0;
        double best = best_distance(dis, POP);

        for (long gen = 0; gen < maxGen; gen++) {
            inver_over_evolve(colony, dis, POP);
            select_elitist(colony, dis, POP);

            best = best_distance(dis, POP);

            // ---- Layer 2: Hierarchical elite exchange ----
            if (gen > 0 && gen % MIG_INTERVAL == 0) {
                int bi = best_idx(dis, POP);

                if (rk == 0) {
                    int recvTour[CITY];
                    double recvD;
                    for (int i = 1; i < np; i++) {
                        MPI_Recv(recvTour, 1, MPI_TOUR, i, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&recvD, 1, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        int wi = worst_idx(dis, POP);
                        if (recvD < dis[wi]) {
                            dis[wi] = recvD;
                            memcpy(colony[wi], recvTour, CITY * sizeof(int));
                        }
                    }
                } else {
                    MPI_Send(colony[bi], 1, MPI_TOUR, 0, 10, MPI_COMM_WORLD);
                    MPI_Send(&dis[bi], 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
                }

                int gTour[CITY];
                double gBest;
                if (rk == 0) {
                    int gbi = best_idx(dis, POP);
                    memcpy(gTour, colony[gbi], CITY * sizeof(int));
                    gBest = dis[gbi];
                }
                MPI_Bcast(gTour, 1, MPI_TOUR, 0, MPI_COMM_WORLD);
                MPI_Bcast(&gBest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                int wi = worst_idx(dis, POP);
                if (gBest < dis[wi]) {
                    dis[wi] = gBest;
                    memcpy(colony[wi], gTour, CITY * sizeof(int));
                }
            }

            // ---- Layer 3: Moving Colony (ring migration) ----
            if (gen > 0 && gen % MOVING_INTERVAL == 0) {
                sort_by_fitness(colony, dis, POP);

                int sendTo = (rk + 1) % np;
                int recvFrom = (rk - 1 + np) % np;

                int packSize = MIG * (CITY + 1);
                int *sendPacked = (int *)malloc(packSize * sizeof(int));
                int *recvPacked = (int *)malloc(packSize * sizeof(int));

                // Pack MIG (25) worst individuals = colony[CORE..CORE+MIG-1]
                for (int m = 0; m < MIG; m++) {
                    int idx = CORE + m;
                    memcpy(&sendPacked[m * CITY], colony[idx], CITY * sizeof(int));
                    sendPacked[MIG * CITY + m] = (int)dis[idx];
                }

                MPI_Sendrecv(sendPacked, packSize, MPI_INT, sendTo, 100,
                             recvPacked, packSize, MPI_INT, recvFrom, 100,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Replace migrant section with received individuals
                for (int m = 0; m < MIG; m++) {
                    int idx = CORE + m;
                    double recvD = (double)recvPacked[MIG * CITY + m];
                    dis[idx] = recvD;
                    memcpy(colony[idx], &recvPacked[m * CITY], CITY * sizeof(int));
                }

                free(sendPacked);
                free(recvPacked);
            }

            // ---- critical velocity ----
            if (gen % 2000 == 0) {
                double currBest = best;
                double elapsed = MPI_Wtime() - t0;
                if (prevCvGen == 0) {
                    prevCvGen = gen;
                    prevCvBest = currBest;
                } else {
                    double velocity = (prevCvBest - currBest) / (elapsed > 0 ? elapsed : 1);
                    if (velocity < 5000 && (rand() / 32768.0) < 0.05)
                        mapping_operator(colony, dis, POP);
                    prevCvGen = gen;
                    prevCvBest = currBest;
                }
            }

            // ---- convergence ----
            if (rk == 0 && gen % 2000 == 0 && convFile)
                log_convergence(convFile, gen, MPI_Wtime() - t0, best);
        }

        double lb = best;
        double gb;
        MPI_Reduce(&lb, &gb, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        if (rk == 0) {
            double t = MPI_Wtime() - run0;
            printf("Run %d: best=%.0f  time=%.2f\\n", run + 1, gb, t);
            if (resFile) fprintf(resFile, "%d\\t%.0f\\t%.2f\\n", run + 1, gb, t);
        }
    }

    if (rk == 0) {
        printf("Total: %.2f s (%d runs)\\n", MPI_Wtime() - t0, runs);
        if (resFile) fclose(resFile);
        if (convFile) fclose(convFile);
    }

    MPI_Type_free(&MPI_TOUR);
    MPI_Finalize();
    return 0;
}
