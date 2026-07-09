// Algorithm 3: 移动种群分层分布式进化算法 (HdEA with Migration Population)
//
// Extends Algorithm 2 (HdEA) with peer-to-peer migration between islands
// in a ring topology. In addition to the hierarchical elite exchange:
//   - Every MIGRATION_INTERVAL, each rank sends its M = 3 best individuals
//     to the next rank in a ring (rank i -> rank (i+1)%np)
//   - Each rank receives M individuals from the previous rank
//   - Received migrants replace the worst individuals in the island
//
// The hierarchy ensures global convergence pressure; migration ensures
// diversity flow between islands, preventing premature convergence.
//
// Compile: cl /EHsc /I"%MSMPI_INC%" TSP_HdEA_MP.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
// Run: mpiexec -n 4 TSP_HdEA_MP.exe

#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define CITY      442
#define POP       50        // individuals per island
#define ELITE     50        // elite population size on master
#define MIGRATE   50        // migration interval (generations)
#define MASTER_GENS 20      // master evolution gens per interval
#define MIG_M     3         // number of migrants per exchange

static int xCity = CITY;
static double probab1 = 0.02;
static long maxGen = 200000;
static double city_dis[CITY][CITY];
static int temp[CITY];

static int  position(int *t, int c);
static void invert(int a, int b);
static double dist(int *t);
static void mutate(int *parent, int *child, double pDis, double *cDis,
                   int (*colony)[CITY], double *dps, int ps, int idx);
static void evolve_island(int (*colony)[CITY], double *dps, int ps);
static int  best_idx(double *dps, int ps);
static int  worst_idx(double *dps, int ps);
static void sort_by_fitness(int (*colony)[CITY], double *dps, int ps);

int main(int argc, char **argv)
{
    int rk, np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int runs = 10;
    if (argc > 1) runs = atoi(argv[1]);
    if (argc > 2) maxGen = atol(argv[2]);
    srand((unsigned)time(NULL) + rk);

    /* ---- read city file on root ---- */
    double xy[CITY][2];
    if (rk == 0) {
        FILE *f = fopen("pcb442.tsp", "r");
        if (!f) { printf("cannot open pcb442.tsp\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        fscanf(f, "%d", &xCity);
        for (int i = 0; i < xCity; i++) fscanf(f, "%*d%lf%lf", &xy[i][0], &xy[i][1]);
        fclose(f);
    }
    MPI_Bcast(&xCity, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(xy, CITY*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* ---- distance matrix ---- */
    for (int i = 0; i < xCity; i++)
        for (int j = 0; j < xCity; j++) {
            if (j > i) {
                double d = (xy[i][0]-xy[j][0])*(xy[i][0]-xy[j][0])
                         + (xy[i][1]-xy[j][1])*(xy[i][1]-xy[j][1]);
                city_dis[i][j] = (int)(sqrt(d)+0.5);
            } else if (j == i) city_dis[i][j] = 0;
            else city_dis[i][j] = city_dis[j][i];
        }

    /* populations */
    int colony[POP][CITY];
    double dis[POP];
    int elite[ELITE][CITY];
    double eliteDis[ELITE];

    /* results */
    FILE *log = NULL;
    if (rk == 0) {
        log = fopen("results_HdEA_MP.txt", "w");
        if (log) fprintf(log, "# HdEA-MP  run\\tbest\\ttime\\n");
    }

    double t0 = MPI_Wtime();
    MPI_Datatype MPI_TOUR;
    MPI_Type_contiguous(CITY, MPI_INT, &MPI_TOUR);
    MPI_Type_commit(&MPI_TOUR);

    /* ---- multi-run ---- */
    for (int run = 0; run < runs; run++) {
        srand((unsigned)time(NULL) + rk + run*137);
        double run0 = MPI_Wtime();

        /* init island & elite with valid random tours */
        for (int i = 0; i < POP; i++) {
            int arr[CITY];
            for (int j = 0; j < xCity; j++) arr[j] = j;
            int m = xCity;
            for (int j = 0; j < xCity; j++) {
                int s = rand() % m;
                colony[i][j] = arr[s];
                arr[s] = arr[m-1];
                if (--m == 1) { colony[i][++j] = arr[0]; break; }
            }
            dis[i] = dist(colony[i]);
        }
        for (int i = 0; i < ELITE; i++) {
            int arr[CITY];
            for (int j = 0; j < xCity; j++) arr[j] = j;
            int m = xCity;
            for (int j = 0; j < xCity; j++) {
                int s = rand() % m;
                elite[i][j] = arr[s];
                arr[s] = arr[m-1];
                if (--m == 1) { elite[i][++j] = arr[0]; break; }
            }
            eliteDis[i] = dist(elite[i]);
        }

        /* evolve */
        for (long gen = 0; gen < maxGen; gen++) {
            evolve_island(colony, dis, POP);

            if (gen > 0 && gen % MIGRATE == 0) {
                /* ===== Step 1: Hierarchical elite exchange ===== */
                int bi = best_idx(dis, POP);

                /* Point-to-point: send best to master */
                if (rk != 0) {
                    MPI_Send(colony[bi], 1, MPI_TOUR, 0, 10, MPI_COMM_WORLD);
                    MPI_Send(&dis[bi], 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
                } else {
                    int recvTour[CITY];
                    double recvDis;
                    for (int i = 1; i < np; i++) {
                        MPI_Recv(recvTour, 1, MPI_TOUR, i, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&recvDis, 1, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        int wi = worst_idx(eliteDis, ELITE);
                        if (recvDis < eliteDis[wi]) {
                            eliteDis[wi] = recvDis;
                            memcpy(elite[wi], recvTour, CITY*sizeof(int));
                        }
                    }
                    int mbi = best_idx(dis, POP);
                    int wi = worst_idx(eliteDis, ELITE);
                    if (dis[mbi] < eliteDis[wi]) {
                        eliteDis[wi] = dis[mbi];
                        memcpy(elite[wi], colony[mbi], CITY*sizeof(int));
                    }

                    for (int mg = 0; mg < MASTER_GENS; mg++)
                        evolve_island(elite, eliteDis, ELITE);
                }

                /* Collective: master broadcasts global best elite */
                int gTour[CITY];
                double gBest;
                if (rk == 0) {
                    int gbi = best_idx(eliteDis, ELITE);
                    memcpy(gTour, elite[gbi], CITY*sizeof(int));
                    gBest = eliteDis[gbi];
                }
                MPI_Bcast(gTour, 1, MPI_TOUR, 0, MPI_COMM_WORLD);
                MPI_Bcast(&gBest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                int wi = worst_idx(dis, POP);
                if (gBest < dis[wi]) {
                    dis[wi] = gBest;
                    memcpy(colony[wi], gTour, CITY*sizeof(int));
                }

                /* ===== Step 2: Peer-to-peer ring migration (NEW in Alg 3) ===== */
                sort_by_fitness(colony, dis, POP);

                int sendTo = (rk + 1) % np;
                int recvFrom = (rk - 1 + np) % np;

                /* Pack all migrants into flat buffer: [MIG_M tours][MIG_M doubles] */
                int packSize = MIG_M * (CITY + 1);
                int *sendPacked = (int *)malloc(packSize * sizeof(int));
                int *recvPacked = (int *)malloc(packSize * sizeof(int));
                for (int m = 0; m < MIG_M && m < POP; m++) {
                    memcpy(&sendPacked[m * CITY], colony[m], CITY * sizeof(int));
                    sendPacked[MIG_M * CITY + m] = (int)dis[m];
                }

                /* Ring migration via MPI_Sendrecv (avoids buffering issues) */
                MPI_Sendrecv(sendPacked, packSize, MPI_INT, sendTo, 100,
                             recvPacked, packSize, MPI_INT, recvFrom, 100,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int m = 0; m < MIG_M; m++) {
                    int wii = worst_idx(dis, POP);
                    double recvD = (double)recvPacked[MIG_M * CITY + m];
                    if (recvD < dis[wii]) {
                        dis[wii] = recvD;
                        memcpy(colony[wii], &recvPacked[m * CITY], CITY * sizeof(int));
                    }
                }

                free(sendPacked);
                free(recvPacked);
            }

            if (gen % 20000 == 0 && rk == 0) {
                double bi = dis[best_idx(dis, POP)];
                double gi = eliteDis[best_idx(eliteDis, ELITE)];
                double b = bi < gi ? bi : gi;
                printf("gen %ld  best=%.0f\n", gen, b);
            }
        }

        /* collect global best */
        double lb = dis[best_idx(dis, POP)];
        double gb;
        MPI_Reduce(&lb, &gb, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        if (rk == 0) {
            double eb = eliteDis[best_idx(eliteDis, ELITE)];
            if (eb < gb) gb = eb;
            double t = MPI_Wtime() - run0;
            printf("Run %d: best=%.0f  time=%.2f\n", run+1, gb, t);
            if (log) fprintf(log, "%d\t%.0f\t%.2f\n", run+1, gb, t);
        }
    }

    if (rk == 0) {
        printf("Total time: %.2f s\n", MPI_Wtime() - t0);
        if (log) fclose(log);
    }
    MPI_Type_free(&MPI_TOUR);
    MPI_Finalize();
    return 0;
}

/* ---- helpers ---- */
static double dist(int *t) {
    double d = 0;
    for (int j = 0; j < xCity-1; j++) d += city_dis[t[j]][t[j+1]];
    return d + city_dis[t[0]][t[xCity-1]];
}

static int position(int *t, int c) {
    for (int j = 0; j < xCity; j++) if (t[j] == c) return j;
    return -1;
}

static void invert(int a, int b) {
    if (a < b) {
        for (int i = a+1, j = b; i <= j; i++, j--)
            { int t = temp[i]; temp[i] = temp[j]; temp[j] = t; }
    } else if (xCity-1-a <= b+1) {
        int i = b, j = a+1;
        for (; j < xCity; i--, j++) { int t = temp[i]; temp[i] = temp[j]; temp[j] = t; }
        for (j = 0; j <= i; j++, i--) { int t = temp[i]; temp[i] = temp[j]; temp[j] = t; }
    } else {
        int i = b, j = a+1;
        for (; i >= 0; i--, j++) { int t = temp[i]; temp[i] = temp[j]; temp[j] = t; }
        for (i = xCity-1; j <= i; j++, i--) { int t = temp[i]; temp[i] = temp[j]; temp[j] = t; }
    }
}

static void mutate(int *parent, int *child, double pDis, double *cDis,
                   int (*colony)[CITY], double *dps, int ps, int idx)
{
    memcpy(temp, parent, CITY*sizeof(int));
    double ch = 0;
    int n = 0, pc = rand() % xCity;
    for (;;) {
        int pc1, c1, j, k1, k2, l1, l2;
        if ((rand()/32768.0) < probab1) {
            do pc1 = rand() % xCity; while (pc1 == pc);
            c1 = colony[idx][pc1];
        } else {
            do j = rand() % ps; while (j == idx);
            int k = position(colony[j], temp[pc]);
            if (k < 0) break;
            c1 = colony[j][(k+1)%xCity];
            pc1 = position(temp, c1);
            if (pc1 < 0) break;
        }
        if ((pc+1)%xCity == pc1 || (pc-1+xCity)%xCity == pc1) break;
        k1 = temp[pc]; k2 = temp[(pc+1)%xCity];
        l1 = temp[pc1]; l2 = temp[(pc1+1)%xCity];
        ch += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
        invert(pc, pc1);
        if (++n > xCity-1) break;
        if (++pc >= xCity) pc = 0;
    }
    *cDis = pDis + ch;
    memcpy(child, temp, CITY*sizeof(int));
}

static void evolve_island(int (*colony)[CITY], double *dps, int ps) {
    int off[POP][CITY];
    double offDis[POP];
    for (int i = 0; i < ps; i++)
        mutate(colony[i], off[i], dps[i], &offDis[i], colony, dps, ps, i);
    for (int i = 0; i < ps; i++)
        if (offDis[i] < dps[i]) {
            dps[i] = offDis[i];
            memcpy(colony[i], off[i], CITY*sizeof(int));
        }
}

static int best_idx(double *dps, int ps) {
    int bi = 0;
    for (int i = 1; i < ps; i++) if (dps[i] < dps[bi]) bi = i;
    return bi;
}

static int worst_idx(double *dps, int ps) {
    int wi = 0;
    for (int i = 1; i < ps; i++) if (dps[i] > dps[wi]) wi = i;
    return wi;
}

static void sort_by_fitness(int (*colony)[CITY], double *dps, int ps) {
    /* Simple insertion sort: best (lowest distance) first */
    for (int i = 1; i < ps; i++) {
        double key = dps[i];
        int keyTour[CITY];
        memcpy(keyTour, colony[i], CITY*sizeof(int));
        int j = i - 1;
        while (j >= 0 && dps[j] > key) {
            dps[j+1] = dps[j];
            memcpy(colony[j+1], colony[j], CITY*sizeof(int));
            j--;
        }
        dps[j+1] = key;
        memcpy(colony[j+1], keyTour, CITY*sizeof(int));
    }
}
