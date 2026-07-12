// HdEA Ring Individual — Hierarchical Distributed EA
// Two-level migration: local ring (within group) + global ring individual (between groups)
// Based on Li & Hu 2014
//
// Process layout: M groups × N subpopulations per group = total np processes
// Compile:
//   cl /EHsc /I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" TSP_HdEA_ring.cpp
//      /link /LIBPATH:"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" msmpi.lib
//      /OUT:TSP_HdEA_ring.exe
// Run:
//   mpiexec -n 4 TSP_HdEA_ring.exe <tspfile> <M> <N> <runs> <migInterval>

#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_CITY 1500
#define POP_SIZE 100
#define LOCAL_RATIO 20
#define GLOBAL_ROUNDS 100
#define LOG_INTERVAL 10000

double probab1 = 0.02;
double probab1_init = 0.02;
int xCity;
int temp[MAX_CITY];

#define COLONY_TOTAL (POP_SIZE * 2)

static int position(int *tmp, int C)
{
    for (int j = 0; j < xCity; j++)
        if (tmp[j] == C) return j;
    return -1;
}

static void invert(int pos_start, int pos_end)
{
    int j, k, t;
    if (pos_start < pos_end) {
        j = pos_start + 1; k = pos_end;
        for (; j <= k; j++, k--) { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
    } else {
        if (xCity - 1 - pos_start <= pos_end + 1) {
            j = pos_end; k = pos_start + 1;
            for (; k < xCity; j--, k++) { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
            k = 0;
            for (; k <= j; k++, j--) { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
        } else {
            j = pos_end; k = pos_start + 1;
            for (; j >= 0; j--, k++) { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
            j = xCity - 1;
            for (; k <= j; k++, j--) { t = temp[j]; temp[j] = temp[k]; temp[k] = t; }
        }
    }
}

static double tour_distance(int *t, double **city_dis)
{
    double d = 0;
    int i;
    for (i = 0; i < xCity - 1; i++) d += city_dis[t[i]][t[i + 1]];
    return d + city_dis[t[0]][t[xCity - 1]];
}

static void init_population(int (*colony)[MAX_CITY], double dis[],
                            int popSize, double **city_dis)
{
    int array[MAX_CITY];
    int i, j;
    for (j = 0; j < xCity; j++) array[j] = j;
    for (i = 0; i < popSize; i++) {
        int mod = xCity;
        for (j = 0; j < xCity; j++) {
            int sign = rand() % mod;
            colony[i][j] = array[sign];
            int t = array[mod - 1];
            array[mod - 1] = array[sign];
            array[sign] = t;
            mod--;
            if (mod == 1) { j++; colony[i][j] = array[0]; break; }
        }
    }
    for (i = 0; i < popSize; i++)
        dis[i] = tour_distance(colony[i], city_dis);
}

static void inver_over_evolve(int (*colony)[MAX_CITY], double dis[],
                              int popSize, double **city_dis)
{
    register int C1, j, k, pos_C, pos_C1;
    int i, k1, k2, l1, l2, pos_flag;
    register double disChange;

    for (i = 0; i < popSize; i++) {
        for (j = 0; j < xCity; j++) temp[j] = colony[i][j];
        disChange = 0;
        pos_flag = 0;
        pos_C = rand() % xCity;
        for (;;) {
            if ((rand() / 32768.0) < probab1) {
                do pos_C1 = rand() % xCity;
                while (pos_C1 == pos_C);
                C1 = colony[i][pos_C1];
            } else {
                do j = rand() % popSize;
                while (j == i);
                k = position(colony[j], temp[pos_C]);
                C1 = colony[j][(k + 1) % xCity];
                pos_C1 = position(temp, C1);
            }
            if ((pos_C + 1) % xCity == pos_C1 ||
                (pos_C - 1 + xCity) % xCity == pos_C1)
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
        dis[popSize + i] = dis[i] + disChange;
        for (j = 0; j < xCity; j++)
            colony[popSize + i][j] = temp[j];
    }
}

static void select_elitist(int (*colony)[MAX_CITY], double dis[], int popSize)
{
    int j, k;
    for (j = 0; j < popSize; j++)
        if (dis[popSize + j] < dis[j]) {
            dis[j] = dis[popSize + j];
            for (k = 0; k < xCity; k++)
                colony[j][k] = colony[popSize + j][k];
        }
}

static int best_idx(double dis[], int ps)
{
    int bi = 0, i;
    for (i = 1; i < ps; i++)
        if (dis[i] < dis[bi]) bi = i;
    return bi;
}

static int worst_idx(double dis[], int ps)
{
    int wi = 0, i;
    for (i = 1; i < ps; i++)
        if (dis[i] > dis[wi]) wi = i;
    return wi;
}

static void mapping_operator(int (*colony)[MAX_CITY], double dis[],
                             int popSize, double **city_dis)
{
    int i, k;
    int a = rand() % popSize;
    int b = rand() % popSize;
    while (b == a) b = rand() % popSize;
    if (dis[a] < dis[b]) { int t = a; a = b; b = t; }

    int start = rand() % xCity;
    int len = (rand() % (xCity / 2)) + 1;
    int firstCity = colony[a][start];
    int bPos = position(colony[b], firstCity);

    int newTour[MAX_CITY];
    int used[MAX_CITY] = {0};
    int pmxSrc[MAX_CITY], pmxDst[MAX_CITY];

    for (i = 0; i < len; i++) {
        int p = (start + i) % xCity;
        newTour[p] = colony[b][(bPos + i) % xCity];
        used[newTour[p]] = 1;
    }
    for (i = 0; i < len; i++) {
        pmxSrc[i] = colony[a][(start + i) % xCity];
        pmxDst[i] = colony[b][(bPos + i) % xCity];
    }
    for (i = 0; i < xCity; i++) {
        int inSeg = 0;
        for (k = 0; k < len; k++)
            if ((start + k) % xCity == i) { inSeg = 1; break; }
        if (inSeg) continue;
        int city = colony[a][i];
        int dup = used[city] ? 1 : 0;
        if (dup) {
            for (k = 0; k < len; k++) {
                if (city == pmxDst[k]) {
                    city = pmxSrc[k];
                    if (!used[city]) { dup = 0; break; }
                }
            }
            if (dup) {
                for (k = 0; k < xCity; k++)
                    if (!used[k]) { city = k; break; }
            }
        }
        newTour[i] = city;
        used[city] = 1;
    }
    memcpy(colony[a], newTour, xCity * sizeof(int));
    dis[a] = tour_distance(colony[a], city_dis);
}

int main(int argc, char **argv)
{
    int rank, np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int M = 2;
    int N = 2;
    const char *tspFile = "pcb442.tsp";
    int numRuns = 30;
    long migInterval = 1000;

    if (argc > 1) tspFile = argv[1];
    if (argc > 2) M = atoi(argv[2]);
    if (argc > 3) N = atoi(argv[3]);
    if (argc > 4) numRuns = atoi(argv[4]);
    if (argc > 5) migInterval = atol(argv[5]);

    if (M * N != np) {
        if (rank == 0)
            printf("Error: M×N(%d×%d=%d) != np(%d). "
                   "Use: mpiexec -n %d ...\n", M, N, M*N, np, M*N);
        MPI_Finalize();
        return 1;
    }

    int myGroup = rank / N;
    int myPos = rank % N;
    long localMigCount = 0;
    long globalMigCount = 0;
    long totalGen = migInterval * LOCAL_RATIO * GLOBAL_ROUNDS;

    double **city_dis = (double **)malloc(MAX_CITY * sizeof(double *));
    for (int i = 0; i < MAX_CITY; i++)
        city_dis[i] = (double *)malloc(MAX_CITY * sizeof(double));

    double coord_x[MAX_CITY], coord_y[MAX_CITY];

    if (rank == 0) {
        FILE *fp = fopen(tspFile, "r");
        if (!fp) {
            printf("Cannot open %s\n", tspFile);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(fp, "%d", &xCity);
        for (int i = 0; i < xCity; i++)
            fscanf(fp, "%*d %lf %lf", &coord_x[i], &coord_y[i]);
        fclose(fp);
    }

    MPI_Bcast(&xCity, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(coord_x, xCity, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(coord_y, xCity, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < xCity; i++)
        for (int j = 0; j < xCity; j++) {
            double dx = coord_x[i] - coord_x[j];
            double dy = coord_y[i] - coord_y[j];
            city_dis[i][j] = sqrt(dx * dx + dy * dy);
        }

    int (*colony)[MAX_CITY] = (int (*)[MAX_CITY])malloc(
        COLONY_TOTAL * MAX_CITY * sizeof(int));
    double *dis = (double *)malloc(COLONY_TOTAL * sizeof(double));
    if (!colony || !dis) {
        printf("[%d] malloc failed!\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    char resFile[256];
    snprintf(resFile, sizeof(resFile), "results_HdEA_ring.txt");

    if (rank == 0) {
        printf("HdEA Ring Individual | %s | %d cities | %d procs (M=%d,N=%d) "
               "| pop=%d | localMig=%ld | ratio=%d:1 | globalRounds=%d | "
               "totalGen=%ld | %d runs\n",
               tspFile, xCity, np, M, N, POP_SIZE,
               migInterval, LOCAL_RATIO, GLOBAL_ROUNDS, totalGen, numRuns);
    }

    for (int run = 0; run < numRuns; run++) {
        srand((unsigned)(time(NULL) + rank * 137 + run * 1000));

        init_population(colony, dis, POP_SIZE, city_dis);
        localMigCount = 0;
        globalMigCount = 0;

        double prevBest = dis[best_idx(dis, POP_SIZE)];
        double runStart = MPI_Wtime();
        long lastImproveGen = 0;

        int localNext = myGroup * N + (myPos + 1) % N;
        int localPrev = myGroup * N + (myPos - 1 + N) % N;
        int tourBuffer[MAX_CITY];

        for (long gen = 0; gen <= totalGen && globalMigCount < GLOBAL_ROUNDS; gen++) {
            probab1 = probab1_init *
                (1.0 - ((double)gen / (double)totalGen) * 0.01);

            inver_over_evolve(colony, dis, POP_SIZE, city_dis);
            select_elitist(colony, dis, POP_SIZE);

            double currBest = dis[best_idx(dis, POP_SIZE)];
            if (currBest < prevBest) {
                prevBest = currBest;
                lastImproveGen = gen;
            }
            if (gen - lastImproveGen > 2000 &&
                (rand() / 32768.0) < 0.05) {
                mapping_operator(colony, dis, POP_SIZE, city_dis);
                lastImproveGen = gen;
                prevBest = dis[best_idx(dis, POP_SIZE)];
            }

            if (gen > 0 && gen % migInterval == 0) {
                int emigrant = rand() % POP_SIZE;
                MPI_Sendrecv(colony[emigrant], xCity, MPI_INT,
                             localNext, 0,
                             tourBuffer, xCity, MPI_INT,
                             localPrev, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int worst = worst_idx(dis, POP_SIZE);
                memcpy(colony[worst], tourBuffer, xCity * sizeof(int));
                dis[worst] = tour_distance(colony[worst], city_dis);

                localMigCount++;
                if (localMigCount % LOCAL_RATIO == 0) {
                    globalMigCount++;
                    int subpopIdx = (globalMigCount - 1) % N;
                    if (myPos == subpopIdx) {
                        int srcGroup = (myGroup - 1 + M) % M;
                        int dstGroup = (myGroup + 1) % M;
                        int srcRank = srcGroup * N + subpopIdx;
                        int dstRank = dstGroup * N + subpopIdx;
                        int best = best_idx(dis, POP_SIZE);

                        MPI_Sendrecv(colony[best], xCity, MPI_INT,
                                     dstRank, 1,
                                     tourBuffer, xCity, MPI_INT,
                                     srcRank, 1,
                                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        worst = worst_idx(dis, POP_SIZE);
                        memcpy(colony[worst], tourBuffer, xCity * sizeof(int));
                        dis[worst] = tour_distance(colony[worst], city_dis);
                    }
                    if (globalMigCount >= GLOBAL_ROUNDS) break;
                }
            }

            if (gen % LOG_INTERVAL == 0 && rank == 0)
                printf("[%s] run=%d gen=%ld best=%.0f glob=%ld\n",
                       tspFile, run, gen, currBest, globalMigCount);
        }

        double localBest = dis[best_idx(dis, POP_SIZE)];
        double globalBest;
        MPI_Reduce(&localBest, &globalBest, 1, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);

        if (rank == 0) {
            double elapsed = MPI_Wtime() - runStart;
            FILE *fp = fopen(resFile, "a");
            if (run == 0) fprintf(fp, "run\tbest\telapsed\tglobalMigs\n");
            fprintf(fp, "%d\t%.0f\t%.2f\t%ld\n",
                    run, globalBest, elapsed, globalMigCount);
            fclose(fp);
            printf("[%s] run=%d best=%.0f time=%.2fs glob=%ld\n",
                   tspFile, run, globalBest, elapsed, globalMigCount);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(dis);
    free(colony);
    for (int i = 0; i < MAX_CITY; i++) free(city_dis[i]);
    free(city_dis);

    MPI_Finalize();
    return 0;
}
