// Algorithm 2: 分层分布式进化算法 (Hierarchical Distributed EA)
//
// Two-level hierarchy:
//   Level 1 (islands): each rank evolves its own independent subpopulation
//   Level 2 (master):  rank 0 collects best individuals from all islands,
//                       runs additional evolution on the elite pool, and
//                       broadcasts the global best back to all islands.
//
// Compile: cl /EHsc /I"%MSMPI_INC%" TSP_HdEA.cpp /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
// Run: mpiexec -n 4 TSP_HdEA.exe [runs] [maxGen]

#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define CITY    442
#define POP     50
#define ELITE   50
#define MIGRATE 50
#define MASTER_GENS 20

static int xCity = CITY;
static double probab1 = 0.02;
static long maxGen = 200000;
static double city_dis[CITY][CITY];
static int temp[CITY];

static int  position(int *t, int c);
static void invert(int a, int b);
static double dist(int *t);
static void mutate(int *parent, int *child, double pDis, double *cDis,
                   int (*pop)[CITY], double *dps, int ps, int idx);
static void evolve_pop(int (*pop)[CITY], double *dps, int ps);
static int  best_idx(double *dps, int ps);
static int  worst_idx(double *dps, int ps);
static void init_rand_tour(int *tour);

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

    /* read city data */
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

    for (int i = 0; i < xCity; i++)
        for (int j = 0; j < xCity; j++) {
            if (j > i) {
                double d = (xy[i][0]-xy[j][0])*(xy[i][0]-xy[j][0])
                         + (xy[i][1]-xy[j][1])*(xy[i][1]-xy[j][1]);
                city_dis[i][j] = (int)(sqrt(d)+0.5);
            } else if (j == i) city_dis[i][j] = 0;
            else city_dis[i][j] = city_dis[j][i];
        }

    int colony[POP][CITY];
    double dis[POP];
    int elite[ELITE][CITY];
    double eliteDis[ELITE];

    FILE *log = NULL;
    if (rk == 0) {
        log = fopen("results_HdEA.txt", "w");
        if (log) fprintf(log, "# HdEA  run\\tbest\\ttime\\n");
        printf("HdEA: %d procs, %d runs, %ld gens, POP=%d\n", np, runs, maxGen, POP);
    }

    double t0 = MPI_Wtime();
    MPI_Datatype MPI_TOUR;
    MPI_Type_contiguous(CITY, MPI_INT, &MPI_TOUR);
    MPI_Type_commit(&MPI_TOUR);

    for (int run = 0; run < runs; run++) {
        srand((unsigned)time(NULL) + rk + run*137);
        double run0 = MPI_Wtime();

        /* init island with valid random tours */
        for (int i = 0; i < POP; i++) {
            init_rand_tour(colony[i]);
            dis[i] = dist(colony[i]);
        }
        /* init elite with valid random tours too */
        for (int i = 0; i < ELITE; i++) {
            init_rand_tour(elite[i]);
            eliteDis[i] = dist(elite[i]);
        }

        for (long gen = 0; gen < maxGen; gen++) {
            evolve_pop(colony, dis, POP);

            if (gen > 0 && gen % MIGRATE == 0) {
                int bi = best_idx(dis, POP);

                /* workers send best to master */
                if (rk != 0) {
                    MPI_Send(colony[bi], 1, MPI_TOUR, 0, 10, MPI_COMM_WORLD);
                    MPI_Send(&dis[bi], 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
                } else {
                    /* master: receive best from each worker */
                    int recvTour[CITY];
                    double recvD;
                    for (int i = 1; i < np; i++) {
                        MPI_Recv(recvTour, 1, MPI_TOUR, i, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&recvD, 1, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        int wi = worst_idx(eliteDis, ELITE);
                        if (recvD < eliteDis[wi]) {
                            eliteDis[wi] = recvD;
                            memcpy(elite[wi], recvTour, CITY*sizeof(int));
                        }
                    }
                    /* also inject master's own best into elite */
                    if (dis[bi] < eliteDis[worst_idx(eliteDis, ELITE)]) {
                        int wi = worst_idx(eliteDis, ELITE);
                        eliteDis[wi] = dis[bi];
                        memcpy(elite[wi], colony[bi], CITY*sizeof(int));
                    }
                    /* evolve elite on master */
                    for (int mg = 0; mg < MASTER_GENS; mg++)
                        evolve_pop(elite, eliteDis, ELITE);
                }

                /* all ranks: wait for global best broadcast */
                int gTour[CITY];
                double gBest;
                if (rk == 0) {
                    int gbi = best_idx(eliteDis, ELITE);
                    memcpy(gTour, elite[gbi], CITY*sizeof(int));
                    gBest = eliteDis[gbi];
                }
                MPI_Bcast(gTour, 1, MPI_TOUR, 0, MPI_COMM_WORLD);
                MPI_Bcast(&gBest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                /* inject global best into local island */
                int wi = worst_idx(dis, POP);
                if (gBest < dis[wi]) {
                    dis[wi] = gBest;
                    memcpy(colony[wi], gTour, CITY*sizeof(int));
                }
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

static void init_rand_tour(int *tour) {
    int arr[CITY];
    for (int j = 0; j < xCity; j++) arr[j] = j;
    int m = xCity;
    for (int j = 0; j < xCity; j++) {
        int s = rand() % m;
        tour[j] = arr[s];
        arr[s] = arr[m-1];
        if (--m == 1) { tour[++j] = arr[0]; break; }
    }
}

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
                   int (*pop)[CITY], double *dps, int ps, int idx)
{
    memcpy(temp, parent, CITY*sizeof(int));
    double ch = 0;
    int n = 0, pc = rand() % xCity;
    for (;;) {
        int pc1, c1, j;
        if ((rand()/32768.0) < probab1) {
            do pc1 = rand() % xCity; while (pc1 == pc);
            c1 = pop[idx][pc1];
        } else {
            do j = rand() % ps; while (j == idx);
            int k = position(pop[j], temp[pc]);
            if (k < 0) break;
            c1 = pop[j][(k+1)%xCity];
            pc1 = position(temp, c1);
            if (pc1 < 0) break;
        }
        if ((pc+1)%xCity == pc1 || (pc-1+xCity)%xCity == pc1) break;
        int k1 = temp[pc], k2 = temp[(pc+1)%xCity];
        int l1 = temp[pc1], l2 = temp[(pc1+1)%xCity];
        ch += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
        invert(pc, pc1);
        if (++n > xCity-1) break;
        if (++pc >= xCity) pc = 0;
    }
    *cDis = pDis + ch;
    memcpy(child, temp, CITY*sizeof(int));
}

static void evolve_pop(int (*pop)[CITY], double *dps, int ps) {
    int off[POP][CITY];
    double offDis[POP];
    for (int i = 0; i < ps; i++)
        mutate(pop[i], off[i], dps[i], &offDis[i], pop, dps, ps, i);
    for (int i = 0; i < ps; i++)
        if (offDis[i] < dps[i]) {
            dps[i] = offDis[i];
            memcpy(pop[i], off[i], CITY*sizeof(int));
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
