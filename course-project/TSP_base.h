#ifndef TSP_BASE_H
#define TSP_BASE_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define CITY 442

extern double city_dis[CITY][CITY];
extern int xCity;
extern double probab1;
extern int temp[CITY];

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

static double tour_distance(int *t)
{
    double d = 0;
    for (int j = 0; j < xCity - 1; j++) d += city_dis[t[j]][t[j + 1]];
    return d + city_dis[t[0]][t[xCity - 1]];
}

static void init_population(int colony[][CITY], double dis[], int popSize)
{
    int array[CITY];
    for (int j = 0; j < xCity; j++) array[j] = j;
    for (int i = 0; i < popSize; i++) {
        int mod = xCity;
        for (int j = 0; j < xCity; j++) {
            int sign = rand() % mod;
            colony[i][j] = array[sign];
            int t = array[mod - 1];
            array[mod - 1] = array[sign];
            array[sign] = t;
            mod--;
            if (mod == 1) { j++; colony[i][j] = array[0]; break; }
        }
    }
    for (int i = 0; i < popSize; i++)
        dis[i] = tour_distance(colony[i]);
}

static void inver_over_evolve(int colony[][CITY], double dis[], int popSize)
{
    register int C1, j, k, pos_C, pos_C1;
    int k1, k2, l1, l2, pos_flag;
    register double disChange;

    for (int i = 0; i < popSize; i++) {
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
        dis[popSize + i] = dis[i] + disChange;
        for (j = 0; j < xCity; j++)
            colony[popSize + i][j] = temp[j];
    }
}

static void select_elitist(int colony[][CITY], double dis[], int popSize)
{
    for (int j = 0; j < popSize; j++)
        if (dis[popSize + j] < dis[j]) {
            dis[j] = dis[popSize + j];
            for (int k = 0; k < CITY; k++)
                colony[j][k] = colony[popSize + j][k];
        }
}

static int best_idx(double dis[], int ps)
{
    int bi = 0;
    for (int i = 1; i < ps; i++)
        if (dis[i] < dis[bi]) bi = i;
    return bi;
}

static int worst_idx(double dis[], int ps)
{
    int wi = 0;
    for (int i = 1; i < ps; i++)
        if (dis[i] > dis[wi]) wi = i;
    return wi;
}

static double best_distance(double dis[], int ps)
{
    int bi = best_idx(dis, ps);
    return dis[bi];
}

static void mapping_operator(int colony[][CITY], double dis[], int popSize)
{
    int a = rand() % popSize;
    int b = rand() % popSize;
    while (b == a) b = rand() % popSize;

    if (dis[a] < dis[b]) { int t = a; a = b; b = t; }

    int start = rand() % xCity;
    int len = (rand() % (xCity / 2)) + 1;

    int firstCity = colony[a][start];
    int bPos = position(colony[b], firstCity);

    int newTour[CITY];
    int used[CITY] = {0};

    for (int i = 0; i < len; i++) {
        int p = (start + i) % xCity;
        newTour[p] = colony[b][(bPos + i) % xCity];
        used[newTour[p]] = 1;
    }

    int pmxSrc[CITY], pmxDst[CITY];
    for (int i = 0; i < len; i++) {
        pmxSrc[i] = colony[a][(start + i) % xCity];
        pmxDst[i] = colony[b][(bPos + i) % xCity];
    }

    for (int i = 0; i < xCity; i++) {
        int inSeg = 0;
        for (int j = 0; j < len; j++) {
            if ((start + j) % xCity == i) { inSeg = 1; break; }
        }
        if (inSeg) continue;

        int city = colony[a][i];
        int dup = used[city] ? 1 : 0;

        if (dup) {
            for (int k = 0; k < len; k++) {
                if (city == pmxDst[k]) {
                    city = pmxSrc[k];
                    if (!used[city]) { dup = 0; break; }
                }
            }
            if (dup) {
                for (int j = 0; j < xCity; j++) {
                    if (!used[j]) { city = j; break; }
                }
            }
        }

        newTour[i] = city;
        used[city] = 1;
    }

    memcpy(colony[a], newTour, xCity * sizeof(int));
    dis[a] = tour_distance(colony[a]);
}

static void sort_by_fitness(int colony[][CITY], double dis[], int ps)
{
    for (int i = 1; i < ps; i++) {
        double key = dis[i];
        int keyTour[CITY];
        memcpy(keyTour, colony[i], CITY * sizeof(int));
        int j = i - 1;
        while (j >= 0 && dis[j] > key) {
            dis[j + 1] = dis[j];
            memcpy(colony[j + 1], colony[j], CITY * sizeof(int));
            j--;
        }
        dis[j + 1] = key;
        memcpy(colony[j + 1], keyTour, CITY * sizeof(int));
    }
}

static void log_convergence(FILE *fp, long gen, double elapsed, double best)
{
    fprintf(fp, "%ld\t%.4f\t%.0f\n", gen, elapsed, best);
}

static double update_mutation_rate(double p, long gen, long maxGen)
{
    double factor = 1.0 - ((double)gen / (double)maxGen) * 0.01;
    return p * factor;
}

static void check_critical_velocity(double prevBest, double currBest,
    double elapsed, int colony[][CITY], double dis[], int popSize,
    double *outPrevBest)
{
    double velocity = (prevBest - currBest) / (elapsed > 0 ? elapsed : 1.0);
    if (velocity < 5000.0 && (rand() / 32768.0) < 0.05)
        mapping_operator(colony, dis, popSize);
    *outPrevBest = currBest;
}

#endif
