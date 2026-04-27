/* dent_check.c — compute dent index (min link turning) of OBJ polyhedra.
 *
 * No args:  reads one OBJ from stdin, prints turning to stdout.
 * With args: each arg is an OBJ file, prints "turning name" per line.
 *
 * Link turning at vertex v: sum of signed spherical exterior angles
 * of the link polygon (neighbors of v on the unit sphere centered at v).
 * Positive = undented.
 *
 * Compile: cc -O3 -o dent_check dent_check.c -lm
 *
 * Usage:
 *   ./dent_check < foo.obj                       # one file, stdin
 *   ./dent_check dir/foo.obj dir/bar.obj                        # batch, one process
 *   find dir -name '*.obj' | xargs -n1000 ./dent_check   # large batch
 *   find dir -name '*.obj' | parallel -j80 -n1000 ./dent_check  # parallel
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libgen.h>

#define MAXV 400
#define MAXF (2*MAXV+4)
#define MAXRING 12

static double VX[MAXV+1], VY[MAXV+1], VZ[MAXV+1];
static int FA[MAXF], FB[MAXF], FC[MAXF];
static int NV, NF;

static int EM[MAXV+1][MAXV+1];
static int NBR[MAXV+1][MAXRING], NNBR[MAXV+1];

static void nbr_add(int u, int w) {
    for (int i = 0; i < NNBR[u]; i++) if (NBR[u][i] == w) return;
    NBR[u][NNBR[u]++] = w;
}

static int parse_obj_fp(FILE *fp) {
    char line[4096];
    NV = 0; NF = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == 'v' && line[1] == ' ') {
            NV++;
            sscanf(line + 2, "%lf %lf %lf", &VX[NV], &VY[NV], &VZ[NV]);
        } else if (line[0] == 'f' && line[1] == ' ') {
            int a, b, c;
            if (sscanf(line + 2, "%d %d %d", &a, &b, &c) == 3) {
                FA[NF] = a; FB[NF] = b; FC[NF] = c; NF++;
            }
        }
    }
    return NV > 0 && NF > 0;
}

static void build(void) {
    memset(EM, 0, sizeof(EM));
    memset(NNBR, 0, sizeof(NNBR));
    for (int i = 0; i < NF; i++) {
        int a = FA[i], b = FB[i], c = FC[i];
        EM[a][b] = c; EM[b][c] = a; EM[c][a] = b;
        nbr_add(a, b); nbr_add(b, a);
        nbr_add(b, c); nbr_add(c, b);
        nbr_add(a, c); nbr_add(c, a);
    }
}

static int cyclic_nbrs(int v, int ring[]) {
    if (!NNBR[v]) return 0;
    int start = NBR[v][0];
    ring[0] = start; int k = 1, cur = start;
    for (;;) {
        int nxt = EM[v][cur];
        if (nxt == start || nxt == 0) break;
        ring[k++] = nxt; cur = nxt;
    }
    return k;
}

static double dent_index(void) {
    double min_turn = 1e30;
    int ring[MAXRING];

    for (int v = 1; v <= NV; v++) {
        int k = cyclic_nbrs(v, ring);
        if (k < 3) continue;

        double dirs[MAXRING][3];
        for (int i = 0; i < k; i++) {
            int nb = ring[i];
            double dx = VX[nb] - VX[v], dy = VY[nb] - VY[v], dz = VZ[nb] - VZ[v];
            double len = sqrt(dx*dx + dy*dy + dz*dz);
            if (len < 1e-15) len = 1e-15;
            dirs[i][0] = dx/len; dirs[i][1] = dy/len; dirs[i][2] = dz/len;
        }

        double turn = 0;
        for (int i = 0; i < k; i++) {
            double *A = dirs[(i-1+k)%k], *B = dirs[i], *C = dirs[(i+1)%k];
            double bc0 = A[1]*C[2] - A[2]*C[1];
            double bc1 = A[2]*C[0] - A[0]*C[2];
            double bc2 = A[0]*C[1] - A[1]*C[0];
            double num = B[0]*bc0 + B[1]*bc1 + B[2]*bc2;
            double den = (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
                        *(B[0]*C[0]+B[1]*C[1]+B[2]*C[2])
                        -(A[0]*C[0]+A[1]*C[1]+A[2]*C[2]);
            turn += atan2(num, den);
        }

        if (turn < min_turn) min_turn = turn;
    }
    return min_turn;
}

/* strip directory and .obj extension from path */
static const char *stem(const char *path) {
    static char buf[4096];
    const char *p = strrchr(path, '/');
    p = p ? p + 1 : path;
    strncpy(buf, p, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    char *dot = strrchr(buf, '.');
    if (dot && strcmp(dot, ".obj") == 0) *dot = '\0';
    return buf;
}

static int process(FILE *fp, const char *name) {
    if (!parse_obj_fp(fp)) {
        fprintf(stderr, "bad input: %s\n", name);
        return 1;
    }
    build();
    if (name)
        printf("%.6e %s\n", dent_index(), name);
    else
        printf("%.6e\n", dent_index());
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        /* stdin mode: one OBJ, print just the number */
        return process(stdin, NULL);
    }

    /* batch mode: each arg is an OBJ file */
    int errors = 0;
    for (int i = 1; i < argc; i++) {
        FILE *fp = fopen(argv[i], "r");
        if (!fp) {
            fprintf(stderr, "can't open: %s\n", argv[i]);
            errors++;
            continue;
        }
        errors += process(fp, stem(argv[i]));
        fclose(fp);
    }
    return errors ? 1 : 0;
}
