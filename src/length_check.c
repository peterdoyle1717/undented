/* length_check.c — report max edge length deviation from 1.
 *
 * No args: one OBJ from stdin, prints max |length - 1| to stdout.
 * With args: each arg is an OBJ file, prints "deviation name" per line.
 *
 * Compile: cc -O3 -o length_check length_check.c -lm
 * Usage:
 *   ./length_check < foo.obj
 *   ./length_check dir/foo.obj dir/bar.obj
 *   find dir -name '*.obj' | xargs -n5000 ./length_check
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV 400
#define MAXF (2*MAXV+4)

static double VX[MAXV+1], VY[MAXV+1], VZ[MAXV+1];
static int FA[MAXF], FB[MAXF], FC[MAXF];
static int NV, NF;

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

static double edge_len(int a, int b) {
    double dx = VX[a]-VX[b], dy = VY[a]-VY[b], dz = VZ[a]-VZ[b];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

static double max_deviation(void) {
    double maxdev = 0;
    /* check each edge of each face */
    for (int i = 0; i < NF; i++) {
        double d;
        d = fabs(edge_len(FA[i], FB[i]) - 1.0); if (d > maxdev) maxdev = d;
        d = fabs(edge_len(FB[i], FC[i]) - 1.0); if (d > maxdev) maxdev = d;
        d = fabs(edge_len(FA[i], FC[i]) - 1.0); if (d > maxdev) maxdev = d;
    }
    return maxdev;
}

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
        fprintf(stderr, "bad input: %s\n", name ? name : "(stdin)");
        return 1;
    }
    double dev = max_deviation();
    if (name)
        printf("%.6e %s\n", dev, name);
    else
        printf("%.6e\n", dev);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2)
        return process(stdin, NULL);

    int errors = 0;
    for (int i = 1; i < argc; i++) {
        FILE *fp = fopen(argv[i], "r");
        if (!fp) { fprintf(stderr, "can't open: %s\n", argv[i]); errors++; continue; }
        errors += process(fp, stem(argv[i]));
        fclose(fp);
    }
    return errors ? 1 : 0;
}
