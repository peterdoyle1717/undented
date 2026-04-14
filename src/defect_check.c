#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { double x, y, z; } Vec3;
typedef struct { int a, b, c; } Face;

typedef struct {
    Vec3 *data;
    size_t n, cap;
} Vec3Array;

typedef struct {
    Face *data;
    size_t n, cap;
} FaceArray;

static void die(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(1);
}

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n);
    if (!q) die("out of memory");
    return q;
}

static void vec3_push(Vec3Array *a, Vec3 v) {
    if (a->n == a->cap) {
        a->cap = a->cap ? 2 * a->cap : 1024;
        a->data = xrealloc(a->data, a->cap * sizeof(*a->data));
    }
    a->data[a->n++] = v;
}

static void face_push(FaceArray *a, Face f) {
    if (a->n == a->cap) {
        a->cap = a->cap ? 2 * a->cap : 1024;
        a->data = xrealloc(a->data, a->cap * sizeof(*a->data));
    }
    a->data[a->n++] = f;
}

static Vec3 vsub(Vec3 a, Vec3 b) {
    Vec3 r = {a.x - b.x, a.y - b.y, a.z - b.z};
    return r;
}

static double vdot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static double vnorm(Vec3 a) {
    return sqrt(vdot(a, a));
}

static double clamp1(double x) {
    if (x < -1.0) return -1.0;
    if (x > 1.0) return 1.0;
    return x;
}

static int parse_obj_index(const char *tok, int nverts) {
    char *end = NULL;
    errno = 0;
    long idx = strtol(tok, &end, 10);
    if (tok == end || errno) die("bad OBJ face index");
    if (idx > 0) {
        idx -= 1;
    } else if (idx < 0) {
        idx = (long)nverts + idx;
    } else {
        die("OBJ indices are 1-based; got 0");
    }
    if (idx < 0 || idx >= nverts) die("OBJ face index out of range");
    return (int)idx;
}

static void read_obj(FILE *fp, Vec3Array *verts, FaceArray *faces) {
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;

    while ((len = getline(&line, &cap, fp)) != -1) {
        if (len == 0) continue;
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '#') continue;

        if (p[0] == 'v' && isspace((unsigned char)p[1])) {
            double x, y, z;
            if (sscanf(p + 1, "%lf %lf %lf", &x, &y, &z) != 3) {
                die("bad OBJ vertex line");
            }
            Vec3 v = {x, y, z};
            vec3_push(verts, v);
            continue;
        }

        if (p[0] == 'f' && isspace((unsigned char)p[1])) {
            int idx[3];
            int count = 0;
            char *q = p + 1;
            while (*q) {
                while (*q && isspace((unsigned char)*q)) q++;
                if (!*q || *q == '\n' || *q == '\r') break;
                char *start = q;
                while (*q && !isspace((unsigned char)*q)) q++;
                char saved = *q;
                *q = '\0';
                if (count >= 3) die("non-triangular OBJ face encountered");
                idx[count++] = parse_obj_index(start, (int)verts->n);
                *q = saved;
            }
            if (count != 3) die("expected triangular OBJ faces");
            Face f = {idx[0], idx[1], idx[2]};
            face_push(faces, f);
            continue;
        }
    }

    free(line);
    if (verts->n == 0) die("no vertices found");
    if (faces->n == 0) die("no faces found");
}

static double corner_angle(Vec3 p, Vec3 q, Vec3 r) {
    Vec3 u = vsub(q, p);
    Vec3 v = vsub(r, p);
    double nu = vnorm(u);
    double nv = vnorm(v);
    if (!(nu > 0.0) || !(nv > 0.0)) die("zero-length edge encountered");
    return acos(clamp1(vdot(u, v) / (nu * nv)));
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
    Vec3Array verts = {0};
    FaceArray faces = {0};

    read_obj(fp, &verts, &faces);

    double *anglesum = calloc(verts.n, sizeof(*anglesum));
    if (!anglesum) die("out of memory");

    for (size_t i = 0; i < faces.n; i++) {
        Face f = faces.data[i];
        Vec3 a = verts.data[f.a];
        Vec3 b = verts.data[f.b];
        Vec3 c = verts.data[f.c];
        anglesum[f.a] += corner_angle(a, b, c);
        anglesum[f.b] += corner_angle(b, c, a);
        anglesum[f.c] += corner_angle(c, a, b);
    }

    double best = HUGE_VAL;
    for (size_t v = 0; v < verts.n; v++) {
        double defect = 6.28318530717958647693 - anglesum[v];
        if (defect < best) best = defect;
    }

    if (name)
        printf("%.6e %s\n", best, name);
    else
        printf("%.6e\n", best);

    free(anglesum);
    free(verts.data);
    free(faces.data);
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
