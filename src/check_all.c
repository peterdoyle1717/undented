/* check_all.c — all five embedding checks in one pass.
 *
 * Output per net: dent embed defect length sep [name]
 *   dent:   min link turning (positive = undented)
 *   embed:  1 = no self-intersection, 0 = self-intersecting
 *   defect: min angular defect (2π - angle_sum)
 *   length: max |edge_length - 1|
 *   sep:    min distance between non-adjacent vertices
 *
 * No args:  one OBJ from stdin.
 * With args: each arg is an OBJ file.
 *
 * Compile: cc -O3 -o check_all check_all.c -lm
 * Usage:
 *   ./check_all < foo.obj
 *   ./check_all dir/foo.obj dir/bar.obj
 *   find dir -name '*.obj' | xargs -n5000 ./check_all
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/* ── dent index ──────────────────────────────────────────────────────────── */
static double calc_dent(void) {
    double min_turn = 1e30;
    int ring[MAXRING];
    for (int v = 1; v <= NV; v++) {
        int k = cyclic_nbrs(v, ring);
        if (k < 3) continue;
        double dirs[MAXRING][3];
        for (int i = 0; i < k; i++) {
            int nb = ring[i];
            double dx = VX[nb]-VX[v], dy = VY[nb]-VY[v], dz = VZ[nb]-VZ[v];
            double len = sqrt(dx*dx+dy*dy+dz*dz);
            if (len < 1e-15) len = 1e-15;
            dirs[i][0]=dx/len; dirs[i][1]=dy/len; dirs[i][2]=dz/len;
        }
        double turn = 0;
        for (int i = 0; i < k; i++) {
            double *A=dirs[(i-1+k)%k], *B=dirs[i], *C=dirs[(i+1)%k];
            double bc0=A[1]*C[2]-A[2]*C[1], bc1=A[2]*C[0]-A[0]*C[2], bc2=A[0]*C[1]-A[1]*C[0];
            double num=B[0]*bc0+B[1]*bc1+B[2]*bc2;
            double den=(A[0]*B[0]+A[1]*B[1]+A[2]*B[2])*(B[0]*C[0]+B[1]*C[1]+B[2]*C[2])
                      -(A[0]*C[0]+A[1]*C[1]+A[2]*C[2]);
            turn += atan2(num, den);
        }
        if (turn < min_turn) min_turn = turn;
    }
    return min_turn;
}

/* ── self-intersection (triangle-vs-triangle) ────────────────────────────── */
/* Möller–Trumbore style: check all pairs of non-adjacent triangles. */

static void cross3(const double a[3], const double b[3], double out[3]) {
    out[0]=a[1]*b[2]-a[2]*b[1]; out[1]=a[2]*b[0]-a[0]*b[2]; out[2]=a[0]*b[1]-a[1]*b[0];
}
static double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static void sub3(const double a[3], const double b[3], double out[3]) {
    out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2];
}

/* Does segment P0→P1 intersect triangle V0,V1,V2? */
static int seg_tri(const double P0[3], const double P1[3],
                   const double V0[3], const double V1[3], const double V2[3]) {
    double dir[3], e1[3], e2[3], h[3], s[3], q[3];
    sub3(P1, P0, dir); sub3(V1, V0, e1); sub3(V2, V0, e2);
    cross3(dir, e2, h);
    double a = dot3(e1, h);
    if (fabs(a) < 1e-15) return 0;
    double f = 1.0/a;
    sub3(P0, V0, s);
    double u = f * dot3(s, h);
    if (u < 0.0 || u > 1.0) return 0;
    cross3(s, e1, q);
    double v = f * dot3(dir, q);
    if (v < 0.0 || u + v > 1.0) return 0;
    double t = f * dot3(e2, q);
    return (t > 1e-10 && t < 1.0 - 1e-10);
}

static int shares_vertex(int f1, int f2) {
    int a1=FA[f1],b1=FB[f1],c1=FC[f1], a2=FA[f2],b2=FB[f2],c2=FC[f2];
    return a1==a2||a1==b2||a1==c2||b1==a2||b1==b2||b1==c2||c1==a2||c1==b2||c1==c2;
}

static int calc_embed(void) {
    for (int i = 0; i < NF; i++) {
        double A[3]={VX[FA[i]],VY[FA[i]],VZ[FA[i]]};
        double B[3]={VX[FB[i]],VY[FB[i]],VZ[FB[i]]};
        double C[3]={VX[FC[i]],VY[FC[i]],VZ[FC[i]]};
        for (int j = i+1; j < NF; j++) {
            if (shares_vertex(i, j)) continue;
            double D[3]={VX[FA[j]],VY[FA[j]],VZ[FA[j]]};
            double E[3]={VX[FB[j]],VY[FB[j]],VZ[FB[j]]};
            double F[3]={VX[FC[j]],VY[FC[j]],VZ[FC[j]]};
            /* check all 6 edge-triangle combinations */
            if (seg_tri(A,B,D,E,F)||seg_tri(B,C,D,E,F)||seg_tri(A,C,D,E,F)||
                seg_tri(D,E,A,B,C)||seg_tri(E,F,A,B,C)||seg_tri(D,F,A,B,C))
                return 0;
        }
    }
    return 1;
}

/* ── angular defect ──────────────────────────────────────────────────────── */
static double calc_defect(void) {
    static double anglesum[MAXV+1];
    memset(anglesum, 0, (NV+1)*sizeof(double));
    for (int i = 0; i < NF; i++) {
        int a=FA[i], b=FB[i], c=FC[i];
        double ax=VX[a],ay=VY[a],az=VZ[a];
        double bx=VX[b],by=VY[b],bz=VZ[b];
        double cx=VX[c],cy=VY[c],cz=VZ[c];
        /* angle at a */
        double ux=bx-ax,uy=by-ay,uz=bz-az, vx=cx-ax,vy=cy-ay,vz=cz-az;
        double nu=sqrt(ux*ux+uy*uy+uz*uz), nv=sqrt(vx*vx+vy*vy+vz*vz);
        double d=(ux*vx+uy*vy+uz*vz)/(nu*nv); if(d>1)d=1; if(d<-1)d=-1;
        anglesum[a] += acos(d);
        /* angle at b */
        ux=ax-bx;uy=ay-by;uz=az-bz; vx=cx-bx;vy=cy-by;vz=cz-bz;
        nu=sqrt(ux*ux+uy*uy+uz*uz); nv=sqrt(vx*vx+vy*vy+vz*vz);
        d=(ux*vx+uy*vy+uz*vz)/(nu*nv); if(d>1)d=1; if(d<-1)d=-1;
        anglesum[b] += acos(d);
        /* angle at c */
        ux=ax-cx;uy=ay-cy;uz=az-cz; vx=bx-cx;vy=by-cy;vz=bz-cz;
        nu=sqrt(ux*ux+uy*uy+uz*uz); nv=sqrt(vx*vx+vy*vy+vz*vz);
        d=(ux*vx+uy*vy+uz*vz)/(nu*nv); if(d>1)d=1; if(d<-1)d=-1;
        anglesum[c] += acos(d);
    }
    double best = 1e30;
    for (int v = 1; v <= NV; v++) {
        double def = 2.0*M_PI - anglesum[v];
        if (def < best) best = def;
    }
    return best;
}

/* ── edge length ─────────────────────────────────────────────────────────── */
static double calc_length(void) {
    double maxdev = 0;
    for (int i = 0; i < NF; i++) {
        int a=FA[i],b=FB[i],c=FC[i];
        double d;
        d=sqrt((VX[a]-VX[b])*(VX[a]-VX[b])+(VY[a]-VY[b])*(VY[a]-VY[b])+(VZ[a]-VZ[b])*(VZ[a]-VZ[b]));
        d=fabs(d-1.0); if(d>maxdev) maxdev=d;
        d=sqrt((VX[b]-VX[c])*(VX[b]-VX[c])+(VY[b]-VY[c])*(VY[b]-VY[c])+(VZ[b]-VZ[c])*(VZ[b]-VZ[c]));
        d=fabs(d-1.0); if(d>maxdev) maxdev=d;
        d=sqrt((VX[a]-VX[c])*(VX[a]-VX[c])+(VY[a]-VY[c])*(VY[a]-VY[c])+(VZ[a]-VZ[c])*(VZ[a]-VZ[c]));
        d=fabs(d-1.0); if(d>maxdev) maxdev=d;
    }
    return maxdev;
}

/* ── vertex separation ──────────────────────────────────────────────────── */
static double calc_sep(void) {
    double minsep = 1e30;
    for (int i = 1; i <= NV; i++) {
        for (int j = i+1; j <= NV; j++) {
            /* Skip adjacent vertices (connected by an edge) */
            if (EM[i][j] || EM[j][i]) continue;
            double dx=VX[i]-VX[j], dy=VY[i]-VY[j], dz=VZ[i]-VZ[j];
            double d = sqrt(dx*dx+dy*dy+dz*dz);
            if (d < minsep) minsep = d;
        }
    }
    return minsep;
}

/* ── main ────────────────────────────────────────────────────────────────── */
static const char *stem(const char *path) {
    static char buf[4096];
    const char *p = strrchr(path, '/');
    p = p ? p + 1 : path;
    strncpy(buf, p, sizeof(buf)-1); buf[sizeof(buf)-1]='\0';
    char *dot = strrchr(buf, '.'); if(dot && strcmp(dot,".obj")==0) *dot='\0';
    return buf;
}

static int process(FILE *fp, const char *name) {
    if (!parse_obj_fp(fp)) {
        fprintf(stderr, "bad input: %s\n", name ? name : "(stdin)");
        return 1;
    }
    build();

    double dent   = calc_dent();
    int    embed  = calc_embed();
    double defect = calc_defect();
    double length = calc_length();
    double sep    = calc_sep();

    if (name)
        printf("%.6e %d %.6e %.6e %.6e %s\n", dent, embed, defect, length, sep, name);
    else
        printf("%.6e %d %.6e %.6e %.6e\n", dent, embed, defect, length, sep);
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
