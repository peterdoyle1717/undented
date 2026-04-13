/*
 * clers_name.c — canonical CLERS name from face list.
 *
 * Input:  face lists, one graph per line: "a,b,c;d,e,f;..."  (1-indexed vertices)
 * Output: canonical CLERS name, one per line.
 *
 * Compile: cc -O3 -o clers_name_c clers_name.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define MAXV    400
#define MAXF    (2*MAXV + 4)
#define MAXLINE (MAXF * 12 + 4)
#define MAXCODE (4*MAXV + 8)
#define STKDEP  (MAXF * 2 + 4)

typedef struct { int a, b, c; } Face;

static int   NV, NF;
static Face  F[MAXF];
static int   EM[MAXV+1][MAXV+1];
static int   DEG[MAXV+1];
static int   NBR[MAXV+1][7];
static int   NNBR[MAXV+1];
static short DU[MAXF*3], DW[MAXF*3]; static int ND;

static inline void nbr_add(int u, int w) {
    for (int i = 0; i < NNBR[u]; i++) if (NBR[u][i] == w) return;
    NBR[u][NNBR[u]++] = w;
}
static void build(void) {
    memset(DEG,  0, (NV+2)*sizeof(int));
    memset(NNBR, 0, (NV+2)*sizeof(int));
    ND = 0;
    for (int i = 0; i < NF; i++) {
        int a=F[i].a, b=F[i].b, c=F[i].c;
        EM[a][b]=c; DU[ND]=a; DW[ND]=b; ND++;
        EM[b][c]=a; DU[ND]=b; DW[ND]=c; ND++;
        EM[c][a]=b; DU[ND]=c; DW[ND]=a; ND++;
        DEG[a]++; DEG[b]++; DEG[c]++;
        nbr_add(a,b); nbr_add(b,a);
        nbr_add(b,c); nbr_add(c,b);
        nbr_add(a,c); nbr_add(c,a);
    }
}
static void build_clear(void) {
    for (int i = 0; i < ND; i++) EM[DU[i]][DW[i]] = 0;
    ND = 0;
}

static uint8_t VSET[MAXV+1];
static uint8_t TREE[MAXV+1][MAXV+1];
static short TSU[MAXV], TSW[MAXV];
static int   SKX[STKDEP], SKY[STKDEP];
static char  CUR[MAXCODE], BEST[MAXCODE];

static void encode(int a0, int b0, int rev) {
    int top=0, len=0; TNS:;
    int TNS = 0;
    memset(VSET+1, 0, NV);
    VSET[a0]=VSET[b0]=1;
    TREE[b0][a0]=1; TSU[TNS]=b0; TSW[TNS]=a0; TNS++;
    SKX[top]=a0; SKY[top]=b0; top++;
    while (top > 0) {
        --top;
        int x=SKX[top], y=SKY[top];
        int z = rev ? EM[y][x] : EM[x][y];
        if (!VSET[z]) {
            CUR[len++]='C';
            TREE[z][y]=1; TSU[TNS]=z; TSW[TNS]=y; TNS++;
            VSET[z]=1; SKX[top]=x; SKY[top]=z; top++;
        } else if (TREE[y][z] && TREE[z][x]) {
            CUR[len++]='E';
        } else if (TREE[z][x]) {
            CUR[len++]='A'; SKX[top]=z; SKY[top]=y; top++;
        } else if (TREE[y][z]) {
            CUR[len++]='B'; SKX[top]=x; SKY[top]=z; top++;
        } else {
            CUR[len++]='D';
            SKX[top]=z; SKY[top]=y; top++;
            SKX[top]=x; SKY[top]=z; top++;
        }
    }
    CUR[len]='\0';
    if (!BEST[0] || strcmp(CUR,BEST)<0) memcpy(BEST,CUR,len+1);
    for (int i=0; i<TNS; i++) TREE[TSU[i]][TSW[i]]=0;
    return; (void)TNS;
}

static void canon_name(char *out) {
    BEST[0]='\0';
    int md=999;
    for (int v=1; v<=NV; v++) if (DEG[v]>0 && DEG[v]<md) md=DEG[v];
    for (int u=1; u<=NV; u++) {
        if (DEG[u]!=md) continue;
        for (int i=0; i<NNBR[u]; i++) {
            int w=NBR[u][i];
            encode(u,w,0); encode(u,w,1);
        }
    }
    strcpy(out,BEST);
}

/* Parse "a,b,c;d,e,f;..." into F[]/NF/NV. Returns NF, or 0 on error. */
static int parse_faces(char *line) {
    NF = 0; NV = 0;
    char *p = line;
    while (*p) {
        int a, b, c;
        if (sscanf(p, "%d,%d,%d", &a, &b, &c) != 3) return 0;
        F[NF++] = (Face){a, b, c};
        if (a > NV) NV = a;
        if (b > NV) NV = b;
        if (c > NV) NV = c;
        /* skip to next face */
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return NF;
}

int main(void) {
    static char line[MAXLINE];
    static char name[MAXCODE];

    while (fgets(line, sizeof(line), stdin)) {
        int ll = strlen(line);
        while (ll > 0 && (line[ll-1] == '\n' || line[ll-1] == '\r')) line[--ll] = '\0';
        if (!ll) continue;
        if (!parse_faces(line)) continue;
        build();
        canon_name(name);
        puts(name);
        build_clear();
    }
    return 0;
}
