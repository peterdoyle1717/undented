/*
 * clers.c — CLERS decode, encode, and canonical naming for triangulated spheres.
 *
 * Three modes:
 *   clers decode   < names.txt    →  face lists (a,b,c;d,e,f;...)
 *   clers encode   < faces.txt    →  canonical CLERS names
 *   clers name     < names.txt    →  canonical CLERS names (decode then encode)
 *
 * One name/face-list per line. Batch stdin, no per-line overhead.
 *
 * Compile: cc -O3 -o clers clers.c
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

/* ── shared graph state ─────────────────────────────────────────────────── */

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

/* ── encoder: face list → canonical CLERS name ──────────────────────────── */

static uint8_t VSET[MAXV+1];
static uint8_t TREE[MAXV+1][MAXV+1];
static short TSU[MAXV], TSW[MAXV];
static int   SKX[STKDEP], SKY[STKDEP];
static char  CUR[MAXCODE], BEST[MAXCODE];

static void encode(int a0, int b0, int rev) {
    int top=0, len=0, TNS=0;
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

/* ── decoder: CLERS name → face list ────────────────────────────────────── */

static int UF[MAXV*2+4];
static int uf_find(int x) {
    while (UF[x]!=x) { UF[x]=UF[UF[x]]; x=UF[x]; }
    return x;
}
static void uf_link(int x, int y) {
    x=uf_find(x); y=uf_find(y);
    if (x!=y) UF[x]=y;
}

#define MAXDQ (MAXV*2+8)
typedef struct { int data[MAXDQ]; int h, t; } Deque;
static void dq_init(Deque *d)                  { d->h=d->t=MAXDQ/2; }
static void dq_push_front(Deque *d, int v)     { d->data[--d->h]=v; }
static void dq_push_back(Deque *d, int v)      { d->data[d->t++]=v; }
static int  dq_pop_back(Deque *d)              { return d->data[--d->t]; }

typedef struct { char tile; int a,b,c,phase; } Frame;
static Frame FSTK[MAXCODE];
static Deque DSTK[MAXV];
static Face  DTRIS[MAXF];

static int decode(const char *code) {
    int n = strlen(code);
    if (!n) return 0;

    for (int i = 0; i < n; i++) {
        char c = code[i];
        if (c != 'E' && c != 'A' && c != 'B' && c != 'C' && c != 'D') {
            fprintf(stderr, "clers decode: invalid character '%c' at position %d (alphabet: ABCDE)\n", c, i);
            return 0;
        }
    }

    int maxv = n + 4;
    for (int i=0; i<=maxv; i++) UF[i]=i;

    int ptr=0, V=2, ntris=0, nstk=0, ndq=0;
    char tile = code[ptr++]; V++;
    FSTK[nstk++] = (Frame){tile, 1, 2, V, 0};

    while (nstk > 0) {
        Frame *f = &FSTK[nstk-1];
        char t=f->tile; int a=f->a, b=f->b, c=f->c, ph=f->phase;

        if (t=='E') {
            nstk--;
            dq_init(&DSTK[ndq]);
            dq_push_back(&DSTK[ndq],a);
            dq_push_back(&DSTK[ndq],c);
            dq_push_back(&DSTK[ndq],b);
            ndq++;
            DTRIS[ntris++]=(Face){a,b,c};
        } else if (t=='A') {
            if (!ph) { f->phase=1; char nt=code[ptr++]; V++;
                       FSTK[nstk++]=(Frame){nt,c,b,V,0}; }
            else { nstk--; dq_push_front(&DSTK[ndq-1],a);
                   DTRIS[ntris++]=(Face){a,b,c}; }
        } else if (t=='B') {
            if (!ph) { f->phase=1; char nt=code[ptr++]; V++;
                       FSTK[nstk++]=(Frame){nt,a,c,V,0}; }
            else { nstk--; dq_push_back(&DSTK[ndq-1],b);
                   DTRIS[ntris++]=(Face){a,b,c}; }
        } else if (t=='C') {
            if (!ph) { f->phase=1; char nt=code[ptr++]; V++;
                       FSTK[nstk++]=(Frame){nt,a,c,V,0}; }
            else { nstk--; Deque *d=&DSTK[ndq-1];
                   dq_pop_back(d); int dd=dq_pop_back(d);
                   dq_push_back(d,b); uf_link(dd,b);
                   DTRIS[ntris++]=(Face){a,b,c}; }
        } else { /* D */
            if (ph==0)    { f->phase=1; char nt=code[ptr++]; V++;
                            FSTK[nstk++]=(Frame){nt,a,c,V,0}; }
            else if (ph==1){ f->phase=2; char nt=code[ptr++]; V++;
                            FSTK[nstk++]=(Frame){nt,c,b,V,0}; }
            else { nstk--; Deque *dq2=&DSTK[--ndq]; Deque *dq=&DSTK[ndq-1];
                   dq_pop_back(dq);
                   for (int i=dq2->h; i<dq2->t; i++) dq_push_back(dq,dq2->data[i]);
                   DTRIS[ntris++]=(Face){a,b,c}; }
        }
    }

    /* relabel */
    static int lbl[MAXV*2+4];
    int maxv2 = n + 4;
    memset(lbl, 0, (maxv2+1)*sizeof(int));
    int k=0;
    for (int i=ntris-1; i>=0; i--) {
        int a=uf_find(DTRIS[i].a), b=uf_find(DTRIS[i].b), c=uf_find(DTRIS[i].c);
        if (!lbl[a]) lbl[a]=++k;
        if (!lbl[b]) lbl[b]=++k;
        if (!lbl[c]) lbl[c]=++k;
        F[ntris-1-i]=(Face){lbl[a],lbl[b],lbl[c]};
    }
    NF=ntris; NV=k;
    return ntris;
}

/* ── face list parser ───────────────────────────────────────────────────── */

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
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return NF;
}

/* ── face list formatter ────────────────────────────────────────────────── */

static void print_faces(void) {
    for (int i = 0; i < NF; i++) {
        if (i) putchar(';');
        printf("%d,%d,%d", F[i].a, F[i].b, F[i].c);
    }
    putchar('\n');
}

/* ── main ───────────────────────────────────────────────────────────────── */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: clers {decode|encode|name|canonical}\n"
                        "  decode:    CLERS name → face list\n"
                        "  encode:    face list  → CLERS name (one traversal)\n"
                        "  name:      face list  → canonical CLERS name\n"
                        "  canonical: CLERS name → canonical CLERS name\n");
        return 1;
    }

    static char line[MAXLINE];
    static char name[MAXCODE];
    int mode; /* 0=decode, 1=encode, 2=name, 3=canonical */

    if      (!strcmp(argv[1], "decode"))    mode = 0;
    else if (!strcmp(argv[1], "encode"))    mode = 1;
    else if (!strcmp(argv[1], "name"))      mode = 2;
    else if (!strcmp(argv[1], "canonical")) mode = 3;
    else { fprintf(stderr, "unknown mode: %s\n", argv[1]); return 1; }

    while (fgets(line, sizeof(line), stdin)) {
        int ll = strlen(line);
        while (ll > 0 && (line[ll-1]=='\n' || line[ll-1]=='\r')) line[--ll] = '\0';
        if (!ll) continue;

        if (mode == 0) {
            /* decode: name → faces */
            if (!decode(line)) continue;
            print_faces();
        } else if (mode == 1) {
            /* encode: faces → one CLERS name (first min-degree edge) */
            if (!parse_faces(line)) continue;
            build();
            int md=999;
            for (int v=1; v<=NV; v++) if (DEG[v]>0 && DEG[v]<md) md=DEG[v];
            for (int u=1; u<=NV; u++) if (DEG[u]==md) {
                encode(u, NBR[u][0], 0);
                break;
            }
            puts(CUR);
            build_clear();
        } else if (mode == 2) {
            /* name: faces → canonical CLERS name */
            if (!parse_faces(line)) continue;
            build();
            canon_name(name);
            puts(name);
            build_clear();
        } else {
            /* canonical: name → canonical name (decode | name) */
            if (!decode(line)) continue;
            build();
            canon_name(name);
            puts(name);
            build_clear();
        }
    }
    return 0;
}
