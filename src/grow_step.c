/*
 * grow_step.c — grow prime 6-nets, names only.
 *
 * Input:  CLERS names, one per line.
 * Output: CLERS names of all offspring (may contain duplicates; caller deduplicates).
 *
 * Pipeline per line:
 *   decode(name) → face list → legal insertions → insert → canon_name → print
 *
 * Compile: cc -O3 -o grow_step grow_step.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

/* ── constants ────────────────────────────────────────────────────────────── */
#define MAXV    400
#define MAXF    (2*MAXV + 4)
#define MAXCODE (4*MAXV + 8)
#define STKDEP  (MAXF * 2 + 4)

typedef struct { int a, b, c; } Face;

/* ── graph state ──────────────────────────────────────────────────────────── */
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

/* ── CLERS encoding ───────────────────────────────────────────────────────── */
static uint8_t VSET[MAXV+1];
static uint8_t TREE[MAXV+1][MAXV+1];
static short TSU[MAXV], TSW[MAXV];
static int   SKX[STKDEP], SKY[STKDEP];
static char  CUR[MAXCODE], BEST[MAXCODE];

static void encode(int a0, int b0, int rev) {
    int top=0, len=0; TNS:;
    /* TREE is clean on entry */
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

/* suppress unused-label warning */
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

/* ── insertion ────────────────────────────────────────────────────────────── */
static int has_de(Face f, int u, int w) {
    return (f.a==u&&f.b==w)||(f.b==u&&f.c==w)||(f.c==u&&f.a==w);
}
static void make_child(Face *par, int pnf, int pnv, int u, int w, int x, int y) {
    int n=pnv+1; NF=0; NV=n;
    for (int i=0; i<pnf; i++) {
        if (has_de(par[i],u,w)||has_de(par[i],w,u)) continue;
        F[NF++]=par[i];
    }
    F[NF++]=(Face){n,x,u}; F[NF++]=(Face){n,w,x};
    F[NF++]=(Face){n,y,w}; F[NF++]=(Face){n,u,y};
}

/* ── CLERS decode ─────────────────────────────────────────────────────────── */

/* union-find */
static int UF[MAXV*2+4];
static int uf_find(int x) {
    while (UF[x]!=x) { UF[x]=UF[UF[x]]; x=UF[x]; }
    return x;
}
static void uf_link(int x, int y) {
    x=uf_find(x); y=uf_find(y);
    if (x!=y) UF[x]=y;
}

/* deque (fixed array, centered so we can grow both ways) */
#define MAXDQ (MAXV*2+8)
typedef struct { int data[MAXDQ]; int h, t; } Deque;
static void dq_init(Deque *d)         { d->h=d->t=MAXDQ/2; }
static void dq_push_front(Deque *d, int v) { d->data[--d->h]=v; }
static void dq_push_back(Deque *d, int v)  { d->data[d->t++]=v; }
static int  dq_pop_back(Deque *d)          { return d->data[--d->t]; }

/* decode stack */
typedef struct { char tile; int a,b,c,phase; } Frame;
static Frame FSTK[MAXCODE];   /* DFS stack */
static Deque DSTK[MAXV];      /* deque stack */
static Face  DTRIS[MAXF];     /* accumulated triangles */

/* Decode a CLERS name into global F[]/NF/NV.
 * Returns NF, or 0 on error. */
static int decode(const char *code) {
    int n = strlen(code);
    if (!n) return 0;

    /* init union-find for up to 2*n+4 vertices */
    int maxv = n + 4;
    for (int i=0; i<=maxv; i++) UF[i]=i;

    int ptr=0, V=2, ntris=0, nstk=0, ndq=0;

    /* push first frame */
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
            if (!ph) {
                f->phase=1; char nt=code[ptr++]; V++;
                FSTK[nstk++]=(Frame){nt,c,b,V,0};
            } else {
                nstk--;
                dq_push_front(&DSTK[ndq-1], a);
                DTRIS[ntris++]=(Face){a,b,c};
            }

        } else if (t=='B') {
            if (!ph) {
                f->phase=1; char nt=code[ptr++]; V++;
                FSTK[nstk++]=(Frame){nt,a,c,V,0};
            } else {
                nstk--;
                dq_push_back(&DSTK[ndq-1], b);
                DTRIS[ntris++]=(Face){a,b,c};
            }

        } else if (t=='C') {
            if (!ph) {
                f->phase=1; char nt=code[ptr++]; V++;
                FSTK[nstk++]=(Frame){nt,a,c,V,0};
            } else {
                nstk--;
                Deque *d=&DSTK[ndq-1];
                dq_pop_back(d);              /* drop trailing c */
                int dd=dq_pop_back(d);        /* penultimate */
                dq_push_back(d, b);
                uf_link(dd, b);
                DTRIS[ntris++]=(Face){a,b,c};
            }

        } else { /* D */
            if (ph==0) {
                f->phase=1; char nt=code[ptr++]; V++;
                FSTK[nstk++]=(Frame){nt,a,c,V,0};
            } else if (ph==1) {
                f->phase=2; char nt=code[ptr++]; V++;
                FSTK[nstk++]=(Frame){nt,c,b,V,0};
            } else {
                nstk--;
                Deque *dq2=&DSTK[--ndq];
                Deque *dq =&DSTK[ndq-1];
                dq_pop_back(dq);             /* drop shared c */
                for (int i=dq2->h; i<dq2->t; i++) dq_push_back(dq, dq2->data[i]);
                DTRIS[ntris++]=(Face){a,b,c};
            }
        }
    }

    /* reverse triangles, apply find(), relabel 1..NV */
    static int lbl[MAXV*2+4];
    memset(lbl, 0, (maxv+1)*sizeof(int));
    int k=0;
    for (int i=ntris-1; i>=0; i--) {
        int a=uf_find(DTRIS[i].a);
        int b=uf_find(DTRIS[i].b);
        int c=uf_find(DTRIS[i].c);
        if (!lbl[a]) lbl[a]=++k;
        if (!lbl[b]) lbl[b]=++k;
        if (!lbl[c]) lbl[c]=++k;
        F[ntris-1-i]=(Face){lbl[a],lbl[b],lbl[c]};
    }
    NF=ntris; NV=k;
    return ntris;
}

/* ── main ─────────────────────────────────────────────────────────────────── */
int main(void) {
    static char line[MAXCODE];
    static Face par[MAXF];
    static int  ins_u[64], ins_w[64], ins_x[64], ins_y[64];
    static char name[MAXCODE];

    while (fgets(line, sizeof(line), stdin)) {
        int ll=strlen(line);
        while (ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r')) line[--ll]='\0';
        if (!ll) continue;

        /* decode name → face list */
        if (!decode(line)) continue;
        int pnf=NF, pnv=NV;
        memcpy(par, F, pnf*sizeof(Face));

        /* build EM/DEG for parent */
        build();

        /* find legal insertions */
        int n_ins=0;
        for (int u=1; u<=pnv; u++) {
            for (int i=0; i<NNBR[u]; i++) {
                int w=NBR[u][i];
                if (u>=w) continue;
                int x=EM[u][w], y=EM[w][u];
                if (!x||!y) continue;
                if (DEG[x]<=5 && DEG[y]<=5) {
                    ins_u[n_ins]=u; ins_w[n_ins]=w;
                    ins_x[n_ins]=x; ins_y[n_ins]=y;
                    n_ins++;
                }
            }
        }
        build_clear();

        /* process each insertion */
        for (int k=0; k<n_ins; k++) {
            make_child(par,pnf,pnv,ins_u[k],ins_w[k],ins_x[k],ins_y[k]);
            build();
            canon_name(name);
            puts(name);
            build_clear();
        }
    }
    return 0;
}
