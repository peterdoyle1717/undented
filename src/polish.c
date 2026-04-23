/* polish.c — length-polish horodump's Klein-expanded .obj to unit edges.
 *
 * For each CLERS on stdin:
 *   1. Read indir/NAME.obj (from horodump's embedded+undented gate).
 *   2. polishA: Newton on edge-length system with residual backtracking.
 *   3. If polishA stalls (final max|r|>1e-10 after 50 iters), detect flat
 *      deg-6 vertices (|turn|<1e-6) and run polishB: augmented Newton with
 *      per-flat 4 coplanarities replacing 4 spokes (per-edge dedup so the
 *      total equation count stays 3V-6).
 *   4. Test embedded (check_all-style tri-tri).
 *   5. Pass → outdir/NAME.obj. Fail → outdir/NAME.fail (same OBJ data for
 *      inspection; reason on first line after an optional '#' marker).
 *
 * The prover (Plan A vs Plan B) dispatch happens downstream by reading
 * the min-turn from the polished OBJ: ~zero → flopper, else rigid.
 *
 * Usage:   ./polish indir outdir < names.txt
 * Compile: cc -O3 -o polish polish.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ── graph state (shared with horosolve.c-style build; lightweight here) ── */
#define MAXV    400
#define MAXF    (2*MAXV+4)
#define MAXCODE (4*MAXV+8)
#define MAXRING 12
#define MAXE    (3*MAXV)
#define MAXFN   (3*MAXV)
#define STKDEP  MAXCODE
#define MAXDQ   (MAXV*2+8)

typedef struct { int a,b,c; } Face;
static int  NV, NF;
static Face F[MAXF];
static int  EM[MAXV+1][MAXV+1];
static int  NBR[MAXV+1][MAXRING], NNBR[MAXV+1];
static short DU[MAXF*3], DW[MAXF*3]; static int ND;
static int eu[MAXE], ev[MAXE], n_edges;

static void nbr_add(int u, int w){for(int i=0;i<NNBR[u];i++)if(NBR[u][i]==w)return;NBR[u][NNBR[u]++]=w;}
static void build_clear(void){for(int i=0;i<ND;i++)EM[DU[i]][DW[i]]=0;ND=0;}
static void build(void){
    memset(EM,0,(NV+2)*(NV+2)*sizeof(int)); memset(NNBR,0,(NV+2)*sizeof(int)); ND=0;
    for(int i=0;i<NF;i++){int a=F[i].a,b=F[i].b,c=F[i].c;
        EM[a][b]=c; DU[ND]=a; DW[ND]=b; ND++;
        EM[b][c]=a; DU[ND]=b; DW[ND]=c; ND++;
        EM[c][a]=b; DU[ND]=c; DW[ND]=a; ND++;
        nbr_add(a,b); nbr_add(b,a); nbr_add(b,c); nbr_add(c,b); nbr_add(a,c); nbr_add(c,a);}
}
static void collect_edges(void){
    /* skip edges among the gauge-fixed v1..v3 — their residual is constant
       and the Jacobian rows are all-zero, which would make LU singular */
    n_edges=0;
    for(int i=1;i<=NV;i++)for(int j=i+1;j<=NV;j++){
        if(!EM[i][j]) continue; if(i<=3 && j<=3) continue;
        eu[n_edges]=i; ev[n_edges]=j; n_edges++;}
}

/* ── CLERS decoder (same as horosolve.c / euclidsolve.c) ─────────────── */
static int UF[MAXV*2+4];
static int uf_find(int x){while(UF[x]!=x){UF[x]=UF[UF[x]];x=UF[x];}return x;}
static void uf_link(int x,int y){x=uf_find(x);y=uf_find(y);if(x!=y)UF[x]=y;}
typedef struct{int data[MAXDQ];int h,t;}Deque;
static void dq_init(Deque*d){d->h=d->t=MAXDQ/2;}
static void dq_push_front(Deque*d,int v){d->data[--d->h]=v;}
static void dq_push_back(Deque*d,int v){d->data[d->t++]=v;}
static int dq_pop_back(Deque*d){return d->data[--d->t];}
typedef struct{char tile;int a,b,c,phase;}Frame;
static Frame FSTK[STKDEP];static Deque DSTK[MAXV];static Face DTRIS[MAXF];
static int decode(const char*code){
    int n=strlen(code);if(!n)return 0;int maxv=n+4;for(int i=0;i<=maxv;i++)UF[i]=i;
    int ptr=0,V=2,ntris=0,nstk=0,ndq=0;char tile=code[ptr++];V++;
    FSTK[nstk++]=(Frame){tile,1,2,V,0};
    while(nstk>0){Frame*f=&FSTK[nstk-1];char t=f->tile;int a=f->a,b=f->b,c=f->c,ph=f->phase;
        if(t=='E'){nstk--;dq_init(&DSTK[ndq]);dq_push_back(&DSTK[ndq],a);dq_push_back(&DSTK[ndq],c);dq_push_back(&DSTK[ndq],b);ndq++;DTRIS[ntris++]=(Face){a,b,c};}
        else if(t=='A'){if(!ph){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,c,b,V,0};}else{nstk--;dq_push_front(&DSTK[ndq-1],a);DTRIS[ntris++]=(Face){a,b,c};}}
        else if(t=='B'){if(!ph){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,a,c,V,0};}else{nstk--;dq_push_back(&DSTK[ndq-1],b);DTRIS[ntris++]=(Face){a,b,c};}}
        else if(t=='C'){if(!ph){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,a,c,V,0};}else{nstk--;Deque*d=&DSTK[ndq-1];dq_pop_back(d);int dd=dq_pop_back(d);dq_push_back(d,b);uf_link(dd,b);DTRIS[ntris++]=(Face){a,b,c};}}
        else{if(ph==0){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,a,c,V,0};}
            else if(ph==1){f->phase=2;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,c,b,V,0};}
            else{nstk--;Deque*dq2=&DSTK[--ndq],*dq1=&DSTK[ndq-1];dq_pop_back(dq1);for(int i=dq2->h;i<dq2->t;i++)dq_push_back(dq1,dq2->data[i]);DTRIS[ntris++]=(Face){a,b,c};}}}
    static int lbl[MAXV*2+4];memset(lbl,0,(maxv+1)*sizeof(int));int k=0;
    for(int i=ntris-1;i>=0;i--){int a=uf_find(DTRIS[i].a),b=uf_find(DTRIS[i].b),c=uf_find(DTRIS[i].c);
        if(!lbl[a])lbl[a]=++k;if(!lbl[b])lbl[b]=++k;if(!lbl[c])lbl[c]=++k;F[ntris-1-i]=(Face){lbl[a],lbl[b],lbl[c]};}
    NF=ntris;NV=k;return ntris;
}

/* ── coords + gauge ──────────────────────────────────────────────────── */
static double pos[MAXV+1][3];        /* v1..vNV Euclidean coords */
static double e_xvec[MAXFN];         /* free coords for v=4..NV (3*(NV-3)) */
static double EC[4][3];              /* gauge anchors after reading OBJ: v1,v2,v3 */
#define EX(v) ((v)<=3?EC[v][0]:e_xvec[3*((v)-4)])
#define EY(v) ((v)<=3?EC[v][1]:e_xvec[3*((v)-4)+1])
#define EZ(v) ((v)<=3?EC[v][2]:e_xvec[3*((v)-4)+2])

/* ── OBJ reader: populates pos[] from file; returns 1 on success ─────── */
static int read_obj_coords(const char *path, int expect_nv){
    FILE *fp=fopen(path,"r");
    if(!fp) return 0;
    char ln[1024]; int vi=0;
    while(fgets(ln,sizeof(ln),fp)){
        if(ln[0]=='v' && (ln[1]==' '||ln[1]=='\t')){
            vi++;
            if(vi>expect_nv){fclose(fp); return 0;}
            if(sscanf(ln+2,"%lf %lf %lf",&pos[vi][0],&pos[vi][1],&pos[vi][2])!=3){fclose(fp); return 0;}
        }
    }
    fclose(fp);
    return vi==expect_nv;
}

/* ── Newton linear solver (LU with partial pivoting) ─────────────────── */
static double NJmat[MAXFN][MAXFN], NFvec[MAXFN];
static int lu_solve_n(int n){
    for(int col=0;col<n;col++){int piv=col;double best=fabs(NJmat[col][col]);
        for(int row=col+1;row<n;row++)if(fabs(NJmat[row][col])>best){best=fabs(NJmat[row][col]);piv=row;}
        if(best<1e-14) return -1;
        if(piv!=col){for(int k=col;k<n;k++){double t=NJmat[col][k];NJmat[col][k]=NJmat[piv][k];NJmat[piv][k]=t;}
                     {double t=NFvec[col];NFvec[col]=NFvec[piv];NFvec[piv]=t;}}
        double inv=1.0/NJmat[col][col];
        for(int row=col+1;row<n;row++){double fac=NJmat[row][col]*inv;
            for(int k=col;k<n;k++)NJmat[row][k]-=fac*NJmat[col][k];
            NFvec[row]-=fac*NFvec[col];}}
    for(int i=n-1;i>=0;i--){double s=NFvec[i];for(int j=i+1;j<n;j++)s-=NJmat[i][j]*NFvec[j];NFvec[i]=s/NJmat[i][i];}
    return 0;
}

/* ── Turning at vertex v (ordered link via EM) ────────────────────────── */
static void cyclic_link(int v, int *link, int *klen_out){
    int k=0;
    if(!NNBR[v]){*klen_out=0; return;}
    int start=NBR[v][0];
    link[k++]=start;
    int cur=start;
    while(k<NNBR[v]){int nxt=EM[v][cur]; if(!nxt||nxt==start) break; link[k++]=nxt; cur=nxt;}
    *klen_out=k;
}
static double vertex_turn(int v){
    int link[MAXRING], k; cyclic_link(v, link, &k);
    if(k<3) return 0.0;
    double dirs[MAXRING][3];
    for(int i=0;i<k;i++){
        double dx=pos[link[i]][0]-pos[v][0];
        double dy=pos[link[i]][1]-pos[v][1];
        double dz=pos[link[i]][2]-pos[v][2];
        double n=sqrt(dx*dx+dy*dy+dz*dz);
        dirs[i][0]=dx/n; dirs[i][1]=dy/n; dirs[i][2]=dz/n;
    }
    double turn=0;
    for(int i=0;i<k;i++){
        int im=(i+k-1)%k, ip=(i+1)%k;
        double *dA=dirs[im], *dB=dirs[i], *dC=dirs[ip];
        double cAC[3]={dA[1]*dC[2]-dA[2]*dC[1],
                      dA[2]*dC[0]-dA[0]*dC[2],
                      dA[0]*dC[1]-dA[1]*dC[0]};
        double num=dB[0]*cAC[0]+dB[1]*cAC[1]+dB[2]*cAC[2];
        double dAB=dA[0]*dB[0]+dA[1]*dB[1]+dA[2]*dB[2];
        double dBC=dB[0]*dC[0]+dB[1]*dC[1]+dB[2]*dC[2];
        double dAC=dA[0]*dC[0]+dA[1]*dC[1]+dA[2]*dC[2];
        double den=dAB*dBC-dAC;
        turn+=atan2(num,den);
    }
    return turn;
}
static double min_turn_all(void){
    double mn=1e30;
    for(int v=1;v<=NV;v++){double t=vertex_turn(v); if(t<mn) mn=t;}
    return mn;
}

/* ── embedded check (tri-tri intersection) ───────────────────────────── */
static int seg_tri(const double P0[3], const double P1[3],
                   const double V0[3], const double V1[3], const double V2[3]){
    double dir[3]={P1[0]-P0[0],P1[1]-P0[1],P1[2]-P0[2]};
    double e1[3]={V1[0]-V0[0],V1[1]-V0[1],V1[2]-V0[2]};
    double e2[3]={V2[0]-V0[0],V2[1]-V0[1],V2[2]-V0[2]};
    double h[3]={dir[1]*e2[2]-dir[2]*e2[1],
                 dir[2]*e2[0]-dir[0]*e2[2],
                 dir[0]*e2[1]-dir[1]*e2[0]};
    double a=e1[0]*h[0]+e1[1]*h[1]+e1[2]*h[2];
    if(fabs(a)<1e-15) return 0;
    double f=1.0/a;
    double s[3]={P0[0]-V0[0],P0[1]-V0[1],P0[2]-V0[2]};
    double u=f*(s[0]*h[0]+s[1]*h[1]+s[2]*h[2]);
    if(u<0.0||u>1.0) return 0;
    double q[3]={s[1]*e1[2]-s[2]*e1[1],
                 s[2]*e1[0]-s[0]*e1[2],
                 s[0]*e1[1]-s[1]*e1[0]};
    double v=f*(dir[0]*q[0]+dir[1]*q[1]+dir[2]*q[2]);
    if(v<0.0||u+v>1.0) return 0;
    double t=f*(e2[0]*q[0]+e2[1]*q[1]+e2[2]*q[2]);
    return (t>1e-10 && t<1.0-1e-10);
}
static int is_embedded(void){
    for(int i=0;i<NF;i++){
        int ai=F[i].a,bi=F[i].b,ci=F[i].c;
        for(int j=i+1;j<NF;j++){
            int aj=F[j].a,bj=F[j].b,cj=F[j].c;
            if(ai==aj||ai==bj||ai==cj||bi==aj||bi==bj||bi==cj||
               ci==aj||ci==bj||ci==cj) continue;
            if(seg_tri(pos[ai],pos[bi],pos[aj],pos[bj],pos[cj])||
               seg_tri(pos[bi],pos[ci],pos[aj],pos[bj],pos[cj])||
               seg_tri(pos[ai],pos[ci],pos[aj],pos[bj],pos[cj])||
               seg_tri(pos[aj],pos[bj],pos[ai],pos[bi],pos[ci])||
               seg_tri(pos[bj],pos[cj],pos[ai],pos[bi],pos[ci])||
               seg_tri(pos[aj],pos[cj],pos[ai],pos[bi],pos[ci])) return 0;
        }
    }
    return 1;
}

/* ── polishA: Newton on edge lengths. Returns final max|r|. ────────────
   Sets *stall_out to 1 if Newton showed stall signals (linear-at-best
   convergence typical of a near-singular Jacobian = flopper territory):
     - >= 15 iterations used (rigids converge in 3-6)
     - any iteration where final backtrack step < 1/32 (heavy backtracking)
     - residual ratio between consecutive iters > 0.5 for >= 3 iters
*/
static double polishA(int max_iter, int *stall_out){
    int n=3*(NV-3);
    double res=1e30, prev_res=1e30;
    int stall_hits=0, iters_used=0, heavy_bt=0;
    for(int iter=0;iter<max_iter;iter++){
        res=0;
        for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            NFvec[k]=dx*dx+dy*dy+dz*dz-1.0;
            double af=fabs(NFvec[k]);if(af>res)res=af;}
        iters_used=iter+1;
        if(iter>0 && res>0.5*prev_res) stall_hits++;
        prev_res=res;
        if(res<1e-14) break;
        memset(NJmat,0,sizeof(double)*(size_t)n*MAXFN);
        for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            if(i>=4){int ci=3*(i-4);NJmat[k][ci]+=2*dx;NJmat[k][ci+1]+=2*dy;NJmat[k][ci+2]+=2*dz;}
            if(j>=4){int cj=3*(j-4);NJmat[k][cj]-=2*dx;NJmat[k][cj+1]-=2*dy;NJmat[k][cj+2]-=2*dz;}}
        for(int k=0;k<n;k++)NFvec[k]=-NFvec[k];
        if(lu_solve_n(n)<0) break;
        double step=1.0;
        for(int bt=0;bt<40;bt++,step*=0.5){
            double res2=0;
            for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
                double dx=(i>=4?EX(i)+step*NFvec[3*(i-4)]:EX(i))-(j>=4?EX(j)+step*NFvec[3*(j-4)]:EX(j));
                double dy=(i>=4?EY(i)+step*NFvec[3*(i-4)+1]:EY(i))-(j>=4?EY(j)+step*NFvec[3*(j-4)+1]:EY(j));
                double dz=(i>=4?EZ(i)+step*NFvec[3*(i-4)+2]:EZ(i))-(j>=4?EZ(j)+step*NFvec[3*(j-4)+2]:EZ(j));
                double af=fabs(dx*dx+dy*dy+dz*dz-1.0);if(af>res2)res2=af;}
            if(res2<res) break;
        }
        if(step<1.0/32.0) heavy_bt++;
        for(int v=4;v<=NV;v++){
            e_xvec[3*(v-4)  ]+=step*NFvec[3*(v-4)  ];
            e_xvec[3*(v-4)+1]+=step*NFvec[3*(v-4)+1];
            e_xvec[3*(v-4)+2]+=step*NFvec[3*(v-4)+2];}
    }
    if(stall_out){
        *stall_out = (iters_used>=15) || (heavy_bt>=1) || (stall_hits>=3);
    }
    return res;
}

/* ── flat detection + polishB (augmented: 4 coplanarities per flat, ─── */
/*    4 spokes dropped per flat, edge-dedup for shared spokes)          */

static int flat_set[MAXV+1];
static int keep1[MAXV+1], keep2[MAXV+1];

static int identify_flats(double thresh){
    /* Requires pos[] current (copy from e_xvec first if needed). */
    memset(flat_set,0,(NV+2)*sizeof(int));
    int count=0;
    for(int v=1;v<=NV;v++){
        if(NNBR[v]!=6) continue;
        double t=vertex_turn(v);
        if(fabs(t)<thresh){flat_set[v]=1; count++;}
    }
    return count;
}

/* Pick 2 link verts to keep at flat v. Prefer flat neighbors. */
static void pick_keep(int v){
    int link[MAXRING], k; cyclic_link(v, link, &k);
    int flat_n[MAXRING], nf=0, non_n[MAXRING], nn=0;
    for(int i=0;i<k;i++){
        if(flat_set[link[i]]) flat_n[nf++]=link[i];
        else non_n[nn++]=link[i];
    }
    if(nf>=2){keep1[v]=flat_n[0]; keep2[v]=flat_n[1];}
    else if(nf==1){keep1[v]=flat_n[0]; keep2[v]=non_n[0];}
    else {keep1[v]=non_n[0]; keep2[v]=non_n[1];}
}

/* Augmented system. Equations:
 *   (A) length-1 for every edge (i,j) NOT replaced by a coplanarity.
 *   (B) for each replaced edge (v,o) where v flat and o∈link(v)\{keep1,keep2}:
 *       det(pos[o]-pos[v], pos[k1]-pos[v], pos[k2]-pos[v]) = 0
 *
 * An edge (v,o) with both endpoints flat is dedup'd: if v's coplanarity
 * claims it, u's coplanarity can't (whichever was first wins).
 */

typedef struct { int v, k1, k2, o; } Coplan;
static Coplan coplans[MAXE];
static int n_coplan;
static int kept_eu[MAXE], kept_ev[MAXE], n_kept_edges;
static int edge_replaced[MAXV+1][MAXV+1];   /* 1 if edge replaced by some coplanarity */

static void build_augmented(void){
    memset(edge_replaced,0,(NV+2)*(NV+2)*sizeof(int));
    n_coplan=0;
    for(int v=1;v<=NV;v++){
        if(!flat_set[v]) continue;
        int link[MAXRING], k; cyclic_link(v, link, &k);
        int k1=keep1[v], k2=keep2[v];
        for(int i=0;i<k;i++){
            int o=link[i];
            if(o==k1||o==k2) continue;
            /* edge (v,o) — only replace if not already replaced by u<v's coplan */
            int a=v<o?v:o, b=v<o?o:v;
            if(edge_replaced[a][b]) continue;
            edge_replaced[a][b]=1;
            coplans[n_coplan++]=(Coplan){v,k1,k2,o};
        }
    }
    n_kept_edges=0;
    for(int k=0;k<n_edges;k++){
        int i=eu[k], j=ev[k];
        int a=i<j?i:j, b=i<j?j:i;
        if(edge_replaced[a][b]) continue;
        kept_eu[n_kept_edges]=i; kept_ev[n_kept_edges]=j; n_kept_edges++;
    }
}

static double polishB(int max_iter){
    int n=3*(NV-3);
    int neq=n_kept_edges+n_coplan;
    double res=1e30;
    for(int iter=0;iter<max_iter;iter++){
        res=0;
        /* edge residuals */
        for(int k=0;k<n_kept_edges;k++){int i=kept_eu[k], j=kept_ev[k];
            double dx=EX(i)-EX(j), dy=EY(i)-EY(j), dz=EZ(i)-EZ(j);
            NFvec[k]=dx*dx+dy*dy+dz*dz-1.0;
            double af=fabs(NFvec[k]); if(af>res) res=af;}
        /* coplan residuals: det([o-v, k1-v, k2-v]) */
        for(int k=0;k<n_coplan;k++){
            int v=coplans[k].v, k1=coplans[k].k1, k2=coplans[k].k2, o=coplans[k].o;
            double ox=EX(o)-EX(v), oy=EY(o)-EY(v), oz=EZ(o)-EZ(v);
            double ax=EX(k1)-EX(v), ay=EY(k1)-EY(v), az=EZ(k1)-EZ(v);
            double bx=EX(k2)-EX(v), by=EY(k2)-EY(v), bz=EZ(k2)-EZ(v);
            double d=ox*(ay*bz-az*by)-oy*(ax*bz-az*bx)+oz*(ax*by-ay*bx);
            NFvec[n_kept_edges+k]=d;
            double af=fabs(d); if(af>res) res=af;
        }
        if(res<1e-14) break;
        memset(NJmat,0,sizeof(double)*(size_t)neq*MAXFN);
        /* edge rows */
        for(int k=0;k<n_kept_edges;k++){int i=kept_eu[k], j=kept_ev[k];
            double dx=EX(i)-EX(j), dy=EY(i)-EY(j), dz=EZ(i)-EZ(j);
            if(i>=4){int ci=3*(i-4); NJmat[k][ci]+=2*dx; NJmat[k][ci+1]+=2*dy; NJmat[k][ci+2]+=2*dz;}
            if(j>=4){int cj=3*(j-4); NJmat[k][cj]-=2*dx; NJmat[k][cj+1]-=2*dy; NJmat[k][cj+2]-=2*dz;}}
        /* coplan rows: Det[o-v, k1-v, k2-v]. Let a=k1-v, b=k2-v, c=o-v.
         * det = c . (a × b)
         * d/do = a × b
         * d/d(k1) = b × c
         * d/d(k2) = c × a
         * d/dv = -(above sum)
         */
        for(int k=0;k<n_coplan;k++){
            int v=coplans[k].v, k1=coplans[k].k1, k2=coplans[k].k2, o=coplans[k].o;
            double ax=EX(k1)-EX(v), ay=EY(k1)-EY(v), az=EZ(k1)-EZ(v);
            double bx=EX(k2)-EX(v), by=EY(k2)-EY(v), bz=EZ(k2)-EZ(v);
            double cx=EX(o)-EX(v),  cy=EY(o)-EY(v),  cz=EZ(o)-EZ(v);
            double gO[3]={ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx};
            double gA[3]={by*cz-bz*cy, bz*cx-bx*cz, bx*cy-by*cx};
            double gB[3]={cy*az-cz*ay, cz*ax-cx*az, cx*ay-cy*ax};
            int row=n_kept_edges+k;
            if(o>=4){int co=3*(o-4); NJmat[row][co]+=gO[0]; NJmat[row][co+1]+=gO[1]; NJmat[row][co+2]+=gO[2];}
            if(k1>=4){int ck=3*(k1-4); NJmat[row][ck]+=gA[0]; NJmat[row][ck+1]+=gA[1]; NJmat[row][ck+2]+=gA[2];}
            if(k2>=4){int cm=3*(k2-4); NJmat[row][cm]+=gB[0]; NJmat[row][cm+1]+=gB[1]; NJmat[row][cm+2]+=gB[2];}
            if(v>=4){int cv=3*(v-4);
                NJmat[row][cv]-=(gO[0]+gA[0]+gB[0]);
                NJmat[row][cv+1]-=(gO[1]+gA[1]+gB[1]);
                NJmat[row][cv+2]-=(gO[2]+gA[2]+gB[2]);}
        }
        for(int k=0;k<neq;k++)NFvec[k]=-NFvec[k];
        /* Square system? (neq should equal n for a well-posed flopper). If
         * neq > n, use the first n rows (over-determined case with linearly
         * dependent extras, which should be OK for valid floppers). */
        int solve_n = neq<=n ? neq : n;
        if(neq<n) break;     /* under-determined — can't proceed */
        if(lu_solve_n(solve_n)<0) break;
        double step=1.0;
        for(int bt=0;bt<40;bt++,step*=0.5){
            double res2=0;
            for(int k=0;k<n_kept_edges;k++){int i=kept_eu[k], j=kept_ev[k];
                double dx=(i>=4?EX(i)+step*NFvec[3*(i-4)]:EX(i))-(j>=4?EX(j)+step*NFvec[3*(j-4)]:EX(j));
                double dy=(i>=4?EY(i)+step*NFvec[3*(i-4)+1]:EY(i))-(j>=4?EY(j)+step*NFvec[3*(j-4)+1]:EY(j));
                double dz=(i>=4?EZ(i)+step*NFvec[3*(i-4)+2]:EZ(i))-(j>=4?EZ(j)+step*NFvec[3*(j-4)+2]:EZ(j));
                double af=fabs(dx*dx+dy*dy+dz*dz-1.0); if(af>res2) res2=af;}
            if(res2<res) break;
        }
        for(int v=4;v<=NV;v++){
            e_xvec[3*(v-4)  ]+=step*NFvec[3*(v-4)  ];
            e_xvec[3*(v-4)+1]+=step*NFvec[3*(v-4)+1];
            e_xvec[3*(v-4)+2]+=step*NFvec[3*(v-4)+2];}
    }
    return res;
}

/* Pull pos[] from current free coords + gauge. */
static void sync_pos_from_state(void){
    pos[1][0]=EC[1][0]; pos[1][1]=EC[1][1]; pos[1][2]=EC[1][2];
    pos[2][0]=EC[2][0]; pos[2][1]=EC[2][1]; pos[2][2]=EC[2][2];
    pos[3][0]=EC[3][0]; pos[3][1]=EC[3][1]; pos[3][2]=EC[3][2];
    for(int v=4;v<=NV;v++){
        pos[v][0]=e_xvec[3*(v-4)];
        pos[v][1]=e_xvec[3*(v-4)+1];
        pos[v][2]=e_xvec[3*(v-4)+2];
    }
}

/* Gauge-align pos[] to (v1=origin, v2=+x, v3=xy upper-half-plane). */
static void gauge_align(void){
    double shift[3]={pos[1][0],pos[1][1],pos[1][2]};
    for(int v=1;v<=NV;v++){pos[v][0]-=shift[0]; pos[v][1]-=shift[1]; pos[v][2]-=shift[2];}
    double ex[3]={pos[2][0],pos[2][1],pos[2][2]};
    double nx=sqrt(ex[0]*ex[0]+ex[1]*ex[1]+ex[2]*ex[2]);
    ex[0]/=nx; ex[1]/=nx; ex[2]/=nx;
    double ey[3]={pos[3][0],pos[3][1],pos[3][2]};
    double dp=ex[0]*ey[0]+ex[1]*ey[1]+ex[2]*ey[2];
    ey[0]-=dp*ex[0]; ey[1]-=dp*ex[1]; ey[2]-=dp*ex[2];
    double ny=sqrt(ey[0]*ey[0]+ey[1]*ey[1]+ey[2]*ey[2]);
    ey[0]/=ny; ey[1]/=ny; ey[2]/=ny;
    double ez[3]={ex[1]*ey[2]-ex[2]*ey[1], ex[2]*ey[0]-ex[0]*ey[2], ex[0]*ey[1]-ex[1]*ey[0]};
    for(int v=1;v<=NV;v++){
        double x=pos[v][0]*ex[0]+pos[v][1]*ex[1]+pos[v][2]*ex[2];
        double y=pos[v][0]*ey[0]+pos[v][1]*ey[1]+pos[v][2]*ey[2];
        double z=pos[v][0]*ez[0]+pos[v][1]*ez[1]+pos[v][2]*ez[2];
        pos[v][0]=x; pos[v][1]=y; pos[v][2]=z;
    }
}

/* ── main ────────────────────────────────────────────────────────────── */
int main(int argc, char **argv){
    if(argc<3){fprintf(stderr,"usage: polish indir outdir < names.txt\n"); return 1;}
    const char *indir=argv[1], *outdir=argv[2];
    static char line[MAXCODE];
    long nets=0, ok_count=0, fail_count=0;
    long pA_ok=0, pB_ok=0, both_fail=0, no_input=0, not_embedded=0;

    long skipped=0;
    while(fgets(line,sizeof(line),stdin)){
        int ll=strlen(line);
        while(ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r')) line[--ll]='\0';
        if(!ll) continue;

        /* Incremental: if outdir/NAME.obj already exists from a prior run,
           skip. Don't retouch successful polishA work. To force a rerun
           on some subset, delete those outputs first. */
        char path[4096];
        snprintf(path,sizeof(path),"%s/%s.obj",outdir,line);
        FILE *check=fopen(path,"r");
        if(check){fclose(check); skipped++; nets++; continue;}

        if(!decode(line)){fprintf(stderr,"decode failed: %s\n",line); continue;}
        build(); collect_edges();

        snprintf(path,sizeof(path),"%s/%s.obj",indir,line);
        if(!read_obj_coords(path, NV)){
            snprintf(path,sizeof(path),"%s/%s.fail",outdir,line);
            FILE *fp=fopen(path,"w"); if(fp){fprintf(fp,"# no-input\n"); fclose(fp);}
            no_input++; fail_count++; nets++; build_clear(); continue;
        }

        /* Snap v1,v2,v3 to EXACT horodump-gauge positions — no rotation, no
           rescaling of the bulk. Input pos[1..3] are approximate (O(1e-10)
           off from these); their imprecision would otherwise propagate as
           unsatisfiable constant edge residuals in Newton and floor the
           polish at ~3e-9 per edge. v4..NV stay at their input coords. */
        EC[1][0] = 0.0;               EC[1][1] = 0.0; EC[1][2] =  0.5;
        EC[2][0] = 0.0;               EC[2][1] = 0.0; EC[2][2] = -0.5;
        EC[3][0] = 0.86602540378443864676; /* √3/2 */
        EC[3][1] = 0.0;               EC[3][2] = 0.0;
        for(int v=4;v<=NV;v++){
            e_xvec[3*(v-4)  ]=pos[v][0];
            e_xvec[3*(v-4)+1]=pos[v][1];
            e_xvec[3*(v-4)+2]=pos[v][2];
        }

        /* polishA */
        int stallA=0;
        double rA=polishA(50, &stallA);
        /* Check the polished geometry too — polishA can drive edge residual
           to 1e-14 but land on a near-flat (crushed pancake) config where
           min_turn ~ 1e-8. That's the case where polishB should retry from
           the Klein starting point. Rigid polyhedra have min_turn ~ O(1),
           well clear of this threshold. */
        sync_pos_from_state();
        double mtA = min_turn_all();
        int pA_success = (rA<1e-10) && (mtA > 1e-4);

        double rB=1e30;
        int pB_success=0;
        if(!pA_success){
            /* Reset free coords to the Klein starting config. polishA may
               have slid the pancake along its flop family toward the flat
               limit, leaving a near-degenerate state that polishB can't
               open back up. Restart polishB from the un-crushed Klein
               input. */
            for(int v=4;v<=NV;v++){
                e_xvec[3*(v-4)  ]=pos[v][0];
                e_xvec[3*(v-4)+1]=pos[v][1];
                e_xvec[3*(v-4)+2]=pos[v][2];
            }
            /* pos[] was overwritten by sync_pos_from_state; reread .obj to
               restore Klein input. */
            char reset_path[4096];
            snprintf(reset_path,sizeof(reset_path),"%s/%s.obj",indir,line);
            read_obj_coords(reset_path, NV);
            for(int v=4;v<=NV;v++){
                e_xvec[3*(v-4)  ]=pos[v][0];
                e_xvec[3*(v-4)+1]=pos[v][1];
                e_xvec[3*(v-4)+2]=pos[v][2];
            }
            /* Need current pos to identify flats — sync from e_xvec state. */
            sync_pos_from_state();
            /* Flat threshold scales to the current residual. Position noise
               is O(rA); turn-noise ~ sqrt(position*eps) per catastrophic-
               cancellation rule ~ sqrt(rA). Keep floor at 1e-6 for tight
               residuals; cap at 1e-3 so we don't grab regular vertices. */
            double thresh=10.0*sqrt(rA);
            if(thresh<1e-6) thresh=1e-6;
            if(thresh>1e-3) thresh=1e-3;
            int nflat=identify_flats(thresh);
            if(nflat>0){
                for(int v=1;v<=NV;v++) if(flat_set[v]) pick_keep(v);
                build_augmented();
                rB=polishB(50);
                pB_success=(rB<1e-10);
            }
        }

        /* Final state */
        sync_pos_from_state();
        int pass=0;
        const char *reason="";
        if(pA_success){pass=is_embedded(); if(!pass) reason="A:not-embedded"; else {pA_ok++;}}
        else if(pB_success){pass=is_embedded(); if(!pass) reason="B:not-embedded"; else {pB_ok++;}}
        else {reason="polish-stall"; both_fail++;}
        if(!pass && *reason==0) reason="unknown";

        snprintf(path,sizeof(path),"%s/%s.%s",outdir,line,pass?"obj":"fail");
        FILE *fp=fopen(path,"w");
        if(fp){
            if(!pass) fprintf(fp,"# %s rA=%.3e rB=%.3e\n",reason,rA,rB);
            for(int v=1;v<=NV;v++)
                fprintf(fp,"v %.17g %.17g %.17g\n",pos[v][0],pos[v][1],pos[v][2]);
            for(int i=0;i<NF;i++)
                fprintf(fp,"f %d %d %d\n",F[i].a,F[i].b,F[i].c);
            fclose(fp);
        }
        if(pass) ok_count++; else {fail_count++; if(pA_success||pB_success) not_embedded++;}

        build_clear(); nets++;
    }
    fprintf(stderr,"polish: nets=%ld ok=%ld fail=%ld skipped=%ld "
            "(pA_ok=%ld pB_ok=%ld both_fail=%ld not_embedded=%ld no_input=%ld)\n",
            nets,ok_count,fail_count,skipped,pA_ok,pB_ok,both_fail,not_embedded,no_input);
    return 0;
}
