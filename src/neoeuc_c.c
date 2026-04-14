/* neoeuc_c.c — Euclidean equal-sphere packing solver for prime 6-nets
 *
 * Pipeline: horoball weights → horosphere positions → finite UHS homotopy
 * (n=12 rho-uniform steps, rho 3→0.25) → hyperbolic bisection until
 * blow-up residual ≤ 0.10 → Euclidean Newton with undented backtracking.
 *
 * Gauge: v1=(0,0,1/2), v2=(0,0,-1/2), v3=(√3/2,0,0).
 *
 * Output to stdout: per net, NV×3 float64 (doubles):
 *   v1..v3: exact reference triangle
 *   v4..NV: Euclidean Newton solution
 *   all NaN if net failed
 *
 * Usage:   ./neoeuc_c < prime/N.txt > euc/N.bin
 * Compile: cc -O3 -o neoeuc_c neoeuc_c.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV    400
#define MAXF    (2*MAXV+4)
#define MAXCODE (4*MAXV+8)
#define STKDEP  MAXCODE
#define MAXDQ   (MAXV*2+8)
#define MAXRING 12
#define MAXE    (3*MAXV)
#define MAXFN   (3*MAXV)

/* ── graph state ─────────────────────────────────────────────────────────── */
typedef struct { int a,b,c; } Face;
static int  NV, NF;
static Face F[MAXF];
static int  EM[MAXV+1][MAXV+1];
static int  DEG[MAXV+1];
static int  NBR[MAXV+1][MAXRING];
static int  NNBR[MAXV+1];
static short DU[MAXF*3], DW[MAXF*3]; static int ND;

static void nbr_add(int u, int w) {
    for (int i=0;i<NNBR[u];i++) if(NBR[u][i]==w) return;
    NBR[u][NNBR[u]++]=w;
}
static void build_clear(void) { for(int i=0;i<ND;i++) EM[DU[i]][DW[i]]=0; ND=0; }
static void build(void) {
    memset(DEG,0,(NV+2)*sizeof(int)); memset(NNBR,0,(NV+2)*sizeof(int)); ND=0;
    for(int i=0;i<NF;i++) {
        int a=F[i].a,b=F[i].b,c=F[i].c;
        EM[a][b]=c; DU[ND]=a; DW[ND]=b; ND++;
        EM[b][c]=a; DU[ND]=b; DW[ND]=c; ND++;
        EM[c][a]=b; DU[ND]=c; DW[ND]=a; ND++;
        DEG[a]++; DEG[b]++; DEG[c]++;
        nbr_add(a,b); nbr_add(b,a); nbr_add(b,c);
        nbr_add(c,b); nbr_add(a,c); nbr_add(c,a);
    }
}

/* ── CLERS decoder ───────────────────────────────────────────────────────── */
static int UF[MAXV*2+4];
static int uf_find(int x){while(UF[x]!=x){UF[x]=UF[UF[x]];x=UF[x];}return x;}
static void uf_link(int x,int y){x=uf_find(x);y=uf_find(y);if(x!=y)UF[x]=y;}
typedef struct{int data[MAXDQ];int h,t;}Deque;
static void dq_init(Deque*d){d->h=d->t=MAXDQ/2;}
static void dq_push_front(Deque*d,int v){d->data[--d->h]=v;}
static void dq_push_back(Deque*d,int v){d->data[d->t++]=v;}
static int  dq_pop_back(Deque*d){return d->data[--d->t];}
typedef struct{char tile;int a,b,c,phase;}Frame;
static Frame FSTK[STKDEP]; static Deque DSTK[MAXV]; static Face DTRIS[MAXF];

static int decode(const char*code){
    int n=strlen(code); if(!n)return 0;
    int maxv=n+4; for(int i=0;i<=maxv;i++) UF[i]=i;
    int ptr=0,V=2,ntris=0,nstk=0,ndq=0;
    char tile=code[ptr++]; V++;
    FSTK[nstk++]=(Frame){tile,1,2,V,0};
    while(nstk>0){
        Frame*f=&FSTK[nstk-1];
        char t=f->tile; int a=f->a,b=f->b,c=f->c,ph=f->phase;
        if(t=='E'){
            nstk--; dq_init(&DSTK[ndq]);
            dq_push_back(&DSTK[ndq],a); dq_push_back(&DSTK[ndq],c); dq_push_back(&DSTK[ndq],b);
            ndq++; DTRIS[ntris++]=(Face){a,b,c};
        }else if(t=='A'){
            if(!ph){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,c,b,V,0};}
            else{nstk--;dq_push_front(&DSTK[ndq-1],a);DTRIS[ntris++]=(Face){a,b,c};}
        }else if(t=='B'){
            if(!ph){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,a,c,V,0};}
            else{nstk--;dq_push_back(&DSTK[ndq-1],b);DTRIS[ntris++]=(Face){a,b,c};}
        }else if(t=='C'){
            if(!ph){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,a,c,V,0};}
            else{
                nstk--;Deque*d=&DSTK[ndq-1];
                dq_pop_back(d); int dd=dq_pop_back(d);
                dq_push_back(d,b); uf_link(dd,b); DTRIS[ntris++]=(Face){a,b,c};
            }
        }else{
            if(ph==0){f->phase=1;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,a,c,V,0};}
            else if(ph==1){f->phase=2;char nt=code[ptr++];V++;FSTK[nstk++]=(Frame){nt,c,b,V,0};}
            else{
                nstk--;Deque*dq2=&DSTK[--ndq],*dq1=&DSTK[ndq-1];
                dq_pop_back(dq1);
                for(int i=dq2->h;i<dq2->t;i++) dq_push_back(dq1,dq2->data[i]);
                DTRIS[ntris++]=(Face){a,b,c};
            }
        }
    }
    static int lbl[MAXV*2+4]; memset(lbl,0,(maxv+1)*sizeof(int)); int k=0;
    for(int i=ntris-1;i>=0;i--){
        int a=uf_find(DTRIS[i].a),b=uf_find(DTRIS[i].b),c=uf_find(DTRIS[i].c);
        if(!lbl[a])lbl[a]=++k; if(!lbl[b])lbl[b]=++k; if(!lbl[c])lbl[c]=++k;
        F[ntris-1-i]=(Face){lbl[a],lbl[b],lbl[c]};
    }
    NF=ntris; NV=k; return ntris;
}

/* ── cyclic neighbors ────────────────────────────────────────────────────── */
static int cyclic_nbrs(int v, int ring[]) {
    if(!NNBR[v]) return 0;
    int start=NBR[v][0]; ring[0]=start; int k=1,cur=start;
    for(;;){ int nxt=EM[v][cur]; if(nxt==start) break; ring[k++]=nxt; cur=nxt; }
    return k;
}

/* ── petal / petal_grad ──────────────────────────────────────────────────── */
static double petal(double ui,double uj,double uk){
    double a=ui*uj,b=ui*uk,c=uj*uk;
    double ct=(a*a+b*b-c*c)/(2.0*a*b);
    if(ct>1)ct=1; if(ct<-1)ct=-1; return acos(ct);
}
static void petal_grad(double ui,double uj,double uk,double*dui,double*duj,double*duk){
    double a=ui*uj,b=ui*uk,c=uj*uk;
    double ct=(a*a+b*b-c*c)/(2.0*a*b);
    if(ct>1)ct=1; if(ct<-1)ct=-1;
    double st=sqrt(1.0-ct*ct); if(st<1e-15){*dui=*duj=*duk=0;return;}
    double s=-1.0/st;
    *dui=s*(uj*uk/(ui*ui*ui));
    *duj=s*(1.0/(2*uk)-uk/(2*uj*uj)-uk/(2*ui*ui));
    *duk=s*(1.0/(2*uj)-uj/(2*uk*uk)-uj/(2*ui*ui));
}

/* ── horou Newton solver ─────────────────────────────────────────────────── */
static double Jmat[MAXV][MAXV], Fvec[MAXV], dxvec[MAXV];
static int lu_solve(double J[][MAXV],double b[],int n){
    for(int col=0;col<n;col++){
        int piv=col; double best=fabs(J[col][col]);
        for(int row=col+1;row<n;row++) if(fabs(J[row][col])>best){best=fabs(J[row][col]);piv=row;}
        if(best<1e-14) return -1;
        if(piv!=col){
            for(int k=col;k<n;k++){double t=J[col][k];J[col][k]=J[piv][k];J[piv][k]=t;}
            {double t=b[col];b[col]=b[piv];b[piv]=t;}
        }
        double inv=1.0/J[col][col];
        for(int row=col+1;row<n;row++){
            double fac=J[row][col]*inv;
            for(int k=col;k<n;k++) J[row][k]-=fac*J[col][k];
            b[row]-=fac*b[col];
        }
    }
    for(int i=n-1;i>=0;i--){double s=b[i];for(int j=i+1;j<n;j++)s-=J[i][j]*b[j];b[i]=s/J[i][i];}
    return 0;
}

static void horou(double u_out[]){
    static int bndry[MAXV+1],int_idx[MAXV+1],interior[MAXV];
    static int ring[MAXV+1][MAXRING],ringlen[MAXV+1];
    static double xvec[MAXV];
    static int ff_a[MAXF],ff_b[MAXF],ff_c[MAXF];
    memset(bndry,0,(NV+2)*sizeof(int)); memset(int_idx,-1,(NV+2)*sizeof(int));
    int bndry_ring[MAXV],nb=cyclic_nbrs(1,bndry_ring);
    for(int i=0;i<nb;i++) bndry[bndry_ring[i]]=1;
    int n_int=0;
    for(int v=2;v<=NV;v++) if(!bndry[v]){int_idx[v]=n_int;interior[n_int++]=v;}
    #define U(v) (bndry[v]?1.0:xvec[int_idx[v]])
    for(int i=0;i<n_int;i++){int v=interior[i];ringlen[v]=cyclic_nbrs(v,ring[v]);}
    int nff=0;
    for(int i=0;i<NF;i++){int a=F[i].a,b=F[i].b,c=F[i].c;
        if(a!=1&&b!=1&&c!=1){ff_a[nff]=a;ff_b[nff]=b;ff_c[nff]=c;nff++;}}
    if(!n_int) goto done;
    for(int i=0;i<n_int;i++) xvec[i]=1.0;
    for(int iter=0;iter<200;iter++){
        double res=0;
        for(int i=0;i<n_int;i++){
            int v=interior[i];double s=0;
            for(int j=0;j<ringlen[v];j++) s+=petal(xvec[i],U(ring[v][j]),U(ring[v][(j+1)%ringlen[v]]));
            Fvec[i]=s-2*M_PI; double af=fabs(Fvec[i]); if(af>res)res=af;
        }
        if(res<1e-10) break;
        memset(Jmat,0,sizeof(double)*n_int*MAXV);
        for(int i=0;i<n_int;i++){
            int v=interior[i];
            for(int j=0;j<ringlen[v];j++){
                int vj=ring[v][j],vk=ring[v][(j+1)%ringlen[v]];
                double dui,duj,duk; petal_grad(xvec[i],U(vj),U(vk),&dui,&duj,&duk);
                Jmat[i][i]+=dui;
                if(int_idx[vj]>=0) Jmat[i][int_idx[vj]]+=duj;
                if(int_idx[vk]>=0) Jmat[i][int_idx[vk]]+=duk;
            }
        }
        for(int i=0;i<n_int;i++) dxvec[i]=-Fvec[i];
        if(lu_solve(Jmat,dxvec,n_int)<0) break;
        double step=1.0;
        for(int bt=0;bt<60;bt++,step*=0.5){
            int ok=1; for(int i=0;i<n_int;i++) if(xvec[i]+step*dxvec[i]<=0){ok=0;break;}
            if(!ok) continue;
            double res2=0;
            for(int i=0;i<n_int;i++){
                int v=interior[i];double s=0;
                for(int j=0;j<ringlen[v];j++){
                    double uj=U(ring[v][j]),uk2=U(ring[v][(j+1)%ringlen[v]]);
                    double ui2=xvec[i]+step*dxvec[i];
                    if(int_idx[ring[v][j]]>=0) uj=xvec[int_idx[ring[v][j]]]+step*dxvec[int_idx[ring[v][j]]];
                    if(int_idx[ring[v][(j+1)%ringlen[v]]]>=0) uk2=xvec[int_idx[ring[v][(j+1)%ringlen[v]]]]+step*dxvec[int_idx[ring[v][(j+1)%ringlen[v]]]];
                    s+=petal(ui2,uj,uk2);
                }
                double af=fabs(s-2*M_PI); if(af>res2)res2=af;
            }
            if(res2<res) break;
        }
        for(int i=0;i<n_int;i++) xvec[i]+=step*dxvec[i];
    }
done:
    u_out[0]=NAN;
    for(int v=2;v<=NV;v++) u_out[v-1]=bndry[v]?1.0:xvec[int_idx[v]];
    #undef U
}

/* ── thirdpoint / horoz ──────────────────────────────────────────────────── */
static void thirdpoint(double xa,double ya,double xb,double yb,
                       double dA,double dB,double*xc,double*yc){
    double dx=xb-xa,dy=yb-ya,L=sqrt(dx*dx+dy*dy);
    double d0=dA/L,d1=dB/L,xn=(1+d0*d0-d1*d1)/2.0;
    double yn2=d0*d0-xn*xn,yn=sqrt(yn2>0?yn2:0);
    *xc=xa+dx*xn+dy*yn; *yc=ya+dy*xn-dx*yn;
}
static void horoz(double u[],double out[]){
    static double hx[MAXV+1],hy[MAXV+1]; static int placed[MAXV+1];
    static int qa[MAXV*2+8],qb[MAXV*2+8];
    memset(placed,0,(NV+2)*sizeof(int));
    placed[1]=placed[2]=placed[3]=1; hx[2]=0;hy[2]=0;hx[3]=1;hy[3]=0;
    int qh=0,qt=0; qa[qt]=3;qb[qt]=2;qt++;
    while(qh<qt){
        int a=qa[qh],b=qb[qh];qh++;
        int c=EM[a][b]; if(!c||c==1||placed[c]) continue;
        double dA=(double)u[a-1]*(double)u[c-1],dB=(double)u[b-1]*(double)u[c-1];
        thirdpoint(hx[a],hy[a],hx[b],hy[b],dA,dB,&hx[c],&hy[c]);
        placed[c]=1; qa[qt]=c;qb[qt]=b;qt++; qa[qt]=a;qb[qt]=c;qt++;
    }
    {int chg=1; while(chg){chg=0;
        for(int i=0;i<NF;i++){int a=F[i].a,b=F[i].b,c=F[i].c;
            if(a==1||b==1||c==1) continue;
            if(placed[a]&&placed[b]&&!placed[c]){thirdpoint(hx[a],hy[a],hx[b],hy[b],(double)u[a-1]*(double)u[c-1],(double)u[b-1]*(double)u[c-1],&hx[c],&hy[c]);placed[c]=1;chg=1;}
            else if(placed[b]&&placed[c]&&!placed[a]){thirdpoint(hx[b],hy[b],hx[c],hy[c],(double)u[b-1]*(double)u[a-1],(double)u[c-1]*(double)u[a-1],&hx[a],&hy[a]);placed[a]=1;chg=1;}
            else if(placed[a]&&placed[c]&&!placed[b]){thirdpoint(hx[c],hy[c],hx[a],hy[a],(double)u[c-1]*(double)u[b-1],(double)u[a-1]*(double)u[b-1],&hx[b],&hy[b]);placed[b]=1;chg=1;}
        }
    }}
    out[0]=out[1]=out[2]=NAN;
    for(int v=2;v<=NV;v++){out[3*(v-1)]=u[v-1];out[3*(v-1)+1]=hx[v];out[3*(v-1)+2]=hy[v];}
}

/* ── reference triangle & edge list ─────────────────────────────────────── */
static double fp_x[4],fp_y[4],fp_t[4];
static void ref_triangle(double A){
    double A2=A*A,A4=A2*A2,A6=A4*A2;
    double Delta=A6-A4+A2-1.0;
    fp_x[1]=0;fp_y[1]=0;fp_t[1]=A;
    fp_x[2]=0;fp_y[2]=0;fp_t[2]=1.0/A;
    fp_t[3]=A*(A4-1.0)/Delta;
    double prod=(A2-1.0)*(A2-1.0)*(A2-1.0)*(A6-1.0);
    fp_x[3]=sqrt(prod>0?prod:0)/Delta; fp_y[3]=0;
}

static int eu[MAXE],ev[MAXE],n_edges;
static void collect_edges(void){
    n_edges=0;
    for(int i=1;i<=NV;i++) for(int j=i+1;j<=NV;j++){
        if(!EM[i][j]) continue;
        if(i<=3&&j<=3) continue;
        eu[n_edges]=i;ev[n_edges]=j;n_edges++;
    }
}

/* ── shared Newton matrix ────────────────────────────────────────────────── */
static double NJmat[MAXFN][MAXFN],NFvec[MAXFN];
static int lu_solve_n(int n){
    for(int col=0;col<n;col++){
        int piv=col;double best=fabs(NJmat[col][col]);
        for(int row=col+1;row<n;row++) if(fabs(NJmat[row][col])>best){best=fabs(NJmat[row][col]);piv=row;}
        if(best<1e-14) return -1;
        if(piv!=col){
            for(int k=col;k<n;k++){double t=NJmat[col][k];NJmat[col][k]=NJmat[piv][k];NJmat[piv][k]=t;}
            {double t=NFvec[col];NFvec[col]=NFvec[piv];NFvec[piv]=t;}
        }
        double inv=1.0/NJmat[col][col];
        for(int row=col+1;row<n;row++){
            double fac=NJmat[row][col]*inv;
            for(int k=col;k<n;k++) NJmat[row][k]-=fac*NJmat[col][k];
            NFvec[row]-=fac*NFvec[col];
        }
    }
    for(int i=n-1;i>=0;i--){double s=NFvec[i];for(int j=i+1;j<n;j++)s-=NJmat[i][j]*NFvec[j];NFvec[i]=s/NJmat[i][i];}
    return 0;
}

/* ── finite UHS Newton ───────────────────────────────────────────────────── */
/* unknowns g_xvec[3*(v-4)+{0,1,2}] = x,y,t for v=4..NV */
static double g_xvec[MAXFN];
#define GX(v) ((v)<=3?fp_x[v]:g_xvec[3*((v)-4)])
#define GY(v) ((v)<=3?fp_y[v]:g_xvec[3*((v)-4)+1])
#define GT(v) ((v)<=3?fp_t[v]:g_xvec[3*((v)-4)+2])

/* ── geodesic predictor ──────────────────────────────────────────────────── */
/* Radial contraction in the Klein model about basepoint j=(0,0,1) in UHS.
   c = log(b)/log(a) is the contraction exponent. Updates g_xvec in-place. */
static void predict_step(double a_old, double a_new){
    if(a_old<=1.0||a_new<=1.0||a_old==a_new) return;
    double c=log(a_new)/log(a_old);
    for(int v=4;v<=NV;v++){
        double x=g_xvec[3*(v-4)],y=g_xvec[3*(v-4)+1],t=g_xvec[3*(v-4)+2];
        /* UHS -> Klein */
        double r2=x*x+y*y, sum2=r2+t*t;
        double den=sum2+1.0;
        double kx=2.0*x/den, ky=2.0*y/den, kz=(sum2-1.0)/den;
        /* radial rescaling */
        double rho2=kx*kx+ky*ky+kz*kz, rho=sqrt(rho2);
        if(rho<1e-15) continue;
        double rho_new=tanh(c*atanh(rho));
        double scale=rho_new/rho;
        double kxn=scale*kx, kyn=scale*ky, kzn=scale*kz;
        /* Klein -> UHS */
        double s2=1.0-(kxn*kxn+kyn*kyn+kzn*kzn);
        double d=1.0-kzn;
        g_xvec[3*(v-4)  ]=kxn/d;
        g_xvec[3*(v-4)+1]=kyn/d;
        g_xvec[3*(v-4)+2]=sqrt(s2>0?s2:0)/d;
    }
}

static int finite_newton(double a, int *iters_out, double *res_out){
    ref_triangle(a);
    double kf=(a-1.0/a)*(a-1.0/a);
    int n=3*(NV-3); double res=1e30; int iter; int lufail=0;
    for(iter=0;iter<100;iter++){
        res=0;
        for(int k=0;k<n_edges;k++){
            int i=eu[k],j=ev[k];
            double xi=GX(i),yi=GY(i),ti=GT(i),xj=GX(j),yj=GY(j),tj=GT(j);
            double dx=xi-xj,dy=yi-yj,dt=ti-tj;
            NFvec[k]=dx*dx+dy*dy+dt*dt-ti*tj*kf;
            double af=fabs(NFvec[k]); if(af>res) res=af;
        }
        if(res<1e-7) break;
        memset(NJmat,0,sizeof(double)*(size_t)n*MAXFN);
        for(int k=0;k<n_edges;k++){
            int i=eu[k],j=ev[k];
            double xi=GX(i),yi=GY(i),ti=GT(i),xj=GX(j),yj=GY(j),tj=GT(j);
            if(i>=4){int ci=3*(i-4);NJmat[k][ci]+=2*(xi-xj);NJmat[k][ci+1]+=2*(yi-yj);NJmat[k][ci+2]+=2*(ti-tj)-tj*kf;}
            if(j>=4){int cj=3*(j-4);NJmat[k][cj]-=2*(xi-xj);NJmat[k][cj+1]-=2*(yi-yj);NJmat[k][cj+2]+=-2*(ti-tj)-ti*kf;}
        }
        for(int k=0;k<n;k++) NFvec[k]=-NFvec[k];
        if(lu_solve_n(n)<0){lufail=1;break;}
        double step=1.0;
        for(int bt=0;bt<40;bt++,step*=0.5){
            int t_ok=1;
            for(int v=4;v<=NV;v++) if(GT(v)+step*NFvec[3*(v-4)+2]<=0){t_ok=0;break;}
            if(!t_ok) continue;
            double res2=0;
            for(int k=0;k<n_edges;k++){
                int i=eu[k],j=ev[k];
                double xi=(i>=4?GX(i)+step*NFvec[3*(i-4)]:GX(i));
                double yi=(i>=4?GY(i)+step*NFvec[3*(i-4)+1]:GY(i));
                double ti=(i>=4?GT(i)+step*NFvec[3*(i-4)+2]:GT(i));
                double xj=(j>=4?GX(j)+step*NFvec[3*(j-4)]:GX(j));
                double yj=(j>=4?GY(j)+step*NFvec[3*(j-4)+1]:GY(j));
                double tj=(j>=4?GT(j)+step*NFvec[3*(j-4)+2]:GT(j));
                double dx=xi-xj,dy=yi-yj,dt=ti-tj,f=dx*dx+dy*dy+dt*dt-ti*tj*kf;
                double af=fabs(f); if(af>res2) res2=af;
            }
            if(res2<res) break;
        }
        for(int v=4;v<=NV;v++){g_xvec[3*(v-4)]+=step*NFvec[3*(v-4)];g_xvec[3*(v-4)+1]+=step*NFvec[3*(v-4)+1];g_xvec[3*(v-4)+2]+=step*NFvec[3*(v-4)+2];}
    }
    *iters_out=iter; *res_out=res;
    return (!lufail && res<1e-5)?1:0;
}

/* ── Euclidean Newton ────────────────────────────────────────────────────── */
static double e_xvec[MAXFN];  /* unknowns: X,Y,Z for v=4..NV */
static const double EC[4][3]={{0,0,0},{0,0,0.5},{0,0,-0.5},{0.8660254037844387,0,0}};
#define EX(v) ((v)<=3?EC[v][0]:e_xvec[3*((v)-4)])
#define EY(v) ((v)<=3?EC[v][1]:e_xvec[3*((v)-4)+1])
#define EZ(v) ((v)<=3?EC[v][2]:e_xvec[3*((v)-4)+2])

static int undented_check(double coords[][3]);  /* forward decl */
static int euclid_newton(int *iters_out, double *res_out, int *bt_total_out){
    int n=3*(NV-3); double res=1e30; int iter; int lufail=0; int bt_total=0;
    for(iter=0;iter<50;iter++){
        res=0;
        for(int k=0;k<n_edges;k++){
            int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            NFvec[k]=dx*dx+dy*dy+dz*dz-1.0;
            double af=fabs(NFvec[k]); if(af>res) res=af;
        }
        if(res<1e-10) break;
        memset(NJmat,0,sizeof(double)*(size_t)n*MAXFN);
        for(int k=0;k<n_edges;k++){
            int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            if(i>=4){int ci=3*(i-4);NJmat[k][ci]+=2*dx;NJmat[k][ci+1]+=2*dy;NJmat[k][ci+2]+=2*dz;}
            if(j>=4){int cj=3*(j-4);NJmat[k][cj]-=2*dx;NJmat[k][cj+1]-=2*dy;NJmat[k][cj+2]-=2*dz;}
        }
        for(int k=0;k<n;k++) NFvec[k]=-NFvec[k];
        if(lu_solve_n(n)<0){lufail=1;break;}
        /* backtrack: halve step if residual doesn't improve or dent introduced */
        double step=1.0;
        static double ctmp[MAXV+1][3];
        for(int bt=0;bt<60;bt++,step*=0.5,bt_total++){
            double res2=0;
            for(int k=0;k<n_edges;k++){
                int i=eu[k],j=ev[k];
                double dx=(i>=4?EX(i)+step*NFvec[3*(i-4)]:EX(i))-(j>=4?EX(j)+step*NFvec[3*(j-4)]:EX(j));
                double dy=(i>=4?EY(i)+step*NFvec[3*(i-4)+1]:EY(i))-(j>=4?EY(j)+step*NFvec[3*(j-4)+1]:EY(j));
                double dz=(i>=4?EZ(i)+step*NFvec[3*(i-4)+2]:EZ(i))-(j>=4?EZ(j)+step*NFvec[3*(j-4)+2]:EZ(j));
                double af=fabs(dx*dx+dy*dy+dz*dz-1.0); if(af>res2) res2=af;
            }
            if(res2>=res) continue;
            /* residual improved — check undentedness of tentative positions */
            for(int v=1;v<=3;v++){ctmp[v][0]=EC[v][0];ctmp[v][1]=EC[v][1];ctmp[v][2]=EC[v][2];}
            for(int v=4;v<=NV;v++){
                ctmp[v][0]=e_xvec[3*(v-4)  ]+step*NFvec[3*(v-4)  ];
                ctmp[v][1]=e_xvec[3*(v-4)+1]+step*NFvec[3*(v-4)+1];
                ctmp[v][2]=e_xvec[3*(v-4)+2]+step*NFvec[3*(v-4)+2];
            }
            if(undented_check(ctmp)) break;
        }
        for(int v=4;v<=NV;v++){e_xvec[3*(v-4)]+=step*NFvec[3*(v-4)];e_xvec[3*(v-4)+1]+=step*NFvec[3*(v-4)+1];e_xvec[3*(v-4)+2]+=step*NFvec[3*(v-4)+2];}
    }
    *iters_out=iter; *res_out=res; *bt_total_out=bt_total;
    return (!lufail && res<1e-8)?1:0;
}

/* ── Klein model & undentedness ──────────────────────────────────────────── */
static void uhs_to_klein(double x,double y,double t,double*kx,double*ky,double*kz){
    double r2=x*x+y*y+t*t,d=r2+1;
    *kx=2*x/d; *ky=2*y/d; *kz=(r2-1)/d;
}

static double signed_sph_angle(double ax,double ay,double az,
                                double bx,double by,double bz,
                                double cx,double cy,double cz){
    double acx=ay*cz-az*cy,acy=az*cx-ax*cz,acz=ax*cy-ay*cx;
    double num=bx*acx+by*acy+bz*acz;
    double den=(ax*bx+ay*by+az*bz)*(bx*cx+by*cy+bz*cz)-(ax*cx+ay*cy+az*cz);
    return atan2(num,den);
}

/* coords[v][0..2] = 3D position of vertex v (1..NV).
   Returns 1 if undented, 0 if any vertex has negative turning. */
static int undented_check(double coords[][3]){
    int ring[MAXRING]; double min_t=1e30;
    for(int v=1;v<=NV;v++){
        int k=cyclic_nbrs(v,ring); if(k<3) continue;
        double dx[MAXRING],dy[MAXRING],dz[MAXRING];
        for(int i=0;i<k;i++){
            int w=ring[i];
            double ex=coords[w][0]-coords[v][0],ey=coords[w][1]-coords[v][1],ez=coords[w][2]-coords[v][2];
            double L=sqrt(ex*ex+ey*ey+ez*ez); if(L<1e-30) L=1e-30;
            dx[i]=ex/L;dy[i]=ey/L;dz[i]=ez/L;
        }
        double tv=0;
        for(int i=0;i<k;i++)
            tv+=signed_sph_angle(dx[(i-1+k)%k],dy[(i-1+k)%k],dz[(i-1+k)%k],
                                  dx[i],dy[i],dz[i],
                                  dx[(i+1)%k],dy[(i+1)%k],dz[(i+1)%k]);
        if(tv<min_t) min_t=tv;
    }
    return (min_t>=0.0)?1:0;
}

/* ── angle → a ───────────────────────────────────────────────────────────── */
static double a_from_alpha(double alpha){
    if(alpha<1e-15) return 1e15;
    double c=cos(alpha),disc=2*c-1;
    if(disc<0) disc=0;
    return sqrt((c+sqrt(disc))/(1.0-c));
}

/* ── output helpers ──────────────────────────────────────────────────────── */
/* Write NV×3 float32 blow-up coords to stdout, fill coords[][] for undented check */
static void write_frame_hyp(double a, float *out_buf, double coords[][3]){
    double rho=log(a);
    /* v1,v2,v3 from reference triangle */
    ref_triangle(a);
    double pts[4][3];
    pts[1][0]=pts[1][1]=0; pts[1][2]=0.5;
    pts[2][0]=pts[2][1]=0; pts[2][2]=-0.5;
    pts[3][0]=fp_x[3]/(2*rho); pts[3][1]=0; pts[3][2]=log(fp_t[3])/(2*rho);
    for(int v=1;v<=3;v++){out_buf[3*(v-1)  ]=(float)pts[v][0];out_buf[3*(v-1)+1]=(float)pts[v][1];out_buf[3*(v-1)+2]=(float)pts[v][2];}
    if(coords){for(int v=1;v<=3;v++){coords[v][0]=pts[v][0];coords[v][1]=pts[v][1];coords[v][2]=pts[v][2];}}
    for(int v=4;v<=NV;v++){
        double x=g_xvec[3*(v-4)],y=g_xvec[3*(v-4)+1],t=g_xvec[3*(v-4)+2];
        double X=x/(2*rho),Y=y/(2*rho),Z=(t>0?log(t)/(2*rho):0);
        out_buf[3*(v-1)  ]=(float)X; out_buf[3*(v-1)+1]=(float)Y; out_buf[3*(v-1)+2]=(float)Z;
        if(coords){coords[v][0]=X;coords[v][1]=Y;coords[v][2]=Z;}
    }
}

static void write_frame_euc(float *out_buf, double coords[][3]){
    for(int v=1;v<=NV;v++){
        double x=EX(v),y=EY(v),z=EZ(v);
        if(out_buf){out_buf[3*(v-1)  ]=(float)x; out_buf[3*(v-1)+1]=(float)y; out_buf[3*(v-1)+2]=(float)z;}
        if(coords){coords[v][0]=x;coords[v][1]=y;coords[v][2]=z;}
    }
}

/* ── main ────────────────────────────────────────────────────────────────── */

int main(int argc, char **argv){
    static char line[MAXCODE];
    static double u[MAXV], hz[MAXV*3];
    static double out_buf[MAXV*3];
    static double coords[MAXV+1][3];
    static char outdir[4096];

    outdir[0]='\0';
    if(argc>=2) strncpy(outdir,argv[1],sizeof(outdir)-1);

    int n = 12;   /* rho schedule steps (bisection does the rest) */
    long nets=0, net_failed=0;

    while(fgets(line,sizeof(line),stdin)){
        int ll=strlen(line);
        while(ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r')) line[--ll]='\0';
        if(!ll) continue;
        if(!decode(line)){fprintf(stderr,"decode failed: %s\n",line);continue;}
        build(); collect_edges();

        /* horoball → horosphere */
        horou(u); horoz(u,hz);

        /* lift to finite UHS at rho=3 */
        double rho1=3.0, a1=exp(rho1);
        ref_triangle(a1);
        for(int v=4;v<=NV;v++){
            g_xvec[3*(v-4)  ]=(double)hz[3*(v-1)+1];
            g_xvec[3*(v-4)+1]=(double)hz[3*(v-1)+2];
            g_xvec[3*(v-4)+2]=(double)(u[v-1]*u[v-1])/a1;
        }

        int net_ok=1;

        /* finite UHS homotopy: rho = 3*(n+1-fr)/n, fr=1..n */
        for(int fr=1;fr<=n;fr++){
            double rho=3.0*(n+1-fr)/n;
            double a=exp(rho);
            if(fr>=2){ double rho_prev=3.0*(n+2-fr)/n; predict_step(exp(rho_prev),a); }
            int iters; double res;
            int ok=finite_newton(a,&iters,&res);
            if(!ok) net_ok=0;
        }

        /* bisect rho gap toward 0 until blow-up residual <= 0.10 */
        {
            double rho_cur=3.0/n;
            double rho_gap=rho_cur;
            double init_res=1e30; int n_bisect=0;
            for(int ref=0;ref<40;ref++){
                double a_cur=exp(rho_cur);
                ref_triangle(a_cur);
                for(int v=4;v<=NV;v++){
                    double x=g_xvec[3*(v-4)],y=g_xvec[3*(v-4)+1],t=g_xvec[3*(v-4)+2];
                    e_xvec[3*(v-4)  ]=x/(2*rho_cur);
                    e_xvec[3*(v-4)+1]=y/(2*rho_cur);
                    e_xvec[3*(v-4)+2]=(t>0?log(t)/(2*rho_cur):0);
                }
                init_res=0;
                for(int k=0;k<n_edges;k++){
                    int i=eu[k],j=ev[k];
                    double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
                    double af=fabs(dx*dx+dy*dy+dz*dz-1.0); if(af>init_res) init_res=af;
                }
                if(init_res<=0.10) break;
                double a_prev=a_cur;
                rho_gap/=2.0; rho_cur-=rho_gap;
                predict_step(a_prev,exp(rho_cur));
                int iw; double rw;
                int okw=finite_newton(exp(rho_cur),&iw,&rw);
                n_bisect++;
                if(!okw) break;
            }
            int iters; double res; int bt_total;
            int ok=euclid_newton(&iters,&res,&bt_total);
            /* polish: Newton with backtracking to avoid worsening residual */
            if(ok){
                for(int pi=0;pi<20&&res>1e-14;pi++){
                    res=0;
                    for(int k=0;k<n_edges;k++){
                        int i=eu[k],j=ev[k];
                        double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
                        NFvec[k]=dx*dx+dy*dy+dz*dz-1.0;
                        double af=fabs(NFvec[k]);if(af>res)res=af;
                    }
                    if(res<1e-14) break;
                    memset(NJmat,0,sizeof(double)*(size_t)n*MAXFN);
                    for(int k=0;k<n_edges;k++){
                        int i=eu[k],j=ev[k];
                        double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
                        if(i>=4){int ci=3*(i-4);NJmat[k][ci]+=2*dx;NJmat[k][ci+1]+=2*dy;NJmat[k][ci+2]+=2*dz;}
                        if(j>=4){int cj=3*(j-4);NJmat[k][cj]-=2*dx;NJmat[k][cj+1]-=2*dy;NJmat[k][cj+2]-=2*dz;}
                    }
                    for(int k=0;k<n;k++) NFvec[k]=-NFvec[k];
                    if(lu_solve_n(n)<0) break;
                    double step=1.0;
                    for(int bt=0;bt<20;bt++,step*=0.5){
                        double res2=0;
                        for(int k=0;k<n_edges;k++){
                            int i=eu[k],j=ev[k];
                            double dx=(i>=4?EX(i)+step*NFvec[3*(i-4)]:EX(i))-(j>=4?EX(j)+step*NFvec[3*(j-4)]:EX(j));
                            double dy=(i>=4?EY(i)+step*NFvec[3*(i-4)+1]:EY(i))-(j>=4?EY(j)+step*NFvec[3*(j-4)+1]:EY(j));
                            double dz=(i>=4?EZ(i)+step*NFvec[3*(i-4)+2]:EZ(i))-(j>=4?EZ(j)+step*NFvec[3*(j-4)+2]:EZ(j));
                            double af=fabs(dx*dx+dy*dy+dz*dz-1.0);if(af>res2)res2=af;
                        }
                        if(res2<res) break;
                    }
                    for(int v=4;v<=NV;v++){e_xvec[3*(v-4)]+=step*NFvec[3*(v-4)];e_xvec[3*(v-4)+1]+=step*NFvec[3*(v-4)+1];e_xvec[3*(v-4)+2]+=step*NFvec[3*(v-4)+2];}
                }
            }
            write_frame_euc(NULL, coords);
            int und=undented_check(coords);
            if(!ok||!und) net_ok=0;
            if(!ok||!und) fprintf(stderr,"FAIL net %ld  bisect=%d  init_res=%.2e  iters=%d  bt=%d  res=%.2e  %s  %s\n",
                    nets, n_bisect, init_res, iters, bt_total, res,
                    ok?"CONVERGED":"FAILED", und?"UNDENTED":"DENTED");
        }

        /* write OBJ or .failed */
        if(outdir[0]){
            char path[4096];
            if(net_ok){
                snprintf(path,sizeof(path),"%s/%s.obj",outdir,line);
                FILE *fp=fopen(path,"w");
                if(fp){
                    fprintf(fp,"v %.17g %.17g %.17g\n",EC[1][0],EC[1][1],EC[1][2]);
                    fprintf(fp,"v %.17g %.17g %.17g\n",EC[2][0],EC[2][1],EC[2][2]);
                    fprintf(fp,"v %.17g %.17g %.17g\n",EC[3][0],EC[3][1],EC[3][2]);
                    for(int v=4;v<=NV;v++)
                        fprintf(fp,"v %.17g %.17g %.17g\n",e_xvec[3*(v-4)],e_xvec[3*(v-4)+1],e_xvec[3*(v-4)+2]);
                    for(int i=0;i<NF;i++)
                        fprintf(fp,"f %d %d %d\n",F[i].a,F[i].b,F[i].c);
                    fclose(fp);
                }
            } else {
                snprintf(path,sizeof(path),"%s/%s.failed",outdir,line);
                fclose(fopen(path,"w"));
                net_failed++;
            }
        } else {
            /* legacy binary mode: NV*3 doubles to stdout */
            if(net_ok){
                out_buf[0]=EC[1][0]; out_buf[1]=EC[1][1]; out_buf[2]=EC[1][2];
                out_buf[3]=EC[2][0]; out_buf[4]=EC[2][1]; out_buf[5]=EC[2][2];
                out_buf[6]=EC[3][0]; out_buf[7]=EC[3][1]; out_buf[8]=EC[3][2];
                for(int v=4;v<=NV;v++){
                    out_buf[3*(v-1)  ]=e_xvec[3*(v-4)];
                    out_buf[3*(v-1)+1]=e_xvec[3*(v-4)+1];
                    out_buf[3*(v-1)+2]=e_xvec[3*(v-4)+2];
                }
            } else {
                for(int i=0;i<NV*3;i++) out_buf[i]=NAN;
                net_failed++;
            }
            fwrite(out_buf,sizeof(double),NV*3,stdout);
        }

        build_clear();
        nets++;
    }
    fprintf(stderr,"nets=%ld  net_failed=%ld\n",nets,net_failed);
    return 0;
}
