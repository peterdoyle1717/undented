/* ideal.c — Ideal horoball solver: CLERS name → ideal OBJ
 *
 * Computes horoball weights (horou) and flat layout (horoz).
 * Writes OBJ with:
 *   v1: (0, 0, 1)              — vertex at infinity
 *   v2: (0, 0, 1)              — boundary, u=1
 *   v3+: (hx, hy, u^2)         — horosphere positions + weight squared
 *
 * The hypersolver reads this OBJ and lifts to finite rho by:
 *   moving v1 to (0, 0, a=e^rho), dividing all z by a.
 *
 * Usage:   ./ideal outdir < prime/N.txt
 * Compile: cc -O3 -o ideal ideal.c -lm
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
    memset(bndry,0,(NV+2)*sizeof(int)); memset(int_idx,-1,(NV+2)*sizeof(int));
    int bndry_ring[MAXV],nb=cyclic_nbrs(1,bndry_ring);
    for(int i=0;i<nb;i++) bndry[bndry_ring[i]]=1;
    int n_int=0;
    for(int v=2;v<=NV;v++) if(!bndry[v]){int_idx[v]=n_int;interior[n_int++]=v;}
    #define U(v) (bndry[v]?1.0:xvec[int_idx[v]])
    for(int i=0;i<n_int;i++){int v=interior[i];ringlen[v]=cyclic_nbrs(v,ring[v]);}
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
static void horoz(double u[],double hx_out[],double hy_out[]){
    static double hx[MAXV+1],hy[MAXV+1]; static int placed[MAXV+1];
    static int qa[MAXV*2+8],qb[MAXV*2+8];
    memset(placed,0,(NV+2)*sizeof(int));
    placed[1]=placed[2]=placed[3]=1; hx[2]=0;hy[2]=0;hx[3]=1;hy[3]=0;
    int qh=0,qt=0; qa[qt]=3;qb[qt]=2;qt++;
    while(qh<qt){
        int a=qa[qh],b=qb[qh];qh++;
        int c=EM[a][b]; if(!c||c==1||placed[c]) continue;
        double dA=u[a-1]*u[c-1],dB=u[b-1]*u[c-1];
        thirdpoint(hx[a],hy[a],hx[b],hy[b],dA,dB,&hx[c],&hy[c]);
        placed[c]=1; qa[qt]=c;qb[qt]=b;qt++; qa[qt]=a;qb[qt]=c;qt++;
    }
    {int chg=1; while(chg){chg=0;
        for(int i=0;i<NF;i++){int a=F[i].a,b=F[i].b,c=F[i].c;
            if(a==1||b==1||c==1) continue;
            if(placed[a]&&placed[b]&&!placed[c]){thirdpoint(hx[a],hy[a],hx[b],hy[b],u[a-1]*u[c-1],u[b-1]*u[c-1],&hx[c],&hy[c]);placed[c]=1;chg=1;}
            else if(placed[b]&&placed[c]&&!placed[a]){thirdpoint(hx[b],hy[b],hx[c],hy[c],u[b-1]*u[a-1],u[c-1]*u[a-1],&hx[a],&hy[a]);placed[a]=1;chg=1;}
            else if(placed[a]&&placed[c]&&!placed[b]){thirdpoint(hx[c],hy[c],hx[a],hy[a],u[c-1]*u[b-1],u[a-1]*u[b-1],&hx[b],&hy[b]);placed[b]=1;chg=1;}
        }
    }}
    for(int v=1;v<=NV;v++){hx_out[v]=hx[v]; hy_out[v]=hy[v];}
}

/* ── main ────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv){
    if(argc<2){fprintf(stderr,"usage: ideal outdir < names.txt\n");return 1;}

    char outdir[4096];
    strncpy(outdir,argv[1],sizeof(outdir)-1); outdir[sizeof(outdir)-1]='\0';

    static char line[MAXCODE];
    static double u[MAXV], hx[MAXV+1], hy[MAXV+1];
    long nets=0;

    while(fgets(line,sizeof(line),stdin)){
        int ll=strlen(line);
        while(ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r')) line[--ll]='\0';
        if(!ll) continue;
        if(!decode(line)){fprintf(stderr,"decode failed: %s\n",line);continue;}
        build();

        horou(u); horoz(u, hx, hy);

        char path[4096];
        snprintf(path,sizeof(path),"%s/%s.obj",outdir,line);
        FILE *fp=fopen(path,"w");
        if(fp){
            /* v1: vertex at infinity */
            fprintf(fp,"v 0 0 1\n");
            /* v2: boundary vertex at origin, u=1 */
            fprintf(fp,"v 0 0 1\n");
            /* v3..NV: (hx, hy, u^2) */
            for(int v=3;v<=NV;v++)
                fprintf(fp,"v %.17g %.17g %.17g\n", hx[v], hy[v], u[v-1]*u[v-1]);
            for(int i=0;i<NF;i++)
                fprintf(fp,"f %d %d %d\n",F[i].a,F[i].b,F[i].c);
            fclose(fp);
        }
        build_clear();
        nets++;
    }
    fprintf(stderr,"ideal: %ld nets\n",nets);
    return 0;
}
