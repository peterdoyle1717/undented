/* euclidsolve.c — UHS last frame → Euclidean OBJ
 *
 * Reads .uhs files (from horodump), takes the last frame,
 * blows up to Euclidean coordinates, optionally polishes
 * with Newton (undented backtracking).
 *
 * Blow-up: x/(2ρ), y/(2ρ), log(t)/(2ρ)
 * Gauge: v1=(0,0,½), v2=(0,0,-½), v3=(√3/2,0,0)
 *
 * Usage:   ./euclidsolve [-polish] indir outdir < names.txt
 * Compile: cc -O3 -o euclidsolve euclidsolve.c -lm
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
static int  NBR[MAXV+1][MAXRING], NNBR[MAXV+1];
static short DU[MAXF*3], DW[MAXF*3]; static int ND;
static int eu[MAXE], ev[MAXE], n_edges;

static void nbr_add(int u, int w){for(int i=0;i<NNBR[u];i++)if(NBR[u][i]==w)return;NBR[u][NNBR[u]++]=w;}
static void build_clear(void){for(int i=0;i<ND;i++)EM[DU[i]][DW[i]]=0;ND=0;}
static void build(void){
    memset(EM,0,(NV+2)*(NV+2)*sizeof(int));memset(NNBR,0,(NV+2)*sizeof(int));ND=0;
    for(int i=0;i<NF;i++){int a=F[i].a,b=F[i].b,c=F[i].c;
        EM[a][b]=c;DU[ND]=a;DW[ND]=b;ND++;EM[b][c]=a;DU[ND]=b;DW[ND]=c;ND++;
        EM[c][a]=b;DU[ND]=c;DW[ND]=a;ND++;
        nbr_add(a,b);nbr_add(b,a);nbr_add(b,c);nbr_add(c,b);nbr_add(a,c);nbr_add(c,a);}
}
static void collect_edges(void){
    n_edges=0;
    for(int i=1;i<=NV;i++)for(int j=i+1;j<=NV;j++){
        if(!EM[i][j])continue;if(i<=3&&j<=3)continue;
        eu[n_edges]=i;ev[n_edges]=j;n_edges++;}
}

/* ── CLERS decoder ───────────────────────────────────────────────────────── */
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
static int cyclic_nbrs(int v,int ring[]){if(!NNBR[v])return 0;int start=NBR[v][0];ring[0]=start;int k=1,cur=start;for(;;){int nxt=EM[v][cur];if(nxt==start)break;ring[k++]=nxt;cur=nxt;}return k;}

/* ── Euclidean coords ────────────────────────────────────────────────────── */
static double e_xvec[MAXFN];
static const double EC[4][3]={{0,0,0},{0,0,0.5},{0,0,-0.5},{0.8660254037844387,0,0}};
#define EX(v) ((v)<=3?EC[v][0]:e_xvec[3*((v)-4)])
#define EY(v) ((v)<=3?EC[v][1]:e_xvec[3*((v)-4)+1])
#define EZ(v) ((v)<=3?EC[v][2]:e_xvec[3*((v)-4)+2])

/* ── undented check ──────────────────────────────────────────────────────── */
static double signed_sph_angle(double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz){
    double acx=ay*cz-az*cy,acy=az*cx-ax*cz,acz=ax*cy-ay*cx;
    double num=bx*acx+by*acy+bz*acz;
    double den=(ax*bx+ay*by+az*bz)*(bx*cx+by*cy+bz*cz)-(ax*cx+ay*cy+az*cz);
    return atan2(num,den);
}
static int undented_check(void){
    int ring[MAXRING]; double min_t=1e30;
    for(int v=1;v<=NV;v++){
        int k=cyclic_nbrs(v,ring);if(k<3)continue;
        double dx[MAXRING],dy[MAXRING],dz[MAXRING];
        for(int i=0;i<k;i++){int w=ring[i];
            double ex=EX(w)-EX(v),ey=EY(w)-EY(v),ez=EZ(w)-EZ(v);
            double L=sqrt(ex*ex+ey*ey+ez*ez);if(L<1e-30)L=1e-30;
            dx[i]=ex/L;dy[i]=ey/L;dz[i]=ez/L;}
        double tv=0;
        for(int i=0;i<k;i++)
            tv+=signed_sph_angle(dx[(i-1+k)%k],dy[(i-1+k)%k],dz[(i-1+k)%k],
                                  dx[i],dy[i],dz[i],dx[(i+1)%k],dy[(i+1)%k],dz[(i+1)%k]);
        if(tv<min_t)min_t=tv;
    }
    return (min_t>=0.0)?1:0;
}

/* ── Newton solver ───────────────────────────────────────────────────────── */
static double NJmat[MAXFN][MAXFN],NFvec[MAXFN];
static int lu_solve_n(int n){
    for(int col=0;col<n;col++){int piv=col;double best=fabs(NJmat[col][col]);
        for(int row=col+1;row<n;row++)if(fabs(NJmat[row][col])>best){best=fabs(NJmat[row][col]);piv=row;}
        if(best<1e-14)return -1;if(piv!=col){for(int k=col;k<n;k++){double t=NJmat[col][k];NJmat[col][k]=NJmat[piv][k];NJmat[piv][k]=t;}{double t=NFvec[col];NFvec[col]=NFvec[piv];NFvec[piv]=t;}}
        double inv=1.0/NJmat[col][col];for(int row=col+1;row<n;row++){double fac=NJmat[row][col]*inv;for(int k=col;k<n;k++)NJmat[row][k]-=fac*NJmat[col][k];NFvec[row]-=fac*NFvec[col];}}
    for(int i=n-1;i>=0;i--){double s=NFvec[i];for(int j=i+1;j<n;j++)s-=NJmat[i][j]*NFvec[j];NFvec[i]=s/NJmat[i][i];}return 0;}

static int euclid_newton(int *iters_out, double *res_out){
    int n=3*(NV-3); double res=1e30; int iter; int lufail=0;
    for(iter=0;iter<50;iter++){
        res=0;
        for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            NFvec[k]=dx*dx+dy*dy+dz*dz-1.0;
            double af=fabs(NFvec[k]);if(af>res)res=af;}
        if(res<1e-10) break;
        memset(NJmat,0,sizeof(double)*(size_t)n*MAXFN);
        for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            if(i>=4){int ci=3*(i-4);NJmat[k][ci]+=2*dx;NJmat[k][ci+1]+=2*dy;NJmat[k][ci+2]+=2*dz;}
            if(j>=4){int cj=3*(j-4);NJmat[k][cj]-=2*dx;NJmat[k][cj+1]-=2*dy;NJmat[k][cj+2]-=2*dz;}}
        for(int k=0;k<n;k++)NFvec[k]=-NFvec[k];
        if(lu_solve_n(n)<0){lufail=1;break;}
        /* backtrack: reject steps that worsen residual or dent */
        double step=1.0;
        double save[MAXFN]; memcpy(save,e_xvec,sizeof(double)*n);
        int accepted=0;
        for(int bt=0;bt<60;bt++,step*=0.5){
            for(int k=0;k<n;k++)e_xvec[k]=save[k]+step*NFvec[k];
            double res2=0;
            for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
                double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
                double af=fabs(dx*dx+dy*dy+dz*dz-1.0);if(af>res2)res2=af;}
            if(res2>=res)continue;
            if(!undented_check())continue;
            accepted=1;break;
        }
        if(!accepted)memcpy(e_xvec,save,sizeof(double)*n);
    }
    *iters_out=iter; *res_out=res;
    return (!lufail && res<1e-8)?1:0;
}

/* ── main ────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv){
    int do_polish=0;
    int argi=1;
    if(argi<argc && !strcmp(argv[argi],"-polish")){do_polish=1;argi++;}
    if(argc-argi<2){fprintf(stderr,"usage: euclidsolve [-polish] indir outdir < names.txt\n");return 1;}
    char *indir=argv[argi], *outdir=argv[argi+1];

    static char line[MAXCODE];
    long nets=0, ok_count=0, fail_count=0, dented=0;

    while(fgets(line,sizeof(line),stdin)){
        int ll=strlen(line);
        while(ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r'))line[--ll]='\0';
        if(!ll)continue;
        if(!decode(line)){fprintf(stderr,"decode failed: %s\n",line);continue;}
        build(); collect_edges();

        /* read last frame from .uhs file */
        char path[4096];
        snprintf(path,sizeof(path),"%s/%s.uhs",indir,line);
        FILE *fp=fopen(path,"r");
        if(!fp){
            /* no .uhs — propagate failure */
            snprintf(path,sizeof(path),"%s/%s.failed",outdir,line);
            FILE *ff=fopen(path,"w");if(ff)fclose(ff);
            fail_count++;nets++;build_clear();continue;
        }
        /* read last line */
        char lastline[65536]; lastline[0]='\0';
        char buf[65536];
        while(fgets(buf,sizeof(buf),fp)) strcpy(lastline,buf);
        fclose(fp);

        /* parse: rho x4 y4 t4 x5 y5 t5 ... */
        double rho_f;
        int pos=0; int nr;
        sscanf(lastline,"%lf%n",&rho_f,&nr); pos+=nr;
        double inv2rho=1.0/(2.0*rho_f);
        for(int v=4;v<=NV;v++){
            double x,y,t;
            sscanf(lastline+pos," %lf %lf %lf%n",&x,&y,&t,&nr); pos+=nr;
            e_xvec[3*(v-4)  ]=x*inv2rho;
            e_xvec[3*(v-4)+1]=y*inv2rho;
            e_xvec[3*(v-4)+2]=(t>0?log(t)*inv2rho:0);
        }

        /* check undented before polish */
        int und=undented_check();
        if(!und) dented++;

        /* compute unpolished residual */
        double res_before=0;
        for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            double af=fabs(dx*dx+dy*dy+dz*dz-1.0);if(af>res_before)res_before=af;}

        /* polish if requested */
        double res_after=res_before;
        int ok=1;
        if(do_polish && und){
            int iters;
            ok=euclid_newton(&iters,&res_after);
            und=undented_check();
        }

        /* write OBJ */
        if(und){
            snprintf(path,sizeof(path),"%s/%s.obj",outdir,line);
            fp=fopen(path,"w");
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
            ok_count++;
        } else {
            snprintf(path,sizeof(path),"%s/%s.failed",outdir,line);
            fp=fopen(path,"w");if(fp)fclose(fp);
            fail_count++;
        }

        build_clear();
        nets++;
    }
    fprintf(stderr,"euclidsolve: nets=%ld ok=%ld fail=%ld dented_before_polish=%ld polish=%s\n",
            nets,ok_count,fail_count,dented,do_polish?"yes":"no");
    return 0;
}
