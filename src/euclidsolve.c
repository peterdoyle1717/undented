/* euclidsolve.c — UHS last frame → Euclidean OBJ
 *
 * Reads binary .uhs files (from horodump), takes the last frame,
 * maps UHS→Klein, scales by 1/(2ρ) for Euclidean initial guess,
 * optionally polishes with Newton (residual backtracking, no dent gating).
 * Always writes .obj — the prover decides if it passes.
 *
 * Gauge: v1=(0,0,½), v2=(0,0,-½), v3=(√3/2,0,0)
 *   applied by pinning v1,v2,v3 to gauge positions.
 *
 * Usage:   ./euclidsolve [-polish] indir outdir < names.txt
 * Compile: cc -O3 -o euclidsolve euclidsolve.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

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

/* ── Euclidean coords ────────────────────────────────────────────────────── */
static double e_xvec[MAXFN];
static const double EC[4][3]={{0,0,0},{0,0,0.5},{0,0,-0.5},{0.8660254037844387,0,0}};
#define EX(v) ((v)<=3?EC[v][0]:e_xvec[3*((v)-4)])
#define EY(v) ((v)<=3?EC[v][1]:e_xvec[3*((v)-4)+1])
#define EZ(v) ((v)<=3?EC[v][2]:e_xvec[3*((v)-4)+2])

/* ── Newton solver (residual backtracking only) ──────────────────────────── */
static double NJmat[MAXFN][MAXFN],NFvec[MAXFN];
static int lu_solve_n(int n){
    for(int col=0;col<n;col++){int piv=col;double best=fabs(NJmat[col][col]);
        for(int row=col+1;row<n;row++)if(fabs(NJmat[row][col])>best){best=fabs(NJmat[row][col]);piv=row;}
        if(best<1e-14)return -1;if(piv!=col){for(int k=col;k<n;k++){double t=NJmat[col][k];NJmat[col][k]=NJmat[piv][k];NJmat[piv][k]=t;}{double t=NFvec[col];NFvec[col]=NFvec[piv];NFvec[piv]=t;}}
        double inv=1.0/NJmat[col][col];for(int row=col+1;row<n;row++){double fac=NJmat[row][col]*inv;for(int k=col;k<n;k++)NJmat[row][k]-=fac*NJmat[col][k];NFvec[row]-=fac*NFvec[col];}}
    for(int i=n-1;i>=0;i--){double s=NFvec[i];for(int j=i+1;j<n;j++)s-=NJmat[i][j]*NFvec[j];NFvec[i]=s/NJmat[i][i];}return 0;}

static void euclid_polish(void){
    int n=3*(NV-3);
    for(int iter=0;iter<50;iter++){
        double res=0;
        for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            NFvec[k]=dx*dx+dy*dy+dz*dz-1.0;
            double af=fabs(NFvec[k]);if(af>res)res=af;}
        if(res<1e-14) break;
        memset(NJmat,0,sizeof(double)*(size_t)n*MAXFN);
        for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            if(i>=4){int ci=3*(i-4);NJmat[k][ci]+=2*dx;NJmat[k][ci+1]+=2*dy;NJmat[k][ci+2]+=2*dz;}
            if(j>=4){int cj=3*(j-4);NJmat[k][cj]-=2*dx;NJmat[k][cj+1]-=2*dy;NJmat[k][cj+2]-=2*dz;}}
        for(int k=0;k<n;k++)NFvec[k]=-NFvec[k];
        if(lu_solve_n(n)<0)break;
        double step=1.0;
        for(int bt=0;bt<40;bt++,step*=0.5){
            double res2=0;
            for(int k=0;k<n_edges;k++){int i=eu[k],j=ev[k];
                double dx=(i>=4?EX(i)+step*NFvec[3*(i-4)]:EX(i))-(j>=4?EX(j)+step*NFvec[3*(j-4)]:EX(j));
                double dy=(i>=4?EY(i)+step*NFvec[3*(i-4)+1]:EY(i))-(j>=4?EY(j)+step*NFvec[3*(j-4)+1]:EY(j));
                double dz=(i>=4?EZ(i)+step*NFvec[3*(i-4)+2]:EZ(i))-(j>=4?EZ(j)+step*NFvec[3*(j-4)+2]:EZ(j));
                double af=fabs(dx*dx+dy*dy+dz*dz-1.0);if(af>res2)res2=af;}
            if(res2<res)break;
        }
        for(int v=4;v<=NV;v++){
            e_xvec[3*(v-4)  ]+=step*NFvec[3*(v-4)  ];
            e_xvec[3*(v-4)+1]+=step*NFvec[3*(v-4)+1];
            e_xvec[3*(v-4)+2]+=step*NFvec[3*(v-4)+2];}
    }
}

/* ── main ────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv){
    int do_polish=0;
    int argi=1;
    if(argi<argc && !strcmp(argv[argi],"-polish")){do_polish=1;argi++;}
    if(argc-argi<2){fprintf(stderr,"usage: euclidsolve [-polish] indir outdir < names.txt\n");return 1;}
    char *indir=argv[argi], *outdir=argv[argi+1];

    static char line[MAXCODE];
    long nets=0, ok_count=0, fail_count=0;

    while(fgets(line,sizeof(line),stdin)){
        int ll=strlen(line);
        while(ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r'))line[--ll]='\0';
        if(!ll)continue;
        if(!decode(line)){fprintf(stderr,"decode failed: %s\n",line);continue;}
        build(); collect_edges();

        int n=3*(NV-3);

        /* read binary .uhs file */
        char path[4096];
        snprintf(path,sizeof(path),"%s/%s.uhs",indir,line);
        FILE *fp=fopen(path,"rb");
        if(!fp){
            snprintf(path,sizeof(path),"%s/%s.failed",outdir,line);
            FILE *ff=fopen(path,"w");if(ff)fclose(ff);
            fail_count++;nets++;build_clear();continue;
        }
        int nv_file,nframes;
        fread(&nv_file,sizeof(int),1,fp);
        fread(&nframes,sizeof(int),1,fp);
        int cpf=1+n; /* doubles per frame: rho + 3*(NV-3) */
        /* seek to last frame */
        fseek(fp,sizeof(int)*2+(long)(nframes-1)*cpf*sizeof(double),SEEK_SET);
        static double frame[1+3*MAXV];
        fread(frame,sizeof(double),cpf,fp);
        fclose(fp);

        double rho_f=frame[0];
        double a=exp(rho_f);
        double inv2rho=1.0/(2.0*rho_f);

        /* ref triangle in UHS */
        double A2=a*a,A4=A2*A2,A6=A4*A2,Delta=A6-A4+A2-1.0;
        double ux[4],uy[4],ut[4];
        ux[1]=0;uy[1]=0;ut[1]=a;
        ux[2]=0;uy[2]=0;ut[2]=1.0/a;
        double prod=(A2-1)*(A2-1)*(A2-1)*(A6-1);
        ux[3]=sqrt(prod>0?prod:0)/Delta; uy[3]=0; ut[3]=a*(A4-1)/Delta;

        /* UHS→Klein→scale for ALL vertices (including gauge v1,v2,v3) */
        static double pos[MAXV+1][3];
        for(int v=1;v<=NV;v++){
            double x,y,t;
            if(v<=3){x=ux[v];y=uy[v];t=ut[v];}
            else{x=frame[1+3*(v-4)];y=frame[1+3*(v-4)+1];t=frame[1+3*(v-4)+2];}
            double r2=x*x+y*y+t*t, d=r2+1.0;
            pos[v][0]=2*x/d*inv2rho;
            pos[v][1]=2*y/d*inv2rho;
            pos[v][2]=(r2-1)/d*inv2rho;
        }

        /* write unpolished OBJ */
        if(!do_polish){
            snprintf(path,sizeof(path),"%s/%s.obj",outdir,line);
            fp=fopen(path,"w");
            if(fp){
                for(int v=1;v<=NV;v++)
                    fprintf(fp,"v %.17g %.17g %.17g\n",pos[v][0],pos[v][1],pos[v][2]);
                for(int i=0;i<NF;i++)
                    fprintf(fp,"f %d %d %d\n",F[i].a,F[i].b,F[i].c);
                fclose(fp);
            }
        } else {
            /* polish: pin gauge, Newton on v4+ */
            for(int v=4;v<=NV;v++){
                e_xvec[3*(v-4)  ]=pos[v][0];
                e_xvec[3*(v-4)+1]=pos[v][1];
                e_xvec[3*(v-4)+2]=pos[v][2];
            }
            euclid_polish();
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
        }
        ok_count++;

        build_clear();
        nets++;
    }
    fprintf(stderr,"euclidsolve: nets=%ld ok=%ld fail=%ld polish=%s\n",
            nets,ok_count,fail_count,do_polish?"yes":"no");
    return 0;
}
