/* euclid.c — Euclidean solver: Klein OBJ → unit-edge Euclidean OBJ
 *
 * Reads Klein model OBJ from hyper solver, converts to Euclidean
 * initial guess by scaling by 1/(2·rho), then runs Newton iteration
 * to solve |e|² = 1 for all edges.
 *
 * Gauge: v1=(0,0,½), v2=(0,0,-½), v3=(√3/2,0,0).
 *
 * Usage:   ./euclid indir outdir [rho] < prime/N.txt
 *          For each CLERS name, reads indir/NAME.obj, writes outdir/NAME.obj
 * Compile: cc -O3 -o euclid euclid.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV    400
#define MAXF    (2*MAXV+4)
#define MAXCODE (4*MAXV+8)
#define MAXRING 12
#define MAXE    (3*MAXV)
#define MAXFN   (3*MAXV)

/* ── geometry ────────────────────────────────────────────────────────────── */
typedef struct { int a,b,c; } Face;
static int NV, NF;
static Face F[MAXF];
static int EM[MAXV+1][MAXV+1];
static int NBR[MAXV+1][MAXRING], NNBR[MAXV+1];
static int eu[MAXE], ev[MAXE], n_edges;

static void nbr_add(int u, int w){for(int i=0;i<NNBR[u];i++)if(NBR[u][i]==w)return;NBR[u][NNBR[u]++]=w;}
static int cyclic_nbrs(int v, int ring[]){
    if(!NNBR[v])return 0; int start=NBR[v][0]; ring[0]=start; int k=1,cur=start;
    for(;;){int nxt=EM[v][cur];if(nxt==start)break;ring[k++]=nxt;cur=nxt;} return k;
}

/* ── read Klein OBJ ──────────────────────────────────────────────────────── */
static double kx[MAXV+1], ky[MAXV+1], kz[MAXV+1];

static int read_obj(const char *path){
    FILE *fp=fopen(path,"r"); if(!fp) return 0;
    char line[1024]; NV=0; NF=0;
    memset(EM,0,sizeof(EM)); memset(NNBR,0,sizeof(NNBR));
    while(fgets(line,sizeof(line),fp)){
        if(line[0]=='v'&&line[1]==' '){
            NV++; sscanf(line+2,"%lf%lf%lf",&kx[NV],&ky[NV],&kz[NV]);
        } else if(line[0]=='f'){
            int a,b,c; sscanf(line+2,"%d%d%d",&a,&b,&c);
            F[NF++]=(Face){a,b,c};
            EM[a][b]=c; EM[b][c]=a; EM[c][a]=b;
            nbr_add(a,b);nbr_add(b,a);nbr_add(b,c);nbr_add(c,b);nbr_add(a,c);nbr_add(c,a);
        }
    }
    fclose(fp);
    /* collect edges (skip gauge triangle edges) */
    n_edges=0;
    for(int i=1;i<=NV;i++) for(int j=i+1;j<=NV;j++){
        if(!EM[i][j]) continue;
        if(i<=3&&j<=3) continue;
        eu[n_edges]=i; ev[n_edges]=j; n_edges++;
    }
    return NV;
}

/* ── Newton solver ───────────────────────────────────────────────────────── */
static double NJmat[MAXFN][MAXFN], NFvec[MAXFN];
static double e_xvec[MAXFN];
static const double EC[4][3]={{0,0,0},{0,0,0.5},{0,0,-0.5},{0.8660254037844387,0,0}};
#define EX(v) ((v)<=3?EC[v][0]:e_xvec[3*((v)-4)])
#define EY(v) ((v)<=3?EC[v][1]:e_xvec[3*((v)-4)+1])
#define EZ(v) ((v)<=3?EC[v][2]:e_xvec[3*((v)-4)+2])

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

/* ── undented check ──────────────────────────────────────────────────────── */
static double signed_sph_angle(double ax,double ay,double az,
                                double bx,double by,double bz,
                                double cx,double cy,double cz){
    double acx=ay*cz-az*cy,acy=az*cx-ax*cz,acz=ax*cy-ay*cx;
    double num=bx*acx+by*acy+bz*acz;
    double den=(ax*bx+ay*by+az*bz)*(bx*cx+by*cy+bz*cz)-(ax*cx+ay*cy+az*cz);
    return atan2(num,den);
}

static int undented_check(void){
    int ring[MAXRING]; double min_t=1e30;
    for(int v=1;v<=NV;v++){
        int k=cyclic_nbrs(v,ring); if(k<3) continue;
        double dx[MAXRING],dy[MAXRING],dz[MAXRING];
        for(int i=0;i<k;i++){
            int w=ring[i];
            double ex=EX(w)-EX(v),ey=EY(w)-EY(v),ez=EZ(w)-EZ(v);
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

/* ── Euclidean Newton with undented backtracking ─────────────────────────── */
static int euclid_newton(int *iters_out, double *res_out){
    int n=3*(NV-3); double res=1e30; int iter; int lufail=0;
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
        /* backtrack: halve step if residual worsens or dent introduced */
        double step=1.0;
        double save[MAXFN]; memcpy(save,e_xvec,sizeof(double)*n);
        int accepted=0;
        for(int bt=0;bt<60;bt++,step*=0.5){
            for(int k=0;k<n;k++) e_xvec[k]=save[k]+step*NFvec[k];
            double res2=0;
            for(int k=0;k<n_edges;k++){
                int i=eu[k],j=ev[k];
                double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
                double af=fabs(dx*dx+dy*dy+dz*dz-1.0); if(af>res2) res2=af;
            }
            if(res2>=res) continue;
            if(!undented_check()) continue;
            accepted=1; break;
        }
        if(!accepted) memcpy(e_xvec,save,sizeof(double)*n);
    }
    *iters_out=iter; *res_out=res;
    return (!lufail && res<1e-8)?1:0;
}

/* ── polish ──────────────────────────────────────────────────────────────── */
static void polish(void){
    int n=3*(NV-3);
    for(int pi=0;pi<20;pi++){
        double res=0;
        for(int k=0;k<n_edges;k++){
            int i=eu[k],j=ev[k];
            double dx=EX(i)-EX(j),dy=EY(i)-EY(j),dz=EZ(i)-EZ(j);
            NFvec[k]=dx*dx+dy*dy+dz*dz-1.0;
            double af=fabs(NFvec[k]); if(af>res) res=af;
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

/* ── main ────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv){
    if(argc<3){fprintf(stderr,"usage: euclid indir outdir [rho] < names.txt\n");return 1;}

    char indir[4096], outdir[4096];
    strncpy(indir,argv[1],sizeof(indir)-1); indir[sizeof(indir)-1]='\0';
    strncpy(outdir,argv[2],sizeof(outdir)-1); outdir[sizeof(outdir)-1]='\0';

    double rho = 0.01;
    if(argc>=4) rho = atof(argv[3]);
    double scale = 1.0/(2.0*rho);

    static char line[MAXCODE];
    long nets=0, ok_count=0, fail_count=0;

    while(fgets(line,sizeof(line),stdin)){
        int ll=strlen(line);
        while(ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r')) line[--ll]='\0';
        if(!ll) continue;

        /* read Klein OBJ */
        char path[4096];
        snprintf(path,sizeof(path),"%s/%s.obj",indir,line);
        if(!read_obj(path)){
            /* no Klein OBJ (hyper failed) — propagate failure */
            snprintf(path,sizeof(path),"%s/%s.failed",outdir,line);
            FILE *fp=fopen(path,"w"); if(fp) fclose(fp);
            fail_count++; nets++; continue;
        }

        /* initial guess: scale Klein coords, pin gauge */
        for(int v=4;v<=NV;v++){
            e_xvec[3*(v-4)  ]=kx[v]*scale;
            e_xvec[3*(v-4)+1]=ky[v]*scale;
            e_xvec[3*(v-4)+2]=kz[v]*scale;
        }

        /* Newton */
        int iters; double res;
        int ok=euclid_newton(&iters,&res);
        if(ok) polish();
        int und=undented_check();
        if(!ok||!und){
            fprintf(stderr,"FAIL %s iters=%d res=%.2e %s %s\n",line,iters,res,
                    ok?"CONVERGED":"FAILED",und?"UNDENTED":"DENTED");
        }

        /* write output */
        snprintf(path,sizeof(path),"%s/%s.%s",outdir,line,(ok&&und)?"obj":"failed");
        FILE *fp=fopen(path,"w");
        if(fp){
            if(ok&&und){
                fprintf(fp,"v %.17g %.17g %.17g\n",EC[1][0],EC[1][1],EC[1][2]);
                fprintf(fp,"v %.17g %.17g %.17g\n",EC[2][0],EC[2][1],EC[2][2]);
                fprintf(fp,"v %.17g %.17g %.17g\n",EC[3][0],EC[3][1],EC[3][2]);
                for(int v=4;v<=NV;v++)
                    fprintf(fp,"v %.17g %.17g %.17g\n",e_xvec[3*(v-4)],e_xvec[3*(v-4)+1],e_xvec[3*(v-4)+2]);
                for(int i=0;i<NF;i++)
                    fprintf(fp,"f %d %d %d\n",F[i].a,F[i].b,F[i].c);
            }
            fclose(fp);
        }
        if(ok&&und) ok_count++; else fail_count++;

        /* clear EM for next net */
        memset(EM,0,sizeof(EM)); memset(NNBR,0,sizeof(NNBR));
        nets++;
    }
    fprintf(stderr,"euclid: nets=%ld ok=%ld fail=%ld\n",nets,ok_count,fail_count);
    return 0;
}
