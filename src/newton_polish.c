/* newton_polish.c — Newton-polish OBJs to machine-precision edge lengths.
 *
 * Reads OBJ from stdin, writes polished OBJ to stdout.
 * Gauge: pins v1, v2.y, v2.z, v3.z (6 DOFs removed).
 * Uses LAPACK dgels for least-squares solve.
 *
 * Compile: cc -O3 -o newton_polish newton_polish.c -lm -llapack
 * Usage:   ./newton_polish < in.obj > out.obj
 *          or: ls dir/*.obj | parallel './newton_polish < {} > outdir/{/}'
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV 200
#define MAXE 600
#define MAXITER 30

static double verts[MAXV][3];
static int faces[MAXE][3], nv, nf;
static int edges[MAXE][2], ne;

extern void dgels_(char*,int*,int*,int*,double*,int*,double*,int*,double*,int*,int*);

static void read_obj(void) {
    char line[1024]; nv=nf=0;
    while(fgets(line,sizeof(line),stdin)) {
        if(line[0]=='v'&&line[1]==' ')
            sscanf(line+2,"%lf%lf%lf",&verts[nv][0],&verts[nv][1],&verts[nv][2]), nv++;
        else if(line[0]=='f')
            sscanf(line+2,"%d%d%d",&faces[nf][0],&faces[nf][1],&faces[nf][2]), nf++;
    }
}

static void build_edges(void) {
    ne=0;
    for(int f=0;f<nf;f++) for(int k=0;k<3;k++) {
        int a=faces[f][k]-1, b=faces[f][(k+1)%3]-1;
        if(a>b){int t=a;a=b;b=t;}
        int found=0;
        for(int e=0;e<ne;e++) if(edges[e][0]==a&&edges[e][1]==b){found=1;break;}
        if(!found){edges[ne][0]=a;edges[ne][1]=b;ne++;}
    }
}

static void polish(void) {
    /* Free DOFs: all except v0(0,1,2), v1.y(4), v1.z(5), v2.z(8) */
    int pinned[] = {0,1,2,4,5,8};
    int npin=6, nfree=3*nv-npin;
    int *free_idx = malloc(nfree*sizeof(int));
    int fi=0;
    for(int i=0;i<3*nv;i++) {
        int p=0; for(int j=0;j<npin;j++) if(pinned[j]==i) p=1;
        if(!p) free_idx[fi++]=i;
    }

    double *J = malloc(ne*nfree*sizeof(double));
    double *r = malloc((ne>nfree?ne:nfree)*sizeof(double));
    double *work = malloc(10*ne*sizeof(double));
    int lwork=10*ne;

    for(int iter=0;iter<MAXITER;iter++) {
        double rho2=0;
        for(int e=0;e<ne;e++) {
            int i=edges[e][0], j=edges[e][1];
            double dx=verts[i][0]-verts[j][0];
            double dy=verts[i][1]-verts[j][1];
            double dz=verts[i][2]-verts[j][2];
            r[e] = dx*dx+dy*dy+dz*dz - 1.0;
            rho2 += r[e]*r[e];
            for(int k=0;k<nfree;k++) {
                int idx=free_idx[k];
                double deriv=0;
                if(idx/3==i) deriv = 2*(verts[i][idx%3]-verts[j][idx%3]);
                else if(idx/3==j) deriv = -2*(verts[i][idx%3]-verts[j][idx%3]);
                J[k*ne+e] = deriv; /* column-major for LAPACK */
            }
        }
        if(rho2 < 1e-30) break;

        /* Solve J^T dx = -r via dgels (least squares, column-major) */
        for(int e=0;e<ne;e++) r[e] = -r[e];
        char trans='N'; int m=ne, n=nfree, nrhs=1, lda=ne, ldb=ne>nfree?ne:nfree, info;
        dgels_(&trans,&m,&n,&nrhs,J,&lda,r,&ldb,work,&lwork,&info);

        double *x = (double*)verts;
        for(int k=0;k<nfree;k++) x[free_idx[k]] += r[k];
    }
    free(J); free(r); free(work); free(free_idx);
}

int main(void) {
    read_obj();
    build_edges();
    polish();
    for(int i=0;i<nv;i++)
        printf("v %.17g %.17g %.17g\n",verts[i][0],verts[i][1],verts[i][2]);
    for(int i=0;i<nf;i++)
        printf("f %d %d %d\n",faces[i][0],faces[i][1],faces[i][2]);
    return 0;
}
