/* horodump.c — UHS homotopy solver with frame dump (binary).
 *
 * Reads CLERS names from stdin, one per line.
 * For each net, writes outdir/NAME.uhs in binary:
 *   header: int32 NV, int32 nframes
 *   frame 0 (ideal): double rho=inf(DBL_MAX), double[3*(NV-3)] {x,y,u²}
 *   frame 1..n:      double rho,               double[3*(NV-3)] {x,y,t}
 * Vertices 1-3 omitted (recoverable from ref_triangle at each rho).
 * Writes outdir/NAME.failed on failure.
 *
 * Uses adaptive restart from horosolve.c.
 *
 * Usage:   ./horodump outdir < prime/N.txt
 * Compile: cc -O3 -o horodump horodump.c -lm
 */
#include "horosolve.c"

#include <float.h>

int main(int argc, char **argv){
    if(argc<2){fprintf(stderr,"usage: horodump outdir < names.txt\n");return 1;}
    char outdir[4096];
    strncpy(outdir,argv[1],sizeof(outdir)-1); outdir[sizeof(outdir)-1]='\0';

    static char line[MAXCODE];
    static double u[MAXV], hz[MAXV*3];
    long nets=0, ok_count=0, fail_count=0;

    while(fgets(line,sizeof(line),stdin)){
        int ll=strlen(line);
        while(ll>0&&(line[ll-1]=='\n'||line[ll-1]=='\r')) line[--ll]='\0';
        if(!ll) continue;
        if(!decode(line)){fprintf(stderr,"decode failed: %s\n",line);continue;}
        build(); collect_edges();
        horou(u); horoz(u,hz);

        int n=3*(NV-3);
        double rho_start=4.0;
        double rho_target=ldexp(1.0,-14);
        double a0=exp(rho_start);
        ref_triangle(a0);
        for(int v=4;v<=NV;v++){
            g_xvec[3*(v-4)  ]=hz[3*(v-1)+1];
            g_xvec[3*(v-4)+1]=hz[3*(v-1)+2];
            g_xvec[3*(v-4)+2]=(u[v-1]*u[v-1])/a0;
        }

        /* Converge at rho_start */
        int ok=1;
        double at_start[MAXFN];
        { int iters; double res;
          ok=finite_newton(exp(rho_start),&iters,&res);
          if(!ok) goto fail;
          memcpy(at_start,g_xvec,sizeof(double)*n);
        }

        /* Collect frames in memory */
        /* Max frames: ideal + start + up to 500 steps = ~502 */
        static double frames[502][1+3*MAXV]; /* rho + coords */
        int nframes=0;

        /* Ideal frame: rho=DBL_MAX, coords = hz positions + u² */
        frames[nframes][0]=DBL_MAX;
        for(int v=4;v<=NV;v++){
            frames[nframes][1+3*(v-4)  ]=hz[3*(v-1)+1];
            frames[nframes][1+3*(v-4)+1]=hz[3*(v-1)+2];
            frames[nframes][1+3*(v-4)+2]=u[v-1]*u[v-1];
        }
        nframes++;

        /* Starting frame */
        {
            double a_cur=exp(rho_start);
            ref_triangle(a_cur);
            frames[nframes][0]=rho_start;
            for(int v=4;v<=NV;v++){
                frames[nframes][1+3*(v-4)  ]=GX(v);
                frames[nframes][1+3*(v-4)+1]=GY(v);
                frames[nframes][1+3*(v-4)+2]=GT(v);
            }
            nframes++;
        }

        /* Adaptive restart loop */
        {
            double mult=2.0;
            int nrestart=0;
            for(;;){
                int failed_step=0;
                double rho=rho_start, a_cur=exp(rho_start);
                for(int step=0;step<500 && rho>rho_target;step++){
                    double rho_next=rho/mult;
                    if(rho_next<rho_target) rho_next=rho_target;
                    double a_next=exp(rho_next);
                    predict_step(a_cur,a_next);
                    int iters; double res;
                    int ok_k=finite_newton(a_next,&iters,&res);
                    if(!ok_k || iters>=10){failed_step=1; break;}
                    a_cur=a_next; rho=rho_next;

                    ref_triangle(a_cur);
                    frames[nframes][0]=rho;
                    for(int v=4;v<=NV;v++){
                        frames[nframes][1+3*(v-4)  ]=GX(v);
                        frames[nframes][1+3*(v-4)+1]=GY(v);
                        frames[nframes][1+3*(v-4)+2]=GT(v);
                    }
                    nframes++;
                }
                if(!failed_step) break;
                mult=sqrt(mult);
                nrestart++;
                if(nrestart>4){ok=0; break;}
                memcpy(g_xvec,at_start,sizeof(double)*n);
                /* Reset to just ideal + start */
                nframes=2;
            }
        }

        if(ok){
            char path[4096];
            snprintf(path,sizeof(path),"%s/%s.uhs",outdir,line);
            FILE *fp=fopen(path,"wb");
            if(fp){
                int nv32=NV, nf32=nframes;
                fwrite(&nv32,sizeof(int),1,fp);
                fwrite(&nf32,sizeof(int),1,fp);
                int coords_per_frame=1+n; /* rho + 3*(NV-3) doubles */
                for(int f=0;f<nframes;f++)
                    fwrite(frames[f],sizeof(double),coords_per_frame,fp);
                fclose(fp);
            }
            ok_count++;
        }

        fail:
        if(!ok){
            char path[4096];
            snprintf(path,sizeof(path),"%s/%s.failed",outdir,line);
            FILE *fp=fopen(path,"w");
            if(fp) fclose(fp);
            fail_count++;
        }

        build_clear();
        nets++;
    }
    fprintf(stderr,"horodump: nets=%ld ok=%ld fail=%ld\n",nets,ok_count,fail_count);
    return 0;
}
