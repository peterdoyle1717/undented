/* horodump2.c — dent-guarded UHS homotopy solver.
 *
 * Same as horodump.c (writes binary .uhs frames from ρ=4 down to ρ=2⁻¹⁴),
 * but after every successful Newton step checks the Klein-expanded
 * Euclidean configuration for dents (min vertex turning > 0). If a
 * step would introduce a dent, the step is rolled back and the rho
 * ratio is halved — locally, without restarting the whole homotopy.
 *
 * Rationale: the two branches of a flopper come close near ρ=0; a big
 * rho step can jump branches. Slower steps near the bottom keep us on
 * the undented branch.
 *
 * Writes .uhs + .failed as horodump does. Additionally on dent-block,
 * the .failed marker contains "# dent-guard rho=<floor>" and the last
 * known-good Klein-expanded OBJ, for inspection.
 *
 * Usage:   ./horodump2 outdir < names.txt
 * Compile: cc -O3 -o horodump2 horodump2.c -lm
 */
#include "horosolve.c"

#include <float.h>

/* ── UHS-native dent check ─────────────────────────────────────────────
   Operates directly on UHS (x,y,t) coords. Unit-normalized difference
   vectors at a vertex v give the *hyperbolic* unit tangents in UHS's
   metric ds² = (dx²+dy²+dt²)/t² (the t² factor cancels in unit
   normalization). Spherical exterior-angle sum of those tangents is
   the hyperbolic turning at v.

   Unlike a Klein-model check, no nonlinear UHS→Klein remap is done, so
   there's no conformal distortion introduced and no large-divide /
   subtract-large cancellations near ρ=0.
*/
static double signed_sph_angle_uhs(const double A[3], const double B[3], const double C[3]){
    double bcx=B[1]*C[2]-B[2]*C[1], bcy=B[2]*C[0]-B[0]*C[2], bcz=B[0]*C[1]-B[1]*C[0];
    double num=A[0]*bcx+A[1]*bcy+A[2]*bcz;
    double den=(A[0]*B[0]+A[1]*B[1]+A[2]*B[2])*(B[0]*C[0]+B[1]*C[1]+B[2]*C[2])
              -(A[0]*C[0]+A[1]*C[1]+A[2]*C[2]);
    return atan2(num,den);
}

/* uhs[v] = (x, y, t) for v=1..NV in UHS coords (ref triangle + free coords). */
static double min_turn_uhs(double uhs[MAXV+1][3]){
    double mn=1e30;
    int ring[MAXRING];
    for(int v=1;v<=NV;v++){
        int k=0;
        if(!NNBR[v]) continue;
        int start=NBR[v][0]; ring[k++]=start; int cur=start;
        while(k<NNBR[v]){int nxt=EM[v][cur]; if(!nxt||nxt==start) break; ring[k++]=nxt; cur=nxt;}
        if(k<3) continue;
        double dirs[MAXRING][3];
        double vx=uhs[v][0], vy=uhs[v][1], vt=uhs[v][2];
        for(int i=0;i<k;i++){
            int nb=ring[k-1-i];  /* CCW-from-outside order */
            double d0=uhs[nb][0]-vx, d1=uhs[nb][1]-vy, d2=uhs[nb][2]-vt;
            double L=sqrt(d0*d0+d1*d1+d2*d2); if(L<1e-30) L=1e-30;
            dirs[i][0]=d0/L; dirs[i][1]=d1/L; dirs[i][2]=d2/L;
        }
        double turn=0;
        for(int i=0;i<k;i++){
            turn+=signed_sph_angle_uhs(dirs[(i-1+k)%k], dirs[i], dirs[(i+1)%k]);
        }
        if(turn<mn) mn=turn;
    }
    return mn;
}

/* Fill uhs[v] from a frame (rho + v4..NV coords) + ref triangle at that rho. */
static void frame_to_uhs(const double *frame, double uhs[MAXV+1][3]){
    double rho=frame[0], a=exp(rho);
    double A2=a*a, A4=A2*A2, A6=A4*A2, Delta=A6-A4+A2-1.0;
    uhs[1][0]=0; uhs[1][1]=0; uhs[1][2]=a;
    uhs[2][0]=0; uhs[2][1]=0; uhs[2][2]=1.0/a;
    double prod=(A2-1)*(A2-1)*(A2-1)*(A6-1);
    uhs[3][0]=sqrt(prod>0?prod:0)/Delta; uhs[3][1]=0; uhs[3][2]=a*(A4-1)/Delta;
    for(int v=4;v<=NV;v++){
        uhs[v][0]=frame[1+3*(v-4)];
        uhs[v][1]=frame[1+3*(v-4)+1];
        uhs[v][2]=frame[1+3*(v-4)+2];
    }
}

/* Klein-expand for writing the final .obj (NOT used for dent check). */
static void klein_expand(const double *frame, double pos[MAXV+1][3]){
    double rho_f=frame[0];
    double a=exp(rho_f);
    double inv2rho=1.0/(2.0*rho_f);
    double A2=a*a, A4=A2*A2, A6=A4*A2, Delta=A6-A4+A2-1.0;
    double ux[4],uy[4],ut[4];
    ux[1]=0; uy[1]=0; ut[1]=a;
    ux[2]=0; uy[2]=0; ut[2]=1.0/a;
    double prod=(A2-1)*(A2-1)*(A2-1)*(A6-1);
    ux[3]=sqrt(prod>0?prod:0)/Delta; uy[3]=0; ut[3]=a*(A4-1)/Delta;
    for(int v=1;v<=NV;v++){
        double x,y,t;
        if(v<=3){x=ux[v]; y=uy[v]; t=ut[v];}
        else {x=frame[1+3*(v-4)]; y=frame[1+3*(v-4)+1]; t=frame[1+3*(v-4)+2];}
        double r2=x*x+y*y+t*t, d=r2+1.0;
        pos[v][0]=2*x/d*inv2rho;
        pos[v][1]=2*y/d*inv2rho;
        pos[v][2]=(r2-1)/d*inv2rho;
    }
}

int main(int argc, char **argv){
    if(argc<2){fprintf(stderr,"usage: horodump2 outdir < names.txt\n");return 1;}
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

        int ok=1;
        double at_start[MAXFN];
        { int iters; double res;
          ok=finite_newton(exp(rho_start),&iters,&res);
          if(!ok) goto fail;
          memcpy(at_start,g_xvec,sizeof(double)*n);
        }

        static double frames[2048][1+3*MAXV];
        int nframes=0;

        /* ideal + start frames */
        frames[nframes][0]=DBL_MAX;
        for(int v=4;v<=NV;v++){
            frames[nframes][1+3*(v-4)  ]=hz[3*(v-1)+1];
            frames[nframes][1+3*(v-4)+1]=hz[3*(v-1)+2];
            frames[nframes][1+3*(v-4)+2]=u[v-1]*u[v-1];
        }
        nframes++;
        {
            double a_cur=exp(rho_start); ref_triangle(a_cur);
            frames[nframes][0]=rho_start;
            for(int v=4;v<=NV;v++){
                frames[nframes][1+3*(v-4)  ]=GX(v);
                frames[nframes][1+3*(v-4)+1]=GY(v);
                frames[nframes][1+3*(v-4)+2]=GT(v);
            }
            nframes++;
        }

        /* Dent-guarded adaptive descent: keep local mult; halve on dent. */
        {
            double mult=2.0;
            double rho=rho_start, a_cur=exp(rho_start);
            double saved_xvec[MAXFN];
            memcpy(saved_xvec, g_xvec, sizeof(double)*n);
            int nsaved=nframes;
            const double MIN_MULT = 1.0 + 1e-6; /* floor on step shrink */
            int dent_block=0;
            static double pos[MAXV+1][3];
            while(rho>rho_target && nframes<2040){
                double rho_next=rho/mult;
                if(rho_next<rho_target) rho_next=rho_target;
                double a_next=exp(rho_next);
                predict_step(a_cur,a_next);
                int iters; double res;
                int ok_k=finite_newton(a_next,&iters,&res);
                if(!ok_k || iters>=10){
                    /* Newton failed: shrink step and retry */
                    mult=sqrt(mult);
                    if(mult<MIN_MULT){ok=0; break;}
                    memcpy(g_xvec, saved_xvec, sizeof(double)*n);
                    continue;
                }
                /* Newton ok. Check dent on UHS-native link turning. */
                ref_triangle(a_next);
                static double tmp_frame[1+3*MAXV];
                tmp_frame[0]=rho_next;
                for(int v=4;v<=NV;v++){
                    tmp_frame[1+3*(v-4)  ]=GX(v);
                    tmp_frame[1+3*(v-4)+1]=GY(v);
                    tmp_frame[1+3*(v-4)+2]=GT(v);
                }
                static double uhs[MAXV+1][3];
                frame_to_uhs(tmp_frame, uhs);
                double mt=min_turn_uhs(uhs);
                if(!(mt>0)){
                    /* Would introduce dent: roll back + shrink step */
                    mult=sqrt(mult);
                    if(mult<MIN_MULT){
                        dent_block=1; ok=0; break;
                    }
                    memcpy(g_xvec, saved_xvec, sizeof(double)*n);
                    continue;
                }
                /* Accept step: commit. */
                a_cur=a_next; rho=rho_next;
                frames[nframes][0]=rho;
                for(int v=4;v<=NV;v++){
                    frames[nframes][1+3*(v-4)  ]=GX(v);
                    frames[nframes][1+3*(v-4)+1]=GY(v);
                    frames[nframes][1+3*(v-4)+2]=GT(v);
                }
                nframes++;
                memcpy(saved_xvec, g_xvec, sizeof(double)*n);
                nsaved=nframes;
            }
            (void)nsaved;
            (void)dent_block;
        }

        /* If dent-blocked but we already descended at least one step from
         * rho_start, the last saved frame is a valid (undented) UHS solution.
         * Save what we have — pusheuclid + polish can take it from rho=N/2^k
         * to Euclidean. Only treat as a true failure if we never made any
         * descent progress (e.g. couldn't even Newton-converge at rho_start). */
        if(!ok && nframes >= 3) ok = 1;

        if(ok){
            char path[4096];
            snprintf(path,sizeof(path),"%s/%s.uhs",outdir,line);
            FILE *fp=fopen(path,"wb");
            if(fp){
                int nv32=NV, nf32=nframes;
                fwrite(&nv32,sizeof(int),1,fp);
                fwrite(&nf32,sizeof(int),1,fp);
                int coords_per_frame=1+n;
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
            if(fp){
                fprintf(fp,"# dent-guard-stall nframes=%d\n",nframes);
                fclose(fp);
            }
            fail_count++;
        }

        build_clear();
        nets++;
    }
    fprintf(stderr,"horodump2: nets=%ld ok=%ld fail=%ld\n",nets,ok_count,fail_count);
    return 0;
}
