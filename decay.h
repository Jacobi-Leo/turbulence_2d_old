        include 'mpif.h'
        double precision rand
        integer my,mx,mmy,mmx2,mx2,nproc,nsteps,nsteps0
        common /data/ my,mx,mmy,mmx2,mx2,nproc,nsteps,nsteps0
        real time,timestar,dt,dt_h,rnu,alpha,ak0,e0
        common /cqd_dat1/ time,timestar,dt,dt_h,rnu,alpha,ak0,e0
        real pi,Re_lmda,lmda,turb_u
        common /cqd_dat2/ pi,Re_lmda,lmda,turb_u        
