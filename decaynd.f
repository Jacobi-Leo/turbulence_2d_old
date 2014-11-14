!*************************************
      program spectral_mpi

      include 'decay.h'

      complex,allocatable,dimension(:, :)::vx,vy,wz,owz
      complex,allocatable,dimension(:, :)::wt, uxt, wt2
      real,allocatable,dimension(:, :)::kx,ky,k2,k2e,tmp,tmp1
      real, allocatable,dimension(:)::wx,wy,ek,e_t
      integer,allocatable,dimension(:)::ipx,ipy         
      integer,allocatable,dimension(:)::iseed

      integer new,id,numb,resultlen
      integer nstep,itout,ieout,iseed_m0
      integer jj,jjj,i,j,k,istep,idp,idp1,iii
      real xk,xx,yy,cy,sx1,sx2,ekr,e_tr
      real umax,rmax1,rmax2,umax_t,xmax,ymax,dt_max
      real flatness,vor(2),vort(2),sum1
      real*8 num,br1,timer2 

      character iopath*50,name*60,fin*80
      common /iopathc/ iopath
      common /timeseq/ rnv,vormax,vorave,gamave,radave,areave,epsave,
     1                 vort,flatness

!  setup MPI environment

      call mpi_init(ierror)
      call mpi_comm_size(mpi_comm_world,nproc,ierror)
      call mpi_comm_rank(mpi_comm_world,id,ierror)
      call mpi_barrier(mpi_comm_world,ierror)
      nallgrp = mpi_comm_world
      call mpi_get_processor_name(name,resultlen,ierror)
      write(*,98) id,name,resultlen
98    format(1x,'myid=',i4,4x,'name=',a20,4x,'resultlen=',i2)

      if (id.eq.0)  then              
        open(1,file='decay.in',status='old') 
           write(*,*)'nproc=', nproc
           read(1,*) iopath
           read(1, *) my
           write(*, *) '  my=', my
           read(1,*) iseed_m0
           write(*,*) '   seed=', iseed_m0
           read(1,*) e0
           write(*,*) '   e0=', e0
           read(1,*) ak0
           write(*,*) '   ak0=', ak0
           read(1,*)  rnu
           write(*,*) '   rnu=', rnu
           read(1,*)  alpha
           write(*,*) '   alpha=', alpha
           read(1,*) dt
           write(*,*) '   dt=', dt
           read(1,*) nstep
           write(*,*) '   nstep=', nstep
           read(1,*) itout
           write(*,*) '   itout=', itout
           read(1,*) ieout
           write(*,*) '   ieout=', ieout
           read(1,*) new
           write(*,*) '   new=', new
           read(1,*) idp
           write(*,*) '   idp(for rerun)=', idp
          close(1)
      endif

! Broadcast inputs across processors, since just read into id=0.
!      call mpi_bcast(iopath,50,MPI_CHARACTER,0,nallgrp,ierror)
      call mpi_bcast(my,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(iseed_m0,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(ak0,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(e0,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(nstep,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(itout,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(ieout,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(dt,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(rnu,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(alpha,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(new,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(idp,1,MPI_INTEGER,0,nallgrp,ierror)
!-------------------------------
      mx2=my/2
      mx=my
      mmx2=mx2/nproc
      mmy=my/nproc
      pi = 4.0 * atan(1.0)
      nek = int(sqrt(2.0)*my/3.0)  
      scale = 1.0/my/mx        
!---------------------------------
      jj = 1
      jjj = 1
	  
      iii = 1
      idp1 = idp        
!----------------------------------
!... allocate memory...................................................

        allocate (vx(my,mmx2) )
        allocate (vy(my,mmx2) )
        allocate (wz(my,mmx2) )
        allocate (owz(my,mmx2) )
        allocate (kx(my,mmx2) )
        allocate (ky(my,mmx2) )
        allocate (k2(my,mmx2) )
        allocate (k2e(my,mmx2) )
        allocate (tmp(my,mmx2) )
        allocate (tmp1(my,mmx2) )
        allocate (wt(mx2,mmy) )
		allocate (wt2(mx2/2, mx2))
        allocate (uxt(mx2,mmy) )       
        allocate (wx(0:(mx2-1)) )
        allocate (wy(0:(my-1)) )
        allocate (ek(nek) )
        allocate (e_t(nek) )
        allocate (ipx(0:mx2) )  
        allocate (ipy(0:my) )
        allocate (iseed(nproc) )

!***********************************
!===output info.txt=========
       if (id.eq.0) then

         write(fin,109) iopath,'info.txt'
         call movespa(fin,80)
         open(70,file=fin,access='append',status='unknown')
             write(70,701) nproc,my,iseed_m0,e0,ak0,rnu,alpha,dt,nstep,
     1                     itout,ieout,new,idp,nek,iopath
         close(70)
701      format(1x,'nproc=',i3,/1x,'my=',i4,4x,'iseed_m0=',i10,
     1         /1x,'e0=',f10.2,4x,'ak0=',f10.2,
     1         /1x,'rnu=',e15.6,4x,'alpha=',e15.6,4x,'dt=',e15.6,
     1         /1x,'nstep=',i10,4x,'itout=',i5,4x,'ieout=',i5,
     1         /1x,'new=',i2,4x,'idp=',i3,5x,'nek=',i4,
     1         /1x,'iopath=',a60)
       endif
109    format(1x,a50,a)

!---Prepare for call mpifft----------

         ipx(0) = 0  
         ipy(0) = 0

!... If the computation is from the initial valve, then at the beginning,
!    we choose a small time step.
        if(new.ne.0) then
          dt = 0.1*dt
          nsteps0 = 0
          time = 0.0
        endif
        dt_h = 0.5*dt

        call wavenumber(kx,ky,k2,k2e,id,0)

!...INITIAL CONDITIONS
!  new = 1:  random vorticity field with specified spectrum 
!             vor(k)=u0*k**6/((1+k/ak0)**18)
!  new = 2:  uniform disstributed vortex:
!             psi(x,y) = sin(ak0*x)cos(ak0*y),               stream function
!             vor(x,y) = -ak0*ak0*sin(ak0*x)cos(ak0*y)

      if(new.ne.0) then
      if (new.eq.1) then
        if (id.eq.0) then
          numb = irand (iseed_m0)
          do i = 1,nproc
            iseed(i) = irand(0)
          enddo
        endif
        call mpi_bcast(iseed,nproc,MPI_INTEGER,0,nallgrp,ierror)
        i = iseed(id+1)
        !num =  drand(i)
		call RANDOM_NUMBER(num)

        call gaussian(vx)
        call gaussian(vy)

!... projection:
         tmp = (kx*real(vx) + ky*real(vy))/k2
         vx = cmplx(real(vx) - kx*tmp, aimag(vx))
         vy = cmplx(real(vy) - ky*tmp, aimag(vy))
         tmp = (kx*aimag(vx) + ky*aimag(vy))/k2
         vx = cmplx(real(vx), aimag(vx) - kx*tmp)
         vy = cmplx(real(vy), aimag(vy) - ky*tmp)

         if(id.eq.0) then
           do i=1,nek
             xk = i
             ek(i) = xk**6/(1.0+xk/ak0)**18
           enddo
           sum1 = sum(ek)
           c0 = e0/sum1
         endif
         call mpi_bcast(c0,1,MPI_REAL,0,nallgrp,ierror) 
         c0 = sqrt(c0)
         tmp1 = sqrt(k2)+0.5
         tmp = c0*tmp1**3/(1.0+tmp1/ak0)**9
         vx = vx*tmp
         vy = vy*tmp
         call symmetrize(vx,id)
         call symmetrize(vy,id)
!... vorticity:
         wz = (0.,1.) * (kx*vy - ky*vx)
         call symmetrize(wz,id)

       elseif(new.eq.2) then
!... uniform distributed vortex: Tayor-Green vortes in two-dimensional:
         xk = 2.*pi/float(my)
         do j=1,mmy
            yy = (id*mmy+j-1)*xk
            cy = cos(ak0*yy)
         do i=1,mx2
            xx = 2*(i-1)*xk
            sx1 = sin(ak0*xx)
            xx = xx+xk
            sx2 = sin(ak0*xx)
            if(sx1*cy.ge.0.1) then
               sx1 = 1.0
            elseif(sx1*cy.le.-0.1) then
               sx1 = -1.0
            endif
            if(sx2*cy.ge.0.1) then
               sx2 = 1.0
            elseif(sx2*cy.le.-0.1) then
               sx2 = -1.0
            endif
            uxt(i,j)= e0*cmplx(sx1,sx2)
         enddo
         enddo
!...TRANSFORM BACK TO K-SPACE
         isign = 1       
         call newfft (wz,uxt, isign,ipx,ipy,wx,wy,id,nallgrp)  
         call symmetrize (wz, id)
         call dealiasing (wz, k2) 
       endif

!--calculating initial spectrum------
        tmp= wz*conjg(wz) / k2
        if (id.eq.0)  tmp(:, 1)=0.50*tmp(:, 1)  
        do i=1,nek
            ek(i)=sum(tmp,mask=(abs(sqrt(k2)-i-0.499999).lt.0.5) )
        enddo
        call mpi_reduce(ek,e_t,nek,MPI_REAL,MPI_SUM,
     +                        0,mpi_comm_world,ierror)
        if (id.eq.0) then
            write(fin,109) iopath,'initial.sp'
            call movespa(fin,80)
            open(70,file=fin,access='append',status='unknown')
            do i=1,nek
              write(70,1002) i,  e_t(i)    
            enddo
            close(70)
1002        format(1x,i6,2x,e15.6)
        endif

      else
         call input (wz,idp,id,nallgrp)
         if (id.eq.0) wz(1,1)=(0.0, 0.0)     
      endif

      call symmetrize (wz, id)
      if(nsteps0.ge.ieout) call pickvor0_input(id)

!...*********** MAIN LOOP *************

      timer1 = mpi_wtime()
      do istep = 0, nstep

        nsteps = nsteps0 + istep

!....If the computation is started from initial value, then the first
!    10 steps is used small time step
        if(nsteps.eq.10) then
          dt = 2.0*dt
          dt_h = 0.5*dt
          call wavenumber(kx,ky,k2,k2e,id,0)
        endif
        if(nsteps.eq.100) then
          dt = 5.0*dt
          dt_h = 0.5*dt
          call wavenumber(kx,ky,k2,k2e,id,0)
        endif

        if (mod(nsteps,100).eq.0.or.nsteps.le.100) then
            timer2 = mpi_wtime()
            if(id.eq.0) write(*,19) nsteps,time,timer2-timer1
            timer1 = timer2
        endif
19      format(/2x,'nsteps =',i6,3x,'Time =',F8.4,
     1           3x,'WallTime =',F12.0,' msec')

        call wavenumber(kx,ky,k2,k2e,id,1)

!!**********************************************************
!---Some extra job need to do here to produce the correct spectrum
!---to do 1----
!...WRITE ENERGY SPECTRUM
         if (mod(nsteps,ieout).eq.0) then
            tmp= wz*conjg(wz) / k2
            if (id.eq.0)  tmp(:, 1)=0.50*tmp(:, 1)  
            do i=1,nek           
                 ek(i)=sum(tmp,mask=(abs(sqrt(k2)-i-0.499999).lt.0.5) )
            enddo
            call mpi_reduce(ek,e_t,nek,MPI_REAL,MPI_SUM,0,
     +                                       nallgrp,ierror)
            if (id.eq.0) then
             write(fin,109) iopath,'spectrum.d'
             call movespa(fin,80)
             open(20,file=fin,access='append',status='unknown')
             write(20,201) jjj-1,nsteps,time,sum(e_t)
             do i=1,nek
                 write(20, 1002) i,  e_t(i)
             enddo
             close(20)
            endif    
             jjj = jjj + 1
         endif
201      format(1x,'#k',3x,'j=',i3,3x,'nsteps=',i6,3x,'Time=',f10.4,3x,
     +                'Te=',e15.6)
   
!...STORE VORTICITY TEMPORARILY
         tmp = real(wz)
         tmp1 = aimag(wz)

!...output K-space vorticity if it's time
         if ( (mod(nsteps,itout).eq.0).and.(istep.ne.0) ) then
             idp1=idp1+1
             call output (wz,idp1,id,nallgrp) 
         endif

!-------------NONLINEAR TERM----------

!...TRANSFER VELOCITY AND VORTICITY TO X-SPACE
         vx = (0.0,1.0) * ky * wz / k2
         call symmetrize (vx, id)
         isign= -1      
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
         call newfft (wz,wt, isign,ipx,ipy,wx,wy,id,nallgrp) 

		 !===== to do 3=============
		 !===DONE===
		 forall (i=1:mx2/2, j=1:mx2)
			 wt2(i,j) = wt(i,j)
		 end forall

         if(nsteps.eq.ieout) call pickvor0(wt2,mx2,id)
         if(mod(nsteps,ieout).eq.0) then
          call pickvor(wt2,mx2)
         endif

!...Calculate flatness of vorticity.................................
!..  vor(1): < vor**2 >         vor(2): < vor**4 >
      if (mod(nsteps,ieout).eq.0) then
            vor(1)=sum( real(wt*scale)**2 + aimag(wt*scale)**2 ) 
            vor(2)=sum( real(wt*scale)**4 + aimag(wt*scale)**4 )
            call mpi_reduce(vor,vort,2,MPI_REAL,MPI_SUM,0,
     +                                   nallgrp,ierror)	
        if (id.eq.0) then
            flatness = vort(2)/vort(1)/vort(1)
            write(*,202) flatness
        endif
      endif 
202   format(18x,'Flatness of voriticity =',e15.6)
!1.....< ux^2 >: only one componet is computed here .......................
        if (mod(nsteps,ieout).eq.0) then   
         ekr = 0.5*sum( real(uxt)**2 + aimag(uxt)**2 )
        endif
!...Find max x-velocity (for CFL condition )--------------------       
!         if ( mod(istep-1, ieout).eq.0 ) then
!              rmax1= maxval (real(uxt)) 
!              rmax2= maxval (aimag(uxt))
!              xmax= amax1(rmax1,rmax2) 
!         endif 
!...FORM THE PRODUCT V*W IN X-SPACE
         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )

!...TRANSFORM BACK TO K-SPACE
         isign = 1       
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp)  
         call symmetrize ( vx, id)
         call dealiasing (vx, k2) 

         vy = vx * (0.0,-0.50)*kx

!===============================================
!     now do on y-component of velocity
!     vort. in x-space is already saved in wt
!     wz also changed 

         vx=-(0.0,1.0)* kx*cmplx(tmp,tmp1)/k2

         isign = -1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp)

!...Find max y-velocity (for CFL condition)-------------------------------
!          if ( mod(istep-1, ieout).eq.0 ) then
!              rmax1= maxval (real(uxt)) 
!              rmax2= maxval (aimag(uxt))
!              ymax= amax1(rmax1,rmax2) 
!              umax= sqrt (xmax*xmax + ymax*ymax) 
!              call mpi_reduce(umax,umax_t,1,MPI_REAL,MPI_MAX,
!     +                                     0,nallgrp,ierror)
!              call mpi_bcast( umax_t, 1, MPI_REAL, 0, nallgrp, ierror) 
!          endif

!...Find kinetic energy in x-space, E = 0.5 * < u^2 > .............
!2.....< u^2 > the second component ..................
       if (mod(nsteps,ieout).eq.0) then          
         ekr = ekr + 0.5*sum( real(uxt)**2 + aimag(uxt)**2 )
         call mpi_reduce(ekr,e_tr,1,MPI_REAL,MPI_SUM,0,nallgrp,ier)
         if (id.eq.0) then
            write(fin,109) iopath,'vorstat.d'
            call movespa(fin,80)
            open(83,file=fin,access='append',status='unknown')
            write(83,89) time,rnv,vormax,vorave,gamave,radave,areave,
     1                   epsave,e_tr,vort(1),vort(2),flatness
            close(83)
         endif                                                  
       endif
89     format(1x,12e15.6)
!-------------------------------------------------------
         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )
         
         isign = 1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
         call symmetrize (vx, id)
         call dealiasing (vx, k2)
   
         vy = vy + (0.0,-0.50) * ky * vx

!-Recover vort.,do phase shift dealiasing on x-comp of velo.
         wz = cmplx (tmp, tmp1)
         vx = (0.0,1.0) * ky * wz / k2

!...PHASE SHIFT
         vx=vx*cmplx(cos(pi/my*(kx+ky)),sin(pi/my*(kx+ky)))
         wz=wz*cmplx(cos(pi/my*(kx+ky)),sin(pi/my*(kx+ky)))

         call symmetrize(vx,id)
         call symmetrize(wz,id)

!...TRANSFER VELOCITY AND VORTICITY TO X-SPACE
         isign = -1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
         call newfft (wz,wt, isign,ipx,ipy,wx,wy,id,nallgrp)  

!...FORM THE PRODUCT V*W IN X-SPACE
         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )
        
!...TRANSFORM TO K-SPACE
         isign = 1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
!...PHASE SHIFT
         vx=vx*cmplx(cos(pi/my*(kx+ky)),-sin(pi/my*(kx+ky)) )
         call symmetrize(vx, id)

         vy = vy + vx*(0.0,-0.50)*kx

!--- do phase shift dealiasing for y-comp
!--- phase-shifted wz in real space already saved in wt

         vx=-(0.0,1.0)* kx*cmplx(tmp, tmp1)/k2

         vx=vx*cmplx(cos(pi/my*(kx+ky)),sin(pi/my*(kx+ky)))

         call symmetrize(vx,id)

         isign = -1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 

         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )

         isign = 1
        call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
        vx=vx*cmplx(cos(pi/my*(kx+ky)),-sin(pi/my*(kx+ky)))

         call symmetrize(vx, id)

         vy = vy + vx*(0.0,-0.50)*ky

!-----------END CONVOLUTION PLUS DEALISING-----------
!-----now, vy = -i k.fft(Vw)
!...RECOVER VORTICITY
           wz = cmplx( tmp,tmp1)

!!---Commented out by LIU on 20141110 for simplicity of code
!!
!!...CALCULATE ENSTROPHY FLUX (it must be done HERE)
!!...enstrophy flux is
!!   2 Imag [ w^*(k) k.F(uw) ]  (F(.) means fourier transform)
!         if (mod(nsteps,ieout).eq.0) then
!            vx =  2.* (0.0, 1.0)*vy
!            kx = real(wz)*aimag(vx) - aimag(wz)*real(vx)
!            if (id.eq.0)  kx (:,1)=0.50*kx (:,1)  
!            sum1=0.0 
!            do i = 1,nek          
!               a = i - 1 + 0.5
!               ek(i) =sum (kx,mask=(sqrt(k2).ge.a.and.sqrt(k2).lt.a+1) )
!            enddo
!            call mpi_reduce(ek, e_t, nek, MPI_REAL,MPI_SUM,0,
!     +           nallgrp,ierror)
!            if (id.eq.0) then
!               write(fin,109) iopath,'enstflx.d'
!               call movespa(fin,80)
!               open(82,file=fin,access='append',status='unknown')
!               write(82, 821) iii,nsteps,time
!               do i=1,nek
!                   sum1=sum1+e_t(i)
!                   write(82, 1002) i, -sum1
!               end do
!               close(82)
!            endif
!            iii = iii + 1
!
!         endif
!821      format(1x,'#k',3x,'i=',i3,3x,'nsteps=',i6,3x,'time=',f10.4) 

!-------------END NONLINEAR TERM----------
         call dealiasing (vy, k2)
         call symmetrize (vy, id)

!... do CFL condition with umax, h=1./my
!             if ( (id.eq.0) .and. (mod(istep-1,100).eq.0) ) then
!                  dt_max=1./ umax_t /float(my)
!c                  if (dt .ge. (0.50*dt_max) ) then
!                     dt = 0.5 * dt_max         
!                     dt_h = 0.5 * dt
!                     call mpi_bcast(dt,1,MPI_REAL,0,nallgrp,ierror)
!                     call mpi_bcast(dt_h,1,MPI_REAL,0,nallgrp,ierror)
!                  endif
!             endif
!---If first-step then use modified Euler method
         if (istep.eq.0) then
             owz = wz
             wz = sqrt(k2e)*(wz+dt_h*vy)
             call output(vy,0,id,nallgrp)
             time=time+dt_h
         else if (istep.eq.1) then
             wz = k2e*(owz+dt*vy/sqrt(k2e))
             time=time+dt_h
             call input(owz,0,id,nallgrp)
!---Adams-banshford--------
         else
            wz = wz + dt_h * (3.0*vy - k2e * owz)
            wz = wz * k2e
            owz=vy
            time=time+dt
         end if

!---Add the update proccegure here.----
!---to do 2 ===DONE===
!
         isign = -1
         call newfft (wz, wt, isign,ipx,ipy,wx,wy,id,nallgrp)  
		 call update (wt, my)
		 isign = 1
         call newfft (wz, wt, isign,ipx,ipy,wx,wy,id,nallgrp)  
		 call symmetrize (wz, id)
		 !call dealiasing (wz, k2)
        enddo                     

      write(*,991) id
991   format(1x,'processor id=',i3,3x,'finished.')
      deallocate(vx,vy,wz,owz,kx,ky,k2,k2e,tmp,tmp1,wt,uxt,wx,wy,ek,e_t)
      deallocate(ipx,ipy,iseed)

      call MPI_FINALIZE(ierror)
      stop
      end

!--------------------------------------------

      subroutine dealiasing(vx,k2)
!...8/9 rule for dealiasing convolution in FFT

      include 'decay.h'
      complex, dimension(my, mmx2)::vx
      real,dimension (my, mmx2)::k2
      real ass

      ass = 2.0/9.0*(my*my)
      where(k2.ge.ass)
        vx =  (0.0, 0.0)
      endwhere
      return
      end
!----------------------------------------------
      subroutine wavenumber(kx,ky,k2,k2e,id,kc)

      include 'decay.h'
      real,dimension (my, mmx2)::kx,ky,k2,k2e
      integer id
      integer i, j, j1
      

      do j=1,mmx2
         j1 = id * mmx2 + j
         kx (:, j) = float( j1-1 )
      enddo
      if(kc.eq.1) return

      do i=1,my
         ky (i, :) = float(mod( i-1+my/2, my)-my/2)
      enddo
      if(kc.eq.2) return
    

      k2 = kx*kx + ky*ky
      k2e = exp(-(rnu*k2+alpha)*dt)
      if (id.eq.0) then
         k2(1,1) = 0.5         
      endif

      return
      end

!------------------------------------------------
      subroutine gaussian (u)
      include 'decay.h'
      complex, dimension(my, mmx2) :: u
      real t1,t2

      u  = (0.0, 0.0)
      do i = 1,my
         do j = 1,mmx2
            call RANDOM_NUMBER(t1)
            call RANDOM_NUMBER(t2)
            !t1 = drand(0)
            !t2 = drand(0)
            if (t1.le.1.e-10) t1 = 1.e-10
            if (t2.le.1.e-10) t2 = 1.e-10
            t2 = 2.0*pi*t2
            u (i, j) = sqrt(-2.0*log(t1))*cmplx(cos(t2),sin(t2)) 
         enddo
      enddo

      return
      end

!-----------------------------------------------------

      subroutine symmetrize(c,id)
      include 'decay.h'
      complex,dimension(my, mmx2)::c

        c(my/2+1,:) = (0.0, 0.0)

        if (id.eq.0) then
          c(1,1) = (0.0, 0.0)
          do iy = 2,my/2-1
            c(iy,1) = .50*(c(iy,1)+ conjg(c(mod(my+1-iy,my)+1,1)))  
          enddo

           do iy = 2, my/2-1
             iy2 = mod(my+1-iy,my) + 1
             c(iy2, 1) = conjg (c(iy, 1) )
            enddo
           
         endif

      return
      end
! ---------------------------------------------------
      subroutine outputp (uxt,idp,id,nallgrp)
      include 'decay.h'
      complex,dimension(mx2,mmy)::uxt
      character iopath*50,fin*80
      common /iopathc/ iopath

      write(fin,1) iopath,idp,id
1     format(1x,a50,'pvort',i3.3,'.',i3.3)
      call movespa(fin,80)
      open(10,file=fin,status='unknown')
         write(10,*) 2*mx2,mmy
         write(10,*) (( real(uxt(i,j)),aimag(uxt(i,j)),i=1,mx2),j=1,mmy)
         write(10,*) nsteps,time
      close(10)

      return
      end
! --------------------------------------------------
      subroutine output (ux,idp,id,nallgrp)
      include 'decay.h'
      complex,dimension(my, mmx2)::ux
      character iopath*50,fin*80
      common /iopathc/ iopath

      write(fin,1) iopath,idp,id
      write(*,*) 'idp,id=',idp,id
      write(*,*) fin 
1     format(1x,a50,'vort',i3.3,'.',i3.3)
      call movespa(fin,80)
      open(10,file=fin,status='unknown', form='unformatted')
         write(10) ux
         write(10) nsteps,time
      close(10)

      return
      end
! --------------------------------------------------
      subroutine input (ux,idp,id,nallgrp)

      include 'decay.h'
      complex,dimension(my, mmx2)::ux
      character iopath*50,fin*80
      common /iopathc/ iopath

      write(fin,1) iopath,idp,id
1     format(1x,a50,'vort',i3.3,'.',i3.3)
      call movespa(fin,80)
      open(10,file=fin,status='unknown',form='unformatted')
         read (10) ux
         read(10) nsteps0,time
      close(10)

      return
      end
! --------------------------------------------------
        subroutine newfft (ux,uxt, isign,ipx,ipy,wx,wy,id,nallgrp)
               
      include 'decay.h'
        complex, dimension (mx2, mmy):: uxt
        complex, dimension (my, mmx2):: ux
        real  ux1(0:mx-1)
        complex  uy1(0:my-1)
       
        real  wx (0:(mx2-1) ) 
        real  wy (0:(my-1) )
        integer ipx (0:mx2)         
        integer ipy (0:my)
        integer isign,i, j

       if (isign.eq.1) then     
!-- rc fft in x-dir-------------
        do j=1, mmy
           do i=1, mx2
              ux1(2*i-2) = real(uxt (i, j) )
              ux1(2*i-1) = aimag(uxt (i, j) )
           enddo

           call rdft (mx, isign, ux1, ipx, wx)
           ux1=ux1/float(mx2)
           ux1(1)=0.0

           do i=1, mx2
              uxt(i, j)=cmplx( ux1(2*i-2), ux1(2*i-1) )
           enddo
        enddo
        
        call transpose_xtok (uxt, ux, id, nallgrp)

!-- cc fft in y-dir-------------       
        do i=1, mmx2
            uy1=ux(:, i)
            call cdft (my*2, isign, uy1, ipy, wy)
            ux(:, i)= uy1/float(my)
         enddo    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
        else if (isign.eq. -1)  then  
         do i=1, mmx2
            uy1=ux(:, i)
            call cdft (my*2, isign, uy1, ipy, wy)
            ux(:, i)= uy1
         enddo

        call transpose_ktox (ux, uxt, id, nallgrp)
        do j=1, mmy
           do i=1, mx2
              ux1(2*i-2) = real(uxt (i, j) )
              ux1(2*i-1) = aimag(uxt (i, j) )
            enddo
            ux1(1)=0.0
            call rdft (mx, isign, ux1, ipx, wx)
            do i=1, mx2
               uxt(i, j)=cmplx( ux1(2*i-2), ux1(2*i-1))
            enddo
        enddo     
       endif

       return
       end
! -------------------------------------------------
      subroutine transpose_ktox (ux,uxt,id,nallgrp)
!     transpose ux to uxt so can do x-dir fft on uxt

      include 'decay.h'
      complex,dimension(my,mmx2)::ux
      complex,dimension(mx2,mmy)::uxt
      complex,dimension(mmy,mmx2)::tmp1,tmp2  
      complex,dimension(mmy,mx2)::tmp
      integer isize,nzm,nzp,status(MPI_STATUS_SIZE,2),req(2)
      integer i,j,k,js,j1,ks,k1, l
 
      if (nproc.gt.1) then
!        isize = mmy*mmx2
!        do i = 1,nproc-1
!         nzp=mod(id+i,nproc)
!         nzm=mod(id-i+nproc,nproc)
!         js = nzp*mmy
!
!         do j = 1,mmy
!              j1=js+j
!              tmp1( j, : ) = ux ( j1, :)
!         enddo
!         call mpi_isend(tmp1, isize, MPI_COMPLEX, nzp, i,
!     +        nallgrp,req(1),ierror)
!         call mpi_irecv(tmp2, isize, MPI_COMPLEX, nzm, i,
!     +        nallgrp,req(2),ierror)
!         call mpi_waitall (2,req,status,ierror)
!       
!         ks = nzm*mmx2
!         do k = 1, mmx2
!              k1 = ks+k
!              tmp( :, k1)= tmp2 ( :, k)
!         enddo
!       enddo
!
!!     does (id,id) spot from ux to tmp so can transpose to uxt
!
!       ks = id*mmx2
!       js = id*mmy
!         k1 = ks + k
!         do j = 1,mmy
!            j1 = js + j
!            tmp( j, k1 ) = ux ( j1, k )
!         enddo
!       enddo
!
!         do k=1, mmy
!            do j=1, mx2
!                uxt ( j, k ) = tmp(k, j )
!            enddo
!         enddo
       
       else if (nproc.eq.1) then
         do j=1, mx2
           do k=1, my
                uxt ( j, k)= ux (k, j ) 
           enddo
        enddo
       end if
       
      return
      end

!------------------------------------------------
      subroutine transpose_xtok (uxt,ux,id,nallgrp)
!     transpose uxt to ux so can do y-dir ifft on ux

      include 'decay.h'
      complex,dimension(my,mmx2)::ux
      complex,dimension(mx2,mmy)::uxt
      complex,dimension(mmx2, mmy)::tmp1,tmp2
      complex,dimension(mmx2,my)::tmp
      integer isize,nzm,nzp,status(MPI_STATUS_SIZE,2),req(2)
      integer i,j,k,js,j1,ks,k1, m

      if (nproc.gt.1) then
       isize = mmy*mmx2
       do i = 1,nproc-1
         nzp=mod(id+i,nproc)
         nzm=mod(id-i+nproc,nproc)
         js = nzp*mmx2
         do j = 1,mmx2
              j1 = js+j
              do m=1, mmy
                    tmp1(j, m) = uxt( j1, m) 
              enddo
         enddo
          call mpi_isend(tmp1, isize, MPI_COMPLEX, nzp, i,
     +        nallgrp,req(1),ierror)
         call mpi_irecv(tmp2, isize, MPI_COMPLEX, nzm, i,
     +        nallgrp,req(2),ierror)
         call mpi_waitall(2,req,status,ierror)

         ks = nzm*mmy
         do k = 1,mmy
              k1 = ks+k
              tmp( :, k1) = tmp2( :, k)
         enddo
       enddo

!     does the (id,id) spot from uxt to tmp so can transpose to ux
        ks = id*mmy
        js = id*mmx2
        do k = 1,mmy
            k1 = ks+k
           do j = 1,mmx2
              j1 = js+j
              tmp( j, k1) = uxt ( j1, k)
           enddo
       enddo
!---important  Transpose here!!
       do k = 1,mmx2
         do j = 1,my
            ux( j, k) = tmp( k, j)
         enddo
       enddo

      else if (nproc.eq.1) then
        do j=1, my  
          do k=1, mmx2  
             ux (j, k)= uxt (k, j) 
          enddo
        enddo
      end if

      return
      end
!--------------------------------------------
      subroutine update(c, my)
      complex,dimension(my, my/2)::c
	  real,dimension(my, my) :: x
	  
	  !do i=1,my
	  !    do j=1,my/2
	  !  	  x(i,2*j-1)=real(c(i,j))
	  !  	  x(i,2*j) = aimag(c(i,j))
	  !    enddo
	  !enddo

	  forall (i=1:my, j=1:my/2)
		  x(i,2*j-1) = real(c(i,j))
		  x(i,2*j) = aimag(c(i,j))
	  end forall

	  !=============================
	  ! Do the update job here
	  !	
	  do i=1,my/2
	      do j=1,my/2
	    	  x(i+my/2,j) = -x(i, my/2-j)
	      enddo
	  enddo

	  do i=my/2+1,my
	      do j=1,my
	    	  x(i, j) = x(i-my/2, j)
	      enddo
	  enddo
	  !
	  ! End of update job
	  !==============================

	  !do i=1,my
	  !    do j=1,my/2
	  !  	  c(i,j) = cmplx(x(2*i-1,j), x(2*i,j))
	  !    enddo
	  !enddo

	  forall (i=1:my, j=1:my/2)
		  c(i,j) = cmplx(x(2*i-1,j), x(2*i,j))
	  end forall
	  
	  end subroutine update
