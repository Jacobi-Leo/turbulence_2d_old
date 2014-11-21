           subroutine mpi_init(ierr)
           ierr = 0
           return
           end

           subroutine mpi_finalize(ierr)
           ierr = 0
           return
           end

           subroutine mpi_comm_size(mpi,npr,ierr)
           ierr = mpi
           ierr = 0
           npr = 1
           return
           end

           subroutine mpi_comm_rank(mpi,id,ierr)
           ierr = mpi
           ierr = 0
           id = 0
           return
           end

           subroutine mpi_send(x,n,mp,idst,itag,mpi,req,ierr)
           integer x(1),req
           ierr = x(1)
           req = x(1)
           ierr = n
           ierr = mp
           ierr = idst
           ierr = itag
           ierr = mpi
           ierr = 0
           return
           end

           subroutine mpi_recv(x,n,mp,istart,itag,mpi,req,ierr)
           integer x(1),req
           ierr = x(1)
           req = x(1)
           ierr = n
           ierr = mp
           ierr = itag
           ierr = istart
           ierr = mpi
           ierr = 0
           return
           end

           subroutine mpi_isend(x,n,mp,idst,itag,mpi,req,ierr)
           integer x(1),req
           ierr = x(1)
           req = x(1)
           ierr = n
           ierr = mp
           ierr = idst
           ierr = itag
           ierr = mpi
           ierr = 0
           return
           end

           subroutine mpi_irecv(x,n,mp,istart,itag,mpi,req,ierr)
           integer x(1),req
           ierr = x(1)
           req = x(1)
           ierr = n
           ierr = mp
           ierr = itag
           ierr = istart
           ierr = mpi
           ierr = 0
           return
           end

           subroutine mpi_barrier(mpi,ierr)
           ierr = mpi
           ierr = 0
           return
           end

           subroutine mpi_waitall(n,req,status,ierr)
           integer req(1),status(1)
           ierr = n
           ierr = req(1)
           ierr = status(1)
           return
           end

           subroutine mpi_get_processor_name(name,len,ierr)
           character*(*) name
           name = 'only one processor'
           len = 18
           ierr = 0
           return
           end

           subroutine mpi_bcast(x,n,mp,id,mpi,ierr)
           integer x(1)
           ierr = x(1)
           ierr = n
           ierr = mp
           ierr = id
           ierr = mpi
           ierr = 0
           return
           end

           subroutine mpi_reduce(x,y,n,mp,mo,id,mpi,ierr)
           real x(1),y(1)
           do i=1,n
             y(i) = x(i)
           enddo
           ierr = n
           ierr = mp
           ierr = mo
           ierr = id
           ierr = mpi
           ierr = 0
           return
           end

           subroutine mpi_allreduce(x,y,n,mp,mo,mpi,ierr)
           real x(1),y(1)
           do i=1,n
             y(i) = x(i)
           enddo
           ierr = n
           ierr = mp
           ierr = mo
           ierr = mpi
           ierr = 0
           return
           end

c----------------------------------------------------
           real*8 function mpi_wtime()
           INTEGER II
           call run_time_msec(ii)
           mpi_wtime = DBLE(ii)
           return
           end



