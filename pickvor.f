!--------------------------------------------------------------
        subroutine pickvor0(vor,n,idn)
        real,dimension(n,n):: vor
        common /McWillimas/delt,crmax,armax,rmin,rmax,deltmax,epsmax,
     1                     vormin
        character iopath*50,fin*80
        common /iopathc/ iopath

!... 搜索全局极大值，并据此给出旋涡的最小涡量极值
        vormax = maxval(abs(vor))

!...  所有小于vormin的局部极值都将被去掉
!...  (McWillimas 建议的取值:约为全最大值的1/20，且不随时间改变)
        vormin = 0.05*vormax

!...  确定旋涡边界的参数
        delt  = 0.2

!...  旋涡的周长最大值为6.28, 即区域的一个边的长度
        crmax = amax1(float(n),256.0)

!...  旋涡的最大面积为 (0.1*2Pi)^2, 即整个区域的1%.
        armax = amax1(float(n*n)/100.0,1000.0)

!...  半径的最小值为网格尺度h, 除了在测试的初期, 此条件很少通不过. 
        rmin = 1.0

!...  轴对称情况测试一: 周长除面积 R=0.5*C/sqrt(Pi A) < Rmax
        rmax = 1.75

!...  轴对称情况测试二: 旋涡重心对中心点(涡量极值点)的偏离
        deltmax = 0.35

!...  轴对称情况测试三: 旋涡椭圆度最大值
        epsmax  = 2.5

!...  上面给出的是McWillimas所建议的, 有时针对具体情况, 可能需要修改.
        crmax = amax1(float(n),256.0)
        armax = 0.25*crmax*crmax/3.14
        rmax = 2.0
        deltmax = 1.0
        epsmax  = 2.5
        delt  = 0.3

10      format(1x,a50,'McWillimas',i3.3,'.dat')
        write(fin,10) iopath,idn
        call movespa(fin,80)
        open(18,file=fin)
          write(18,*) delt,crmax,armax,rmin,rmax,deltmax,epsmax,vormin
        close(18)
20      continue

        return
        end
!--------------------------------------------------------------
        subroutine pickvor0_input(idn)
        common /McWillimas/delt,crmax,armax,rmin,rmax,deltmax,epsmax,
     1                     vormin
        character iopath*50,fin*80
        common /iopathc/ iopath

10      format(1x,a50,'McWillimas',i3.3,'.dat')
        write(fin,10) iopath,idn
        call movespa(fin,80)
        open(18,file=fin,status='old',err=20)
          read(18,*,err=20,end=20) delt,crmax,armax,rmin,
     1                             rmax,deltmax,epsmax,vormin
        close(18)
        return

20      continue

        return
        end
!--------------------------------------------------------------
        subroutine pickvor(vor,n)
        real,dimension(n,n):: vor
        common /vortexdat/ ii0,jj0,gamv,arev,radv,epsv
        common /McWillimas/delt,crmax,armax,rmin,rmax,deltmax,epsmax,
     1                     vormin
        common /timeseq/ rnv,vormax,vorave,gamave,radave,areave,epsave,
     1                   vort(2),flatness

        vormax = maxval(abs(vor))

!... 搜索所有的局部极值, 搜索宽度为ld(一般取ld=2), 然后，判断这些局部极值
!... 是否是一个旋涡结构，如果是，则统计如下量：平均涡量极值，平均环量，平
!... 均半径，平均面积，平均椭圆度，涡量的数目为nv。

        vorave = 0.0
        gamave = 0.0
        radave = 0.0
        areave = 0.0
        epsave = 0.0

        nv = 0
        ld = 2
        do j=1,n
        do i=1,n
          vorij = vor(i,j)
          if(abs(vorij).lt.vormin) goto 10
          if(vorij.gt.0.0) then
!... 正的最大值
             do jj=j-ld,j+ld
               jj1 = mod(jj-1+n,n)+1
             do ii=i-ld,i+ld
               ii1 = mod(ii-1+n,n)+1
               if(i.ne.ii1.or.j.ne.jj1) then
                 if(vor(ii1,jj1).gt.vorij) goto 10
               endif
             enddo
             enddo
          else
!... 负的最小值
             do jj=j-ld,j+ld
               jj1 = mod(jj-1+n,n)+1
             do ii=i-ld,i+ld
               ii1 = mod(ii-1+n,n)+1
               if(i.ne.ii1.or.j.ne.jj1) then
                 if(vor(ii1,jj1).lt.vorij) goto 10
               endif
             enddo
             enddo
          endif

!... 程序运行到此处，说明(i,j)处为一个局部极值
          ii0 = i
          jj0 = j
!... 调用测试程序，判断此局部极值处是否有一个旋涡结构，极值位置是通过公共
!... 块传递的。
          call test_vor(vor,n,key)
          if(key.eq.1) then 
            nv = nv+1
            vorave = vorave + abs(vor(i,j))
            gamave = gamave + abs(gamv)
            radave = radave + radv
            areave = areave + arev
            epsave = epsave + epsv
          endif
10        continue
        enddo
        enddo
        if(nv.eq.0) nv=1
        rnv = float(nv)
        vorave = vorave/rnv
        gamave = gamave/rnv
        radave = radave/rnv
        epsave = epsave/rnv

!... 计算涡量的二阶矩和四阶矩以及平坦因子。
        sca = 1./float(n*n)
        vort(1) = sum((vor*sca)*(vor*sca))
        vort(2) = sum((vor*sca)**4)
        flatness = vort(2)/vort(1)/vort(1)

        return
        end
!-------------------------------------------------------------------
        subroutine test_vor(vor,n,key)
!... McWillimas 建议的取值，这些值将用来判断一个旋涡的合理性。
        common /McWillimas/delt,crmax,armax,rmin,rmax,deltmax,epsmax,
     1                     vormin
        integer,dimension(n*n)::iv,jv
        real vor(n,n)
        common /vortexdat/ ii0,jj0,gamv,arev,radv,epsv
        common /vortexbox/ ii,jj,ii1,jj1
		integer :: itest

!... 首先假设此极值点不是旋涡结构，如果程序中途返回，则表明不是旋涡。
        key = 0
        vormax = vor(ii0,jj0)

        ii = ii0
        jj = jj0

		itest = 0

!... 沿+x方向搜索第一个边界点 
10      continue           
        ii = mod(ii,n)+1                      ! This is infact "ii=ii+1"
        vorijm = vor(ii,jj)/vormax
!.... 如果出现新的极值，则忽略此极值
        if(vorijm.gt.1.0) return
!.... 此点仍然处于有效区间内，返回到10继续搜索下一点
        if(vorijm.gt.delt) goto 10

!..... 找到边界的第一个点(ii1,jj0)，搜索方向改为+y方向
        ii = mod(ii-2+n,n)+1                  ! this is ii=ii-1
!....如果只有一个孤立的点满足条件,则删除
!...  ii=ii0时，右点已经不属于旋涡，依次再判断左点，下点和上点。
        if(ii.eq.ii0) then
          nc = 0
          ii1 = mod(ii-2+n,n)+1
          if(vor(ii1,jj0)/vormax.lt.delt) nc=nc+1
          ii1 = mod(jj0-2+n,n)+1
          if(vor(ii0,ii1)/vormax.lt.delt) nc=nc+1
          ii1 = mod(jj0,n)+1
          if(vor(ii0,ii1)/vormax.lt.delt) nc=nc+1
          if(nc.eq.3) return
        endif

!....记录初始点(ii1,jj0)和初始搜索方向；
        ii1 = ii
        idx = 0
        idy = 1
!... 面积和周长都以网格尺度h为单位，
!... 于是得到的周长乘以h就是实际长度，面积则需要乘以h^2才得到
!... 实际面积。从找到的第一点开始，记录周长为0，点数为1。
        cir = 0
        nc = 1
        iv(nc) = ii
        jv(nc) = jj

21      continue           
!... 1. 判断右转90度的点是否属于此旋涡
        iit = mod(ii+idy-1+n,n)+1
        jjt = mod(jj-idx-1+n,n)+1
		!itest = itest + 1
        vorijm = vor(iit,jjt)/vormax
        if(vorijm.gt.delt) then
!......属于：
            cir = cir+1.0
!......回到了第一点，此旋涡的搜索过程结束，此为肯定出口。
            if(iit.eq.ii1.and.jjt.eq.jj0) goto 30
!......此旋涡的周长过长，不是一个正常的旋涡
            if(cir.gt.crmax) return
!......记录该点，当前方向右旋90度，继续搜索
            ii = iit
            jj = jjt
            nc = nc+1
            iv(nc) = ii
            jv(nc) = jj
            idxyt = idx
            idx = idy
            idy = -idxyt
            goto 21
        endif

22      continue           
!... 右旋90度的点不属于此旋涡，判断当前方向的下一个点
        iit = mod(ii+idx-1+n,n)+1
        jjt = mod(jj+idy-1+n,n)+1
		itest = itest + 1
        vorijm = vor(iit,jjt)/vormax
        if(vorijm.gt.delt) then
			itest = 0
!............. 属于：记录此点，并沿同方向继续搜索
               ii = iit
               jj = jjt
               cir = cir+1.0
               if(cir.gt.crmax) return
!......回到了第一点，此旋涡的搜索过程结束，此为肯定出口。
               if(ii.eq.ii1.and.jj.eq.jj0) goto 30
               nc = nc+1
               iv(nc) = ii
               jv(nc) = jj
               goto 21
        endif
!... 当前方向的下一个点不属于此旋涡，当前方向左转90度，继续搜索
        if(itest .gt. n*4) then
		    !write(*,*) 
			!write(*,*) 'There is a bug...'
			return
		endif
        idxyt = idx
        idx = -idy
        idy = idxyt
        goto 22

!... 搜索回到了初始点，搜索过程结束
30      continue

!... 调整数组iv,jv，使得所有点的标记在同一个连续的周期内(以ii0,jj0为准)
!..... 第一点必在中心极值点的右边
        if(iv(1).lt.ii0) iv(1)=iv(1)+n
!..... 后面的每一个点都以它的前一个点为参照，看是否需要修正
        do i=2,nc
              idx = iv(i)-iv(i-1)
              idy = jv(i)-jv(i-1)
              if(idx.gt. 1) iv(i) = iv(i)-n
              if(idx.lt.-1) iv(i) = iv(i)+n
              if(idy.gt. 1) jv(i) = jv(i)-n
              if(idy.lt.-1) jv(i) = jv(i)+n
        enddo

!... 计算一阶空间矩(重心位移)和二阶矩矩阵(这里要对整个旋涡的面积积分)
!... 这里对属于此集合的所有点求和(假设网格的尺度为1，网格面积也是1)
        ii = minval(iv(1:nc))
        jj = minval(jv(1:nc))
        ii1= maxval(iv(1:nc)) 
        jj1= maxval(jv(1:nc))
        arev = 0.0
        dx1 = 0.0
        dx2 = 0.0
        rm11 = 0.0
        rm22 = 0.0
        rm12 = 0.0
        vormax = vormax*1.01
        do j=jj,jj1
              idy = mod(j-1+n,n)+1
        do i=ii,ii1
              idx = mod(i-1+n,n)+1 
!... 测试(i,j)是否在区域内(包括边界上)，在idxyt=1，否则idxyt=0
              idxyt = 1
!... 1. 点(i,j)处于边界上，这样的点在统计环量时考虑，统计其他量时不
!...    考虑，可以避免把拖一个很细的长尾巴的涡排除。
!...    如果用goto 50一句，则是边界点也全部考虑。
              do jjt = 1,nc
                if(iv(jjt).eq.i.and.jv(jjt).eq.j) then
!                   gamv = gamv + vor(idx,idy)
!                   goto 51
                   goto 50
                endif
              enddo
!... 2. 点(i,j)处于边界点包围的区域内部
              call POLYGON_CP2DI(iv,jv,nc,real(i),real(j),INK)
              if(ink.gt.nc) goto 50
51            continue
              idxyt = 0
50            continue
              if(idxyt.eq.1) then
!....如果内部有更强的极值，则忽略此较小的极值
                 if(abs(vor(idx,idy)).gt.abs(vormax)) return
!....如果内部有不属于此旋涡的值，则忽略此较小的极值
                 if(vor(idx,idy)/vormax .le. 0.0) return
                 gamv = gamv + vor(idx,idy)
                 arev = arev + 1.0
                 di =  float(i-ii0)
                 dj =  float(j-jj0)
                 dx1 = dx1 + di
                 dx2 = dx2 + dj
                 rm11 = rm11 + di*di
                 rm22 = rm22 + dj*dj
                 rm12 = rm12 + di*dj
              endif
        enddo
        enddo
        dx1 = dx1/arev
        dx2 = dx2/arev
        rm11 = rm11/arev
        rm12 = rm12/arev
        rm22 = rm22/arev
        radv = sqrt(arev/3.1416)

        call plot_vortex(n,iv,jv,nc,14)


!....此旋涡的面积过大，不是一个正常的旋涡
        if(arev.gt.armax) return

!....找到的旋涡半径过小，不能构成一个合理的旋涡
        if(radv.lt.rmin) return

!... 周长和面积不相配，说明旋涡的形状不圆
        if(0.5*cir/sqrt(3.1416*arev).gt.rmax) return

!... 极值点和重心的偏移过大
        if( sqrt(dx1*dx1+dx2*dx2)/radv.gt.deltmax) return

!... 椭圆度过大

        cir = sqrt((rm11-rm22)*(rm11-rm22)+4.0*rm12*rm12)
        dx1 = rm11+rm22+cir
        dx2 = rm11+rm22-cir
        if(dx2.le.0.0) return
        epsv = sqrt(dx1/dx2-1.0)
        if(epsv.gt.epsmax) return

!... 通过了所有的测试，表明此极值所在位置是一个旋涡结构

        call plot_vortex(n,iv,jv,nc,15)

        key = 1

        return
        end
!-------------------------------------------------------------------------
        SUBROUTINE POLYGON_CP2DI(X,Y,N,XX,YY,INK)
        PARAMETER(PI2=0.7853981633974483,COE=0.6180339887498949)
        PARAMETER(DA=PI2*COE,ZERO=1.0E-6,UNIN=1.0-ZERO,ZEROB=1.2*ZERO)
        common /vortexbox/ ii,jj,ii1,jj1
        integer X(N),Y(N)

        INK = 0
        XMAX = ii1
        YMAX = jj1
        XMIN = ii
        YMIN = jj
        RR = AMAX1(XMAX-XMIN, YMAX-YMIN)
        X0 = .5*(XMIN+XMAX)
        Y0 = .5*(YMIN+YMAX)

        CT = 0.0
        IC = 0
100     CONTINUE
        CT = CT+DA
        IC = IC+1
        IF(MOD(IC,10).EQ.0) WRITE(*,101) IC
        if(ic.gt.100) then
          write(*,*) 'xx,yy=',xx,yy
	    return
        endif
101     FORMAT(1X,'IN POLYGON_CP2DI IC=',I3)
        XP = X0+RR*COS(CT)
        YP = Y0+RR*SIN(CT)
        XMIN = AMIN1(XX,XP)
        XMAX = AMAX1(XX,XP)
        YMIN = AMIN1(YY,YP)
        YMAX = AMAX1(YY,YP)
        IPCROSS = 0

        X1 = X(N)
        Y1 = Y(N)
        DO I=1,N
            IB = I-1
            IF(IB.EQ.0) IB=N
            X2 = X(I)
            Y2 = Y(I)
            IF( AMAX1(X1,X2).LT.XMIN .OR.
     &          AMIN1(X1,X2).GT.XMAX .OR.
     &          AMAX1(Y1,Y2).LT.YMIN .OR.
     &          AMIN1(Y1,Y2).GT.YMAX     ) GOTO 201

            DD = (X2-X1)*(YY-YP)-(Y2-Y1)*(XX-XP)
            IF(ABS(DD).LE.ZERO) THEN
!...  此时，1-2和(XX,YY)-(XP,YP)平行，需要计算两个平行线之间的距离。
              DD = (X1-XX)*(YP-YY)-(Y1-YY)*(XP-XX)
!...  注意，这里所用的判据是判断1-2边上的点1到(XX,YY)-(XP,YP)的距离，而不是
!...  计算到1-2的距离，这样做的好处是在1－2边的长度接近于零时，也不会有什么
!.... 障碍。根据这里所用的计算方法，线段(XX,YY)-(XP,YP)的长度总不会是零。
              IF(ABS(DD).LT.ZERO) THEN 
!.... 1点到(XX,YY)-(XP,YP)的距离是零，再判断1-2长度，只要它不是零，总能计算出来，
                DD = ABS(X2-X1)+ABS(Y2-Y1)
                IF(ABS(DD).LT.ZERO) THEN
!.... 如果1-2的距离也是零，说明1－2边几乎是一个点，此时只要被判断的点不是在此边
!.... 附近，那么通过移动P点可以避免这种不确定性，
                  DD = ABS(XX-X1)+ABS(YY-Y1)
                  IF(DD.LT.ZERO) THEN
                    INK = IB
                    RETURN
                  ELSE
                    GOTO 100
                  ENDIF
                ELSE
                  GOTO 100
                ENDIF
              ENDIF
            ELSE
               A1 = ((XX-X1)*(YY-YP)-(YY-Y1)*(XX-XP))/DD
               A2 = ((X2-X1)*(YY-Y1)-(Y2-Y1)*(XX-X1))/DD
              IF(A1.GE.ZERO.AND.A1.LE.UNIN) THEN
!1. (X1,Y1)--(X2,Y2)和(XX,YY)--(XP,YP)或其延长线相交。
!     1.1 (X1,Y1)--(X2,Y2)和(XX,YY)--(XP,YP)相交
                IF(A2.GE.ZERO.AND.A2.LE.UNIN) THEN
                  IPCROSS=IPCROSS+1
                ELSEIF(A2.LE.(-ZERO).OR.A2.GE.1.0+ZERO) THEN
!     1.2 (X1,Y1)--(X2,Y2)和(XX,YY)--(XP,YP)的延长线相交。

                ELSEIF(ABS(A2).LT.ZEROB) THEN 
!     1.3 (XX,YY)在(X1,Y1)--(X2,Y2)上。
                  INK = IB
                  RETURN
                ELSEIF(ABS(A2-1.0).LT.ZEROB) THEN 
!     1.4 (XP,YP)在(X1,Y1)--(X2,Y2)上，不应有此情况发生。
                  STOP 'ERROR HAPPEDED IN POLYGON_CP2DI'
                ELSE
                  STOP 'ERROR 1 IN POLYGON_CP2DI'                
                ENDIF
              ELSEIF(ABS(A1).LT.ZEROB) THEN
!2. (X1,Y1)在(XX,YY)--(XP,YP)或其延长线附近。
                IF(ABS(A2).LT.ZEROB) THEN
!     2.1 (X1,Y1)在(XX,YY)附近
                   INK = -IB
                   RETURN
                ELSEIF(A2.GT.ZERO.AND.A2.LT.UNIN) THEN
!     2.2 (X1,Y1)在(XX,YY)--(XP,YP)附近 
                   GOTO 100
                ELSEIF(A2.LE.-ZERO) THEN
!     2.3 (X1,Y1)--(X2,Y2)和(XX,YY)--(XP,YP)不相交 

                ELSE
                   STOP 'ERROR 2 IN POLYGON_CP2DI' 
                ENDIF
              ELSEIF(ABS(A1-1.0).LT.ZEROB) THEN
!3. (X2,Y2)在(XX,YY)--(XP,YP)或其延长线附近。
                IF(ABS(A2).LT.ZEROB) THEN
!     3.1 (X2,Y2)在(XX,YY)附近
                   INK = -I
                   RETURN
                ELSEIF(A2.GT.ZERO.AND.A2.LT.UNIN) THEN
!     3.2 (X2,Y2)在(XX,YY)--(XP,YP)附近 
                   GOTO 100
                ELSEIF(A2.LE.-ZERO) THEN
!     3.3 (X1,Y1)--(X2,Y2)和(XX,YY)--(XP,YP)不相交 

                ELSE
                   STOP 'ERROR 3 IN POLYGON_CP2DI' 
                ENDIF
              ELSEIF(A1.LE.(-ZERO).OR.A1.GE.1.0+ZERO) THEN
!4. (X1,Y1)--(X2,Y2)和(XX,YY)--(XP,YP)或其延长线不相交。
              ELSE
                 write(*,*) 'a1,a2=',a1,a2
                 STOP 'ERROR Last IN POLYGON_CP2DI'
              ENDIF

            ENDIF
201         CONTINUE
            X1 = X2
            Y1 = Y2
        ENDDO

        IF(MOD(IPCROSS,2).EQ.0) THEN
             INK = 0
        ELSE
             INK = N+1
        ENDIF

        RETURN
        END
