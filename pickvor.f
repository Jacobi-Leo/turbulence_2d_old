c--------------------------------------------------------------
        subroutine pickvor0(vor,n,idn)
        real,dimension(n,n):: vor
        common /McWillimas/delt,crmax,armax,rmin,rmax,deltmax,epsmax,
     1                     vormin
        character iopath*50,fin*80
        common /iopathc/ iopath

c... ����ȫ�ּ���ֵ�����ݴ˸������е���С������ֵ
        vormax = maxval(abs(vor))

c...  ����С��vormin�ľֲ���ֵ������ȥ��
c...  (McWillimas �����ȡֵ:ԼΪȫ���ֵ��1/20���Ҳ���ʱ��ı�)
        vormin = 0.05*vormax

c...  ȷ�����б߽�Ĳ���
        delt  = 0.2

c...  ���е��ܳ����ֵΪ6.28, �������һ���ߵĳ���
        crmax = amax1(float(n),256.0)

c...  ���е�������Ϊ (0.1*2Pi)^2, �����������1%.
        armax = amax1(float(n*n)/100.0,1000.0)

c...  �뾶����СֵΪ����߶�h, �����ڲ��Եĳ���, ����������ͨ����. 
        rmin = 1.0

c...  ��Գ��������һ: �ܳ������ R=0.5*C/sqrt(Pi A) < Rmax
        rmax = 1.75

c...  ��Գ�������Զ�: �������Ķ����ĵ�(������ֵ��)��ƫ��
        deltmax = 0.35

c...  ��Գ����������: ������Բ�����ֵ
        epsmax  = 2.5

c...  �����������McWillimas�������, ��ʱ��Ծ������, ������Ҫ�޸�.
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
c--------------------------------------------------------------
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
c--------------------------------------------------------------
        subroutine pickvor(vor,n)
        real,dimension(n,n):: vor
        common /vortexdat/ ii0,jj0,gamv,arev,radv,epsv
        common /McWillimas/delt,crmax,armax,rmin,rmax,deltmax,epsmax,
     1                     vormin
        common /timeseq/ rnv,vormax,vorave,gamave,radave,areave,epsave,
     1                   vort(2),flatness

        vormax = maxval(abs(vor))

c... �������еľֲ���ֵ, �������Ϊld(һ��ȡld=2), Ȼ���ж���Щ�ֲ���ֵ
c... �Ƿ���һ�����нṹ������ǣ���ͳ����������ƽ��������ֵ��ƽ��������ƽ
c... ���뾶��ƽ�������ƽ����Բ�ȣ���������ĿΪnv��

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
c... �������ֵ
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
c... ������Сֵ
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

c... �������е��˴���˵��(i,j)��Ϊһ���ֲ���ֵ
          ii0 = i
          jj0 = j
c... ���ò��Գ����жϴ˾ֲ���ֵ���Ƿ���һ�����нṹ����ֵλ����ͨ������
c... �鴫�ݵġ�
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

c... ���������Ķ��׾غ��Ľ׾��Լ�ƽ̹���ӡ�
        sca = 1./float(n*n)
        vort(1) = sum((vor*sca)*(vor*sca))
        vort(2) = sum((vor*sca)**4)
        flatness = vort(2)/vort(1)/vort(1)

        return
        end
c-------------------------------------------------------------------
        subroutine test_vor(vor,n,key)
c... McWillimas �����ȡֵ����Щֵ�������ж�һ�����еĺ����ԡ�
        common /McWillimas/delt,crmax,armax,rmin,rmax,deltmax,epsmax,
     1                     vormin
        integer,dimension(n*n)::iv,jv
        real vor(n,n)
        common /vortexdat/ ii0,jj0,gamv,arev,radv,epsv
        common /vortexbox/ ii,jj,ii1,jj1

c... ���ȼ���˼�ֵ�㲻�����нṹ�����������;���أ�������������С�
        key = 0
        vormax = vor(ii0,jj0)

        ii = ii0
        jj = jj0

c... ��+x����������һ���߽�� 
10      continue           
        ii = mod(ii,n)+1                      ! This is infact "ii=ii+1"
        vorijm = vor(ii,jj)/vormax
c.... ��������µļ�ֵ������Դ˼�ֵ
        if(vorijm.gt.1.0) return
c.... �˵���Ȼ������Ч�����ڣ����ص�10����������һ��
        if(vorijm.gt.delt) goto 10

c..... �ҵ��߽�ĵ�һ����(ii1,jj0)�����������Ϊ+y����
        ii = mod(ii-2+n,n)+1                  ! this is ii=ii-1
c....���ֻ��һ�������ĵ���������,��ɾ��
c...  ii=ii0ʱ���ҵ��Ѿ����������У��������ж���㣬�µ���ϵ㡣
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

c....��¼��ʼ��(ii1,jj0)�ͳ�ʼ��������
        ii1 = ii
        idx = 0
        idy = 1
c... ������ܳ���������߶�hΪ��λ��
c... ���ǵõ����ܳ�����h����ʵ�ʳ��ȣ��������Ҫ����h^2�ŵõ�
c... ʵ����������ҵ��ĵ�һ�㿪ʼ����¼�ܳ�Ϊ0������Ϊ1��
        cir = 0
        nc = 1
        iv(nc) = ii
        jv(nc) = jj

21      continue           
c... 1. �ж���ת90�ȵĵ��Ƿ����ڴ�����
        iit = mod(ii+idy-1+n,n)+1
        jjt = mod(jj-idx-1+n,n)+1
        vorijm = vor(iit,jjt)/vormax
        if(vorijm.gt.delt) then
c......���ڣ�
            cir = cir+1.0
c......�ص��˵�һ�㣬�����е��������̽�������Ϊ�϶����ڡ�
            if(iit.eq.ii1.and.jjt.eq.jj0) goto 30
c......�����е��ܳ�����������һ������������
            if(cir.gt.crmax) return
c......��¼�õ㣬��ǰ��������90�ȣ���������
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
c... ����90�ȵĵ㲻���ڴ����У��жϵ�ǰ�������һ����
        iit = mod(ii+idx-1+n,n)+1
        jjt = mod(jj+idy-1+n,n)+1
        vorijm = vor(iit,jjt)/vormax
        if(vorijm.gt.delt) then
c............. ���ڣ���¼�˵㣬����ͬ�����������
               ii = iit
               jj = jjt
               cir = cir+1.0
               if(cir.gt.crmax) return
c......�ص��˵�һ�㣬�����е��������̽�������Ϊ�϶����ڡ�
               if(ii.eq.ii1.and.jj.eq.jj0) goto 30
               nc = nc+1
               iv(nc) = ii
               jv(nc) = jj
               goto 21
        endif
c... ��ǰ�������һ���㲻���ڴ����У���ǰ������ת90�ȣ���������
        idxyt = idx
        idx = -idy
        idy = idxyt
        goto 22

c... �����ص��˳�ʼ�㣬�������̽���
30      continue

c... ��������iv,jv��ʹ�����е�ı����ͬһ��������������(��ii0,jj0Ϊ׼)
c..... ��һ��������ļ�ֵ����ұ�
        if(iv(1).lt.ii0) iv(1)=iv(1)+n
c..... �����ÿһ���㶼������ǰһ����Ϊ���գ����Ƿ���Ҫ����
        do i=2,nc
              idx = iv(i)-iv(i-1)
              idy = jv(i)-jv(i-1)
              if(idx.gt. 1) iv(i) = iv(i)-n
              if(idx.lt.-1) iv(i) = iv(i)+n
              if(idy.gt. 1) jv(i) = jv(i)-n
              if(idy.lt.-1) jv(i) = jv(i)+n
        enddo

c... ����һ�׿ռ��(����λ��)�Ͷ��׾ؾ���(����Ҫ���������е��������)
c... ��������ڴ˼��ϵ����е����(��������ĳ߶�Ϊ1���������Ҳ��1)
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
c... ����(i,j)�Ƿ���������(�����߽���)����idxyt=1������idxyt=0
              idxyt = 1
c... 1. ��(i,j)���ڱ߽��ϣ������ĵ���ͳ�ƻ���ʱ���ǣ�ͳ��������ʱ��
c...    ���ǣ����Ա������һ����ϸ�ĳ�β�͵����ų���
c...    �����goto 50һ�䣬���Ǳ߽��Ҳȫ�����ǡ�
              do jjt = 1,nc
                if(iv(jjt).eq.i.and.jv(jjt).eq.j) then
c                   gamv = gamv + vor(idx,idy)
c                   goto 51
                   goto 50
                endif
              enddo
c... 2. ��(i,j)���ڱ߽���Χ�������ڲ�
              call POLYGON_CP2DI(iv,jv,nc,real(i),real(j),INK)
              if(ink.gt.nc) goto 50
51            continue
              idxyt = 0
50            continue
              if(idxyt.eq.1) then
c....����ڲ��и�ǿ�ļ�ֵ������Դ˽�С�ļ�ֵ
                 if(abs(vor(idx,idy)).gt.abs(vormax)) return
c....����ڲ��в����ڴ����е�ֵ������Դ˽�С�ļ�ֵ
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


c....�����е�������󣬲���һ������������
        if(arev.gt.armax) return

c....�ҵ������а뾶��С�����ܹ���һ�����������
        if(radv.lt.rmin) return

c... �ܳ�����������䣬˵�����е���״��Բ
        if(0.5*cir/sqrt(3.1416*arev).gt.rmax) return

c... ��ֵ������ĵ�ƫ�ƹ���
        if( sqrt(dx1*dx1+dx2*dx2)/radv.gt.deltmax) return

c... ��Բ�ȹ���

        cir = sqrt((rm11-rm22)*(rm11-rm22)+4.0*rm12*rm12)
        dx1 = rm11+rm22+cir
        dx2 = rm11+rm22-cir
        if(dx2.le.0.0) return
        epsv = sqrt(dx1/dx2-1.0)
        if(epsv.gt.epsmax) return

c... ͨ�������еĲ��ԣ������˼�ֵ����λ����һ�����нṹ

        call plot_vortex(n,iv,jv,nc,15)

        key = 1

        return
        end
C-------------------------------------------------------------------------
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
C...  ��ʱ��1-2��(XX,YY)-(XP,YP)ƽ�У���Ҫ��������ƽ����֮��ľ��롣
              DD = (X1-XX)*(YP-YY)-(Y1-YY)*(XP-XX)
C...  ע�⣬�������õ��о����ж�1-2���ϵĵ�1��(XX,YY)-(XP,YP)�ľ��룬������
c...  ���㵽1-2�ľ��룬�������ĺô�����1��2�ߵĳ��Ƚӽ�����ʱ��Ҳ������ʲô
c.... �ϰ��������������õļ��㷽�����߶�(XX,YY)-(XP,YP)�ĳ����ܲ������㡣
              IF(ABS(DD).LT.ZERO) THEN 
C.... 1�㵽(XX,YY)-(XP,YP)�ľ������㣬���ж�1-2���ȣ�ֻҪ�������㣬���ܼ��������
                DD = ABS(X2-X1)+ABS(Y2-Y1)
                IF(ABS(DD).LT.ZERO) THEN
C.... ���1-2�ľ���Ҳ���㣬˵��1��2�߼�����һ���㣬��ʱֻҪ���жϵĵ㲻���ڴ˱�
c.... ��������ôͨ���ƶ�P����Ա������ֲ�ȷ���ԣ�
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
C1. (X1,Y1)--(X2,Y2)��(XX,YY)--(XP,YP)�����ӳ����ཻ��
C     1.1 (X1,Y1)--(X2,Y2)��(XX,YY)--(XP,YP)�ཻ
                IF(A2.GE.ZERO.AND.A2.LE.UNIN) THEN
                  IPCROSS=IPCROSS+1
                ELSEIF(A2.LE.(-ZERO).OR.A2.GE.1.0+ZERO) THEN
C     1.2 (X1,Y1)--(X2,Y2)��(XX,YY)--(XP,YP)���ӳ����ཻ��

                ELSEIF(ABS(A2).LT.ZEROB) THEN 
C     1.3 (XX,YY)��(X1,Y1)--(X2,Y2)�ϡ�
                  INK = IB
                  RETURN
                ELSEIF(ABS(A2-1.0).LT.ZEROB) THEN 
C     1.4 (XP,YP)��(X1,Y1)--(X2,Y2)�ϣ���Ӧ�д����������
                  STOP 'ERROR HAPPEDED IN POLYGON_CP2DI'
                ELSE
                  STOP 'ERROR 1 IN POLYGON_CP2DI'                
                ENDIF
              ELSEIF(ABS(A1).LT.ZEROB) THEN
C2. (X1,Y1)��(XX,YY)--(XP,YP)�����ӳ��߸�����
                IF(ABS(A2).LT.ZEROB) THEN
C     2.1 (X1,Y1)��(XX,YY)����
                   INK = -IB
                   RETURN
                ELSEIF(A2.GT.ZERO.AND.A2.LT.UNIN) THEN
C     2.2 (X1,Y1)��(XX,YY)--(XP,YP)���� 
                   GOTO 100
                ELSEIF(A2.LE.-ZERO) THEN
C     2.3 (X1,Y1)--(X2,Y2)��(XX,YY)--(XP,YP)���ཻ 

                ELSE
                   STOP 'ERROR 2 IN POLYGON_CP2DI' 
                ENDIF
              ELSEIF(ABS(A1-1.0).LT.ZEROB) THEN
C3. (X2,Y2)��(XX,YY)--(XP,YP)�����ӳ��߸�����
                IF(ABS(A2).LT.ZEROB) THEN
C     3.1 (X2,Y2)��(XX,YY)����
                   INK = -I
                   RETURN
                ELSEIF(A2.GT.ZERO.AND.A2.LT.UNIN) THEN
C     3.2 (X2,Y2)��(XX,YY)--(XP,YP)���� 
                   GOTO 100
                ELSEIF(A2.LE.-ZERO) THEN
C     3.3 (X1,Y1)--(X2,Y2)��(XX,YY)--(XP,YP)���ཻ 

                ELSE
                   STOP 'ERROR 3 IN POLYGON_CP2DI' 
                ENDIF
              ELSEIF(A1.LE.(-ZERO).OR.A1.GE.1.0+ZERO) THEN
C4. (X1,Y1)--(X2,Y2)��(XX,YY)--(XP,YP)�����ӳ��߲��ཻ��
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
