	subroutine ptbp(modinput,delta,nbp,pft,inbiop,ibiomp,b4bp,sout,ldelta)
! ibiomp = liste des numeros de biomes pollen existant
! inbiop = numero du biome pollen de l échantillon étudié
! b4bp = matrice de transfert
	parameter (pi=3.1415926541,nbparam=4,nb4=28)
	dimension atemp(12),aprec(12),aclou(12),temp(12),prec(12)
	integer js(1),jw(1),h,inbiop,ibiomp(nbp)
	double precision modinput(45),delta(nbparam),ldelta,sout(11),b4bp(nb4,nbp),pft(nbp)
	double precision vars_in(50),output(500)
    real vth(nbp),soil(5)
	INTEGER OT,OL,OD
	COMMON IT,IL,OT,OL,ID,OD,ITST,IAUTO,ITST2
!	common /climat/atemp,aprec,aclou,soil,co2,xlon,xlat,sigma,zalti
!-----------------------------------------------------------------------------
!   real Anntemp,Annperc,Biotype,NPP,YLAT,CDveg,CDlit,CDsoil
!	real Area,CDtal,atotal,atoveg,atolit,atosoil
!-----------------------------------------------------------------------------
	js=maxloc(atemp(1:12))
	jw=minloc(atemp(1:12))
	sw=(maxval(atemp)-minval(atemp))/2

	it=8
	atemp(1:12)=modinput(1:12)
	aprec(1:12)=modinput(13:24)
	aclou(1:12)=modinput(25:36)
	soil(1:5)=modinput(37:41)
	co2=modinput(42)
	xlon=modinput(43)
	xlat=modinput(44)
	alti=modinput(45)
	js=maxloc(atemp(1:12))
	jw=minloc(atemp(1:12))
	sw=(maxval(atemp)-minval(atemp))/2
	vars_in(4)=soil(5)+1.1*min(delta(1),delta(2))  ! absolute min temp
	vars_in(2)=co2
	vars_in(3)=1e5-8.8*alti
	vars_in(49)=xlon
	vars_in(50)=0
	vars_in(1)=xlat
	vars_in(41:44)=soil(1:4)
	vars_in(45:48)=0

	gdd0=0
	do j=0,11
		vars_in(j+5)=atemp(j+1)+sin(0.262*J)*(delta(2)-delta(1))+delta(1)
		if(vars_in(42).eq.-9)then
			vars_in(j+5)=max(5.0,vars_in(j+5))
		endif
		gdd0=gdd0+max(0.,vars_in(j+5))*30
	enddo
	temp(1:12)=vars_in(5:16)

	if(gdd0.lt.150.or.vars_in(41).eq.-9)then
		ldelta=-999.0
		return
	endif
	vars_in(4)=vars_in(5)+(soil(5)-atemp(1))

!	boucle sur les precip
		do j=0,11
			prec(j+1)=aprec(j+1)*(sin(0.262*j)*(delta(4)-delta(3))+delta(3)+100)/100.0
			if(prec(j+1).lt.0.5)prec(j+1)=0.5
		enddo
		annperc=sum(prec(1:12))
		vars_in(17:28)=prec(1:12)
		
! calcul des ensoleillements
		if(maxval(delta(1:4)).lt.1.0.and.minval(delta(1:4)).gt.-1.0)then
			vars_in(29:40)=aclou(1:12)
		else
			vars_in(29)=49.6-0.083*prec(1)+0.68*temp(1)
			vars_in(30)=52.4-0.091*prec(2)+0.66*temp(2)
			vars_in(31)=50.9-0.074*prec(3)+0.69*temp(3)
			vars_in(32)=45.8-0.068*prec(4)+0.95*temp(4)
			vars_in(33)=38.3-0.067*prec(5)+1.3*temp(5)
			vars_in(34)=34.9-0.064*prec(6)+1.3*temp(6)
			vars_in(35)=36.3-0.063*prec(7)+1.2*temp(7)
			vars_in(36)=33.8-0.074*prec(8)+1.4*temp(8)
			vars_in(37)=33.4-0.090*prec(9)+1.6*temp(9)
			vars_in(38)=39.1-0.089*prec(10)+1.4*temp(10)
			vars_in(39)=44.8-0.092*prec(11)+1.04*temp(11)
			vars_in(40)=46.9-0.091*prec(12)+0.79*temp(12)
			vars_in(29:40)=max(vars_in(29:40),0.0)
			vars_in(29:40)=min(vars_in(29:40),100.0)
		endif

! calculation of the NPP
		output(1:500)=0
		call biome4(vars_in,output)
!---------------------------------------------------------------------------
!     Anntemp=output(458)
!	  Annperc=output(13)
!	   Bitype =output(1)
!	   NPP=output(3)*2.0   !carbon content is 0.5
!	   YLAT=output(49)

!	  CALL CARBON(Anntemp,Annperc,NPP,Bitype,YLAT,AREA,CDtal,CDveg,CDlit,&
!                	 CDsoil,atotal,atoveg,atolit,atosoil)
!---------------------------------------------------------------------------
! output
		ib4=output(1)
		sout(4)=output(3)   ! npp
!---------------------------------------------------------------------------
!		sout(5)=CDveg       ! 5 c veg
!		sout(6)=CDlit       ! 6 c lit
!		sout(7)=CDsoil      ! 7 c soil
!---------------------------------------------------------------------------
		sout(5:11)=output(331:337)  ! 'MTCO','MTWA','E/PE','P-E','GDD5','TANN','PANN'


		sout(1)=inbiop   ! biomee pollen obs
		sout(3)=output(1)   ! biome4 simul 
		js=maxloc(b4bp(ib4,1:nbp))
		sout(2)=js(1)    ! biome pollen simulé 
! calcul de la vraisemblance
		ds=0
		sup=maxval(b4bp(ib4,:))
		if(sup.gt.0.0)then
			do i=1,nbp
				dd1=b4bp(ib4,i)
				dd2=pft(i)
				ds=ds+(dd1-dd2)**2
			enddo
		else
			ds=999
		endif
		ldelta=-ds
	return
	end



	subroutine pt_c13(modinput,delta,sout,inbiom,c13,ldelta)
	parameter (pi=3.1415926541,nbparam=4)
	dimension atemp(12),aprec(12),aclou(12),temp(12),prec(12)
	double precision delta(nbparam),c13,sout(11),modinput(45),ldelta
	integer js(1),jw(1),h,inbiom(2),tbiom(28)
	dimension clou(12),c13(3)
	character li*50,fil*80,cident(100)*8,BIOM(28)*8
	real soil(5)
	double precision vars_in(50),output(500)
	INTEGER OT,OL,OD
	COMMON IT,IL,OT,OL,ID,OD,ITST,IAUTO,ITST2
!	common /climat/atemp,aprec,aclou,soil,co2,xlon,xlat,sigma,alti
	data biom/'TrEgFo','TrSeDeFo','TrDeFo',&
	'TeDeFo','TeCoFo','WaMxFo','CoMxFo','CoCoFo',&
	'ClMxFo','EgTaig','DeTaig',&
	'TrSav','TrXeShl',&
	'TeXsShl','TeXsWol','TeBlSa',&
	'OpCoWol','BoPrkl',&
	'TrGrl',&
	'TeGrl',&
	'Desert',&
	'StTund','ShTund','DShTund','PsShTund','CuTund',&
	'Barren','LIce'/
	data tbiom/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28/
!	data tbiom/1,1,1,  4,4,4,4,4,  5,5,5,  2,2,3,3,3,  5,6,  2,3,2,  6,6,6,6,6,  7,7/	
!	ci-dessus "temperate": 1:trop. forest, 2:trop. grass, 3:temp. grass, 4:temp. forest, 5: boreal forest, 6: boreal grass, 7: none

	it=8
	atemp(1:12)=modinput(1:12)
	aprec(1:12)=modinput(13:24)
	aclou(1:12)=modinput(25:36)
	soil(1:5)=modinput(37:41)
	co2=modinput(42)
	xlon=modinput(43)
	xlat=modinput(44)
	alti=modinput(45)

	js=maxloc(atemp(1:12))
	jw=minloc(atemp(1:12))
	sw=(maxval(atemp)-minval(atemp))/2
	vars_in(4)=soil(5)+1.1*min(delta(1),delta(2))  ! absolute min temp
	vars_in(2)=co2
	vars_in(3)=1e5-8.8*alti
	vars_in(49)=xlon
	vars_in(50)=0
	vars_in(1)=xlat
	vars_in(41:44)=soil(1:4)
	vars_in(45:48)=0

1 continue
	li=' '
	gdd0=0
	do j=0,11
		vars_in(j+5)=atemp(j+1)+sin(0.262*J)*(delta(2)-delta(1))+delta(1)
		if(vars_in(42).eq.-9)then
			vars_in(j+5)=max(5.0,vars_in(j+5))
		endif
		gdd0=gdd0+max(0.,vars_in(j+5))*30
	enddo
	temp(1:12)=vars_in(5:16)
	if(gdd0.lt.150.or.vars_in(41).eq.-9)then
		ldelta=-999.0
		return
	endif
	vars_in(4)=vars_in(5)+(soil(5)-atemp(1))

!	boucle sur les precip
	do j=0,11
		vars_in(j+17)=aprec(j+1)*(sin(0.262*j)*(delta(4)-delta(3))+delta(3)+100)/100.0
		if(vars_in(j+17).lt.2)vars_in(j+17)=2.0
	enddo
	prec(1:12)=vars_in(17:28)

! calcul des ensoleillements
	if(maxval(delta(1:4)).lt.1.0.and.minval(delta(1:4)).gt.-1.0)then
		vars_in(29:40)=aclou(1:12)
	else
		vars_in(29)=49.6-0.083*prec(1)+0.68*temp(1)
		vars_in(30)=52.4-0.091*prec(2)+0.66*temp(2)
		vars_in(31)=50.9-0.074*prec(3)+0.69*temp(3)
		vars_in(32)=45.8-0.068*prec(4)+0.95*temp(4)
		vars_in(33)=38.3-0.067*prec(5)+1.3*temp(5)
		vars_in(34)=34.9-0.064*prec(6)+1.3*temp(6)
		vars_in(35)=36.3-0.063*prec(7)+1.2*temp(7)
		vars_in(36)=33.8-0.074*prec(8)+1.4*temp(8)
		vars_in(37)=33.4-0.090*prec(9)+1.6*temp(9)
		vars_in(38)=39.1-0.089*prec(10)+1.4*temp(10)
		vars_in(39)=44.8-0.092*prec(11)+1.04*temp(11)
		vars_in(40)=46.9-0.091*prec(12)+0.79*temp(12)
		vars_in(29:40)=max(vars_in(29:40),0.0)
		vars_in(29:40)=min(vars_in(29:40),100.0)
	endif

! calculation of the C13
	output(1:500)=0
	call biome4(vars_in,output)
	sout(1)=output(1)
!	output(54) = C3+C4 / 55 = C3  / 56 = C4
	sout(2)=c13(3)-output(55)*(1+c13(3)/1000)
	sout(3)=amax1(output(57),0.0)   !  %NPPC4 
	sout(4)=output(3)   ! NPP
	sout(5:11)=output(331:337)

! calcul de la vraisemblance (on divise par l'mplitude max au carré)
	ds3=min(abs(inbiom(1)-tbiom(int(sout(1)))),abs(inbiom(2)-tbiom(int(sout(1))))) 
!	if(ds3.gt.0.1)then
!		delta(1)=delta(1)+boxmu(bidon)
!		delta(2)=delta(2)+boxmu(bidon)
!		delta(3)=delta(3)*(1+boxmu(bidon)*0.10)
!		delta(4)=delta(4)*(1+boxmu(bidon)*0.10)
!		write(8,*)delta(1:4)
!		goto 1
!	endif
	ds=(c13(1)-sout(2))**2 
	ds2=0
	if(c13(2).ne.-999.0)ds2=0.005*(sout(3)-c13(2))**2
	ldelta=-ds-ds2-ds3**2
!	else
!		ldelta=-999.0
!	endif
return
end

      FUNCTION BOXMU(BIDON)
!
!  Lebart et al., 1979. Dunod, Paris
!
!  VARIABLE PSEUDO ALEATOIRE NORMALE (0,1) --- BOX ET MULLER, 1

      DATA CC/6.28318501/
      DATA Z1/0.0/,Z2/0.0/,KK/1/
      IF(KK.EQ.0)GOTO 10
      U1=SEN3A(BIDON)
      U2=SEN3A(BIDON)
      P=SQRT(ABS(-2.0*ALOG(U1)))
      Q=CC*U2
      Z1=P*SIN(Q)
      Z2=P*COS(Q)
      BOXMU=Z1
      GOTO 20
10    BOXMU=Z2
20    KK=1-KK
      RETURN
      END


      FUNCTION SEN3A(BIDON)
!  TIRAGE PSEUDO-ALEATOIRE UNIFORME SUR (0,1)
!  Lebart et al., 1979. Dunod, Paris
      DATA M12/4096/
      DATA F1/2.44140625E-04/,F2/5.96046448E-08/,F3/1.45519152E-11/
      DATA J1/3823/,J2/4006/,J3/2903/
      DATA I1/3823/,I2/4006/,I3/2903/
      K3=I3*J3
      L3=K3/M12
      K2=I2*J3+I3*J2+L3
      L2=K2/M12
      K1=I1*J3+I2*J2+I3*J1+L2
      L1=K1/M12
      I1=K1-L1*M12
      I2=K2-L2*M12
      I3=K3-L3*M12
      SEN3A=F1*FLOAT(I1)+F2*FLOAT(I2)+F3*FLOAT(I3)
      RETURN
      END


! borns
      SUBROUTINE BORNS(N,X,XMIN,XMAX)
      DIMENSION X(1)
      XMIN=X(1)
      XMAX=X(1)
      IF(N.EQ.1)RETURN
      DO I=2,N
      IF(X(I).LT.XMIN)XMIN=X(I)
      IF(X(I).GT.XMAX)XMAX=X(I)
      ENDDO
      RETURN
      END


! quant
      SUBROUTINE QUANT(ICARD,V,QT,NB,KV)
!
!  CALCUL DANS QT(NB) DES 1/(NB+1) QUANTILES DE V(ICARD)
!    EXEMPLE POUR NB=3 QUANTILES 0.25 0.50 0.75
!            POUR NB=19          0.05 0.10 0.15 ... 0.90 0.95
!    LEBART, MORINEAU, FENELON, 1979, DUNOD, PARIS
      DIMENSION V(1),KV(1),QT(1)
      CALL SHELK(ICARD,V,KV)
      NCLAS=NB+1
      DO K=1,NB
      I1=0.9999+FLOAT(K*ICARD)/NCLAS
      I2=1.0000+FLOAT(K*ICARD)/NCLAS
      QT(K)=(V(I1)+V(I2))/2.
      enddo
      RETURN
      END


!fuzzyk
      subroutine fuzzyk(x,n,xf,kc,xb)
      dimension x(n),xf(n,1),xb(1)
!  dimension xf(n,kc),xb(kc)
!  abs(kc) = nb categories
!  if kc>0 the categories limits are calculated with this dataset
!  if kc<0 the categories limits must be provided in xb(1) --> xb(abs(kc))
      kcat=iabs(kc)
      if(kcat.gt.10)stop'Too many categories'
      if(kcat.lt.2)stop'Too few categories'
      if(kc.gt.0)call borns(n,x,xb(1),xb(kcat))
      xs=(xb(kcat)-xb(1))/(kcat-1)
      do i=2,kcat-1
        xb(i)=xb(i-1)+xs
      enddo
      do i=1,n
!        sum=0
        do k=1,kcat
          if(xs.gt.0.001)then
            xf(i,k)=exp(-((x(i)-xb(k))/xs)**2)
!            sum=sum+xf(i,k)
          else
            xf(i,k)=0
          endif
        enddo
        if(x(i).lt.xb(1))then
          do k=1,kcat
            xf(i,k)=0
          enddo
          xf(i,1)=1
        endif
        if(x(i).gt.xb(kcat))then
          do k=1,kcat
            xf(i,k)=0
          enddo
          xf(i,kcat)=1
        endif
      enddo
      return
      end subroutine fuzzyk


!
!  nnapp3
!
	subroutine nnapp3(nom,dd,ident,mt,doo,mp)
	PARAMETER (IM3=1000,IM8=10000)
	double precision DOO(mt),DD(mp)
	CHARACTER nom*255,ident(mt)*8,li*80,nom2*6
	allocatable:: dd2(:),w1m(:),dd1(:),y(:),w1(:),w(:),x0(:),y0(:),nv(:)
	character,allocatable:: dent(:)*8
	COMMON IT,IL,OT,OL,ID,OD,ITST,IAUTO,ITST2
	INTEGER OT,OL,OD
    common /tp/me,errg,rl,rli,rld,cmoment,er,nstd,itransfer
	data nom2/'nn.tmp'/
!
! initialisation
!
	Allocate(DD2(IM3),W1M(IM3),DD1(IM3),Y(IM8),W1(IM3),W(IM3),x0(IM3),y0(IM3),nv(IM3),dent(im3))
	open(98,file=nom,status='old',iostat=iosl)
	open(it,file=nom2,status='replace',iostat=k)
	if(iosl.ne.0)write(*,*)' NN file (',nom2(1:len_trim(nom2)),')  NOT FOUND'
	nnnet=0
	nli=0
	dmin=1.e+10
	imin=0
	do while (iosl>=0)
		read(98,'(a)',end=98)li
		nli=nli+1
		if(li(1:8).eq.'X:Centre')then
			nnnet=nnnet+1
			if(nnnet.gt.im3)Stop 'Too many NN sets'
			nv(nnnet)=nli+1
			backspace 98
			read(98,*)li,x0(nnnet),xmin,xmax 
			read(98,*)li,y0(nnnet),ymin,ymax
			if(li(1:8).ne.'Y:Centre')Stop 'Bad NN file'
			nli=nli+1
			dst=(dd(1)-x0(nnnet))**2+(dd(2)-y0(nnnet))**2
			if(dd(1).lt.xmin.or.dd(1).gt.xmax)cycle
			if(dd(2).lt.ymin.or.dd(2).gt.ymax)cycle
			if(dst.lt.dmin)then
				dmin=dst
				imin=nnnet
			endif
		endif
	enddo
98 	if(nnnet.gt.0)then
		nv(nnnet+1)=nli+1
		rewind 98
		read(98,*)m,mo,k1,ncas,nstd,itransfer
		read(98,*)npx,npy
		write(it,*)m,mo,k1,ncas,nstd,itransfer
		write(it,*)npx,npy
		mp=0
		if(imin.eq.0)Stop 'No NN set found'
		do i=3,nv(imin+1)
			read(98,'(a)',end=1)li
			if(i.gt.nv(imin))write(it,'(a)')li(:len_trim(li))
		enddo
1		close(98)
	else
		rewind 98
		do while (iosl>=0)
			read(98,'(a)',end=2)li	
			write(it,'(a)')li(:len_trim(li))
		enddo
2		close(98)
      endif
	close(it)
	open(it,file=nom2,status='old')
	read(it,*)m,mo,k1,ncas,nstd,itransfer
	read(it,*)npx,npy
	do j=1,m+mo
		read(it,'(5X,A8,25F10.3)')dent(j),dd1(j),dd2(j)
	enddo
	close(it)
!
!    recherche des predicteurs
!
	do j=1,m
	  kk=0
	  w1m(j)=0
	  kl=len_trim(dent(j))
	  do k=1,mt
			if(dent(j)(1:kl).eq.ident(k)(1:len_trim(ident(k))))kk=k
	  enddo
	  if(kk.eq.0)then
!			write(ol,*)'Var. ',dent(j),' is not in the input data'
!			write(ol,*)'We consider it as constantly 0'
!			write(*,*)'Var. ',dent(j),' is not in the input data'
			kk=j
	  endif
	  w1m(j)=dd(kk)
	enddo
!
!  standardisation normale des regresseurs
!
	if(nstd.le.1)then
		w1(1:m)=(w1m(1:m)-dd1(1:m))/dd2(1:m)
	else
		write(*,*)'standardization not allowed'
	endif
!
!  estimation des predictands
!
	if(ncas.gt.1)then
		ncas=1
		write(*,*)'Bootstrap estimate not allowed - ncas set to 1'
	endif
	call app_nnet(w1,1,m,mo,ncas,k1,nom2)
	do j=1,mo
		ident(mt+j)=dent(m+j)
		doo(j)=w1(j)*dd2(m+j)+dd1(m+j)
		if(nstd.eq.-1)doo(j)=10.0**doo(j)-1.0
	enddo
	mp=mo
	DeAllocate(DD2,W1M,DD1,Y,W1,W,X0,Y0,nv,dent,stat=ierror)
	return
	end subroutine nnapp3


!
! app_nnet
!
      subroutine app_nnet(dd,n,m,mo,ncas,k1,nom)
      dimension dd(n,1),qt(20)
      real, allocatable:: s0(:),s1(:),B1(:),B2(:),W1(:,:),W2(:,:),y(:,:),in(:),cp(:,:)
      character nom*(*),ent*80
      integer ot,ol,od
      COMMON IT,IL,OT,OL,ID,OD,ITST
      common /tp/me,errg,rl,rli,rld,cmoment,er,nstd,itransfer
!
!  initialisation
!
      mrg=m
      if(nstd.gt.2)mrg=m*nstd
      m3=max0(mo,k1,mrg)
      Allocate(s0(m3),y(ncas,mo),in(ncas),s1(m3),w1(mrg,k1),w2(k1,mo),b1(k1),b2(mo),cp(mo,mo),stat=ierror)    
      mvd=mo
      npy=1
      if(nstd.eq.2)then
         open(it,file=nom,status='old')
         read(it,'(a)')ent
         read(it,*)npx,npy
         do j=1,m
           read(it,'(a)')ent
         enddo
         do j=1,mo
           read(it,'(29X,23F10.3)')(cp(j,k),k=1,npy)
         enddo
         mrg=npx
         mo=npy
         close(it)
      endif
      nq=19
      if(ncas.le.20)nq=9
      if(n.gt.1)write(*,*)
!
!  boucle sur les observations
!
      do i=1,n
        if(n.gt.1)write(*,'(1H+,I6)')I
        open(it,file=nom,status='old')
        do j=1,m+mvd+2
          read(it,'(a)')ent
        enddo
!
!   boucle sur les simulations bootstrap
!
        do l=1,ncas
           read(it,*)((w1(j,k),k=1,k1),j=1,mrg)
           read(it,*)(b1(k),k=1,k1)
!   mise a zero des connexions bloquees
           do k=1,k1
             do j=1,mrg
               if(w1(j,k).eq.-99.0)w1(j,k)=0
             enddo
             if(b1(k).eq.-99.0)b1(k)=0
           enddo
           read(it,*)((w2(j,k),k=1,mo),j=1,k1)
           read(it,*)(b2(k),k=1,mo)
!   mise a zero des connexions bloquees
           do k=1,mo
             do j=1,k1
               if(w2(j,k).eq.-99.0)w2(j,k)=0
             enddo
             if(b2(k).eq.-99.0)b2(k)=0
           enddo
!
!  calcul des outputs
!
           do j=1,mrg
             s0(j)=dd(i,j)
           enddo
           call slogsig(w1,s0,b1,1,mrg,k1,s1)
           if(itransfer.eq.0)then
             call spurlin(w2,s1,b2,1,k1,mo,s0)
           else
             call slogsig(w2,s1,b2,1,k1,mo,s0)
           endif
           if(nstd.ne.2)then
             do k=1,mo
               y(l,k)=s0(k)
             enddo
           else
             do j=1,mvd
               y(l,j)=0
               do k=1,mo
                 y(l,j)=y(l,j)+cp(j,k)*s0(k)
               enddo
             enddo
          endif
        enddo        
        close(it)
!
!  sorties: calcul des triplets/singletons
!
        do k=1,mvd
          if(ncas.gt.3)then
            call quant(ncas,y(1,k),qt,nq,in)
            dd(i,(k-1)*3+1)=qt(1)
            dd(i,(k-1)*3+2)=qt((nq+1)/2)
            dd(i,(k-1)*3+3)=qt(nq)
          else
            if(ncas.le.1)then
              dd(i,k)=y(1,k)
            else
              call shelk(3,y(1,k),in)
              dd(i,(k-1)*3+1)=y(1,k)
              dd(i,(k-1)*3+2)=y(2,k)
              dd(i,(k-1)*3+3)=y(3,k)
            endif
          endif
        enddo
      enddo
      mo=mvd
      DeAllocate(s0,y,in,s1,w1,w2,b1,b2,cp,stat=ierror)    
      return
      end subroutine app_nnet
!
! shelk
!
      SUBROUTINE SHELK(N,X,KX)
!
!  RANGE LE VECTEUR X(N) EN ORDRE CROISSANT DANS X(N)
!  KX(J) DONNE LES POSITIONS INITIALES
!    J. BOOTHROYD/SHELLSORT ALGORITHM.201/COMM.ACM/VOL 6 (1963)
!    NO 8, 445/
!    D.A. SHELL/A HIGH-SPEED SORTING PROCEDURE/COMM.ACM/VOL2, 1959, 30-3
!  See  Lebart et al., 1979. Dunod, Paris
!
      DIMENSION X(N),KX(N)
      DO J=1,N
      	KX(J)=J
      enddo
      I=1
   20 I=I+I
      IF(I.LE.N)GOTO 20
      M=I-1
   30 M=M/2
      IF(M.EQ.0)GOTO 70
      K=N-M
      DO 60 J=1,K
      JM=J+M
   40 JM=JM-M
      IF(JM.LE.0)GOTO 60
   50 L=JM+M
      IF(X(L).GE.X(JM))GOTO 60
      PIV=X(JM)
      X(JM)=X(L)
      X(L)=PIV
      KPIV=KX(JM)
      KX(JM)=KX(L)
      KX(L)=KPIV
      GOTO 40
   60 CONTINUE
      GOTO 30
   70 RETURN
      END 

!
! slogsig
!
      subroutine slogsig(w1,d,b1,n,m,k1,a1)
!
!   fonction sigmoid de k combili de m inputs
!    done on n observations
!
      dimension d(n,m),w1(m,k1),b1(k1),a1(n,k1)
	real*8 x
      do i=1,n
        do k=1,k1
          x=sum(d(i,:)*w1(:,k))+b1(k)
          if(dabs(x).lt.20.)then
            x=1./(1+exp(-x))
          else
            x=0.
	    if(x.gt.0.)x=1.
          endif
          a1(i,k)=x
        enddo
      enddo
      return
      end

!
! spurlin
!
      subroutine spurlin(w1,d,b1,n,m,k1,a1)
!   fonction tansigmoid de k combili de m inputs
!    done on n observations
      dimension d(n,m),w1(m,k1),b1(k1),a1(n,k1)
      do i=1,n
        do k=1,k1
          a1(i,k)=sum(d(i,:)*w1(:,k))+b1(k)
        enddo
      enddo
      return
      end
!
! stansig
!
      subroutine stansig(w1,d,b1,n,m,k1,a1)
!
!   fonction tansigmoid de k combili de m inputs
!    done on n observations
!
      dimension d(n,m),w1(m,k1),b1(k1),a1(n,k1)
	real*8 x
      do i=1,n
        do k=1,k1
          x=sum(d(i,:)*w1(:,k))+b1(k)
          if(dabs(x).lt.10.)then
            x=2./(1+exp(-2*x))-1
          else
            x=-1
				if(x.gt.0)x=1
          endif
          a1(i,k)=x
        enddo
      enddo
      return
      end





!    The BIOME4-system:	biome4.f	4.2b1	22.10.99
!
!	Copyright (c) 1999 by Jed O. Kaplan
!     
!	See COPYING file for copying and redistribution conditions.
!
!	This program is free software; you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation; version 2 of the License.
!
!	This program is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	Contact info: jkaplan@bgc-jena.mpg.de
!---------------------------------------------------------------------------
!
!				B I O M E 4 . F
!
!---------------------------------------------------------------------------
!      This is BIOME4, based on BIOME3 (Haxeltine and Prentice 1996).
!
!      All argguments to and from this subroutine are passed through two arrays: 
!
!              "vars_in" and "output"
!
!      You can customize the contents of these arrays to suit your own needs
!      based on the assignments in the beginning of the program and the 
!      output assignments in the subroutine "competition2".  I suggest each
!      user write his own driver software.
!
!      The basic input variables required by this model are:
!      Monthly mean- temperature, precipitation, %cloudiness (or sunshine)
!      Absolute minimum temperature
!      Water holding capacity (mm/m) in the top 30cm of the soil
!      Water holding capacity (mm/m) in the rest of the soil
!      Conductivity index of water through the soil column (see the soil
!      and hydrology subroutines for more information).
!
!      Caveat: compared to the original BIOME3 there is a lot of added
!      functionality in this model, and in some cases complexity.  The code
!      may neither compile nor run on certain systems (it works well now
!      compiled with g77 and running on an Intel Linux machine).
!
!      Please direct your questions to me at jed.kaplan@bgc-jena.mpg.de
!
!      By Jed O. Kaplan 1994-1999.
!
!      For further reference see:
!
!	Haxeltine and Prentice, I.C. 1996. BIOME3: An equilibrium 
!	terrestrial biosphere model based on ecophysiological 
!	constraints, resource availibility and competition among
!  	plant functional types, Global Biogeochemical Cycles, 10(4) 693-709.
!
!      Please remove this header and change the name of your new model
!      if you make ANY alterations to the code!
!
!      Author:	Jed O. Kaplan
!      Date:	22 October 1999
!      Version:	v4.2b1
!      Revised:
!
!------------------------------------------------------------------------

      subroutine biome4(vars_in,output)

      implicit none

      integer numofpfts

      parameter(numofpfts=13)

      logical diagmode

      integer optdata(0:numofpfts,500)
      integer biome,pft
      integer pfts(numofpfts),iopt
      integer i

      real lon,lat,co2,p,tdif,temp(12),prec(12),clou(12),soil(5)
      real dtemp(365),dprec(365),dclou(365),dprecin(365)
      real dpet(365),dphen(365,2),dmelt(365),realout(0:numofpfts,200)
      real optnpp(0:numofpfts),optlai(0:numofpfts)
      real pftpar(25,25),k(12),dayl(12),sun(12),wetness(0:numofpfts)
      real gdd5,tcm,twm,tprec,gdd0,rad0,tmin,tminin
      real tsoil(12),ddayl(365),maxdepth
      real radanom(12),alttmin

      double precision vars_in(50),output(500)

!----------------------------
!     assign the variables that arrived in the array vars_in
      lon=vars_in(49)
      lat=vars_in(1) !vars_in(49)-(vars_in(1)/vars_in(50))
      co2=vars_in(2)
      p=vars_in(3)
      tminin=vars_in(4)
      
      do i=1,12
       temp(i)=vars_in(4+i)
       prec(i)=vars_in(16+i)
       clou(i)=vars_in(28+i)
      end do
       
      do i=1,4
       soil(i)=vars_in(40+i)
      end do
      
      iopt=nint(vars_in(46))

      if (iopt.eq.1.) then
       diagmode=.true.
      else
       diagmode=.false.
      end if
      
!----------------------------
 
!      set a dummy rad anomaly (not used in this version)
       do i=1,12
        radanom(i)=1.0
       end do

!-------------------------------------------------------------------------
!      Reset the output matrix
       do pft=0,numofpfts
        do i=1,500
         optdata(pft,i)=0
        end do
       end do
!-------------------------------------------------------------------------
!      Initialize soil texture specific parameters
       k(1)=soil(1)
       k(2)=soil(2)
       k(5)=soil(3)
       k(6)=soil(4)
!      call soildata(k,soil)
!-------------------------------------------------------------------------
!      Linearly interpolate mid-month values to quasi-daily values:
       call daily(temp,dtemp)
       call daily(clou,dclou)   
       call daily(prec,dprecin)   
!--------------------------------------------------------------------------
!      Initialize parameters derived from climate data:
       call climdata(tcm,twm,gdd5,gdd0,tprec,temp,prec,dtemp,alttmin)
       call soiltemp(temp,soil,tsoil)
!--------------------------------------------------------------------------
!      Calculate mid-month values for pet,sun & dayl from temp,cloud & lat:
       call ppeett(lat,dtemp,dclou,dpet,temp,sun,dayl,rad0,ddayl,radanom)  
!-------------------------------------------------------------------------
!      Run snow model:
       call snow(dtemp,dprec,dmelt,dprecin,maxdepth)   
!-------------------------------------------------------------------------

!      Initialize the evergreen phenology
       do i=1,365
        dphen(i,1)=1.0
        dphen(i,2)=1.0
       end do

!      Initialize pft specific parameters
       call pftdata(pftpar)   

!--------------------------------------------------------------------------
!      Rulebase of absolute constraints to select potentially presents pfts:  
        call constraints(tcm,twm,tminin,gdd5,rad0,pfts,tmin,maxdepth,gdd0)

!       The tropical evergreen pft is not used in this version
!       of the model.  This is because the tropical deciduous tree
!       will be evergreen if it is not subject to water stress.
!       Otherwise the two pft's are parameterized in the same way,
!       so not using the pft saves computation time.

!        pfts(1)=0

!--------------------------------------------------------------------------

      if (diagmode) then
      
       write(*,'(A,F6.1,A,F6.1,A)')'Tcm and Tmin are:',tcm,' and',tmin,' degrees C respectively.'
       write(*,'(A,F8.1,A,F8.1,A)')'GDD5 is:',gdd5,' and total annual precip is:',tprec,' mm.'
       write(*,'(A,F8.1,A)')'Maximum snowdepth is:',maxdepth*10.,' mm.'
       write(*,'(A,2F7.2,3F7.1)')'The soil parameters are:',soil
       write(*,*)'The following PFTs will be computed:'

      end if

!--------------------------------------------------------------------------
!     Calculate optimal LAI & NPP for the selected pfts:
        
      do pft=1,numofpfts

       optlai(pft)=0.0
       optnpp(pft)=0.0

       if (pfts(pft).ne.0) then

        if (pftpar(pft,1).ge.2) then
!        Initialize the generic summergreen phenology
         call phenology(dphen,dtemp,temp,tcm,tdif,tmin,pft,ddayl,pftpar)
        end if
   call findnpp(pfts,pft,optlai(pft),optnpp(pft),tprec,dtemp,sun,temp,dprec,dmelt,dpet,dayl, &
  k,pftpar,optdata,dphen,co2,p,tsoil,realout,numofpfts)
       end if
      end do
!------------------------------------------------------------------------------
!      Select dominant plant type/s on the basis of modelled optimal NPP & LAI:
       
       call competition2(optnpp,optlai,wetness,tmin,tprec,pfts,optdata,output,diagmode, &
     biome,numofpfts,gdd0,gdd5,tcm,pftpar,soil)

!------------------------------------------------------------------------------
!      Final output biome is given by the integer biome:
 
       output(1)=biome
	   output(453)=tcm
	   output(454)=twm
	   output(455)=tmin
	   output(456)=gdd0
	   output(457)=gdd5
	   output(458)=sum(temp(1:12))/12.0
	   output(331:332)=output(453:454)
	   output(333:334)=output(10:11)
	   output(335:336)=output(457:458)
	   output(337)=output(13)

       output(48)=lon
       output(49)=lat
       
       do pft=1,numofpfts
       
        output(300+pft)=nint(optnpp(pft))
        output(300+numofpfts+pft)=nint(optlai(pft)*100.0)
       
!-----------------------
        if (pfts(pft).ne.0) then
        if (diagmode) then
!       type some diagnostic output here
         write(*,10)pft,optlai(pft),optnpp(pft),wetness(pft),((optdata(pft,i)/10.),i=37,48) 
10       format(I3,F5.2,F7.1,F6.1,12F6.1)
        end if
        end if

!-----------------------

       end do
       
      if (diagmode) then
       write(*,*)'press return to continue'
       read(*,*)
      end if
     
      return
      end

!*******************************************************************************
!a la fin, modif Hatte

      subroutine competition2 &
     (optnpp,optlai,wetness,tmin,tprec,pfts,optdata,output,diagmode, &
      biome,numofpfts,gdd0,gdd5,tcm,pftpar,soil)

      implicit none

      integer numofpfts

      real optnpp(0:numofpfts),optlai(0:numofpfts),maxnpp,maxlai
      real maxdiffnpp,woodylai,grasslai,driest(0:14)  !numofpfts
      real tmin,tcm,pftpar(25,25),soil(5),mwet,wettest(0:14)
      real temperatenpp,tprec,wetness(0:numofpfts),lai,npp
      real woodnpp,grassnpp,subnpp,gdd0,gdd5,ratio,lairatio,nppdif
      real wetlayer(0:numofpfts,2),wetratio(0:numofpfts)
      
      double precision output(500)
      real somnpp, somfractt, meanfractt
      real somC3npp,somC4npp,somC3fractt,somC4fractt,meanC3fractt,meanC4fractt

	   
!      real woodypercent,grasspercent

      integer pft,pftmaxlai,pftmaxnpp,optpft,subpft
      integer pfts(numofpfts),optdata(0:numofpfts,500)
      integer dom,wdom,grasspft,m,pos,month,wetpft,firedays
      integer biome,greendays,drymonth(0:14),subfiredays !numofpfts

      logical grass(14),present(14),flop  !numofpfts
      logical diagmode

      do pft=1,numofpfts
       if (pft.ge.8) then
        grass(pft)=.true.
       else
        grass(pft)=.false.
       end if

       grass(10)=.false.

       if (optnpp(pft).gt.0.0) then
        present(pft)=.true.
       else
        present(pft)=.false.
       end if
      end do
      
      present(12)=.true.

!----Initialize all of the variables that index an array---
      optpft=0
      subpft=0
      grasspft=0
      pftmaxnpp=0
      pftmaxlai=0
      dom=0
      wdom=0
      wetpft=0

      maxnpp=0.0
      maxlai=0.0
      temperatenpp=0.0
      maxdiffnpp=0.0
      grassnpp=0.0

!------------------------------------------------------------------------

!---- choose the dominant woody PFT on the basis of NPP:

!------------------------------------------------------------------------------
!-----Find the PFTs with the highest npp and lai-------------------------

      do pft=1,12

       if (grass(pft)) then                 !grass PFT's
        if (optnpp(pft).gt.grassnpp) then
         grassnpp=optnpp(pft) 
         grasspft=pft
        end if

       else

        if (optnpp(pft).gt.maxnpp) then     !woody PFT's
         maxnpp=optnpp(pft)
         pftmaxnpp=pft
        end if
        if (optlai(pft).gt.maxlai) then
         maxlai=optlai(pft)
         pftmaxlai=pft
        else if (optlai(pft).eq.maxlai) then
         maxlai=optlai(pftmaxnpp)
         pftmaxlai=pftmaxnpp
        end if
       end if

      end do
 
!-----Find average annual soil moisture value for all PFTs:----------


      do pft=1,numofpfts

       wetness(pft)=0.
       wetlayer(pft,1)=0.0
       wetlayer(pft,2)=0.0
       drymonth(pft)=0
       wettest(pft)=-1.0
       driest(pft)=101.0

       do m=1,12
!        mwet=optdata(pft,m+12)
        mwet=optdata(pft,m+412)
        wetness(pft)=wetness(pft)+mwet/12.
        wetlayer(pft,1)=wetlayer(pft,1)+optdata(pft,m+412)/12.  !top
        wetlayer(pft,2)=wetlayer(pft,2)+optdata(pft,m+424)/12.  !bottom
        if (mwet.gt.wettest(pft)) wettest(pft)=mwet
	if (mwet.lt.driest(pft)) then 
	 drymonth(pft)=m
	 driest(pft)=mwet
	end if
       end do

      end do
!----------------------------------------

      optpft=pftmaxnpp
      wdom=optpft

!     find the subdominant woody pft (2nd in NPP)

      subnpp=0.0
      subpft=0

      do pft=1,7
       if (pft.ne.wdom) then
        if (optnpp(pft).gt.subnpp) then
         subnpp=optnpp(pft)
         subpft=pft
        end if
       end if
      end do

!-----------------------------------------------------------

      flop=.false.

1     continue

      woodylai= optlai(wdom)
      woodnpp = optnpp(wdom)
      grasslai= optlai(grasspft)

      if (wdom.ne.0) then
       firedays=optdata(wdom,199)
       subfiredays=optdata(subpft,199)
       greendays=optdata(wdom,200)
      else 
       firedays=0
       subfiredays=0
       greendays=0
      end if

      nppdif=optnpp(wdom)-optnpp(grasspft)

      ratio=0.0

!------------------------------------------------------------

      if ((wdom.eq.3.or.wdom.eq.5).and.tmin.gt.0.0) then
       if (gdd5.gt.5000.0) then
        wdom=2
        goto 1
       end if	
      end if

!------------------------------------------------------------
!     Under certain conditions grass will be the dominant or co-dominant PFT:

      if (wdom.eq.1) then
       if (optnpp(wdom).lt.2000.0) then
        wdom=2
	subpft=1
	goto 1
       end if
      end if  	

      if (wdom.eq.2) then
       if (woodylai.lt.2.0) then
        optpft=grasspft
       else if (grasspft.eq.9.and.woodylai.lt.3.6) then
        optpft=14
       else if(greendays.lt.270.and.tcm.gt.21.0.and.tprec.lt.1700.0)then
        optpft=14
       else
        optpft=wdom
       end if
      end if

      if (wdom.eq.3) then
       if (optnpp(wdom).lt.140.0) then
        optpft=grasspft
       else if (woodylai.lt.1.0) then
        optpft=grasspft
       else if (woodylai.lt.2.0) then
        optpft=14
       else
        optpft=wdom
       end if
      end if 

      if (wdom.eq.4) then
       if (woodylai.lt.2.0) then
        optpft=grasspft
       else if (firedays.gt.210.and.nppdif.lt.0.0) then
	if (.not.flop.and.subpft.ne.0) then
	 wdom=subpft
	 subpft=4
	 flop=.true.
	 goto 1
        else
	 optpft=grasspft
	end if
       else if (woodylai.lt.3.0.or.firedays.gt.180) then
        if (nppdif.lt.0.0) then
         optpft=14
	else if (.not.flop.and.subpft.ne.0) then
	 wdom=subpft
	 subpft=4
	 flop=.true.
	 goto 1
	end if
       else
	optpft=wdom
       end if
      end if

      if (wdom.eq.5) then
       if (present(3)) then
        wdom=3
	subpft=5
	goto 1
!       else if (nppdif.lt.0.0) then
!	if (.not.flop.and.subpft.ne.0) then
!	 wdom=subpft
!	 subpft=5
!	 flop=.true.
!	 goto 1
!	end if  
       else if (optnpp(wdom).lt.140.0) then
        optpft=grasspft
       else if (woodylai.lt.1.2) then
        optpft=14
       else
        optpft=wdom
       end if
      end if 


!     add npp limits on other PFT's (tropical mountain story)

      if (wdom.eq.6) then
       if (optnpp(wdom).lt.140.0) then
        optpft=grasspft
       else if (firedays.gt.90) then
	if (.not.flop.and.subpft.ne.0) then
	 wdom=subpft
	 subpft=6
	 flop=.true.
	 goto 1
	end if  
       else
        optpft=wdom
       end if
      end if

      if (wdom.eq.7) then
       if (optnpp(wdom).lt.120.0) then
        optpft=grasspft
!       else if (optnpp(wdom).lt.120.0) then
!        optpft=14
       else if (wetness(wdom).lt.30.0.and.nppdif.lt.0.0) then
        optpft=grasspft
       else	
        optpft=wdom
       end if
      end if 		  

      if (wdom.eq.0) then 
       if (grasspft.ne.0) then
        optpft=grasspft
       else if (optnpp(13).ne.0.0) then
        optpft=13
       else
        optpft=0
       end if
      end if

      if (optpft.eq.0.and.present(10)) optpft=10

      if (optpft.eq.10) then
       if (grasspft.ne.9.and.(optnpp(grasspft).gt.optnpp(10))) then
        optpft=grasspft
       else
        optpft=10
       end if
      end if 	

      if (optpft.eq.grasspft) then
       if (optlai(grasspft).lt.1.8.and.present(10)) then
        optpft=10
       else
        optpft=grasspft
       end if
      end if

      if (optpft.eq.11) then
       if (wetness(optpft).le.25.0.and.(present(12))) then
        optpft=12
       end if
      end if

!----------------------------------------------------------------------

!----------------------------------------------------------------------

!     output some diagnostic results

      if (diagmode) then
       do pft=1,numofpfts
        if (pfts(pft).ne.0) then
         write(*,4)pft,drymonth(pft),driest(pft),wetness(pft),optdata(pft,199),optdata(pft,200)
4        format(2I5,2F6.2,2I5)
        end if
       end do

       write(*,*)' wpft  woodynpp   woodylai gpft grassnpp subpft phi'
       write(*,5)wdom,woodnpp,woodylai,grasspft,grassnpp,subpft,optdata(8,52)/100.
5      format(I5,F10.2,F10.2,I5,F10.2,I5,F8.2,F8.2)
      
      end if
!------------------------------------------------------
!     put some variables into format for output

      dom=optpft


! fractionnement isotopique Hatte
      somnpp=0
      somC3npp=0
      somC4npp=0
      somfractt=0
      somC3fractt=0
      somC4fractt=0
      meanfractt=0
      meanC3fractt=0
      meanC4fractt=0
      do pft=1,numofpfts
        somnpp = somnpp + float(optdata (pft,53))/10.
        somfractt=somfractt +float(optdata (pft,53))*float(optdata (pft,59))/100.0
        somC3npp = somC3npp + float(optdata (pft,55))/10.
        somC3fractt=somC3fractt +float(optdata (pft,55))*float(optdata (pft,60))/100.0
        somC4npp = somC4npp + float(optdata (pft,57))/10.
        somC4fractt=somC4fractt +float(optdata (pft,57))*float(optdata (pft,61))/100.0
      end do

!      if (optpft.eq.14) then
!         somnpp = somnpp + optdata (grasspft,1)
!         somfractt=somfractt +optdata (grasspft,1)*optdata (grasspft,53)/100.0
!      endif
      if (somnpp.gt.0.0)meanfractt = somfractt/somnpp
      output(54)=meanfractt
!	  write(8,*)optpft,somfractt,somnpp
      output(57)=somC4npp/somnpp*100
	  if(somC3npp.gt.0.0) meanC3fractt = somC3fractt/somC3npp
      output(55)=meanC3fractt
	  if(somC4npp.gt.0.0) meanC4fractt = somC4fractt/somC4npp
      output(56)=meanC4fractt





!-----------------------------------------------------savanne
      if (optpft.eq.14) then
       dom=wdom
       lai=(woodylai+(2.*grasslai))/3.0
       npp=(woodnpp+(2.*grassnpp))/3.0
! NEP
       do pos=137,148
        optdata(dom,pos)=(optdata(wdom,pos)+(2.*(optdata(grasspft,pos))))/3.0
       end do
! DeltaA
       do pos=80,91
        optdata(dom,pos)=(optdata(wdom,pos)+(2.*(optdata(grasspft,pos))))/3.0
       end do
! NPP
       do pos=37,48
        optdata(dom,pos)=(optdata(wdom,pos)+(2.*(optdata(grasspft,pos))))/3.0
       end do
! Rh
       do pos=113,124
        optdata(dom,pos)=(optdata(wdom,pos)+(2.*(optdata(grasspft,pos))))/3.0
       end do
! DeltaE
       optdata(dom,101)=(optdata(wdom,101)+(2.*(optdata(grasspft,101))))/3.0
! %C4 NPP
       optdata(dom,98)=(optdata(wdom,98)+(2.*(optdata(grasspft,98))))/3.0
      end if
!------------------------------------------------------

      if (optlai(dom).eq.0.0) optpft=0

!      npp=optnpp(dom)
!      lai=optlai(dom)

      npp=optnpp(wdom)
      lai=optlai(wdom)
      grasslai=optlai(grasspft)

!      npp=optnpp(grasspft)
!      lai=optlai(grasspft)

!      lai=woodylai
!      lai=grasslai

      call newassignbiome(optpft,wdom,grasspft,subpft,npp,woodnpp,grassnpp,subnpp, &
      greendays,biome,gdd0,gdd5,tcm,present,woodylai,grasslai,tmin)

!-----------------------------------------------------------------------
!       The values of all output variables, except the actual biome type
!       [output(1)] are assigned here:
 
        output(2)=lai
        output(3)=npp  
        output(4)=optlai(wdom)
        output(5)=optnpp(wdom)
        output(6)=optlai(grasspft)
        output(7)=optnpp(grasspft)

!       Annual APAR / annual PAR expressed as a percentage:
        output(8) = optdata(dom,8)

!       Respiration costs (for dom plant type, wood or grass):
        output(9) = optdata(dom,9)

!       Soil moisture for dominant pft:
        output(10)= wetness(dom)

!       Soil moisture for dominant pft:
!        output(10)= nint(wetness(wdom)*10.)

!       Predicted runoff (for dom plant type, wood or grass):
        output(11)= optdata(wdom,6)

!	annual aet  (Guiot 23/09/2019)
	      output(19)=optdata(wdom,3)

!       Number of the dominant (woody) pft:  
        output(12)=optpft
!        output(12)=wdom

!       Total annual precipitation (hopefully<9999mm):
        output(13)=tprec
        if(output(13).gt.9999) output(13)=9999

!       Total annual PAR MJ.m-2.yr-1
        output(14)=optdata(dom,7)

        output(15)=lairatio

        output(16)=nppdif

        if (lai.lt.2.0) then               !FVC
         output(17)=nint(lai/2.0*100.)
        else
         output(17)=100
        end if

        if (dom.eq.0) then 
         output(18)=0.0
        else
         output(18)=nint(pftpar(dom,6)*100.)  !root percent
        end if      
!       Store monthly fpar values in positions 25-36:
        do 12 month=1,12
         output(24+month)=optdata(dom,24+month)
 12     continue

!      Annual mean delta 13C stored in position 50-51
!! modified by C. Hatte and J. Guiot
       output(50)=optdata(dom,50)/10.0  !C3 photosynthesis
       output(51)=optdata(dom,51)/10.0  !C4 photosynthesis
       output(52)=optdata(8,52)    !phi value
       output(53)=optdata(dom,53)/10.0  !C3+C4 photosynthesis
!!!	write(8,'(a,f10.2)')'C3+C4 photosynth (avant)',output(53)
!!!	write(8,'(a,f10.2)')'meanfractt',output(54)
!!!	write(8,'(a,f10.2)')'meanC3fractt',output(55)
!!!	write(8,'(a,f10.2)')'meanC4fractt',output(56)
!!!	write(8,'(10a9)')'pft n°','misnpp','misgpp','c3npp','c3gpp','c4npp','c4gpp','isoC','isoc3','isoc4'
!!!	do pft=1,numofpfts
!!!		write(8,'(I9,9f9.1)')pft,float(optdata(pft,53:61))/10.0
!!!	enddo

!------------------------------------
!     changed to NPP for temporary output

       output(60:70)=optnpp(1:11)   !optimized NPP for all PFTs

!------------------------------------

!       output(60)=optdata(1,50)  !mean C3 photosynthesis discrimination
!       output(61)=optdata(2,50)  !for all PFTs
!       output(62)=optdata(3,50)  
!       output(63)=optdata(4,50)
!       output(64)=optdata(5,50)
!       output(65)=optdata(6,50)
!       output(66)=optdata(7,50)
!       output(67)=optdata(8,50)
!       output(68)=optdata(9,50)
!        
!       output(69)=optdata(8,51)  !mean C4 photosynthesis discrimination
!       output(70)=optdata(9,51)  !for grasses

!------------------------------------

       do pos=80,91
        output(pos)=optdata(dom,pos)  !monthly discrimination for one pft
       end do

       do pos=37,48
        output(pos)=optdata(dom,pos)  !monthly npp, one pft
       end do

       do pos=101,112
        output(pos)=optdata(dom,pos)  !monthly delta e, dominant PFT
       end do

       do pos=113,124
        output(pos)=optdata(dom,pos)  !monthly het resp, dom pft
       end do

       do pos=125,136
        output(pos)=optdata(dom,pos)  !monthly isoresp (product)
       end do

       do pos=137,148
        output(pos)=optdata(dom,pos)  !monthly net C flux (npp-resp)
       end do

       do pos=160,172
        output(pos)=optdata(dom,pos)  !monthly mean gc
       end do

       do pos=173,184
        output(pos)=optdata(dom,pos)  !monthly LAI
       end do

       output(149)=optdata(dom,149)   !annual NEP
       output(150)=optdata(dom,150)   !annual mean A/g

       output(97)=optdata(dom,97)  !Mean annual hetresp scalar
       output(98)=optdata(dom,98)  !pct of NPP that is C4
       output(99)=optdata(dom,99)  !annual het resp

       output(199)=optdata(wdom,199) !firedays
       output(200)=optdata(dom,200) !greendays

       do pos=201,241
        output(pos)=optdata(dom,pos) !ten-day lai*100
       end do

!      Monthly soil moisture, mean, top, and bottom layers *100
       
       do pos=389,400
        output(pos)=optdata(dom,pos-376)  !mean
       end do

       do pos=413,424
        output(pos)=optdata(dom,pos)      !top
       end do

       do pos=425,436
        output(pos)=optdata(dom,pos)      !bottom
       end do
       
       output(425)=wetlayer(dom,1)
       output(426)=wetlayer(dom,2)
       output(427)=wetratio(dom)*100.
       
       output(450)=optdata(dom,450)    !meanKlit
       output(451)=optdata(dom,451)    !meanKsoil

!---------------------------------------------------------------------

      return

      end

!*******************************************************************************
      subroutine newassignbiome &
     (optpft,woodpft,grasspft,subpft,optnpp,woodnpp,grassnpp,subnpp, &
      greendays,biome,gdd0,gdd5,tcm,present,woodylai,grasslai,tmin)

!     this is a new subroutine for assigning biomes in BIOME3.5
!     according to a new scheme of biomes.
!     Jed Kaplan 3/1998

      implicit none

      integer biome,greendays
      integer optpft,woodpft,grasspft,subpft

      real optnpp,woodnpp,grassnpp,subnpp,nppdif,gdd0,gdd5,tcm
      real woodylai,grasslai,tmin
      
      logical present(14)

!     list of the 28 biomes assigned, including land ice
!-------------------------------------------------------
!     1  Tropical evergreen forest
!     2  Tropical semi-deciduous forest
!     3  Tropical deciduous forest/woodland
!     4  Temperate deciduous forest
!     5  Temperate conifer forest
!     6  Warm mixed forest
!     7  Cool mixed forest
!     8  Cool conifer forest
!     9  Cold mixed forest
!     10 Evegreen taiga/montane forest
!     11 Deciduous taiga/montane forest
!     12 Tropical savanna
!     13 Tropical xerophytic shrubland
!     14 Temperate xerophytic shrubland
!     15 Temperate sclerophyll woodland
!     16 Temperate broadleaved savanna
!     17 Open conifer woodland
!     18 Boreal parkland
!     19 Tropical grassland
!     20 Temperate grassland
!     21 Desert
!     22 Steppe tundra
!     23 Shrub tundra
!     24 Dwarf shrub tundra
!     25 Prostrate shrub tundra
!     26 Cushion forb lichen moss tundra
!     27 Barren
!     28 Land ice
!-------------------------------------------------------

!-----barren-------------------------
      if (optpft.eq.0) then
       biome=27
       goto 200
      end if
!------------------------------------

!-----arctic/alpine biomes-----------
      if (optpft.eq.13) then
       biome=26    !cushion-forb tundra 
       goto 200
      end if

      if (optpft.eq.11) then
        if (gdd0.lt.200.0) then
         biome=25  !prostrate shrub tundra
         goto 200
        else if (gdd0.lt.500.0) then
         biome=24  !dwarf shrub tundra 
        else
         biome=23  !shrub tundra
        end if 
       goto 200
      else if (optpft.eq.12) then
       biome=22    !steppe-tundra
       goto 200
      end if
!------------------------------------
!-----desert-------------------------
	if (optpft.eq.10) then
		if (grasslai.gt.1.0) then
			if (tmin.ge.0.0) then
				biome=13
				goto 200
			else
				biome=14
				goto 200
			end if
		else
			biome=21
			goto 200
		end if
	else if (optnpp.le.100.0) then
		if (optpft.le.5.or.optpft.eq.9.or.optpft.eq.10) then
			biome=21
			goto 200
		else if (optpft.eq.8) then
			if (subpft.ne.6.or.subpft.ne.7) then
				biome=21
				goto 200
			end if
		end if
	end if
!------------------------------------

!-----boreal biomes------------------
	if (optpft.eq.6) then
 		if (gdd5.gt.900.0.and.tcm.gt.-19.0) then
 			if (present(4)) then
				biome=7
			else
				biome=8
			end if 
 		else
 			if (present(4)) then
				biome=9
			else
				biome=10
			end if 
		end if
		goto 200
	end if

	if (optpft.eq.7) then
		if (subpft.eq.4) then
			biome=4
			goto 200
		else if (subpft.eq.5) then  !.or.subpft.eq.6
			biome=9
			goto 200
 		else if (gdd5.gt.900.0.and.tcm.gt.-19.0) then
 			biome=9
 			goto 200
		else
			biome=11
			goto 200
		end if
	end if       
!------------------------------------

!-----temperate biomes---------------
      nppdif=optnpp-subnpp

	if (optpft.eq.8) then
		if (gdd0.ge.800.0) then
			biome=20 
 			goto 200
		else 
			biome=22
			goto 200
		end if
	end if

	if (optpft.eq.3) then
		biome=6
		goto 200
	end if

	if (optpft.eq.4) then
		if (present(6)) then
			if (tcm.lt.-15.) then
				biome=9 !cold mixed
			else
				biome=7 !cool mixed
			end if
			goto 200
		else if (present(3).or.(present(5).and.gdd5.gt.3000.0.and.tcm.gt.3.0))then
			biome=6
			goto 200
		else
			biome=4  !TeDeFo
			goto 200
		end if
 	end if

      if (optpft.eq.5) then
       if (present(3)) then 
        biome=6
        goto 200
       else if (subpft.eq.4.and.nppdif.lt.50.) then
        biome=5
        goto 200
       else if (subpft.eq.7) then
        biome=9
        goto 200
       else
        biome=5
        goto 200
       end if
      end if

!------------------------------------
!-----savanna and woodland-----------
      if (optpft.eq.14) then
       if (woodpft.le.2) then
        if (woodylai.gt.4.0) then
         biome=12  !tropical savanna
        else	
         biome=13  !tropical xero scrub 
	end if
        goto 200
       else if (woodpft.eq.3) then
        biome=15   !sclerophyll woodland
        goto 200
       else if (woodpft.eq.4) then
        biome=16   !temp brdlf savanna
        goto 200
       else if (woodpft.eq.5) then
        biome=17   !open conifer
        goto 200
       else if (woodpft.eq.7.or.woodpft.eq.6) then
        biome=18   !boreal parkland
        goto 200
       end if
      end if
!------------------------------------

!-----tropical biomes----------------
      if (optpft.le.2.or.optpft.eq.9) then

       if (optpft.eq.1) then
        biome=1
        goto 200
       end if

       if (optpft.eq.2) then
        if (greendays.gt.300) then
         biome=1
         goto 200
        else if (greendays.gt.250) then
         biome=2
         goto 200
        else
         biome=3
         goto 200
        end if
       end if

       if (optpft.eq.9) then
        biome=19
        goto 200
       end if

      end if
!------------------------------------
      
      biome=0

200   return

      end

!*************************************************************************
!      Run npp optimization model for one pft:
   
       subroutine findnpp(pfts,pft,optlai,optnpp,wst,dtemp,sun,temp,dprec,dmelt,dpet,dayl,k, &
       pftpar,optdata,dphen,co2,p,tsoil,realout,numofpfts)
       implicit none                                   
       integer numofpfts
       integer pft,i,pfts(numofpfts),optdata(0:numofpfts,500)
       integer inv(500),iterate

       real temp(12),sun(12),dayl(12),tsoil(12)
       real dprec(365),dphen(365,2)
       real dpet(365),co2,p,realout(0:numofpfts,200),realin(200)
       real optnpp,optlai,dtemp(365)
       real wst,k(12),pftpar(25,25)
       real npp,dmelt(365)
       real lowbound,range,alai(2)
	   real meanC3,meanC4,partc3,partc4
	   common /izzo/meanC3,meanC4,partc3,partc4

!      Set output values to zero:
       do i=1,500
        optdata(pft,i)=0
       end do

       optnpp = 0.
       optlai = 0.

!      If pft.ne.0 this is a dummy call of subroutine, return zero values:
       if(pfts(pft).ne.1) then
        return
       end if

!       print*,'entered findnpp',pft

!---------------------------------------------------------------------      
     
!      Calculate NPP at a range of different leaf areas by iteration:

      lowbound=0.01
      range=8.0

      do iterate=1,8

       alai(1)=lowbound+(1./4.)*range
       alai(2)=lowbound+(3./4.)*range

       call growth(npp,alai(1),wst,sun,temp,dprec,dmelt,dpet, &
      k,pftpar,pft,dayl,dtemp,inv,dphen,co2,p,tsoil,realin)

        if (npp.ge.optnpp) then
         optlai = alai(1)
         optnpp = npp
         do i=1,500
          optdata(pft,i)=inv(i)
!          realout(pft,i)=realin(i)
         end do
        end if

       call growth(npp,alai(2),wst,sun,temp,dprec,dmelt,dpet, &
      k,pftpar,pft,dayl,dtemp,inv,dphen,co2,p,tsoil,realin)

!       Find the leaf area which gives the highest NPP:
        if (npp.ge.optnpp) then
         optlai = alai(2)
         optnpp = npp
         do i=1,500
          optdata(pft,i)=inv(i)
!          realout(pft,i)=realin(i)
         end do
        end if
       range=range/2.0
       lowbound=optlai-range/2.0
       if (lowbound.le.0.0) lowbound=0.01

      end do

!---------------------------------------------------------------------      

20     return

       end
!***************************************************************************
!      Subroutine growth calculates NPP of one PFT
 
       subroutine growth(npp,maxlai,annp,sun,temp,dprec,dmelt, &
       dpet,k,pftpar,pft,dayl,dtemp,outv,dphen,co2,p,tsoil,realin)

       implicit none

       integer pft,m,month,j,days(12),phentype!,midday(12)
       integer outv(500),grass,i,c4months,greendays,day,mcount

       logical c4,c4month(12),wilt

       real dprec(365),dpet(365),annp,npp,p
       real gpp,dayl(12),pftpar(25,25),maxfvc
       real sun(12),temp(12),root,c4pot
       real k(12),age,emax
       real meangc(12),realin(200)
       real mgpp(12),gphot
       real rainscalar,wst,doptgc(365),optgc(12)
       real dtemp(365),dphen(365,2)
       real xmid,aday,dx,fmid
       real igphot,gt,meanfvc(12),dmelt(365)
       real ap,fpar,tsecs,meanwr(12,3)
       real rtbis,x1,x2,annaet,co2,ca
       real pgphot,leafresp,mlresp(12),leaftime,maxgc
       real alresp,mgmin,lresp(12),fr,runoff
       real optratio,kk(13),stemresp,maxlai,optratioa(13)
       real annualapar,annualparr,annualfpar
       real monthlyfpar(12),monthlyparr(12),monthlyapar(12)
       real backleafresp(12)

 !! modified by C. Hatte & J. Guiot
	   real isoC
       real CCratio(12),isoresp(12),phi,isoC3,isoC4,dayfvc(365)
       real C3DA(12),c4gpp(12),c4pct,annresp,annnep,C4DA(12)
       real maintresp(12),mstemresp(12),mrootresp(12),mgrowresp(12)
       real mnpp(12),meanaet(12),tsoil(12),cflux(12)
       real rlit(12),rfst(12),rslo(12),riso(12),rtot(12),riflux(12)
       real wet(365),firedays,annc4npp,c4ccratio(12),c4leafresp(12)
       real c4fpar(12),c4parr(12),c4apar(12),monthlylai(12)
       real tendaylai(40),totnpp,c4mnpp(12),nppsum,anngasum,Rmean

       real meanKlit,meanKsoil
	   Logical isc4month(12)
	   integer IsC4months
	   real C4npp(12),c3npp(12),c3gpp(12),c3ccratio(12),isnpp(12),isgpp(12)
	   real misnpp,misgpp,misc3npp,misc3gpp,misc4npp,misc4gpp   	

!       data (midday(m),m=1,12)
!     *   / 16,44,75,105,136,166,197,228,258,289,319,350 /

       data (days(month),month=1,12) &
        /  31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.  /

!      This array defines pft-specific maximum Ci/Ca ratios.
       data (optratioa(i),i=1,13) &
     /0.95,0.9, 0.8,0.8,0.9, 0.8,0.9, 0.65,0.65, 0.70, 0.90,0.75,0.80/

       data (kk(i),i=1,13) &
       / 0.7,0.7,0.6,0.6,0.5,0.5,0.4,0.4,0.4,0.3,0.5,0.3,0.6 /

!---------------------------------------------------------------------------
       ca=co2*1e-6
!---------------------------------------------------------------------------
!      Initialize day one value for soil moisture
!      COULD MAYBE DO THIS AS A FUNCTION OF AET/PET! 
       rainscalar=1000.
       wst = annp / rainscalar
       if(wst.ge.1.)  wst=1.
!---------------------------------------------------------------------------
!      Assign pft specific parameters for photosynthesis model 

       phentype  = nint(pftpar(pft,1))
       mgmin     = pftpar(pft,2)
       root      = pftpar(pft,6)
       age       = pftpar(pft,7)
       c4pot     = pftpar(pft,11) 
       grass     = nint(pftpar(pft,10))
       emax      = pftpar(pft,3)
!----------------------------------------------------------------------------
!      Calculate the maxfvc from the maxlai and (fixed) k value
!       kk=0.5
       maxfvc=1.-exp(-kk(pft)*maxlai)

!------------------------------------------------------------------------    
!      Set the value of optratio depending on whether c4 plant or not.

       if (pft.eq.8.or.pft.eq.9.or.pft.eq.10) then
		c4=.true.
		
       else
        c4=.false.
       endif

 100   continue

       if (c4) then
        optratio=0.4
       else
        optratio=optratioa(pft)  !assigns the pft-specific ratio
       endif

!----------------------------------------------------------------------------
!      Calculate monthly values for the optimum non-water-stressed gc (optgc)
       
      maxgc=0.
      do 120 month=1,12
      m=month
      tsecs=3600.*dayl(m)

!------------------------------------------------------------------------
!      First find the gc value resulting from the max ci/ca ratio
    
!      Find mid-monthly-day daily tstressed photosynthesis value
       fpar = 1.-exp(-kk(pft)*maxlai)


       if (c4) then
        call c4photo(optratio,sun(m),dayl(m),temp(m), &
       age,lresp(m),pgphot,aday,fpar,p,ca,pft)
       else
        call photosynthesis(optratio,sun(m),dayl(m),temp(m), &
       age,lresp(m),pgphot,aday,fpar,p,ca,pft)
       end if

!     changed to include possible aday zero
!     Calculate gt using the physical eqn from aday and ci/ca ratio
      if (tsecs.gt.0.and.aday.gt.0.0) then
       gt = mgmin + ( (1.6*aday) / (ca*(1.-optratio)) ) / tsecs
      else
       gt=0.0
      end if
   
!      This gives us the final non-water-stressed gc value
       optgc(m) = gt

!      Store output values:
       if(maxgc.le.optgc(m)) maxgc=optgc(m)

 120   continue
!-------------------------------------------------------------------------
!      Calculate water balance and phenology for this pft/s and fvc/s
!      Subroutine hydrology returns monthly mean gc & summed fvc value
   
!      Linearly interpolate the mid-month optgc & ga values to daily values
       call daily(optgc,doptgc)
    
       call hydrology &
      (dprec,dmelt,dpet,root,k,maxfvc,pft,phentype,wst, &
       doptgc,meangc,meanfvc,meanwr,meanaet,annaet,mgmin,dphen,dtemp, &
       grass,runoff,wet,greendays,dayfvc,emax,wilt,pftpar)      

!-------------------------------------------------------------------------
!     Now use the monthly values of fvc & meangc to calculate net & gross
!     photosynthesis for an "average" day in the month and multiply by the
!     number of days in the month to get total monthly photosynthesis.

      alresp   = 0.
      gpp      = 0.
      leaftime = 0.
      annualparr=0.
      annualapar=0.

!-------------------------------------------------------------------------      
      do month=1,12                                   !begin monthly loop here
      m=month

!      c4gpp(m)=0.0

!     If meangc is zero then photosynthesis must also be zero
      if(meangc(m).eq.0.)then
       gphot=0.
       rtbis=0.
       leafresp=lresp(m)*(meanfvc(m)/maxfvc)

      else   
!     Iterate to a solution for gphot given this meangc value!
!....................................................................
!      This is a tailored implementation of the bisection method
!      with a fixed 8 bisections and assuming root is bracketed and 
!      that f(x1)<0 and f(x2)>0   

       if(C4)then 
		x1=0.02 
	   else 
	    x1=0.5 
	   endif
       x2=optratio+0.05
       rtbis=x1
       dx=x2-x1
       do 30 j=1,10
        dx=dx*0.5
        xmid=rtbis+dx
!...................................................
!      Evaluate fmid=Anetdt-Ap at the point ci/ca = xmid     
              
       fpar = meanfvc(m)


       if (c4) then 
	    call c4photo(xmid,sun(m),dayl(m),temp(m),age,leafresp,igphot,aday,fpar,p,ca,pft)
       else
        call photosynthesis(xmid,sun(m),dayl(m),temp(m),age,leafresp,igphot,aday,fpar,p,ca,pft)
       end if
  
       gt = 3600.*dayl(m)*meangc(m)

       if (gt.eq.0.0) then
        ap=0.0
       else
        ap = mgmin + (gt/1.6)*(ca*(1.-xmid))
       end if

       fmid = aday - ap
!....................................................
!      If fmid is closer to the root store new values
       if(fmid.le.0.)then
        rtbis=xmid    
        gphot = igphot
       endif
 30    continue
!....................................................................        
       endif
      
!     We already include the albedo in the calculation of the net 
!     short wave radiation so apar here really is the absorbed PAR:
!     Idea here is that annualfpar should be the total amount of PAR
!     ansorbed during the year divided by the total amount of PAR per
!     unit area.
!     Get PAR in units of MJ.m-2.month-1 (hence *1e-6)

      monthlyfpar(m) = meanfvc(m)
      monthlyparr(m) = sun(m)*days(m)*1e-6
      monthlyapar(m) = monthlyparr(m)*monthlyfpar(m)
      annualapar  = annualapar + monthlyapar(m) 
      annualparr  = annualparr + monthlyparr(m)

!     Monthly gross photosynthesis (=numdays*average-daily-photosynthesis)
      mgpp(m) = days(m)*gphot
      gpp     = gpp + mgpp(m)

!     Calculate monthly leaf respiration (=numdays*average-daily-leafresp)
      mlresp(m) = days(m)*leafresp
      alresp    = alresp + mlresp(m)    

!     store monthly values of Ca/Cst and leaf resp.

      CCratio(m) = rtbis
      isoresp(m) = mlresp(m)

      if (c4) then
       c4gpp(m)=mgpp(m)            !store monthly c4 values here
       c4fpar(m)=monthlyfpar(m)
       c4parr(m)=monthlyparr(m)
       c4apar(m)=monthlyapar(m)
       c4ccratio(m)=CCratio(m)
       c4leafresp(m)=isoresp(m)
      end if

      end do                                              !the loop ends

!--------------------------------------------------------------------------
!     Calculate monthly LAI

      do m=1,12
       monthlylai(m)=(log(amax1(1-monthlyfpar(m),1e-8)))/(-1.0*kk(pft))
      end do

!     And ten-day lai
      i=1
      do day=1,365,10
       tendaylai(i)=(log(amax1(1-dayfvc(day),1e-8)))/(-1.*kk(pft))
       i=i+1
      end do

!---------------------------------
!     Calculate annual FPAR (%) from annual totals of APAR and PAR

      if (annualapar.eq.0.) then
       annualfpar=0.0
      else
       annualfpar=100.*annualapar/annualparr
      end if

!---------------------------------------

!     Calculate annual respiration costs to find annual npp:
      call respiration &
     (npp,gpp,alresp,temp,grass,maxlai,stemresp,fr, &
     mstemresp,mrootresp,pft,mlresp,monthlyfpar,backleafresp)

!     If the plant wilted during the year, it is dead for this LAI.

      if (wilt) npp=-9999.0

!---------------------------------------
!     calculate monthly NPP

      nppsum=0.0

      do m=1,12
       mnpp(m)=0.0
      end do

      do m=1,11
       maintresp(m)=mlresp(m)+backleafresp(m)+mstemresp(m)+mrootresp(m)
       mgrowresp(m) = (0.02*(mgpp(m+1)-maintresp(m+1)))
        if (mgrowresp(m).lt.0.0) mgrowresp(m)=0.0
       mnpp(m) = mgpp(m)-(maintresp(m)+mgrowresp(m))
      end do
      
      maintresp(12)=mlresp(12)+backleafresp(12)+mstemresp(12)+mrootresp(12)
      mgrowresp(12) = (0.02*(mgpp(1)-maintresp(1)))
       if (mgrowresp(12).lt.0.0) mgrowresp(12)=0.0
      mnpp(12) = mgpp(12)-(maintresp(12)+mgrowresp(12))

      do m=1,12
       if (c4) c4mnpp(m) = mnpp(m)
       nppsum=nppsum+mnpp(m)
      end do

! C. Hatte
! à partir d'ici je vais faire un stockage des valeurs qui me serviront pour les calculs 
! isotopiques. Je ne peux pas dans le calcul du fractionnement ne considérer les choses 
! que si elles sont majoritaires. il me faut tout a condition bien entendu que le cycle de vie
! soit respecté. Je ne préfère pas transformer le reste du programme alors il est certain 
! qu'il va y avoir des doublons avec ce qui est déjà stocké par ailleurs. On verra si on veut 
! l'intégrer par la suite. .. masi as-t-on le droit de triturer comme ça le boulot des autres?

	do m=1,12
		if (c4)then
		 	C4npp(m) = mnpp(m)
			C4ccratio(m) = ccratio(m)
			C4gpp(m) = mgpp(m)
		else
			C3npp(m) = mnpp(m)
			C3ccratio(m) = ccratio(m)
			C3gpp(m) = mgpp(m)
		endif
	end do


!----------------------------------------------------------------------

!     If this is a temperate grass or desert shrub pft, compare both
!     C3 and C4 NPP and choose the more productive one on a monthly basis. 
!     However the C4 advantage period must be for at least two months 
!     (ie. long enough to complete a life cycle).

      if (pft.eq.8.or.pft.eq.10) then
       if (c4) then
        c4=.false.
        goto 100       !compute everything again, with C3 pathway 
       end if
      end if
      c4months=0
      annc4npp=0.0
      do m=1,12
       if (pft.eq.9) then
        c4month(m)=.true.
        IsC4month(m)=.true.        !C. Hatté
       else
        c4month(m)=.false.
       end if
      end do

      if (pft.eq.8.or.pft.eq.10) then
		 IsC4months=0
         do m=1,12
           if (c4mnpp(m).gt.mnpp(m))c4months=c4months+1 
! C. Hatte
		   if (c4mnpp(m).gt.0.and.mnpp(m).gt.0.0) IsC4months=IsC4months+1
!
         end do

         if (c4months.ge.3) then
           do m=1,12
             if (c4mnpp(m).gt.mnpp(m)) then
               c4month(m)=.true.
             else
               c4month(m)=.false.
             end if
           end do
! c hatte
         else
           do m=1,12
             c4month(m)=.false.
           enddo
         end if
! C. Hatte
         if (IsC4months.ge.3) then
          do m=1,12
            if (C4mnpp(m).gt.0.and.mnpp(m).gt.0.0) then
              IsC4month(m)=.true.
            else
              IsC4month(m)=.false.
            end if
          enddo
		endif
!
      end if

      totnpp=0.0

      do m=1,12
       if (c4month(m)) then

        mnpp(m)=c4mnpp(m)
        annc4npp=annc4npp+c4mnpp(m)
        totnpp=totnpp+mnpp(m)

        monthlyfpar(m)=c4fpar(m)
        monthlyparr(m)=c4parr(m)
        monthlyapar(m)=c4apar(m)
        CCratio(m)=c4ccratio(m)
        isoresp(m)=c4leafresp(m)
       else
        totnpp=totnpp+mnpp(m)
       end if
      end do

      if (c4months.ge.2) nppsum=totnpp

!     Calculate % of annual npp that is C4
      c4pct=annc4npp/npp

!---------------------------------------------------------------------------

      if (npp.ne.nppsum) npp=nppsum

      if (npp.le.0.0) then
       return
      else
       continue
      end if

!---------------------------------------------------------------------------
      if (gpp.gt.0.0) then

!     calculate the phi term that is used in the C4 13C fractionation
!     routines
      if (pft.ge.8) then
       call calcphi(mgpp,phi)
      end if
!     calculate carbon isotope fractionation in plants
!! C. Hatte
		misnpp=0 ; misgpp=0
		misc3npp=0 ;misc3gpp=0; misc4npp=0; misc4gpp=0   
      do m=1,12
	    if (IsC4month(m)) then
		  Isnpp(m) = C4npp(m) + C3npp(m)
		  Isgpp(m) = C4gpp(m) + C3gpp(m)
	    else
		  Isnpp(m) = C3npp(m)
		  Isgpp(m) = C3gpp(m)
	    endif
        mIsnpp = mIsnpp + Isnpp(m)
        mIsgpp = mIsgpp + Isgpp(m)
		misc3npp=misc3npp+c3npp(m)
		misc3gpp=misc3gpp+c3gpp(m)
		misc4npp=misc4npp+c4npp(m)
		misc4gpp=misc4gpp+c4gpp(m)

      end do
      call isotope(C3CCratio,C4CCratio,ca,temp,isoresp,IsC4month,C3npp,C4npp,phi,isoC3,isoC4,isoC,C3DA,C4DA,mIsnpp)
!      call isotope(CCratio,ca,temp,isoresp,c4month,mgpp,phi,isoC3,isoC4,isoC,C3DA,C4DA,gpp)
 !! modified by C. Hatte & J. Guiot

      end if

!----------------------------------------------------------------------
!     Call subroutine to calculate heterotrophic respiration
      call hetresp  &
     (pft,npp,temp,tsoil,meanaet,meanwr,rlit,rfst,rslo,rtot,isoC3,riso, &
     riflux,Rmean,meanKlit,meanKsoil)

      annresp=0.0   !zero the annual hetresp calculation
      do m=1,12    !sum up monthlies to get an ann hetresp = npp
       annresp=annresp+rtot(m)
      end do
!      write(*,*)annresp

      annnep=0.0
!     calculate monthly ecosystem carbon flux NPP-Hetresp
      do m=1,12
       cflux(m)=mnpp(m)-rtot(m)
       annnep=annnep+cflux(m)
      end do

!-----------------------------------------------
!     Call fire subroutine
      call fire (wet,pft,maxlai,npp,firedays)
      
!--------------------------------------------------------------------------
!     Outputs section:
     
!     Output the monthly soil moisture value for this pft and lai:   
      do m=1,12

       outv(12+m)=nint(100.*meanwr(m,1))
       outv(412+m)=nint(100.*meanwr(m,2))
       outv(424+m)=nint(100.*meanwr(m,3))

       outv(24+m)=nint(100.*monthlyfpar(m))
      end do

!     Record values of output variables:
      outv(1) = nint(npp)                
      outv(3) = nint(annaet)
      outv(4) = nint(maxgc)
      outv(5) = nint(stemresp)
      outv(6) = nint(runoff)
      outv(7) = nint(annualparr)
      outv(8) = nint(annualfpar)
      outv(9) = nint(fr)
 !! modified by C. Hatte & J. Guiot
		outv(53)=nint(mIsnpp*10.)
		outv(54)=nint(mIsgpp*10.)
		outv(55)=nint(mIsC3npp*10.)
		outv(56)=nint(mIsC3gpp*10.)
		outv(57)=nint(mIsC4npp*10.)
		outv(58)=nint(mIsC4gpp*10.)
		outv(59)=nint(isoC*10.)
		outv(60)= nint(isoC3*10.)
		outv(61)= nint(isoC4*10.)
      outv(52)= nint(phi*100.)
      outv(97)= nint(Rmean*100.)
      outv(98)= nint(c4pct*100.)
      outv(99)= nint(annresp*10.)
      anngasum=0.0
      mcount=0

      do m=1,12
       realin(36+m)=mnpp(m)
       outv(79+m)=nint(C3DA(m)*100.)    !includes c3 and c4
       outv(36+m)=nint(mnpp(m)*10.)
!       outv(100+m)=nint(mlresp(m)*10.)
       outv(100+m)=nint(riso(m)*10.)
       outv(112+m)=nint(rtot(m)*10.)
       outv(124+m)=nint(riflux(m)*10.)
       outv(136+m)=nint(cflux(m)*10.)
       outv(160+m)=nint(meangc(m))
       outv(172+m)=nint(monthlylai(m)*100.)
       
       if (meangc(m).ne.0) then
        mcount=mcount+1
        anngasum=anngasum+(mgpp(m)/meangc(m))  !new line for A/gc
       end if
       
      end do
      
      outv(150)=nint((anngasum/mcount)*100.)            !new line for A/gc

      outv(149)=nint(annnep*10.)

      outv(199)=nint(firedays)
      outv(200)=greendays

      do i=1,40
       outv(200+i)=nint(tendaylai(i)*100.)
      end do
      
      outv(450)=nint(meanKlit*100.)
      outv(451)=nint(meanKsoil*100.)
      
!------------------------------------------------------------------------
      return

      end

!***************************************************************************
!      C3 Photosynthesis subroutine:
   
       subroutine photosynthesis(ratio,dsun,daytime,temp,age,leafresp,grossphot,aday,fpar,p,ca,pft)

       implicit none      

       integer pft,n

       real dsun,temp,daytime
       real c1,c2,s,teta,drespc3,qeffc3
       real ko,kc,ts,abs1,tao,vmax
       real ko25,kc25,tao25,koq10,kcq10,taoq10
       real p,pi,o2,jtoe,cmass,kk
       real tstress,ca,t0(13)
       real z,aday,grossphot,fpar
       real oc,ratio,age,leafresp,tune,optratio
       real je,jc,wif,adaygc,twigloss
       real slo2,mfo2
       real leafcost,tcurve(13),mintemp

       parameter(qeffc3=0.08)
       parameter(drespc3=0.015)
       parameter(abs1=1.,teta=0.7)
       parameter(slo2=20.9*1e3,jtoe=2.3*1e-6,optratio=0.95)  !figure this out
       parameter(ko25=30.*1e3,kc25=30.,tao25=2600.,cmass=12.)
       parameter(kcq10=2.1,koq10=1.2,taoq10=0.57,twigloss=1.)
                         
!      t0 defines the minimum mean monthly temperature
!      at which photosynthesis takes place

       data (t0(n),n=1,13) &
        /10.,10., 5.,4.,3., 0.,0., 4.5,10., 5., -7.,-7.,-12./

       data (tcurve(n),n=1,13) &
        /1.0,1.0, 1.0,1.0,0.9, 0.8,0.8, 1.0,1.0, 1.0, 0.6,0.6,0.5 /

!      Twigloss assumes some % of FPAR is lost to absorption by stems
!      A tune of 50% accounts for all losses!
       parameter(tune=1.0)

!      The leafcost parameter is related to expected leaf longevity

       leafcost =(age/12.)**0.25

!      Need to change the partial pressure of o2 also:

       mfo2=slo2/1e5
       o2=p*mfo2 

!      If daytime<=1 hour set to one hour to avoid /0.
       if(daytime.le.4.) daytime=4.

!---------------------------------------------------------------------------
!      tstress limits quantum efficiency as a function of temperature
!      Because of temperature optimization (and maybe nitrogen limitations)
!      this value is PFT specific.

        mintemp=t0(pft)
        if (temp.gt.mintemp+0.1) then   !the +0.1 is here to prevent underflows
         tstress=tcurve(pft)*exp(-10.0/(temp-mintemp))
        else
         tstress=0.0
        endif

!--------------------------------------------------------------------
!      work out temperature adjusted values of parameters
       ko  =    ko25*koq10**((temp-25.)/10.)
       kc  =    kc25*kcq10**((temp-25.)/10.)
       tao =  tao25*taoq10**((temp-25.)/10.)

!      Set non-co2-dependent parameters
       s  = drespc3*(24./daytime)
       ts = o2 / (2.*tao)
       kk = kc*(1. + (o2/ko)) 
       z  = cmass*jtoe*dsun*fpar*twigloss*tune

!--------------------------------------------------------------------
!      First calculate the vm value based on a ratio=0.95

       pi = optratio*ca*p
       c1 = tstress*qeffc3*( (pi-ts)/(pi+2.*ts) )       
       c2 = (pi - ts) / (pi + kk) 
       oc = ( (s-teta*s)/(c2-teta*s) )**0.5      

!      Estimate the optimal value of Vm at ratio=0.95 g(C).m-2.day-1

       if (z.eq.0.0) then
        vmax=0.0
       else
        vmax = (z/drespc3)*(c1/c2)*((2.*teta-1.)*s-(2.*teta*s-c2)*oc)
       end if

!.........................................................
!      Now use this vm value to calculate actual photosynthesis

!      Calculate inter-cellular co2 partial pressure
       pi = ratio*ca*p 

!      If pi is less than the compensation point then grossphot=0.

       if(pi.le.ts)then
        grossphot=0.0
       else
!      Otherwise calculate grossphot
        c1 = tstress*qeffc3*( (pi-ts)/(pi+2.*ts) )       
        c2 = (pi - ts) / (pi + kk)

        if (z.eq.0.0) then
         je=0.0
        else
         je=c1*z    / daytime
        end if

        if (vmax.eq.0.0) then
         jc=0.0
        else
         jc=c2*vmax / 24.
        end if

        wif = daytime/(2.*teta)

        if (je.eq.0.0.and.jc.eq.0.0) then
         grossphot=0.0
        else
         grossphot=wif*(je+jc-((je+jc)**2.-4.*teta*je*jc)**0.5)
        end if
       endif

!      Calculate net-daytime-photosynthesis and daily dark respiration
       adaygc = grossphot - (daytime/24.)*drespc3*vmax

       leafresp = drespc3*vmax*leafcost
       if (leafresp.lt.0.0) leafresp=0.0

!      Change Aday from gC.m-2.day-1 to umol.day-1
       if (adaygc.eq.0.0) then
        aday=0.0
       else
        aday = (adaygc/cmass)*(8.314*(temp+273.3)/p)*1000.
       end if

       return

       end

!**************************************************************************
!     C4 photosynthesis subroutine:
 
      subroutine c4photo(ratio,dsun,daytime,temp,age,leafresp,grossphot,aday,fpar,p,ca,pft)

       implicit none      

       integer pft

       real dsun,temp,daytime
       real c1,c2,s,teta
       real ts,abs1
       real ko25,kc25,tao25,koq10,kcq10,taoq10,tao
       real p,pi,jtoe,cmass
       real tstress,ca,t0(10)
       real z,aday,grossphot,fpar
       real oc,ratio,age,leafresp,optratio
       real je,jc,wif,twigloss
       real qeffc4,adaygcc4,grossphotc4,drespc4
       real adayc4,leafrespc4,vmaxc4
       real damage,leafcost,mintemp,maxtemp,tune
       real o2,mfo2,slo2

       parameter(drespc4=0.03)
       parameter(abs1=1.,teta=0.7)
       parameter(slo2=20.9*1e3, jtoe=2.3*1e-6, optratio=0.95)
       parameter(ko25=30.*1e3,  kc25=30.,tao25=2600., cmass=12.)
       parameter(kcq10=2.1, koq10=1.2, taoq10=0.57, twigloss=1.)

!      This leafcost parameter is a way of making the respiration costs
!      of evergreen leaves largger than deciduous leaves.

       leafcost=(age/12.0)**0.25

       t0(8)  = 10.0
       t0(9)  = 10.0
       t0(10) = 10.0

!      C4 quantum efficiency is a function of actual photosynthetic pathway
!      (NAD-me or NADP-me) and plant functional type (dicot or monocot).
!      Here we use a mean between the pathways for both woody and grass types.

       if (pft.eq.8.or.pft.eq.9) then
        qeffc4=0.0633
        tune=1.0
       else if (pft.eq.10) then
        qeffc4=0.0565
        tune=0.75     !for widely spread leaves and stems
       end if
!---------------------------------------------------------------------------
!      tstress limits quantum efficiency as a function of temperature
!      Because of temperature optimization (and maybe nitrogen limitations)
!      this value is PFT specific.

        mintemp=t0(pft)
        maxtemp=55.0
        if (temp.gt.mintemp+0.1.and.temp.lt.maxtemp) then
         tstress=exp(-10.0/(temp-mintemp))
        else
         tstress=0.0
        endif

        if (tstress.gt.1.0) tstress=1.0

!--------------------------------------------------------------------
!      work out temperature adjusted values of parameters
       tao =  tao25*taoq10**((temp-25.)/10.)

!      Set non-co2-dependent parameters
       ts = o2 / (2.*tao)
       z  = cmass*jtoe*dsun*fpar*twigloss*tune

!      Need to change the partial pressure of o2 also:
       mfo2=slo2/1e5
       o2=p*mfo2 

!--------------------------------------------------------------------
       pi = optratio*ca*p
       s  = drespc4*(24./daytime)
       c1 = qeffc4*tstress
       c2 = 1.
       oc = ( (s-teta*s)/(c2-teta*s) )**0.5

!      Estimate the optimal value of Vm at ratio=0.8 g(C).m-2.day-1
       if (z.eq.0.0) then
        vmaxc4 = 0.0
       else
        vmaxc4 = (z/drespc4)*(c1/c2)*((2.*teta-1.)*s-(2.*teta*s-c2)*oc)
       end if
!.........................................................
!      Now use this vm value to calculate actual photosynthesis

!      If pi is less than the compensation point then grossphot=0.
       if(pi.le.ts)then
        grossphotc4=0.0
       else
!      Otherwise calculate grossphot
        c1 = qeffc4*tstress
        c2 = 1.

        if (z.eq.0.0) then
         je=0.0
        else
         je=c1*z / daytime
        end if
     
        if (vmaxc4.eq.0.0) then
         jc=0.0
        else
         jc=c2*vmaxc4 / 24.
        end if

!       Damage gives the limitation of c4 photosynthesis by pi
        if(ratio.lt.0.4)then
         damage = ratio/0.4
        else
         damage = 1.
        endif
        wif = damage*daytime/(2.*teta)
  
        if (je.eq.0.0.and.jc.eq.0.0) then
         grossphotc4=0.0
        else
         grossphotc4=wif*(je+jc-((je+jc)**2.-4.*teta*je*jc)**0.5)
        end if

       endif

!      Calculate net-daytime-photosynthesis and daily dark respiration
       adaygcc4 = grossphotc4 - (daytime/24.)*drespc4*vmaxc4
       leafrespc4 = drespc4*vmaxc4*leafcost

!      Change Aday from gC.m-2.day-1 to mm.day-1
       if (grossphotc4.eq.0..and.vmaxc4.eq.0.)then
        adayc4=0.0
       else
        adayc4 = (adaygcc4/cmass)*(8.314*(temp+273.3)/p)*1000.
       end if

       aday=adayc4
       leafresp=leafrespc4
       grossphot=grossphotc4

       return
       end

!**************************************************************************
!      Subroutine to calculate annual respiration costs
 
       subroutine respiration(npp,gpp,alresp,temp,grass,lai, &
        stemresp,percentcost,mstemresp,mrootresp,pft,mlresp, &
        fpar,backleafresp)

       implicit none

       integer m,grass,pft

       real npp,gpp,alresp,temp(12),lai,mlresp(12)
       real growthresp,stemresp,leafresp,finerootresp
       real minallocation,ln,percentcost!,days(12)
       real litterfall,stemcarbon,y,p1,e0,tref,m10
       real mstemresp(12),mrootresp(12),t0,respfact(13)
       real allocfact(13)
       real fpar(12),backleafresp(12),leafmaint

!       data (days(m),m=1,12) 
!     *   /  31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.  /

!      Grass defines if there is sapwood respiration

!      Ln  = Leaf litterfall per unit Leaf area index
!      m10 = Maintenance respiration of sapwood at 10oC, gC.month-1.KgC-1
!      k   = Extinction coefficient
!      y   = Efficiency with which carbon is turned into biomass
!      p1  = Fine root respiration to litterfall ratio

!      stemcarbon = sapwood mass as KgC.m-2(leaf area).m-2(ground area)
       parameter(Ln=50.,y=0.8,m10=1.6,p1=0.25,stemcarbon=0.5) 

!      e0, t0 and tref are parameters from Lloyd & Taylor 1995
       parameter(e0=308.56, tref=10.0,t0=46.02)

!      t0 can be used as a pft specific parameter to modify the shape of
!      the curve describing the temperature dependence of respiration

       data (respfact(m),m=1,13) &
       /0.8,0.8, 1.4,1.6,0.8, 4.0,4.0, 1.6,0.8, 1.4, 4.0,4.0,4.0/

       data (allocfact(m),m=1,13) &
       /1.0,1.0, 1.2,1.2,1.2, 1.2,1.2, 1.0,1.0, 1.0, 1.0,1.0,1.5/

!------------------------------------------------------------------------
!      Calculate leafmass (gC.m-2) and litterfall (gC.m-2.year-1)
       litterfall=lai*Ln*allocfact(pft)

!      Calculate stem maintenace respiration costs in gC.year-1 

       stemresp = 0.0

       do m=1,12
        if (temp(m).le.-46.02) then
         mstemresp(m)=0.0
        else
         mstemresp(m) = lai*stemcarbon*respfact(pft)*  &
         exp( e0*(1./(tref+t0) - 1./(temp(m)+t0)) )
        end if
        stemresp = mstemresp(m) + stemresp
       end do

!      Calculate belowground maintenance respiration costs in gC.year-1     

       leafmaint=0.0
       finerootresp = p1*litterfall
       do m=1,12
        mrootresp(m)=(mstemresp(m)/stemresp)*finerootresp
        backleafresp(m)=mrootresp(m)*fpar(m)*4.0
        leafmaint=backleafresp(m)+leafmaint
       end do
   
!      Assume leaf respiration costs supplied are in units of gC.year-1
       leafresp = alresp+leafmaint

!      if this is a grass clear out the stem respiration

       if (grass.eq.2) then
        stemresp=0.0
        do m=1,12
         mstemresp(m)=0.0
        end do
       end if

!      20% of whats left goes to construction respiration, gC.year-1
       growthresp = (1.-y)*(gpp-stemresp-leafresp-finerootresp)

!      Finally calculate the resulting annual NPP, gC.year-1
       npp = gpp -stemresp-leafresp-finerootresp-growthresp     

!-------------------------------------------------------------------------

!      Now calculate the minumum allocation requirement. If NPP<litterfall
!      this lai is not sustainable so set NPP=-9999. for output purposes
       minallocation = 1.*litterfall

       if(npp.lt.minallocation) then
        npp=-9999.0
       end if

!-------------------------------------------------------------------------

!      Find respiration costs as a percentage of GPP
       if(gpp.gt.0.and.npp.ne.-9999.)then
        percentcost = 100.*(gpp-npp)/gpp
       else
        percentcost = 0.
       end if   

       return
       end     

!*************************************************************************

!      Subroutine hydrology
!      calculates the actual values of gc and soil moisture

       subroutine hydrology &
       (dprec,dmelt,deq,root,k,maxfvc,pft &
       ,phentype,wst,gcopt,meangc,meanfvc,meanwr,meanaet,annaet,mgmin &
       ,dphen,dtemp,grass,sumoff,wet,greendays,dayfvc,emax,wilt,pftpar)

       implicit none

       integer month,twice,d,days(12),dayofmonth,phentype,grass
       integer greendays,pft

       real w(2),root,fvc,k(12),aet,r1(2),emax
       real perc,runnoff,wr,meanfvc(12),meangc(12)
       real dprec(365),deq(365),dphen(365,2),gcopt(365)
       real maxfvc,wst,gc,dtemp(365),wet(365),dayfvc(365)
       real annaet,drainage,meanwr(12,3),mgmin,meanaet(12)
       real dmelt(365),supply,alfa,alfam,gm,gsurf,gmin
       real onnw,offw,sumoff,demand,waste,wetphytomass,a

       real pftpar(25,25)

       logical wilt

       parameter(alfam=1.4,gm=5.)
!       parameter(onnw=0.4,offw=0.3)

       onnw=pftpar(pft,4)
       offw=pftpar(pft,4)

       data days / 31,28,31,30,31,30,31,31,30,31,30,31  /
           
!      Initialize soil moisture stores for day one of the "spin-up" year 
       w(1) = wst
       w(2) = wst
!----------------------------------------------------------------------------
!      Run the hydrology and phenology models!
     
!      Run for one "spin-up" year & then use output from the 2nd year 
       do 100 twice = 1,2
  
       d=0   
       annaet=0.
       sumoff=0.
       greendays=0
       wilt=.false.

!      Do the daily calculations to find the monthly output values
       do 110 month = 1,12

        meanfvc(month) = 0.
        meangc(month)  = 0.    
        meanwr(month,1)= 0.
        meanwr(month,2)= 0.
        meanwr(month,3)= 0.
	    meanaet(month) = 0.0    

        do 120 dayofmonth = 1,days(month)
          d = d+1

!      Calculate effective soil moisture in rooting zone 
          wr =  root*w(1) + (1.-root)*w(2)

!      deq is the daily PET
!      maxfvc is foliar vegetation cover (related to LAI).
!      phentype is phenological type (1 evergreen, 2 summergreen, 3 watergreen)
!      offw is soil moisture threshold for leaf drop


!      Calculate vegetation phenology for today     

          if(phentype.eq.1)then                   !evergreen    
            fvc = maxfvc                                        

          else if(phentype.eq.2)then              !cold deciduous
            fvc = maxfvc*dphen(d,grass)  

          else if(grass.eq.2)then                 !cold deciduous
            fvc = maxfvc*dphen(d,grass)  

!------new code-02.05.99----
            if(fvc.gt.0.01.and.wr.gt.offw)then  !drought deciduous
              fvc = fvc
	        else if(fvc.lt.0.01.and.wr.gt.onnw)then
	          fvc = fvc
	        else
	          fvc = 0.0
	        end if
!------end new code---------

       elseif(fvc.gt.0.01.and.wr.gt.offw)then  !drought deciduous
        fvc = maxfvc
       elseif(fvc.lt.0.01.and.wr.gt.onnw)then    
        fvc = maxfvc
       else
        fvc = 0.0
       endif

       if (fvc.gt.0.0) greendays=greendays+1

!------------------------------------------------------------------------
       if(dtemp(d).le.-10.0) then
        gc=0.0
        aet=0.0
        perc=0.0
       else    !begin a loop here to calculate under non-freezing conditions

!---------------------------
!      If the fvc (ie. phenology) indicates vegetation bare of leaves,
!      there can still be some evaporation and loss of water from
!      the soil and stems (25% after Larcher 1995).

       if (fvc.eq.0.0) then
        aet=0.25*deq(d)
       else
!---------------------------

!      Calculate the optimal conductance for today
       gmin = mgmin*fvc
       gc   = gcopt(d)*(fvc/maxfvc)   
       gsurf = gc + gmin 
      
!      Calculate aet from gc & Eq (=deq) using eqn from Monteith 1995
       if (gsurf.gt.0.) then           
        alfa  = alfam*(1.-exp(-gsurf/gm) )  
        aet   = alfa*deq(d)
       else
        alfa  = 0.
        aet   = 0.
       end if

       end if                 !end phenology sensitive loop


!      Not all soil water, or precip. is effectively used, and some water
!      goes into the new phytomass

       wetphytomass = 0.01*aet
       waste = 0.01*aet

       demand = aet + wetphytomass + waste   

!      Calculate daily supply function in mm/day
       supply = emax*wr

!      If this optimal value of aet is greater than the supply function
!      Then limit aet to the supply function rate.

!      If the gmin cannot be satisfied then the plant wilts (dies) and this
!      LAI is not sustainable.

       if(demand.gt.supply)then
        a = (1. - supply/(deq(d)*alfam))
        if (a.lt.0.0) a=0.0000001
        gsurf    =  -gm*log(a)
        aet      =  supply
        gc = gsurf - gmin
!       Constrain gc value to zero!
        if (gc.le.0.0) then
         gc=0.0
         wilt=.true.
        end if
       endif
      
!      Calculate daily percolation from layer 1 to 2
       if(w(1).gt.1e-4)then
		perc = k(1)*w(1)**k(2)
	   else
	    perc=0.0
	   endif

       end if                 !the cold temp sensitive loop ends here
!-----------------------------------------------------------------------  
       
!      Exrati give rates of extraction from upper and lower soil layers

       if (wr.gt.1e-8) then   
        r1(1) =      root  * (w(1)/wr)
        r1(2) =  (1.-root) * (w(2)/wr)       
       else
        r1(1) = 0.
        r1(2) = 0.
       endif  

!....................................................

!      Carry out daily water balance accounting!
       if (abs(k(5)).gt.1e-8) then
        w(1)=w(1)+(dprec(d)+dmelt(d)-perc-r1(1)*aet)/k(5)
       end if

       if (abs(k(6)).gt.1e-8) then
        w(2)=w(2)+(perc-r1(2)*aet)/k(6)
       end if

!.....................................................

!      To prevent errors, set moisture back to wp

       if (w(1).le.0.) w(1)=0. 
       if (w(2).le.0.) w(2)=0.

!      If w2 is filled beyond fc then get drainage from layer two
       drainage=0.

       if(w(2).ge.1.)then
        drainage = (w(2)-1.)*k(6)
        w(2) = 1.
       endif

!      If w1 is filled above fc then get surface or sub-surface runoff
       runnoff=0.

       if(w(1).ge.1.)then
        runnoff = (w(1)-1.)*k(5)
        w(1) = 1.
       endif
!---------------------------------------------------------------------------
!      Sum the daily aet values:
       annaet = annaet + aet 

!      Sum the total runoff:
       sumoff=sumoff+runnoff+drainage

!      Sum the daily fvc value and average the daily gc value

       meanwr(month,1) = meanwr(month,1) + wr   / days(month)
       meanwr(month,2) = meanwr(month,2) + w(1) / days(month)
       meanwr(month,3) = meanwr(month,3) + w(2) / days(month)

       if (gc.ne.0.0) then
        meangc(month) =  meangc(month) +   gc / days(month)
       end if
       if (fvc.ne.0.0) then
        meanfvc(month)= meanfvc(month) +  fvc / days(month)
       end if
       meanaet(month)= meanaet(month) +  aet / days(month)

       wet(d)=wr
       dayfvc(d)=fvc

 120   continue
 110   continue
 100   continue
    
       return
       end
!***************************************************************************
!      Calculate the a generic phenology for any summergreen pft
!      A three month period centred around the coldest month is
!      defined as the minimum period during which foliage is not
!      present. Plants then start growing leaves and the end of this 
!      3 month period or when the temperature gos above 5oC if this
!      occurs later. Plants take 200 gdd5 to grow a full leaf canopy:

       subroutine phenology(dphen,dtemp,temp,tcm,tdif,tmin,pft,ddayl,pftpar)

       implicit none

       integer day,spinup,m,dayofmonth,ncm,daysinmonth(12)
       integer coldm(3),phencase,winter,flip,pft,hotm

       real dphen(365,2),dtemp(365),ramp(2),tcm,gdd,temp(12),ont
       real today,tdif,tmin,ddayl(365),warm,pftpar(25,25)

       data (daysinmonth(m),m=1,12) &
        / 31,28,31,30,31,30,31,31,30,31,30,31 /

!-----------------------------------------------------------------------

       ramp(1)=pftpar(pft,8)   

       if(pft.eq.7)then
        ont=0.0                  !this sets the minimum temp for growth
       else
        ont=5.0
       endif


!      Set grass ramp value:
       ramp(2)=pftpar(pft,9)

!      Find the months with the warmest and coldest temperatures (cm):
       warm=tcm
       do m=1,12
        if(temp(m).eq.tcm) ncm=m
        if(temp(m).gt.warm) then
         warm=temp(m)
         hotm=m
        end if
       end do

       do 10 phencase=1,2

       coldm(2)=ncm
       coldm(1)=coldm(2)-1
       coldm(3)=coldm(2)+1

       if (coldm(1).eq.0) coldm(1)=12
       if (coldm(3).eq.13) coldm(3)=1
       if (hotm.eq.12) hotm=0

       gdd = 0.
       winter=0

       do 20 spinup=1,2
        day = 0

        do 30 m=1,12
         do 40 dayofmonth=1,daysinmonth(m)            
          day=day+1

          if (dtemp(day).gt.ont) then
          if (m.ne.coldm(1).and.m.ne.coldm(2).and.m.ne.coldm(3)) then
           today=(dtemp(day))
           if(today.le.0.) today=0.
           gdd = gdd + today
           if (gdd.eq.0.0) then
            dphen(day,phencase) = 0.0
           else
            dphen(day,phencase) = gdd/ramp(phencase)
           end if
           if (gdd.ge.ramp(phencase)) dphen(day,phencase)=1.
           flip=1
          else
           if (flip.eq.1) winter=0
           winter=winter+1
           dphen(day,phencase) = 0.
           gdd  = 0.
           flip = 0
          end if
          end if

!         remove deciduous leaves in the fall when the temp or 
!         photoperiod reaches a threshold.

          if (phencase.eq.1) then
           if (m.ge.hotm) then
            if (dtemp(day).lt.-10.0.or.ddayl(day).lt.10.0) then
             dphen(day,phencase)=0.
            end if
           else if (m.eq.coldm(1)) then
             dphen(day,phencase)=0.
           end if
          else if (phencase.eq.2) then
           if (dtemp(day).lt.-5.0) dphen(day,phencase)=0. 
          end if
       
 40      continue   !daily loop ends
 30     continue    !montlhy loop ends
 20    continue     !annual loop ends
 10    continue     !case (grass,tree) loop ends

 100   continue
       
       return
       end
!*************************************************************************
!      subroutine to calc GDDs, TCM, wrin, and total precipitation

       subroutine climdata(cold,warm,gdd5,gdd0,rain,temp,prec,dtemp,alttmin)
       implicit none
       integer m,day
       real cold,gdd5,warm,gdd0,rain
       real temp(12),prec(12),dtemp(365)
       real minus0,minus5,gdd10,minus10
       real alttmin

       cold = 100.
       warm = -100.
       rain = 0.

       do m = 1,12
        if(temp(m).lt.cold) cold=temp(m)
        if(temp(m).gt.warm) warm=temp(m)
        rain = rain + prec(m)
       end do

       gdd10= 0.
       gdd5 = 0.
       gdd0 = 0.

       do day = 1,365
        minus10= dtemp(day) - 10.
        minus5 = dtemp(day) - 5.
        minus0 = dtemp(day)
        if(minus10.le.0.) minus10 = 0.
        if(minus5.le.0.)  minus5  = 0.
        if(minus0.le.0.)  minus0  = 0.
        gdd10 = gdd10 + minus10
        gdd5  = gdd5  + minus5
        gdd0  = gdd0  + minus0
       end do

       alttmin=(0.006*cold**2)+(1.316*cold)-21.9

       return
       end
!--------------------------------------------------------------------------
!      subroutine provides PFT specific parameters stored within subroutine

       subroutine pftdata(pftpar)

       implicit none

       integer iv,ip,npft,npar

       real var(25,25),pftpar(25,25)

       parameter(npar=11,npft=13)

!      Define all PFT specific parameters
!          1  = Phenological type 1=evergreen,2=summergreen,3=raingreen
!          2  = maximal value for minimum canopy conductance
!          3  = value for Emax, maximum daily transpiration rate 
!          4  = value of sw below which raingreen leaves drop
!          5  = value of sw above which raingreen leaves appear
!          6  = fraction of roots in top soil layer, 30 cm from Jackson et al.
!          7  = Expected leaf longevity in months
!          8  = Number of GDD5 required for full leaf out 
!          9  = Number of GDD0 required for full leaf out
!         10  = presence of sapwood respiration
!         11  = c4 plant or not
!    LIST OF ALL PLANT TYPES
!    1 = tet Tropical Evergreen Trees   
!    2 = trt Tropical Drought-deciduous Trees (raingreens) 
!    3 = tbe Temperate Broadleaved Evergreen Trees 
!    4 = tst Temperate Deciduous Trees
!    5 = ctc Cool Conifer Trees
!    6 = bec Boreal Evergreen Trees     
!    7 = bst Boreal Deciduous Trees   
!    8 = C3/C4 temperate grass plant type
!    9 = C4 tropical grass plant type
!    10= C3/C4 woody desert plant type
!    11= Tundra shrub type
!    12= cold herbaceous type
!    13= Lichen/forb type

       data ((var(iv,ip),ip=1,npar),iv=1,npft)     /   &
       1., 0.5, 10.0, -99.,-99. ,0.69 ,18. , -99.,-99., 1., 0., &
       3., 0.5, 10.0, 0.5, 0.6  ,0.70 , 9. , -99.,-99., 1., 0., &
       1., 0.2,  4.8, -99.,-99. ,0.67 ,18. , -99.,-99., 1., 0., &
       2., 0.8, 10.0, -99.,-99. ,0.65 , 7. , 200.,-99., 1., 0., &
       1., 0.2,  4.8, -99.,-99. ,0.52 ,30. , -99.,-99., 1., 0., &
       1., 0.5,  4.5, -99.,-99. ,0.83 ,24. , -99.,-99., 1., 0., &
       2., 0.8, 10.0, -99.,-99. ,0.83 ,24. , 200.,-99., 1., 0., &
       3., 0.8,  6.5,  0.2, 0.3 ,0.83 , 8. , -99.,100., 2., 1., &
       3., 0.8,  8.0,  0.2, 0.3 ,0.57 ,10. , -99.,-99., 2., 1., &
       1., 0.1,  1.0, -99.,-99. ,0.53 ,12. , -99.,-99., 1., 1., &
       1., 0.8,  1.0, -99.,-99. ,0.93 , 8. , -99.,-99., 1., 0., &
       2., 0.8,  1.0, -99.,-99. ,0.93 , 8. , -99., 25., 2., 0., &
       1., 0.8,  1.0, -99.,-99. ,0.93 , 8. , -99.,-99., 1., 0.   /

       do 10 iv = 1,npft
       do 10 ip = 1,npar
   10  pftpar(iv,ip) = var(iv,ip)

       return
       end
!**********************************************************************
!      Subroutine to get soil parameters given soil class number

       subroutine soildata(k,soil)
       implicit none
       integer i,j,stype
       real k(12),whc,d1,d2,store(3,12),soil(*)
       parameter(d1=300.,d2=1200.)
       data ((store(i,j),i=1,2),j=1,9) / &
!-------------------------------------------------------------
!       FAO soil texture data set values: 
!       1 coarse, 2 medium, 3 fine
!       4 medium-coarse, 
! 	5 fine-coarse, 
!	6 fine-medium, 
!	7 fine-medium-coarse,
!       8 Organic, 9 Ice
!       Selected FAO soil types: 
!       1 Kastozems,2 Chernozems,3 Vertisols, 4 Phaeozems
!
!       k1      whc  
       5., 	0.11,  &
       4.,	0.15,   &
       3.,	0.12,   &
       4.5,	0.13,   &
       4.,	0.115,  &
       3.5,	0.135,   &
       4.,	0.127, &
!       Organic soils:
       9.,	0.30,  &
!       Extra soil type:
       0.2,	0.10 /

!     *  5., 	0.11,  
!     *  4.,	0.15,   
!     *  3.,	0.12,
!     *  4.5,	0.13,   
!     *  4.,	0.115,  
!     *  3.5,	0.135,   
!     *  4.,	0.127,
!-------------------------------------------------------------       
!      Model should not be called for grid cells mapped as ice (soil=9)
       if(soil(1).ge.9.0) stop 'Value of soil type is not allowed!'

!      Assign the model soil type based on the FAO soil data base:
!      if zobler soil type is fine AND soil is a vertisol the define
!      heavy clay soil textural type:
       if(soil(2).eq.1.0.and.soil(1).eq.3.0)then
       stype=9
       else
       stype=soil(1)
       endif

!      Define water holding capacity for both soil layers (mm)
       whc  = store(2,stype) !real(soil(2))

       k(5) = whc*d1   !call this the top 30cm of a 1.5m deep soil
       k(6) = whc*d2   !and this the rest

!      Define k1 value (texture-dependent)
       k(1) = store(1,stype)
     
!      Define k2 value (not texture-dependent) 
       k(2) = 4.


       return
       end
!**********************************************************************
!     subroutine constraints provides environmental sieve

      subroutine constraints(tcm,twm,tminin,gdd5,rad0,pfts,tmin,maxdepth,gdd0)

      implicit none        

      integer npft,nclin,ip,iv,il
      integer pfts(13)

      parameter(npft=13,nclin=6)

      real limits(npft,nclin,2),clindex(nclin),undef
      real tmin,tcm,twm,ts,gdd5,rad0,gdd0,tminin,maxdepth

      parameter(undef=-99.9)

!------------------------------------------------------
!    LIST OF THE THIRTEEN PLANT FUNCTIONAL TYPES
!    1 = tet = Tropical Evergreen 
!    2 = trt = Tropical Raingreen
!    3 = wte = Temperate Broadleaved Evergreen   
!    4 = tst = Temperate Summergreen 
!    5 = ctc = Temperate Evergreen Conifer        
!    6 = bec = Boreal Evergreen       
!    7 = bst = Boreal Deciduous  
!    8 = temperate grass
!    9 = tropical/warm-temperate grass
!   10 = Desert woody plant type C3, C4
!   11 = Tundra shrub type
!   12 = Cold herbaceous type
!   13 = Lichen/forb type

!     Define and initialize the limits of the climatic indices 

      data(((limits(ip,iv,il),il=1,2),iv=1,5),ip=1,npft) / &
     -99.9,-99.9,   0.0,-99.9, -99.9,-99.9, -99.9,-99.9,  10.0,-99.9,  & !1
     -99.9,-99.9,   0.0,-99.9, -99.9,-99.9, -99.9,-99.9,  10.0,-99.9,  & !2
     -99.9,-99.9,  -8.0,  5.0,1200.0,-99.9, -99.9,-99.9,  10.0,-99.9,  & !3
     -15.0,-99.9, -99.9, -8.0,1200.0,-99.9, -99.9,-99.9, -99.9,-99.9,  & !4
      -2.0,-99.9, -99.9, 10.0, 900.0,-99.9, -99.9,-99.9,  10.0,-99.9,  & !5
     -32.5, -2.0, -99.9,-99.9, -99.9,-99.9, -99.9,-99.9, -99.9, 21.0,  & !6
     -99.9,  5.0, -99.9,-10.0, -99.9,-99.9, -99.9,-99.9, -99.9, 21.0,  & !7
     -99.9,-99.9, -99.9,  0.0, 550.0,-99.9, -99.9,-99.9, -99.9,-99.9,  & !8
     -99.9,-99.9,  -3.0,-99.9, -99.9,-99.9, -99.9,-99.9,  10.0,-99.9,  & !9
     -99.9,-99.9, -45.0,-99.9, 500.0,-99.9, -99.9,-99.9,  10.0,-99.9,  & !10
     -99.9,-99.9, -99.9,-99.9, -99.9,-99.9,  50.0,-99.9, -99.9, 15.0,  & !11
     -99.9,-99.9, -99.9,-99.9, -99.9,-99.9,  50.0,-99.9, -99.9, 15.0,  & !12
     -99.9,-99.9, -99.9,-99.9, -99.9,-99.9, -99.9,-99.9, -99.9, 15.0  / !13
!      ltcm, utcm,  lmin, umin,  lgdd, ugdd, lgdd0,ugdd0,  ltwm, utwm 

      data((limits(ip,6,il),il=1,2),ip=1,npft) / &
     -99.9,-99.9, & !1
     -99.9,-99.9, & !2
     -99.9,-99.9, & !3
     -99.9,-99.9, & !4
     -99.9,-99.9, & !5
     -99.9,-99.9, & !6
     -99.9,-99.9, & !7
     -99.9,-99.9, & !8
     -99.9,-99.9, & !9
     -99.9,-99.9, & !10
      15.0,-99.9, & !11
     -99.9,-99.9, & !12
     -99.9,-99.9 / !13
!     lsnow,usnow

!     Assign tmin value
      if (tminin.le.tcm) then
       tmin=tminin
      else
       tmin=tcm-5.0
      end if

!     Assign ts value
      ts=twm-tcm

!     set up climate indices array:
      clindex(1)=tcm      !temperature of the coldest month
      clindex(2)=tmin     !absolute minimum temperature
      clindex(3)=gdd5     !GDDays above 5 deg C
      clindex(4)=gdd0     !gdd0
      clindex(5)=twm
      clindex(6)=maxdepth !snow depth

! Determines the values of the climatic indices are within the climatic limits 

      do 100 ip=1,npft
         do 101 iv=1,nclin
!
! both limits given, value inside - PRESENT:
!
            if(limits(ip,iv,1).le.clindex(iv).and. &
              limits(ip,iv,1).ne.undef.and. &
              limits(ip,iv,2).ne.undef.and. &
              limits(ip,iv,2).gt.clindex(iv))then
               pfts(ip)=1
!
! or lower limit missing, value below upper level - PRESENT:
!
            elseif(limits(ip,iv,1).eq.undef.and.  &
                  limits(ip,iv,2).ne.undef.and. &
                  limits(ip,iv,2).gt.clindex(iv))then
               pfts(ip)=1
!
! or upper limit missing, value above lower level - PRESENT:
!
            elseif(limits(ip,iv,1).le.clindex(iv).and. &
                  limits(ip,iv,1).ne.undef.and. &
                  limits(ip,iv,2).eq.undef)then
               pfts(ip)=1
!
! both limits missing - PRESENT:
!
            elseif(limits(ip,iv,1).eq.undef.and. &
                  limits(ip,iv,2).eq.undef)then
               pfts(ip)=1
!
! none of these - ABSENT:
!
            else
               pfts(ip)=0

               goto 100
            endif
101         continue
100      continue

       return
       end
!******************************************************************************
!     subroutine snow masks precip to account for effects of snow

      subroutine snow(dtemp,dprec,dmelt,dprecin,maxdepth)
      implicit none
      integer day,it
      real dtemp(1:365),dprec(1:365),dprecin(1:365)
      real tsnow,km,snowpack,snowmelt,newsnow,drain
      real dmelt(365),sum1,sum2,maxdepth
      parameter(tsnow=-1.,km=0.7)
      
      snowpack = 0.0
      maxdepth = 0.0

      do it =1,2
       sum1=0.
       sum2=0.

       do day=1,365
   
        drain = dprecin(day) / (365./12.)

!       Calculate snow melt and new snow for today
        if(dtemp(day).lt.tsnow)then
         newsnow  = drain
         snowmelt = 0.
        else
         newsnow  = 0.
         snowmelt = km*(dtemp(day)-tsnow)
        endif

!       Reduce snowmelt if greater than total snow remaining
        if (snowmelt.gt.snowpack) snowmelt = snowpack
    
!       Update snowpack store
        snowpack = snowpack + newsnow - snowmelt
        if (snowpack.gt.maxdepth) maxdepth=snowpack

!       Calculate effective water supply (as daily values in mm/day)
        dprec(day) = drain - newsnow 
        dmelt(day) = snowmelt

        sum1=sum1+dprec(day)+dmelt(day)
        sum2=sum2+drain

       end do
      end do

      return
      end
!******************************************************************************
!      Calculates insolation and PET for each month
 
       subroutine  ppeett (lat,dtemp,dclou,dpet,temp,sun,dayl,rad0,ddayl,radanom)

       implicit none

       integer month,day,dayofm,midday(12),daysinmonth(12)

       real dtemp(365),dclou(365),lat
       real dpet(365),ddayl(365),temp(12)
       real sun(12),dayl(12),rad0
       real dip,pie,a,sat,cla,sla,ho,rl,fd,qo,rs
       real psi,l,b,radup,qoo,c,d,albedo,hos,u,v,us,vs
       real radanom(12)

       parameter(b=0.2,radup=107.,qoo=1360.,d=0.5,c=0.25)
       parameter(albedo=0.17)

       data (midday(month),month=1,12) &
        / 16,44,75,105,136,166,197,228,258,289,319,350 /
       data (daysinmonth(month),month=1,12) &
        / 31,28,31,30,31,30,31,31,30,31,30,31          /

       pie = 4.*atan(1.)
       dip = pie/180.
    
!      Daily loop
       day=0
       rad0=0.
       do 10 month  = 1,12
       do 20 dayofm = 1,daysinmonth(month)
       day=day+1

!      Find psi and l for this temperature from lookup table
!      psychrometer constant (pa/oc), latent heat lamba (mj/kg) 
       call table(dtemp(day),psi,l)
    
!      Calculation of longwave radiation       
       rl = (b + (1-b)*(dclou(day)/100.))*(radup- dtemp(day))

!      Since changes in radiation (short or long) will mainly be due
!      to changes in cloudiness, apply the (short wave) anomaly here too.
!      Per B. Smith 1998
       
       rl=rl*radanom(month)

!      c=0.29*cos(lat) to emphasize the effect of clouds at high latitude
!      c=0.29*cos(lat*dip)

!      Calculation of short wave radiation
       qo  =  qoo*(1.+2.*0.01675*cos(dip*(360.*real(day))/365.))
       rs  =  qo*(c+d*(dclou(day)/100.))*(1.-albedo)

       rs=rs*radanom(month)

       a   = -dip*23.4*cos(dip*360.*(real(day)+10.)/365.)
       cla =  cos(lat*dip)*cos(a)
       sla =  sin(lat*dip)*sin(a)
       u = rs*sla - rl
       v = rs*cla      

!      Check for polar day and polar night
       if(u.ge.v)then
!      polar day:
       ho = pie
       elseif(u.le.(0.-v))then
!      polar night:
       ho = 0.
       else
!      normal day and night: (find ho the time of dawn)
       ho =  acos(-u/v)
       endif
     
!      Equations for demand function
       sat=(2.5*10**6.*exp((17.27*dtemp(day))/(237.3+dtemp(day)))) &
               /((237.3+dtemp(day))**2.)
!      Multiply l by e6 to convert from mj/kg to j/kg 
       fd = (3600./(l*1e6))*(sat/(sat+psi))

!      Store total daily equilibrium transpiration rate as dpet
       dpet(day)=fd*2.*((rs*sla-rl)*ho+rs*cla*sin(ho))/(pie/12.)

!      Calculate daylength in hours
       if (ho.eq.0.0) then
        ddayl(day)=0.0
       else
        ddayl(day) = 24.*(ho/pie)
       end if

!      If at a mid-month day then record mid-month daily sun and dayl
       if(day.eq.midday(month))then


!        First record the day length 
         dayl(month)=ddayl(day)

!        Now calculate daily total irradiance (j/m2) & record in sun
         us = rs*sla
         vs = rs*cla
!        check for polar day and polar night
         if(us.ge.vs)then
!        polar day:
         hos = pie
         elseif(us.le.(0.-vs))then
!        polar night (also h1=0. for polar night)
         hos = 0.
         else
!        normal day and night, find hos the time of dawn
         hos =  acos(-us/vs)
         endif

!        Find total insolation for this day, units are j/m2
         sun(month)= 2.*(rs*sla*hos+rs*cla*sin(hos))*(3600.*12./pie)
!        Do not allow negative values for insolation 
         if(sun(month).le.0.) sun(month)=0.

!        Sum total annual radiation for months with t>0oC (GJs PAR year-1)
!        (assuming 50% of short wave radiation is PAR)
         if(temp(month).gt.0.)then
         rad0=rad0+real(daysinmonth(month))*sun(month)*1e-9*0.5
         endif

       endif

 20    continue
 10    continue

       return
       end
!***************************************************************************
!      subroutine table from bucket subroutine

      subroutine table(tc,gamma,lambda)

! looks up gamma and lambda from table (essential part of EVAPO.F)

! Author: Wolfgang Cramer, Dept. of Geography, Trondheim University-AVH,
! N-7055 Dragvoll, Norway.

! latest revisions 14/2-1991

! enable this when you run on a compiler allowing for it:
      implicit none

! on UNIX, please compile with "f77 -u"

      integer ir,il
      real gbase(2,11),lbase(2,11)
      real gamma,lambda,tc

      data ((gbase(ir,il),ir=1,2),il=1,11) &
      /-5.,64.6, 0.,64.9, 5.,65.2,10.,65.6,15.,65.9,20.,66.1, &
       25.,66.5,30.,66.8,35.,67.2,40.,67.5,45.,67.8/
      data ((lbase(ir,il),ir=1,2),il=1,11) &
      /-5.,2.513, 0.,2.501, 5.,2.489,10.,2.477,15.,2.465,20.,2.454, &
       25.,2.442,30.,2.430,35.,2.418,40.,2.406,45.,2.394/

! temperature above highest value - set highest gamma and lambda and return

      if(tc.gt.gbase(1,11)) then
         gamma=gbase(2,11)
         lambda=lbase(2,11)
         return
      endif

! temperature at or below value - set gamma and lambda

      do 100 il=1,11
         if(tc.le.gbase(1,il)) then
            gamma=gbase(2,il)
            lambda=lbase(2,il)
            return
         endif
100   continue

      end
!***************************************************************************
      subroutine daily(mly,dly)
      implicit none
      real mly(12),dly(365),midday(12),vinc
      integer im,id
      data (midday(im),im=1,12)/16., 44., 75.,105.,136.,166., &
                              197.,228.,258.,289.,319.,350./

      vinc=(mly(1)-mly(12))/31.0
      dly(350)=mly(12)
      do 100 id=351,365
         dly(id)=dly(id-1)+vinc
100      continue
      dly(1)=dly(365)+vinc
      do 101 id=2,15
         dly(id)=dly(id-1)+vinc
101      continue
      do 103 im=1,11
         vinc=(mly(im+1)-mly(im))/(midday(im+1)-midday(im))
         dly(int(midday(im)))=mly(im)
         do 104 id=int(midday(im))+1,int(midday(im+1))-1
            dly(id)=dly(id-1)+vinc
104         continue
103      continue
      return
      end   
!**************************************************************************
!      subroutine isotope(Cratio,Ca,temp,Rd,c4month,mgpp,phi,meanC3,meanC4,meanC,C3DA,C4DA,gpp)
       subroutine isotope(C3ratio,C4ratio,Ca,temp,Rd,IsC4month,C3npp,C4npp,phi,meanC3,meanC4,meanC,C3DA,C4DA,mIsnpp)
! je préfère changer aussi les argguments c4month et gpp parce que IsC4month veut dire presence de C4 alors que C4month 
! veut dire présence majoritairement de C4, et Isgpp veut dire gpp de tout alors que gpp veut dire gpp des plantes qui 
! sont mjoritaires par mois.
! je trouve plus logique d'utiliser les npp plutot que les gpp: le gpp correspond à tout le CO2 entrant, le npp est 
! plutôt le carbone fixé après utilisation des métabolites intermédiaires.
!     This subroutine is for calculating the total fractionation of 13C
!     as it goes from free air (as 13CO2) to fixed carbon in the leaf.
!     For use with the BIOME3 model of A. Haxeltine (1996).
!     There are separate routines for calculating fractionation by both
!     C3 and C4 plants.  This program is based upon the model used by
!     Lloyd and Farquhar (1994).

      implicit none

      logical Isc4month(12)

      integer m

!! C. Hatte
      real C3ratio(12),C4ratio(12),Ca,temp(12),Rd(12),C4npp(12),C3npp(12), mIsnpp
!      real Cratio(12),Ca,temp(12),Rd(12),C4gpp(12), mgpp(12), gpp
      real C3DA(12),C4DA(12),meanC3,meanC4,meanC
      real wtC3,wtC4,mC3,mc4
	  real SC4npp,SC3npp, partC3, partC4 
      real delC3,delC4,phi
		common /izzo/mC3,mC4,partc3,partc4
!      open (unit=1,file='del13C3.out',status='unknown')
!      open (unit=2,file='del13C4.out',status='unknown')



      wtC3=0.0
      wtC4=0.0
      SC4npp=0.0 
      SC3npp=0.0 
      do m=1,12

!      if (mgpp(m).gt.0.0) then
	   if (C3npp(m).gt.0.0.or.C4npp(m).gt.0.0) then

       if (C3ratio(m).lt.0.05) C3ratio(m)=0.05

       if (Isc4month(m)) then
        call isoC4(C4ratio(m),phi,temp(m),delC4)
        C4DA(m)=delC4
!		wtc4 = wtc4+delC4*mgpp(m)
		wtc4 = wtc4+delC4*C4npp(m)
!        SC4gpp=SC4gpp+mgpp(m) 
        SC4npp=SC4npp+C4npp(m) 
        call isoC3(C3ratio(m),Ca,temp(m),Rd(m),delC3)
        C3DA(m)=delC3
!        wtC3 = wtC3+delC3*mgpp(m)
        wtC3 = wtC3+delC3*C3npp(m)
!		SC3gpp=SC3gpp+mgpp(m) 
		SC3npp=SC3npp+C3npp(m) 
   
       else
        call isoC3(C3ratio(m),Ca,temp(m),Rd(m),delC3)
        C3DA(m)=delC3
        wtC3 = wtC3+delC3*C3npp(m)
		SC3npp=SC3npp+C3npp(m) 
       end if
      else
       C3DA(m)=0.0
       C4DA(m)=0.0 
      end if
      end do

      meanC3 = wtC3/amax1(SC3npp,1.0)
      meanC4 = wtC4/amax1(SC4npp,1.0)
      meanC = (wtc3+wtc4)/mIsnpp
	  mC3=meanC3
	  mC4=meanC4
      partC3=SC3npp/mIsnpp
      partC4=SC4npp/mIsnpp
      return
      end

!----------------------------------------------------------------
!     This part calculates fractionation for C3 photosynthesis.
      subroutine isoC3(Cratio,Ca,temp,Rd,delC3)
  
      implicit none
      real DeltaA,delC3
      real a,es,a1,b,e,k,f,gamma,Catm
      real Cratio,Ca,Rd,temp,leaftemp,co2
      real q,r,s,t

!     define fractionation parameters

	co2=ca*1.e+6
       a= 4.4
      es= 1.1
      a1= 0.7
       b=27.5
       e= 0.0
!       f= 8.0
! correction C. Hatteco2
       b=27.5+(25-temp)/2
       e=1.0
	   f=(-8.0+(350-Co2)/7)
      Catm= 0.0 

      if (Rd.le.0) Rd=0.01
      leaftemp = 1.05*(temp+2.5)
      gamma = 1.54*leaftemp
      Rd = Rd/(86400.0*12.0)
      Catm = Ca*1.0e6
      k = Rd/11.0            !From Farquhar et al. 1982 p. 126

!     calculate the fractionation

      q = a*(1-Cratio+0.025)
      r = 0.075*(es+a1)
      s = b*(Cratio-0.1)
      t = (e*Rd/k+f*gamma)/Catm    
!		t=0.0

      DeltaA = q+r+s-t

      delC3 = DeltaA
      return
      end
!---------------------------------------------------------------
!     This part calculates fractionation for C4 photosynthesis.

      subroutine isoc4(Cratio,phi,temp,delC4)

      implicit none
      real DeltaA,delC4
      real a,es,a1,b4,b3,phi,temp
      real Cratio

         a= 4.4
        es= 1.1
        a1= 0.7
        b3=30.0
!       phi= 0.2

      b4=(26.19-(9483/(273.2+temp)))
 
      DeltaA=a*(1-(Cratio)+0.0125)+0.0375*(es+a1)+(b4+(b3-es-a1)*phi)*((Cratio)-0.05)

      delC4 = DeltaA
      
      return
      end
!--------------------------------------------------------------------
!     This subroutine is for calculating the phi variable used in 
!     C4 photosynethsis isotope fractionation calculations

      subroutine calcphi(gpp,phi)

      implicit none
      real gpp(12),totgpp,meangpp,normgpp(12)
      real snormavg(4),svar(4)
      real avar,phi,a
      integer m,s

      totgpp=0.0       !initialize a few variables
      do s=1,4
       svar(s)=0.0
      end do

!     This first part of the subroutine estimates annual variability of
!     GPP first by normalizing and then summing seasonal variability
!     which compensates for amplitude and seasonal variation in GPP.

      do m=1,12
       totgpp=totgpp+gpp(m)
      end do

      meangpp=totgpp/12.0

      do m=1,12
       normgpp(m)=gpp(m)/meangpp
      end do 

      snormavg(1)=(normgpp(1)+normgpp(2)+normgpp(3))/3.0
      snormavg(2)=(normgpp(4)+normgpp(5)+normgpp(6))/3.0
      snormavg(3)=(normgpp(7)+normgpp(8)+normgpp(9))/3.0
      snormavg(4)=(normgpp(10)+normgpp(11)+normgpp(12))/3.0

!     calculate the population variances by season

      do m=1,3
       a=((normgpp(m)-snormavg(1))**2)/3
       svar(1)=svar(1)+a
      end do

      do m=4,6
       a=((normgpp(m)-snormavg(2))**2)/3
       svar(2)=svar(2)+a
      end do

      do m=7,9
       a=((normgpp(m)-snormavg(3))**2)/3
       svar(3)=svar(3)+a
      end do

      do m=10,12
       a=((normgpp(m)-snormavg(4))**2)/3
       svar(4)=svar(4)+a
      end do

      avar=svar(1)+svar(2)+svar(3)+svar(4)
!------------------------------------------------

!     This part sets the phi value based upon the annual variability.
!     The equation is a simple regresion based upon hypothetical extreme
!     scenarios of phi.

      phi=0.3518717*avar+0.2552359

      if (phi.ge.1.0) phi=phi/10.0

      return

      end

!---------------------------------------------------------------------------
!
!     This subroutine models heterotrophic respiration of litter and soíl
!     organic carbon in both a fast and a slow pool.  It assumes equilibrium 
!     and so decays all of a given year's NPP.  The 13C composition of respired
!     CO2 is also modelled.  Models are based on the work of Foley, Lloyd and
!     Taylor, and Sitch.

      subroutine hetresp(pft,nppann,tair,tsoil,aet,moist,Rlit,Rfst,Rslo,Rtot, &
      isoveg,isoR,isoflux,Rmean,meanKlit,meanKsoil)

      implicit none

      real tair(12),moist(12),tsoil(12)
      real Plit,Pfst,Pslo,Rtot(12),aet(12)
      real Klit(12),Kfst(12),Kslo(12)
      real Rlit(12),Rfst(12),Rslo(12),isoR(12),Rmean
      real klitsum,kfstsum,kslosum
      real Rten,mfact,nppann,isoveg,isoatm
      real isolit(12),isofst(12),isoslo(12),isoflux(12)
      integer m,pft
      
      real meanKlit,meanKsoil

      parameter(isoatm=-8.0)
     
!     P is pool sizes for partitioning, R is respired CO2

!     the soil temp subroutine must have been called by now 

!     partition annual npp into pools according to Foley strategy

      if (nppann.le.0.0) then
       do m=1,12
        Rlit(m)=0.0
        Rfst(m)=0.0
        Rslo(m)=0.0
        Rtot(m)=0.0
        isoR(m)=0.0
        isoflux(m)=0.0
       end do
       isoveg=0.0
       return

      else               !begin the real routine here

       if (pft.eq.1.or.pft.eq.2) then
         Plit=0.650*nppann
         Pfst=0.980*0.350*nppann
         Pslo=0.020*0.350*nppann
       else
         Plit=0.700*nppann
         Pfst=0.985*0.300*nppann
         Pslo=0.015*0.300*nppann
       end if

!     Calculate respiration for each pool with an R10 base resp.
!     Litter needs to decay according to a basic temp and moist function.
!     Soil decay can be calculated according to temp. response of
!     Lloyd and moisture of Foley with a turnover time built into the Rten

!     Two ways to decay NPP, one based on surface temp and AET for litter 
!     (Foley).  The other is for soil decay and is based on soil 
!     temperature and moisture.

      Rten=1.0
      Klitsum=0.0
      Kfstsum=0.0
      Kslosum=0.0

      do m=1,12
       mfact=0.25+0.75*moist(m)

       Klit(m)=1.0*10.0**(-1.4553+0.0014175*aet(m))
       Klitsum=Klitsum+Klit(m)

       Kfst(m)=mfact*Rten*EXP(308.56*((1/56.02)-(1/(tsoil(m)+273.-227.13))))
       Kfstsum=Kfstsum+kfst(m)

       Kslo(m)=mfact*Rten*EXP(308.56*((1/56.02)-(1/(tsoil(m)+273.-227.13)))) 
       Kslosum=Kslosum+Kslo(m)
      end do

!-----added 23.07.99---------------

      meanKlit=Klitsum/12.
      meanKsoil=Kfstsum/12.

!----------------------------------

      Rmean=0

      do m=1,12
       Rlit(m)=Plit*(Klit(m)/Klitsum)
       Rfst(m)=Pfst*(Kfst(m)/Kfstsum)
       Rslo(m)=Pslo*(Kslo(m)/Kslosum)
       Rtot(m)=Rlit(m)+Rfst(m)+Rslo(m)
       Rmean=Rmean+(Rtot(m)/12.)
      end do

!     calculate the isotope ratio of respired CO2 based on
!     the NPP weighted mean 13C in the vegetation
!     Since 13C is enriched in organic matter over time add factors

      do m=1,12
       isolit(m)=isoveg-0.75  !these factors represent enrichment
       isofst(m)=isoveg-1.5
       isoslo(m)=isoveg-2.25
       isoR(m)=((Plit/nppann)*isolit(m))+ &
        ((Pfst/nppann)*isofst(m))+((Pslo/nppann)*isoslo(m))
       isoflux(m)=(isoatm-isoR(m))*Rtot(m)
      end do

      end if
      return

      end 

!-----------------------------------------------------------------
!     This subroutine calculates monthly mean soil temperature 
!     based on monthly mean air temperature assuming a thermal
!     conductivity of the soil and a time lag between soil and air 
!     temperatures. Based on work by S. Sitch.

      subroutine soiltemp(tair,soiltext,tsoil)

      implicit none

      real tair(12),tsoil(12)
      real diffus,damp,amp,lag,pie
      real therm(9),soiltext(4)
      real sumtemp,meantemp
      integer m,i

      pie = 4.*ATAN(1.)    

      data (therm(i),i=1,9)  &
     / 8.0,4.5,1.0,5.25,4.5,2.75,1.0,1.0,8.0 / !check value for soil 8

      sumtemp=0.
!----------------------------------
!     calculate a soil-texture based thermal cond. and lag time

      diffus = therm(2)
      damp=0.25/(sqrt(diffus))
      lag=damp*(6/pie)
      amp=exp(-damp)

!     calculate mean annual air temperature

      do m=1,12
       sumtemp=sumtemp+tair(m)
      end do
      meantemp=sumtemp/12.
      
!     calculate soil temperature
 
      tsoil(1) = (1.-amp)*meantemp+amp*(tair(12)+(1.-lag)*(tair(1)-tair(12)))
     
      do m=2,12
       tsoil(m) = (1.-amp)*meantemp+amp*(tair(m-1)+(1.-lag)*(tair(m)-tair(m-1)))
      end do

!     due to snow cover don't allow soil temp < -10

      do m=1,12
       if (tsoil(m).lt.-10.) tsoil(m)=-10.
      end do

      return

      end

!---------------------------------------------------------------------
      subroutine fire (wet,pft,lai,npp,firedays)

!     This subroutine calculates the number of potential fire days
!     in a year based on threshold values for soil moisture, which are
!     PFT dependent.  Jed Kaplan 1998.

!     May 1998 now includes a parameter to scale firedays in terms of
!     annual NPP so that firedays are reduced linearly to 0 below 1000gC/m2

      implicit none

      real wet(365),threshold(13),firedays,burn(365)
      real wetday,dryday,lai,npp,litter
      real firefraction,burnfraction

      integer day,pft,i

      data (threshold(i),i=1,13) &
        /0.25,0.20, 0.40,0.33,0.40, 0.33,0.33, 0.40,0.40, 0.33, &
         0.33,0.33,0.33/

      firedays =   0.0
      wetday   =   0.0
      dryday   = 100.0

      do day=1,365
    
       if (wet(day).lt.threshold(pft)) then
        burn(day)=1.0
       else if (wet(day).gt.threshold(pft)+0.05) then
        burn(day)=0.0
       else
        burn(day)=1/exp(wet(day)-threshold(pft))
       end if

       if (wet(day).gt.wetday) wetday=wet(day)
       if (wet(day).lt.dryday) dryday=wet(day)

       firedays=firedays+burn(day)

      end do

      firefraction=firedays/365.0

      litter=(lai/5.0)*npp

      burnfraction=litter*(1.-(EXP(-0.2*firefraction**1.5))**1.5)

      if (npp.lt.1000.0) then
       firedays=firedays*(npp/1000.0)
      end if

      return

      end


