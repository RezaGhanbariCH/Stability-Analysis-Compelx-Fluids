
c     storage of varaibles
c=====================================================================================================================       
c  p | ur | ut | uz | t | trr1 | ...  | trrn | trt1 | ... | trtn | trz1 | ... | trzn | ttt1 | ... | tttn | ttz1 |....
c=====================================================================================================================    
c	Linear stability analysis for dean flow problem
c	Started March 26,2002.
c	Completed April 3,2002
c	Giesekus fluid
c
c	---------------MAIN PROGRAM-------------------------

	subroutine lsataylor(k,q)
	implicit double precision (a-h,o-z)
	parameter(iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
	parameter(nvars=(6*nmodes+5)*n,nrhs=1,nn=2*nvars)
	parameter(lwork=33*nvars)
	common/variables/ eigvar(nvars),deigvar(nvars),
     &	ddeigvar(nvars),y(n),r(n),der(n,n),dder(n,n),dfirst(n,n),
     &	dsecond(n,n)
      common/steadyvar/uthss(n),temss(n),trrss(n,nrows),
     & 	trthss(n,nrows), tththss(n,nrows),duthss(n),dtemss(n),
     &	dduthss(n),ddtemss(n),dtrrss(n,nrows),ddtrrss(n,nrows),
     &	dtrthss(n,nrows),ddtrthss(n,nrows),dtththss(n,nrows),
     &	ddtththss(n,nrows)
	common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &	epsilonlam(nrows),alp(nrows),delham(nrows),
     &	delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond
	common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &    gammadot
c	double precision L
	common/chain/L
	common/wavenum/si,alpha
	common/mat/AA(nvars,nvars),BB(nvars,nvars)
	common/remat/AAA(nn,nn),BBB(nn,nn)
	double precision sigmar(nn),sigmai(nn),beta(nn),z(nn,nn),
     &	axial(60),sum,realeig(nvars),mageig(nvars),eigr,eigi
	complex*16 AA, BB,eig(nn)
	complex*16 alphaa(nvars),betaa(nvars),vl(1,nvars),
     &	vr(nvars,nvars),work(lwork)
	double precision rwork(8*nvars),eigmaxr,eigmaxi,q(10)
	integer i,j,k,matz,col(nn),pivnum
	
	data pi/3.1415926535897930d0/
c	data L/1000.d0/

      write(*,*) ' L', L
	
       

	matz=1.d0	
	call inputvalue
	
c	write(*,*) 'UINVEL=',uin,'UOUTVEL=',uout
	call meshpoint(y,r,nmodes,n,delta)
	call derivative(der,y,n)
	call squarematrix(der,n,dder)
	call steadystate
c			write(*,*) 'L', L
	call constructmatA(AA,q)	
	call constructmatB(BB,q)
	call realmatrix

	
c	write(*,*) nn	

c	do i=nvars+1,nn
c	do j=1,nn
c	if (BBB(i,j).NE.0.d0) then
c	write(*,*) i,j,BBB(i,j)
c	endif
	
c	enddo
c	enddo

	call zggev('N','V',nvars,AA,nvars,BB,nvars,alphaa,betaa,vl,1,vr
     $           ,nvars,work,lwork,rwork,info)

      
	eigmaxr = -100.d0
      k=0
      do i=1,nvars
        if(zabs(betaa(i)).GE.1.d-14) then
          k=k+1
          col(k)=i
          eig(i)=alphaa(i)/betaa(i)
          x1=dreal(eig(i))
          x2=dimag(eig(i))
          
		if (x1.gt.eigmaxr) then
            eigmaxr = x1
            eigmaxi = x2
            pivnum = i
          endif
          
		write(14,'(12x,2e26.12)') x1,x2
        endif
      enddo
      close(unit=14)
      
	
	do i = 1,nvars
        write(50,*) dreal(vr(i,pivnum)),dimag(vr(i,pivnum))
      enddo
      close(unit=50)

c	do i=1,nvars
c       read(50,*) realeig(i),mageig(i)
c     enddo
c     close(unit=50)

c	-------------rgg----------
c	call rgg(nn,nn,AAA,BBB,sigmar,sigmai,beta,matz,z,ierr)
	
c	k=0
c	do i=1,nn
c	  if(dabs(beta(i)).GE.1.d-14) then
c	  k=k+1
c	  col(k)=i
c	  eig(i)=dcmplx(sigmar(i),sigmai(i))/beta(i)
c	  x1=dreal(eig(i))
c	  x2=dimag(eig(i))

c	  write(*,*) x1,x2,i

c	  write(14,'(12x,2e26.12)') x1,x2 
c	  endif
c	enddo

c	write(*,*) k
c	close(unit=14)
c	call iterate(k,pivnum)

c	----------rgg end----------------

c	write(*,*) 'UinROT=',uin,'UoutROT=',uout


c	do i=1,k
c       write(15,*) i,col(i)
c     enddo
c     write(*,*) pivnum,col(pivnum)
c     write(*,*) sigmar(col(pivnum))/beta(col(pivnum))
c     write(*,*) sigmai(col(pivnum))

c     do i=1,nn
c       write(49,*) z(i,col(pivnum))
c     enddo


c     do i=1,51
c       axial(i)=2.d0*(i-1)*pi/(alpha*50.d0)
c     enddo

c     do i=1,n
c       write(50,*) r(i),z(i,col(pivnum))
c     enddo

c     do i=1,n
c       sum=0.d0
c       do j=1,n
c       sum=sum+dfirst(i,j)*z(j,col(87))
c     enddo
c     write(51,*) r(i),sum
c     enddo


c	do i=1,nvars
c       read(50,*) realeig(i),mageig(i)
c     enddo
c     close(unit=50)
	
c     do i=1,n
c       do j=1,51
c         write(55,*) r(i), axial(j),z(n+i,col(87))*dcos(alpha*axial(j))-
c     &   z(nvars+n+i,col(87))*dsin(alpha*axial(j))

c	    write(55,*) r(i),axial(j),realeig(i+n)*dcos(alpha*
c     &   axial(j)+0.d0*2.d0*pi*1.d0/2.d0)-mageig(i+n)*
c     &   dsin(alpha*axial(j)+0.d0*2.d0*pi*1.d0/2.d0)
c       enddo
c     enddo

c	do i=1,n
c	  do j=1,51
c	    j=1
c	    write(50,*) r(i), -z(i,col(pivnum))*dcos(alpha*axial(j))+z(nvars
c     &	+n+i,col(pivnum))*dsin(alpha*axial(j))+2.d0*(betas/r(i))*(
c    &	z(n+i,col(pivnum))*dcos(alpha*axial(j))-z(nvars+n+i,col(pivnum))
c     &	*dsin(alpha*axial(j)))+(1.d0-betas)*(z(8*n+i,col(pivnum))*dcos(
c     &	alpha*axial(j))-z(nvars+8*n+i,col(pivnum))*dsin(alpha*axial(j)))
cc	  enddo
c	enddo

c     do i=1,n
c       write(50,*) r(i),z(i,col(41)),z(nvars+i,col(41))
c       write(51,*) r(i),z(n+i,col(41)),z(nvars+n+i,col(41))
c       write(52,*) r(i),z(2*n+i,col(41)),z(nvars+2*n+i,col(41))
c       write(53,*) r(i),z(3*n+i,col(41)),z(nvars+3*n+i,col(41))
c       write(54,*) r(i),z(4*n+i,col(41)),z(nvars+4*n+i,col(41))
c     enddo





	return
	end
c	---------------END OF MAIN PROGRAM-----------------

	subroutine iterate(eigennum,pivnum)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
      parameter(nvars=(6*nmodes+5)*n,nrhs=1,nn=2*nvars)

      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &            epsilonlam(nrows),alp(nrows),delham(nrows),
     &            delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond    
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &           gammadot
c	double precision L
	common/chain/L
      common/wavenum/si,alpha

      double precision eigenx(nn), eigeny(nn),temp,pivot
      integer i,j,eigennum,pivnum

c     write(*,*) eigennum
      do i=1,eigennum
        read(14,*) eigenx(i),eigeny(i)
      enddo
	close(unit=14)
      
	j=1
      pivot=eigenx(1)
      temp=pivot
      do i=2,eigennum
        if (eigenx(i).GE.pivot) then
          pivot=eigenx(i)
          j=i
        else
          temp=pivot
          pivot=temp
        endif
      enddo

c     write(*,*) pivot,j, eigeny(j)
      pivnum=j
      
	return
      end






c	--------------SUBROUTINE REALMATRIX---------------
	subroutine realmatrix
	implicit double precision (a-h,o-z)
      parameter(iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
      parameter(nvars=(6*nmodes+5)*n,nrhs=1,nn=2*nvars)
      common/variables/ eigvar(nvars),deigvar(nvars),
     &      ddeigvar(nvars),y(n),r(n),der(n,n),dder(n,n),dfirst(n,n),
     &	  dsecond(n,n)
      common/steadyvar/uthss(n),temss(n),trrss(n,nrows),
     &  trthss(n,nrows), tththss(n,nrows),duthss(n),dtemss(n),
     &  dduthss(n),ddtemss(n),dtrrss(n,nrows),ddtrrss(n,nrows),
     &  dtrthss(n,nrows),ddtrthss(n,nrows),dtththss(n,nrows),
     &  ddtththss(n,nrows)
      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &  epsilonlam(nrows),alp(nrows),delham(nrows),
     &  delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &  gammadot
c	double precision L
	common/chain/L
      common/wavenum/si,alpha
	common/mat/AA(nvars,nvars),BB(nvars,nvars)
	common/remat/AAA(nn,nn),BBB(nn,nn)
	complex*16 AA,BB
      integer i,j,k

	data pi/3.1415926535897930d0/
	
	do i=1,nn/2
	  do j=1,nn/2
	    AAA(i,j)=dreal(AA(i,j))
	    BBB(i,j)=dreal(BB(i,j))
	    AAA(i,j+nn/2)=-dimag(AA(i,j))
	    BBB(i,j+nn/2)=-dimag(BB(i,j))
	    AAA(i+nn/2,j)=dimag(AA(i,j))
	    BBB(i+nn/2,j)=dimag(BB(i,j))
	    AAA(i+nn/2,j+nn/2)=dreal(AA(i,j))
	    BBB(i+nn/2,j+nn/2)=dreal(BB(i,j))
        enddo
	enddo
		
	end
		
c	------------END OF SUBROUTINE REALMATRIX---------
c      ----------------SUBROUTINE INPUTVALUE ----------------

        
	subroutine inputvalue

      implicit double precision (a-h,o-z)
      parameter (iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
      parameter (nvars = (6*nmodes+5)*n,nrhs=1,nn=2*nvars)
	common/variables/ eigvar(nvars),deigvar(nvars),
     &	ddeigvar(nvars),der(n,n),dder(n,n),dfirst(n,n),dsecond(n,n)
      common/poly/lamdap(nrows)
      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &  epsilonlam(nrows),alp(nrows),delhlam(nrows),
     &  delhp(nrows),delhs,epsilons,betas
   	common/physical/visc,dens,cond 
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uout
     &      gammadot
c	double precision L
	common/chain/L
	common/wavenum/si,alpha
      double precision lamdap
      data pi/3.1415926535897930d0/

      integer i,j,k
	
      read(10,*) Br
      read(10,*) T0
      read(10,*) T1
      read(10,*) T2
      read(10,*) delta
      read(10,*) r21
      read(10,*) betas
      read(10,*) delhs
	read(10,*) si
	read(10,*) alpha
	read(10,*) Re
	read(10,*) epsilons
	read(10,*) Pe

c	read(10,*) uin
c	read(10,*) uout
c	read(10,*) rot_p
c	read(10,*) visc
c	read(10,*) dens
c	read(10,*) cond 
c	read(10,*) gammadot

	if(nmodes.NE.0) then
	  do i=1,nmodes
	    read(9,*) alp(i)
	    read(9,*) epsilonp(i)
	    read(9,*) epsilonlam(i)
	    read(9,*) betap(i)
	    read(9,*) de(i)
	  enddo
	endif

	
c     mode parameters
c     do i=1,nmodes
c       de(i)=50.32d0
c       alp(i)=0.d0
c	  delhlam(i)=61210.d0
c	  delhp(i)=61210.d0
c 	  epsilonp(i)=delhp(i)/(8.314d0*T0)
c       epsilonlam(i)=delhlam(i)/(8.314d0*T0) 
c	  lamdap(i)=0.d0
c       betap(i)=1-betas
         
c     enddo

	close(unit=10)
	close(unit=9)


	return
      end

c       -------------END OF SUBROUTINE INPUTVALUE----------

c       --------------SUBROUTINE MESHPOINT------------------

      subroutine meshpoint(cheb,rpt,modes,ncheb,dr1)
	implicit double precision (a-h,o-z)
      integer modes,ncheb,i,j,k
      double precision cheb(ncheb),rpt(ncheb),dr1

      data pi/3.1415926535897930d0/


      do 1 j=1,ncheb

	  cheb(j) = dcos(pi*dfloat(j-1)/dfloat(ncheb-1))
1     continue

      do 2 j=1,ncheb
        rpt(j) =(1.d0-cheb(j))/2.d0+(1.d0/dr1)
2     continue
      
	
	return
      end
c       -----------END OF SUBROUTINE MESHPOINT--------------

c       -----------SUBROUTINE DERIVATIVE---------------

      subroutine derivative(dmatrix,cheb,ncheb)
	implicit double precision (a-h,o-z)
      integer i,j,ncheb
      double precision c(ncheb),dmatrix(ncheb,ncheb),cheb(ncheb)

      data pi/3.1415926535897930d0/

      do  i =1,ncheb
        if((i.EQ.1).or.(i.EQ.ncheb)) then
          c(i)=2.d0
        else
  	    c(i)=1.d0
        end if
      end do

      do i=1,ncheb
        do j=1,ncheb
          if (i.NE.j) then
            dmatrix(i,j)=c(i)*(-1.d0)**(i+j)/(c(j)*
     &                   (cheb(i)-cheb(j)))

          elseif((i.GE.2).and.(i.EQ.j).and.(j.LE.(ncheb-1))) then
            dmatrix(i,j)= -1.d0*cheb(j)/(2.d0*(1.d0-cheb(j)**2))

          elseif ((i.EQ.1).and.(j.EQ.1)) then
            dmatrix(i,j)=(2.d0*dfloat((ncheb-1)**2)+1.d0)/6.d0

          elseif ((i.EQ.ncheb).and.(j.EQ.ncheb)) then
            dmatrix(i,j)=(-1.d0)*(2.d0*dfloat((ncheb-1)**2)
     &           +1.d0)/6.d0
          endif
        enddo
      enddo


	return
      end
c       ----------------END OF SUBROUTINE DERIVATIVE-----------------

c       --------------subroutine square matrix----------------------

      subroutine squarematrix(mat,ncheb,sqmat)
	implicit double precision (a-h,o-z)
      integer i,j,k,ncheb
      double precision sum,mat(ncheb,ncheb),sqmat(ncheb,ncheb)
      data pi/3.1415926535897930d0/


      do 4 i=1,ncheb
        do 5  k=1,ncheb
          sum=0.d0
          do 6 j=1,ncheb
            sum=sum+mat(i,j)*mat(j,k)
6         continue
          
		sqmat(i,k)=sum
5       continue
4     continue

c     do i=1,ncheb

c       do j=1,ncheb
c         write (*,*) mat(i,j), sqmat(i,j)
c       enddo
c     enddo
      
	
	return
      end

c       ---------------end of squarematrix subroutine------------------
c	---------SUBROUTINE STEADYSTATE ------------------------

	subroutine steadystate

      implicit double precision (a-h,o-z)
      parameter(iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
      parameter(nvars=(6*nmodes+5)*n,nrhs=1,nn=2*nvars)
      common/variables/ eigvar(nvars),deigvar(nvars),
     &    ddeigvar(nvars),y(n),r(n),der(n,n),dder(n,n),dfirst(n,n),
     &	dsecond(n,n)
      common/steadyvar/uthss(n),temss(n),trrss(n,nrows),
     & 	trthss(n,nrows), tththss(n,nrows),duthss(n),dtemss(n),
     &	dduthss(n),ddtemss(n),dtrrss(n,nrows),ddtrrss(n,nrows),
     &	dtrthss(n,nrows),ddtrthss(n,nrows),dtththss(n,nrows),
     &	ddtththss(n,nrows)
      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &  epsilonlam(nrows),alp(nrows),delham(nrows),
     &  delhp(nrows),delhs,epsilons,betas
     	common/physical/visc,dens,cond 
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &  gammadot
c	double precision L
c	common/chain/L
      common/wavenum/si,alpha
	integer i,j,k
      double precision sum(6)
	data pi/3.1415926535897930d0/


	do i=1,n
	  uthss(i)=0.d0
	  temss(i)=0.d0
	  if (nmodes.NE.0) then
	    do k=1,nmodes
	      trrss(i,k)=0.d0
	      trthss(i,k)=0.d0
	      tththss(i,k)=0.d0
	    enddo
	  endif
	enddo

	read(11,*)  (uthss(i),i=1,n)                          ! Uth
      read(11,*)  (temss(i),i=1,n)                          ! Tem

	if (nmodes.NE.0) then
	  do k=1,nmodes
	    read(11,*) (trrss(i,k),i=1,n)                     ! Trr
	  enddo
	  
	  do k=1,nmodes
	    read(11,*) (trthss(i,k),i=1,n)                    ! Trt
	  enddo
	  
	  do k=1,nmodes
	    read(11,*) (tththss(i,k),i=1,n)                   ! Ttt
	  enddo
     	endif 
	close(unit=11)


	do 10 i=1,n
	  do 11 j=1,n
	    dfirst(i,j)=der(i,j)*(-2.d0)
		dsecond(i,j)=dder(i,j)*4.d0
11	  continue
10	continue
	
	
	k=1
	do 12 i=1,n
	  sum(k)=0.d0
	  sum(k+1)=0.d0
	  sum(k+2)=0.d0
	  sum(k+3)=0.d0
	  
	  do 13 j=1,n
	    sum(k)=sum(k)+dfirst(i,j)*uthss(j)
	    sum(k+1)=sum(k+1)+dsecond(i,j)*uthss(j)
	    sum(k+2)=sum(k+2)+dfirst(i,j)*temss(j)
	    sum(k+3)=sum(k+3)+dsecond(i,j)*temss(j)
13	  continue
	  
	  duthss(i)=sum(k)                                   ! dUth/dr
	  dduthss(i)=sum(k+1)                                ! d2Uth/dr2
	  dtemss(i)=sum(k+2)                                 ! dTem/dr
	  ddtemss(i)=sum(k+3)                                ! d2Tem/dr2
12	continue
	
	
	if (nmodes.NE.0) then
	  do 14 l=1,nmodes
	    k=1
	
	    do 15 i=1,n
	      sum(k)=0.d0
	      sum(k+1)=0.d0
	      sum(k+2)=0.d0
	      sum(k+3)=0.d0
	      sum(k+4)=0.d0
	      sum(k+5)=0.d0
	      
		  do 16 j=1,n
	        sum(k)=sum(k)+dfirst(i,j)*trrss(j,l)
	        sum(k+1)=sum(k+1)+dsecond(i,j)*trrss(j,l)
	        sum(k+2)=sum(k+2)+dfirst(i,j)*trthss(j,l)
	        sum(k+3)=sum(k+3)+dsecond(i,j)*trthss(j,l)
	        sum(k+4)=sum(k+4)+dfirst(i,j)*tththss(j,l)
	        sum(k+5)=sum(k+5)+dsecond(i,j)*tththss(j,l)
16	      continue
	      
		  dtrrss(i,l)=sum(k)                               ! dTrr/dr
	      ddtrrss(i,l)=sum(k+1)                            ! d2Trr/dr
	      dtrthss(i,l)=sum(k+2)                            ! dTrt/dr
	      ddtrthss(i,l)=sum(k+3)                           ! d2Trt/dr
	      dtththss(i,l)=sum(k+4)                           ! dTtt/dr
	      ddtththss(i,l)=sum(k+5)                          ! d2Ttt/dr
15	    continue
14	  continue
	endif
	

	return
	end 
c	---------END OF SUBROUTINE STEADYSTATE-----------------		
c	-----------SUBROUTINE CONSTRUCTMATRIX-------------------
	subroutine constructmatA(AA,q)
	
	implicit double precision (a-h,o-z)
      parameter(iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
      parameter(nvars=(6*nmodes+5)*n,nrhs=1,nn=2*nvars)
      common/variables/ eigvar(nvars),deigvar(nvars),
     &    ddeigvar(nvars),y(n),r(n),der(n,n),dder(n,n),dfirst(n,n),
     &	dsecond(n,n)
      common/steadyvar/uthss(n),temss(n),trrss(n,nrows),
     & 	trthss(n,nrows), tththss(n,nrows),duthss(n),dtemss(n),
     &	dduthss(n),ddtemss(n),dtrrss(n,nrows),ddtrrss(n,nrows),
     &	dtrthss(n,nrows),ddtrthss(n,nrows),dtththss(n,nrows),
     &	ddtththss(n,nrows)
      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &    epsilonlam(nrows),alp(nrows),delham(nrows),
     &    delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &    gammadot
c	double precision L
	common/chain/L
     	common/wavenum/si,alpha
   	complex*16 AA(nvars,nvars), BB(nvars,nvars), csi,ca
	double precision fun1(n),fun2(n),funp(n,nrows),
     &	funlam(n,nrows),funp2(n,nrows),funli(n,nrows),funp3(n,nrows),
     &	funp4(n,nrows)
	double precision strrss(n),strthss(n),stththss(n),sum(3),q(10)  
      integer i,j,k
   
	data pi/3.1415926535897930d0/
c       data L/1000.d0/
c      f=(L*L-3.d0)*(1/(L*L-trrss(i,k)-tththss(i,k)))

c      write(*,*) 'peterlin function is' , f;
c	Initializing matrix AA

	do i=1,nvars
	  do j=1,nvars
	    AA(i,j)=dcmplx(0.d0,0.d0)
	  enddo
	enddo
	
	csi=dcmplx(0.d0,si)
	ca=dcmplx(0.d0,alpha)
	
	do i=1,n

	  fun1(i)=dexp(epsilons*((1.d0/temss(i))-1.d0))
	  fun2(i)=-epsilons/(temss(i)**2)
c	  fun2(i)=0.d0
	  strrss(i)=0.d0
	  stththss(i)=0.d0
	  strthss(i)=0.d0
	enddo
	
	if (nmodes.NE.0) then
	  do i=1,n
	    sum(1)=0.d0
	    sum(2)=0.d0
	    sum(3)=0.d0
	    
		do k=1,nmodes
	      sum(1)=sum(1)+trrss(i,k)*betap(k)
	      sum(2)=sum(2)+tththss(i,k)*betap(k)
	      sum(3)=sum(3)+trthss(i,k)*betap(k)
	    enddo
	    
		strrss(i)=sum(1)
	    stththss(i)=sum(2)
	    strthss(i)=sum(3)
	  enddo
	endif


	do 16 i=1,n
	  do 17 j=1,n
c           write(*,*) ' deborah is', de(k)
	    AA(i,j+n)=AA(i,j+n)+dfirst(i,j)                                        ! for Continuity, dur/dr
	

	    AA(i+n,j)=AA(i+n,j)-dfirst(i,j)	                                       ! for ur, -dp/dr
	    
		AA(i+n,j+n)=AA(i+n,j+n)+betas*fun1(i)*(dsecond(i,j)+(2.d0*             ! for ur, betas*exp(eps_s*(1/Tem-1))*(d2ur/dr2+
     &	            dfirst(i,j)*fun2(i)*dtemss(i))+(dfirst(i,j)/r(i)))         !            2*(-eps_s/Tem**2)*dTem/dr*dur/dr+dur/dr/r) 


	    AA(i+2*n,j+2*n)=AA(i+2*n,j+2*n)+betas*fun1(i)*((fun2(i)*               ! for ut, betas*exp(eps_s*(1/Tem-1))*(
     &	     dtemss(i)*dfirst(i,j))+dsecond(i,j)+(dfirst(i,j)/r(i)))           !            (-eps_s/Tem**2)*dTem/dr*dut/dr+d2ut/dr2+dut/dr/r)
          		
	    AA(i+2*n,j+4*n)=AA(i+2*n,j+4*n)+(duthss(i)-uthss(i)/r(i))              ! for ut, betas*exp(eps_s*(1/Tem-1))*(-eps_s/Tem**2)*
     &	               *fun2(i)*dfirst(i,j)*betas*fun1(i)                      !            (dUth/dr-Uth/r)
          	
	    AA(i+3*n,j+3*n)=AA(i+3*n,j+3*n)+betas*fun1(i)*((fun2(i)*               ! for uz, betas*exp(eps_s*(1/Tem-1))*(
     &         dtemss(i)*dfirst(i,j))+(dfirst(i,j)/r(i))+(dsecond(i,j)))         !            (-eps_s/Tem**2)*dTem/dr*duz/dr+duz/dr+d2uz/dr2)

c		AA(i+4*n,j+n)=AA(i+4*n,j+n)+Br*strrss(i)*dfirst(i,j)                   ! for tem, Br*(Sigma betap*Trr)*dur/dr
c	    	  	    	write(*,*) 'AA(i+4*n,j+n)is =', 
c     &strrss(i)/de        
	    
c		AA(i+4*n,j+2*n)=AA(i+4*n,j+2*n)+2.d0*Br*betas*fun1(i)*(                ! for tem, 2*Br*betas*exp(eps_s*(1/Tem-1))*(dUth/dr-Uth/r)*dut/dr
c     &	duthss(i)-uthss(i)/r(i))*dfirst(i,j)+Br*strthss(i)*dfirst(i,j)         !              +Br*(Sigma betap*Trt)*dut/dr
	  
c         AA(i+4*n,j+4*n)=AA(i+4*n,j+4*n)+(1.d0/r(i))*dfirst(i,j)+                  ! for tem, dtem/dr/r+d2tem/dr2 
c     &	dsecond(i,j)                   
        

	    if(nmodes.NE.0) then 			
	      do k=1,nmodes

	    AA(i+n,j+4*n+n*k)=AA(i+n,j+4*n+n*k)+betap(k)*                          ! for ur, ( betap*dcrr/dr )_pi
     &    ((L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
     &    *trrss(i,k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-
     &	2.d0*trrss(i,k))))*(1.d0/de(k))*dfirst(i,j)

                                                                															   
          AA(i+n,j+4*n+3*n*nmodes+n*k)=AA(i+n,j+4*n+3*n*nmodes+n*k)              ! for ur, ( betap*dctt/dr )_pi
     &	+betap(k)*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)
     &    )**2)*trrss(i,k)*(1.d0/de(k))*dfirst(i,j)


	    AA(i+n,j+4*n+5*n*nmodes+n*k)=AA(i+n,j+4*n+5*n*nmodes+n*k)              ! for ur, ( betap*dczz/dr )_pi
     &	+betap(k)*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0
     &    *trrss(i,k))**2)*trrss(i,k)*(1.d0/de(k))*dfirst(i,j)

            
c	        AA(i+n,j+4*n+n*k)=AA(i+n,j+4*n+n*k)+dfirst(i,j)*betap(k)           ! for ur, ( betap*dtrr/dr )_pi

          AA(i+2*n,j+4*n+n*k)=AA(i+2*n,j+4*n+n*k)+betap(k)*(1.d0/de(k))          ! for ut, ( betap*dcrr/dr )_pi
     &	*(L*L-3.d0)*(1.0d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
     &	*trthss(i,k)*dfirst(i,j)

 
             
	    AA(i+2*n,j+4*n+n*nmodes+n*k)=AA(i+2*n,j+4*n+n*nmodes+n*k)+             ! for ut, ( betap*dcrt/dr )_pi
     &	betap(k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
     & 	-2.d0*trrss(i,k)))*dfirst(i,j)


	    AA(i+2*n,j+4*n+3*n*nmodes+n*k)=AA(i+2*n,j+4*n+3*n*nmodes+n             ! for ut, ( betap*dctt/dr )_pi
     &	*k)+betap(k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k) 
     &	-2.d0*trrss(i,k))**2)*trthss(i,k)*dfirst(i,j)

         
     	    AA(i+2*n,j+4*n+5*n*nmodes+n*k)=AA(i+2*n,j+4*n+5*n*nmodes+n             ! for ut, ( betap*dczz/dr )_pi
     &	*k)+betap(k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
     &	-2.d0*trrss(i,k))**2)*trthss(i,k)*dfirst(i,j)  
 
	
                 	        
c			AA(i+2*n,j+4*n+n*nmodes+n*k)=AA(i+2*n,j+4*n+n*nmodes+n*k)+         ! for ut, ( betap*dtrt/dr )_pi
c     &                  betap(k)*dfirst(i,j)

          AA(i+3*n,j+4*n+2*n*nmodes+n*k)=AA(i+3*n,j+4*n+2*n*nmodes+              ! for uz, ( betap*dcrz/dr )_pi
     &	n*k)+betap(k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k
     &    )-2.d0*trrss(i,k)))*dfirst(i,j)

	        
c			AA(i+3*n,j+4*n+2*n*nmodes+n*k)=                                    ! for uz, ( betap*dtrz/dr )_pi
c     &			AA(i+3*n,j+4*n+2*n*nmodes+n*k)+dfirst(i,j)*betap(k)


            
         AA(i+4*n,j+n)=AA(i+4*n,j+n)-Br*(1.d0/de(k))*betap(k)                    ! for tem, ( betap*dur/dr )_pi
     &   *dfirst(i,j)+Br*betap(k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L
     &   -tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)*dfirst(i,j) 

	
         AA(i+4*n,j+2*n)=AA(i+4*n,j+2*n)+Br*(1.d0/de(k))*betap(k)*               ! for tem, ( betap*dut/dr )_pi
     &   (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))  
     &   *trthss(i,k)*dfirst(i,j)+2.d0*Br*betas*fun1(i)
     &   *(duthss(i)-uthss(i)/r(i))*dfirst(i,j)

  	           
		AA(i+4*n,j+4*n)=AA(i+4*n,j+4*n)+dfirst(i,j)/r(i)+dsecond(i,j)          ! for tem, dtem/dr/r+d2tem/dr2



	      enddo
	    endif
17	  continue	
16	continue
 

	do i=1,n
	
	  AA(i,i+n)=AA(i,i+n)+(1.d0/r(i))                                          ! for continuity, ur/r
	    
	  AA(i,i+2*n)=AA(i,i+2*n)+csi/r(i)                                         ! for continuity, i*si*ut/r

	  AA(i,i+3*n)=AA(i,i+3*n)+ca	                                           ! for continuity, i*al*uz

	  AA(i+n,i+n)=AA(i+n,i+n)-Re*uthss(i)*csi/r(i)+betas*(-(si**2)*            ! for ur, -i*si*Re*Uth*ur/r+betas*exp(eps_s*(1/Tem-1))*
     &	    fun1(i)/(r(i)**2)-(fun1(i)/r(i)**2)-(alpha**2)*fun1(i))            !           (-si**2*ur/r**2-ur/r**2-aL*L*ur)

	  AA(i+n,i+2*n)=AA(i+n,i+2*n)+(2.d0*uthss(i)*Re/r(i))-(2.d0*csi*           ! for ur, 2*Re*Uth*/r*ut-i*2*si*betas*exp(eps_s*(1/Tem-1))/r**2*ut
     &	    fun1(i)*betas/r(i)**2)	
	  
	  AA(i+n,i+4*n)=AA(i+n,i+4*n)+betas*fun2(i)*fun1(i)*(duthss(i)-            ! for ur, i*si*betas*exp(eps_s*(1/Tem-1))*(-eps_s/Tem**2)*
     &	    uthss(i)/r(i))*csi/r(i)                                            !            *(dUth/dr-Uth/r)/r*tem

	  AA(i+2*n,i)=AA(i+2*n,i)-csi/r(i)                                         ! for ut, -i*si/r*p
	  
	  AA(i+2*n,i+n)=AA(i+2*n,i+n)-Re*duthss(i)-Re*uthss(i)/r(i)+               ! for ut, -Re*(dUth/dr+Uth/r)*ur+i*2*si*betas*exp(eps_s*(1/Tem-1))/r**2*ur
     &        betas*fun1(i)*2.d0*csi/r(i)**2
     &       +betas*fun1(i)*fun2(i)*dtemss(i)*csi/r(i) 
	  
	  AA(i+2*n,i+2*n)=AA(i+2*n,i+2*n)-Re*csi*uthss(i)/r(i)+betas*              ! for ut, -i*si*Re*Uth/r*ut+betas*exp(eps_s*(1/Tem-1))*(
     &	    fun1(i)*((-(si**2)/(r(i)**2))-(fun2(i)*                            !            -si**2/r**2*ut-ut/r**2-aL*L*ut-(-eps_s/Tem**2)*dTem/dr/r*ut)
     &	    dtemss(i)/r(i))-(1.d0/r(i)**2)-alpha**2)
        
	  AA(i+2*n,i+4*n)=AA(i+2*n,i+4*n)+betas*fun1(i)*((-2.d0*fun2(i)*           ! for ut, betas*exp(eps_s*(1/Tem-1))*(-eps_s/Tem**2)*dTem/dr*
     &	(duthss(i)-uthss(i)/r(i))*dtemss(i)/temss(i))+(fun2(i)**2*             !            *(dUth/dr-Uth/r)(-2/Tem-(-eps_s/Tem**2))*tem
     &	dtemss(i)*(duthss(i)-uthss(i)/r(i))))+betas*fun1(i)*fun2(i)*(          !            +betas*exp(eps_s*(1/Tem-1))*(-eps_s/Tem**2)*
     &	(-uthss(i)/r(i)**2)+(duthss(i)/r(i))+dduthss(i))                       !               *(d2Uth/dr2+dUth/dr/r-Uth/r**2)*tem
	
	  AA(i+3*n,i)=AA(i+3*n,i)-ca                                               ! for uz, -i*al*p
	  
	  AA(i+3*n,i+n)=AA(i+3*n,i+n)+betas*fun1(i)*fun2(i)*                       ! for uz, i*al*betas*exp(eps_s*(1/Tem-1))*(-eps_s/Tem**2)*dTem/dr*ur
     &        dtemss(i)*ca
	
	  AA(i+3*n,i+3*n)=AA(i+3*n,i+3*n)-Re*uthss(i)*csi/r(i)+betas*              ! for uz, -i*si*Re*Uth/r*uz+betas*exp(eps_s*(1/Tem-1))*
     &        fun1(i)*((-si**2)/(r(i)**2)-alpha**2)                              !            (-si**2/r**2-aL*L)*uz


		
c	  AA(i+4*n,i+n)=AA(i+4*n,i+n)+(Br*stththss(i)/r(i))+(Br*strthss(i)         ! for tem, Br*(Sigma betap*Ttt)/r*ur+i*si*Br*(Sigma betap*Trt)/r*ur
c     &	*csi/r(i))-Pe*dtemss(i)+2.d0*csi*Br*betas*fun1(i)*(duthss(i)-          !            -Pe*dTem/dr*ur+i*si*2*Br*betas*exp(eps_s*(1/Tem-1))*
c     &	uthss(i)/r(i))/r(i)                                                    !                                (dUth/dr-Uth/r)/r*ur
	  
c	  AA(i+4*n,i+2*n)=AA(i+4*n,i+2*n)-(2.d0*Br*betas*fun1(i)*                  ! for tem, -2*Br*betas*exp(eps_s*(1/Tem-1))*(dUth/dr-Uth/r)/r*ut
c     &	(duthss(i)-uthss(i)/r(i))/r(i))+(Br*stththss(i)*csi/r(i))              !             +i*si*Br*(Sigma betap*Ttt)/r*ut-Br*(Sigma betap*Trt)/r*ut
c     &	-(Br*strthss(i)/r(i))


        AA(i+4*n,i+4*n)=AA(i+4*n,i+4*n)-Pe*(1.d0/r(i))*uthss(i)*csi              ! for tem, tem
     &  +(1.d0/r(i))*(1.d0/r(i))*(-1.d0)*(si**2)-(alpha**2)+
     &  Br*betas*fun1(i)*fun2(i)*(duthss(i)**2+(-uthss(i)/r(i))**2+
     &  2.d0*duthss(i)*(-uthss(i)/r(i)))

	        
c	  AA(i+4*n,i+4*n)=AA(i+4*n,i+4*n)-(Pe*uthss(i)*csi/r(i))-(si**2/           ! for tem, -i*si*Pe*Uth/r*tem-si**2/r**2*tem-aL*L*tem
c     &	r(i)**2)-(alpha**2)+Br*betas*fun1(i)*fun2(i)*(duthss(i)**2+            !            +Br*betas*exp(eps_s*(1/Tem-1))*(-eps_s/Tem**2)*
c     &	(-uthss(i)/r(i))**2+2.d0*duthss(i)*(-uthss(i)/r(i)))                   !              (dUth/dr-Uth/r)**2*tem
	
	  if (nmodes.NE.0) then
	    do 18 k=1,nmodes
	      
	  AA(i+n,i+4*n+n*k)=AA(i+n,i+4*n+n*k)+((L*L-3.d0)*(1.d0/r(i))*            ! for ur, +( betap/crr )_pi
     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))+(L*L-3.d0)*(1.d0/
     &  r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trrss(i,k)
     &  +(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*
     &  (2.d0*dtrrss(i,k)+dtththss(i,k))+(L*L-3.d0)*(-2.d0)*
     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**3)*(2.d0*dtrrss(i,k)
     &  +dtththss(i,k))*trrss(i,k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-
     &  2.d0*trrss(i,k))**2)*dtrrss(i,k)+(L*L-3.d0)*(1.d0/r(i))*csi
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trthss(i,k)-(L*L
     &  -3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
     &  *tththss(i,k))*betap(k)*(1.d0/de(k))



                                   ! for ur, +( betap/crr )_pi

	
c     &-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
c     &*tththss(i,k)*betap(k)*(1.d0/de(k))
c		  AA(i+n,i+4*n+n*k)=AA(i+n,i+4*n+n*k)+betap(k)/r(i)                    ! for ur, +( betap/r*trr )_pi


        AA(i+n,i+4*n+n*nmodes+n*k)=AA(i+n,i+4*n+n*nmodes+n*k)+                   ! for ur, +(betap*crt )_pi	
     &  (L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)
     &  ))*csi*betap(k)*(1.d0/de(k))


	      
c		  AA(i+n,i+4*n+n*nmodes+n*k)=AA(i+n,i+4*n+n*nmodes+n*k)+csi*           ! for ur, +i*si*( betap/r*trt)_pi
c     &	       betap(k)/r(i)

        AA(i+n,i+4*n+3*n*nmodes+n*k)=AA(i+n,i+4*n+3*n*nmodes+n*k)                ! for ur, +(betap/ctt)_pi
     &  +betap(k)*((L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0
     &  *trrss(i,k))**2)*trrss(i,k)+(L*L-3.d0)*-2.d0*(1.d0/(L*L-tthths
     &  s(i,k)-2.d0*trrss(i,k))**3)*(2.d0*dtrrss(i,k)+dtththss(i,k))*
     &  trrss(i,k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))
     &  **2)*dtrrss(i,k)+(L*L-3.d0)*(1.d0/r(i))*csi*(1.d0/(L*L-tththss
     &  (i,k)-2.d0*trrss(i,k))**2)*trthss(i,k)-(L*L-3.d0)*(1.d0/r(i))
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))-(L*L-3.d0)
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*tththss(i,k)*
     &  (1.d0/r(i)))*(1.d0/de(k))


	 
c		  AA(i+n,i+4*n+2*n*nmodes+n*k)=AA(i+n,i+4*n+2*n*nmodes+n*k)+           ! for ur, i*al*( betap*trz )_pi
c     &	       betap(k)*ca

        AA(i+n,i+4*n+2*n*nmodes+n*k)=AA(i+n,i+4*n+2*n*nmodes+n*k)+               ! for ur, +( betap*czr)_pi
     &  betap(k)*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &  *ca*(1.d0/de(k))


c	      AA(i+n,i+4*n+3*n*nmodes+n*k)=AA(i+n,i+4*n+3*n*nmodes+n*k)-           ! for ur, -( betap/r*ttt )_pi
c     &           betap(k)/r(i)
        
        AA(i+n,i+4*n+5*n*nmodes+n*k)=AA(i+n,i+4*n+5*n*nmodes+n*k)                ! for ur, +( betap*czz)_pi
     &  +betap(k)*((L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0
     &  *trrss(i,k))**2)*trrss(i,k)+(L*L-3.d0)*-2.d0*(1.d0/(L*L-
     &  tththss(i,k)-2.d0*trrss(i,k))**3)*(2.d0*dtrrss(i,k)+
     &  dtththss(i,k))*trrss(i,k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
     &  -2.d0*trrss(i,k))**2)*dtrrss(i,k)-(L*L-3.d0)*(1.d0/r(i))
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*tththss(i,k)+
     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
     &  *(1.d0/r(i))*csi*trthss(i,k))*(1.d0/de(k))
  


	  AA(i+2*n,i+4*n+n*k)=AA(i+2*n,i+4*n+n*k)+betap(k)*(1.d0/de(k))           ! for ut, +( betap*crr)_pi
     &  *((L*L-3.d0)*(-2.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))*
     &  *3)*(2.d0*dtrrss(i,k)+dtththss(i,k))*trthss(i,k)+(L*L-3.d0)*
     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*dtrthss(i,k)+2.d0
     &  *(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0
     &  *trrss(i,k))**2)*trthss(i,k)+(L*L-3.d0)*(1.d0/r(i))*(1.d0/
     &  (L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*csi*tththss(i,k))



c	      AA(i+2*n,i+4*n+n*nmodes+n*k)=AA(i+2*n,i+4*n+n*nmodes+n*k)+           ! for ut, +2*( betap/r*trt )_pi
c     &	       2.d0*betap(k)/r(i)

        AA(i+2*n,i+4*n+n*nmodes+n*k)=AA(i+2*n,i+4*n+n*nmodes+n*k)                ! for ut, +( betap*crt )_pi
     &  +betap(k)*(1.d0/de(k))*(L*L-3.d0)*(2.d0/r(i))*(1.d0/(L*L-
     &  tththss(i,k)-2.d0*trrss(i,k)))

	      
c		  AA(i+2*n,i+4*n+3*n*nmodes+n*k)=                                      ! for ut, +i*si*( betap/r*ttt)_pi
c    &	       AA(i+2*n,i+4*n+3*n*nmodes+n*k)+betap(k)*csi/r(i)

        AA(i+2*n,i+4*n+3*n*nmodes+n*k)=AA(i+2*n,i+4*n+3*n*nmodes+n*k             ! for ut, +( betap*ctt )_pi
     &  )+betap(k)*(1.d0/de(k))*((L*L-3.d0)*(-2.d0)*(1.d0/(L*L-tththss
     &  (i,k)-2.d0*trrss(i,k))**3)*(2.d0*dtrrss(i,k)+dtththss(i,k))
     &  *trthss(i,k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
     &  trrss(i,k))**2)*dtrthss(i,k)+2.d0*(L*L-3.d0)*(1.d0/r(i))*(1.d0/
     &  (L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trthss(i,k)+(L*L-3.d0)
     &  *(1.d0/r(i))*csi*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
     &  *tththss(i,k)+(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-
     &  2.d0*trrss(i,k)))*csi)


	      
c		  AA(i+2*n,i+4*n+4*n*nmodes+n*k)=                                      ! for ut, +i*al*( betap*ttz)_pi
c     &	       AA(i+2*n,i+4*n+4*n*nmodes+n*k)+ca*betap(k)

        AA(i+2*n,i+4*n+4*n*nmodes+n*k)=AA(i+2*n,i+4*n+4*n*nmodes+n*k            ! for ut, +( betap*czt )_pi
     &  )+betap(k)*(1.d0/de(k))*((L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-
     &  2.d0*trrss(i,k)))*ca)


 
        AA(i+2*n,i+4*n+5*n*nmodes+n*k)=AA(i+2*n,i+4*n+5*n*nmodes+n*k            ! for ut, +( betap*czz )_pi
     &  )+betap(k)*(1.d0/de(k))*((L*L-3.d0)*(-2.d0)*(1.d0/(L*L-
     &  tththss(i,k)-2.d0*trrss(i,k))**3)*(2.d0*dtrrss(i,k)+dtththss
     &  (i,k))*trthss(i,k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0
     &  *trrss(i,k))**2)*dtrthss(i,k)+2.d0*(L*L-3.d0)*(1.d0/r(i))
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trthss(i,k)+
     &  (L*L-3.d0)*(1.d0/r(i))*csi*(1.d0/(L*L-tththss(i,k)-2.d0
     &  *trrss(i,k))**2)*tththss(i,k))

c	      AA(i+3*n,i+4*n+n*k)=AA(i+3*n,i+4*n+n*k)+betap(k)*(1/de(k))*
c	      ((L*L-3)*(1/(L*L-tththss(i,k)-2*trrss(i,k))**3)*(


        AA(i+3*n,i+4*n+n*k)=AA(i+3*n,i+4*n+n*k)+betap(k)*(1.d0/de(k))           ! for uz, +( betap*crr)_pi
     &  *(L*L-3.d0)*ca*trrss(i,k)*(1.d0/(L*L-tththss(i,k)
     &  -2.d0*trrss(i,k))**2)  
     
      
	       
	  AA(i+3*n,i+4*n+3*n*nmodes+n*k)=AA(i+3*n,i+4*n+3*n*nmodes+n*k)           ! for uz, +( betap*ctt)_pi
     &  +betap(k)*(1.d0/de(k))*(L*L-3.d0)*ca*trrss(i,k)
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)


     
	  AA(i+3*n,i+4*n+5*n*nmodes+n*k)=AA(i+3*n,i+4*n+5*n*nmodes+n*k)           ! for uz, +( betap*czz)_pi
     &  +betap(k)*(1.d0/de(k))*((L*L-3.d0)*ca*trrss(i,k)
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)+(L*L-3.d0)*ca
     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))))

	  AA(i+3*n,i+4*n+2*n*nmodes+n*k)=AA(i+3*n,i+4*n+2*n*nmodes+n              ! for uz, +( betap*crz)_pi
     &  *k)+betap(k)*(1.d0/de(k))*((L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-
     &  2.d0*trrss(i,k))**2)*(2.d0*dtrrss(i,k)+dtththss(i,k))+(L*L-3.d0
     &  )*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))))

 	 
c	      AA(i+3*n,i+4*n+2*n*nmodes+n*k)=                                      ! for uz, +( betap/r*trz)_pi
c     &	       AA(i+3*n,i+4*n+2*n*nmodes+n*k)+betap(k)/r(i)

        AA(i+3*n,i+4*n+4*n*nmodes+n*k)=AA(i+3*n,i+4*n+4*n*nmodes+n*             ! for uz, +( betap*ctz)_pi
     &  k)+betap(k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/r(i))*csi*
     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))

 
c		  AA(i+3*n,i+4*n+4*n*nmodes+n*k)=                                      ! for uz, +i*si*( betap/r*ttz)_pi
c     &	       AA(i+3*n,i+4*n+4*n*nmodes+n*k)+betap(k)*csi/r(i)

	      
c		  AA(i+3*n,i+4*n+5*n*nmodes+n*k)=                                      ! for uz, +i*al*( betap*tzz)_pi


        AA(i+4*n,i+n)=AA(i+4*n,i+n)-Pe*dtemss(i)+2.d0*Br*betas*fun1(i)*          ! for tem, ur
     & (-uthss(i)/r(i)+duthss(i))*csi*(1.d0/r(i))+Br*(1.d0/de(k))
     & *betap(k)*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     & *trthss(i,k)*(1.d0/r(i))*csi+Br*(1.d0/de(k))*betap(k)*(L*L-3.d0)
     & *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*tththss(i,k)
     & *(1.d0/r(i))-Br*(1.d0/de(k))*betap(k)*(1.d0/r(i))
           

	  AA(i+4*n,i+2*n)=AA(i+4*n,i+2*n)+2.d0*Br*betas*fun1(i)*                   ! for tem, ut
     &  (-uthss(i)/r(i)+duthss(i))*(-1.d0/r(i))+Br*(1.d0/de(k))
     &  *(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*
     &  trthss(i,k)*betap(k)*(-1.d0/r(i))+Br*(1.d0/de(k))*betap(k)
     &  *(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*
     &  tththss(i,k)*(1.d0/r(i))*csi-Br*(1.d0/de(k))*betap(k)
     &  *(1.d0/r(i))*csi

c		      	write(*,*) ' AA(i+4*n,i+2*n)is =', 
   
c     &AA(i+4*n,i+2*n)
    

	  AA(i+4*n,i+3*n)=AA(i+4*n,i+3*n)+Br*(1.d0/de(k))*betap(k)*                ! for tem, uz
     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &  *trrss(i,k)*ca-Br*(1.d0/de(k))*betap(k)*ca
     


c     &	       AA(i+3*n,i+4*n+5*n*nmodes+n*k)+ca*betap(k)


        AA(i+4*n,i+4*n+n*k)=AA(i+4*n,i+4*n+n*k)+Br*(1.d0/de(k))*                 ! for tem, crr
     &  betap(k)*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-
     &  2.d0*trrss(i,k))**2)*trthss(i,k)*(-uthss(i)/r(i)+duthss(i))

	  AA(i+4*n,i+4*n+n*nmodes+n*k)=AA(i+4*n,i+4*n+n*nmodes+n*k)                ! for tem, crt
     &  +Br*betap(k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-
     &  2.d0*trrss(i,k)))*(-uthss(i)/r(i)+duthss(i))

	



	  AA(i+4*n,i+4*n+3*n*nmodes+n*k)=                                          ! for tem, ctt
     &  AA(i+4*n,i+4*n+3*n*nmodes+n*k)+Br*(1.d0/de(k))*betap(k)*
     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
     &  *trthss(i,k)*(-uthss(i)/r(i)+duthss(i))




        AA(i+4*n,i+4*n+5*n*nmodes+n*k)=                                          ! for tem, czz
     &  AA(i+4*n,i+4*n+5*n*nmodes+n*k)+Br*(1.d0/de(k))*betap(k)*
     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
     &  *trthss(i,k)*(-uthss(i)/r(i)+duthss(i)) 
   			

	
c	      AA(i+4*n,i+4*n+n*nmodes+n*k)=AA(i+4*n,i+4*n+n*nmodes+n*k)+           ! for tem, +Br*( betap*(dUth/dr-Uth/r)*trt )_pi
c     &	       Br*betap(k)*(duthss(i)-uthss(i)/r(i))

18	    continue
	  endif
	enddo
	

	if (nmodes.NE.0) then
	  do k=1,nmodes
	    do i=1,n
	      funlam(i,k)=dexp(-epsilonlam(k)*(1.d0/temss(i)-1.d0))                ! exp(-eps_li*(1/Tem-1))
	      funli(i,k)=epsilonlam(k)/temss(i)**2                                 ! eps_li/Tem**2
c	      funli(i,k)=0.d0
	      funp(i,k)=dexp(-epsilonp(k)*(1.d0/temss(i)-1.d0))                    ! exp(-eps_pi*(1/Tem-1))
	      funp2(i,k)=epsilonp(k)/temss(i)**2                                   ! eps_pi/Tem**2
c	      funp2(i,k)=0.d0
	      funp3(i,k)=dexp(epsilonp(k)*(1.d0/temss(i)-1.d0))                    ! exp(eps_pi*(1/Temp-1))
	      funp4(i,k)=-epsilonp(k)/temss(i)**2                                  ! eps_pi/Tem**2
c	      funp4(i,k)=0.d0
	    enddo
	  enddo
	endif
	
	
	if(nmodes.NE.0) then	
	  do k=1,nmodes
	    do 19 i=1,n
	      do 20 j=1,n


	        AA(i+4*n+n*k,j+n)=AA(i+4*n+n*k,j+n)+2.d0*trrss(i,k)                      ! for crr, dur/dr
     & *dfirst(i,j) 

   

c        AA(i+4*n+n*k,j+n)=AA(i+4*n+n*k,j+n)+2.d0*((L*L-3.d0)*                   ! for crr, dur/dr
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)-1.d0)
c     &  *dfirst(i,j)+2.d0*funlam(i,k)*funp3(i,k)*dfirst(i,j)



c	        AA(i+4*n+n*k,j+n)=AA(i+4*n+n*k,j+n)+de(k)*                         ! for trr, 2*De*Trr*dur/dr+2*exp(-eps_li*(1/Tem-1))
c     &	      (2.d0*trrss(i,k)*dfirst(i,j))+2.d0*funlam(i,k)*                                 *exp(eps_pi*(1/Tem-1))*dur/dr
c     &	      funp3(i,k)*dfirst(i,j)

        AA(i+4*n+n*nmodes+n*k,j+n)=AA(i+4*n+n*nmodes+n*k,j+n)+                     ! for crt, dur/dr
     &  trthss(i,k)*dfirst(i,j)

c        AA(i+4*n+n*nmodes+n*k,j+n)=AA(i+4*n+n*nmodes+n*k,j+n)+                     ! for crt, dur/dr
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *trthss(i,k)*dfirst(i,j)


	 	  

c	        AA(i+4*n+n*nmodes+n*k,j+n)=AA(i+4*n+n*nmodes+n*k,j+n)              ! for trt, De*Trt*dur/dr
c     &	      +de(k)*(trthss(i,k)*dfirst(i,j))


         
        AA(i+4*n+n*nmodes+n*k,j+2*n)=AA(i+4*n+n*nmodes+n*k,j+2*n)+               ! for crt, dut/dr
     &  trrss(i,k)*dfirst(i,j)
 

c        AA(i+4*n+n*nmodes+n*k,j+2*n)=AA(i+4*n+n*nmodes+n*k,j+2*n)+               ! for crt, dut/dr
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *trrss(i,k)*dfirst(i,j)-1.d0*dfirst(i,j)+funlam(i,k)*funp3(i,k)*
c     &  dfirst(i,j)



c			AA(i+4*n+n*nmodes+n*k,j+2*n)=AA(i+4*n+n*nmodes+n*k,j+2*n)          ! for trt, De*Trr*dut/dr+
c     &	       +de(k)*trrss(i,k)*dfirst(i,j)                                   !          exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Tem-1))*dut/dr
c     &           +funlam(i,k)*funp3(i,k)*dfirst(i,j)

         AA(i+4*n+2*n*nmodes+n*k,j+3*n)=                                          ! for crz, duz/dr
     &  AA(i+4*n+2*n*nmodes+n*k,j+3*n)+trrss(i,k)*dfirst(i,j)


c        AA(i+4*n+2*n*nmodes+n*k,j+3*n)=                                          ! for crz, duz/dr
c     &  AA(i+4*n+2*n*nmodes+n*k,j+3*n)+((L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)-1.d0+
c     &  funlam(i,k)*funp3(i,k))*dfirst(i,j)

c	        AA(i+4*n+2*n*nmodes+n*k,j+3*n)=                                    ! for trz, De*Trr*duz/dr+
c     &	       AA(i+4*n+2*n*nmodes+n*k,j+3*n)+                                 !          exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Tem-1))*duz/dr
c     &           de(k)*trrss(i,k)*dfirst(i,j)+
c     &           funlam(i,k)*funp3(i,k)*dfirst(i,j)

         AA(i+4*n+3*n*nmodes+n*k,j+2*n)=                                          ! for ctt, dut/dr
     &  AA(i+4*n+3*n*nmodes+n*k,j+2*n)+2.d0*trthss(i,k)*dfirst(i,j)


c        AA(i+4*n+3*n*nmodes+n*k,j+2*n)=                                          ! for ctt, dut/dr
c     &  AA(i+4*n+3*n*nmodes+n*k,j+2*n)+2.d0*((L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k))*dfirst(i,j)


        	  AA(i+4*n+4*n*nmodes+n*k,j+3*n)=                                          ! for ctz, duz/dr
     &  AA(i+4*n+4*n*nmodes+n*k,j+3*n)+trthss(i,k)*dfirst(i,j)


c	  AA(i+4*n+4*n*nmodes+n*k,j+3*n)=                                          ! for ctz, duz/dr
c     &  AA(i+4*n+4*n*nmodes+n*k,j+3*n)+(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)*dfirst(i,j)


c	        AA(i+4*n+3*n*nmodes+n*k,j+2*n)=                                    ! for ttt, 2*De*Trt*dut/dr
c     &		AA(i+4*n+3*n*nmodes+n*k,j+2*n)
c     &	    +2.d0*de(k)*trthss(i,k)*dfirst(i,j)
	
c	        AA(i+4*n+4*n*nmodes+n*k,j+3*n)=                                    ! for ttz, De*Trt*duz/dr
c     &		AA(i+4*n+4*n*nmodes+n*k,j+3*n)
c     &	    +de(k)*trthss(i,k)*dfirst(i,j)
	
20	      continue
19	    continue
	  enddo
	endif


	if(nmodes.NE.0) then
	  do k=1,nmodes
	    do i=1,n

c       iso
c	AA(i+4*n+n*k,i+n)=AA(i+4*n+n*k,i+n)-dtrrss(i,k)                            ! for crr, ur
c     &+2.d0*trthss(i,k)*(1.d0/r(i))*csi-(2.d0*dtrrss(i,k)+dtththss(i,k))
c     &*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*
c     &(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)

c       noniso

	AA(i+4*n+n*k,i+n)=AA(i+4*n+n*k,i+n)-dtrrss(i,k)                            ! for crr, ur
     &+2.d0*trthss(i,k)*(1.d0/r(i))*csi-(2.d0*dtrrss(i,k)+dtththss(i,k))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)
     &+dtemss(i)*(1.d0/temss(i))*trrss(i,k)-(1.d0/(L*L-3.d0))*
     &(L*L-tththss(i,k)-2.d0*trrss(i,k))*dtemss(i)*(1.d0/temss(i))


c        AA(i+4*n+n*k,i+n)=AA(i+4*n+n*k,i+n)-(L*L-3.d0)*(1.d0/(L*L-              ! for crr, ur
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)*(dtththss(i,k)+2.d0*dtrrss                   
c     &  (i,k))*trrss(i,k)-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k)))*dtrrss(i,k)+2.d0*(1.d0/r(i))*(L*L-3.d0)
c     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)*csi
c     &  +(L*L-3.d0)*(dtemss(i)/temss(i))*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k)))*trrss(i,k)-dtemss(i)/temss(i)

c	      AA(i+4*n+n*k,i+n)=AA(i+4*n+n*k,i+n)+de(k)*                          ! for trr, -De*(dTrr/dr*ur+i*si*2*Trt/r*ur
c     &	    (-dtrrss(i,k)+2.d0*trthss(i,k)*csi/r(i)                           !           +dTem/dr*Trr/Tem*ur)
c     &	    +dtemss(i)*trrss(i,k)/temss(i))
	       
c      iso
c       AA(i+4*n+n*k,i+4*n)=AA(i+4*n+n*k,i+4*n)+(L*L-3.d0)*(1.d0/r(i))*          ! for crr, tem
c     &  uthss(i)*csi*(1.d0/temss(i))*(1.d0/(L*L-tththss(i,k)-2.d0 
c     &  *trrss(i,k)))*trrss(i,k)-(1.d0/r(i))*uthss(i)*csi*(1.d0/temss(i)
c     &  )

c      noniso
       AA(i+4*n+n*k,i+4*n)=AA(i+4*n+n*k,i+4*n)+(1.d0/r(i))*uthss(i)              ! for crr, tem
     &*csi*(1.d0/temss(i))*trrss(i,k)-(1.d0/(L*L-3.d0))*
     &(L*L-2.d0*trrss(i,k)-tththss(i,k))*(1/r(i))
     &*uthss(i)*csi*(1.0/temss(i))


c     +funli(i,k)*((L*L-3.d0)*(1.d0/r(i))*uthss(i)*(1.d0/(L*L-
c     &  tththss(i,k)-2*trrss(i,k)))*dtrrss(i,k)-(L*L-3.d0)*2.d0*(1.d0/
c     &  r(i))*uthss(i)*(1.d0/(L*L-tththss(i,k)-2*trrss(i,k)))*trthss
c     &  (i,k)+2*(L*L-3.d0)*(1.d0/r(i))*uthss(i)*(1.d0/(L*L-tththss(i,k
c     &  )-2.d0*trrss(i,k)))*trthss(i,k))
	      


c		  AA(i+4*n+n*k,i+4*n)=AA(i+4*n+n*k,i+4*n)+de(k)*                      ! for trr, De*(i*si*Uth*Trr/r*tem-exp(-eps_pi*(1/Tem-1))*eps_pi/Tem**2*
c     &	    ((csi*uthss(i)*trrss(i,k)/temss(i)/r(i))-funp2(i,k)*              !              (Trr**2+Trt**2)*alpha*tem)
c     &	    (trrss(i,k)**2+trthss(i,k)**2)*alp(k)*funp(i,k))-                 !            -exp(-eps_li*(1/Tem-1))*eps_li/Tem**2*Trr*tem
c     &	    funlam(i,k)*funli(i,k)*trrss(i,k)

c      iso
c        AA(i+4*n+n*k,i+4*n+n*k)=AA(i+4*n+n*k,i+4*n+n*k)                         ! for crr, crr
c     &-(1.d0/r(i))*uthss(i)*csi-(1.d0/de(k))*trrss(i,k)*(L*L)*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &-(1.d0/de(k))*(L*L)*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c      noniso

        AA(i+4*n+n*k,i+4*n+n*k)=AA(i+4*n+n*k,i+4*n+n*k)                         ! for crr, crr
     &-(1.d0/r(i))*uthss(i)*csi-
     &funlam(i,k)*((1.d0/de(k))*trrss(i,k)*(L*L)*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &+(1.d0/de(k))*(L*L)*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k))))
     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c        AA(i+4*n+n*k,i+4*n+n*k)=AA(i+4*n+n*k,i+4*n+n*k)-funlam(i,k)*             ! for crr, crr
c     &  ((1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss
c     &  (i,k)))-(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0
c     &  *trrss(i,k))**2)*trrss(i,k))-(L*L-3.d0)*(1.d0/r(i))*(1.d0
c     &  /(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*uthss(i)*trrss(i,k)*csi
c     &  -(L*L-3.d0)*(1.d0/r(i))*uthss(i)*csi*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k)))

     	  
c            +2*(1/r(i))*(L*L-3)*uthss(i)*trthss(i,k)*(1/(L*L-
c     &	  tththss(i,k)-2*trrss(i,k))**2)
     
c            -2*(L*L-3)*(1/r(i))*uthss(i)
c     &	  *trthss(i,k)*(1/(L*L-tththss(i,k)-2*trrss(i,k))**2)
	      
c		  AA(i+4*n+n*k,i+4*n+n*k)=AA(i+4*n+n*k,i+4*n+n*k)-funlam(i,k)         ! for trr, -exp(-eps_li*(1/Tem-1))*trr-
c     &	    -(uthss(i)*csi*de(k)/r(i))                                        !             i*si*Uth*De/r*trr-exp(-eps_pi*(1/Tem-1))*alpha*2*Trr*De*trr
c     &	    -(funp(i,k)*alp(k)*2.d0*trrss(i,k)*
c     &	    de(k))

c        AA(i+4*n+n*k,i+4*n+n*nmodes+n*k)=AA(i+4*n+n*k,i+4*n+n*nmodes            ! for crr, crt
c     &  +n*k)+2.d0*(L*L-3.d0)*uthss(i)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k)))


c	  AA(i+4*n+n*k,i+4*n+2*n*nmodes+n*k)=                                     ! for crr, crz
c     &  AA(i+4*n+n*k,i+4*n+2*n*nmodes+n*k)-2.d0*(L*L-3.d0)*(1.d0/r(i))*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*uthss(i)

c      iso
c        AA(i+4*n+n*k,i+4*n+3*n*nmodes+n*k)=AA(i+4*n+n*k,i+4*n                   ! for crr, ctt
c     &+3*n*nmodes+n*k)-(1.d0/de(k))*trrss(i,k)*(L*L)*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c      noniso
        AA(i+4*n+n*k,i+4*n+3*n*nmodes+n*k)=AA(i+4*n+n*k,i+4*n                   ! for crr, ctt
     &+3*n*nmodes+n*k)-funlam(i,k)*(1.d0/de(k))*trrss(i,k)*(L*L)*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c	  AA(i+4*n+n*k,i+4*n+3*n*nmodes+n*k)=AA(i+4*n+n*k,i+4*n                   ! for crr, ctt
c     &  +3*n*nmodes+n*k)-funlam(i,k)*(1.d0/de(k))*(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trrss(i,k)-
c     &  (L*L-3.d0)*(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*(1.d0/(L*L
c     &  -tththss(i,k)-2.d0*trrss(i,k))**2)

	 
c           +2*(L*L-3)*(1/r(i))*uthss(i)*trthss(i,k)*(1/(L*L-
c     &	  tththss(i,k)-2*trrss(i,k))**2)-2*(L*L-3)*(1/r(i))*uthss(i)* 
c     &	  trthss(i,k)*(1/(L*L-tththss(i,k)-2*trrss(i,k))**2)

c      iso
c        AA(i+4*n+n*k,i+4*n+5*n*nmodes+n*k)=AA(i+4*n+n*k,i+4*n+5*n*              ! for crr, czz
c     &  nmodes+n*k)-(1.d0/de(k))*trrss(i,k)*(L*L)*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c      noniso
	  AA(i+4*n+n*k,i+4*n+5*n*nmodes+n*k)=AA(i+4*n+n*k,i+4*n+5*n*              ! for crr, czz
     &nmodes+n*k)-funlam(i,k)*(1.d0/de(k))*trrss(i,k)*(L*L)*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c	  AA(i+4*n+n*k,i+4*n+5*n*nmodes+n*k)=AA(i+4*n+n*k,i+4*n+5*n*              ! for crr, czz
c     &  nmodes+n*k)-funlam(i,k)*(1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)*trrss(i,k)-(L*L-3.d0)
c     &  *(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k))**2)


     
c            +2*(L*L-3)*uthss(i)*trrss(i,k)*csi*(1/(L*L-tththss(i,k)
c     &	  -2*trrss(i,k))**2)-2*(L*L-3)*(1/r(i))*uthss(i)*trthss(i,k)*
c     &	  (1/(L*L-tththss(i,k)-2*trrss(i,k))**2)
	  
	      
c		  AA(i+4*n+n*k,i+4*n+n*nmodes+n*k)=                                   ! for trr, -exp(-eps_pi*(1/Tem-1))*alpha*2*Trt*De*trt
c     &	    AA(i+4*n+n*k,i+4*n+n*nmodes+n*k)
c     &        -funp(i,k)*alp(k)*2.d0*trthss(i,k)*de(k)


c       iso
c        AA(i+4*n+n*nmodes+n*k,i+n)=AA(i+4*n+n*nmodes+n*k,i+n)-                  ! for crt, ur
c     &dtrthss(i,k)+tththss(i,k)*(1.d0/r(i))*csi+trthss(i,k)*(1.0/r(i))
c     &-(2.d0*dtrrss(i,k)+dtththss(i,k))
c     &*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*
c     &(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)

c      noniso
	  AA(i+4*n+n*nmodes+n*k,i+n)=AA(i+4*n+n*nmodes+n*k,i+n)-                  ! for crt, ur
     &dtrthss(i,k)+tththss(i,k)*(1.d0/r(i))*csi+trthss(i,k)*(1.0/r(i))
     &-(2.d0*dtrrss(i,k)+dtththss(i,k))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)
     &+dtemss(i)*(1.d0/temss(i))*trthss(i,k)



c        AA(i+4*n+n*nmodes+n*k,i+n)=AA(i+4*n+n*nmodes+n*k,i+n)-                  ! for crt, ur
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*
c     &  (dtththss(i,k)+2.d0*dtrrss(i,k))*trthss(i,k)-(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*dtrthss(i,k)+
c     &  (L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)
c     &  ))*tththss(i,k)*csi-csi-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k)))*trthss(i,k)*(1.d0/r(i))-(L*L-3.d0)
c     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)
c     &  *dtemss(i)/temss(i)+funlam(i,k)*funp3(i,k)*csi*(1.d0/r(i))



c	      AA(i+4*n+n*nmodes+n*k,i+n)=AA(i+4*n+n*nmodes+n*k,i+n)+de(k)*        ! for trt, De*(-dTrt/dr+Trt/r+i*si*Ttt/r+dTem/dr*Trt/Tem)*ur
c     &	    (-dtrthss(i,k)+(trthss(i,k)/r(i))+(tththss(i,k)*csi               !         +i*si*exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Temp-1))/r*ur
c     &	    /r(i))+q(2)*(dtemss(i)*
c     &	    trthss(i,k)/temss(i)))+(funlam(i,k)*funp3(i,k)*csi/r(i))

c       iso
c        AA(i+4*n+n*nmodes+n*k,i+2*n)=AA(i+4*n+n*nmodes+n*k,i+2*n)-              ! for crt, ut
c     &  (1.d0/r(i))*trrss(i,k)+trthss(i,k)*(1.d0/r(i))*csi

c       noniso
        AA(i+4*n+n*nmodes+n*k,i+2*n)=AA(i+4*n+n*nmodes+n*k,i+2*n)-              ! for crt, ut
     &  (1.d0/r(i))*trrss(i,k)+trthss(i,k)*(1.d0/r(i))*csi



c        AA(i+4*n+n*nmodes+n*k,i+2*n)=AA(i+4*n+n*nmodes+n*k,i+2*n)-              ! for crt, ut
c     &  (L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k)))*trrss(i,k)+(1.d0/r(i))+(L*L-3.d0)*(1.d0/(L*L
c     &  -tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)*csi*(1.d0/r(i))-
c     &  funlam(i,k)*funp3(i,k)*(1.d0/r(i))


c	      AA(i+4*n+n*nmodes+n*k,i+2*n)=AA(i+4*n+n*nmodes+n*k,i+2*n)           ! for trt, De*(-Trr/r+i*si*Trt/r)*ut
c     &	    +de(k)*(-trrss(i,k)/r(i)+(trthss(i,k)*csi/r(i)))-                 !             -exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Temp-1))/r*ut
c     &	    funlam(i,k)*funp3(i,k)/r(i)

c       iso
c        AA(i+4*n+n*nmodes+n*k,i+4*n)=AA(i+4*n+n*nmodes+n*k,i+4*n)               ! for crt, tem 
c     &  +(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *trthss(i,k)*(1.d0/r(i))*uthss(i)*csi*(1.d0/temss(i))+funli(i,k)
c     &  *((L*L-3.d0)*uthss(i)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0
c     &  *trrss(i,k)))*trrss(i,k)-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k)))*trrss(i,k)*duthss(i)+duthss(i))-funp2(i,k)
c     &  *funlam(i,k)*funp3(i,k)*(-uthss(i)/r(i)+duthss(i))

c       noniso
        AA(i+4*n+n*nmodes+n*k,i+4*n)=AA(i+4*n+n*nmodes+n*k,i+4*n)               ! for crt, tem 
     &+(1.d0/r(i))*uthss(i)*csi*(1.d0/temss(i))*trthss(i,k)
     &+funli(i,k)*((1.0/r(i))*uthss(i)*trrss(i,k)-(1.d0/r(i))
     &*uthss(i)*tththss(i,k)+(1.d0/r(i))*uthss(i)*tththss(i,k)-
     &trrss(i,k)*duthss(i))

              	  
   

c	      AA(i+4*n+n*nmodes+n*k,i+4*n)=AA(i+4*n+n*nmodes+n*k,i+4*n)+          ! for trt, De*(i*si*Uth*Trt/Tem/r-eps_pi/Tem**2*(Trr*Trt+Trt*Ttt)*
c     &	    de(k)*((uthss(i)*csi*trthss(i,k)/(r(i)*temss(i)))-                !             exp(-eps_pi*(1/Tem-1))*alpha)*tem+exp(-eps_li*(1/Tem-1))*
c     &	    funp2(i,k)*(trrss(i,k)*trthss(i,k)+trthss(i,k)*                   !               (eps_li/Tem**2-eps_pi/Tem**2)*exp(eps_pi*(1/Temp-1))*
c     &	    tththss(i,k))*funp(i,k)*alp(k))+funlam(i,k)*(funli(i,k)-          !                  (dUth/dr-Uth/r)*tem
c     &	    funp2(i,k))*funp3(i,k)*(duthss(i)-(uthss(i)/r(i)))                !           -exp(-eps_li*(1/Tem-1))*eps_li/Tem**2*Trt*tem
c     &	    -funli(i,k)*funlam(i,k)*trthss(i,k)	

c???	      AA(i+4*n+n*nmodes+n*k,i+4*n)=                                       ! for trt, ?
c     &		q(1)* AA(i+4*n+n*nmodes+n*k,i+4*n)

c      iso
c          AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)=                                       ! for crt, crr
c     &AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)-(1.d0/r(i))*uthss(i)+duthss(i)
c     &-trthss(i,k)*(1.d0/de(k))*(L*L)
c     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trthss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))

c      noniso
          AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)=                                       ! for crt, crr
     &AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)-(1.d0/r(i))*uthss(i)+duthss(i)
     &-trthss(i,k)*funlam(i,k)*(1.d0/de(k))*(L*L)
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &-(1.d0/r(i))*uthss(i)*trthss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))



c        AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)=                                       ! for crt, crr
c     &  AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)-funlam(i,k)*(1.d0/de(k))*
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
c     &  *trthss(i,k)-(L*L-3.d0)*uthss(i)*(1.d0/r(i))*trthss(i,k)*csi*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)-(L*L-3.d0)
c     &  *(1.d0/r(i))*uthss(i)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  -(L*L-3.d0)*(1.d0/r(i))*uthss(i)*trrss(i,k)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)+(L*L-3.d0)*(1.d0/
c     &  (L*L-tththss(i,k)-2.d0*trrss(i,k)))*duthss(i)+(L*L-3.d0)
c     &  *trrss(i,k)*duthss(i)*(1.d0/(L*L-tththss(i,k)-
c     &  2.d0*trrss(i,k))**2)





c	      AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)=                                   ! for trt, De*(-Uth/r+dUth/dr-exp(-eps_pi*(1/Tem-1))*
c     &	    AA(i+4*n+n*nmodes+n*k,i+4*n+n*k)+                                 !            alpha*Trt)*trr
c     &        de(k)*((-uthss(i)/r(i))+duthss(i)-(funp(i,k)*
c     &	    alp(k)*trthss(i,k)))


c      iso
c        AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                              ! for crt,crt
c     &  AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)-
c     &  (1.d0/r(i))*uthss(i)*csi-(1.d0/de(k))*(L*L)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))


c      noniso
        AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                               ! for crt,crt
     &  AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)-
     &  (1.d0/r(i))*uthss(i)*csi-funlam(i,k)*(1.d0/de(k))*(L*L)*
     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))

	


c        AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                              ! for crt,crt
c     &  AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)-funlam(i,k)*  
c     &  (1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0
c     &  *trrss(i,k)))-(L*L-3.d0)*(1.d0/r(i))*uthss(i)*(1.d0/(L*L
c     &  -tththss(i,k)-2.d0*trrss(i,k)))*csi


c      iso
c	  AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                            ! for crt,ctt
c    &  AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-(1.d0/de(k))
c     & *trthss(i,k)*(L*L)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     & *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     & -(1.d0/r(i))*uthss(i)*trthss(i,k)*csi*
c     & (1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))

c      noniso
	  AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                            ! for crt,ctt
     &AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)
     &-funlam(i,k)*(1.d0/de(k))*trthss(i,k)*(L*L)
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &-(1.d0/r(i))*uthss(i)*trthss(i,k)*csi
     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     

c	  AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                            ! for crt,ctt
c     &  AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-funlam(i,k)*
c     &  (1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*trthss(i,k)-(L*L-3.d0)*(1.d0/r(i))*uthss(i)*
c     &  trthss(i,k)*csi*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
c     &  -(L*L-3.d0)*(1.d0/r(i))*uthss(i)*trrss(i,k)*(1.d0/(L*L-tththss
c     &  (i,k)-2.d0*trrss(i,k))**2)+(L*L-3.d0)*trrss(i,k)*duthss(i)
c     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2) 
c      iso


c	  AA(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                                ! for crt,czz
c     &  AA(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)-(1.d0/de(k))
c     & *trthss(i,k)*(L*L)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     & *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trthss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c      noniso

	  AA(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                                ! for crt,czz
     &AA(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)-
     &funlam(i,k)*(1.d0/de(k))*trthss(i,k)*(L*L)
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &-(1.d0/r(i))*uthss(i)*trthss(i,k)*csi
     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c	  AA(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                                ! for crt,czz
c     &  AA(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)-funlam(i,k)*
c     &  (1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*trthss(i,k)-(L*L-3.d0)*(1.d0/r(i))*uthss(i)*
c     &  trthss(i,k)*csi*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)-
c     &  (L*L-3.d0)*(1.d0/r(i))*uthss(i)*trrss(i,k)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)+(L*L-3.d0)*trrss(i,k)
c     &  *duthss(i)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)

		    
c		  AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                          ! for trt, -exp(-eps_li*(1/Tem-1))*trt
c     &	    AA(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)                         !            +De*(-i*si*Uth/r-exp(-eps_pi*(1/Tem-1))*alpha*
c     &	    -funlam(i,k)+de(k)*((-uthss(i)*csi/                               !                 (Trr+Ttt))*trt
c     &	    r(i))-(funp(i,k)*alp(k)*(trrss(i,k)+tththss(i,k))))
	     
c		  AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                        ! for trt, De*(-exp(-eps_pi*(1/Tem-1))*alpha*Trt)*ttt
c     &	    AA(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)
c     &	    +de(k)*(-funp(i,k)*alp(k)*trthss(i,k))

c        Omit

        AA(i+4*n+2*n*nmodes+n*k,i+n)=AA(i+4*n+2*n*nmodes+n*k,i+n)               ! for crz, ur
     &  +ca*trrss(i,k)


           
c        AA(i+4*n+2*n*nmodes+n*k,i+n)=AA(i+4*n+2*n*nmodes+n*k,i+n)               ! for crz, ur
c     &  +ca*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *trrss(i,k)-ca+ca*funp3(i,k)*funlam(i,k)


c	  AA(i+4*n+2*n*nmodes+n*k,i+n)=AA(i+4*n+2*n*nmodes+n*k,i+n)               ! for crz, ur
c     &  +ca*funp3(i,k)*funlam(i,k)

c	      AA(i+4*n+2*n*nmodes+n*k,i+n)=AA(i+4*n+2*n*nmodes+n*k,i+n)           ! for trz, +i*al*exp(eps_pi*(1/Temp-1))*exp(-eps_li*(1/Tem-1))*ur
c     &	    +funp3(i,k)*ca*funlam(i,k)

                AA(i+4*n+2*n*nmodes+n*k,i+3*n)=                                    ! for crz, uz
     &  AA(i+4*n+2*n*nmodes+n*k,i+3*n)+trthss(i,k)*(1.d0/r(i))*csi
 	
c        AA(i+4*n+2*n*nmodes+n*k,i+3*n)=                                         ! for crz, uz
c     &  AA(i+4*n+2*n*nmodes+n*k,i+3*n)+(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)*(1.d0/r(i))*csi


	    
c		  AA(i+4*n+2*n*nmodes+n*k,i+3*n)=                                     ! for trz, +i*si*De*Trt/r*uz
c    &	    AA(i+4*n+2*n*nmodes+n*k,i+3*n)+de(k)*trthss(i,k)*csi/r(i)


c       iso
c         AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                          ! for crz, crz
c     &  AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)-
c     &  (1.d0/r(i))*uthss(i)*csi-(1.d0/de(k))*(L*L)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))

c       noniso
         AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                          ! for crz, crz
     &  AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)-
     &  (1.d0/r(i))*uthss(i)*csi-funlam(i,k)*(1.d0/de(k))*(L*L)*
     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))



c        AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                          ! for crz, crz
c     &  AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)-funlam(i,k)*
c     &  (1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k)))-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k)))*uthss(i)*(1.d0/r(i))*csi


c		  AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                      ! for trz, -exp(-eps_li*(1/Tem-1))*trz-i*si*De*Uth/r*trz
c     &	    AA(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)                     !           -De*exp(-eps_pi*(1/Tem-1))*al*Trr*trz
c     &        -funlam(i,k)-(uthss(i)*csi*de(k)
c     &	    /r(i))-de(k)*trrss(i,k)*alp(k)*funp(i,k)	

	     
c		  AA(i+4*n+2*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                      ! for trz, -De*exp(-eps_pi*(1/Tem-1))*al*Trt*ttz
c     &	    AA(i+4*n+2*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)
c     &	    -de(k)*funp(i,k)*alp(k)*trthss(i,k)


c      iso
c	  AA(i+4*n+3*n*nmodes+n*k,i+n)=AA(i+4*n+3*n*nmodes+n*k,i+n)               ! for ctt, ur
c     &-dtththss(i,k)+2.d0*tththss(i,k)*(1.d0/r(i))
c     &-(2.d0*dtrrss(i,k)+dtththss(i,k))
c     &*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*
c     &(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*tththss(i,k)

c      noniso
	  AA(i+4*n+3*n*nmodes+n*k,i+n)=AA(i+4*n+3*n*nmodes+n*k,i+n)               ! for ctt, ur
     &-dtththss(i,k)+2.d0*tththss(i,k)*(1.d0/r(i))
     &-(2.d0*dtrrss(i,k)+dtththss(i,k))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*tththss(i,k)
     &+dtemss(i)*(1.d0/temss(i))*tththss(i,k)-(1.d0/(L*L-3.d0))*
     &(L*L-2.d0*trrss(i,k)-tththss(i,k))*dtemss(i)*(1.d0/temss(i))



c	  AA(i+4*n+3*n*nmodes+n*k,i+n)=AA(i+4*n+3*n*nmodes+n*k,i+n)               ! for ctt, ur
c     &  -(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*
c     &  (dtththss(i,k)+2*dtrrss(i,k))*tththss(i,k)-(L*L-3.d0)*
c     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*dtththss(i,k)
c     &  +2.d0*(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0
c     &  *trrss(i,k)))*tththss(i,k)-2.d0*(1.d0/r(i))+(L*L-3.d0)*(1.d0/
c     &  (L*L-tththss(i,k)-2.d0*trrss(i,k)))*tththss(i,k)*dtemss(i)*
c     &  (1.d0/temss(i))-dtemss(i)*(1.d0/temss(i))+funlam(i,k)*funp3(i,k)
c     &  *2.d0*(1.d0/r(i))

	     
       
c	      AA(i+4*n+3*n*nmodes+n*k,i+n)=AA(i+4*n+3*n*nmodes+n*k,i+n)           ! for ttt, +De*(-dTtt/dr+2*Ttt/r+dTem/dr*Ttt/Tem)*ur
c     &	    +de(k)*(-dtththss(i,k)+(2.d0*tththss(i,k)/r(i))                   !           +2*exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Temp-1))/r*ur
c     &        +q(3)*(dtemss(i)*tththss(i,k)/temss(i)))
c     &        +funlam(i,k)*funp3(i,k)*2.d0/r(i)


c       iso
c        AA(i+4*n+3*n*nmodes+n*k,i+2*n)=                                         ! for ctt, ut   
c     &  AA(i+4*n+3*n*nmodes+n*k,i+2*n)-2.d0*(1.d0/r(i))*trthss(i,k)
c     &  +2.d0*tththss(i,k)*(1.d0/r(i))*csi

c        AA(i+4*n+3*n*nmodes+n*k,i+2*n)=                                         ! for ctt, ut   
c     &  AA(i+4*n+3*n*nmodes+n*k,i+2*n)-(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*trthss(i,k)*2.d0*(1.d0/r(i))
c     &  +2.d0*(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0
c     &  *trrss(i,k)))*tththss(i,k)*csi-(2.d0/r(i))*csi+funlam(i,k)
c     &  *funp3(i,k)*2.d0*csi*(1.d0/r(i))

c      noniso

        AA(i+4*n+3*n*nmodes+n*k,i+2*n)=                                         ! for ctt, ut   
     &  AA(i+4*n+3*n*nmodes+n*k,i+2*n)-2.d0*(1.d0/r(i))*trthss(i,k)
     &  +2.d0*tththss(i,k)*(1.d0/r(i))*csi

	
 
c		  AA(i+4*n+3*n*nmodes+n*k,i+2*n)=                                     ! for ttt, +De*(-2*Trt/r+i*si*2*Ttt/r)*ut
c     &		AA(i+4*n+3*n*nmodes+n*k,i+2*n)                                    !           +i*si*2*exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Temp-1))/r*ut
c     &	    +de(k)*((-2.d0*trthss(i,k)/r(i))+(2.d0*tththss(i,k)*
c     &	    csi/r(i)))+funlam(i,k)*funp3(i,k)*2.d0*csi/r(i)

c      iso
c        AA(i+4*n+3*n*nmodes+n*k,i+4*n)=                                         ! for ctt, tem
c     &  AA(i+4*n+3*n*nmodes+n*k,i+4*n)+(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*(1.d0/r(i))*uthss(i)*csi*
c     &  (1.d0/temss(i))*tththss(i,k)-(1.d0/r(i))*uthss(i)*csi*(1.d0/
c     &  temss(i))+funli(i,k)*((L*L-3)*2*(1.d0/r(i))*(1/(L*L-tththss
c     &  (i,k)-2.d0*trrss(i,k)))*uthss(i)*trthss(i,k)-(L*L-3.d0)*2.d0
c     &  *(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *duthss(i)*trthss(i,k))

c      noniso
        AA(i+4*n+3*n*nmodes+n*k,i+4*n)=                                         ! for ctt, tem
     &  AA(i+4*n+3*n*nmodes+n*k,i+4*n)+(1.d0/r(i))*uthss(i)*csi
     &*(1.d0/temss(i))*tththss(i,k)-(1.d0/(L*L-3.d0))*
     &(L*L-2.d0*trrss(i,k)-tththss(i,k))*(1.d0/r(i))*uthss(i)*
     &csi*(1.d0/temss(i))
     &+funli(i,k)*((2.d0/r(i))*uthss(i)*trthss(i,k)-2.d0*
     &trthss(i,k)*duthss(i))


c	      AA(i+4*n+3*n*nmodes+n*k,i+4*n)=                                     ! for ttt, +De*(i*si*Uth*Ttt/r/Tem-eps_pi/Tem**2*(Trt**2+Ttt*2)*
c     &		AA(i+4*n+3*n*nmodes+n*k,i+4*n)                                    !              alpha*exp(-eps_pi*(1/Tem-1)))*tem-
c     &	    +de(k)*(uthss(i)*csi*tththss(i,k)/(r(i)*temss(i))-                !            exp(-eps_li*(1/Tem-1))*eps_li/Tem**2*Ttt*tem  
c     &	    funp2(i,k)*(trthss(i,k)**2+tththss(i,k)**2)*alp(k)
c     &	    *funp(i,k))-funlam(i,k)*funli(i,k)*tththss(i,k)
c	
c	  AA(i+4*n+3*n*nmodes+n*k,i+4*n) =                                        ! for ttt, ?
c     &  q(1)*AA(i+4*n+3*n*nmodes+n*k,i+4*n)

c      iso
c	  AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)=                                     ! for ctt,crr  
c     &AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)-(1.d0/de(k))*tththss(i,k)
c     &*(L*L)*(1.d0/(L*L-tththss(i,k)-2.0*trrss(i,k)))
c     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*tththss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))

c      noniso
	  AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)=                                     ! for ctt,crr  
     &AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)-funlam(i,k)
     &*(1.d0/de(k))*tththss(i,k)*(L*L)*
     &(1.d0/(L*L-tththss(i,k)-2.0*trrss(i,k)))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &-(1.d0/r(i))*uthss(i)*tththss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))



c	  AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)=                                     ! for ctt,crr  
c     &  AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)-funlam(i,k)*(1.d0/de(k))*
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
c     &  *tththss(i,k)-(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k))**2)*uthss(i)*csi*tththss(i,k)-2.d0*(L*L-3.d0)
c     &  *(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)
c     &  *trthss(i,k)*uthss(i)+2.d0*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k))**2)*duthss(i)*trthss(i,k)

c      iso and noniso                   
          AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                            ! for ctt, crt
     &  AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*nmodes+n*k)-
     &  2.d0*(1.d0/r(i))*uthss(i)+2.d0*duthss(i)   



c        AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                            ! for ctt, crt
c     &  AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*nmodes+n*k)-2.d0*(L*L-3.d0)
c     &  *(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*uthss(i)
c     &  +2.d0*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *duthss(i)



c	      AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                        ! for ttt, De*(-2*Uth/r+2*dUth/dr-2*exp(-eps_pi*(1/Tem-1))*alpha*Trt)*trt
c     &        AA(i+4*n+3*n*nmodes+n*k,i+4*n+n*nmodes+n*k)
c     &        +de(k)*((-2.d0*uthss(i)/r(i))
c     &	    +(2.d0*duthss(i))-(2.d0*funp(i,k)*alp(k)*trthss(i,k)))

c       iso
        
c        AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for ctt, ctt
c     &AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-(1.d0/r(i))
c     &*uthss(i)*csi-(1.d0/de(k))*tththss(i,k)
c     &*(L*L)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))-(1.d0/de(k)) 
c     &*(L*L)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))) 
c     &-(1.d0/r(i))*uthss(i)*tththss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))

c       noniso
        AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for ctt, ctt
     &AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-(1.d0/r(i))
     &*uthss(i)*csi-funlam(i,k)*((1.d0/de(k))*tththss(i,k)
     &*(L*L)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))+(1.d0/de(k)) 
     &*(L*L)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))) 
     &-(1.d0/r(i))*uthss(i)*tththss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))

   




c        AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for ctt, ctt
c     &  AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-funlam(i,k)* 
c     &  (1.d0/de(k))*((L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k)))+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0
c     &  *trrss(i,k))**2)*tththss(i,k))-(L*L-3.d0)*(1.d0/r(i))*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*uthss(i)*csi*
c     &  tththss(i,k)-(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-
c     &  2.d0*trrss(i,k)))*uthss(i)*csi-2.d0*(L*L-3.d0)*(1.d0/r(i))*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trthss(i,k)
c     &  *uthss(i)+2.d0*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k))**2)*duthss(i)*trthss(i,k)


	  
c		  AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                      ! for ttt, -exp(-eps_li*(1/Tem-1))*ttt+
c     &	    AA(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)                     !            De*(-i*si*Uth/r-2*exp(-eps_pi*(1/Tem-1))*alpha*Ttt)*ttt
c     &        -funlam(i,k)+de(k)*((-uthss(i)*csi
c     &	    /r(i))-(funp(i,k)*alp(k)*2.d0*tththss(i,k)))

c       iso
c        AA(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for ctt, czz
c     &AA(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)
c     &-(1.d0/de(k))*tththss(i,k)*(L*L)
c     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*tththss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c      noniso
       AA(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                            ! for ctt, czz
     &AA(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)
     &-funlam(i,k)*(1.d0/de(k))*tththss(i,k)*(L*L)
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
     &-(1.d0/r(i))*uthss(i)*tththss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c        AA(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for ctt, czz
c     &  AA(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)-funlam(i,k)*
c     &  (1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*tththss(i,k)-(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L
c     &  -tththss(i,k)-2.d0*trrss(i,k))**2)*uthss(i)*csi*tththss(i,k)-
c     &  2.d0*(L*L-3.d0)*(1.d0/r(i))*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*uthss(i)*trthss(i,k)+2.d0*(L*L-3.d0)*(1.d0/
c     &  (L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*duthss(i)*trthss(i,k)

	

c        Omit


c      iso
c        AA(i+4*n+5*n*nmodes+n*k,i+n)=                                           ! for czz, ur
c     &  AA(i+4*n+5*n*nmodes+n*k,i+n)-dtrrss(i,k)
c     &-(2.d0*dtrrss(i,k)+dtththss(i,k))
c     &*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*
c     &(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)

c      noniso	
        AA(i+4*n+5*n*nmodes+n*k,i+n)=                                           ! for czz, ur
     &  AA(i+4*n+5*n*nmodes+n*k,i+n)-dtrrss(i,k)
     &-(2.d0*dtrrss(i,k)+dtththss(i,k))
     &*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)
     &+dtemss(i)*(1.d0/temss(i))*trrss(i,k)-(1.d0/(L*L-3.d0))*
     &(L*L-2.d0*trrss(i,k)-tththss(i,k))*dtemss(i)*(1.d0/temss(i))




c        AA(i+4*n+5*n*nmodes+n*k,i+n)=                                           ! for czz, ur
c     &  AA(i+4*n+5*n*nmodes+n*k,i+n)-(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)*(dtththss(i,k)+2.d0*
c     &  dtrrss(i,k))*trrss(i,k)-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k)))*dtrrss(i,k)+(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*(dtemss(i)/temss(i))*trrss(i,k)
c     &  -dtemss(i)/temss(i)


c      iso and noniso
c        AA(i+4*n+5*n*nmodes+n*k,i+3*n)=                                         ! for czz, uz
c     &  AA(i+4*n+5*n*nmodes+n*k,i+3*n)+2.d0*trrss(i,k)
c     &+funp3(i,k)*funlam(i,k)*2.d0*ca



        AA(i+4*n+5*n*nmodes+n*k,i+3*n)=                                         ! for czz, uz
     &  AA(i+4*n+5*n*nmodes+n*k,i+3*n)+2.d0*(L*L-3.d0)*(1.d0/(L*L-
     &  tththss(i,k)-2.d0*trrss(i,k)))*trrss(i,k)*ca-2*ca+funlam(i,k)*
     &  funp3(i,k)*2.d0*ca

          
c        Omit

c      iso
c        AA(i+4*n+5*n*nmodes+n*k,i+4*n)=                                         ! for czz, tem
c     &  AA(i+4*n+5*n*nmodes+n*k,i+4*n)+(L*L-3.d0)*(1.d0/(L*L-
c    &  tththss(i,k)-2.d0*trrss(i,k)))*uthss(i)*(1.d0/r(i))*csi*
c     &  trrss(i,k)-(1.d0/r(i))*uthss(i)*csi*(1.d0/temss(i))

c      nonsio	
	    AA(i+4*n+5*n*nmodes+n*k,i+4*n)=                                         ! for czz, tem
     & AA(i+4*n+5*n*nmodes+n*k,i+4*n)+(1.d0/r(i))*uthss(i)*csi*
     & (1.d0/temss(i))*trrss(i,k)-(1.d0/(L*L-3.d0))*
     & (L*L-2.d0*trrss(i,k)-tththss(i,k))*(1.d0/r(i))*uthss(i)*csi*
     & (1.d0/temss(i))

c      iso          
c	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)=                                     ! for czz, crr
c     &AA(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)-(1.d0/de(k))*trrss(i,k)
c     &*(L*L)*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c      noniso
		  AA(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)=                                     ! for czz, crr
     &AA(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)-
     &funlam(i,k)*(1.d0/de(k))*trrss(i,k)*(L*L)
     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))




          
c	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)=                                     ! for czz, crr
c    &  AA(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)-funlam(i,k)*(1.d0/de(k))*
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*
c     &  trrss(i,k)-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*(1.d0/r(i))*uthss(i)*csi*trrss(i,k)

c      iso

c	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for czz, ctt
c     &AA(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-(1.d0/de(k))
c     &*trrss(i,k)*(L*L)*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


c      noniso
	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for czz, ctt
     &AA(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-
     &funlam(i,k)*(1.d0/de(k))*trrss(i,k)*(L*L)
     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))


  
c	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for czz, ctt
c     &  AA(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)-funlam(i,k)*
c     &  (1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*trrss(i,k)-(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)*(1.d0/r(i))*uthss(i)*csi
c     &  *trrss(i,k)

c      iso
c	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for czz, czz
c     &AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)-(1.d0/r(i))*uthss(i)
c     &*csi-(1.d0/de(k))*(L*L)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))  
c     &-(1.d0/de(k))*trrss(i,k)*(L*L)
c     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
c     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
c     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))

c      noniso
	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for czz, czz
     &AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)-(1.d0/r(i))*uthss(i)
     &*csi-funlam(i,k)*((1.d0/de(k))*(L*L)*  
     &(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))+(1.d0/de(k))*trrss(i,k)
     &*(L*L)*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))
     &*(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k))))
     &-(1.d0/r(i))*uthss(i)*trrss(i,k)*csi*
     &(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))




c	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for czz, czz
c     &  AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)-funlam(i,k)*
c     &  (1.d0/de(k))*((L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k)))+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*trrss(i,k))-(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)*(1.d0/r(i))*uthss(i)*csi*
c     &  trrss(i,k)-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-
c     &  2.d0*trrss(i,k)))*uthss(i)*csi*(1/r(i))

	  AA(i+4*n+4*n*nmodes+n*k,i+2*n)=                                         ! for ctz, ut
     &  AA(i+4*n+4*n*nmodes+n*k,i+2*n)+trrss(i,k)*ca



c	  AA(i+4*n+4*n*nmodes+n*k,i+2*n)=                                         ! for ctz, ut
c     &  AA(i+4*n+4*n*nmodes+n*k,i+2*n)+funlam(i,k)*funp3(i,k)*ca
c     &  +(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *trrss(i,k)*ca-ca



c	      AA(i+4*n+4*n*nmodes+n*k,i+2*n)=                                     ! for ttz, +i*al*exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Temp-1))*ut
c     &	    AA(i+4*n+4*n*nmodes+n*k,i+2*n)+funp3(i,k)*funlam(i,k)*ca


	  AA(i+4*n+4*n*nmodes+n*k,i+3*n)=                                         ! for ctz, uz
     &  AA(i+4*n+4*n*nmodes+n*k,i+3*n)+tththss(i,k)*csi*(1.d0/r(i))
       
c	  AA(i+4*n+4*n*nmodes+n*k,i+3*n)=                                         ! for ctz, uz
c     &  AA(i+4*n+4*n*nmodes+n*k,i+3*n)+(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k)))*tththss(i,k)*csi*(1.d0/r(i))
c     &  -csi*(1.d0/r(i))+funlam(i,k)*funp3(i,k)*csi*(1.d0/r(i))



	
c	      AA(i+4*n+4*n*nmodes+n*k,i+3*n)=                                     ! for ttz, +i*si*De*Ttt/r*uz+
c     &		AA(i+4*n+4*n*nmodes+n*k,i+3*n)                                    !           i*si*exp(-eps_li*(1/Tem-1))*exp(eps_pi*(1/Temp-1))/r*uz 
c     &	    +(de(k)*tththss(i,k)*csi/r(i))
c     &        +(funlam(i,k)*funp3(i,k)*csi/r(i))


      AA(i+4*n+4*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                          ! for ctz, crz
     &AA(i+4*n+4*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)-(1.d0/r(i))*uthss(i)
     &+duthss(i)
 

c        AA(i+4*n+4*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                          ! for ctz, crz
c     &  AA(i+4*n+4*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)-(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*(1.d0/r(i))*uthss(i)+
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*duthss(i)

c      iso
c         AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                          ! for ctz, ctz
c     &  AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)-(1.d0/r(i))*
c     &  uthss(i)*csi-(1.d0/de(k))*(L*L)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))


c      noniso
	         AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                          ! for ctz, ctz
     &  AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)-(1.d0/r(i))*
     &  uthss(i)*csi-funlam(i,k)*(1.d0/de(k))*(L*L)*
     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))

 
c        AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                          ! for ctz, ctz
c     &  AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)-funlam(i,k)*
c     &  (1.d0/de(k))*(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k
c     &  )))-(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *uthss(i)*(1.d0/r(i))*csi
			 
       	      	
		  	     
c		  AA(i+4*n+4*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                      ! for ttz, +De*(-Uth/r+dUth/dr-exp(-eps_pi*(1/Tem-1))*alpha*Trt)*trz
c     &	  AA(i+4*n+4*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)
c     &      +de(k)*((-uthss(i)/
c     &	  r(i))+duthss(i)-(funp(i,k)*alp(k)*trthss(i,k)))
	     
c		  AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                      ! for ttz, -exp(-eps_li*(1/Tem-1))*ttz
c     &	  AA(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)                       !           +De*(-i*si*Uth/r-exp(-eps_pi*(1/Tem-1))*alpha*Ttt)*ttz
c     &      -funlam(i,k)+de(k)
c     &	  *((-uthss(i)*csi/r(i))-(funp(i,k)*alp(k)*tththss(i,k)))

c	      AA(i+4*n+5*n*nmodes+n*k,i+3*n)=                                     ! for tzz, +i*al*2*exp(eps_pi*(1/Temp-1))*exp(-eps_li*(1/Tem-1)*uz
c    &		AA(i+4*n+5*n*nmodes+n*k,i+3*n)
c     &	    +funp3(i,k)*funlam(i,k)*2.d0*ca
	
c	      AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                      ! for tzz, -exp(-eps_li*(1/Tem-1)*tzz-i*si*De*Uth/r*tzz
c     & 	  AA(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)
c     &      -funlam(i,k)-(de(k)*uthss(i)*csi/r(i))

	    enddo
	  enddo
	endif


c	boundary conditions	
	
	do i=1,nvars
	  AA(n+1,i)=dcmplx(0.d0,0.d0)                                             ! BC for ur at r=r1 or r2, ur=0 
	  AA(i,n+1)=dcmplx(0.d0,0.d0)                                              
	  AA(2*n,i)=dcmplx(0.d0,0.d0)                                             
	  AA(i,2*n)=dcmplx(0.d0,0.d0)                                             
	
	  AA(2*n+1,i)=dcmplx(0.d0,0.d0)                                           ! BC for ut at r=r1 or r2, ur=0                                           
	  AA(i,2*n+1)=dcmplx(0.d0,0.d0)                                           
	  AA(3*n,i)=dcmplx(0.d0,0.d0)                                             
	  AA(i,3*n)=dcmplx(0.d0,0.d0)                                              
	
	  AA(3*n+1,i)=dcmplx(0.d0,0.d0)                                           ! BC for ur at r=r1 or r2, uz=0                                           
	  AA(i,3*n+1)=dcmplx(0.d0,0.d0)                                           
	  AA(4*n,i)=dcmplx(0.d0,0.d0)                                             
	  AA(i,4*n)=dcmplx(0.d0,0.d0)                                             

	  AA(4*n+1,i)=dcmplx(0.d0,0.d0)                                           ! BC for tempature at r=r1, tem=0    

c	  turn off below term for Neumann
	  AA(i,4*n+1)=dcmplx(0.d0,0.d0)
	  AA(5*n,i)=dcmplx(0.d0,0.d0)                                             ! BC for tempature at r=r2, tem=0

c	  turn off below term for Neumann
	  AA(i,5*n)=dcmplx(0.d0,0.d0)
	enddo

	AA(n+1,n+1)=dcmplx(-20.d0,0.d0)
	AA(2*n,2*n)=dcmplx(-20.d0,0.d0)
	AA(2*n+1,2*n+1)=dcmplx(-20.d0,0.d0)
	AA(3*n,3*n)=dcmplx(-20.d0,0.d0)
	AA(3*n+1,3*n+1)=dcmplx(-20.d0,0.d0)
	AA(4*n,4*n)=dcmplx(-20.d0,0.d0)

c	turn off both terms below for neumann
	AA(4*n+1,4*n+1)=dcmplx(-20.d0,0.d0)
	AA(5*n,5*n)=dcmplx(-20.d0,0.d0)

c	include below do loop and the next for neumann
c	do i=1,n
c	  AA(4*n+1,4*n+i)=dfirst(1,i)
c        AA(5*n,4*n+i)=-dfirst(n,i)
c         enddo
c        AA(5*n,5*n)=AA(5*n,5*n)-10.d0

	return	
	end
	
c	------------END OF SUBROUTINE CONSTRUCTMATA-------------------

c	----------SUBROUTINE CONSTRUCTMATB------------------------

	subroutine constructmatB(BB,q)

	implicit double precision (a-h,o-z)
      parameter(iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
      parameter(nvars=(6*nmodes+5)*n,nrhs=1,nn=2*nvars)
      common/variables/ eigvar(nvars),deigvar(nvars),
     &  ddeigvar(nvars),y(n),r(n),der(n,n),dder(n,n),dfirst(n,n),
     &  dsecond(n,n)
      common/steadyvar/uthss(n),temss(n),trrss(n,nrows),
     & 	trthss(n,nrows), tththss(n,nrows),duthss(n),dtemss(n),
     &	dduthss(n),ddtemss(n),dtrrss(n,nrows),ddtrrss(n,nrows),
     &	dtrthss(n,nrows),ddtrthss(n,nrows),dtththss(n,nrows),
     &	ddtththss(n,nrows)
      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &    epsilonlam(nrows),alp(nrows),delham(nrows),
     &    delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &    gammadot
c	double precision L
	common/chain/L
     	common/wavenum/si,alpha
 	complex*16 AA(nvars,nvars), BB(nvars,nvars), csi,ca
	double precision fun1(n),fun2(n),funp(n,nrows),
     &	funlam(n,nrows),funp2(n,nrows)
	double precision strrss(n),strthss(n),stththss(n),sum(3),q(10)
      integer i,j,k
   
	data pi/3.1415926535897930d0/
c	data L/1000.d0/
	do i=1,nvars
	  do j=1,nvars	
	    BB(i,j)=dcmplx(0.d0,0.d0)	
	  enddo
	enddo
c

	if (nmodes.NE.0) then
	  do k=1,nmodes
	    do i=1,n
	      funlam(i,k)=dexp(epsilonlam(k)*(1.d0/temss(i)-1.d0))
	    enddo
	  enddo	
	endif
c		
c	write(*,*) 'Re=',Re, 'Pe=',Pe, 'T0=',T0, 'T1=',T1, 'T2=',T2
c     write(*,*) 'epsilons=',epsilons, 'betas=',betas
c     write(*,*) 'epsilons=',epsilons,'si=',si, 'alpha=',alpha


	do i=1,n
	  BB(i+n,i+n)=BB(i+n,i+n)+Re                                               ! for ur, sig*Re*ur
	  BB(i+2*n,i+2*n)=BB(i+2*n,i+2*n)+Re                                       ! for ut, sig*Re*ut
	  BB(i+3*n,i+3*n)=BB(i+3*n,i+3*n)+Re                                       ! for uz, sig*Re*tz
	  BB(i+4*n,i+4*n)=BB(i+4*n,i+4*n)+Pe                                       ! for tem, sig*Pe*tem
	enddo


	if (nmodes.NE.0) then
	  do k=1,nmodes
	    do i=1,n
c       omit sigma-T

c      iso
c        BB(i+4*n+n*k,i+4*n)=BB(i+4*n+n*k,i+4*n)-(L*L-3.d0)*(1.d0/               ! for crr, tem
c     &  temss(i))*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))
c     &  *trrss(i,k)+1.d0/temss(i)

c      noniso

        BB(i+4*n+n*k,i+4*n)=BB(i+4*n+n*k,i+4*n)-                                 ! for crr, tem
     &(1.d0/temss(i))*trrss(i,k)+(1.d0/(L*L-3.d0))*(L*L-2.d0*trrss(i,k)
     &-tththss(i,k))*(1.d0/temss(i))

c	      BB(i+4*n+n*k,i+4*n)=BB(i+4*n+n*k,i+4*n)-trrss(i,k)*de(k)             ! for trr, -sig*De*Trr/Tem*tem
c     &	       /temss(i)


c       iso and noniso
       BB(i+4*n+n*k,i+4*n+n*k)=BB(i+4*n+n*k,i+4*n+n*k)+1.d0                     ! for crr, crr
     & +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trrss(i,k)



c        BB(i+4*n+n*k,i+4*n+n*k)=BB(i+4*n+n*k,i+4*n+n*k)+(L*L-3.d0)*             ! for crr, crr
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trrss(i,k)+ 
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))

	      
c		  BB(i+4*n+n*k,i+4*n+n*k)=BB(i+4*n+n*k,i+4*n+n*k)+de(k)                ! for trr, sig*De*trr

c        BB(i+4*n+n*k,i+4*n+3*n*nmodes+n*k)=BB(i+4*n+n*k,i+4*n+3*n*               ! for crr, ctt
c     &  nmodes+n*k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0
c     &  *trrss(i,k))**2)*trrss(i,k)

c       iso and noniso

       BB(i+4*n+n*k,i+4*n+3*n*nmodes+n*k)=                                        ! for crr,ctt
     & BB(i+4*n+n*k,i+4*n+3*n*nmodes+n*k)
     & +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trrss(i,k)


c        BB(i+4*n+n*k,i+4*n+5*n*nmodes+n*k)=BB(i+4*n+n*k,i+4*n+5*n*               ! for crr, czz
c     &  nmodes+n*k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*
c     &  trrss(i,k))**2)*trrss(i,k)


c       iso and noniso
       BB(i+4*n+n*k,i+4*n+5*n*nmodes+n*k)=                                        ! for crr,czz
     & BB(i+4*n+n*k,i+4*n+5*n*nmodes+n*k)
     & +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trrss(i,k)

c       omit sigma-T
c       iso	      
c       BB(i+4*n+n*nmodes+n*k,i+4*n)=BB(i+4*n+n*nmodes+n*k,i+4*n)                ! for crt, tem
c     &  -(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))*
c     &  trthss(i,k)*(1.d0/temss(i))

c       nonsio
	   BB(i+4*n+n*nmodes+n*k,i+4*n)=BB(i+4*n+n*nmodes+n*k,i+4*n)                ! for crt, tem
     &  -(1.d0/temss(i))*trthss(i,k)



c		  BB(i+4*n+n*nmodes+n*k,i+4*n)=BB(i+4*n+n*nmodes+n*k,i+4*n)            ! for trt, -sig*De*Trt/Tem*tem
c     &	       -q(1)*trthss(i,k)*de(k)/temss(i)


c        BB(i+4*n+n*nmodes+n*k,i+4*n+n*k)=BB(i+4*n+n*nmodes+n*k,i+4*n             ! for crt, crr
c     &  +n*k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*
c     &  trthss(i,k)


c       iso and noniso
        BB(i+4*n+n*nmodes+n*k,i+4*n+n*k)=                                       ! for crt, crr
     &  BB(i+4*n+n*nmodes+n*k,i+4*n+n*k)
     &  +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trthss(i,k)

c       iso and noniso
        BB(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)=BB(i+4*n+n*nmodes             ! for crt, crt
     &  +n*k,i+4*n+n*nmodes+n*k)+1.d0 

c       iso and noniso	
        BB(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                              ! for crt, ctt
     &  BB(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)
     &  +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trthss(i,k)

c       iso and noniso	
        BB(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                              ! for crt, czz
     &  BB(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)
     &  +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trthss(i,k)

c        BB(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)=BB(i+4*n+n*nmodes             ! for crt, crt
c     &  +n*k,i+4*n+n*nmodes+n*k)+(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)
c     &  -2.d0*trrss(i,k)))                                    

            
c	  BB(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                            ! for crt, ctt
c     &  BB(i+4*n+n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)+(L*L-3.d0)*(1.d0/
c     &  (L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trthss(i,k)


c        BB(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                                ! for crt, czz
c     &  BB(i+4*n+n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)+(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trthss(i,k)

c       iso and noniso
	  BB(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                          ! for crz, crz
     &  BB(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)+1.d0


c	  BB(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                          ! for crz, crz
c     &  BB(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)+
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))


c        BB(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)=                                     ! for ctt, crr
c     &  BB(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)+(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2.d0*trrss(i,k))**2)*tththss(i,k)

c       iso and noniso
       BB(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)=                                      ! for ctt, crr
     & BB(i+4*n+3*n*nmodes+n*k,i+4*n+n*k)
     &  +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*tththss(i,k)
c       iso and noniso
	  BB(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for ctt, ctt
     &  BB(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)+1.d0
     &  +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*tththss(i,k)
c       iso and noniso
	  BB(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for ctt, czz
     &  BB(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)
     &  +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*tththss(i,k)

c	  BB(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for ctt, ctt
c     &  BB(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)+(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*tththss(i,k)+
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))


c	  BB(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for ctt, czz
c     &  BB(i+4*n+3*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)+(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*tththss(i,k)
            
c       omit sigma-T

c       iso		   
c	  BB(i+4*n+3*n*nmodes+n*k,i+4*n)=                                         ! for ctt, tem
c     &  BB(i+4*n+3*n*nmodes+n*k,i+4*n)-(L*L-3.d0)*(1.d0/(L*L-
c     &  tththss(i,k)-2*trrss(i,k)))*tththss(i,k)*(1.d0/temss(i))
c     &  +1.d0/temss(i)


c       noniso		   
	  BB(i+4*n+3*n*nmodes+n*k,i+4*n)=                                         ! for ctt, tem
     &  BB(i+4*n+3*n*nmodes+n*k,i+4*n)-(1.d0/temss(i))*tththss(i,k)+
     &  (1.d0/(L*L-3.d0))
     &  *(L*L-2.d0*trrss(i,k)-tththss(i,k))*(1.d0/temss(i))


c       iso and noniso
	  BB(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                          ! for ctz, ctz
     &  BB(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)+1.d0


c	  BB(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                          ! for ctz, ctz
c     &  BB(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)+
c     &  (L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))


c	  BB(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)=                                     ! for czz, crr
c     &  BB(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)+(L*L-3.d0)*(1.d0/(L*L-
c     &	  tththss(i,k)-2.d0*trrss(i,k))**2)*trrss(i,k)
		 	  		

			  
c        BB(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                          ! for czz, ctt
c     &  BB(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)+(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trrss(i,k)

c       iso and noniso
       BB(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)=                                      ! for czz, crr
     & BB(i+4*n+5*n*nmodes+n*k,i+4*n+n*k)
     & +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trrss(i,k)

c       iso and noniso
       BB(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                           ! for czz, ctt
     & BB(i+4*n+5*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)
     & +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trrss(i,k)

c       iso and noniso
	  BB(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for czz, czz
     &  BB(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)+1.d0
     &  +(1.d0/(L*L-2.d0*trrss(i,k)-tththss(i,k)))*trrss(i,k)
		  
c	  BB(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                          ! for czz, czz
c     &  BB(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)+(L*L-3.d0)*
c     &  (1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k))**2)*trrss(i,k)
c     &  +(L*L-3.d0)*(1.d0/(L*L-tththss(i,k)-2.d0*trrss(i,k)))

c       omit sigma-T
c       iso 
c	      BB(i+4*n+5*n*nmodes+n*k,i+4*n)=                                      ! for czz, tem
c     &	  BB(i+4*n+5*n*nmodes+n*k,i+4*n)-(L*L-3.d0)*(1.d0/(L*L-
c     &	  tththss(i,k)-2.d0*trrss(i,k)))*(1.d0/temss(i))+1.d0/temss(i)



c      noniso

	 BB(i+4*n+5*n*nmodes+n*k,i+4*n)=                                          ! for czz, tem
     & BB(i+4*n+5*n*nmodes+n*k,i+4*n)-(1.d0/temss(i))*trrss(i,k)+
     & (1.d0/(L*L-3.d0))
     & *(L*L-2.d0*trrss(i,k)-tththss(i,k))*(1.d0/temss(i))



c		  BB(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)=                           ! for trt, sig*De*trt
c     &	       BB(i+4*n+n*nmodes+n*k,i+4*n+n*nmodes+n*k)+de(k)
	      
c		  BB(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)=                       ! for trz, sig*De*trz
c     &	  BB(i+4*n+2*n*nmodes+n*k,i+4*n+2*n*nmodes+n*k)+de(k)
	
c	      BB(i+4*n+3*n*nmodes+n*k,i+4*n)=                                      ! for ttt, -sig*De*Ttt/Tem*tem
c     &	  BB(i+4*n+3*n*nmodes+n*k,i+4*n)
c     &      -q(1)*tththss(i,k)*de(k)/temss(i)
	      
c		  BB(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)=                       ! for ttt, sig*De*ttt
c     &	  BB(i+4*n+3*n*nmodes+n*k,i+4*n+3*n*nmodes+n*k)+de(k)
	      
c		  BB(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)=                       ! for ttz, sig*De*ttz
c     &	  BB(i+4*n+4*n*nmodes+n*k,i+4*n+4*n*nmodes+n*k)+de(k)
	      
c		  BB(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)=                       ! for tzz, sig*De*tzz
c     &      BB(i+4*n+5*n*nmodes+n*k,i+4*n+5*n*nmodes+n*k)+de(k)
	    enddo
	  enddo
	endif


c	boundary conditions
	
	BB(n+1,n+1)=dcmplx(1.d0,0.d0)
	BB(2*n,2*n)=dcmplx(1.d0,0.d0)
	BB(2*n+1,2*n+1)=dcmplx(1.d0,0.d0)
	BB(3*n,3*n)=dcmplx(1.d0,0.d0)
	BB(3*n+1,3*n+1)=dcmplx(1.d0,0.d0)
	BB(4*n,4*n)=dcmplx(1.d0,0.d0)
c	turn off both terms below for neumann
	BB(4*n+1,4*n+1)=dcmplx(1.d0,0.d0)
	BB(5*n,5*n)=dcmplx(1.d0,0.d0)

c	include both 2 terms below for neumann
c	 BB(4*n+1,4*n+1)=dcmplx(0.d0,0.d0)
c       BB(5*n,5*n)=dcmplx(0.d0,0.d0)
	
	return
	
	
	end
	
c	-------------END OF SUBROUTINE CONSTRUCTMATB-----------------

	
