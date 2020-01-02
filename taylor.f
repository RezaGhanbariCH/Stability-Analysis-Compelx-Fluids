c     storage of varaibles
c1      n       2n       3n       4n ... (2+modes)*n  (3+modes)*n ... (2+2*modes)*n (2+2*modes)*n  ...  (2+3*modes)*n
c=====================================================================================================================       
c  Uth  |  Temp  | Trr_1  | Trr_2  | ... | Trt_1      |  Trt_2  | ... | Ttt_1       | Ttt_2      | ...  | Ttt_modes |
c=====================================================================================================================    
c	main program for taylor flow problem
c	modified April 24,2002  completed on April 24, 2002
c	Giesekus fluid 
c	-----------------------MAIN PROGRAM----------------
	subroutine taylor(velmean)
	implicit double precision (a-h,o-z)
	parameter (iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
	parameter (nvars = (3*nmodes+2)*n,nrhs=1)
	common/var/xvar(nvars),gre,y(n),r(n),der(n,n),dder(n,n),
     &    dxvar(nvars),ddxvar(nvars)
	common/poly/lamdap(nrows)
	common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &	epsilonlam(nrows),alpha(nrows),delhlam(nrows),
     &	delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond
	common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &    gammadot
c	double precision L
	common/chain/L
	common/matrix/resid(nvars),sumstress(nmodes),a(nvars,nvars),
     &    djacob(nvars,nvars),djacobinv(nvars,nvars),b(nvars,nrhs)
	double precision lamdap
	integer i,j,k,ipiv(nvars),info,iter
	double precision sum,idmatrix(nvars,nvars),u(n),trt(n,1),
     &	tth(n,1),ss(n),velmean	
	data pi/3.1415926535897930d0/
c	data L/1000.d0/


	uin=1.d0                   ! inner cylinder velocity 
	uout=0.d0                  ! outer cylinder velocity
	rot_p=-1.d0	               ! co-rotation parameter
	iter=0

	call input
	call meshpt(y,r,nmodes,n,delta)
	call derivmat(der,y,n)
	call squaremat(der,n,dder)
	call initialvalue(xvar,nvars,nmodes,r,T1,T2,T0,delta,u,de(1)
     &             ,uout,uin)                                                                                                                    
	
	
	do i=1,n
	  ss(i)=0.d0
	  do j=1,n
	    ss(i)=ss(i)+der(i,j)*(-2.d0)*u(j)/r(j)          ! dr(uth/r)

	  enddo 
	enddo
	



	do i=1,n
c	  trt(i,1)=r(i)*ss(i)                               ! T_rt for isothermal case
c	  tth(i,1)=2.d0*de(1)*trt(i,1)*r(i)*ss(i)           ! T_tt for isothermal case



        trt(i,1)=de(1)*r(i)*ss(i)*
     &  (L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))*
     &  (L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))/
     &  (L*L)/(L*L-3.d0)

	  tth(i,1)=2.d0*de(1)*de(1)*r(i)*ss(i)*r(i)*ss(i)*
     &(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))*
     &(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))*
     &(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))/(L*L)/(L*L)/
     &(L*L-3.d0)+(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))/
     &(L*L-3.d0)


	enddo

	
	do i=1,n
	  write(32,*) r(i),tth(i,1)
	enddo
	close(unit=32)	
	

100	call gradient

      
c     write(*,*) 'Deborah for the problem is', de(1)
	call residual
	call jacobian
      
	iter=iter+1
	do i=1,nvars
	  do j=1,nvars
	    a(i,j)=djacob(i,j)
	  enddo
	  
	  b(i,nrhs)=resid(i)
	enddo
c      write(*,*) nvars
	call dgesv(nvars,nrhs,a,nvars,ipiv,b,nvars,info)

c	write(*,*) gre
	
	call newvalues(b,xvar,nvars,nrhs)
c	write(*,*) iter

	call iteration(resid,nvars,gre)
      
            
c	  write(*,*) gre 


	if(gre.LE.(0.0000001d0)) then 
	  go to 200
	else
	  go to 100	

	endif
                  
      
200	do i=1,n
	  write(20,*) r(i),xvar(i)
	  write(30,*) r(i),xvar(i+n)
	
	  if (nmodes.EQ.1) then
	    write(40,*) r(i),xvar(i+2*n)
	    write(41,*) r(i),xvar(i+3*n)
c	    write(42,*) r(i),xvar(i+4*n)*de(1)/10.d0+1.d0/de(1)/10.d0) 
            write(42,*) r(i),xvar(i+4*n)
		write(43,*) r(i),(L*L-3.d0)/
     & (L*L-2.d0*xvar(i+2*n)-xvar(i+4*n))
	  endif
	enddo
		
	
	do i=1,nvars
	  write(11,*) xvar(i)                                                                                                                                                                                                             
	enddo	
	close(unit=20)
      close(unit=30)
      close(unit=40)
      close(unit=41)
      close(unit=42)
	close(unit=43)
      close(unit=11)
c	call meanvelocity(r,xvar,n,nvars,Re,delta,velmean)
c	write(*,*) 'delta,Re(d/R1)^1/2', delta, Re*dsqrt(delta)
c	write(*,*) 'w.r.t mean vel, Re(d/R1)^1/2',
c     &	Re*velmean*dsqrt(delta)	
c	write(*,*) 'uin=',uin,'uoutvelocity=',uout
	return
	end 
		

c	---------------END OF MAIN PROGRAM-------------------

c       --------------SUBROUTINE MEANVELOCITY---------------

      subroutine meanvelocity(rpt,x,ncheb,nv,Re,del,meanvel)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer i,j,k,ncheb,nv

      double precision rpt(ncheb),x(nv),velcoeff(ncheb),
     &  CKJ(ncheb,ncheb),c(ncheb),meanvel,sum,Re,del


      data pi/3.1415926535897930d0/

      
	c(1)=2.d0;
      c(ncheb)=2.d0;

      do i=2,ncheb-1
        c(i)=1.d0
      enddo

      do k=1,ncheb
        do j=1,ncheb
          CKJ(k,j)=2.d0*dcos(pi*dfloat(j-1)*dfloat(k-1)
     &            /dfloat(ncheb-1))/(dfloat(ncheb-1)*c(j)*c(k))
        enddo
      enddo
      
	
	do k=1,ncheb
        sum=0.d0
        do j=1,ncheb
          sum=sum+CKJ(k,j)*x(j)
        enddo
        velcoeff(k)=sum
      enddo
	sum=4.d0*velcoeff(1)
      
	
	do k=3,ncheb
        sum=sum+velcoeff(k)*(-dcos(dfloat(k)*pi)/dfloat(k)+
     &  dcos(dfloat(k-2)*pi)/dfloat(k-2)+dfloat(1/k)-dfloat(1/(k-2)))

      enddo
      meanvel=sum/4.d0
c	  write(*,*) 'mean velocity=',meanvel
c        write(*,*) 'Re w.r.t mean velocity=',Re*meanvel
c        write(*,*) 'Re(d/R1)^0.5=',del,Re*meanvel*dsqrt(del)
      
	
	return
      end
c       ------------END OF SUBROUTINE MEANVELOCITY-----------------



c	----------------SUBROUTINE INPUTVALUE ----------------

	subroutine input

	implicit double precision (a-h,o-z)
	parameter (iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
	parameter (nvars = (3*nmodes+2)*n)
	common/var/xvar(nvars),y(n),r(n),der(n,n),dder(n,n),
     &    dxvar(nvars),ddxvar(nvars)
	common/poly/lamdap(nrows)
	common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &	epsilonlam(nrows),alpha(nrows),delhlam(nrows),
     &	delhp(nrows),delhs,epsilons,betas,alph,si
	common/physical/visc,dens,cond
	common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &    gammadot
c	double precision L
	common/chain/L
	double precision lamdap,E
	data pi/3.1415926535897930d0/
		
	integer i,j,k
		
c	flow parameters----
c	Pe=0.d0
c	Pe = 20000.d0
c	Br = 0.000004019d0
c 	Br=0.0244d0
c	Re = 101.90d0
c	T0 = 298.15d0
c	T1 = 298.15d0
c	T2 = 298.15d0

c	delta=.0965
c	r21=1.0965
c	delta = 0.052631579d0  
c	r21 = 1.052631579d0    

c	taylor newtonian
c	delta=.209189842d0
c	r21=1.209189843d0
c	Pe=7975.157895d0*Re
c	Pe=11000.d0*Re

c	rheological parameters-------
c	betas=1.d0
c	betas = 0.79d0
c	delhs = 61790.d0
c	delhs=61664.940
c	epsilons = delhs/(8.314d0*T0)
c	epsilons=0.d0
	open(unit=10)
	read(10,*) Br
	read(10,*) T0
	read(10,*) T1
	read(10,*) T2
	read(10,*) delta
	read(10,*) r21
	read(10,*) betas
	read(10,*) delhs


c	mode parameters
	if (nmodes.NE.0) then
        open(unit=9)
	  do i=1,nmodes
c	    de(i)=38.d0
c	    de(i)=19.35d0 
c	    alpha(i)=0.d0
c	    delhlam(i)=61210.d0
c	    delhp(i)=61210.d0
c	    epsilonp(i)=delhp(i)/(8.314d0*T0)
c	    lamdap(i)=0.d0
c	    betap(i)=1-betas
c	    epsilonlam(i)=delhlam(i)/(8.314d0*T0)
c	    epsilonlam(i)=0.d0
c	    epsilonp(i)=0.d0
	    read(9,*) alpha(i)
	    read (9,*) epsilonp(i)
	    read(9,*) epsilonlam(i)
	    read(9,*) betap(i)
	    read(9,*) de(i)
	  enddo		
	endif
	
	

c	Pe=20000.d0
c	Pe=7975.157895d0*Re
c	Pe=0.d0
c	epsilons=delhs/(8.314d0*T0)
c	epsilons=0.d0
	read(10,*) si
	read(10,*) alph
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
	close(unit=10)
	close(unit=9)
c	write(*,*) Br, Re, Pe, betas,T0, T1, T2,epsilons
	
	
	return
	end

c	-------------END OF SUBROUTINE INPUTVALUE-----------

c	----------------SUBROUTINE INITIALGUESS-------------
	
	subroutine initialvalue(x,n,modes,rpt,temp1,temp2,ref,del,
     &                        u,de,uout,uin)
      
	implicit double precision (a-h,o-z)
	integer i,j,k,n,modes,nos
	double precision x(n),temp1,temp2,ref,rpt(n),del,A,B,eta
     &                ,ss,uout,uin
	common/chain/L
      
	
c      data L/1000.d0/
     
    
     
	double precision u(n),c11,c22,trt(n,1),tth(n,1),r1,r2,de

      La = 20.d0
	eta=1.d0/(1.d0+del)                           ! R1/R2
c	write(*,*),'remember to change boundary conditions for velocity'
c	rw=0.5d0
	rw=0.d0                                       ! Omig1/Omig2
	nos=n/(modes*3+2)
	
	A=(rw/eta-eta)/(1.d0+eta)
	B=(rw-1.d0)/(1.d0-eta)/(eta-1.d0/eta)
	
	do i=1,nos
	  x(i)=A*rpt(i)+B/rpt(i)                                                    ! Taylor-Couette flow, velocity 
	  u(i)=x(i)	
	  

	

	  
	  if(modes.EQ.1) then
                                                               ! Taylor-Couette flow, non-zero polymer stress, Oldroy-B model
        trt(i,1)= De*(-2.d0*B/rpt(i)**2)
     &*(1.d0/((L*L-3.d0)*sqrt(6.d0)*(1.d0/(2.d0*L*L*L))*De
     &*(-2.d0*B/rpt(i)**2)*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)
     &*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*(-2.d0*B/rpt(i)**2)+
     &dsqrt(1.d0+(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))
     &*(-2.d0*B/rpt(i)**2))**2))))))
	  
	  tth(i,1) = De*De*(8.d0*B**2/rpt(i)**4)
     &*(1.d0/((L*L-3.d0)*sqrt(6.d0)*(1.d0/(2.d0*L*L*L))*De
     &*(-2.d0*B/rpt(i)**2)*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)
     &*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*(-2.d0*B/rpt(i)**2)+
     &dsqrt(1.d0+(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))
     &*(-2.d0*B/rpt(i)**2))**2))))))
     &+(1.d0/((L*L-3.d0)*sqrt(6.d0)*(1.d0/(2.d0*L*L*L))*De
     &*(-2.d0*B/rpt(i)**2)*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)
     &*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*(-2.d0*B/rpt(i)**2)+
     &dsqrt(1.d0+(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))
     &*(-2.d0*B/rpt(i)**2))**2))))))

	  endif


	enddo
c	--------------isothermal solution-------------
c	r1=1/del
c	r2=r1+1.d0

c	c11=(1.d0/del)*(r1**2*dlog(r1)-r2**2*dlog(r2))/(r1**2-r2**2)
c	c22=-(0.5d0/del)*(1-r1**2*dlog(r1))-c11*r1**2/2.d0
	
c	do i=1,nos
c	  u(i)=(0.5d0/del)*(1.d0/rpt(i)-rpt(i)*dlog(rpt(i)))+
c    &	   c11*rpt(i)*0.5d0+c22/rpt(i)
c	enddo

	do i=1,nos
	  write(31,*) rpt(i),u(i), tth(i,1)
	enddo
c	---------------end of isothermal-------------

c	close(unit=31)

	do i=2,nos-1 
	  x(i+nos)=1.00d0                              ! Taylor-Couette flow, temperature 
	enddo
	 
c     read(1,*) (x(i+nos),i=2,nos-1)
c	 do k=1,nmodes
c	    read(2,*) (x(i),i=(k+1)*nos+1,(k+1)*nos+nos)
c		read(3,*) (x(i+modes*nos),i=(k+1)*nos+1,(k+1)*nos+nos)
c		read(4,*) (x(i+2*modes*nos),i=(k+1)*nos+1,(k+1)*nos+nos)                    
c	  enddo

        if (modes.NE.0) then                  ! Taylor-Couette flow, polymer stress
	  do 50 k=1,modes
		do 51 i=(k+1)*nos+1,(k+1)*nos+nos
		  x(i)=1.d0                ! T_rr
            x(i+modes*nos)=1.0d0                             ! T_rt
 		  x(i+2*modes*nos)=1.d0               ! T_tt
51		continue
		  
		  
50	  continue
	endif
c		if (i.eq.1) then
c	write(*,*) x((k+1)*nos+3*nos*modes-1)
c	write(*,*) n
c	endif



c	if(modes.EQ.1) then
c	  do i=1,nos
c		x(i+2*nos) =0.d0
c		x(i+3*nos)=trt(i,1)
c		x(i+4*nos)=tth(i,1)
c	  enddo
c	endif


c	---setting boundary values------
	x(1)=uin
c	x(nos)=0.5263d0
	x(nos)=uout
	x(nos+1)=temp1/ref
	x(nos+nos)=temp2/ref

c	do i=1,n
c	  read(60,*) x(i)
c	enddo


c	do i=1,n
c	  write(60,*) x(i)
c	enddo
	close(unit=60)
	
	
	return
	end
c	--------------END OF SUBROUTINE INITIALGUESS-----------

c	--------------SUBROUTINE NEWVALUES---------------------

	subroutine newvalues(rmat,x,dim,nrhs)
	implicit double precision (a-h,o-z)
	
	integer i,j,dim,nrhs
	double precision rmat(dim,nrhs),x(dim),sum
	
	
	do i=1,dim
	  x(i)=x(i)-rmat(i,nrhs)
	enddo
	
	
	return
	end
c	--------------END OF SUBROUTINE NEWVALUES--------------
	

c	-------------SUBROUTINE ITERATION-------------------	

	subroutine iteration(residue,numvar,gre)
	implicit double precision (a-h,o-z)

	integer i,j,numvar
	double precision residue(numvar),gre,temp
						

	gre=dabs(residue(1))
	temp=gre
	
	do i=2,numvar
	  if (dabs(residue(i)).GE.gre) then
	    gre=dabs(residue(i))
	  else
	    temp=gre
	    gre=temp
	  endif
	enddo
c	 write(*,*) 'gre for the problem is', gre
	
	return
	end

c	-------------END OF SUBROUTINE ITERATION------------

c	--------------SUBROUTINE MESHPOINT------------------

	subroutine meshpt(cheb,rpt,modes,ncheb,dr1)
	implicit double precision (a-h,o-z)

	integer modes,ncheb,i,j,k
	double precision cheb(ncheb),rpt(ncheb),dr1
	
	data pi/3.1415926535897930d0/
		

	do 1 j=1,ncheb
	  cheb(j) = dcos(pi*dfloat(j-1)/dfloat(ncheb-1))
1	continue
	
	do 2 j=1,ncheb
	  rpt(j) =(1.d0-cheb(j))/2.d0+(1.d0/dr1)
2	continue
	
	
	return
	end
c	-----------END OF SUBROUTINE MESHPOINT--------------
	
c	-----------SUBROUTINE DERIVATIVE---------------
	
	subroutine derivmat(dmatrix,cheb,ncheb)
	
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
     &		           (cheb(i)-cheb(j)))
	
		elseif((i.GE.2).and.(i.EQ.j).and.(j.LE.(ncheb-1))) then
		  dmatrix(i,j)= -1.d0*cheb(j)/(2.d0*(1.d0-cheb(j)**2))
		
		elseif ((i.EQ.1).and.(j.EQ.1)) then
		  dmatrix(i,j)=(2.d0*dfloat((ncheb-1)**2)+1.d0)/6.d0
		
		elseif ((i.EQ.ncheb).and.(j.EQ.ncheb)) then
		  dmatrix(i,j)=(-1.d0)*(2.d0*dfloat((ncheb-1)**2)
     &		           +1.d0)/6.d0
		endif
	  enddo
	enddo


	return
	end



c	----------------END OF SUBROUTINE DERIVATIVE-----------------

c	--------------subroutine square matrix----------------------
	
	subroutine squaremat(mat,ncheb,sqmat)
      implicit double precision (a-h,o-z)

	integer i,j,k,ncheb
	double precision sum,mat(ncheb,ncheb),sqmat(ncheb,ncheb)
	data pi/3.1415926535897930d0/	

	do 4 i=1,ncheb
	  do 5  k=1,ncheb
	    sum=0.d0
		
		do 6 j=1,ncheb
		  sum=sum+mat(i,j)*mat(j,k)
6		continue
		
		sqmat(i,k)=sum
5	  continue
4	continue

c	do i=1,ncheb
c	  do j=1,ncheb
c	    write (*,*) mat(i,j), sqmat(i,j)
c	  enddo
c	enddo
	
	
	return
	end	

c	---------------end of squarematrix subroutine------------------

	subroutine gradient
c	
	implicit double precision (a-h,o-z)
	parameter (iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
	parameter (nvars = (3*nmodes+2)*n)
	common/var/xvar(nvars),gre,y(n),r(n),der(n,n),dder(n,n),
     &           dxvar(nvars),ddxvar(nvars)
	common/poly/lamdap(nrows)
	common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &	epsilonlam(nrows),alpha(nrows),delhlam(nrows),
     &	delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond
	common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &	gammadot
c	double precision L
c	common/chain/L
	real*8 lamdap
	integer i,j,k,l,m	
	common/grad/du(n),dsqu(n),dtemp(n),ddtemp(n),
     &    dtrr(nmodes,n),dtrth(nmodes,n),dthth(nmodes,n)
	common/matrix/resid(nvars),sumstress(nmodes),a(nvars,nvars),
     &    djacob(nvars,nvars),djacobinv(nvars,nvars),b(nvars,nvars)
	double precision sumu1,sumu2,sumt1,sumt2,dfirst(n,n),
     &    dsecond(n,n)
	double precision sum(nmodes)	
      data pi/3.1415926535897930d0/
		
	
c	chain rule for differentiation from physical space to chebyshev space
c	 r(j) =(1.d0-y(j))/2.d0+(1.d0/delta)
c 	multiply first derivative matrix by -2,
c 	multiply second derivative matrix by 4
	
	do 11 i=1,n
	  do 12 j=1,n
	    dfirst(i,j)=der(i,j)*(-2.d0)
	    dsecond(i,j)=dder(i,j)*4.d0
12	  continue
11	continue

100	do 13 i=1,n
	  sumu1=0.d0
	  sumu2=0.d0
	  sumt1=0.d0
	  sumt2=0.d0
c	  sumtrh=0.d0
	  
	  do 14 j=1,n
		sumu1=sumu1+dfirst(i,j)*xvar(j) 
		sumu2=sumu2+dsecond(i,j)*xvar(j)
		sumt1=sumt1+dfirst(i,j)*xvar(j+n)
c		sumtrh=sumtrh+dfirst(i,j)*xvar(j+3*n)
		sumt2=sumt2+dsecond(i,j)*xvar(j+n)		
14	  continue
	  
	  du(i)=sumu1                          ! first derivative of Uth
	  dsqu(i)=sumu2                        ! second derivative of Uth
	  dtemp(i)=sumt1                       ! first derivative of Temperature
	  ddtemp(i)=sumt2                      ! second derivative of Temperature
c	  dtrth(i)=sumtrh
13	continue
	
	do 15 i=1,n
	  dxvar(i)=du(i) 
	  dxvar(i+n)=dtemp(i)
	  ddxvar(i)=dsqu(i)
	  ddxvar(i+n)=ddtemp(i)
15	continue

c 	----1st and 2nd derivative of trr,trth,tthth---------------
	
	if (nmodes.NE.0) then 
	  do 16 k=1,nmodes
		do 17 i=1,n
		  do 18 l=1,3
		    sumstress(k)=0.d0
		 	sum(k)=0.d0
			
			do 19 j=1,n
	 		  sumstress(k)=sumstress(k)+dfirst(i,j)*
     &	 	 	    xvar((k+1)*n+(l-1)*nmodes*n+j)	
			  sum(k)=sum(k)+dsecond(i,j)*xvar((k+1)*
     &				n+(l-1)*nmodes*n+j)
19			continue	
			
			dxvar(i+(k+1)*n+(l-1)*nmodes*n)=sumstress(k)
			ddxvar(i+(k+1)*n+(l-1)*nmodes*n)=sum(k)
			
18		  continue
17		continue
16	  continue
	endif

      return
	end
c	-------------end of subroutine gradient---------------------
c	-------------subroutine residual--------------------------
	subroutine residual

	implicit double precision (a-h,o-z)
	parameter (iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
	parameter (nvars = (3*nmodes+2)*n)
	common/var/xvar(nvars),gre,y(n),r(n),der(n,n),dder(n,n),
     &    dxvar(nvars),ddxvar(nvars)
	common/poly/lamdap(nrows)
	common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &	epsilonlam(nrows),alpha(nrows),delhlam(nrows),
     &	delhp(nrows),delhs,epsilons,betas
	common/physical/visc,dens,cond
	common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2,uout,rot_p,uin,
     &    gammadot
c	double precision L
	common/chain/L
	double precision lamdap
	integer i,j,k,l,m	
	common/grad/du(n),dsqu(n),dtemp(n),ddtemp(n),
     &    dtrr(nmodes,n),dtrth(nmodes,n),dthth(nmodes,n)
	common/matrix/resid(nvars),sumstress(nmodes),a(nvars,nvars),
     &    djacob(nvars,nvars),djacobinv(nvars,nvars),b(nvars,nvars)
	double precision sterm1(n),sum1,sum2,sterm2(n)
	
      data pi/3.1415926535897930d0/
	
c	data L/1000.d0/

c	---boundary conditions
	resid(1)=xvar(1)-uin
c	resid(n)=xvar(n)-0.5d0/0.95d0      
	resid(n)=xvar(n)-uout

	resid(n+1)=xvar(n+1)-T1/T0
	resid(n+n)=xvar(n+n)-T2/T0	

c	Neumann
c	resid(n+1)=dxvar(n+1)-0.d0
c	resid(n+n)=-dxvar(n+n)-10.d0*(xvar(n+n)-(T0-1.d0)/T0)-0.d0

	if (nmodes.EQ.0) then
	  do 20 i=1,n	
	  	sterm1(i)=0.d0
	  	sterm2(i)=0.d0
20	  continue
	else	
	  do 21 i=1,n
		sum1=0.d0
		sum2=0.d0
		
		do 22 j=1,nmodes


c       a)	                                
	 sum1=sum1+(betap(j)/de(j))*(2.d0*(1.d0/((1.d0-y(i))/2.d0+                                         ! SIGMA betap*(dTrt/dr+2Trt/r)      for Uth- equation
     &(1.d0/delta)))*(L*L-3.d0)/(L*L-2.d0*xvar(i+(j+1)*n)-
     &xvar(i+(j+1)*n+2*n*nmodes))*xvar(i+(j+1)*n+n*nmodes)+
     &xvar(i+(j+1)*n+n*nmodes)*(2.d0*dxvar(i+(j+1)*n)+
     &dxvar(i+(j+1)*n+2*n*nmodes))*(L*L-3.d0)/(L*L-2.d0*xvar(i+(j+1)*n)
     &-xvar(i+(j+1)*n+2*n*nmodes))/(L*L-2.d0*xvar(i+(j+1)*n)-
     &xvar(i+(j+1)*n+2*n*nmodes))+dxvar(i+(j+1)*n+n*nmodes)*(L*L-3.d0)/
     &(L*L-2.d0*xvar(i+(j+1)*n)-xvar(i+(j+1)*n+2*n*nmodes)))  
     
c       sum1 = 0.d0

c     							  	write(*,*) 'sum1=', 
c     &2.d0*(1.d0/((1.d0-y(i))/2.d0+(1.d0/delta)))*(L*L-3.d0)/
c     &(L*L-2.d0*xvar(i+(j+1)*n)-xvar(i+(j+1)*n+2*n*nmodes))
c     &*(betap(j)/de(j))*xvar(i+(j+1)*n+n*nmodes)

c      b) (in terms of new f(c))
c      sum1=sum1+(betap(j)/de(j))*2.d0*(1.d0/((1.d0-y(i))/2.d0+
c     &(1.d0/delta)))*
c     &(L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))
c     &+(betap(j)/de(j))*xvar(i+(j+1)*n+
c     &n*nmodes)*(2.d0*dxvar(i+(j+1)*n)+dxvar(i+(j+1)*n+2*n*nmodes))
c     &*(1.d0/(L*L-3.d0))*
c     &(L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))
c     &*(L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))
c     &+(betap(j)/de(j))*dxvar(i+(j+1)*n+n*nmodes)*
c     &(L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))  

c						      			write(*,*) 'ss(i)**2=', 
c     &xvar(i)
     
 	
c			           write(*,*) 'xvar(i+(k+1)*n+2*n*nmodes)=', 
c     & dxvar(i)

c      f(r) = 1
  	 
c		 sum1=sum1+(betap(j)/de(j))*(2.d0*xvar((j+1)*n                                                      ! SIGMA betap*(dTrt/dr+2Trt/r)      for Uth- equation
c     & +nmodes*n+i)/((1.d0-y(i))/2.d0+ 
c     &(1.d0/delta)))*1.d0+(betap(j)/de(j))*1.d0*
c     & dxvar((j+1)*n+nmodes*n+i)


c       a)
        sum2=sum2+(betap(j)/de(j))*(L*L-3.d0)/(L*L-2.d0*xvar((j+1)*n+i)
     &-xvar((j+1)*n+2*n*nmodes+i))*(dxvar(i)
     &-xvar(i)/((1.d0-y(i))/2.d0+(1.d0/delta)))*xvar((j+1)*n+
     &n*nmodes+i)



c      b)
	
c	  sum2=sum2+(betap(j)/de(j))*     
c     &(L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(j)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))*
c     &(dxvar(i)
c     &-xvar(i)/((1.d0-y(i))/2.d0+(1.d0/delta)))*xvar((j+1)*n+
c     &n*nmodes+i)
	      
c		  sum2=sum2+betap(j)*(xvar((j+1)*n+nmodes*n+i)*               ! SIGMA betap*Trt*(dUth/dr-Uth/r)   for Temp- equation
c     &             dxvar(i)-xvar((j+1)*n+nmodes*n+i)*
c     & 			 xvar(i)/((1.d0-y(i))/2.d0+(1.d0/delta)))
22		continue
		
		sterm1(i)=sum1
		sterm2(i)=sum2
21	  continue
	endif
	
	
	do 23 i=2,n-1
	
	  resid(i) =                                                      ! betas*exp(eps_s*(1/Tem-1))*( dUth/dr/r-Uth/r**2+d2Uth/dr2                                                                                                                                               
     &	 betas*dexp(epsilons*(1.d0/xvar(i+n)-1.d0))*(dxvar(i)/        !                         -eps_s/Tem**2*(dUth/dr-Uth/r)*dTem/dr )
     &     ((1.d0-y(i))/2.d0+(1.d0/delta))-xvar(i)/((1.d0-y(i))/2.d0+   ! + SIGMA betap*(dTrt/dr+2Trt/r)
     &     (1.d0/delta))**2+ddxvar(i)-(epsilons/(xvar(i+n))**2)*
     &	 (dxvar(i)-(xvar(i)/((1.d0-y(i))/2.d0+(1.d0/delta))))*
     &	 dxvar(i+n))+sterm1(i)	


	
	  resid(i+n)=ddxvar(i+n)+(dxvar(i+n)/((1.d0-y(i))/2.d0+           ! d2Tem/dr2 + dTem/dr/r + 
     &	 (1.d0/delta)))+Br*(betas*dexp(epsilons*                      ! Br*(exp(eps_s*(1/Tem-1))*(dUth/dr-Uth/r)**2 +
     &     ((1.d0/xvar(i+n))-1.d0))*((dxvar(i))**2-(2.d0*xvar(i)*       !     SIGMA betap*Trt*(dUth/dr-Uth/r) )  
     &	 dxvar(i)/((1.d0-y(i))/2.d0+(1.d0/delta)))+(xvar(i)**2/
     &	 ((1.d0-y(i))/2.d0+(1.d0/delta))**2))+sterm2(i)) 
          
23	continue


	
	if (nmodes.NE.0) then
	  do 24 k=1,nmodes
		do 25 i=1,n

c     a)      

        resid(i+(k+1)*n)=xvar(i+(k+1)*n)-(L*L-(2.d0*xvar(i+(k+1)*n)+
     &  xvar(i+(k+1)*n+2*n*nmodes)))*(1.d0/(L*L-3.d0))

c     b)
c         resid(i+(k+1)*n)=xvar(i+(k+1)*n)-
c     &(1.d0/((L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))))
c
c		  resid(i+(k+1)*n)=xvar(i+(k+1)*n)+                           ! Trr + De*exp(eps_l*(1/Tem-1))*exp(-eps_p*(1/Tem-1))*
c     & 		 de(k)*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))*         !       alpha*(Trr*Trr+Trt*trt)
c     &		 (dexp(-epsilonp(k)*(1.d0/xvar(i+n)-1.d0))*alpha(k)*
c     &       	 (xvar(i+(k+1)*n)**2+xvar(i+(k+1)*n+nmodes*n)**2))

                               

c		     resid(i+(k+1)*n+nmodes*n)=xvar(i+(k+1)*n+nmodes*n)-         ! Trt - De*exp(eps_l*(1/Tem-1))*(Trr*dUth/dr-Trr*Uth/r)
c     &         de(k)*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))*         !     + De*exp(eps_l*(1/Tem-1))*exp(-eps_p*(1/Tem-1))*alpha*
c     &         (xvar(i+(k+1)*n)*dxvar(i)-(xvar(i+(k+1)*n)*              !           (Trt*Trr+Trt*Ttt)   
c     &		 xvar(i)/((1.d0-y(i))/2.d0+(1.d0/delta))))+               !     - exp(-eps_p*(1/Tem-1))*(dUth/dr-Uth/r)
c     &		 de(k)*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))*
c     & 		 dexp(-epsilonp(k)*(1.d0/xvar(i+n)-1.d0))*
c     &		 alpha(k)*(xvar(i+(k+1)*n+nmodes*n)*
c     & 		 xvar(i+(k+1)*n)+xvar(i+(k+1)*n+nmodes*n)*
c     &		 xvar(i+(k+1)*n+2*nmodes*n))-dexp(epsilonp(k)*
c     &		 (1.d0/xvar(i+n)-1.d0))*(dxvar(i)-xvar(i)/
c     &		 ((1.d0-y(i))/2.d0+(1.d0/delta)))

c      a)

      resid(i+(k+1)*n+nmodes*n)=xvar(i+(k+1)*n+n*nmodes)-
     &dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))*de(k)
     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
     &*xvar(i+(k+1)*n)
     &*(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
     &/(L*L)   

c      b)

c      resid(i+(k+1)*n+nmodes*n)=xvar(i+(k+1)*n+n*nmodes)-
c     &dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))*de(k)
c     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
c     &*(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
c     &/(L*L)/(L*L-3.d0)

c      c)

c      resid(i+(k+1)*n+nmodes*n)=xvar(i+(k+1)*n+n*nmodes)-
c     &dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))*de(k)
c     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(L*L-3.d0)*(L*L-3.d0)/(L*L)/(L*L)
c     &*(1.d0/((L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))))
c     &*(1.d0/((L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))))

c	write(*,*) ' pet', 
      
c     &(L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))
     &

c        write(*,*) 'peterlin function is', (L*L-3.d0)/
c     & (L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
		
c		    resid(i+(k+1)*n+2*nmodes*n)= xvar(i+(k+1)*n+2*             ! Ttt - 2*De*exp(eps_l*(1/Tem-1))*(Trt*dUth/r-Trt*Uth/r)
c     &		 nmodes*n)-2.d0*de(k)*dexp(epsilonlam(k)*(1.d0/           !     + De*exp(eps_l*(1/Tem-1))*exp(-eps_p*(1/Tem-1))*alpha*
c     &		 xvar(i+n)-1.d0))*(xvar(i+(k+1)*n+nmodes*n)*                           (Trt*Trt+Ttt*Ttt)
c     &		 dxvar(i)-(xvar(i+(k+1)*n+nmodes*n)*xvar(i)/
c     &		 ((1.d0-y(i))/2.d0+(1.d0/delta))))
c     &		 +de(k)*
c     &		 dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))*
c     &		 dexp(-epsilonp(k)*(1.d0/xvar(i+n)-1.d0))*alpha(k)*(
c     &		 xvar(i+(k+1)*n+nmodes*n)**2+xvar(i+(k+1)*n+2*
c     &		 nmodes*n)**2)


c      a)
      resid(i+(k+1)*n+2*nmodes*n)=xvar(i+(k+1)*n+2*n*nmodes)-
     &2.d0*de(k)*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))
     &*xvar(i+(k+1)*n+n*nmodes)
     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
     &*(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
     &/(L*L)
     &-(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
     &/(L*L-3.d0)



c      b)

c     resid(i+(k+1)*n+2*nmodes*n)=xvar(i+(k+1)*n+2*n*nmodes)-
c    &2.d0*de(k)*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))
c     &*(1)*de(k)
c     &*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))
c     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
c     &*(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
c     &*(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
c     &/(L*L)/(L*L)/(L*L-3.d0)
c     &-(L*L-2.d0*xvar(i+(k+1)*n)-xvar(i+(k+1)*n+2*n*nmodes))
c     &/(L*L-3.d0)


c      c)
c      resid(i+(k+1)*n+2*nmodes*n)=xvar(i+(k+1)*n+2*n*nmodes)-
c     &2.d0*de(k)*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))
c     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(L*L-3.d0)*(1.d0/(L*L))*
c     &*(1.d0/((L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))))
c     &*(1.d0/((L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))))
c     &*(1.d0/((L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))))
c     &*(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(L*L-3.d0)*(L*L-3.d0)/(L*L)/(L*L)
c     &*de(k)*dexp(epsilonlam(k)*(1.d0/xvar(i+n)-1.d0))-
c     &(1.d0/((L*L-3.d0)*dsqrt(6.d0)*(1.d0/(2.d0*L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))
c     &*(1.d0/sinh((1.d0/3.d0)*dlog(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*
c     &(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i))+dsqrt(1.d0+
c     &(3.d0*dsqrt(6.d0)*(1.d0/2.d0)*(L*L-3.d0)*(1.d0/(L*L*L))*de(k)*
c     &(-xvar(i)/((1.d0-y(i))/2.d0+1.d0/delta)+dxvar(i)))**2))))))

25		continue
24	  continue

	endif
        

130 	continue
      
c	do i=1,n+(n+1)*n+2*n*nmodes
c      if (resid(i) .NE. resid(i)) then
c	write(*,*) resid(i)
c	endif
c	enddo
	return
	end

c	--------------end of subroutine residual------------------
c	-------------subroutine jacobian matrix----------------

	

	subroutine jacobian
		
	implicit double precision (a-h,o-z)
	parameter (iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
	parameter (nvars = (3*nmodes+2)*n)
	common/var/xvar(nvars)
c	double precision L
c	common/chain/L
	common/matrix/resid(nvars),sumstress(nmodes),a(nvars,nvars),
     &    djacob(nvars,nvars),djacobinv(nvars,nvars),b(nvars,nvars)
	double precision tempx(nvars),incr,temresid(nvars)
	integer i,j,k
	
	
	incr=0.0000001d0
	do i=1,nvars
	  tempx(i)=xvar(i)
	  temresid(i)=resid(i)
	enddo
	
	do 26 k=1,nvars
	  xvar(k)=xvar(k)+incr
	  call gradient
	  call residual
	  
	  do 27 i=1,nvars
	    djacob(i,k)=(resid(i)-temresid(i))/incr
27	  continue
	  
	  xvar(k)=xvar(k)-incr	
	  do i=1,nvars
	    resid(i)=temresid(i)
	  enddo
26	continue
	
	
	return
	end
c	------------end of subroutine jacobian------------


