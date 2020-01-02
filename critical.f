
c     ----------------------------------------------

	subroutine critical(alpha1,alpha2,recrit,rey,pivot,
     &                    velmean,sigmar,sigmai,q,elas)

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(iy=8,ny=2**iy+1,n=ny,nmodes=1,nrows=nmodes+1)
      parameter(nvars=(6*nmodes+5)*n,nrhs=1,nn=2*nvars)

      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &            epsilonlam(nrows),alp(nrows),delhlam(nrows),
     &            delhp(nrows),delhs,epsilons,betas
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2	
      common/wavenum/si,alpha
c	double precision L
	common/chain/L

      double precision eigenx(nn), eigeny(nn),temp,pivot(4),reynold,
     &   rey(4),slope,recrit,incre,decre,sigmar,sigmai,velmean,temp1(4),
     &   renext(40),alpha1,alpha2,rt1(12),q(10),elas
      integer i,j,eigennum,iter,step,num

c      write(*,*) ' Deborah is =', de(1)

20    do iter =1,2

        reynold=rey(iter)

        call inputval(reynold,si,alpha,elas)
	write(*,*) 'Deborah is', de(1)
        call taylor(velmean)
        call lsataylor(eigennum,q)
        
	  do i=1,eigennum
          read(14,*) eigenx(i),eigeny(i)
        enddo
        close(unit=14)
	
c	write(*,*) 'Reynolds number = ',reynold 
c       write(*,*)  'CHECK after writing fort.14 in critical.f'
c	stop

c       j=1
c       pivot(iter)=eigenx(1)
	  j=3
	  pivot(iter)=eigenx(3)

        temp=pivot(iter)

c       do i=2,eigennum
	  do i=4,eigennum

c	    if (eigenx(i).GE.10000) then
c	      pivot(iter)=eigenx(i+1)
c 	    endif
	
          if (eigenx(i).GE.pivot(iter)) then
	      pivot(iter)=eigenx(i)
   	      j=i
          else
            temp=pivot(iter)
            pivot(iter)=temp
          endif
	
          if (eigenx(i).GE.10000) then
            pivot(iter)=temp
          endif

        enddo
      enddo

c     do iter=1,2
c       write(*,*)
c       write(*,*) rey(iter),pivot(iter),j, eigeny(j)
c       write(*,*)
c     enddo

      if(pivot(1).LE.1.d-6.AND.pivot(2).GE.1.5d-6) then
        recrit=(rey(1)+rey(2))/2.d0
        goto 10
      endif


      if(pivot(1).GE.1.d-6) then
c       rey(1)=rey(1)-5.d0
c	  rey(1)=rey(1)-0.5d0
 	  rey(1)=rey(1)-2.d0

      endif
      
	
	if(pivot(2).LE.1.5d-6) then
c	  rey(2)=rey(2)+5.d0
c	  rey(2)=rey(2)+0.5d0
	  rey(2)=rey(2)+2.d0

      endif
      goto 20


10    call inputval(recrit,si,alpha,elas)
      call taylor(velmean)
      call lsataylor(eigennum,q)


      do i=1,eigennum
        read(14,*) eigenx(i),eigeny(i)
      enddo
      close(unit=14)
      
	

c     j=1
c     pivot(3)=eigenx(1)
	j=3	
	pivot(3)=eigenx(3)
      temp=pivot(3)
c     do i=2,eigennum

	do i=4,eigennum

        if (eigenx(i).GE.pivot(3)) then
          pivot(3)=eigenx(i)
          j=i
        else
          temp=pivot(3)
          pivot(3)=temp
          sigmar=eigenx(j)
          sigmai=eigeny(j)
        endif
      enddo


      if(pivot(3).LE.1.d-6) then
        rey(1)=recrit
        recrit=(rey(1)+rey(2))/2.d0
        goto 10
      endif

      
	if(pivot(3).GE.1.5d-6) then
        rey(2)=recrit
        recrit=(rey(1)+rey(2))/2.d0
        goto 10
      endif

	
	if(pivot(3).GE.1.d-6.AND.pivot(3).LE.1.5d-6) then
	  goto 90
	endif
	
90	return

	end






c	-----------------------------------------------
      subroutine inputval(Reynold,si,alph,elas)

      implicit double precision (a-h,o-z)
      parameter (iy=8,ny=257,n=ny,nmodes=1,nrows=nmodes+1)
      parameter (nvars = (3*nmodes+2)*n)
      common/var/xvar(nvars),y(n),r(n),der(n,n),dder(n,n),
     &           dxvar(nvars),ddxvar(nvars)
      common/poly/lamdap(nrows)
      common/rheo/epsilonp(nrows),betap(nrows),de(nrows),
     &            epsilonlam(nrows),alp(nrows),delhlam(nrows),
     &            delhp(nrows),delhs,epsilons,betas
      common/geo/Pe,Br,Re,r21,r1,d,delta,T0,T1,T2
c	common/chain/L
      double precision lamdap,alph,si,Reynold,elas
      data pi/3.1415926535897930d0/

      integer i,j,k
	if(Re.eq.0.d0) then
	if(elas.lt.1.d0)then	
	Re=Reynold                        ! Reynolds number
	else
	Re = Reynold/elas	! Reynold variable is the Deborah number
	endif
	else
	de(1) = Reynold
	endif
	
c     flow parameters----

c     rheological parameters-------

c       thermal parameters ------
      open(unit=10)
      write(10,*) Br                                ! Brinkman number
      write(10,*) T0                                ! reference temperature
      write(10,*) T1                                ! inner cylinder temperature
      write(10,*) T2                                ! outer cylinder temperature
      write(10,*) delta
      write(10,*) r21
      write(10,*) betas                             ! ratio of solvent viscosity to total viscosity
      write(10,*) delhs
      write(10,*) si                                ! azimuthal wavenumber
      write(10,*) alph                              ! axial wavenumber
      write(10,*) Re
      write(10,*) epsilons
      write(10,*) Pe

        close(unit=10)

c     mode parameters
        if(nmodes.EQ.1)then

        open(unit=9)
c	if(elas.lt.1.d0) then
c         de(1) = Reynold*elas
c	else
c	de(1) = Reynold
c	endif
        alp(1) = 0.d0
        delhp(1) = 0.d0
        delhlam(1) =0.d0
        lamdap(1) = 0.793d0
        betap(1) = 1.d0-betas
        epsilonp(1) = delhp(1)/(8.314d0*T0)  ! activation energy for polymer viscosity
        epsilonlam(1)=delhlam(1)/(8.314*T0) ! activation energy for polymer relaxation time
          write(9,*) alp(1)
          write(9,*) epsilonp(1)
          write(9,*) epsilonlam(1)
          write(9,*) betap(1)
          write(9,*) de(1)
        close(unit=9)
        endif

      if (nmodes.GT.1) then
        open(unit=11,FILE='polyparam.dat',STATUS='OLD')
        open(unit=9)
          do i=1,nmodes

         read(11,*) de(i)
        read(11,*) alp(i)
        read(11,*) delhp(i)
        read(11,*) delhlam(i)
        read(11,*) lamdap(i)
        read(11,*) betap(i)
        epsilonp(i) = delhp(i)/(8.314d0*T0)
        epsilonlam(i) = delhlam(i)/(8.314d0*T0)

          write(9,*) alp(i)
          write(9,*) epsilonp(i)
          write(9,*) epsilonlam(i)
          write(9,*) betap(i)
          write(9,*) de(i)
        enddo
        close(unit=11)
        close(unit=9)
      endif

	return
	end


