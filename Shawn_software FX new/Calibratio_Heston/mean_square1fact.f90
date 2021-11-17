	program mean_square
!! This program compute the minimization of the objective function defined by
!! the sum of the difference between the datum and the price computes with the
!!  formula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	include 'mpif.h'
	integer irank,nproc, ierr
!	use IMSL      
	integer ibtype,iprint,N,Nm,M,Max,ME,Nmpi,k0,k10
	integer NNEW,nobs,nobsp1,ini1
	parameter(NNEW=8)
	parameter(Nobsp1=400)
!!!!!!!!!!!!!!!!!1
!! Me: numereo di vincoli saturati
!! M: numeri di vincoli oltre upper bound e lower bound
!! in particulare 2\chi\theta/eps^2>1
    parameter (ibtype=0, iprint=3, M=1, Max=10, ME=0,N=14,nxdata=2**NNEW)
     
	   real*8 ertheta,ermax,xold(N),grad(N),alp,E,s0, rtax,Tm,aux
        real*8 hderiv,hstep,hs(N),xvero(N),vvv,tobs(Nobsp1),enormv
		real*8 option(Nobsp1), optobs(Nobsp1), ivett(2**NNEW),opt
		real*8 pl,vol, tod
	   integer nvar,jg
      real*8 ones(2**NNEW),kprinc(2**NNEW),erfun
    real*8 pi,erv1,erv2,erx,erz,X,Z,er,er1,er2,vvm,riskl(2)
	real*8 xguess(N),diff1,dtt,ddtt,s00
	real*8 xnew(N),flike,flikeold,fk,fk1,diag(N)
	real*8 xprim(30), xprim2,zprim2
	real*8 xscale(N),v0(2),z0,xlb(N),xub(N)
	real*8 mu,eps(2),theta(2),kkp(2),vinc(M)
	real*8 tau(Nobsp1), kprice(Nobsp1), rvalueop(Nobsp1)
	real*8  funzionale, valint
	complex*16  transl(2**NNEW), transl1(2**NNEW)
	real*8 val1(2**NNEW)
!!	complex*16 Ad(2**NNEW,2**NNEW), Bd(2**NNEW,2**NNEW)
	real*8 xdata1(2**NNEW),xdata2(2**NNEW),ydata(2**NNEW)
	!! New variables
      real*8 rho12(2),V1,V2
	  real*8 dgamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Aggiunte da me
	real*8 vals(N)
	real*8 t(Nobsp1),ror2

	integer ilinx,iliny,ier
	
	real a1,b1,c1,d1,ror1




        common/proce/nproc, irank,ierr
        common/var/k0,k10,vm0,vni0,vgei0,vm10,vni10,vgei10
      
!	call mpi_init(ierr)
!	call mpi_comm_rank(mpi_comm_world,irank,ierr)
!        call mpi_comm_size(mpi_comm_world,nproc,ierr)
!        call mpi_barrier(mpi_comm_world,ierr)
         ierr=0
         nproc=1
         irank=0
		 open(10,file='like1fact031105.txt',status='unknown') 
         open(15,file='rec1fact031105.txt',status='unknown')
!		 open(16,file='probcondleggi.txt',status='unknown')
	Nmpi=2**NNEW
	Nm=Nmpi
	pi=dacos(-1.d0)
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! initial parameters
        v0(1)=0.5d0
		v0(2)=0.005d0

    	mu=5.6d-2

		mu=0.3
        eps(1)=0.306d0
        eps(2)=0.05d0

        rho12(1)=-5.76d-1
        rho12(2)=-0.1d0


        kkp(1)=9.94d0
        kkp(2)=0.05d0

        theta(1)=0.5159d0
        theta(2)=0.8d0

		riskl(1)=-0.1d0
		riskl(2)=-0.5d0



    	xguess(1)=eps(1)
	    xguess(2)=theta(1)
	    xguess(3)=rho12(1)

!! kkp e' la chi degli appunti
	    xguess(4)=kkp(1)
	    xguess(5)=v0(1)

        xguess(6)=mu
		xguess(7)=riskl(1)
		
		
		xguess(8)=eps(2)
		xguess(9)=theta(2)
		xguess(10)=rho12(2)
		xguess(11)=kkp(2)
		xguess(12)=v0(2)
		xguess(13)=riskl(2)
			
         nvar=7


        enormv=0.d0
	    do j=1, nvar
	       xvero(j)=xguess(j)
		   enormv=enormv+xvero(j)**2
	       xguess(j)=xguess(j)*(1.d0+(-1.d0+rand()*2.d0)*0.5)
!	  xguess(j)=xvero(j)
	    enddo
		enormv=dsqrt(enormv)

         
!		   do j=1, nvar
!		     write(10,*) j,xguess(j)
!			 enddo
!         write(*,*) 'prima pausa'

	write(10,*) 'j=1-eps_1(vol of vol)'
	write(10,*) 'j=2-theta_1'
	write(10,*) 'j=3-rho12_1'
	write(10,*)'j=4- chi_1'
	write(10,*) 'j=5-v0_1'
	write(10,*)'j=6-mu (drift)'
	write(10,*)'j=7-lambda_1'
	
	diff1=2.d0*kkp(1)*theta(1)/eps(1)**2-1.d0
	
	
	 if(diff1.lt.0.d0) then
        write(*,*) 'negative diff1=',diff1
		stop
	 endif

	
           xguess(1)=  0.151760029191969d0     
           xguess(2)=  0.350687523477068d0     
           xguess(3)= -0.634077992295314d0     
           xguess(4)=   13.8323847581696d0             
            xguess(5)=  0.743014064941567d0       
           xguess(6)=  0.216770867118359d0 
		   xguess(7)= -0.102152598628821d0     
!           xguess(8)=  4.511707416188458d-2
!            xguess(9)=  0.609524599905211d0
              xguess(8)=10.d0
			  xguess(9)=0.0001d0
			  xguess(11)=90.d0     
           xguess(10)= -0.124448095518608d0     
!           xguess(11)=9.521945336786091d-04
           xguess(12)=  6.266012310524148d-3 
		   xguess(13)= -0.539862111347517d0
		   
!!!!!!!!!!!!!!!!! optimal point
!            xguess(1)=  0.162075039973550     
!            xguess(2)=  1.243207326437537E-002
!            xguess(3)= -0.635992120641391     
!            xguess(4)=   13.8467130933956     
!            xguess(5)=  4.657581898049529E-002
!            xguess(6)=  1.124933267467420E-002
!            xguess(7)= -8.305477783607002E-002

!           xguess(1)=  0.180413870501784d0
!		   xguess(2)=  2.037230512396974d-2
!		   xguess(3)= -0.640888190849594d0
!		   xguess(4)= 13.8465164919938d0
!		   xguess(5)=  1.215112091475418d-2
!		   xguess(6)=  1.667326436215557d-2
!		   xguess(7)= -8.326359747085277d-2
		   

!!! griglie di punti ydata griglia per x, xdata griglia per v, zdata griglia per z
 
	xdata1(1:2**NNEW)=0.d0
	xdata2(1:2**NNEW)=0.d0
	ydata(1:2**NNEW)=0.d0
     
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! operazione di filtraggio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Numero osservazioni: nobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Attenzione Nobsp1=nobs+ini1       ATTENZIONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!! mese novembre nobs=303	  			 
!!!! Reading the data set November 3, 2005
!!! S&P500 index value=1214.76
		    open(43,file='november031105.txt')
			nobs=303
!             nobs=281
!! mese giugno reading 3, June 2005
!! dati 281 S&P500 value=1196.02
!					    open(43,file='giugno030605.txt')    
!! dati 281 S&P500 value=1196.02
!					    open(43,file='gennaio030105.txt')


         
		 rewind(43)
!		 s0=1214.76d0

		 rtax=0.03d0
		 

		    do j=1, nobs
               read(43,*) tod, s0,kprice(j),tau(j),rvalueop(j),pl,vol
			enddo




       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! parametri per 2^8 punti
!         erx=300.d0
!		 erv1=600.d0

!         erx=300.d0
!		 erv=600.d0
		 
		 erx=200.d0
         erv1=300.d0
		 erv2=300.d0

	  NN=Nm-1

!    V2=pi/(erv2)*(NN+1)
!	V1=pi/(erv1)*(NN+1)
!	X=pi/(2.d0*erx)*(NN+1)
	
!	do j=1,Nm
!       xdata1(j)=dfloat(j-1)*V/dfloat(NN+1)
!       ydata(j)=-X+2.d0*X*dfloat(j-1)/dfloat(NN+1)
!       xdata2(j)=-Z+2.d0*Z/(NN+1)*(j-1)
!       write(*,*)'v',j,xdata(j)    
!	end do
!    do j=1,Nm
!	  write(*,*) 'x',ydata(j)
!	  enddo

	  	do j=1,NN+1
	      ones(j)=dfloat(j-1)
	    end do	


      do j=1,NN+1	   
       	ivett(j)=pi*dfloat(NN+1)*(-0.5+dfloat(j-1)/dfloat(NN+1))/erx
	  enddo
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Ciclo di minimizzazione
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   	   ermax=5.d-02
	  

	   hderiv=0.00001d0
	  ! hstep=0.01d0
         hstep=0.01d0
        
		flikeold=1000000.d0
		do j=1,nvar
		  hs(j)=hstep
		  enddo
	   do j=1,nvar
	      xold(j)=0.d0
	   enddo


             	   do 2322 iter=1,1000


        if(iter.lt.50) hstep=0.5
		if(iter.lt.200) hstep=0.1
		if(iter.lt.1000) hstep=0.001

		do j=1,nvar
		  hs(j)=hstep
		  enddo
              ivar1=1
			  ivar2=nvar	   

      CALL CPU_TIME ( time2_begin )


!!!! inizio calcolo funzione
23232        continue

            ! do j=1, nvar
			 !  write(*,*) j,xguess(j)
			  ! enddo
			   
			                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
                           
      CALL CPU_TIME ( time2_begin )

              
           flike=0.d0
           do jj=1, nobs
    		   E=kprice(jj)/s0
	      	   Tm=tau(jj)
               s00=1.d0
         rtax=xguess(6)  
		    
 
            call Hestonv0(nvar,xguess,E,rtax,Tm,opt)
            option(jj)=opt
        !       write(*,*) jj,option(jj),rvalueop(jj)/s0
     
			flike=flike+(option(jj)-rvalueop(jj)/s0)**2
			enddo
			flike=dsqrt(flike)

           !  write(*,*) 'funzione=',flike

      CALL CPU_TIME ( time2_end )

		write(*,*) 'time=',time2_end-time2_begin
	       


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Calcolata flike
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             write(*,*) 'funzione nel punto=',flike
	      
           !! erfun=dabs(flike-flikeold)/dabs(flike)
		    erfun=dabs(flike)
		if(erfun.le.ermax) go to 2323
          write(*,*) 'errfun=',erfun
		  write(10,*) 'erfun=',erfun
          write(*,*)'iter=',iter, flike,flikeold
          write(10,*)'iter=',iter,'funzione nuova e vecchia=',flike,flikeold

           do j=1,nvar
		     write(*,*) j,xguess(j),xold(j)
			 write(10,*) j,xguess(j),xold(j)
		   enddo
                
         
              if(flike.gt.flikeold) then
			      icont=icont+1
			     write(*,*) 'errato passo'

                do j=1,nvar
				   hs(j)=hs(j)*0.5
				   xguess(j)=xold(j)+hs(j)*grad(j)*diag(j)
				   enddo
				   flike=flikeold
				   go to 23232
                   endif
                 flikeold=flike
                 do j=1,nvar
				  xold(j)=xguess(j)
				  enddo


		   
           flikeold=flike




             write(*,*) 'inizio calcolo gradiente'
!			                 stop

!!! inizio calcolo gradiente
              
             do 2224 jg=ivar1,ivar2

               do j=1,nvar
			   xnew(j)=xguess(j)
			   enddo
			   xnew(jg)=xnew(jg)+hderiv


              
           fk=0.d0
           do jj=1, nobs
    		   E=kprice(jj)/s0
	      	   Tm=tau(jj)
               s00=1.d0
            rtax=xnew(6)
            call Hestonv0(nvar,xnew,E,rtax,Tm,opt)
            option(jj)=opt

			fk=fk+(option(jj)-rvalueop(jj)/s0)**2
			enddo
			
             fk=dsqrt(fk)




                   grad(jg)=-(fk-flike)/hderiv
				   write(*,*) jg,grad(jg)

2224          continue

               write(10,*) 'minus gradient'
			   do jg=1,nvar
			     write(10,*) jg,grad(jg)
               enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! matrice vincoli bound
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 	 alp=1.d0
					 diag(1)=xguess(1)**2
					 diag(2)=xguess(2)**2
					 diag(3)=(1.d0-xguess(3))**2*(xguess(3)+1.d0)**2
					 diag(4)=xguess(4)**2
					 diag(5)=xguess(5)**2
	                 diag(6)=1.d0
                     diag(7)=1.d0
					  do j=1,nvar
				      diag(j)=1.d0
					  enddo 

                     write(*,*) 'entro step minimizzazione'
!!!!! opero step di massimizzazione
23200               continue

                    do  2326 j=1,2
					
2325                if(hs(j).lt.1.d-18) go to 2326
                     

					   xnew(j)=xguess(j)+hs(j)*diag(j)*grad(j)
                          
!!! opero controllo sulle variabili
                    
					 if(xnew(j).lt.0.d0) then
					    hs(j)=hs(j)*0.5d0
						go to 2325
					 endif
2326                 continue
!! variabile rho12
                     do  2327 j=3,3
					  
23255                if(hs(j).lt.1.d-18) go to 2327
                     

					   xnew(j)=xguess(j)+hs(j)*diag(j)*grad(j)
                         
					 if(dabs(xnew(j)).gt.1.d0) then
					    hs(j)=hs(j)*0.5d0
						go to 23255
					 endif
2327                 continue                     


!!!! variabili chi, v0,eps (gamma),beta 
                    do  2328 j=4,5
					  
23256                if(hs(j).lt.1.d-18) go to 2328
                       xnew(j)=xguess(j)+hs(j)*diag(j)*grad(j)
                          
					 if(xnew(j).lt.0.d0) then
					    hs(j)=hs(j)*0.5d0
						go to 23256
					 endif
2328                 continue
                
                      xnew(6)=xguess(6)+hs(6)*diag(6)*grad(6)
					  xnew(7)=xguess(7)+hs(7)*diag(7)*grad(7)
                  	diff1=2*xguess(4)*xguess(2)/xguess(1)**2-1.d0 
!                    if(diff1.lt.0.d0) then
!					 do j=1,nvar
!					   hs(j)=hs(j)*0.5d0
!					 enddo
!   					  write(*,*)'violato vincolo non lineare'
!					  go to 23200
!					endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! passo effettuato
	                  					
			do j=1,nvar
			  xold(j)=xguess(j)
			  xguess(j)=xnew(j)
			enddo		
									                
2322  continue

2323  write(*,*) 'successful relative error function',erfun

      write(*,*) 'optimal point'
	  do j=1,nvar
	     write(*,*) j,xguess(j)
	  enddo

	  
      write(10,*) 'errore funzione=',erfun

	  write(10,*) 'optimal point'
	  do j=1,nvar
	     write(10,*) j,xguess(j)
	  enddo
!!!! printing the value in the optimal point
           do jj=1, nobs
    		   E=kprice(jj)/s0
	      	   Tm=tau(jj)
               s00=1.d0
         rtax=xguess(6)   
            call Hestonv0(nvar,xguess,E,rtax,Tm,opt)
            option(jj)=opt*s0
            write(15,925) Tm,s0, kprice(jj),rvalueop(jj),option(jj)
			enddo

925         format(5(1x,d12.6))

       stop
	end




      subroutine mmm(x,n,sol)
      integer n,j,k
        real*8 a(0:2*n-1)
	complex*16 pip
	complex*16 i
	complex*16 x(0:n-1),sol(0:n-1)
	real*8 pi,as
	i=(0,1)
	as=-1.d0
	pi=dacos(as)
	do j=0,n-1
	a(2*j) = Real(x(j))
	a(2*j+1) = Imag(x(j))
      end do
	call cdft(2*n, cos(pi/real(n)), sin(pi/real(n)), a)
	do j=0,n-1
	sol(j)= a(2*j)+i*a(2*j+1)
      end do
	return
	end

	!
      subroutine cdft(n, wr, wi, a)
      integer n, i, j, k, l, m
      real*8 wr, wi, a(0 : n - 1), wmr, wmi, wkr, wki, wdr, wdi, &
         ss, xr, xi
      wmr = wr
      wmi = wi
      m = n
      do while (m .gt. 4)
          l = m / 2
          wkr = 1
          wki = 0
          wdr = 1 - 2 * wmi * wmi
          wdi = 2 * wmi * wmr
          ss = 2 * wdi
          wmr = wdr
          wmi = wdi
          do j = 0, n - m, m
              i = j + l
              xr = a(j) - a(i)
              xi = a(j + 1) - a(i + 1)
              a(j) = a(j) + a(i)
              a(j + 1) = a(j + 1) + a(i + 1)
              a(i) = xr
              a(i + 1) = xi
              xr = a(j + 2) - a(i + 2)
              xi = a(j + 3) - a(i + 3)
              a(j + 2) = a(j + 2) + a(i + 2)
              a(j + 3) = a(j + 3) + a(i + 3)
              a(i + 2) = wdr * xr - wdi * xi
              a(i + 3) = wdr * xi + wdi * xr
          end do
          do k = 4, l - 4, 4
              wkr = wkr - ss * wdi
              wki = wki + ss * wdr
              wdr = wdr - ss * wki
              wdi = wdi + ss * wkr
              do j = k, n - m + k, m
                  i = j + l
                  xr = a(j) - a(i)
                  xi = a(j + 1) - a(i + 1)
                  a(j) = a(j) + a(i)
                  a(j + 1) = a(j + 1) + a(i + 1)
                  a(i) = wkr * xr - wki * xi
                  a(i + 1) = wkr * xi + wki * xr
                  xr = a(j + 2) - a(i + 2)
                  xi = a(j + 3) - a(i + 3)
                  a(j + 2) = a(j + 2) + a(i + 2)
                  a(j + 3) = a(j + 3) + a(i + 3)
                  a(i + 2) = wdr * xr - wdi * xi
                  a(i + 3) = wdr * xi + wdi * xr
              end do
          end do
          m = l
      end do
      if (m .gt. 2) then
          do j = 0, n - 4, 4
              xr = a(j) - a(j + 2)
              xi = a(j + 1) - a(j + 3)
              a(j) = a(j) + a(j + 2)
              a(j + 1) = a(j + 1) + a(j + 3)
              a(j + 2) = xr
              a(j + 3) = xi
          end do
      end if
      if (n .gt. 4) call bitrv2(n, a)
	end
!      

      subroutine bitrv2(n, a)
      integer n, j, j1, k, k1, l, m, m2, n2
      real*8 a(0 : n - 1), xr, xi
      m = n / 4
      m2 = 2 * m
      n2 = n - 2
      k = 0
      do j = 0, m2 - 4, 4
          if (j .lt. k) then
              xr = a(j)
              xi = a(j + 1)
              a(j) = a(k)
              a(j + 1) = a(k + 1)
              a(k) = xr
              a(k + 1) = xi
          else if (j .gt. k) then
              j1 = n2 - j
              k1 = n2 - k
              xr = a(j1)
              xi = a(j1 + 1)
              a(j1) = a(k1)
              a(j1 + 1) = a(k1 + 1)
              a(k1) = xr
              a(k1 + 1) = xi
          end if
          k1 = m2 + k
          xr = a(j + 2)
          xi = a(j + 3)
          a(j + 2) = a(k1)
          a(j + 3) = a(k1 + 1)
          a(k1) = xr
          a(k1 + 1) = xi
          l = m
          do while (k .ge. l)
              k = k - l
              l = l / 2
          end do
          k = k + l
      end do
      end




      subroutine funzHeston(tau,k,epsil,theta,rho12,k1,v0,NN,A)
	integer NN
	real*8 tau,pi,arg,k(1:NN+1),xip,nupu,v0(2)
	real*8 epsil(2),rho12(2),k1(2),theta(2),epsq
       complex*16 cden, expression,kapp
	complex*16 i,A(1:NN+1),mu,zita,lambda,tildev
	complex*16 sg,sd,sb,ep,atmt, cdif,tildem
	arg=-1.d0
	pi=dacos(arg)
	i=(0.d0,1.d0)
    
	
	do j=1,NN+1
  	  kapp=k(j)
        cden=2.d0-3.d0*i*kapp-kapp**2
       cdif=-cden*0.5d0
     A(j)=(1.d0,0.d0)

	 do nb=1,1

	       epsq=epsil(nb)**2
  	 nupu=2.d0*k1(nb)*theta(nb)/epsq

        mu=-0.5d0*(k1(nb)+epsil(nb)*rho12(nb)*(i*kapp-2.d0))
        
      zita=0.5d0*(4.d0*mu**2+epsq*(kapp**2+3.d0*i*kapp-2.d0))**0.5
      

	   ep=exp(-2.d0*tau*zita)
       sb=-mu+zita+(mu+zita)*ep
	   sg=1.d0-ep
	  
	   sd=(mu+zita)+(zita-mu)*ep

	   	  tildem=2.d0/(epsil(nb)**2*sg)

	   tildev=4.d0*zita**2*ep*v0(nb)/sb

        expression=nupu*(log(sb/(2.d0*zita))+(mu+zita)*tau)+&
         v0(nb)*sg*cdif/sb

                  
	   if(dreal(expression).lt.400.d0) then 
	   
	   A(j)=A(j)*(1.d0,0.d0)/exp(expression)
        
		else
		
		A(j)=(0.d0,0.d0)
		
		endif
        enddo
		A(j)=A(j)/cden
	enddo
      return
	end






       subroutine Hestonv0(nvar,xguesss,E,rtax,Tm,val)
         integer N1,nvar
		 parameter(N1=2**12)
	   real*8 val
	   real*8  v0(2)
	   real*8  E, s00,rtax,Tm, elog, pert, aux
	   
	    real*8 thetan(2),k1n(2),riskl(2),xguesss(nvar)
	    
		real*8 eps(2),t,V,Y,Z,arg,pi,mu1,fk1,xprim,vprim
	    integer j,jj,M1,j1,jxi,jin
	    integer N,NN
	    real*8 epq
	    real*8 er,er1,M,mu
	    complex*16 i,x2
	    complex*16 x1(0:N1-1),x3(1:N1)
	    complex*16 A(1:N1),solu2(0:N1-1),solu1(1:N1)
	    real*8 kp(1:N1)
	    real*8 ones(1:N1)
        real*8 rho12(2),theta(2),k1(2)
	    complex*16 caux,csum
                 
        arg=-1.d0
   	    pi=dacos(arg)
	    i=(0.d0,1.d0)

           eps(1)=xguesss(1)
		   thetan(1)=xguesss(2)
		   rho12(1)=xguesss(3)
           k1n(1)=xguesss(4)
		   v0(1)=xguesss(5)
		   mu1=xguesss(6)
		   riskl(1)=xguesss(7)
		   
		  ! eps(2)=xguesss(8)
		  ! thetan(2)=xguesss(9)
		  ! rho12(2)=xguesss(10)
          ! k1n(2)=xguesss(11)
		  ! v0(2)=xguesss(12)
		  ! riskl(2)=xguesss(13)
		   
		   	
	       k1(1)=k1n(1)+riskl(1)
		   theta(1)=thetan(1)*k1n(1)/(k1n(1)+riskl(1))
		   
	       k1(2)=k1n(2)+riskl(2)
		   theta(2)=thetan(2)*k1n(2)/(k1n(2)+riskl(2))
	NN=N1-1
	do j=1,NN+1
	ones(j)=dfloat(j-1)
	end do
		
	er=500.d0
	do j=1,NN+1
	kp(j)=-er+2.d0*(er/(NN+1))*(j-1)
		end do

!! Attenzione ora calcolo integrali al variare di v'

       
       
	  


     do iob=1, 1

!!! xprim =log(S/S_0) when S=S_0

        xprim=0.d0
	          
          call funzHeston(Tm,kp,eps,theta,rho12,k1,v0,NN,A)
	     
!!!       
!!! calcolo integrale con la formula dei rettangoli	
!!!!	   
	        csum=(0.d0,0.d0)

	  do jj=1,NN
         csum=csum+E**(i*kp(jj)-1)*cdexp(i*kp(jj)*(-xprim-mu1*Tm))*A(jj)*(kp(jj+1)-kp(jj))  
	  end do

      
         csum=csum*dexp(-rtax*Tm)*dexp(mu1*2.d0*Tm+2.d0*xprim)/(2.d0*pi)

       val=csum

!!! close loop iop
    enddo




       return
	  end

