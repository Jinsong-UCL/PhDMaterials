!  Calibratio.f90 
!
!  FUNCTIONS:
!  Calibratio - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Calibratio
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Calibratio

   
          implicit none
      integer np,npmat, nv

! np numero opzioni call and put l'indice varia da 0 a np
	parameter(np=5)
	parameter(npmat=18)

! nv e' il numero delle variabili
	parameter(nv=5+5+1+1)

      integer nptime, jj,j,istrike, neta, maxdist,L, itmax,ikk,ik,ip
	integer ij, iv,ii
    
    real*8 v0_1,eps_1,chi_1,theta_1
    real*8 v0_2,eps_2,chi_2,theta_2
    
	real*8 tau, pi, futprice  
	real*8  texpire_c_p(1:npmat,1:np), tempo(1:np)
	real*8 diag(nv), diagold(nv)
	
	real*8 rdS2(nv), rdS1(nv), eta(1) 
	real*8 xc_1, xc_1new, xc_2,xc2new,tt
	real*8  vol_1,vol_1new

	real*8   var1media
	real*8  alphaS1, beta
    real*8 param1(nv), param2(nv), alp
    real*8 rparam(nv),vparam(nv), delta
	
	real*8 dt, Tmax,Kstrike(1:np),t, S0, auxop

	real*8 optS1(1,npmat,1:np), optS2(1,npmat,1:np)

	real*8  callval(1:npmat,1:np),putval(1:npmat,1:np)
	real*8 optvalp(1:npmat,1:np), optvalc(1:npmat,1:np),  auxpp
      real*8 future_value(1:npmat), crncy_value(1:npmat)
      integer idx(1:npmat)
       integer iseed
	real*8  param(nv), paramnew(nv),somma,sommanew
	real*8 hstep,hderiv,aux, ermax
	real*8 grad(nv), paramold(nv), sommaold
    real*8 omega
	
      
      integer jalp  
      
    open(1,file='Call_March_2013_out_of_sample_67_84.txt',status='unknown')
       open(2,file='SP500_out_of_sample_67_84.txt',status='unknown')
          open(3,file='Put_March_2013_out_of_sample_67_84.txt',status='unknown')

       open(7,file='param_recH0.txt',status='unknown')

              open(5,file='out_sampleH0_call_67_84.txt',status='unknown')
          open(6,file='out_sampleH0_put_67_84.txt',status='unknown')  
      
!!
!! Call on SP&500 Time to maturity March 2013

!! np number of strike np=5
      

!! idx is an integer vector that contains the index of the data that we want to use in the calibration

!! npmat is the number of option daily observations npmat =212	
	

!! nptime=1
!! nptime is the number of days used in the calibration
!! idx(1) corresponds to set September 27
      
	
!! idx(1)=1 September 19, 2010
!      idx(1)=1
!! idx(1)=21 October 12, 2010
 !     idx(1)=21
!! idx(1)=36 Novembre 15,2010
!      idx(1)=36
!! 
!      idx(1)=1
!! SaBR
!      neta=1
!	eta(1)=0.d0	
             
!! the first line of the file contains the strike prices 
!! Note that call and put options have the same strike prices
        rewind(1)
	   read(1,*) (Kstrike(ii), ii=1,np)

    do j=1, npmat
	   read(1,*) (callval(j,ii), ii=1,np)
	enddo
	
       write(*,*) 'end of the first file'
      rewind(3)
	  read(3,*) 
	do j=1, npmat
	   read(3,*) (putval(j,ii), ii=1,np)
	enddo
      write(*,*) 'end of third file'   

      rewind(2)
	do j=1,npmat
	  read(2,*) future_value(j)
	enddo
     write(*,*) 'end of the second file'   
  !    rewind(11)
!	do j=1,npmat
!	  read(11,*)  crncy_value(j)
!	enddo
      



!! defining the day to the maturity
!! Our valuation starts April 2, 2012 and ends July 27, 2012. From March 19, 2012 to July 27, 2013 we 92 observations
!! A month is composed of 22 days. So that a year is 22*12=264

!!  maxdist e' il numero di giorni che intercorre dalla prima osservazione alla scadenza contato dalla serie storica messa
!! a disposizione da Bloomberg.

!!!!!!!!!! attenzione nel lognormal la variabile e' ln(S_t/S_0)

      maxdist=264-22
      do j=1, npmat
         do ii=1, np
	     texpire_c_p(j,ii)=dfloat(maxdist-j+1)/264.d0
        enddo
	enddo

!! Note that the unit is a day rember that a year is made by 261 daily observations on Bloomberg.


	
	L=2
	pi=dacos(-1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! parametri opzione
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
     
      ermax=1.d-05

      ermax=5.d-06
	       

! param(1)=delta   		 
! param(2) = eps_1 (vol of vol)
! param (3)= v0_1 (volatility)
! param(4)=theta_1
!param(5)=chi_1
! param (6)= rho_1 (correlation volatility - asset)
! param(7)= eps_2 (vol of vol)
! param(8) =v0_2 (volatility)
! param(9)=theta_2
!param(10)=chi_2
! param(11)= rho_2 (correlation coefficient)
!param(12)=Omega
!! Starting point

     
      
     nptime=1
      rewind(7)
      
      do jalp=1, 60
          
          do j=1, nv
              read(7,*) jj, param(jj)
          enddo
      enddo
      
      
      do 2011 jalp=1,npmat

     
         nptime=1
	     idx(1)=jalp
!      idx(2)=jalp+1
!      idx(3)=jalp+2
        
!! ciclo di minimizzazione
      
!cc calcolo funzione teorica
 144        call calcolocall(nv,np,npmat,nptime,idx,param,texpire_c_p,future_value,Kstrike,optvalc,optvalp)
!cc calcolo funzione
            somma=0.d0
            do j=1,nptime 
               do ik = 1, np
                  tau=texpire_c_p(j,ik)
  	              auxop=optvalc(j,ik)
	          !    somma=somma+(auxop-callval(idx(j),ik))**2
                           somma=somma+(auxop-callval(idx(j),ik))**2/callval(idx(j),ik)**2
                  auxpp=optvalp(j,ik)
	              somma=somma+(auxpp-putval(idx(j),ik))**2/putval(idx(j),ik)**2
                  futprice=future_value(idx(j))
        	      write(*,945) tau, future_value(idx(j)),auxop,callval(idx(j),ik),auxpp,putval(idx(j),ik)
               enddo	

                             write(5,947) (optvalc(j,ik),callval(idx(j),ik),ik=1,5)
                             write(6,947) (optvalp(j,ik),putval(idx(j),ik),ik=1,5)

            enddo
      somma=dsqrt(somma)
	write(*,*) 'vecchia=',sommaold,' nuova=',somma


945 format(6(1x,d12.6))
947 format(10(1x,f9.3))
    

112   continue
     
2011  continue
	
	    end program Calibratio


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! la routine che segue implementa la formula corretta ma numericamente
!!  e' meglio la formula implementata in L2
      subroutine calcolocall(nvv,npp,nppmat,npptime,idxx,par,tmp_c_p,fut,KKstrike,optc,optp)

! nvv: number of variables
! npp: number of call and put options
! par: vettore dei parametri
!!
! tmp_c_p: time to maturity  of call options
!! KKstrike: exercise prices of call options. 
!!           Note that call options must have the same exercise prices

      implicit none
      integer nvv,npp, nppmat,npptime, neta
      

      integer  ii,jj,j,istrike, itemp,  maxdist,L, itmax,ikk,ik,ip
      integer ij, iv, iteta
      real*8 tau, pi,vvv, delta  

      integer idxx(1:npptime)
	  real*8 par(nvv),tmp_c_p(1:nppmat,1:npp),fut(1:nppmat)
	  real*8 xx0, prod
	  real*8 KKstrike(1:npp),optc(1:nppmat,1:npp),optp(1:nppmat,1:npp)
      
      
         integer N1
		 parameter(N1=2**14)
	   real*8 val
	   real*8  v0(2)
	   real*8  E, s00,rtax,Tm, elog, pert, aux
	   
	    real*8 thetan(2),k1n(2),riskl(2)
	    real*8 texpire_c_p(nppmat,npp)
		real*8 eps(2),t,V,Y,Z,arg,mu1,fk1,xprim,vprim
	    integer M1,j1,jxi,jin
	    integer N,NN
	    real*8 epq
	    real*8 er,er1,M,mu
	    complex*16 i,x2
	    complex*16 x1(0:N1-1),x3(1:N1)
	    complex*16 Ac(1:N1),Ap(1:N1),solu2(0:N1-1),solu1(1:N1)
	    real*8 kp(1:N1)
	    real*8 ones(1:N1)
        real*8 rho12(2),theta(2),k1(2)
	    complex*16 caux,csum, csump
        real*8 omega
        complex*16 phiq
       
! param(1)=delta   		 
! param(2) = eps_1 (vol of vol)
! param (3)= v0_1 (volatility)
! param(4)=theta_1
!param(5)=chi_1
! param (6)= rho_1 (correlation volatility - asset)
! param(7)= eps_2 (vol of vol)
! param(8) =v0_2 (volatility)
! param(9)=theta_2
!param(10)=chi_2
! param(11)= rho_2 (correlation coefficient)
! param(12)=omega

          
        arg=-1.d0
   	    pi=dacos(arg)
	    i=(0.d0,1.d0)

           delta=par(1)
           eps(1)=par(2)
		   v0(1)=par(3)
           thetan(1)=par(4)
           k1n(1)=par(5)
		   
           rho12(1)=par(6)
           
		  
		   eps(2)=par(7)
           v0(2)=par(8)
		   thetan(2)=par(9)
           k1n(2)=par(10)
		   rho12(2)=par(11)
           omega=par(12)
	NN=N1-1
	do j=1,NN+1
	ones(j)=dfloat(j-1)
	end do
		
	er=500d0
	do j=1,NN+1
	kp(j)=-er+2.d0*(er/dfloat(NN))*(j-1)
		end do

!! Attenzione ora calcolo integrali al variare di v'

   do j=1,npptime
               do ik = 1, npp
         Tm=tmp_c_p(j,ik)
         
         E=KKStrike(ik)/fut(idxx(j))
         
!!! xprim =log(S/S_0) when S=S_0
        xprim=0.d0
	          
          call funzHeston(Tm,kp,eps,thetan,rho12,k1n,v0,delta,omega,NN,Ac, Ap)
	     
!!!       
!!! calcolo integrale con la formula dei rettangoli	
!!!!	   

      csum=(0.d0,0.d0)
      csump=(0.d0,0.d0)

	  do jj=1,NN
         csum=csum+E**(i*kp(jj))*cdexp(i*kp(jj)*v0(2)*Tm)*Ac(jj)*(kp(jj+1)-kp(jj))
           csump=csump+E**(i*kp(jj))*cdexp(i*kp(jj)*v0(2)*Tm)*Ap(jj)*(kp(jj+1)-kp(jj))  
	  end do

      
            csum=csum*dexp(v0(2)*Tm*2.d0)*dexp(-v0(2)*Tm)*fut(idxx(j))**2/(2.d0*pi*KKStrike(ik))
       !     csum=csum+(fut(idxx(j))-KKStrike(ik))*dexp(-v0(2)*Tm*0.5d0)

            !! Here put price formula  differs from the call price formula 
            csump=csump*dexp(-v0(2)*Tm*2.d0)*dexp(-v0(2)*Tm)*KKStrike(ik)**3/(2.d0*pi*fut(idxx(j))**2)
 !           csump=csump+(KKStrike(ik)-fut(idxx(j)))*dexp(-v0(2)*Tm*0.5d0)

       val=csum
       

	    optc(j,ik)=val 
         optp(j,ik)=csump 
       ! write(*,*) val, fut(idxx(j))-KKStrike(ik)

       enddo
   enddo
   
       return
    end

    
      subroutine funzHeston(tau,k,epsil,theta,rho12,k1,v0,deltaa,omeg,NN,Acc,App)
	  implicit none
      integer NN, nb
	  real*8 tau,pi,arg,k(1:NN+1),xip,nupu,v0(2),omeg
	  real*8 epsil(2),rho12(2),k1(2),theta(2),epsq, deltaa,psii
       complex*16 cden, expression,kapp
	  complex*16 i,Acc(1:NN+1),mu,zita,lambda,tildev
      complex*16 App(1:NN+1)
	   complex*16 sg,sd,sb,ep,atmt, cdif(2),tildem
    complex*16 cdenp, cdenc, phiq
       real*8 qq
       integer j
	   arg=-1.d0
	   pi=dacos(arg)
	   i=(0.d0,1.d0)
    
	psii=deltaa**2+2.d0*deltaa*rho12(1)+1.d0
    
	do j=1,NN+1
  	    
        kapp=k(j)
        
        
         cdenc=2.d0-3.d0*i*kapp-kapp**2
        cdenp=6.d0+5.d0*i*kapp-kapp**2
        
     
     !!! start call price computation   
        qq=2.d0
        
       phiq=0.5d0*(kapp**2+i*kapp*(2.0d0*qq-1.d0)-(qq**2-qq))
        
     
        
        cdif(1)=phiq*psii
        cdif(2)=phiq*omeg**2-qq+i*kapp
     
        !! imitialization for a product
        Acc(j)=(1.d0,0.d0)
        App(j)=(1.d0,0.d0)
        
	    do nb=1,1

	       epsq=epsil(nb)**2
           
    	    nupu=2.d0*k1(nb)*theta(nb)/epsq
        
        !! here define mu for both models
       
            if(nb.eq.1) then
        !! here define mu for both models
               mu=-0.5d0*(k1(nb)+epsil(nb)*(rho12(nb)+deltaa)*(i*kapp-qq))
            else
                   mu=-0.5d0*(k1(nb)+epsil(nb)*rho12(nb)*omeg*(i*kapp-qq))
            endif
            
        !! here we define zeta
        zita=0.5d0*(4.d0*mu**2+2.d0*epsq*cdif(nb))**0.5
      
	    ep=cdexp(-2.d0*tau*zita)
        sb=-mu+zita+(mu+zita)*ep
	    sg=1.d0-ep
	    sd=(mu+zita)+(zita-mu)*ep

	   	!!! Remark: tildem does not contains sg. We have rewritten the equations in order to put in evidence sg
        !! this is relevant to avoid explotions when time to maturity is small.
        
         tildem=2.d0*sb/(epsil(nb)**2)

	     tildev=4.d0*zita**2*ep*v0(nb)/(sb)**2

        expression=nupu*(cdlog(sb/(2.d0*zita))+(mu+zita)*tau)+v0(nb)*sg*cdif(nb)/sb

         
        if(nb.eq.2) then
         
         expression=expression+nupu*cdlog((tildem+tau*0.5d0*sg)/tildem)+0.5d0*tau*tildem*tildev/(tildem+0.5d0*tau*sg)
         
         endif
         
	   if(dreal(expression).lt.400.d0) then 
	   
	   Acc(j)=Acc(j)*(1.d0,0.d0)/cdexp(expression)
        
		else
		
		Acc(j)=(0.d0,0.d0)
		
		endif
        enddo
     
        
        
        
        !!! start call price computation   
        qq=-2.d0
        
       phiq=0.5d0*(kapp**2+i*kapp*(2.0d0*qq-1.d0)-(qq**2-qq))
        
     
        
        cdif(1)=phiq*psii
        cdif(2)=phiq*omeg**2-qq+i*kapp
     
        !! imitialization for a product
        App(j)=(1.d0,0.d0)
        
	    do nb=1,1

	       epsq=epsil(nb)**2
           
    	    nupu=2.d0*k1(nb)*theta(nb)/epsq
        
        !! here define mu for both models
        
            if(nb.eq.1) then
        !! here define mu for both models
               mu=-0.5d0*(k1(nb)+epsil(nb)*(rho12(nb)+deltaa)*(i*kapp-qq))
            else
                   mu=-0.5d0*(k1(nb)+epsil(nb)*rho12(nb)*omeg*(i*kapp-qq))
            endif
            
        !! here we define zeta
        zita=0.5d0*(4.d0*mu**2+2.d0*epsq*cdif(nb))**0.5
      
	    ep=cdexp(-2.d0*tau*zita)
        sb=-mu+zita+(mu+zita)*ep
	    sg=1.d0-ep
	    sd=(mu+zita)+(zita-mu)*ep

	   	!!! Remark: tildem does not contains sg. We have rewritten the equations in order to put in evidence sg
        !! this is relevant to avoid explotions when time to maturity is small.
        
         tildem=2.d0*sb/(epsil(nb)**2)

	     tildev=4.d0*zita**2*ep*v0(nb)/(sb)**2

        expression=nupu*(cdlog(sb/(2.d0*zita))+(mu+zita)*tau)+v0(nb)*sg*cdif(nb)/sb

         
        if(nb.eq.2) then
         
         expression=expression+nupu*cdlog((tildem+tau*0.5d0*sg)/tildem)+0.5d0*tau*tildem*tildev/(tildem+0.5d0*tau*sg)
         
         endif
         
	   if(dreal(expression).lt.400.d0) then 
	   
	   App(j)=App(j)*(1.d0,0.d0)/cdexp(expression)
        
		else
		
		App(j)=(0.d0,0.d0)
		
		endif
 enddo
     
     
!! modifica inserita 7/12/07 devi integrare analiticamente il resto
!		App(j)=(App(j)-(1.d0,0.d0))/cdenp
!        Acc(j)=(Acc(j)-(1.d0,0.d0))/cdenc

		App(j)=App(j)/cdenp
        Acc(j)=Acc(j)/cdenc

    enddo
    
    
    
    
      return
	end