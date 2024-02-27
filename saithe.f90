!|--------------------------------------------------------------------------------------|
!| This program is based on Lattice Boltzmann Method v2.0                               |
!| with Immersed-Boundary of "bounce back" Boundary condition                           |
!| written by vandine, 2010/9/23                                                        |
!|--------------------------------------------------------------------------------------|
program main
	implicit none
	!============= Parameters setting Area ===============================!	
	! Reynolds number,U_infinity,Reduced frequency,Plunge Amplitude
	real*8,parameter::Re=1000d0,Mach=0.1d0,u_inf=Mach/sqrt(3d0),lambda=1.0d0,freq=2.0d0
	! Airfoil Thickness,Chord,and outer spacing 
	real*8,parameter::thick=0.12d0,c=200d0,epsilon=0.2d0
	real*8,parameter::t2=thick*c+2d0*epsilon,Dc2=1d0/(c+2d0*epsilon)
	! Wind tunnel Domain size and frame
	real*8,parameter::H_c=3d0,Ly=H_c*c,Lx=2d0*Ly
	integer,parameter::mx=int(Lx),my=int(Ly),q=9,frame=10,frame2=50,mchord=int(c),TOTAL=int(Lx/u_inf)
	integer,parameter::mx0=my/2-int(c)/4,mx2=my/2+3*int(c)/4
	integer,parameter::my0=my/2-int(c)/4,my2=my/2+int(c)/4
	!=====================================================================!
	real*8,parameter::pi=4d0*atan(1d0),D6=1d0/6d0,D4=0.25d0
	real*8,dimension(1:q)::w,f1,uex,uey,windvec
	integer,dimension(1:q)::ex,ey,phi2,lbb
	real*8,dimension(-1:Mx+1)::x
	real*8,dimension(-1:My+1)::y
	integer,dimension(-1:Mx+1,-1:My+1)::phi,wall,downwind
	real*8,dimension(-1:Mx+1,-1:My+1)::rho,u,v
	real*8,dimension(1:q,-1:Mx+1,-1:My+1)::f,f0
	real*8 tau,omega,nu,q_inf,usqr,vec,err,mu,A,frequency,feq
	integer i,j,k,l,t,ld
	real*8 Px,Py,x2,y2,x_c,sqr,r,distance,omega2,a_angle,temp
	real*8 Fx,Fy,FxP,FyP,FxV,FyV
	real*8 Fx2,Fy2,Fx2P,Fx2V,Fy2P,Fy2V
	real*8 tempxP,tempyP,tempxV,tempyV
	real*8 moment,moment2,Drhox,Drhoy
	real*8,dimension(0:int(c))::amp,ym,vm,yb
	real tbegin,tend
	open(UNIT= 9,FILE='Cd.dat',FORM='FORMATTED',STATUS='REPLACE',ACCESS='APPEND')
	open(UNIT=10,FILE='Cl.dat',FORM='FORMATTED',STATUS='REPLACE',ACCESS='APPEND')
	open(UNIT=11,FILE='Mo.dat',FORM='FORMATTED',STATUS='REPLACE',ACCESS='APPEND')
	close(9); close(10); close(11)
	! D2Q9
	w(1)=4d0/9d0; w(3)=1d0/9d0;  w(5)=1d0/9d0;  w(7)=1d0/9d0;  w(9)=1d0/9d0
	              w(2)=1d0/36d0; w(4)=1d0/36d0; w(6)=1d0/36d0; w(8)=1d0/36d0
	ex(1)= 0; ex(2)=-1; ex(3)=-1; ex(4)=-1; ex(5)= 0; ex(6)= 1; ex(7)= 1; ex(8)= 1; ex(9)= 0
	ey(1)= 0; ey(2)= 1; ey(3)= 0; ey(4)=-1; ey(5)=-1; ey(6)=-1; ey(7)= 0; ey(8)= 1; ey(9)= 1
	uex(1)= 0d0; uex(3)=-1d0; uex(5)= 0d0; uex(7)= 1d0; uex(9)= 0d0
	             uex(2)=-1d0/sqrt(2d0); uex(4)=-1d0/sqrt(2d0); uex(6)= 1d0/sqrt(2d0); uex(8)= 1d0/sqrt(2d0)
	uey(1)= 0d0; uey(3)= 0d0; uey(5)=-1d0; uey(7)= 0d0; uey(9)= 1d0
	             uey(2)= 1d0/sqrt(2d0); uey(4)=-1d0/sqrt(2d0); uey(6)=-1d0/sqrt(2d0); uey(8)= 1d0/sqrt(2d0)	
	lbb(1)=1; lbb(2)=6; lbb(3)=7; lbb(4)=8; lbb(5)=9; lbb(6)=2; lbb(7)=3; lbb(8)=4; lbb(9)=5	             
	! fluid properties	
	nu=u_inf*c/Re
	tau=0.5d0+3d0*nu
	omega=1d0/tau
	q_inf=2d0/(u_inf*u_inf*c)
	frequency=freq*2d0*pi*u_inf/c
	! grid	
	do i=-1,Mx+1; x(i)=dble(i); enddo
	do j=-1,My+1; y(j)=dble(j); enddo
	! initial data ( steady state at zero AOA )	
	t=0
	do j=0,My; do i=0,Mx
	   rho(i,j)=1d0
	   u(i,j)=u_inf; v(i,j)=0d0
	enddo; enddo
	Py=0.5d0*Ly
	Px=0.5d0*Ly
	do k=0,mchord
	   x2=dble(k)/c
	   amp(k)=(0.1625d0*x2*x2-0.0825d0*x2+0.02d0)*c
	   x2=(dble(k)+epsilon)*Dc2
	   sqr=x2*x2
	   yb(k)=5d0*t2*(0.2969d0*sqrt(x2)-0.1260d0*x2-0.3516d0*sqr+0.2843d0*sqr*x2-0.1015d0*sqr*sqr)
	enddo	
	! initial density distribution
	do j=0,My; do i=0,Mx
	   usqr=u(i,j)*u(i,j)+v(i,j)*v(i,j)
	   do l=1,q
	      vec=u(i,j)*ex(l)+v(i,j)*ey(l)
	      f(l,i,j)= w(l)*rho(i,j)*( 1d0+3d0*vec+4.5d0*vec*vec-1.5d0*usqr )
	   enddo
	enddo; enddo
	tbegin=secnds(0.0)		
100	t=t+1
	f0=f
!=========================================================================================================!
!================== Main Algorithm starts from here ======================================================!
! 1.Immersed Boundary phase ( NACA-0012 airfoil )
!   ******Set prefered Immersed Boundary phase and motion in the following section**********   !
	do k=0,mchord
	   x2=dble(k)/c
	   ym(k)=amp(k)*sin(2d0*pi*x2/lambda-frequency*dble(t))
	   vm(k)=-frequency*amp(k)*cos(2d0*pi*x2/lambda-frequency*dble(t))
	enddo	
	do j=my0,my2; do i=mx0,mx2
	   y2=y(j)-Py-ym(i-mx0)
	   if (  y2<=yb(i-mx0) .AND. y2>=-yb(i-mx0) ) then
	      phi(i,j)=1
	      u(i,j)=0d0
	      v(i,j)=vm(i-mx0)
	   else
	      phi(i,j)=0
	   endif
	enddo; enddo
!   ****************************************************************************************   !
! 2.wall determination
	do j=my0,my2; do i=mx0,mx2
	   if (phi(i,j)==1) then
	      do l=1,q
	         phi2(l)=phi(i+ex(l),j+ey(l))
	      enddo
	      wall(i,j)=1-minval(phi2)
	   else
	      wall(i,j)=0
	   endif
	enddo; enddo
! 3.flow field collision and streaming
	do j=0,My; do i=0,Mx
	   if (phi(i,j)==0) then
	      usqr=u(i,j)*u(i,j)+v(i,j)*v(i,j)
	      do l=1,q
	         vec=u(i,j)*ex(l)+v(i,j)*ey(l)
	         feq = w(l)*rho(i,j)*( 1d0+3d0*vec+4.5d0*vec*vec-1.5d0*usqr )
	         f(l,i+ex(l),j+ey(l))=f0(l,i,j)-omega*(f0(l,i,j)-feq)	                         
	      enddo
	   endif
	enddo; enddo
! 4.IB wall collision,streaming, bounce back,momentum exchange,and downwind-direction determination  
	Fx=0d0; Fy=0d0; moment=0d0 
	do j=my0,my2; do i=mx0,mx2	    
	   if (wall(i,j)==1) then
	      f1(:)=f(:,i,j)
	      usqr=u(i,j)*u(i,j)+v(i,j)*v(i,j)   
	      do l=1,q
	         vec=u(i,j)*ex(l)+v(i,j)*ey(l)
	         feq = w(l)*rho(i,j)*( 1d0+3d0*vec+4.5d0*vec*vec-1.5d0*usqr )
	         if (phi(i+ex(l),j+ey(l))==0) then ! bounce back and momentum exchange
	            f(l,i+ex(l),j+ey(l))=f0(l,i,j)-omega*(f0(l,i,j)-feq)        
	            f(l,i,j)=f1(lbb(l))+6d0*w(l)*rho(i,j)*vec
	            temp=f1(lbb(l))+f(l,i+ex(l),j+ey(l))
	            Drhox=temp*ex(lbb(l))
	            Drhoy=temp*ey(lbb(l))
	            temp=Fx+Drhox; Fx=temp
	            temp=Fy+Drhoy; Fy=temp
	            temp=moment+Drhox*(y(j)-Py)-Drhoy*(x(i)-Px); moment=temp
	         elseif (phi(i-ex(l),j-ey(l))==1) then ! not streaming out (colision only)
	            f(l,i,j)=f0(l,i,j)-omega*(f0(l,i,j)-feq)
	            temp=f0(l,i,j)-f(l,i,j)
	            Drhox=temp*ex(l)
	            Drhoy=temp*ey(l)
	            temp=Fx+Drhox; Fx=temp
	            temp=Fy+Drhoy; Fy=temp
	            temp=moment+Drhox*(y(j)-Py)-Drhoy*(x(i)-Px); moment=temp
	         endif
	         windvec(l)=uex(l)*u(i,j)+uey(l)*v(i,j) ! compute u dot unit vector of e
	      enddo
	      downwind(i,j)=maxloc(windvec,DIM=1) ! downwind-side direction
	   endif
	enddo; enddo
! 5.Macroscopic variables evaluation(rho,u,v) 
	do j=0,My; do i=0,Mx
	   rho(i,j)=f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)+f(9,i,j)
	   if (phi(i,j)==0) then
	      temp=1d0/rho(i,j)
	      u(i,j)=(f(6,i,j)+f(7,i,j)+f(8,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j))*temp
	      v(i,j)=(f(2,i,j)+f(9,i,j)+f(8,i,j)-f(6,i,j)-f(5,i,j)-f(4,i,j))*temp
	   endif
	enddo; enddo
! 6.Boundary condition	
        do i=1,Mx-1
           ! bottom boundary / free stream
	   j=0    
	   rho(i,j)=rho(i,j+1)
	   u(i,j)=u(i,j+1); v(i,j)=v(i,j+1)
	   vec=-u(i,j)+v(i,j)
	   f(2,i,j)=f(6,i,j)+6d0*w(6)*rho(i,j)*vec
	   vec=       +v(i,j)
	   f(9,i,j)=f(5,i,j)+6d0*w(5)*rho(i,j)*vec
	   vec= u(i,j)+v(i,j)
	   f(8,i,j)=f(4,i,j)+6d0*w(4)*rho(i,j)*vec
	   ! top    boundary / free stream
	   j=my   
	   rho(i,j)=rho(i,j-1)
	   u(i,j)=u(i,j-1); v(i,j)=v(i,j-1)
	   vec=-u(i,j)-v(i,j)
	   f(4,i,j)=f(8,i,j)+6d0*w(8)*rho(i,j)*vec
	   vec=       -v(i,j)
	   f(5,i,j)=f(9,i,j)+6d0*w(9)*rho(i,j)*vec
	   vec= u(i,j)-v(i,j)
	   f(6,i,j)=f(2,i,j)+6d0*w(2)*rho(i,j)*vec 
	enddo
	do j=1,My-1
	   ! left   boundary / forced inlet
	   i=0    
	   rho(i,j)=1d0
	   u(i,j)=u_inf; v(i,j)=0d0
	   vec= u(i,j)+v(i,j)	   
	   f(8,i,j)=f(4,i,j)+6d0*w(4)*rho(i,j)*vec
	   vec= u(i,j)	    
	   f(7,i,j)=f(3,i,j)+6d0*w(3)*rho(i,j)*vec
	   vec= u(i,j)-v(i,j)
	   f(6,i,j)=f(2,i,j)+6d0*w(2)*rho(i,j)*vec
	   ! right  boundary / free-stream outlet
	   i=mx   
	   rho(i,j)=1d0
	   u(i,j)=u(i-1,j); v(i,j)=v(i-1,j)
	   vec=-u(i,j)-v(i,j)	   
	   f(4,i,j)=f(8,i,j)+6d0*w(8)*rho(i,j)*vec
	   vec=-u(i,j)
	   f(3,i,j)=f(7,i,j)+6d0*w(7)*rho(i,j)*vec
	   vec=-u(i,j)+v(i,j)
	   f(2,i,j)=f(6,i,j)+6d0*w(6)*rho(i,j)*vec 
	enddo
	! left  bottom corner
	i=0; j=0        
	rho(i,j)=1d0
	u(i,j)=u_inf; v(i,j)=0d0
	vec= u(i,j)	    
	f(7,i,j)=f(3,i,j)+6d0*w(3)*rho(i,j)*vec
	vec= u(i,j)+v(i,j)	   
	f(8,i,j)=f(4,i,j)+6d0*w(4)*rho(i,j)*vec	
	vec=       +v(i,j)
	f(9,i,j)=f(5,i,j)+6d0*w(5)*rho(i,j)*vec
	! left  top    corner
	i=0; j=my         
	rho(i,j)=1d0
	u(i,j)=u_inf; v(i,j)=0d0
	vec=       -v(i,j)
	f(5,i,j)=f(9,i,j)+6d0*w(9)*rho(i,j)*vec
	vec= u(i,j)-v(i,j)
	f(6,i,j)=f(2,i,j)+6d0*w(2)*rho(i,j)*vec
	vec= u(i,j)	    
	f(7,i,j)=f(3,i,j)+6d0*w(3)*rho(i,j)*vec
	! right bottom corner 
	i=mx; j=0       
	rho(i,j)=1d0
	u(i,j)=u(i-1,j+1); v(i,j)=v(i-1,j+1)
	vec=       +v(i,j)
	f(9,i,j)=f(5,i,j)+6d0*w(5)*rho(i,j)*vec	
	vec=-u(i,j)
	f(3,i,j)=f(7,i,j)+6d0*w(7)*rho(i,j)*vec
	vec=-u(i,j)+v(i,j)
	f(2,i,j)=f(6,i,j)+6d0*w(6)*rho(i,j)*vec
	! right top    corner 
	i=mx; j=my        
	rho(i,j)=1d0
	u(i,j)=u(i-1,j-1); v(i,j)=v(i-1,j-1)
	vec=-u(i,j)
	f(3,i,j)=f(7,i,j)+6d0*w(7)*rho(i,j)*vec
	vec=-u(i,j)-v(i,j)
	f(4,i,j)=f(8,i,j)+6d0*w(8)*rho(i,j)*vec
	vec=       -v(i,j)
	f(5,i,j)=f(9,i,j)+6d0*w(9)*rho(i,j)*vec	
! 7.Correct density inside IB and reset Distribution	 
	do j=my0,my2; do i=mx0,mx2
	   if (wall(i,j)==1 ) then
	      ld=downwind(i,j)
	      if (phi(i+ex(ld),j+ey(ld))==1 .AND. wall(i+ex(ld),j+ey(ld))==0) then 
	         rho(i+ex(ld),j+ey(ld))=rho(i,j)
	         usqr=u(i,j)*u(i,j)+v(i,j)*v(i,j)	     
	         do l=1,q
	            vec=u(i,j)*ex(l)+v(i,j)*ey(l)
	            f(l,i+ex(ld),j+ey(ld))=w(l)*rho(i,j)*( 1d0+3d0*vec+4.5d0*vec*vec-1.5d0*usqr )
	         enddo
	      endif
	   endif
	enddo; enddo	
!================== Main Algorithm ends here =============================================================!
!=========================================================================================================!	
	! terminal output at each frame stage  	    	
	if (t/frame*frame==t) then
	   write(*,*)t,'Cd=',Fx*q_inf,'Cl=',Fy*q_inf   	   
	   ! Hydrodynamics force ( volume integral method)
	   Fx2=0d0; Fy2=0d0; moment2=0d0
	   do j=my0,my2; do i=mx0,mx2	    
	      if (phi(i,j)==1) then
	         mu=rho(i,j)*nu
	         tempxP=-(rho(i+1,j)-rho(i-1,j))*D6+2d0*mu*(u(i+1,j)-2d0*u(i,j)+u(i-1,j)) 
	         tempyP=-(rho(i,j+1)-rho(i,j-1))*D6+2d0*mu*(v(i,j+1)-2d0*v(i,j)+v(i,j-1)) 
	         tempxV=mu*( u(i,j+1)+u(i,j-1)-2d0*u(i,j)+D4*(v(i+1,j+1)+v(i-1,j-1)-v(i+1,j-1)-v(i-1,j+1)) ) 
	         tempyV=mu*( v(i+1,j)+v(i-1,j)-2d0*v(i,j)+D4*(u(i+1,j+1)+u(i-1,j-1)-u(i+1,j-1)-u(i-1,j+1)) ) 
	         temp=Fx2+tempxP+tempxV; Fx2=temp
	         temp=Fy2+tempyP+tempyV; Fy2=temp
	         temp=moment2+(tempxP+tempxV)*(y(j)-Py)-(tempyP+tempyV)*(x(i)-Px); moment2=temp
	      endif
	   enddo; enddo
	   open(UNIT= 9,FILE='Cd.dat',FORM='FORMATTED',STATUS='OLD',ACCESS='APPEND')
	   open(UNIT=10,FILE='Cl.dat',FORM='FORMATTED',STATUS='OLD',ACCESS='APPEND')
	   open(UNIT=11,FILE='Mo.dat',FORM='FORMATTED',STATUS='OLD',ACCESS='APPEND')
	   write( 9,*)t,Fx*q_inf,Fx2*q_inf,ym(mchord) 
	   write(10,*)t,Fy*q_inf,Fy2*q_inf
	   write(11,*)t,moment*q_inf/c,moment2*q_inf/c
	   close(9); close(10); close(11)
!	   if (t/frame2*frame2==t)call output(mx,my,x,y,phi,rho,u,v,t/frame2,u_inf)  
	endif 
	! time iteration
	if (t<total) goto 100
	tend=secnds(tbegin)
	write(*,*)'duration=',tend
	call output(mx,my,x,y,phi,rho,u,v,0,u_inf) 


end program	
!|-------------------------------------------------------------------------------|
!| The main code is end                                                          |
!|-------------------------------------------------------------------------------|

	


			
               	
! =========================================================
! subroutine : output file for post processor (Paraview)
subroutine output(mx,my,x,y,phi,rho,u,v,count,u_inf)
        implicit none
	integer mx,my,count
	real*8 u_inf,q_inf,Du_inf
	real*8,dimension(-1:Mx+1)::x
	real*8,dimension(-1:My+1)::y
	integer,dimension(-1:Mx+1,-1:My+1)::phi                           
	real*8,dimension(-1:Mx+1,-1:My+1)::rho,u,v
        integer i,j
        character(len=12) filename
        q_inf=2d0/(3d0*u_inf*u_inf)
        Du_inf=1d0/u_inf	
	write(filename,'(A5,I3.3,A4)') 'post_',count,'.vtk'
	open(UNIT=8,FILE=filename,FORM='FORMATTED',STATUS='REPLACE')	
	write(8,'(A)')'# vtk DataFile Version 3.0'
        write(8,'(A10,I4,A1,I5)')'pitching acceleration'
        write(8,'(A)')'ASCII'
        write(8,'(A)')'DATASET STRUCTURED_POINTS'
        write(8,'(A,I5,I5,I5)')'DIMENSIONS',(Mx+1),(My+1),1
        write(8,'(A,F12.8,F12.8,F12.8)')'ORIGIN',x(0),y(0),0d0  
        write(8,'(A,I3,I3,I3)')'SPACING',1,1,1
	write(8,'(A,I8)')'POINT_DATA', (Mx+1)*(My+1)
        write(8,'(A)') 'SCALARS phi float 1'
	write(8,'(A)') 'LOOKUP_TABLE default'
        do j=0,My;do i=0,Mx
           write(8,'(I1)') phi(i,j)
        enddo; enddo 
        write(8,'(A)') 'SCALARS Cp float 1'
	write(8,'(A)') 'LOOKUP_TABLE default'
        do j=0,My;do i=0,Mx
           if (phi(i,j)==1) then
              write(8,'(I1)') 0
           else
              write(8,'(E11.4E2)') (rho(i,j)-1d0)*q_inf
           endif
        enddo; enddo 
        write(8,'(A)') 'SCALARS vorticity float 1'
	write(8,'(A)') 'LOOKUP_TABLE default'
        do j=0,My;do i=0,Mx
           if (i==0 .OR. i==mx .OR. j==0 .OR. j==my .OR. phi(i,j)==1) then
              write(8,'(I1)') 0
           else
              write(8,'(E11.4E2)') 0.5d0*(v(i+1,j)-v(i-1,j)-u(i,j+1)+u(i,j-1))*Du_inf
           endif
        enddo; enddo
	write(8,'(A)') 'VECTORS velocity float' 
        do j=0,My;do i=0,Mx
           write(8,'(E11.4E2,E12.4E2,I2)') u(i,j)*Du_inf,v(i,j)*Du_inf,0
        enddo; enddo   
	close(8)
end subroutine  
