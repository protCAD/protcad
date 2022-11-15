C	bestfit.f		Version 1 5/19/1999
c
c	This subroutine performs the best fit of two structures
c	(i.e. moves structure 2 such as to provide the best superposition
c	of the two structure)
c	The new structure 2 is called structure 3
c	It is based on the algorithm of Mclachlan
c
c	Input:
c		-coord1 (size 3xnat1) : vector containing the coordinates
c					of all atoms of molecule 1
c					arranged as x1,y1,z1,x2,y2,z2...
c					Molecule 1 is considered the
c					"reference" structure
c		-nat1		      :	number of atoms in molecule 1
c		-coord2 (size 3xnat2) : coordinates of all atoms of
c					molecule 2
c					Molecule 2 is considered the "test"
c					structure; this is the one
c					that will be moved 
c		-nat2		      : number of atoms in molecule 2
c		-nat		      : number of atoms picked in both
c					structures for superposition
c		list1 (size nat)      :	position of the nat atoms for
c					superposition in molecule 1
c		list2 (size nat)      :	position of the nat atoms for
c					superposition in molecule 2
c	Output:
c		-coord3 (size 3xnat2) : coordinates of all atoms of molecule
c					2 after superposition. coord2
c					remains unchanged
c		-rmsd		      : coordinate RMS between the 2
c					structure (caculated over the
c					nat atoms picked)
c		-ierr		      : flag: 0 if computation of rmsd
c					was OK, 1 otherwise
c
c	Subroutine used:
c		-deter		: compute a 3x3 determinant
c		-dsvdc		: compute the SVD of a matrix
c
c	All computations in double precision
c
	subroutine bestfit(coord1,nat1,coord2,nat2,nat,coord3,
     1		list1,list2,rmsd,ierr,r,xc1,xc2,rmsdat)
c
	integer	i,j,k,nat1,nat2,nat,job,info,ierr,natmax1,natmax2
	integer	list1(nat),list2(nat)
c       
	real*8	deter,tiny
	real*8	rmsd,sign,det,norm
	real*8	a(3,3),u(3,3),v(3,3),r(3,3),d(3),work(3)
	real*8	coord1(3,nat1),coord2(3,nat2)
	real*8	coord3(3,nat2),rmsdat(nat)
	real*8	xc1(3),xc2(3),c(3)
c
	ierr = 0
	tiny = 1.d-14
c
c 1	format ('we are in Fortran')
c	write (6,1)
c
c	Find center of mass of the two sets of atoms which are used
c	for the bestfit
c
	do 100 i = 1,3
		xc1(i) = 0
		xc2(i) = 0 
100	continue
c
	do 300 i = 1,nat
		do 200 j = 1,3
			xc1(j) = xc1(j) + coord1(j,list1(i))
			xc2(j) = xc2(j) + coord2(j,list2(i))
200		continue
300	continue
c
	do 400 i = 1,3
		xc1(i) = xc1(i)/nat
		xc2(i) = xc2(i)/nat
400	continue
c
c	Calculate Covariance matrix :
c
	do 700 i = 1,3
		do 600 j = 1,3
			a(i,j) = 0
			do 500 k = 1,nat
				a(i,j) = a(i,j) + (coord1(i,list1(k))
     1			-xc1(i))*(coord2(j,list2(k))-xc2(j))
500			continue
600		continue
700	continue
c
c	calculate determinant of covariance matrix
c
	det = deter(a)
c
	if(det.eq.0) then
c		write(6,*) 'error in best fit : det = 0 !'
		ierr = 1
c		return
	endif
c
	sign = 1.0
	if(det.le.0.0) sign = -1.0
c
c	Perform SVD on covariance matrix
c
	call dsvdc(a,3,3,3,3,d,u,3,v,3,work,info)
c
	if(info.ne.0) then
c		write(6,*) 'Error in bestfit : SVD failed !'
		ierr = 1
		return
	endif
c
c	Calculate bestfit rotational matrix r
c
	if(d(2).gt.tiny) then
		if(d(3).le.tiny) then
			sign = 1.0
			u(1,3) = u(2,1)*u(3,2) - u(3,1)*u(2,2)
			u(2,3) = u(3,1)*u(1,2) - u(1,1)*u(3,2)
			u(3,3) = u(1,1)*u(2,2) - u(2,1)*u(1,2)
			v(1,3) = v(2,1)*v(3,2) - v(3,1)*v(2,2)
			v(2,3) = v(3,1)*v(1,2) - v(1,1)*v(3,2)
			v(3,3) = v(1,1)*v(2,2) - v(2,1)*v(1,2)
		endif
		do 900 i = 1,3
			do 800 j = 1,3
				r(i,j) = u(i,1)*v(j,1)+u(i,2)*v(j,2)
     1				+ sign*u(i,3)*v(j,3)
800			continue
900		continue
	else
		work(1) = u(2,1)*v(3,1)-u(3,1)*v(2,1)
		work(2) = u(3,1)*v(1,1)-u(1,1)*v(3,1)
		work(3) = u(1,1)*v(2,1)-u(2,1)*v(1,1)
		norm = 0
		do 1000 i = 1,3
			norm = norm + work(i)*work(i)
1000		continue
		if(norm.ne.0) then
			do 1100 i = 1,3
				work(i) = u(i,1) + v(i,1)
1100			continue
		else
			work(1) = -u(2,1)
			work(2) = u(1,1)
			work(3) = 0
		endif
		norm = 0
		do 1200 i = 1,3
			norm = norm + work(i)*work(i)
1200		continue
		norm = sqrt(norm)
		if(norm.eq.0) then
			ierr = 1
			return
		endif
		do 1300 i = 1,3
			work(i) = work(i)/norm
1300		continue
		do 1500 i = 1,3
			do 1400 j = 1,3
				r(i,j) = 2*work(j)*work(i)
1400			continue
			r(i,i) = r(i,i) - 1
1500		continue
	endif

	det = deter(r)
c
c	Apply rotation on molecule 2 and store results in coord3 :
c
	do 1900 i = 1,nat2
		do 1700 j = 1,3
			c(j) = 0
			do 1600 k = 1,3
				c(j) = c(j) + r(j,k)*(coord2(k,i)-
     1				xc2(k))
1600			continue
1700		continue
		do 1800 j = 1,3
			coord3(j,i) = c(j) + xc1(j)
1800		continue
1900	continue
c
c	Calculate rmsd
c
	rmsd = 0
	do 2100 i = 1,nat
        rmsdat(i) = 0
		do 2000 j = 1,3
			rmsdat(i) = rmsdat(i) + (coord3(j,list2(i)) -
     1			coord1(j,list1(i)))**2
2000	continue
        rmsd = rmsd + rmsdat(i)
2100    continue
c
	rmsd = sqrt(rmsd/nat)
c
	return
	end
c
c	Deter.for	Version 1 23/11/1987		Patrice Koehl
c
c	This function calculates a 3 x 3 determinant by developping along
c	the first column.
c	det(i,j), i = 1,3 ; j = 1,3 are the coefficients of the determinant
c
	function deter(det)
c
	real*8	a1,a2,a3,deter
	real*8	det(3,3)
c
	a1	= det(2,2)*det(3,3) - det(3,2)*det(2,3)
	a2	= det(3,2)*det(1,3) - det(1,2)*det(3,3)
	a3	= det(1,2)*det(2,3) - det(2,2)*det(1,3)
c
	deter   = a1*det(1,1) + a2*det(2,1) + a3*det(3,1)
c
	return
	end
c
C////////////////////////////////////////////////////////////////////   
C                                                                       
C THIS SUBROUTINE PERFORMS THE SINGULAR VALUE DECOMPOSITION             
C                                                                       
C#NUMPAC#DSVDC                REVISED ON 1984-11-30                      
C     AMAX1--> DMAX1, ABS--> DABS, SIGN--> DSIGN, SQRT--> DSQRT         
C////////////////////////////////////////////////////////////////////   
C                                                                       
      SUBROUTINE DSVDC(A,KA,M,N,ISW,Q,U,KU,V,KV,W,IND)                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(KA,N),U(KU,N),V(KV,N),Q(N),W(N)                       
      IND=30000                                                         
      MN=MIN0(M,N)                                                      
      IF(MN.LT.1.OR.M.GT.KA.OR.M.GT.KU) GO TO 490                       
      MU=ISW/2                                                          
      MV=MOD(ISW,2)                                                     
      IF(MU.LT.0.OR.MU.GT.1.OR.MV.LT.0.OR.MV.GT.1) GO TO 490            
      IF(MV.EQ.1.AND.N.GT.KV) GO TO 490                                 
      M1N=MIN0(M+1,N)                                                   
      DO 10 J=1,N                                                       
      DO 10 I=1,M                                                       
   10 U(I,J)=A(I,J)                                                     
      ANORM=0.                                                          
      G=0.                                                              
      DO 100 I=1,M1N                                                    
      Q(I)=0.                                                           
      W(I)=G                                                            
      IF(I.GT.M) GO TO 100                                              
      IP1=I+1                                                           
      G=U(I,I)                                                          
      IF(I.EQ.M) GO TO 30                                               
      SUM=0.                                                            
      DO 20 K=I,M                                                       
   20 SUM=U(K,I)*U(K,I)+SUM                                             
      S=SUM                                                             
      G=-DSIGN(DSQRT(S),G)                                              
      H=U(I,I)*G-S                                                      
      U(I,I)=U(I,I)-G                                                   
   30 Q(I)=G                                                            
      IF(I.EQ.N) GO TO 100                                              
      IF(S.EQ.0..OR.I.EQ.M) GO TO 60                                    
      DO 50 J=IP1,N                                                     
      SUM=0.                                                            
      DO 40 K=I,M                                                       
   40 SUM=U(K,I)*U(K,J)+SUM                                             
      F=SUM/H                                                           
      DO 50 K=I,M                                                       
   50 U(K,J)=U(K,I)*F+U(K,J)                                            
   60 G=U(I,IP1)                                                        
      IF(IP1.EQ.N) GO TO 100                                            
      SUM=0.                                                            
      DO 70 K=IP1,N                                                     
   70 SUM=U(I,K)*U(I,K)+SUM                                             
      S=SUM                                                             
      G=-DSIGN(DSQRT(S),G)                                              
      H=U(I,IP1)*G-S                                                    
      U(I,IP1)=U(I,IP1)-G                                               
      IF(S.EQ.0..OR.I.EQ.M) GO TO 100                                   
      DO 90 J=IP1,M                                                     
      SUM=0.                                                            
      DO 80 K=IP1,N                                                     
   80 SUM=U(I,K)*U(J,K)+SUM                                             
      F=SUM/H                                                           
      DO 90 K=IP1,N                                                     
   90 U(J,K)=U(I,K)*F+U(J,K)                                            
  100 ANORM=DMAX1(DABS(Q(I))+DABS(W(I)),ANORM)                          
      TOL=DMACH(TOL)*ANORM                                              
      IF(MV.EQ.0) GO TO 180                                             
      DO 170 II=1,M1N                                                   
      I=M1N+1-II                                                        
      IF(I.EQ.N) GO TO 170                                              
      IP1=I+1                                                           
      IF(I.EQ.M1N) GO TO 150                                            
      IF(IP1.EQ.N.OR.W(IP1).EQ.0.) GO TO 130                            
      H=U(I,IP1)*W(IP1)                                                 
      DO 120 J=IP1,M1N                                                  
      SUM=0.                                                            
      DO 110 K=IP1,N                                                    
  110 SUM=U(I,K)*V(K,J)+SUM                                             
      F=SUM/H                                                           
      DO 120 K=IP1,N                                                    
  120 V(K,J)=U(I,K)*F+V(K,J)                                            
  130 DO 140 J=IP1,M1N                                                  
  140 V(I,J)=0.                                                         
  150 DO 160 J=IP1,N                                                    
  160 V(J,I)=0.                                                         
  170 V(I,I)=1.                                                         
  180 IF(MU.EQ.0) GO TO 260                                             
      DO 250 II=1,MN                                                    
      I=MN+1-II                                                         
      IF(I.EQ.MN) GO TO 200                                             
      IP1=I+1                                                           
      DO 190 J=IP1,MN                                                   
  190 U(I,J)=0.                                                         
  200 IF(Q(I).EQ.0.) GO TO 250                                          
      IF(I.EQ.MN) GO TO 230                                             
      H=U(I,I)*Q(I)                                                     
      DO 220 J=IP1,MN                                                   
      SUM=0.                                                            
      DO 210 K=IP1,M                                                    
  210 SUM=U(K,I)*U(K,J)+SUM                                             
      F=SUM/H                                                           
      DO 220 K=I,M                                                      
  220 U(K,J)=U(K,I)*F+U(K,J)                                            
  230 DO 240 K=I,M                                                      
  240 U(K,I)=U(K,I)/Q(I)                                                
  250 IF(I.LT.M.OR.Q(I).EQ.0.) U(I,I)=U(I,I)+1.                         
  260 IF(ANORM.EQ.0.) GO TO 470                                         
      DO 390 KK=1,M1N                                                   
      K=M1N+1-KK                                                        
      DO 360 IT=1,30                                                    
      DO 270 LL=1,K                                                     
      L=K+1-LL                                                          
      IF(DABS(W(L)).LT.TOL) GO TO 310                                   
      IF(DABS(Q(L)).LT.TOL) GO TO 280                                   
  270 CONTINUE                                                          
  280 C=0.                                                              
      S=-1.                                                             
      DO 300 II=2,L                                                     
      I=L+1-II                                                          
      F=-W(I+1)*S                                                       
      W(I+1)=W(I+1)*C                                                   
      IF(DABS(F).LT.TOL) GO TO 310                                      
      G=Q(I)                                                            
      Q(I)=DSQRT(G*G+F*F)                                               
      C=G/Q(I)                                                          
      S=F/Q(I)                                                          
      IF(MV.EQ.0) GO TO 300                                             
      DO 290 J=1,N                                                      
      X=V(J,I)                                                          
      V(J,I)=V(J,L)*S+X*C                                               
  290 V(J,L)=V(J,L)*C-X*S                                               
  300 CONTINUE                                                          
  310 IF(L.EQ.K) GO TO 370                                              
      G=W(K-1)                                                          
      H=W(K)                                                            
      X=Q(L)                                                            
      Y=Q(K-1)                                                          
      Z=Q(K)                                                            
      F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(H*Y*2.)                              
      F=((X-Z)*(X+Z)+H*(Y/(DSIGN(DSQRT(F*F+1.),F)+F)-H))/X              
      C=1.                                                              
      S=1.                                                              
      LP1=L+1                                                           
      DO 350 I=LP1,K                                                    
      H=W(I)*S                                                          
      G=W(I)*C                                                          
      W(I-1)=DSQRT(F*F+H*H)                                             
      C=F/W(I-1)                                                        
      S=H/W(I-1)                                                        
      F=X*C+G*S                                                         
      G=G*C-X*S                                                         
      H=Q(I)*S                                                          
      Y=Q(I)*C                                                          
      IF(MV.EQ.0) GO TO 330                                             
      DO 320 J=1,N                                                      
      X=V(J,I-1)                                                        
      V(J,I-1)=V(J,I)*S+X*C                                             
  320 V(J,I)=V(J,I)*C-X*S                                               
  330 Q(I-1)=DSQRT(F*F+H*H)                                             
      C=F/Q(I-1)                                                        
      S=H/Q(I-1)                                                        
      F=G*C+Y*S                                                         
      X=Y*C-G*S                                                         
      IF(MU.EQ.0) GO TO 350                                             
      DO 340 J=1,M                                                      
      Y=U(J,I-1)                                                        
      U(J,I-1)=U(J,I)*S+Y*C                                             
  340 U(J,I)=U(J,I)*C-Y*S                                               
  350 CONTINUE                                                          
      W(L)=0.                                                           
      W(K)=F                                                            
  360 Q(K)=X                                                            
      GO TO 480                                                         
  370 IF(Q(K).GE.0.) GO TO 390                                          
      Q(K)=-Q(K)                                                        
      IF(MV.EQ.0) GO TO 390                                             
      DO 380 J=1,N                                                      
  380 V(J,K)=-V(J,K)                                                    
  390 CONTINUE                                                          
      IF(N.EQ.1) GO TO 470                                              
      K=MN                                                              
  400 L=1                                                               
      II=1                                                              
      LL=1                                                              
      DO 420 I=2,K                                                      
      IF(Q(I).GT.Q(L)) GO TO 410                                        
      L=I                                                               
      GO TO 420                                                         
  410 II=I                                                              
      LL=L                                                              
  420 CONTINUE                                                          
      IF(II.EQ.LL) GO TO 460                                            
      S=Q(II)                                                           
      Q(II)=Q(LL)                                                       
      Q(LL)=S                                                           
      IF(MV.EQ.0) GO TO 440                                             
      DO 430 J=1,N                                                      
      S=V(J,II)                                                         
      V(J,II)=V(J,LL)                                                   
  430 V(J,LL)=S                                                         
  440 IF(MU.EQ.0) GO TO 460                                             
      DO 450 J=1,M                                                      
      S=U(J,II)                                                         
      U(J,II)=U(J,LL)                                                   
  450 U(J,LL)=S                                                         
  460 K=II-1                                                            
      IF(K.GE.2) GO TO 400                                              
  470 IND=0                                                             
      RETURN                                                            
  480 IND=20000                                                         
  490 RETURN                                                            
      END                                                               
C                                                                       
C////////////////////////////////////////////////////////////////////   
C                                                                       
C#NUMPAC#DMACH               REVISED ON 1984-11-30                      
C			     PK:Revised on 11/5/1997: "save" eps (required
C				when using f2c which does not assume variable as being static)
C                                                                       
C////////////////////////////////////////////////////////////////////   
C                                                                       
      FUNCTION DMACH(X)                                                 
      REAL*8 DMACH,X,ONE,EPS                                            
	SAVE EPS
      DATA ONE,IFIRST/1.0,1/                                            
      IF(IFIRST.EQ.0) GO TO 20                                          
      IFIRST=0                                                          
      EPS=ONE                                                           
   10 EPS=EPS*0.5D0                                                     
      IF(EPS+ONE.NE.ONE) GO TO 10                                       
      EPS=EPS+EPS                                                       
   20 DMACH=EPS                                                         
      RETURN                                                            
      END
c	Vector		Version 1 8/2/1990		Patrice Koehl
c
c	This file contains several subroutine that can be used for any
c	vector operations (vector in 3D cartesian space)
c
c	This includes :	crossvect	: cross vector of two vectors
c			dotvect		: dot product of two vectors
c			normvect	: norm of a vector
c			detvect		: determinant of three vectors
c			diffvect	: substract two vectors
c			addvect		: add two vectors
c
c	For each subroutine : u is for vectors (arrays of size 3)
c			      all other value are scalar
c			      calculations are done in double precision
c
c	1 . crossvect :
c
	subroutine crossvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	u3(1) = u1(2)*u2(3) - u1(3)*u2(2)
	u3(2) = -u1(1)*u2(3) + u1(3)*u2(1)
	u3(3) = u1(1)*u2(2) - u1(2)*u2(1)
c
	return
	end
c
c	2. dotvect :
c
	subroutine dotvect(u1,u2,dot)
c
	integer	i
c
	real*8	u1(3),u2(3),dot
c
	dot = 0.d0
	do 100 i = 1,3
		dot = dot + u1(i)*u2(i)
100	continue
c
	return
	end
c
c	3. normvect :
c
	subroutine normvect(u1,norm)
c
	real*8	u1(3),norm
c
	call dotvect(u1,u1,norm)
	norm = dsqrt(norm)
c
	return
	end
c
c	4. detvect :
c
	subroutine detvect(u1,u2,u3,det)
c
	real*8	u1(3),u2(3),u3(3),det,u4(3)
c
	call crossvect(u2,u3,u4)
	call dotvect(u1,u4,det)
c
	return
	end
c
c	5. diffvect :
c
	subroutine diffvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u2(i) - u1(i)
100	continue
c
	return
	end
c
c	6. addvect :
c
	subroutine addvect(u1,u2,u3)
c
	real*8	u1(3),u2(3),u3(3)
c
	integer i
c
	do 100 i = 1,3
		u3(i) = u1(i) + u2(i)
100	continue
c
	return
	end
c
c	7. Normalise a vector : given a vector u1, output u1/norm(u1) :
c
	subroutine unitvector(u1,u2)
c
	real*8  u1(3),u2(3),norm
c
	integer i
c
	call normvect(u1,norm)
c
	do 100 i = 1,3
		u2(i) = u1(i)/norm
100	continue
c
	return
	end
