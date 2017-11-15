Program main
	
	Implicit None
	
	Integer :: I
    Integer :: J
    Integer :: M
    Integer :: NN,G
	Integer :: NL
    Integer,Parameter :: IS1=4
    Integer,Parameter :: JS1=4
    Integer,Parameter :: IE1=35
    Integer,Parameter :: JE1=35
    Double Precision,Parameter :: R1=16.0d0
    Double Precision,Parameter :: DX1=1.0d0/R1
    Double Precision,Parameter :: DY1=DX1
    Double Precision,Parameter :: CFL1=0.05d0
	Double Precision,Parameter :: CFL2=0.8d0
	Double Precision :: CFL
!    Double Precision,Parameter :: CIGUMA=9.09375d-05
    Double Precision :: DT1
	Double Precision :: DT1_A
	Double Precision :: DT1_B
    Double Precision :: A1
    Double Precision :: B1
    Double Precision :: VEL1
    Double Precision :: VEL2
    Double Precision :: PRE1
    Double Precision :: PRE2
    Double Precision :: DEN1
    Double Precision :: DEN2
    Double Precision :: CZ1
    Double Precision :: SDT
	Double Precision :: t1
	Double Precision :: t2

	Double Precision :: CV_R
	Double Precision :: CP_R
    
    Integer,Parameter :: IS2=4
    Integer,Parameter :: JS2=4
    Integer,Parameter :: IE2=255
    Integer,Parameter :: JE2=255
    Integer,Parameter :: IS1_2=4
    Integer,Parameter :: JS1_2=4
    Integer,Parameter :: IE1_2=IS2+(IE2-IS2+1)/2
    Integer,Parameter :: JE1_2=JS2+(JE2-JS2+1)/2
    Double Precision :: DT2
    Double Precision,Parameter :: DX2=0.5d0*DX1
    Double Precision,Parameter :: DY2=DX2
    
    Integer,Parameter :: IS3=4
    Integer,Parameter :: JS3=4
    Integer,Parameter :: IE3=485
    Integer,Parameter :: JE3=485
    Integer,Parameter :: IS2_3=4
    Integer,Parameter :: JS2_3=4
    Integer,Parameter :: IE2_3=IS3+(IE3-IS3+1)/2
    Integer,Parameter :: JE2_3=JS3+(JE3-JS3+1)/2
    Double Precision :: DT3
    Double Precision,Parameter :: DX3=0.5d0*DX2
    Double Precision,Parameter :: DY3=DX3
 
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxMAT=4
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(MaxM) :: GAM
	Double Precision,Dimension(MaxM) :: PA
	Double Precision,Dimension(IE1+3) :: RR1
    Double Precision,Dimension(IE1+3,JE1+3) :: F1
    Double Precision,Dimension(MaxM,IE1+3,JE1+3) :: P1
	Double Precision,Dimension(MaxM,IE1+3,JE1+3) :: U1
	Double Precision,Dimension(MaxM,IE1+3,JE1+3) :: V1
	Double Precision,Dimension(MaxM,IE1+3,JE1+3) :: C1
	Double Precision,Dimension(MaxM,MaxN,MaxL,IE1+3,JE1+3) :: Q1
	
	Double Precision,Dimension(IE2+3) :: RR2
	Double Precision,Dimension(IE2+3,JE2+3) :: F2
    Double Precision,Dimension(MaxM,IE2+3,JE2+3) :: P2
	Double Precision,Dimension(MaxM,IE2+3,JE2+3) :: U2
	Double Precision,Dimension(MaxM,IE2+3,JE2+3) :: V2
	Double Precision,Dimension(MaxM,IE2+3,JE2+3) :: C2
	Double Precision,Dimension(MaxM,MaxN,MaxL,IE2+3,JE2+3) :: Q2
	
	Double Precision,Dimension(IE3+3) :: RR3
	Double Precision,Dimension(IE3+3,JE3+3) :: F3
    Double Precision,Dimension(MaxM,IE3+3,JE3+3) :: P3
	Double Precision,Dimension(MaxM,IE3+3,JE3+3) :: U3
	Double Precision,Dimension(MaxM,IE3+3,JE3+3) :: V3
	Double Precision,Dimension(MaxM,IE3+3,JE3+3) :: C3
	Double Precision,Dimension(MaxM,MaxN,MaxL,IE3+3,JE3+3) :: Q3
	
	Do J=JS1-3,JE1+3
    Do I=IS1-3,IE1+3
    
        F1(I,J)=0.0d0
        
    End Do 
    End Do

	Call cpu_time(t2)
    
 !   A1=IS1+(IE1-IS1)/2
 !   B1=JS1+(JE1-JS1)/2
    
    A1=Dble(IS1)+Dble(0.5d0)
!	A1=Dble(JS1+R1+60.0d0)
!	B1=Dble(IS1)+Dble(0.5d0)
	B1=Dble(JS1)+R1+1.0d0
!    B1=Dble(JS1+R1+60.0d0)
 !      B1=Dble(JS1+R1+300.0d0)
    
    Do I=IS1+1,IE1-1
        RR1(I)=DX1*(Dble(I-IS1)-0.5d0)
    End Do
    
    Do I=IS2+1,IE2-1
        RR2(I)=DX2*(Dble(I-IS2)-0.5d0)
    End Do
    
    Do I=IS3+1,IE3-1
        RR3(I)=DX3*(Dble(I-IS3)-0.5d0)
    End Do
    
    Do J=JS1+1,JE1
    Do I=IS1+1,IE1
    
        F1(I,J)=DSqrt(((Dble(I)-A1)*DX1)**2.0d0+(Dble(J-B1)*DY1)**2.0d0)-R1*DX1
!	F1(I,J)=Dble(J-20)*DX1
        Write(1,*)I,J,F1(I,J)

    End Do
    End Do

    
	Do I=IS1+1,IE1-1

!		Write(1,*)I,F1(I,JS1+1)
	End Do
    
    VEL1=100.0d0
    VEL2=0.0d0
    PRE1=101300.0d0
    PRE2=101300.0d0

    DEN1=1000.0d0
    DEN2=1.2d0
    
    GAM(1)=4.4d0 !5.0
    GAM(2)=1.4d0
    PA(1)=(6.0d08)/(DEN1*VEL1*VEL1)
    PA(2)=0.0d0
    
    Do J=JS1-3,JE1+3
	Do I=IS1-3,IE1+3
	    P1(1,I,J)=PRE1/(DEN1*VEL1*VEL1) !+2.0d0*CIGUMA
	    Q1(1,1,1,I,J)=DEN1/DEN1
	    P1(2,I,J)=PRE2/(DEN1*VEL1*VEL1)
	    Q1(2,1,1,I,J)=DEN2/DEN1

!	    If(F1(I,J).LE.0.0d0)Then
	        U1(1,I,J)=0.0d0 !-VEL1/VEL1
	        V1(1,I,J)=-VEL1/VEL1
!	    ElseIf(F1(I,J).GT.0.0d0)Then
	        U1(2,I,J)=0.0d0 !-VEL1/VEL1
	        V1(2,I,J)=0.0d0 !-VEL1/VEL1
!	    EndIf
	    
	End Do
	End Do

	Do I=IS1+1,IS1+10
	Do J=JS1+1,JE1-1

!		P1(2,I,J)=PRE2/(DEN1*VEL1*VEL1)*2.0d0

	End Do
	End Do
	
	Do M=1,2
	Do J=JS1-3,JE1+3
	Do I=IS1-3,IE1+3
	 
        Q1(M,1,2,I,J)=Q1(M,1,1,I,J)*U1(M,I,J)
	    Q1(M,1,3,I,J)=Q1(M,1,1,I,J)*V1(M,I,J)
	    Q1(M,1,4,I,J)=(P1(M,I,J)+GAM(M)*PA(M))/(GAM(M)-1.0d0)+0.5d0*Q1(M,1,1,I,J)*(U1(M,I,J)*U1(M,I,J)+V1(M,I,J)*V1(M,I,J))
	    
	End Do
	End Do
	End Do

	
!	Call Boundary(Q1,IS1,IE1,JS1,JE1)
	
!	Call Boundaryf3(F1,IS1,IE1,JS1,JE1)
	
!	Call BoundaryLayerA_B(IS1_2,JS1_2,IE1_2,JE1_2,IS2,JS2,IE2,JE2,Q1,Q2,F1,F2)

!	Call Boundary2(Q2,IS2,IE2,JS2,JE2)

!	Call Boundaryf3(F2,IS2,IE2,JS2,JE2)
	
!	Call BoundaryLayerA_B(IS2_3,JS2_3,IE2_3,JE2_3,IS3,JS3,IE3,JE3,Q2,Q3,F2,F3)

!	A2=Dble(IS2)+0.5d0
!    B2=Dble(JS2)+2.0d0*R1+5.5d0
    
!    Do J=JS2+1,JE2-1
!    Do I=IS2+1,IE2-1
    
!       F2(I,J)=DSqrt(((Dble(I)-A2)*DX2)**2.0d0+(Dble(J-B2)*DY2)**2.0d0)-2.0d0*R1*DX2
!    Write(1,*)I,J,F2(I,J)
        
!    End Do
!    End Do

	NL=2

	NN=1

	Do While(NN.GE.0)
	
	
	Do M=1,2
    Do J=JS1-3,JE1+3
    Do I=IS1-3,IE1+3
    
        C1(M,I,J)=dsqrt(GAM(M)*(P1(M,I,J)+PA(M))/Q1(M,1,1,I,J))
    
    End Do
    End Do
    End Do
    
    CZ1=0.0
    Do M=1,2
    Do J=JS1+1,JE1-1
    Do I=IS1+1,IE1-1
    
        If(Abs(U1(M,I,J)+C1(M,I,J)).GT.CZ1) CZ1=Abs(U1(M,I,J)+C1(M,I,J))
        If(Abs(V1(M,I,J)+C1(M,I,J)).GT.CZ1) CZ1=Abs(V1(M,I,J)+C1(M,I,J))
        If(Abs(U1(M,I,J)-C1(M,I,J)).GT.CZ1) CZ1=Abs(U1(M,I,J)-C1(M,I,J))
        If(Abs(V1(M,I,J)-C1(M,I,J)).GT.CZ1) CZ1=Abs(V1(M,I,J)-C1(M,I,J))
        
    End Do
    End Do
    End Do

	Do I=IS1+1,IE1-1

		If((F1(I,JS1+1).LE.0.0d0).AND.(F1(I+1,JS1+1).GT.0.0d0))CP_R=dble(I-IS1)+dabs(F1(I,JS1+1))/(dabs(F1(I,JS1+1))+dabs(F1(I+1,JS1+1)))-0.5d0
		If(CP_R.EQ.0)CP_R=0.5d0

	End Do

	CV_R=dsqrt(1.0d0-(CP_R/R1)*(CP_R/R1))/(CP_R/R1)

	DT1_A=CFL2*DX1/CV_R
	DT1_B=CFL1*DX1/CZ1


	If(DT1_A.GT.DT1_B)Then
		DT1=DT1_B

	Else


		DT1=DT1_A
	End If

	
    SDT=SDT+DT1
    Write(*,*)NN,DT1,SDT,CFL,CV_R,CZ1,CP_R
    
    Call Calculation(DT1,DX1,DY1,U1,V1,P1,Q1,F1,IS1,IE1,JS1,JE1,GAM,PA,RR1,NN,1,1,SDT,NL,R1,A1,B1)

    Call Boundary(Q1,IS1,IE1,JS1,JE1,F1,DX1)

	Call Boundaryf(F1,IS1,IE1,JS1,JE1,DX1)
	
!    Call LayerBoundaryA_B(IS1_2,JS1_2,IE1_2,JE1_2,IS2,JS2,IE2,JE2,Q1,Q2,F1,F2)

    
    DT2=0.5d0*DT1
    
!    Call Calculation(DT2,DX2,DY2,U2,V2,P2,Q2,F2,IS2,IE2,JS2,JE2,GAM,PA,RR2,NN,2,0,SDT)
    
!    Call Calculation(DT2,DX2,DY2,U2,V2,P2,Q2,F2,IS2,IE2,JS2,JE2,GAM,PA,RR2,NN,2,1,SDT)


!	Call Boundary2(Q2,IS2,IE2,JS2,JE2)

!	Call Boundaryf2(F2,IS2,IE2,JS2,JE2)
	
!    Call LayerBoundaryA_B(IS2_3,JS2_3,IE2_3,JE2_3,IS3,JS3,IE3,JE3,Q2,Q3,F2,F3)

	DT3=0.5d0*DT2
    
!    Call Calculation(DT3,DX3,DY3,U3,V3,P3,Q3,F3,IS3,IE3,JS3,JE3,GAM,PA,RR3,NN,3,0,SDT)
    
!    Call Calculation(DT3,DX3,DY3,U3,V3,P3,Q3,F3,IS3,IE3,JS3,JE3,GAM,PA,RR3,NN,3,0,SDT)

!	Call Calculation(DT3,DX3,DY3,U3,V3,P3,Q3,F3,IS3,IE3,JS3,JE3,GAM,PA,RR3,NN,3,0,SDT)
    
!    Call Calculation(DT3,DX3,DY3,U3,V3,P3,Q3,F3,IS3,IE3,JS3,JE3,GAM,PA,RR3,NN,3,1,SDT)


!	Call Boundary2(Q3,IS3,IE3,JS3,JE3)

!	Call Boundaryf2(F3,IS3,IE3,JS3,JE3)

!    Call BoundaryLayerB_A(IS2_3,JS2_3,IE2_3,JE2_3,IS3,JS3,IE3,JE3,Q2,Q3,F2,F3)


!	Call Boundary2(Q2,IS2,IE2,JS2,JE2)

!	Call Boundaryf2(F2,IS2,IE2,JS2,JE2)

!    Call BoundaryLayerB_A(IS1_2,JS1_2,IE1_2,JE1_2,IS2,JS2,IE2,JE2,Q1,Q2,F1,F2)

    
!    If(F1(IE1-1,JS1+1).LT.0.0d0)EXIT 
 !    If(F1(180,JS1+1).LT.0.0)EXIT
 
!    If(F1(60,JS1+1).LT.0.0)EXIT

!    If(F1(IS1+20,JS1+1).LE.0.0)EXIT

!	If(SDT.GE.0.5)EXIT

	      
	NN=NN+1

	If(NN.EQ.101)Then

		Call cpu_time(t1)

		Write(2,*)t2-t1

	EndIf
    
    End Do
    
    contains
    
	Subroutine Level(DX,DY,U,V,F,LF,IS,IE,JS,JE)
    
    Implicit None
    
    Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Double Precision :: DX
    Double Precision :: DY
    Double Precision :: VA
    Double Precision :: VB
    Double Precision :: VC
    Double Precision :: VD
    Double Precision :: VE
    Double Precision :: FM
    Double Precision :: FP
    Double Precision :: DFDX
    Double Precision :: DFDY
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
    Double Precision,Dimension(:,:) :: F
    Double Precision,Dimension(:,:) :: LF
    Double Precision,Dimension(:,:,:) :: U
	Double Precision,Dimension(:,:,:) :: V
	
         
    Do J=JS+1,JE-1
    Do I=IS+1,IE-1
    
        VA=(F(I-2,J)-F(I-3,J))/DX
        VB=(F(I-1,J)-F(I-2,J))/DX
        VC=(F(I,J)-F(I-1,J))/DX
        VD=(F(I+1,J)-F(I,J))/DX
        VE=(F(I+2,J)-F(I+1,J))/DX 
        Call Weno(VA,VB,VC,VD,VE,FM)
        
        VA=(F(I+3,J)-F(I+2,J))/DX
        VB=(F(I+2,J)-F(I+1,J))/DX
        VC=(F(I+1,J)-F(I,J))/DX
        VD=(F(I,J)-F(I-1,J))/DX
        VE=(F(I-1,J)-F(I-2,J))/DX
        Call Weno(VA,VB,VC,VD,VE,FP)
        
        If(F(I,J).LE.0.0)Then 
        
	        If(U(1,I,J).GE.0.0)Then
	            DFDX=FM
            Else
                DFDX=FP
            EndIf
            
        Else
        
            If(U(2,I,J).GE.0.0)Then
	            DFDX=FM
            Else
                DFDX=FP
            EndIf
            
        EndIf 
        
        VA=(F(I,J-2)-F(I,J-3))/DY
        VB=(F(I,J-1)-F(I,J-2))/DY
        VC=(F(I,J)-F(I,J-1))/DY
        VD=(F(I,J+1)-F(I,J))/DY
        VE=(F(I,J+2)-F(I,J+1))/DY
        Call Weno(VA,VB,VC,VD,VE,FM)
        VA=(F(I,J+3)-F(I,J+2))/DY
        VB=(F(I,J+2)-F(I,J+1))/DY
        VC=(F(I,J+1)-F(I,J))/DY
        VD=(F(I,J)-F(I,J-1))/DY
        VE=(F(I,J-1)-F(I,J-2))/DY
        Call Weno(VA,VB,VC,VD,VE,FP)

        If(F(I,J).LE.0.0)Then
        
	        If(V(1,I,J).GE.0.0)Then
	            DFDY=FM
            Else
                DFDY=FP
            EndIf
            
        Else
        
            If(V(2,I,J).GE.0.0)Then
	            DFDY=FM
            Else
                DFDY=FP
            EndIf
            
        EndIf
        
        If(F(I,J).LE.0.0)Then
        
            LF(I,J)=-U(1,I,J)*DFDX-V(1,I,J)*DFDY
        Else
        
            LF(I,J)=-U(2,I,J)*DFDX-V(2,I,J)*DFDY
            
        EndIf
    
    End Do
    End Do
	
    End Subroutine Level
    
    
    Subroutine Reinitial(DX,DY,F,LF,IS,IE,JS,JE,SF0)
    
    Implicit None
    
    Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Double Precision :: DX
    Double Precision :: DY
    Double Precision :: VA
    Double Precision :: VB
    Double Precision :: VC
    Double Precision :: VD
    Double Precision :: VE
    Double Precision :: FP
    Double Precision :: FM
    Double Precision :: DFDX
    Double Precision :: DFDY
    
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
    Double Precision,Dimension(:,:) :: F
    Double Precision,Dimension(:,:) :: SF0
    Double Precision,Dimension(:,:) :: LF
    
    Do J=JS+1,JE-1
    Do I=IS+1,IE-1
    
        VA=(F(I-2,J)-F(I-3,J))/DX
        VB=(F(I-1,J)-F(I-2,J))/DX
        VC=(F(I,J)-F(I-1,J))/DX
        VD=(F(I+1,J)-F(I,J))/DX
        VE=(F(I+2,J)-F(I+1,J))/DX 
        Call Weno(VA,VB,VC,VD,VE,FM)
        
        VA=(F(I+3,J)-F(I+2,J))/DX
        VB=(F(I+2,J)-F(I+1,J))/DX
        VC=(F(I+1,J)-F(I,J))/DX
        VD=(F(I,J)-F(I-1,J))/DX
        VE=(F(I-1,J)-F(I-2,J))/DX
        Call Weno(VA,VB,VC,VD,VE,FP)
        
        If((FP*SF0(I,J).LE.0.0).AND.(-FP*SF0(I,J).GE.FM*SF0(I,J)))Then
            DFDX=FP
        ElseIf((FM*SF0(I,J).GE.0.0).AND.(FP*SF0(I,J).GE.-FM*SF0(I,J)))Then
            DFDX=FM
        ElseIf((FM*SF0(I,J).LE.0.0).AND.(FP*SF0(I,J).GE.0.0))Then
            DFDX=0.5*(FP+FM)
        EndIf
        
        VA=(F(I,J-2)-F(I,J-3))/DY
        VB=(F(I,J-1)-F(I,J-2))/DY
        VC=(F(I,J)-F(I,J-1))/DY
        VD=(F(I,J+1)-F(I,J))/DY
        VE=(F(I,J+2)-F(I,J+1))/DY
        Call Weno(VA,VB,VC,VD,VE,FM)
        VA=(F(I,J+3)-F(I,J+2))/DY
        VB=(F(I,J+2)-F(I,J+1))/DY
        VC=(F(I,J+1)-F(I,J))/DY
        VD=(F(I,J)-F(I,J-1))/DY
        VE=(F(I,J-1)-F(I,J-2))/DY
        Call Weno(VA,VB,VC,VD,VE,FP)


        If((FP*SF0(I,J).LE.0.0).AND.(-FP*SF0(I,J).GE.FM*SF0(I,J)))Then
            DFDY=FP
        ElseIf((FM*SF0(I,J).GE.0.0).AND.(FP*SF0(I,J).GE.-FM*SF0(I,J)))Then
            DFDY=FM
        ElseIf((FM*SF0(I,J).LE.0.0).AND.(FP*SF0(I,J).GE.0.0))Then
            DFDY=0.5*(FP+FM)
        EndIf
        
        LF(I,J)=(-SF0(I,J)*(SQRT(DFDX*DFDX+DFDY*DFDY)-1.0d0))
        
    End Do
    End Do
    
    End Subroutine Reinitial
    
    
    Subroutine Boundaryf(F,IS,IE,JS,JE,DX)
    
    Implicit None
    
    Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
	Double Precision :: DX
   
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
    Double Precision,Dimension(:,:) :: F
	
	Do I=IS+1,IE-1
!		If(F(I,JS+1).LE.-0.5*DX)Then
!	    F(I,JS)=F(I,JS+1)
!	    F(I,JS-1)=F(I,JS+2)
!	    F(I,JS-2)=F(I,JS+3)
!	    F(I,JS-3)=F(I,JS+4)
		F(I,JS)=2.0d0*F(I,JS+1)-F(I,JS+2)
	    F(I,JS-1)=2.0d0*F(I,JS)-F(I,JS+1)
	    F(I,JS-2)=2.0d0*F(I,JS-1)-F(I,JS)
	    F(I,JS-3)=2.0d0*F(I,JS-2)-F(I,JS-1)
!		Else
!		F(I,JS)=F(I,JS+1)
!	    F(I,JS-1)=F(I,JS+1)
!	    F(I,JS-2)=F(I,JS+1)
!	    F(I,JS-3)=F(I,JS+1)
!		EndIf
	    F(I,JE)=F(I,JE-1)
	    F(I,JE+1)=F(I,JE-1)
	    F(I,JE+2)=F(I,JE-1)
	    F(I,JE+3)=F(I,JE-1)
!        F(I,JE)=2.0*F(I,JE-1)-F(I,JE-2)
!	    F(I,JE+1)=2.0*F(I,JE)-F(I,JE-1)
!	    F(I,JE+2)=2.0*F(I,JE+1)-F(I,JE)
!	    F(I,JE+3)=2.0*F(I,JE+2)-F(I,JE+1)
	End Do
	
	Do J=JS-3,JE+3
	    F(IS,J)=F(IS+1,J)
	    F(IS-1,J)=F(IS+1,J)
	    F(IS-2,J)=F(IS+1,J)
	    F(IS-3,J)=F(IS+1,J)
	    F(IE,J)=F(IE-1,J)
	    F(IE+1,J)=F(IE-1,J)
	    F(IE+2,J)=F(IE-1,J)
	    F(IE+3,J)=F(IE-1,J)
!	    F(IE,J)=2.0*F(IE-1,J)-F(IE-2,J)
!	    F(IE+1,J)=2.0*F(IE,J)-F(IE-1,J)
!	    F(IE+2,J)=2.0*F(IE+1,J)-F(IE,J)
!	    F(IE+3,J)=2.0*F(IE+2,J)-F(IE+1,J)
	    F(IS,J)=F(IS+1,J)
	    F(IS-1,J)=F(IS+2,J)
	    F(IS-2,J)=F(IS+3,J)
	    F(IS-3,J)=F(IS+4,J)
	End Do
	
	End Subroutine Boundaryf


	Subroutine Boundaryf2(F,IS,IE,JS,JE)
    
    Implicit None
    
    Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
   
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
    Double Precision,Dimension(:,:) :: F
	
	Do I=IS+1,IE-1
	    F(I,JS)=F(I,JS+1)
	    F(I,JS-1)=F(I,JS+2)
	    F(I,JS-2)=F(I,JS+3)
	    F(I,JS-3)=F(I,JS+4)
!	    F(I,JE)=2.0*F(I,JE-1)-F(I,JE-2)
!	    F(I,JE+1)=2.0*F(I,JE)-F(I,JE-1)
!	    F(I,JE+2)=2.0*F(I,JE+1)-F(I,JE)
!	    F(I,JE+3)=2.0*F(I,JE+2)-F(I,JE+1)
	End Do
	
	Do J=JS-3,JE+3
	    F(IS,J)=F(IS+1,J)
	    F(IS-1,J)=F(IS+1,J)
	    F(IS-2,J)=F(IS+1,J)
	    F(IS-3,J)=F(IS+1,J)
	    F(IS,J)=F(IS+1,J)
	    F(IS-1,J)=F(IS+2,J)
	    F(IS-2,J)=F(IS+3,J)
	    F(IS-3,J)=F(IS+4,J)
!	    F(IE,J)=2.0*F(IE-1,J)-F(IE-2,J)
!	    F(IE+1,J)=2.0*F(IE,J)-F(IE-1,J)
!	    F(IE+2,J)=2.0*F(IE+1,J)-F(IE,J)
!	    F(IE+3,J)=2.0*F(IE+2,J)-F(IE+1,J)
	End Do
	
	End Subroutine Boundaryf2


	Subroutine Boundaryf3(F,IS,IE,JS,JE)
    
    Implicit None
    
    Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
   
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
    Double Precision,Dimension(:,:) :: F
	
	Do I=IS+1,IE-1
	    F(I,JS)=2.0*F(I,JS+1)-F(I,JS+2)
	    F(I,JS-1)=2.0*F(I,JS)-F(I,JS+1)
	    F(I,JS-2)=2.0*F(I,JS-1)-F(I,JS)
	    F(I,JS-3)=2.0*F(I,JS-2)-F(I,JS-1)
!	    F(I,JE)=F(I,JE-1)
!	    F(I,JE+1)=F(I,JE-1)
!	    F(I,JE+2)=F(I,JE-1)
!	    F(I,JE+3)=F(I,JE-1)
        F(I,JE)=2.0*F(I,JE-1)-F(I,JE-2)
	    F(I,JE+1)=2.0*F(I,JE)-F(I,JE-1)
	    F(I,JE+2)=2.0*F(I,JE+1)-F(I,JE)
	    F(I,JE+3)=2.0*F(I,JE+2)-F(I,JE+1)
	End Do
	
	Do J=JS-3,JE+3
	    F(IS,J)=F(IS+1,J)
	    F(IS-1,J)=F(IS+1,J)
	    F(IS-2,J)=F(IS+1,J)
	    F(IS-3,J)=F(IS+1,J)
!	    F(IE,J)=F(IE-1,J)
!	    F(IE+1,J)=F(IE-1,J)
!	    F(IE+2,J)=F(IE-1,J)
!	    F(IE+3,J)=F(IE-1,J)
	    F(IE,J)=2.0*F(IE-1,J)-F(IE-2,J)
	    F(IE+1,J)=2.0*F(IE,J)-F(IE-1,J)
	    F(IE+2,J)=2.0*F(IE+1,J)-F(IE,J)
	    F(IE+3,J)=2.0*F(IE+2,J)-F(IE+1,J)
	    F(IS,J)=F(IS+1,J)
	    F(IS-1,J)=F(IS+2,J)
	    F(IS-2,J)=F(IS+3,J)
	    F(IS-3,J)=F(IS+4,J)
	End Do
	
	End Subroutine Boundaryf3
    
    
    Subroutine Weno(VA,VB,VC,VD,VE,F)
    
    Implicit None
    
    Double Precision :: VA
    Double Precision :: VB
    Double Precision :: VC
    Double Precision :: VD
    Double Precision :: VE
    Double Precision :: SA
    Double Precision :: SB
    Double Precision :: SC
    Double Precision :: AA
    Double Precision :: AB
    Double Precision :: AC
    Double Precision :: WAA
    Double Precision :: WAB
    Double Precision :: WAC
    Double Precision :: F
    Double Precision,Parameter :: EPU=0.000001
    
        SA=13.0*(VA-2.0*VB+VC)*(VA-2.0*VB+VC)/12.0+(VA-4.0*VB+3.0*VC)*(VA-4.0*VB+3.0*VC)/4.0
        SB=13.0*(VB-2.0*VC+VD)*(VB-2.0*VC+VD)/12.0+(VB-VD)*(VB-VD)/4.0
        SC=13.0*(VC-2.0*VD+VE)*(VC-2.0*VD+VE)/12.0+(3.0*VC-4.0*VD+VE)*(3.0*VC-4.0*VD+VE)/4.0
        AA=1.0/((EPU+SA)*(EPU+SA)*10.0)
        AB=6.0/((EPU+SB)*(EPU+SB)*10.0)
        AC=3.0/((EPU+SC)*(EPU+SC)*10.0)
        WAA=AA/(AA+AB+AC)
        WAB=AB/(AA+AB+AC)
        WAC=AC/(AA+AB+AC)
        F=WAA*(VA/3.0-7.0*VB/6.0+11.0*VC/6.0)+WAB*(-VB/6.0+5.0*VC/6.0+VD/3.0)+WAC*(VC/3.0+5.0*VD/6.0-VE/6.0)
        
    End Subroutine Weno
    
    
    Subroutine Levelset(DT,DX,DY,U,V,Q,F,IS,IE,JS,JE,LL)

	Implicit None
	
	Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    Integer :: G
    Integer :: LL
    Double Precision :: DX
    Double Precision :: DY
    Double Precision :: DT
    Double Precision :: FTO
    Double Precision :: DTA
      
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(:,:) :: F
	Double Precision,allocatable,Dimension(:,:) :: FN
	Double Precision,allocatable,Dimension(:,:) :: LF
	Double Precision,allocatable,Dimension(:,:) :: SF0
    Double Precision,Dimension(:,:,:) :: U
	Double Precision,Dimension(:,:,:) :: V
	Double Precision,Dimension(:,:,:,:,:) :: Q
	
	Allocate(LF(IE+3,JE+3))
	Allocate(FN(IE+3,JE+3))
	Allocate(SF0(IE+3,JE+3))
	
	If(LL.EQ.1)Then
		Call Boundary(Q,IS,IE,JS,JE,F,DX)
	Else
		Call Boundary2(Q,IS,IE,JS,JE)
	EndIf
	
	Do M=1,2
	Do J=JS-3,JE+3
	Do I=IS-3,IE+3
	
!	    U(M,I,J)=Q(M,1,2,I,J)/Q(M,1,1,I,J)
!	    V(M,I,J)=Q(M,1,3,I,J)/Q(M,1,1,I,J)
	
	End Do
	End Do
    End Do
    
    If(LL.EQ.1) Then
		Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf

    Call Level(DX,DY,U,V,F,LF,IS,IE,JS,JE)
    
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
	    FN(I,J)=F(I,J)
	    
	End Do
	End Do 
    
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        F(I,J)=FN(I,J)+DT*LF(I,J)
        
	End Do
	End Do
	
	If(LL.EQ.1) Then
		Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf

    Call Level(DX,DY,U,V,F,LF,IS,IE,JS,JE)
	
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        F(I,J)=(3.0d0*FN(I,J)+F(I,J)+DT*LF(I,J))/4.0d0
            
	End Do
	End Do
	
	If(LL.EQ.1) Then
		Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf

    Call Level(DX,DY,U,V,F,LF,IS,IE,JS,JE)
	
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        F(I,J)=(FN(I,J)+2.0d0*F(I,J)+2.0d0*DT*LF(I,J))/3.0d0
            
	End Do
	End Do
	
	If(LL.EQ.1) Then
		Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf
	
	Do J=JS-3,JE+3
	Do I=IS-3,IE+3
	
        SF0(I,J)=F(I,J)/(dSQRT(F(I,J)*F(I,J)+DX*DX))
        
    End Do
    End Do
    
	
	FTO=1.0
    G=0
    DTA=DT
    Do While(FTO.GT.0.000001)
    
    G=G+1
	
	If(LL.EQ.1) Then
		Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf
	
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        FN(I,J)=F(I,J)
            
	End Do
	End Do
	
    Call Reinitial(DX,DY,F,LF,IS,IE,JS,JE,SF0)
    
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        F(I,J)=FN(I,J)+DTA*LF(I,J)
            
	End Do
	End Do
	
	If(LL.EQ.1) Then
		Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf
	
	Call Reinitial(DX,DY,F,LF,IS,IE,JS,JE,SF0)
	
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        F(I,J)=(3.0d0*FN(I,J)+F(I,J)+DTA*LF(I,J))/4.0d0
            
	End Do
	End Do
	
	If(LL.EQ.1) Then
		Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf
	
	Call Reinitial(DX,DY,F,LF,IS,IE,JS,JE,SF0)
	
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        F(I,J)=(FN(I,J)+2.0d0*F(I,J)+2.0d0*DTA*LF(I,J))/3.0d0
            
	End Do
	End Do
	
	FTO=0.0
    
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
        FTO=FTO+(F(I,J)-FN(I,J))*(F(I,J)-FN(I,J))
        
    End Do
    End Do
    
    FTO=SQRT(FTO/((IE-IS+1)*(JE-JS+1)))
    
    If(G.GE.20)EXIT
    
    End Do    

    End Subroutine Levelset
    
    
	Subroutine heapsort(IS,IE,JS,JE,DX,F,FI,FJ,L,gb,ph)
    
    Implicit None
    
    Integer :: IS
    Integer :: IE
    Integer :: JS
    Integer :: JE
    Integer :: L
    Integer :: max
    Integer :: n
    Integer :: mm
	Integer :: ph
    
    Double precision :: DX
    Double precision :: m
    Double precision :: gb
    
    Integer,Dimension(:) :: FI
    Integer,Dimension(:) :: FJ
    Double Precision,Dimension(:,:) :: F
    
    Double Precision,allocatable,Dimension(:) :: A
    
    Allocate(A(IE*JE))

	L=0

	If(ph.EQ.1)Then

	Do J=JS+1,JE-1
	Do I=IS+1,IE-1

!		If(F(I,J).LE.gb*DX.AND.F(I,J).GE.-DX)Then

		If(F(I,J).LE.gb*DX.AND.F(I,J).GE.0.0d0)Then
		
		    L=L+1

			A(L)=F(I,J)
			FI(L)=I
			FJ(L)=J
	
		End If

	End Do
	End Do

	ElseIf(ph.EQ.0)Then

	Do J=JS+1,JE-1
	Do I=IS+1,IE-1

!		If(F(I,J).GE.-gb*DX.AND.F(I,J).LE.DX)Then

		If(F(I,J).GE.-gb*DX.AND.F(I,J).LE.0.0d0)Then
		
		    L=L+1

			A(L)=-F(I,J)
			FI(L)=I
			FJ(L)=J
	
		End If

	End Do
	End Do

	EndIf

	max=L

	Do While(max.GT.0)

		n=max/2

		Do While(n.GT.0)

			j=j+1

			If(2*n.LT.max)Then
			If(A(2*n).LT.A(2*n+1))Then
				m=A(2*n)
				A(2*n)=A(2*n+1)
				A(2*n+1)=m
				mm=FI(2*n)
				FI(2*n)=FI(2*n+1)
				FI(2*n+1)=mm
				mm=FJ(2*n)
				FJ(2*n)=FJ(2*n+1)
				FJ(2*n+1)=mm
			EndIf
		EndIf
			If(A(n).LT.A(2*n))Then
				m=A(n)
				A(n)=A(2*n)
				A(2*n)=m
				mm=FI(n)
				FI(n)=FI(2*n)
				FI(2*n)=mm
				mm=FJ(n)
				FJ(n)=FJ(2*n)
				FJ(2*n)=mm
			EndIf

		n=n-1

	End Do

	m=A(1)
	A(1)=A(max)
	A(max)=m
	mm=FI(1)
	FI(1)=FI(max)
	FI(max)=mm
	mm=FJ(1)
	FJ(1)=FJ(max)
	FJ(max)=mm

	max=max-1

	End Do

	End Subroutine heapsort


    Subroutine Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI,FJ,NMAX,NN,ph)
    
    Implicit None
    
    Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    Integer :: LL
	Integer :: NMAX
	Integer :: N
	Integer :: NN
	Integer :: ph
    Double Precision :: DX
    Double Precision :: DY
    Double Precision :: DFDX
    Double Precision :: DFDY
    Double Precision :: SITA
    Double Precision,Parameter :: PAI=3.141592653589793
	Double Precision :: iso
    Double Precision :: UVT
    Double Precision :: VVT
    Double Precision :: VNX
    Double Precision :: VNY
    Double Precision :: VTX
    Double Precision :: VTY
	Double Precision :: VTX1
    Double Precision :: VTY1
	Double Precision :: VTX2
    Double Precision :: VTY2
	Double Precision :: VNX1
    Double Precision :: VNY1
	Double Precision :: VNX2
    Double Precision :: VNY2
	Double Precision :: VN1
    Double Precision :: VN2
	Double Precision :: VNI1
    Double Precision :: VNI2
	Double Precision :: PI1
    Double Precision :: PI2
	Double Precision :: LO1
    Double Precision :: LO2
	Double Precision :: VNI
	Double Precision :: PI
    Double Precision :: LO
	Double Precision,Parameter :: CIGUMA=0.0d0 !9.09375d-05
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(:) :: GAM
	Double Precision,Dimension(:) :: PA
	Double Precision,Dimension(:) :: RR
	Double Precision,Dimension(:,:) :: F
	Integer,Dimension(:) :: FI
	Integer,Dimension(:) :: FJ
	Double Precision,allocatable,Dimension(:,:) :: NX
	Double Precision,allocatable,Dimension(:,:) :: NY
	Double Precision,allocatable,Dimension(:,:) :: K
    Double Precision,Dimension(:,:,:) :: U
	Double Precision,Dimension(:,:,:) :: V
	Double Precision,Dimension(:,:,:) :: P
	Double Precision,Dimension(:,:,:) :: S
	Double Precision,Dimension(:,:,:,:,:) :: Q
	Double Precision,allocatable,Dimension(:,:,:) :: C
	
	Allocate(NX(IE+3,JE+3))
	Allocate(NY(IE+3,JE+3))
	Allocate(K(IE+3,JE+3))
	Allocate(C(MaxM,IE+3,JE+3))

	iso=0.0d0

	If(F(IS+1,JS+1).LE.0.0d0)iso=1.0d0
    
    If(LL.EQ.1)Then
		 Call Boundaryf(F,IS,IE,JS,JE,DX)
	Else
		Call Boundaryf2(F,IS,IE,JS,JE)
	EndIf
    
    Do J=JS-2,JE+2
	Do I=IS-2,IE+2
        DFDX=(0.5d0*(F(I+1,J)+F(I,J))-0.5d0*(F(I,J)+F(I-1,J)))/DX
        DFDY=(0.5d0*(F(I,J+1)+F(I,J))-0.5d0*(F(I,J)+F(I,J-1)))/DY
        NX(I,J)=DFDX/(SQRT(DFDX*DFDX+DFDY*DFDY))
        NY(I,J)=DFDY/(SQRT(DFDX*DFDX+DFDY*DFDY))
        If((DFDX.EQ.0.0).AND.(DFDY.EQ.0.0))NX(I,J)=0.0d0
        If((DFDX.EQ.0.0).AND.(DFDY.EQ.0.0))NY(I,J)=0.0d0
    End Do
    End Do

	Do I=IS+1,IE-1

!		NX(I,JS+1)=1.0d0
!		NY(I,JS+1)=0.0d0

	End Do
    
    Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	    If(Abs(F(I,J)).LE.(5.0d0*DX))Then
	    K(I,J)=(0.5d0*(NX(I+1,J)+NX(I,J))-0.5d0*(NX(I,J)+NX(I-1,J)))/DX+(0.5d0*(NY(I,J+1)+NY(I,J))-0.5d0*(NY(I,J)+NY(I,J-1)))/DY
	    K(I,J)=(0.5d0*(NX(I+1,J)+NX(I,J))-0.5d0*(NX(I,J)+NX(I-1,J)))/DX+(0.5d0*(NY(I,J+1)+NY(I,J))-0.5d0*(NY(I,J)+NY(I,J-1)))/DY+NX(I,J)/RR(I)
	    Else
	    K(I,J)=0.0d0
	    EndIf 
	End Do
	End Do

	
	Do M=1,2
    Do J=JS-3,JE+3
	Do I=IS-3,IE+3

        U(M,I,J)=Q(M,1,2,I,J)/Q(M,1,1,I,J)
        V(M,I,J)=Q(M,1,3,I,J)/Q(M,1,1,I,J)
        P(M,I,J)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J)-0.5d0*Q(M,1,1,I,J)*(U(M,I,J)*U(M,I,J)+V(M,I,J)*V(M,I,J)))-GAM(M)*PA(M)
        S(M,I,J)=(P(M,I,J)+PA(M))/(Q(M,1,1,I,J)**GAM(M))
	C(M,I,J)=dsqrt(GAM(M)*(P(M,I,J)+PA(M))/Q(M,1,1,I,J))

	End Do
	End Do
	End Do
	


    	LL=0

	Do N=1,NMAX

    	LL=LL+1
   

	If(NX(FI(N),FJ(N)).GE.0.0.AND.NY(FI(N),FJ(N)).GE.0.0)Then
	
!	If(F(FI(N),FJ(N)).GT.-iso*DX)Then

	If(ph.EQ.1)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0
!			If(F(FI(N),FJ(N)).GE.0.0)Then

!		    	U(1,FI(N),FJ(N))=((2.0d0*U(1,FI(N)-1,FJ(N))-U(1,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(1,FI(N),FJ(N)-1)-U(1,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=((2.0d0*V(1,FI(N)-1,FJ(N))-V(1,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(1,FI(N),FJ(N)-1)-V(1,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			U(1,FI(N),FJ(N))=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))

!			EndIf

!		    	S(1,FI(N),FJ(N))=((2.0d0*S(1,FI(N)-1,FJ(N))-S(1,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(1,FI(N),FJ(N)-1)-S(1,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!            		P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N)) !+CIGUMA*K(FI(N),FJ(N))

			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(1,FI(N)-2,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(2,FI(N)+1,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(1,FI(N)-2,FJ(N))-P(2,FI(N)+1,FJ(N))+Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))*VN1+Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))*VN2)/(Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))+Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N)))
			PI1=(Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))*Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))*(VN1-VN2)+Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))*P(1,FI(N)-2,FJ(N))+Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))*P(2,FI(N)+1,FJ(N)))/(Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))+Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N)))
			LO1=((PI1+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI1-0.5d0*(PI1-P(1,FI(N)-2,FJ(N))))*Q(1,1,1,FI(N)-2,FJ(N))/((P(1,FI(N)-2,FJ(N))+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N)-2,FJ(N))+0.5d0*(PI1-P(1,FI(N)-2,FJ(N))))
			
			VN1=V(1,FI(N),FJ(N)-2) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(2,FI(N),FJ(N)+1) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(1,FI(N),FJ(N)-2)-P(2,FI(N),FJ(N)+1)+Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)*VN1+Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)*VN2)/(Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)+Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1))
			PI2=(Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)*Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)*(VN1-VN2)+Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)*P(1,FI(N),FJ(N)-2)+Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)*P(2,FI(N),FJ(N)+1))/(Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)+Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2))
			LO2=((PI2+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI2-0.5d0*(PI2-P(1,FI(N),FJ(N)-2)))*Q(1,1,1,FI(N),FJ(N)-2)/((P(1,FI(N),FJ(N)-2)+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N),FJ(N)-2)+0.5d0*(PI2-P(1,FI(N),FJ(N)-2)))

			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(1,FI(N),FJ(N))=VNX+VTX
			V(1,FI(N),FJ(N))=VNY+VTY
			P(1,FI(N),FJ(N))=PI
			S(1,FI(N),FJ(N))=(PI+PA(1))/(LO**GAM(1))

			Else

!			VNX=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(1,FI(N),FJ(N))=VNX+VTX
!			V(1,FI(N),FJ(N))=VNY+VTY
!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))
!			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			U(1,FI(N),FJ(N))=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	V(1,FI(N),FJ(N))=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			If(F(FI(N),FJ(N)).GT.0.0d0)Then
			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))+CIGUMA*K(FI(N),FJ(N))
			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			EndIf

			EndIf
			

!		ElseIf(F(FI(N),FJ(N)).LT.0.0)Then
		ElseIf(ph.EQ.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0
!			If(F(FI(N),FJ(N)).LT.0.0)Then

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

!		    	UVT=((2.0d0*U(2,FI(N)+1,FJ(N))-U(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(2,FI(N),FJ(N)+1)-U(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!		    	VVT=((2.0d0*V(2,FI(N)+1,FJ(N))-V(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(2,FI(N),FJ(N)+1)-V(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		    	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		    	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
		
!	       		U(2,FI(N),FJ(N))=VNX+VTX
!		    	V(2,FI(N),FJ(N))=VNY+VTY

!			P(2,FI(N),FJ(N))=((2.0d0*P(2,FI(N)+1,FJ(N))-P(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*P(2,FI(N),FJ(N)+1)-P(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=(P(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
!			EndIf
!		    	S(2,FI(N),FJ(N))=((2.0d0*S(2,FI(N)+1,FJ(N))-S(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(2,FI(N),FJ(N)+1)-S(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(1,FI(N)-1,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(2,FI(N)+2,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(1,FI(N)-1,FJ(N))-P(2,FI(N)+2,FJ(N))+Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))*VN1+Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))*VN2)/(Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))+Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N)))
			PI1=(Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))*Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))*(VN1-VN2)+Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))*P(1,FI(N)-1,FJ(N))+Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))*P(2,FI(N)+2,FJ(N)))/(Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))+Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N)))
			LO1=((PI1+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI1-0.5d0*(PI1-P(2,FI(N)+2,FJ(N))))*Q(2,1,1,FI(N)+2,FJ(N))/((P(2,FI(N)+2,FJ(N))+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N)+2,FJ(N))+0.5d0*(PI1-P(2,FI(N)+2,FJ(N))))
		
			VN1=V(1,FI(N),FJ(N)-1) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(2,FI(N),FJ(N)+2) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(1,FI(N),FJ(N)-1)-P(2,FI(N),FJ(N)+2)+Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)*VN1+Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)*VN2)/(Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)+Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2))
			PI2=(Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)*Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)*(VN1-VN2)+Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)*P(1,FI(N),FJ(N)-1)+Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)*P(2,FI(N),FJ(N)+2))/(Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)+Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1))
			LO2=((PI2+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI2-0.5d0*(PI2-P(2,FI(N),FJ(N)+2)))*Q(2,1,1,FI(N),FJ(N)+2)/((P(2,FI(N),FJ(N)+2)+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N),FJ(N)+2)+0.5d0*(PI2-P(2,FI(N),FJ(N)+2)))
	
			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=PI
			S(2,FI(N),FJ(N))=(PI+PA(2))/(LO**GAM(2))

			Else

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(2,FI(N),FJ(N))=VNX+VTX
!			V(2,FI(N),FJ(N))=VNY+VTY
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
!			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=(P(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

			EndIf

		EndIf
		
	ElseIf(NX(FI(N),FJ(N)).LE.0.0.AND.NY(FI(N),FJ(N)).GE.0.0)Then
	
!	    	If(F(FI(N),FJ(N)).GT.-iso*DX)Then
			If(ph.EQ.1)Then
			SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		    	If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0
!			If(F(FI(N),FJ(N)).GE.0.0)Then

!			U(1,FI(N),FJ(N))=((2.0d0*U(1,FI(N)+1,FJ(N))-U(1,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(1,FI(N),FJ(N)-1)-U(1,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=((2.0d0*V(1,FI(N)+1,FJ(N))-V(1,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(1,FI(N),FJ(N)-1)-V(1,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			U(1,FI(N),FJ(N))=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))
!			EndIf
!		    	S(1,FI(N),FJ(N))=((2.0d0*S(1,FI(N)+1,FJ(N))-S(1,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(1,FI(N),FJ(N)-1)-S(1,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

!            		P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N)) !+CIGUMA*K(FI(N),FJ(N))

			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(2,FI(N)-1,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(1,FI(N)+2,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(2,FI(N)-1,FJ(N))-P(1,FI(N)+1,FJ(N))+Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))*VN1+Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))*VN2)/(Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))+Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N)))
			PI1=(Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))*Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))*(VN1-VN2)+Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))*P(2,FI(N)-1,FJ(N))+Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))*P(1,FI(N)+2,FJ(N)))/(Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))+Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N)))
			LO1=((PI1+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI1-0.5d0*(PI1-P(1,FI(N)+2,FJ(N))))*Q(1,1,1,FI(N)+2,FJ(N))/((P(1,FI(N)+2,FJ(N))+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N)+2,FJ(N))+0.5d0*(PI1-P(1,FI(N)+2,FJ(N))))
			
			VN1=V(1,FI(N),FJ(N)-2) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(2,FI(N),FJ(N)+1) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(1,FI(N),FJ(N)-2)-P(2,FI(N),FJ(N)+1)+Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)*VN1+Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)*VN2)/(Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)+Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1))
			PI2=(Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)*Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)*(VN1-VN2)+Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)*P(1,FI(N),FJ(N)-2)+Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2)*P(2,FI(N),FJ(N)+1))/(Q(2,1,1,FI(N),FJ(N)+1)*C(2,FI(N),FJ(N)+1)+Q(1,1,1,FI(N),FJ(N)-2)*C(1,FI(N),FJ(N)-2))
			LO2=((PI2+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI2-0.5d0*(PI2-P(1,FI(N),FJ(N)-2)))*Q(1,1,1,FI(N),FJ(N)-2)/((P(1,FI(N),FJ(N)-2)+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N),FJ(N)-2)+0.5d0*(PI2-P(1,FI(N),FJ(N)-2)))
	
			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(1,FI(N),FJ(N))=VNX+VTX
			V(1,FI(N),FJ(N))=VNY+VTY
			P(1,FI(N),FJ(N))=PI
			S(1,FI(N),FJ(N))=(PI+PA(1))/(LO**GAM(1))

			Else

!			VNX=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(1,FI(N),FJ(N))=VNX+VTX
!			V(1,FI(N),FJ(N))=VNY+VTY
!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))
!			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			U(1,FI(N),FJ(N))=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	V(1,FI(N),FJ(N))=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			If(F(FI(N),FJ(N)).GT.0.0d0)Then
			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))+CIGUMA*K(FI(N),FJ(N))
			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			EndIf

			EndIf

!		ElseIf(F(FI(N),FJ(N)).LT.0.0)Then
		ElseIf(ph.EQ.0)Then
			SITA=dACOS(ABS(NX(FI(N),FJ(N))))
			If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0
			If(F(FI(N),FJ(N)).LT.0.0)Then

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

!	        	UVT=((2.0d0*U(2,FI(N)-1,FJ(N))-U(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(2,FI(N),FJ(N)+1)-U(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!		    	VVT=((2.0d0*V(2,FI(N)-1,FJ(N))-V(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(2,FI(N),FJ(N)+1)-V(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		    	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		    	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
		
!	       		U(2,FI(N),FJ(N))=VNX+VTX
!		    	V(2,FI(N),FJ(N))=VNY+VTY

!			P(2,FI(N),FJ(N))=((2.0d0*P(2,FI(N)-1,FJ(N))-P(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*P(2,FI(N),FJ(N)+1)-P(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=(P(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
			EndIf
!		    	S(2,FI(N),FJ(N))=((2.0d0*S(2,FI(N)-1,FJ(N))-S(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(2,FI(N),FJ(N)+1)-S(2,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)


			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(2,FI(N)-2,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(1,FI(N)+1,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(2,FI(N)-2,FJ(N))-P(1,FI(N)+1,FJ(N))+Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))*VN1+Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))*VN2)/(Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))+Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N)))
			PI1=(Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))*Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))*(VN1-VN2)+Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))*P(2,FI(N)-2,FJ(N))+Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))*P(1,FI(N)+1,FJ(N)))/(Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))+Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N)))
			LO1=((PI1+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI1-0.5d0*(PI1-P(2,FI(N)-2,FJ(N))))*Q(2,1,1,FI(N)-2,FJ(N))/((P(2,FI(N)-2,FJ(N))+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N)-2,FJ(N))+0.5d0*(PI1-P(2,FI(N)-2,FJ(N))))
			
			VN1=V(1,FI(N),FJ(N)-1) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(2,FI(N),FJ(N)+2) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(1,FI(N),FJ(N)-1)-P(2,FI(N),FJ(N)+2)+Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)*VN1+Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)*VN2)/(Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)+Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2))
			PI2=(Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)*Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)*(VN1-VN2)+Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)*P(1,FI(N),FJ(N)-1)+Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1)*P(2,FI(N),FJ(N)+2))/(Q(2,1,1,FI(N),FJ(N)+2)*C(2,FI(N),FJ(N)+2)+Q(1,1,1,FI(N),FJ(N)-1)*C(1,FI(N),FJ(N)-1))
			LO2=((PI2+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI2-0.5d0*(PI2-P(2,FI(N),FJ(N)+2)))*Q(2,1,1,FI(N),FJ(N)+2)/((P(2,FI(N),FJ(N)+2)+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N),FJ(N)+2)+0.5d0*(PI2-P(2,FI(N),FJ(N)+2)))
	
			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=PI
			S(2,FI(N),FJ(N))=(PI+PA(2))/(LO**GAM(2))

			Else

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(2,FI(N),FJ(N))=VNX+VTX
!			V(2,FI(N),FJ(N))=VNY+VTY
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
!			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=(P(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

			EndIf

		EndIf
		
	ElseIf(NX(FI(N),FJ(N)).LE.0.0.AND.NY(FI(N),FJ(N)).LE.0.0)Then
	
!	    If(F(FI(N),FJ(N)).GT.-iso*DX)Then
		If(ph.EQ.1)Then
		    SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		    If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0
			If(F(FI(N),FJ(N)).GE.0.0)Then

!			U(1,FI(N),FJ(N))=((2.0d0*U(1,FI(N)+1,FJ(N))-U(1,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(1,FI(N),FJ(N)+1)-U(1,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=((2.0d0*V(1,FI(N)+1,FJ(N))-V(1,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(1,FI(N),FJ(N)+1)-V(1,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			U(1,FI(N),FJ(N))=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))
			EndIf
!		    	S(1,FI(N),FJ(N))=((2.0d0*S(1,FI(N)+1,FJ(N))-S(1,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(1,FI(N),FJ(N)+1)-S(1,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

 !           		P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N)) !+CIGUMA*K(FI(N),FJ(N))

			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(2,FI(N)-1,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(1,FI(N)+2,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(2,FI(N)-1,FJ(N))-P(1,FI(N)+1,FJ(N))+Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))*VN1+Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))*VN2)/(Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))+Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N)))
			PI1=(Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))*Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))*(VN1-VN2)+Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))*P(2,FI(N)-1,FJ(N))+Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N))*P(1,FI(N)+2,FJ(N)))/(Q(1,1,1,FI(N)+2,FJ(N))*C(1,FI(N)+2,FJ(N))+Q(2,1,1,FI(N)-1,FJ(N))*C(2,FI(N)-1,FJ(N)))
			LO1=((PI1+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI1-0.5d0*(PI1-P(1,FI(N)+2,FJ(N))))*Q(1,1,1,FI(N)+2,FJ(N))/((P(1,FI(N)+2,FJ(N))+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N)+2,FJ(N))+0.5d0*(PI1-P(1,FI(N)+2,FJ(N))))
			
			VN1=V(2,FI(N),FJ(N)-1) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(1,FI(N),FJ(N)+2) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(2,FI(N),FJ(N)-1)-P(1,FI(N),FJ(N)+2)+Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)*VN1+Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)*VN2)/(Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)+Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2))
			PI2=(Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)*Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)*(VN1-VN2)+Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)*P(2,FI(N),FJ(N)-1)+Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)*P(1,FI(N),FJ(N)+2))/(Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)+Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1))
			LO2=((PI2+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI2-0.5d0*(PI2-P(1,FI(N),FJ(N)+2)))*Q(1,1,1,FI(N),FJ(N)+2)/((P(1,FI(N),FJ(N)+2)+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N),FJ(N)+2)+0.5d0*(PI2-P(1,FI(N),FJ(N)+2)))
	
			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(1,FI(N),FJ(N))=VNX+VTX
			V(1,FI(N),FJ(N))=VNY+VTY
			P(1,FI(N),FJ(N))=PI
			S(1,FI(N),FJ(N))=(PI+PA(1))/(LO**GAM(1))

			Else

!			VNX=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(1,FI(N),FJ(N))=VNX+VTX
!			V(1,FI(N),FJ(N))=VNY+VTY
!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))
!			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

			U(1,FI(N),FJ(N))=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	V(1,FI(N),FJ(N))=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			If(F(FI(N),FJ(N)).GT.0.0d0)Then
			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))+CIGUMA*K(FI(N),FJ(N))
			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			EndIf

			EndIf
			
!		ElseIf(F(FI(N),FJ(N)).LT.0.0)Then
		ElseIf(ph.EQ.0)Then
			SITA=dACOS(ABS(NX(FI(N),FJ(N))))
     			If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0
			If(F(FI(N),FJ(N)).LT.0.0)Then

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

!			UVT=((2.0d0*U(2,FI(N)-1,FJ(N))-U(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(2,FI(N),FJ(N)-1)-U(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!		    	VVT=((2.0d0*V(2,FI(N)-1,FJ(N))-V(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(2,FI(N),FJ(N)-1)-V(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		    	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

!		    	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
		
!	       		U(2,FI(N),FJ(N))=VNX+VTX
!		    	V(2,FI(N),FJ(N))=VNY+VTY

!			P(2,FI(N),FJ(N))=((2.0d0*P(2,FI(N)-1,FJ(N))-P(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*P(2,FI(N),FJ(N)-1)-P(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=(P(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
			EndIf
!		    	S(2,FI(N),FJ(N))=((2.0d0*S(2,FI(N)-1,FJ(N))-S(2,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(2,FI(N),FJ(N)-1)-S(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(2,FI(N)-2,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(1,FI(N)+1,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(2,FI(N)-2,FJ(N))-P(1,FI(N)+1,FJ(N))+Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))*VN1+Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))*VN2)/(Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))+Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N)))
			PI1=(Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))*Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))*(VN1-VN2)+Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))*P(2,FI(N)-2,FJ(N))+Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N))*P(1,FI(N)+1,FJ(N)))/(Q(1,1,1,FI(N)+1,FJ(N))*C(1,FI(N)+1,FJ(N))+Q(2,1,1,FI(N)-2,FJ(N))*C(2,FI(N)-2,FJ(N)))
			LO1=((PI1+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI1-0.5d0*(PI1-P(2,FI(N)-2,FJ(N))))*Q(2,1,1,FI(N)-2,FJ(N))/((P(2,FI(N)-2,FJ(N))+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N)-2,FJ(N))+0.5d0*(PI1-P(2,FI(N)-2,FJ(N))))
			
			VN1=V(2,FI(N),FJ(N)-2) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(1,FI(N),FJ(N)+1) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(2,FI(N),FJ(N)-2)-P(1,FI(N),FJ(N)+1)+Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)*VN1+Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)*VN2)/(Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)+Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1))
			PI2=(Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)*Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)*(VN1-VN2)+Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)*P(2,FI(N),FJ(N)-2)+Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)*P(1,FI(N),FJ(N)+1))/(Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)+Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2))
			LO2=((PI2+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI2-0.5d0*(PI2-P(2,FI(N),FJ(N)-2)))*Q(2,1,1,FI(N),FJ(N)-2)/((P(2,FI(N),FJ(N)-2)+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N),FJ(N)-2)+0.5d0*(PI2-P(2,FI(N),FJ(N)-2)))
	
			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=PI
			S(2,FI(N),FJ(N))=(PI+PA(2))/(LO**GAM(2))

			Else

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(2,FI(N),FJ(N))=VNX+VTX
!			V(2,FI(N),FJ(N))=VNY+VTY
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
!			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			UVT=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=(P(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			EndIf

	EndIf
		
	ElseIf(NX(FI(N),FJ(N)).GE.0.0.AND.NY(FI(N),FJ(N)).LE.0.0)Then
	
!	If(F(FI(N),FJ(N)).GT.-iso*DX)Then
	If(ph.EQ.1)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0
			If(F(FI(N),FJ(N)).GE.0.0)Then
		  
!			U(1,FI(N),FJ(N))=((2.0d0*U(1,FI(N)-1,FJ(N))-U(1,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(1,FI(N),FJ(N)+1)-U(1,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=((2.0d0*V(1,FI(N)-1,FJ(N))-V(1,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(1,FI(N),FJ(N)+1)-V(1,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			U(1,FI(N),FJ(N))=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		    	V(1,FI(N),FJ(N))=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))
			EndIf
!		    	S(1,FI(N),FJ(N))=((2.0d0*S(1,FI(N)-1,FJ(N))-S(1,FI(N)-2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(1,FI(N),FJ(N)+1)-S(1,FI(N),FJ(N)+2))*SITA)/(0.5d0*PAI)
!			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

!            		P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N)) !+CIGUMA*K(FI(N),FJ(N))

			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(1,FI(N)-2,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(2,FI(N)+1,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(1,FI(N)-2,FJ(N))-P(2,FI(N)+1,FJ(N))+Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))*VN1+Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))*VN2)/(Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))+Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N)))
			PI1=(Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))*Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))*(VN1-VN2)+Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))*P(1,FI(N)-2,FJ(N))+Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N))*P(2,FI(N)+1,FJ(N)))/(Q(2,1,1,FI(N)+1,FJ(N))*C(2,FI(N)+1,FJ(N))+Q(1,1,1,FI(N)-2,FJ(N))*C(1,FI(N)-2,FJ(N)))
			LO1=((PI1+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI1-0.5d0*(PI1-P(1,FI(N)-2,FJ(N))))*Q(1,1,1,FI(N)-2,FJ(N))/((P(1,FI(N)-2,FJ(N))+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N)-2,FJ(N))+0.5d0*(PI1-P(1,FI(N)-2,FJ(N))))
		
			VN1=V(2,FI(N),FJ(N)-1) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(1,FI(N),FJ(N)+2) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(2,FI(N),FJ(N)-1)-P(1,FI(N),FJ(N)+2)+Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)*VN1+Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)*VN2)/(Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)+Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2))
			PI2=(Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)*Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)*(VN1-VN2)+Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)*P(2,FI(N),FJ(N)-1)+Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1)*P(1,FI(N),FJ(N)+2))/(Q(1,1,1,FI(N),FJ(N)+2)*C(1,FI(N),FJ(N)+2)+Q(2,1,1,FI(N),FJ(N)-1)*C(2,FI(N),FJ(N)-1))
			LO2=((PI2+GAM(1)*PA(1))/(GAM(1)-1.0d0)+PI2-0.5d0*(PI2-P(1,FI(N),FJ(N)+2)))*Q(1,1,1,FI(N),FJ(N)+2)/((P(1,FI(N),FJ(N)+2)+GAM(1)*PA(1))/(GAM(1)-1.0d0)+P(1,FI(N),FJ(N)+2)+0.5d0*(PI2-P(1,FI(N),FJ(N)+2)))
	
			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(1,FI(N),FJ(N))=VNX+VTX
			V(1,FI(N),FJ(N))=VNY+VTY
			P(1,FI(N),FJ(N))=PI
			S(1,FI(N),FJ(N))=(PI+PA(1))/(LO**GAM(1))

			Else

!			VNX=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(2,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(2,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(1,FI(N),FJ(N))=VNX+VTX
!			V(1,FI(N),FJ(N))=VNY+VTY
!			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))
!			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

			U(1,FI(N),FJ(N))=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
		   	V(1,FI(N),FJ(N))=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			If(F(FI(N),FJ(N)).GT.0.0d0)Then
			P(1,FI(N),FJ(N))=P(2,FI(N),FJ(N))+CIGUMA*K(FI(N),FJ(N))
			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			EndIf

			EndIf
			
!	ElseIf(F(FI(N),FJ(N)).LT.0.0)Then
	ElseIf(ph.EQ.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0
			If(F(FI(N),FJ(N)).LT.0.0)Then

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

!			UVT=((2.0d0*U(2,FI(N)+1,FJ(N))-U(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*U(2,FI(N),FJ(N)-1)-U(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!		    	VVT=((2.0d0*V(2,FI(N)+1,FJ(N))-V(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*V(2,FI(N),FJ(N)-1)-V(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		    	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		    	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
		
!	       		U(2,FI(N),FJ(N))=VNX+VTX
!		    	V(2,FI(N),FJ(N))=VNY+VTY

!			P(2,FI(N),FJ(N))=((2.0d0*P(2,FI(N)+1,FJ(N))-P(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*P(2,FI(N),FJ(N)-1)-P(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=(P(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
			EndIf
!		    	S(2,FI(N),FJ(N))=((2.0d0*S(2,FI(N)+1,FJ(N))-S(2,FI(N)+2,FJ(N)))*(0.5d0*PAI-SITA)+(2.0d0*S(2,FI(N),FJ(N)-1)-S(2,FI(N),FJ(N)-2))*SITA)/(0.5d0*PAI)
!			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)


			If(NN.EQ.1)Then !If(Abs(F(FI(N),FJ(N))).LE.DX)Then

			VN1=U(1,FI(N)-1,FJ(N)) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=U(2,FI(N)+2,FJ(N)) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI1=(P(1,FI(N)-1,FJ(N))-P(2,FI(N)+2,FJ(N))+Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))*VN1+Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))*VN2)/(Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))+Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N)))
			PI1=(Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))*Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))*(VN1-VN2)+Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))*P(1,FI(N)-1,FJ(N))+Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N))*P(2,FI(N)+2,FJ(N)))/(Q(2,1,1,FI(N)+2,FJ(N))*C(2,FI(N)+2,FJ(N))+Q(1,1,1,FI(N)-1,FJ(N))*C(1,FI(N)-1,FJ(N)))
			LO1=((PI1+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI1-0.5d0*(PI1-P(2,FI(N)+2,FJ(N))))*Q(2,1,1,FI(N)+2,FJ(N))/((P(2,FI(N)+2,FJ(N))+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N)+2,FJ(N))+0.5d0*(PI1-P(2,FI(N)+2,FJ(N))))
		
			VN1=V(2,FI(N),FJ(N)-2) !dsqrt(VNX1*VNX1+VNY1*VNY1)
			VN2=V(1,FI(N),FJ(N)+1) !dsqrt(VNX2*VNX2+VNY2*VNY2)
			VNI2=(P(2,FI(N),FJ(N)-2)-P(1,FI(N),FJ(N)+1)+Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)*VN1+Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)*VN2)/(Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)+Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1))
			PI2=(Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)*Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)*(VN1-VN2)+Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)*P(2,FI(N),FJ(N)-2)+Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2)*P(1,FI(N),FJ(N)+1))/(Q(1,1,1,FI(N),FJ(N)+1)*C(1,FI(N),FJ(N)+1)+Q(2,1,1,FI(N),FJ(N)-2)*C(2,FI(N),FJ(N)-2))
			LO2=((PI2+GAM(2)*PA(2))/(GAM(2)-1.0d0)+PI2-0.5d0*(PI2-P(2,FI(N),FJ(N)-2)))*Q(2,1,1,FI(N),FJ(N)-2)/((P(2,FI(N),FJ(N)-2)+GAM(2)*PA(2))/(GAM(2)-1.0d0)+P(2,FI(N),FJ(N)-2)+0.5d0*(PI2-P(2,FI(N),FJ(N)-2)))

			PI=(PI1*(0.5d0*PAI-SITA)+PI2*SITA)/(0.5d0*PAI)
			LO=(LO1*(0.5d0*PAI-SITA)+LO2*SITA)/(0.5d0*PAI)
			VNX=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
			VNY=(VNI1*NX(FI(N),FJ(N))+VNI2*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))

			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=PI
			S(2,FI(N),FJ(N))=(PI+PA(2))/(LO**GAM(2))

			Else

!			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
!		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
!		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
!			U(2,FI(N),FJ(N))=VNX+VTX
!			V(2,FI(N),FJ(N))=VNY+VTY
!			P(2,FI(N),FJ(N))=P(1,FI(N),FJ(N))
!			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			VNX=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VNY=(U(1,FI(N),FJ(N))*NX(FI(N),FJ(N))+V(1,FI(N),FJ(N))*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			UVT=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VVT=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
		   	VTX=UVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NX(FI(N),FJ(N))
		    	VTY=VVT-(UVT*NX(FI(N),FJ(N))+VVT*NY(FI(N),FJ(N)))*NY(FI(N),FJ(N))
			U(2,FI(N),FJ(N))=VNX+VTX
			V(2,FI(N),FJ(N))=VNY+VTY
			P(2,FI(N),FJ(N))=(P(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+P(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

			EndIf

	EndIf
		
	EndIf
	
	End Do

	If(NN.EQ.1)Then

!	Call heapsort(IS,IE,JS,JE,DX,F,FI,FJ,NMAX,1.0d0)

	LL=0

	Do N=1,NMAX

    	LL=LL+1
   

	If(NX(FI(N),FJ(N)).GE.0.0.AND.NY(FI(N),FJ(N)).GE.0.0)Then
	
		If(F(FI(N),FJ(N)).LE.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0

			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			U(1,FI(N),FJ(N))=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			V(1,FI(N),FJ(N))=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)			

		ElseIf(F(FI(N),FJ(N)).GT.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0

			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			U(2,FI(N),FJ(N))=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			V(2,FI(N),FJ(N))=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

		EndIf
		
	ElseIf(NX(FI(N),FJ(N)).LE.0.0.AND.NY(FI(N),FJ(N)).GE.0.0)Then
	
	    	If(F(FI(N),FJ(N)).LE.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0

			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			U(1,FI(N),FJ(N))=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			V(1,FI(N),FJ(N))=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

		ElseIf(F(FI(N),FJ(N)).GT.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0

			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			U(2,FI(N),FJ(N))=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			V(2,FI(N),FJ(N))=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)

		EndIf
		
	ElseIf(NX(FI(N),FJ(N)).LE.0.0.AND.NY(FI(N),FJ(N)).LE.0.0)Then
	
		If(F(FI(N),FJ(N)).LE.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0

			S(1,FI(N),FJ(N))=(S(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			U(1,FI(N),FJ(N))=(U(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			V(1,FI(N),FJ(N))=(V(1,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			
		ElseIf(F(FI(N),FJ(N)).GT.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
     		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0

			S(2,FI(N),FJ(N))=(S(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			U(2,FI(N),FJ(N))=(U(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			V(2,FI(N),FJ(N))=(V(2,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

		EndIf
		
	ElseIf(NX(FI(N),FJ(N)).GE.0.0.AND.NY(FI(N),FJ(N)).LE.0.0)Then
	
		If(F(FI(N),FJ(N)).LE.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0

			S(1,FI(N),FJ(N))=(S(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+S(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			U(1,FI(N),FJ(N))=(U(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+U(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			V(1,FI(N),FJ(N))=(V(1,FI(N)+1,FJ(N))*(0.5d0*PAI-SITA)+V(1,FI(N),FJ(N)-1)*SITA)/(0.5d0*PAI)
			
		ElseIf(F(FI(N),FJ(N)).GT.0.0)Then
		SITA=dACOS(ABS(NX(FI(N),FJ(N))))
		If(ABS(NX(FI(N),FJ(N))).GE.1.0)SITA=0.0d0

			S(2,FI(N),FJ(N))=(S(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+S(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			U(2,FI(N),FJ(N))=(U(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+U(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)
			V(2,FI(N),FJ(N))=(V(2,FI(N)-1,FJ(N))*(0.5d0*PAI-SITA)+V(2,FI(N),FJ(N)+1)*SITA)/(0.5d0*PAI)

		EndIf
		
	EndIf
	
	End Do

	EndIf
	
	Do I=IS+1,IE-1

		If(F(I,JS+1).GE.0.0d0)Then

!			U(1,I,JS+1)=U(1,I-1,JS+1)
!			V(1,I,JS+1)=V(1,I-1,JS+1)

		End If

	End Do
	    
	Do M=1,1
    	Do J=JS-3,JE+3
	Do I=IS-3,IE+3

        	Q(M,1,1,I,J)=((P(M,I,J)+PA(M))/S(M,I,J))**(1.0/GAM(M))
        	Q(M,1,2,I,J)=Q(M,1,1,I,J)*U(M,I,J)
        	Q(M,1,3,I,J)=Q(M,1,1,I,J)*V(M,I,J)
        	Q(M,1,4,I,J)=(P(M,I,J)+GAM(M)*PA(M))/(GAM(M)-1.0d0)+0.5d0*Q(M,1,1,I,J)*(U(M,I,J)*U(M,I,J)+V(M,I,J)*V(M,I,J))

	End Do
	End Do
	End Do
	
	End Subroutine Gf
	
	
	Subroutine Lxx(C,Q,IS,IE,JS,JE,LX,DX,GAM,PA)

	Implicit None
	
	Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    Integer :: N
    Integer :: L
    Integer :: II
    Integer :: III
    Double Precision :: DX
    Double Precision :: DY
    Double Precision :: Q2
    Double Precision :: A2
    Double Precision :: A3
    Double Precision :: B1
    Double Precision :: B2
    Double Precision :: B3
    Double Precision :: C1
    Double Precision :: C2
    Double Precision :: C3
    Double Precision :: PM
    Double Precision :: PM2
    Double Precision :: HH
    Double Precision :: HP
    Double Precision :: HM
    
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Integer,Parameter :: MaxMAT=4
	Double Precision,Dimension(:) :: GAM
	Double Precision,Dimension(:) :: PA
	Double Precision,Dimension(MaxMat) :: VL
	Double Precision,Dimension(MaxMat) :: VR
	Double Precision,Dimension(MaxMat) :: VV
	Double Precision,Dimension(MaxMat,MaxMat) :: LXMAT
	Double Precision,Dimension(MaxMat,MaxMat) :: RXMAT
    Double Precision,allocatable,Dimension(:,:,:) :: U
	Double Precision,allocatable,Dimension(:,:,:) :: V
	Double Precision,allocatable,Dimension(:,:,:) :: P
	Double Precision,allocatable,Dimension(:,:,:) :: S
	Double Precision,Dimension(:,:,:) :: C
	Double Precision,allocatable,Dimension(:,:,:) :: H
	Double Precision,allocatable,Dimension(:,:,:) :: LOAVE
	Double Precision,allocatable,Dimension(:,:,:) :: UAVE
	Double Precision,allocatable,Dimension(:,:,:) :: VAVE
	Double Precision,allocatable,Dimension(:,:,:) :: HAVE
	Double Precision,allocatable,Dimension(:,:,:) :: CAVE
	Double Precision,allocatable,Dimension(:,:,:,:) :: FLX
	Double Precision,allocatable,Dimension(:,:,:,:) :: WX
	Double Precision,allocatable,Dimension(:,:,:,:) :: FWX
	Double Precision,allocatable,Dimension(:,:,:,:) :: FLUX
	Double Precision,allocatable,Dimension(:,:,:,:) :: FLUXX
	Double Precision,Dimension(:,:,:,:) :: LX
	Double Precision,Dimension(:,:,:,:,:) :: Q
	
	Allocate(U(MaxM,IE+3,JE+3))
	Allocate(V(MaxM,IE+3,JE+3))
	Allocate(P(MaxM,IE+3,JE+3))
	Allocate(S(MaxM,IE+3,JE+3))
	Allocate(H(MaxM,IE+3,JE+3))
	Allocate(LOAVE(MaxM,IE+3,JE+3))
	Allocate(UAVE(MaxM,IE+3,JE+3))
	Allocate(VAVE(MaxM,IE+3,JE+3))
	Allocate(HAVE(MaxM,IE+3,JE+3))
	Allocate(CAVE(MaxM,IE+3,JE+3))
	Allocate(FLX(MaxM,MaxL,IE+3,JE+3))
	Allocate(WX(MaxM,MaxL,IE+3,JE+3))
	Allocate(FWX(MaxM,MaxL,IE+3,JE+3))
	Allocate(FLUX(MaxM,MaxL,IE+3,JE+3))
	Allocate(FLUXX(MaxM,MaxL,IE+3,JE+3))
	    
	Do M=1,1
	Do I=IS,IE-1
	Do J=JS,JE-1
	
	    U(M,I,J)=Q(M,1,2,I,J)/Q(M,1,1,I,J)
	    V(M,I,J)=Q(M,1,3,I,J)/Q(M,1,1,I,J)
	    P(M,I,J)=(GAM(M)-1.0)*(Q(M,1,4,I,J)-0.5d0*Q(M,1,1,I,J)*(U(M,I,J)*U(M,I,J)+V(M,I,J)*V(M,I,J)))-GAM(M)*PA(M)
	    H(M,I,J)=Q(M,1,4,I,J)/Q(M,1,1,I,J)+P(M,I,J)/Q(M,1,1,I,J)
	    U(M,I+1,J)=Q(M,1,2,I+1,J)/Q(M,1,1,I+1,J)
	    V(M,I+1,J)=Q(M,1,3,I+1,J)/Q(M,1,1,I+1,J)
	    P(M,I+1,J)=(GAM(M)-1.0)*(Q(M,1,4,I+1,J)-0.5d0*Q(M,1,1,I+1,J)*(U(M,I+1,J)*U(M,I+1,J)+V(M,I+1,J)*V(M,I+1,J)))-GAM(M)*PA(M)
	    H(M,I+1,J)=Q(M,1,4,I+1,J)/Q(M,1,1,I+1,J)+P(M,I+1,J)/Q(M,1,1,I+1,J)
	    LOAVE(M,I,J)=dSQRT(Q(M,1,1,I+1,J)*Q(M,1,1,I,J))
	    UAVE(M,I,J)=(dSQRT(Q(M,1,1,I+1,J))*U(M,I+1,J)+dSQRT(Q(M,1,1,I,J))*U(M,I,J))/(dSQRT(Q(M,1,1,I+1,J))+dSQRT(Q(M,1,1,I,J)))
	    VAVE(M,I,J)=(dSQRT(Q(M,1,1,I+1,J))*V(M,I+1,J)+dSQRT(Q(M,1,1,I,J))*V(M,I,J))/(dSQRT(Q(M,1,1,I+1,J))+dSQRT(Q(M,1,1,I,J)))
	    HAVE(M,I,J)=(dSQRT(Q(M,1,1,I+1,J))*H(M,I+1,J)+dSQRT(Q(M,1,1,I,J))*H(M,I,J))/(dSQRT(Q(M,1,1,I+1,J))+dSQRT(Q(M,1,1,I,J)))
	    CAVE(M,I,J)=dSQRT((GAM(M)-1.0)*(HAVE(M,I,J)-0.5d0*(UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J))))
	
        Q2=UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J)
		B1=(GAM(M)-1.0d0)/(CAVE(M,I,J)*CAVE(M,I,J))
		B2=0.5d0*Q2*B1
		LXMAT(1,1)=0.5d0*(B2+UAVE(M,I,J)/CAVE(M,I,J))
		LXMAT(1,2)=-0.5d0*(1.0d0/CAVE(M,I,J)+B1*UAVE(M,I,J))
		LXMAT(1,3)=-0.5d0*B1*VAVE(M,I,J)
		LXMAT(1,4)=0.5d0*B1
		LXMAT(2,1)=1.0d0-B2
		LXMAT(2,2)=B1*UAVE(M,I,J)
		LXMAT(2,3)=B1*VAVE(M,I,J)
		LXMAT(2,4)=-B1
		LXMAT(3,1)=0.5d0*(B2-UAVE(M,I,J)/CAVE(M,I,J))
		LXMAT(3,2)=0.5d0*(1.0d0/CAVE(M,I,J)-B1*UAVE(M,I,J))
		LXMAT(3,3)=-0.5d0*B1*VAVE(M,I,J)
		LXMAT(3,4)=0.5d0*B1
		LXMAT(4,1)=-VAVE(M,I,J)
		LXMAT(4,2)=0.0d0
		LXMAT(4,3)=1.0d0
		LXMAT(4,4)=0.0d0
	
	Do N=1,6

	    U(M,I+N-3,J)=Q(M,1,2,I+N-3,J)/Q(M,1,1,I+N-3,J)
	    V(M,I+N-3,J)=Q(M,1,3,I+N-3,J)/Q(M,1,1,I+N-3,J)
	    P(M,I+N-3,J)=(GAM(M)-1.0d0)*(Q(M,1,4,I+N-3,J)-0.5d0*Q(M,1,1,I+N-3,J)*(U(M,I+N-3,J)*U(M,I+N-3,J)+V(M,I+N-3,J)*V(M,I+N-3,J)))-GAM(M)*PA(M)
	    FLX(M,1,I+N-3,J)=Q(M,1,2,I+N-3,J)
	    FLX(M,2,I+N-3,J)=P(M,I+N-3,J)+Q(M,1,2,I+N-3,J)*U(M,I+N-3,J)
        FLX(M,3,I+N-3,J)=Q(M,1,1,I+N-3,J)*U(M,I+N-3,J)*V(M,I+N-3,J)
        FLX(M,4,I+N-3,J)=(P(M,I+N-3,J)+Q(M,1,4,I+N-3,J))*U(M,I+N-3,J)
        
        Do L=1,4
        
        WX(M,L,I+N-3,J)=LXMAT(L,1)*Q(M,1,1,I+N-3,J)+LXMAT(L,2)*Q(M,1,2,I+N-3,J)+LXMAT(L,3)*Q(M,1,3,I+N-3,J)+LXMAT(L,4)*Q(M,1,4,I+N-3,J)
        FWX(M,L,I+N-3,J)=LXMAT(L,1)*FLX(M,1,I+N-3,J)+LXMAT(L,2)*FLX(M,2,I+N-3,J)+LXMAT(L,3)*FLX(M,3,I+N-3,J)+LXMAT(L,4)*FLX(M,4,I+N-3,J)
        
        End Do   
    
    End Do
	
	    VL(1)=U(M,I,J)-C(M,I,J)
	    VL(2)=U(M,I,J)
	    VL(3)=U(M,I,J)+C(M,I,J)
	    VL(4)=U(M,I,J)
	    VR(1)=U(M,I+1,J)-C(M,I+1,J)
	    VR(2)=U(M,I+1,J)
	    VR(3)=U(M,I+1,J)+C(M,I+1,J)
	    VR(4)=U(M,I+1,J)
	
	Do L=1,4
	
	    VV(L)=MAX(ABS(VR(L)),ABS(VL(L)))
	
	End Do
	
	Do L=1,4
	
	    PM=1.0
	    III=0
	    C1=0.5*(FWX(M,L,I+III,J)+PM*VV(L)*WX(M,L,I+III,J))
        A2=0.25*((FWX(M,L,I+III+1,J)-FWX(M,L,I+III,J))+PM*VV(L)*(WX(M,L,I+III+1,J)-WX(M,L,I+III,J)))
        B2=0.25*((FWX(M,L,I+III,J)-FWX(M,L,I+III-1,J))+PM*VV(L)*(WX(M,L,I+III,J)-WX(M,L,I+III-1,J)))
        
        If(ABS(A2).GE.ABS(B2))Then
            C2=B2
            II=-1
		PM2=0.5d0+PM*1.5d0
        Else
            C2=A2
            II=0
		PM2=0.5d0-PM*1.5d0
        EndIf
        
        A3=((FWX(M,L,I+III+2+II,J)-2.0d0*FWX(M,L,I+III+1+II,J)+FWX(M,L,I+III+II,J))+PM*VV(L)*(WX(M,L,I+III+2+II,J)-2.0d0*WX(M,L,I+III+1+II,J)+WX(M,L,I+III+II,J)))/12.0d0
        B3=((FWX(M,L,I+III+1+II,J)-2.0d0*FWX(M,L,I+III+II,J)+FWX(M,L,I+III-1+II,J))+PM*VV(L)*(WX(M,L,I+III+1+II,J)-2.0d0*WX(M,L,I+III+II,J)+WX(M,L,I+III-1+II,J)))/12.0d0
    	
	    If(ABS(A3).GE.ABS(B3))Then
	        C3=B3
!	        PM2=0.5d0+PM*1.5d0
	    Else
	        C3=A3
!	        PM2=0.5d0-PM*1.5d0
	    EndIf
    	
	    HH=C1+PM*C2+PM2*C3
	    HP=HH
	
	    PM=-1.0d0
	    III=1
	    C1=0.5d0*(FWX(M,L,I+III,J)+PM*VV(L)*WX(M,L,I+III,J))
        A2=0.25d0*((FWX(M,L,I+III+1,J)-FWX(M,L,I+III,J))+PM*VV(L)*(WX(M,L,I+III+1,J)-WX(M,L,I+III,J)))
        B2=0.25d0*((FWX(M,L,I+III,J)-FWX(M,L,I+III-1,J))+PM*VV(L)*(WX(M,L,I+III,J)-WX(M,L,I+III-1,J)))
        
        If(ABS(A2).GE.ABS(B2))Then
            C2=B2
            II=-1
	PM2=0.5d0+PM*1.5d0
        Else
            C2=A2
            II=0
	PM2=0.5d0-PM*1.5d0
        EndIf
        
        A3=((FWX(M,L,I+III+2+II,J)-2.0d0*FWX(M,L,I+III+1+II,J)+FWX(M,L,I+III+II,J))+PM*VV(L)*(WX(M,L,I+III+2+II,J)-2.0d0*WX(M,L,I+III+1+II,J)+WX(M,L,I+III+II,J)))/12.0d0
        B3=((FWX(M,L,I+III+1+II,J)-2.0d0*FWX(M,L,I+III+II,J)+FWX(M,L,I+III-1+II,J))+PM*VV(L)*(WX(M,L,I+III+1+II,J)-2.0d0*WX(M,L,I+III+II,J)+WX(M,L,I+III-1+II,J)))/12.0d0
    	
	    If(ABS(A3).GE.ABS(B3))Then
	        C3=B3
!	        PM2=0.5d0+PM*1.5d0
	    Else
	        C3=A3
!	        PM2=0.5d0-PM*1.5d0
	    EndIf
    	
	    HH=C1+PM*C2+PM2*C3
	    HM=HH
	
	    FLUXX(M,L,I,J)=HP+HM
	
	End Do
	
!	    U(M,I,J)=Q(M,1,2,I,J)/Q(M,1,1,I,J)
!	    V(M,I,J)=Q(M,1,3,I,J)/Q(M,1,1,I,J)
!	    P(M,I,J)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J)-0.5d0*Q(M,1,1,I,J)*(U(M,I,J)*U(M,I,J)+V(M,I,J)*V(M,I,J)))-GAM(M)*PA(M)
!	    H(M,I,J)=Q(M,1,4,I,J)/Q(M,1,1,I,J)+P(M,I,J)/Q(M,1,1,I,J)
!	    U(M,I+1,J)=Q(M,1,2,I+1,J)/Q(M,1,1,I+1,J)
!	    V(M,I+1,J)=Q(M,1,3,I+1,J)/Q(M,1,1,I+1,J)
!	    P(M,I+1,J)=(GAM(M)-1.0d0)*(Q(M,1,4,I+1,J)-0.5d0*Q(M,1,1,I+1,J)*(U(M,I+1,J)*U(M,I+1,J)+V(M,I+1,J)*V(M,I+1,J)))-GAM(M)*PA(M)
!	    H(M,I+1,J)=Q(M,1,4,I+1,J)/Q(M,1,1,I+1,J)+P(M,I+1,J)/Q(M,1,1,I+1,J)
!	    LOAVE(M,I,J)=dSQRT(Q(M,1,1,I+1,J)*Q(M,1,1,I,J))
!	    UAVE(M,I,J)=(dSQRT(Q(M,1,1,I+1,J))*U(M,I+1,J)+dSQRT(Q(M,1,1,I,J))*U(M,I,J))/(dSQRT(Q(M,1,1,I+1,J))+dSQRT(Q(M,1,1,I,J)))
!	    VAVE(M,I,J)=(dSQRT(Q(M,1,1,I+1,J))*V(M,I+1,J)+dSQRT(Q(M,1,1,I,J))*V(M,I,J))/(dSQRT(Q(M,1,1,I+1,J))+dSQRT(Q(M,1,1,I,J)))
!	    HAVE(M,I,J)=(dSQRT(Q(M,1,1,I+1,J))*H(M,I+1,J)+dSQRT(Q(M,1,1,I,J))*H(M,I,J))/(dSQRT(Q(M,1,1,I+1,J))+dSQRT(Q(M,1,1,I,J)))
!	    CAVE(M,I,J)=dSQRT((GAM(M)-1.0d0)*(HAVE(M,I,J)-0.5d0*(UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J))))
	
	    Q2=UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J)
		RXMAT(1,1)=1.0d0
		RXMAT(1,2)=1.0d0
		RXMAT(1,3)=1.0d0
		RXMAT(1,4)=0.0d0
		RXMAT(2,1)=UAVE(M,I,J)-CAVE(M,I,J)
		RXMAT(2,2)=UAVE(M,I,J)
		RXMAT(2,3)=UAVE(M,I,J)+CAVE(M,I,J)
		RXMAT(2,4)=0.0d0
		RXMAT(3,1)=VAVE(M,I,J)
		RXMAT(3,2)=VAVE(M,I,J)
		RXMAT(3,3)=VAVE(M,I,J)
		RXMAT(3,4)=1.0d0
		RXMAT(4,1)=HAVE(M,I,J)-CAVE(M,I,J)*UAVE(M,I,J)
		RXMAT(4,2)=0.5d0*Q2
		RXMAT(4,3)=HAVE(M,I,J)+CAVE(M,I,J)*UAVE(M,I,J)
		RXMAT(4,4)=VAVE(M,I,J)
	
	Do L=1,4
    
	    FLUX(M,L,I,J)=RXMAT(L,1)*FLUXX(M,1,I,J)+RXMAT(L,2)*FLUXX(M,2,I,J)+RXMAT(L,3)*FLUXX(M,3,I,J)+RXMAT(L,4)*FLUXX(M,4,I,J)
	    
	End Do
	
	End Do
	End Do
	End Do
	
    Do M=1,1
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	Do L=1,4
	
	    LX(M,L,I,J)=-(FLUX(M,L,I,J)-FLUX(M,L,I-1,J))/DX
	    
	End Do
	End Do
	End Do
	End Do

	End Subroutine Lxx
	
	
	Subroutine Lyy(C,Q,IS,IE,JS,JE,LY,DY,GAM,PA)

	Implicit None
	
	Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    Integer :: N
    Integer :: L
    Integer :: JJ
    Integer :: JJJ
    Double Precision :: DX
    Double Precision :: DY
    Double Precision :: Q2
    Double Precision :: A2
    Double Precision :: A3
    Double Precision :: B1
    Double Precision :: B2
    Double Precision :: B3
    Double Precision :: C1
    Double Precision :: C2
    Double Precision :: C3
    Double Precision :: PM
    Double Precision :: PM2
    Double Precision :: HH
    Double Precision :: HP
    Double Precision :: HM
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Integer,Parameter :: MaxMAT=4
	Double Precision,Dimension(:) :: GAM
	Double Precision,Dimension(:) :: PA
	Double Precision,Dimension(MaxMat) :: VL
	Double Precision,Dimension(MaxMat) :: VR
	Double Precision,Dimension(MaxMat) :: VV
	Double Precision,Dimension(MaxMat,MaxMat) :: LYMAT
	Double Precision,Dimension(MaxMat,MaxMat) :: RYMAT
    Double Precision,allocatable,Dimension(:,:,:) :: U
	Double Precision,allocatable,Dimension(:,:,:) :: V
	Double Precision,allocatable,Dimension(:,:,:) :: P
	Double Precision,allocatable,Dimension(:,:,:) :: S
	Double Precision,Dimension(:,:,:) :: C
	Double Precision,allocatable,Dimension(:,:,:) :: H
	Double Precision,allocatable,Dimension(:,:,:) :: LOAVE
	Double Precision,allocatable,Dimension(:,:,:) :: UAVE
	Double Precision,allocatable,Dimension(:,:,:) :: VAVE
	Double Precision,allocatable,Dimension(:,:,:) :: HAVE
	Double Precision,allocatable,Dimension(:,:,:) :: CAVE
	Double Precision,allocatable,Dimension(:,:,:,:) :: FLY
	Double Precision,allocatable,Dimension(:,:,:,:) :: WY
	Double Precision,allocatable,Dimension(:,:,:,:) :: FWY
	Double Precision,allocatable,Dimension(:,:,:,:) :: FLUX
	Double Precision,allocatable,Dimension(:,:,:,:) :: FLUXY
	Double Precision,Dimension(:,:,:,:) :: LY
	Double Precision,Dimension(:,:,:,:,:) :: Q
	
	Allocate(U(MaxM,IE+3,JE+3))
	Allocate(V(MaxM,IE+3,JE+3))
	Allocate(P(MaxM,IE+3,JE+3))
	Allocate(S(MaxM,IE+3,JE+3))
	Allocate(H(MaxM,IE+3,JE+3))
	Allocate(LOAVE(MaxM,IE+3,JE+3))
	Allocate(UAVE(MaxM,IE+3,JE+3))
	Allocate(VAVE(MaxM,IE+3,JE+3))
	Allocate(HAVE(MaxM,IE+3,JE+3))
	Allocate(CAVE(MaxM,IE+3,JE+3))
	Allocate(FLY(MaxM,MaxL,IE+3,JE+3))
	Allocate(WY(MaxM,MaxL,IE+3,JE+3))
	Allocate(FWY(MaxM,MaxL,IE+3,JE+3))
	Allocate(FLUX(MaxM,MaxL,IE+3,JE+3))
	Allocate(FLUXY(MaxM,MaxL,IE+3,JE+3))
	
	Do M=1,1
	Do I=IS,IE-1
	Do J=JS,JE-1
	
	    U(M,I,J)=Q(M,1,2,I,J)/Q(M,1,1,I,J)
	    V(M,I,J)=Q(M,1,3,I,J)/Q(M,1,1,I,J)
	    P(M,I,J)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J)-0.5d0*Q(M,1,1,I,J)*(U(M,I,J)*U(M,I,J)+V(M,I,J)*V(M,I,J)))-GAM(M)*PA(M)
	    H(M,I,J)=Q(M,1,4,I,J)/Q(M,1,1,I,J)+P(M,I,J)/Q(M,1,1,I,J)
	    U(M,I,J+1)=Q(M,1,2,I,J+1)/Q(M,1,1,I,J+1)
	    V(M,I,J+1)=Q(M,1,3,I,J+1)/Q(M,1,1,I,J+1)
	    P(M,I,J+1)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J+1)-0.5d0*Q(M,1,1,I,J+1)*(U(M,I,J+1)*U(M,I,J+1)+V(M,I,J+1)*V(M,I,J+1)))-GAM(M)*PA(M)
	    H(M,I,J+1)=Q(M,1,4,I,J+1)/Q(M,1,1,I,J+1)+P(M,I,J+1)/Q(M,1,1,I,J+1)
	    LOAVE(M,I,J)=dSQRT(Q(M,1,1,I,J+1)*Q(M,1,1,I,J))
	    UAVE(M,I,J)=(dSQRT(Q(M,1,1,I,J+1))*U(M,I,J+1)+dSQRT(Q(M,1,1,I,J))*U(M,I,J))/(dSQRT(Q(M,1,1,I,J+1))+dSQRT(Q(M,1,1,I,J)))
	    VAVE(M,I,J)=(dSQRT(Q(M,1,1,I,J+1))*V(M,I,J+1)+dSQRT(Q(M,1,1,I,J))*V(M,I,J))/(dSQRT(Q(M,1,1,I,J+1))+dSQRT(Q(M,1,1,I,J)))
	    HAVE(M,I,J)=(dSQRT(Q(M,1,1,I,J+1))*H(M,I,J+1)+dSQRT(Q(M,1,1,I,J))*H(M,I,J))/(dSQRT(Q(M,1,1,I,J+1))+dSQRT(Q(M,1,1,I,J)))
	    CAVE(M,I,J)=dSQRT((GAM(M)-1.0d0)*(HAVE(M,I,J)-0.5d0*(UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J))))

	
        Q2=UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J)
		B1=(GAM(M)-1.0d0)/(CAVE(M,I,J)*CAVE(M,I,J))
		B2=0.5d0*Q2*B1
		LYMAT(1,1)=0.5d0*(B2+VAVE(M,I,J)/CAVE(M,I,J))
		LYMAT(1,2)=-0.5d0*B1*UAVE(M,I,J)
		LYMAT(1,3)=-0.5d0*(1.0d0/CAVE(M,I,J)+B1*VAVE(M,I,J))
		LYMAT(1,4)=0.5d0*B1
		LYMAT(2,1)=1.0d0-B2
		LYMAT(2,2)=B1*UAVE(M,I,J)
		LYMAT(2,3)=B1*VAVE(M,I,J)
		LYMAT(2,4)=-B1
		LYMAT(3,1)=0.5d0*(B2-VAVE(M,I,J)/CAVE(M,I,J))
		LYMAT(3,2)=-0.5d0*B1*UAVE(M,I,J)
		LYMAT(3,3)=0.5d0*(1.0d0/CAVE(M,I,J)-B1*VAVE(M,I,J))
		LYMAT(3,4)=0.5d0*B1
		LYMAT(4,1)=-UAVE(M,I,J)
		LYMAT(4,2)=1.0d0
		LYMAT(4,3)=0.0d0
		LYMAT(4,4)=0.0d0
    
    Do N=1,6

	    U(M,I,J+N-3)=Q(M,1,2,I,J+N-3)/Q(M,1,1,I,J+N-3)
	    V(M,I,J+N-3)=Q(M,1,3,I,J+N-3)/Q(M,1,1,I,J+N-3)
	    P(M,I,J+N-3)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J+N-3)-0.5d0*Q(M,1,1,I,J+N-3)*(U(M,I,J+N-3)*U(M,I,J+N-3)+V(M,I,J+N-3)*V(M,I,J+N-3)))-GAM(M)*PA(M)
	    FLY(M,1,I,J+N-3)=Q(M,1,3,I,J+N-3)
	    FLY(M,2,I,J+N-3)=Q(M,1,1,I,J+N-3)*U(M,I,J+N-3)*V(M,I,J+N-3)
        FLY(M,3,I,J+N-3)=P(M,I,J+N-3)+Q(M,1,3,I,J+N-3)*V(M,I,J+N-3)
        FLY(M,4,I,J+N-3)=(P(M,I,J+N-3)+Q(M,1,4,I,J+N-3))*V(M,I,J+N-3)
        
        Do L=1,4
        
        WY(M,L,I,J+N-3)=LYMAT(L,1)*Q(M,1,1,I,J+N-3)+LYMAT(L,2)*Q(M,1,2,I,J+N-3)+LYMAT(L,3)*Q(M,1,3,I,J+N-3)+LYMAT(L,4)*Q(M,1,4,I,J+N-3)
        FWY(M,L,I,J+N-3)=LYMAT(L,1)*FLY(M,1,I,J+N-3)+LYMAT(L,2)*FLY(M,2,I,J+N-3)+LYMAT(L,3)*FLY(M,3,I,J+N-3)+LYMAT(L,4)*FLY(M,4,I,J+N-3)
        
        End Do   
    
    End Do
	
	    VL(1)=V(M,I,J)-C(M,I,J)
	    VL(2)=V(M,I,J)
	    VL(3)=V(M,I,J)+C(M,I,J)
	    VL(4)=V(M,I,J)
	    VR(1)=V(M,I,J+1)-C(M,I,J+1)
	    VR(2)=V(M,I,J+1)
	    VR(3)=V(M,I,J+1)+C(M,I,J+1)
	    VR(4)=V(M,I,J+1)
	
	Do L=1,4
	
	    VV(L)=MAX(ABS(VR(L)),ABS(VL(L)))
	
	End Do
	
	Do L=1,4
	
	    PM=1.0d0
	    JJJ=0
	    C1=0.5*(FWY(M,L,I,J+JJJ)+PM*VV(L)*WY(M,L,I,J+JJJ))
        A2=0.25*((FWY(M,L,I,J+JJJ+1)-FWY(M,L,I,J+JJJ))+PM*VV(L)*(WY(M,L,I,J+JJJ+1)-WY(M,L,I,J+JJJ)))
        B2=0.25*((FWY(M,L,I,J+JJJ)-FWY(M,L,I,J+JJJ-1))+PM*VV(L)*(WY(M,L,I,J+JJJ)-WY(M,L,I,J+JJJ-1)))
        
        If(ABS(A2).GE.ABS(B2))Then
            C2=B2
            JJ=-1
		PM2=0.5d0+PM*1.5d0
        Else
            C2=A2
            JJ=0
		PM2=0.5d0-PM*1.5d0
        EndIf
        
        A3=((FWY(M,L,I,J+JJJ+2+JJ)-2.0d0*FWY(M,L,I,J+JJJ+1+JJ)+FWY(M,L,I,J+JJJ+JJ))+PM*VV(L)*(WY(M,L,I,J+JJJ+2+JJ)-2.0d0*WY(M,L,I,J+JJJ+1+JJ)+WY(M,L,I,J+JJJ+JJ)))/12.0d0
        B3=((FWY(M,L,I,J+JJJ+1+JJ)-2.0d0*FWY(M,L,I,J+JJJ+JJ)+FWY(M,L,I,J+JJJ-1+JJ))+PM*VV(L)*(WY(M,L,I,J+JJJ+1+JJ)-2.0d0*WY(M,L,I,J+JJJ+JJ)+WY(M,L,I,J+JJJ-1+JJ)))/12.0d0
    	
	    If(ABS(A3).GE.ABS(B3))Then
	        C3=B3
!	        PM2=0.5d0+PM*1.5d0
	    Else
	        C3=A3
!	        PM2=0.5d0-PM*1.5d0
	    EndIf
    	
	    HH=C1+PM*C2+PM2*C3
	    HP=HH
	
	    PM=-1.0d0
	    JJJ=1
	    C1=0.5d0*(FWY(M,L,I,J+JJJ)+PM*VV(L)*WY(M,L,I,J+JJJ))
        A2=0.25d0*((FWY(M,L,I,J+JJJ+1)-FWY(M,L,I,J+JJJ))+PM*VV(L)*(WY(M,L,I,J+JJJ+1)-WY(M,L,I,J+JJJ)))
        B2=0.25d0*((FWY(M,L,I,J+JJJ)-FWY(M,L,I,J+JJJ-1))+PM*VV(L)*(WY(M,L,I,J+JJJ)-WY(M,L,I,J+JJJ-1)))
        
        If(ABS(A2).GE.ABS(B2))Then
            C2=B2
            JJ=-1
		PM2=0.5d0+PM*1.5d0
        Else
            C2=A2
            JJ=0
		PM2=0.5d0-PM*1.5d0
        EndIf
        
        A3=((FWY(M,L,I,J+JJJ+2+JJ)-2.0d0*FWY(M,L,I,J+JJJ+1+JJ)+FWY(M,L,I,J+JJJ+JJ))+PM*VV(L)*(WY(M,L,I,J+JJJ+2+JJ)-2.0d0*WY(M,L,I,J+JJJ+1+JJ)+WY(M,L,I,J+JJJ+JJ)))/12.0d0
        B3=((FWY(M,L,I,J+JJJ+1+JJ)-2.0d0*FWY(M,L,I,J+JJJ+JJ)+FWY(M,L,I,J+JJJ-1+JJ))+PM*VV(L)*(WY(M,L,I,J+JJJ+1+JJ)-2.0d0*WY(M,L,I,J+JJJ+JJ)+WY(M,L,I,J+JJJ-1+JJ)))/12.0d0
    	
	    If(ABS(A3).GE.ABS(B3))Then
	        C3=B3
!	        PM2=0.5d0+PM*1.5d0
	    Else
	        C3=A3
!	        PM2=0.5d0-PM*1.5d0
	    EndIf
    	
	    HH=C1+PM*C2+PM2*C3
	    HM=HH
	
	    FLUXY(M,L,I,J)=HP+HM
	
	End Do
	
!	    U(M,I,J)=Q(M,1,2,I,J)/Q(M,1,1,I,J)
!	    V(M,I,J)=Q(M,1,3,I,J)/Q(M,1,1,I,J)
!	    P(M,I,J)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J)-0.5d0*Q(M,1,1,I,J)*(U(M,I,J)*U(M,I,J)+V(M,I,J)*V(M,I,J)))-GAM(M)*PA(M)
!	    H(M,I,J)=Q(M,1,4,I,J)/Q(M,1,1,I,J)+P(M,I,J)/Q(M,1,1,I,J)
!	    U(M,I,J+1)=Q(M,1,2,I,J+1)/Q(M,1,1,I,J+1)
!	    V(M,I,J+1)=Q(M,1,3,I,J+1)/Q(M,1,1,I,J+1)
!	    P(M,I,J+1)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J+1)-0.5d0*Q(M,1,1,I,J+1)*(U(M,I,J+1)*U(M,I,J+1)+V(M,I,J+1)*V(M,I,J+1)))-GAM(M)*PA(M)
!	    H(M,I,J+1)=Q(M,1,4,I,J+1)/Q(M,1,1,I,J+1)+P(M,I,J+1)/Q(M,1,1,I,J+1)
!	    LOAVE(M,I,J)=dSQRT(Q(M,1,1,I,J+1)*Q(M,1,1,I,J))
!	    UAVE(M,I,J)=(dSQRT(Q(M,1,1,I,J+1))*U(M,I,J+1)+dSQRT(Q(M,1,1,I,J))*U(M,I,J))/(dSQRT(Q(M,1,1,I,J+1))+dSQRT(Q(M,1,1,I,J)))
!	    VAVE(M,I,J)=(dSQRT(Q(M,1,1,I,J+1))*V(M,I,J+1)+dSQRT(Q(M,1,1,I,J))*V(M,I,J))/(dSQRT(Q(M,1,1,I,J+1))+dSQRT(Q(M,1,1,I,J)))
!	    HAVE(M,I,J)=(dSQRT(Q(M,1,1,I,J+1))*H(M,I,J+1)+dSQRT(Q(M,1,1,I,J))*H(M,I,J))/(dSQRT(Q(M,1,1,I,J+1))+dSQRT(Q(M,1,1,I,J)))
!	    CAVE(M,I,J)=dSQRT((GAM(M)-1.0d0)*(HAVE(M,I,J)-0.5d0*(UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J))))
	
	    Q2=UAVE(M,I,J)*UAVE(M,I,J)+VAVE(M,I,J)*VAVE(M,I,J)
		RYMAT(1,1)=1.0d0
		RYMAT(1,2)=1.0d0
		RYMAT(1,3)=1.0d0
		RYMAT(1,4)=0.0d0
		RYMAT(2,1)=UAVE(M,I,J)
		RYMAT(2,2)=UAVE(M,I,J)
		RYMAT(2,3)=UAVE(M,I,J)
		RYMAT(2,4)=1.0d0
		RYMAT(3,1)=VAVE(M,I,J)-CAVE(M,I,J)
		RYMAT(3,2)=VAVE(M,I,J)
		RYMAT(3,3)=VAVE(M,I,J)+CAVE(M,I,J)
		RYMAT(3,4)=0.0d0
		RYMAT(4,1)=HAVE(M,I,J)-CAVE(M,I,J)*VAVE(M,I,J)
		RYMAT(4,2)=0.5d0*Q2
		RYMAT(4,3)=HAVE(M,I,J)+CAVE(M,I,J)*VAVE(M,I,J)
		RYMAT(4,4)=UAVE(M,I,J)
	
	Do L=1,4 
    
	    FLUX(M,L,I,J)=RYMAT(L,1)*FLUXY(M,1,I,J)+RYMAT(L,2)*FLUXY(M,2,I,J)+RYMAT(L,3)*FLUXY(M,3,I,J)+RYMAT(L,4)*FLUXY(M,4,I,J)
	    
	End Do
	
	End Do
	End Do
	End Do
	
	Do M=1,1
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	Do L=1,4
	
	    LY(M,L,I,J)=-(FLUX(M,L,I,J)-FLUX(M,L,I,J-1))/DY
	    
	End Do
	End Do
	End Do
	End Do

	End Subroutine Lyy
	
	
	Subroutine Lss(Q,U,V,P,RR,LS,IS,IE,JS,JE)

	Implicit None
	
	Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Integer,Parameter :: MaxMAT=4
	Double Precision,Dimension(:) :: RR
    Double Precision,Dimension(:,:,:) :: U
	Double Precision,Dimension(:,:,:) :: V
	Double Precision,Dimension(:,:,:) :: P
	Double Precision,Dimension(:,:,:,:) :: LS
	Double Precision,Dimension(:,:,:,:,:) :: Q
	
	Do M=1,1
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
	    LS(M,1,I,J)=-Q(M,1,2,I,J)/RR(I)
	    LS(M,2,I,J)=-Q(M,1,2,I,J)*U(M,I,J)/RR(I)
	    LS(M,3,I,J)=-Q(M,1,2,I,J)*V(M,I,J)/RR(I)
	    LS(M,4,I,J)=-(Q(M,1,4,I,J)+P(M,I,J))*U(M,I,J)/RR(I)
	    
	End Do
	End Do
	End Do

	End Subroutine Lss
	
	
	Subroutine Boundary(Q,IS,IE,JS,JE,F,DX)

	Implicit None
	
	Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    Integer :: L
	Double Precision :: DX
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
    Double Precision,Dimension(:,:,:,:,:) :: Q
	Double Precision,Dimension(:,:) :: F
	
	Do M=1,1
	Do L=1,4
	
	Do I=IS+1,IE-1
	    
	    Q(M,1,L,I,JE)=Q(M,1,L,I,JE-1)
	    Q(M,1,L,I,JE+1)=Q(M,1,L,I,JE-1)
	    Q(M,1,L,I,JE+2)=Q(M,1,L,I,JE-1)
	    Q(M,1,L,I,JE+3)=Q(M,1,L,I,JE-1)
!		Q(M,1,1,I,JE)=0.001d0
!	    Q(M,1,1,I,JE+1)=0.001d0
!	    Q(M,1,1,I,JE+2)=0.001d0
!	    Q(M,1,1,I,JE+3)=0.001d0
!		Q(M,1,3,I,JE)=-0.001d0
!	    Q(M,1,3,I,JE+1)=-0.001d0
!	    Q(M,1,3,I,JE+2)=-0.001d0
!	    Q(M,1,3,I,JE+3)=-0.001d0
	    
!		If(F(I,JS+1).LE.-0.5*DX)Then
	    Q(M,1,L,I,JS)=Q(M,1,L,I,JS+1)
	    Q(M,1,L,I,JS-1)=Q(M,1,L,I,JS+2)
	    Q(M,1,L,I,JS-2)=Q(M,1,L,I,JS+3)
	    Q(M,1,L,I,JS-3)=Q(M,1,L,I,JS+4)
		Q(M,1,3,I,JS)=-Q(M,1,3,I,JS+1)
	    Q(M,1,3,I,JS-1)=-Q(M,1,3,I,JS+2)
	    Q(M,1,3,I,JS-2)=-Q(M,1,3,I,JS+3)
	    Q(M,1,3,I,JS-3)=-Q(M,1,3,I,JS+4)
!		Else
!		Q(M,1,L,I,JS)=Q(M,1,L,I,JS+1)
!	    Q(M,1,L,I,JS-1)=Q(M,1,L,I,JS+1)
!	    Q(M,1,L,I,JS-2)=Q(M,1,L,I,JS+1)
!	    Q(M,1,L,I,JS-3)=Q(M,1,L,I,JS+1)
!		EndIf
	    
	End Do
	
	Do J=JS-3,JE+3
	
!	    Q(M,1,L,IS,J)=Q(M,1,L,IS+1,J)
!	    Q(M,1,L,IS-1,J)=Q(M,1,L,IS+1,J) 
!	    Q(M,1,L,IS-2,J)=Q(M,1,L,IS+1,J) 
!	    Q(M,1,L,IS-3,J)=Q(M,1,L,IS+1,J)
	    Q(M,1,L,IE,J)=Q(M,1,L,IE-1,J)
	    Q(M,1,L,IE+1,J)=Q(M,1,L,IE-1,J)
	    Q(M,1,L,IE+2,J)=Q(M,1,L,IE-1,J)
	    Q(M,1,L,IE+3,J)=Q(M,1,L,IE-1,J)
	   
	    Q(M,1,L,IS,J)=Q(M,1,L,IS+1,J)
	    Q(M,1,L,IS-1,J)=Q(M,1,L,IS+2,J)
	    Q(M,1,L,IS-2,J)=Q(M,1,L,IS+3,J)
	    Q(M,1,L,IS-3,J)=Q(M,1,L,IS+4,J)
	    Q(M,1,2,IS,J)=-Q(M,1,2,IS+1,J)
	    Q(M,1,2,IS-1,J)=-Q(M,1,2,IS+2,J)
	    Q(M,1,2,IS-2,J)=-Q(M,1,2,IS+3,J)
	    Q(M,1,2,IS-3,J)=-Q(M,1,2,IS+4,J)
	    
	End Do
	
	End Do
	End Do
    
    End Subroutine Boundary


	Subroutine Boundary2(Q,IS,IE,JS,JE)

	Implicit None
	
	Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    Integer :: L
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
    Double Precision,Dimension(:,:,:,:,:) :: Q
	
	Do M=1,1
	Do L=1,4
	
	Do I=IS+1,IE-1
	    
	    Q(M,1,L,I,JS)=Q(M,1,L,I,JS+1)
	    Q(M,1,L,I,JS-1)=Q(M,1,L,I,JS+2)
	    Q(M,1,L,I,JS-2)=Q(M,1,L,I,JS+3)
	    Q(M,1,L,I,JS-3)=Q(M,1,L,I,JS+4)
	
	    Q(M,1,3,I,JS)=-Q(M,1,3,I,JS+1)
	    Q(M,1,3,I,JS-1)=-Q(M,1,3,I,JS+2)
	    Q(M,1,3,I,JS-2)=-Q(M,1,3,I,JS+3)
	    Q(M,1,3,I,JS-3)=-Q(M,1,3,I,JS+4)
	    
	End Do
	
	Do J=JS-3,JE+3
	
!	    Q(M,1,L,IS,J)=Q(M,1,L,IS+1,J)
!	    Q(M,1,L,IS-1,J)=Q(M,1,L,IS+1,J) 
!	    Q(M,1,L,IS-2,J)=Q(M,1,L,IS+1,J) 
!	    Q(M,1,L,IS-3,J)=Q(M,1,L,IS+1,J)
	   
	    Q(M,1,L,IS,J)=Q(M,1,L,IS+1,J)
	    Q(M,1,L,IS-1,J)=Q(M,1,L,IS+2,J)
	    Q(M,1,L,IS-2,J)=Q(M,1,L,IS+3,J)
	    Q(M,1,L,IS-3,J)=Q(M,1,L,IS+4,J)
	    Q(M,1,2,IS,J)=-Q(M,1,2,IS+1,J)
	    Q(M,1,2,IS-1,J)=-Q(M,1,2,IS+2,J)
	    Q(M,1,2,IS-2,J)=-Q(M,1,2,IS+3,J)
	    Q(M,1,2,IS-3,J)=-Q(M,1,2,IS+4,J)
	    
	End Do
	
	End Do
	End Do
    
    End Subroutine Boundary2
	

    Subroutine Eluer(DT,DX,DY,U,V,P,Q,F,IS,IE,JS,JE,GAM,PA,RR,LL)

	Implicit None

    Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
    Integer :: L
    Integer :: LL
	Integer :: NMAX
    Double Precision :: DT
    Double Precision :: DX
    Double Precision :: DY 
	
    
    Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(MaxM) :: GAM
	Double Precision,Dimension(MaxM) :: PA
	Double Precision,Dimension(:) :: RR
	Integer,Allocatable,Dimension(:) :: FI1
	Integer,Allocatable,Dimension(:) :: FJ1
	Integer,Allocatable,Dimension(:) :: FI0
	Integer,Allocatable,Dimension(:) :: FJ0
	Double Precision,Dimension(:,:) :: F
	Double Precision,Dimension(:,:,:) :: U
	Double Precision,Dimension(:,:,:) :: V
	Double Precision,Allocatable,Dimension(:,:,:) :: C
	Double Precision,Dimension(:,:,:) :: P
	Double Precision,Allocatable,Dimension(:,:,:) :: S
	Double Precision,Allocatable,Dimension(:,:,:,:) :: LX
	Double Precision,Allocatable,Dimension(:,:,:,:) :: LY
	Double Precision,Allocatable,Dimension(:,:,:,:) :: LS
	Double Precision,Dimension(:,:,:,:,:) :: Q
	
	Allocate(LX(MaxM,MaxL,IE+3,JE+3))
	Allocate(LY(MaxM,MaxL,IE+3,JE+3))
	Allocate(LS(MaxM,MaxL,IE+3,JE+3))
	Allocate(C(MaxM,IE+3,JE+3))
	Allocate(S(MaxM,IE+3,JE+3))
	Allocate(FI1(IE*JE))
	Allocate(FJ1(IE*JE))
	Allocate(FI0(IE*JE))
	Allocate(FJ0(IE*JE))
	
    Do M=1,1
	Do L=1,4
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
       Q(M,2,L,I,J)=Q(M,1,L,I,J)
   
	End Do
	End Do
	End Do
	End Do

	LX=0.0d0
	LY=0.0d0
	LS=0.0d0

	Call heapsort(IS,IE,JS,JE,DX,F,FI1,FJ1,NMAX,5.0d0,1)

!	Call heapsort(IS,IE,JS,JE,DX,F,FI0,FJ0,NMAX,5.0d0,0)

    If(LL.EQ.1) Then
		Call Boundary(Q,IS,IE,JS,JE,F,DX)
	Else
		Call Boundary2(Q,IS,IE,JS,JE)
	EndIf
    Call Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI1,FJ1,NMAX,0,1)
!	Call Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI0,FJ0,NMAX,0,0)

    If(LL.EQ.1) Then
		Call Boundary(Q,IS,IE,JS,JE,F,DX)
	Else
		Call Boundary2(Q,IS,IE,JS,JE)
	EndIf
    
    Do M=1,1
    Do J=JS-3,JE+3
    Do I=IS-3,IE+3
    
        C(M,I,J)=dSQRT(GAM(M)*(P(M,I,J)+PA(M))/Q(M,1,1,I,J))
    
    End Do
    End Do
    End Do
    
    Call Lxx(C,Q,IS,IE,JS,JE,LX,DX,GAM,PA)
    Call Lyy(C,Q,IS,IE,JS,JE,LY,DY,GAM,PA)
    Call Lss(Q,U,V,P,RR,LS,IS,IE,JS,JE) 
	    
    Do M=1,1
	Do L=1,4
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
	    Q(M,1,L,I,J)=Q(M,2,L,I,J)+DT*(LX(M,L,I,J)+LY(M,L,I,J)+LS(M,L,I,J))
		    
	End Do
	End Do
	End Do
	End Do

	If(LL.EQ.1) Then
		Call Boundary(Q,IS,IE,JS,JE,F,DX)
	Else
		Call Boundary2(Q,IS,IE,JS,JE)
	EndIf
    Call Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI1,FJ1,NMAX,0,1)
!	Call Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI0,FJ0,NMAX,0,0)
    If(LL.EQ.1) Then
		Call Boundary(Q,IS,IE,JS,JE,F,DX)
	Else
		Call Boundary2(Q,IS,IE,JS,JE)
	EndIf
    
    Do M=1,1
    Do J=JS-3,JE+3
    Do I=IS-3,IE+3
    
        C(M,I,J)=dSQRT(GAM(M)*(P(M,I,J)+PA(M))/Q(M,1,1,I,J))
    
    End Do
    End Do
    End Do
    
    Call Lxx(C,Q,IS,IE,JS,JE,LX,DX,GAM,PA)
    Call Lyy(C,Q,IS,IE,JS,JE,LY,DY,GAM,PA)
    Call Lss(Q,U,V,P,RR,LS,IS,IE,JS,JE)

    Do M=1,1
	Do L=1,4
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
		Q(M,1,L,I,J)=(3.0d0*Q(M,2,L,I,J)+Q(M,1,L,I,J)+DT*(LX(M,L,I,J)+LY(M,L,I,J)+LS(M,L,I,J)))/4.0d0
    
	End Do	
	End Do
	End Do
	End Do
	
	If(LL.EQ.1) Then
		Call Boundary(Q,IS,IE,JS,JE,F,DX)
	Else
		Call Boundary2(Q,IS,IE,JS,JE)
	EndIf
    Call Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI1,FJ1,NMAX,0,1)
!	Call Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI0,FJ0,NMAX,0,0)
    If(LL.EQ.1) Then
		Call Boundary(Q,IS,IE,JS,JE,F,DX)
	Else
		Call Boundary2(Q,IS,IE,JS,JE)
	EndIf
    
    Do M=1,1
    Do J=JS-3,JE+3
    Do I=IS-3,IE+3
    
        C(M,I,J)=dSQRT(GAM(M)*(P(M,I,J)+PA(M))/Q(M,1,1,I,J))
    
    End Do
    End Do
    End Do
    
    Call Lxx(C,Q,IS,IE,JS,JE,LX,DX,GAM,PA)
    Call Lyy(C,Q,IS,IE,JS,JE,LY,DY,GAM,PA)
    Call Lss(Q,U,V,P,RR,LS,IS,IE,JS,JE)

    Do M=1,1
	Do L=1,4
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
	
		Q(M,1,L,I,J)=(Q(M,2,L,I,J)+2.0d0*Q(M,1,L,I,J)+2.0d0*DT*(LX(M,L,I,J)+LY(M,L,I,J)+LS(M,L,I,J)))/3.0d0
		
	End Do
	End Do
	End Do
	End Do

!	Call Gf(DX,DY,F,U,V,S,P,Q,IS,IE,JS,JE,GAM,PA,RR,FI,FJ,NMAX,1)

	End Subroutine Eluer
	
	Subroutine Calculation(DT,DX,DY,U,V,P,Q,F,IS,IE,JS,JE,GAM,PA,RR,MM,LL,L,SDT,NL,R2,A2,B2)

	Implicit None
	
	Integer :: I
    Integer :: IS
    Integer :: IE
    Integer :: J
    Integer :: JS
    Integer :: JE
    Integer :: M
	Integer :: MM
	Integer :: LL
	Integer :: L
	Integer :: EPSJ
	Integer :: N
	Integer::NN
	Integer :: NL
	Integer :: ML
	Integer :: II1
	Integer :: II2

	Integer,Parameter :: MA=1
	Double Precision :: DT
	Double Precision :: DX
	Double Precision :: DY
	Double Precision :: SDT
	Double Precision :: VMAX
	Double Precision :: XMAX
	Double Precision :: YMAX
	Double Precision :: EPS
	Double precision :: aaa
	Double precision :: bbb
	Double precision :: ccc
	Double precision :: pgrad
	Double precision :: pgradX
	Double precision :: pgradY
	Double precision :: pgrad1
	Double precision :: I1
	Double precision :: PP
	Double precision :: ca
	Double precision :: sumvv
	Double precision :: sumvvv
	Double precision,Parameter :: pai=3.141592653589793
	Double precision :: ff
	Double precision :: A2
	Double precision :: B2
	Double precision :: R2
	
	Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxMAT=4
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(:) :: GAM
	Double Precision,Dimension(:) :: PA
	Double Precision,Dimension(:) :: RR
	Double Precision,Dimension(:,:) :: F
	Double Precision,Allocatable,Dimension(:,:) :: K
	Double Precision,Allocatable,Dimension(:) :: a
	Double Precision,Allocatable,Dimension(:) :: b
	Double Precision,Allocatable,Dimension(:) :: d
	Double Precision,Allocatable,Dimension(:) :: e
	Double Precision,Allocatable,Dimension(:) :: g
	Double Precision,Dimension(:,:,:) :: U
	Double Precision,Dimension(:,:,:) :: V
	Double Precision,Dimension(:,:,:) :: P
	Double Precision,Allocatable,Dimension(:,:,:) :: S
	Double Precision,Allocatable,Dimension(:,:,:) :: C
	Double Precision,Dimension(:,:,:,:,:) :: Q
	Double Precision,Allocatable,Dimension(:) :: Jb
	Double Precision,Allocatable,Dimension(:) :: nz
	Double Precision,Allocatable,Dimension(:) :: nr

	Double Precision,Allocatable,Dimension(:) :: X
	Double Precision,Allocatable,Dimension(:) :: Y
	Double Precision,Allocatable,Dimension(:,:) :: Q1A
	Double Precision,Allocatable,Dimension(:,:) :: Q2A
	Double Precision,Allocatable,Dimension(:,:) :: Q3A
	Double Precision,Allocatable,Dimension(:,:) :: Q4A

	Double Precision,Allocatable,Dimension(:,:,:) :: UA
	Double Precision,Allocatable,Dimension(:,:,:) :: VA

	Double Precision,Allocatable,Dimension(:,:) :: FP
	
	Allocate(K(IE+3,JE+3))
	Allocate(a(IE+3))
	Allocate(b(IE+3))
	Allocate(d(IE+3))
	Allocate(e(IE+3))
	Allocate(S(MaxM,IE+3,JE+3))
	Allocate(C(MaxM,IE+3,JE+3))
	Allocate(X(IE+3))
	Allocate(Y(JE+3))
	Allocate(Q1A(IE+3,JE+3))
	Allocate(Q2A(IE+3,JE+3))
	Allocate(Q3A(IE+3,JE+3))
	Allocate(Q4A(IE+3,JE+3))
	Allocate(UA(MaxM,IE+3,JE+3))
	Allocate(VA(MaxM,IE+3,JE+3))
	Allocate(Jb(IE+3))
	Allocate(g(IE+3))
	Allocate(FP(IE+3,JE+3))
	Allocate(nz(IE+3))
	Allocate(nr(IE+3))
	
   
    Call Eluer(DT,DX,DY,U,V,P,Q,F,IS,IE,JS,JE,GAM,PA,RR,LL)

	Do M=1,2
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
		UA(M,I,J)=0.0d0
		VA(M,I,J)=-1.0d0
	End Do
	End Do
	End Do

!	If(F(IS+1,JS+1).GT.-0.4d0*DX)Then
!    Do J=JS+1,JE-1
!	Do I=IS+1,IE-1
!		F(I,J)=F(I,J)-DT
!	End Do
!	End Do
!	Else
 
!	Call Levelset(DT,DX,DY,UA,VA,Q,F,IS,IE,JS,JE,LL)
!	EndIf

	B2=B2-DT/DX

	Do J=JS+1,JE
    	Do I=IS+1,IE
    
        	F(I,J)=DSqrt(((Dble(I)-A2)*DX)**2.0d0+((Dble(J)-B2)*DY)**2.0d0)-R2*DX

    	End Do
    	End Do
    
    Do M=1,2
    Do J=JS-3,JE+3
	Do I=IS-3,IE+3

        U(M,I,J)=Q(M,1,2,I,J)/Q(M,1,1,I,J)
        V(M,I,J)=Q(M,1,3,I,J)/Q(M,1,1,I,J)
        P(M,I,J)=(GAM(M)-1.0d0)*(Q(M,1,4,I,J)-0.5d0*Q(M,1,1,I,J)*(U(M,I,J)*U(M,I,J)+V(M,I,J)*V(M,I,J)))-GAM(M)*PA(M)
        S(M,I,J)=(P(M,I,J)+GAM(M)*PA(M))/(Q(M,1,1,I,J)**GAM(M))

	End Do
	End Do
	End Do

	If(F(IS+1,JS+1).LE.0.0d0)Then

	Do I=IS+1,IE-1
	
	If(F(I,JS+1).LE.0.0)Then

	b(I)=0.0d0
	Jb(I)=0.0d0
	pgrad=0.0d0

	Do J=JS+1,JE-1
		
		If(F(I,J).LE.0.0d0)Then

		pgradX=0.5d0*(P(1,I+1,J)-P(1,I-1,J))/DX
		pgradY=0.5d0*(P(1,I,J+1)-P(1,I,J-1))/DY

		pgrad1=dsqrt(pgradX*pgradX+pgradY*pgradY)

		If(pgrad1.GE.pgrad)Then

			pgrad=pgrad1
			Jb(I)=J
			nz(I)=pgradY/pgrad1
			nr(I)=pgradX/pgrad1

		EndIf
		
		Else

		Jb(I)=0.0d0
		nz(I)=0.0d0
		nr(I)=0.0d0

		EndIf
		

	End Do
	
	Else

	EndIf

	End Do

	b(I)=0.0d0
	Do I=IS+1,IE-1

	PP=0.0d0

		ca=ABS(nz(I)/nr(I))


!		If(F(I,Jb(I)).LE.0.0d0)Then

			If(Abs(nz(I)).LE.0.01d0)Then

				Do M=IS+1,I

					PP=P(1,M,Jb(I))

					If(PP.GE.b(I))b(I)=PP

					!If(I.EQ.IE/2)Write(MM*10+6,*)PP

				End Do


			ElseIf(Abs(nz(I)).GE.0.99d0)Then

				Do M=JS+1,Jb(I)

					PP=P(1,I,M)

					If(PP.GE.b(I))b(I)=PP

					!If(I.EQ.IE/2)Write(MM*10+6,*)PP

				End Do

			Else

				Do M=JS+1,Jb(I)

				I1=Dble(I)+(Dble(M)-Jb(I))/ca

				If(I1.GT.IS)Then

				II1=Dble(Int(I1))
				II2=II1+1

				PP=Abs(Dble(II1)-I1)*P(1,II2,M)+Abs(I1-Dble(II2))*P(1,II1,M)

				If(PP.GE.b(I))b(I)=PP

				EndIf

				!If(I.EQ.IE/2)Write(MM*10+6,*)PP

				End Do

			EndIf


!		EndIf


	End Do

	Close(MM*10+6)

	Do I=IS+1,IE-1

		If(F(I,JS+1).LE.0.0d0)Then

			b(I)=P(1,I,JS+1)

		Else

			b(I)=0.0d0

		End If

	End Do


	open (20, file='pre.csv') !, status='new')

	write (20,100) (b(I),I=IS+1,IE-1,MA)
	100 Format(f7.4,2000(',',1x,E13.6))

	Do I=IS+1,IE-1

		If(F(I,JS+1).LE.0.0d0)Then

			b(I)=U(1,I,JS+1)

		Else

			b(I)=100.0d0

		End If

	End Do


	open (21, file='uvel.csv') !, status='new')

	write (21,101) (b(I),I=IS+1,IE-1,MA)
	101 Format(f7.4,2000(',',1x,E13.6))


	Do I=IS+1,IE-1

		If(F(I,JS+1).LE.0.0d0)Then

			b(I)=V(1,I,JS+1)

		Else

			b(I)=100.0d0

		End If

	End Do


	open (22, file='vvel.csv') !, status='new')

	write (22,102) (b(I),I=IS+1,IE-1,MA)
	102 Format(f7.4,2000(',',1x,E13.6))

	Do I=IS+1,IE-1

		If(F(I,JS+1).LE.0.0d0)Then

			b(I)=F(I,JS+1)

		Else

			b(I)=100.0d0

		End If

	End Do


	open (23, file='f.csv') !, status='new')

	write (23,103) (b(I),I=IS+1,IE-1,MA)
	103 Format(f7.4,2000(',',1x,E13.6))




	Do I=IS+1,IE-1

		If(F(I,JS+1).LE.0.0d0)Then

			b(I)=P(1,I,JS+1)

		Else

			b(I)=100.0d0

		End If

	End Do


	open (27, file='pre2.csv') !, status='new')

	write (27,107) (b(I),I=IS+1,IE-1,MA)
	107 Format(f7.4,2000(',',1x,E13.6))

	Do I=IS+1,IE-1

		If(F(I,JS+2).LE.0.0d0)Then

			b(I)=U(1,I,JS+2)

		Else

			b(I)=100.0d0

		End If

	End Do


	open (28, file='uvel2.csv') !, status='new')

	write (28,108) (b(I),I=IS+1,IE-1,MA)
	108 Format(f7.4,2000(',',1x,E13.6))


	Do I=IS+1,IE-1

		If(F(I,JS+2).LE.0.0d0)Then

			b(I)=V(1,I,JS+2)

		Else

			b(I)=100.0d0

		End If

	End Do


	open (29, file='vvel2.csv') !, status='new')

	write (29,109) (b(I),I=IS+1,IE-1,MA)
	109 Format(f7.4,2000(',',1x,E13.6))

	Do I=IS+1,IE-1

		If(F(I,JS+2).LE.0.0d0)Then

			b(I)=F(I,JS+2)

		Else

			b(I)=100.0d0

		End If

	End Do


	open (30, file='f2.csv') !, status='new')

	write (30,110) (b(I),I=IS+1,IE-1,MA)
	110 Format(f7.4,2000(',',1x,E13.6))

	Do I=IS+1,IE-1

		If(P(1,I,JS+1).GE.d(I))Then

			d(I)=P(1,I,JS+1)
			e(I)=U(1,I,JS+1)
			g(I)=V(1,I,JS+1)

		EndIf

	End Do





    If(MM-MM/100*100.EQ.0.0)Then
	    Do J=JS+1,JE-1,4
	    Do I=IS+1,IE-1,4
	    If(F(I,J).LT.0.0)Then
	        Write(MM*10+2,'(i5,i5,F15.5,F15.5)')I,J,P(1,I,J),F(I,J)
	        Write(MM*10+3,'(i5,i5,F15.5,F15.5)')I,J,U(1,I,J),V(1,I,J)
!		    Write(MM*10+4,'(F15.5,F15.5)')F(I,J),Q(1,1,1,I,J)
	    Else
	        Write(MM*10+2,'(i5,i5,F15.5,F15.5)')I,J,P(2,I,J),F(I,J)
	        Write(MM*10+3,'(i5,i5,F15.5,F15.5)')I,J,U(2,I,J),V(2,I,J)
!		    Write(MM*10+4,'(F15.5,F15.5)')F(I,J),Q(2,1,1,I,J)
	    EndIf
	    End Do
	    End Do

		Close(MM*10+2)
	    Close(MM*10+3)
	    Close(MM*10+4)


	Do J=JS+1,JE-1
	Do I=IS+1,IE-1
		
		X(I)=I
		Y(J)=J
	
		If(F(I,J).LE.0.0)Then
			Q1A(I,J)=Q(1,1,1,I,J)
			Q2A(I,J)=Q(1,1,2,I,J)
			Q3A(I,J)=Q(1,1,3,I,J)
			Q4A(I,J)=Q(1,1,4,I,J)
		Else
			Q1A(I,J)=Q(2,1,1,I,J)
			Q2A(I,J)=Q(2,1,2,I,J)
			Q3A(I,J)=Q(2,1,3,I,J)
			Q4A(I,J)=Q(2,1,4,I,J)
		EndIf
	
	End Do
	End Do

	Open(MM*10,status='unknown',form='binary')
     
	Do J=JS+1,JE-1
	Do I=IS+1,IE-1

		If(LL.EQ.1) Write(MM*10) X(I),Y(J),F(I,J),Q1A(I,J),Q2A(I,J),Q3A(I,J),Q4A(I,J)
	
	End Do
	End Do
      
	Close(MM*10)
	

	EndIf



    
    If(F((IE-1),JS+1).LE.0.0d0)Then
	open (24, file='p.csv') !, status='new')

	write (24,104) (d(I),I=IS+1,IE-1)
	104 Format(f7.4,2000(',',1x,E13.6))

	open (25, file='u.csv') !, status='new')

	write (25,105) (e(I),I=IS+1,IE-1)
	105 Format(f7.4,2000(',',1x,E13.6))

	open (26, file='v.csv') !, status='new')

	write (26,106) (g(I),I=IS+1,IE-1)
	106 Format(f7.4,2000(',',1x,E13.6))

	stop

	EndIf

	Write(3,*)SDT


!	Close(MM*10+1)


	EndIf
	
    
    End Subroutine Calculation
    
    
    Subroutine LayerBoundaryA_B(ISA,JSA,IEA,JEA,ISB,JSB,IEB,JEB,QA,QB,FA,FB)

	Implicit None

	Integer :: I
	Integer :: J
	Integer :: M
	Integer :: L
	Integer :: N
	Integer :: NN
	Integer :: ISA
	Integer :: JSA
	Integer :: IEA
	Integer :: JEA
	Integer :: ISB
	Integer :: JSB
	Integer :: IEB
	Integer :: JEB
	
	Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(:,:) :: FA
	Double Precision,Dimension(:,:) :: FB
	Double Precision,Dimension(:,:,:,:,:) :: QA
	Double Precision,Dimension(:,:,:,:,:) :: QB

    Do M=1,2
	
	Do L=1,4
	
	Do N=1,2

	NN=0
    
	Do I=ISB-3,IEB+1,2
	
		QB(M,1,L,I,JSB-N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JSA-(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JSA+1-(N-1))+QA(M,1,L,ISA-2+NN,JSA-(N-1)))+QA(M,1,L,ISA-2+NN,JSA+1-(N-1)))/16.0d0
		QB(M,1,L,I,JSB-1-N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JSA-(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JSA-1-(N-1))+QA(M,1,L,ISA-2+NN,JSA-(N-1)))+QA(M,1,L,ISA-2+NN,JSA-1-(N-1)))/16.0d0
		QB(M,1,L,I+1,JSB-N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JSA-(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JSA+1-(N-1))+QA(M,1,L,ISA+NN,JSA-(N-1)))+QA(M,1,L,ISA+NN,JSA+1-(N-1)))/16.0d0
		QB(M,1,L,I+1,JSB-1-N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JSA-(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JSA-1-(N-1))+QA(M,1,L,ISA+NN,JSA-(N-1)))+QA(M,1,L,ISA+NN,JSA-1-(N-1)))/16.0d0	
		QB(M,1,L,I,JEB+N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JEA+(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JEA-1+(N-1))+QA(M,1,L,ISA-2+NN,JEA+(N-1)))+QA(M,1,L,ISA-2+NN,JEA-1+(N-1)))/16.0d0		
		QB(M,1,L,I,JEB+1+N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JEA+(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JEA+1+(N-1))+QA(M,1,L,ISA-2+NN,JEA+(N-1)))+QA(M,1,L,ISA-2+NN,JEA+1+(N-1)))/16.0d0
		QB(M,1,L,I+1,JEB+N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JEA+(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JEA-1+(N-1))+QA(M,1,L,ISA+NN,JEA+(N-1)))+QA(M,1,L,ISA+NN,JEA-1+(N-1)))/16.0d0
		QB(M,1,L,I+1,JEB+1+N*(N-1))=(9.0d0*QA(M,1,L,ISA-1+NN,JEA+(N-1))+3.0d0*(QA(M,1,L,ISA-1+NN,JEA+1+(N-1))+QA(M,1,L,ISA+NN,JEA+(N-1)))+QA(M,1,L,ISA+NN,JEA+1+(N-1)))/16.0d0
        FB(I,JSB-N*(N-1))=(9.0d0*FA(ISA-1+NN,JSA-(N-1))+3.0d0*(FA(ISA-1+NN,JSA+1-(N-1))+FA(ISA-2+NN,JSA-(N-1)))+FA(ISA-2+NN,JSA+1-(N-1)))/16.0d0
		FB(I,JSB-1-N*(N-1))=(9.0d0*FA(ISA-1+NN,JSA-(N-1))+3.0d0*(FA(ISA-1+NN,JSA-1-(N-1))+FA(ISA-2+NN,JSA-(N-1)))+FA(ISA-2+NN,JSA-1-(N-1)))/16.0d0
		FB(I+1,JSB-N*(N-1))=(9.0d0*FA(ISA-1+NN,JSA-(N-1))+3.0d0*(FA(ISA-1+NN,JSA+1-(N-1))+FA(ISA+NN,JSA-(N-1)))+FA(ISA+NN,JSA+1-(N-1)))/16.0d0
		FB(I+1,JSB-1-N*(N-1))=(9.0d0*FA(ISA-1+NN,JSA-(N-1))+3.0d0*(FA(ISA-1+NN,JSA-1-(N-1))+FA(ISA+NN,JSA-(N-1)))+FA(ISA+NN,JSA-1-(N-1)))/16.0d0
		FB(I,JEB+N*(N-1))=(9.0d0*FA(ISA-1+NN,JEA+(N-1))+3.0d0*(FA(ISA-1+NN,JEA-1+(N-1))+FA(ISA-2+NN,JEA+(N-1)))+FA(ISA-2+NN,JEA-1+(N-1)))/16.0d0
		FB(I,JEB+1+N*(N-1))=(9.0d0*FA(ISA-1+NN,JEA+(N-1))+3.0d0*(FA(ISA-1+NN,JEA+1+(N-1))+FA(ISA-2+NN,JEA+(N-1)))+FA(ISA-2+NN,JEA+1+(N-1)))/16.0d0
		FB(I+1,JEB+N*(N-1))=(9.0d0*FA(ISA-1+NN,JEA+(N-1))+3.0d0*(FA(ISA-1+NN,JEA-1+(N-1))+FA(ISA+NN,JEA+(N-1)))+FA(ISA+NN,JEA-1+(N-1)))/16.0d0
		FB(I+1,JEB+1+N*(N-1))=(9.0d0*FA(ISA-1+NN,JEA+(N-1))+3.0d0*(FA(ISA-1+NN,JEA+1+(N-1))+FA(ISA+NN,JEA+(N-1)))+FA(ISA+NN,JEA+1+(N-1)))/16.0d0

    NN=NN+1
    
	End Do

	NN=0
    
	Do J=JSB-3,JEB+1,2
	
		QB(M,1,L,ISB-N*(N-1),J)=(9.0d0*QA(M,1,L,ISA-(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,ISA+1-(N-1),JSA-1+NN)+QA(M,1,L,ISA-(N-1),JSA-2+NN))+QA(M,1,L,ISA+1-(N-1),JSA-2+NN))/16.0d0
		QB(M,1,L,ISB-1-N*(N-1),J)=(9.0d0*QA(M,1,L,ISA-(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,ISA-1-(N-1),JSA-1+NN)+QA(M,1,L,ISA-(N-1),JSA-2+NN))+QA(M,1,L,ISA-1-(N-1),JSA-2+NN))/16.0d0
		QB(M,1,L,ISB-N*(N-1),J+1)=(9.0d0*QA(M,1,L,ISA-(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,ISA+1-(N-1),JSA-1+NN)+QA(M,1,L,ISA-(N-1),JSA+NN))+QA(M,1,L,ISA+1-(N-1),JSA+NN))/16.0d0
		QB(M,1,L,ISB-1-N*(N-1),J+1)=(9.0d0*QA(M,1,L,ISA-(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,ISA-1-(N-1),JSA-1+NN)+QA(M,1,L,ISA-(N-1),JSA+NN))+QA(M,1,L,ISA-1-(N-1),JSA+NN))/16.0d0
		QB(M,1,L,IEB+N*(N-1),J)=(9.0d0*QA(M,1,L,IEA+(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,IEA-1+(N-1),JSA-1+NN)+QA(M,1,L,IEA+(N-1),JSA-2+NN))+QA(M,1,L,IEA-1+(N-1),JSA-2+NN))/16.0d0
		QB(M,1,L,IEB+1+N*(N-1),J)=(9.0d0*QA(M,1,L,IEA+(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,IEA+1+(N-1),JSA-1+NN)+QA(M,1,L,IEA+(N-1),JSA-2+NN))+QA(M,1,L,IEA+1+(N-1),JSA-2+NN))/16.0d0
		QB(M,1,L,IEB+N*(N-1),J+1)=(9.0d0*QA(M,1,L,IEA+(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,IEA-1+(N-1),JSA-1+NN)+QA(M,1,L,IEA+(N-1),JSA+NN))+QA(M,1,L,IEA-1+(N-1),JSA+NN))/16.0d0
		QB(M,1,L,IEB+1+N*(N-1),J+1)=(9.0d0*QA(M,1,L,IEA+(N-1),JSA-1+NN)+3.0d0*(QA(M,1,L,IEA+1+(N-1),JSA-1+NN)+QA(M,1,L,IEA+(N-1),JSA+NN))+QA(M,1,L,IEA+1+(N-1),JSA+NN))/16.0d0    
        FB(ISB-N*(N-1),J)=(9.0d0*FA(ISA-(N-1),JSA-1+NN)+3.0d0*(FA(ISA+1-(N-1),JSA-1+NN)+FA(ISA-(N-1),JSA-2+NN))+FA(ISA+1-(N-1),JSA-2+NN))/16.0d0
		FB(ISB-1-N*(N-1),J)=(9.0d0*FA(ISA-(N-1),JSA-1+NN)+3.0d0*(FA(ISA-1-(N-1),JSA-1+NN)+FA(ISA-(N-1),JSA-2+NN))+FA(ISA-1-(N-1),JSA-2+NN))/16.0d0
		FB(ISB-N*(N-1),J+1)=(9.0d0*FA(ISA-(N-1),JSA-1+NN)+3.0d0*(FA(ISA+1-(N-1),JSA-1+NN)+FA(ISA-(N-1),JSA+NN))+FA(ISA+1-(N-1),JSA+NN))/16.0d0
		FB(ISB-1-N*(N-1),J+1)=(9.0d0*FA(ISA-(N-1),JSA-1+NN)+3.0d0*(FA(ISA-1-(N-1),JSA-1+NN)+FA(ISA-(N-1),JSA+NN))+FA(ISA-1-(N-1),JSA+NN))/16.0d0
		FB(IEB+N*(N-1),J)=(9.0d0*FA(IEA+(N-1),JSA-1+NN)+3.0d0*(FA(IEA-1+(N-1),JSA-1+NN)+FA(IEA+(N-1),JSA-2+NN))+FA(IEA-1+(N-1),JSA-2+NN))/16.0d0
		FB(IEB+1+N*(N-1),J)=(9.0d0*FA(IEA+(N-1),JSA-1+NN)+3.0d0*(FA(IEA+1+(N-1),JSA-1+NN)+FA(IEA+(N-1),JSA-2+NN))+FA(IEA+1+(N-1),JSA-2+NN))/16.0d0
		FB(IEB+N*(N-1),J+1)=(9.0d0*FA(IEA+(N-1),JSA-1+NN)+3.0d0*(FA(IEA-1+(N-1),JSA-1+NN)+FA(IEA+(N-1),JSA+NN))+FA(IEA-1+(N-1),JSA+NN))/16.0d0
		FB(IEB+1+N*(N-1),J+1)=(9.0d0*FA(IEA+(N-1),JSA-1+NN)+3.0d0*(FA(IEA+1+(N-1),JSA-1+NN)+FA(IEA+(N-1),JSA+NN))+FA(IEA+1+(N-1),JSA+NN))/16.0d0

    NN=NN+1
    
	End Do

	End Do
	
	End Do
	
	End Do
	
	End Subroutine LayerBoundaryA_B


	Subroutine BoundaryLayerB_A(ISA,JSA,IEA,JEA,ISB,JSB,IEB,JEB,QA,QB,FA,FB)

	Implicit None

	Integer :: I
	Integer :: J
	Integer :: M
	Integer :: L
	Integer :: NN
	Integer :: MM

	Integer :: ISA
	Integer :: IEA
	Integer :: JSA
	Integer :: JEA
	Integer :: ISB
	Integer :: IEB
	Integer :: JSB
	Integer :: JEB

	Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(:,:) :: FA
	Double Precision,Dimension(:,:) :: FB
	Double Precision,Dimension(:,:,:,:,:) :: QA
	Double Precision,Dimension(:,:,:,:,:) :: QB

	NN=1

	Do I=ISA+1,IEA-1

	MM=1

	Do J=JSA+1,JEA-1
	
	Do M=1,2
	
	Do L=1,4
	
	
		QA(M,1,L,I,J)=0.25d0*(QB(M,1,L,ISB+NN,JSB+MM)+QB(M,1,L,ISB+1+NN,JSB+MM)+QB(M,1,L,ISB+NN,JSB+1+MM)+QB(M,1,L,ISB+1+NN,JSB+1+MM))
    
!        QA(M,1,L,I,J)=(9.0*(QB(M,1,L,ISB+NN,JSB+MM)+QB(M,1,L,ISB+1+NN,JSB+MM)+QB(M,1,L,ISB+NN,JSB+1+MM)+QB(M,1,L,ISB+1+NN,JSB+1+MM))+3.0*(QB(M,1,L,ISB-1+NN,JSB+MM)+QB(M,1,L,ISB+NN,JSB-1+MM)+QB(M,1,L,ISB+2+NN,JSB+MM)+QB(M,1,L,ISB+1+NN,JSB-1+MM)+QB(M,1,L,ISB-1+NN,JSB+1+MM)+QB(M,1,L,ISB+NN,JSB+2+MM)+QB(M,1,L,ISB+2+NN,JSB+1+MM)+QB(M,1,L,ISB+1+NN,JSB+2+MM))+QB(M,1,L,ISB-1+NN,JSB-1+MM)+QB(M,1,L,ISB+2+NN,JSB-1+MM)+QB(M,1,L,ISB-1+NN,JSB+2+MM)+QB(M,1,L,ISB+2+NN,JSB+2+MM))/64.0
	
	End Do
	
	End Do
	
	    FA(I,J)=0.25d0*(FB(ISB+NN,JSB+MM)+FB(ISB+1+NN,JSB+MM)+FB(ISB+NN,JSB+1+MM)+FB(ISB+1+NN,JSB+1+MM))
	    
!        FA(I,J)=(9.0*(FB(ISB+NN,JSB+MM)+FB(ISB+1+NN,JSB+MM)+FB(ISB+NN,JSB+1+MM)+FB(ISB+1+NN,JSB+1+MM))+3.0*(FB(ISB-1+NN,JSB+MM)+FB(ISB+NN,JSB-1+MM)+FB(ISB+2+NN,JSB+MM)+FB(ISB+1+NN,JSB-1+MM)+FB(ISB-1+NN,JSB+1+MM)+FB(ISB+NN,JSB+2+MM)+FB(ISB+2+NN,JSB+1+MM)+FB(ISB+1+NN,JSB+2+MM))+FB(ISB-1+NN,JSB-1+MM)+FB(ISB+2+NN,JSB-1+MM)+FB(ISB-1+NN,JSB+2+MM)+FB(ISB+2+NN,JSB+2+MM))/64.0
	
		MM=MM+2
	
	End Do

	NN=NN+2

	End Do


	End Subroutine BoundaryLayerB_A
	
	
    Subroutine BoundaryLayerA_B(ISA,JSA,IEA,JEA,ISB,JSB,IEB,JEB,QA,QB,FA,FB)
	
	Implicit None
	
	Integer :: I
	Integer :: J
	Integer :: M
	Integer :: L
	Integer :: NN
	Integer :: MM

	Integer :: ISA
	Integer :: IEA
	Integer :: JSA
	Integer :: JEA
	Integer :: ISB
	Integer :: IEB
	Integer :: JSB
	Integer :: JEB
	
	Integer,Parameter :: MaxM=2
!	Integer,Parameter :: MaxX=800
!	Integer,Parameter :: MaxY=800
	Integer,Parameter :: MaxN=2
	Integer,Parameter :: MaxL=4
	Double Precision,Dimension(:,:) :: FA
	Double Precision,Dimension(:,:) :: FB
	Double Precision,Dimension(:,:,:,:,:) :: QA
	Double Precision,Dimension(:,:,:,:,:) :: QB
	
	MM=0
	
	Do J=JSB+1,JEB-1,2
	
	NN=0
	
	Do I=ISB+1,IEB-1,2
	
	Do M=1,2
	
	Do L=1,4
	
	    
		QB(M,1,L,I,J)=(9.0d0*QA(M,1,L,ISA+1+NN,JSA+1+MM)+3.0d0*(QA(M,1,L,ISA+NN,JSA+1+MM)+QA(M,1,L,ISA+1+NN,JSA+MM))+QA(M,1,L,ISA+NN,JSA+MM))/16.0d0
		QB(M,1,L,I,J+1)=(9.0d0*QA(M,1,L,ISA+1+NN,JSA+1+MM)+3.0d0*(QA(M,1,L,ISA+NN,JSA+1+MM)+QA(M,1,L,ISA+1+NN,JSA+2+MM))+QA(M,1,L,ISA+NN,JSA+2+MM))/16.0d0
		QB(M,1,L,I+1,J)=(9.0d0*QA(M,1,L,ISA+1+NN,JSA+1+MM)+3.0d0*(QA(M,1,L,ISA+2+NN,JSA+1+MM)+QA(M,1,L,ISA+1+NN,JSA+MM))+QA(M,1,L,ISA+2+NN,JSA+MM))/16.0d0
		QB(M,1,L,I+1,J+1)=(9.0d0*QA(M,1,L,ISA+1+NN,JSA+1+MM)+3.0d0*(QA(M,1,L,ISA+2+NN,JSA+1+MM)+QA(M,1,L,ISA+1+NN,JSA+2+MM))+QA(M,1,L,ISA+2+NN,JSA+2+MM))/16.0d0
	
	End Do
	
	End Do
		
!		FB(I,J)=(0.769750166235208000d0*FA(ISA+1+NN,JSA+1+MM)+0.206482365662665000d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+0.055388058619298400d0*FA(ISA+NN,JSA+MM)-0.074333651638559300d0*(FA(ISA+1+NN,JSA+2+MM)+FA(ISA+2+NN,JSA+1+MM))-0.019939701102947400d0*(FA(ISA+NN,JSA+2+MM)+FA(ISA+2+NN,JSA+MM))-0.015709187066024700d0*(FA(ISA-1+NN,JSA+1+MM)+FA(ISA+1+NN,JSA-1+MM))-0.004213925829850300d0*(FA(ISA+NN,JSA-1+MM)+FA(ISA-1+NN,JSA+MM))+0.001517013298746110d0*(FA(ISA-1+NN,JSA+2+MM)+FA(ISA+2+NN,JSA-1+MM))+0.007178292397061080d0*FA(ISA+2+NN,JSA+2+MM)+0.000320595654408667d0*FA(ISA-1+NN,JSA-1+MM))/1.02024294d0
!		FB(I,J+1)=(0.769750166235208000d0*FA(ISA+1+NN,JSA+1+MM)+0.206482365662665000d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+0.055388058619298400d0*FA(ISA+NN,JSA+2+MM)-0.074333651638559300d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))-0.019939701102947400d0*(FA(ISA+2+NN,JSA+2+MM)+FA(ISA+NN,JSA+MM))-0.015709187066024700d0*(FA(ISA+1+NN,JSA+3+MM)+FA(ISA-1+NN,JSA+1+MM))-0.004213925829850300d0*(FA(ISA+NN,JSA+3+MM)+FA(ISA-1+NN,JSA+2+MM))+0.001517013298746110d0*(FA(ISA+2+NN,JSA+3+MM)+FA(ISA-1+NN,JSA+MM))+0.007178292397061080d0*FA(ISA+2+NN,JSA+MM)+0.000320595654408667d0*FA(ISA-1+NN,JSA+3+MM))/1.02024294d0
!		FB(I+1,J)=(0.769750166235208000d0*FA(ISA+1+NN,JSA+1+MM)+0.206482365662665000d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+0.055388058619298400d0*FA(ISA+2+NN,JSA+MM)-0.074333651638559300d0*(FA(ISA+1+NN,JSA+2+MM)+FA(ISA+NN,JSA+1+MM))-0.019939701102947400d0*(FA(ISA+2+NN,JSA+2+MM)+FA(ISA+NN,JSA+MM))-0.015709187066024700d0*(FA(ISA+3+NN,JSA+1+MM)+FA(ISA+1+NN,JSA-1+MM))-0.004213925829850300d0*(FA(ISA+3+NN,JSA+MM)+FA(ISA+2+NN,JSA-1+MM))+0.001517013298746110d0*(FA(ISA+3+NN,JSA+2+MM)+FA(ISA+NN,JSA-1+MM))+0.007178292397061080d0*FA(ISA+NN,JSA+2+MM)+0.000320595654408667d0*FA(ISA+3+NN,JSA-1+MM))/1.02024294d0
!		FB(I+1,J+1)=(0.769750166235208000d0*FA(ISA+1+NN,JSA+1+MM)+0.206482365662665000d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+0.055388058619298400d0*FA(ISA+2+NN,JSA+2+MM)-0.074333651638559300d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))-0.019939701102947400d0*(FA(ISA+NN,JSA+2+MM)+FA(ISA+2+NN,JSA+MM))-0.015709187066024700d0*(FA(ISA+1+NN,JSA+3+MM)+FA(ISA+3+NN,JSA+1+MM))-0.004213925829850300d0*(FA(ISA+2+NN,JSA+3+MM)+FA(ISA+3+NN,JSA+2+MM))+0.001517013298746110d0*(FA(ISA+NN,JSA+3+MM)+FA(ISA+3+NN,JSA+MM))+0.007178292397061080d0*FA(ISA+NN,JSA+MM)+0.000320595654408667d0*FA(ISA+3+NN,JSA+3+MM))/1.02024294d0
		
!		FB(I,J)=(3249.0d0*FA(ISA+1+NN,JSA+1+MM)+1083.0d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+361.0d0*FA(ISA+NN,JSA+MM)-1323.0d0*(FA(ISA+1+NN,JSA+2+MM)+FA(ISA+2+NN,JSA+1+MM))-225.0d0*(FA(ISA+NN,JSA+2+MM)+FA(ISA+2+NN,JSA+MM))-441.0d0*(FA(ISA-1+NN,JSA+1+MM)+FA(ISA+1+NN,JSA-1+MM))-75.0d0*(FA(ISA+NN,JSA-1+MM)+FA(ISA-1+NN,JSA+MM))+27.0d0*(FA(ISA-1+NN,JSA+2+MM)+FA(ISA+2+NN,JSA-1+MM))+81.0d0*FA(ISA+2+NN,JSA+2+MM)+9.0d0*FA(ISA-1+NN,JSA-1+MM))/1792.0d0
!		FB(I,J+1)=(3249.0d0*FA(ISA+1+NN,JSA+1+MM)+1083.0d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+361.0d0*FA(ISA+NN,JSA+2+MM)-1323.0d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))-225.0d0*(FA(ISA+2+NN,JSA+2+MM)+FA(ISA+NN,JSA+MM))-441.0d0*(FA(ISA+1+NN,JSA+3+MM)+FA(ISA-1+NN,JSA+1+MM))-75.0d0*(FA(ISA+NN,JSA+3+MM)+FA(ISA-1+NN,JSA+2+MM))+27.0d0*(FA(ISA+2+NN,JSA+3+MM)+FA(ISA-1+NN,JSA+MM))+81.0d0*FA(ISA+2+NN,JSA+MM)+9.0d0*FA(ISA-1+NN,JSA+3+MM))/1792.0d0
!		FB(I+1,J)=(3249.0d0*FA(ISA+1+NN,JSA+1+MM)+1083.0d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+361.0d0*FA(ISA+2+NN,JSA+MM)-1323.0d0*(FA(ISA+1+NN,JSA+2+MM)+FA(ISA+NN,JSA+1+MM))-225.0d0*(FA(ISA+2+NN,JSA+2+MM)+FA(ISA+NN,JSA+MM))-441.0d0*(FA(ISA+3+NN,JSA+1+MM)+FA(ISA+1+NN,JSA-1+MM))-75.0d0*(FA(ISA+3+NN,JSA+MM)+FA(ISA+2+NN,JSA-1+MM))+27.0d0*(FA(ISA+3+NN,JSA+2+MM)+FA(ISA+NN,JSA-1+MM))+81.0d0*FA(ISA+NN,JSA+2+MM)+9.0d0*FA(ISA+3+NN,JSA-1+MM))/1792.0d0
!		FB(I+1,J+1)=(3249.0d0*FA(ISA+1+NN,JSA+1+MM)+1083.0d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+361.0d0*FA(ISA+2+NN,JSA+2+MM)-1323.0d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))-225.0d0*(FA(ISA+NN,JSA+2+MM)+FA(ISA+2+NN,JSA+MM))-441.0d0*(FA(ISA+1+NN,JSA+3+MM)+FA(ISA+3+NN,JSA+1+MM))-75.0d0*(FA(ISA+2+NN,JSA+3+MM)+FA(ISA+3+NN,JSA+2+MM))+27.0d0*(FA(ISA+NN,JSA+3+MM)+FA(ISA+3+NN,JSA+MM))+81.0d0*FA(ISA+NN,JSA+MM)+9.0d0*FA(ISA+3+NN,JSA+3+MM))/1792.0d0
		FB(I,J)=(9.0d0*FA(ISA+1+NN,JSA+1+MM)+3.0d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+FA(ISA+NN,JSA+MM))/16.0d0
		FB(I,J+1)=(9.0d0*FA(ISA+1+NN,JSA+1+MM)+3.0d0*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+FA(ISA+NN,JSA+2+MM))/16.0d0
		FB(I+1,J)=(9.0d0*FA(ISA+1+NN,JSA+1+MM)+3.0d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+FA(ISA+2+NN,JSA+MM))/16.0d0
		FB(I+1,J+1)=(9.0d0*FA(ISA+1+NN,JSA+1+MM)+3.0d0*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+FA(ISA+2+NN,JSA+2+MM))/16.0d0
!		FB(I,J)=(dsqrt(18.0d0)*FA(ISA+1+NN,JSA+1+MM)+dsqrt(10.0d0)*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+dsqrt(2.0d0)*FA(ISA+NN,JSA+MM))/(dsqrt(18.0d0)+2.0d0*dsqrt(10.0d0)+dsqrt(2.0d0))
!		FB(I,J+1)=(dsqrt(18.0d0)*FA(ISA+1+NN,JSA+1+MM)+dsqrt(10.0d0)*(FA(ISA+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+dsqrt(2.0d0)*FA(ISA+NN,JSA+2+MM))/(dsqrt(18.0d0)+2.0d0*dsqrt(10.0d0)+dsqrt(2.0d0))
!		FB(I+1,J)=(dsqrt(18.0d0)*FA(ISA+1+NN,JSA+1+MM)+dsqrt(10.0d0)*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+MM))+dsqrt(2.0d0)*FA(ISA+2+NN,JSA+MM))/(dsqrt(18.0d0)+2.0d0*dsqrt(10.0d0)+dsqrt(2.0d0))
!		FB(I+1,J+1)=(dsqrt(18.0d0)*FA(ISA+1+NN,JSA+1+MM)+dsqrt(10.0d0)*(FA(ISA+2+NN,JSA+1+MM)+FA(ISA+1+NN,JSA+2+MM))+dsqrt(2.0d0)*FA(ISA+2+NN,JSA+2+MM))/(dsqrt(18.0d0)+2.0d0*dsqrt(10.0d0)+dsqrt(2.0d0))
	
	NN=NN+1
	
	End Do
	
	MM=MM+1
	
	End Do
	
	End Subroutine BoundaryLayerA_B
    End Program