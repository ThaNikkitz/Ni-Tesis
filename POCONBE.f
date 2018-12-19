	Program POCONBE
	Implicit none

	Character*10 FILEIN, FILEOUT

	Real*8 :: X(101), Y(101), XM(100), YM(100), F1(100), DF1(100)
	Real*8 :: KODE(100), CX(20), CY(20), POT(20), FLUX1(20), FLUX2(20)

	Real*8 :: G(100,100), H(100,100)
	Real*8 :: D, A, B, CH, DQ1, DQ2, DU1, DU2, XP, YP, X1, Y1, X2, Y2

	Integer :: NX, I, J, K, KK
	Integer :: N, L, INP, IPR

	NX = 100

	INP = 5
	IPR = 6

	Write(*, ' (A) ') ' Name of Input File (Max. 10 Chart.)'
	Read(*, ' (A) ') FILEIN
	Open(INP, File = FILEIN, Status = 'OLD')
	Write(*, ' (A) ') 'Name of the OUTPUT file (Max. 10 Chart.)'
	Read(*, ' (A) ') FILEOUT
	Open(IPR, File = FILEOUT, Status = 'NEW')

	Call INPUTPC(CX, CY, X, Y, KODE, F1)

	Call GHMATPC(X, Y, XM, YM, G, H, F1, DF1, KODE, NX)

	Call SLNPD(G, DF1, D, N, NX)

	Call INTERPC(F1, DF1, KODE, CX, CY, X, Y, POT, FLUX1, FLUX2)

	Call OUTPTPC(XM, YM, F1, DF1, CX, CY, POT, FLUX1, FLUX2)

	Close (IPR)
	Close (INP)
	Stop
	End



	!___________________________________________________________________



	Subroutine INPUTPC(CX, CY, X, Y, KODE, F1)
	Implicit None


	Real*8 :: XM(100), YM(100), DF1(100)
	Real*8 :: POT(20), FLUX1(20), FLUX2(20)

	Real*8 :: G(100,100), H(100,100)
	Real*8 :: D, A, B, CH, DQ1, DQ2, DU1, DU2, XP, YP, X1, Y1, X2, Y2

	Integer :: NX, I, J, K, KK
	Integer :: N, L, INP, IPR

	Character*80 Title
	Common N, L, INP, IPR
	Real*8 :: CX(1), CY(1), X(1), Y(1), KODE(1), F1(1)

	Write(IPR, 100)
100	Format(' ',79('*'))

	Read(INP, ' (A) ') Title
	Write(IPR, ' (A) ') Title

	Read(INP, *)N, L
	Write(IPR, 300)N, L
300	Format (//' DATA'//2X, 'NUMBER OF BOUNDARY ELEMENTS =', I3/2X,
	1'NUMBER OF INTERNAL POINTS WHERE THE FUNCTION IS'
	1'CALCULATED =', I3)

	Write(IPR, 500)
500	Format(2X, 'COORDINATES OF THE EXTREME POINTS OF THE'
	1'BOUNDARY ELEMENTS', 1X, 'POINT', 7X, 'X', 15X, 'Y')

	Read(INP, *) (X(I), Y(I), I=1,N)
	Do 10 I=1,N
10	Write(IPR, 700)I, X(I), Y(I)
700	Format(2X, I3, 2(2X,E14.5))

	Write(IPR, 800)
800	Format(2X, 'BOUNDARY CONDITIONS', 2X, 'NODE', 6X, 'CO'
	1'DE', 7X, 'PRESCRIBED VALUE')

	Do 20 I=1,N
	Read(INP,*) KODE(I), F1(I)
20	Write(IPR, 950)I, KODE(I), F1(I)
950	Format(2X, I3, 9X, I1, 8X, E14.5)

	If(L .EQ. 0) GO TO 30
	Read(INP,*) (CX(I), CY(I), I=1,L)
30	Return
	End



	!___________________________________________________________________




	Subroutine GHMATPC(X, Y, XM, YM, G, H, F1, DF1, KODE, NX)
	Implicit None

	Real*8 :: CX(20), CY(20), POT(20), FLUX1(20), FLUX2(20)
	Real*8 :: D, A, B, CH, DQ1, DQ2, DU1, DU2, XP, YP, X1, Y1, X2, Y2

	Integer :: NX, I, J, K, KK
	Integer :: N, L, INP, IPR

	Common N, L, INP, IPR
	Real*8 :: X(1), Y(1), XM(1), YM(1), F1(1), KODE(1)
	Real*8 :: DF1(1), G(NX,NX), H(NX,NX)

	X(N+1) = X(1)
	Y(N+1) = Y(1)
	Do 10 I=1,N
	XM(I) = (X(I) + X(I+1))/2
10	YM(I) = (Y(I) + Y(I+1))/2

	Do 30 I = 1,N
	Do 30 J = 1,N
	KK = J+1
	If (I-J)20, 25, 20
20	Call EXTINPC(XM(I), YM(I), X(J), Y(J), X(KK), Y(KK), H(I,J),
	1 G(I, J), DQ1, DQ2, DU1, DU2, 0)
	Go to 30
25	Call LOCINPC(X(J), Y(J), X(KK), Y(KK), G(I,J))
	H(I,J) = 3.1415926
30	Continue

	Do 55 J = 1,N
	IF(KODE(J))55, 55, 40
40	Do 50 I = 1,N
	CH = G(I,J)
	G(I,J) = -H(I,J)
	H(I,J) = -CH
50	Continue
55	Continue

	Do 60 I = 1,N
	DF1(I) = 0.
	Do 60 J = 1,N
	DF1(I) = DF1(I) + H(I,J)*F1(J)
60	Continue
	Return
	End



	!___________________________________________________________________




	Subroutine EXTINPC(XP, YP, X1, Y1, X2, Y2, H, G, DQ1, DQ2,
	1 DU1, DU2, K)
	Implicit None

	Real*8 :: X(101), Y(101), XM(100), YM(100), F1(100), DF1(100)
	Real*8 :: KODE(100), CX(20), CY(20), POT(20), FLUX1(20), FLUX2(20)

	Real*8 :: G, H
	Real*8 :: D, A, B, CH, DQ1, DQ2, DU1, DU2, XP, YP, X1, Y1, X2, Y2
	Real*8 :: AX, BX, AY, BY, SL, ETA1, ETA2, RD1, RD2, RDN
	Real*4 :: RA

	Integer :: NX, I, J, K, KK
	Integer :: N, L, INP, IPR

	Real*8 :: XCO(4), YCO(4), GI(4), OME(4)
	Data GI/0.86113631, -0.86113631, 0.33998104, -0.33998104/
	Data OME/0.34785485, 0.34785185, 0.65214515, 0.65214515/

	AX = (X2 - X1)/2.
	BX = (X2 + X1)/2.
	AY = (Y2 - Y1)/2.
	BY = (Y2 + Y1)/2.
	SL = SQRT(AX**2 + AY**2)
	ETA1 = AY/SL
	ETA2 = -AX/SL
	G = 0.
	H = 0.
	DU1 = 0.
	DU2 = 0.
	DQ1 = 0.
	DQ2 = 0.

	Do 40 I = 1,4
	XCO(I) = AX*GI(I) + BX
	YCO(I) = AY*GI(I) + BY
	RA = SQRT((XP - XCO(I))**2 + (YP - YCO(I))**2)
	RD1 = (XCO(I) - XP)/RA
	RD2 = (YCO(I) - YP)/RA
	RDN = RD1*ETA1 + RD2*ETA2

	If(K) 30, 30, 10
10	DU1 = DU1 + RD1*OME(I)*SL/RA
	DU2 = DU2 + RD2*OME(I)*SL/RA
	DQ1 = DQ1 - ((2.*RD1**2 - 1.)*ETA1 + 2.*RD1*RD2*ETA2)*OME(I)*SL/RA**2
	DQ2 = DQ2 - ((2.*RD2**2 - 1.)*ETA2 + 2.*RD1*RD2*ETA1)*OME(I)*SL/RA**2
30	G = G + Alog(1/RA)*OME(I)*SL
40	H = H - RDN*OME(I)*SL/RA

	Return
	End



	!___________________________________________________________________




	Subroutine LOCINPC(X1, Y1, X2, Y2, G)
	Implicit None

	Real*8 :: X1, Y1, X2, Y2, G, AX, AY
	Real*4 :: SR

	AX = (X2 - X1)/2.
	AY = (Y2 - Y1)/2.
	SR = Sqrt(AX**2 + AY**2)
	G = 2*SR*(1. - Alog(SR))

	Return
	End



	!___________________________________________________________________



	Subroutine SLNPD(A, B, D, N, NX)
	Implicit None

	Real*8 :: X(101), Y(101), XM(100), YM(100), F1(100), DF1(100)
	Real*8 :: KODE(100), CX(20), CY(20), POT(20), FLUX1(20), FLUX2(20)

	Real*8 :: G(100,100), H(100,100)
	Real*8 :: D, CH, DQ1, DQ2, DU1, DU2, XP, YP, X1, Y1, X2, Y2
	Real*8 :: C, TOL

	Integer :: NX, I, J, K, KK, K1, N1
	Integer :: N, L, INP, IPR

	Real*8 :: B(NX), A(NX,NX)
	TOL = 1.E-6
	N1 = N - 1
	Do 100 K = 1, N1
	K1 = K + 1
	C = A(K,K)

	If(Abs(C) - TOL)1, 1, 3
1	Do 7 J = K1,N

	If(Abs(A(J,K)) - TOL)7, 7, 5
5	Do 6 L = K, N
	C = A(K,L)
	A(K,L) = A(J,L)
6	A(J,L) = C
	C = B(K)
	B(K) = B(J)
	B(J) = C
	C = A(K,K)
	Go to 3
7	Continue
	Go to 8

3	C = A(K,K)
	Do 4 J = K1, N
4	A(K,J) = A(K,J)/C
	B(K) = B(K)/C

	Do 10 I = K1, N
	C = A(I, K)
	Do 9 J = K1, N
9	A(I,J) = A(I,J) - C*A(K,J)
10	B(I) = B(I) - C*B(K)
100	Continue

	If(Abs(A(N,N)) - TOL)8, 8, 101
101	B(N) = B(N)/A(N,N)

	Do 200 L = 1, N1
	K = N - L
	K1 = K + 1
	Do 200 J = K1, N
200	B(K) = B(K) - A(K,J)*B(J)

	D = 1.

	Do 250 I = 1, N
250	D = D*A(I,I)
	Go to 300
8	Write(*,2) K
2	Format(' *** Singularity In Row', I5)
	D = 0.
300	Return
	End



	!___________________________________________________________________



	Subroutine INTERPC(F1, DF1, KODE, CX, CY, X, Y, POT, FLUX1, FLUX2)
	Implicit None

	Real*8 :: G(100,100), H(100,100)
	Real*8 :: D, CH, DQ1, DQ2, DU1, DU2, XP, YP, X1, Y1, X2, Y2
	Real*8 :: A, B

	Integer :: NX, I, J, K, KK
	Integer :: N, L, INP, IPR

	Common N, L, INP, IPR
	Real*8 :: F1(1), DF1(1), KODE(1), CX(1), CY(1), X(1), Y(1)
	Real*8 :: POT(1), FLUX1(1), FLUX2(1)

	Do 20 I = 1, N
	If(KODE(I)) 20, 20, 10
10	CH = F1(I)
	F1(I) = DF1(I)
	DF1(I) = CH
20	Continue

	If(L.EQ.0) Go to 50
	Do 40 K = 1, L
	POT(K) = 0.
	FLUX1(K) = 0.
	FLUX2(k) = 0.

	Do 30 J = 1, N
	KK = J + 1

	Call EXTINPC(CX(K), CY(K), X(J), Y(J), X(KK), Y(KK), A, B
	1 , DQ1, DQ2, DU1, DU2, 1)
	POT(K) = POT(K) + DF1(J)*B - F1(J)*A
	FLUX1(K) = FLUX1(K) + DF1(J)*DU1 - F1(J)*DQ1
30	FLUX2(K) = FLUX2(K) + DF1(J)*DU2 - F1(J)*DQ2
	POT(K) = POT(K)/(2.*3.1415926)
	FLUX1(K) = FLUX1(K)/(2.*3.1415926)
40	FLUX2(K) = FLUX2(K)/(2.*3.1415926)
50	Return
	End



	!___________________________________________________________________



	Subroutine OUTPTPC(XM, YM, F1, DF1, CX, CY, POT, FLUX1, FLUX2)
	Implicit None

	Real*8 :: X(101), Y(101)
	Real*8 :: KODE(100)

	Real*8 :: G(100,100), H(100,100)
	Real*8 :: D, A, B, CH, DQ1, DQ2, DU1, DU2, XP, YP, X1, Y1, X2, Y2

	Integer :: NX, I, J, K, KK
	Integer :: N, L, INP, IPR

	Common N, L, INP, IPR
	Real*8 :: XM(1), YM(1), F1(1), DF1(1), CX(1), CY(1)
	Real*8 :: POT(1), FLUX1(1), FLUX2(1)

	Write(IPR, 100)
100	Format(' ', 79('*')//1X, 'Results'//2X, 'Boundary Nodes'//8X,
	1 'X', 15X,'Y', 13X, 'Potential', 3X, 'Potential Derivative'/)

	Do 10 I = 1, N
10	Write(IPR, 200) XM(I), YM(I), F1(I), DF1(I)
200	Format(4(2X, E14.5))

	If(L.EQ.0) Go to 30
	Write(IPR, 300)
300	Format(//, 2X, 'Internal Points', //8X, 'X', 15X, 'Y', 13X,
	1 'Potential', 9X, 'Flux X', 10X, 'Flux Y'/)

	Do 20 K = 1, L
20	Write(IPR, 400) CX(K), CY(K), POT(K), FLUX1(K), FLUX2(K)
400	Format(5(2X, E14.5))
30	Write(IPR, 500)
500	Format(' ',79('*'))
	Return
	End














