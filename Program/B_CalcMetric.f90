Subroutine B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
REAL X(NI,NJ),Y(NI,NJ),& !координаты узлов сетки
CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),& !центры ячеек и объемы ячеек
IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! центры и векторы граней для I-граней
JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2) ! центры и векторы граней для J-граней
REAL r(2) ! вектор с компонентами r(1) - X, r(2) - Y

!=== FACE CENTERS AND FACE VECTORS ===
! I-DIRECTION (расчет для граней в направлении I)
DO J = 1,NJ-1
	DO I = 1,NI
		r(1) = X(I,J+1) - X(I,J) ! r = вектор от одного узла к другому
		r(2) = Y(I,J+1) - Y(I,J)
		IFaceVector(I,J,1) = r(2) ! IFaceVector = r поверн-го на 90 градусов
		IFaceVector(I,J,2) =-r(1) ! IFaceVector направлен на увеличение I-индекса
		IFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i,j+1)) ! центры грани для I-граней
		IFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i,j+1))
	ENDDO
ENDDO

! J-DIRECTION (расчет для граней в направлении J)
DO J = 1,NJ
	DO I = 1,NI-1
		r(1) = X(I+1,J) - X(I,J) ! r = вектор от одного узла к другому
		r(2) = Y(I+1,J) - Y(I,J)
		JFaceVector(I,J,1) =-r(2) ! JFaceVector = r поверн-го на -90 градусов
		JFaceVector(I,J,2) = r(1) ! JFaceVector направлен на увеличение J-индекса
		JFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i+1,j))
		JFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i+1,j))
	ENDDO
ENDDO

!=== CELL VOLUMES ===
! Расчет объемов ячеек
DO J = 1,NJ-1
	DO I = 1,NI-1
		! диагональный вектор
		r(1)=X(I+1,J+1) - X(I,J)
		r(2)=Y(I+1,J+1) - Y(I,J)
		! сумма объемов 2 треуг, раздел-х диаг-ю
		CellVolume(I,J) = 0.5*DOT_PRODUCT(IFaceVector(I,J,:),r)& !скалярно площадь грани x диаг
		+ 0.5*DOT_PRODUCT(JFaceVector(I,J,:),r)
	ENDDO
ENDDO

!=== CELL CENTERS ===
! FOR INNER CELLS: CENTER OF CONTOUR (sum of FaceCenter*FaceLength/Perimeter)
DO J = 1,NJ-1
	DO I = 1,NI-1
		! центр ячейки для внутренних ячеек
		CellCenter(I,J,:) = ( IFaceCenter(I ,J,:)*Norm2(IFaceVector(I ,J,:))+&
		IFaceCenter(I+1,J,:)*Norm2(IFaceVector(I+1,J,:))+&
		JFaceCenter(I,J ,:)*Norm2(JFaceVector(I,J ,:))+&
		JFaceCenter(I,J+1,:)*Norm2(JFaceVector(I,J+1,:)) )&
		/( Norm2(IFaceVector(I,J,:))+Norm2(IFaceVector(I+1,J,:))+&
		Norm2(JFaceVector(I,J,:))+Norm2(JFaceVector(I,J+1,:)) )
	ENDDO
ENDDO

! FOR DUMMY CELLS ON BOUNDARIES: CELL CENTER = FACE CENTER
! I-BOUNDARIES -----------------------------------------------------

DO NBOUND = 1,2

	IF (NBOUND.EQ.1) THEN
		IBOUND = 1; IOUT = 0
	ELSE
		IBOUND = NI; IOUT = NI
	ENDIF
	
	DO J = 1,NJ-1
		! граничная ячейка: центр ячейки = центр соотв грани
		CellCenter(IOUT,J,:) = IFaceCenter(IBOUND,J,:)
	ENDDO
	
ENDDO

! J-BOUNDARIES -----------------------------------------------------
DO NBOUND = 1,2

	IF (NBOUND.EQ.1) THEN
		JBOUND = 1; JOUT = 0
	ELSE
		JBOUND = NJ; JOUT = NJ
	ENDIF
	
	DO I = 1,NI-1
		CellCenter(I,JOUT,:) = JFaceCenter(I,JBOUND,:)
	ENDDO
	
ENDDO

END SUBROUTINE

