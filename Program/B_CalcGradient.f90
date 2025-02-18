Subroutine B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
IGrad,P,GradP)
REAL X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),& ! центры и объемы ячеек
IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! центры и нормали I-граней
JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали J-граней
P(0:NI,0:NJ),& ! поле давления
GradP_Exact(0:NI,0:NJ,2), GradP_Error(0:NI,0:NJ,2), & !точное и ошибка для GradP
GradP(0:NI,0:NJ,2),GradP_tmp(0:NI,0:NJ,2) ! GradP
integer IGrad !метод расчета
real Sf(4,2),rf(4,2),Pf(4),GPE(2),rE(2)
integer NeighCell(4,2),GrG_iter
real d,dn, p_dum
integer i,j,inc,jnc,iface
select case (IGrad)

!метод Грина-Гаусса
case(1)
!$OMP DO private(i,j,Sf,rf,NeighCell,iface,inc,jnc,d,dn,Pf)
do i = 1,NI-1
	do j = 1,NJ-1
		! вектор площади грани с учетом внешней нормали
		Sf(1,:) = -IFaceVector(i,j,:)
		Sf(2,:) = IFaceVector(i+1,j,:)
		Sf(3,:) = -JFaceVector(i,j,:)
		Sf(4,:) = JFaceVector(i,j+1,:)
		
		! координаты центров граней
		rf(1,:) = IFaceCenter(i,j,:)
		rf(2,:) = IFaceCenter(i+1,j,:)
		rf(3,:) = JFaceCenter(i,j,:)
		rf(4,:) = JFaceCenter(i,j+1,:)
		
		! соседние ячейки
		NeighCell(1,:) = [i-1,j]
		NeighCell(2,:) = [i+1,j]
		NeighCell(3,:) = [i,j-1]
		NeighCell(4,:) = [i,j+1]
		
		! инициализация GradP
		GradP(i,j,:) = 0
		
		! цикл по всем граням
		do iface = 1,4
			inc = NeighCell(iface,1)
			jnc = NeighCell(iface,2)
			d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! расстояние от центра ячейки до грани
			dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! расст до центра соседней ячейк
			Pf(iface) = RLinearInterp(d,dn,P(i,j),P(inc,jnc)) ! интерполяция P на грань
			
			!расчет в приграничной ячейке
			if (dn < 1e-7) then
				! значение в центре заграничной ячейке
				p_dum = Pressure(2*rf(iface,1) - CellCenter(i,j,1),2*rf(iface,2) - CellCenter(i,j,2))
				Pf(iface) = 0.5*(p_dum + P(i,j))
			end if
			
			! обновляем GradP
			GradP(i,j,:) = GradP(i,j,:) + Pf(iface)*Sf(iface,:)
		end do
		GradP(i,j,:) = GradP(i,j,:)/CellVolume(i,j)
		! вычисление точного значения GradP
		call Calc_GradP_Exact(CellCenter(I,J,1),CellCenter(i,j,2),GradP_Exact(I,J,:))
	end do
end do
!$OMP END DO

!расчет ошибки GradP
GradP_Error = ABS((GradP_Exact-GradP)/GradP_Exact)
write(*,*) '1', maxval(GradP_Error(1:NI-1,1:NJ-1,:))

!метод Грина-Гаусса с итерациями
case(2)
GrG_iter = 10 !кол-во итераций

do k = 1,GrG_iter

	!$OMP DO private(i,j,Sf,rf,NeighCell,iface,inc,jnc,d,dn,PE,rE,GPE,Pf)
	do i = 1,NI-1
		do j = 1,NJ-1
			! инициализация GradP
			GradP_tmp(i,j,:) = 0
			! вектор площади грани с учетом внешней нормали
			Sf(1,:) = -IFaceVector(i,j,:)
			Sf(2,:) = IFaceVector(i+1,j,:)
			Sf(3,:) = -JFaceVector(i,j,:)
			Sf(4,:) = JFaceVector(i,j+1,:)
			
			! координаты центров граней
			rf(1,:) = IFaceCenter(i,j,:)
			rf(2,:) = IFaceCenter(i+1,j,:)
			rf(3,:) = JFaceCenter(i,j,:)
			rf(4,:) = JFaceCenter(i,j+1,:)
			
			! соседние ячейки
			NeighCell(1,:) = [i-1,j]
			NeighCell(2,:) = [i+1,j]
			NeighCell(3,:) = [i,j-1]
			NeighCell(4,:) = [i,j+1]
			
			do iface = 1,4
				inc = NeighCell(iface,1)
				jnc = NeighCell(iface,2)
				d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! расстояние от центра яч до грани
				dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! расст до соседней ячейки
				!Pf(iface) = RLinearInterp(d,dn,P(i,j),P(inc,jnc))
				PE = RLinearInterp(d,dn,P(i,j),P(inc,jnc)) ! интерп P для получ знач в (.) E
				!if (dn < 1e-7) then
				! p_dum = Pressure(2*rf(iface,1) - CellCenter(i,j,1),2*rf(iface,2) - CellCenter(i,j,2))
				! PE = 0.5*(p_dum + P(i,j))
				!end if
				
				! координаты точки E
				rE(1) = RLinearInterp(d,dn,CellCenter(i,j,1),CellCenter(inc,jnc,1))
				rE(2) = RLinearInterp(d,dn,CellCenter(i,j,2),CellCenter(inc,jnc,2))
				
				! компоненты градиента в точке E
				GPE(1) = RLinearInterp(d,dn,GradP(i,j,1),GradP(inc,jnc,1))
				GPE(2) = RLinearInterp(d,dn,GradP(i,j,2),GradP(inc,jnc,2))
				
				! вычисление P на грани c поправки
				Pf(iface) = PE + dot_product(GPE(:),rf(iface,:) - rE(:))
				
				! обновление GradP с учетом итерации
				GradP_tmp(i,j,:) = GradP_tmp(i,j,:) + Pf(iface)*Sf(iface,:)
				
			end do
			GradP_tmp(i,j,:) = GradP_tmp(i,j,:)/CellVolume(i,j)
			GradP(i,j,:) = GradP_tmp(i,j,:)
			call Calc_GradP_Exact(CellCenter(I,J,1),CellCenter(i,j,2),GradP_Exact(I,J,:))
		end do
	end do
	!$OMP END DO
	! вычисление ошибки GradP	
	!write(*,*) k, GradP_Exact((NI-1)/2,(NJ-1)/2,:)
	GradP_Error = ABS((GradP_Exact-GradP)/GradP_Exact)
	!write(*,*) k, maxval(GradP_Error(1:NI-1,1:NJ-1,:))
end do
end select
End Subroutine