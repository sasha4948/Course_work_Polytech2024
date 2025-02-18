Subroutine B_CalcLapl(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
P,GradP,LaplP)
REAL X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),& ! центры и объемы ячеек
IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! центры и нормали I-граней
JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали J-граней
P(0:NI,0:NJ),GradP(0:NI,0:NJ,2),& ! P и массив для GradP
LaplP(0:NI,0:NJ) ! LaplP
real Sf(4,2),rf(4,2),ri(4,2),GPE(2),rE(2),Pf(4),rnc(2),Nf(2)
integer NeighCell(4,2)
real d,dn,dnc,dpdn,dpdn_c
integer i,j,inc,jnc,iface,icell
! инициализация LaplP
LaplP = 0
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
		
		! цикл по всем граням
		do iface = 1,4
			inc = NeighCell(iface,1)
			jnc = NeighCell(iface,2)
			d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! расстояние от центра ячейки до грани
			dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! расстояние до центра соседней ячейки
			Pf = RLinearInterp(d,dn,P(i,j),P(inc,jnc))
			! расст между 2 соседн ячейками (модуль вектотра)
			dnc = norm2(CellCenter(inc,jnc,:) - CellCenter(i,j,:))
			dpdn = (P(inc,jnc) - P(i,j))/dnc ! производная по направлению к соседней ячеке
			! единичный вектор к соседн ячейке (l_e)
			rnc(:) = (CellCenter(inc,jnc,:) - CellCenter(i,j,:))/dnc
			rE(1) = RLinearInterp(d,dn,CellCenter(i,j,1),CellCenter(inc,jnc,1))
			rE(2) = RLinearInterp(d,dn,CellCenter(i,j,2),CellCenter(inc,jnc,2))
			! комп градиента в точке E
			GPE(1) = RLinearInterp(d,dn,GradP(i,j,1),GradP(inc,jnc,1))
			GPE(2) = RLinearInterp(d,dn,GradP(i,j,2),GradP(inc,jnc,2))
			! нормированная внешняя нормаль к грани
			Nf(:) = Sf(iface,:)/norm2(Sf(iface,:))
			! для приграничной ячейки
			if (dn.lt.1e-7) then
				dpdn_c = dot_product(GradP(i,j,:),Nf(:)) ! (gradP)_P * n
				!dpdn = dpdn + (dpdn - dpdn_c) ! 1 - order
				dpdn = 5./3.*dpdn - 2./3.*dpdn_c ! 2 - order
				GPE(:) = GradP(i,j,:)
			end if
			! обновляем производную
			dpdn = dpdn + dot_product(GPE(:),Nf(:) - rnc(:))
			! обновляем LaplP
			LaplP(i,j) = LaplP(i,j) + dpdn*norm2(Sf(iface,:))
		end do
		LaplP(i,j) = LaplP(i,j)/CellVolume(i,j)
	end do
end do
End Subroutine