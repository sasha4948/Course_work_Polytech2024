Subroutine B_CalcRot(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
V,RotV)
REAL X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),& ! центры и объемны ячеек
IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! центры и нормали I-граней
JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали J-граней
V(0:NI,0:NJ,2),RotV(0:NI,0:NJ) ! V и массив для rotV
real Sf(4,2),rf(4,2),Vf(2),GPE(2),rE(2)
integer NeighCell(4,2)
real d,dn
integer i,j,inc,jnc,iface,icell
do i = 1,NI-1
	do j = 1,NJ-1
		! вектор прощади грани с учетом внешней нормали
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
		! инициализация rotV
		RotV(i,j) = 0
		do iface = 1,4
			inc = NeighCell(iface,1)
			jnc = NeighCell(iface,2)
			d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! расстояние от центра ячейки до грани
			dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! расстояние до центра соседней ячейки
			! интерполяция скорости на грани
			Vf(1) = RLinearInterp(d,dn,V(i,j,1),V(inc,jnc,1))
			Vf(2) = RLinearInterp(d,dn,V(i,j,2),V(inc,jnc,2))
			! обновление rotV
			RotV(i,j) = RotV(i,j) + (Sf(iface,1)*Vf(2) - Sf(iface,2)*Vf(1))
		end do
		RotV(i,j) = RotV(i,j)/CellVolume(i,j)
	end do
end do
End Subroutine