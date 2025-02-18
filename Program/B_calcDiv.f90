Subroutine B_CalcDiv(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
V,DivV,P,GradP,scheme)
REAL X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),& ! центры и объемны ячеек
IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! центры и нормали I-граней
JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали J-граней
V(0:NI,0:NJ,2),DivV(0:NI,0:NJ),& ! V и массив divV
P(0:NI,0:NJ),GradP(0:NI,0:NJ,2) ! P и массив GradP
real Sf(4,2),rf(4,2),Vf(2),A_ls(2,2),b_ls(2),ri(4,2),GPE(2),rE(2),Pf(4)
integer NeighCell(4,2)
real d,dn,det_A
integer i,j,inc,jnc,iface,icell,scheme
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
		
		! инициализация divV
		DivV(i,j) = 0
		
		! цикл по всем граням
		do iface = 1,4
		inc = NeighCell(iface,1)
		jnc = NeighCell(iface,2)
		d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! расстояние от центра ячейки до грани
		dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! расстояние до центра соседней ячейки
		
		! интерполяция скорости на грани
		Vf(1) = RLinearInterp(d,dn,V(i,j,1),V(inc,jnc,1))
		Vf(2) = RLinearInterp(d,dn,V(i,j,2),V(inc,jnc,2))
		
		select case(scheme)
		!DivV - простая дивергенция
		case(0)
		Pf(iface) = 1
		!Central - центральная схема (линейная интерполяция)
		case(1)
		Pf(iface) = RLinearInterp(d,dn,P(i,j),P(inc,jnc))
		!FOU - противопоточная схема первого порядка
		case(2)
		! если G (расход) > 0 - P в тек ячейке (слева от грани)
		if (dot_product(Vf(:),Sf(iface,:)) > 0.0) then
			Pf(iface) = P(i,j)
		else
			Pf(iface) = P(inc,jnc)
			! для приграничной ячейки
			if (dn < 1e-7) then
				Pf(iface) = Pf(iface) + P(inc,jnc) - P(i,j)
			end if
		end if
		
		!SOU - противопоточная схема второго порядка
		case(3)
		! если G (расход) > 0 - P в тек ячейке (слева от грани) с учетом поправки
		if (dot_product(Vf(:),Sf(iface,:)) > 0.0) then
			Pf(iface) = P(i,j) + dot_product(GradP(i,j,:),rf(iface,:) - CellCenter(i,j,:))
		else
			Pf(iface) = P(inc,jnc) + dot_product(GradP(inc,jnc,:),rf(iface,:) - &
			CellCenter(inc,jnc,:))
			! учет приграничных ячеек
			if (dn < 1e-7) then
				Pf(iface) = 2*P(inc,jnc) - P(i,j) -&
				4*P(inc,jnc) + 4*P(i,j) +&
				3*dot_product(GradP(i,j,:),rf(iface,:)-CellCenter(i,j,:))
				!Pf(iface) = P(i,j) + dot_product(GradP(i,j,:),rf(iface,:) - CellCenter(i,j,:))
			end if
		end if
		end select
		! обновление divV
		DivV(i,j) = DivV(i,j) + dot_product(Pf(iface)*Vf(:),Sf(iface,:))
		end do
		DivV(i,j) = DivV(i,j)/CellVolume(i,j)
	end do
end do
End Subroutine