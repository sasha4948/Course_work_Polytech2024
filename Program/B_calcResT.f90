Subroutine B_CalcResT(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
Reyn,Pr,VFlos,T,GradT,ResT,scheme_res,CFL,VNM,a_diff,dtau)
REAL X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),& ! центры и объемы ячеек
IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! центры и нормали I-граней
JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали J-граней
T(0:NI,0:NJ),GradT(0:NI,0:NJ,2),& ! T и массив для GradT
VFlos(0:NI,0:NJ,2),& ! Поле V из Flos
ResT(0:NI,0:NJ),dtau(0:NI,0:NJ) ! ResT и шаг по псевдоврем dtau
real Sf(4,2),rf(4,2),ri(4,2),GTE(2),rE(2),Tf,rnc(2),Nf(2),Vf(2), VTf(2)
integer NeighCell(4,2)
real d,dn,dnc,dTdn,dTdn_c,Reyn,Pr,CFL,VNM,a_diff,dtau_c,dtau_d
integer i,j,inc,jnc,iface,icell,scheme_res
!$OMP DO private(i,j,Sf,rf,NeighCell,dtau_c,dtau_d,iface,inc,jnc,d,dn,Vf,Tf,VTf,&
!$OMP& dnc,dTdn,dTdn_c,rnc,rE,GTE,Nf)
do j = 1,NJ-1
	do i = 1,NI-1
		! инициализация невязки ResT
		ResT(i,j) = 0.0
		
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
		
		! инициализация шагов по псевдоврем
		dtau_c = 0.0 ! конвективный
		dtau_d = 0.0 ! диффузионный
		dtau(i,j) = 0.0
		
		! цикл по всем граням
		do iface = 1,4
		inc = NeighCell(iface,1)
		jnc = NeighCell(iface,2)
		d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! расстояние от центра ячейки до грани
		dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! расстояние до центра соседней ячейки
		
		!расчет конвективного слагаемого
		Vf(1) = RLinearInterp(d,dn,VFlos(i,j,1),VFlos(inc,jnc,1)) ! линейная интерполяция скорости на грани
		Vf(2) = RLinearInterp(d,dn,VFlos(i,j,2),VFlos(inc,jnc,2))
		
		! Central - центральная схема расчета конв слагаемого
		Tf = RLinearInterp(d,dn,T(i,j),T(inc,jnc))
		
		select case(scheme_res)
		! случай отсутсвия конвективного слагаемого
		case(0)
		Vf(:) = 0.0
		
		!Upwind gradient based
		case(2)
		
		! если G (расход) > 0 - на тек (L) ячейке
		if (dot_product(Vf(:),Sf(iface,:)) > 0.0) then
			 Tf = 2*(T(i,j) + dot_product(GradT(i,j,:),rf(iface,:) - CellCenter(i,j,:))) - Tf
		else
			Tf = 2*(T(inc,jnc) + dot_product(GradT(inc,jnc,:),rf(iface,:) - CellCenter(inc,jnc,:))) - Tf
			! для приграничной ячейки
			if (dn < 1e-7) then
				Tf = 0.5*RLinearInterp(d,dn,T(i,j),T(inc,jnc)) + &
				0.5*(2*T(inc,jnc) - T(i,j) - 4*T(inc,jnc) + 4*T(i,j) +&
				3*dot_product(GradT(i,j,:),rf(iface,:)-CellCenter(i,j,:)))
			end if
		end if
		
		end select
		! итоговое конвективное слагаемое TV
		VTf(:) = Vf(:)*Tf
		! расчет диффузионного слагаемого (по аналогии с lapl)
		! расст между 2 соседн ячейками (модуль вект)
		dnc = norm2(CellCenter(inc,jnc,:) - CellCenter(i,j,:))
		dTdn = (T(inc,jnc) - T(i,j))/dnc ! произв по напр к соседн
		! единичный вектор к соседн ячейке (l_e)
		rnc = (CellCenter(inc,jnc,:) - CellCenter(i,j,:))/dnc
		rE(1) = RLinearInterp(d,dn,CellCenter(i,j,1),CellCenter(inc,jnc,1))
		rE(2) = RLinearInterp(d,dn,CellCenter(i,j,2),CellCenter(inc,jnc,2))
		! комп градиента в (.) E
		GTE(1) = RLinearInterp(d,dn,GradT(i,j,1),GradT(inc,jnc,1))
		GTE(2) = RLinearInterp(d,dn,GradT(i,j,2),GradT(inc,jnc,2))
		! нормир-я внешняя нормаль к грани
		Nf = Sf(iface,:)/norm2(Sf(iface,:))
		! для приграничной ячейки
		if (dn.lt.1e-7) then
		dTdn_c = dot_product(GradT(i,j,:),Nf(:))
		!dTdn = dTdn + (dTdn - dTdn_c) ! 1 - order
		dTdn = 5./3.*dTdn - 2./3.*dTdn_c ! 2 - order
		GTE(:) = GradT(i,j,:)
		end if
		! обновляем производную
		dTdn = dTdn + dot_product(GTE(:),Nf(:) - rnc(:))
		! граничное условие (адиабатические стенки)
		if ((inc == 0) .or. (inc == NI)) then
		dTdn = 0.0
		T(inc,jnc) = T(i,j) + dnc*3./5.*(2./3.*dot_product(GradT(i,j,:),Nf(:)) -&
		dot_product(GradT(i,j,:),Nf(:) - rnc(:)))
		end if
		! обновление невязки с учетом конвктивного и диффузионного слагаемого
		ResT(i,j) = ResT(i,j) + dot_product(VTf(:),Sf(iface,:)) - &
		dTdn*norm2(Sf(iface,:))/(Reyn*Pr)
		! обновление шагов по псевдовремени
		dtau_c = dtau_c + abs(dot_product(Vf(:),Sf(iface,:)))
		dtau_d = dtau_d + a_diff/norm2(Sf(iface,:))
		end do
		ResT(i,j) = ResT(i,j)/CellVolume(i,j)
		! обновление общего шага по псевдо времени
		! учитываем, дает ли вклад, чтобы малая величина не делилась
		if (dtau_c > 1e-7) then
		dtau(i,j) = dtau(i,j) + dtau_c/(CFL*CellVolume(i,j))
		end if
		if (dtau_d > 1e-7) then
		dtau(i,j) = dtau(i,j) + 2*dtau_d/VNM
		end if
		dtau(i,j) = 1.0/dtau(i,j)
	end do
end do
!$OMP END DO
End Subroutine