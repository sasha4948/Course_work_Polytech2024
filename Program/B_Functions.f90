! Функция для вычисления P
real function Pressure(X,Y)
!Pressure = 2*x + 2*y
!Pressure =2*x**2 + 4*y**2 + 7*y + 3*x
Pressure = x**3 + y**3 + 2*x + 3*y
!Pressure =x**4 + y**4 + y**2 + x**2
End Function

!Функция для вычисления точного значения GradP
Subroutine Calc_GradP_Exact(x,y,GradP)
real :: x, y, GradP(2)
!GradP(1) = 2.0
!GradP(2) = 2.0
!GradP(1) = 4*x + 3
!GradP(2) = 8*y + 7
GradP(1) = 3*x**2 + 2.0
GradP(2) = 3*y**2 + 3.0
!GradP(1) = 4*x**3 + 2.0*x
!GradP(2) = 4*y**3 + 2.0*y
End Subroutine

! Функция для вычисления поля скорости
subroutine Velocity(X,Y,V)
real :: x, y, V(2)
!V(1) = 1.0
!V(2) = 1.0
!V(1) = - 3*y
!V(2) = x
V(1) = -3*y**2 - x*y - y
V(2) = x**2 + x*y + x
!V(1) = -3*y**3 - x*y - y
!V(2) = x**3 + x*y + x
End subroutine

! Функция для вычисления точного значения divV
real function Calc_DivV_Exact(X,Y)
!Calc_DivV_Exact = 4.0
!Calc_DivV_Exact = 4*x + 8*y + 10.0
Calc_DivV_Exact = 3*x**2 + 3*y**2 + 5.0
End Function

! Функция для вычисления точного значения LaplP
real function Calc_LaplP_Exact(X,Y)
!Calc_LaplP_Exact = 12.0
Calc_LaplP_Exact = 6*x + 6*y
!Calc_LaplP_Exact = 12*x**2 + 12*y**2 + 4.0
End Function

! Функция для вычисления точного значения rotV
real function Calc_RotV_Exact(X,Y)
!Calc_RotV_Exact = 4.0
Calc_RotV_Exact = 3*x + 7*y + 2
!Calc_RotV_Exact = 3*x**2 + 9*y**2 + x + y + 2.0
End Function

! Линейная интерполяция: выч знач в (.) между 2 изв знач
! d1, d2 - расстояния от (.) интерп до 2 изв (.)
! p1, p2 - знач в 2 изв (.)
real function RLinearInterp(d1,d2,p1,p2)
real :: d1, d2, p1, p2
RLinearInterp = (d1*p2 + d2*p1)/(d1+d2)
End Function