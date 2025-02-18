Program Main
USE OMP_LIB
character(*), parameter:: InputFile='input.txt',OutputFile='data.plt' !файлы на ввод и вывод
character MeshFile*30 !файл с сеткой
integer, parameter:: IO = 12 
integer :: IGrad,scheme,iter,niter,i,j,NI,NJ,cavity, scheme_res
real,allocatable,dimension(:,:):: X,Y,P,CellVolume ! скалярные величины
real,allocatable,dimension(:,:,:):: CellCenter, IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! геометрические характеристики
real,allocatable,dimension(:,:,:):: GradP, GradP_Exact, GradP_Error ! градиент давления и ошибка
real,allocatable,dimension(:,:):: DivV, DivV_Exact, DivV_Error ! дивергенция скорости
real,allocatable,dimension(:,:,:):: V !поле скорости
real,allocatable,dimension(:,:):: LaplP, LaplP_Exact, LaplP_Error !лапласиан давления
real,allocatable,dimension(:,:):: RotV, RotV_Exact, RotV_Error !ротор скорости
real,allocatable,dimension(:,:,:):: VFlos, GradT !V и T из Flos
real,allocatable,dimension(:,:):: T, ResT, TFlos, T_Err, dtau ! Поля температуры
real :: Reyn, Pr, Vs, Ls, CFL, VNM, a_diff, t1,t2 !параметры для задачи о каверне

!=== READ INPUT FILE ===
WRITE(*,*) 'Чтение входного файла: ', InputFile
OPEN(IO,FILE=InputFile)
READ(IO,*) MeshFile ! Расчетная сетка
READ(IO,*) IGrad ! Выбор метода для градиента
READ(IO,*) scheme ! схема для divV
READ(IO,*) scheme_res ! схема для resT
READ(IO,*) Vs ! масштаб скорости
READ(IO,*) Ls ! масштаб длины
READ(IO,*) Reyn ! Re
READ(IO,*) Pr ! Pr
READ(IO,*) CFL ! CFL
READ(IO,*) VNM ! VNM
READ(IO,*) niter ! число итераций
READ(IO,*) cavity ! тип задачи: 0 - операторы; 1 - каверна
CLOSE(IO)

!=== READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
WRITE(*,*) 'Чтение размеров сетки из файла: ', MeshFile
OPEN(IO,FILE = MeshFile)
READ(IO,*) NI,NJ
WRITE(*,*) 'NI, NJ = ',NI,NJ

!=== ALLOCATE ALL ARRAYS ===
WRITE(*,*) 'Выделение памяти под массивы'
allocate(X(NI,NJ)) ! координаты X узлов сетки
allocate(Y(NI,NJ)) ! координаты Y узлов сетки
allocate(P(0:NI,0:NJ)) ! давление
allocate(CellVolume(NI-1,NJ-1)) ! объем ячеек
allocate(CellCenter(0:NI,0:NJ,2)) ! центры ячеек
allocate(IFaceCenter( NI,NJ-1,2)) ! центры граней для I-граней
allocate(IFaceVector( NI,NJ-1,2)) ! вектора для I-граней
allocate(JFaceCenter( NI-1,NJ,2)) ! центры граней для J-граней
allocate(JFaceVector( NI-1,NJ,2)) ! вектора для J-граней
allocate(GradP(0:NI,0:NJ,2)) ! GradP
allocate(GradP_Exact(0:NI,0:NJ,2)) ! точный GradP
allocate(GradP_Error(0:NI,0:NJ,2)) ! ошибка GradP
allocate(DivV(0:NI,0:NJ)) ! дивергенция скорости
allocate(DivV_Exact(0:NI,0:NJ)) ! точная дивергениция
allocate(DivV_Error(0:NI,0:NJ)) ! ошибка дивергениции
allocate(LaplP(0:NI,0:NJ)) ! Лапласиан давления
allocate(LaplP_Exact(0:NI,0:NJ)) ! точный Лапласиан
allocate(LaplP_Error(0:NI,0:NJ)) ! ошибка Лапласиана
allocate(RotV(0:NI,0:NJ)) ! ротор
allocate(RotV_Exact(0:NI,0:NJ)) ! точный ротор
allocate(RotV_Error(0:NI,0:NJ)) ! ошибка ротора
allocate(GradT(0:NI,0:NJ,2)) ! градиент температуры
allocate(VFlos(0:NI,0:NJ,2)) ! Поле скорости из Flos
allocate(T(0:NI,0:NJ)) ! температура
allocate(ResT(0:NI,0:NJ)) ! невязка T
allocate(TFlos(0:NI,0:NJ)) ! T из Flos
allocate(T_Err(0:NI,0:NJ)) ! ошибка T
allocate(dtau(0:NI,0:NJ)) ! шаг по времени
allocate(V(0:NI,0:NJ,2)) ! Поле скорости

!=== READ GRID ===
WRITE(*,*) 'Чтение сетки из файла: ', MeshFile
READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
CLOSE(IO)

!Случай расчета задачи о каверне
if (cavity == 1) then

!=== READ FLOS VELOCITY FIELD ===
OPEN(IO,FILE = 'VelocityField.txt')
READ(IO,*) ((VFlos(I,J,1),VFlos(I,J,2),I=0,NI),J=0,NJ)
CLOSE(IO)

!=== READ FLOS TEMPERATURE FIELD ===
OPEN(IO,FILE = 'Temperature.txt')
READ(IO,*) ((TFlos(I,J),I=0,NI),J=0,NJ)
CLOSE(IO)
end if

!=== CALCULATE METRIC ===
WRITE(*,*) 'Вычисление метрик'
Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

!=== INITIATE FIELDS ===
WRITE(*,*) 'Инициализация полей'
DO J = 0,NJ
	DO I = 0,NI
		P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(i,j,2))
		call Calc_GradP_Exact(CellCenter(I,J,1),CellCenter(i,j,2),GradP_Exact(I,J,:))
		call Velocity(CellCenter(I,J,1),CellCenter(i,j,2),V(I,J,:))
		DivV_Exact(I,J) = Calc_DivV_Exact(CellCenter(I,J,1),CellCenter(i,j,2))
		LaplP_Exact(I,J) = Calc_LaplP_Exact(CellCenter(I,J,1),CellCenter(i,j,2))
		RotV_Exact(I,J) = Calc_RotV_Exact(CellCenter(I,J,1),CellCenter(i,j,2))
	ENDDO
ENDDO

GradP = 0

!=== CALCULATE GRADIENT ===
WRITE(*,*) 'Вычисление операторов'
Call B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
IGrad,P,GradP)
GradP_Error = ABS((GradP_Exact-GradP)/GradP_Exact)
write(*,*) 'Max error GradP:', maxval(GradP_Error(1:NI-1,1:NJ-1,:))
write(*,*) 'Error (0.5, 0.5)', GradP_Error((NI-1)/2,(NJ-1)/2,:) 

!=== CALCULATE DIVERGENCE ===
Call B_CalcDiv(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
V,DivV,P,GradP,scheme)
DivV_Error = ABS((DivV_Exact-DivV)/DivV_Exact)
write(*,*) 'Max error divV:', maxval(DivV_Error(1:NI-1,1:NJ-1))

!=== CALCULATE LAPLACIAN ===
Call B_CalcLapl(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
P,GradP,LaplP)
LaplP_Error = ABS((LaplP_Exact-LaplP)/LaplP_Exact)
write(*,*) 'Max error LaplP:', maxval(LaplP_Error(1:NI-1,1:NJ-1))

!=== CALCULATE ROTOR ===
Call B_CalcRot(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
V,RotV)
RotV_Error = ABS((RotV_Exact-RotV)/RotV_Exact)
write(*,*) 'Max error rotV:', maxval(RotV_Error(1:NI-1,1:NJ-1))

!=== CALCULATE RESIDUAL T ===

! инициализация поля
T(:,:) = 1.0
! стенки при постоянной T
T(:,0) = 1.0 !низ
T(:,NJ) = 2.0 !верх
dtau = 0.01 !псевдо время
a_diff = Vs*Ls/(Reyn*Pr) !коэффициент температуропроводности

OPEN(IO,FILE = 'Residual.plt')
write(io,*) 'Variables = "iterations", "Res T"'

t1 = OMP_GET_WTIME()
!$OMP parallel
PRINT *, "THREAD_NUM: ", OMP_GET_THREAD_NUM()
do iter = 1,niter
	! вычисление GradT
	Call B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,&
	JFaceCenter,JFaceVector,IGrad,T,GradT)
	! вычисление невязки для GradT
	Call B_CalcResT(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,&
	JFaceCenter,JFaceVector,Reyn,Pr,VFlos,T,GradT,ResT,scheme_res,&
	CFL,VNM,a_diff,dtau)
	!$OMP single
	!write(*,*) iter, maxval(abs(ResT(1:NI-1,1:NJ-1)))
	write(io,*) iter, maxval(abs(ResT(1:NI-1,1:NJ-1)))
	!$OMP end single
	!$OMP DO private(i,j)
	! итерационно обновляем T c учетом Res и dtau
	do j = 1, NJ-1
		do i = 1, NI-1
		T(i,j) = T(i,j) - ResT(i,j)*dtau(i,j)
		end do
	end do
	!$OMP END DO
end do
!$OMP end parallel
t2 = OMP_GET_WTIME()
print *, 'time = ',t2 - t1
close(io)
!вычисление ошибки T

if (cavity == 1) then
T_Err = ABS((TFlos-T)/TFlos)
end if

write(*,*) 'Max error T:',maxval(abs(T_Err(0:NI,0:NJ)))

!=== OUTPUT FIELDS ===
WRITE(*,*) 'Запись полей в файл: ', OutputFile
Open(IO,FILE=OutputFile)
Call B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradP_Error,V,DivV,DivV_Exact,DivV_Error,&
LaplP,LaplP_Exact,LaplP_Error,RotV,RotV_Exact,RotV_Error,VFlos,T,GradT,ResT,TFlos,T_Err)
Close(IO)
END PROGRAM Main