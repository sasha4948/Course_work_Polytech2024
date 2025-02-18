Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradP_Error,V,DivV,DivV_Exact,DivV_Error, &
LaplP,LaplP_Exact,LaplP_Error,RotV,RotV_Exact,RotV_Error,VFlos,T,GradT,ResT,TFlos,T_Err)
Real,Dimension(NI,NJ):: X,Y
Real,Dimension(0:NI,0:NJ)::P,DivV, DivV_Exact, DivV_Error,LaplP,LaplP_Exact,LaplP_Error,&
RotV,RotV_Exact,RotV_Error,T,ResT,TFlos,T_Err
Real,Dimension(0:NI,0:NJ,2)::GradP,GradP_Error,V,GradT,VFlos
Write(IO,*) 'VARIABLES = "X","Y","P","GradP_X","GradP_Y","GradP_Error_X","GradP_Error_Y","V_x","V_y", &
"DivV","DivV_Exact","DivV_Error","LaplP","LaplP_Exact","LaplP_Error","RotV","RotV_Exact","RotV_Error",&
"VFlos_X","VFlos_Y","T","GradT_X","GradT_Y","ResT","TFlos","T_Error"'
Write(IO,*) 'ZONE I=',NI,',J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
Write(IO,'(100F14.7)') X(1:NI,1:NJ)
Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1)
Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2)
Write(IO,'(100F14.7)') GradP_Error(1:NI-1,1:NJ-1,1)
Write(IO,'(100F14.7)') GradP_Error(1:NI-1,1:NJ-1,2)
Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)
Write(IO,'(100F14.7)') DivV(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') DivV_Exact(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') DivV_Error(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') LaplP(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') LaplP_Exact(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') LaplP_Error(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') RotV(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') RotV_Exact(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') RotV_Error(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') VFlos(1:NI-1,1:NJ-1,1)
Write(IO,'(100F14.7)') VFlos(1:NI-1,1:NJ-1,2)
Write(IO,'(100F14.7)') T(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') GradT(1:NI-1,1:NJ-1,1)
Write(IO,'(100F14.7)') GradT(1:NI-1,1:NJ-1,2)
Write(IO,'(100F14.7)') ResT(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') TFlos(1:NI-1,1:NJ-1)
Write(IO,'(100F14.7)') T_Err(1:NI-1,1:NJ-1)
End Subroutine