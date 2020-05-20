 SUBROUTINE PLAG_RFLU_FindCellsHardCode(pRegion)
 IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: icg,iCol,iDep,iPcl,iRow
  REAL(RFREAL) :: drad,dthe,dz,radLoc,theLoc,xLoc,yLoc,zLoc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsHardCode',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

 DO iPcl = 1,pPlag%nPcls

! ==============================================================================  
!   Get particle position
! ==============================================================================  

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)

    radLoc = SQRT(xLoc**2.0_RFREAL + yLoc**2.0_RFREAL)
    theLoc = ATAN2(yLoc,xLoc)
    drad = global%drad
    dthe = global%dthe
    dz   = global%dz
    iRow = FLOOR((radLoc - (global%radMin-0.5_RFREAL*drad))/drad + 1)
    iCol = FLOOR((theLoc - (global%theMin-0.5_RFREAL*dthe))/dthe + 1)
    iDep = FLOOR((radLoc - (global%zMin-0.5_RFREAL*dz))/dz + 1)
    icg  = (iRow-1)*(global%Imax*global%Jmax) + (iDep-1)*(global%Imax) + iCol

    pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg

 END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************
 
  CALL DeregisterFunction(global)

 END SUBROUTINE PLAG_RFLU_FindCellsHardCode
