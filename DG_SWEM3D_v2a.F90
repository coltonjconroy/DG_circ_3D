!*************************************************************************************
! Declaration of Precision of Variables
!*************************************************************************************

MODULE precisions

IMPLICIT NONE
SAVE

INTEGER, PARAMETER :: int_p  = SELECTED_INT_KIND(16)
INTEGER, PARAMETER :: real_p = SELECTED_REAL_KIND(16)

END MODULE precisions
!*************************************************************************************
MODULE global
USE precisions
IMPLICIT NONE

!*************************************************************************************
!....Global variables
!*************************************************************************************

INTEGER(int_p) :: i, j, ii, kk, NDOF, nelemsX, nelemsZ, nnodesX, nnodesZ, hx, hy, hz, p, pu, plot_n
INTEGER(int_p) :: nelemsY, nnodesY, HEXelems, HEXnnodes, L2_DEG_3D, irkW, LNrk, NDOF_2D, ielems
INTEGER(int_p) :: nelemsXF, nnodesXF, DIM, irk, LE, coupled, Qelems, Qnnodes, NPTS_2D, L2_2D
INTEGER(int_p) :: P_method, NDOF_3D, NPTS_3D, L2_3D, meshX, meshY, meshZ, L2TP, LRK, TP_2D, TP_3D
INTEGER(int_p) :: nrkW, KERNEL, ielems_2D, stencil_2D, Z_level, var, NPTS_K, NPTS_K2D, pw, NDOF_3Dw
INTEGER(int_p) :: FULL2D, w_method, ph, NCONST, NEfacesX, NEfacesY, NEedgesX, NEedgesY, IRAMP, NRAMP
INTEGER(int_p), PARAMETER :: debug = 15, optimal = 1, star = 1, VpH = 1, L2int = 10, NON_LINEAR = 0
INTEGER(int_p), PARAMETER :: p2 = 1, NbPTS = 2, VID = 1, HOT_START = 0, SHIFTint = 6
INTEGER(int_p), PARAMETER :: PROB = 0, NUM_ERROR = 1, BASIS_OUTPUT = 0, QUAD_TYPE = 1, FULL3D = 1
INTEGER(int_p), ALLOCATABLE :: HEXzxy(:,:), HEXsigma(:,:), HEXgrid(:,:,:),  zrkCONN(:,:),HEXconn(:,:)
INTEGER(int_p), ALLOCATABLE :: QedgeX(:,:), QedgeY(:,:), HEXfaceX(:,:), HEXfaceY(:,:), yEDGE(:,:)
INTEGER(int_p), ALLOCATABLE :: zrkzCONN(:), RK_HEXconn(:,:), CONNzx(:,:), CONNgrid(:,:),CONNsigma(:,:)
INTEGER(int_p), ALLOCATABLE :: zEDGE(:,:), xEDGE(:,:), xFACE(:,:), yFACE(:,:), zFACE(:,:), Qstencil(:,:)
INTEGER(int_p), ALLOCATABLE :: CONNielem(:,:), KCelems(:), HEXfaceZ(:,:), topLAYER(:),ExEDGE(:), EyEDGE(:)
REAL(real_p) :: L2_u, L2_zeta, Linf_ubar, Linf_gse, dx, dxf, dy, dz, dt, advFLUX, difFLUX
REAL(real_p) :: CFL,t, Tend, SS, fc, DRAMP, RAMP
REAL(real_p), ALLOCATABLE :: U(:,:,:,:,:), ZETA(:,:,:), Uz(:), Ubar(:,:), Q(:,:), PHIzeta(:,:)
REAL(real_p), ALLOCATABLE :: RHSu(:,:,:,:,:), RHSzeta(:,:,:), dZETA(:,:), RHSubar(:,:,:), Hexact(:,:)
REAL(real_p), ALLOCATABLE :: dZETA_COUPLE(:,:), dZETA_MOM(:,:), USE(:,:,:), Uh(:,:), NX(:,:), KX(:,:)
REAL(real_p), ALLOCATABLE :: X(:), Y(:), Z(:), XF(:), U_3D(:,:,:), V_3D(:,:,:), RHSu_3D(:,:,:)
REAL(real_p), ALLOCATABLE :: xELEM_jac(:), yELEM_jac(:), zELEM_jac(:), xfELEM_jac(:), Hh(:,:)
REAL(real_p), ALLOCATABLE :: xELEM_nodes(:,:), yELEM_nodes(:,:), zELEM_nodes(:,:), xfELEM_nodes(:,:)
REAL(real_p), ALLOCATABLE :: xNODE_elems(:,:), yNODE_elems(:,:), zNODE_elems(:,:), xfNODE_elems(:,:)
REAL(real_p), ALLOCATABLE :: bPTSz(:,:,:), bWTSz(:), C_3D(:,:), B_3D(:,:,:), DPHIzeta_3D(:,:,:)
REAL(real_p), ALLOCATABLE :: PHI(:,:),DPHI(:,:),L2PHI(:,:),L2DPHI(:,:),bPHI(:,:),bDPHI(:,:)
REAL(real_p), ALLOCATABLE :: PHIz(:,:), DPHIz(:,:), L2PHIz(:,:), L2DPHIz(:,:), bPHIz(:,:), bDPHIz(:,:)
REAL(real_p), ALLOCATABLE :: L3PHI(:,:), L3DPHI(:,:), Qnode(:,:), RHSu_2D(:,:,:), RHO(:,:)
REAL(real_p), ALLOCATABLE :: CONN(:,:), CONN_jac(:,:), CONN_jacZ(:,:), CONN_jacX(:,:), Rh(:,:)
REAL(real_p), ALLOCATABLE :: A(:,:), B1(:), B2(:), C(:,:), L2(:,:), L3(:,:), L3U(:,:), DPHI_2D(:,:,:)
REAL(real_p), ALLOCATABLE :: AU(:,:), B1U(:), B2U(:), CU(:,:), L2U(:,:), L2DPHI_2D(:,:,:), L2WTSz(:)
REAL(real_p), ALLOCATABLE :: alpha(:,:),beta(:,:),ct(:), BASIS_2D(:), DBASIS_2D(:,:), L2PTSz(:,:)
REAL(real_p), ALLOCATABLE :: bPHI_2D(:,:,:), bDPHI_2D(:,:,:), A_2D(:,:,:), BU(:,:,:), U_2D(:,:,:)
REAL(real_p), ALLOCATABLE :: DPHIzeta(:,:), CONN_jacF(:,:), UintPTS(:,:,:), PHIint(:,:,:), PHIzeta2(:,:)
REAL(real_p), ALLOCATABLE :: elemPTS(:,:),elemWTS(:),L2PTS(:,:),L2WTS(:),bPTS(:), A2(:,:)
REAL(real_p), ALLOCATABLE :: elemPTSz(:,:),elemWTSz(:),L3PTS(:,:),L3WTS(:), BASIS_3D(:), DBASIS_3D(:,:)
REAL(real_p), ALLOCATABLE :: PHI_3D(:,:), DPHI_3D(:,:,:), L2PHI_3D(:,:), L2DPHI_3D(:,:,:)
REAL(real_p), ALLOCATABLE :: BASIS(:), DBASIS(:), BASISz(:), DBASISz(:), BX(:,:,:), bPHI_3D(:,:,:)
REAL(real_p), ALLOCATABLE :: bPTSx(:,:,:), bPHI_2Dx(:,:,:), CONN_jacS(:,:), bPTS_3D(:,:,:)
REAL(real_p), ALLOCATABLE :: V(:,:,:,:,:), V_2D(:,:,:), Vbar(:,:), VintPTS(:,:,:), VPHIint(:,:,:)
REAL(real_p), ALLOCATABLE :: L2ptsNEW(:,:,:), L2PTS_3D(:,:), L2WTS_3D(:), elemPTS_3D(:,:), elemWTS_3D(:)
REAL(real_p), ALLOCATABLE :: HEXnode(:,:), HEX_jac(:,:), xk(:), RHSv_3D(:,:,:), PHIzeta_2D3D(:,:)
REAL(real_p), ALLOCATABLE :: HEX_jacZ(:,:), HEX_jacX(:,:), HEX_jacY(:,:), HEX_jacS(:,:)
REAL(real_p), ALLOCATABLE :: HEX_jacF(:,:), A_3D(:,:,:), L2_U3D(:,:), bPHI_3D2D(:,:)
REAL(real_p), ALLOCATABLE :: UintPTS_3D(:,:,:), PHIint_3D(:,:,:), PHIzeta_3D(:,:), bPTS_3D2D(:,:)
REAL(real_p), ALLOCATABLE :: Hu(:,:,:), Hv(:,:,:), W_3D(:,:), VSE(:,:,:), b_PHI_2D3Dx(:,:,:)
REAL(real_p), ALLOCATABLE :: AwT(:,:,:),L2wT(:,:),CwT(:,:),BwT(:,:,:), b_PHI_2D3Dy(:,:,:)
REAL(real_p), ALLOCATABLE :: AwB(:,:,:),L2wB(:,:),CwB(:,:),BwB(:,:,:), HEXcolumn(:,:), CONN_jacY(:,:)
REAL(real_p), ALLOCATABLE :: HEX_jacZw(:,:), HEX_jacXw(:,:), HEX_jacYw(:,:), RHSvbar(:,:,:)
REAL(real_p), ALLOCATABLE :: L2_pts1D(:,:), L2_wts1D(:), ELEM_pts1D(:,:), ELEM_wts1D(:)
REAL(real_p), ALLOCATABLE :: QW(:,:), DQW(:,:), RHSw(:,:), nodesZRK(:), dzRK(:), uavgENHANCED_2D(:,:) 
REAL(real_p), ALLOCATABLE :: TPTS_2D(:,:), TWTS_2D(:), TPTS_3D(:,:), TWTS_3D(:), KPTS_1D(:,:), KWTS_1D(:)
REAL(real_p), ALLOCATABLE :: L3w(:,:), L3PHIw(:,:), alphaW(:,:), betaW(:,:), ctW(:), PHI_3Dw(:,:)
REAL(real_p), ALLOCATABLE :: bPTS_3DTP(:,:,:), TPTS_1D(:,:), TWTS_1D(:), W_2D(:,:), bPHI_3Dw(:,:,:)
REAL(real_p), ALLOCATABLE :: bPTS_2DTP(:,:,:), PHI_2D3D(:,:), PHI_TP2D(:,:), DPHI_TP2D(:,:,:)
REAL(real_p), ALLOCATABLE :: bPHI_3DTP(:,:,:), bPHI_2DTP(:,:,:), PHIw(:,:,:), bPHIw(:,:,:,:)
REAL(real_p), ALLOCATABLE :: KC(:,:,:,:,:), DKC(:,:,:,:,:), KR(:,:,:,:,:), DKR(:,:,:,:,:)
REAL(real_p), ALLOCATABLE :: DKL(:,:,:,:,:), DKL2(:,:,:,:,:), KL(:,:,:,:,:), KL2(:,:,:,:,:)
REAL(real_p), ALLOCATABLE :: KC2(:,:,:,:,:), DKC2(:,:,:,:,:), KR2(:,:,:,:,:), DKR2(:,:,:,:,:)
REAL(real_p), ALLOCATABLE :: ZETAstar(:,:), DZETAstar(:,:), UBARstar(:,:), DUBARstar(:,:), bPTSy(:,:,:)
REAL(real_p), ALLOCATABLE :: shiftPTS(:,:), shiftWTS(:), PHIs(:,:), DPHIs(:,:), bPHIs(:,:), bDPHIs(:,:)
REAL(real_p), ALLOCATABLE :: CC(:,:), CD(:,:), ZETA_USE(:,:), zetaENHANCED(:), zetaENHANCED2(:)
REAL(real_p), ALLOCATABLE :: KC_2D(:,:,:,:), zetaENHANCED_2D(:,:), PHIkernel(:,:), wENHANCED_2D(:,:)
REAL(real_p), ALLOCATABLE :: KPTS_2D(:,:), KWTS_2D(:), L2K(:,:), ZETA_2D(:,:,:), SFw(:,:), bPHI_2Dy(:,:,:)
REAL(real_p), ALLOCATABLE :: PHI_TP2Dh(:,:), bPHI_2DTPh(:,:,:), Hu_bar(:,:,:), Hv_bar(:,:,:), Vh(:,:)
REAL(real_p), ALLOCATABLE :: BASIS_2Dw(:), DBASIS_2Dw(:,:), BASIS_3Dw(:), DBASIS_3Dw(:,:), L2PHI_3Dw(:,:)
REAL(real_p), ALLOCATABLE :: L2PHI_h(:,:), PHI_h(:,:), DPHI_h(:,:,:), PHIL2_h(:,:), b_PHI_2D3Dxh(:,:,:)
REAL(real_p), ALLOCATABLE :: b_PHI_2D3Dyh(:,:,:), bPHI_2Dxh(:,:,:), bPHI_2Dyh(:,:,:), PHI_2D3Dh(:,:), L2Bh(:,:)
REAL(real_p), ALLOCATABLE :: PHIr(:,:,:), PHIs_z(:,:,:), PHIq(:,:,:), PHIs_1(:), BASISs(:), DBASISs(:)
REAL(real_p), ALLOCATABLE :: GSExh(:,:), ExFACE(:), EyFACE(:), EAx(:,:,:), EPx(:,:,:), EAy(:,:,:), EPy(:,:,:)
REAL(real_p), ALLOCATABLE :: DPHIzeta_3Dh(:,:,:), hPHI_TP2Dh(:,:), hbPHI_2DTPh(:,:,:), AMIG(:)
REAL(real_p), ALLOCATABLE :: EAx2D(:,:,:), EPx2D(:,:,:), EAy2D(:,:,:), EPy2D(:,:,:)
REAL(real_p), PARAMETER :: z0 = -1.0d0, zN = 0.0d0, PI = 3.1415926535897932d0
COMPLEX(16) :: tt

!****************
!....Data types
!****************

TYPE LYNCH
    REAL(real_p) :: amp
    REAL(real_p) :: freq 
    REAL(real_p) :: phase
    REAL(real_p) :: k
    REAL(real_p) :: g
    REAL(real_p) :: N0
    REAL(real_p) :: h0 
    REAL(real_p) :: x0 
    REAL(real_p) :: xN
    REAL(real_p) :: y0
    REAL(real_p) :: yN
    REAL(real_p) :: Hslope
    REAL(real_p) :: ampU
    REAL(real_p) :: r0
    REAL(real_p) :: rS
    REAL(real_p) :: S
    REAL(real_p) :: f
    REAL(real_p) :: E 
    REAL(real_p) :: gR
    REAL(real_p) :: Lr
    REAL(real_p) :: lambda
    REAL(real_p) :: Lx
END TYPE LYNCH

TYPE(LYNCH) :: PROB_DATA
 
END MODULE global
!*************************************************************************************
PROGRAM DG_SWEM3D
!*************************************************************************************
! A 3D discontinuous Galerkin shallow water equation model (3D DGSWEM)
!             
!                         Developed by:
!
!                        Colton J. Conroy
!                              and
!                        Ethan J. Kubatko
! 
!     @ The Isaac Newton Institute for Mathematical Sciences and
!     The Computational Hydrodynamics and Informatics Lab (C.H.I.L)
!                  @ The Ohio State University
!                            2011-20**
!
!*************************************************************************************
!
! Input Parameters, i.e., Solver Options
!
! coupled    == Method used to couple 3D momentum to PCE (i.e. calc of depth int. vel.)
!            == 1 (Eqns decoupled, i.e., exact solution used for coupling - debugging)
!            == 2 (Depth int. vel. calculated via numerical integration)
!            == 3 (Depth int. mom. equation solved, i.e., mode-splitting)
!            == 4 (Discrete depth int. vel. evaluated exactly - default)
! 
! w_method   == Method used to calculate the transformed vertical velocity
!            == 1 (Calculates vert. vel. via 2D slices using a RKDG method) 
!            == 2 (Calculates vert. vel. using a modified basis - default)
!
! QUADTYPE   == Type of numerical integration points
!            == 1 (Optimal integration points)
!            == 2 (Tensor product points)
!
! FULL3D     == Type of basis used for 3D variables
!            == 0 (Uses full tensor basis - includes cross terms, requires higher int. rules)
!            == 1 (Uses basis without the cross terms, similar to basis used for trianglular prisms)
!
! FULL2D     == Type of basis used for 2D variables
!            == 0 (Uses full tensor basis)
!            == 1 (No cross terms)
!
! NON_LINEAR == Flag to solve either the linear or nonlinear eqns
!            == 0 (Solves the linear eqns)
!            == 1 (Solves the full nonlinear eqns)
! 
! P_Method   == Method to discretize pressure forcing
!            == 1 (Uses derivative of basis functions, i.e., treatment as a source)
!            == 2 (Treat pressure as a flux - default) 
!
! NEW        == RK time steppers flag
!            == 0 (Traditional RK time steppers optimized for TVD stability)
!            == 1 (RK time steppers optimized for L2 stability, i.e., the more stringent condition)
!
! KERNEL     == Flag to turn on kernel convolution (enhances accuracy of 2D variables)
!            == 0 (Does not utilize kernel)
!            == 1 (Enhances accuracy of 2D variables - future work extend to 3D)
!
! HOT_START  == Flag for hot start (ie Tinit /= 0) or cold start (ie Tinit == 0)
!            == 0 (Cold start)
!            == 1 (Hot start - requires input for U, V, W, ZETA, etc.)
!
! VID        == Flag for video output, i.e., files that can be read into MatLab to create a movie
!            == 0 (No video ouput)
!            == 1 (Video output)
!
! PROB       == Problem to be solved (1-5 are analytic test cases with 1-4 being barotropic tests)
!            == 0 (Application/Evaluation)
!            == 1 (Lynch tidal inlet w/quadratic bathymetry)
!            == 2 (Vert. Vel. test case)
!            == 3 (Smooth nonlinear test case with forcing functions)
!            == 4 (Dam break)
!            == 5 (Lynch and Loder shallow stratified sea - baroclinic test)
!
!*************************************************************************************

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ALLSTAT, ALLSTAT2, nrk, stages, NT, n, DEG, REGION, nframes
INTEGER(int_p) :: DEG_3D, RKorder, kind, WRKorder, Wstages, Qe, sumDOF, NEW
INTEGER(int_p), PARAMETER :: REGIONs = 1, timing = 1
REAL(real_p) :: TIME, T1, T2, cpuT, Tinput, tol_check, tol_checkU, tol_checkV
REAL(real_p) :: tol_checkN
REAL(real_p), ALLOCATABLE :: RHS_L2(:), BU_L2(:)

!********************
!....User Input
!********************

WRITE(*,*)'*************************************************************'
WRITE(*,*)'*                                                           *'
WRITE(*,*)'* This program solves the 3D shallow water equations using  *'
WRITE(*,*)'* discontinuous Galerkin methods with SSP RK methods.       *'
WRITE(*,*)'*                                                           *'
WRITE(*,*)'*************************************************************'

WRITE(*,*)'Semi coupled? (1) Integrate Vel? (2) Depth Int. Momentum? (3) or Exact? (4)'
READ(*,*)coupled

WRITE(*,*)'Method to calculate W? (1) RK method (2) Exact method'
READ(*,*)w_method

WRITE(*,*)'Time to run simulation? (secs)'
READ(*,*)Tend

IF (PROB /= 0) THEN

    WRITE(*,*)'Enter x-spacing hi --> hi = 1 => 4elem, hi = 2 => 8elem, etc.'
    READ(*,*)hx
    
    WRITE(*,*)'Enter y-spacing hi --> hi = 1 => 4elem, hi = 2 => 8 elem, etc'
    READ(*,*)hy
    
    WRITE(*,*)'Enter z-spacing hi --> hi = 1 => 4elem, hi = 2 => 8elem, etc.'
    READ(*,*)hz
    
ELSE

    WRITE(*,*)'Would you like to use a ramp for the initial conditions?'
    WRITE(*,*)'1 == yes, 0 == no'
    READ(*,*)NRAMP
    IF (NRAMP == 1) THEN
        WRITE(*,*)'How many days for the ramp?'
        READ(*,*)DRAMP
    ENDIF
    
ENDIF

IF (coupled == 4 .and. w_method == 2) THEN

    WRITE(*,*)'Enter polynomial approximation'
    READ(*,*)p
    pu = p
    pw = p
 
ELSEIF (coupled == 4) THEN

    WRITE(*,*)'Enter polynomial approx. for the Surface Elevation and horizontal velocity'
    READ(*,*)p
    pu = p
    WRITE(*,*)'Enter polynomial approximation for w'
    READ(*,*)pw

ELSE
      
    WRITE(*,*)'Enter polynomial approximation to use for the Surface Elevation'
    READ(*,*)p

    WRITE(*,*)'Enter polynomial approximation to use for the horizontal velocity'
    READ(*,*)pu
    
    WRITE(*,*)'Enter polynomial approximation for w'
    READ(*,*)pw

ENDIF

WRITE(*,*)'Enter polynomial approximation for the Bathymetry'
READ(*,*)ph

WRITE(*,*)'Apply post-processing kernel?'
READ(*,*)KERNEL

P_method = 2
FULL2D   = 1

IF (KERNEL == 1 .and. FULL2D == 0) THEN
    WRITE(*,*)'************************************'
    WRITE(*,*)'* WARNING: Kernel can only be used *'
    WRITE(*,*)'* with full 2D tensor basis -->    *'
    WRITE(*,*)'* Execution will continue BUT with *'
    WRITE(*,*)'* the full 2D tensor basis.        *'
    WRITE(*,*)'************************************'
    FULL2D = 1
ENDIF 
 
IF (timing == 1) THEN
    CALL CPU_TIME(T1)
ENDIF

!************************************
! # Degrees of Freedom/Preliminaries
!************************************

L2TP     = 2*(pu+1)                            !....Number of points used for L2 tensor product
LE       = p+1                                 !....Number of lines per horizontal element    
LRK      = pw+4                                !....Number of 1D int pts for Vertical Vel. Calc.
NPTS_K   = p+4                                 !....Number of Gauss-Radau points for kernel
NPTS_K2D = (NPTS_K)**2                         !....# of 2D G-R points for kernel     
TP_2D    = (LRK)**2                            !....Number of lines per horizontal for vertical vel. calc 
TP_3D    = (LRK)**3
NDOF_2D  = (pw+1)**2                           !....Number of degrees of freedom for 2D calcs for w (RK method)

IF (FULL3D == 0) THEN
    sumDOF = 0
    DO i = 1, pw+1
        DO j = 1, i
            sumDOF = sumDOF + j
        ENDDO
    ENDDO
    NDOF_3Dw = sumDOF 
ELSEIF (FULL3D == 1) THEN
    NDOF_3Dw = (pw+1)**3
ENDIF   
 
REGION = 2  
       
IF (FULL2D == 0) THEN                          !....2D degrees of freedom (Quads)
    NDOF = (p+2)*(p+1)/2
ELSE       
    NDOF = (p+1)**2
ENDIF
                                
IF (p == 0) THEN
    DEG = 1
ELSEIF (p == 1) THEN
    DEG = (p+1)**2 + 1    
ELSEIF (p == 2) THEN                           !....Degree of area integral
    DEG = 10
ELSEIF (p >= 3) THEN
    DEG = 11
ENDIF    
CALL NUM_INT_PTS(DEG,NPTS_2D,REGION)           !....Number of Area integration pts
CALL NUM_INT_PTS(DEG*2,L2_2D,REGION)           !....Number of L2 integration pts

REGION = 3

IF (REGION == 1 .or. REGION == 2) THEN
    IF (coupled == 2) THEN
        ALLOCATE(UintPTS(L2int,2,LE))
        ALLOCATE(PHIint(L2int,NDOF,LE))
        ALLOCATE(VintPTS(L2int,2,L2_2D))
        ALLOCATE(VPHIint(L2int,NDOF,L2_2D))
    ENDIF 
ENDIF    

IF (FULL3D == 0) THEN
    sumDOF = 0
    DO i = 1, pu+1
        DO j = 1, i
            sumDOF = sumDOF + j
        ENDDO
    ENDDO
    NDOF_3D = sumDOF 
ELSEIF (FULL3D == 1) THEN
    NDOF_3D = (pu+1)**3
ENDIF
    
IF (QUAD_TYPE == 1) THEN !....Optimal quadrature rules
    IF (pu == 0) THEN
        DEG_3D = 1
    ELSEIF (pu >= 3) THEN
        DEG_3D  = (pu**2) + 1
    ELSEIF (pu == 1) THEN
        IF (FULL3D == 1) THEN
            DEG_3D  = (pu+1)**2 + pu
        ELSE     
            DEG_3D = 2*pu + 1
        ENDIF 
    ELSE
        DEG_3D  = (pu+1)**2 + pu 
    ENDIF
    CALL NUM_INT_PTS(DEG_3D,NPTS_3D,REGION)
    IF (pu < 2) THEN
        L2_DEG_3D = 10
        CALL NUM_INT_PTS(L2_DEG_3D,L2_3D,REGION) 
    ELSEIF (pu >= 2) THEN
        L2_DEG_3D = 11
        CALL NUM_INT_PTS(L2_DEG_3D,L2_3D,REGION)
    ENDIF 
ELSE  !....Tensor product quadrature
    NPTS_3D = (pu+1)**3
    IF (pu > 2) THEN
        L2TP = 6
    ENDIF    
    L2_3D = (L2TP)**3        
ENDIF
           
IF (coupled == 2) THEN
    ALLOCATE(UintPTS_3D(L2int,3,L2_2D))
    ALLOCATE(PHIint_3D(L2int,NDOF_3D,L2_2D))
    ALLOCATE(VintPTS(L2int,2,L2_2D))
    ALLOCATE(VPHIint(L2int,NDOF_3D,L2_2D))
ENDIF 

IF (pu > 3) THEN
    WRKorder = 4
ELSEIF (pu == 0) THEN
    WRKorder = 2
ELSE
    WRKorder = pw+1
ENDIF

Wstages = 8

CALL WRK_COEFFICIENTS(WRKorder,Wstages)

nrkW = Wstages

!************************
!....Allocate matricies
!************************

REGION = 3
CALL MATRIX_ALLOCATION(REGION)  

!***********************
!....Get problem data
!***********************

CALL PROBLEM_DATA()

!**********************
!....Time step
!**********************
IF (PROB /= 0) THEN
    IF (NON_LINEAR == 1) THEN
        IF (hz == 0) THEN
            dt = 3.0d0/500.0d0
        ELSEIF (hz == 1) THEN
            IF (pu == 1) THEN
                dt = 3.0d0/(500.0d0*hz**3) 
            ELSE
                dt = 3.0d0/(500.0d0*hz**3)   
            ENDIF     
        ELSE 
            dt = 3.0d0/(200.0d0*hz**3)  
        ENDIF
    ELSE
        IF (hz == 0) THEN
            dt = 3.0d0/300.0d0
        ELSEIF (hz == 1) THEN
            IF (pu == 1) THEN
                dt = 3.0d0/(100.0d0*hz**3) 
            ELSE
                dt = 3.0d0/(100.0d0*hz**3)   
            ENDIF     
        ELSE 
            dt = 30.0d0/(300.0d0*hz**3)  
        ENDIF  
    ENDIF       
ELSE
    IF (p == 0) THEN
        dt = 0.500d0
    ELSEIF (p == 1) THEN
        dt = 0.15d0
    ELSE
        dt = 0.020d0    
    ENDIF    
ENDIF
  
write(*,*)'dt=',dt

IF (PROB == 1) THEN
    advFLUX = PROB_DATA%h0  !....Advective flux
    difFLUX = 6.50d-5       !....Diffusive flux
ELSEIF (PROB == 5) THEN
    advFLUX = 0.0d0         !....Problem dependent need to fix
    difFLUX = PROB_DATA%E 
ELSE
    advFLUX = 0.0d0         
    difFLUX = 0.0d0          
ENDIF    

!*********************
!....Discretize domain
!*********************

IF (PROB == 0) THEN

    CALL READ_FORT_1D_MESH()

    CALL READ_FORT_2D_MESH()

    CALL READ_FORT_3D_MESH()
    
    CALL READ_FORT_15()
    
ELSE

    meshX = 1  !....(1) Structured (2) Chebyshev 
    meshY = 1
    meshZ = 1

    CALL LINE_MESH()

    IF (REGION == 2) THEN
        CALL QUAD_MESH(REGION)
    ELSEIF (REGION == 3) THEN
        CALL QUAD_MESHy(REGION)
        CALL HEX_MESH()
    ENDIF    
    
ENDIF

!***********************************
!....Quadrature points and weights
!***********************************

REGION = 3
CALL DG_QUADRATURE(REGION,REGIONs,DEG,DEG_3D)

!**********************************
!....DG Basis
!**********************************

REGION = 3
CALL DG_BASIS(REGION,REGIONs)
REGION = 3

!**********************************
!....Kernel setup
!**********************************

IF (KERNEL == 1) THEN
    CALL KERNEL_PREPROCESS()
ENDIF

!**********************************
!....Mesh for vert. vel. (RK method)
!**********************************

CALL RK_MESH() 

!**********************************
!....DG matricies
!**********************************

!....Continuity eqn

IF (REGION == 1 .OR. REGION == 2) THEN

    CALL DG_MATRICIES(elemWTS,L2WTS,L3WTS,p+1,A,L2,C,B1,B2,L2PHI,PHI,DPHI,bPHI,L2int,NbPTS,L3,LE,L3PHI)
    
ELSEIF (REGION == 3) THEN    

    CALL DG_MATRICIES_2D(p,elemWTSz,L2WTSz,bWTSz,NPTS_2D,L2_2D,p+1,NDOF,A_2D,L2U,CU,BU,BX,L2PHIz,PHIz,DPHI_2D,bPHI_2D,bPHI_2Dx,A2)
    
ENDIF    


!....Momentum eqn

!--------------------------------------------------------------------------------------------------------------------
IF (REGION == 1) THEN       !....Lines
!--------------------------------------------------------------------------------------------------------------------

    CALL DG_MATRICIES(elemWTSz,L2WTS,L3WTS,pu+1,AU,L2U,CU,B1U,B2U,L2PHIz,PHIz,DPHIz,bPHIz,L2int,NbPTS,L3U,LE,L3PHI)
    
!--------------------------------------------------------------------------------------------------------------------
ELSEIF (REGION == 2) THEN   !....Quads
!--------------------------------------------------------------------------------------------------------------------

    CALL DG_MATRICIES_2D(pu,elemWTSz,L2WTSz,bWTSz,NPTS_2D,L2_2D,pu+1,NDOF,A_2D,L2U,CU,BU,BX,L2PHIz,PHIz,DPHI_2D,bPHI_2D,bPHI_2Dx)
    
!--------------------------------------------------------------------------------------------------------------------    
ELSEIF (REGION == 3) THEN   !....Hexs
!--------------------------------------------------------------------------------------------------------------------

    CALL DG_MATRICIES_3D(pu,elemWTS_3D,L2WTS_3D,elemWTSz,NPTS_3D,L2_3D,NPTS_2D,NDOF_3D,A_3D,L2_U3D,C_3D,B_3D,L2PHI_3D,PHI_3D,DPHI_3D,bPHI_3D)

ENDIF

!....Vertical velocity W

CALL DG_W_MATRICIES(pw,p,TWTS_3D,TWTS_2D,TWTS_1D,KWTS_2D,TP_3D,TP_2D,LRK,NPTS_K2D,NDOF_3Dw,NDOF_2D,NDOF,L2K,L3w,L3PHIw,PHIkernel,PHI_TP2D,DPHI_TP2D,bPHI_2DTP,AwT,BwT,KERNEL)

!*********************
!....RK coefficients
!*********************

IF (p > 3) THEN
    RKorder = 4
    stages = 5
ELSEIF (p == 0) THEN
    RKorder = p+1
    stages  = p+1
ELSE
    RKorder = p+1
    stages  = p+2
ENDIF
NEW = 1 !....New optimal time steppers? (1 = yes)
CALL RK_COEFFICIENTS(RKorder,stages,NEW)

nrk = stages

!****************************
!....Allocate main variables
!****************************

!....Depth-independent variables (zeta,ubar,vbar,etc)

IF (REGION == 1 .OR. REGION == 2) THEN

    ALLOCATE(ZETA(p+1,nelemsX,nrk+1))
    ALLOCATE(RHSzeta(p+1,nelemsX,nrk))
    ALLOCATE(USE(p+1,nelemsX,nrk+1))
    ALLOCATE(RHSubar(p+1,nelemsX,nrk))
    ALLOCATE(RHSvbar(p+1,nelemsX,nrk))
    ALLOCATE(Uh(L2int,nelemsX))   
    ALLOCATE(Vh(L2int,nelemsY)) 
    ALLOCATE(NX(p+1,nelemsX))
    ALLOCATE(KX(p+1,nelemsX))
    ALLOCATE(Hh(p+1,nelemsX))
    ALLOCATE(Hexact(L2int,nelemsX))
    
ELSEIF (REGION == 3) THEN

    ALLOCATE(ZETA(NDOF,Qelems,nrk+1))    
    ALLOCATE(RHSzeta(NDOF,Qelems,nrk+1))
    ALLOCATE(USE(NDOF,Qelems,nrk+1))
    ALLOCATE(VSE(NDOF,Qelems,nrk+1))
    ALLOCATE(Hu_bar(NDOF,Qelems,nrk+1))
    ALLOCATE(Hv_bar(NDOF,Qelems,nrk+1))
    ALLOCATE(RHSubar(NDOF,Qelems,nrk))
    ALLOCATE(RHSvbar(NDOF,Qelems,nrk))
    ALLOCATE(Uh(L2_2D,Qelems))
    ALLOCATE(Vh(L2_2D,Qelems))
    ALLOCATE(NX(NDOF,Qelems))
    ALLOCATE(KX(NDOF,Qelems))
    !ALLOCATE(Hh(NDOF,Qelems))
    ALLOCATE(Hexact(L2_2D,Qelems))
    IF (PROB == 5) THEN
        ALLOCATE(GSExh(NDOF,Qelems))
        GSExh = 0.0d0
    ENDIF
ENDIF    

!....Depth-dependent variables (u,v,w,etc)

IF (REGION == 1) THEN

    ALLOCATE(U(pu+1,nelemsZ,nelemsX,LE,nrk+1))
    ALLOCATE(V(pu+1,nelemsZ,nelemsX,LE,nrk+1))
    ALLOCATE(RHSu(pu+1,nelemsZ,nelemsX,LE,nrk))
    ALLOCATE(Q(pu+1,nelemsZ))
    ALLOCATE(Ubar(LE,nelemsX))
    ALLOCATE(Uz(nelemsZ))
    ALLOCATE(dZETA(LE,nelemsX))

ELSEIF (REGION == 2) THEN

    ALLOCATE(U_2D(NDOF,Qelems,nrk+1))
    ALLOCATE(V_2D(NDOF,Qelems,nrk+1))
    ALLOCATE(Hu(NDOF,Qelems,nrk+1))
    ALLOCATE(Hv(NDOF,Qelems,nrk+1))    
    ALLOCATE(RHSu_2D(NDOF,Qelems,nrk))
    ALLOCATE(Q(NDOF,Qelems))
    ALLOCATE(dZETA(L2_2D,nelemsX))
    IF (coupled == 2) THEN
        ALLOCATE(Ubar(LE,nelemsX))
        ALLOCATE(Vbar(LE,nelemsX))
        Ubar    = 0.0d0
        Vbar    = 0.0d0
    ENDIF
    
ELSEIF (REGION == 3) THEN

    ALLOCATE(U_3D(NDOF_3D,HEXelems,nrk+1))
    ALLOCATE(V_3D(NDOF_3D,HEXelems,nrk+1))
    ALLOCATE(W_3D(NDOF_3Dw,HEXelems))
    ALLOCATE(W_2D(NDOF_2D,Qelems))
    ALLOCATE(Hu(NDOF_3D,HEXelems,nrk+1))
    ALLOCATE(Hv(NDOF_3D,HEXelems,nrk+1))
    ALLOCATE(RHSu_3D(NDOF_3D,HEXelems,nrk))
    ALLOCATE(RHSv_3D(NDOF_3D,HEXelems,nrk))    
    ALLOCATE(RHSw(NDOF_2D,nrkW))
    ALLOCATE(Q(NDOF_3D,HEXelems))
    ALLOCATE(QW(NDOF_3D,HEXelems))
    ALLOCATE(DQW(NDOF_3D,HEXelems))
    ALLOCATE(dZETA(L2_3D,Qelems))
    ALLOCATE(RHO(NDOF_3D,HEXelems))
    ALLOCATE(Rh(NDOF_3D,HEXelems))
    IF (coupled == 2) THEN
        ALLOCATE(Ubar(L2_2D,Qelems))
        ALLOCATE(Vbar(L2_2D,Qelems))
        Ubar    = 0.0d0
        Vbar    = 0.0d0
    ENDIF       
        
ENDIF

!***************************
!....Initialize matricies
!***************************

ZETA    = 0.0d0
dZETA   = 0.0d0
RHSzeta = 0.0d0

IF (REGION == 1) THEN
    U       = 0.0d0
    V       = 0.0d0
    Uz      = 0.0d0
    RHSu    = 0.0d0
ELSEIF (REGION == 2) THEN
    U_2D    = 0.0d0
    V_2D    = 0.0d0
    RHSu_2D = 0.0d0
ELSEIF (REGION == 3) THEN
    U_3D    = 0.0d0
    V_3D    = 0.0d0
    W_3D    = 0.0d0
    QW      = 0.0d0
    DQW     = 0.0d0
    Hu      = 0.0d0
    Hv      = 0.0d0
    RHO     = 0.0d0
    Rh      = 0.0d0
    RHSu_3D = 0.0d0   
    RHSv_3D = 0.0d0 
    RHSw    = 0.0d0
ENDIF

IF (NON_LINEAR == 1) THEN
    ALLOCATE(ZETA_2D(NDOF,Qelems,nrk+1))
    ZETA_2D = 0.0d0
ENDIF    

Q = 0.0d0

USE     = 0.0d0
VSE     = 0.0d0
Hu_bar  = 0.0d0
Hv_bar  = 0.0d0
RHSubar = 0.0d0
RHSvbar = 0.0d0
Uh      = 0.0d0
Vh      = 0.0d0
!write(*,*)'hi'
!*******************************************
!....L2 projection of the initial condition
!*******************************************

IF (HOT_START == 1) THEN
    REGION = 3
    Tinput = 255.9522222222220d0
    Tend = Tend - Tinput
    tt = dcmplx(Tinput,0.0d0)
    CALL H_START(tt,REGION)
ELSE
    CALL L2_PROJECTION(REGION)
ENDIF
    
CALL TIDAL_CONST()

!*************************
!....Time stepping setup
!*************************

NT  = floor(Tend/dt)-1
write(*,*)'Number of time steps = ',NT

IF (NRAMP == 1) THEN
    IRAMP = floor((DRAMP*86400.0d0)/dt)-1
ENDIF    

!....Hot or cold start?

IF (HOT_START == 1) THEN
    TIME = Tinput
ELSE
    TIME = 0.0d0
ENDIF 

!*****************
!....Video setup
!*****************

plot_n  = 0
nframes = 0
IF (VID == 1) THEN
    OPEN(17,FILE='time.1d')
ENDIF

!****************
!....If T = 0
!****************

IF ( NT <= 0) THEN
    REGION = 3
    irk = 1
    IF (NON_LINEAR == 1) THEN
        CALL CALC_SE(REGION)
        CALL CALC_U(REGION)
        CALL CALC_V(REGION)
        !CALL CALC_DQW()
    ENDIF
    IF (w_method == 1) THEN
        CALL CALC_W(TIME) 
    ELSE
        CALL W_EXACT()
    ENDIF           
    CALL RHS_dZETA(REGION)
    CALL CALC_R()
    IF (coupled == 2) THEN   
        CALL DEPTH_VEL(REGION)
        DO i = 1, Qelems
            USE(:,i,irk) = MATMUL(L2U,Ubar(:,i))  
            VSE(:,i,irk) = MATMUL(L2U,Vbar(:,i))
        ENDDO 
        CALL CALC_HUbar(REGION)
        CALL CALC_HVbar(REGION)
    ELSEIF (coupled == 3) THEN
        CALL CALC_HUbar(REGION)
        CALL CALC_HVbar(REGION)    
    ELSEIF (coupled == 4) THEN
        CALL DEPTH_VEL_EXACT(REGION)   
        CALL CALC_HUbar(REGION)  
        CALL CALC_HVbar(REGION)         
    ENDIF    
ELSE
    CALL CALC_R()
    IF (coupled == 3) THEN
        irk = 1
        CALL CALC_HUbar(REGION)
        CALL CALC_HVbar(REGION)
    ENDIF
ENDIF

REGION = 3

!**********************
!....Main time loop
!**********************

IF (PROB /= 5) THEN

DO n = 0, NT

    DO irk = 1, nrk

        t  = TIME + ct(irk)*dt
        tt = dcmplx(t,0.0d0)
        
        IF ((NRAMP == 1).AND.(n <= IRAMP)) THEN
            RAMP = TANH((2.0D0*(n + ct(irk))*dt/86400.0D0)/DRAMP)
        ELSE
            RAMP = 1.0d0
        ENDIF       
        
        IF (NON_LINEAR == 1) THEN
	        CALL CALC_SE(REGION)
	        CALL CALC_U(REGION)
            CALL CALC_V(REGION)   
            IF (coupled == 3) THEN
                CALL CALC_Ubar(REGION)
                CALL CALC_Vbar(REGION)
            ENDIF
	    ENDIF
        
        IF (coupled == 2) THEN      !....Depth integrated velocities via summation in vertical  
	        CALL DEPTH_VEL(REGION)
	    ELSEIF (coupled == 4) THEN
	        CALL DEPTH_VEL_EXACT(REGION)
	    ENDIF
	    
	    IF (PROB == 0) THEN         !....Calculate vertical eddy viscosity
	       ! IF (p == 0) THEN
	       !     CALL CALC_Nz(REGION)
	       ! ELSEIF ((p > 0).AND.(n > IRAMP)) THEN
	       !     CALL CALC_Nz(REGION)
	       ! ENDIF      
	    ENDIF    
	    
        IF (NON_LINEAR == 1) THEN
            IF (w_method == 1) THEN
	            CALL CALC_W(TIME)
	        ELSE
	            CALL W_EXACT()
	        ENDIF    
        ENDIF
  
	    CALL RHS_dZETA(REGION)       !....Gradient of surface elevation

        IF (REGION == 1 .OR. REGION == 2) THEN
	    
            IF (NON_LINEAR == 1) THEN         !....DG for surface elevation
                CALL RHS_ZETA_NL(tt)	       
            ELSEIF (NON_LINEAR == 0) THEN    
                CALL RHS_ZETA(tt)  
            ENDIF               
            
        ELSEIF (REGION == 3) THEN

            IF (NON_LINEAR == 1) THEN
                CALL RHS_ZETA2D_NL(tt)
            ELSEIF (NON_LINEAR == 0) THEN
                CALL RHS_ZETA2D(tt)
            ENDIF    

        ENDIF    

        IF (coupled == 1 .or. coupled == 3) THEN      !....Depth integrated velocities via solution of depth int. momentum eqn
        
            IF (REGION == 1 .or. REGION == 2) THEN
                CALL RHS_Ubar(tt,REGION)
            ELSEIF (REGION == 3 .and. NON_LINEAR == 0) THEN
                CALL RHS_Ubar2D(tt,REGION)
                !CALL RHS_Vbar2D(tt,REGION) (need to add)
            ELSEIF (REGION == 3 .and. NON_LINEAR == 1) THEN
                CALL RHS_Ubar2D_NL(tt,REGION)   
                CALL RHS_Vbar2D_NL(tt,REGION) 
            ENDIF    
                
        ENDIF
        
        IF (REGION == 1 .and. NON_LINEAR == 0) THEN       !....DG for horizontal velocity via lines
        
	        DO ii = 1, nelemsX

	            DO kk = 1, LE

		            CALL RHS_U()
		            Q = 0.0d0

	            ENDDO

	        ENDDO
	        
	    ELSEIF(REGION == 2 .and. NON_LINEAR == 0) THEN    !....DG for horizontal velocity via quads
	        
	        CALL RHS_U2D(tt)
	        Q = 0.0d0
	        
	    ELSEIF (REGION == 3 .and. NON_LINEAR == 0) THEN
	    
	        CALL RHS_U3D(tt)
	        Q = 0.0d0  
	        IF (PROB /= 1) THEN 
	            CALL RHS_V3D(tt)
	            Q = 0.0d0 	         
	        ENDIF
	        
	    ELSEIF (REGION == 2 .and. NON_LINEAR == 1) THEN
	    
	        CALL RHS_HU2D(tt)
	        Q = 0.0d0
	    
	    ELSEIF (REGION == 3 .and. NON_LINEAR == 1) THEN
	    
	        CALL RHS_HU3D(tt)
	        Q = 0.0d0
	        IF (PROB /= 4)THEN
	            CALL RHS_HV3D(tt)
	            Q = 0.0d0
	        ENDIF    

	    ENDIF
	    
        DO i = 1, irk
        
            IF (REGION == 1 .and. NON_LINEAR == 0) THEN

                U(:,:,:,:,irk+1) = U(:,:,:,:,irk+1)+alpha(irk,i)*U(:,:,:,:,i)  &
                                 + beta(irk,i)*dt*RHSu(:,:,:,:,i)
                                 
            ELSEIF (REGION == 2 .and. NON_LINEAR == 0) THEN
            
                U_2D(:,:,irk+1) = U_2D(:,:,irk+1)+alpha(irk,i)*U_2D(:,:,i)     &
                                + beta(irk,i)*dt*RHSu_2D(:,:,i)
                                               
                                
            ELSEIF (REGION == 3 .and. NON_LINEAR == 0) THEN
            
                U_3D(:,:,irk+1) = U_3D(:,:,irk+1)+alpha(irk,i)*U_3D(:,:,i)     &
                                + beta(irk,i)*dt*RHSu_3D(:,:,i)      
                                
                V_3D(:,:,irk+1) = V_3D(:,:,irk+1)+alpha(irk,i)*V_3D(:,:,i)     &
                                + beta(irk,i)*dt*RHSv_3D(:,:,i)                                                          
                            
            ELSEIF (REGION == 2 .and. NON_LINEAR == 1) THEN
            
                Hu(:,:,irk+1) = Hu(:,:,irk+1)+alpha(irk,i)*Hu(:,:,i)     &
                                + beta(irk,i)*dt*RHSu_2D(:,:,i) 
            
            ELSEIF (REGION == 3 .and. NON_LINEAR == 1) THEN
            
                Hu(:,:,irk+1) = Hu(:,:,irk+1)+alpha(irk,i)*Hu(:,:,i)     &
                              + beta(irk,i)*dt*RHSu_3D(:,:,i)    
                              
                IF (PROB /= 4) THEN
                    Hv(:,:,irk+1) = Hv(:,:,irk+1)+alpha(irk,i)*Hv(:,:,i)     &
                                  + beta(irk,i)*dt*RHSv_3D(:,:,i) 
                ENDIF                                               
                             
            ENDIF            

	        ZETA(:,:,irk+1) = ZETA(:,:,irk+1)+alpha(irk,i)*ZETA(:,:,i)     &
                            + beta(irk,i)*dt*RHSzeta(:,:,i)                   
                                  
            IF (coupled == 3 .and. NON_LINEAR == 0) THEN
                            
                USE(:,:,irk+1)  = USE(:,:,irk+1)+alpha(irk,i)*USE(:,:,i)   &
                                + beta(irk,i)*dt*RHSubar(:,:,i)       
                                
                VSE(:,:,irk+1)  = VSE(:,:,irk+1)+alpha(irk,i)*VSE(:,:,i)   &
                                + beta(irk,i)*dt*RHSvbar(:,:,i) 
                                
            ELSEIF (coupled == 3 .and. NON_LINEAR == 1) THEN
            
                Hu_bar(:,:,irk+1) = Hu_bar(:,:,irk+1)+alpha(irk,i)*Hu_bar(:,:,i)   &
                                  + beta(irk,i)*dt*RHSubar(:,:,i)       
                                
                Hv_bar(:,:,irk+1) = Hv_bar(:,:,irk+1)+alpha(irk,i)*Hv_bar(:,:,i)   &
                                  + beta(irk,i)*dt*RHSvbar(:,:,i)                                                                           
                                         
            ENDIF
            
        ENDDO
        
        
     ENDDO
       
     ZETA(:,:,1)       = ZETA(:,:,nrk+1)
     ZETA(:,:,2:nrk+1) = 0.0000000000000d0
     RHSzeta(:,:,:)    = 0.0000000000000d0
     
     IF (REGION == 1) THEN

        U(:,:,:,:,1)       = U(:,:,:,:,nrk+1)
        U(:,:,:,:,2:nrk+1) = 0.000000000000000d0
        RHSu(:,:,:,:,:)    = 0.000000000000000d0
        
     ELSEIF (REGION == 2) THEN
     
        IF (NON_LINEAR == 0) THEN
     
            U_2D(:,:,1)       = U_2D(:,:,nrk+1)
            U_2D(:,:,2:nrk+1) = 0.000000000000000d0
            RHSu_2D(:,:,:)    = 0.000000000000000d0
            
        ELSEIF (NON_LINEAR == 1) THEN
        
            Hu(:,:,1)         = Hu(:,:,nrk+1)
            Hu(:,:,2:nrk+1)   = 0.000000000000000d0
            RHSu_2D(:,:,:)    = 0.000000000000000d0     
            
        ENDIF           
        
     ELSEIF (REGION == 3) THEN
     
        IF (NON_LINEAR == 0) THEN
     
            U_3D(:,:,1)       = U_3D(:,:,nrk+1)
            U_3D(:,:,2:nrk+1) = 0.000000000000000d0
            RHSu_3D(:,:,:)    = 0.000000000000000d0   
            
            V_3D(:,:,1)       = V_3D(:,:,nrk+1)
            V_3D(:,:,2:nrk+1) = 0.000000000000000d0
            RHSv_3D(:,:,:)    = 0.000000000000000d0              
            
        ELSEIF (NON_LINEAR == 1) THEN

            Hu(:,:,1)         = Hu(:,:,nrk+1)
            Hu(:,:,2:nrk+1)   = 0.000000000000000d0
            RHSu_3D(:,:,:)    = 0.000000000000000d0      
            
            Hv(:,:,1)         = Hv(:,:,nrk+1)
            Hv(:,:,2:nrk+1)   = 0.000000000000000d0
            RHSv_3D(:,:,:)    = 0.000000000000000d0             
              
        ENDIF
                 
     ENDIF
     
     IF (coupled == 3 .and. NON_LINEAR == 0) THEN
     
        USE(:,:,1)       = USE(:,:,nrk+1)
        USE(:,:,2:nrk+1) = 0.000000000000000d0
        RHSubar(:,:,:)   = 0.000000000000000d0
        
        VSE(:,:,1)       = VSE(:,:,nrk+1)
        VSE(:,:,2:nrk+1) = 0.000000000000000d0
        RHSvbar(:,:,:)   = 0.000000000000000d0   
             
     ELSEIF (coupled == 3 .and. NON_LINEAR == 1) THEN
     
        Hu_bar(:,:,1)       = Hu_bar(:,:,nrk+1)
        Hu_bar(:,:,2:nrk+1) = 0.000000000000000d0
        RHSubar(:,:,:)      = 0.000000000000000d0
        
        Hv_bar(:,:,1)       = Hv_bar(:,:,nrk+1)
        Hv_bar(:,:,2:nrk+1) = 0.000000000000000d0
        RHSvbar(:,:,:)      = 0.000000000000000d0      
     
     ENDIF
     
     IF (VID == 1) THEN
        CALL VIDEO(n,REGION,nframes,TIME)
     ENDIF     

     TIME = TIME + dt
     CALL CPU_TIME(T2)
     cpuT = T2 - T1
     IF (cpuT >= 590000.0d0) THEN
         CALL GLOBAL_OUTPUT(REGION)
         tt = dcmplx(TIME,0.0d0)
         CALL COMPUTE_ERROR(tt,REGION)
         write(*,*)'TIME=',TIME
         write(*,*)n,'time steps complete of',NT
         STOP
     ENDIF
       
ENDDO

ELSE 
    
    irk = 1
    CALL RHS_BAROCLINIC_U()
    U_3D(:,:,1) = RHSv_3D(:,:,1)
!    ALLOCATE(RHS_L2(L2_3D))
!    ALLOCATE(BU_L2(L2_3D))
!    RHS_L2 = 0.0d0
!    BU_L2  = 0.0d0
!    DO i = 1, HEXelems
!        RHS_L2 = MATMUL(L2PHI_3D,RHSv_3D(:,i,1))
!        DO j = 1, L2_3D
!            BU_L2(j) = RHS_L2(j)/PROB_DATA%f
!        ENDDO
!        U_3D(:,i,1) = MATMUL(L2_U3D,BU_L2)
!    ENDDO


    !CALL RHS_BAROCLINIC_V

ENDIF

IF (NON_LINEAR == 0) THEN
    IF (w_method == 1) THEN
	    CALL CALC_W(TIME)
	ELSE
        CALL W_EXACT()
    ENDIF  
ENDIF

!**************************
!....Application of KERNEL
!**************************

IF (KERNEL == 1) THEN
    var = 1
    CALL HYPERBOLIC_KERNEL_2D(KC_2D,ZETA(:,:,1),p,Qelems,ielems_2D,zetaENHANCED_2D,stencil_2D,NDOF,Qstencil,nelemsX,nelemsY,CONNielem,NPTS_K)    
    CALL HYPERBOLIC_KERNEL_2D(KC_2D,USE(:,:,1),p,Qelems,ielems_2D,uavgENHANCED_2D,stencil_2D,NDOF,Qstencil,nelemsX,nelemsY,CONNielem,NPTS_K)
    CALL KERNEL_2D_ERROR(TIME)
ENDIF

IF (timing == 1) THEN
    CALL CPU_TIME(T2)
    WRITE(*,*)'CPU TIME=',T2-T1
ENDIF

IF (VID == 1) THEN
    CLOSE(17)
ENDIF

WRITE(*,*)'TIME=',TIME

IF (VID == 1) THEN
    write(*,*)'# frames=',nframes
ENDIF

REGION = 3
!*******************
!....Global output
!*******************

IF (coupled == 4) THEN
    irk = 1
    CALL DEPTH_VEL_EXACT(REGION)
ENDIF

CALL GLOBAL_OUTPUT(REGION)

!*********************
!....Compute Errors
!*********************

REGION = 3
tt = dcmplx(TIME,0.0d0)

IF (PROB == 4 .and. hx <=2 ) THEN
    DO i = 1, Qelems
        write(*,*)'Zeta=',ZETA(1,i,1)
    ENDDO
    DO i = 1, Qelems
        write(*,*)'Ubar=',USE(1,i,1)
    ENDDO
ELSEIF (PROB /= 0) THEN
    CALL COMPUTE_ERROR(tt,REGION)    
ENDIF        


STOP
END PROGRAM DG_SWEM3D
!-------------------------------------------------------------------
!
!                 DG RHS Continuity EQN Subroutine
!                   Written by Colton J. Conroy
!                   @ the Isaac Newton Institute
!                           12.14.12
!
!--------------------------------------------------------------------
SUBROUTINE RHS_ZETA(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt, UbarL, UbarR
REAL(real_p) :: UL, UR, FL, FR, Fhat, nn, k, hpt
REAL(real_p), DIMENSION(L2int) :: L2Vavg
REAL(real_p), DIMENSION(p+1) :: Uaj, Fj
REAL(real_p), PARAMETER :: half = 0.500000000000000000d0


IF (coupled == 1) THEN

!....Exact solution for depth averaged velocities

    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, nelemsX                                                     !....Loop over elem.

        INVJAC  = 1.0d0/xELEM_jac(i)

        DO j = 1, L2int                                                 !....Loop over L2pts

            xpt = X(i)  + INVJAC*(1.0d0+L2pts(j,1))                          !....Transform pts
            hpt = PROB_DATA%h0 + PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real
                                                                        !....variables to complex
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)   !....Compute exact soln
                                                                        !....@ the L2pts
            L2Vavg(j)    = Vavg
        
        ENDDO

        USE(:,i,irk)    = MATMUL(L2,L2Vavg)                                 !....L2 projection Ubar        

    ENDDO

ELSEIF (coupled == 2) THEN

!....Couple vertical to horizontal

    DO i = 1, nelemsX

	    USE(:,i,irk) = MATMUL(L3,Ubar(:,i))
        
    ENDDO

ENDIF

!....Loop over Elements for area calculations

DO i = 1, nelemsX

    Uaj = MATMUL(PHI,USE(:,i,irk)) 
    Fj  = advE(Uaj)
    RHSzeta(:,i,irk) = xELEM_jac(i)*MATMUL(A,Fj)

ENDDO

IF (debug == 9 .and. irk == 1) THEN
    write(*,*)'Fj=',Fj,'--'
    write(*,*)PROB_DATA%h0*Uaj
ENDIF

!....Loop over element boundaries

UbarL = 0.0000000000000000d0
Fhat = advB(UbarL)
RHSzeta(:,1,irk) = RHSzeta(:,1,irk)+xELEM_jac(1)*B1*Fhat

DO i = 2, nnodesX-1
    
    jL      = xNODE_elems(i,1)
    jR      = xNODE_elems(i,2)
    UL      = DOT_PRODUCT(bPHI(2,:),ZETA(:,jL,irk))
    UR      = DOT_PRODUCT(bPHI(1,:),ZETA(:,jR,irk))
    UbarL   = DOT_PRODUCT(bPHI(2,:),USE(:,jL,1))
    UbarR   = DOT_PRODUCT(bPHI(1,:),USE(:,jR,1))
    FL      = advB(UbarL)
    FR      = advB(UbarR)
    Fhat    = half*(FL + FR - sqrt(PROB_DATA%g*PROB_DATA%h0)*(UR-UL))
    RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) + xELEM_jac(jL)*B2*Fhat
    RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + xELEM_jac(jR)*B1*Fhat

    IF (debug == 9 .and. i ==  nnodesX-1) THEN
        WRITE(*,*)'jL=',jL,'--'
	    WRITE(*,*)'jR=',jR,'--'
	    WRITE(*,*)'UL=',UL,'--'
	    WRITE(*,*)'UR=',UR,'--'
	    WRITE(*,*)'UbarL=',UbarL,'--'
	    WRITE(*,*)'UbarR=',UbarR,'--'
	    WRITE(*,*)'FL=',FL,'--'
	    WRITE(*,*)'FR=',FR,'--'
	    WRITE(*,*)'FHAT=',Fhat,'--'
    ENDIF

ENDDO

UL     = DOT_PRODUCT(bPHI(2,:),ZETA(:,nelemsX,irk))
UbarL  = DOT_PRODUCT(bPHI(2,:),USE(:,nelemsX,irk))
FL     = advB(UbarL) 

!....Dirichlet BC

zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)
xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
hpt = PROB_DATA%h0 + PROB_DATA%Hslope*PROB_DATA%xN**2
CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
UR     = SE
UbarR  = Vavg
FR     = advB(UbarR)
Fhat   = half*(FL + FR - sqrt(PROB_DATA%g*PROB_DATA%h0)*(UR-UL))
RHSzeta(:,nelemsX,irk) = RHSzeta(:,nelemsX,irk) + xELEM_jac(nelemsX)*B2*Fhat

CONTAINS

    FUNCTION advE(Y)
      REAL(real_p), DIMENSION(p+1) :: advE, Y
      advE = advFLUX*Y
    END FUNCTION advE

    FUNCTION advB(Y)
      REAL(real_p) :: advB, Y
      advB = advFLUX*Y
    END FUNCTION advB

END SUBROUTINE RHS_ZETA
!-------------------------------------------------------------------
!
!                DG NL RHS Continuity EQN Subroutine
!                   Written by Colton J. Conroy
!                        @ the C.H.I.L
!                            12.4.13
!
!--------------------------------------------------------------------
SUBROUTINE RHS_ZETA_NL(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt, UbarL, UbarR
REAL(real_p) :: HL, HR, FL, FR, Fhat, g, EIGmax
REAL(real_p), DIMENSION(L2int) :: L2Vavg
REAL(real_p), DIMENSION(p+1) :: Uaj, Fj, Haj
REAL(real_p), PARAMETER :: half = 0.500000000000000000d0

g = PROB_DATA%g

IF (coupled == 1) THEN

!....Exact solution for depth averaged velocities

    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, nelemsX                                                   !....Loop over elem.

        INVJAC = 1/xELEM_jac(i)

        DO j = 1, L2int                                                 !....Loop over L2pts

            xpt = X(i) + INVJAC*(1+L2pts(j,1))                          !....Transform pts
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real
                                                                        !....variables to complex
            IF (PROB == 1) THEN
                !CALL LG3D_SOLUTIONS(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg)   !....Compute exact soln
            ELSE
                WRITE(*,*)'Cannot use exact couple for this problem '   !....@ the L2pts
                STOP    
            ENDIF    
            L2Vavg(j)    = Vavg
        
        ENDDO

        USE(:,i,irk)    = MATMUL(L2,L2Vavg)                                 !....L2 projection Ubar        

    ENDDO

ELSEIF (coupled == 2) THEN

!....Couple vertical to horizontal

    DO i = 1, nelemsX

	    USE(:,i,irk) = MATMUL(L3,Ubar(:,i))
        
    ENDDO
    
ENDIF

!....Loop over Elements for area calculations

DO i = 1, nelemsX

    Uaj = MATMUL(PHI,USE(:,i,irk)) 
    Haj = MATMUL(PHI,ZETA(:,i,irk))
    DO j = 1, p+1
        Fj(j)  = Haj(j)*Uaj(j)
    ENDDO    
    RHSzeta(:,i,irk) = xELEM_jac(i)*MATMUL(A,Fj)

ENDDO

IF (debug == 9 .and. irk == 1) THEN
    write(*,*)'Fj=',Fj,'--'
    write(*,*)PROB_DATA%h0*Uaj
ENDIF

!....Loop over element boundaries

IF (PROB == 1) THEN
    UbarL = 0.0000000000000000d0
    Fhat = advB(UbarL)
ELSEIF (PROB == 2) THEN
    UbarL = 0.0000000000000000d0
    HL    = 10.000000000000000d0
    Fhat  = HL*UbarL
ENDIF    
RHSzeta(:,1,irk) = RHSzeta(:,1,irk)+xELEM_jac(1)*B1*Fhat

DO i = 2, nnodesX-1
    
    jL      = xNODE_elems(i,1)
    jR      = xNODE_elems(i,2)
    HL      = DOT_PRODUCT(bPHI(2,:),ZETA(:,jL,irk))
    HR      = DOT_PRODUCT(bPHI(1,:),ZETA(:,jR,irk))
    UbarL   = DOT_PRODUCT(bPHI(2,:),USE(:,jL,irk))
    UbarR   = DOT_PRODUCT(bPHI(1,:),USE(:,jR,irk))
    EIGmax  = EIG(HL,HR,UbarL,UbarR)
    Fhat    = LLF(HL,HR,UbarL,UbarR,EIGmax)
    RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) + xELEM_jac(jL)*B2*Fhat
    RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + xELEM_jac(jR)*B1*Fhat

    IF (debug == 9 .and. i ==  nnodesX-1) THEN
        WRITE(*,*)'jL=',jL,'--'
	    WRITE(*,*)'jR=',jR,'--'
	    WRITE(*,*)'HL=',HL,'--'
	    WRITE(*,*)'HR=',HR,'--'
	    WRITE(*,*)'UbarL=',UbarL,'--'
	    WRITE(*,*)'UbarR=',UbarR,'--'
	    WRITE(*,*)'FHAT=',Fhat,'--'
    ENDIF

ENDDO

IF (PROB == 1) THEN
    HL     = DOT_PRODUCT(bPHI(2,:),ZETA(:,nelemsX,irk))
    UbarL  = DOT_PRODUCT(bPHI(2,:),USE(:,nelemsX,irk))
    FL     = advB(UbarL) 

    !....Dirichlet BC

    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)
    xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
    !CALL LG3D_SOLUTIONS(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg)
    HR     = SE
    UbarR  = Vavg
    FR     = advB(UbarR)
    Fhat   = half*(FL + FR - sqrt(PROB_DATA%g*PROB_DATA%h0)*(HR-HL))
ELSEIF(PROB == 3) THEN
    UbarR = 0.0d0
    HR = 5.0d0
    Fhat = HR*UbarR
ENDIF
    
RHSzeta(:,nelemsX,irk) = RHSzeta(:,nelemsX,irk) + xELEM_jac(nelemsX)*B2*Fhat

CONTAINS

    FUNCTION LLF(HL,HR,UbarL,UbarR,EIGmax)
        REAL(real_p) :: HL, HR, UbarL, UbarR, EIGmax
        REAL(real_p) :: FL, FR, LLF
        FL = HL*UbarL
        FR = HR*UbarR  
        LLF = half*(FL + FR - abs(EIGmax)*(HR-HL))  
    END FUNCTION LLF      

    FUNCTION EIG(HL,HR,UbarL,UbarR)
        REAL(real_p) :: HL, HR, UbarL,UbarR
        REAL(real_p) :: LEp, LEm, REp, REm, EIG
            LEp = UbarL + sqrt(g*HL)
            LEm = UbarL - sqrt(g*HL)
            REp = UbarR + sqrt(g*HR)
            REm = UbarR - sqrt(g*HR)
            EIG = max(abs(LEp),abs(LEm))
            EIG = max(abs(EIG),abs(REp))
            EIG = max(abs(EIG),abs(REm))
    END FUNCTION EIG        

    FUNCTION advE(Y)
      REAL(real_p), DIMENSION(p+1) :: advE, Y
      advE = advFLUX*Y
    END FUNCTION advE

    FUNCTION advB(Y)
      REAL(real_p) :: advB, Y
      advB = advFLUX*Y
    END FUNCTION advB

END SUBROUTINE RHS_ZETA_NL
!-----------------------------------------------------------------
!
!          2D DG RHS for Depth Integrated Continuity Eqn
!                    Written by Colton J. Conroy
!                         @ the C.H.I.L
!                             8.13.13
!
!------------------------------------------------------------------
SUBROUTINE RHS_ZETA2D(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD, yCOORD
INTEGER(int_p) :: jL, jR, nedgesX, nedgesY, ExELEM, EyELEM, jj, NCYC
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt, ypt, nn, k, hpt
REAL(real_p) :: PER, ARGJ, RFF, ARG
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Ynode
REAL(real_p), DIMENSION(p+1) :: UL, UR, FL, FR, Fhat, UbarL, UbarR
REAL(real_p), DIMENSION(p+1) :: VbarL, VbarR, Huse, HL, HR
REAL(real_p), DIMENSION(L2_2D) :: L2Vavg
REAL(real_p), DIMENSION(NPTS_2D) :: Uaj, Fxj, Fyj, Haj, Vaj
REAL(real_p), PARAMETER :: half = 0.500000000000000000d0

nedgesX = size(xEDGE,1)
nedgesY = size(yEDGE,1)

IF (coupled == 1) THEN

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
        Ynode(1) = Qnode(globalNODES(1),2)
        Ynode(2) = Qnode(globalNODES(2),2)
        
        DO j = 1, L2_2D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
            Hpt = PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt,0.0d0)
            yCOORD = dcmplx(ypt,0.0d0)
            
            IF (PROB == 1) THEN
            
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
                L2Vavg(j)    = Vavg
                    
            ENDIF                     
        ENDDO
        
        USE(:,i,irk)  = MATMUL(L2U,L2Vavg)                                 !....L2 projection Ubar      

    ENDDO

ELSEIF (coupled == 2) THEN

        DO i = 1, Qelems
        
            USE(:,i,irk) = MATMUL(L2U,Ubar(:,i))  
            VSE(:,i,irk) = MATMUL(L2U,Vbar(:,i))
            
        ENDDO 

ENDIF      

!....Loop over elements for area integrals

DO i = 1, Qelems

    Uaj = MATMUL(PHIz,USE(:,i,irk))
    Vaj = MATMUL(PHIz,VSE(:,i,irk))
    !Haj = MATMUL(PHIz,Hh(:,i))
    Haj = MATMUL(PHI_h,Hh(:,i))
    DO j = 1, NPTS_2D
        Fxj(j)  = Haj(j)*Uaj(j) 
        Fyj(j)  = Haj(j)*Vaj(j)
    ENDDO

    RHSzeta(:,i,irk) = MATMUL(A_2D(:,:,1),Fxj)*CONN_jacX(i,1) + MATMUL(A_2D(:,:,2),Fyj)*CONN_jacY(i,1)
    
ENDDO

!....Loop over edges for flux calculations  

DO i = 1, nedgesX

    jL = xEDGE(i,1)
    jR = xEDGE(i,2)
    
    IF (jL == -1) THEN
    
        Fhat = 0.0d0
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)
        
    ELSEIF (jR == 0) THEN
       
        IF (PROB == 1) THEN
        
            xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
            hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
            zCOORD = dcmplx(hpt*Z(nelemsZ), 0.0d0)
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            UR     = SE
            UbarR  = Vavg
            HR     = hpt
            UL     = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
            UbarL  = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
            !HL     = MATMUL(bPHI_2Dx(:,:,2),Hh(:,jL))
            HL     = MATMUL(bPHI_2Dxh(:,:,2),Hh(:,jL))
            DO j = 1, p+1
                FL(j) = HL(j)*UbarL(j)
                FR(j) = HR(j)*UbarR(j)
                Huse(j) = MAX(HL(j),HR(j))
                Fhat(j) = half*(FL(j) + FR(j) - sqrt(PROB_DATA%g*Huse(j))*(UR(j)-UL(j)))
            ENDDO
            
        ELSE
        
            Fhat = 0.0d0
            
        ENDIF
                
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
        
    ELSEIF (jR == -2) THEN  !....Elevation specified B.C.
                
        !....Left and right states
            
        UL    = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))  
        UR    = 0.0d0  !....Initialize UR
        HL    = MATMUL(bPHI_2Dxh(:,:,2),Hh(:,jL))
        Huse  = HL
        UbarL  = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
        UbarR  = UbarL
            
        !....Get edge-element index
            
        DO j = 1, NEedgesX
            IF (ExEDGE(j) == jL) THEN
                ExELEM = j
                exit
            ENDIF
        ENDDO    
            
        !....Loop over int. points
                
        DO j = 1, p+1
                
            DO jj = 1, NCONST  !....# of tidal constituents
                
                IF (AMIG(jj) < 1.0d-12) THEN
                    NCYC = 0
                    PER  = 0.0d0
                ELSE
                    PER  = 2.0d0*pi/AMIG(jj)
                    NCYC = INT(real(TIME)/PER)
                ENDIF

                ARGJ  = AMIG(jj)*(real(TIME) - NCYC*PER)
                RFF   = RAMP
                ARG   = ARGJ  - EPx2D(j,ExELEM,jj)
                UR(j) = UR(j) + EAx2D(j,ExELEM,jj)*RFF*COS(ARG)
                    
            ENDDO
                
            FL(j) = HL(j)*UbarL(j)
            FR(j) = HL(j)*UbarR(j)                
            Fhat(j) = half*(FL(j) + FR(j) - SQRT(PROB_DATA%g*Huse(j))*(UR(j) - UL(j)))
                
        ENDDO  
            
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)              
    
    ELSE
    
        UL    = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
        UbarL = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
        !HL    = MATMUL(bPHI_2Dx(:,:,2),Hh(:,jL))
        HL    = MATMUL(bPHI_2Dxh(:,:,2),Hh(:,jL))
        UR    = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
        UbarR = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
        !HR    = MATMUL(bPHI_2Dx(:,:,1),Hh(:,jR))
        HR    = MATMUL(bPHI_2Dxh(:,:,1),Hh(:,jR))
        DO j = 1, p+1
            FL(j)   = HL(j)*UbarL(j)
            FR(j)   = HR(j)*UbarR(j)
            Huse(j) = MAX(HL(j),HR(j)) 
            Fhat(j) = half*(FL(j) + FR(j) - SQRT(PROB_DATA%g*Huse(j))*(UR(j)-UL(j)))
        ENDDO   
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)
        
    ENDIF
    
ENDDO    

IF (PROB /= 1) THEN

DO i = 1, nedgesY

    jL = yEDGE(i,1)
    jR = yEDGE(i,2)
    
    IF (jL == -1) THEN

        Fhat = 0.0d0
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)
        
    ELSEIF (jL == -2) THEN  !....Elevation specified B.C.
                
        !....Left and right states
            
        UL = 0.0d0  !....Initialize UL
        UR    = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))  
        HR    = MATMUL(bPHI_2Dyh(:,:,1),Hh(:,jR))
        Huse  = HR
        UbarR = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
        UbarL = UbarR
            
        !....Get face-element index
            
        DO j = 1, NEedgesY
            IF (EyEDGE(j) == jR) THEN
                EyELEM = j
                exit
            ENDIF
        ENDDO    
            
        !....Loop over int. points
                
        DO j = 1, p+1
                
            DO jj = 1, NCONST  !....# of tidal constituents
                
                IF (AMIG(jj) < 1.0d-12) THEN
                    NCYC = 0.0d0
                    PER  = 0.0d0
                ELSE
                    PER  = 2.0d0*pi/AMIG(jj)
                    NCYC = INT(TIME/PER)
                ENDIF
        
                ARGJ  = AMIG(jj)*(TIME - NCYC*PER)
                RFF   = RAMP
                ARG   = ARGJ  - EPy2D(j,EyELEM,jj)
                UL(j) = UL(j) + EAy2D(j,EyELEM,jj)*RFF*COS(ARG)
                    
            ENDDO
                
            FL(j)   = HR(j)*UbarL(j)
            FR(j)   = HR(j)*UbarR(j)                
            Fhat(j) = half*(FL(j) + FR(j) - sqrt(PROB_DATA%g*Huse(j))*(UR(j)-UL(j)))
                
        ENDDO  
            
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)        
        
    ELSEIF (jR == 0) THEN
       
        Fhat = 0.0d0
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
    
    ELSE
    
        UL    = MATMUL(bPHI_2Dy(:,:,2),ZETA(:,jL,irk))
        UbarL = MATMUL(bPHI_2Dy(:,:,2),VSE(:,jL,irk))
        HL    = MATMUL(bPHI_2Dyh(:,:,2),Hh(:,jL))
        UR    = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))
        UbarR = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
        HR    = MATMUL(bPHI_2Dyh(:,:,1),Hh(:,jR))
        DO j = 1, p+1
            FL(j)   = HL(j)*UbarL(j)
            FR(j)   = HR(j)*UbarR(j)
            Huse(j) = MAX(HL(j),HR(j)) 
            Fhat(j) = half*(FL(j) + FR(j) - sqrt(PROB_DATA%g*Huse(j))*(UR(j)-UL(j)))
        ENDDO   
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)
        
    ENDIF
    
ENDDO  

ENDIF

CONTAINS

    FUNCTION advE(Y)
      REAL(real_p), DIMENSION(NPTS_2D) :: advE, Y
      advE = advFLUX*Y
    END FUNCTION advE

    FUNCTION advB(Y)
      REAL(real_p), DIMENSION(p+1) :: advB, Y
      advB = advFLUX*Y
    END FUNCTION advB

END SUBROUTINE RHS_ZETA2D
!-----------------------------------------------------------------
!
!          2D DG RHS for Depth Integrated Continuity Eqn (NL)
!                    Written by Colton J. Conroy
!                         @ the C.H.I.L
!                             2.18.14
!
!------------------------------------------------------------------
SUBROUTINE RHS_ZETA2D_NL(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD, yCOORD
INTEGER(int_p) :: jL, jR, nedgesX, nedgesY
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt, ypt, zpt, nn, k, hpt, g, FxmA
REAL(real_p) :: solnT, Vv, Vw, Huse, Zb, DSEx, DSEy, Uavg, Fpce, Fxm, Fym, Fw, FymA
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Ynode
REAL(real_p), DIMENSION(p+1) :: UL, UR, FL, FR, Fhat, UbarL, UbarR, HL, HR, EIGMAX
REAL(real_p), DIMENSION(L2_2D) :: L2Vavg, L2Uavg
REAL(real_p), DIMENSION(NPTS_2D) :: Uaj, Fxj, Fyj, Haj, Vaj, Saj
REAL(real_p), DIMENSION(L2_2D,Qelems) :: L2Fpce
REAL(real_p), PARAMETER :: half = 0.50d0

nedgesX = size(xEDGE,1)
nedgesY = size(yEDGE,1)
solnT = REAL(TIME)
g = PROB_DATA%g

UL  = 0.0d0 
UR  = 0.0d0 
FL  = 0.0d0 
FR  = 0.0d0 
HL  = 0.0d0 
HR  = 0.0d0 
Uaj = 0.0d0 
Fxj = 0.0d0 
Fyj = 0.0d0 
Haj = 0.0d0 
Vaj = 0.0d0 
Saj = 0.0d0
L2Fpce = 0.0d0
Fhat   = 0.0d0 
L2Vavg = 0.0d0
L2Uavg = 0.0d0
UbarL  = 0.0d0 
UbarR  = 0.0d0  
EIGMAX = 0.0d0

IF (PROB == 3) THEN

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
        Ynode(1) = Qnode(globalNODES(1),2)
        Ynode(2) = Qnode(globalNODES(2),2)
        
        DO j = 1, L2_2D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
            Hpt = PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt,0.0d0)
            yCOORD = dcmplx(ypt,0.0d0)
            zCOORD = dcmplx(Hpt*Z(1),0.0d0)
            IF (PROB == 1) THEN
            
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
                L2Vavg(j)    = Vavg
                
            ELSEIF (PROB == 3) THEN
            
                zpt = Z(1)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2Uavg(j) = Uavg
                L2Vavg(j) = Vavg
                L2Fpce(j,i) = Fpce        
                    
            ENDIF                     
        ENDDO
        
        IF (coupled == 1) THEN
            USE(:,i,irk)  = MATMUL(L2U,L2Uavg)                                 !....L2 projection Ubar      
            VSE(:,i,irk)  = MATMUL(L2U,L2Vavg)
        !   FS(:,i)       = MATMUL(L2U,L2Fpce)
        ENDIF

    ENDDO
    
ELSE

    L2Fpce = 0.0d0    

ENDIF      

IF (coupled == 2) THEN
    DO i = 1, Qelems
        
        USE(:,i,irk) = MATMUL(L2U,Ubar(:,i))  
        VSE(:,i,irk) = MATMUL(L2U,Vbar(:,i))
            
    ENDDO 
ENDIF        

!....Loop over elements for area integrals

DO i = 1, Qelems

    Uaj = MATMUL(PHIz,USE(:,i,irk))
    Vaj = MATMUL(PHIz,VSE(:,i,irk))
    Haj = MATMUL(PHIz,ZETA(:,i,irk))
    !Saj = MATMUL(L2PHIz,FS(:,i))
    DO j = 1, NPTS_2D
        Fxj(j)  = Haj(j)*Uaj(j) 
        Fyj(j)  = Haj(j)*Vaj(j)
    ENDDO
    RHSzeta(:,i,irk) = MATMUL(A_2D(:,:,1),Fxj)*CONN_jacX(i,1) + MATMUL(A_2D(:,:,2),Fyj)*CONN_jacY(i,1) &
                     + MATMUL(CU,L2Fpce(:,i))
ENDDO

!....Loop over edges for flux calculations  

DO i = 1, nedgesX

    jL = xEDGE(i,1)
    jR = xEDGE(i,2)
    
    IF (jL == -1) THEN
        IF (PROB == 1) THEN
            UbarL = 0.0d0
            UbarR = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
            UR    = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
            UL    = UR
            HR    = MATMUL(bPHI_2Dx(:,:,1),Hh(:,jR))
            hpt = PROB_DATA%Hslope*PROB_DATA%x0**2
            DO j = 1, p+1
                FL(j)   = 0.0d0
                FR(j)   = HR(j)*UbarR(j)
                Huse    = hpt
                Fhat(j) = half*(FL(j) + FR(j) - sqrt(PROB_DATA%g*Huse)*(UR(j)-UL(j)))
            ENDDO
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UR = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
            HR = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,1,1))*Xnode(1)+(1.0d0+bPTSx(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,2,1))*Ynode(1)+(1.0d0+bPTSx(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HL(j)   = Huse
                UL(j)   = Uavg
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat = LLF(UL,UR,HL,HR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF        
            
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)
        
    ELSEIF (jR == 0) THEN
       
        IF (PROB == 1) THEN
            xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
            hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
            zCOORD = dcmplx(hpt*Z(nelemsZ), 0.0d0)
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            UR     = SE
            UbarR  = Vavg
            HR     = hpt
            UL     = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
            UbarL  = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
            HL     = MATMUL(bPHI_2Dx(:,:,2),Hh(:,jL))
            DO j = 1, p+1
                FL(j) = HL(j)*UbarL(j)
                FR(j) = HR(j)*UbarR(j)
                Huse  = MAX(HL(j),HR(j))
                Fhat(j) = half*(FL(j) + FR(j) - sqrt(PROB_DATA%g*Huse)*(UR(j)-UL(j)))
            ENDDO
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jL,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UL = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
            HL = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
            DO j = 1, p+1
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,1,2))*Xnode(1)+(1.0d0+bPTSx(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,2,2))*Ynode(1)+(1.0d0+bPTSx(j,2,2))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HR(j)   = Huse
                UR(j)   = Uavg
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat = LLF(UL,UR,HL,HR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF             
        
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
    
    ELSE
    
        HL = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
        UL = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
        HR = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
        UR = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
        EIGmax = EIG(HL,HR,UL,UR)
        Fhat = LLF(UL,UR,HL,HR,EIGmax)   
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)
        
    ENDIF
    
ENDDO  

!....Loop over y-edges    

DO i = 1, nedgesY

    jL = yEDGE(i,1)
    jR = yEDGE(i,2)
    
    IF (jL == -1) THEN
        IF (PROB == 1) THEN
            Fhat = 0.0d0
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UR = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
            HR = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,1,1))*Xnode(1)+(1.0d0+bPTSy(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,2,1))*Ynode(1)+(1.0d0+bPTSy(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HL(j)   = Huse
                UL(j)   = Vavg
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat = LLF(UL,UR,HL,HR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF        
            
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)
        
    ELSEIF (jR == 0) THEN
       
        IF (PROB == 1) THEN
            Fhat = 0.0d0
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jL,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UL = MATMUL(bPHI_2Dy(:,:,2),VSE(:,jL,irk))
            HL = MATMUL(bPHI_2Dy(:,:,2),ZETA(:,jL,irk))
            DO j = 1, p+1
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,1,2))*Xnode(1)+(1.0d0+bPTSy(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,2,2))*Ynode(1)+(1.0d0+bPTSy(j,2,2))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HR(j)   = Huse
                UR(j)   = Vavg
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat = LLF(UL,UR,HL,HR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF             
        
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
    
    ELSE
    
        HL = MATMUL(bPHI_2Dy(:,:,2),ZETA(:,jL,irk))
        UL = MATMUL(bPHI_2Dy(:,:,2),VSE(:,jL,irk))
        HR = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))
        UR = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
        EIGmax = EIG(HL,HR,UL,UR)
        Fhat = LLF(UL,UR,HL,HR,EIGmax)   
        RHSzeta(:,jL,irk) = RHSzeta(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
        RHSzeta(:,jR,irk) = RHSzeta(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)

    ENDIF
    
ENDDO  

CONTAINS

    FUNCTION LLF(UL,UR,HL,HR,EIGmax)
    INTEGER(int_p) :: j
    REAL(real_p), DIMENSION(p+1) :: FL, FR, UL, UR, LLF
    REAL(real_p), DIMENSION(p+1) :: HL, HR, EIGmax
        DO j = 1, p+1
            FL(j)  = HL(j)*UL(j) 
            FR(j)  = HR(j)*UR(j)  
            LLF(j) = (1.0d0/2.0d0)*(FL(j) + FR(j) - abs(EIGmax(j))*(HR(j)-HL(j))) 
        ENDDO
    END FUNCTION LLF
    
    FUNCTION EIG(HL,HR,UL,UR)
    REAL(real_p), DIMENSION(p+1) :: UL,UR, HL, HR
    REAL(real_p), DIMENSION(p+1) :: LEp, LEm, REp, REm, EIG
    DO j = 1, p+1
        LEp(j) = UL(j) + sqrt(g*HL(j))
        LEm(j) = UL(j) - sqrt(g*HL(j))
        REp(j) = UR(j) + sqrt(g*HR(j))
        REm(j) = UR(j) - sqrt(g*HR(j))
        EIG(j) = max(abs(LEp(j)),abs(LEm(j)))
        EIG(j) = max(abs(EIG(j)),abs(REp(j)))
        EIG(j) = max(abs(EIG(j)),abs(REm(j)))
        EIG(j) = max(abs(EIG(j)),abs(UL(j)))
        EIG(j) = max(abs(EIG(j)),abs(UR(j)))
    ENDDO    
    END FUNCTION EIG     

    FUNCTION advE(Y)
      REAL(real_p), DIMENSION(NPTS_2D) :: advE, Y
      advE = advFLUX*Y
    END FUNCTION advE

    FUNCTION advB(Y)
      REAL(real_p), DIMENSION(p+1) :: advB, Y
      advB = advFLUX*Y
    END FUNCTION advB

END SUBROUTINE RHS_ZETA2D_NL
!------------------------------------------------------------------
!
!                LDG RHS for Momentum Eqn Subroutine
!                    Written by Colton J. Conroy
!		     @ the Isaac Newton Institute
!                           12.18.12
!
!------------------------------------------------------------------
SUBROUTINE RHS_U()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: jL, jR
REAL(real_p) :: UL, UR, QL, QR, Ghat, Fhat
REAL(real_p), DIMENSION(pu+1) :: Uj, Gj, Qj, Rj, Sj, Bj
REAL(real_p), DIMENSION(L2int) :: dZETAj, RRj

!....Loop over boundaries for Q

UL     = dot_product(bPHIz(1,:),U(:,1,ii,kk,irk))
Ghat   = gfB(UL)
Q(:,1) = Q(:,1) + zELEM_jac(1)*B1U*Ghat

DO i = 2, nnodesZ-1

    jL      = zNODE_elems(i,1)
    jR      = zNODE_elems(i,2)
    UL      = dot_product(bPHIz(2,:),U(:,jL,ii,kk,irk))
    Ghat    = gfB(UL)
    Q(:,jL) = Q(:,jL) + zELEM_jac(jL)*B2U*Ghat
    Q(:,jR) = Q(:,jR) + zELEM_jac(jR)*B1U*Ghat

ENDDO

UL           = dot_product(bPHIz(2,:),U(:,nelemsZ,ii,kk,irk))
Ghat         = gfB(UL)
Q(:,nelemsZ) = Q(:,nelemsZ) + zELEM_jac(nelemsZ)*B2U*Ghat

!....Surface Gradient (Pressure)

dZETAj = -PROB_DATA%g*dZETA(kk,ii)

!....Loop over Elements for Q and U

DO i = 1, nelemsZ

    Uj     = MATMUL(PHIz,U(:,i,ii,kk,irk))
    Gj     = gfE(Uj)
    Q(:,i) = Q(:,i) + zELEM_jac(i)*MATMUL(AU,Gj)
    Qj                  = MATMUL(PHIz,Q(:,i))*sqrt(difFLUX)
    Rj                  = MATMUL(AU,Qj)*zELEM_jac(i)
    Sj                  = MATMUL(CU,dZETAj)
   
    RHSu(:,i,ii,kk,irk) = Rj(:) + Sj(:) 

ENDDO

!....Loop over boundaries for U

UL   = dot_product(bPHIz(1,:),U(:,1,ii,kk,irk))
QL   = -PROB_DATA%k/sqrt(PROB_DATA%N0)*UL
Fhat = SQRT(difFLUX)*QL
RHSu(:,1,ii,kk,irk) = RHSu(:,1,ii,kk,irk) + zELEM_jac(1)*B1U*Fhat

DO i = 2, nnodesZ-1

    jL                   = zNODE_elems(i,1)
    jR                   = zNODE_elems(i,2)
    QR                   = dot_product(bPHIz(1,:),Q(:,jR))
    Fhat                 = SQRT(difFLUX)*QR
    RHSu(:,jL,ii,kk,irk) = RHSu(:,jL,ii,kk,irk) + zELEM_jac(jL)*B2U*Fhat
    RHSu(:,jR,ii,kk,irk) = RHSu(:,jR,ii,kk,irk) + zELEM_jac(jR)*B1U*Fhat

ENDDO

QR = 0.0000000000000d0
Fhat = SQRT(difFLUX)*QR
RHSu(:,nelemsZ,ii,kk,irk) = RHSu(:,nelemsZ,ii,kk,irk) + zELEM_jac(nelemsZ)*B2U*Fhat


IF (debug == 13) THEN
    DO i = 1, nelemsZ
        write(*,*)'RHSu=',RHSu(:,i,ii,kk,irk)
    ENDDO
ENDIF


CONTAINS

    FUNCTION gfB(Y)
	REAL(real_p) :: Y, gfB
	 gfB = sqrt(difFLUX)*Y
    END FUNCTION gfB

    FUNCTION gfE(Y)
	REAL(real_p), DIMENSION(pu+1) :: Y, gfE
	 gfE = sqrt(difFLUX)*Y
    END FUNCTION gfE

END SUBROUTINE RHS_U
!-----------------------------------------------------------------
!
!          RHS for 2D Horizontal Velocities
!              Written by Colton J. Conroy
!                     @ The C.H.I.L
!                          2.1.13 
!
!-----------------------------------------------------------------
SUBROUTINE RHS_U2D(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR, nedgesZ, H_Elem, nedgesX, hL, hR
INTEGER(int_p), PARAMETER :: cheat = 2
REAL(real_p), DIMENSION(pu+1) :: UL, UR, QL, QR, Ghat, Fhat, FL, FR
REAL(real_p), DIMENSION(NPTS_2D) :: Uj, Gj, Qj,RHOj, Eaj
REAL(real_p), DIMENSION(NDOF) :: Rj, Sj, BRj
REAL(real_p), DIMENSION(L2_2D) :: dZETAj, Bj, Vj, RRj
REAL(real_p), DIMENSION(NDOF,Qelems) :: RHOint
REAL(real_p) :: INVJAC, xpt, SE, grad_SE, Vu, Vavg, g, h0, rhoS, C1
REAL(real_p) :: EL, ER, s, hpt, nn, k
REAL(real_p), PARAMETER :: half = 0.5000000000000000000000000000d0

nedgesZ = size(zEDGE,1)
nedgesX = size(xEDGE,1)

IF (cheat == 1) THEN
    
    dZETA = 0.0d0
    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, nelemsX                                                   !....Loop over elem.

        INVJAC = 1.0d0/xELEM_jac(i)

        DO j = 1, L2_2D                                                 !....Loop over L2pts

            xpt = X(i) + INVJAC*(1.0d0+L2PTSz(j,1))                          !....Transform pts
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real
                                                                        !....variables to complex
            CALL LG3D_SOLUTIONS(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg)   !....Compute exact soln
                                                                        !....@ the L2pts
            dZETA(j,i) = grad_SE
        
        ENDDO     

    ENDDO
    
ENDIF
!....Loop over edges for Q

Do i = 1, nedgesZ
    
    jL = zEDGE(i,1)
    jR = zEDGE(i,2)
    
    IF (jL == -1) THEN
        
        UL      = MATMUL(bPHI_2D(:,:,1),U_2D(:,jR,irk))
        Ghat    = gfB(UL)
        Q(:,jR) = Q(:,jR)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Ghat)
        
    ELSEIF (jR == 0) THEN
    
        UL      = MATMUL(bPHI_2D(:,:,2),U_2D(:,jL,irk))
        Ghat    = gfB(UL)
        Q(:,jL) = Q(:,jL)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Ghat)
        
    ELSE
    
        UL   = MATMUL(bPHI_2D(:,:,2),U_2D(:,jL,irk))
        Ghat = gfB(UL)
        Q(:,jL) = Q(:,jL)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Ghat)
        Q(:,jR) = Q(:,jR)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Ghat)
        
    ENDIF
    
ENDDO

!....Loop over elements for Q and U 

g      = PROB_DATA%g
h0     = PROB_DATA%h0

DO i = 1, Qelems

    H_Elem           = CONNzx(i,1)    
    Uj               = MATMUL(PHIz,U_2D(:,i,irk))
    Gj               = gfE(Uj)
    Q(:,i)           = Q(:,i) + CONN_jacF(i,1)*MATMUL(A_2D(:,:,2),Gj)
    Qj               = MATMUL(PHIz,Q(:,i))*SQRT(difFLUX)
    Rj               = MATMUL(A_2D(:,:,2),Qj)*CONN_jacF(i,1)
    IF (P_method == 1) THEN
        dZETAj           = -PROB_DATA%g*dZETA(:,H_elem)
        Sj               = MATMUL(CU,dZETAj)
    ELSEIF (P_method == 2) THEN
        Eaj              = MATMUL(PHIzeta,ZETA(:,H_Elem,irk))
        Sj               = MATMUL(A_2D(:,:,1),Eaj)*CONN_jacX(i,1)*g
    ENDIF
  
    RHSu_2D(:,i,irk) = Rj(:) + Sj(:) 
    
ENDDO

!....Loop over edges for U

DO i = 1, nedgesZ

    jL = zEDGE(i,1)
    jR = zEDGE(i,2)
    
    IF (jL == -1) THEN

        UL   = MATMUL(bPHI_2D(:,:,1),U_2D(:,jR,irk))
        QL   = -PROB_DATA%k/SQRT(PROB_DATA%N0)*UL
        Fhat = SQRT(difFLUX)*QL
        RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Fhat)
        
    ELSEIF (jR == 0) THEN 
        
        QR   = 0
        Fhat = SQRT(difFLUX)*QR
        RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Fhat)
        
    ELSE
     
        QR                = MATMUL(bPHI_2D(:,:,1),Q(:,jR))
        Fhat              = SQRT(difFLUX)*QR
        RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Fhat)
        RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Fhat)
        
    ENDIF
    
ENDDO  

IF (P_method == 2) THEN
    DO i = 1, nedgesX

        jL = xEDGE(i,1)
        jR = xEDGE(i,2)
        hL = CONNzx(jL,1)
        hR = CONNzx(jR,1)

        IF (jL == -1) THEN
    
            ER   = DOT_PRODUCT(bPHI(1,:),ZETA(:,hR,irk))*g
            FR   = ER
            FL   = FR
            UR   = MATMUL(bPHI_2Dx(:,:,1),U_2D(:,jR,irk))
            UL   = UR
            Fhat = FR
            RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacX(jR,1)*MATMUL(BX(:,:,1),Fhat)
            
        ELSEIF (jR == 0) THEN

            EL     = DOT_PRODUCT(bPHI(2,:),ZETA(:,hL,irk))*g
            FL     = EL
            zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)
            xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
            hpt = PROB_DATA%h0 + PROB_DATA%Hslope*PROB_DATA%xN**2
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            FR   = SE
            UL   = MATMUL(bPHI_2Dx(:,:,2),U_2D(:,jL,irk))
            UR   = UL
            Fhat = FR
            RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacX(jL,1)*MATMUL(BX(:,:,2),Fhat)

        ELSE
    
            EL = DOT_PRODUCT(bPHI(2,:),ZETA(:,hL,irk))
            FL = EL*g
            ER = DOT_PRODUCT(bPHI(1,:),ZETA(:,hR,irk))
            FR = ER*g
            UL = MATMUL(bPHI_2Dx(:,:,2),U_2D(:,jL,irk))
            UR = MATMUL(bPHI_2Dx(:,:,1),U_2D(:,jR,irk))
            s = (FR(1) - FL(1))/(UR(1) - UL(1))
            IF (s > 0.0d0) THEN
                Fhat = FL
            ELSE    
                Fhat = FR
            ENDIF
            RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacX(jR,1)*MATMUL(BX(:,:,1),Fhat)
            RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacX(jL,1)*MATMUL(BX(:,:,2),Fhat)        

        ENDIF  
         
    ENDDO 
    
    
ENDIF 

CONTAINS

    FUNCTION gfB(Y)
	REAL(real_p), DIMENSION(pu+1) :: Y, gfB
	    gfB = sqrt(difFLUX)*Y
    END FUNCTION gfB
    
    FUNCTION gfE(Y)
    REAL(real_p), DIMENSION(NPTS_2D) :: Y, gfE
        gfE = sqrt(difFLUX)*Y
    END FUNCTION gfE

END SUBROUTINE RHS_U2D
!-----------------------------------------------------------------
!
!          Non Linear RHS for 2D Horizontal Velocities
!              Written by Colton J. Conroy
!                     @ The C.H.I.L
!                          2.1.13 
!
!-----------------------------------------------------------------
SUBROUTINE RHS_HU2D(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR, nedgesZ, H_Elem, nedgesX, hL, hR
INTEGER(int_p), PARAMETER :: cheat = 2, flux_type = 1
REAL(real_p), DIMENSION(pu+1) :: UL, UR, QL, QR, Ghat, Fhat, FL, FR
REAL(real_p), DIMENSION(pu+1) :: HuL, HuR, EIGmax
REAL(real_p), DIMENSION(NPTS_2D) :: Uj, Gj, Qj, Eaj, Faj, Huj
REAL(real_p), DIMENSION(NDOF) :: Rj, Sj, BRj
REAL(real_p), DIMENSION(L2_2D) :: dZETAj, Bj, Vj, RRj
REAL(real_p) :: INVJAC, xpt, SE, grad_SE, Vu, Vavg, g, h0, rhoS, C1
REAL(real_p) :: EL, ER, s
REAL(real_p), PARAMETER :: half = 0.5000000000000000000000000000d0

nedgesZ = size(zEDGE,1)
nedgesX = size(xEDGE,1)

IF (cheat == 1) THEN
    
    dZETA = 0.0d0
    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, nelemsX                                                   !....Loop over elem.

        INVJAC = 1/xELEM_jac(i)

        DO j = 1, L2_2D                                                 !....Loop over L2pts

            xpt = X(i) + INVJAC*(1+L2PTSz(j,1))                          !....Transform pts
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real
                                                                        !....variables to complex
            !CALL LG3D_SOLUTIONS(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg)   !....Compute exact soln
                                                                        !....@ the L2pts
            dZETA(j,i) = grad_SE
        
        ENDDO     

    ENDDO
    
ENDIF
!....Loop over edges for Q

Do i = 1, nedgesZ
    
    jL = zEDGE(i,1)
    jR = zEDGE(i,2)
    
    IF (jL == -1) THEN
        
        UL      = MATMUL(bPHI_2D(:,:,1),U_2D(:,jR,irk))
        Ghat    = gfB(UL)
        Q(:,jR) = Q(:,jR)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Ghat)
        
    ELSEIF (jR == 0) THEN
    
        UL      = MATMUL(bPHI_2D(:,:,2),U_2D(:,jL,irk))
        Ghat    = gfB(UL)
        Q(:,jL) = Q(:,jL)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Ghat)
        
    ELSE
    
        UL   = MATMUL(bPHI_2D(:,:,2),U_2D(:,jL,irk))
        Ghat = gfB(UL)
        Q(:,jL) = Q(:,jL)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Ghat)
        Q(:,jR) = Q(:,jR)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Ghat)
        
    ENDIF
    
ENDDO

!....Loop over elements for Q and U 

g      = PROB_DATA%g
h0     = PROB_DATA%h0

DO i = 1, Qelems

    H_Elem  = CONNzx(i,1)    
    Uj      = MATMUL(PHIz,U_2D(:,i,irk))
    Huj     = MATMUL(PHIz,Hu(:,i,irk))
    Gj      = gfE(Uj)
    Q(:,i)  = Q(:,i) + CONN_jacF(i,1)*MATMUL(A_2D(:,:,2),Gj)
    Qj      = MATMUL(PHIz,Q(:,i))*SQRT(difFLUX)
    Rj      = MATMUL(A_2D(:,:,2),Qj)*CONN_jacF(i,1)
    IF (P_method == 1) THEN
        dZETAj = -PROB_DATA%g*dZETA(:,H_elem)
        DO j = 1, NPTS_2D
            Faj(j) = Huj(j)*Uj(j)
        ENDDO    
        Sj     = MATMUL(CU,dZETAj) + MATMUL(A_2D(:,:,1),Faj)*CONN_jacX(i,1)
    ELSEIF (P_method == 2) THEN
        Eaj = MATMUL(PHIzeta,ZETA(:,H_Elem,irk))
        DO j = 1, NPTS_2D
            Faj(j) = Huj(j)*Uj(j) + (1.0d0/2.0d0)*g*Eaj(j)**2 ! For Dam break h^2 = 0
        ENDDO    
        Sj = MATMUL(A_2D(:,:,1),Faj)*CONN_jacX(i,1)
    ENDIF
   
    RHSu_2D(:,i,irk) = Rj(:) + Sj(:) 
    
ENDDO

!....Loop over edges for U

DO i = 1, nedgesZ

    jL = zEDGE(i,1)
    jR = zEDGE(i,2)
    
    IF (jL == -1) THEN

        UL   = MATMUL(bPHI_2D(:,:,1),U_2D(:,jR,irk))
        QL   = -PROB_DATA%k/SQRT(PROB_DATA%N0)*UL
        Fhat = SQRT(difFLUX)*QL
        RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Fhat)
        
    ELSEIF (jR == 0) THEN 
        
        QR   = 0
        Fhat = SQRT(difFLUX)*QR
        RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Fhat)
        
    ELSE
     
        QR                = MATMUL(bPHI_2D(:,:,1),Q(:,jR))
        Fhat              = SQRT(difFLUX)*QR
        RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacF(jR,1)*MATMUL(BU(:,:,1),Fhat)
        RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacF(jL,1)*MATMUL(BU(:,:,2),Fhat)
        
    ENDIF
    
ENDDO  

IF (P_method == 2) THEN
    DO i = 1, nedgesX

        jL = xEDGE(i,1)
        jR = xEDGE(i,2)
        hL = CONNzx(jL,1)
        hR = CONNzx(jR,1)

        IF (jL == -1) THEN
            
            IF (PROB == 1) THEN
                ER   = DOT_PRODUCT(bPHI(1,:),ZETA(:,hR,irk))
                FR   = g*ER
                FL   = FR
                UR   = MATMUL(bPHI_2Dx(:,:,1),U_2D(:,jR,irk))
                UL   = UR
                Fhat = FR
            ELSEIF (PROB == 3) THEN
                EL   = 10.0d0
                UL   = 0.0d0
                HuL  = 0.0d0
                Fhat = (1.0d0/2.0d0)*g*EL**2 
            ENDIF    
            RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacX(jR,1)*MATMUL(BX(:,:,1),Fhat)
            
        ELSEIF (jR == 0) THEN
            
            IF (PROB == 1) THEN
                EL     = DOT_PRODUCT(bPHI(2,:),ZETA(:,hL,irk))*g
                FL     = EL
                zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)
                xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
                !CALL LG3D_SOLUTIONS(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg)
                FR   = SE
                UL   = MATMUL(bPHI_2Dx(:,:,2),U_2D(:,jL,irk))
                UR   = UL
                Fhat = g*FR
            ELSEIF (PROB == 3) THEN
                ER   = 5.0d0
                UR   = 0.0d0
                HuR  = 0.0d0
                Fhat = (1.0d0/2.0d0)*g*ER**2
            ENDIF        
            RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacX(jL,1)*MATMUL(BX(:,:,2),Fhat)

        ELSE
    
            EL   = DOT_PRODUCT(bPHI(2,:),ZETA(:,hL,irk))
            ER   = DOT_PRODUCT(bPHI(1,:),ZETA(:,hR,irk))
            UL   = MATMUL(bPHI_2Dx(:,:,2),U_2D(:,jL,irk))
            UR   = MATMUL(bPHI_2Dx(:,:,1),U_2D(:,jR,irk))
            HuL  = MATMUL(bPHI_2Dx(:,:,2),Hu(:,jL,irk))
            HuR  = MATMUL(bPHI_2Dx(:,:,1),Hu(:,jR,irk))
            IF (flux_type == 1) THEN
                EIGmax = EIG(EL,ER,UL,UR)
                Fhat = LLF(EL,ER,UL,UR,HuL,HuR,EIGmax)
            ELSEIF (flux_type == 2) THEN    
                Fhat = Godunov(EL,ER,UL,UR,HuL,HuR)
            ENDIF    
            RHSu_2D(:,jR,irk) = RHSu_2D(:,jR,irk)+CONN_jacX(jR,1)*MATMUL(BX(:,:,1),Fhat)
            RHSu_2D(:,jL,irk) = RHSu_2D(:,jL,irk)-CONN_jacX(jL,1)*MATMUL(BX(:,:,2),Fhat)        

        ENDIF  
         
    ENDDO 
    
    
ENDIF 

CONTAINS

    FUNCTION LLF(EL,ER,UL,UR,HuL,HuR,EIGmax)
    INTEGER(int_p) :: j
    REAL(real_p), DIMENSION(pu+1) :: FL, FR, UL, UR, LLF
    REAL(real_p), DIMENSION(pu+1) :: HuL, HuR, EIGmax
    REAL(real_p) :: EL, ER
        DO j = 1, pu + 1
            FL(j)  = HuL(j)*UL(j) + (1.0d0/2.0d0)*g*EL**2 
            FR(j)  = HuR(j)*UR(j) + (1.0d0/2.0d0)*g*ER**2 
            LLF(j) = (1.0d0/2.0d0)*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLF

    FUNCTION Godunov(EL,ER,UL,UR,HuL,HuR)
    INTEGER(int_p) :: j
    REAL(real_p), DIMENSION(pu+1) :: FL, FR, UL, UR, Godunov
    REAL(real_p), DIMENSION(pu+1) :: HuL, HuR
    REAL(real_p) :: s, DEN, EL, ER
        DO j = 1, pu + 1 
            FL(j) = HuL(j)*UL(j) + (1.0d0/2.0d0)*g*EL**2 ! for dam break h**2 = 0 b/c flat bottom
            FR(j) = HuR(j)*UR(j) + (1.0d0/2.0d0)*g*ER**2 !
            DEN = (HuR(j) - HuL(j))
            IF (abs(DEN) < 1d-16) THEN
                s = 0.0d0
            ELSE    
                s = (FR(j) - FL(j))/DEN
            ENDIF    
            IF (s > 0.0d0) THEN
                Godunov(j) = FL(j)
                IF (UL(j) < 0 .and. UR(j) > 0) THEN  ! Transonic case
                    Godunov(j) = (1.0d0/2.0d0)*g*EL**2
                ENDIF 
            ELSE    
                Godunov(j) = FR(j)
                IF (UL(j) < 0 .and. UR(j) > 0) THEN  ! Transonic case
                    Godunov(j) = (1.0d0/2.0d0)*g*ER**2
                ENDIF 
            ENDIF       
    
        ENDDO
    END FUNCTION Godunov
    
    FUNCTION EIG(EL,ER,UL,UR)
    REAL(real_p), DIMENSION(pu+1) :: UL,UR
    REAL(real_p), DIMENSION(pu+1) :: LEp, LEm, REp, REm, EIG
    REAL(real_p) :: EL, ER
    DO j = 1, pu+1
        LEp(j) = UL(j) + sqrt(g*EL)
        LEm(j) = UL(j) - sqrt(g*EL)
        REp(j) = UR(j) + sqrt(g*ER)
        REm(j) = UR(j) - sqrt(g*ER)
        EIG(j) = max(abs(LEp(j)),abs(LEm(j)))
        EIG(j) = max(abs(EIG(j)),abs(REp(j)))
        EIG(j) = max(abs(EIG(j)),abs(REm(j)))
    ENDDO    
    END FUNCTION EIG  
    
    FUNCTION gfB(Y)
	REAL(real_p), DIMENSION(pu+1) :: Y, gfB
	    gfB = sqrt(difFLUX)*Y
    END FUNCTION gfB
    
    FUNCTION gfE(Y)
    REAL(real_p), DIMENSION(NPTS_2D) :: Y, gfE
        gfE = sqrt(difFLUX)*Y
    END FUNCTION gfE

END SUBROUTINE RHS_HU2D
!------------------------------------------------------------------
!
!             RHS for 3D Horizontal Velocity (Linear Eqns)
!               Written by Colton J. Conroy
!                     @ The C.H.I.L
!                         8.13.13
!
!------------------------------------------------------------------
SUBROUTINE RHS_U3D(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR, nfacesZ, H_Elem, nfacesX, hL, hR
INTEGER(int_p) :: NCYC, jj, ExELEM
INTEGER(int_p), PARAMETER :: flux_type = 2, LDG_FLUX = 0
REAL(real_p), DIMENSION(NPTS_2D) :: EL, ER, ZL, ZR, RR, RL, Hj
REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, QL, QR, Ghat, Nj
REAL(real_p), DIMENSION(NPTS_2D) :: Fhat, FL, FR, Kh, h0, GL, GR
REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Gj, Qj, Eaj, Raj, Naj, Haj
REAL(real_p), DIMENSION(NDOF_3D) :: Rj, Sj, BRj
REAL(real_p), DIMENSION(L2_3D) :: dZETAj, Bj, Vj
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode
REAL(real_p) :: INVJAC, xpt, SE, grad_SE, Vu, Vavg, g, s, hpt, nn, k
REAL(real_p) :: f, PER, ARG, ARGJ, RFF, dzH, dzL, dzR
REAL(real_p), PARAMETER :: half = 0.5000000000000000000000000000d0

!....# of Faces

nfacesZ = size(zFACE,1)
nfacesX = size(xFACE,1)

!....Initialize matricies

EL = 0.0d0 
ER = 0.0d0 
ZL = 0.0d0 
ZR = 0.0d0
UL = 0.0d0
UR = 0.0d0 
QL = 0.0d0 
QR = 0.0d0 
FL = 0.0d0 
FR = 0.0d0 
Kh = 0.0d0 
Nj = 0.0d0
Hj = 0.0d0 
h0 = 0.0d0
Uj = 0.0d0 
Gj = 0.0d0 
Qj = 0.0d0 
Rj = 0.0d0
Sj = 0.0d0
Eaj    = 0.0d0
Naj    = 0.0d0
Haj    = 0.0d0
Ghat   = 0.0d0
Fhat   = 0.0d0
dZETAj = 0.0d0 
Xnode  = 0.0d0
globalNODES = 0.0d0

!....Debugging purposes

IF (coupled == 1) THEN

    dZETA = 0.0d0
    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
       
        DO j = 1, L2_3D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,1))*Xnode(1)+(1.0d0+L2PTS_3D(j,1))*Xnode(2))
            xCOORD = dcmplx(xpt,0.0d0)
            hpt = PROB_DATA%Hslope*xpt**2
                        
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            dZETA(j,i) = grad_SE
                 
        ENDDO
             
    ENDDO
    
ENDIF


!....Loop over faces for Q

DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        H_Elem  = HEXzxy(jR,1)
        Nj      = MATMUL(PHIz,NX(:,H_Elem))
        Hj      = MATMUL(PHI_h,Hh(:,H_elem))
        UL      = MATMUL(bPHI_3D(:,:,3),U_3D(:,jR,irk))
        Ghat    = gfB(UL,Nj,Hj)
        Q(:,jR) = Q(:,jR) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Ghat)
        
    ELSEIF (jR == 0) THEN

        H_Elem  = HEXzxy(jL,1)
        Nj      = MATMUL(PHIz,NX(:,H_Elem))
        Hj      = MATMUL(PHI_h,Hh(:,H_elem))    
        UL      = MATMUL(bPHI_3D(:,:,4),U_3D(:,jL,irk))
        Ghat    = gfB(UL,Nj,Hj)
        Q(:,jL) = Q(:,jL) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Ghat)
        
    ELSE 
    
        H_Elem  = HEXzxy(jL,1)
        Nj      = MATMUL(PHIz,NX(:,H_Elem))
        Hj      = MATMUL(PHI_h,Hh(:,H_elem))        
        UL      = MATMUL(bPHI_3D(:,:,4),U_3D(:,jL,irk))
        IF (LDG_FLUX == 0) THEN
            Ghat = gfB(UL,Nj,Hj)
        ELSE
            UR   = MATMUL(bPHI_3D(:,:,3),U_3D(:,jR,irk)) 
            GR   = gfB(UR,Nj,Hj)
            GL   = gfB(UL,Nj,Hj)
            Ghat = half*(GL + GR) - half*(GR - GL) 
        ENDIF        
        Q(:,jL) = Q(:,jL) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Ghat)
        Q(:,jR) = Q(:,jR) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Ghat)

    ENDIF
        
ENDDO

!....Loop over elements for Q and U          

g = PROB_DATA%g
f = PROB_DATA%f

DO i = 1, HEXelems

    H_Elem = HEXzxy(i,1)
    Naj    = MATMUL(PHI_2D3D,NX(:,H_Elem))
    Haj    = MATMUL(PHI_2D3Dh,Hh(:,H_Elem))
    Uj     = MATMUL(PHI_3D,U_3D(:,i,irk))
    Gj     = gfE(Uj,Naj,Haj)
    Q(:,i) = Q(:,i) + HEX_jacZ(i,1)*MATMUL(A_3D(:,:,3),Gj)
    Qj     = MATMUL(PHI_3D,Q(:,i))
    Gj     = gfE(Qj,Naj,Haj)
    Rj     = MATMUL(A_3D(:,:,3),Gj)*HEX_jacZ(i,1)
    IF (P_method == 1) THEN
        dZETAj = -g*dZETA(:,H_elem)
        Sj     = MATMUL(C_3D,dZETAj)
    ELSEIF (P_method == 2) THEN
        Eaj    = MATMUL(PHI_2D3D,ZETA(:,H_Elem,irk))
        IF (PROB == 5) THEN   !....Baroclinic --> need to fix all of this not necessarily correct!
            Raj    = MATMUL(PHI_3D,Rh(:,i))
            Eaj    = Eaj - g*Raj
            dZETAj = MATMUL(PHIzeta_3D,GSExh(:,H_Elem))
            dZETAj = g*dZETAj
            Vj     = MATMUL(L2PHI_3D,V_3D(:,i,irk))
            Vj     = f*Vj
            BRj    = MATMUL(C_3D,dZETAj) + MATMUL(C_3D,Vj)
        ELSEIF (PROB == 0) THEN
            Vj     = MATMUL(L2PHI_3D,V_3D(:,i,irk))
            Vj     = f*Vj
            BRj    = MATMUL(C_3D,Vj)            
        ELSE
            BRj    = 0.0d0
        ENDIF
        Sj     = MATMUL(A_3D(:,:,1),Eaj)*HEX_jacX(i,1)*g
    ENDIF         
    RHSu_3D(:,i,irk) = Rj(:) + Sj(:) + BRj(:)
    
ENDDO    

!....Loop over faces for U

!....Z-faces

DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        H_Elem = HEXzxy(jR,1)
        UL = MATMUL(bPHI_3D(:,:,3),U_3D(:,jR,irk))
        Kh = MATMUL(PHIz,KX(:,H_Elem))
        Nj = MATMUL(PHIz,NX(:,H_Elem))
        Hj = MATMUL(PHI_h,Hh(:,H_elem)) 

        DO j = 1, NPTS_2D
            QL(j)   = -Kh(j)/SQRT(Nj(j))*UL(j)
            difFLUX = Nj(j)/Hj(j)**2.0d0
            Fhat(j) = SQRT(difFLUX)*QL(j)
        ENDDO    

        RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)
        
    ELSEIF (jR == 0) THEN
    
        QR   = 0.0d0
        Fhat = SQRT(difFLUX)*QR
        RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)
        
    ELSE 
    
        H_Elem = HEXzxy(jR,1)
        Nj   = MATMUL(PHIz,NX(:,H_Elem))
        Hj   = MATMUL(PHI_h,Hh(:,H_elem))
        QR   = MATMUL(bPHI_3D(:,:,3),Q(:,jR))
        IF (LDG_FLUX == 0) THEN
            Fhat = gfB(QR,Nj,Hj)
        ELSE
            UR   = MATMUL(bPHI_3D(:,:,3),U_3D(:,jR,irk))
            UL   = MATMUL(bPHI_3D(:,:,4),U_3D(:,jL,irk))
            QL   = MATMUL(bPHI_3D(:,:,4),Q(:,jL))  
            FL   = gfB(QL,Nj,Hj)   
            FR   = gfB(QR,Nj,Hj)
            dzL  = (2.0d0)*(1.0d0/HEX_jacZ(jL,1))
            dzR  = (2.0d0)*(1.0d0/HEX_jacZ(jR,1))
            dzH  = MAX(dzL,dzR)
            Fhat = half*(FL + FR) - dzH*(UR - UL)
        ENDIF    
        RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)     
        RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)
        
    ENDIF    
        
ENDDO       

!....X-faces
  
IF (P_method == 2) THEN 
                       
    DO i = 1, nfacesX
    
        jL = xFACE(i,1)
        jR = xFACE(i,2)
        hL = HEXzxy(jL,1)
        hR = HEXzxy(jR,1)
        
        IF (jL == -1) THEN
        
            ER = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
            RR = MATMUL(bPHI_3D(:,:,1),Rh(:,jR))            
            DO j = 1, NPTS_2D
                IF (PROB == 5) THEN
                    Fhat(j) = g*(ER(j) - RR(j))
                ELSE
                    Fhat(j) = g*ER(j)
                ENDIF    
            ENDDO    
            RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)           
        
        ELSEIF (jR == 0) THEN
          
            IF (PROB == 1) THEN
                xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
                hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
                zCOORD = dcmplx(hpt*Z(nelemsZ), 0.0d0)
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
                FR   = g*SE
                Fhat = FR   
            ELSEIF (PROB == 5) THEN
                RL = MATMUL(bPHI_3D(:,:,2),Rh(:,jL))
                EL = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
                DO j = 1, NPTS_2D
                    Fhat(j) = g*(EL(j) - RL(j))
                ENDDO    
            
            ELSE
                EL   = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
                Fhat = g*EL                       
            ENDIF
            RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)
        
        ELSEIF (jR == -2) THEN  !....Elevation specified B.C.
                
            !....Left and right states
            
            EL = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))  
            ER = 0.0d0  !....Initialize ER
            ZL = MATMUL(b_PHI_2D3Dxh(:,:,2),Hh(:,hL))
            h0 = ZL
            UL = MATMUL(bPHI_3D(:,:,2),U_3D(:,jL,irk))
            UR = UL
            
            !....Get face-element index
            
            DO j = 1, NEfacesX
                IF (ExFACE(j) == jL) THEN
                    ExELEM = j
                    exit
                ENDIF
            ENDDO    
            
            !....Loop over int. points
                
            DO j = 1, NPTS_2D
                
                DO jj = 1, NCONST  !....# of tidal constituents
                
                    IF (AMIG(jj) < 1.0d-12) THEN
                        NCYC = 0.0d0
                        PER  = 0.0d0
                    ELSE
                        PER  = 2.0d0*pi/AMIG(jj)
                        NCYC = INT(TIME/PER)
                    ENDIF
        
                    ARGJ  = AMIG(jj)*(TIME - NCYC*PER)
                    RFF   = RAMP
                    ARG   = ARGJ - EPx(j,ExELEM,jj)
                    ER(j) = ER(j) + EAx(j,ExELEM,jj)*RFF*COS(ARG)
                    
                ENDDO
                
                FL(j) = g*EL(j)
                FR(j) = g*ER(j)                
                !Fhat(j) = half*(FL(j) + FR(j) - SQRT(g*h0(j))*(UR(j) - UL(j)))
                Fhat(j) = FR(j)
                
            ENDDO        
            
        ELSE
        
            
            EL = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
            !ZL = MATMUL(b_PHI_2D3Dx(:,:,2),Hh(:,hL))
            ZL = MATMUL(b_PHI_2D3Dxh(:,:,2),Hh(:,hL))
            ER = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
            !ZR = MATMUL(b_PHI_2D3Dx(:,:,1),Hh(:,hR))
            ZR = MATMUL(b_PHI_2D3Dxh(:,:,1),Hh(:,hR))
            UL = MATMUL(bPHI_3D(:,:,2),U_3D(:,jL,irk))
            UR = MATMUL(bPHI_3D(:,:,1),U_3D(:,jR,irk))
            DO j = 1, NPTS_2D
                FL(j) = g*EL(j)
                FR(j) = g*ER(j)
                h0(j) = half*(ZL(j)+ZR(j))
            ENDDO    
            IF (flux_type == 1) THEN  !....Godunov
                DO j = 1, NPTS_2D
                    s  = (FR(j)-FL(j))/(UR(j)-UL(j))
                    IF (s > 0.0d0) THEN
                        Fhat(j) = FL(j)
                    ELSE
                        Fhat(j) = FR(j)
                    ENDIF
                ENDDO    
            ELSEIF (flux_type == 2) THEN  !....LLF
                DO j = 1, NPTS_2D
                    Fhat(j) = half*(FL(j) + FR(j) - SQRT(g*h0(j))*(UR(j) - UL(j)))
                ENDDO
            ENDIF
            IF (PROB == 5) THEN
                RL = MATMUL(bPHI_3D(:,:,2),Rh(:,jL))
                RR = MATMUL(bPHI_3D(:,:,1),Rh(:,jR))   
                DO j = 1, NPTS_2D
                    Fhat(j) = Fhat(j) - half*g*(RL(j) + RR(j))
                ENDDO 
            ENDIF                
                  
            RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
            RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)
            
        ENDIF
        
    ENDDO    
        
ENDIF

!....No y-faces for linear problem        

CONTAINS

    FUNCTION gfB(Uz,Nz,Hz)
	REAL(real_p), DIMENSION(NPTS_2D) :: Uz, gfB
	REAL(real_p), DIMENSION(NPTS_2D) :: Nz, Hz
	    DO j = 1, NPTS_2D
            difFLUX = Nz(j)/Hz(j)**2.0d0
            gfB(j)  = sqrt(difFLUX)*Uz(j)
        ENDDO
    END FUNCTION gfB
    
    FUNCTION gfE(Uez,Nez,Hez)
    REAL(real_p), DIMENSION(NPTS_3D) :: Uez, gfE
    REAL(real_p), DIMENSION(NPTS_3D) :: Nez, Hez
        DO j = 1, NPTS_3D
            difFLUX = Nez(j)/Hez(j)**2.0d0
            gfE(j)  = sqrt(difFLUX)*Uez(j)
        ENDDO    
    END FUNCTION gfE

END SUBROUTINE RHS_U3D
!------------------------------------------------------------------
!
!             RHS for 3D Horizontal Velocity (Linear Eqns)
!               Written by Colton J. Conroy
!                     @ The C.H.I.L
!                         5.12.14
!
!------------------------------------------------------------------
SUBROUTINE RHS_V3D(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR, nfacesZ, H_Elem, nfacesY, hL, hR
INTEGER(int_p) :: EyELEM, NCYC, jj
INTEGER(int_p), PARAMETER :: flux_type = 2, LDG_FLUX = 0
REAL(real_p), DIMENSION(NPTS_2D) :: EL, ER, ZL, ZR, Hj
REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, QL, QR, Ghat, GL, GR
REAL(real_p), DIMENSION(NPTS_2D) :: Fhat, FL, FR, Kh, Nj, h0
REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Gj, Qj, Eaj, Haj, Naj
REAL(real_p), DIMENSION(NDOF_3D) :: Rj, Sj, BRj 
REAL(real_p), DIMENSION(L2_3D) :: dZETAj, Vj
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode
REAL(real_p) :: INVJAC, xpt, SE, grad_SE, Vu, Vavg, g, s, hpt, nn, k
REAL(real_p) :: f, L2_RESID, PER, ARG, ARGJ, RFF, dzH, dzL, dzR
REAL(real_p), PARAMETER :: half = 0.5000000000000000000000000000d0

!....# of Faces

nfacesZ = size(zFACE,1)
nfacesY = size(yFACE,1)

!....Initialize matricies

EL = 0.0d0 
ER = 0.0d0 
ZL = 0.0d0 
ZR = 0.0d0
UL = 0.0d0
UR = 0.0d0 
QL = 0.0d0 
QR = 0.0d0 
FL = 0.0d0 
FR = 0.0d0 
Kh = 0.0d0 
Nj = 0.0d0 
Hj = 0.0d0
h0 = 0.0d0
Uj = 0.0d0 
Gj = 0.0d0 
Qj = 0.0d0 
Rj = 0.0d0
Sj = 0.0d0
Eaj    = 0.0d0
Naj    = 0.0d0
Haj    = 0.0d0
Ghat   = 0.0d0
Fhat   = 0.0d0
dZETAj = 0.0d0 
Xnode  = 0.0d0
globalNODES = 0.0d0

!....Debugging purposes

IF (coupled == 1) THEN

    dZETA = 0.0d0
    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
       
        DO j = 1, L2_3D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,1))*Xnode(1)+(1.0d0+L2PTS_3D(j,1))*Xnode(2))
            xCOORD = dcmplx(xpt,0.0d0)
            hpt = PROB_DATA%Hslope*xpt**2
                        
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            dZETA(j,i) = grad_SE
                 
        ENDDO
             
    ENDDO
    
ENDIF

!....Loop over faces for Q

DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        H_Elem  = HEXzxy(jR,1)
        Nj      = MATMUL(PHIz,NX(:,H_Elem))
        Hj      = MATMUL(PHI_h,Hh(:,H_elem))        
        UL      = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk))
        Ghat    = gfB(UL,Nj,Hj)
        Q(:,jR) = Q(:,jR) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Ghat)

    ELSEIF (jR == 0) THEN
    
        H_Elem  = HEXzxy(jL,1)
        Nj      = MATMUL(PHIz,NX(:,H_Elem))
        Hj      = MATMUL(PHI_h,Hh(:,H_elem))        
        UL      = MATMUL(bPHI_3D(:,:,4),V_3D(:,jL,irk))
        Ghat    = gfB(UL,Nj,Hj)
        Q(:,jL) = Q(:,jL) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Ghat)
        
    ELSE 
    
        H_Elem  = HEXzxy(jL,1)
        Nj      = MATMUL(PHIz,NX(:,H_Elem))
        Hj      = MATMUL(PHI_h,Hh(:,H_elem))        
        UL      = MATMUL(bPHI_3D(:,:,4),V_3D(:,jL,irk))
        IF (LDG_FLUX == 0) THEN
            Ghat = gfB(UL,Nj,Hj)
        ELSE
            UR   = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk)) 
            GR   = gfB(UR,Nj,Hj)
            GL   = gfB(UL,Nj,Hj)
            Ghat = half*(GL + GR) - half*(GR - GL) 
        ENDIF
        Q(:,jL) = Q(:,jL) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Ghat)
        Q(:,jR) = Q(:,jR) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Ghat)
        
    ENDIF
    
ENDDO

!....Loop over elements for Q and U          

g = PROB_DATA%g
f = PROB_DATA%f

DO i = 1, HEXelems

    H_Elem = HEXzxy(i,1)
    Naj    = MATMUL(PHI_2D3D,NX(:,H_Elem))
    Haj    = MATMUL(PHI_2D3Dh,Hh(:,H_Elem))    
    Uj     = MATMUL(PHI_3D,V_3D(:,i,irk))
    Gj     = gfE(Uj,Naj,Haj)
    Q(:,i) = Q(:,i) + HEX_jacZ(i,1)*MATMUL(A_3D(:,:,3),Gj)
    Qj     = MATMUL(PHI_3D,Q(:,i))
    Gj     = gfE(Qj,Naj,Haj)
    Rj     = MATMUL(A_3D(:,:,3),Qj)*HEX_jacZ(i,1)
    IF (P_method == 1) THEN
        dZETAj           = -g*dZETA(:,H_elem)
        Sj               = MATMUL(C_3D,dZETAj)
    ELSEIF (P_method == 2) THEN
        Eaj              = MATMUL(PHI_2D3D,ZETA(:,H_Elem,irk))
        IF (PROB == 5) THEN
            !Raj    = MATMUL(PHI_3D,Rh(:,i))
            !Eaj    = Eaj + Raj
            !dZETAj = MATMUL(PHIzeta_3D,GSExh(:,H_Elem))
            !dZETAj = -g*dZETAj
            Vj  = MATMUL(L2PHI_3D,U_3D(:,i,irk))
            Vj  = -f*Vj
            !BRj    = MATMUL(C_3D,dZETAj) + MATMUL(C_3D,Vj)
            BRj = MATMUL(C_3D,Vj)
        ELSEIF (PROB == 0) THEN
            Vj  = MATMUL(L2PHI_3D,U_3D(:,i,irk))
            Vj  = -f*Vj
            BRj = MATMUL(C_3D,Vj)
        ELSE
            BRj = 0.0d0    
        ENDIF        
        Sj = MATMUL(A_3D(:,:,2),Eaj)*HEX_jacY(i,1)*g
    ENDIF         
    RHSv_3D(:,i,irk) = Rj(:) + Sj(:) + BRj(:)

ENDDO    

!....Loop over faces for U

!....Z-faces

DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        H_Elem = HEXzxy(jR,1)
        UL   = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk))
        Kh   = MATMUL(PHIz,KX(:,H_Elem))
        Nj   = MATMUL(PHIz,NX(:,H_Elem))
        Hj   = MATMUL(PHI_h,Hh(:,H_elem)) 

        DO j = 1, NPTS_2D
            QL(j)   = -Kh(j)/SQRT(Nj(j))*UL(j)
            difFLUX = Nj(j)/Hj(j)**2.0d0
            Fhat(j) = SQRT(difFLUX)*QL(j)
        ENDDO    

        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)

    ELSEIF (jR == 0) THEN
    
        QR   = 0.0d0
        Fhat = SQRT(difFLUX)*QR
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)

    ELSE 
    
        H_Elem = HEXzxy(jR,1)
        Nj   = MATMUL(PHIz,NX(:,H_Elem))
        Hj   = MATMUL(PHI_h,Hh(:,H_elem))        
        QR   = MATMUL(bPHI_3D(:,:,3),Q(:,jR))
        IF (LDG_FLUX == 0) THEN !....upwind flux
            Fhat = gfB(QR,Nj,Hj)
        ELSE !....LDG flux 
            UR   = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk))
            UL   = MATMUL(bPHI_3D(:,:,4),V_3D(:,jL,irk))
            QL   = MATMUL(bPHI_3D(:,:,4),Q(:,jL))  
            FL   = gfB(QL,Nj,Hj)   
            FR   = gfB(QR,Nj,Hj)
            dzL  = (2.0d0)*(1.0d0/HEX_jacZ(jL,1))
            dzR  = (2.0d0)*(1.0d0/HEX_jacZ(jR,1))
            dzH  = MAX(dzL,dzR)
            Fhat = half*(FL + FR) - dzH*(UR - UL)
        ENDIF 
        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)     
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)

    ENDIF    
        
ENDDO      

!....Y-faces

IF (P_method == 2) THEN 
                        
    DO i = 1, nfacesY
    
        jL = yFACE(i,1)
        jR = yFACE(i,2)
        hL = HEXzxy(jL,1)
        hR = HEXzxy(jR,1)
        
        IF (jL == -1) THEN
        
            ER   = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
            DO j = 1, NPTS_2D
                Fhat(j) = g*ER(j)
            ENDDO
                
            RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
            
        ELSEIF (jL == -2) THEN  !....Elevation specified B.C.
                
            !....Left and right states
            
            EL = 0.0d0  !....Initialize EL
            ER = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))  
            ZR = MATMUL(b_PHI_2D3Dyh(:,:,1),Hh(:,hR))
            h0 = ZR
            UR = MATMUL(bPHI_3D(:,:,5),V_3D(:,jR,irk))
            UL = UR
            
            !....Get face-element index
            
            DO j = 1, NEfacesY
                IF (EyFACE(j) == jR) THEN
                    EyELEM = j
                    exit
                ENDIF
            ENDDO    
            
            !....Loop over int. points
                
            DO j = 1, NPTS_2D
                
                DO jj = 1, NCONST  !....# of tidal constituents
                
                    IF (AMIG(jj) < 1.0d-12) THEN
                        NCYC = 0.0d0
                        PER  = 0.0d0
                    ELSE
                        PER  = 2.0d0*pi/AMIG(jj)
                        NCYC = INT(TIME/PER)
                    ENDIF
        
                    ARGJ  = AMIG(jj)*(TIME - NCYC*PER)
                    RFF   = RAMP
                    ARG   = ARGJ - EPy(j,EyELEM,jj)
                    EL(j) = EL(j) + EAy(j,EyELEM,jj)*RFF*COS(ARG)
                    
                ENDDO
                
                FL(j) = g*EL(j)
                FR(j) = g*ER(j)                
                Fhat(j) = half*(FL(j) + FR(j) - SQRT(g*h0(j))*(UR(j) - UL(j)))
                !Fhat(j) = FL(j)
                
            ENDDO  
            
            RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)           
        
        ELSEIF (jR == 0) THEN
          
            IF (PROB == 1) THEN
                xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
                hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
                zCOORD = dcmplx(hpt*Z(nelemsZ), 0.0d0)
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
                FR   = g*SE
                Fhat = FR  
            ELSEIF (PROB == 5) THEN
                Fhat = 0.0d0
            ELSE
                EL = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
                Fhat = g*EL    
            ENDIF
                 
            RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
        ELSE
        
            EL = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
            !ZL = MATMUL(b_PHI_2D3Dy(:,:,2),Hh(:,hL))
            ZL = MATMUL(b_PHI_2D3Dyh(:,:,2),Hh(:,hL))
            ER = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
            !ZR = MATMUL(b_PHI_2D3Dy(:,:,1),Hh(:,hR))
            ZR = MATMUL(b_PHI_2D3Dyh(:,:,1),Hh(:,hR))
            UL = MATMUL(bPHI_3D(:,:,6),V_3D(:,jL,irk))
            UR = MATMUL(bPHI_3D(:,:,5),V_3D(:,jR,irk))
            DO j = 1, NPTS_2D
                FL(j) = g*EL(j)
                FR(j) = g*ER(j)
                h0(j) = half*(ZL(j)+ZR(j))
            ENDDO    
            IF (flux_type == 1) THEN  !....Godunov
                DO j = 1, NPTS_2D
                    s  = (FR(j)-FL(j))/(UR(j)-UL(j))
                    IF (s > 0.0d0) THEN
                        Fhat(j) = FL(j)
                    ELSE
                        Fhat(j) = FR(j)
                    ENDIF
                ENDDO    
            ELSEIF (flux_type == 2) THEN  !....LLF
                DO j = 1, NPTS_2D
                    Fhat(j) = half*(FL(j) + FR(j) - SQRT(g*h0(j))*(UR(j) - UL(j)))
                ENDDO
            ENDIF
                  
            RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
            RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
        ENDIF
        
    ENDDO    
     
ENDIF    

!....No x-faces for linear eqns    

CONTAINS

    FUNCTION gfB(Uz,Nz,Hz)
	REAL(real_p), DIMENSION(NPTS_2D) :: Uz, gfB
	REAL(real_p), DIMENSION(NPTS_2D) :: Nz, Hz
	    DO j = 1, NPTS_2D
            difFLUX = Nz(j)/Hz(j)**2.0d0
            gfB(j)  = sqrt(difFLUX)*Uz(j)
        ENDDO
    END FUNCTION gfB
    
    FUNCTION gfE(Uez,Nez,Hez)
    REAL(real_p), DIMENSION(NPTS_3D) :: Uez, gfE
    REAL(real_p), DIMENSION(NPTS_3D) :: Nez, Hez
        DO j = 1, NPTS_3D
            difFLUX = Nez(j)/Hez(j)**2.0d0
            gfE(j)  = sqrt(difFLUX)*Uez(j)
        ENDDO    
    END FUNCTION gfE

END SUBROUTINE RHS_V3D
!------------------------------------------------------------------
!
!   RHS for 3D Horizontal Velocity (Linear Eqns - BAROCLINIC)
!               Written by Colton J. Conroy
!                     @ The C.H.I.L
!                         8.12.14
!
!------------------------------------------------------------------
SUBROUTINE RHS_BAROCLINIC_U()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: jL, jR, nfacesZ, H_Elem, nfacesY, hL, hR
INTEGER(int_p), PARAMETER :: flux_type = 2
REAL(real_p), DIMENSION(NPTS_2D) :: EL, ER, ZL, ZR, h0
REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, QL, QR, Ghat
REAL(real_p), DIMENSION(NPTS_2D) :: Fhat, FL, FR, Kh, Nh
REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Gj, Qj, Eaj, Vj 
REAL(real_p), DIMENSION(L2_3D) :: dZETAj
REAL(real_p) :: g, f, s
REAL(real_p), PARAMETER :: half = 0.5000000000000d0

!....# of Faces

nfacesZ = size(zFACE,1)
nfacesY = size(yFACE,1)

!....Initialize matricies

EL = 0.0d0 
ER = 0.0d0 
ZL = 0.0d0 
ZR = 0.0d0
UL = 0.0d0
UR = 0.0d0 
QL = 0.0d0 
QR = 0.0d0 
FL = 0.0d0 
FR = 0.0d0 
h0 = 0.0d0
Kh = 0.0d0 
Nh = 0.0d0 
Uj = 0.0d0 
Gj = 0.0d0 
Qj = 0.0d0 
Eaj    = 0.0d0
Ghat   = 0.0d0
Fhat   = 0.0d0
dZETAj = 0.0d0 

!....Loop over faces for Q

DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        UL      = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk))
        Ghat    = gfB(UL)
        Q(:,jR) = Q(:,jR) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Ghat)
    
    ELSEIF (jR == 0) THEN
    
        UL      = MATMUL(bPHI_3D(:,:,4),V_3D(:,jL,irk))
        Ghat    = gfB(UL)
        Q(:,jL) = Q(:,jL) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Ghat)
        
    ELSE 
    
        UL   = MATMUL(bPHI_3D(:,:,4),V_3D(:,jL,irk))
        Ghat = gfB(UL)
        Q(:,jL) = Q(:,jL) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Ghat)
        Q(:,jR) = Q(:,jR) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Ghat)
        
    ENDIF
    
ENDDO

!....Loop over elements for Q and U          

g = PROB_DATA%g
f = PROB_DATA%f

DO i = 1, HEXelems

    Uj     = MATMUL(PHI_3D,V_3D(:,i,irk))
    Gj     = SQRT(difFLUX)*Uj
    Q(:,i) = Q(:,i) + HEX_jacZ(i,1)*MATMUL(A_3D(:,:,3),Gj)
    
    Qj               = MATMUL(PHI_3D,Q(:,i))
    Qj               = SQRT(difFLUX)*Qj    
    RHSv_3D(:,i,irk) = MATMUL(A_3D(:,:,3),Qj)*HEX_jacZ(i,1)

ENDDO    

!....Loop over faces for U

!....Z-faces

DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        H_Elem = HEXzxy(jR,1)
        UL   = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk))
        Kh   = MATMUL(PHIz,KX(:,H_Elem))
        Nh   = MATMUL(PHIz,NX(:,H_Elem))

        DO j = 1, NPTS_2D
            QL(j)   = -Kh(j)/SQRT(f*Nh(j))*UL(j)
            Fhat(j) = SQRT(difFLUX)*QL(j)
        ENDDO    

        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)

    ELSEIF (jR == 0) THEN
    
        QR   = 0.0d0
        Fhat = SQRT(difFLUX)*QR
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)

    ELSE 
    
        QR   = MATMUL(bPHI_3D(:,:,3),Q(:,jR))
        Fhat = SQRT(difFLUX)*QR
        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)     
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)

    ENDIF    
        
ENDDO      

!....Y-faces (Currently not taking into account b/c dzeta/dy = 0 and dRdy = 0)

!IF (P_method == 2) THEN 
!                        
!    DO i = 1, nfacesY
!    
!        jL = yFACE(i,1)
!        jR = yFACE(i,2)
!        hL = HEXzxy(jL,1)
!        hR = HEXzxy(jR,1)
!        
!        IF (jL == -1) THEN
!        
!            ER   = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
!            DO j = 1, NPTS_2D
!                Fhat(j) = g*ER(j)
!            ENDDO
!                
!            RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
!        
!        ELSEIF (jR == 0) THEN
!          
!            IF (PROB == 5) THEN
!                Fhat = 0.0d0
!            ENDIF
!                 
!            RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
!            
!        ELSE
!        
!            EL = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
!            ZL = MATMUL(b_PHI_2D3Dy(:,:,2),Hh(:,hL))
!            ER = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
!            ZR = MATMUL(b_PHI_2D3Dy(:,:,1),Hh(:,hR))
!            UL = MATMUL(bPHI_3D(:,:,6),V_3D(:,jL,irk))
!            UR = MATMUL(bPHI_3D(:,:,5),V_3D(:,jR,irk))
!            DO j = 1, NPTS_2D
!                FL(j) = g*EL(j)
!                FR(j) = g*ER(j)
!                h0(j) = half*(ZL(j)+ZR(j))
!            ENDDO    
!            IF (flux_type == 1) THEN  !....Godunov
!                DO j = 1, NPTS_2D
!                    s  = (FR(j)-FL(j))/(UR(j)-UL(j))
!                    IF (s > 0.0d0) THEN
!                        Fhat(j) = FL(j)
!                    ELSE
!                        Fhat(j) = FR(j)
!                    ENDIF
!                ENDDO    
!            ELSEIF (flux_type == 2) THEN  !....LLF
!                DO j = 1, NPTS_2D
!                    Fhat(j) = half*(FL(j) + FR(j) - SQRT(g*h0(j))*(UR(j) - UL(j)))
!                ENDDO
!            ENDIF
!                  
!            RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
!            RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
!            
!        ENDIF
!        
!    ENDDO    
!     
!ENDIF    

!....No x-faces for linear eqns    

CONTAINS

    FUNCTION gfB(Y)
	REAL(real_p), DIMENSION(NPTS_2D) :: Y, gfB
	    gfB = sqrt(difFLUX)*Y
    END FUNCTION gfB
    
    FUNCTION gfE(Y)
    REAL(real_p), DIMENSION(NPTS_3D) :: Y, gfE
        gfE = sqrt(difFLUX)*Y
    END FUNCTION gfE

END SUBROUTINE RHS_BAROCLINIC_U
!------------------------------------------------------------------
!
!             RHS for 3D Horizontal Velocity
!               Written by Colton J. Conroy
!                     @ The C.H.I.L
!                         2.18.14
!
!------------------------------------------------------------------
SUBROUTINE RHS_HU3D(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR, nfacesZ, H_Elem, nfacesX, nfacesY, hL, hR
INTEGER(int_p) :: REGION, check
INTEGER(int_p), PARAMETER :: flux_type = 1
REAL(real_p) :: INVJAC, xpt, ypt, zpt, SE, grad_SE, Vu, Vavg, g, h0, s, hpt, nn, k
REAL(real_p) :: Vv, Vw, DVw, Huse, Zb, DSEx, DSEy, Uavg, Fpce, Fxm, Fym, Fw, solnT
REAL(real_p) :: freq, Lx, FxmA, FymA, dxh, dHdX, ampU
REAL(real_p), DIMENSION(NPTS_2D) :: EL, ER, ZL, ZR, WL, WR
REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, QL, QR, Ghat, Fhat, FL, FR, Kh, Nh
REAL(real_p), DIMENSION(NPTS_2D) :: ZWL, ZWR, HUL, HUR, EIGMAX, Hhj, FQhat
REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Gj, Qj, Saj, HUj, HVj, Ubj, Vbj, Vj, Zj
REAL(real_p), DIMENSION(NPTS_3D) :: HUxj, HVyj, Hj, Nj, Kj, Wj, Fxj, Fyj, Fzj
REAL(real_p), DIMENSION(NDOF_3D) :: Sj, Rzj, Rxj, Ryj
REAL(real_p), DIMENSION(L2_2D) :: L2SE, Kuse, Nuse, L2H
REAL(real_p), DIMENSION(L2_3D) :: dZETAj, L2V_3D, L2W_3D, HL2j, FL2j, Ej, DZj
REAL(real_p), DIMENSION(L2_3D,HEXelems) :: L2Fx
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode
REAL(real_p), PARAMETER :: half = 0.50d0

nfacesZ = size(zFACE,1)
nfacesX = size(xFACE,1)
nfacesY = size(yFACE,1)
solnT   = REAL(TIME)
freq    = PROB_DATA%freq
Lx      = PROB_DATA%xN - PROB_DATA%x0
check   = 0
ampU    = PROB_DATA%ampU

EL = 0.0d0
ER = 0.0d0
ZL = 0.0d0
ZR = 0.0d0
WL = 0.0d0
WR = 0.0d0
UL = 0.0d0
UR = 0.0d0
QL = 0.0d0
QR = 0.0d0
FL = 0.0d0
FR = 0.0d0
Kh = 0.0d0
Nh = 0.0d0
ZWL = 0.0d0
ZWR = 0.0d0
HUL = 0.0d0
HUR = 0.0d0
Hhj = 0.0d0
Uj  = 0.0d0
Gj  = 0.0d0
Qj  = 0.0d0
Saj = 0.0d0
HUj = 0.0d0
HVj = 0.0d0
Ubj = 0.0d0
Vbj = 0.0d0
Vj  = 0.0d0
Zj  = 0.0d0
HUxj = 0.0d0
HVyj = 0.0d0
Hj   = 0.0d0
Nj   = 0.0d0
Kj   = 0.0d0
Wj   = 0.0d0
Fxj  = 0.0d0
Fyj  = 0.0d0
Fzj  = 0.0d0
Sj   = 0.0d0 
Rzj  = 0.0d0 
Rxj  = 0.0d0 
Ryj  = 0.0d0
L2SE = 0.0d0 
Kuse = 0.0d0 
Nuse = 0.0d0 
L2H  = 0.0d0
dZETAj = 0.0d0 
L2V_3D = 0.0d0 
L2W_3D = 0.0d0 
HL2j   = 0.0d0 
FL2j   = 0.0d0 
Ej     = 0.0d0 
DZj    = 0.0d0
Ghat   = 0.0d0
Fhat   = 0.0d0
FQhat  = 0.0d0
EIGMAX = 0.0d0
L2Fx   = 0.0d0
Xnode  = 0.0d0
Znode  = 0.0d0
Ynode  = 0.0d0
globalNODES    = 0.0d0
globalNODES_3D = 0.0d0

!....Debug check

IF (PROB == 3) THEN

    DO i = 1, HEXelems
        
        globalNODES_3D = HEXconn(i,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, L2_3D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,1))*Xnode(1)+(1.0d0+L2PTS_3D(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,2))*Ynode(1)+(1.0d0+L2PTS_3D(j,2))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,3))*Znode(1)+(1.0d0+L2PTS_3D(j,3))*Znode(2))
            CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
            L2V_3D(j) = Huse*Vv
            L2W_3D(j) = Vw
            L2Fx(j,i) = Fxm
                                         
        ENDDO       
        
        IF (coupled == 1) THEN
            Hv(:,i,irk) = MATMUL(L2_U3D,L2V_3D)  
            W_3D(:,i)   = MATMUL(L2_U3D,L2W_3D)  
        ENDIF    
        
    ENDDO
    
    IF (coupled == 1) THEN
        DO i = 1, Qelems
    
            globalNODES = CONN(i,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
       
            DO j = 1, L2_2D
        
                xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
                zpt = Z(nelemsZ)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2SE(j)   = Huse
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2H(j)    = Zb     
                             
            ENDDO
                                          
            Hh(:,i)       = MATMUL(L2Bh,L2H)
            NX(:,i)       = MATMUL(L2U,Nuse)
            KX(:,i)       = MATMUL(L2U,Kuse)   
            ZETA(:,i,irk) = MATMUL(L2U,L2SE)
            
        ENDDO     
    ENDIF    
          
ENDIF

IF (coupled == 1) THEN
    REGION = 3
    CALL CALC_SE(REGION)
    CALL CALC_U(REGION)
    CALL CALC_V(REGION)
ENDIF

!....Elements for U          

g      = PROB_DATA%g
h0     = PROB_DATA%h0

DO i = 1, HEXelems
        
    !....Elem calc for Hu
    
    H_Elem  = HEXzxy(i,1)
    Nj      = MATMUL(PHI_2D3D,NX(:,H_Elem))
    Uj      = MATMUL(PHI_3D,U_3D(:,i,irk))
    HUj     = MATMUL(PHI_3D,Hu(:,i,irk))
    Vj      = MATMUL(PHI_3D,V_3D(:,i,irk))
    Hj      = MATMUL(PHI_2D3D,ZETA(:,H_Elem,irk))
    Ej      = MATMUL(PHIzeta_3D,ZETA_2D(:,H_Elem,irk))
    Wj      = MATMUL(PHI_3Dw,W_3D(:,i))
    !Zj      = MATMUL(PHI_2D3D,Hh(:,H_Elem))
    Zj      = MATMUL(PHI_2D3Dh,Hh(:,H_Elem))
    !Qj      = MATMUL(PHI_3D,Q(:,i))
    
    IF (p == 0) THEN
        CALL dHdX_P0(H_Elem,dHdX,solnT)
        DZj = dHdX
    ELSE
        DZj = MATMUL(DPHIzeta_3Dh(:,:,1),Hh(:,H_Elem))*CONN_jacX(H_Elem,1)   
    ENDIF     
    
    !....Sigma Contribution
    
    Fzj   = sigE(Uj,Wj)
    Rzj   = MATMUL(A_3D(:,:,3),Fzj)*HEX_jacZ(i,1)
    
    !....X Contribution
  
    Fxj = xmE(HUj,Uj,Hj,Zj,g)
    Rxj  = MATMUL(A_3D(:,:,1),Fxj)*HEX_jacX(i,1)
        
    !....Y Contribution
    
    Fyj = xymE(Huj,Vj)
    Ryj = MATMUL(A_3D(:,:,2),Fyj)*HEX_jacY(i,1)
    
    !....Source Contribution
    
    Saj  = xmSE(DZj,Ej,g)
    Sj   = MATMUL(C_3D,Saj) + MATMUL(C_3D,L2Fx(:,i))
    
    !....Add it all up
    
    IF (PROB == 4) THEN
        Sj = 0.0d0
    ENDIF
    
    RHSu_3D(:,i,irk) = Rxj(:) + Ryj(:) + Rzj(:) + Sj(:)
        
    !....Calculations for DQW
!
!    HUj  = MATMUL(PHI_3D,Hu(:,i,irk))
!    HVj  = MATMUL(PHI_3D,Hv(:,i,irk))
!    Ubj  = MATMUL(PHI_2D3D,USE(:,H_Elem,irk)) 
!    Vbj  = MATMUL(PHI_2D3D,VSE(:,H_Elem,irk))  
!    Hj   = MATMUL(PHI_2D3D,ZETA(:,H_Elem,irk))
!    DO j = 1, NPTS_3D
!        HUxj(j) = -Hj(j)*Ubj(j) + Huj(j)
!        HVyj(j) = -Hj(j)*Vbj(j) + Hvj(j)
!    ENDDO
!    Rxj = MATMUL(A_3D(:,:,1),HUxj)*HEX_jacX(i,1)
!    Ryj = MATMUL(A_3D(:,:,2),HVyj)*HEX_jacY(i,1)
!    DQW(:,i) = Rxj + Ryj
 
ENDDO    


DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        H_Elem = HEXzxy(jR,1)
        UL   = MATMUL(bPHI_3D(:,:,3),U_3D(:,jR,irk))
        Kh   = MATMUL(PHIz,KX(:,H_Elem))
        Nh   = MATMUL(PHIz,NX(:,H_Elem))
        Hhj  = MATMUL(PHIz,Hh(:,H_Elem))

        DO j = 1, NPTS_2D
            QL(j)   = -Kh(j)/SQRT(Nh(j))*UL(j)
        ENDDO  
         
        Fhat = sigB(QL,Nh,Hhj) 
        
        IF (PROB == 4) THEN
            Fhat = 0.0d0
        ENDIF

        RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)
        
    ELSEIF (jR == 0) THEN  
    
        globalNODES_3D = HEXconn(jL,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
    
        IF (PROB == 3) THEN
        
            H_Elem = HEXzxy(jL,1)
            Nh     = MATMUL(PHIz,NX(:,H_Elem))
            Hhj    = MATMUL(PHIz,Hh(:,H_Elem)) 
            Kh     = MATMUL(PHIz,KX(:,H_Elem))
            DO j = 1, NPTS_2D
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,4))*Xnode(1)+(1.0d0+bPTS_3D(j,1,4))*Xnode(2))
                QR(j) = sqrt(Hhj(j)/Nh(j))*Kh(j)*exp(ampU*sin(2.0d0*pi*xpt/Lx + freq*solnT))
            ENDDO    
            Fhat  = sigB(QR,Nh,Hhj)
            
        ELSEIF (PROB == 4) THEN
        
            Fhat = 0.0d0
            
        ENDIF       
            
        RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)
        
    ELSE 
    
        H_Elem = HEXzxy(jR,1)      
        
        WR     = MATMUL(bPHI_3Dw(:,:,1),W_3D(:,jR)) !....Need to fix need bPHI_3DW ie basis evaluated for pw +1 
        WL     = MATMUL(bPHI_3Dw(:,:,2),W_3D(:,jL))
        HuL    = MATMUL(bPHI_3D(:,:,4),Hu(:,jL,irk))
        HuR    = MATMUL(bPHI_3D(:,:,3),Hu(:,jR,irk))
        EL     = MATMUL(PHIz,ZETA(:,H_Elem,irk))
        ER     = MATMUL(PHIz,ZETA(:,H_Elem,irk))
        UR     = MATMUL(bPHI_3D(:,:,3),U_3D(:,jR,irk)) 
        UL     = MATMUL(bPHI_3D(:,:,4),U_3D(:,jL,irk))
        EIGMAX = EIGz(WR,WL,EL,ER)
        Fhat   = LLFz(UL,UR,WL,WR,HuL,HuR,EIGmax)
        
        RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)     
        RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)
        
    ENDIF    
        
ENDDO        


DO i = 1, nfacesX
    
    jL = xFACE(i,1)
    jR = xFACE(i,2)
    hL = HEXzxy(jL,1)
    hR = HEXzxy(jR,1)
        
    IF (jL == -1) THEN
        
        IF (PROB == 1) THEN
        
            ER   = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
            Fhat = g*ER(1)
            
        ELSEIF (PROB == 3) THEN
            
            DO j = 1, NPTS_2D
        
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,1))*Xnode(1)+(1.0d0+bPTS_3D(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,1))*Ynode(1)+(1.0d0+bPTS_3D(j,2,1))*Ynode(2))
                zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,1))*Znode(1)+(1.0d0+bPTS_3D(j,3,1))*Znode(2))
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                HuL(j) = Huse*Vu
                UL(j)  = Vu
                EL(j)  = Huse
                ZL(j)  = Zb
                Fhat(j) = HuL(j)*UL(j) + half*g*(EL(j)**2.0d0 - ZL(j)**2.0d0)
                                         
            ENDDO          
                
            HuR = MATMUL(bPHI_3D(:,:,1),Hu(:,jR,irk))
            UR  = MATMUL(bPHI_3D(:,:,1),U_3D(:,jR,irk))
            ER  = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
            ZR  = MATMUL(b_PHI_2D3Dx(:,:,1),Hh(:,hR))  
            EIGMAX = EIG(EL,ER,UL,UR)
            !Fhat = LLFx(EL,ER,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
            
        ELSEIF (PROB == 4) THEN
        
            ER = 10.0d0
            DO j = 1, NPTS_2D
                Fhat(j) = half*g*ER(j)**2.0d0
            ENDDO    
            
        ENDIF    
            
        RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
            
        !....Calculations for DQW
            
!        Fhat = 0.0d0 !....Land BC
!        DQW(:,jR) = DQW(:,jR) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)            
        
    ELSEIF (jR == 0) THEN
        
        IF (PROB == 1) THEN
          
            xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
            hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
            zCOORD = dcmplx(hpt*Z(nelemsZ), 0.0d0)
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            FR   = g*SE
            Fhat = FR   
                
        ELSEIF (PROB == 3) THEN
            
            DO j = 1, NPTS_2D
        
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,2))*Xnode(1)+(1.0d0+bPTS_3D(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,2))*Ynode(1)+(1.0d0+bPTS_3D(j,2,2))*Ynode(2))
                zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,2))*Znode(1)+(1.0d0+bPTS_3D(j,3,2))*Znode(2))
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                HuR(j) = Huse*Vu
                UR(j)  = Vu
                ER(j)  = Huse
                ZR(j)  = Zb
                Fhat(j) = HuR(j)*UR(j) + half*g*(ER(j)**2.0d0-ZR(j)**2.0d0)                     
            ENDDO          
                
            HuL = MATMUL(bPHI_3D(:,:,2),Hu(:,jL,irk))
            UL  = MATMUL(bPHI_3D(:,:,2),U_3D(:,jL,irk))
            EL  = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
            ZL  = MATMUL(b_PHI_2D3Dx(:,:,2),Hh(:,hL))  
            EIGMAX = EIG(EL,ER,UL,UR)
            !Fhat = LLFx(EL,ER,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
            
        ELSEIF (PROB == 4) THEN
        
            EL = 5.0d0
            DO j = 1, NPTS_2D
                Fhat(j) = half*g*EL(j)**2.0d0
            ENDDO              
                
        ENDIF          
            
        RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)        
            
            !....Calculations for DQW
            
                !Fhat = 0.0d0  !....Land BC
                !DQW(:,jL) = DQW(:,jL) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)            
            
    ELSE
        
        HuL = MATMUL(bPHI_3D(:,:,2),Hu(:,jL,irk))
        UL  = MATMUL(bPHI_3D(:,:,2),U_3D(:,jL,irk))
        EL  = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
        !ZL  = MATMUL(b_PHI_2D3Dx(:,:,2),Hh(:,hL))
        ZL  = MATMUL(b_PHI_2D3Dxh(:,:,2),Hh(:,hL))
        HuR = MATMUL(bPHI_3D(:,:,1),Hu(:,jR,irk))
        UR  = MATMUL(bPHI_3D(:,:,1),U_3D(:,jR,irk))
        ER  = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
        !ZR  = MATMUL(b_PHI_2D3Dx(:,:,1),Hh(:,hR))
        ZR  = MATMUL(b_PHI_2D3Dxh(:,:,1),Hh(:,hR))
        IF (flux_type == 1) THEN    
            EIGMAX = EIG(EL,ER,UL,UR)
            Fhat   = LLFx(EL,ER,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
        ELSE
            Fhat = Godunov(ZL,ZR,EL,ER,UL,UR,HuL,HuR)    
        ENDIF  
      
        RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
        RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)
            
            !....Calculations for DQW
                
!            Fhat = 0.0d0    
!            ZWL  = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
!            UL   = MATMUL(b_PHI_2D3Dx(:,:,2),USE(:,hL,irk))
!            ZWR  = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
!            UR   = MATMUL(b_PHI_2D3Dx(:,:,1),USE(:,hR,irk))
!            HUL  = MATMUL(bPHI_3D(:,:,2),Hu(:,jL,irk))
!            HUR  = MATMUL(bPHI_3D(:,:,1),Hu(:,jR,irk))
!            DO j = 1, NPTS_2D
!                FR(j)   = -ZWR(j)*UR(j)+HUR(j)
!                FL(j)   = -ZWL(j)*UL(j)+HUL(j)
!                Fhat(j) = half*(FR(j)+FL(j)) 
!            ENDDO    
!            DQW(:,jR) = DQW(:,jR) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
!            DQW(:,jL) = DQW(:,jL) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)            
            
    ENDIF
        
ENDDO    
   
!....Loop over y-faces

DO i = 1, nfacesY    !....The max eigenvalues in this case will be max(v + sqrt(gH),v - sqrt(gH),v) Not u!
    
    jL = yFACE(i,1)
    jR = yFACE(i,2)
    hL = HEXzxy(jL,1)
    hR = HEXzxy(jR,1)
        
    IF (jL == -1) THEN
    
            IF (PROB == 3) THEN
            
                DO j = 1, NPTS_2D
        
                    xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,5))*Xnode(1)+(1.0d0+bPTS_3D(j,1,5))*Xnode(2))
                    ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,5))*Ynode(1)+(1.0d0+bPTS_3D(j,2,5))*Ynode(2))
                    zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,5))*Znode(1)+(1.0d0+bPTS_3D(j,3,5))*Znode(2))
                    CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                    HuL(j) = Huse*Vu
                    UL(j)  = Vv
                    EL(j)  = Huse
                    ZL(j)  = Zb
                    Fhat(j) = HuL(j)*UL(j)                     
                ENDDO          
                
                HuR = MATMUL(bPHI_3D(:,:,5),Hu(:,jR,irk))
                UR  = MATMUL(bPHI_3D(:,:,5),V_3D(:,jR,irk))
                ER  = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
                ZR  = MATMUL(b_PHI_2D3Dy(:,:,1),Hh(:,hR))  
                EIGmax = EIG(EL,ER,UL,UR)
                !Fhat = LLFy(UL,UR,HuL,HuR,EIGmax)
                
            ELSEIF (PROB == 4) THEN
            
                Fhat = 0.0d0    
                
            ENDIF    
            
            RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)    
    
        !....DQW
        
!        Fhat = 0.0d0  !....Land BC       
!        DQW(:,jR) = DQW(:,jR) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
        
    ELSEIF (jR == 0) THEN
    
            IF (PROB == 3) THEN
            
                DO j = 1, NPTS_2D
        
                    xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,6))*Xnode(1)+(1.0d0+bPTS_3D(j,1,6))*Xnode(2))
                    ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,6))*Ynode(1)+(1.0d0+bPTS_3D(j,2,6))*Ynode(2))
                    zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,6))*Znode(1)+(1.0d0+bPTS_3D(j,3,6))*Znode(2))
                    CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                    HuR(j) = Huse*Vu
                    UR(j)  = Vv
                    ER(j)  = Huse
                    ZR(j)  = Zb
                    Fhat(j) = HuR(j)*UR(j)                     
                ENDDO          
                
                HuL = MATMUL(bPHI_3D(:,:,6),Hu(:,jL,irk))
                UL  = MATMUL(bPHI_3D(:,:,6),V_3D(:,jL,irk))
                EL  = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
                ZL  = MATMUL(b_PHI_2D3Dy(:,:,2),Hh(:,hL))  
                EIGmax = EIG(EL,ER,UL,UR)
                !Fhat = LLFy(UL,UR,HuL,HuR,EIGmax)
                
            ELSEIF (PROB == 4) THEN
            
                Fhat = 0.0d0    
                
            ENDIF    
            
            RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat) 
    
        !....DQW 
             
!        Fhat = 0.0d0  !....Land BC
!        DQW(:,jL) = DQW(:,jL) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
    ELSE
    
            HuL = MATMUL(bPHI_3D(:,:,6),Hu(:,jL,irk))
            UL  = MATMUL(bPHI_3D(:,:,6),V_3D(:,jL,irk))
            EL  = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
            HuR = MATMUL(bPHI_3D(:,:,5),Hu(:,jR,irk))
            UR  = MATMUL(bPHI_3D(:,:,5),V_3D(:,jR,irk))
            ER  = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
            
            EIGmax = EIG(EL,ER,UL,UR)
            Fhat   = LLFy(UL,UR,HuL,HuR,EIGmax)
            
            RHSu_3D(:,jR,irk) = RHSu_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
            RHSu_3D(:,jL,irk) = RHSu_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)    

        !....DQW
        
!        ZWL  = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
!        UL   = MATMUL(b_PHI_2D3Dy(:,:,2),VSE(:,hL,irk))
!        ZWR  = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
!        UR   = MATMUL(b_PHI_2D3Dy(:,:,1),VSE(:,hR,irk))
!        HUL  = MATMUL(bPHI_3D(:,:,6),Hv(:,jL,irk))
!        HUR  = MATMUL(bPHI_3D(:,:,5),Hv(:,jR,irk))
!        DO j = 1, NPTS_2D
!            FR(j)   = -ZWR(j)*UR(j)+HUR(j)
!            FL(j)   = -ZWL(j)*UL(j)+HUL(j)
!            Fhat(j) = half*(FR(j)+FL(j)) 
!        ENDDO 
!        DQW(:,jR) = DQW(:,jR) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
!        DQW(:,jL) = DQW(:,jL) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
    ENDIF
    
ENDDO 

CONTAINS

    FUNCTION Godunov(ZL,ZR,EL,ER,UL,UR,HuL,HuR)
    INTEGER(int_p) :: j
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, Godunov
    REAL(real_p), DIMENSION(NPTS_2D) :: HuL, HuR, EL, ER, ZL, ZR
    REAL(real_p) :: s, DEN
        DO j = 1, NPTS_2D
            FL(j) = HuL(j)*UL(j) + (1.0d0/2.0d0)*g*(EL(j)**2 - ZL(j)**2) ! for dam break h**2 = 0 b/c flat bottom
            FR(j) = HuR(j)*UR(j) + (1.0d0/2.0d0)*g*(ER(j)**2 - ZR(j)**2) !
            DEN = (HuR(j) - HuL(j))
            IF (abs(DEN) < 1d-16) THEN
                s = 0.0d0
            ELSE    
                s = (FR(j) - FL(j))/DEN
            ENDIF    
            IF (s > 0.0d0) THEN
                Godunov(j) = FL(j)
                IF (UL(j) < 0 .and. UR(j) > 0) THEN  ! Transonic case
                    Godunov(j) = (1.0d0/2.0d0)*g*(EL(j)**2 - ZL(j)**2)
                ENDIF 
            ELSE    
                Godunov(j) = FR(j)
                IF (UL(j) < 0 .and. UR(j) > 0) THEN  ! Transonic case
                    Godunov(j) = (1.0d0/2.0d0)*g*(ER(j)**2 - ZR(j)**2)
                ENDIF 
            ENDIF       
    
        ENDDO
    END FUNCTION Godunov

    FUNCTION LLFx(HL,HR,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, LLFx
    REAL(real_p), DIMENSION(NPTS_2D) :: HuL, HuR, EIGmax
    REAL(real_p), DIMENSION(NPTS_2D) :: HL, HR, ZL, ZR
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_2D
            FL(j)   = HuL(j)*UL(j) + half*g*(HL(j)**2.0d0 - ZL(j)**2.0d0) 
            FR(j)   = HuR(j)*UR(j) + half*g*(HR(j)**2.0d0 - ZR(j)**2.0d0) 
            LLFx(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFx
    
    FUNCTION LLFy(UL,UR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, LLFy
    REAL(real_p), DIMENSION(NPTS_2D) :: HuL, HuR, EIGmax
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_2D
            FL(j)   = HuL(j)*UL(j) 
            FR(j)   = HuR(j)*UR(j)  
            LLFy(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFy    
    
    FUNCTION LLFz(UL,UR,WL,WR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, LLFz
    REAL(real_p), DIMENSION(NPTS_2D) :: WL, WR, HuL, HuR, EIGmax
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_2D
            FL(j)   = -WL(j)*UL(j) 
            FR(j)   = -WR(j)*UR(j)  
            LLFz(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFz 

    FUNCTION EIG(HL,HR,UL,UR)
    REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, HL, HR, EIG
    REAL(real_p) :: LEp, LEm, REp, REm
        DO j = 1, NPTS_2D
            LEp    = UL(j) + sqrt(g*HL(j))
            LEm    = UL(j) - sqrt(g*HL(j))
            REp    = UR(j) + sqrt(g*HR(j))
            REm    = UR(j) - sqrt(g*HR(j))
            EIG(j) = max(abs(LEp),abs(LEm))
            EIG(j) = max(abs(EIG(j)),abs(REp))
            EIG(j) = max(abs(EIG(j)),abs(REm))
            EIG(j) = max(abs(EIG(j)),abs(UL(j)))
            EIG(j) = max(abs(EIG(j)),abs(UR(j)))
        ENDDO    
    END FUNCTION EIG  
    
    FUNCTION EIGz(WR,WL,EL,ER)
    REAL(real_p), DIMENSION(NPTS_2D) :: WR, WL, EL, ER, EIGz
    REAL(real_p) :: LE, RE
        DO j = 1, NPTS_2D
            LE = WL(j)/EL(j)
            RE = WR(j)/ER(j)
            EIGz(j) = max(abs(LE),abs(RE))
        ENDDO    
    END FUNCTION EIGz

    FUNCTION gfB(Uj,Nj,Hj)
	REAL(real_p), DIMENSION(NPTS_2D) :: Uj, Nj, Hj, gfB
	    DO j = 1, NPTS_2D
	        gfB = sqrt(Nj(j)/Hj(j))*Uj(j)
	    ENDDO    
    END FUNCTION gfB
    
    FUNCTION gfE(Uj,Nj,Hj)
    REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Nj, Hj, gfE
        DO j = 1, NPTS_3D
            gfE(j) = sqrt(Nj(j)/Hj(j))*Uj(j)
        ENDDO    
    END FUNCTION gfE
    
    FUNCTION sigE(Uj,Wj) 
    REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Wj, sigE
        DO j = 1, NPTS_3D
            sigE(j) = -Uj(j)*Wj(j)
        ENDDO    
    END FUNCTION sigE
    
    FUNCTION sigB(Qj,Nj,Hj) 
    REAL(real_p), DIMENSION(NPTS_2D) :: Qj, Nj, Hj, sigB
        DO j = 1, NPTS_2D
            sigB(j) = sqrt(Nj(j)/Hj(j))*Qj(j) 
        ENDDO    
    END FUNCTION sigB    
    
    FUNCTION xmE(HUj,Uj,Hj,Zj,g)
    REAL(real_p) :: g
    REAL(real_p), DIMENSION(NPTS_3D) :: Huj, Uj, Hj, Zj, xmE
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_3D
            xmE(j) = Huj(j)*Uj(j) + half*g*(Hj(j)**2.0d0 - Zj(j)**2.0d0)
        ENDDO
    END FUNCTION xmE   
    
    FUNCTION xymE(Huj,Vj)
    REAL(real_p), DIMENSION(NPTS_3D) :: Huj, Vj, xymE     
        DO j = 1, NPTS_3D
            xymE(j) = Huj(j)*Vj(j)
        ENDDO
    END FUNCTION xymE   
    
    FUNCTION xmSE(DZj,Ej,g)
    REAL(real_p) :: g
    REAL(real_p), DIMENSION(L2_3D) :: xmSE, DZj, Ej 
        DO j = 1, L2_3D
            xmSE(j) = g*Ej(j)*DZj(j)
        ENDDO
    END FUNCTION xmSE    
            
END SUBROUTINE RHS_HU3D
!------------------------------------------------------------------
!
!             RHS for 3D Horizontal Velocity (y-mom)
!               Written by Colton J. Conroy
!                     @ The C.H.I.L
!                         2.18.14
!
!------------------------------------------------------------------
SUBROUTINE RHS_HV3D(TIME)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR, nfacesZ, H_Elem, nfacesX, nfacesY, hL, hR
INTEGER(int_p) :: REGION, check
INTEGER(int_p), PARAMETER :: flux_type = 2
REAL(real_p) :: INVJAC, xpt, ypt, zpt, SE, grad_SE, Vu, Vavg, g, h0, s, hpt, nn, k
REAL(real_p) :: Vv, Vw, DVw, Huse, Zb, DSEx, DSEy, Uavg, Fpce, Fxm, Fym, Fw, solnT
REAL(real_p) :: freq, Lx, FxmA, FymA, dxh, dHdX, ampU, dHdY, Ly
REAL(real_p), DIMENSION(NPTS_2D) :: EL, ER, ZL, ZR, WL, WR
REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, QL, QR, Ghat, Fhat, FL, FR, Kh, Nh
REAL(real_p), DIMENSION(NPTS_2D) :: ZWL, ZWR, HUL, HUR, EIGMAX, Hhj, FQhat
REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Gj, Qj, Saj, HUj, HVj, Ubj, Vbj, Vj, Zj
REAL(real_p), DIMENSION(NPTS_3D) :: HUxj, HVyj, Hj, Nj, Kj, Wj, Fxj, Fyj, Fzj
REAL(real_p), DIMENSION(NDOF_3D) :: Sj, Rzj, Rxj, Ryj
REAL(real_p), DIMENSION(L2_2D) :: L2SE, Kuse, Nuse, L2H
REAL(real_p), DIMENSION(L2_3D) :: dZETAj, L2V_3D, L2W_3D, HL2j, FL2j, Ej, DZj
REAL(real_p), DIMENSION(L2_3D,HEXelems) :: L2Fy
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode
REAL(real_p), PARAMETER :: half = 0.50d0

nfacesZ = size(zFACE,1)
nfacesX = size(xFACE,1)
nfacesY = size(yFACE,1)
solnT   = REAL(TIME)
freq    = PROB_DATA%freq
Ly      = PROB_DATA%yN - PROB_DATA%y0
check   = 0
ampU    = PROB_DATA%ampU

EL = 0.0d0
ER = 0.0d0
ZL = 0.0d0
ZR = 0.0d0
WL = 0.0d0
WR = 0.0d0
UL = 0.0d0
UR = 0.0d0
QL = 0.0d0
QR = 0.0d0
FL = 0.0d0
FR = 0.0d0
Kh = 0.0d0
Nh = 0.0d0
ZWL = 0.0d0
ZWR = 0.0d0
HUL = 0.0d0
HUR = 0.0d0
Hhj = 0.0d0
Uj  = 0.0d0
Gj  = 0.0d0
Qj  = 0.0d0
Saj = 0.0d0
HUj = 0.0d0
HVj = 0.0d0
Ubj = 0.0d0
Vbj = 0.0d0
Vj  = 0.0d0
Zj  = 0.0d0
HUxj = 0.0d0
HVyj = 0.0d0
Hj   = 0.0d0
Nj   = 0.0d0
Kj   = 0.0d0
Wj   = 0.0d0
Fxj  = 0.0d0
Fyj  = 0.0d0
Fzj  = 0.0d0
Sj   = 0.0d0 
Rzj  = 0.0d0 
Rxj  = 0.0d0 
Ryj  = 0.0d0
L2SE = 0.0d0 
Kuse = 0.0d0 
Nuse = 0.0d0 
L2H  = 0.0d0
dZETAj = 0.0d0 
L2V_3D = 0.0d0 
L2W_3D = 0.0d0 
HL2j   = 0.0d0 
FL2j   = 0.0d0 
Ej     = 0.0d0 
DZj    = 0.0d0
Ghat   = 0.0d0
Fhat   = 0.0d0
FQhat  = 0.0d0
EIGMAX = 0.0d0
L2Fy   = 0.0d0
Xnode  = 0.0d0
Znode  = 0.0d0
Ynode  = 0.0d0
globalNODES    = 0.0d0
globalNODES_3D = 0.0d0

!....Debug check

IF (PROB == 3) THEN

    DO i = 1, HEXelems
        
        globalNODES_3D = HEXconn(i,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, L2_3D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,1))*Xnode(1)+(1.0d0+L2PTS_3D(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,2))*Ynode(1)+(1.0d0+L2PTS_3D(j,2))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,3))*Znode(1)+(1.0d0+L2PTS_3D(j,3))*Znode(2))
            CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
            L2V_3D(j) = Huse*Vu
            L2W_3D(j) = Vw
            L2Fy(j,i) = Fym
                                         
        ENDDO       
        
        IF (coupled == 1) THEN
            Hu(:,i,irk) = MATMUL(L2_U3D,L2V_3D)  
            W_3D(:,i)   = MATMUL(L2_U3D,L2W_3D)  !....would need to change if pw =/ pu
        ENDIF    
        
    ENDDO
    
    IF (coupled == 1) THEN
        DO i = 1, Qelems
    
            globalNODES = CONN(i,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
       
            DO j = 1, L2_2D
        
                xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
                zpt = Z(nelemsZ)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2SE(j)   = Huse
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2H(j)    = Zb     
                             
            ENDDO
                                          
            Hh(:,i)       = MATMUL(L2Bh,L2H)
            NX(:,i)       = MATMUL(L2U,Nuse)
            KX(:,i)       = MATMUL(L2U,Kuse)   
            ZETA(:,i,irk) = MATMUL(L2U,L2SE)
            
        ENDDO     
    ENDIF
          
ENDIF

IF (coupled == 1) THEN
    REGION = 3
    CALL CALC_SE(REGION)
    CALL CALC_U(REGION)
    CALL CALC_V(REGION)
ENDIF

!....Elements for U          

g      = PROB_DATA%g
h0     = PROB_DATA%h0

DO i = 1, HEXelems
        
    !....Elem calc for Hu
    
    H_Elem  = HEXzxy(i,1)
    Nj      = MATMUL(PHI_2D3D,NX(:,H_Elem))
    Uj      = MATMUL(PHI_3D,U_3D(:,i,irk))
    HVj     = MATMUL(PHI_3D,Hv(:,i,irk))
    Vj      = MATMUL(PHI_3D,V_3D(:,i,irk))
    Hj      = MATMUL(PHI_2D3D,ZETA(:,H_Elem,irk))
    Ej      = MATMUL(PHIzeta_3D,ZETA_2D(:,H_Elem,irk))
    Wj      = MATMUL(PHI_3Dw,W_3D(:,i))
    Zj      = MATMUL(PHI_2D3Dh,Hh(:,H_Elem))
    !Qj      = MATMUL(PHI_3D,Q(:,i))
    
    IF (p == 0) THEN
        CALL dHdY_P0(H_Elem,dHdY,solnT) 
        DZj = dHdY
    ELSE
        DZj = MATMUL(DPHIzeta_3Dh(:,:,2),Hh(:,H_Elem))*CONN_jacY(H_Elem,1)   
    ENDIF     
    
    !....Sigma Contribution
    
    Fzj   = sigE(Vj,Wj)
    Rzj   = MATMUL(A_3D(:,:,3),Fzj)*HEX_jacZ(i,1)
    
    !....Y Contribution
    
    Fyj = ymE(HVj,Vj,Hj,Zj,g)
    Ryj  = MATMUL(A_3D(:,:,2),Fyj)*HEX_jacY(i,1)
        
    !....X Contribution
    
    Fxj = xymE(Hvj,Uj)
    Rxj = MATMUL(A_3D(:,:,1),Fxj)*HEX_jacX(i,1)
    
    !....Source Contribution
    
    Saj  = ymSE(DZj,Ej,g)
    Sj   = MATMUL(C_3D,Saj) + MATMUL(C_3D,L2Fy(:,i))
    
    !....Add it all up
    
    IF (PROB == 4) THEN
        Sj = 0.0d0
    ENDIF
    
    RHSv_3D(:,i,irk) = Rxj(:) + Ryj(:) + Rzj(:) + Sj(:)
    
    !....Calculations for DQW
!
!    HUj  = MATMUL(PHI_3D,Hu(:,i,irk))
!    HVj  = MATMUL(PHI_3D,Hv(:,i,irk))
!    Ubj  = MATMUL(PHI_2D3D,USE(:,H_Elem,irk)) 
!    Vbj  = MATMUL(PHI_2D3D,VSE(:,H_Elem,irk))  
!    Hj   = MATMUL(PHI_2D3D,ZETA(:,H_Elem,irk))
!    DO j = 1, NPTS_3D
!        HUxj(j) = -Hj(j)*Ubj(j) + Huj(j)
!        HVyj(j) = -Hj(j)*Vbj(j) + Hvj(j)
!    ENDDO
!    Rxj = MATMUL(A_3D(:,:,1),HUxj)*HEX_jacX(i,1)
!    Ryj = MATMUL(A_3D(:,:,2),HVyj)*HEX_jacY(i,1)
!    DQW(:,i) = Rxj + Ryj
 
ENDDO    


DO i = 1, nfacesZ

    jL = zFACE(i,1)
    jR = zFACE(i,2)
    
    IF (jL == -1) THEN
        
        H_Elem = HEXzxy(jR,1)
        UL   = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk))
        Kh   = MATMUL(PHIz,KX(:,H_Elem))
        Nh   = MATMUL(PHIz,NX(:,H_Elem))
        Hhj  = MATMUL(PHI_h,Hh(:,H_Elem))

        DO j = 1, NPTS_2D
            QL(j)   = -Kh(j)/SQRT(Nh(j))*UL(j)
        ENDDO  
         
        Fhat = sigB(QL,Nh,Hhj) 
        
        IF (PROB == 4) THEN
            Fhat = 0.0d0
        ENDIF

        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)
        
    ELSEIF (jR == 0) THEN  
    
        globalNODES_3D = HEXconn(jL,:)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
    
        IF (PROB == 3) THEN
        
            H_Elem = HEXzxy(jL,1)
            Nh     = MATMUL(PHIz,NX(:,H_Elem))
            Hhj    = MATMUL(PHI_h,Hh(:,H_Elem)) 
            Kh     = MATMUL(PHIz,KX(:,H_Elem))
            DO j = 1, NPTS_2D
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,4))*Ynode(1)+(1.0d0+bPTS_3D(j,2,4))*Ynode(2))
                QR(j) = sqrt(Hhj(j)/Nh(j))*Kh(j)*exp(ampU*sin(2.0d0*pi*ypt/Ly + freq*solnT))
            ENDDO    
            Fhat  = sigB(QR,Nh,Hhj)
            
        ELSEIF (PROB == 4) THEN
        
            Fhat = 0.0d0
            
        ENDIF       
            
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)
        
    ELSE 
    
        H_Elem = HEXzxy(jR,1)      
        
        WR     = MATMUL(bPHI_3Dw(:,:,1),W_3D(:,jR)) 
        WL     = MATMUL(bPHI_3Dw(:,:,2),W_3D(:,jL))
        HuL    = MATMUL(bPHI_3D(:,:,4),Hv(:,jL,irk))
        HuR    = MATMUL(bPHI_3D(:,:,3),Hv(:,jR,irk))
        EL     = MATMUL(PHIz,ZETA(:,H_Elem,irk))
        ER     = MATMUL(PHIz,ZETA(:,H_Elem,irk))
        UR     = MATMUL(bPHI_3D(:,:,3),V_3D(:,jR,irk)) 
        UL     = MATMUL(bPHI_3D(:,:,4),V_3D(:,jL,irk))
        EIGMAX = EIGz(WR,WL,EL,ER)
        Fhat   = LLFz(UL,UR,WL,WR,HuL,HuR,EIGmax)
        
        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacZ(jR,1)*MATMUL(B_3D(:,:,3),Fhat)     
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacZ(jL,1)*MATMUL(B_3D(:,:,4),Fhat)
        
    ENDIF    
        
ENDDO        


DO i = 1, nfacesX
    
    jL = xFACE(i,1)
    jR = xFACE(i,2)
    hL = HEXzxy(jL,1)
    hR = HEXzxy(jR,1)
        
    IF (jL == -1) THEN
        
        IF (PROB == 1) THEN
        
            ER   = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
            Fhat = g*ER(1)
            
        ELSEIF (PROB == 3) THEN
            
            DO j = 1, NPTS_2D
        
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,1))*Xnode(1)+(1.0d0+bPTS_3D(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,1))*Ynode(1)+(1.0d0+bPTS_3D(j,2,1))*Ynode(2))
                zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,1))*Znode(1)+(1.0d0+bPTS_3D(j,3,1))*Znode(2))
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                HuL(j) = Huse*Vv
                UL(j)  = Vu
                EL(j)  = Huse
                ZL(j)  = Zb
                Fhat(j) = HuL(j)*UL(j) 
                                         
            ENDDO          
                
            HuR = MATMUL(bPHI_3D(:,:,1),Hv(:,jR,irk))
            UR  = MATMUL(bPHI_3D(:,:,1),U_3D(:,jR,irk))
            ER  = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
            ZR  = MATMUL(b_PHI_2D3Dxh(:,:,1),Hh(:,hR))  
            EIGMAX = EIG(EL,ER,UL,UR)
            !Fhat = LLFx(EL,ER,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
            
        ELSEIF (PROB == 4) THEN
        
            Fhat = 0.0d0   
            
        ENDIF    
            
        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
            
        !....Calculations for DQW
            
!        Fhat = 0.0d0 !....Land BC
!        DQW(:,jR) = DQW(:,jR) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)            
        
    ELSEIF (jR == 0) THEN
        
        IF (PROB == 1) THEN
          
            xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
            hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
            zCOORD = dcmplx(hpt*Z(nelemsZ), 0.0d0)
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            FR   = g*SE
            Fhat = FR   
                
        ELSEIF (PROB == 3) THEN
            
            DO j = 1, NPTS_2D
        
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,2))*Xnode(1)+(1.0d0+bPTS_3D(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,2))*Ynode(1)+(1.0d0+bPTS_3D(j,2,2))*Ynode(2))
                zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,2))*Znode(1)+(1.0d0+bPTS_3D(j,3,2))*Znode(2))
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                HuR(j) = Huse*Vv
                UR(j)  = Vu
                ER(j)  = Huse
                ZR(j)  = Zb
                Fhat(j) = HuR(j)*UR(j)                  
            ENDDO          
                
            HuL = MATMUL(bPHI_3D(:,:,2),Hv(:,jL,irk))
            UL  = MATMUL(bPHI_3D(:,:,2),U_3D(:,jL,irk))
            EL  = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
            ZL  = MATMUL(b_PHI_2D3Dxh(:,:,2),Hh(:,hL))  
            EIGMAX = EIG(EL,ER,UL,UR)
            !Fhat = LLFx(EL,ER,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
            
        ELSEIF (PROB == 4) THEN
        
            Fhat = 0.0d0              
                
        ENDIF          
            
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)        
            
            !....Calculations for DQW
            
                !Fhat = 0.0d0  !....Land BC
                !DQW(:,jL) = DQW(:,jL) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)            
            
    ELSE
        
        HuL = MATMUL(bPHI_3D(:,:,2),Hv(:,jL,irk))
        UL  = MATMUL(bPHI_3D(:,:,2),U_3D(:,jL,irk))
        EL  = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
        !ZL  = MATMUL(b_PHI_2D3Dx(:,:,2),Hh(:,hL))
        HuR = MATMUL(bPHI_3D(:,:,1),Hv(:,jR,irk))
        UR  = MATMUL(bPHI_3D(:,:,1),U_3D(:,jR,irk))
        ER  = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
        !ZR  = MATMUL(b_PHI_2D3Dx(:,:,1),Hh(:,hR))
        !IF (flux_type == 1) THEN    
            EIGMAX = EIG(EL,ER,UL,UR)
            Fhat = LLFx(UL,UR,HuL,HuR,EIGmax)
        !ELSE
            !Fhat = Godunov(EL,ER,UL,UR,HuL,HuR)    
        !ENDIF    
               
        RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
        RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)
            
            !....Calculations for DQW
                
!            Fhat = 0.0d0    
!            ZWL  = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
!            UL   = MATMUL(b_PHI_2D3Dx(:,:,2),USE(:,hL,irk))
!            ZWR  = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
!            UR   = MATMUL(b_PHI_2D3Dx(:,:,1),USE(:,hR,irk))
!            HUL  = MATMUL(bPHI_3D(:,:,2),Hu(:,jL,irk))
!            HUR  = MATMUL(bPHI_3D(:,:,1),Hu(:,jR,irk))
!            DO j = 1, NPTS_2D
!                FR(j)   = -ZWR(j)*UR(j)+HUR(j)
!                FL(j)   = -ZWL(j)*UL(j)+HUL(j)
!                Fhat(j) = half*(FR(j)+FL(j)) 
!            ENDDO    
!            DQW(:,jR) = DQW(:,jR) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
!            DQW(:,jL) = DQW(:,jL) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)            
            
    ENDIF
        
ENDDO    
   
!....Loop over y-faces

DO i = 1, nfacesY    !....The max eigenvalues in this case will be max(v + sqrt(gH),v - sqrt(gH),v) Not u!
    
    jL = yFACE(i,1)
    jR = yFACE(i,2)
    hL = HEXzxy(jL,1)
    hR = HEXzxy(jR,1)
        
    IF (jL == -1) THEN
    
            IF (PROB == 3) THEN
            
                DO j = 1, NPTS_2D
        
                    xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,5))*Xnode(1)+(1.0d0+bPTS_3D(j,1,5))*Xnode(2))
                    ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,5))*Ynode(1)+(1.0d0+bPTS_3D(j,2,5))*Ynode(2))
                    zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,5))*Znode(1)+(1.0d0+bPTS_3D(j,3,5))*Znode(2))
                    CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                    HuL(j) = Huse*Vv
                    UL(j)  = Vv
                    EL(j)  = Huse
                    ZL(j)  = Zb
                    Fhat(j) = HuL(j)*UL(j) + half*g*(EL(j)**2.0d0 - ZL(j)**2.0d0)                    
                ENDDO          
                
                HuR = MATMUL(bPHI_3D(:,:,5),Hv(:,jR,irk))
                UR  = MATMUL(bPHI_3D(:,:,5),V_3D(:,jR,irk))
                ER  = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
                ZR  = MATMUL(b_PHI_2D3Dyh(:,:,1),Hh(:,hR))  
                EIGmax = EIG(EL,ER,UL,UR)
                !Fhat = LLFy(UL,UR,HuL,HuR,EIGmax)
                
            ELSEIF (PROB == 4) THEN
            
                Fhat = 0.0d0    
                
            ENDIF    
            
            RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)    
    
        !....DQW
        
!        Fhat = 0.0d0  !....Land BC       
!        DQW(:,jR) = DQW(:,jR) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
        
    ELSEIF (jR == 0) THEN
    
            IF (PROB == 3) THEN
            
                DO j = 1, NPTS_2D
        
                    xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,6))*Xnode(1)+(1.0d0+bPTS_3D(j,1,6))*Xnode(2))
                    ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,6))*Ynode(1)+(1.0d0+bPTS_3D(j,2,6))*Ynode(2))
                    zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,6))*Znode(1)+(1.0d0+bPTS_3D(j,3,6))*Znode(2))
                    CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                    HuR(j) = Huse*Vv
                    UR(j)  = Vv
                    ER(j)  = Huse
                    ZR(j)  = Zb
                    Fhat(j) = HuR(j)*UR(j) + half*g*(ER(j)**2.0d0 - ZR(j)**2.0d0)                     
                ENDDO          
                
                HuL = MATMUL(bPHI_3D(:,:,6),Hv(:,jL,irk))
                UL  = MATMUL(bPHI_3D(:,:,6),V_3D(:,jL,irk))
                EL  = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
                ZL  = MATMUL(b_PHI_2D3Dyh(:,:,2),Hh(:,hL))  
                EIGmax = EIG(EL,ER,UL,UR)
                !Fhat = LLFy(UL,UR,HuL,HuR,EIGmax)
                
            ELSEIF (PROB == 4) THEN
            
                Fhat = 0.0d0    
                
            ENDIF    
            
            RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat) 
    
        !....DQW 
             
!        Fhat = 0.0d0  !....Land BC
!        DQW(:,jL) = DQW(:,jL) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
    ELSE
    
            HuL = MATMUL(bPHI_3D(:,:,6),Hv(:,jL,irk))
            UL  = MATMUL(bPHI_3D(:,:,6),V_3D(:,jL,irk))
            EL  = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
            ZL  = MATMUL(b_PHI_2D3Dyh(:,:,2),Hh(:,hL))
            HuR = MATMUL(bPHI_3D(:,:,5),Hv(:,jR,irk))
            UR  = MATMUL(bPHI_3D(:,:,5),V_3D(:,jR,irk))
            ER  = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
            ZR  = MATMUL(b_PHI_2D3Dyh(:,:,1),Hh(:,hR))
            
            EIGmax = EIG(EL,ER,UL,UR)
            Fhat   = LLFy(EL,ER,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
            
            RHSv_3D(:,jR,irk) = RHSv_3D(:,jR,irk) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
            RHSv_3D(:,jL,irk) = RHSv_3D(:,jL,irk) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)    

        !....DQW
        
!        ZWL  = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
!        UL   = MATMUL(b_PHI_2D3Dy(:,:,2),VSE(:,hL,irk))
!        ZWR  = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
!        UR   = MATMUL(b_PHI_2D3Dy(:,:,1),VSE(:,hR,irk))
!        HUL  = MATMUL(bPHI_3D(:,:,6),Hv(:,jL,irk))
!        HUR  = MATMUL(bPHI_3D(:,:,5),Hv(:,jR,irk))
!        DO j = 1, NPTS_2D
!            FR(j)   = -ZWR(j)*UR(j)+HUR(j)
!            FL(j)   = -ZWL(j)*UL(j)+HUL(j)
!            Fhat(j) = half*(FR(j)+FL(j)) 
!        ENDDO 
!        DQW(:,jR) = DQW(:,jR) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
!        DQW(:,jL) = DQW(:,jL) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
    ENDIF
    
ENDDO 

CONTAINS

    FUNCTION Godunov(ZL,ZR,EL,ER,UL,UR,HuL,HuR)
    INTEGER(int_p) :: j
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, Godunov
    REAL(real_p), DIMENSION(NPTS_2D) :: HuL, HuR, EL, ER, ZL, ZR
    REAL(real_p) :: s, DEN
        DO j = 1, NPTS_2D
            FL(j) = HuL(j)*UL(j) + (1.0d0/2.0d0)*g*(EL(j)**2-ZL(j)**2) ! for dam break h**2 = 0 b/c flat bottom
            FR(j) = HuR(j)*UR(j) + (1.0d0/2.0d0)*g*(ER(j)**2-ZR(j)**2) !
            DEN = (HuR(j) - HuL(j))
            IF (abs(DEN) < 1d-16) THEN
                s = 0.0d0
            ELSE    
                s = (FR(j) - FL(j))/DEN
            ENDIF    
            IF (s > 0.0d0) THEN
                Godunov(j) = FL(j)
                IF (UL(j) < 0 .and. UR(j) > 0) THEN  ! Transonic case
                    Godunov(j) = (1.0d0/2.0d0)*g*(EL(j)**2-ZL(j)**2)
                ENDIF 
            ELSE    
                Godunov(j) = FR(j)
                IF (UL(j) < 0 .and. UR(j) > 0) THEN  ! Transonic case
                    Godunov(j) = (1.0d0/2.0d0)*g*(ER(j)**2-ZR(j)**2)
                ENDIF 
            ENDIF       
    
        ENDDO
    END FUNCTION Godunov

    FUNCTION LLFy(HL,HR,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, LLFy
    REAL(real_p), DIMENSION(NPTS_2D) :: HuL, HuR, EIGmax
    REAL(real_p), DIMENSION(NPTS_2D) :: HL, HR, ZL, ZR
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_2D
            FL(j)   = HuL(j)*UL(j) + half*g*(HL(j)**2.0d0 - ZL(j)**2.0d0) 
            FR(j)   = HuR(j)*UR(j) + half*g*(HR(j)**2.0d0 - ZR(j)**2.0d0) 
            LLFy(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFy
    
    FUNCTION LLFx(UL,UR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, LLFx
    REAL(real_p), DIMENSION(NPTS_2D) :: HuL, HuR, EIGmax
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_2D
            FL(j)   = HuL(j)*UL(j) 
            FR(j)   = HuR(j)*UR(j)  
            LLFx(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFx    
    
    FUNCTION LLFz(UL,UR,WL,WR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(NPTS_2D) :: FL, FR, UL, UR, LLFz
    REAL(real_p), DIMENSION(NPTS_2D) :: WL, WR, HuL, HuR, EIGmax
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_2D
            FL(j)   = -WL(j)*UL(j) 
            FR(j)   = -WR(j)*UR(j)  
            LLFz(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFz 

    FUNCTION EIG(HL,HR,UL,UR)
    REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, HL, HR, EIG
    REAL(real_p) :: LEp, LEm, REp, REm
        DO j = 1, NPTS_2D
            LEp    = UL(j) + sqrt(g*HL(j))
            LEm    = UL(j) - sqrt(g*HL(j))
            REp    = UR(j) + sqrt(g*HR(j))
            REm    = UR(j) - sqrt(g*HR(j))
            EIG(j) = max(abs(LEp),abs(LEm))
            EIG(j) = max(abs(EIG(j)),abs(REp))
            EIG(j) = max(abs(EIG(j)),abs(REm))
            EIG(j) = max(abs(EIG(j)),abs(UL(j)))
            EIG(j) = max(abs(EIG(j)),abs(UR(j)))
        ENDDO    
    END FUNCTION EIG  
    
    FUNCTION EIGz(WR,WL,EL,ER)
    REAL(real_p), DIMENSION(NPTS_2D) :: WR, WL, EL, ER, EIGz
    REAL(real_p) :: LE, RE
        DO j = 1, NPTS_2D
            LE = WL(j)/EL(j)
            RE = WR(j)/ER(j)
            EIGz(j) = max(abs(LE),abs(RE))
        ENDDO    
    END FUNCTION EIGz

    FUNCTION gfB(Uj,Nj,Hj)
	REAL(real_p), DIMENSION(NPTS_2D) :: Uj, Nj, Hj, gfB
	    DO j = 1, NPTS_2D
	        gfB = sqrt(Nj(j)/Hj(j))*Uj(j)
	    ENDDO    
    END FUNCTION gfB
    
    FUNCTION gfE(Uj,Nj,Hj)
    REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Nj, Hj, gfE
        DO j = 1, NPTS_3D
            gfE(j) = sqrt(Nj(j)/Hj(j))*Uj(j)
        ENDDO    
    END FUNCTION gfE
    
    FUNCTION sigE(Vj,Wj) 
    REAL(real_p), DIMENSION(NPTS_3D) :: Vj, Wj, sigE
        DO j = 1, NPTS_3D
            sigE(j) = -Vj(j)*Wj(j)
        ENDDO    
    END FUNCTION sigE
    
    FUNCTION sigB(Qj,Nj,Hj) 
    REAL(real_p), DIMENSION(NPTS_2D) :: Qj, Nj, Hj, sigB
        DO j = 1, NPTS_2D
            sigB(j) = sqrt(Nj(j)/Hj(j))*Qj(j) 
        ENDDO    
    END FUNCTION sigB    
    
    FUNCTION ymE(Hvj,Vj,Hj,Zj,g)
    REAL(real_p) :: g
    REAL(real_p), DIMENSION(NPTS_3D) :: Hvj, Vj, Hj, Zj, ymE
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, NPTS_3D
            ymE(j) = Hvj(j)*Vj(j) + half*g*(Hj(j)**2.0d0 - Zj(j)**2.0d0)
        ENDDO
    END FUNCTION ymE   
    
    FUNCTION xymE(Hvj,Vj)
    REAL(real_p), DIMENSION(NPTS_3D) :: Hvj, Vj, xymE     
        DO j = 1, NPTS_3D
            xymE(j) = Hvj(j)*Vj(j)
        ENDDO
    END FUNCTION xymE   
    
    FUNCTION ymSE(DZj,Ej,g)
    REAL(real_p) :: g
    REAL(real_p), DIMENSION(L2_3D) :: ymSE, DZj, Ej 
        DO j = 1, L2_3D
            ymSE(j) = g*Ej(j)*DZj(j)
        ENDDO
    END FUNCTION ymSE    
            
END SUBROUTINE RHS_HV3D
!------------------------------------------------------------------
!
!           RHS for Depth Integrated Velocities
!               Written by Colton J. Conroy
!                      @ The C.H.I.L
!                          1.17.13
!
!-----------------------------------------------------------------
SUBROUTINE RHS_Ubar(TIME,ELEM)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: jL, jR, ELEM, H_Elem
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt
REAL(real_p) :: UL, UR, EL, ER, Ghat, GL, GR, hpt, nn, k
REAL(real_p), DIMENSION(L2int) :: L2Vavg
REAL(real_p), DIMENSION(p+1) :: Uaj, Gaj
REAL(real_p), DIMENSION(L2int) :: Saj 
REAL(real_p), PARAMETER :: half = 0.500000000000000000d0

!....Get horizontal velocity @ bottom boundary

IF (ELEM == 1) THEN

    DO i = 1, nelemsX

        DO j = 1, LE
        
            Uh(j,i) = DOT_PRODUCT(bPHIz(1,:),U(:,1,i,j,irk))
        
        ENDDO
    
    ENDDO
    
ELSEIF (ELEM == 2) THEN

    DO i = 1, Qelems
    
        IF (zEDGE(i,1) == -1) THEN
            
            H_Elem  = CONNzx(i,1)
            Uh(:,H_Elem) = MATMUL(bPHI_2D(:,:,1),U_2D(:,i,irk))
            
        ENDIF
        
    ENDDO
    
ELSE 

    zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, nelemsX

        INVJAC = 1/xELEM_jac(i)

        DO j = 1, L2int

            xpt = X(i) + INVJAC*(1+L2pts(j,1))                          !....Transform pts
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real
            hpt = PROB_DATA%h0 + PROB_DATA%Hslope*xpt**2                                                        !....variables to complex
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
        
            Uh(j,i) = Vu
        
        ENDDO
    
    ENDDO

ENDIF

!....Loop over Elements

DO i = 1, nelemsX

    Uaj              = MATMUL(PHI,ZETA(:,i,irk))
    Saj              = -PROB_DATA%k/PROB_DATA%h0*Uh(:,i)
    Gaj              = gE(Uaj)
    RHSubar(:,i,irk) = xELEM_jac(i)*MATMUL(A,Gaj) + MATMUL(C,Saj)
    
ENDDO

!....Loop over Boundaries

UL               = 0.000000000000d0
UR               = DOT_PRODUCT(bPHI(1,:),USE(:,1,irk))
EL               = DOT_PRODUCT(bPHI(1,:),ZETA(:,1,irk))
GL               = gB(EL)
Ghat             = half*(GL + GL - sqrt(PROB_DATA%g*PROB_DATA%h0)*(UR-UL))
RHSubar(:,1,irk) = RHSubar(:,1,irk) + xELEM_jac(1)*B1*Ghat

DO i = 2, nnodesX-1

    jL                = xNODE_elems(i,1)
    jR                = xNODE_elems(i,2)
    UL                = DOT_PRODUCT(bPHI(2,:),USE(:,jL,irk))
    UR                = DOT_PRODUCT(bPHI(1,:),USE(:,jR,irk))
    EL                = DOT_PRODUCT(bPHI(2,:),ZETA(:,jL,irk))
    ER                = DOT_PRODUCT(bPHI(1,:),ZETA(:,jR,irk))
    GL                = gB(EL)
    GR                = gB(ER)
    Ghat              = half*(GL + GR - sqrt(PROB_DATA%g*PROB_DATA%h0)*(UR-UL))
    RHSubar(:,jL,irk) = RHSubar(:,jL,irk) + xELEM_jac(jL)*B2*Ghat
    RHSubar(:,jR,irk) = RHSubar(:,jR,irk) + xELEM_jac(jR)*B1*Ghat
    
ENDDO

UL = DOT_PRODUCT(bPHI(2,:),USE(:,nelemsX,irk))
EL = DOT_PRODUCT(bPHI(2,:),ZETA(:,nelemsX,irk))
GL = gB(EL)

zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)
xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
hpt = PROB_DATA%h0 + PROB_DATA%Hslope*PROB_DATA%xN**2
CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)

UR   = Vavg
ER   = SE
GR   = gB(SE)
Ghat = half*(GL + GR - sqrt(PROB_DATA%g*PROB_DATA%h0)*(UR-UL))

RHSubar(:,nelemsX,irk) = RHSubar(:,nelemsX,irk) + xELEM_jac(nelemsX)*B2*Ghat

CONTAINS

    FUNCTION gE(Y)
    REAL(real_p), DIMENSION(p+1) :: gE, Y
     gE = PROB_DATA%g*Y
    END FUNCTION gE
    
    FUNCTION gB(Y)
    REAL(real_p) :: gB, Y
     gB = PROB_DATA%g*Y
    END FUNCTION gB

END SUBROUTINE RHS_Ubar
!-----------------------------------------------------------------
!
!           RHS for 2D Depth Integrated Velocity
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                            8.15.13
!
!-----------------------------------------------------------------
SUBROUTINE RHS_Ubar2D(TIME,ELEM)  !....Need to UPdate to take into account d/dy(zeta)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, yCOORD, zCOORD
INTEGER(int_p) :: jL, jR, ELEM, H_Elem, nedgesX, nfacesZ
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt, ypt, hpt, nn, k, EIG
REAL(real_p), DIMENSION(p+1) :: UL, UR, EL, ER, Ghat, GL, GR, HR, HL
REAL(real_p), DIMENSION(NPTS_2D) :: Uaj, Gaj
REAL(real_p), DIMENSION(L2_2D) :: Saj, Kh, Haj
REAL(real_p), DIMENSION(2) :: Xnode, Ynode
REAL(real_p), DIMENSION(4) :: globalnodes
REAL(real_p), PARAMETER :: half = 0.500000000000000000d0

nedgesX = size(xEDGE,1)

UL   = 0.0d0
UR   = 0.0d0
EL   = 0.0d0
ER   = 0.0d0
Ghat = 0.0d0
GL   = 0.0d0
GR   = 0.0d0
HR   = 0.0d0
HL   = 0.0d0
Uaj  = 0.0d0
Gaj  = 0.0d0
Saj  = 0.0d0
Kh   = 0.0d0
Haj  = 0.0d0
Xnode = 0.0d0
Ynode = 0.0d0
globalnodes = 0.0d0

!....Get horizontal velocity @ bottom boundary

IF (ELEM == 1) THEN

    DO i = 1, nelemsX

        DO j = 1, LE
        
            Uh(j,i) = DOT_PRODUCT(bPHIz(1,:),U(:,1,i,j,irk))
        
        ENDDO
    
    ENDDO
    
ELSEIF (ELEM == 2) THEN

    DO i = 1, Qelems
    
        IF (zEDGE(i,1) == -1) THEN ! fix this
            
            H_Elem  = CONNzx(i,1)
            Uh(:,H_Elem) = MATMUL(bPHI_2D(:,:,1),U_2D(:,i,irk))
            
        ENDIF
        
    ENDDO
    
ELSEIF (ELEM == 3) THEN

    nfacesZ = size(zFACE,1)
    
    DO i = 1, nfacesZ
        IF (zFACE(i,1) == -1) THEN
            jR = zFACE(i,2)
            H_Elem = HEXzxy(jR,1)
            Uh(:,H_Elem) = MATMUL(bPHI_3D2D,U_3D(:,jR,irk))
        ENDIF                    
    ENDDO
    
ENDIF   
 

IF (coupled == 1) THEN

    !zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
        Ynode(1) = Qnode(globalNODES(1),2)
        Ynode(2) = Qnode(globalNODES(2),2)
       
        DO j = 1, L2_2D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
            xCOORD = dcmplx(xpt,0.0d0)
            yCOORD = dcmplx(ypt,0.0d0)
            !hpt = PROB_DATA%h0 + PROB_DATA%Hslope*xpt**2
            hpt = PROB_DATA%Hslope*xpt**2
            zCOORD = dcmplx(hpt*Z(1), 0.0d0)
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
            Uh(j,i) = Vu
                    
                   
        ENDDO
        

    ENDDO
    
ENDIF    

!....Loop over elements

DO i = 1, Qelems

    Uaj = MATMUL(PHIz,ZETA(:,i,irk))
    !Haj = MATMUL(L2PHIz,Hh(:,i))
    Haj =  MATMUL(PHIL2_h,Hh(:,i))
    Kh  = MATMUL(L2PHIz,KX(:,i))
    DO j = 1, L2_2D
        Saj(j) = -Kh(j)/Haj(j)*Uh(j,i)
    ENDDO
    !Saj = -PROB_DATA%k/PROB_DATA%h0*Uh(:,i)
    Gaj = gE(Uaj)
    RHSubar(:,i,irk) = CONN_jacX(i,1)*MATMUL(A_2D(:,:,1),Gaj) + MATMUL(CU,Saj)
    
ENDDO

!....Loop over boundaries

DO i = 1, nedgesX

    jL = xEDGE(i,1)
    jR = xEDGE(i,2)
    
    IF (jL == -1) THEN
        
        UL = 0.0d0
        UR = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
        ER = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
        GL = gB(ER)
        HR = MATMUL(bPHI_2Dxh(:,:,1),Hh(:,jR))
        Ghat = half*(GL + GL - SQRT(PROB_DATA%g*HR)*(UR-UL))
        RHSubar(:,jR,irk) = RHSubar(:,jR,irk) + MATMUL(BX(:,:,1),Ghat)*CONN_jacX(jR,1)
        
    ELSEIF (jR == 0) THEN
    
        UL = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
        EL = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
        GL = gB(EL)

        !zCOORD = dcmplx(PROB_DATA%h0*Z(1), 0.0d0)
        xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
        !hpt = PROB_DATA%h0 + PROB_DATA%Hslope*PROB_DATA%xN**2
        hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
        zCOORD = dcmplx(hpt*Z(1), 0.0d0)
        CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)

        UR   = Vavg
        ER   = SE
        GR   = gB(ER)
        Ghat = half*(GL + GR - sqrt(PROB_DATA%g*hpt)*(UR-UL))
        RHSubar(:,jL,irk) = RHSubar(:,jL,irk) - MATMUL(BX(:,:,2),Ghat)*CONN_jacX(jL,1)
        
    ELSE
    
        UL = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
        EL = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
        GL = gB(EL)     
        UR = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
        ER = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
        GR = gB(ER) 
        HR    = MATMUL(bPHI_2Dxh(:,:,1),Hh(:,jR))
        HL    = MATMUL(bPHI_2Dxh(:,:,2),Hh(:,jL))     
        EIG   = max(abs(HL(1)),abs(HR(1)))     
        Ghat = half*(GL + GR - sqrt(PROB_DATA%g*EIG)*(UR-UL))
        RHSubar(:,jR,irk) = RHSubar(:,jR,irk) + MATMUL(BX(:,:,1),Ghat)*CONN_jacX(jR,1)
        RHSubar(:,jL,irk) = RHSubar(:,jL,irk) - MATMUL(BX(:,:,2),Ghat)*CONN_jacX(jL,1)
        
    ENDIF
   
ENDDO    



CONTAINS

    FUNCTION gE(Y)
    REAL(real_p), DIMENSION(NPTS_2D) :: gE, Y
     gE = PROB_DATA%g*Y
    END FUNCTION gE
    
    FUNCTION gB(Y)
    REAL(real_p), DIMENSION(p+1) :: gB, Y
     gB = PROB_DATA%g*Y
    END FUNCTION gB

END SUBROUTINE RHS_Ubar2D    
!-----------------------------------------------------------------
!
!          RHS for Non-linear 2D Depth Integrated Velocity
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                            3.18.14
!
!-----------------------------------------------------------------
SUBROUTINE RHS_Ubar2D_NL(TIME,ELEM)  

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD, yCOORD
INTEGER(int_p) :: jL, jR, nedgesX, nedgesY, ELEM, H_Elem, nfacesZ
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt, ypt, zpt, nn, k, hpt, g, freq, Lx
REAL(real_p) :: solnT, Vv, Vw, Huse, Zb, DSEx, DSEy, Uavg, Fpce, Fxm, Fym, Fw,FxmA, FymA
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Ynode
REAL(real_p), DIMENSION(p+1) :: UL, UR, FL, FR, Fhat, HL, HR, EIGMAX
REAL(real_p), DIMENSION(p+1) :: HUL, HUR, ZL, ZR
REAL(real_p), DIMENSION(L2_2D) :: L2Vavg, L2Hse, L2Hz, L2Kx, Saj, Khj, DZj, SEj, Sxj
REAL(real_p), DIMENSION(NPTS_2D) :: Uaj, Fxj, Fyj, Haj, Vaj, Hhj, Huj
REAL(real_p), DIMENSION(L2_2D,Qelems) :: L2FxmA, L2wind
REAL(real_p), PARAMETER :: half = 0.500000000000000000d0

nedgesX = size(xEDGE,1)
nedgesY = size(yEDGE,1)
solnT = REAL(TIME)
g     = PROB_DATA%g
freq  = PROB_DATA%freq
Lx    = PROB_DATA%xN - PROB_DATA%x0

UL   = 0.0d0
UR   = 0.0d0
FL   = 0.0d0
FR   = 0.0d0
Fhat = 0.0d0
ZL   = 0.0d0
ZR   = 0.0d0
HR   = 0.0d0
HL   = 0.0d0
HUL  = 0.0d0
HUR  = 0.0d0
Uaj  = 0.0d0
Fxj  = 0.0d0
Fyj  = 0.0d0
Vaj  = 0.0d0
Hhj  = 0.0d0
Huj  = 0.0d0
Saj  = 0.0d0
DZj  = 0.0d0
L2Kx = 0.0d0
L2Hz = 0.0d0
L2Vavg = 0.0d0
L2HSe = 0.0d0
SEj   = 0.0d0
Sxj   = 0.0d0
Khj   = 0.0d0
Haj   = 0.0d0
L2FxmA = 0.0d0
L2wind = 0.0d0
EIGMAX = 0.0d0
Xnode = 0.0d0
Ynode = 0.0d0
globalnodes = 0.0d0

IF (PROB == 3) THEN  

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
        Ynode(1) = Qnode(globalNODES(1),2)
        Ynode(2) = Qnode(globalNODES(2),2)
        
        DO j = 1, L2_2D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
            Hpt = PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt,0.0d0)
            yCOORD = dcmplx(ypt,0.0d0)
            zCOORD = dcmplx(Hpt*Z(1),0.0d0)
            IF (PROB == 1) THEN
            
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
                L2Vavg(j)    = Vavg
                
            ELSEIF (PROB == 3) THEN
            
                zpt = Z(1)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2Hse(j)    = Huse
                L2Hz(j)     = Zb
                L2Vavg(j)   = Vavg
                L2FxmA(j,i) = FxmA 
                Uh(j,i)     = Vu  
                L2Kx(j)     = PROB_DATA%k     
                L2wind(j,i) = PROB_DATA%k*exp(sin((2.0d0*PI*xpt/Lx) + solnT*freq))  
                  
            ENDIF                     
        ENDDO
        IF (coupled == 1) THEN
            ZETA(:,i,irk)  = MATMUL(L2U,L2Hse)               
            VSE(:,i,irk)   = MATMUL(L2U,L2Vavg)
            Hh(:,i)        = MATMUL(L2Bh,L2Hz)
            KX(:,i)        = MATMUL(L2U,L2Kx)
        ENDIF
    ENDDO
    IF (coupled == 1) THEN
        CALL CALC_SE(ELEM)
        CALL CALC_Ubar(ELEM)
    ENDIF
ENDIF  

IF (coupled == 3) THEN

    nfacesZ = size(zFACE,1)
    
    DO i = 1, nfacesZ
        IF (zFACE(i,1) == -1) THEN
            jR = zFACE(i,2)
            H_Elem = HEXzxy(jR,1)
            Uh(:,H_Elem) = MATMUL(bPHI_3D2D,U_3D(:,jR,irk))
        ENDIF                    
    ENDDO
    
ENDIF

IF (PROB == 4) THEN
    L2wind = 0.0d0
    L2FxmA = 0.0d0
    Uh     = 0.0d0
ENDIF    
  
!....Loop over elements for area integrals

DO i = 1, Qelems  

    Uaj = MATMUL(PHIz,USE(:,i,irk))
    Vaj = MATMUL(PHIz,VSE(:,i,irk))
    Haj = MATMUL(PHIz,ZETA(:,i,irk))
    Huj = MATMUL(PHIz,Hu_bar(:,i,irk))
    Hhj = MATMUL(PHI_h,Hh(:,i))

    DO j = 1, NPTS_2D
        Fxj(j)  = Huj(j)*Uaj(j) + half*g*(Haj(j)**2.0d0 - Hhj(j)**2.0d0)
        Fyj(j)  = Huj(j)*Vaj(j)
    ENDDO
    
    DZj = MATMUL(DPHI_h(:,:,1),Hh(:,i))*CONN_jacX(i,1)
    SEj = MATMUL(L2PHIz,ZETA_2D(:,i,irk))
    Khj = MATMUL(L2PHIz,KX(:,i))    
    
    DO j = 1, L2_2D
        Sxj(j) = g*SEj(j)*DZj(j)
        Saj(j) = -Khj(j)*Uh(j,i)
    ENDDO
    
    RHSubar(:,i,irk) = MATMUL(A_2D(:,:,1),Fxj)*CONN_jacX(i,1) + MATMUL(A_2D(:,:,2),Fyj)*CONN_jacY(i,1) &
                     + MATMUL(CU,L2FxmA(:,i)) + MATMUL(CU,Sxj) + MATMUL(CU,Saj) + MATMUL(CU,L2wind(:,i))
ENDDO

!....Loop over edges for flux calculations  

DO i = 1, nedgesX

    jL = xEDGE(i,1)
    jR = xEDGE(i,2)
    
    IF (jL == -1) THEN
        IF (PROB == 1) THEN

        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            HuR = MATMUL(bPHI_2Dx(:,:,1),Hu_bar(:,jR,irk))
            UR  = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
            HR  = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
            ZR  = MATMUL(bPHI_2Dxh(:,:,1),Hh(:,jR))
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,1,1))*Xnode(1)+(1.0d0+bPTSx(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,2,1))*Ynode(1)+(1.0d0+bPTSx(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HL(j)   = Huse
                UL(j)   = Uavg
                HuL(j)  = Huse*Uavg
                ZL(j)   = Zb
                Fhat(j) = Huse*Uavg*Uavg + half*g*(Huse**2.0d0 - Zb**2.0d0)
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            !Fhat = LLFx(HL,HR,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
        ELSEIF (PROB == 4) THEN
            UL = 0.0d0
            HL = 10.0d0
            DO j = 1, p+1
                Fhat(j) = half*g*HL(j)**2.0d0
            ENDDO    
        ENDIF        
            
        RHSubar(:,jR,irk) = RHSubar(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)
        
    ELSEIF (jR == 0) THEN
       
        IF (PROB == 1) THEN

        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jL,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            HuL = MATMUL(bPHI_2Dx(:,:,2),Hu_bar(:,jL,irk))
            UL  = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
            HL  = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
            ZL  = MATMUL(bPHI_2Dxh(:,:,2),Hh(:,jL))
            DO j = 1, p+1
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,1,2))*Xnode(1)+(1.0d0+bPTSx(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,2,2))*Ynode(1)+(1.0d0+bPTSx(j,2,2))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HR(j)   = Huse
                UR(j)   = Uavg
                HuR(j)  = Huse*Uavg
                ZR(j)   = Zb
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat = LLFx(HL,HR,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
        ELSEIF (PROB == 4) THEN
            UR = 0.0d0
            HR = 5.0d0
            DO j = 1, p+1
                Fhat(j) = half*g*HR(j)**2.0d0
            ENDDO   
        ENDIF             
        
        RHSubar(:,jL,irk) = RHSubar(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
    
    ELSE
    
        HL  = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
        UL  = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
        HuL = MATMUL(bPHI_2Dx(:,:,2),Hu_bar(:,jL,irk))
        ZL  = MATMUL(bPHI_2Dxh(:,:,2),Hh(:,jL))
        HR  = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
        UR  = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
        HuR = MATMUL(bPHI_2Dx(:,:,1),Hu_bar(:,jR,irk))
        ZR  = MATMUL(bPHI_2Dxh(:,:,1),Hh(:,jR))
        EIGmax = EIG(HL,HR,UL,UR)
        Fhat   = LLFx(HL,HR,ZL,ZR,UL,UR,HuL,HuR,EIGmax)  
        RHSubar(:,jL,irk) = RHSubar(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
        RHSubar(:,jR,irk) = RHSubar(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)

    ENDIF
    
ENDDO  

!....Loop over y-edges    

DO i = 1, nedgesY

    jL = yEDGE(i,1)
    jR = yEDGE(i,2)
    
    IF (jL == -1) THEN
        IF (PROB == 1) THEN
            Fhat = 0.0d0
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UR  = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
            HuR = MATMUL(bPHI_2Dy(:,:,1),Hu_bar(:,jR,irk))
            HR  = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,1,1))*Xnode(1)+(1.0d0+bPTSy(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,2,1))*Ynode(1)+(1.0d0+bPTSy(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                UL(j)   = Vavg
                HuL(j)  = Huse*Uavg
                HL(j)   = Huse
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat   = LLFy(UL,UR,HuL,HuR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF        
            
        RHSubar(:,jR,irk) = RHSubar(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)
        
    ELSEIF (jR == 0) THEN
       
        IF (PROB == 1) THEN
            Fhat = 0.0d0
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jL,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UL  = MATMUL(bPHI_2Dy(:,:,2),VSE(:,jL,irk))
            HL  = MATMUL(bPHI_2Dy(:,:,2),ZETA(:,jL,irk))
            HuL = MATMUL(bPHI_2Dy(:,:,2),Hu_bar(:,jL,irk))
            DO j = 1, p+1
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,1,2))*Xnode(1)+(1.0d0+bPTSy(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,2,2))*Ynode(1)+(1.0d0+bPTSy(j,2,2))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HR(j)   = Huse
                UR(j)   = Vavg
                HuR(j)  = Huse*Uavg
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat   = LLFy(UL,UR,HuL,HuR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF             
        
        RHSubar(:,jL,irk) = RHSubar(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
    
    ELSE
    
        HL  = MATMUL(bPHI_2Dy(:,:,2),ZETA(:,jL,irk))
        UL  = MATMUL(bPHI_2Dy(:,:,2),VSE(:,jL,irk))
        HuL = MATMUL(bPHI_2Dy(:,:,2),Hu_bar(:,jL,irk))
        HR  = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))
        UR  = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
        HuR = MATMUL(bPHI_2Dy(:,:,1),Hu_bar(:,jR,irk))
        EIGmax = EIG(HL,HR,UL,UR)
        Fhat   = LLFy(UL,UR,HuL,HuR,EIGmax)   
        RHSubar(:,jL,irk) = RHSubar(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
        RHSubar(:,jR,irk) = RHSubar(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)
        
    ENDIF
    
ENDDO  

CONTAINS

    FUNCTION LLFx(HL,HR,ZL,ZR,UL,UR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(p+1) :: FL, FR, UL, UR, LLFx
    REAL(real_p), DIMENSION(p+1) :: HuL, HuR, EIGmax
    REAL(real_p), DIMENSION(p+1) :: HL, HR, ZL, ZR
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, p+1
            FL(j)   = HuL(j)*UL(j) + half*g*(HL(j)**2.0d0 - ZL(j)**2.0d0) 
            FR(j)   = HuR(j)*UR(j) + half*g*(HR(j)**2.0d0 - ZR(j)**2.0d0) 
            LLFx(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFx
    
    FUNCTION LLFy(UL,UR,HuL,HuR,EIGmax)
    REAL(real_p), DIMENSION(p+1) :: FL, FR, UL, UR, LLFy
    REAL(real_p), DIMENSION(p+1) :: HuL, HuR, EIGmax
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, p+1
            FL(j)   = HuL(j)*UL(j) 
            FR(j)   = HuR(j)*UR(j)  
            LLFy(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HuR(j)-HuL(j))) 
        ENDDO
    END FUNCTION LLFy  
    
    FUNCTION EIG(HL,HR,UL,UR)
    REAL(real_p), DIMENSION(p+1) :: UL,UR, HL, HR
    REAL(real_p), DIMENSION(p+1) :: LEp, LEm, REp, REm, EIG
    DO j = 1, p+1
        LEp(j) = UL(j) + sqrt(g*HL(j))
        LEm(j) = UL(j) - sqrt(g*HL(j))
        REp(j) = UR(j) + sqrt(g*HR(j))
        REm(j) = UR(j) - sqrt(g*HR(j))
        EIG(j) = max(abs(LEp(j)),abs(LEm(j)))
        EIG(j) = max(abs(EIG(j)),abs(REp(j)))
        EIG(j) = max(abs(EIG(j)),abs(REm(j)))
        EIG(j) = max(abs(EIG(j)),abs(UL(j)))
        EIG(j) = max(abs(EIG(j)),abs(UR(j)))
    ENDDO    
    END FUNCTION EIG     

END SUBROUTINE RHS_Ubar2D_NL 
!-----------------------------------------------------------------
!
!          RHS for Non-linear 2D Depth Integrated Velocity
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                            3.18.14
!
!-----------------------------------------------------------------
SUBROUTINE RHS_Vbar2D_NL(TIME,ELEM)  

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD, yCOORD
INTEGER(int_p) :: jL, jR, nedgesX, nedgesY, ELEM, H_Elem, nfacesZ
REAL(real_p) :: SE, grad_SE, Vu, Vavg, INVJAC, xpt, ypt, zpt, nn, k, hpt, g, freq, Ly
REAL(real_p) :: solnT, Vv, Vw, Huse, Zb, DSEx, DSEy, Uavg, Fpce, Fxm, Fym, Fw, FxmA, FymA
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Ynode
REAL(real_p), DIMENSION(p+1) :: UL, UR, FL, FR, Fhat, HL, HR, EIGMAX
REAL(real_p), DIMENSION(p+1) :: HvL, HvR, ZL, ZR
REAL(real_p), DIMENSION(L2_2D) :: L2Uavg, L2Hse, L2Hz, L2Kx, Saj, Khj, DZj, SEj, Syj
REAL(real_p), DIMENSION(NPTS_2D) :: Uaj, Fxj, Fyj, Haj, Vaj, Hhj, Hvj
REAL(real_p), DIMENSION(L2_2D,Qelems) :: L2FymA, L2wind
REAL(real_p), PARAMETER :: half = 0.500000000000000000d0

nedgesX = size(xEDGE,1)
nedgesY = size(yEDGE,1)
solnT = REAL(TIME)
g     = PROB_DATA%g
freq  = PROB_DATA%freq
Ly    = PROB_DATA%yN - PROB_DATA%y0

UL   = 0.0d0
UR   = 0.0d0
FL   = 0.0d0
FR   = 0.0d0
Fhat = 0.0d0
ZL   = 0.0d0
ZR   = 0.0d0
HR   = 0.0d0
HL   = 0.0d0
HvL  = 0.0d0
HvR  = 0.0d0
Uaj  = 0.0d0
Fxj  = 0.0d0
Fyj  = 0.0d0
Vaj  = 0.0d0
Hhj  = 0.0d0
Hvj  = 0.0d0
Saj  = 0.0d0
DZj  = 0.0d0
L2Kx = 0.0d0
L2Hz = 0.0d0
EIGMAX = 0.0d0
L2Uavg = 0.0d0
L2HSe  = 0.0d0
SEj    = 0.0d0
Syj    = 0.0d0
Khj    = 0.0d0
Haj    = 0.0d0
L2FymA = 0.0d0
L2wind = 0.0d0
Xnode  = 0.0d0
Ynode  = 0.0d0
globalnodes = 0.0d0

IF (PROB == 3) THEN 

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
        Ynode(1) = Qnode(globalNODES(1),2)
        Ynode(2) = Qnode(globalNODES(2),2)
        
        DO j = 1, L2_2D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
            Hpt = PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt,0.0d0)
            yCOORD = dcmplx(ypt,0.0d0)
            zCOORD = dcmplx(Hpt*Z(1),0.0d0)
            IF (PROB == 1) THEN
            
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
                L2Uavg(j)    = Vavg
                
            ELSEIF (PROB == 3) THEN
            
                zpt = Z(1)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2Hse(j)    = Huse
                L2Hz(j)     = Zb
                L2Uavg(j)   = Uavg
                L2FymA(j,i) = FymA 
                Vh(j,i)     = Vv  
                L2Kx(j)     = PROB_DATA%k     
                L2wind(j,i) = PROB_DATA%k*exp(sin((2.0d0*PI*ypt/Ly) + solnT*freq))  
                  
            ENDIF                     
        ENDDO
        IF (coupled == 1) THEN
            ZETA(:,i,irk)  = MATMUL(L2U,L2Hse)               
            USE(:,i,irk)   = MATMUL(L2U,L2Uavg)
            Hh(:,i)        = MATMUL(L2Bh,L2Hz)
            KX(:,i)        = MATMUL(L2U,L2Kx)
        ENDIF
    ENDDO
    IF (PROB == 1) THEN
        CALL CALC_SE(ELEM)
        CALL CALC_Vbar(ELEM)
    ENDIF
ENDIF  

IF (coupled == 3) THEN

    nfacesZ = size(zFACE,1)
    
    DO i = 1, nfacesZ
        IF (zFACE(i,1) == -1) THEN
            jR = zFACE(i,2)
            H_Elem = HEXzxy(jR,1)
            Vh(:,H_Elem) = MATMUL(bPHI_3D2D,V_3D(:,jR,irk))
        ENDIF                    
    ENDDO
    
ENDIF

IF (PROB == 4) THEN
    L2wind = 0.0d0
    L2FymA = 0.0d0
    Vh     = 0.0d0
ENDIF
  
!....Loop over elements for area integrals

DO i = 1, Qelems  

    !....Element Integrals

    Uaj = MATMUL(PHIz,USE(:,i,irk))
    Vaj = MATMUL(PHIz,VSE(:,i,irk))
    Haj = MATMUL(PHIz,ZETA(:,i,irk))
    Hvj = MATMUL(PHIz,Hv_bar(:,i,irk))
    Hhj = MATMUL(PHI_h,Hh(:,i))

    DO j = 1, NPTS_2D
        Fxj(j)  = Hvj(j)*Uaj(j) 
        Fyj(j)  = Hvj(j)*Vaj(j) + half*g*(Haj(j)**2.0d0 - Hhj(j)**2.0d0)
    ENDDO
    
    !....Source Integrals
    
    DZj = MATMUL(DPHI_h(:,:,2),Hh(:,i))*CONN_jacY(i,1)
    SEj = MATMUL(L2PHIz,ZETA_2D(:,i,irk))
    Khj = MATMUL(L2PHIz,KX(:,i))    
    
    DO j = 1, L2_2D
        Syj(j) = g*SEj(j)*DZj(j)
        Saj(j) = -Khj(j)*Vh(j,i)
    ENDDO
    
    RHSvbar(:,i,irk) = MATMUL(A_2D(:,:,1),Fxj)*CONN_jacX(i,1) + MATMUL(A_2D(:,:,2),Fyj)*CONN_jacY(i,1) &
                     + MATMUL(CU,L2FymA(:,i)) + MATMUL(CU,Syj) + MATMUL(CU,Saj) + MATMUL(CU,L2wind(:,i))
ENDDO

!....Loop over edges for flux calculations  

DO i = 1, nedgesX

    jL = xEDGE(i,1)
    jR = xEDGE(i,2)
    
    IF (jL == -1) THEN
        IF (PROB == 1) THEN

        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            HvR = MATMUL(bPHI_2Dx(:,:,1),Hv_bar(:,jR,irk))
            UR  = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
            HR  = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,1,1))*Xnode(1)+(1.0d0+bPTSx(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,2,1))*Ynode(1)+(1.0d0+bPTSx(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HL(j)   = Huse
                UL(j)   = Uavg
                HvL(j)  = Huse*Vavg
                !Fhat(j) = Huse*Vavg*Uavg 
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat   = LLFx(UL,UR,HvL,HvR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF        
            
        RHSvbar(:,jR,irk) = RHSvbar(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)
        
    ELSEIF (jR == 0) THEN
       
        IF (PROB == 1) THEN

        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jL,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            HvL = MATMUL(bPHI_2Dx(:,:,2),Hv_bar(:,jL,irk))
            UL  = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
            HL  = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
            DO j = 1, p+1
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,1,2))*Xnode(1)+(1.0d0+bPTSx(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,2,2))*Ynode(1)+(1.0d0+bPTSx(j,2,2))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HR(j)   = Huse
                UR(j)   = Uavg
                HvR(j)  = Huse*Vavg
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat   = LLFx(UL,UR,HvL,HvR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF             
        
        RHSvbar(:,jL,irk) = RHSvbar(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
    
    ELSE
    
        HL  = MATMUL(bPHI_2Dx(:,:,2),ZETA(:,jL,irk))
        UL  = MATMUL(bPHI_2Dx(:,:,2),USE(:,jL,irk))
        HvL = MATMUL(bPHI_2Dx(:,:,2),Hv_bar(:,jL,irk))
        HR  = MATMUL(bPHI_2Dx(:,:,1),ZETA(:,jR,irk))
        UR  = MATMUL(bPHI_2Dx(:,:,1),USE(:,jR,irk))
        HvR = MATMUL(bPHI_2Dx(:,:,1),Hv_bar(:,jR,irk))
        EIGmax = EIG(HL,HR,UL,UR)
        Fhat   = LLFx(UL,UR,HvL,HvR,EIGmax)
        RHSvbar(:,jL,irk) = RHSvbar(:,jL,irk) - MATMUL(BX(:,:,2),Fhat)*CONN_jacX(jL,1)
        RHSvbar(:,jR,irk) = RHSvbar(:,jR,irk) + MATMUL(BX(:,:,1),Fhat)*CONN_jacX(jR,1)
        
    ENDIF
    
ENDDO  

!....Loop over y-edges    

DO i = 1, nedgesY

    jL = yEDGE(i,1)
    jR = yEDGE(i,2)
    
    IF (jL == -1) THEN
        IF (PROB == 1) THEN
            
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UR  = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
            HvR = MATMUL(bPHI_2Dy(:,:,1),Hv_bar(:,jR,irk))
            HR  = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))
            ZR  = MATMUL(bPHI_2Dyh(:,:,1),Hh(:,jR))
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,1,1))*Xnode(1)+(1.0d0+bPTSy(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,2,1))*Ynode(1)+(1.0d0+bPTSy(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                UL(j)   = Vavg
                HvL(j)  = Huse*Vavg
                HL(j)   = Huse
                ZL(j)   = Zb
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat   = LLFy(HL,HR,ZL,ZR,UL,UR,HvL,HvR,EIGmax)
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF        
            
        RHSvbar(:,jR,irk) = RHSvbar(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)
        
    ELSEIF (jR == 0) THEN
       
        IF (PROB == 1) THEN
            
        ELSEIF (PROB == 3) THEN
            globalNODES = CONN(jL,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            UL  = MATMUL(bPHI_2Dy(:,:,2),VSE(:,jL,irk))
            HL  = MATMUL(bPHI_2Dy(:,:,2),ZETA(:,jL,irk))
            ZL  = MATMUL(bPHI_2Dyh(:,:,2),Hh(:,jL))
            HvL = MATMUL(bPHI_2Dy(:,:,2),Hv_bar(:,jL,irk))
            DO j = 1, p+1
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,1,2))*Xnode(1)+(1.0d0+bPTSy(j,1,2))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,2,2))*Ynode(1)+(1.0d0+bPTSy(j,2,2))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                HR(j)   = Huse
                ZR(j)   = Zb
                UR(j)   = Vavg
                HvR(j)  = Huse*Vavg
            ENDDO
            EIGmax = EIG(HL,HR,UL,UR)
            Fhat   = LLFy(HL,HR,ZL,ZR,UL,UR,HvL,HvR,EIGmax) 
        ELSEIF (PROB == 4) THEN
            Fhat = 0.0d0    
        ENDIF             
        
        RHSvbar(:,jL,irk) = RHSvbar(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
    
    ELSE
    
        HL  = MATMUL(bPHI_2Dy(:,:,2),ZETA(:,jL,irk))
        UL  = MATMUL(bPHI_2Dy(:,:,2),VSE(:,jL,irk))
        HvL = MATMUL(bPHI_2Dy(:,:,2),Hv_bar(:,jL,irk))
        ZL  = MATMUL(bPHI_2Dy(:,:,2),Hh(:,jL))
        HR  = MATMUL(bPHI_2Dy(:,:,1),ZETA(:,jR,irk))
        ZR  = MATMUL(bPHI_2Dy(:,:,1),Hh(:,jR))
        UR  = MATMUL(bPHI_2Dy(:,:,1),VSE(:,jR,irk))
        HvR = MATMUL(bPHI_2Dy(:,:,1),Hv_bar(:,jR,irk))
        EIGmax = EIG(HL,HR,UL,UR)
        Fhat   = LLFy(HL,HR,ZL,ZR,UL,UR,HvL,HvR,EIGmax)   
        RHSvbar(:,jL,irk) = RHSvbar(:,jL,irk) - MATMUL(BU(:,:,2),Fhat)*CONN_jacY(jL,1)
        RHSvbar(:,jR,irk) = RHSvbar(:,jR,irk) + MATMUL(BU(:,:,1),Fhat)*CONN_jacY(jR,1)

    ENDIF
    
ENDDO  

CONTAINS

    FUNCTION LLFy(HL,HR,ZL,ZR,UL,UR,HvL,HvR,EIGmax)
    REAL(real_p), DIMENSION(p+1) :: FL, FR, UL, UR, LLFy
    REAL(real_p), DIMENSION(p+1) :: HvL, HvR, EIGmax
    REAL(real_p), DIMENSION(p+1) :: HL, HR, ZL, ZR
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, p+1
            FL(j)   = HvL(j)*UL(j) + half*g*(HL(j)**2.0d0 - ZL(j)**2.0d0) 
            FR(j)   = HvR(j)*UR(j) + half*g*(HR(j)**2.0d0 - ZR(j)**2.0d0) 
            LLFy(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HvR(j)-HvL(j))) 
        ENDDO
    END FUNCTION LLFy
    
    FUNCTION LLFx(UL,UR,HvL,HvR,EIGmax)
    REAL(real_p), DIMENSION(p+1) :: FL, FR, UL, UR, LLFx
    REAL(real_p), DIMENSION(p+1) :: HvL, HvR, EIGmax
    REAL(real_p), PARAMETER :: half = 0.50d0
        DO j = 1, p+1
            FL(j)   = HvL(j)*UL(j) 
            FR(j)   = HvR(j)*UR(j)  
            LLFx(j) = half*(FL(j) + FR(j) - abs(EIGmax(j))*(HvR(j)-HvL(j))) 
        ENDDO
    END FUNCTION LLFx  
    
    FUNCTION EIG(HL,HR,UL,UR)
    REAL(real_p), DIMENSION(p+1) :: UL,UR, HL, HR
    REAL(real_p), DIMENSION(p+1) :: LEp, LEm, REp, REm, EIG
    DO j = 1, p+1
        LEp(j) = UL(j) + sqrt(g*HL(j))
        LEm(j) = UL(j) - sqrt(g*HL(j))
        REp(j) = UR(j) + sqrt(g*HR(j))
        REm(j) = UR(j) - sqrt(g*HR(j))
        EIG(j) = max(abs(LEp(j)),abs(LEm(j)))
        EIG(j) = max(abs(EIG(j)),abs(REp(j)))
        EIG(j) = max(abs(EIG(j)),abs(REm(j)))
        EIG(j) = max(abs(EIG(j)),abs(UL(j)))
        EIG(j) = max(abs(EIG(j)),abs(UR(j)))
    ENDDO    
    END FUNCTION EIG     

END SUBROUTINE RHS_Vbar2D_NL
!-----------------------------------------------------------------
!
!                      Calculate DW/Dsigma  
!                  Written by Colton J. Conroy  
!                         @ the C.H.I.L
!                            12.17.13
!
!-----------------------------------------------------------------
SUBROUTINE CALC_DQW()

USE precisions
USE global
IMPLICIT NONE


INTEGER(int_p) :: jL, jR, nfacesX, nfacesY, nfacesZ, hL, hR, H_Elem
REAL(real_p), DIMENSION(NPTS_2D) :: UL, UR, Fhat
REAL(real_p), DIMENSION(NPTS_2D) :: HUL, HUR, ZL, ZR
REAL(real_p), DIMENSION(NPTS_3D) :: Uj, Hj, HUj, HUxj, HVyj, HVj, Vj
REAL(real_p), DIMENSION(NDOF_3D) :: Rxj, Ryj, Rj
REAL(real_p) :: INVJAC, g, s, h0, xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb
REAL(real_p) :: FL, FR, Havg
REAL(real_p), PARAMETER :: half = 0.50d0
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode

nfacesX  = size(xFACE,1)
nfacesY  = size(yFACE,1)
nfacesZ  = size(zFACE,1)
DQW      = 0.0d0
g        = PROB_DATA%g

UL   = 0.0d0
UR   = 0.0d0
Fhat = 0.0d0
HUL  = 0.0d0
HUR  = 0.0d0
ZL   = 0.0d0
ZR   = 0.0d0
Uj   = 0.0d0
Hj   = 0.0d0
HUj  = 0.0d0
HUxj = 0.0d0
HVyj = 0.0d0
HVj  = 0.0d0
Vj   = 0.0d0
Rxj  = 0.0d0
Ryj  = 0.0d0
Rj   = 0.0d0
globalNODES_3D = 0.0d0
Xnode = 0.0d0
Znode = 0.0d0
Ynode = 0.0d0

!....DQw

!....Element integration

DO i = 1, HEXelems
    
    H_Elem = HEXzxy(i,1)
    HUj    = MATMUL(PHI_3D,Hu(:,i,irk))
    HVj    = MATMUL(PHI_3D,Hv(:,i,irk))
    Uj     = MATMUL(PHI_2D3D,USE(:,H_Elem,irk)) 
    Vj     = MATMUL(PHI_2D3D,VSE(:,H_Elem,irk)) 
    Hj     = MATMUL(PHI_2D3D,ZETA(:,H_Elem,irk)) 
    DO j = 1, NPTS_3D
        HUxj(j) = -Hj(j)*Uj(j) + Huj(j)
        HVyj(j) = -Hj(j)*Vj(j) + Hvj(j)
    ENDDO
    Rxj = MATMUL(A_3D(:,:,1),HUxj)*HEX_jacX(i,1)
    Ryj = MATMUL(A_3D(:,:,2),HVyj)*HEX_jacY(i,1)
    DQW(:,i) = Rxj + Ryj
 
ENDDO   

!....Loop over x-faces

DO i = 1, nfacesX
    
    jL = xFACE(i,1)
    jR = xFACE(i,2)
    hL = HEXzxy(jL,1)
    hR = HEXzxy(jR,1)
        
    IF (jL == -1) THEN
    
        globalNODES_3D = HEXconn(jR,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, NPTS_2D
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,1))*Xnode(1)+(1.0d0+bPTS_3D(j,1,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,1))*Ynode(1)+(1.0d0+bPTS_3D(j,2,1))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,1))*Znode(1)+(1.0d0+bPTS_3D(j,3,1))*Znode(2))
            CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
            Fhat(j) = Huse*Vu
        ENDDO
            
        DQW(:,jR) = DQW(:,jR) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
        
    ELSEIF (jR == 0) THEN
    
        globalNODES_3D = HEXconn(jL,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, NPTS_2D
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,2))*Xnode(1)+(1.0d0+bPTS_3D(j,1,2))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,2))*Ynode(1)+(1.0d0+bPTS_3D(j,2,2))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,2))*Znode(1)+(1.0d0+bPTS_3D(j,3,2))*Znode(2))
            CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
            Fhat(j) = Huse*Vu
        ENDDO   
             
        DQW(:,jL) = DQW(:,jL) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)
            
    ELSE
        
        ZL  = MATMUL(b_PHI_2D3Dx(:,:,2),ZETA(:,hL,irk))
        UL  = MATMUL(b_PHI_2D3Dx(:,:,2),USE(:,hL,irk))
        ZR  = MATMUL(b_PHI_2D3Dx(:,:,1),ZETA(:,hR,irk))
        UR  = MATMUL(b_PHI_2D3Dx(:,:,1),USE(:,hR,irk))
        HUL = MATMUL(bPHI_3D(:,:,2),Hu(:,jL,irk))
        HUR = MATMUL(bPHI_3D(:,:,1),Hu(:,jR,irk))
        DO j = 1, NPTS_2D
            FR      = -ZR(j)*UR(j)+HUR(j)
            FL      = -ZL(j)*UL(j)+HUL(j)
            Havg    = half*(ZL(j)+ZR(j))
            Fhat(j) = half*(FL + FR)
        ENDDO    
        DQW(:,jR) = DQW(:,jR) + HEX_jacX(jR,1)*MATMUL(B_3D(:,:,1),Fhat)
        DQW(:,jL) = DQW(:,jL) - HEX_jacX(jL,1)*MATMUL(B_3D(:,:,2),Fhat)
            
    ENDIF
        
ENDDO 

!....Loop over y-faces

DO i = 1, nfacesY
    
    jL = yFACE(i,1)
    jR = yFACE(i,2)
    hL = HEXzxy(jL,1)
    hR = HEXzxy(jR,1)
        
    IF (jL == -1) THEN
    
        globalNODES_3D = HEXconn(jR,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, NPTS_2D
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,5))*Xnode(1)+(1.0d0+bPTS_3D(j,1,5))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,5))*Ynode(1)+(1.0d0+bPTS_3D(j,2,5))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,5))*Znode(1)+(1.0d0+bPTS_3D(j,3,5))*Znode(2))
            CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
            Fhat(j) = Huse*Vv
        ENDDO    
        
        DQW(:,jR) = DQW(:,jR) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
        
    ELSEIF (jR == 0) THEN
    
        globalNODES_3D = HEXconn(jL,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, NPTS_2D
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,6))*Xnode(1)+(1.0d0+bPTS_3D(j,1,6))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,6))*Ynode(1)+(1.0d0+bPTS_3D(j,2,6))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,6))*Znode(1)+(1.0d0+bPTS_3D(j,3,6))*Znode(2))
            CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
            Fhat(j) = Huse*Vv
        ENDDO    
             
        DQW(:,jL) = DQW(:,jL) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
    ELSE
        
        ZL  = MATMUL(b_PHI_2D3Dy(:,:,2),ZETA(:,hL,irk))
        UL  = MATMUL(b_PHI_2D3Dy(:,:,2),VSE(:,hL,irk))
        ZR  = MATMUL(b_PHI_2D3Dy(:,:,1),ZETA(:,hR,irk))
        UR  = MATMUL(b_PHI_2D3Dy(:,:,1),VSE(:,hR,irk))
        HUL = MATMUL(bPHI_3D(:,:,6),Hv(:,jL,irk))
        HUR = MATMUL(bPHI_3D(:,:,5),Hv(:,jR,irk))
        DO j = 1, NPTS_2D
            FR      = -ZR(j)*UR(j)+HUR(j)
            FL      = -ZL(j)*UL(j)+HUL(j)
            Havg    = half*(ZL(j)+ZR(j))
            Fhat(j) = half*(FL + FR)
        ENDDO 
        DQW(:,jR) = DQW(:,jR) + HEX_jacY(jR,1)*MATMUL(B_3D(:,:,5),Fhat)
        DQW(:,jL) = DQW(:,jL) - HEX_jacY(jL,1)*MATMUL(B_3D(:,:,6),Fhat)
            
    ENDIF
        
ENDDO

RETURN
END SUBROUTINE CALC_DQW
!------------------------------------------------------------------
!
!                    Calculate U from Hu
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                          11.20.13
!
!------------------------------------------------------------------
SUBROUTINE CALC_U(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: H_elem, ELEM
REAL(real_p), ALLOCATABLE :: Hu_L2(:), H_L2(:), U_L2(:)

IF (ELEM == 2) THEN
    ALLOCATE(Hu_L2(L2_2D))
    ALLOCATE(H_L2(L2_2D))
    ALLOCATE(U_L2(L2_2D))
    Hu_L2 = 0.0d0
    H_L2  = 0.0d0
    U_L2  = 0.0d0
    DO i = 1, Qelems
        H_elem = CONNzx(i,1) 
        Hu_L2  = MATMUL(L2PHIz,Hu(:,i,irk))
        H_L2  = MATMUL(PHIzeta2,ZETA(:,H_elem,irk))
        DO j = 1, L2_2D
            U_L2(j) = Hu_L2(j)/H_L2(j)
        ENDDO
        U_2D(:,i,irk) = MATMUL(L2U,U_L2)
    ENDDO    
    DEALLOCATE(Hu_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(U_L2)
ELSEIF (ELEM == 3) THEN
    ALLOCATE(Hu_L2(L2_3D))
    ALLOCATE(H_L2(L2_3D))
    ALLOCATE(U_L2(L2_3D))
    Hu_L2 = 0.0d0
    H_L2  = 0.0d0
    U_L2  = 0.0d0
    DO i = 1, HEXelems
        H_elem = HEXzxy(i,1)
        Hu_L2  = MATMUL(L2PHI_3D,Hu(:,i,irk))
        H_L2 = MATMUL(PHIzeta_3D,ZETA(:,H_elem,irk))
        DO j = 1, L2_3D
            U_L2(j) = Hu_L2(j)/H_L2(j)
        ENDDO
        U_3D(:,i,irk) = MATMUL(L2_U3D,U_L2)
    ENDDO 
    DEALLOCATE(Hu_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(U_L2)         
ENDIF    
 
RETURN
END SUBROUTINE CALC_U
!------------------------------------------------------------------
!
!                   Calculate Uavg from HUavg
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                           3.20.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_Ubar(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM
REAL(real_p), ALLOCATABLE :: Hu_L2(:), H_L2(:), U_L2(:)

IF (ELEM == 3) THEN
    ALLOCATE(Hu_L2(L2_2D))
    ALLOCATE(H_L2(L2_2D))
    ALLOCATE(U_L2(L2_2D))
    Hu_L2 = 0.0d0
    H_L2  = 0.0d0
    U_L2  = 0.0d0    
    DO i = 1, Qelems 
        Hu_L2 = MATMUL(L2PHIz,Hu_bar(:,i,irk))
        H_L2  = MATMUL(L2PHIz,ZETA(:,i,irk))
        DO j = 1, L2_2D
            U_L2(j) = Hu_L2(j)/H_L2(j)
        ENDDO
        USE(:,i,irk) = MATMUL(L2U,U_L2)
    ENDDO  
    DEALLOCATE(Hu_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(U_L2)             
ENDIF    
 
RETURN
END SUBROUTINE CALC_Ubar
!------------------------------------------------------------------
!
!                   Calculate Vavg from HVavg
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                           3.20.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_Vbar(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM
REAL(real_p), ALLOCATABLE :: Hu_L2(:), H_L2(:), U_L2(:)

IF (ELEM == 3) THEN
    ALLOCATE(Hu_L2(L2_2D))
    ALLOCATE(H_L2(L2_2D))
    ALLOCATE(U_L2(L2_2D))
    Hu_L2 = 0.0d0
    H_L2  = 0.0d0
    U_L2  = 0.0d0     
    DO i = 1, Qelems 
        Hu_L2 = MATMUL(L2PHIz,Hv_bar(:,i,irk))
        H_L2  = MATMUL(L2PHIz,ZETA(:,i,irk))
        DO j = 1, L2_2D
            U_L2(j) = Hu_L2(j)/H_L2(j)
        ENDDO
        VSE(:,i,irk) = MATMUL(L2U,U_L2)
    ENDDO   
    DEALLOCATE(Hu_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(U_L2)            
ENDIF    
 
RETURN
END SUBROUTINE CALC_Vbar
!------------------------------------------------------------------
!
!                    Calculate V from Hv
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                          11.20.13
!
!------------------------------------------------------------------
SUBROUTINE CALC_V(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: H_elem, ELEM
REAL(real_p), ALLOCATABLE :: Hv_L2(:), H_L2(:), V_L2(:)

IF (ELEM == 2) THEN
    ALLOCATE(Hv_L2(L2_2D))
    ALLOCATE(H_L2(L2_2D))
    ALLOCATE(V_L2(L2_2D))
    HV_L2 = 0.0d0
    H_L2  = 0.0d0
    V_L2  = 0.0d0
    DO i = 1, Qelems
        H_elem = CONNzx(i,1) 
        Hv_L2  = MATMUL(L2PHIz,Hv(:,i,irk))
        H_L2  = MATMUL(PHIzeta2,ZETA(:,H_elem,irk))
        DO j = 1, L2_2D
            V_L2(j) = Hv_L2(j)/H_L2(j)
        ENDDO
        V_2D(:,i,irk) = MATMUL(L2U,V_L2)
    ENDDO    
    DEALLOCATE(HV_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(V_L2)
ELSEIF (ELEM == 3) THEN
    ALLOCATE(Hv_L2(L2_3D))
    ALLOCATE(H_L2(L2_3D))
    ALLOCATE(V_L2(L2_3D))
    HV_L2 = 0.0d0
    H_L2  = 0.0d0
    V_L2  = 0.0d0    
    DO i = 1, HEXelems
        H_elem = HEXzxy(i,1)
        Hv_L2  = MATMUL(L2PHI_3D,Hv(:,i,irk))
        H_L2 = MATMUL(PHIzeta_3D,ZETA(:,H_elem,irk))
        DO j = 1, L2_3D
            V_L2(j) = Hv_L2(j)/H_L2(j)
        ENDDO
        V_3D(:,i,irk) = MATMUL(L2_U3D,V_L2)
    ENDDO 
    DEALLOCATE(HV_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(V_L2)               
ENDIF    
     

RETURN
END SUBROUTINE CALC_V
!------------------------------------------------------------------
!
!                    Calculate SE from H
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                          2.14.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_SE(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM
REAL(real_p), ALLOCATABLE :: Hh_L2(:), H_L2(:), SE_L2(:)

   
IF (ELEM == 3) THEN
    ALLOCATE(Hh_L2(L2_2D))
    ALLOCATE(H_L2(L2_2D))
    ALLOCATE(SE_L2(L2_2D))
    Hh_L2 = 0.0d0
    H_L2  = 0.0d0
    SE_L2 = 0.0d0
    DO i = 1, Qelems
        Hh_L2  = MATMUL(PHIL2_h,Hh(:,i))
        H_L2 = MATMUL(L2PHIz,ZETA(:,i,irk))
        DO j = 1, L2_2D
            SE_L2(j) = H_L2(j) - Hh_L2(j)
        ENDDO
        ZETA_2D(:,i,irk) = MATMUL(L2U,SE_L2)
    ENDDO  
    DEALLOCATE(Hh_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(SE_L2)          
ENDIF    
     
RETURN
END SUBROUTINE CALC_SE
!------------------------------------------------------------------
!
!                      Calculate HUavg
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                          3.18.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_HUbar(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM
REAL(real_p), DIMENSION(L2_2D) :: Hu_L2, H_L2, u_L2
   
IF (ELEM == 3) THEN
    Hu_L2 = 0.0d0
    H_L2 = 0.0d0
    u_L2 = 0.0d0
    DO i = 1, Qelems
        u_L2 = MATMUL(L2PHIz,USE(:,i,irk))
        H_L2 = MATMUL(L2PHIz,ZETA(:,i,irk))
        DO j = 1, L2_2D
            Hu_L2(j) = H_L2(j)*u_L2(j)
        ENDDO
        Hu_bar(:,i,irk) = MATMUL(L2U,Hu_L2)
    ENDDO            
ENDIF    
     
RETURN
END SUBROUTINE CALC_HUbar
!------------------------------------------------------------------
!
!                      Calculate HVavg
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                          3.18.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_HVbar(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM
REAL(real_p), ALLOCATABLE :: Hu_L2(:), H_L2(:), u_L2(:)

   
IF (ELEM == 3) THEN
    ALLOCATE(Hu_L2(L2_2D))
    ALLOCATE(H_L2(L2_2D))
    ALLOCATE(u_L2(L2_2D))
    Hu_L2 = 0.0d0
    H_L2 = 0.0d0
    u_L2 = 0.0d0
    DO i = 1, Qelems
        u_L2  = MATMUL(L2PHIz,VSE(:,i,irk))
        H_L2 = MATMUL(L2PHIz,ZETA(:,i,irk))
        DO j = 1, L2_2D
            Hu_L2(j) = H_L2(j)*u_L2(j)
        ENDDO
        Hv_bar(:,i,irk) = MATMUL(L2U,Hu_L2)
    ENDDO   
    DEALLOCATE(Hu_L2)
    DEALLOCATE(H_L2)
    DEALLOCATE(u_L2)         
ENDIF    
     
RETURN
END SUBROUTINE CALC_HVbar
!------------------------------------------------------------------
!
!               Calculate Vertical Eddy Viscosity
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                          9.19.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_Nz(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM
REAL(real_p), ALLOCATABLE :: u_L2(:), N_L2(:), v_L2(:)
REAL(real_p) :: MAG, Wa, Kt

Kt = 2.0d-5
Wa = 1.0d-4
   
IF (ELEM == 3) THEN
    ALLOCATE(u_L2(L2_2D))
    ALLOCATE(N_L2(L2_2D))
    ALLOCATE(v_L2(L2_2D))
    u_L2 = 0.0d0
    N_L2 = 0.0d0
    v_L2 = 0.0d0
    DO i = 1, Qelems
        u_L2  = MATMUL(L2PHIz,USE(:,i,irk))
        v_L2  = MATMUL(L2PHIz,VSE(:,i,irk))
        DO j = 1, L2_2D
            MAG = u_L2(j)**2.0d0 + v_L2(j)**2.0d0
            MAG = SQRT(MAG)
            N_L2(j) = Kt*(MAG/Wa)
        ENDDO
        NX(:,i) = MATMUL(L2U,N_L2)
    ENDDO   
    DEALLOCATE(u_L2)
    DEALLOCATE(N_L2)
    DEALLOCATE(v_L2)         
ENDIF    
     
RETURN
END SUBROUTINE CALC_Nz
!------------------------------------------------------------------
!
!                        Calculate W
!                  Written by Colton J. Conroy
!                       @ the C.H.I.L
!                          1.17.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_W(solnT)

USE precisions
USE global
USE IFPORT
IMPLICIT NONE

INTEGER(int_p) :: Ldz, iii, k, HEX, n, REGION, Rcount
INTEGER(int_p) :: jj, V_elem, Locx, Locy, Locz, Witer
INTEGER(int_p), DIMENSION(LRK) :: zpts
REAL(real_p) :: zCOORD, xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb,L_inf,L_infOLD
REAL(real_p) :: cpuT, T1, T2, TIME1, RANval, solnT
!REAL(real_p), DIMENSION(TP_2D,LNrk) :: Wsoln
!REAL(real_p), DIMENSION(TP_2D,nrkW+1) :: WRK
REAL(real_p), DIMENSION(NDOF_2D,LNrk) :: Wsoln
REAL(real_p), DIMENSION(NDOF_2D,nrkW+1) :: WRK
REAL(real_p), DIMENSION(TP_3D) :: Wpt
REAL(real_p), DIMENSION(TP_2D) :: Horizontal_soln
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode

W_2D = 0.0d0
Ldz = nelemsZ*LRK+1
L_infOLD = 0.0d0
L_inf = 1e-12
REGION = 3
Z_level = FLOOR(Ldz/2.0d0)

Horizontal_soln = 0.0d0
globalNODES_3D  = 0.0d0
TIME1 = 0.0d0
Xnode = 0.0d0
Znode = 0.0d0
Ynode = 0.0d0


DO i = 1, Qelems

    zCOORD = z0
    Wsoln  = 0.0d0
    WRK    = 0.0d0
    
    DO n = 0, Ldz-1
        HEX = RK_HEXconn(n+1,i)
        DO irkW = 1, nrkW
            zpt = zCOORD + ctW(irkW)*dzRK(n+1)
            CALL RHS_W(HEX,zpt,n,solnT)
            DO j = 1, irkW
                WRK(:,irkW+1) = WRK(:,irkW+1) + alphaW(irkW,j)*WRK(:,j) &
                              + betaW(irkW,j)*dzRK(n+1)*RHSw(:,j)
            ENDDO
        ENDDO
        Wsoln(:,n+2)    = WRK(:,nrkW+1)
        WRK(:,1)        = WRK(:,nrkW+1)
        IF (n == Z_level) THEN
            W_2D(:,i) = WRK(:,nrkW+1)
        ENDIF
        WRK(:,2:nrkW+1) = 0.0d0
        RHSw(:,:)       = 0.0d0
        zCOORD = zCOORD + dzRK(n+1)
    ENDDO    

    DO j = 1, nelemsZ

        V_elem = HEXcolumn(j,i)
        zpts   = zrkCONN(j,:)   
        Wpt    = 0.0d0
        kk     = 0
        globalNODES_3D = HEXconn(V_elem,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)

        DO k = 1, LRK   
            Horizontal_soln = MATMUL(PHI_TP2D,Wsoln(:,zpts(k)))
            DO ii = 1, TP_2D
                kk  = kk + 1
                Wpt(kk) = Horizontal_soln(ii)
                !xpt = (1.0d0/2.0d0)*((1.0d0-TPTS_3D(kk,1))*Xnode(1)+(1.0d0+TPTS_3D(kk,1))*Xnode(2))
                !ypt = (1.0d0/2.0d0)*((1.0d0-TPTS_3D(kk,2))*Ynode(1)+(1.0d0+TPTS_3D(kk,2))*Ynode(2))
                !zpt = (1.0d0/2.0d0)*((1.0d0-TPTS_3D(kk,3))*Znode(1)+(1.0d0+TPTS_3D(kk,3))*Znode(2))
                !CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                !L_inf = abs(Wpt(kk)-Vw)
                !IF (L_inf > L_infOLD) THEN
                !    L_infOLD = L_inf
                !    Locx = xpt
                !    Locy = ypt
                !    Locz = zpt
                !ENDIF
            ENDDO                    
        ENDDO
        W_3D(:,V_elem) = MATMUL(L3w,Wpt)
    ENDDO    
    
!     CALL CPU_TIME(T2)
!     cpuT = T2 - T1
!     IF (cpuT >= 86400.0d0) THEN
!         CALL GLOBAL_OUTPUT(REGION)
!         tt = dcmplx(TIME1,0.0d0)
!         CALL COMPUTE_ERROR(tt,REGION)
!         STOP
!     ENDIF
    
ENDDO  

IF (KERNEL == 1) THEN
    var = 2
    CALL HYPERBOLIC_KERNEL_2D(KC_2D,W_2D,pw,Qelems,ielems_2D,wENHANCED_2D,stencil_2D,NDOF_2D,Qstencil,nelemsX,nelemsY,CONNielem,NPTS_K)    
    CALL KERNEL_2D_ERROR(solnT)
ENDIF

!write(*,*)'L_infOLD=',L_infOLD  
!write(*,*)'LocX=',Locx
!write(*,*)'LocY=',Locy
!write(*,*)'LocZ=',Locz
                    
RETURN
END SUBROUTINE CALC_W
!------------------------------------------------------------------
!
!                       RHS W Subroutine
!                  Written by Colton J. Conroy
!                        @ the C.H.I.L
!                            1.17.14
!
!------------------------------------------------------------------
SUBROUTINE RHS_W(HEX,zpt,n,solnT)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD
INTEGER(int_p) :: iii, HEX, REGION, n, jL, jR, hL, hR, Witer
INTEGER(int_p) :: hLL, hRR, jLL, jRR, switch, flux_type, peek
INTEGER(int_p), DIMENSION(2) :: edge, face
REAL(real_p) :: Lpt, xpt, ypt, zpt, Vu, Vv, Vw, DVw, Huse, Zb
REAL(real_p) :: FL, FR, g, CF, FRR, FLL, SE, solnT, FxmA, FymA
REAL(real_p) :: DSEx, DSEy, Uavg, Fpce, Fxm, Fym, Fw
REAL(real_p) :: grad_SE, Vavg, hpt, nn, k
REAL(real_p), DIMENSION(LRK) :: ZLL, ZRR, HULL, HURR, ULL, URR
REAL(real_p), DIMENSION(LRK) :: EIGpce, EIGcon, EIGmax, UbL, UbR, s
REAL(real_p), DIMENSION(TP_2D) :: HUj, HVj, Uj, Vj, Hj, HUxj, HVyj, Zj
REAL(real_p), DIMENSION(NDOF_2D) :: Rxj, Ryj
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode
REAL(real_p), DIMENSION(LRK) :: Fhat, UL, UR, ZL, ZR, HUL, HUR, WUL, WUR
REAL(real_p), PARAMETER :: half = 0.50d0


Lpt       = 2.0d0/dzRK(n+1)*(zpt-nodesZRK(n+1)) - 1.0d0
g         = PROB_DATA%g
CF        = 1.0d0/3.0d0
flux_type = 3
switch    = flux_type
  
ZLL  = 0.0d0
ZRR  = 0.0d0
HULL = 0.0d0 
HURR = 0.0d0
ULL  = 0.0d0
URR  = 0.0d0
EIGpce = 0.0d0
EIGcon = 0.0d0
EIGmax = 0.0d0
UbL = 0.0d0
UbR = 0.0d0
s   = 0.0d0
HUj = 0.0d0
HVj = 0.0d0
Uj  = 0.0d0
Zj  = 0.0d0
Vj  = 0.0d0
Hj  = 0.0d0
HUxj = 0.0d0
HVyj = 0.0d0
Rxj  = 0.0d0
Ryj  = 0.0d0
globalNODES_3D = 0.0d0
Xnode = 0.0d0
Ynode = 0.0d0
Znode = 0.0d0
Fhat  = 0.0d0
UL    = 0.0d0
UR    = 0.0d0
ZL    = 0.0d0
ZR    = 0.0d0
HUL   = 0.0d0
HUR   = 0.0d0
WUL   = 0.0d0
WUR   = 0.0d0 
  
!....Element calculations

IF (NON_LINEAR == 1) THEN
    HUj = MATMUL(PHIw(:,:,irkW),Hu(:,HEX,irk))
    HVj = MATMUL(PHIw(:,:,irkW),Hv(:,HEX,irk))
    Uj  = MATMUL(PHI_TP2Dh,USE(:,i,irk))
    Vj  = MATMUL(PHI_TP2Dh,VSE(:,i,irk))
    Hj  = MATMUL(PHI_TP2Dh,ZETA(:,i,irk))
    DO j = 1, TP_2D
        HUxj(j) = -Hj(j)*Uj(j) + HUj(j)
        HVyj(j) = -Hj(j)*Vj(j) + HVj(j)
    ENDDO
ELSE
    Zj  = MATMUL(hPHI_TP2Dh,Hh(:,i))
    Uj = MATMUL(PHIw(:,:,irkW),U_3D(:,HEX,irk))
    Vj = MATMUL(PHIw(:,:,irkW),V_3D(:,HEX,irk))    
    DO j = 1, TP_2D
        HUxj(j) =  Zj(j)*Uj(j)
        HVyj(j) =  Zj(j)*Vj(j)
    ENDDO
ENDIF           

Rxj = MATMUL(AwT(:,:,1),HUxj)*CONN_jacX(i,1)
Ryj = MATMUL(AwT(:,:,2),HVyj)*CONN_jacY(i,1)

RHSw(:,irkW) = Rxj + Ryj

!....Flux Calculations

!....x- edges and faces

edge(1) = QedgeX(i,1)
edge(2) = QedgeX(i,2)
face(1) = HEXfaceX(HEX,1)
face(2) = HEXfaceX(HEX,2)

DO j = 1, 2
    
    jL = xEDGE(edge(j),1)
    jR = xEDGE(edge(j),2)
    hL = xFACE(face(j),1)
    hR = xFACE(face(j),2)
    
    IF (j == 1) THEN
        jRR = xEDGE(edge(j+1),2)
        hRR = xFACE(face(j+1),2)
    ELSEIF (j == 2) THEN
        jLL = xEDGE(edge(j-1),1)
        hLL = xFACE(face(j-1),1)
    ENDIF
    
    IF (hL == -1) THEN
    
        globalNODES_3D = HEXconn(hR,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO ii = 1, LRK
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,1,1))*Xnode(1)+(1.0d0+bPTS_2DTP(ii,1,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,2,1))*Ynode(1)+(1.0d0+bPTS_2DTP(ii,2,1))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-Lpt)*nodesZRK(n+1)+(1.0d0+Lpt)*nodesZRK(n+2))
            IF (PROB == 1) THEN
                Fhat(ii) = 0.0d0
            ELSEIF (PROB == 2) THEN
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                Fhat(ii) = Huse*Vu
            ELSEIF (PROB == 3) THEN
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                Fhat(ii) = -Huse*Uavg + Huse*Vu
            ENDIF  
        ENDDO
            
        RHSw(:,irkW) = RHSw(:,irkW) + CONN_jacX(i,1)*MATMUL(BwT(:,:,1),Fhat)    
        
    ELSEIF (hR == 0 .or. hR == -2) THEN
    
        globalNODES_3D = HEXconn(hL,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        ZL  = MATMUL(hbPHI_2DTPh(:,:,2),Hh(:,jL))                        
        UL  = MATMUL(bPHIw(:,:,irkW,2),U_3D(:,hL,irk))         
        
        DO ii = 1, LRK
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,1,2))*Xnode(1)+(1.0d0+bPTS_2DTP(ii,1,2))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,2,2))*Ynode(1)+(1.0d0+bPTS_2DTP(ii,2,2))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-Lpt)*nodesZRK(n+1)+(1.0d0+Lpt)*nodesZRK(n+2))
            IF (PROB == 1) THEN
                xCOORD = dcmplx(PROB_DATA%xN, 0.0d0)
                TIME   = dcmplx(solnT,0.0d0)
                hpt = PROB_DATA%Hslope*PROB_DATA%xN**2
                zCOORD = dcmplx(hpt*Z(nelemsZ), 0.0d0)
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,hpt,nn,k)
                Fhat = hpt*Vu           
            ELSEIF (PROB == 2) THEN
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                Fhat(ii) = Huse*Vu
            ELSEIF (PROB == 3) THEN
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                Fhat(ii) = -Huse*Uavg + Huse*Vu
            ELSEIF (PROB == 0) THEN
                Fhat(ii) = ZL(ii)*UL(ii)     
            ENDIF 
        ENDDO
            
        RHSw(:,irkW) = RHSw(:,irkW) - CONN_jacX(i,1)*MATMUL(BwT(:,:,2),Fhat)  
        
    ELSE
        
        IF (NON_LINEAR == 1) THEN
        
            IF (hLL == -1 .and. flux_type == 2) THEN
                flux_type = 1
                switch = 2
            ELSEIF (hRR == 0 .and. flux_type == 2) THEN
                flux_type = 1
                switch = 2
            ENDIF

            IF (flux_type == 1) THEN
                UL  = MATMUL(bPHI_2DTPh(:,:,2),USE(:,jL,irk))
                UR  = MATMUL(bPHI_2DTPh(:,:,1),USE(:,jR,irk))
                ZL  = MATMUL(bPHI_2DTPh(:,:,2),ZETA(:,jL,irk))            
                ZR  = MATMUL(bPHI_2DTPh(:,:,1),ZETA(:,jR,irk))
                HUL = MATMUL(bPHIw(:,:,irkW,2),Hu(:,hL,irk))
                HUR = MATMUL(bPHIw(:,:,irkW,1),Hu(:,hR,irk))
            ELSEIF (flux_type == 2) THEN 
                UL   = MATMUL(bPHI_2DTPh(:,:,1),USE(:,jL,irk))
                ZL   = MATMUL(bPHI_2DTPh(:,:,1),ZETA(:,jL,irk)) 
                HUL  = MATMUL(bPHIw(:,:,irkW,1),Hu(:,hL,irk)) 
                UR   = MATMUL(bPHI_2DTPh(:,:,2),USE(:,jR,irk)) 
                ZR   = MATMUL(bPHI_2DTPh(:,:,2),ZETA(:,jR,irk))
                HUR  = MATMUL(bPHIw(:,:,irkW,2),Hu(:,hR,irk))           
                IF (j == 1) THEN
                    URR  = MATMUL(bPHI_2DTPh(:,:,2),USE(:,jRR,irk))   
                    ZRR  = MATMUL(bPHI_2DTPh(:,:,2),ZETA(:,jRR,irk))
                    HURR = MATMUL(bPHIw(:,:,irkW,2),Hu(:,hRR,irk))
                ELSEIF (j == 2) THEN
                    ULL  = MATMUL(bPHI_2DTPh(:,:,1),USE(:,jLL,irk))   
                    ZLL  = MATMUL(bPHI_2DTPh(:,:,1),ZETA(:,jLL,irk))
                    HULL = MATMUL(bPHIw(:,:,irkW,1),Hu(:,hLL,irk))    
                ENDIF     
            ELSEIF (flux_type == 3) THEN
                UbL = MATMUL(bPHI_2DTPh(:,:,2),USE(:,jL,irk))
                UbR = MATMUL(bPHI_2DTPh(:,:,1),USE(:,jR,irk))
                ZL  = MATMUL(bPHI_2DTPh(:,:,2),ZETA(:,jL,irk))            
                ZR  = MATMUL(bPHI_2DTPh(:,:,1),ZETA(:,jR,irk))
                HUL = MATMUL(bPHIw(:,:,irkW,2),Hu(:,hL,irk))
                HUR = MATMUL(bPHIw(:,:,irkW,1),Hu(:,hR,irk))
                UL  = MATMUL(bPHIw(:,:,irkW,2),U_3D(:,hL,irk))
                UR  = MATMUL(bPHIw(:,:,irkW,1),U_3D(:,hR,irk))                        
            ENDIF  
        
            IF (flux_type == 1) THEN
                EIGmax = 0.0d0
                Fhat = LLF(ZL,ZR,UL,UR,HuL,HuR,EIGmax) 
            ELSEIF (flux_type == 2) THEN
                Fhat = CF3(ZL,ZR,ZLL,ZRR,UL,UR,ULL,URR,HuL,HuR,HuLL,HuRR,CF,j)
            ELSEIF (flux_type == 3) THEN
                EIGpce = EIG(ZL,ZR,UbL,UbR)    
                EIGcon = EIG(ZL,ZR,UL,UR)    
                EIGmax = abs(EIGpce-EIGcon)
                Fhat   = LLF(ZL,ZR,UbL,UbR,HuL,HuR,EIGmax)
            ENDIF      
              
         ELSE
            
            ZL  = MATMUL(hbPHI_2DTPh(:,:,2),Hh(:,jL))            
            ZR  = MATMUL(hbPHI_2DTPh(:,:,1),Hh(:,jR))            
            UL  = MATMUL(bPHIw(:,:,irkW,2),U_3D(:,hL,irk))  
            UR  = MATMUL(bPHIw(:,:,irkW,1),U_3D(:,hR,irk))  ! PROBLEM IS HERE (could be hR or bPHIw)
            DO ii = 1, LRK
                Fhat(ii) = half*(ZL(ii)*UL(ii) + ZR(ii)*UR(ii))
            ENDDO
            
         ENDIF                  
!        DO ii = 1, LRK
!            IF (flux_type == 1) THEN
!                FL       = -ZL(ii)*UL(ii) + HUL(ii)                    
!                FR       = -ZR(ii)*UR(ii) + HUR(ii)
!                Fhat(ii) = half*(FL+FR) 
!            ELSEIF (flux_type == 2) THEN         
!                IF (j == 1) THEN                                  
!                    FL       = -ZL(ii)*UL(ii) + HUL(ii)
!                    FR       = -ZR(ii)*UR(ii) + HUR(ii)
!                    FRR      = -ZRR(ii)*URR(ii) + HURR(ii) 
!                    Fhat(ii) = CF*FL + FR - CF*FRR
!                ELSEIF (j == 2) THEN
!                    FLL      = -ZLL(ii)*ULL(ii) + HULL(ii)                
!                    FL       = -ZL(ii)*UL(ii) + HUL(ii)
!                    FR       = -ZR(ii)*UR(ii) + HUR(ii)                      
!                    Fhat(ii) = CF*FR + FL - CF*FLL
!                ENDIF    
!            ENDIF    
!        ENDDO         
        
        IF (j == 1) THEN    
            RHSw(:,irkW) = RHSw(:,irkW) + CONN_jacX(i,1)*MATMUL(BwT(:,:,1),Fhat)
        ELSEIF (j == 2) THEN
            RHSw(:,irkW) = RHSw(:,irkW) - CONN_jacX(i,1)*MATMUL(BwT(:,:,2),Fhat)   
        ENDIF
        
        IF (switch == 2) THEN
            flux_type = 2
        ENDIF
        
    ENDIF    
         
    
ENDDO    

!....y- edges and faces

edge(1) = QedgeY(i,1)
edge(2) = QedgeY(i,2)
face(1) = HEXfaceY(HEX,1)
face(2) = HEXfaceY(HEX,2)

DO j = 1, 2
    
    jL = yEDGE(edge(j),1)
    jR = yEDGE(edge(j),2)
    hL = yFACE(face(j),1)
    hR = yFACE(face(j),2)
    
    IF (j == 1) THEN
        jRR = yEDGE(edge(j+1),2)
        hRR = yFACE(face(j+1),2)
    ELSEIF (j == 2) THEN
        jLL = yEDGE(edge(j-1),1)
        hLL = yFACE(face(j-1),1)
    ENDIF    
    
    IF (hL == -1 .or. hL == -2) THEN
    
        globalNODES_3D = HEXconn(hR,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
                    
        ZR  = MATMUL(hbPHI_2DTPh(:,:,3),Hh(:,jR))             
        UR  = MATMUL(bPHIw(:,:,irkW,3),V_3D(:,hR,irk))        
        
        DO ii = 1, LRK
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,1,3))*Xnode(1)+(1.0d0+bPTS_2DTP(ii,1,3))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,2,3))*Ynode(1)+(1.0d0+bPTS_2DTP(ii,2,3))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-Lpt)*nodesZRK(n+1)+(1.0d0+Lpt)*nodesZRK(n+2))
            IF (PROB == 1) THEN
                Fhat(ii) = 0.0d0
            ELSEIF (PROB == 2) THEN
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                Fhat(ii) = Huse*Vv
            ELSEIF (PROB == 3) THEN
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                Fhat(ii) = -Huse*Vavg + Huse*Vv
            ELSEIF (PROB == 0) THEN
                Fhat(ii) = ZR(ii)*UR(ii)    
            ENDIF 
        ENDDO
            
        RHSw(:,irkW) = RHSw(:,irkW) + CONN_jacY(i,1)*MATMUL(BwT(:,:,3),Fhat)    
        
    ELSEIF (hR == 0) THEN
    
        globalNODES_3D = HEXconn(hL,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO ii = 1, LRK
            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,1,4))*Xnode(1)+(1.0d0+bPTS_2DTP(ii,1,4))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_2DTP(ii,2,4))*Ynode(1)+(1.0d0+bPTS_2DTP(ii,2,4))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-Lpt)*nodesZRK(n+1)+(1.0d0+Lpt)*nodesZRK(n+2))
            IF (PROB == 1) THEN
                Fhat = 0.0d0
            ELSEIF (PROB == 2) THEN
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                Fhat(ii) = Huse*Vv
            ELSEIF (PROB == 3) THEN
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                Fhat(ii) = -Huse*Vavg + Huse*Vv
            ENDIF
        ENDDO
            
        RHSw(:,irkW) = RHSw(:,irkW) - CONN_jacY(i,1)*MATMUL(BwT(:,:,4),Fhat)  
        
    ELSE
    
        IF (NON_LINEAR == 1) THEN
    
            IF (hLL == -1 .and. flux_type == 2) THEN
                flux_type = 1
                switch = 2
            ELSEIF (hRR == 0 .and. flux_type == 2) THEN
                flux_type = 1
                switch = 2
            ENDIF    
        
            IF (flux_type == 1) THEN
                UL  = MATMUL(bPHI_2DTPh(:,:,4),VSE(:,jL,irk))
                UR  = MATMUL(bPHI_2DTPh(:,:,3),VSE(:,jR,irk))
                ZL  = MATMUL(bPHI_2DTPh(:,:,4),ZETA(:,jL,irk))            
                ZR  = MATMUL(bPHI_2DTPh(:,:,3),ZETA(:,jR,irk))
                HUL = MATMUL(bPHIw(:,:,irkW,4),Hv(:,hL,irk))
                HUR = MATMUL(bPHIw(:,:,irkW,3),Hv(:,hR,irk))
            ELSEIF (flux_type == 2) THEN    
                UL   = MATMUL(bPHI_2DTPh(:,:,3),VSE(:,jL,irk))
                ZL   = MATMUL(bPHI_2DTPh(:,:,3),ZETA(:,jL,irk)) 
                HUL  = MATMUL(bPHIw(:,:,irkW,3),Hv(:,hL,irk)) 
                UR   = MATMUL(bPHI_2DTPh(:,:,4),VSE(:,jR,irk)) 
                ZR   = MATMUL(bPHI_2DTPh(:,:,4),ZETA(:,jR,irk))
                HUR  = MATMUL(bPHIw(:,:,irkW,4),Hv(:,hR,irk))           
                IF (j == 1) THEN
                    URR  = MATMUL(bPHI_2DTPh(:,:,4),VSE(:,jRR,irk))   
                    ZRR  = MATMUL(bPHI_2DTPh(:,:,4),ZETA(:,jRR,irk))
                    HURR = MATMUL(bPHIw(:,:,irkW,4),Hv(:,hRR,irk))
                ELSEIF (j == 2) THEN
                    ULL  = MATMUL(bPHI_2DTPh(:,:,3),VSE(:,jLL,irk))   
                    ZLL  = MATMUL(bPHI_2DTPh(:,:,3),ZETA(:,jLL,irk))
                    HULL = MATMUL(bPHIw(:,:,irkW,3),Hv(:,hLL,irk))    
                ENDIF 
            ELSEIF (flux_type == 3) THEN
                UbL  = MATMUL(bPHI_2DTPh(:,:,4),VSE(:,jL,irk))
                UbR  = MATMUL(bPHI_2DTPh(:,:,3),VSE(:,jR,irk))
                ZL   = MATMUL(bPHI_2DTPh(:,:,4),ZETA(:,jL,irk))            
                ZR   = MATMUL(bPHI_2DTPh(:,:,3),ZETA(:,jR,irk))
                HUL  = MATMUL(bPHIw(:,:,irkW,4),Hv(:,hL,irk))
                HUR  = MATMUL(bPHIw(:,:,irkW,3),Hv(:,hR,irk))  
                UL   = MATMUL(bPHIw(:,:,irkW,4),V_3D(:,hL,irk))
                UR   = MATMUL(bPHIw(:,:,irkW,3),V_3D(:,hR,irk))                  
            ENDIF
                    
!        DO ii = 1, LRK
!            IF (flux_type == 1) THEN
!                FL       = -ZL(ii)*UL(ii) + HUL(ii)                    
!                FR       = -ZR(ii)*UR(ii) + HUR(ii)
!                Fhat(ii) = half*(FL+FR) 
!            ELSEIF (flux_type == 2) THEN         
!                IF (j == 1) THEN                                  
!                    FL       = -ZL(ii)*UL(ii) + HUL(ii)
!                    FR       = -ZR(ii)*UR(ii) + HUR(ii)
!                    FRR      = -ZRR(ii)*URR(ii) + HURR(ii) 
!                    Fhat(ii) = CF*FL + FR - CF*FRR
!                ELSEIF (j == 2) THEN
!                    FLL      = -ZLL(ii)*ULL(ii) + HULL(ii)                
!                    FL       = -ZL(ii)*UL(ii) + HUL(ii)
!                    FR       = -ZR(ii)*UR(ii) + HUR(ii)                      
!                    Fhat(ii) = CF*FR + FL - CF*FLL
!                ENDIF 
!            ENDIF
!        ENDDO      

            IF (flux_type == 1) THEN
                EIGmax = 0.0d0
                Fhat = LLF(ZL,ZR,UL,UR,HuL,HuR,EIGmax) 
            ELSEIF (flux_type == 2) THEN
                Fhat = CF3(ZL,ZR,ZLL,ZRR,UL,UR,ULL,URR,HuL,HuR,HuLL,HuRR,CF,j)
            ELSEIF (flux_type == 3) THEN
                EIGpce = EIG(ZL,ZR,UbL,UbR)    
                EIGcon = EIG(ZL,ZR,UL,UR)    
                EIGmax = abs(EIGpce - EIGcon)
                Fhat   = LLF(ZL,ZR,UbL,UbR,HuL,HuR,EIGmax)
            ENDIF
            
        ELSE
        
            ZL  = MATMUL(hbPHI_2DTPh(:,:,4),Hh(:,jL))            
            ZR  = MATMUL(hbPHI_2DTPh(:,:,3),Hh(:,jR))            
            UL  = MATMUL(bPHIw(:,:,irkW,4),V_3D(:,hL,irk))  
            UR  = MATMUL(bPHIw(:,:,irkW,3),V_3D(:,hR,irk))
            DO ii = 1, LRK
                Fhat(ii) = half*(ZL(ii)*UL(ii) + ZR(ii)*UR(ii))
            ENDDO        

        ENDIF
        
        peek = HEXelems/4

!        IF (HEX == peek) THEN
!            write(*,*)Fhat
!        ENDIF             
            
        IF (j == 1) THEN    
            RHSw(:,irkW) = RHSw(:,irkW) + CONN_jacY(i,1)*MATMUL(BwT(:,:,3),Fhat)
        ELSEIF (j == 2) THEN
            RHSw(:,irkW) = RHSw(:,irkW) - CONN_jacY(i,1)*MATMUL(BwT(:,:,4),Fhat)   
        ENDIF
        
        IF (switch == 2) THEN
            flux_type = 2
        ENDIF
        
    ENDIF    
         
    
ENDDO  

contains

   FUNCTION LLF(EL,ER,UbL,UbR,HuL,HuR,EIGmax)
    INTEGER(int_p) :: ii
    REAL(real_p), DIMENSION(LRK) :: EL, ER, UbL, UbR, LLF
    REAL(real_p), DIMENSION(LRK) :: HuL, HuR, EIGmax
    REAL(real_p) :: FL, FR
        DO ii = 1, LRK
            FL      = -EL(ii)*UbL(ii) + HuL(ii) 
            FR      = -ER(ii)*UbR(ii) + HuR(ii) 
            LLF(ii) = (1.0d0/2.0d0)*(FL + FR - abs(EIGmax(ii))*(ER(ii)-EL(ii))) 
        ENDDO                                                                    
    END FUNCTION LLF

    FUNCTION EIG(ZL,ZR,UL,UR)
    REAL(real_p), DIMENSION(LRK) :: UL, UR, ZL, ZR, EIG
    REAL(real_p) ::  LEp, LEm, REp, REm
    INTEGER(int_p) :: ii
    DO ii = 1, LRK
        LEp    = UL(ii) + sqrt(g*abs(ZL(ii)))
        LEm    = UL(ii) - sqrt(g*abs(ZL(ii)))
        REp    = UR(ii) + sqrt(g*abs(ZR(ii)))
        REm    = UR(ii) - sqrt(g*abs(ZR(ii)))
        EIG(ii) = max(abs(LEp),abs(LEm))
        EIG(ii) = max(abs(EIG(ii)),abs(REp))
        EIG(ii) = max(abs(EIG(ii)),abs(REm))
        EIG(ii) = max(abs(EIG(ii)),abs(UL(ii)))
        EIG(ii) = max(abs(EIG(ii)),abs(UR(ii)))
    ENDDO    
    END FUNCTION EIG
    
   FUNCTION CF3(EL,ER,ZLL,ZRR,UbL,UbR,ULL,URR,HuL,HuR,HuLL,HuRR,CF,j)
    INTEGER(int_p) :: j, ii
    REAL(real_p), DIMENSION(LRK) :: EL, ER, ZLL, ZRR, UbL, UbR, CF3
    REAL(real_p), DIMENSION(LRK) :: HuL, HuR, HuLL, HuRR, EIGmax, ULL, URR
    REAL(real_p) :: FL, FR, CF
        DO ii = 1, LRK
            IF (j == 1) THEN                                  
                FL       = -EL(ii)*UL(ii) + HUL(ii)
                FR       = -ER(ii)*UR(ii) + HUR(ii)
                FRR      = -ZRR(ii)*URR(ii) + HURR(ii) 
                CF3(ii) = CF*FL + FR - CF*FRR
            ELSEIF (j == 2) THEN
                FLL      = -ZLL(ii)*ULL(ii) + HULL(ii)                
                FL       = -EL(ii)*UL(ii) + HUL(ii)
                FR       = -ER(ii)*UR(ii) + HUR(ii)                      
                CF3(ii) = CF*FR + FL - CF*FLL
            ENDIF  
        ENDDO
    END FUNCTION CF3   
       
END SUBROUTINE RHS_W
!------------------------------------------------------------------
!        
!           Calculate Vertical Velocity Exactly
!              Written by Colton J. Conroy  
!                    @ the C.H.I.L
!                       7.16.14
!
!------------------------------------------------------------------
SUBROUTINE W_EXACT()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM_3Da, ELEM_3Dv, jj, xy, Sp1, Sm1
INTEGER(int_p) :: r, qq, s
REAL(real_p) :: xpt, ypt, zpt, WhL2, Psp1, Psm1, CS1x, CS2x
REAL(real_p) :: DENs, DENr, CR1x, CR2x, CH1x, CH2x
REAL(real_p) :: CS1y, CS2y, CR1y, CR2y, CH1y, CH2y
REAL(real_p) :: JACx, JACy, JACz
REAL(real_p), DIMENSION(L2_3D) :: WL2
!REAL(real_p), DIMENSION(8) :: globalNODES_3D
!REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode

JACx = 2.0d0/dx
JACy = 2.0d0/dy
JACz = dz/2.0d0
        
DO i = 1, Qelems

    DO j = 1, nelemsZ
    
        ELEM_3Da = HEXcolumn(j,i)        
        
        !globalNODES_3D = HEXconn(ELEM_3Da,:)
        !Xnode(1) = HEXnode(globalNODES_3D(1),1)
        !Xnode(2) = HEXnode(globalNODES_3D(3),1)
        !Znode(1) = HEXnode(globalNODES_3D(1),3)
        !Znode(2) = HEXnode(globalNODES_3D(2),3)
        !Ynode(1) = HEXnode(globalNODES_3D(1),2)
        !Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        WL2 = 0.0d0
        
        DO ii = 1, L2_3D
        
            !xpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(ii,1))*Xnode(1)+(1.0d0+L2PTS_3D(ii,1))*Xnode(2))
            !ypt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(ii,2))*Ynode(1)+(1.0d0+L2PTS_3D(ii,2))*Ynode(2))
            !zpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(ii,3))*Znode(1)+(1.0d0+L2PTS_3D(ii,3))*Znode(2))
            
            WhL2 = 0.0d0
            DENr = (L2PTS_3D(ii,3)+1.0d0)
            
            DO jj = 1, j
            
                ELEM_3Dv = HEXcolumn(jj,i)
                kk = 0                           !....Counter for 3D DOF
                
                DO s = 1, pu+1                   !....Sigma DOF
                    
                    xy   = 0                     !....Counter for 2D DOF        
                    Sp1  = s+1
                    Sm1  = s-1           
                    DENs = 1.0d0/(2.0d0*(Sm1)+1.0d0)   
                    
                    DO qq = 1, p+1               !....Y DOF
                    
                        DO r = 1, p+1            !....X DOF
                        
                            kk = kk + 1
                            xy = xy + 1
                         
                            IF (jj == j) THEN
                                
                                IF (Sm1 == 0) THEN
                                    Psm1 = -1.0d0
                                ELSE
                                    Psm1 = PHIs_z(ii,Sm1,1)
                                ENDIF
                                
                                Psp1  = PHIs_z(ii,Sp1,1)    
                                CS1x  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,1)*U_3D(kk,ELEM_3Dv,irk)
                                CS2x  = DENs*(Psp1 - Psm1)*PHIr(ii,r,2)*PHIq(ii,qq,1)*U_3D(kk,ELEM_3Dv,irk)*JACx
                                CR1x  = DENr*PHIr(ii,r,1)*PHIq(ii,qq,1)*USE(xy,i,irk)
                                CR2x  = DENr*PHIr(ii,r,2)*PHIq(ii,qq,1)*USE(xy,i,irk)*JACx
                                CH1x  = PHIr(ii,r,2)*PHIq(ii,qq,1)*ZETA(xy,i,irk)*JACx*JACz
                                CH2x  = PHIr(ii,r,1)*PHIq(ii,qq,1)*ZETA(xy,i,irk)*JACz
                                    
                                CS1y  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,1)*V_3D(kk,ELEM_3Dv,irk)
                                CS2y  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,2)*V_3D(kk,ELEM_3Dv,irk)*JACy
                                CR1y  = DENr*PHIr(ii,r,1)*PHIq(ii,qq,1)*VSE(xy,i,irk)
                                CR2y  = DENr*PHIr(ii,r,1)*PHIq(ii,qq,2)*VSE(xy,i,irk)*JACy
                                CH1y  = PHIr(ii,r,1)*PHIq(ii,qq,2)*ZETA(xy,i,irk)*JACy*JACz
                                CH2y  = PHIr(ii,r,1)*PHIq(ii,qq,1)*ZETA(xy,i,irk)*JACz   
                                
                            ELSE
                            
                                IF (Sm1 == 0) THEN
                                    Psm1 = -1.0d0
                                ELSE
                                    Psm1 = PHIs_1(Sm1)                                 
                                ENDIF  
                                
                                Psp1  = PHIs_1(Sp1)      
                                CS1x  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,1)*U_3D(kk,ELEM_3Dv,irk)
                                CS2x  = DENs*(Psp1 - Psm1)*PHIr(ii,r,2)*PHIq(ii,qq,1)*U_3D(kk,ELEM_3Dv,irk)*JACx
                                CR1x  = 2.0d0*PHIr(ii,r,1)*PHIq(ii,qq,1)*USE(xy,i,irk)
                                CR2x  = 2.0d0*PHIr(ii,r,2)*PHIq(ii,qq,1)*USE(xy,i,irk)*JACx
                                CH1x  = PHIr(ii,r,2)*PHIq(ii,qq,1)*ZETA(xy,i,irk)*JACx*JACz
                                CH2x  = PHIr(ii,r,1)*PHIq(ii,qq,1)*ZETA(xy,i,irk)*JACz
                                    
                                CS1y  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,1)*V_3D(kk,ELEM_3Dv,irk)
                                CS2y  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,2)*V_3D(kk,ELEM_3Dv,irk)*JACy
                                CR1y  = 2.0d0*PHIr(ii,r,1)*PHIq(ii,qq,1)*VSE(xy,i,irk)
                                CR2y  = 2.0d0*PHIr(ii,r,1)*PHIq(ii,qq,2)*VSE(xy,i,irk)*JACy
                                CH1y  = PHIr(ii,r,1)*PHIq(ii,qq,2)*ZETA(xy,i,irk)*JACy*JACz
                                CH2y  = PHIr(ii,r,1)*PHIq(ii,qq,1)*ZETA(xy,i,irk)*JACz   
                                           
                            ENDIF
                            
                            WhL2 = WhL2 + CH1x*(CR1x - CS1x) + CH2x*(CR2x - CS2x)  &
                                        + CH1y*(CR1y - CS1y) + CH2y*(CR2y - CS2y)
                                                                    
                        ENDDO                
                        
                    ENDDO    
                    
                ENDDO
                
            ENDDO    
            
            WL2(ii) = WhL2
                
        ENDDO
        
        W_3D(:,ELEM_3Da) = MATMUL(L2_U3D,WL2)
            
    ENDDO
    
ENDDO

    
RETURN
END SUBROUTINE W_EXACT
!------------------------------------------------------------------
!        
!              Calculates Baroclinic Head
!              Written by Colton J. Conroy  
!                    @ the C.H.I.L
!                       8.18.14
!
!------------------------------------------------------------------
SUBROUTINE CALC_R()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM_3Da, ELEM_3Dv, jj, xy, Sp1, Sm1
INTEGER(int_p) :: r, qq, s
REAL(real_p) :: xpt, RhL2, Psp1, Psm1, rSc
REAL(real_p) :: DENs, CR1, JACz, BC, r0, rS, Sr, Lx
REAL(real_p), DIMENSION(L2_3D) :: RL2
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(2) :: Xnode

globalNODES_3D = 0.0d0
Xnode = 0.0d0
JACz  = dz/2.0d0
r0    = PROB_DATA%r0
rS    = PROB_DATA%rS
Sr    = PROB_DATA%S
Lx    = PROB_DATA%Lx
rSc   = (rS/4.0d0)    !....For Lynch/Loder test case Boundary Condition
        
DO i = 1, Qelems

    DO j = 1, nelemsZ
    
        ELEM_3Da = HEXcolumn(j,i)        
        
        globalNODES_3D = HEXconn(ELEM_3Da,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        
        RL2 = 0.0d0

        DO ii = 1, L2_3D
        
            xpt  = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(ii,1))*Xnode(1)+(1.0d0+L2PTS_3D(ii,1))*Xnode(2))
            BC   = Sr*(rSc*(1.0d0 + cos(PI*xpt/Lx)))/r0 
            RhL2 = 0.0d0
            
            DO jj = 1, j

                ELEM_3Dv = HEXcolumn(jj,i)
                kk = 0                           !....Counter for 3D DOF

               DO s = 1, pu+1                   !....Sigma DOF
                  
                    xy   = 0                     !....Counter for 2D DOF        
                    Sp1  = s+1
                    Sm1  = s-1           
                    DENs = 1.0d0/(2.0d0*(Sm1)+1.0d0)   

                    DO qq = 1, p+1               !....Y DOF
 
                        DO r = 1, p+1            !....X DOF

                            kk = kk + 1
                            xy = xy + 1

                            IF (jj == j) THEN
                                
                                IF (Sm1 == 0) THEN
                                    Psm1 = -1.0d0
                                ELSE
                                    Psm1 = PHIs_z(ii,Sm1,1)
                                ENDIF
                                
                                Psp1 = PHIs_z(ii,Sp1,1)    
                                CR1  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,1)*JACz*PROB_DATA%h0  !....Check h0 part
           
                            ELSE
      
                                IF (Sm1 == 0) THEN
                                    Psm1 = -1.0d0
                                ELSE
                                    Psm1 = PHIs_1(Sm1)                                 
                                ENDIF  
                                
                                Psp1 = PHIs_1(Sp1)      
                                CR1  = DENs*(Psp1 - Psm1)*PHIr(ii,r,1)*PHIq(ii,qq,1)*JACz*PROB_DATA%h0 !....Check h0 part
         
                            ENDIF
                            
                            RhL2 = RhL2 - CR1*RHO(kk,ELEM_3Dv) 
                                                                    
                        ENDDO                
                        
                    ENDDO    
                    
                ENDDO
                
            ENDDO    
            
            RL2(ii) = RhL2 + BC
                
        ENDDO
        
        Rh(:,ELEM_3Da) = MATMUL(L2_U3D,RL2)
            
    ENDDO
    
ENDDO

    
RETURN
END SUBROUTINE CALC_R
!------------------------------------------------------------------
!        
!           Calculate Depth Averaged Velocity Exactly
!                  Written by Colton J. Conroy  
!                          @ the C.H.I.L
!                              9.16.13
!
!------------------------------------------------------------------
SUBROUTINE DEPTH_VEL_EXACT(REGION)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: REGION, H_Elem, L
REAL(real_p) :: INVJAC

USE(:,:,irk) = 0.0d0
VSE(:,:,irk) = 0.0d0
INVJAC = 1.0d0/nelemsZ

IF (REGION == 2) THEN

    DO L = 1, p+1

        DO i = 1, Qelems
        
            H_Elem = CONNzx(i,1)
            USE(L,H_Elem,irk) = USE(L,H_Elem,irk) + U_2D(L,i,irk)*INVJAC
            VSE(L,H_Elem,irk) = VSE(L,H_Elem,irk) + V_2D(L,i,irk)*INVJAC
        
        ENDDO
    
    ENDDO
    
ELSEIF (REGION == 3) THEN

    DO L = 1, NDOF
    
        DO i = 1, HEXelems
        
            H_Elem = HEXzxy(i,1)
            USE(L,H_Elem,irk) = USE(L,H_Elem,irk) + U_3D(L,i,irk)*INVJAC
            VSE(L,H_Elem,irk) = VSE(L,H_Elem,irk) + V_3D(L,i,irk)*INVJAC
            
        ENDDO
        
    ENDDO
    
ENDIF    

RETURN
END SUBROUTINE DEPTH_VEL_EXACT             
!------------------------------------------------------------------
!
!           Calculate Depth Averaged Velocities Subroutine
!                  Written by Colton J. Conroy
!                  @ the Isaac Newton Institute
!                             12.19.12
!
!------------------------------------------------------------------
SUBROUTINE DEPTH_VEL(REGION)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: k, REGION, H_Elem
REAL(real_p) :: INVJAC, Uint, Vint
REAL(real_p), DIMENSION(L2int) :: INT 

INT = 0.0d0

IF (REGION == 1) THEN

    DO i = 1, nelemsX

        DO j = 1, LE

	        DO k = 1, nelemsZ

	            INVJAC = 1.0d0/zELEM_jac(k)
	            INT    = MATMUL(L2PHIz,U(:,k,i,j,irk))
                Uz(k)  = dot_product(L2wts,int)*INVJAC

	        ENDDO

	        Ubar(j,i) = SUM(Uz)
        
        ENDDO

    ENDDO
    
ELSEIF (REGION == 2) THEN
    
    Ubar = 0.0d0
    DO i = 1, Qelems
    
        INVJAC = 1.0d0/CONN_jacZ(i,1)
        H_Elem = CONNzx(i,1)
        
        DO j = 1, LE
        
            INT  = MATMUL(PHIint(:,:,j),U_2D(:,i,irk))
            Uint = DOT_PRODUCT(L2wts,INT)*INVJAC 
            Ubar(j,H_Elem) = Ubar(j,H_Elem) + Uint

        ENDDO
    
    ENDDO    
    

    Vbar = 0.0d0
    DO i = 1, Qelems
        
        INVJAC = 1.0d0/CONN_jacZ(i,1)
        H_Elem = CONNzx(i,1)
            
        DO j = 1, L2_2D
            
            INT  = MATMUL(VPHIint(:,:,j),V_2D(:,i,irk)) 
            Vint = DOT_PRODUCT(L2wts,INT)*INVJAC 
            Vbar(j,H_Elem) = Vbar(j,H_Elem) + Vint
                
        ENDDO    
               
    ENDDO   
         
ELSEIF (REGION == 3) THEN

    Ubar = 0.0d0
    DO i = 1, HEXelems
    
        INVJAC = 1.0d0/HEX_jacZ(i,1)
        H_Elem = HEXzxy(i,1)
        
        DO j = 1, L2_2D
        
            INT            = MATMUL(PHIint_3D(:,:,j),U_3D(:,i,irk))
            Uint           = DOT_PRODUCT(L2wts,INT)*INVJAC
            Ubar(j,H_Elem) = Ubar(j,H_Elem) + Uint
            
        ENDDO
        
    ENDDO  
    
    Vbar = 0.0d0
    DO i = 1, HEXelems
    
        INVJAC = 1.0d0/HEX_jacZ(i,1)
        H_Elem = HEXzxy(i,1)
        
        DO j = 1, L2_2D
        
            INT            = MATMUL(PHIint_3D(:,:,j),V_3D(:,i,irk))
            Vint           = DOT_PRODUCT(L2wts,INT)*INVJAC
            Vbar(j,H_Elem) = Vbar(j,H_Elem) + Vint
            
        ENDDO
        
    ENDDO        
    
ENDIF    



IF (debug == 14) THEN
    DO i = 1, Qelems
	    write(*,*)'Ubar',Ubar(:,i)
    ENDDO
ENDIF

RETURN
END SUBROUTINE DEPTH_VEL 
!-------------------------------------------------------------------
!
!                   Surface Gradient Subroutine
!                   Written by Colton J. Conroy
!                   @ the Isaac Newton Institute
!                            12.17.12
!
!------------------------------------------------------------------
SUBROUTINE RHS_dZETA(ELEM)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ELEM

IF (ELEM == 1) THEN

    DO i = 1, nelemsX

        dZETA(:,i) = MATMUL(L3DPHI,ZETA(:,i,irk))*xELEM_jac(i)

    ENDDO
    
ELSEIF (ELEM == 2) THEN

    DO i = 1, nelemsX
    
        dZETA(:,i) = MATMUL(DPHIzeta,ZETA(:,i,irk))*xELEM_jac(i)
      
    ENDDO

ELSEIF (ELEM == 3) THEN

    IF (NON_LINEAR == 1) THEN 

        DO i = 1, Qelems
        
            dZETA(:,i) = MATMUL(DPHIzeta_3D(:,:,1),ZETA_2D(:,i,irk))*CONN_jacX(i,1)
        
        ENDDO
    
    ELSE

        DO i = 1, Qelems
        
            dZETA(:,i) = MATMUL(DPHIzeta_3D(:,:,1),ZETA(:,i,irk))*CONN_jacX(i,1)
        
        ENDDO
      
    ENDIF    
        
ENDIF

   
IF (debug == 11) THEN
    IF (ELEM == 1 .or. ELEM == 2) THEN
        DO i = 1, nelemsX
            WRITE(*,*)dZETA(:,i)
        ENDDO
    ELSEIF (ELEM == 3) THEN
        DO i = 1, Qelems
            WRITE(*,*)dZETA(:,i)
        ENDDO
    ENDIF
ENDIF

END SUBROUTINE RHS_dZETA
!-------------------------------------------------------------------
!
!                       dHdX P^0 Subroutine
!                   Written by Colton J. Conroy
!                   @ the Isaac Newton Institute
!                            4.3.14
!
!------------------------------------------------------------------
SUBROUTINE dHdX_P0(H_Elem,dHdX,solnT)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: jL, jR, H_Elem
INTEGER(int_p), DIMENSION(2) :: edge, face
REAL(real_p), DIMENSION(2) :: Zuse
REAL(real_p), DIMENSION(p+1) :: ZL
REAL(real_p) :: dxh, dHdX,Fw,FxmA,FymA
REAL(real_p) :: SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym
REAL(real_p) :: xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode

Zuse  = 0.0d0
ZL    = 0.0d0
edge  = 0.0d0
face  = 0.0d0
Xnode = 0.0d0
Ynode = 0.0d0
Znode = 0.0d0
globalNODES = 0.0d0

edge(1) = QedgeX(H_Elem,1)
edge(2) = QedgeX(H_Elem,2)

DO ii = 1, 2
    
    jL = xEDGE(edge(ii),1)
    jR = xEDGE(edge(ii),2)
  
    IF (jL == -1) THEN
        
        IF (PROB == 1) THEN
        
            
        ELSEIF (PROB == 3) THEN
            
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,1,1))*Xnode(1)+(1.0d0+bPTSx(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSx(j,2,1))*Ynode(1)+(1.0d0+bPTSx(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                Zuse(ii)   = Zb
            ENDDO
                                                      
        ENDIF    
            
!    ELSEIF (jR == 0) THEN
!        
!        IF (PROB == 1) THEN
!          
!  
!                
!        ELSEIF (PROB == 3) THEN
!            
!            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,2))*Xnode(1)+(1.0d0+bPTS_3D(j,1,2))*Xnode(2))
!            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,2))*Ynode(1)+(1.0d0+bPTS_3D(j,2,2))*Ynode(2))
!            zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,2))*Znode(1)+(1.0d0+bPTS_3D(j,3,2))*Znode(2))
!            CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
!
!            Zuse(ii)  = Zb
!                                               
!        ENDIF          
                 
    ELSE
        
        ZL       = MATMUL(bPHI_2Dxh(:,:,2),Hh(:,jL))
        Zuse(ii) = ZL(1)
                    
    ENDIF
        
ENDDO

dxh  = (2.0d0)*(1.0d0/CONN_jacX(H_Elem,1))
dHdX = (Zuse(2)-Zuse(1))/dxh

RETURN
END SUBROUTINE dHdX_P0    
!-------------------------------------------------------------------
!
!                       dHdY P^0 Subroutine
!                   Written by Colton J. Conroy
!                   @ the Isaac Newton Institute
!                            4.14.14
!
!------------------------------------------------------------------
SUBROUTINE dHdY_P0(H_Elem,dHdY,solnT)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: jL, jR, H_Elem
INTEGER(int_p), DIMENSION(2) :: edge, face
REAL(real_p), DIMENSION(2) :: Zuse
REAL(real_p), DIMENSION(p+1) :: ZL
REAL(real_p) :: dyh, dHdY,Fw,FxmA,FymA
REAL(real_p) :: SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym
REAL(real_p) :: xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode

Zuse  = 0.0d0
ZL    = 0.0d0
edge  = 0.0d0
face  = 0.0d0
Xnode = 0.0d0
Ynode = 0.0d0
Znode = 0.0d0
globalNODES = 0.0d0

edge(1) = QedgeY(H_Elem,1)
edge(2) = QedgeY(H_Elem,2)

DO ii = 1, 2
    
    jL = yEDGE(edge(ii),1)
    jR = yEDGE(edge(ii),2)
  
    IF (jL == -1) THEN
        
        IF (PROB == 1) THEN
        
            
        ELSEIF (PROB == 3) THEN
            
            globalNODES = CONN(jR,:)
            Xnode(1) = Qnode(globalNODES(1),1)
            Xnode(2) = Qnode(globalNODES(3),1)
            Ynode(1) = Qnode(globalNODES(1),2)
            Ynode(2) = Qnode(globalNODES(2),2)
            DO j = 1, p+1    
                xpt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,1,1))*Xnode(1)+(1.0d0+bPTSy(j,1,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-bPTSy(j,2,1))*Ynode(1)+(1.0d0+bPTSy(j,2,1))*Ynode(2))  
                zpt = Z(1) 
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                Zuse(ii)   = Zb
            ENDDO
                                                      
        ENDIF    
            
!    ELSEIF (jR == 0) THEN
!        
!        IF (PROB == 1) THEN
!          
!  
!                
!        ELSEIF (PROB == 3) THEN
!            
!            xpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,1,2))*Xnode(1)+(1.0d0+bPTS_3D(j,1,2))*Xnode(2))
!            ypt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,2,2))*Ynode(1)+(1.0d0+bPTS_3D(j,2,2))*Ynode(2))
!            zpt = (1.0d0/2.0d0)*((1.0d0-bPTS_3D(j,3,2))*Znode(1)+(1.0d0+bPTS_3D(j,3,2))*Znode(2))
!            CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
!
!            Zuse(ii)  = Zb
!                                               
!        ENDIF          
                 
    ELSE
        
        ZL       = MATMUL(bPHI_2Dyh(:,:,2),Hh(:,jL))
        Zuse(ii) = ZL(1)
                    
    ENDIF
        
ENDDO

dyh  = (2.0d0)*(1.0d0/CONN_jacY(H_Elem,1))
dHdY = (Zuse(2)-Zuse(1))/dyh

RETURN
END SUBROUTINE dHdY_P0  
!--------------------------------------------------------------------
!
!                     L2 Projection Subroutine
!                   Written by Colton J. Conroy
!                   @ the Isaac Newton Institute
!                           12.13.12
!
!--------------------------------------------------------------------
SUBROUTINE L2_PROJECTION(ELEM)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD, yCOORD
INTEGER(int_p) :: ELEM
REAL(real_p) :: SE, grad_SE, Vv, Vavg, INVJAC, xpt, INVJACz, zpt
REAL(real_p) :: Vu, Vw, ypt, Huse, nn , k, Zb, DVw, solnT, Uavg
REAL(real_p) :: Fpce, Fxm, Fym, Fw, FxmA, FymA, DSEx, DSEy
REAL(real_p) :: Ax1, Rpt, DRDXpt, RHOhpt, RHOpt
REAL(real_p), ALLOCATABLE :: L2SE(:), L2Vavg(:), L2UU(:), L2VV(:)
REAL(real_p), ALLOCATABLE :: L2Uavg(:), Nuse(:), Kuse(:)
REAL(real_p), DIMENSION(L2_2D) :: L2U_2D, L2V_2D, L2H
REAL(real_p), DIMENSION(L2_3D) :: L2U_3D, L2V_3D, L2W_3D, L2RHO_3D
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode

!....L2 projection of bathymetry

IF (ELEM == 3) THEN
    CALL L2_BATHY()
ENDIF    

!....L2 projection for depth-averaged continuity and momentum equations

TIME   = dcmplx(0.0d0, 0.0d0)
solnT  = 0.0d0
L2U_2D = 0.0d0
L2V_2D = 0.0d0
L2H    = 0.0d0
L2U_3D = 0.0d0
L2V_3D = 0.0d0
L2W_3D = 0.0d0
Xnode  = 0.0d0
Ynode  = 0.0d0
Znode  = 0.0d0
globalNODES    = 0.0d0
globalNODES_3D = 0.0d0

!-----------------------------------------------------------------------------------------------------
IF (ELEM == 1 .OR. ELEM == 2) THEN  !....Lines or Quads for Momentum so 1D for Depth Int. Eqns
!-----------------------------------------------------------------------------------------------------

    ALLOCATE(L2SE(L2int))
    ALLOCATE(L2Vavg(L2int))
    ALLOCATE(Kuse(L2int))
    ALLOCATE(Nuse(L2int))
    L2SE   = 0.0d0
    L2Vavg = 0.0d0
    Kuse   = 0.0d0
    Nuse   = 0.0d0
    
    DO i = 1, nelemsX                                                   !....Loop over elem.

        INVJAC = 1/xELEM_jac(i)
                                               
        DO j = 1, L2int                                                 !....Loop over L2pts
                                     
            xpt = X(i) + INVJAC*(1.0d0+L2pts(j,1))                          !....Transform pts
            Hexact(j,i) = PROB_DATA%Hslope*xpt**2
            zCOORD = dcmplx(Hexact(j,i)*Z(1), 0.0d0)
            Huse = Hexact(j,i)
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real 
                                                                    !....variables to complex
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)   !....Compute exact soln
                                                                    !....@ the L2pts
            L2SE(j)      = SE
            L2Vavg(j)    = Vavg
            Nuse(j)      = nn
            Kuse(j)      = k
            
        ENDDO

        ZETA(:,i,1) = MATMUL(L2,L2SE)                                   !....L2 projection Zeta
        USE(:,i,1)  = MATMUL(L2,L2Vavg)                                 !....L2 projection Ubar    
        Hh(:,i)     = MATMUL(L2,Hexact(:,i))
        NX(:,i)     = MATMUL(L2,Nuse)
        KX(:,i)     = MATMUL(L2,Kuse)        
        
    ENDDO
    DEALLOCATE(L2SE)
    DEALLOCATE(L2Vavg)
    DEALLOCATE(Kuse)
    DEALLOCATE(Nuse)   
     
!--------------------------------------------------------------------------------------------------------    
ELSEIF (ELEM == 3) THEN  !....Hexs for Momentum so Quads for Depth Int. Eqns
!--------------------------------------------------------------------------------------------------------
    
    ALLOCATE(L2SE(L2_2D))
    ALLOCATE(L2Vavg(L2_2D))    
    ALLOCATE(L2Uavg(L2_2D))  
    ALLOCATE(Kuse(L2_2D))
    ALLOCATE(Nuse(L2_2D))  
    L2SE   = 0.0d0
    L2Vavg = 0.0d0
    L2Uavg = 0.0d0
    Kuse   = 0.0d0
    Nuse   = 0.0d0

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
        Ynode(1) = Qnode(globalNODES(1),2)
        Ynode(2) = Qnode(globalNODES(2),2)
       
        DO j = 1, L2_2D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
            xCOORD = dcmplx(xpt,0.0d0)
            yCOORD = dcmplx(ypt,0.0d0)
            
            IF (PROB == 0) THEN
            
                L2SE(j) = 0.0d0
                L2Vavg(j) = 0.0d0
                L2Uavg(j) = 0.0d0                
                Kuse(j) = PROB_DATA%k
                Nuse(j) = 0.1d-8

            ELSEIF (PROB == 1) THEN
            
                Hexact(j,i) = PROB_DATA%Hslope*xpt**2
                zCOORD = dcmplx(Hexact(j,i)*Z(1), 0.0d0)
                Huse = Hexact(j,i)            
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                L2SE(j)      = SE
                L2Uavg(j)    = Vavg
                L2Vavg(j)    = 0.0d0
                Kuse(j)      = k
                Nuse(j)      = nn
                L2H(j)       = Huse
                
            ELSEIF (PROB == 2) THEN
                
                zpt = Z(1)
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                L2SE(j)   = Huse
                L2Vavg(j) = 0.0d0
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2H(j)    = Zb
                
            ELSEIF (PROB == 3) THEN
            
                zpt = Z(nelemsZ)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2SE(j)   = Huse
                L2Vavg(j) = Vavg
                L2Uavg(j) = Uavg
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2H(j)    = Zb     
                
            ELSEIF (PROB == 4) THEN
            
                IF (xpt <= 0.0d0) THEN
                    L2SE(j) = 10.0d0
                ELSE
                    L2SE(j) = 5.0d0
                ENDIF
                L2Vavg(j) = 0.0d0
                L2Uavg(j) = 0.0d0   
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2H(j)    = 0.0d0  
                
            ELSEIF (PROB == 5) THEN
            
                IF (i == 1 .and. j == 1) THEN          
                    CALL LL3D_CLOSURE(xCOORD,Ax1,Vavg)
                ENDIF
                zCOORD = dcmplx(Z(1), 0.0d0)
                CALL LL3D_SOLUTIONS(xCOORD,zCOORD,Vu,Vv,DSEx,Ax1,Rpt,DRDXpt,RHOhpt,RHOpt)
                L2SE(j)   = DSEx  
                L2Vavg(j) = 0.0d0 
                L2Uavg(j) = 0.0d0
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2H(j)    = PROB_DATA%h0                
                         
            ENDIF
                                
        ENDDO
        
        IF (PROB == 5) THEN
            GSExh(:,i)   = MATMUL(L2U,L2SE)
            ZETA(:,i,1) = 0.0d0
        ELSE
            ZETA(:,i,1) = MATMUL(L2U,L2SE)                               !....L2 projection Zeta
        ENDIF
        USE(:,i,1)  = MATMUL(L2U,L2Uavg)                                 !....L2 projection Ubar 
        VSE(:,i,1)  = MATMUL(L2U,L2Vavg)
        IF (PROB /= 0) THEN
            Hh(:,i)     = MATMUL(L2Bh,L2H)
        ENDIF    
        NX(:,i)     = MATMUL(L2U,Nuse)
        KX(:,i)     = MATMUL(L2U,Kuse)   
          
        
    ENDDO
    
    DEALLOCATE(L2SE)
    DEALLOCATE(L2Vavg)    
    DEALLOCATE(L2Uavg)  
    DEALLOCATE(Kuse)
    DEALLOCATE(Nuse)    

ENDIF

!....L2 projection for momentum eqn

!-------------------------------------------------------------------------------------------------------

IF (ELEM == 1) THEN  !....LINES

!-------------------------------------------------------------------------------------------------------

    ALLOCATE(L2UU(L2int))
    ALLOCATE(L2VV(L2int))
    L2UU = 0.0d0
    L2VV = 0.0d0
    
    DO i = 1, nelemsX

        INVJAC = 1/xELEM_jac(i)
  
        DO j = 1, LE
        
            xpt    = X(i) + INVJAC*(1+L3PTS(j,1))
            xCOORD = dcmplx(xpt,0.0d0) 
            Huse   = PROB_DATA%Hslope*xpt**2
            
	        DO kk = 1, nelemsZ

	            INVJACz = 1/zELEM_jac(kk)

	            DO ii = 1, L2int

	                zpt    = Z(kk) + INVJACz*(1.0d0+L2pts(ii,1))
		            zCOORD = dcmplx(Huse*zpt,0.0d0)
		            
		            IF (PROB == 1) THEN
	        
		                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)   !....Compute exact soln
                        L2UU(ii) = Vu
                        

                    ENDIF                         
                
                ENDDO

                U(:,kk,i,j,1) = MATMUL(L2U,L2UU)
                IF (PROB == 2) THEN
                     V(:,kk,i,j,1) = MATMUL(L2U,L2VV)
                ENDIF                 
                
	        ENDDO

        ENDDO

    ENDDO
    DEALLOCATE(L2UU)
    DEALLOCATE(L2VV)
    
!-----------------------------------------------------------------------------------------------------------

ELSEIF (ELEM == 2) THEN  !....Quads

!-----------------------------------------------------------------------------------------------------------
    
    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1) = Qnode(globalNODES(1),1)
        Xnode(2) = Qnode(globalNODES(3),1)
        Znode(1) = Qnode(globalNODES(1),2)
        Znode(2) = Qnode(globalNODES(2),2)
       
        DO j = 1, L2_2D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Znode(1)+(1.0d0+L2PTSz(j,2))*Znode(2))
            Huse   = PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt,0.0d0)
            zCOORD = dcmplx(Huse*zpt,0.0d0)
            
            IF (PROB == 1) THEN
            
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                L2U_2D(j) = Vu
                    
            ENDIF                     
        ENDDO
        
        IF (NON_LINEAR == 0) THEN
        
            U_2D(:,i,1) = MATMUL(L2U,L2U_2D)  
            
        ELSEIF (NON_LINEAR == 1) THEN  
        
            Hu(:,i,1)   = MATMUL(L2U,L2U_2D)
            
        ENDIF    
                  

    ENDDO

    IF(debug == 7) THEN
        write(*,*)'U=',U_2D(:,:,1)
    ENDIF

!------------------------------------------------------------------------------------------------------------
ELSEIF (ELEM == 3) THEN
!------------------------------------------------------------------------------------------------------------

    DO i = 1, HEXelems
        
        globalNODES_3D = HEXconn(i,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, L2_3D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,1))*Xnode(1)+(1.0d0+L2PTS_3D(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,2))*Ynode(1)+(1.0d0+L2PTS_3D(j,2))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,3))*Znode(1)+(1.0d0+L2PTS_3D(j,3))*Znode(2))
            
            IF (PROB == 0) THEN
            
                L2U_3D(j) = 1.0d-7
                L2V_3D(j) = 1.0d-7
            
            ELSEIF (PROB == 1) THEN
            
                Huse   = PROB_DATA%Hslope*xpt**2
                xCOORD = dcmplx(xpt,0.0d0)
                yCOORD = dcmplx(ypt,0.0d0)
                zCOORD = dcmplx(Huse*zpt,0.0d0)              
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                L2U_3D(j) = Vu
                L2V_3D(j) = 0.0d0
                
            ELSEIF (PROB == 2) THEN
            
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                L2U_3D(j) = Huse*Vu
                L2V_3D(j) = Huse*Vv
                L2W_3D(j) = 0.0d0
                
            ELSEIF (PROB == 3) THEN
            
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2U_3D(j) = Huse*Vu
                L2V_3D(j) = Huse*Vv
                L2W_3D(j) = 0.0d0
                
            ELSEIF (PROB == 4) THEN
            
                L2U_3D(j) = 0.0d0
                L2V_3D(j) = 0.0d0
                L2W_3D(j) = 0.0d0   
                
            ELSEIF (PROB == 5) THEN
            
                xCOORD = dcmplx(xpt,0.0d0)
                zCOORD = dcmplx(zpt,0.0d0)  
                    
                CALL LL3D_SOLUTIONS(xCOORD,zCOORD,Vu,Vv,DSEx,Ax1,Rpt,DRDXpt,RHOhpt,RHOpt)  
                
                L2U_3D(j)   = 0.0d0
                L2V_3D(j)   = Vv
                L2RHO_3D(j) = RHOhpt
               
            ENDIF                     
        ENDDO       
        
        IF (NON_LINEAR == 0) THEN
        
            U_3D(:,i,1) = MATMUL(L2_U3D,L2U_3D)   
            V_3D(:,i,1) = MATMUL(L2_U3D,L2V_3D)
            
            IF (PROB == 5) THEN
                RHO(:,i) = MATMUL(L2_U3D,L2RHO_3D)
            ENDIF
              
        ELSEIF (NON_LINEAR == 1) THEN     
        
            Hu(:,i,1) = MATMUL(L2_U3D,L2U_3D) 
            Hv(:,i,1) = MATMUL(L2_U3D,L2V_3D)  
            !W_3D(:,i) = MATMUL(L2_U3D,L2W_3D) don't need b/c calc W from Hu and Hv    
            
        ENDIF       

    ENDDO
    
    IF(debug == 7) THEN
        write(*,*)'U=',U_3D(:,:,1)
    ENDIF

ENDIF    

RETURN
END SUBROUTINE L2_PROJECTION
!--------------------------------------------------------------------
!
!                  Compute L2 projection of Bathymetry
!                      Written by Colton J. Conroy
!                           @ the C.H.I.L.
!                               5.21.14
!
!-------------------------------------------------------------------
SUBROUTINE L2_BATHY()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: NPTS_BATHY, NPTS_B2D, k, REGION, meshH, dumb, NDOFh
INTEGER(int_p), PARAMETER :: REGIONs = 1, debug_check = 1
REAL(real_p) :: Mx, My
REAL(real_p), ALLOCATABLE :: BPTS_1D(:,:), BWTS_1D(:), BPTS_2D(:,:), BWTS_2D(:) 
REAL(real_p), ALLOCATABLE :: IL2(:,:), L2h(:), M(:,:)
REAL(real_p), ALLOCATABLE :: BASIS_2DH(:), DBASIS_2DH(:,:)


IF (FULL2D == 0) THEN                       
    NDOFh = (ph+2)*(ph+1)/2
ELSE       
    NDOFh = (ph+1)**2
ENDIF

IF (PROB == 0) THEN
    OPEN(10,FILE='fortH.xy')
    READ(10,'(A24)') meshH
    READ(10,*)dumb, NPTS_BATHY
    NPTS_B2D   = (NPTS_BATHY)**2
ELSE 
    NPTS_BATHY = L2int
    NPTS_B2D   = (NPTS_BATHY)**2
ENDIF         

ALLOCATE(Hh(NDOFh,Qelems))
ALLOCATE(BPTS_1D(NPTS_BATHY,1))
ALLOCATE(BWTS_1D(NPTS_BATHY))
ALLOCATE(BPTS_2D(NPTS_B2D,2))
ALLOCATE(BWTS_2D(NPTS_B2D))
ALLOCATE(L2PHI_h(NPTS_B2D,NDOFh))
ALLOCATE(DPHI_h(L2_2D,NDOFh,2))
ALLOCATE(PHIL2_h(L2_2D,NDOFh))
ALLOCATE(PHI_h(NPTS_2D,NDOFh))
ALLOCATE(M(NDOFh,NDOFh))
IF (PROB == 0) THEN
    ALLOCATE(IL2(NDOFh,NPTS_B2D))
    ALLOCATE(L2Bh(NDOFh,NPTS_B2D))
ELSE
    ALLOCATE(IL2(NDOFh,L2_2D))
    ALLOCATE(L2Bh(NDOFh,L2_2D))
ENDIF    
ALLOCATE(L2h(NPTS_B2D))
ALLOCATE(BASIS_2DH(NDOFh))
ALLOCATE(DBASIS_2DH(NDOFh,2))
ALLOCATE(bPHI_2Dyh(p+1,NDOFh,2))
ALLOCATE(bPHI_2Dxh(p+1,NDOFh,2))
ALLOCATE(PHI_2D3Dh(NPTS_3D,NDOFh))
ALLOCATE(b_PHI_2D3Dxh(NPTS_2D,NDOFh,2))
ALLOCATE(b_PHI_2D3Dyh(NPTS_2D,NDOFh,2))
ALLOCATE(DPHIzeta_3Dh(L2_3D,NDOFh,2))
ALLOCATE(hPHI_TP2Dh(TP_2D,NDOFh))
ALLOCATE(hbPHI_2DTPh(LRK,NDOFh,4))

Hh           = 0.0d0
BPTS_1D      = 0.0d0
BWTS_1D      = 0.0d0
BPTS_2D      = 0.0d0
BWTS_2D      = 0.0d0
BASIS_2DH    = 0.0d0
DBASIS_2DH   = 0.0d0
L2PHI_h      = 0.0d0
PHI_2D3Dh    = 0.0d0
DPHI_h       = 0.0d0
PHI_h        = 0.0d0
IL2          = 0.0d0
L2Bh         = 0.0d0
L2h          = 0.0d0
M            = 0.0d0
b_PHI_2D3Dxh = 0.0d0
b_PHI_2D3Dyh = 0.0d0
DPHIzeta_3Dh = 0.0d0
hPHI_TP2Dh   = 0.0d0
hbPHI_2DTPh  = 0.0d0

DIM = 1
CALL QUADRATURE(NPTS_BATHY,REGIONs,BPTS_1D,BWTS_1D,NPTS_BATHY,DIM)
CALL TENSOR_QUAD_2D(BPTS_2D,BWTS_2D,BPTS_1D,BWTS_1D,NPTS_BATHY)

DIM = 2
REGION = 2
DO i = 1, NPTS_B2D
    CALL ORTHOGONAL_BASIS(REGION, BPTS_2D(i,:), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        L2PHI_h(i,j)  = BASIS_2DH(j)
    ENDDO
ENDDO

DO i = 1, NPTS_2D
    CALL ORTHOGONAL_BASIS(REGION, elemPTSz(i,:), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        PHI_h(i,j)  = BASIS_2DH(j)
    ENDDO
ENDDO

DO i = 1, L2_3D
    CALL ORTHOGONAL_BASIS(REGION, L2PTS_3D(i,1:2), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        DPHIzeta_3Dh(i,j,1) = DBASIS_2DH(j,1)     !....Derivative of zeta wrt x @ 3D L2 pts
        DPHIzeta_3Dh(i,j,2) = DBASIS_2DH(j,2)     !....Derivative of zeta wrt y @ 3D L2 pts
    ENDDO
ENDDO

DO i = 1, TP_2D
    CALL ORTHOGONAL_BASIS(REGION, TPTS_2D(i,:), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        hPHI_TP2Dh(i,j) = BASIS_2DH(j)
    ENDDO
ENDDO 

DO i = 1, L2_2D
    CALL ORTHOGONAL_BASIS(REGION, L2PTSz(i,:), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        PHIL2_h(i,j)   = BASIS_2DH(j)
        DPHI_h(i,j,1)  = DBASIS_2DH(j,1)
        DPHI_h(i,j,2)  = DBASIS_2DH(j,2)
    ENDDO
ENDDO

DO i = 1, NPTS_3D
    CALL ORTHOGONAL_BASIS(REGION, elemPTS_3D(i,1:2), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        PHI_2D3Dh(i,j) = BASIS_2DH(j)
    ENDDO
ENDDO 

DO i = 1, p+1
    CALL ORTHOGONAL_BASIS(REGION, bPTSx(i,:,1), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        bPHI_2Dxh(i,j,1) = BASIS_2DH(j)                  !....Left boundary
    ENDDO
    CALL ORTHOGONAL_BASIS(REGION, bPTSx(i,:,2), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)   
    DO j = 1, NDOFh
        bPHI_2Dxh(i,j,2) = BASIS_2DH(j)                  !....Right boundary
    ENDDO         
ENDDO     
    
DO i = 1, p+1
    CALL ORTHOGONAL_BASIS(REGION, bPTSy(i,:,1), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        bPHI_2Dyh(i,j,1) = BASIS_2DH(j)                  !....Left boundary
    ENDDO
    CALL ORTHOGONAL_BASIS(REGION, bPTSy(i,:,2), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)   
    DO j = 1, NDOFh
        bPHI_2Dyh(i,j,2) = BASIS_2DH(j)                  !....Right boundary
    ENDDO         
ENDDO 

DO i = 1, NPTS_2D
    CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,1), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        b_PHI_2D3Dxh(i,j,1) = BASIS_2DH(j)
    ENDDO
    CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,2), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        b_PHI_2D3Dxh(i,j,2) = BASIS_2DH(j)
    ENDDO
ENDDO
    
DO i = 1, NPTS_2D
    CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,5), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        b_PHI_2D3Dyh(i,j,1) = BASIS_2DH(j)
    ENDDO
    CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,6), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        b_PHI_2D3Dyh(i,j,2) = BASIS_2DH(j)
    ENDDO
ENDDO

DO i = 1, LRK
    CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,1), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)
    DO j = 1, NDOFh
        hbPHI_2DTPh(i,j,1) = BASIS_2DH(j)                  !....left x-face
    ENDDO
    CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,2), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)   
    DO j = 1, NDOFh
        hbPHI_2DTPh(i,j,2) = BASIS_2DH(j)                  !....right x-face
    ENDDO    
    CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,3), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)   
    DO j = 1, NDOFh
        hbPHI_2DTPh(i,j,3) = BASIS_2DH(j)                  !....bottom y-face
    ENDDO               
    CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,4), ph+1, DIM, BASISz, DBASISz, BASIS_2DH, DBASIS_2DH)   
    DO j = 1, NDOFh
        hbPHI_2DTPh(i,j,4) = BASIS_2DH(j)                  !....top y-face
    ENDDO          
ENDDO

IF (FULL2D == 1) THEN
    k = 0
    DO i = 0, ph
        DO j = 0, ph
            k = k+1
            Mx = (2.0d0*j+1.0d0)/2.0d0       !....Inverse of Mass Matrix
            My = (2.0d0*i+1.0d0)/2.0d0
            M(k,k) = Mx*My
        ENDDO
    ENDDO
ELSEIF (FULL2D == 0) THEN
    IF (PROB == 0) THEN
        DO i = 1, NDOFh
            DO j = 1, NDOFh
                DO k = 1, NPTS_B2D
                    M(i,j) = M(i,j) + BWTS_2D(k)*L2PHI_h(k,i)*L2PHI_h(k,j)
                ENDDO
                IF (i == j) THEN
                    M(i,j) = 1.0d0/M(i,j)
                ENDIF    
            ENDDO
        ENDDO
    ELSE
        DO i = 1, NDOFh
            DO j = 1, NDOFh
                DO k = 1, L2_2D
                    M(i,j) = M(i,j) + L2WTSz(k)*PHIL2_h(k,i)*PHIL2_h(k,j)
                ENDDO
                IF (i == j) THEN
                    M(i,j) = 1.0d0/M(i,j)
                ENDIF    
            ENDDO
        ENDDO  
    ENDIF      
ENDIF        

!....L2 matrix

IF (PROB == 0) THEN
    DO j = 1, NPTS_B2D
        DO i = 1, NDOFh
            IL2(i,j) = BWTS_2D(j)*L2PHI_h(j,i)
        ENDDO
    ENDDO
ELSE
    DO j = 1, L2_2D
        DO i = 1, NDOFh
            IL2(i,j) = L2WTSz(j)*PHIL2_h(j,i)
        ENDDO
    ENDDO    
ENDIF     

L2Bh = MATMUL(M,IL2)

IF (PROB == 0) THEN
    DO i = 1, Qelems
        DO j = 1, NPTS_B2D
            READ(10,*)dumb, L2h(j)
        ENDDO
        Hh(:,i) = MATMUL(L2Bh,L2h)
    ENDDO  
    CLOSE(10)
ENDIF  

IF (debug_check == 1) THEN
    OPEN(11,FILE='Hh.2D')
    DO i = 1, Qelems
        DO j = 1, NDOFh
            WRITE(11,*)Hh(j,i)
        ENDDO  
    ENDDO
    CLOSE(11)
ENDIF


DEALLOCATE(BPTS_1D)
DEALLOCATE(BWTS_1D)
DEALLOCATE(BPTS_2D)
DEALLOCATE(BWTS_2D)
DEALLOCATE(M)
DEALLOCATE(IL2)
DEALLOCATE(L2h)
DEALLOCATE(BASIS_2DH)
DEALLOCATE(DBASIS_2DH)

RETURN
END SUBROUTINE L2_BATHY
!--------------------------------------------------------------------
!
!                   Compute Numerical Error Subroutine
!                       Written by Colton J. Conroy
!                       @ the Isaac Newton Institute
!                                12.14.12
!
!-------------------------------------------------------------------
SUBROUTINE COMPUTE_ERROR(TIME,ELEM)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD, yCOORD
INTEGER(int_p) :: ELEM, check, Qe
REAL(real_p) :: Ax1, Rpt, DRDXpt, RHOhpt, RHOpt, RHO_L2error, R_L2error, N_L2error
REAL(real_p) :: L2_error, INVJAC, ZETA_L2error, U_L2error, INVJACz, GSE_LINFerror
REAL(real_p) :: grad_SE, Vv, SE, Vavg, xpt, zpt, Vavg_L2error, L2v_error, point1, point2
REAL(real_p) :: g, h0, Vu, L2Vv_error, ypt, maxcheck, nn, k, Huse, L2h_error, H_L2error
REAL(real_p) :: L2w_error, Zb, Vw, DVw, L2dw_error, L2vse_error, FxmA, L2hu_error
REAL(real_p) :: solnT, Uavg, Fpce, Fxm, Fym, Fw, Uavg_L2error, DSEx, DSEy, FymA, L2hv_error
REAL(real_p), ALLOCATABLE :: L2SE(:), L2Vavg(:), L2ZETA(:), L2USE(:), L2DIF(:), L2DIFv(:)
REAL(real_p), ALLOCATABLE :: L2U_2D(:), L2V_2D(:), Ut_2D(:), Vt_2D(:), L2DIF2(:), L2DIFVvV(:)
REAL(real_p), ALLOCATABLE :: L2H(:), L2DIFh(:), Nuse(:), Kuse(:), Wt(:), L2DIFw(:), L2W(:)
REAL(real_p), ALLOCATABLE :: L2GSE2(:), DIF2(:,:), L2Hz(:), L2DW(:), DWt(:), L2DIFdw(:)
REAL(real_p), ALLOCATABLE :: L2VSE(:), L2DIFvse(:), L2Uavg(:), L2Hu(:), L2DIFhu(:), L2HUs(:)
REAL(real_p), ALLOCATABLE :: L2Hv(:), L2DIFhv(:), L2HVs(:), L2Nx(:), L2DIFNX(:)
REAL(real_p), DIMENSION(L2int) :: L2DIFVv, L2VV, L2Vt, L2UU, L2Ut
REAL(real_p), DIMENSION(LE) :: L2GSE
REAL(real_p), DIMENSION(LE,nelemsX) :: DIF
REAL(real_p), DIMENSION(L2_3D) :: L2RHO_3D, RHOt, L2DIFrho, L2R_3D, Rt, L2DIFr
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(8) :: globalNODES_3D
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode

!....L2 Errors

L2_error       = 0.000000000000000d0
L2v_error      = 0.000000000000000d0
L2h_error      = 0.000000000000000d0
L2vse_error    = 0.000000000000000d0
L2hu_error     = 0.000000000000000d0
L2hv_error     = 0.000000000000000d0
solnT          = REAL(TIME)
L2DIFVv        = 0.0d0 
L2VV           = 0.0d0
L2Vt           = 0.0d0
L2UU           = 0.0d0
L2Ut           = 0.0d0
L2GSE          = 0.0d0
DIF            = 0.0d0
globalNODES    = 0.0d0
globalNODES_3D = 0.0d0
Xnode          = 0.0d0
Ynode          = 0.0d0
Znode          = 0.0d0

!-----------------------------------------------------------------------------------------------------------
IF (ELEM == 1 .OR. ELEM == 2) THEN !....Depth Int. for 2D problem
!-----------------------------------------------------------------------------------------------------------

    ALLOCATE(L2SE(L2int))
    ALLOCATE(L2Vavg(L2int))
    ALLOCATE(L2ZETA(L2int))
    ALLOCATE(L2USE(L2int))
    ALLOCATE(L2DIF(L2int))
    ALLOCATE(L2DIFv(L2int))
    ALLOCATE(Nuse(L2int))
    ALLOCATE(Kuse(L2int))
    ALLOCATE(L2H(L2int))
    ALLOCATE(L2DIFh(L2int))
    L2SE   = 0.0d0
    L2Vavg = 0.0d0
    L2ZETA = 0.0d0
    L2USE  = 0.0d0
    L2DIF  = 0.0d0
    L2DIFv = 0.0d0
    Nuse   = 0.0d0
    Kuse   = 0.0d0
    L2H    = 0.0d0
    L2DIFh = 0.0d0

    DO i = 1, nelemsX
    
        INVJAC = 1/xELEM_jac(i)
    
        DO j = 1, L2int

            xpt = X(i) + INVJAC*(1.0d0+L2pts(j,1))                  !....Transform pts
            xCOORD = dcmplx(xpt,0.0d0)                              !....Convert real
            Hexact(j,i) = PROB_DATA%Hslope*xpt**2
            zCOORD = dcmplx(Hexact(j,i)*Z(1), 0.0d0)
            Huse = Hexact(j,i)
                                                                    !....variables to complex
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)

            L2SE(j)   = SE
            L2Vavg(j) = Vavg
            Nuse(j)   = nn
            Kuse(j)   = k

        ENDDO

        L2ZETA = MATMUL(L2PHI,ZETA(:,i,1))
        L2USE  = MATMUL(L2PHI,USE(:,i,1))
        L2H    = MATMUL(L2PHI,Hh(:,i))

        DO j = 1, L2int
   
            L2DIF(j)  = (L2SE(j) - L2ZETA(j))**2
            L2DIFv(j) = (L2Vavg(j) - L2USE(j))**2
            L2DIFh(j) = (Hexact(j,i) - L2H(j))**2

        ENDDO
    
        L2_error  = L2_error + DOT_PRODUCT(L2wts,L2DIF)*INVJAC
        L2v_error = L2v_error + DOT_PRODUCT(L2wts,L2DIFv)*INVJAC
        L2h_error = L2h_error + DOT_PRODUCT(L2wts,L2DIFh)*INVJAC
    
    ENDDO
    DEALLOCATE(L2SE)
    DEALLOCATE(L2Vavg)
    DEALLOCATE(L2ZETA)
    DEALLOCATE(L2USE)
    DEALLOCATE(L2DIF)
    DEALLOCATE(L2DIFv)
    DEALLOCATE(Nuse)
    DEALLOCATE(Kuse)
    DEALLOCATE(L2H)
    DEALLOCATE(L2DIFh)
        
!-----------------------------------------------------------------------------------------------------------    
ELSEIF (ELEM == 3) THEN   
!-----------------------------------------------------------------------------------------------------------
    
    ALLOCATE(L2SE(L2_2D))
    ALLOCATE(L2Vavg(L2_2D))
    ALLOCATE(L2Uavg(L2_2D))
    ALLOCATE(L2ZETA(L2_2D))
    ALLOCATE(L2USE(L2_2D))
    ALLOCATE(L2DIF(L2_2D))
    ALLOCATE(L2DIFv(L2_2D))
    ALLOCATE(Nuse(L2_2D))
    ALLOCATE(Kuse(L2_2D))
    ALLOCATE(L2H(L2_2D))
    ALLOCATE(L2Hz(L2_2D))
    ALLOCATE(L2DIFh(L2_2D))  
    ALLOCATE(L2VSE(L2_2D))  
    ALLOCATE(L2DIFvse(L2_2D))
    ALLOCATE(L2Hu(L2_2D))
    ALLOCATE(L2DIFhu(L2_2D))
    ALLOCATE(L2HUs(L2_2D))
    ALLOCATE(L2Hv(L2_2D))
    ALLOCATE(L2DIFhv(L2_2D))
    ALLOCATE(L2HVs(L2_2D))
    ALLOCATE(L2Nx(L2_2D))
    ALLOCATE(L2DIFNX(L2_2D))
    L2SE     = 0.0d0
    L2Vavg   = 0.0d0
    L2Uavg   = 0.0d0
    L2ZETA   = 0.0d0
    L2USE    = 0.0d0
    L2DIF    = 0.0d0
    L2DIFv   = 0.0d0
    Nuse     = 0.0d0
    Kuse     = 0.0d0
    L2H      = 0.0d0
    L2Hz     = 0.0d0
    L2DIFh   = 0.0d0
    L2VSE    = 0.0d0
    L2DIFvse = 0.0d0
    L2Hu     = 0.0d0
    L2DIFhu  = 0.0d0
    L2HUs    = 0.0d0
    L2Hv     = 0.0d0
    L2DIFhv  = 0.0d0
    L2HVs    = 0.0d0
    L2Nx     = 0.0d0
    L2DIFNX  = 0.0d0
    N_L2error = 0.0d0
    
    DO i = 1, Qelems

        globalNODES = CONN(i,:)
        Xnode(1)    = Qnode(globalNODES(1),1)
        Xnode(2)    = Qnode(globalNODES(3),1)
        Ynode(1)    = Qnode(globalNODES(1),2)
        Ynode(2)    = Qnode(globalNODES(2),2)
        
        DO j = 1, L2_2D
        
            point1 = L2PTSz(j,1)
            point2 = L2PTSz(j,2)
            xpt    = (1.0d0/2.0d0)*((1.0d0-point1)*Xnode(1)+(1.0d0+point1)*Xnode(2))
            ypt    = (1.0d0/2.0d0)*((1.0d0-point2)*Ynode(1)+(1.0d0+point2)*Ynode(2))
            Hexact(j,i) = PROB_DATA%Hslope*xpt**2
            zCOORD = dcmplx(Hexact(j,i)*Z(1), 0.0d0)
            xCOORD = dcmplx(xpt, 0.0d0)
            yCOORD = dcmplx(ypt, 0.0d0)
                        
            IF (PROB == 1) THEN
            
                Huse = Hexact(j,i)
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)

                L2SE(j)   = SE
                L2Vavg(j) = 0.0d0
                L2Uavg(j) = Vavg
                Nuse(j)   = nn
                Kuse(j)   = k
                L2Hz(j)   = Huse
                
            ELSEIF (PROB == 2) THEN
            
                zpt = Z(1)
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                L2SE(j)   = Huse
                L2Vavg(j) = 0.0d0
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2Hz(j)   = Zb
                
            ELSEIF (PROB == 3) THEN
            
                zpt = Z(nelemsZ)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)                
                L2SE(j)   = Huse
                L2Uavg(j) = Uavg
                L2Vavg(j) = Vavg
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2Hz(j)   = Zb 
                L2Hu(j)   = Huse*Uavg      
                L2Hv(j)   = Huse*Vavg  
                
            ELSEIF (PROB == 5) THEN
            
                IF (i == 1 .and. j == 1) THEN          
                    CALL LL3D_CLOSURE(xCOORD,Ax1,Vavg)
                ENDIF
                zCOORD = dcmplx(Z(1), 0.0d0)
                CALL LL3D_SOLUTIONS(xCOORD,zCOORD,Vu,Vv,DSEx,Ax1,Rpt,DRDXpt,RHOhpt,RHOpt)
                L2SE(j)   = 0.0d0   
                L2Vavg(j) = Vavg 
                L2Uavg(j) = 0.0d0
                Kuse(j)   = PROB_DATA%k
                Nuse(j)   = PROB_DATA%N0
                L2Hz(j)   = PROB_DATA%h0                                    
                
            ENDIF          
                         

        ENDDO

        L2ZETA = MATMUL(L2PHIz,ZETA(:,i,1))
        L2USE  = MATMUL(L2PHIz,USE(:,i,1))
        L2VSE  = MATMUL(L2PHIz,VSE(:,i,1))
        L2H    = MATMUL(PHIL2_h,Hh(:,i))
        L2HUs  = MATMUL(L2PHIz,Hu_bar(:,i,1))
        L2HVs  = MATMUL(L2PHIz,Hv_bar(:,i,1))
        L2Nx   = MATMUL(L2PHIz,NX(:,i))
        
        DO j = 1, L2_2D
   
            L2DIF(j)    = (L2SE(j) - L2ZETA(j))**2.0d0
            L2DIFv(j)   = (L2Uavg(j) - L2USE(j))**2.0d0
            L2DIFh(j)   = (L2Hz(j) - L2H(j))**2.0d0
            L2DIFvse(j) = (L2VSE(j) - L2Vavg(j))**2.0d0
            L2DIFhu(j)  = (L2Hu(j) - L2HUs(j))**2.0d0
            L2DIFhv(j)  = (L2Hv(j) - L2HVs(j))**2.0d0
            L2DIFNX(j)  = (L2Nx(j) - Nuse(j))**2.0d0
            
        ENDDO
    
        L2_error    = L2_error + DOT_PRODUCT(L2WTSz,L2DIF)*(1.0d0/CONN_jac(i,1))
        L2v_error   = L2v_error + DOT_PRODUCT(L2WTSz,L2DIFv)*(1.0d0/CONN_jac(i,1))
        L2h_error   = L2h_error + DOT_PRODUCT(L2WTSz,L2DIFh)*(1.0d0/CONN_jac(i,1))
        L2vse_error = L2vse_error + DOT_PRODUCT(L2WTSz,L2DIFvse)*(1.0d0/CONN_jac(i,1))
        L2hu_error  = L2hu_error + DOT_PRODUCT(L2WTSz,L2DIFhu)*(1.0d0/CONN_jac(i,1))
        L2hv_error  = L2hv_error + DOT_PRODUCT(L2WTSz,L2DIFhv)*(1.0d0/CONN_jac(i,1))
        N_L2error   = N_L2error + DOT_PRODUCT(L2WTSz,L2DIFNX)*(1.0d0/CONN_jac(i,1))
        
    ENDDO 
    
    DEALLOCATE(L2SE)
    DEALLOCATE(L2Vavg)
    DEALLOCATE(L2Uavg)
    DEALLOCATE(L2ZETA)
    DEALLOCATE(L2USE)
    DEALLOCATE(Nuse)
    DEALLOCATE(Kuse)
    DEALLOCATE(L2H)
    DEALLOCATE(L2Hz)
    DEALLOCATE(L2DIFh)  
    DEALLOCATE(L2VSE)  
    DEALLOCATE(L2DIFvse)
    DEALLOCATE(L2Hu)
    DEALLOCATE(L2DIFhu)
    DEALLOCATE(L2HUs)
    DEALLOCATE(L2Hv)
    DEALLOCATE(L2DIFhv)
    DEALLOCATE(L2HVs)    
    DEALLOCATE(L2Nx)
    DEALLOCATE(L2DIFNX)

ENDIF

ZETA_L2error = SQRT(L2_error)
Uavg_L2error = SQRT(L2v_error)
H_L2error    = SQRT(L2h_error)
Vavg_L2error = SQRT(L2vse_error)

!....Display numerical error

IF (NUM_ERROR == 1) THEN
    write(*,*)'ZETA L2_error = ', ZETA_L2error
    write(*,*)'Uavg L2_error = ', Uavg_L2error
    write(*,*)'Vavg L2_error = ', Vavg_L2error
    IF (coupled == 3 .and. NON_LINEAR == 1) THEN
        write(*,*)'HUavg L2_error =', SQRT(L2hu_error)
        write(*,*)'HVavg L2_error =', SQRT(L2hv_error)
    ENDIF    
    write(*,*)'H L2_error = ', H_L2error
    IF (PROB == 1) THEN
        write(*,*)'N L2_error=', N_L2error
    ENDIF
ENDIF

!....L2 Errors for Horizontal Velocity

g    = PROB_DATA%g
h0   = PROB_DATA%h0

L2DIF      = 0.000000000000000d0
L2_error   = 0.000000000000000d0
L2DIF      = 0.000000000000000d0
L2DIFv     = 0.000000000000000d0
L2_error   = 0.000000000000000d0
L2v_error  = 0.000000000000000d0
L2Vv_error = 0.00000000000000d0
L2w_error  = 0.00000000000000d0
L2dw_error = 0.00000000000000d0

!----------------------------------------------------------------------------------------------------
IF (ELEM == 1) THEN !....LINES
!----------------------------------------------------------------------------------------------------

    DO i = 1, nelemsX  ! Loop over x elements

        INVJAC = 1/xELEM_jac(i)  ! Transformation matrix

        DO j = 1, LE     ! Loop over Lines beneath int. points

            xpt    = X(i) + INVJAC*(1+L3PTS(j,1))
            xCOORD = dcmplx(xpt,0.0d0)
            Huse   = PROB_DATA%Hslope*xpt**2
	
            DO kk = 1, nelemsZ  ! Loop over Line elements

                INVJACz = 1/zELEM_jac(kk)

                DO ii = 1, L2int

                    zpt    = Z(kk) + INVJACz*(1+L2pts(ii,1))
                    zCOORD = dcmplx(Huse*zpt,0.0d0)
                    
                    IF (PROB == 1) THEN

                        CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)   !....Compute exact soln
                        L2UU(ii) = Vu

                    ENDIF        

                ENDDO
	    
                L2Ut = MATMUL(L2PHIz,U(:,kk,i,j,1))
                IF (PROB == 2) THEN
                    L2Vt = MATMUL(L2PHIz,V(:,kk,i,j,1))
                ENDIF    
	    
	            DO ii = 1, L2int

        	        L2DIF(ii) = (L2UU(ii) - L2Ut(ii))**2
        	        IF (PROB == 2) THEN
        	            L2DIFVv(ii) = (L2VV(ii) - L2Vt(ii))**2
        	        ENDIF    

    	        ENDDO

    	        L2_error = L2_error + DOT_PRODUCT(L2wts,L2DIF)*INVJACz
    	        IF (PROB == 2) THEN
    	            L2Vv_error = L2Vv_error + DOT_PRODUCT(L2wts,L2DIFVv)*INVJACz 
    	        ENDIF    

            ENDDO

        ENDDO

    ENDDO

!---------------------------------------------------------------------------------------------------------
ELSEIF (ELEM == 2) THEN ! Quads
!---------------------------------------------------------------------------------------------------------

    ALLOCATE(L2U_2D(L2_2D))
    ALLOCATE(L2V_2D(L2_2D))
    ALLOCATE(Ut_2D(L2_2D))
    ALLOCATE(Vt_2D(L2_2D))
    ALLOCATE(L2DIF2(L2_2D))
    ALLOCATE(L2DIFVvV(L2_2D))
    L2U_2D   = 0.0d0
    L2V_2D   = 0.0d0
    Ut_2D    = 0.0d0
    Vt_2D    = 0.0d0
    L2DIF2   = 0.0d0
    L2DIFVvV = 0.0d0
    
    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1)    = Qnode(globalNODES(1),1)
        Xnode(2)    = Qnode(globalNODES(3),1)
        Znode(1)    = Qnode(globalNODES(1),2)
        Znode(2)    = Qnode(globalNODES(2),2)

        DO j = 1, L2_2D
        
            point1 = L2PTSz(j,1)
            point2 = L2PTSz(j,2)
            xpt    = (1.0d0/2.0d0)*((1.0d0-point1)*Xnode(1)+(1.0d0+point1)*Xnode(2))
            zpt    = (1.0d0/2.0d0)*((1.0d0-point2)*Znode(1)+(1.0d0+point2)*Znode(2))
            Huse   = PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt, 0.0d0)
            zCOORD = dcmplx(Huse*zpt, 0.0d0)
            
            IF (PROB == 1) THEN

                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                L2U_2D(j) = Vu      
            
            ENDIF
            
        ENDDO
        
        Ut_2D = MATMUL(L2PHIz,U_2D(:,i,1))

        IF (PROB == 2) THEN
            Vt_2D = MATMUL(L2PHIz,V_2D(:,i,1))
        ENDIF
        
        DO ii = 1, L2_2D

            L2DIF2(ii) = (L2U_2D(ii) - Ut_2D(ii))**2
            IF (PROB == 2) THEN
                L2DIFVvV(ii) = (L2V_2D(ii) - Vt_2D(ii))**2
            ENDIF
                           
    	ENDDO
        IF (PROB == 1) THEN
    	    L2_error = L2_error + DOT_PRODUCT(L2WTSz,L2DIF2)*(1.0d0/CONN_jac(i,1))
    	ELSEIF (PROB == 2) THEN
    	    L2_error = L2_error + DOT_PRODUCT(L2WTSz,L2DIF2)*(1.0d0/CONN_jac(i,1))
    	    L2Vv_error = L2Vv_error + DOT_PRODUCT(L2WTSz,L2DIFVvV)*(1.0d0/CONN_jac(i,1))     	
        ENDIF
        
    ENDDO
    
    DEALLOCATE(L2U_2D)
    DEALLOCATE(L2V_2D)
    DEALLOCATE(Ut_2D)
    DEALLOCATE(Vt_2D)
    DEALLOCATE(L2DIF2)
    DEALLOCATE(L2DIFVvV)    
    
!----------------------------------------------------------------------------------------------------------
ELSEIF (ELEM == 3) THEN !....3D errors
!----------------------------------------------------------------------------------------------------------
    
    ALLOCATE(L2U_2D(L2_3D))
    ALLOCATE(L2V_2D(L2_3D))
    ALLOCATE(L2W(L2_3D))
    ALLOCATE(L2DW(L2_3D))
    ALLOCATE(Ut_2D(L2_3D))
    ALLOCATE(Vt_2D(L2_3D))
    ALLOCATE(Wt(L2_3D))
    ALLOCATE(DWt(L2_3D))
    ALLOCATE(L2DIF2(L2_3D))
    ALLOCATE(L2DIFVvV(L2_3D))
    ALLOCATE(L2DIFw(L2_3D))
    ALLOCATE(L2DIFdw(L2_3D))
    L2U_2D   = 0.0d0
    L2V_2D   = 0.0d0
    L2W      = 0.0d0
    L2DW     = 0.0d0
    Ut_2D    = 0.0d0
    Vt_2D    = 0.0d0
    Wt       = 0.0d0
    DWt      = 0.0d0
    L2DIF2   = 0.0d0
    L2DIFVvV = 0.0d0
    L2DIFw   = 0.0d0
    L2DIFdw  = 0.0d0
    L2RHO_3D = 0.0d0
    L2DIFrho = 0.0d0
    RHOt     = 0.0d0
    L2R_3D   = 0.0d0
    Rt       = 0.0d0
    L2DIFr   = 0.0d0
    RHO_L2error = 0.0d0
    R_L2error   = 0.0d0
    
    DO i = 1, HEXelems

        globalNODES_3D = HEXconn(i,:)
        Xnode(1) = HEXnode(globalNODES_3D(1),1)
        Xnode(2) = HEXnode(globalNODES_3D(3),1)
        Znode(1) = HEXnode(globalNODES_3D(1),3)
        Znode(2) = HEXnode(globalNODES_3D(2),3)
        Ynode(1) = HEXnode(globalNODES_3D(1),2)
        Ynode(2) = HEXnode(globalNODES_3D(5),2)
        
        DO j = 1, L2_3D
        
            xpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,1))*Xnode(1)+(1.0d0+L2PTS_3D(j,1))*Xnode(2))
            ypt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,2))*Ynode(1)+(1.0d0+L2PTS_3D(j,2))*Ynode(2))
            zpt = (1.0d0/2.0d0)*((1.0d0-L2PTS_3D(j,3))*Znode(1)+(1.0d0+L2PTS_3D(j,3))*Znode(2))
            Huse   = PROB_DATA%Hslope*xpt**2
            xCOORD = dcmplx(xpt,0.0d0)
            yCOORD = dcmplx(ypt,0.0d0)
            zCOORD = dcmplx(Huse*zpt,0.0d0)  
            
            IF (PROB == 1) THEN
            
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                L2U_2D(j) = Vu
                
            ELSEIF (PROB == 2) THEN 
            
                Huse = 0.0d0
                CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                L2U_2D(j) = Huse*Vu
                L2V_2D(j) = Huse*Vv  
                L2DW(j)   = DVw
                L2W(j)    = Vw    
                
            ELSEIF (PROB == 3) THEN
            
                Huse = 0.0d0
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                L2U_2D(j) = Huse*Vu
                L2V_2D(j) = Huse*Vv  
                L2W(j)    = Vw    
                
            ELSEIF (PROB == 5) THEN
            
                xCOORD = dcmplx(xpt,0.0d0)
                zCOORD = dcmplx(zpt,0.0d0)  
                    
                CALL LL3D_SOLUTIONS(xCOORD,zCOORD,Vu,Vv,DSEx,Ax1,Rpt,DRDXpt,RHOhpt,RHOpt)  
                
                L2U_2D(j)   = Vu
                L2V_2D(j)   = Vv
                L2RHO_3D(j) = RHOhpt     
                L2R_3D(j)   = Rpt           
                                                     
            ENDIF   
                              
        ENDDO       
        
        IF (PROB == 1) THEN
        
            Ut_2D = MATMUL(L2PHI_3D,U_3D(:,i,1))
            
        ELSEIF (PROB == 2) THEN    
            
            Ut_2D = MATMUL(L2PHI_3D,Hu(:,i,1))    
            Vt_2D = MATMUL(L2PHI_3D,Hv(:,i,1))  
            Wt    = MATMUL(L2PHI_3Dw,W_3D(:,i))
            DWt   = MATMUL(L2PHI_3D,DQW(:,i))
            
        ELSEIF (PROB == 3) THEN
        
            Ut_2D = MATMUL(L2PHI_3D,Hu(:,i,1))    
            Vt_2D = MATMUL(L2PHI_3D,Hv(:,i,1))  
            Wt    = MATMUL(L2PHI_3Dw,W_3D(:,i))    
            
        ELSEIF (PROB == 5) THEN
        
            Ut_2D = MATMUL(L2PHI_3D,U_3D(:,i,1))    
            Vt_2D = MATMUL(L2PHI_3D,V_3D(:,i,1))  
            RHOt  = MATMUL(L2PHI_3D,RHO(:,i))  
            Rt    = MATMUL(L2PHI_3D,Rh(:,i))                      
        
        ENDIF       
        
        DO ii = 1, L2_3D
        
            L2DIF2(ii) = (L2U_2D(ii) - Ut_2D(ii))**2.0d0
            IF (PROB == 2) THEN
                L2DIFVvV(ii) = (L2V_2D(ii) - Vt_2D(ii))**2.0d0
                L2DIFw(ii)   = (L2W(ii) - Wt(ii))**2.0d0
                L2DIFdw(ii)  = (L2DW(ii) - DWt(ii))**2.0d0
            ELSEIF (PROB == 3) THEN
                L2DIFVvV(ii) = (L2V_2D(ii) - Vt_2D(ii))**2.0d0
                L2DIFw(ii)   = (L2W(ii) - Wt(ii))**2.0d0         
            ELSEIF (PROB == 5) THEN
                L2DIFVvV(ii) = (L2V_2D(ii) - Vt_2D(ii))**2.0d0
                L2DIFrho(ii) = (L2RHO_3D(ii) - RHOt(ii))**2.0d0  
                L2DIFr(ii)   = (L2R_3D(ii) - Rt(ii))**2.0d0                   
            ENDIF
                   
    	ENDDO
        IF (PROB == 1) THEN
    	    L2_error = L2_error + DOT_PRODUCT(L2WTS_3D,L2DIF2)*(1.0d0/HEX_jac(i,1))
    	ELSEIF (PROB == 2) THEN
    	    L2_error   = L2_error + DOT_PRODUCT(L2WTS_3D,L2DIF2)*(1.0d0/HEX_jac(i,1))
    	    L2Vv_error = L2Vv_error + DOT_PRODUCT(L2WTS_3D,L2DIFVvV)*(1.0d0/HEX_jac(i,1))     	
    	    L2w_error  = L2w_error + DOT_PRODUCT(L2WTS_3D,L2DIFw)*(1.0d0/HEX_jac(i,1))   
    	    L2dw_error = L2dw_error + DOT_PRODUCT(L2WTS_3D,L2DIFdw)*(1.0d0/HEX_jac(i,1)) 
    	ELSEIF (PROB == 3) THEN
    	    L2_error   = L2_error + DOT_PRODUCT(L2WTS_3D,L2DIF2)*(1.0d0/HEX_jac(i,1))
    	    L2Vv_error = L2Vv_error + DOT_PRODUCT(L2WTS_3D,L2DIFVvV)*(1.0d0/HEX_jac(i,1))     	
    	    L2w_error  = L2w_error + DOT_PRODUCT(L2WTS_3D,L2DIFw)*(1.0d0/HEX_jac(i,1))  
    	ELSEIF (PROB == 5) THEN
    	    L2_error    = L2_error    + DOT_PRODUCT(L2WTS_3D,L2DIF2)*(1.0d0/HEX_jac(i,1))
    	    L2Vv_error  = L2Vv_error  + DOT_PRODUCT(L2WTS_3D,L2DIFVvV)*(1.0d0/HEX_jac(i,1))     	
    	    !L2w_error   = L2w_error   + DOT_PRODUCT(L2WTS_3D,L2DIFw)*(1.0d0/HEX_jac(i,1))
    	    RHO_L2error = RHO_L2error + DOT_PRODUCT(L2WTS_3D,L2DIFrho)*(1.0d0/HEX_jac(i,1))   
    	    R_L2error   = R_L2error + DOT_PRODUCT(L2WTS_3D,L2DIFr)*(1.0d0/HEX_jac(i,1))    	       	        
        ENDIF        
        
    ENDDO
    
    DEALLOCATE(L2U_2D)
    DEALLOCATE(L2V_2D)
    DEALLOCATE(L2W)
    DEALLOCATE(L2DW)
    DEALLOCATE(Ut_2D)
    DEALLOCATE(Vt_2D)
    DEALLOCATE(Wt)
    DEALLOCATE(DWt)
    DEALLOCATE(L2DIF2)
    DEALLOCATE(L2DIFVvV)
    DEALLOCATE(L2DIFw)
    DEALLOCATE(L2DIFdw)    
                  
ENDIF

 U_L2error = SQRT(L2_error)
    
 IF (NUM_ERROR == 1) THEN
    write(*,*)'U L2_error = ', U_L2error
    IF (PROB == 2) THEN
        write(*,*)'V L2_error = ', SQRT(L2Vv_error)
        write(*,*)'W L2_error = ', SQRT(L2w_error)
        write(*,*)'DW L2_error =', SQRT(L2dw_error)
    ELSEIF (PROB == 3) THEN
        write(*,*)'V L2_error = ', SQRT(L2Vv_error)
        write(*,*)'W L2_error = ', SQRT(L2w_error)    
    ELSEIF (PROB == 5) THEN
        write(*,*)'V L2_error = ', SQRT(L2Vv_error)
        !write(*,*)'W L2_error = ', SQRT(L2w_error)    
        write(*,*)'RHO L2_error =', SQRT(RHO_L2error)      
        write(*,*)'R L2_error =', SQRT(R_L2error)      
    ENDIF    
 ENDIF

!------------------------------------------------------------------------------------------------------------
!   Surface Gradient
!------------------------------------------------------------------------------------------------------------    

DIF    = 0.000000000000000d0

IF (ELEM == 1) THEN

    DO i = 1, nelemsX
    
        INVJAC = 1/xELEM_jac(i)
    
        DO j = 1, LE

            xpt = X(i) + INVJAC*(1.0d0+L3PTS(j,1))                      !....Transform pts
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real variables to complex
            Huse   = PROB_DATA%Hslope*xpt**2 
            zCOORD = dcmplx(Huse*Z(1), 0.0d0)    
            IF (PROB == 1) THEN
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                 L2GSE(j) = grad_SE
            ELSEIF (PROB == 3) THEN
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                 L2GSE(j) = DSEx   
            ENDIF                   

        ENDDO


       DO j = 1, LE
   
            DIF(j,i)  = SQRT((L2GSE(j) - dZETA(j,i))**2)
        
        ENDDO
    
    ENDDO
    
    GSE_LINFerror = MAXVAL(DIF)
    
ELSEIF (ELEM == 2) THEN

    ALLOCATE(L2GSE2(L2_2D))
    ALLOCATE(DIF2(L2_2D,nelemsX))
    L2GSE2 = 0.0d0
    DIF2   = 0.0d0

    DO i = 1, nelemsX
    
        INVJAC = 1/xELEM_jac(i)
    
        DO j = 1, L2_2D

            xpt = X(i) + INVJAC*(1.0d0+L2PTSz(j,1))                     !....Transform pts
            xCOORD = dcmplx(xpt,0.0d0)                                  !....Convert real
            Huse   = PROB_DATA%Hslope*xpt**2 
            zCOORD = dcmplx(Huse*Z(1), 0.0d0)                                                                     !....variables to complex
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)

            L2GSE2(j) = grad_SE

        ENDDO
        
        DO j = 1, L2_2D
        
            DIF2(j,i) = SQRT((L2GSE2(j)-dZETA(j,i))**2)
            
        ENDDO
        
    ENDDO

    GSE_LINFerror = MAXVAL(DIF2)
    
    DEALLOCATE(L2GSE2)
    DEALLOCATE(DIF2)    
    
ELSEIF (ELEM == 3) THEN

    ALLOCATE(L2GSE2(L2_3D))
    ALLOCATE(DIF2(L2_3D,Qelems))
    L2GSE2 = 0.0d0
    DIF2   = 0.0d0

    DO i = 1, Qelems
    
        globalNODES = CONN(i,:)
        Xnode(1)    = Qnode(globalNODES(1),1)
        Xnode(2)    = Qnode(globalNODES(3),1)
        Ynode(1)    = Qnode(globalNODES(1),2)
        Ynode(2)    = Qnode(globalNODES(2),2)

        DO j = 1, L2_3D
        
            point1 = L2PTS_3D(j,1)
            point2 = L2PTS_3D(j,2)
            xpt    = (1.0d0/2.0d0)*((1.0d0-point1)*Xnode(1)+(1.0d0+point1)*Xnode(2))
            ypt    = (1.0d0/2.0d0)*((1.0d0-point2)*Ynode(1)+(1.0d0+point2)*Ynode(2))
            Huse   = PROB_DATA%Hslope*xpt**2 
            zCOORD = dcmplx(Huse*Z(1), 0.0d0)             
            xCOORD = dcmplx(xpt, 0.0d0)
            yCOORD = dcmplx(ypt, 0.0d0)
            
            IF (PROB == 1) THEN  
                CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                 L2GSE2(j) = grad_SE
            ELSEIF (PROB == 3) THEN
                zpt = Z(nelemsZ)
                CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)  
                L2GSE2(j) = DSEx   
            ELSEIF (PROB == 5) THEN
                xCOORD = dcmplx(xpt,0.0d0)
                zCOORD = dcmplx(Z(1),0.0d0)      
                CALL LL3D_SOLUTIONS(xCOORD,zCOORD,Vu,Vv,DSEx,Ax1,Rpt,DRDXpt,RHOhpt,RHOpt) 
                L2GSE2(j) = DSEx                                 
            ENDIF 
              
        ENDDO
        
        DO j = 1, L2_3D
        
            DIF2(j,i) = SQRT((L2GSE2(j)-dZETA(j,i))**2)

        ENDDO
            
    ENDDO
    
    GSE_LINFerror = MAXVAL(DIF2)
    
    DEALLOCATE(L2GSE2)
    DEALLOCATE(DIF2)
   
ENDIF

IF (NUM_ERROR == 1) THEN
    write(*,*)'dzeta/dx LINF_error = ', GSE_LINFerror
ENDIF

RETURN
END SUBROUTINE COMPUTE_ERROR
!---------------------------------------------------------------------
!
!                  Basis Subroutine
!             Written by Colton J. Conroy
!                   @ the C.H.I.L
!
!--------------------------------------------------------------------
SUBROUTINE DG_BASIS(REGION,REGIONs)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: REGION, REGIONs, jj
REAL(real_p) :: dzeta2, dzetaNEW, zpt
REAL(real_p), DIMENSION(1,3) :: Wpt
REAL(real_p), DIMENSION(1,2) :: Kpt

Wpt = 0.0d0
Kpt = 0.0d0

!IF (w_method == 2) THEN
    DIM = 1
    DO i = 1, L2_3D
        CALL ORTHOGONAL_BASIS(REGIONs, L2PTS_3D(i,1), p+1, DIM, BASIS, DBASIS, BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            PHIr(i,j,1) = BASIS(j)
            PHIr(i,j,2) = DBASIS(j)
        ENDDO    
        CALL ORTHOGONAL_BASIS(REGIONs, L2PTS_3D(i,2), p+1, DIM, BASIS, DBASIS,BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            PHIq(i,j,1) = BASIS(j)
            PHIq(i,j,2) = DBASIS(j)
        ENDDO
        CALL ORTHOGONAL_BASIS(REGIONs, L2PTS_3D(i,3), pu+2, DIM, BASISs, DBASISs,BASIS_2D, DBASIS_2D)
        DO j = 1, pu+2
            PHIs_z(i,j,1) = BASISs(j)
            PHIs_z(i,j,2) = DBASISs(j)
        ENDDO
    ENDDO
    
    zpt = 1.0d0
    CALL ORTHOGONAL_BASIS(REGIONs, zpt, pu+2, DIM, BASISs, DBASISs,BASIS_2D, DBASIS_2D)
    DO j = 1, pu+2
        PHIs_1(j) = BASISs(j)
    ENDDO

!ENDIF        
        
IF (REGION == 1 .OR. REGION == 2) THEN

!....1D basis evaluated at element integration points for continuity eqn

    DIM = 1
    DO i = 1, p+1 
        CALL ORTHOGONAL_BASIS(REGIONs, elemPTS(i,1), p+1, DIM, BASIS, DBASIS,BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            PHI(i,j)  = BASIS(j)
            DPHI(i,j) = DBASIS(j)
        ENDDO 
    ENDDO

!....1D basis evaluated at L2 points

    DO i = 1, L2int
        CALL ORTHOGONAL_BASIS(REGIONs, L2PTS(i,1), p+1, DIM, BASIS, DBASIS,BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            L2PHI(i,j)  = BASIS(j)
            L2DPHI(i,j) = DBASIS(j)
        ENDDO
    ENDDO

!....1D basis evaluated at the boundary points

    bPTS(1) = -1.0d0
    bPTS(2) = 1.0d0
    DO i = 1, NbPTS
        CALL ORTHOGONAL_BASIS(REGIONs, bPTS(i), p+1, DIM, BASIS, DBASIS,BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            bPHI(i,j) = BASIS(j)
            bDPHI(i,j) = DBASIS(j)
        ENDDO
    ENDDO

    IF (debug == 2) THEN
        write(*,*)'L2DPHI=',bDPHI
    ENDIF

ENDIF

!-------------------------------------------------------------------------

IF (REGION == 1) THEN     !....1D basis for momentum eqn

!-------------------------------------------------------------------------
    DO i = 1, pu+1
        CALL ORTHOGONAL_BASIS(REGION, elemPTSz(i,1), pu+1, DIM, BASISz, DBASISz,BASIS_2D, DBASIS_2D)
        DO j = 1, pu+1
            PHIz(i,j)  = BASISz(j)
            DPHIz(i,j) = DBASISz(j)
        ENDDO
    ENDDO

!....1D basis evaluated at L2 points

    DO i = 1, L2int
        CALL ORTHOGONAL_BASIS(REGION, L2PTS(i,1), pu+1, DIM, BASISz, DBASISz,BASIS_2D, DBASIS_2D)
        DO j = 1, pu+1
            L2PHIz(i,j)  = BASISz(j)
            L2DPHIz(i,j) = DBASISz(j)
        ENDDO
    ENDDO

!....1D basis evaluated at the boundary points

    bPTS(1) = -1.0d0
    bPTS(2) = 1.0d0
    DO i = 1, NbPTS
        CALL ORTHOGONAL_BASIS(REGION, bPTS(i), pu+1, DIM, BASISz, DBASISz,BASIS_2D, DBASIS_2D)
        DO j = 1, pu+1
            bPHIz(i,j) = BASISz(j)
            bDPHIz(i,j) = DBASISz(j)
        ENDDO
    ENDDO

!....1D basis evaluated at L2 points

    DO i = 1, LE
        CALL ORTHOGONAL_BASIS(REGION, L3PTS(i,1), p+1, DIM, BASIS, DBASIS,BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            L3PHI(i,j)  = BASIS(j)
            L3DPHI(i,j) = DBASIS(j)
        ENDDO
    ENDDO

    IF (debug == 12) THEN
        write(*,*)'L3DPHI',L3DPHI
    ENDIF
    
!-------------------------------------------------------------------------------------

ELSEIF (REGION == 2) THEN     !....2D basis for momentum equation

!-------------------------------------------------------------------------------------
    DIM = 2

    DO i = 1, NPTS_2D
        CALL ORTHOGONAL_BASIS(REGION, elemPTSz(i,:), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            PHIz(i,j)  = BASIS_2D(j)
            DPHI_2D(i,j,1) = DBASIS_2D(j,1)
            DPHI_2D(i,j,2) = DBASIS_2D(j,2)
        ENDDO
    ENDDO
    
    DO i = 1, L2_2D
        CALL ORTHOGONAL_BASIS(REGION, L2PTSz(i,:), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            L2PHIz(i,j)  = BASIS_2D(j)
            L2DPHI_2D(i,j,1) = DBASIS_2D(j,1)
            L2DPHI_2D(i,j,2) = DBASIS_2D(j,2)
        ENDDO
    ENDDO
    
    IF (coupled == 2) THEN !....If coupling via summation evaluate 2D shape functions along lines below int. points
        DO ii = 1, LE
            DO i = 1, L2int
                CALL ORTHOGONAL_BASIS(REGION, UintPTS(i,:,ii), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
                DO j = 1, NDOF
                    PHIint(i,j,ii)  = BASIS_2D(j)
                ENDDO
            ENDDO
        ENDDO
        DO ii = 1, L2_2D
            DO i = 1, L2int
                CALL ORTHOGONAL_BASIS(REGION, VintPTS(i,:,ii), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
                DO j = 1, NDOF
                    VPHIint(i,j,ii)  = BASIS_2D(j)
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    
    DO i = 1, pu+1
        CALL ORTHOGONAL_BASIS(REGION, bPTSz(i,:,1), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2D(i,j,1)  = BASIS_2D(j)                   !....first page is bottom boundary
            bDPHI_2D(i,j,1) = DBASIS_2D(j,1)                !....derivative w/respect to x
            bDPHI_2D(i,j,2) = DBASIS_2D(j,2)                !....derivative w/respect to y
        ENDDO 
        CALL ORTHOGONAL_BASIS(REGION, bPTSz(i,:,2), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2D(i,j,2)  = BASIS_2D(j)                   !....second page is top boundary
            bDPHI_2D(i,j,3) = DBASIS_2D(j,1)                !....derivative w/respect to x
            bDPHI_2D(i,j,4) = DBASIS_2D(j,2)                !....derivative w/respect to y
        ENDDO 
    ENDDO
    
    DO i = 1, pu+1
        CALL ORTHOGONAL_BASIS(REGION, bPTSx(i,:,1), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2Dx(i,j,1) = BASIS_2D(j)                  !....Left boundary
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTSx(i,:,2), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)   
        DO j = 1, NDOF
            bPHI_2Dx(i,j,2) = BASIS_2D(j)                  !....Right boundary
        ENDDO         
    ENDDO     
        
    !....New Points
    
    DO i = 1, NPTS_2D
        dzeta2 = 1.0d0 + elemPTSz(i,2)
        DO j = 1, L2int
            dzetaNEW      = dzeta2*(L2PTS(j,1) + 1.0d0)/2.0d0
            L2ptsNEW(j,2,i) = -1.0d0 + dzetaNEW
            L2ptsNEW(j,1,i) = elemPTSz(i,1)
        ENDDO
    ENDDO
    

    !....Horizontal basis evaluated at 2D x points for evaluating the gradient of the surface elevation
    
    DIM = 1
    DO i = 1, NPTS_2D
        CALL ORTHOGONAL_BASIS(REGIONs, elemPTSz(i,1), p+1, DIM, BASIS, DBASIS, BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            PHIzeta(i,j) = BASIS(j)
        ENDDO
    ENDDO
    
    DIM = 1
    DO i = 1, L2_2D
        CALL ORTHOGONAL_BASIS(REGIONs, L2PTSz(i,1), p+1, DIM, BASIS, DBASIS, BASIS_2D, DBASIS_2D)
        DO j = 1, p+1
            PHIzeta2(i,j) = BASIS(j)
            DPHIzeta(i,j) = DBASIS(j)
        ENDDO
    ENDDO

    IF (coupled == 2) THEN
        DO i = 1, LE
            CALL ORTHOGONAL_BASIS(REGIONs, L3PTS(i,1), p+1, DIM, BASIS, DBASIS,BASIS_2D, DBASIS_2D)
            DO j = 1, p+1
                L3PHI(i,j)  = BASIS(j)
                L3DPHI(i,j) = DBASIS(j)
            ENDDO
        ENDDO  
    ENDIF     
!-------------------------------------------------------------------------------------------------------------          
ELSEIF (REGION == 3) THEN
!-------------------------------------------------------------------------------------------------------------
    
    DIM = 2
    REGION = 2
    DO i = 1, NPTS_2D
        CALL ORTHOGONAL_BASIS(REGION, elemPTSz(i,:), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            PHIz(i,j)  = BASIS_2D(j)
            DPHI_2D(i,j,1) = DBASIS_2D(j,1)
            DPHI_2D(i,j,2) = DBASIS_2D(j,2)
        ENDDO
    ENDDO
   
    DO i = 1, L2_2D
        CALL ORTHOGONAL_BASIS(REGION, L2PTSz(i,:), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            L2PHIz(i,j)  = BASIS_2D(j)
            L2DPHI_2D(i,j,1) = DBASIS_2D(j,1)
            L2DPHI_2D(i,j,2) = DBASIS_2D(j,2)
        ENDDO
    ENDDO

!    IF (coupled == 2) THEN !....If coupling via summation evaluate 2D shape functions along lines below int. points
!        DO ii = 1, LE
!            DO i = 1, L2int
!                CALL ORTHOGONAL_BASIS(REGION, UintPTS(i,:,ii), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
!                DO j = 1, NDOF
!                    PHIint(i,j,ii)  = BASIS_2D(j)
!                ENDDO
!            ENDDO
!        ENDDO
!        DO ii = 1, L2_2D
!            DO i = 1, L2int
!                CALL ORTHOGONAL_BASIS(REGION, VintPTS(i,:,ii), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
!                DO j = 1, NDOF
!                    VPHIint(i,j,ii)  = BASIS_2D(j)
!                ENDDO
!            ENDDO
!        ENDDO
!    ENDIF
    
    DO i = 1, p+1
        CALL ORTHOGONAL_BASIS(REGION, bPTSz(i,:,1), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2D(i,j,1)  = BASIS_2D(j)                   !....first page is bottom boundary
            bDPHI_2D(i,j,1) = DBASIS_2D(j,1)                !....derivative w/respect to x
            bDPHI_2D(i,j,2) = DBASIS_2D(j,2)                !....derivative w/respect to y
        ENDDO 
        CALL ORTHOGONAL_BASIS(REGION, bPTSz(i,:,2), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2D(i,j,2)  = BASIS_2D(j)                   !....second page is top boundary
            bDPHI_2D(i,j,3) = DBASIS_2D(j,1)                !....derivative w/respect to x
            bDPHI_2D(i,j,4) = DBASIS_2D(j,2)                !....derivative w/respect to y
        ENDDO 
    ENDDO
    
    DO i = 1, p+1
        CALL ORTHOGONAL_BASIS(REGION, bPTSx(i,:,1), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2Dx(i,j,1) = BASIS_2D(j)                  !....Left boundary
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTSx(i,:,2), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)   
        DO j = 1, NDOF
            bPHI_2Dx(i,j,2) = BASIS_2D(j)                  !....Right boundary
        ENDDO         
    ENDDO     
    
    DO i = 1, p+1
        CALL ORTHOGONAL_BASIS(REGION, bPTSy(i,:,1), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2Dy(i,j,1) = BASIS_2D(j)                  !....Left boundary
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTSy(i,:,2), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)   
        DO j = 1, NDOF
            bPHI_2Dy(i,j,2) = BASIS_2D(j)                  !....Right boundary
        ENDDO         
    ENDDO 
    
    !....Tensor Product

    DO i = 1, LRK
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,1), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw)
        DO j = 1, NDOF_2D
            bPHI_2DTP(i,j,1) = BASIS_2Dw(j)                  !....left x-face
        ENDDO
        
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,2), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw)   
        DO j = 1, NDOF_2D
            bPHI_2DTP(i,j,2) = BASIS_2Dw(j)                  !....right x-face
        ENDDO    
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,3), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw)   
        DO j = 1, NDOF_2D
            bPHI_2DTP(i,j,3) = BASIS_2Dw(j)                  !....bottom y-face
        ENDDO               
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,4), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw)   
        DO j = 1, NDOF_2D
            bPHI_2DTP(i,j,4) = BASIS_2Dw(j)                  !....top y-face
        ENDDO          
    ENDDO      
    
    DO i = 1, LRK
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,1), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            bPHI_2DTPh(i,j,1) = BASIS_2D(j)                  !....left x-face
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,2), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)   
        DO j = 1, NDOF
            bPHI_2DTPh(i,j,2) = BASIS_2D(j)                  !....right x-face
        ENDDO    
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,3), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)   
        DO j = 1, NDOF
            bPHI_2DTPh(i,j,3) = BASIS_2D(j)                  !....bottom y-face
        ENDDO               
        CALL ORTHOGONAL_BASIS(REGION, bPTS_2DTP(i,:,4), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)   
        DO j = 1, NDOF
            bPHI_2DTPh(i,j,4) = BASIS_2D(j)                  !....top y-face
        ENDDO          
    ENDDO    
    
    !....2D basis evaluated at face (x,y) points
    
    DO i = 1, NPTS_2D
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,1), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            b_PHI_2D3Dx(i,j,1) = BASIS_2D(j)
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,2), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            b_PHI_2D3Dx(i,j,2) = BASIS_2D(j)
        ENDDO
    ENDDO
    
    DO i = 1, NPTS_2D
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,5), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            b_PHI_2D3Dy(i,j,1) = BASIS_2D(j)
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,1:2,6), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            b_PHI_2D3Dy(i,j,2) = BASIS_2D(j)
        ENDDO
    ENDDO
    
    !....2D basis evaluated at 3D element pts
    
    DO i = 1, NPTS_3D
        CALL ORTHOGONAL_BASIS(REGION, elemPTS_3D(i,1:2), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            PHI_2D3D(i,j) = BASIS_2D(j)
        ENDDO
    ENDDO    
    
    !....Vertical Velocity
    
    DO i = 1, TP_2D
        CALL ORTHOGONAL_BASIS(REGION, TPTS_2D(i,:), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw)
        DO j = 1, NDOF_2D
            PHI_TP2D(i,j)    = BASIS_2Dw(j)
            DPHI_TP2D(i,j,1) = DBASIS_2Dw(j,1)
            DPHI_TP2D(i,j,2) = DBASIS_2Dw(j,2)
        ENDDO
    ENDDO
    
    DO i = 1, TP_2D
        CALL ORTHOGONAL_BASIS(REGION, TPTS_2D(i,:), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            PHI_TP2Dh(i,j) = BASIS_2D(j)
        ENDDO
    ENDDO     
        
    !....New Points
    
!    DO i = 1, NPTS_2D
!        dzeta2 = 1.0d0 + elemPTSz(i,2)
!        DO j = 1, L2int
!            dzetaNEW      = dzeta2*(L2PTS(j,1) + 1.0d0)/2.0d0
!            L2ptsNEW(j,2,i) = -1.0d0 + dzetaNEW
!            L2ptsNEW(j,1,i) = elemPTSz(i,1)
!        ENDDO
!    ENDDO
    

    !....Horizontal basis evaluated at 2D x points for evaluating the gradient of the surface elevation
    
!    DO i = 1, NPTS_2D
!        CALL ORTHOGONAL_BASIS(REGION, elemPTS_3D(i,1), p+1, DIM, BASIS, DBASIS, BASIS_2D, DBASIS_2D)
!        DO j = 1, p+1
!            PHIzeta(i,j) = BASIS(j)
!        ENDDO
!    ENDDO
!    
!    DIM = 1
!    DO i = 1, L2_2D
!        CALL ORTHOGONAL_BASIS(REGIONs, L2PTSz(i,1), p+1, DIM, BASIS, DBASIS, BASIS_2D, DBASIS_2D)
!        DO j = 1, p+1
!            PHIzeta2(i,j) = BASIS(j)
!            DPHIzeta(i,j) = DBASIS(j)
!        ENDDO
!    ENDDO
!
!    IF (coupled == 2) THEN
!        DO i = 1, LE
!            CALL ORTHOGONAL_BASIS(REGIONs, L3PTS(i,1), p+1, DIM, BASIS, DBASIS,BASIS_2D, DBASIS_2D)
!            DO j = 1, p+1
!                L3PHI(i,j)  = BASIS(j)
!                L3DPHI(i,j) = DBASIS(j)
!            ENDDO
!        ENDDO  
!    ENDIF

    DIM = 3
    REGION = 3
    DO i = 1, NPTS_3D
        CALL ORTHOGONAL_BASIS(REGION, elemPTS_3D(i,:), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)
        DO j = 1, NDOF_3D
            PHI_3D(i,j)  = BASIS_3D(j)
            DPHI_3D(i,j,1) = DBASIS_3D(j,1)
            DPHI_3D(i,j,2) = DBASIS_3D(j,2)
            DPHI_3D(i,j,3) = DBASIS_3D(j,3)
        ENDDO
    ENDDO
    !PHI_3D(:,1) = 1.0d0
    OPEN(10,FILE='PHI_3D.2d')
    DO i = 1, NPTS_3D
        WRITE(10,*)PHI_3D(i,:)  
        write(10,*)'i=',i
        write(10,*)'********************'
    ENDDO
    CLOSE(10)
    
    DO i = 1, L2_3D
        CALL ORTHOGONAL_BASIS(REGION, L2PTS_3D(i,:), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)
        DO j = 1, NDOF_3D
            L2PHI_3D(i,j)  = BASIS_3D(j)
            L2DPHI_3D(i,j,1) = DBASIS_3D(j,1)
            L2DPHI_3D(i,j,2) = DBASIS_3D(j,2)
            L2DPHI_3D(i,j,3) = DBASIS_3D(j,3)
        ENDDO
    ENDDO
    
    DO i = 1, L2_3D
        CALL ORTHOGONAL_BASIS(REGION, L2PTS_3D(i,:), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw)
        DO j = 1, NDOF_3Dw
            L2PHI_3Dw(i,j)  = BASIS_3Dw(j)
        ENDDO
    ENDDO     
    
    !....Vertical Velocity L3 projection
    
    DO i = 1, TP_3D
        CALL ORTHOGONAL_BASIS(REGION, TPTS_3D(i,:), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw)
        DO j = 1, NDOF_3Dw
            L3PHIw(i,j) = BASIS_3Dw(j)
        ENDDO
    ENDDO    
       
    !....Boundary Faces optimal

    DO i = 1, NPTS_2D
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,1), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D )
        DO j = 1, NDOF_3D
            bPHI_3D(i,j,1) = BASIS_3D(j)                  !....left x-face
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,2), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D )   
        DO j = 1, NDOF_3D
            bPHI_3D(i,j,2) = BASIS_3D(j)                  !....right x-face
        ENDDO    
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,3), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D )   
        DO j = 1, NDOF_3D
            bPHI_3D(i,j,3) = BASIS_3D(j)                  !....bottom z-face
        ENDDO               
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,4), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D )   
        DO j = 1, NDOF_3D
            bPHI_3D(i,j,4) = BASIS_3D(j)                  !....top z-face
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,5), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D )   
        DO j = 1, NDOF_3D
            bPHI_3D(i,j,5) = BASIS_3D(j)                  !....front y-face
        ENDDO                       
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,6), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D )   
        DO j = 1, NDOF_3D
            bPHI_3D(i,j,6) = BASIS_3D(j)                  !....back y-face
        ENDDO           
    ENDDO  
    
    !....Sigma Faces for W
    
    DO i = 1, NPTS_2D
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,3), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )   
        DO j = 1, NDOF_3Dw
            bPHI_3Dw(i,j,1) = BASIS_3Dw(j)                  !....bottom z-face
        ENDDO               
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D(i,:,4), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )   
        DO j = 1, NDOF_3Dw
            bPHI_3Dw(i,j,2) = BASIS_3Dw(j)                  !....top z-face
        ENDDO
    ENDDO        
    
    !....Elements for W
    
    DO i = 1, NPTS_3D
        CALL ORTHOGONAL_BASIS(REGION, elemPTS_3D(i,:), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw)
        DO j = 1, NDOF_3Dw
            PHI_3Dw(i,j)  = BASIS_3Dw(j)
        ENDDO
    ENDDO    
    
    !....Tensor Product
    
    DO i = 1, TP_2D             
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3DTP(i,:,1), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )
        DO j = 1, NDOF_3D
            bPHI_3DTP(i,j,1) = BASIS_3Dw(j)                  !....left x-face
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3DTP(i,:,2), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )   
        DO j = 1, NDOF_3D
            bPHI_3DTP(i,j,2) = BASIS_3Dw(j)                  !....right x-face
        ENDDO    
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3DTP(i,:,3), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )   
        DO j = 1, NDOF_3D
            bPHI_3DTP(i,j,3) = BASIS_3Dw(j)                  !....bottom z-face
        ENDDO               
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3DTP(i,:,4), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )   
        DO j = 1, NDOF_3D
            bPHI_3DTP(i,j,4) = BASIS_3Dw(j)                  !....top z-face
        ENDDO
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3DTP(i,:,5), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )   
        DO j = 1, NDOF_3D
            bPHI_3DTP(i,j,5) = BASIS_3Dw(j)                  !....front y-face
        ENDDO                       
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3DTP(i,:,6), pw+1, DIM, BASISz, DBASISz, BASIS_2Dw, DBASIS_2Dw, BASIS_3Dw, DBASIS_3Dw )   
        DO j = 1, NDOF_3D
            bPHI_3DTP(i,j,6) = BASIS_3Dw(j)                  !....back y-face
        ENDDO           
    ENDDO    
    
    !....Vertical Velocity RHS
    
    DO jj = 1, nrkW
    
        Wpt(1,3) = 2.0d0*ctW(jj) - 1.0d0
        
        DO ii = 1, TP_2D
            Wpt(1,1) = TPTS_2D(ii,1)
            Wpt(1,2) = TPTS_2D(ii,2)
            CALL ORTHOGONAL_BASIS(REGION, Wpt, pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)
            DO j = 1, NDOF_3D
                PHIw(ii,j,jj)  = BASIS_3D(j)
            ENDDO
        ENDDO      

        DO i = 1, LRK
            Wpt(1,1) = bPTS_2DTP(i,1,1)
            Wpt(1,2) = bPTS_2DTP(i,2,1)
            CALL ORTHOGONAL_BASIS(REGION, Wpt, pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)
            DO j = 1, NDOF_3D
                bPHIw(i,j,jj,1) = BASIS_3D(j)                  !....left x-face
            ENDDO
            Wpt(1,1) = bPTS_2DTP(i,1,2)
            Wpt(1,2) = bPTS_2DTP(i,2,2)    
            CALL ORTHOGONAL_BASIS(REGION, Wpt, pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)   
            DO j = 1, NDOF_3D
                bPHIw(i,j,jj,2) = BASIS_3D(j)                  !....right x-face
            ENDDO    
            Wpt(1,1) = bPTS_2DTP(i,1,3)
            Wpt(1,2) = bPTS_2DTP(i,2,3)    
            CALL ORTHOGONAL_BASIS(REGION, Wpt, pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)   
            DO j = 1, NDOF_3D
                bPHIw(i,j,jj,3) = BASIS_3D(j)                  !....bottom y-face
            ENDDO     
            Wpt(1,1) = bPTS_2DTP(i,1,4)
            Wpt(1,2) = bPTS_2DTP(i,2,4)              
            CALL ORTHOGONAL_BASIS(REGION, Wpt, pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)   
            DO j = 1, NDOF_3D
                bPHIw(i,j,jj,4) = BASIS_3D(j)                  !....top y-face
            ENDDO          
        ENDDO
        
    ENDDO         

    !....Evaluate U at bottom boundary

    DO i = 1, L2_2D
        CALL ORTHOGONAL_BASIS(REGION, bPTS_3D2D(i,:), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D )
        DO j = 1, NDOF_3D
            bPHI_3D2D(i,j) = BASIS_3D(j)
        ENDDO
    ENDDO        
    
    IF (coupled == 2) THEN
    
        DO ii = 1, L2_2D
            DO i = 1, L2int
                CALL ORTHOGONAL_BASIS(REGION, UintPTS_3D(i,:,ii), pu+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D, BASIS_3D, DBASIS_3D)
                DO j = 1, NDOF_3D
                    PHIint_3D(i,j,ii) = BASIS_3D(j)
                ENDDO
            ENDDO
        ENDDO
        
    ENDIF    
    
    DIM = 2
    REGION = 2
    DO i = 1, L2_3D
        CALL ORTHOGONAL_BASIS(REGION, L2PTS_3D(i,1:2), p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
        DO j = 1, NDOF
            PHIzeta_3D(i,j)    = BASIS_2D(j)
            DPHIzeta_3D(i,j,1) = DBASIS_2D(j,1)     !....Derivative of zeta wrt x @ 3D L2 pts
            DPHIzeta_3D(i,j,2) = DBASIS_2D(j,2)     !....Derivative of zeta wrt y @ 3D L2 pts
        ENDDO
    ENDDO    
    
    IF (KERNEL == 1) THEN
        ALLOCATE(PHIkernel((NPTS_K)**2,NDOF))
        ii = 0
        DO i = 1, NPTS_K
            DO jj = 1, NPTS_K
                ii = ii + 1
                Kpt(1,1) = KPTS_1D(jj,1)
                Kpt(1,2) = KPTS_1D(i,1)

                CALL ORTHOGONAL_BASIS(REGION, Kpt, p+1, DIM, BASISz, DBASISz, BASIS_2D, DBASIS_2D)
                DO j = 1, NDOF
                    PHIkernel(ii,j) = BASIS_2D(j)
                ENDDO
            ENDDO
        ENDDO        
    ENDIF    
        
             
ENDIF

RETURN
END SUBROUTINE DG_BASIS
!--------------------------------------------------------------------
!
!                   Quadrature Subroutine
!                 Written by Colton J. Conroy
!                       @ the C.H.I.L
!
!--------------------------------------------------------------------
SUBROUTINE DG_QUADRATURE(REGION,REGIONs,DEG,DEG_3D)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: REGION, REGIONs, DEG, DEG_3D, switch


!....Quadrature points and weights 1D (Continuity Eqn)
DIM = 1
CALL QUADRATURE(p+1,REGIONs,elemPTS,elemWTS,p+1,DIM)             !....Horizontal element integration pts
CALL QUADRATURE(L2int,REGIONs,L2PTS,L2WTS,L2int,DIM)             !....L2 integration pts
CALL QUADRATURE(pu+1,REGIONs,ELEM_pts1D,ELEM_wts1D,pu+1,DIM)
CALL QUADRATURE(L2TP,REGIONs,L2_pts1D,L2_wts1D,L2TP,DIM)
CALL QUADRATURE(LRK,REGIONs,TPTS_1D,TWTS_1D,LRK,DIM)  
IF (KERNEL == 1) THEN
    CALL GR_QUADRATURE(NPTS_K-1,KPTS_1D,KWTS_1D,NPTS_K)
ENDIF
    
IF (REGION == 1 .OR. REGION == 2) THEN    
    IF (coupled == 2) THEN                                           !....Int pts for depth integrated vel. for 
        DO i = 1, LE                                                 !....coupling #2
            DO j = 1, L2int
                UintPTS(j,1,i) = elemPTS(i,1)
                UintPTS(j,2,i) = L2PTS(j,1)
            ENDDO
        ENDDO
        CALL QUADRATURE(LE,REGIONs,L3PTS,L3WTS,L2int,DIM)  
    ENDIF
ENDIF

IF (REGION == 1) THEN

    CALL QUADRATURE(pu+1,REGION,elemPTSz,elemWTSz,pu+1,DIM)          !....Vertical element integration pts
    CALL QUADRATURE(LE,REGION,L3PTS,L3WTS,L2int,DIM)                 !....L3 integration pts
    
ELSEIF (REGION == 2 .OR. REGION == 3) THEN

    DIM = 2
    IF (REGION == 2) THEN
        CALL QUADRATURE(pu+1,REGIONs,bPTSz,bWTSz,pu+1,DIM)       !....Boundary pts
    ELSEIF (REGION == 3) THEN
        CALL QUADRATURE(p+1,REGIONs,bPTSz,bWTSz,p+1,DIM)  
    ENDIF       

    bPTSz(:,2,1) = -1.0d0                                        !....Bottom y coordinate
    bPTSz(:,:,2) = bPTSz(:,:,1)                         
    bPTSz(:,2,2) = 1.0d0                                         !....Top y coordinate
    
    bPTSx(:,2,1) = bPTSz(:,1,1)                                  !....X Boundary pts   
    bPTSx(:,1,1) = -1.0d0                                        !....Left boundary
    bPTSx(:,2,2) = bPTSx(:,2,1)
    bPTSx(:,1,2) = 1.0d0                                         !....Right boundary  
    
    bPTSy(:,1,1) = elemPTS(:,1)
    bPTSy(:,2,1) = -1.0d0
    bPTSy(:,1,2) = elemPTS(:,1)
    bPTSy(:,2,2) = 1.0d0
      
    IF (coupled == 2) THEN
        DO i = 1, L2_2D                                            
            DO j = 1, L2int
                VintPTS(j,1,i) = L2PTSz(i,1)
                VintPTS(j,2,i) = L2PTS(j,1)
            ENDDO
        ENDDO     
    ENDIF
    
    IF (REGION == 3) THEN
        switch = 1
    ELSE
        switch = 0
    ENDIF
    
    REGION = 2
    CALL QUADRATURE(DEG,REGION,elemPTSz,elemWTSz,NPTS_2D,DIM)        !....Element points
    CALL QUADRATURE(DEG*2,REGION,L2PTSz,L2WTSz,L2_2D,DIM)            !....L2 pts

    IF (switch == 1) THEN
        REGION = 3
    ENDIF
    
  
ENDIF

IF (REGION == 3) THEN

    !....Left x-face
    
    bPTS_3D(:,1,1) = -1.000000000000000d0
    bPTS_3D(:,2,1) = elemPTSz(:,1)
    bPTS_3D(:,3,1) = elemPTSz(:,2)
    
    !....Rigth x-face
 
    bPTS_3D(:,1,2) = 1.000000000000000d0
    bPTS_3D(:,2,2) = elemPTSz(:,1)
    bPTS_3D(:,3,2) = elemPTSz(:,2)   
    
    !....Bottom z-face
    
    bPTS_3D(:,1,3) = elemPTSz(:,1)
    bPTS_3D(:,2,3) = elemPTSz(:,2)
    bPTS_3D(:,3,3) = -1.0000000000000000d0
    
    !....Top z-face
    
    bPTS_3D(:,1,4) = elemPTSz(:,1)
    bPTS_3D(:,2,4) = elemPTSz(:,2)
    bPTS_3D(:,3,4) = 1.0000000000000000d0    
    
    !....Front y-face
    
    bPTS_3D(:,1,5) = elemPTSz(:,1)
    bPTS_3D(:,2,5) = -1.0000000000000000d0
    bPTS_3D(:,3,5) = elemPTSz(:,2)
    
    !....Back y-face
    
    bPTS_3D(:,1,6) = elemPTSz(:,1)
    bPTS_3D(:,2,6) = 1.0000000000000000d0
    bPTS_3D(:,3,6) = elemPTSz(:,2)    
    
    !....Bottom boundary pts for evaluation of bottom velocity
    
    bPTS_3D2D(:,1) = L2PTSz(:,1)
    bPTS_3D2D(:,2) = L2PTSz(:,2)
    bPTS_3D2D(:,3) = -1.000000000000000d0
    
    !....Volume points (element pts. and L2 pts.)
    
    DIM = 3
    IF (QUAD_TYPE == 1) THEN    
        CALL QUADRATURE(DEG_3D,REGION,elemPTS_3D,elemWTS_3D,NPTS_3D,DIM)
        IF (pu < 2) THEN
            CALL QUADRATURE(L2_DEG_3D,REGION,L2PTS_3D,L2WTS_3D,L2_3D,DIM)
        ELSEIF (pu >= 2)THEN
            CALL QUADRATURE(L2_DEG_3D,REGION,L2PTS_3D,L2WTS_3D,L2_3D,DIM)
        ENDIF
    ELSE
        CALL TENSOR_QUAD_3D(elemPTS_3D,elemWTS_3D,ELEM_pts1D,ELEM_wts1D,pu+1)
        CALL TENSOR_QUAD_3D(L2PTS_3D,L2WTS_3D,L2_pts1D,L2_wts1D,L2TP)
    ENDIF
    
    !....Points for vertical velocity w
    
    CALL TENSOR_QUAD_2D(TPTS_2D,TWTS_2D,TPTS_1D,TWTS_1D,LRK)
    CALL TENSOR_QUAD_3D(TPTS_3D,TWTS_3D,TPTS_1D,TWTS_1D,LRK)
    
    IF (KERNEL == 1) THEN
        CALL TENSOR_QUAD_2D(KPTS_2D,KWTS_2D,KPTS_1D,KWTS_1D,NPTS_K)
    ENDIF    
    
   !....2D Boundary points
    
   !....Left x-edge
    
    bPTS_2DTP(:,1,1) = -1.000000000000000d0
    bPTS_2DTP(:,2,1) = TPTS_1D(:,1)
    
    !....Right x-edge
 
    bPTS_2DTP(:,1,2) = 1.000000000000000d0
    bPTS_2DTP(:,2,2) = TPTS_1D(:,1)   
    
    !....Bottom y-edge
    
    bPTS_2DTP(:,1,3) = TPTS_1D(:,1)
    bPTS_2DTP(:,2,3) = -1.0000000000000000d0
    
    !....Top y-edge
    
    bPTS_2DTP(:,1,4) = TPTS_1D(:,1)
    bPTS_2DTP(:,2,4) = 1.0000000000000000d0  
    
    !....3D Tensor Product
    
    !....Left x-face
    
    bPTS_3DTP(:,1,1) = -1.000000000000000d0
    bPTS_3DTP(:,2,1) = TPTS_2D(:,1)
    bPTS_3DTP(:,3,1) = TPTS_2D(:,2)
    
    !....Right x-face
 
    bPTS_3DTP(:,1,2) = 1.000000000000000d0
    bPTS_3DTP(:,2,2) = TPTS_2D(:,1)
    bPTS_3DTP(:,3,2) = TPTS_2D(:,2)   
    
    !....Bottom z-face
    
    bPTS_3DTP(:,1,3) = TPTS_2D(:,1)
    bPTS_3DTP(:,2,3) = TPTS_2D(:,2)
    bPTS_3DTP(:,3,3) = -1.0000000000000000d0
    
    !....Top z-face
    
    bPTS_3DTP(:,1,4) = TPTS_2D(:,1)
    bPTS_3DTP(:,2,4) = TPTS_2D(:,2)
    bPTS_3DTP(:,3,4) = 1.0000000000000000d0    
    
    !....Front y-face
    
    bPTS_3DTP(:,1,5) = TPTS_2D(:,1)
    bPTS_3DTP(:,2,5) = -1.0000000000000000d0
    bPTS_3DTP(:,3,5) = TPTS_2D(:,2)
    
    !....Back y-face
    
    bPTS_3DTP(:,1,6) = TPTS_2D(:,1)
    bPTS_3DTP(:,2,6) = 1.0000000000000000d0
    bPTS_3DTP(:,3,6) = TPTS_2D(:,2)    
            
    IF (coupled == 2) THEN
      
        DO i = 1, L2_2D
            DO j = 1, L2int
                UintPTS_3D(j,1,i) = L2PTSz(i,1)
                UintPTS_3D(j,2,i) = L2PTSz(i,2)
                UintPTS_3D(j,3,i) = L2PTS(j,1)
            ENDDO
        ENDDO
        
    ENDIF        
    
ENDIF

RETURN
END SUBROUTINE DG_QUADRATURE
!--------------------------------------------------------------------
!
!                  Allocate Matricies Subroutine
!                   Written by Colton J. Conroy
!                          @ the C.H.I.L
!
!--------------------------------------------------------------------
SUBROUTINE MATRIX_ALLOCATION(REGION)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: REGION, ierr

!....Continuity eqn

ALLOCATE(elemPTS(p+1,1))
ALLOCATE(elemWTS(p+1))
ALLOCATE(L2PTS(L2int,1))
ALLOCATE(L2WTS(L2int))
ALLOCATE(L3PTS(LE,1))
ALLOCATE(L3WTS(LE))
ALLOCATE(bPTS(NbPTS))
ALLOCATE(L2_pts1D(L2TP,1))
ALLOCATE(L2_wts1D(L2TP))
ALLOCATE(ELEM_pts1D(pu+1,1))
ALLOCATE(ELEM_wts1D(pu+1))
ALLOCATE(PHI(p+1,p+1))
ALLOCATE(DPHI(p+1,p+1))
ALLOCATE(L2PHI(L2int,p+1))
ALLOCATE(L2DPHI(L2int,p+1))
ALLOCATE(L3PHI(LE,p+1))
ALLOCATE(L3DPHI(LE,p+1))
ALLOCATE(bPHI(NbPTS,p+1))
ALLOCATE(bDPHI(NbPTS,p+1))
ALLOCATE(BASIS(p+1))
ALLOCATE(DBASIS(p+1))
ALLOCATE(L2(p+1,L2int))
ALLOCATE(A(p+1,p+1))
ALLOCATE(C(p+1,L2int))
ALLOCATE(B1(p+1))
ALLOCATE(B2(p+1))
ALLOCATE(L3(p+1,LE))
ALLOCATE(TPTS_1D(LRK,1))
ALLOCATE(TWTS_1D(LRK))
ALLOCATE(TPTS_2D(TP_2D,2))
ALLOCATE(TWTS_2D(TP_2D))
ALLOCATE(TPTS_3D(TP_3D,3))
ALLOCATE(TWTS_3D(TP_3D))
ALLOCATE(L3w(NDOF_3Dw,TP_3D))
ALLOCATE(L3PHIw(TP_3D,NDOF_3Dw))
ALLOCATE(bPTS_2DTP(LRK,2,4))
ALLOCATE(bPTS_3DTP(TP_2D,3,6))
ALLOCATE(PHI_TP2D(TP_2D,NDOF_2D))
ALLOCATE(DPHI_TP2D(TP_2D,NDOF_2D,2))
ALLOCATE(bPHI_3DTP(TP_2D,NDOF_3D,6))
ALLOCATE(bPHI_2DTP(LRK,NDOF_2D,4))
ALLOCATE(PHI_TP2Dh(TP_2D,NDOF))
ALLOCATE(bPHI_2DTPh(LRK,NDOF,4))
!IF (w_method == 2) THEN
    ALLOCATE(PHIr(L2_3D,p+1,2))
    ALLOCATE(PHIs_z(L2_3D,pu+2,2))
    ALLOCATE(PHIq(L2_3D,p+1,2))
    ALLOCATE(PHIs_1(pu+2))
    ALLOCATE(BASISs(pu+2))
    ALLOCATE(DBASISs(pu+2))
    PHIr    = 0.0d0
    PHIs_z  = 0.0d0
    PHIq    = 0.0d0
    PHIs_1  = 0.0d0
    BASISs  = 0.0d0
    DBASISs = 0.0d0
!ENDIF
IF (REGION == 3) THEN
    ALLOCATE(elemPTSz(NPTS_2D,2))
    ALLOCATE(elemWTSz(NPTS_2D))
    ALLOCATE(L2PTSz(L2_2D,2))
    ALLOCATE(L2WTSz(L2_2D))
    ALLOCATE(bPTSz(p+1,2,2))
    ALLOCATE(bPTSx(p+1,2,2))
    ALLOCATE(bPTSy(p+1,2,2))
    ALLOCATE(bWTSz(p+1))
    ALLOCATE(PHIz(NPTS_2D,NDOF))
    ALLOCATE(DPHI_2D(NPTS_2D,NDOF,2))
    ALLOCATE(L2PHIz(L2_2D,NDOF))
    ALLOCATE(L2DPHI_2D(L2_2D,NDOF,2))
    ALLOCATE(bPHI_2D(p+1,NDOF,2))
    ALLOCATE(bPHI_2Dx(p+1,NDOF,2))    
    ALLOCATE(bPHI_2Dy(p+1,NDOF,2))
    ALLOCATE(bDPHI_2D(p+1,NDOF,4))
    ALLOCATE(PHIzeta(NPTS_2D,p+1))
    ALLOCATE(PHIzeta2(L2_2D,p+1))
    ALLOCATE(DPHIzeta(L2_2D,p+1))
    ALLOCATE(BASIS_2D(NDOF))
    ALLOCATE(DBASIS_2D(NDOF,2))
    ALLOCATE(BASIS_2Dw(NDOF_2D))
    ALLOCATE(DBASIS_2Dw(NDOF_2D,2))     
    ALLOCATE(L2U(NDOF,L2_2D))
    ALLOCATE(A_2D(NDOF,NPTS_2D,2))
    ALLOCATE(CU(NDOF,L2_2D))
    ALLOCATE(BU(NDOF,p+1,2))
    ALLOCATE(BX(NDOF,p+1,2))
    ALLOCATE(L2ptsNEW(L2int,2,NPTS_2D))
ENDIF

!....Momentum Eqn
IF (REGION == 1) THEN
    ALLOCATE(elemPTSz(pu+1,2))
    ALLOCATE(elemWTSz(pu+1))
    ALLOCATE(PHIz(pu+1,pu+1))
    ALLOCATE(DPHIz(pu+1,pu+1))
    ALLOCATE(L2PHIz(L2int,pu+1))
    ALLOCATE(L2DPHIz(L2int,pu+1))
    ALLOCATE(bPHIz(NbPTS,pu+1))
    ALLOCATE(bDPHIz(NbPTS,pu+1))
    ALLOCATE(BASISz(pu+1))
    ALLOCATE(DBASISz(pu+1))
    ALLOCATE(L2U(pu+1,L2int))
    ALLOCATE(AU(pu+1,pu+1))
    ALLOCATE(CU(pu+1,L2int))
    ALLOCATE(B1U(pu+1))
    ALLOCATE(B2U(pu+1))
    ALLOCATE(L3U(p+1,LE))
ELSEIF (REGION == 2) THEN                           
    ALLOCATE(elemPTSz(NPTS_2D,2))
    ALLOCATE(elemWTSz(NPTS_2D))
    ALLOCATE(L2PTSz(L2_2D,2))
    ALLOCATE(L2WTSz(L2_2D))
    ALLOCATE(bPTSz(pu+1,2,2))
    ALLOCATE(bPTSx(pu+1,2,2))
    ALLOCATE(bWTSz(pu+1))
    ALLOCATE(PHIz(NPTS_2D,NDOF))
    ALLOCATE(DPHI_2D(NPTS_2D,NDOF,2))
    ALLOCATE(L2PHIz(L2_2D,NDOF))
    ALLOCATE(L2DPHI_2D(L2_2D,NDOF,2))
    ALLOCATE(bPHI_2D(pu+1,NDOF,2))
    ALLOCATE(bPHI_2Dx(pu+1,NDOF,2))    
    ALLOCATE(bDPHI_2D(pu+1,NDOF,4))
    ALLOCATE(PHIzeta(NPTS_2D,p+1))
    ALLOCATE(PHIzeta2(L2_2D,p+1))
    ALLOCATE(DPHIzeta(L2_2D,p+1))
    ALLOCATE(BASIS_2D(NDOF))
    ALLOCATE(DBASIS_2D(NDOF,2))
    ALLOCATE(L2U(NDOF,L2_2D))
    ALLOCATE(A_2D(NDOF,NPTS_2D,2))
    ALLOCATE(CU(NDOF,L2_2D))
    ALLOCATE(BU(NDOF,pu+1,2))
    ALLOCATE(BX(NDOF,pu+1,2))
    ALLOCATE(L2ptsNEW(L2int,2,NPTS_2D))   
ELSEIF (REGION == 3) THEN
    ALLOCATE(elemPTS_3D(NPTS_3D,3))
    ALLOCATE(elemWTS_3D(NPTS_3D))
    ALLOCATE(L2PTS_3D(L2_3D,3))
    ALLOCATE(L2WTS_3D(L2_3D))
    ALLOCATE(BASIS_3D(NDOF_3D))
    ALLOCATE(DBASIS_3D(NDOF_3D,3))  
    ALLOCATE(BASIS_3Dw(NDOF_3Dw))
    ALLOCATE(DBASIS_3Dw(NDOF_3Dw,3))      
    ALLOCATE(PHI_3D(NPTS_3D,NDOF_3D))
    ALLOCATE(DPHI_3D(NPTS_3D,NDOF_3D,3))     
    ALLOCATE(L2PHI_3D(L2_3D,NDOF_3D))
    ALLOCATE(L2DPHI_3D(L2_3D,NDOF_3D,3))    
    ALLOCATE(bPTS_3D(NPTS_2D,3,6))
    ALLOCATE(bPTS_3D2D(L2_2D,3))
    ALLOCATE(bPHI_3D(NPTS_2D,NDOF_3D,6))
    ALLOCATE(bPHI_3D2D(L2_2D,NDOF_3D))
    ALLOCATE(A_3D(NDOF_3D,NPTS_3D,3))
    ALLOCATE(A2(NDOF,NPTS_2D))
    ALLOCATE(L2_U3D(NDOF_3D,L2_3D))
    ALLOCATE(C_3D(NDOF_3D,L2_3D))
    ALLOCATE(B_3D(NDOF_3D,NPTS_2D,6))
    ALLOCATE(PHIzeta_3D(L2_3D,NDOF))
    ALLOCATE(DPHIzeta_3D(L2_3D,NDOF,2))   
    ALLOCATE(AwT(NDOF_2D,TP_2D,2))
    ALLOCATE(AwB(NDOF_3D,NPTS_3D,3))
    ALLOCATE(L2wT(NDOF_3D,L2_3D))
    ALLOCATE(L2wB(NDOF_3D,L2_3D))
    ALLOCATE(CwT(NDOF_3D,L2_3D))
    ALLOCATE(CwB(NDOF_3D,L2_3D))
    ALLOCATE(BwT(NDOF_2D,LRK,4))
    ALLOCATE(BwB(NDOF_3D,NPTS_2D,6))
    ALLOCATE(b_PHI_2D3Dx(NPTS_2D,NDOF,2)) 
    ALLOCATE(b_PHI_2D3Dy(NPTS_2D,NDOF,2))   
    ALLOCATE(PHI_2D3D(NPTS_3D,NDOF))
    ALLOCATE(PHIw(TP_2D,NDOF_3D,nrkW)) 
    ALLOCATE(PHI_3Dw(NPTS_3D,NDOF_3Dw))
    ALLOCATE(bPHIw(LRK,NDOF_3D,nrkW,4))  
    ALLOCATE(bPHI_3Dw(NPTS_2D,NDOF_3Dw,2))
    ALLOCATE(L2PHI_3Dw(L2_3D,NDOF_3Dw))
ENDIF

ALLOCATE(KPTS_1D(NPTS_K,1))
ALLOCATE(KWTS_1D(NPTS_K))
ALLOCATE(KPTS_2D(NPTS_K2D,2))
ALLOCATE(KWTS_2D(NPTS_K2D))
ALLOCATE(L2K(NDOF,NPTS_K2D))
KPTS_1D = 0.0d0
KWTS_1D = 0.0d0
KPTS_2D = 0.0d0
KWTS_2D = 0.0d0       

RETURN
END SUBROUTINE MATRIX_ALLOCATION
!---------------------------------------------------------------------
!
!                   Problem Data Subroutine
!                 Written by Colton J. Conroy
!                 @ the Isaac Newton Institute
!                          12.13.12
!
!---------------------------------------------------------------------
SUBROUTINE PROBLEM_DATA()

USE global
IMPLICIT NONE

IF (PROB == 1) THEN
    PROB_DATA%amp    = 0.500000000000000d0
    PROB_DATA%freq   = 0.000140750000000d0
    PROB_DATA%phase  = 0.000000000000000d0
    PROB_DATA%k      = 0.010000000000000d0
    PROB_DATA%g      = 9.806650000000000d0
    PROB_DATA%N0     = 0.001500000000000d0
    PROB_DATA%h0     = 1.000000000000000d0
    PROB_DATA%x0     = 60000.000000000000000d0
    PROB_DATA%xN     = 220000.0000000000d0
    PROB_DATA%y0     = 0.000000000000000d0
    PROB_DATA%yN     = 45000.00000000000d0 
    PROB_DATA%Hslope = 0.111111111111111d-8 
    PROB_DATA%ampU   = 0.0d0
ELSEIF (PROB == 2) THEN
    PROB_DATA%amp    = 0.500000000000000d0
    PROB_DATA%freq   = 0.000140750000000d0
    PROB_DATA%phase  = 0.000000000000000d0
    PROB_DATA%k      = 0.010000000000000d0
    PROB_DATA%g      = 9.806650000000000d0
    PROB_DATA%N0     = 0.001500000000000d0
    PROB_DATA%h0     = 1.000000000000000d0
    PROB_DATA%x0     = 0.000000000000000d0
    PROB_DATA%xN     = 160000.0000000000d0
    PROB_DATA%y0     = 0.000000000000000d0
    PROB_DATA%yN     = 160000.0000000000d0 
    PROB_DATA%Hslope = 0.0000000000d0   
    PROB_DATA%ampU   = 0.0d0
ELSEIF (PROB == 3) THEN
    PROB_DATA%amp    = 0.500000000000000d0
    PROB_DATA%freq   = 0.000140750000000d0
    PROB_DATA%phase  = 0.000000000000000d0
    PROB_DATA%k      = 0.010000000000000d0
    PROB_DATA%g      = 9.806650000000000d0
    PROB_DATA%N0     = 0.050000000000000d0
    PROB_DATA%h0     = 8.000000000000000d0
    PROB_DATA%x0     = 0.000000000000000d0
    PROB_DATA%xN     = 160000.0000000000d0
    PROB_DATA%y0     = 0.000000000000000d0
    PROB_DATA%yN     = 160000.0000000000d0 
    PROB_DATA%Hslope = 0.00020d0    
    PROB_DATA%ampU   = 0.5000d0
ELSEIF (PROB == 4) THEN
    PROB_DATA%amp    = 0.500000000000000d0
    PROB_DATA%freq   = 0.000140750000000d0
    PROB_DATA%phase  = 0.000000000000000d0
    PROB_DATA%k      = 0.010000000000000d0
    PROB_DATA%g      = 9.806650000000000d0
    PROB_DATA%N0     = 0.050000000000000d0
    PROB_DATA%h0     = 5.000000000000000d0
    PROB_DATA%x0     = -50.000000000000000d0
    PROB_DATA%xN     = 50.0000000000d0
    PROB_DATA%y0     = 0.000000000000000d0
    PROB_DATA%yN     = 50.0000000000d0 
    PROB_DATA%Hslope = 0.00000d0    
    PROB_DATA%ampU   = 0.000d0   
ELSEIF (PROB == 5) THEN
    PROB_DATA%amp    = 0.100000000000000d0
    PROB_DATA%freq   = 1.000000000000000d-12
    PROB_DATA%phase  = 0.000000000000000d0
    PROB_DATA%g      = 9.806650000000000d0
    PROB_DATA%h0     = 1.000000000000000d0
    PROB_DATA%x0     = 0.000000000000000d0
    PROB_DATA%xN     = 100000.0000000000d0
    PROB_DATA%y0     = 0.000000000000000d0
    PROB_DATA%yN     = 50.0000000000d0 
    PROB_DATA%Hslope = 0.00000d0    
    PROB_DATA%ampU   = 0.00000d0  
    PROB_DATA%r0     = 24.0000d0
    PROB_DATA%rS     =  1.2000d0
    PROB_DATA%S      =  0.5000d0
    PROB_DATA%f      =  0.99460d-4
    PROB_DATA%E      =  0.01
    PROB_DATA%N0     =  PROB_DATA%E*PROB_DATA%f*PROB_DATA%h0**2.0d0
    PROB_DATA%gR     =  PROB_DATA%g*PROB_DATA%rS/PROB_DATA%r0
    PROB_DATA%Lr     =  (PROB_DATA%gR*PROB_DATA%h0/PROB_DATA%f**2.0d0)**(1.0d0/2.0d0)
    PROB_DATA%lambda =  SQRT(PROB_DATA%E)
    PROB_DATA%k      =  PROB_DATA%N0/(PROB_DATA%lambda*PROB_DATA%h0)
    PROB_DATA%Lx     =  PROB_DATA%xN - PROB_DATA%x0 
ELSEIF (PROB == 0) THEN
    PROB_DATA%amp    = 0.000000000000000d0
    PROB_DATA%freq   = 0.000000000000000d0
    PROB_DATA%phase  = 0.000000000000000d0
    PROB_DATA%g      = 9.806650000000000d0
    PROB_DATA%h0     = 0.000000000000000d0
    PROB_DATA%x0     = 0.000000000000000d0
    PROB_DATA%xN     = 0.0000000000d0
    PROB_DATA%y0     = 0.000000000000000d0
    PROB_DATA%yN     = 0.0000000000d0 
    PROB_DATA%Hslope = 0.00000d0    
    PROB_DATA%ampU   = 0.00000d0  
    PROB_DATA%r0     = 0.0000d0
    PROB_DATA%rS     = 0.0000d0
    PROB_DATA%S      = 0.0000d0
    PROB_DATA%f      = 9.3490d-5
    PROB_DATA%E      = 0.00d0
    PROB_DATA%N0     = 100.00d0
    PROB_DATA%gR     = 0.00d0
    PROB_DATA%Lr     = 0.00d0
    PROB_DATA%lambda = 0.00d0
    PROB_DATA%k      = 0.00010d0
    PROB_DATA%Lx     =  PROB_DATA%xN - PROB_DATA%x0     
ENDIF

IF (debug == 5) THEN
    write(*,*)PROB_DATA%N0
    write(*,*)'-----------'
    write(*,*)PROB_DATA%gR
    write(*,*)'-----------'
    write(*,*)PROB_DATA%Lr
    write(*,*)'-----------'
    write(*,*)PROB_DATA%lambda
    write(*,*)'-----------'
    write(*,*)PROB_DATA%k
    write(*,*)'-----------'
    write(*,*)PROB_DATA%Lx
    write(*,*)'-----------'
    write(*,*)PROB_DATA%h0
    write(*,*)'-----------'
    write(*,*)PROB_DATA%x0
    write(*,*)'-----------'
    write(*,*)PROB_DATA%xN
    write(*,*)'-----------'
    write(*,*)PROB_DATA%y0
    write(*,*)'-----------'
    write(*,*)PROB_DATA%yN
    write(*,*)'-----------'
    write(*,*)PROB_DATA%Hslope
ENDIF

RETURN
END SUBROUTINE PROBLEM_DATA
!---------------------------------------------------------------------
!
!      Closure for Lynch/Loder Stratified Shallow Sea Test Case
!                  Written by Colton J. Conroy
!                           @ the C.H.I.L
!                               8.14.14
!
!---------------------------------------------------------------------
SUBROUTINE LL3D_CLOSURE(xCOORD,Ax1,Vavg)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: jj
REAL(real_p) :: Ax1, Ux1, INVJAC, Vavg
REAL(real_p) :: S, r0, rS, xcheck, zcheck, rhoB, h0, Lx
REAL(real_p), DIMENSION(L2int) :: A1val, A2val
REAL(real_p), DIMENSION(nelemsZ) :: UA1, UA2
COMPLEX(16)  :: xCOORD, zCOORD, f, lambda, Unorm, Ax2, Ax3
COMPLEX(16)  :: gamma, E, Ro, d1a, d2a, d1, d2, UIa, UIb, UIc, UId
COMPLEX(16)  :: UI, UII, UIII, UIIa, UIIb, UT
COMPLEX(16)  :: IITopA, IITopB, IIBot, IIITopA, IIITopB, IIIBot
COMPLEX(16)  :: UIIIa, UIIIb, UIIIc, Lr, Lxc, i1, PIc, d1b, d2b
COMPLEX(16), PARAMETER :: m = -1.0d0


i1     = SQRT(m)
PIc    = dcmplx(PI, 0.0d0)
Lxc    = dcmplx(PROB_DATA%Lx, 0.0d0)
f      = dcmplx(PROB_DATA%f, 0.0d0)
lambda = dcmplx(PROB_DATA%lambda, 0.0d0)
E      = dcmplx(PROB_DATA%E, 0.0d0)
Lr     = dcmplx(PROB_DATA%Lr, 0.0d0)
Lx     = PROB_DATA%Lx
h0     = PROB_DATA%h0
r0     = PROB_DATA%r0
rS     = PROB_DATA%rS
S      = PROB_DATA%S
Unorm  = (PIc/4.0d0)*f*Lr*(Lr/Lx)*SIN(PIc*xCOORD/Lx)
Ro     = Unorm/(f*Lx)
Ax2    = (PIc/4.0d0)*(Lr/Lx)**2.0d0*S/(E*Ro)*SIN(PIc*xCOORD/Lx)
Ax3    = (1.0d0/4.0d0)*(Lr/Lx)**2.0d0*1.0d0/(E*Ro)*SIN(PIc*xCOORD/Lx)
gamma  = (2.0d0*E)**(-1.0d0/2.0d0)

!....Contribution from cross frontal sea surface slope

d1a = 1.0d0 - lambda*gamma*(1.0d0 + i1)
d2a = 1.0d0 + lambda*gamma*(1.0d0 + i1)
d1b = i1 + lambda*gamma*(1.0d0 - i1)
d2b = i1 - lambda*gamma*(1.0d0 - i1)

d1  = EXP(-gamma)*(d1a*COS(gamma) - d1b*SIN(gamma))
d2  = EXP(gamma)*(d2a*COS(gamma) + d2b*SIN(gamma))

UIa = (-i1)/(2.0d0*gamma**2.0d0*(d1+d2))
UId = i1/(2.0d0*gamma**2.0d0)

!....Contribution from cross-frontal gradient in depth-avgeraged density

IITopA = (-d2*(1.0d0 + i1) + i1*2.0d0*gamma*(lambda + 1.0d0));
IITopB = (d1*(1.0d0  + i1) + i1*2.0d0*gamma*(lambda + 1.0d0));
IIBot  = 4.0d0*gamma**3.0d0*(d1 + d2);

!....Contribution from cross-frontal gradient in stratification

IIITopA = PIc*(PIc**2.0d0 - i1/E)*(-lambda - d2/(2.0d0*gamma)*(1.0d0 - i1));
IIITopB = PIc*(PIc**2.0d0 - i1/E)*(-lambda - d1/(2.0d0*gamma)*(1.0d0 - i1));
IIIBot  = (d1 + d2)*(PIc**4.0d0 + 1.0d0/E**2.0d0);
UIIIc   = (PIc**2.0d0 - i1/E)/(PIc**4.0d0 + 1.0d0/E**2.0d0);

DO ii = 1, nelemsZ

    INVJAC = 1.0d0/zELEM_jac(ii)

    DO jj = 1, L2int
    
        zCOORD = (1.0d0/2.0d0)*((1.0d0 - L2PTS(jj,1))*Z(ii) + (1.0d0 + L2PTS(jj,1))*Z(ii+1))
        
        !....Cross frontal sea surface slope
        
        UIb = EXP(gamma*zCOORD)*(COS(gamma*zCOORD) + i1*SIN(gamma*zCOORD))
        UIc = EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD))
        UI  = UIa*(UIb + UIc) + UId
        
        !....Cross-frontal gradient in depth-avgeraged density
        
        UIIa = IITopA/IIBot*EXP(gamma*zCOORD)*(COS(gamma*zCOORD)  + i1*SIN(gamma*zCOORD));
        UIIb = IITopB/IIBot*EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD));
        UII  = UIIa + UIIb + i1*(zCOORD/(2.0d0*gamma**2.0d0));
        
        !....Cross-frontal gradient in stratification
        
        UIIIa = IIITopA/IIIBot*(EXP(gamma*zCOORD)*(COS(gamma*zCOORD) + i1*SIN(gamma*zCOORD)));
        UIIIb = IIITopB/IIIBot*(EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD)));
        UIII  = UIIIa + UIIIb + UIIIc*SIN(PIc*zCOORD);
        
        A1val(jj) = -(real(Ax2*UII) + real(Ax3*UIII))    !....Numerator
        A2val(jj) = real(UI)                             !....Denominator
        
    ENDDO
    
    UA1(ii) = DOT_PRODUCT(L2wts,A1val)*INVJAC   !....Integrate over each z element
    UA2(ii) = DOT_PRODUCT(L2wts,A2val)*INVJAC
    
ENDDO    

!....Sum up all the elements

Ax1 = SUM(UA1)
Ux1 = SUM(UA2)

!....Closure coefficient

Ax1 = Ax1/Ux1

!....Calculate depth-integrated velocity in y-dir, i.e., the imag(U)

A1val = 0.0d0
UA1   = 0.0d0

DO ii = 1, nelemsZ

    INVJAC = 1.0d0/zELEM_jac(ii)

    DO jj = 1, L2int
    
        zCOORD = (1.0d0/2.0d0)*((1.0d0 - L2PTS(jj,1))*Z(ii) + (1.0d0 + L2PTS(jj,1))*Z(ii+1))
        
        !....Cross frontal sea surface slope
        
        UIb = EXP(gamma*zCOORD)*(COS(gamma*zCOORD) + i1*SIN(gamma*zCOORD))
        UIc = EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD))
        UI  = UIa*(UIb + UIc) + UId
        
        !....Cross-frontal gradient in depth-avgeraged density
        
        UIIa = IITopA/IIBot*EXP(gamma*zCOORD)*(COS(gamma*zCOORD)  + i1*SIN(gamma*zCOORD));
        UIIb = IITopB/IIBot*EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD));
        UII  = UIIa + UIIb + i1*(zCOORD/(2.0d0*gamma**2.0d0));
        
        !....Cross-frontal gradient in stratification
        
        UIIIa = IIITopA/IIIBot*(EXP(gamma*zCOORD)*(COS(gamma*zCOORD) + i1*SIN(gamma*zCOORD)));
        UIIIb = IIITopB/IIIBot*(EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD)));
        UIII  = UIIIa + UIIIb + UIIIc*SIN(PIc*zCOORD);
        
        A1val(jj) = imag(Ax1*UI) + imag(Ax2*UII) + imag(Ax3*UIII)    
        
    ENDDO
    
    UA1(ii) = DOT_PRODUCT(L2wts,A1val)*INVJAC   !....Integrate over each z element
    
ENDDO

Vavg = SUM(UA1)*Unorm

END SUBROUTINE LL3D_CLOSURE
!---------------------------------------------------------------------
!
!            Lynch/Loder Stratified Shallow Sea Test Case
!                  Written by Colton J. Conroy
!                           @ the C.H.I.L
!                               8.13.14
!
!---------------------------------------------------------------------
SUBROUTINE LL3D_SOLUTIONS(xCOORD,zCOORD,Vu,Vv,DSE,Ax1,Rp,DRDX,DENSITYh,DENSITY)

USE precisions
USE global
IMPLICIT NONE

REAL(real_p) :: Vu, Vv, Rp, DRDX, DENSITYh, DENSITY, Ax1, g
REAL(real_p) :: S, r0, rS, xcheck, zcheck, rhoB, h0, Lx, DSE
COMPLEX(16)  :: xCOORD, zCOORD, f, lambda, Unorm, Ax2, Ax3
COMPLEX(16)  :: gamma, E, Ro, d1a, d2a, d1, d2, UIa, UIb, UIc
COMPLEX(16)  :: UI, UII, UIII, TopA, TopB, Bot, UIIa, UIIb, UT
COMPLEX(16)  :: UIIIa, UIIIb, UIIIc, Lr, Lxc, i1, PIc, d1b, d2b
COMPLEX(16), PARAMETER :: m = -1.0d0

!....Constants

i1     = SQRT(m)
PIc    = dcmplx(PI, 0.0d0)
Lxc    = dcmplx(PROB_DATA%Lx, 0.0d0)
f      = dcmplx(PROB_DATA%f, 0.0d0)
lambda = dcmplx(PROB_DATA%lambda, 0.0d0)
E      = dcmplx(PROB_DATA%E, 0.0d0)
Lr     = dcmplx(PROB_DATA%Lr, 0.0d0)
Lx     = PROB_DATA%Lx
h0     = PROB_DATA%h0
r0     = PROB_DATA%r0
rS     = PROB_DATA%rS
S      = PROB_DATA%S
g      = PROB_DATA%g
Unorm  = (PIc/4.0d0)*f*Lr*(Lr/Lxc)*SIN(PIc*xCOORD/Lxc)
Ro     = Unorm/(f*Lxc)
Ax2    = (PIc/4.0d0)*(Lr/Lxc)**2.0d0*S/(E*Ro)*SIN(PIc*xCOORD/Lxc)
Ax3    = (1.0d0/4.0d0)*(Lr/Lxc)**2.0d0*1.0d0/(E*Ro)*SIN(PIc*xCOORD/Lxc)
gamma  = (2.0d0*E)**(-1.0d0/2.0d0)

!....Density

xcheck = real(xCOORD)
zcheck = real(zCOORD)

IF (xcheck <= 0.0d0) THEN
    rhoB = rS/2.0d0
ELSEIF (xcheck >= Lx) THEN
    rhoB = 0.0d0
ELSE
    rhoB = rS/4.0d0*(1.0d0 + COS(PI*xcheck/Lx))
ENDIF

DENSITY  = r0 + S*rhoB - rhoB*COS(PI*zcheck)
DENSITYh = (DENSITY - r0)/r0
Rp       = rS/(4.0d0*r0)*(1.0d0 + COS(PI*xcheck/Lx))*(h0/PI*SIN(PI*zcheck) - S*zcheck)
DRDX     = (pi*rS)/(4.0d0*Lx*r0)*SIN(PI*xcheck/Lx)*(S*zcheck - h0/PI*SIN(PI*zcheck))    

!....Contribution from cross frontal sea surface slope

d1a = 1.0d0 - lambda*gamma*(1.0d0 + i1)
d2a = 1.0d0 + lambda*gamma*(1.0d0 + i1)
d1b = i1 + lambda*gamma*(1.0d0 - i1)
d2b = i1 - lambda*gamma*(1.0d0 - i1)

d1  = EXP(-gamma)*(d1a*COS(gamma) - d1b*SIN(gamma))
d2  = EXP(gamma)*(d2a*COS(gamma) + d2b*SIN(gamma))

UIa = (-i1)/(2.0d0*gamma**2.0d0*(d1+d2))
UIb = EXP(gamma*zCOORD)*(COS(gamma*zCOORD) + i1*SIN(gamma*zCOORD))
UIc = EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD))

UI  = UIa*(UIb + UIc) + i1/(2.0d0*gamma**2.0d0)

!....Contribution from cross-frontal gradient in depth-avgeraged density

TopA = (-d2*(1.0d0 + i1) + i1*2.0d0*gamma*(lambda + 1.0d0));
TopB = (d1*(1.0d0  + i1) + i1*2.0d0*gamma*(lambda + 1.0d0));
Bot  = 4.0d0*gamma**3.0d0*(d1 + d2);

UIIa = TopA/Bot*EXP(gamma*zCOORD)*(COS(gamma*zCOORD)  + i1*SIN(gamma*zCOORD));
UIIb = TopB/Bot*EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD));

UII  = UIIa + UIIb + i1*(zCOORD/(2.0d0*gamma**2.0d0));

!....Contribution from cross-frontal gradient in stratification

TopA = PIc*(PIc**2.0d0 - i1/E)*(-lambda - d2/(2.0d0*gamma)*(1.0d0 - i1));
TopB = PIc*(PIc**2.0d0 - i1/E)*(-lambda - d1/(2.0d0*gamma)*(1.0d0 - i1));
Bot  = (d1 + d2)*(PIc**4.0d0 + 1.0d0/E**2.0d0);

UIIIa = TopA/Bot*(EXP(gamma*zCOORD)*(COS(gamma*zCOORD) + i1*SIN(gamma*zCOORD)));
UIIIb = TopB/Bot*(EXP(-gamma*zCOORD)*(COS(gamma*zCOORD) - i1*SIN(gamma*zCOORD)));
UIIIc = (PIc**2.0d0 - i1/E)/(PIc**4.0d0 + 1.0d0/E**2.0d0);

UIII  = UIIIa + UIIIb + UIIIc*SIN(PIc*zCOORD);

!....Total Solution

UT  = Ax1*UI + Ax2*UII + Ax3*UIII;

Vu  = real(UT)*Unorm
Vv  = imag(UT)*Unorm
DSE = real((Ax1*Unorm**2.0d0*E)/(g*Ro)) 


RETURN
END SUBROUTINE LL3D_SOLUTIONS
!---------------------------------------------------------------------
!
!              Smooth Solutions Non Linear 3D Test Case
!                   Written by Colton J. Conroy
!                            @ the C.H.I.L
!                                2.13.14
!
!---------------------------------------------------------------------
SUBROUTINE NL_SMOOTH_3D(xCOORD,yCOORD,zCOORD,solnT,Vu,Vv,Vw,H,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)

USE precisions  
IMPLICIT NONE 

REAL(real_p) :: xCOORD, yCOORD, zCOORD, Vu, Vv, Vw, solnT, DSEx, DSEy
REAL(real_p) :: H, Zb, SE, Uavg, Vavg, Fpce, Fxm, Fym, Fw, g, N0, k, FymA
REAL(real_p) :: Lx, Ly, amp, h0, s, freq, B1A, B2A, B2B, B3A, B3B, FxmA, ampU
REAL(real_p), PARAMETER :: PI = 3.1415926535897932d0

TYPE LYNCH
    REAL(real_p) :: amp
    REAL(real_p) :: freq 
    REAL(real_p) :: phase
    REAL(real_p) :: k
    REAL(real_p) :: g
    REAL(real_p) :: N0
    REAL(real_p) :: h0 
    REAL(real_p) :: x0 
    REAL(real_p) :: xN
    REAL(real_p) :: y0
    REAL(real_p) :: yN
    REAL(real_p) :: Hslope
    REAL(real_p) :: ampU
END TYPE LYNCH

TYPE(LYNCH) :: PROB_DATA

amp  = PROB_DATA%amp
s    = PROB_DATA%Hslope
g    = PROB_DATA%g
h0   = PROB_DATA%h0
freq = PROB_DATA%freq
Lx   = PROB_DATA%xN - PROB_DATA%x0
Ly   = PROB_DATA%yN - PROB_DATA%y0
N0   = PROB_DATA%N0
k    = PROB_DATA%k
ampU = PROB_DATA%ampU

!....2D solutions

!....Bathymetry

Zb = h0 + s*(xCOORD + yCOORD)

!....Surface Elevation

SE = amp*(cos(2.0d0*PI*xCOORD/Lx + freq*solnT) + cos(2.0d0*PI*yCOORD/Ly + freq*solnT))

DSEx = - amp*2.0d0*PI/Lx*sin(2.0d0*PI*xCOORD/Lx + freq*solnT)

DSEy = - amp*2.0d0*PI/Ly*sin(2.0d0*PI*yCOORD/Ly + freq*solnT)

!....Total height of water column

H = Zb + SE

!....3D solutions

Vu = exp(k/N0*H*zCOORD + ampU*sin(2*pi*xCOORD/Lx + solnT*freq));
Vv = exp(k/N0*H*zCOORD + ampU*sin(2*pi*yCOORD/Ly + solnT*freq));

!....Depth averaged solutions


Uavg = (N0*(exp(ampU*sin(freq*solnT + (2.0d0*pi*xCOORD)/Lx)) - exp(ampU*sin(freq*solnT + (2.0d0*pi*xCOORD)/Lx) &
     - (h0*k)/N0 - (k*s*xCOORD)/N0 - (k*s*yCOORD)/N0 - (amp*k*cos(freq*solnT + (2.0d0*pi*xCOORD)/Lx))/N0       &
     - (amp*k*cos(freq*solnT + (2.0d0*pi*yCOORD)/Ly))/N0)))/(k*(h0 + s*xCOORD + s*yCOORD + amp*cos(freq*solnT  &
     + (2.0d0*pi*xCOORD)/Lx) + amp*cos(freq*solnT + (2.0d0*pi*yCOORD)/Ly)))

Vavg = (N0*(exp(ampU*sin(freq*solnT + (2.0d0*pi*yCOORD)/Ly)) - exp(ampU*sin(freq*solnT + (2.0d0*pi*yCOORD)/Ly) &
     - (h0*k)/N0 - (k*s*xCOORD)/N0 - (k*s*yCOORD)/N0 - (amp*k*cos(freq*solnT + (2.0d0*pi*xCOORD)/Lx))/N0       &
     - (amp*k*cos(freq*solnT + (2.0d0*pi*yCOORD)/Ly))/N0)))/(k*(h0 + s*xCOORD + s*yCOORD + amp*cos(freq*solnT  &
     + (2.0d0*pi*xCOORD)/Lx) + amp*cos(freq*solnT + (2.0d0*pi*yCOORD)/Ly)))
     
!....Transformed vertical velocity

Vw = (2.0d0*Ly*N0*pi*ampU*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
   + 2.0d0*Lx*N0*pi*ampU*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
   + 2.0d0*Ly*N0*pi*ampU*zCOORD*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx) + 2.0d0*Lx*N0*pi*ampU*zCOORD*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*cos((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly) + Lx*Ly*k*s*zCOORD*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0) &
   + Lx*Ly*k*s*zCOORD*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0) &
   - 2.0d0*Ly*N0*pi*ampU*exp((amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))/N0)*exp((amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))/N0)*exp((k*s*xCOORD*zCOORD)/N0)*exp((k*s*yCOORD*zCOORD)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))*exp((h0*k*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
   - 2.0d0*Lx*N0*pi*ampU*exp((amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))/N0)*exp((amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))/N0)*exp((k*s*xCOORD*zCOORD)/N0)*exp((k*s*yCOORD*zCOORD)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))*exp((h0*k*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - Lx*Ly*k*s*zCOORD*exp((amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp((amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp((k*s*xCOORD*zCOORD)/N0)*exp((k*s*yCOORD*zCOORD)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))*exp((h0*k*zCOORD)/N0) - Lx*Ly*k*s*zCOORD*exp((amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp((amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp((k*s*xCOORD*zCOORD)/N0)*exp((k*s*yCOORD*zCOORD)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))*exp((h0*k*zCOORD)/N0) + 2.0d0*Ly*pi*amp*k*zCOORD*exp((amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp((amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp((k*s*xCOORD*zCOORD)/N0)*exp((k*s*yCOORD*zCOORD)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))*exp((h0*k*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Lx*pi*amp*k*zCOORD*exp((amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp((amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp((k*s*xCOORD*zCOORD)/N0)*exp((k*s*yCOORD*zCOORD)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))*exp((h0*k*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*ampU*zCOORD*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))*exp(-(amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
   - 2.0d0*Lx*N0*pi*ampU*zCOORD*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*cos((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly) - 2.0d0*Ly*pi*amp*k*zCOORD*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))*exp(-(amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
   - 2.0d0*Lx*pi*amp*k*zCOORD*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
   + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*sin((2.0d0*pi*yCOORD &
   + Ly*freq*solnT)/Ly))/(Lx*Ly*k)
   
!....Forcing Functions

!....PCE 

Fpce = -(Lx*Ly*amp*freq*k*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + Lx*Ly*amp*freq*k*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*ampU*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) - 2.0d0*Lx*N0*pi*ampU*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) - Lx*Ly*k*s*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0) &
     - Lx*Ly*k*s*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0) &
     + 2.0d0*Ly*N0*pi*ampU*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + 2.0d0*Lx*N0*pi*ampU*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + 2.0d0*Ly*pi*amp*k*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + 2.0d0*Lx*pi*amp*k*exp(-(h0*k)/N0)*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*exp(-(amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/N0)*exp(-(amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/N0)*exp(-(k*s*xCOORD)/N0)*exp(-(k*s*yCOORD)/N0)*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly))/(Lx*Ly*k)                                           

!....X-Momentum Equation

Fxm = -(Lx*Ly*h0*k**3.0d0*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    + 2.0d0*Ly*N0**2.0d0*pi*ampU*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - 2.0d0*Ly*N0**2.0d0*pi*ampU*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + Lx*Ly*amp*k**3.0d0*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + Lx*Ly*amp*k**3.0d0*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + Lx*Ly*k**3.0d0*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) + Lx*Ly*k**3.0d0*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    - Lx*Ly*N0*k*s*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0) + 2.0d0*Lx*N0**2.0d0*pi*ampU*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - Lx*Ly*N0*k*s*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    - 2.0d0*Lx*N0**2.0d0*pi*ampU*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - Lx*Ly*k**2.0d0*s**2.0d0*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    - Lx*Ly*k**2.0d0*s**2.0d0*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    + Ly*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*xCOORD &
    + Lx*freq*solnT))/Lx) + 2.0d0*Ly*N0*pi*amp*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - Lx*Ly*h0*k**2.0d0*s*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    + Ly*N0*pi*amp**2.0d0*g*k*sin((2.0d0*(2.0d0*pi*xCOORD + Lx*freq*solnT))/Lx) - Lx*Ly*h0*k**2.0d0*s*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0) - Lx*Ly*k**2.0d0*s**2.0d0*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD &
    - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) - Lx*Ly*k**2.0d0*s**2.0d0*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD &
    - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) + Lx*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*yCOORD + Ly*freq*solnT))/Ly) &
    + 2.0d0*Lx*N0*pi*amp*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*amp*ampU*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx)**2.0d0 - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*pi*amp*h0*k**2.0d0*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Lx*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD &
    - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*h0*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*ampU*h0*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD &
    + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*amp*g*h0*k*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Lx*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly)**2.0d0 + 2.0d0*Lx*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD &
    + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)**2.0d0 + (Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*xCOORD &
    + Lx*freq*solnT))/Lx))/2.0d0 + (Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*yCOORD &
    + Ly*freq*solnT))/Ly))/2.0d0 + 2.0d0*Ly*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Ly*N0*pi*ampU*h0*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + Lx*Ly*N0*amp*freq*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + Lx*Ly*N0*amp*freq*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD &
    - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + 2.0d0*Lx*pi*amp*h0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*amp**2.0d0*g*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)**2.0d0 - Lx*Ly*N0*ampU*freq*h0*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + 2.0d0*Ly*N0*pi*ampU*h0*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)**2.0d0 &
    + Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - 2.0d0*Ly*N0*pi*amp*ampU*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*ampU*h0*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Ly*N0*pi*ampU*k*s*xCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Ly*N0*pi*ampU*k*s*yCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Lx*pi*amp*k**2.0d0*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*Lx*pi*amp*k**2.0d0*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - Lx*Ly*N0*amp*ampU*freq*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)**2.0d0 &
    + 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx)**2.0d0 - 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)**2.0d0 &
    + Lx*Ly*amp*freq*h0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + Lx*Ly*amp*freq*h0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + 2.0d0*Lx*N0*pi*ampU*h0*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*pi*amp*k**2.0d0*s*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + 2.0d0*Ly*pi*amp*k**2.0d0*s*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Lx*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD &
    + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - 2.0d0*Lx*N0*pi*ampU*h0*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - 2.0d0*Lx*N0*pi*ampU*k*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + 2.0d0*Lx*N0*pi*ampU*k*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - 2.0d0*Lx*N0*pi*ampU*k*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*ampU*k*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*amp*g*k*s*xCOORD*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Ly*N0*pi*amp*g*k*s*yCOORD*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + Lx*Ly*amp*freq*k**2.0d0*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + Lx*Ly*amp*freq*k**2.0d0*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + Lx*Ly*amp*freq*k**2.0d0*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + Lx*Ly*amp*freq*k**2.0d0*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + 2.0d0*Lx*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - Lx*Ly*N0*amp*ampU*freq*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - Lx*Ly*N0*ampU*freq*k*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - Lx*Ly*N0*ampU*freq*k*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + 2.0d0*Ly*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Ly*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - 2.0d0*Ly*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - 2.0d0*Ly*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx))/(Lx*Ly*N0*k)

FxmA = k*exp(-(h0*k - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) &
     - (g*(2.0d0*s*(h0 + s*xCOORD + s*yCOORD) - (2.0d0*(Lx*s - 2.0d0*pi*amp*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))*(h0 + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + s*xCOORD + s*yCOORD))/Lx))/2.0d0 - k*exp(ampU*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx)) - (N0*((freq*exp(-(h0*k - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
     + k*s*xCOORD + k*s*yCOORD)/N0)*(N0*ampU*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)))/N0 - ampU*freq*exp(ampU*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)))/k - amp*g*s*(cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)) + (2.0d0*N0**2.0d0*((exp(-(h0*k &
     - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0)*(2.0d0*pi*amp*k*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) - Lx*k*s + 2.0d0*N0*pi*ampU*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)))/(Lx*N0) &
     - (2.0d0*pi*ampU*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/Lx)*(exp(-(h0*k - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) & 
     + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))))/(k**2.0d0*(h0 &
     + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + s*xCOORD &
     + s*yCOORD)) + (N0**2.0d0*((exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
     + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
     + k*s*xCOORD + k*s*yCOORD)/N0)*(2.0d0*pi*amp*k*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - Ly*k*s &
     + 2.0d0*N0*pi*ampU*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)))/(Ly*N0) - (2.0d0*pi*ampU*exp(ampU*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly))*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/Ly)*(exp(-(h0*k &
     - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))))/(k**2.0d0*(h0 + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + s*xCOORD + s*yCOORD)) - (N0**2.0d0*(Lx*s - 2.0d0*pi*amp*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))*(exp(-(h0*k - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
     + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)))**2.0d0)/(Lx*k**2.0d0*(h0 &
     + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + s*xCOORD &
     + s*yCOORD)**2.0d0) - (N0**2.0d0*(Ly*s - 2.0d0*pi*amp*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*(exp(-(h0*k &
     - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx)))*(exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) &
     - exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))))/(Ly*k**2.0d0*(h0 + amp*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + s*xCOORD + s*yCOORD)**2.0d0) &
     - (N0*exp(-(h0*k - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0)*(Ly*s &
     - 2.0d0*pi*amp*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*(exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))))/(Ly*k*(h0 &
     + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + s*xCOORD &
     + s*yCOORD))

Fym = -(Lx*Ly*h0*k**3.0d0*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) + 2.0d0*Lx*N0**2.0d0*pi*ampU*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0**2.0d0*pi*ampU*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + Lx*Ly*amp*k**3.0d0*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + Lx*Ly*amp*k**3.0d0*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + Lx*Ly*k**3.0d0*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) - Lx*Ly*N0*k*s*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) + Lx*Ly*k**3.0d0*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) + 2.0d0*Ly*N0**2.0d0*pi*ampU*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - Lx*Ly*N0*k*s*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    - 2.0d0*Ly*N0**2.0d0*pi*ampU*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - Lx*Ly*k**2.0d0*s**2.0d0*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) - Lx*Ly*k**2.0d0*s**2.0d0*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    + Lx*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*yCOORD + Ly*freq*solnT))/Ly) &
    + 2.0d0*Lx*N0*pi*amp*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - Lx*Ly*h0*k**2.0d0*s*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) &
    - Lx*Ly*N0*g*k*s**2.0d0*xCOORD - Lx*Ly*N0*g*k*s**2.0d0*yCOORD + Lx*N0*pi*amp**2.0d0*g*k*sin((2.0d0*(2.0d0*pi*yCOORD &
    + Ly*freq*solnT))/Ly) - Lx*Ly*h0*k**2.0d0*s*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD &
    - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) - Lx*Ly*k**2.0d0*s**2.0d0*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) - Lx*Ly*k**2.0d0*s**2.0d0*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0) + Ly*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*xCOORD &
    + Lx*freq*solnT))/Lx) + 2.0d0*Ly*N0*pi*amp*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - Lx*Ly*N0*g*h0*k*s &
    - 2.0d0*Lx*N0*pi*amp*ampU*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly)**2.0d0 - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + 2.0d0*Lx*pi*amp*h0*k**2.0d0*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*Ly*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - 2.0d0*Ly*N0*pi*ampU*h0*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Ly*N0*pi*ampU*h0*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + 2.0d0*Lx*N0*pi*amp*g*h0*k*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx)**2.0d0 + 2.0d0*Ly*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)**2.0d0 &
    + (Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*xCOORD &
    + Lx*freq*solnT))/Lx))/2.0d0 + (Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*(2.0d0*pi*yCOORD + Ly*freq*solnT))/Ly))/2.0d0 + 2.0d0*Lx*pi*amp**2.0d0*k**2.0d0*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*h0*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + Lx*Ly*N0*amp*freq*k*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + Lx*Ly*N0*amp*freq*k*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - Lx*Ly*amp*k**2.0d0*s*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*pi*amp*h0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + 2.0d0*Lx*N0*pi*amp**2.0d0*g*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)**2.0d0 &
    - Lx*Ly*N0*ampU*freq*h0*k*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*ampU*h0*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)**2.0d0 + Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + Lx*Ly*amp**2.0d0*freq*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - 2.0d0*Lx*N0*pi*amp*ampU*k*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*h0*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*k*s*xCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - 2.0d0*Lx*N0*pi*ampU*k*s*yCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*pi*amp*k**2.0d0*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Ly*pi*amp*k**2.0d0*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - Lx*Ly*N0*amp*ampU*freq*k*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)**2.0d0 + 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)**2.0d0 - 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)**2.0d0 + Lx*Ly*amp*freq*h0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + Lx*Ly*amp*freq*h0*k**2.0d0*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*ampU*h0*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Lx*pi*amp*k**2.0d0*s*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Lx*pi*amp*k**2.0d0*s*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*amp*ampU*k*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD &
    + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*ampU*h0*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - 2.0d0*Ly*N0*pi*ampU*k*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Ly*N0*pi*ampU*k*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD &
    + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - 2.0d0*Ly*N0*pi*ampU*k*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + 2.0d0*Ly*N0*pi*ampU*k*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*h0*k*zCOORD + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + 2.0d0*amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*k*s*xCOORD*zCOORD + 2.0d0*k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Lx*N0*pi*amp*g*k*s*xCOORD*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*amp*g*k*s*yCOORD*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + Lx*Ly*amp*freq*k**2.0d0*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + Lx*Ly*amp*freq*k**2.0d0*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + Lx*Ly*amp*freq*k**2.0d0*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + Lx*Ly*amp*freq*k**2.0d0*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + 2.0d0*Ly*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + 2.0d0*Ly*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - h0*k + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - 2.0d0*Ly*N0*pi*amp*ampU*k*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Ly*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - 2.0d0*Ly*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
    - Lx*Ly*N0*amp*ampU*freq*k*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    - Lx*Ly*N0*ampU*freq*k*s*xCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - Lx*Ly*N0*ampU*freq*k*s*yCOORD*exp((N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD &
    + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + 2.0d0*Lx*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD &
    + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + 2.0d0*Lx*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - h0*k - amp*k*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) - amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD - k*s*xCOORD - k*s*yCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*amp*ampU*k*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
    + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*k*s*xCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD &
    + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) - 2.0d0*Lx*N0*pi*ampU*k*s*yCOORD*zCOORD*exp((2.0d0*N0*ampU*sin((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly) + h0*k*zCOORD + amp*k*zCOORD*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*zCOORD*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD*zCOORD + k*s*yCOORD*zCOORD)/N0)*cos((2.0d0*pi*yCOORD &
    + Ly*freq*solnT)/Ly))/(Lx*Ly*N0*k)

FymA = k*exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) &
     - (g*(2.0d0*s*(h0 + s*xCOORD + s*yCOORD) - (2.0d0*(Ly*s - 2.0d0*pi*amp*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly))*(h0 + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + s*xCOORD + s*yCOORD))/Ly))/2.0d0 - k*exp(ampU*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly)) - (N0*((freq*exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
     + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
     + k*s*xCOORD + k*s*yCOORD)/N0)*(N0*ampU*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)))/N0 - ampU*freq*exp(ampU*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly))*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)))/k - amp*g*s*(cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)) + (N0**2.0d0*((exp(-(h0*k - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD &
     + k*s*yCOORD)/N0)*(2.0d0*pi*amp*k*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) - Lx*k*s + 2.0d0*N0*pi*ampU*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx)))/(Lx*N0) - (2.0d0*pi*ampU*exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx))/Lx)*(exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) &
     - exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))))/(k**2.0d0*(h0 + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + s*xCOORD + s*yCOORD)) + (2.0d0*N0**2.0d0*((exp(-(h0*k &
     - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0)*(2.0d0*pi*amp*k*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) - Ly*k*s + 2.0d0*N0*pi*ampU*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)))/(Ly*N0) &
     - (2.0d0*pi*ampU*exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))/Ly)*(exp(-(h0*k &
     - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) &
     + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly))))/(k**2.0d0*(h0 + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + s*xCOORD + s*yCOORD)) - (N0**2.0d0*(Ly*s - 2.0d0*pi*amp*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))*(exp(-(h0*k &
     - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly)))**2.0d0)/(Ly*k**2.0d0*(h0 &
     + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + s*xCOORD &
     + s*yCOORD)**2.0d0) - (N0**2.0d0*(Lx*s - 2.0d0*pi*amp*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*(exp(-(h0*k &
     - N0*ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx)))*(exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly))))/(Lx*k**2.0d0*(h0 &
     + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + s*xCOORD &
     + s*yCOORD)**2.0d0) - (N0*exp(-(h0*k - N0*ampU*sin((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + amp*k*cos((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0)*(Lx*s &
     - 2.0d0*pi*amp*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))*(exp(-(h0*k - N0*ampU*sin((2.0d0*pi*xCOORD &
     + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*k*cos((2.0d0*pi*yCOORD &
     + Ly*freq*solnT)/Ly) + k*s*xCOORD + k*s*yCOORD)/N0) - exp(ampU*sin((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx))))/(Lx*k*(h0 &
     + amp*cos((2.0d0*pi*xCOORD + Lx*freq*solnT)/Lx) + amp*cos((2.0d0*pi*yCOORD + Ly*freq*solnT)/Ly) &
     + s*xCOORD + s*yCOORD))
     
Fw = 0.0d0
   
RETURN
END SUBROUTINE NL_SMOOTH_3D
!---------------------------------------------------------------------
!
!                   Vertical Velocity Test Case
!                   Written by Colton J. Conroy
!                            @ the C.H.I.L
!                               12.16.13
!
!---------------------------------------------------------------------
SUBROUTINE Vertical_Test(xCOORD,yCOORD,zCOORD,Vu,Vv,Vw,DVw,H,Zb)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p), PARAMETER :: w_test = 1
REAL(real_p) :: xCOORD, yCOORD, zCOORD, solnT, Vu, Vv, Vw
REAL(real_p) :: H, Zb, amp, L, SE, Bxy, DVw

amp = PROB_DATA%amp

L   = PROB_DATA%xN - PROB_DATA%x0

IF (w_test == 1) THEN
    Zb  = -10.0d0 + cos(2.0d0*PI*xCOORD/L) + cos(2.0d0*PI*yCOORD/L)
    SE  = amp*(cos(2.0d0*PI*xCOORD/L) + cos(2.0d0*PI*yCOORD/L))
    H   = Zb + SE
    Vu  = -2.0d0*PI*sin(2.0d0*PI*zCOORD)*amp*sin(2.0d0*PI*xCOORD/L)
    Vv  = -2.0d0*PI*sin(2.0d0*PI*zCOORD)*amp*sin(2.0d0*PI*yCOORD/L)
    DVw = (4.0d0*amp*pi**2.0d0*cos((2.0d0*pi*xCOORD)/L)*sin(2.0d0*pi*zCOORD)*(cos((2.0d0*pi*xCOORD)/L)  &
    + cos((2.0d0*pi*yCOORD)/L) + amp*(cos((2.0d0*pi*xCOORD)/L) + cos((2.0d0*pi*yCOORD)/L)) - 10.0d0))/L &
    - 2.0d0*amp*pi*sin((2.0d0*pi*yCOORD)/L)*sin(2.0d0*pi*zCOORD)*((2.0d0*pi*sin((2.0d0*pi*yCOORD)/L))/L &
    + (2.0d0*amp*pi*sin((2.0d0*pi*yCOORD)/L))/L) - 2.0d0*amp*pi*sin((2.0d0*pi*xCOORD)/L)*sin(2.0d0*pi*zCOORD)*((2.0d0*pi*sin((2.0d0*pi*xCOORD)/L))/L &
    + (2.0d0*amp*pi*sin((2.0d0*pi*xCOORD)/L))/L) + (4.0d0*amp*pi**2.0d0*cos((2.0d0*pi*yCOORD)/L)*sin(2.0d0*pi*zCOORD)*(cos((2.0d0*pi*xCOORD)/L) &
    + cos((2.0d0*pi*yCOORD)/L) + amp*(cos((2.0d0*pi*xCOORD)/L) + cos((2.0d0*pi*yCOORD)/L)) - 10.0d0))/L
    Vw  = (40.0d0*amp*pi*cos((2.0d0*pi*xCOORD)/L)*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L          &
    + (40.0d0*amp*pi*cos((2.0d0*pi*yCOORD)/L)*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L          &
    - (4.0d0*amp**2.0d0*pi*(cos(pi*zCOORD)**2.0d0 - 1.0d0)*(cos((2.0d0*pi*xCOORD)/L)**2.0d0 - 1.0d0))/L &
    - (4.0d0*amp**2.0d0*pi*(cos(pi*zCOORD)**2.0d0 - 1.0d0)*(cos((2.0d0*pi*yCOORD)/L)**2.0d0 - 1.0d0))/L &
    - (4.0d0*amp*pi*cos((2.0d0*pi*xCOORD)/L)**2.0d0*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L             &
    - (4.0d0*amp*pi*cos((2.0d0*pi*yCOORD)/L)**2.0d0*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L             &
    - (4.0d0*amp*pi*(cos(pi*zCOORD)**2.0d0 - 1.0d0)*(cos((2.0d0*pi*xCOORD)/L)**2.0d0 - 1.0d0))/L   &
    - (4.0d0*amp**2.0d0*pi*cos((2.0d0*pi*xCOORD)/L)**2.0d0*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L      &
    - (4.0d0*amp*pi*(cos(pi*zCOORD)**2.0d0 - 1.0d0)*(cos((2.0d0*pi*yCOORD)/L)**2.0d0 - 1.0d0))/L   &
    - (4.0d0*amp**2.0d0*pi*cos((2.0d0*pi*yCOORD)/L)**2.0d0*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L      &
    - (8.0d0*amp*pi*cos((2.0d0*pi*xCOORD)/L)*cos((2.0d0*pi*yCOORD)/L)*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L & 
    - (8.0d0*amp**2.0d0*pi*cos((2.0d0*pi*xCOORD)/L)*cos((2.0d0*pi*yCOORD)/L)*(cos(pi*zCOORD)**2.0d0 - 1.0d0))/L
ELSEIF (w_test == 2) THEN    
    Zb = -10.0d0
    SE = cos(2.0d0*PI*xCOORD/L) + cos(2.0d0*PI*yCOORD/L)
    H  = Zb + SE
    Vu = 2.0d0*PI*sin(2.0d0*PI*zCOORD)*xCOORD/L
    Vv = 2.0d0*PI*sin(2.0d0*PI*zCOORD)*yCOORD/L
    Bxy = 2.0d0*pi*xCOORD/(L*L)*sin(2.0d0*pi*xCOORD/L)+2.0d0*pi*yCOORD/(L*L)*sin(2.0d0*pi*yCOORD/L)+(2.0d0/L)*(10.0d0-cos(2.0d0*PI*xCOORD/L)-cos(2.0d0*PI*yCOORD/L))
    Vw = (-cos(2.0d0*pi*zCOORD)+1)*Bxy
    DVw = (2.0d0*pi*sin(2.0d0*pi*zCOORD))*Bxy    
ENDIF    

RETURN
END SUBROUTINE Vertical_Test
!---------------------------------------------------------------------
!
!          Lynch Test Case Exact Solution (Quadratic Bathymetry)
!                       Written by Colton J. Conroy
!                             @ the C.H.I.L
!                                 8.26.13
! 
!---------------------------------------------------------------------
SUBROUTINE LG3D_SOLUTIONS_H(xCOORD,zCOORD,solnT,grad_SE,Vu,SE,Vavg,hmax,N,k)

USE precisions
USE global
IMPLICIT NONE

REAL(real_p) :: KHN, grad_SE, Vu, SE, Vavg, wt, hmax
REAL(real_p) :: N, k, N0H0, VwT
COMPLEX(16) :: cgrad_SE, CV, CSE, CVavg,V0,s1, s2, C1, C2, C3
COMPLEX(16) :: xCOORD, zCOORD, solnT, VZ, VBARR, zH, x1, x2
COMPLEX(16) :: amplitude, iwt, lambda, tau, betaA, i1, h
COMPLEX(16) :: Nc, KHNc, alpha1, alpha2, gamma, C11, C22, C33
COMPLEX(16) :: Delta, C12, VW
COMPLEX(16), PARAMETER :: m = -1.0d0 

!....Complex amplitude and the produce of i,freq, and t

i1        = SQRT(m)
amplitude = PROB_DATA%amp*EXP(-i1*PROB_DATA%phase)
iwt       = i1*PROB_DATA%freq*solnT
wt        = PROB_DATA%freq*REAL(solnT)

!....Compute normalized depth

h         = dcmplx(hmax, 0.0d0)
zH        = zCOORD/hmax

!....Compute dimensionless slip coefficient and V(z) based on vertical eddy viscosity

N0H0      = difFLUX
N         = N0H0*h**2
KHN       = 100.0d0/3.0d0
k         = (KHN*N)/h

Nc        = dcmplx(N, 0.0d0)
KHNc      = dcmplx(KHN, 0.0d0)

lambda    = sqrt(i1*PROB_DATA%freq*(h**2)/Nc)
VZ        = 1.0d0-CCOSH(lambda*zH)/(CCOSH(lambda)*(1.0d0+lambda/KHNc*CTANH(lambda)))
VBARR      = 1.0d0-(CTANH(lambda))/(lambda*(1.0d0+lambda/KHNc*CTANH(lambda)))

!....Compute constants

tau       = Nc/(h**2)*lambda**2*CTANH(lambda)/(lambda +       &
            (lambda**2/KHNc-1.0d0)*CTANH(lambda))

betaA     = SQRT((PROB_DATA%freq**2 - i1*PROB_DATA%freq*tau)/(PROB_DATA%g*PROB_DATA%Hslope))

!....Compute complex surface elevation zeta and gradient of surface elevation

x1        = dcmplx(PROB_DATA%x0, 0.0)
x2        = dcmplx(PROB_DATA%xN, 0.0)
 
s1       = dcmplx(-(1.0d0/2.0d0),0.0d0) + SQRT(dcmplx(1.0d0/4.0d0,0.0d0) - betaA**2)
s2       = dcmplx(-(1.0d0/2.0d0),0.0d0) - SQRT(dcmplx(1.0d0/4.0d0,0.0d0) - betaA**2)
C1       = s2*x2**(s1)*x1**(s2) - s1*x1**(s1)*x2**(s2)
C2       = (amplitude*s2*x1**(s2))/C1
C3       = -(amplitude*s1*x1**(s1))/C1
CSE      = C2*xCOORD**(s1) + C3*xCOORD**(s2)
cgrad_SE = C2*s1*xCOORD**(s1-1.0d0) + C3*s2*xCOORD**(s2-1.0d0)      

!....Compute vertical velocity

gamma    = (PROB_DATA%g*PROB_DATA%Hslope)/(i1*PROB_DATA%freq)*exp(iwt)
alpha1   = C2*s1*xCOORD**(s1) + C3*s2*xCOORD**(s2)
alpha2   = C2*(s1**2.0d0)*xCOORD**(s1) + C3*(s2**2.0d0)*xCOORD**(s2) 
Delta    = lambda/(CCOSH(lambda)*(dcmplx(1.0d0,0.0d0)+(lambda/KHNc)*CTANH(lambda)))
C11      = (zH*CCOSH(lambda*zH) + CCOSH(lambda))/lambda
C12      = (CSINH(lambda*zH) + CSINH(lambda))/lambda**(2.0d0)
C22      = zH + (Delta*CSINH(lambda*zH) + CSINH(lambda))/lambda**(2.0d0)
C33      = (dcmplx(1.0d0,0.0d0) - (Delta*CCOSH(lambda))/lambda)
VW       = dcmplx(2.0d0,0.0d0)*gamma*alpha1*Delta*(C11 - C12) &
         + gamma*alpha2*C22 + dcmplx(2.0d0,0.0d0)*gamma*alpha1*C33   
 

!....Compute time dependent solutions

V0        = -PROB_DATA%g*cgrad_SE/(i1*PROB_DATA%freq) 

grad_SE   = REAL(cgrad_SE*EXP(iwt))
SE        = REAL(CSE*EXP(iwt))
Vu        = REAL(V0*VZ*EXP(iwt))
Vavg      = REAL(V0*VBARR*EXP(iwt))
VWT       = REAL(VW) - REAL(zH*Vu*2.0d0*xCOORD*PROB_DATA%Hslope)

IF (debug == 6) THEN
    write(*,*)'SE=',SE,'--'
    write(*,*)'grad_SE=',grad_SE,'--'
    write(*,*)'V=',Vu,'--'
    write(*,*)'Vavg=',Vavg,'--'
ENDIF

CONTAINS 

  FUNCTION CTANH(Y)
    COMPLEX(16) :: Y, CTANH
    COMPLEX(16), PARAMETER :: TWO = 2.0d0
    CTANH = (EXP(TWO*Y) - 1)/(EXP(TWO*Y) + 1)
  END FUNCTION CTANH

  FUNCTION CCOSH(Y)
    COMPLEX(16) :: Y, CCOSH
    COMPLEX(16) :: ONE = 1.0000d0, HALF = 0.50000d0	
    CCOSH = HALF*(EXP(Y) + EXP(-Y))
  END FUNCTION
  
  FUNCTION CSINH(Y)
    COMPLEX(16) :: Y, CSINH
    COMPLEX(16) :: ONE = 1.0000d0, HALF = 0.50000d0	
    CSINH = HALF*(EXP(Y) - EXP(-Y))
  END FUNCTION  

END SUBROUTINE LG3D_SOLUTIONS_H
!--------------------------------------------------------------------
!
!                         KERNEL preprocessor
!                     Written by Colton J. Conroy
!                          @ the C.H.I.L
!                              1.28.14
!
!--------------------------------------------------------------------
SUBROUTINE KERNEL_PREPROCESS()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: HK, Ltype, Lx, Ly, Jx, Jy, Qintx, Qinty, Kcount, QQ
INTEGER(int_p) :: stencil_1D, jj, intELEM, startQ, finishQ, Kcount2
INTEGER(int_p) :: check, k, startX, startY, offset
INTEGER(int_p), PARAMETER :: SHIFTED = 11, REGIONs = 1
REAL(real_p) :: xpt, INVJAC

!....1D interior elements

ielems = 0
IF (p == 1) THEN
    startQ = 2
    DO j = 3, nelemsX-2
        ielems = ielems + 1
    ENDDO
    finishQ = nelemsX-1
ELSEIF (p == 2) THEN
    startQ = 4
    DO j = 5, nelemsX-4
        ielems = ielems + 1
    ENDDO
    finishQ = nelemsX-3    
ENDIF

ielems_2D  = ielems**2
stencil_1D = 4*p+1
stencil_2D = (4*p+1)**2

!....Allocate 1D Matricies

ALLOCATE(KC(2**(p+1)+1,p+1,ielems,2,NPTS_K))
ALLOCATE(DKC(2**(p+1)+1,p+1,ielems,2,NPTS_K))  
ALLOCATE(KR(2**(p+1)+1,p+1,2,2,NPTS_K))
ALLOCATE(DKR(2**(p+1)+1,p+1,2,2,NPTS_K))   
ALLOCATE(KL(2**(p+1)+1,p+1,2,2,NPTS_K))
ALLOCATE(DKL(2**(p+1)+1,p+1,2,2,NPTS_K))      
ALLOCATE(KC2(2**(p+1)+1,p+1,ielems,2,NPTS_K))
ALLOCATE(DKC2(2**(p+1)+1,p+1,ielems,2,NPTS_K))  
ALLOCATE(KR2(2**(p+1)+1,p+1,2,2,NPTS_K))
ALLOCATE(DKR2(2**(p+1)+1,p+1,2,2,NPTS_K))   
ALLOCATE(KL2(2**(p+1)+1,p+1,2,2,NPTS_K))
ALLOCATE(DKL2(2**(p+1)+1,p+1,2,2,NPTS_K))        
ALLOCATE(ZETAstar(NPTS_K,ielems))
ALLOCATE(DZETAstar(NPTS_K,ielems))
ALLOCATE(UBARstar(NPTS_K,ielems))
ALLOCATE(DUBARstar(NPTS_K,ielems))
ALLOCATE(shiftPTS(SHIFTint,1))
ALLOCATE(shiftWTS(SHIFTint))   
ALLOCATE(PHIs(SHIFTint,p+1))
ALLOCATE(DPHIs(SHIFTint,p+1))  
ALLOCATE(CC(2**(p+1)+1,p+1))
ALLOCATE(CD(2**(p+1)+1,p+1))
ALLOCATE(ZETA_USE(p+1,nelemsX))
ALLOCATE(zetaENHANCED(nnodesX))
ALLOCATE(zetaENHANCED_2D(NPTS_K2D,Qelems))
ALLOCATE(wENHANCED_2D(NPTS_K2D,Qelems))
ALLOCATE(uavgENHANCED_2D(NPTS_K2D,Qelems))
ALLOCATE(zetaENHANCED2(nnodesX))
ALLOCATE(bPHIs(NbPTS,p+1))
ALLOCATE(bDPHIs(NbPTS,p+1))
ALLOCATE(KC_2D(stencil_2D,NDOF,ielems_2D,(NPTS_K)**2))
ALLOCATE(Qstencil(stencil_2D,Qelems))
ALLOCATE(KCelems(ielems_2D))

Qstencil  = 0
KC        = 0.0d0
DKC       = 0.0d0
KR        = 0.0d0
DKR       = 0.0d0
KL        = 0.0d0
DKL       = 0.0d0
KC2       = 0.0d0
DKC2      = 0.0d0
KR2       = 0.0d0
DKR2      = 0.0d0
KL2       = 0.0d0
DKL2      = 0.0d0
ZETAstar  = 0.0d0
DZETAstar = 0.0d0
UBARstar  = 0.0d0
DUBARstar = 0.0d0
shiftPTS  = 0.0d0
shiftWTS  = 0.0d0
PHIs      = 0.0d0
DPHIs     = 0.0d0
CC        = 0.0d0
CD        = 0.0d0
bPHIs     = 0.0d0
bDPHIs    = 0.0d0
KC_2D     = 0.0d0
zetaENHANCED_2D  = 0.0d0
wENHANCED_2D     = 0.0d0
uavgENHANCED_2D = 0.0d0
 
!....Get shifted points for kernel

CALL QUADRATURE(SHIFTint,SHIFTED,shiftPTS,shiftWTS,SHIFTint)

!....Basis functions at shifted points

DO i = 1, SHIFTint
    CALL ORTHOGONAL_BASIS(SHIFTED, shiftPTS(i,1), p+1, DIM, BASIS, DBASIS)
    DO j = 1, p+1
        PHIs(i,j)  = BASIS(j)
        DPHIs(i,j) = DBASIS(j)
    ENDDO 
ENDDO

bPTS(1) = 0.0d0
bPTS(2) = 1.0d0
DO i = 1, NbPTS
    CALL ORTHOGONAL_BASIS(KERNEL,REGIONs, bPTS(i), p+1, DIM, BASIS, DBASIS)
    DO j = 1, p+1
      bPHIs(i,j) = BASIS(j)
      bDPHIs(i,j) = DBASIS(j)
    ENDDO
ENDDO

!....Kernel coefficients in x-dir.

DO Kcount = 1, NPTS_K

    Ltype = 1

    IF (p == 1 .and. nelemsX >= 5) THEN
        
        DO J = 1, nelemsX
            INVJAC = 1/xELEM_jac(J)
            IF (J < 3) THEN
                HK = -1
            ELSEIF (J > nelemsX-2) THEN
                HK = 1
            ELSE 
                HK = 0
            ENDIF
            xpt = X(J) + INVJAC*(1.0d0 + KPTS_1D(Kcount,1)) 
            CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,X,dx,nelemsX,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
            IF (HK == -1) THEN !....Kernel on left boundary pointing in right direction
                KR(:,:,J,1,Kcount)  = CC(:,:)
                DKR(:,:,J,1,Kcount) = CD(:,:)
            ELSEIF (HK == 0) THEN !....Centered Kernel
                KC(:,:,J-2,1,Kcount)  = CC(:,:)
                DKC(:,:,J-2,1,Kcount) = CD(:,:)
            ELSEIF (HK == 1) THEN !....Kernel on right boundary pointed in left direction
                KL(:,:,J-(nelemsX-2),1,Kcount)  = CC(:,:)
                DKL(:,:,J-(nelemsX-2),1,Kcount) = CD(:,:)
            ENDIF           
        ENDDO
              
    ELSEIF (p == 2 .and. nelemsX >= 9) THEN
        HK = 0
        DO J = 5, nelemsX-4
            INVJAC = 1/xELEM_jac(J)
            xpt    = X(J) + INVJAC*(1.0d0 + KPTS_1D(Kcount,1)) 
            CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,X,dx,nelemsX,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
            KC(:,:,J-4,1,Kcount)  = CC(:,:)
            DKC(:,:,J-4,1,Kcount) = CD(:,:)
        ENDDO
    ELSE       
        PRINT*,'  ******** WARNING!!! # of elements   ********  '
        PRINT*,'  ******** must be >= 2**(p+1)+1 to    ********  '
        PRINT*,'  ******** use the kernel, execution  ********  '
        PRINT*,'  ******** will continue w/out the    ********  '
        PRINT*,'  ******** use of the kernel.         ********  '
        KERNEL = 0
        GOTO 10
    ENDIF
    
!    IF (p == 1 .and. nelemsX >= 5) THEN
!        
!        DO J = 1, nelemsX
!            IF (J < 3) THEN
!                HK = -1
!            ELSEIF (J > nelemsX-2) THEN
!                HK = 1
!            ELSE 
!                HK = 0
!            ENDIF
!            xpt = X(J+1) - (1.0d0 + KPTS_1D(Kcount,1)) 
!            CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,X,dx,nelemsX,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
!            IF (HK == -1) THEN !....Kernel on left boundary pointing in right direction
!                KR2(:,:,J,1,Kcount)  = CC(:,:)
!                DKR2(:,:,J,1) = CD(:,:)
!            ELSEIF (HK == 0) THEN !....Centered Kernel
!                KC2(:,:,J-2,1)  = CC(:,:)
!                DKC2(:,:,J-2,1) = CD(:,:)
!            ELSEIF (HK == 1) THEN !....Kernel on right boundary pointed in left direction
!                KL2(:,:,J-(nelemsX-2),1)  = CC(:,:)
!                DKL2(:,:,J-(nelemsX-2),1) = CD(:,:)
!            ENDIF           
!        ENDDO
!              
!    ELSEIF (p == 2 .and. nelemsX >= 9) THEN
!        HK = 0
!        DO J = 5, nelemsX-4
!            xpt = X(J) 
!            CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,X,dx,nelemsX,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
!            KC2(:,:,J-4,1)  = CC(:,:)
!            DKC2(:,:,J-4,1) = CD(:,:)
!        ENDDO
!    ENDIF 

!....Kernel coefficients in y-dir.

    Ltype = 2

    IF (p == 1 .and. nelemsY >= 5) THEN
        
        DO J = 1, nelemsY
            INVJAC = 1/yELEM_jac(J)
            IF (J < 3) THEN
                HK = -1
            ELSEIF (J > nelemsY-2) THEN
                HK = 1
            ELSE 
                HK = 0
            ENDIF
            xpt    = Y(J) + INVJAC*(1.0d0 + KPTS_1D(Kcount,1)) 
            CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,Y,dy,nelemsY,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
            IF (HK == -1) THEN !....Kernel on left boundary pointing in right direction
                KR(:,:,J,2,Kcount)  = CC(:,:)
                DKR(:,:,J,2,Kcount) = CD(:,:)
            ELSEIF (HK == 0) THEN !....Centered Kernel
                KC(:,:,J-2,2,Kcount)  = CC(:,:)
                DKC(:,:,J-2,2,Kcount) = CD(:,:)
            ELSEIF (HK == 1) THEN !....Kernel on right boundary pointed in left direction
                KL(:,:,J-(nelemsX-2),2,Kcount)  = CC(:,:)
                DKL(:,:,J-(nelemsX-2),2,Kcount) = CD(:,:)
            ENDIF           
        ENDDO
              
    ELSEIF (p == 2 .and. nelemsY >= 9) THEN
        HK = 0
        DO J = 5, nelemsY-4
            INVJAC = 1/yELEM_jac(J)
            xpt    = Y(J) + INVJAC*(1.0d0 + KPTS_1D(Kcount,1)) 
            CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,Y,dy,nelemsY,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
            KC(:,:,J-4,2,Kcount)  = CC(:,:)
            DKC(:,:,J-4,2,Kcount) = CD(:,:)
        ENDDO
    ELSE       
        PRINT*,'  ******** WARNING!!! # of elements   ********  '
        PRINT*,'  ******** must be >= 2**(p+1)+1 to    ********  '
        PRINT*,'  ******** use the kernel, execution  ********  '
        PRINT*,'  ******** will continue w/out the    ********  '
        PRINT*,'  ******** use of the kernel.         ********  '
        KERNEL = 0
        GOTO 10
    ENDIF
    
!IF (p == 1 .and. nelemsY >= 5) THEN
!        
!    DO J = 1, nelemsY
!        IF (J < 3) THEN
!            HK = -1
!        ELSEIF (J > nelemsY-2) THEN
!            HK = 1
!        ELSE 
!            HK = 0
!        ENDIF
!        xpt    = Y(J+1) 
!        CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,Y,dy,nelemsY,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
!        IF (HK == -1) THEN !....Kernel on left boundary pointing in right direction
!            KR2(:,:,J,2)  = CC(:,:)
!            DKR2(:,:,J,2) = CD(:,:)
!        ELSEIF (HK == 0) THEN !....Centered Kernel
!            KC2(:,:,J-2,2)  = CC(:,:)
!            DKC2(:,:,J-2,2) = CD(:,:)
!        ELSEIF (HK == 1) THEN !....Kernel on right boundary pointed in left direction
!            KL2(:,:,J-(nelemsX-2),2)  = CC(:,:)
!            DKL2(:,:,J-(nelemsX-2),2) = CD(:,:)
!        ENDIF           
!    ENDDO
!              
!ELSEIF (p == 2 .and. nelemsY >= 9) THEN
!    HK = 0
!    DO J = 5, nelemsY-4
!        xpt = Y(J) 
!        CALL KERNEL_COEFFICIENTS(HK,xpt,J,p,Y,dy,nelemsY,PHIs,DPHIs,bPHIs,SHIFTint,shiftPTS,shiftWTS,CC,CD,Ltype)
!        KC2(:,:,J-4,2)  = CC(:,:)
!        DKC2(:,:,J-4,2) = CD(:,:)
!    ENDDO
!ENDIF

ENDDO

!....Tensor Product of coefficients (Currently just centered Kernel)

    IF (p == 1) THEN
        offset = 2
    ELSEIF (p == 2) THEN
        offset = 4
    ENDIF
QQ = 0    
DO Kcount = 1, NPTS_K
    DO Kcount2 = 1, NPTS_K 
        QQ = QQ + 1   
        jj = 0
        DO i = 1, Qelems
            Qintx = CONNielem(i,1)
            Qinty = CONNielem(i,2)
            IF (QintX > startQ .and. QintX < finishQ) THEN
                IF (QintY > startQ .and. QintY < finishQ) THEN
                    jj = jj + 1
                    kk = 0
                    DO Ly = 1, p+1
                        DO Lx = 1, p+1
                            ii = 0
                            kk = kk + 1
                            DO Jx = 1, stencil_1D
                                DO Jy = 1, stencil_1D
                                    ii = ii + 1
                                    KC_2D(ii,kk,jj,QQ) = KC(Jx,Lx,Qintx-offset,1,Kcount2)*KC(Jy,Ly,Qinty-offset,2,Kcount)
                                ENDDO
                            ENDDO
                        ENDDO             
                    ENDDO
                ENDIF    
            ENDIF
        ENDDO   
    ENDDO            
ENDDO

!....Element stencil connectivity

intELEM = 0
DO kk = 1, Qelems
    Qintx = CONNielem(kk,1)
    Qinty = CONNielem(kk,2)
    IF (QintX > startQ .and. QintX < finishQ) THEN
        IF (QintY > startQ .and. QintY < finishQ) THEN
            intELEM = intELEM + 1
            KCelems(intELEM) = kk
            startX = QintX - offset
            startY = QintY - offset
            k = 0
            DO i = 0, stencil_1D - 1
                DO j = 0, stencil_1D - 1
                    k = k + 1
                    ii = startX + i
                    jj = startY + j
                    Qstencil(k,kk) = CONNgrid(jj,ii)
                ENDDO
            ENDDO
        ENDIF        
    ENDIF                
ENDDO            

        
10 RETURN
END SUBROUTINE KERNEL_PREPROCESS
!---------------------------------------------------------------------
!
!             Hyperbolic Kernel Coefficients Subroutine
!                     Written by Colton J. Conroy
!                           @ the C.H.I.L
!                               2.20.13
!
!---------------------------------------------------------------------
SUBROUTINE KERNEL_COEFFICIENTS(HK,xi,J,p,X,dx,nelems,PHI,DPHI,bPHI,SHIFTint,shiftPTS,shiftWTS,CC,CD)
         
USE precisions
IMPLICIT NONE

INTEGER(int_p) :: L, i, ii, J, jj, k, kk, p, SHIFTint, start1, finish1, start2, finish2
INTEGER(int_p) :: nelems, pp, HK
REAL(real_p) :: dx, sum1, sum2, sum1d, sum2d, BS, DBS, point, xi
REAL(real_p), DIMENSION(4*p+1,p+1) :: CC, CD
REAL(real_p), DIMENSION(nelems+1) :: X
REAL(real_p), DIMENSION(SHIFTint,1) :: shiftPTS
REAL(real_p), DIMENSION(SHIFTint) :: shiftWTS
REAL(real_p), DIMENSION(SHIFTint,p+1) :: PHI, DPHI
REAL(real_p), DIMENSION(2,p+1) :: bPHI
REAL(real_p), ALLOCATABLE :: KERN1(:,:,:), KERN2(:,:,:), G(:), C(:)

!....Allocate and initialize matricies based on degree of Legendre polynomial

CC      =  0.0d0
CD      =  0.0d0

IF (p == 1) THEN

    ALLOCATE(G(3))
    ALLOCATE(C(3))
    ALLOCATE(KERN1(SHIFTint,p+4,3))
    ALLOCATE(KERN2(SHIFTint,p+4,3))
    KERN1 = 0.0d0
    KERN2 = 0.0d0
    
    !....Type of Kernel
    
    IF (HK == 0) THEN  !....Centered Kernel
    
    !....Gamma matrix
    
        G(1)    = -1.0d0
        G(2)    =  0.0d0
        G(3)    =  1.0d0
        
    
    
    !....Indexing
    
        start1  =  J - 2
        finish1 =  J + 2
        start2  = -1
        finish2 =  1
    
    !....Coefficients for p = 1
    
        C(1)    =  -(1.0d0/12.0d0)
        C(2)    =  7.0d0/6.0d0
        C(3)    =  C(1)
         
    ELSEIF (HK == -1) THEN  !....Kernel for left boundary
    
    !....Gamma matrix
    
        G(1)    =  3.0d0
        G(2)    =  2.0d0
        G(3)    =  1.0d0
        
    
    
    !....Indexing
    
        start1  =  J - 0
        finish1 =  J + 4
        start2  = -1
        finish2 =  1
    
    !....Coefficients for p = 1
    
        C(1)    =   11.0d0/12.0d0
        C(2)    =  -17.0d0/6.0d0
        C(3)    =   35.0d0/12.0d0    
    
    ELSEIF (HK == 1) THEN !....Kernel for right boundary
    
    !....Gamma matrix
    
        G(1)    =  -3.0d0
        G(2)    =  -2.0d0
        G(3)    =  -1.0d0
        
    
    
    !....Indexing
    
        start1  =  J - 4
        finish1 =  J - 0
        start2  = -1
        finish2 =  1
    
    !....Coefficients for p = 1
    
        C(1)    =   11.0d0/12.0d0
        C(2)    =  -17.0d0/6.0d0
        C(3)    =   35.0d0/12.0d0            
   
   ENDIF
        
    pp = p
    
ELSEIF (p == 2) THEN
    
    ALLOCATE(G(5))
    ALLOCATE(C(5))
    ALLOCATE(KERN1(SHIFTint,p+7,5))
    ALLOCATE(KERN2(SHIFTint,p+7,5))  
    KERN1 = 0.0d0
    KERN2 = 0.0d0  
    
    !....Gamma matrix
    
    G(1) = -2.0d0
    G(2) = -1.0d0
    G(3) =  0.0d0
    G(4) =  1.0d0
    G(5) =  2.0d0
    
    !....Indexing
    
    start1  = J - 4
    finish1 = J + 4
    start2  = -2
    finish2 =  2
    
    !....Coefficients for p = 2
    
    C(1) =  37.0d0/1920.0d0
    C(2) =  -97.0d0/480.0d0
    C(3) = 437.0d0/320.0d0
    C(4) = C(2) 
    C(5) = C(1) 
    
    pp = p
        
ENDIF

IF (p == 1) THEN

!....Evaluate B-splines at integration points

    kk = 0
    DO i = start1, finish1  !....J-(p+1):J+(p+1)
        kk = kk + 1         !....counter for i loop
        DO ii = start2, finish2  !....-p:p
            DO jj = 1, SHIFTint
                point = ((shiftPTS(jj,1)*dx+X(i))-xi)/dx - G(ii+2)  
                CALL B_SPLINE(pp,point,BS,DBS)
                KERN1(jj,kk,ii+2) = BS                 !....B-spline evaluated at point
                CALL B_SPLINE(pp,point,BS,DBS)
                KERN2(jj,kk,ii+2) = BS                 !....Derivative of B-spline
            ENDDO    
        ENDDO
    ENDDO

!....Evaluate Coefficients

    DO L = 0, p                   !....loop over degrees of freedom
        kk = 0
        DO i = start1, finish1    !....J-(p+1):J+(p+1)
            kk = kk + 1           !....counter for i loop
            sum2  = 0.0d0
            sum2d = 0.0d0
            DO jj = start2, finish2    !....-p:p
                sum1  = 0.0d0
                sum1d = 0.0d0
                DO ii = 1, SHIFTint   
                    sum1  = sum1  + shiftWTS(ii)*KERN1(ii,kk,jj+2)*PHI(ii,L+1)
                    !sum1d = sum1d + shiftWTS(ii)*KERN2(ii,kk,jj+2)*DPHI(ii,L+1)
                    sum1d = sum1d + shiftWTS(ii)*KERN2(ii,kk,jj+2)*(bPHI(2,L+1)-bPHI(1,L+1))
                ENDDO                
                sum2  = sum2  + C(jj+2)*sum1
                sum2d = sum2d + C(jj+2)*sum1d
            ENDDO
            CC(kk,L+1) = sum2
            CD(kk,L+1) = sum2d
        ENDDO
    ENDDO

ELSEIF (p == 2) THEN

    kk = 0
    DO i = start1, finish1  !....J-(p+1):J+(p+1)
        kk = kk + 1         !....counter for i loop
        DO ii = start2, finish2  !....-p:p
            DO jj = 1, SHIFTint
                point = ((shiftPTS(jj,1)*dx+X(i))-xi)/dx - G(ii+3)  
                CALL B_SPLINE(pp,point,BS,DBS)
                KERN1(jj,kk,ii+3) = BS                 !....B-spline evaluated at point
                CALL B_SPLINE(pp,point,BS,DBS)
                KERN2(jj,kk,ii+3) = BS                 !....Derivative of B-spline
            ENDDO    
        ENDDO
    ENDDO

!....Evaluate Coefficients

    DO L = 0, p                   !....loop over degrees of freedom
        kk = 0
        DO i = start1, finish1    !....J-(p+1):J+(p+1)
            kk = kk + 1           !....counter for i loop
            sum2  = 0.0d0
            sum2d = 0.0d0
            DO jj = start2, finish2    !....-p:p
                sum1  = 0.0d0
                sum1d = 0.0d0
                DO ii = 1, SHIFTint   
                    sum1  = sum1  + shiftWTS(ii)*KERN1(ii,kk,jj+3)*PHI(ii,L+1)
                    sum1d = sum1d + shiftWTS(ii)*KERN2(ii,kk,jj+3)*DPHI(ii,L+1)
                ENDDO                
                sum2  = sum2  + C(jj+3)*sum1
                sum2d = sum2d + C(jj+3)*sum1d
            ENDDO
            CC(kk,L+1) = sum2
            CD(kk,L+1) = sum2d
        ENDDO
    ENDDO

ENDIF

DEALLOCATE(G)
DEALLOCATE(C)
DEALLOCATE(KERN1)
DEALLOCATE(KERN2)

RETURN
END SUBROUTINE KERNEL_COEFFICIENTS
!-------------------------------------------------------------------
!
!                 2D Hyperbolic Kernel Subroutine
!                    Written by Colton J. Conroy
!                          @ the C.H.I.L
!                              1.29.14
!
!--------------------------------------------------------------------
SUBROUTINE HYPERBOLIC_KERNEL_2D(KC,U,p,Qelems,ielems_2D,Ustar,stencil_2D,NDOF,Qstencil,nelemsX,nelemsY,CONNielem,NPTS_K)

USE precisions
IMPLICIT NONE

INTEGER(int_p) :: i, j, k, Q, L, Qelems, ielems_2D, p, nelemsY
INTEGER(int_p) :: stencil_2D, NDOF, QintX, QintY, nelemsX
INTEGER(int_p) :: startQ, finishQ, NPTS_K, Kcount, finishK
INTEGER(int_p), DIMENSION(stencil_2D,Qelems) :: Qstencil
INTEGER(int_p), DIMENSION(Qelems,2) :: CONNielem
REAL(real_p) :: sum1, sum2
REAL(real_p), DIMENSION(stencil_2D,NDOF,ielems_2D,NPTS_K) :: KC
REAL(real_p), DIMENSION((NPTS_K)**2,Qelems) :: Ustar
REAL(real_p), DIMENSION(NDOF,Qelems) :: U

IF (p == 1) THEN
    startQ = 2
    finishQ = nelemsX-1
ELSEIF (p == 2) THEN
    startQ = 4
    finishQ = nelemsX-3    
ENDIF

finishK = NPTS_K**2
Ustar = 0.0d0

DO Kcount = 1, finishK
    k = 0
    DO i = 1, Qelems
        Qintx = CONNielem(i,1)
        Qinty = CONNielem(i,2)
        IF (QintX > startQ .and. QintX < finishQ) THEN
            IF (QintY > startQ .and. QintY < finishQ) THEN
                k = k + 1
                sum2 = 0.0d0   
                DO j = 1, stencil_2D
                    Q = Qstencil(j,i)
                    sum1 = 0.0d0
                    DO L = 1, NDOF
                        sum1 = sum1 + U(L,Q)*KC(j,L,k,Kcount)
                    ENDDO   
                    sum2 = sum2 + sum1
                ENDDO    
                Ustar(Kcount,i) = sum2
            ENDIF
        ENDIF
    ENDDO
ENDDO 
        
RETURN
END SUBROUTINE HYPERBOLIC_KERNEL_2D
!--------------------------------------------------------------------
!
!                           Check KERNEL Error
!                       Written by Colton J. Conroy
!                             @ the C.H.I.L
!                                 1.29.14
!
!--------------------------------------------------------------------
SUBROUTINE KERNEL_2D_ERROR(solnT)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: Qe, Qk1, Qk2, jj, finishK
REAL(real_p) :: xpt, ypt, zpt, Vu, Vv, Vw, DVw, Huse, Zb, L2SE, L2Ub
REAL(real_p) :: L2sum1, L2sum2, solnT, SE, DSEx, DSEy, Uavg, L2sum3, L2sum4
REAL(real_p) :: Vavg, Fpce, Fxm, Fym, Fw, nn, k, grad_SE, FxmA, FymA
REAL(real_p), DIMENSION((NPTS_K)**2) :: L2DG, L2DG2
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Ynode
REAL(real_p), DIMENSION((NPTS_K)**2,ielems_2D) :: LinfS, LinfDG, LinfS2, LinfDG2
COMPLEX(16) :: xCOORD, yCOORD, zCOORD, TIME

TIME    = dcmplx(solnT,0.0d0)
Xnode   = 0.0d0
Ynode   = 0.0d0
LinfS   = 0.0d0
LinfDG  = 0.0d0
LinfS2  = 0.0d0
LinfDG2 = 0.0d0
L2DG    = 0.0d0
L2DG2   = 0.0d0
globalNODES = 0.0d0

IF (var == 1) THEN !....Surface elevation

    DO i = 1, ielems_2D
        jj = 0
        Qe = KCelems(i)
        L2DG  = MATMUL(PHIkernel,ZETA(:,Qe,1))
        L2DG2 = MATMUL(PHIkernel,USE(:,Qe,1))
        DO Qk1 = 1, NPTS_K
            DO Qk2 = 1, NPTS_K
                jj = jj + 1
                globalNODES = CONN(Qe,:)
                Xnode(1) = Qnode(globalNODES(1),1)
                Xnode(2) = Qnode(globalNODES(3),1)
                Ynode(1) = Qnode(globalNODES(1),2)
                Ynode(2) = Qnode(globalNODES(2),2)
                xpt = (1.0d0/2.0d0)*((1.0d0-KPTS_1D(Qk2,1))*Xnode(1)+(1.0d0+KPTS_1D(Qk2,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-KPTS_1D(Qk1,1))*Ynode(1)+(1.0d0+KPTS_1D(Qk1,1))*Ynode(2))       
                zpt = nodesZRK(Z_level)
                IF (PROB == 1) THEN
                    Huse = PROB_DATA%Hslope*xpt**2
                    zCOORD = dcmplx(Huse*zpt, 0.0d0)
                    xCOORD = dcmplx(xpt, 0.0d0)
                    yCOORD = dcmplx(ypt, 0.0d0)
                    CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
                    L2SE = SE
                    L2Ub = Vavg
                ELSEIF (PROB == 2) THEN
                    CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                    L2SE = Huse
                ELSEIF (PROB == 3) THEN
                    CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                    L2SE = Huse
                    L2Ub = Uavg
                ENDIF          
                LinfS(jj,i)   = abs(L2SE-zetaENHANCED_2D(jj,Qe))
                LinfDG(jj,i)  = abs(L2SE-L2DG(jj))
                LinfS2(jj,i)  = abs(L2Ub-uavgENHANCED_2D(jj,Qe))
                LinfDG2(jj,i) = abs(L2Ub-L2DG2(jj)) 
            ENDDO
        ENDDO    
    ENDDO
    finishK = NPTS_K*NPTS_K
    L2sum1 = 0.0d0
    L2sum2 = 0.0d0
    L2sum3 = 0.0d0
    L2sum4 = 0.0d0
    DO i = 1, finishK
        L2sum1 = L2sum1 + MAXVAL(LinfS(i,:))
        L2sum2 = L2sum2 + MAXVAL(LinfDG(i,:))  
        L2sum3 = L2sum3 + MAXVAL(LinfS2(i,:))
        L2sum4 = L2sum4 + MAXVAL(LinfDG2(i,:)) 
    ENDDO   
    write(*,*)'ZETA* error =', L2sum1/finishK
    write(*,*)'ZETAdg error =', L2sum2/finishK
    write(*,*)'Uavg* error =', L2sum3/finishK
    write(*,*)'UavgDG error =', L2sum4/finishK

ELSEIF (var == 2) THEN !....Vertical velocity

    DO i = 1, ielems_2D
        jj = 0
        Qe = KCelems(i)
        L2DG  = MATMUL(PHIkernel,W_2D(:,Qe))
        DO Qk1 = 1, NPTS_K
            DO Qk2 = 1, NPTS_K
                jj = jj + 1
                globalNODES = CONN(Qe,:)
                Xnode(1) = Qnode(globalNODES(1),1)
                Xnode(2) = Qnode(globalNODES(3),1)
                Ynode(1) = Qnode(globalNODES(1),2)
                Ynode(2) = Qnode(globalNODES(2),2)
                xpt = (1.0d0/2.0d0)*((1.0d0-KPTS_1D(Qk2,1))*Xnode(1)+(1.0d0+KPTS_1D(Qk2,1))*Xnode(2))
                ypt = (1.0d0/2.0d0)*((1.0d0-KPTS_1D(Qk1,1))*Ynode(1)+(1.0d0+KPTS_1D(Qk1,1))*Ynode(2))       
                zpt = nodesZRK(Z_level)
                IF (PROB == 2) THEN
                    CALL Vertical_Test(xpt,ypt,zpt,Vu,Vv,Vw,DVw,Huse,Zb)
                ELSEIF (PROB == 3) THEN
                    CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
                ENDIF    
                L2SE         = Vw
                LinfS(jj,i)  = abs(Vw-wENHANCED_2D(jj,Qe))
                LinfDG(jj,i) = abs(Vw-L2DG(jj))
            ENDDO
        ENDDO    
    ENDDO
    finishK = NPTS_K*NPTS_K
    L2sum1 = 0.0d0
    L2sum2 = 0.0d0
    DO i = 1, finishK
        L2sum1 = L2sum1 + MAXVAL(LinfS(i,:))
        L2sum2 = L2sum2 + MAXVAL(LinfDG(i,:))  
    ENDDO   
    IF (PROB == 2 .OR. PROB == 3) THEN
        write(*,*)'W* error =', L2sum1/finishK
        write(*,*)'Wdg error=', L2sum2/finishK 
    ENDIF         
        
ENDIF    
    
RETURN
END SUBROUTINE KERNEL_2D_ERROR
!--------------------------------------------------------------------
!
!                 1D Hyperbolic Kernel Subroutine
!                    Written by Colton J. Conroy
!                           @ the C.H.I.L
!                               2.20.13
!
!--------------------------------------------------------------------
SUBROUTINE HYPERBOLIC_KERNEL_1D(HK,KC,DKC,U,p,nelems,NPTS,ielems,Ustar,DUstar,dx,Ltype)

USE precisions
IMPLICIT NONE

INTEGER(int_p) :: i, j, jj, ii, kk, L, nelems, p, NPTS, start, finish, ielems
INTEGER(int_p) :: start1, finish1, HK, Ltype
REAL(real_p) :: sum1, sum2, sum1d, sum2d, dx
REAL(real_p), DIMENSION(p+1,nelems) :: U
REAL(real_p), DIMENSION(NPTS,ielems) :: Ustar, DUstar
REAL(real_p), DIMENSION(2**(p+1)+1,p+1,NPTS,ielems) :: KC, DKC

IF (Ltype == 1) THEN
    kk = 1
ELSEIF (Ltype == 2) THEN
    kk = 2
ENDIF

IF (p == 1) THEN
    IF (HK == 0) THEN
        start1  = 3
        finish1 = nelems-2 
    ELSEIF (HK == -1) THEN
        start1  = 1
        finish1 = 2 
    ELSEIF (HK == 1) THEN
        start1  = nelems-1
        finish1 = nelems
    ENDIF    
ELSEIF (p == 2) THEN
    start1  = 5
    finish1 = nelems-4
ENDIF
j = 0
DO i = start1, finish1
    j = j + 1
    IF (p == 1) THEN
        IF (HK == 0) THEN
            start  = i - 2
            finish = i + 2    
        ELSEIF (HK == -1) THEN
            start  = i - 0
            finish = i + 4    
        ELSEIF (HK == 1) THEN
            start   = i - 4 
            finish  = i - 0
        ENDIF    
    ELSEIF (p == 2) THEN
        start  = i - 4
        finish = i + 4    
    ENDIF
        ii = 0
        sum2  = 0.0d0
        sum2d = 0.0d0
        DO jj = start, finish
            ii = ii + 1
            sum1  = 0.0d0
            sum1d = 0.0d0
            DO L = 0, p
                sum1  = sum1 + U(L+1,jj)*KC(ii,L+1,j,kk)
                sum1d = sum1d + U(L+1,jj)*DKC(ii,L+1,j,kk)*(1.0d0/dx)
            ENDDO
            sum2  = sum2 + sum1
            sum2d = sum2d + sum1d
        ENDDO

        Ustar(1,j)  = sum2
        DUstar(1,j) = sum2d  
ENDDO

RETURN
END SUBROUTINE HYPERBOLIC_KERNEL_1D
!---------------------------------------------------------------------
!
!                        B-spline Subroutine
!                    Written by Colton J. Conroy
!                           @ the C.H.I.L
!                               2.19.13
!
!----------------------------------------------------------------------
SUBROUTINE B_SPLINE(p,Xi,BS,DBS)

USE precisions
IMPLICIT NONE

INTEGER(int_p) :: i, j, k, p, nknots, nintervals, extra
REAL(real_p) :: dx, term1, term2, term1d, term2d, x0, xN, temp, Xi
REAL(real_p) :: BS, DBS
REAL(real_p), ALLOCATABLE :: x(:), B(:,:), BD(:,:)

!....Constants

x0         = -(p+1)*(1.0d0/2.0d0)             !....initial knot
xN         = -x0                              !....end of interval
nknots     = p + 2                            !....number of knots
nintervals = nknots-1                         !....number of intervals
dx         = (xN-x0)/nintervals               !....spacing of knots

IF (p > 2) THEN                               !....extra knots for p > 2
    extra = p - 2
    nknots = nknots + extra
ENDIF

ALLOCATE(x(nknots))
ALLOCATE(B(p+1,p+1))
ALLOCATE(BD(p,p))

x  = 0.0d0
B  = 0.0d0
BD = 0.0d0

!....knots (constant spacing)

x(1) = x0
DO j = 2, nknots
    x(j) = x(j-1) + dx    
ENDDO 

!....constant B-splines

DO k = 0, p
    B(k+1,1) = B0(Xi,x(k+1),x(k+2))
ENDDO

!....recursion relation from The Art of Scientific Computing pp. 736

i = p
DO k = 1, p
    DO j = 1, i
        term1    = (Xi-x(j))/(x(j+k)-x(j))
        term2    = (x(j+k+1)-Xi)/(x(j+k+1)-x(j+1))
        B(j,k+1) = term1*B(j,k) + term2*B(j+1,k)
        !....derivative
        term1d   = 1.0d0/(x(j+k)-x(j))
        term2d   = 1.0d0/(x(j+k+1)-x(j+1))
        BD(j,k)  = k*(term1d*B(j,k) - term2d*B(j+1,k))
    ENDDO
    i = p - 1
ENDDO

IF (p >= 1) THEN
    BS  = B(1,p+1)
    DBS = BD(1,p)
ELSE
    BS  = B(1,p+1)
    DBS = 0    
ENDIF

CONTAINS

    FUNCTION B0(a,x1,x2)
    REAL(real_p) :: a, x1, x2, B0
        IF (a >= x1 .and. a < x2) THEN
            B0 = 1.0d0
        ELSE
            B0 = 0.0d0
        ENDIF
    END FUNCTION

END SUBROUTINE B_SPLINE
!--------------------------------------------------------------------
!
!                          Output for Video 
!                     Written by Colton J. Conroy
!                           @ the C.H.I.L
!                               3.21.13
!
!---------------------------------------------------------------------
SUBROUTINE VIDEO(n,REGION,nframes,TIME)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: n, REGION, nframes
REAL(real_p) :: TIME
CHARACTER(50) ::  xstring, estring, ystring

IF (plot_n == n) THEN
    nframes = nframes + 1
    !IF (PROB == 1) THEN
    
        IF (REGION == 1 .or. REGION == 2) THEN
            WRITE(estring,'(i10)')plot_n+1
            OPEN(10,FILE=estring)
            DO i = 1, nelemsX
                WRITE(10,*)ZETA(:,i,1)
            ENDDO
            CLOSE(10)
        ELSEIF (REGION == 3) THEN
            WRITE(estring,'(i10)')plot_n+1
            OPEN(10,FILE=estring)
            DO i = 1, Qelems
                WRITE(10,*)ZETA(:,i,1)
            ENDDO
            CLOSE(10)
        ENDIF    

    !ENDIF
    
    WRITE(xstring,'(i10)')plot_n
    WRITE(ystring,'(i10)')plot_n+2

    IF (REGION == 1) THEN
        OPEN(11,FILE=xstring)
        DO i = 1, nelemsX
            DO j = 1, LE
                DO kk = 1, nelemsZ
                    WRITE(11,*)U(:,kk,i,j,1)
                ENDDO
            ENDDO
        ENDDO
        CLOSE(11)
    ELSEIF (REGION == 2) THEN
        OPEN(11,FILE=xstring)
        DO i = 1, Qelems
            WRITE(11,*)U_2D(:,i,1)
        ENDDO
        CLOSE(11)
    ELSEIF (REGION == 3) THEN
        OPEN(11,FILE=xstring)
        DO i = 1, HEXelems
            WRITE(11,*)U_3D(:,i,1)
        ENDDO
        CLOSE(11)   
        !....V-component 
        OPEN(12,FILE=ystring)
        DO i = 1, HEXelems
            WRITE(12,*)V_3D(:,i,1)
        ENDDO
        CLOSE(12)        
    ENDIF
    !write(*,*)TIME
    WRITE(17,*)TIME
    plot_n = plot_n + FLOOR(450/dt)
    !plot_n = plot_n + FLOOR(225/dt)
ENDIF

RETURN
END SUBROUTINE VIDEO
!---------------------------------------------------------------------
!
!                DG W Matricies Subroutine
!              Written by Colton J. Conroy
!                    @ the C.H.I.L
!                        12.17.13
!
!---------------------------------------------------------------------
SUBROUTINE DG_W_MATRICIES(p,pz,WTS1,WTS2,WTS3,WTS4,PTS1,PTS2,PTS3,PTS4,NDOF,NDOF2,NDOF3,L2K,L3,L3PHI,PHIk,PHI,DPHI,bPHI,A,B,KERNEL)

USE precisions
USE global, ONLY: FULL3D, FULL2D
IMPLICIT NONE

INTEGER(int_p) :: i, j, k, ii, p, PTS1, NDOF, NDOF2
INTEGER(int_p) :: PTS2, PTS3, PTS4, pz, NDOF3, KERNEL
REAL(real_p) :: Mx, My, Mz
REAL(real_p), DIMENSION(PTS1) :: WTS1
REAL(real_p), DIMENSION(PTS2) :: WTS2
REAL(real_p), DIMENSION(PTS3) :: WTS3
REAL(real_p), DIMENSION(PTS4) :: WTS4
REAL(real_p), DIMENSION(NDOF,NDOF) :: M
REAL(real_p), DIMENSION(NDOF,PTS1) :: L3, I3
REAL(real_p), DIMENSION(PTS1,NDOF) :: L3PHI
REAL(real_p), DIMENSION(NDOF3,PTS4) :: L2K, IL2
REAL(real_p), DIMENSION(PTS4,NDOF3) :: PHIk
REAL(real_p), DIMENSION(NDOF2,NDOF2) :: M2
REAL(real_p), DIMENSION(NDOF3,NDOF3) :: M3
REAL(real_p), DIMENSION(NDOF2,PTS2,2) :: A, IA
REAL(real_p), DIMENSION(NDOF2,PTS3,4) :: B, IB
REAL(real_p), DIMENSION(PTS2,NDOF2) :: PHI
REAL(real_p), DIMENSION(PTS2,NDOF2,2) :: DPHI
REAL(real_p), DIMENSION(PTS3,NDOF2,4) :: bPHI


!....Initialize matricies

M   = 0.0d0
I3  = 0.0d0
L3  = 0.0d0

!....Mass Matrix

IF (FULL3D == 1) THEN
    ii = 0
    DO i = 0, p
        DO j = 0, p
            DO k = 0, p
                ii = ii + 1
                Mx = (2.0d0*j+1.0d0)/2.0d0       !....Inverse of Mass Matrix
                My = (2.0d0*i+1.0d0)/2.0d0
                Mz = (2.0d0*k+1.0d0)/2.0d0
                M(ii,ii) = Mx*My*Mz
            ENDDO    
        ENDDO
    ENDDO
ELSEIF (FULL3D == 0) THEN
    DO i = 1, NDOF
        DO j = 1, NDOF
            DO k = 1, PTS1
                M(i,j) = M(i,j) + WTS1(k)*L3PHI(k,i)*L3PHI(k,j)
            ENDDO
            IF (i == j) THEN
                M(i,j) = 1.0d0/M(i,j)
            ENDIF    
        ENDDO
    ENDDO
ENDIF 

!....L3 (Map to vertical to horiztonal) integration

DO i = 1, PTS1
    DO j = 1, NDOF
	    I3(j,i) = WTS1(i)*L3PHI(i,j)
    ENDDO
ENDDO

!....L3 matrix

L3 = MATMUL(M,I3)

!....Initialize matricies

M2  = 0.0d0
IA  = 0.0d0
IB  = 0.0d0
A   = 0.0d0
B   = 0.0d0

!....Mass Matrix
k = 0
DO i = 0, p
    DO j = 0, p
        k = k+1
        Mx = (2.0d0*j+1.0d0)/2.0d0       !....Inverse of Mass Matrix
        My = (2.0d0*i+1.0d0)/2.0d0
        M2(k,k) = Mx*My
    ENDDO
ENDDO

!....Advection Matrix

DO i = 1, PTS2
    DO j = 1, NDOF2     
       IA(j,i,1) = WTS2(i)*DPHI(i,j,1)   !....x derivative
       IA(j,i,2) = WTS2(i)*DPHI(i,j,2)   !....y derivative        
    ENDDO
ENDDO

A(:,:,1) = MATMUL(M2,IA(:,:,1))          !....x advection matrix
A(:,:,2) = MATMUL(M2,IA(:,:,2))          !....z advection matrix

!....Boundary Matrix

DO i = 1, PTS3
    DO j = 1, NDOF2
       IB(j,i,1) = WTS3(i)*bPHI(i,j,1)   !....bottom boundary
       IB(j,i,2) = WTS3(i)*bPHI(i,j,2)   !....top boundary
       IB(j,i,3) = WTS3(i)*bPHI(i,j,3)   !....Left boundary
       IB(j,i,4) = WTS3(i)*bPHI(i,j,4)   !....Right boundary       
    ENDDO
ENDDO

B(:,:,1) = MATMUL(M2,IB(:,:,1))           !....bottom boundary matrix
B(:,:,2) = MATMUL(M2,IB(:,:,2))           !....top boundary matrix
B(:,:,3) = MATMUL(M2,IB(:,:,3))           !....Left boundary matrix
B(:,:,4) = MATMUL(M2,IB(:,:,4))           !....Right boundary matrix

IF (KERNEL == 1) THEN !....Change for non-full matrix?

    M3   = 0.0d0
    IL2  = 0.0d0
    L2K  = 0.0d0

    !....Mass Matrix
    
    IF (FULL2D == 1) THEN
        k = 0
        DO i = 0, pz
            DO j = 0, pz
                k = k+1
                Mx = (2.0d0*j+1.0d0)/2.0d0       !....Inverse of Mass Matrix
                My = (2.0d0*i+1.0d0)/2.0d0
                M3(k,k) = Mx*My
            ENDDO
        ENDDO
    ELSEIF (FULL2D == 0) THEN
        DO i = 1, NDOF3
            DO j = 1, NDOF3
                DO k = 1, PTS4
                    M3(i,j) = M3(i,j) + WTS4(k)*PHIk(k,i)*PHIk(k,j)
                ENDDO
                IF (i == j) THEN
                    M3(i,j) = 1.0d0/M3(i,j)
                ENDIF    
            ENDDO
        ENDDO    
    ENDIF    

    !....L2k matrix

    DO j = 1, PTS4
        DO i = 1, NDOF3
            IL2(i,j) = WTS4(j)*PHIk(j,i)
        ENDDO
    ENDDO 

    L2K = MATMUL(M3,IL2)  
    
ENDIF    

RETURN
END SUBROUTINE DG_W_MATRICIES
!---------------------------------------------------------------------
!
!                DG Matricies Subroutine
!              Written by Colton J. Conroy
!                    @ the C.H.I.L
!                        1.29.13
!
!---------------------------------------------------------------------
SUBROUTINE DG_MATRICIES_3D(p,WTS1,WTS2,WTS3,PTS1,PTS2,PTS3,NDOF,A,L2,C,B,L2PHI,PHI,DPHI,bPHI)

USE precisions
USE global, ONLY: FULL3D
IMPLICIT NONE

INTEGER(int_p) :: i, j, k, ii, p, PTS1, PTS2, PTS3, NDOF
REAL(real_p) :: Mx, My, Mz
REAL(real_p), DIMENSION(PTS1) :: WTS1
REAL(real_p), DIMENSION(PTS2) :: WTS2
REAL(real_p), DIMENSION(PTS3) :: WTS3
REAL(real_p), DIMENSION(NDOF,NDOF) :: M
REAL(real_p), DIMENSION(NDOF,PTS1,3) :: A, IA
REAL(real_p), DIMENSION(NDOF,PTS2) :: L2, C, IL2, IC
REAL(real_p), DIMENSION(NDOF,PTS3,6) :: B, IB
REAL(real_p), DIMENSION(PTS1,NDOF) :: PHI
REAL(real_p), DIMENSION(PTS2,NDOF) :: L2PHI
REAL(real_p), DIMENSION(PTS1,NDOF,3) :: DPHI
REAL(real_p), DIMENSION(PTS3,NDOF,6) :: bPHI

!....Initialize matricies

M   = 0.0d0
IA  = 0.0d0
IB  = 0.0d0
IC  = 0.0d0
IL2 = 0.0d0
A   = 0.0d0
B   = 0.0d0
C   = 0.0d0
L2  = 0.0d0

!....Mass Matrix
IF (FULL3D == 1) THEN
    ii = 0
    DO i = 0, p
        DO j = 0, p
            DO k = 0, p
                ii = ii + 1
                Mx = (2.0d0*j+1.0d0)/2.0d0       !....Inverse of Mass Matrix
                My = (2.0d0*i+1.0d0)/2.0d0
                Mz = (2.0d0*k+1.0d0)/2.0d0
                M(ii,ii) = Mx*My*Mz
            ENDDO    
        ENDDO
    ENDDO
ELSEIF (FULL3D == 0) THEN
    DO i = 1, NDOF
        DO j = 1, NDOF
            DO k = 1, PTS2
                M(i,j) = M(i,j) + WTS2(k)*L2PHI(k,i)*L2PHI(k,j)
            ENDDO
            IF (i == j) THEN
                M(i,j) = 1.0d0/M(i,j)
            ENDIF    
        ENDDO
    ENDDO
ENDIF        

!....Advection Matrix

DO i = 1, PTS1
    DO j = 1, NDOF
       IA(j,i,1) = WTS1(i)*DPHI(i,j,1)   !....x derivative
       IA(j,i,2) = WTS1(i)*DPHI(i,j,2)   !....y derivative
       IA(j,i,3) = WTS1(i)*DPHI(i,j,3)   !....z derivative
    ENDDO
ENDDO

A(:,:,1) = MATMUL(M,IA(:,:,1))          !....x advection matrix
A(:,:,2) = MATMUL(M,IA(:,:,2))          !....y advection matrix
A(:,:,3) = MATMUL(M,IA(:,:,3))          !....z advection matrix

!....L2 matrix

DO j = 1, PTS2
    DO i = 1, NDOF
        IL2(i,j) = WTS2(j)*L2PHI(j,i)
    ENDDO
ENDDO 

L2 = MATMUL(M,IL2)                     
C  = L2                                 !....Source matrix

!....Boundary Matrix

DO i = 1, PTS3
    DO j = 1, NDOF
       IB(j,i,1)  = WTS3(i)*bPHI(i,j,1)   !....Left Face   (-x)
       IB(j,i,2)  = WTS3(i)*bPHI(i,j,2)   !....Right Face  (+x)
       IB(j,i,3)  = WTS3(i)*bPHI(i,j,3)   !....Bottom Face (-z)
       IB(j,i,4)  = WTS3(i)*bPHI(i,j,4)   !....Top Face    (+z)
       IB(j,i,5)  = WTS3(i)*bPHI(i,j,5)   !....Front Face  (-y)
       IB(j,i,6)  = WTS3(i)*bPHI(i,j,6)   !....Back Face   (+y)            
    ENDDO
ENDDO

B(:,:,1) = MATMUL(M,IB(:,:,1))           !....Left Face   (-x) matrix
B(:,:,2) = MATMUL(M,IB(:,:,2))           !....Right Face  (+x) matrix
B(:,:,3) = MATMUL(M,IB(:,:,3))           !....Bottom Face (-z) matrix
B(:,:,4) = MATMUL(M,IB(:,:,4))           !....Top Face    (+z) matrix
B(:,:,5) = MATMUL(M,IB(:,:,5))           !....Front Face  (-y) matrix
B(:,:,6) = MATMUL(M,IB(:,:,6))           !....Back Face   (+y) matrix


RETURN
END SUBROUTINE DG_MATRICIES_3D 
!---------------------------------------------------------------------
!
!                DG Matricies Subroutine
!              Written by Colton J. Conroy
!                    @ the C.H.I.L
!                        1.29.13
!
!---------------------------------------------------------------------
SUBROUTINE DG_MATRICIES_2D(p,WTS1,WTS2,WTS3,PTS1,PTS2,PTS3,NDOF,A,L2,C,B,B2,L2PHI,PHI,DPHI,bPHI,bPHI2,A2)

USE precisions
USE global, ONLY: FULL2D 
IMPLICIT NONE

INTEGER(int_p) :: i, j, k, p, PTS1, PTS2, PTS3, NDOF
REAL(real_p) :: Mx, My
REAL(real_p), DIMENSION(PTS1) :: WTS1
REAL(real_p), DIMENSION(PTS2) :: WTS2
REAL(real_p), DIMENSION(PTS3) :: WTS3
REAL(real_p), DIMENSION(NDOF,NDOF) :: M
REAL(real_p), DIMENSION(NDOF,PTS1) :: A2, IA2
REAL(real_p), DIMENSION(NDOF,PTS1,2) :: A, IA
REAL(real_p), DIMENSION(NDOF,PTS2) :: L2, C, IL2, IC
REAL(real_p), DIMENSION(NDOF,PTS3,2) :: B, B2, IB, IB2
REAL(real_p), DIMENSION(PTS1,NDOF) :: PHI
REAL(real_p), DIMENSION(PTS2,NDOF) :: L2PHI
REAL(real_p), DIMENSION(PTS1,NDOF,2) :: DPHI
REAL(real_p), DIMENSION(PTS3,NDOF,4) :: bPHI, bPHI2

!....Initialize matricies

M   = 0.0d0
IA  = 0.0d0
IA2 = 0.0d0
IB  = 0.0d0
IC  = 0.0d0
IL2 = 0.0d0
A   = 0.0d0
A2  = 0.0d0
B   = 0.0d0
B2  = 0.0d0
C   = 0.0d0
L2  = 0.0d0

!....Mass Matrix
IF (FULL2D == 1) THEN
    k = 0
    DO i = 0, p
        DO j = 0, p
            k = k+1
            Mx = (2.0d0*j+1.0d0)/2.0d0       !....Inverse of Mass Matrix
            My = (2.0d0*i+1.0d0)/2.0d0
            M(k,k) = Mx*My
        ENDDO
    ENDDO
ELSEIF (FULL2D == 0) THEN
    DO i = 1, NDOF
        DO j = 1, NDOF
            DO k = 1, PTS2
                M(i,j) = M(i,j) + WTS2(k)*L2PHI(k,i)*L2PHI(k,j)
            ENDDO
            IF (i == j) THEN
                M(i,j) = 1.0d0/M(i,j)
            ENDIF    
        ENDDO
    ENDDO
ENDIF        

!....Advection Matrix

DO i = 1, PTS1
    DO j = 1, NDOF     
       IA(j,i,1) = WTS1(i)*DPHI(i,j,1)   !....x derivative
       IA(j,i,2) = WTS1(i)*DPHI(i,j,2)   !....z derivative
       IA2(j,i)  = WTS1(i)*PHI(i,j)         
    ENDDO
ENDDO

A(:,:,1) = MATMUL(M,IA(:,:,1))          !....x advection matrix
A(:,:,2) = MATMUL(M,IA(:,:,2))          !....z advection matrix
A2       = MATMUL(M,IA2)

!....L2 matrix

DO j = 1, PTS2
    DO i = 1, NDOF
        IL2(i,j) = WTS2(j)*L2PHI(j,i)
    ENDDO
ENDDO 

L2 = MATMUL(M,IL2)                     
C  = L2                                 !....Source matrix

!....Boundary Matrix

DO i = 1, PTS3
    DO j = 1, NDOF
       IB(j,i,1)  = WTS3(i)*bPHI(i,j,1)   !....bottom boundary
       IB(j,i,2)  = WTS3(i)*bPHI(i,j,2)   !....top boundary
       IB2(j,i,1) = WTS3(i)*bPHI2(i,j,1)   !....Left boundary
       IB2(j,i,2) = WTS3(i)*bPHI2(i,j,2)   !....Right boundary       
    ENDDO
ENDDO

B(:,:,1) = MATMUL(M,IB(:,:,1))           !....bottom boundary matrix
B(:,:,2) = MATMUL(M,IB(:,:,2))           !....top boundary matrix
B2(:,:,1) = MATMUL(M,IB2(:,:,1))           !....Left boundary matrix
B2(:,:,2) = MATMUL(M,IB2(:,:,2))           !....Right boundary matrix


RETURN
END SUBROUTINE DG_MATRICIES_2D 
!----------------------------------------------------------------------
!
!                1D DG Matricies Subroutine
!                Written by Colton J. Conroy
!               @ the Isaac Newton Institute
!                       12.11.12
!
!----------------------------------------------------------------------
SUBROUTINE DG_MATRICIES(WTS1,WTS2,WTS3,D,A,L2,C,B1,B2,L2PHI,PHI,DPHI,bPHI,PTS1,PTS2,L3,PTS3,L3PHI)

USE precisions
IMPLICIT NONE

INTEGER(int_p) :: i, j, p, PTS1, PTS2, PTS3, D, debug, ELEM
REAL(real_p), DIMENSION(D) :: WTS1
REAL(real_p), DIMENSION(PTS1) :: WTS2
REAL(real_p), DIMENSION(PTS3) :: WTS3
REAL(real_p), DIMENSION(D,D) :: M, IA, A, DPHI, PHI
REAL(real_p), DIMENSION(D,PTS1) :: IL2, IC, C, L2
REAL(real_p), DIMENSION(D,PTS3) :: I3, L3
REAL(real_p), DIMENSION(D) :: BB1, BB2, B1, B2
REAL(real_p), DIMENSION(PTS1,D) :: L2PHI
REAL(real_p), DIMENSION(PTS2,D) :: bPHI
REAL(real_p), DIMENSION(PTS3,D) :: L3PHI

!....Initialize matricies

M   = 0.0d0
IA  = 0.0d0
A   = 0.0d0
IL2 = 0.0d0
IC  = 0.0d0
L2  = 0.0d0
L3  = 0.0d0
C   = 0.0d0
B1  = 0.0d0
B2  = 0.0d0
BB1 = 0.0d0
BB2 = 0.0d0
I3  = 0.0d0

!....Mass matrix and element integration

p = D-1
DO i = 0, p
       M(i+1,i+1) =  (2.0d0*i+1.0d0)/2.0d0 ! Invert mass matrix
    DO j = 1, p+1
       IA(i+1,j) = WTS1(j)*DPHI(j,i+1)
    ENDDO
ENDDO

!....L2 Numerical integration

DO i = 0, p
    DO j = 1, PTS1
        IL2(i+1,j) = WTS2(j)*L2PHI(j,i+1)
    ENDDO
ENDDO 

!....L3 (Map to vertical to horiztonal) integration

DO i = 0, p
    DO j = 1, PTS3
	    I3(i+1,j) = WTS3(j)*L3PHI(j,i+1)
    ENDDO
ENDDO

!....Advection matrix

A = MATMUL(M,IA)

!....L2 matrix and source matrix

L2 = MATMUL(M,IL2)
C = L2

!....L3 matrix

L3 = MATMUL(M,I3)

!....Boundary matrix

DO j = 1, p+1
    BB1(j) = bPHI(1,j)
    BB2(j) = bPHI(2,j)
ENDDO

!....Left boundary at x = -1

B1 = MATMUL(M,BB1)

!....Right boundary at x = 1

B2 = MATMUL(M,BB2)
B2 = -B2

RETURN
END SUBROUTINE DG_MATRICIES
!-----------------------------------------------------------------------
!
!           RK Mesh subroutine for Vertical Velocity
!                 Written by Colton J. Conroy
!                       @ the C.H.I.L
!                           1.16.14
!
!-----------------------------------------------------------------------
SUBROUTINE RK_MESH()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: Ldz, V_elem, n
REAL(real_p) :: INVJAC

LNrk = nelemsZ*LRK+2
ALLOCATE(nodesZRK(LNrk))
ALLOCATE(zrkCONN(nelemsZ,LRK))

Ldz = nelemsZ*LRK+1

ALLOCATE(dzRK(Ldz))
ALLOCATE(zrkzCONN(Ldz))
ALLOCATE(RK_HEXconn(Ldz,Qelems))

nodesZRK   = 0.0d0
zrkCONN    = 0.0d0
dzRK       = 0.0d0
zrkzCONN   = 0.0d0
RK_HEXconn = 0.0d0

!....1D mesh for RK method for vertical velocity

kk = 1
ii = 0
nodesZRK(1) = z0
DO i = 1, nelemsZ
    INVJAC = 1/zELEM_jac(i)
    DO j = 1, LRK
        kk = kk + 1
        ii = ii + 1
        nodesZRK(kk) = Z(i) + INVJAC*(1.0d0+TPTS_1D(j,1))
        zrkCONN(i,j) = kk
        zrkzCONN(ii) = i
    ENDDO
ENDDO
nodesZRK(Ldz+1) = zN     
zrkzCONN(Ldz)   = nelemsZ   

!....Mesh spacing of RK mesh

DO i = 1, Ldz
    dzRK(i) = nodesZRK(i+1)-nodesZRK(i)
ENDDO

!....Connectivity

DO i = 1, Qelems
    DO n = 0, Ldz-1    
        V_elem            = zrkzCONN(n+1)
        RK_HEXconn(n+1,i) = HEXcolumn(V_elem,i)
    ENDDO    
ENDDO        
   
RETURN
END SUBROUTINE RK_MESH
!-----------------------------------------------------------------------
!
!               1D Mesh Subroutine
!           Written by Colton J. Conroy      
!           @ the Isaac Newton Institute
!                    12.10.12
!
!-----------------------------------------------------------------------
SUBROUTINE LINE_MESH()

USE global
IMPLICIT NONE

REAL(real_p) :: xN, x0, yN, y0

!....Mesh data

xN = PROB_DATA%xN
x0 = PROB_DATA%x0
IF (hx == 0) THEN
    dx = (xN-x0)/2.0d0
ELSE
    dx = (xN-x0)/(4.0d0*2.0d0**(hx-1))    !....continuity spacing (even # of elements)
ENDIF    
dxf = dx/VpH                              !....momentum spacing in horizontal (for quads)
yN = PROB_DATA%yN          
y0 = PROB_DATA%y0
IF (hy == 0) THEN
    dy = (yN-y0)/2.0d0
ELSE
    dy = (yN-y0)/(4.0d0*2.0d0**(hy-1))
ENDIF
IF (hz == 0) THEN    
    dz = (zN-z0)/2.0d0
ELSE    
    dz = (zN-z0)/(4.0d0*2.0d0**(hz-1))    !....momentum spacing in vertical
ENDIF    
nelemsX = 2**(1+hx)           !....number of x elements (continuity)
nnodesX = nelemsX + 1         !....number of x nodes    ("        ")
nelemsZ = 2**(1+hz)           !....etc
nnodesZ = nelemsZ + 1
nelemsY = 2**(1+hy)
nnodesY = nelemsY + 1
nelemsXF = VpH*nelemsX
nnodesXF = nelemsXF + 1

!....Allocate matricies

ALLOCATE(X(nnodesX))
ALLOCATE(Y(nnodesY))
ALLOCATE(Z(nnodesZ))
ALLOCATE(xELEM_nodes(nelemsX,2))
ALLOCATE(yELEM_nodes(nelemsY,2))
ALLOCATE(zELEM_nodes(nelemsZ,2))
ALLOCATE(xNODE_elems(nnodesX,2))
ALLOCATE(yNODE_elems(nnodesY,2))
ALLOCATE(zNODE_elems(nnodesZ,2))
ALLOCATE(xELEM_jac(nelemsX))
ALLOCATE(yELEM_jac(nelemsY))
ALLOCATE(zELEM_jac(nelemsZ))
ALLOCATE(XF(nnodesXF))
ALLOCATE(xfELEM_nodes(nelemsXF,2))
ALLOCATE(xfNODE_elems(nnodesXF,2))
ALLOCATE(xfELEM_jac(nelemsXF))

!....Initialize Matricies

X  = 0.0d0
Y  = 0.0d0
Z  = 0.0d0
XF = 0.0d0
xELEM_nodes  = 0.0d0
yELEM_nodes  = 0.0d0
zELEM_nodes  = 0.0d0
xfELEM_nodes = 0.0d0
xNODE_elems  = 0.0d0
yNODE_elems  = 0.0d0
zNODE_elems  = 0.0d0
xfNODE_elems = 0.0d0
xELEM_jac    = 0.0d0
yELEM_jac    = 0.0d0
zELEM_jac    = 0.0d0
xfELEM_jac   = 0.0d0

!....Continuity mesh

IF (meshX == 1) THEN
    X(1) = x0 
    DO i = 2, nelemsX
        X(i) = X(i-1) + dx
    ENDDO
    X(nnodesX) = xN
ELSEIF (meshX == 2) THEN
    ALLOCATE(xk(nelemsX-1))
    xk = 0.0d0
    CALL CHEBY_ROOTS(nelemsX-2,x0,xN,xk)
    X(1) = x0
    DO i = 2, nelemsX
        X(i) = xk(i-1)
    ENDDO
    X(nnodesX) = xN
    DEALLOCATE(xk)
ENDIF    

IF (meshY == 1) THEN
    Y(1) = y0
    DO i = 2, nelemsY
        Y(i) = Y(i-1) + dy
    ENDDO
    Y(nnodesY) = yN
ELSEIF (meshY == 2) THEN
    ALLOCATE(xk(nelemsY-1))
    xk = 0.0d0
    CALL CHEBY_ROOTS(nelemsY-2,y0,yN,xk)
    Y(1) = y0
    DO i = 2, nelemsY
        Y(i) = xk(i-1)
    ENDDO
    Y(nnodesY) = yN
    DEALLOCATE(xk)
ENDIF    
    
!....Vertical mesh for momentum

IF (meshZ == 1) THEN
    Z(1) = z0 ! z-mesh
    DO i = 2, nelemsZ
        Z(i) = Z(i-1) + dz
    ENDDO
    Z(nnodesZ) = zN
ELSEIF (meshZ == 2) THEN
    ALLOCATE(xk(nelemsZ-1))
    xk = 0.0d0
    CALL CHEBY_ROOTS(nelemsZ-2,z0,zN,xk)
    Z(1) = z0
    DO i = 2, nelemsZ
        Z(i) = xk(i-1)
    ENDDO
    Z(nnodesZ) = zN
ENDIF            

!....Horizontal mesh for momentum

IF (meshX == 1) THEN

    XF(1) = x0 ! fine-mesh
    DO i = 2, nelemsXF
        XF(i) = XF(i-1) + dxf
    ENDDO
    XF(nnodesXF) = xN

ELSEIF (meshX == 2) THEN

    XF = X
    
ENDIF

!....Connectivity and Jacobians

DO i = 1, nelemsX 
    xELEM_nodes(i,1) = i
    xELEM_nodes(i,2) = i+1
    xELEM_jac(i) = 2.0d0/(X(i+1)-X(i)) 
ENDDO

DO i = 1, nelemsXF 
    xfELEM_nodes(i,1) = i
    xfELEM_nodes(i,2) = i+1
    xfELEM_jac(i) = 2.0d0/(XF(i+1)-XF(i))
ENDDO

DO i = 1, nelemsY
    yELEM_nodes(i,1) = i
    yELEM_nodes(i,2) = i+1
    yELEM_jac(i) = 2.0d0/(Y(i+1)-Y(i))
ENDDO    

DO i = 1, nelemsZ
    zELEM_nodes(i,1) = i
    zELEM_nodes(i,2) = i+1
    zELEM_jac(i) = 2.0d0/(Z(i+1)-Z(i))
ENDDO

!....Set up node element connectivity for periodic BCs

xNODE_elems(1,1) = nelemsX
xNODE_elems(1,2) = 1
DO i = 2, nnodesX-1
    xNODE_elems(i,1) = i-1
    xNODE_elems(i,2) = i
ENDDO
xNODE_elems(nnodesX,1) = nelemsX
xNODE_elems(nnodesX,2) = 1

xfNODE_elems(1,1) = nelemsXF 
xfNODE_elems(1,2) = 1
DO i = 2, nnodesXF-1
    xfNODE_elems(i,1) = i-1
    xfNODE_elems(i,2) = i
ENDDO
xfNODE_elems(nnodesXF,1) = nelemsXF
xfNODE_elems(nnodesXF,2) = 1

yNODE_elems(1,1) = nelemsY
yNODE_elems(1,2) = 1
DO i = 2, nnodesY-1
    yNODE_elems(i,1) = i-1
    yNODE_elems(i,2) = i
ENDDO
yNODE_elems(nnodesY,1) = nelemsY
yNODE_elems(nnodesY,2) = 1

zNODE_elems(1,1) = nelemsZ
zNODE_elems(1,2) = 1
DO i = 2, nnodesZ-1
    zNODE_elems(i,1) = i-1
    zNODE_elems(i,2) = i
ENDDO
zNODE_elems(nnodesZ,1) = nelemsZ
zNODE_elems(nnodesZ,2) = 1

IF (debug == 1) THEN
   write(*,*)'y=', yNODE_elems, 'yELEM=', yELEM_nodes, 'y_jac=', yELEM_jac
ENDIF

END SUBROUTINE LINE_MESH
!-----------------------------------------------------------------------
!   
!                         Subroutine Quad Mesh (XY)
!                        Written by Colton J Conroy
!                             @ the C.H.I.L
!                                 1.24.13
!
!-----------------------------------------------------------------------
SUBROUTINE QUAD_MESHy(REGION)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: check, nedgesY, REGION, k
REAL(real_p) :: b, h
REAL(real_p), DIMENSION(4) :: node

!....Constants

Qelems  = nelemsY*nelemsXF                       !....# of Quads
Qnnodes = nnodesXF*nnodesY                       !....# of Quad. nodes 

!....Allocate Matricies

ALLOCATE(Qnode(nnodesY*nnodesXF,2))
ALLOCATE(CONN(Qelems,4))
ALLOCATE(CONN_jac(Qelems,1))
ALLOCATE(CONNzx(Qelems,1))
ALLOCATE(CONNsigma(Qelems,1))
ALLOCATE(CONNgrid(nelemsY,nelemsX))
ALLOCATE(CONNielem(Qelems,2))
ALLOCATE(CONN_jacY(Qelems,1))
ALLOCATE(CONN_jacX(Qelems,1))
ALLOCATE(CONN_jacF(Qelems,1))
ALLOCATE(CONN_jacS(Qelems,1))
ALLOCATE(xEDGE(nnodesXF*nelemsY,3))
ALLOCATE(yEDGE(nelemsXF*nnodesY,3))
     
!....Initialize Matricies

Qnode     = 0.0d0
CONN      = 0.0d0
CONN_jac  = 0.0d0
CONN_jacY = 0.0d0
CONN_jacX = 0.0d0
CONN_jacS = 0.0d0
node      = 0.0d0
xEDGE     = 0.0d0
yEDGE     = 0.0d0

!....Nodes of Quads

kk = 0

DO i = 1, nnodesXF

    DO j = 1, nnodesY
    
        kk = kk + 1
        Qnode(kk,1) = XF(i)
        Qnode(kk,2) = Y(j)
        
    ENDDO
    
ENDDO

!....Node-Element Connectivity

ii = 0
kk = 0

DO i = 1, nelemsXF
    
    DO j = 1, nelemsY
    
        kk    = kk + 1
        ii    = ii + 1
        check = (i-1)*nnodesY
         
        IF (kk == check) THEN
        
            CONN(ii,1) = kk + 1
            CONN(ii,2) = kk + 2
            CONN(ii,3) = kk + nnodesY + 2
            CONN(ii,4) = kk + nnodesY + 1
            kk = kk + 1
            
        ELSE

            CONN(ii,1) = kk 
            CONN(ii,2) = kk + 1
            CONN(ii,3) = kk + nnodesY + 1
            CONN(ii,4) = kk + nnodesY
            
        ENDIF
        
   ENDDO
   
ENDDO

!....xz connectivity

ii = 0

DO i = 1, nelemsX
    
    DO j = 1, VpH
    
        DO kk = 1, nelemsY
    
            ii         = ii + 1
            CONNzx(ii,1) = i
            
        ENDDO
        
    ENDDO
    
ENDDO    

!....Sigma Layer info

ii = 0

DO i = 1, nelemsX

    DO j = 1, nelemsY
        
        ii = ii + 1
        CONNsigma(ii,1) = j
        CONNgrid(j,i)   = ii
        CONNielem(ii,1) = i
        CONNielem(ii,2) = j
        
    ENDDO
    
ENDDO    

!....Element Jacobians

kk = 0

DO i = 1, nelemsXF

    DO j = 1, nelemsY
    
        kk        = kk + 1
        node      = CONN(kk,:)
        b         = Qnode(node(4),1) - Qnode(node(1),1)
        h         = Qnode(node(2),2) - Qnode(node(1),2)

        CONN_jac(kk,1)  = 4.0d0/(b*h)
        CONN_jacY(kk,1) = 2.0d0/h
        CONN_jacX(kk,1) = 2.0d0/b
        CONN_jacF(kk,1) = CONN_jacY(kk,1)
        CONN_jacS(kk,1) = h/b

    ENDDO
    
ENDDO

!....Edge Element connectivity

kk = 0

DO j = 1, nelemsXF

    DO i = 1, nnodesY
    
        kk = kk + 1
        
        IF (i == 1) THEN
        
            yEDGE(kk,1) = -1
            yEDGE(kk,2) = kk - j + 1
            yEDGE(kk,3) = yEDGE(kk,2) + (nelemsY - 1)
            
        ELSEIF (i == nnodesY) THEN
        
            yEDGE(kk,1) = kk - j
            yEDGE(kk,2) = 0
            yEDGE(kk,3) = yEDGE(kk,1) - (nelemsY - 1)
            
        ELSE
        
            yEDGE(kk,1) = kk - j
            yEDGE(kk,2) = kk - j + 1
            
        ENDIF

    ENDDO
    
ENDDO

kk = 0

DO j = 1, nnodesXF

    DO i = 1, nelemsY
    
        kk = kk + 1
        
        IF (j == 1) THEN
        
            xEDGE(kk,1) = -1
            xEDGE(kk,2) = kk
            xEDGE(kk,3) = xEDGE(kk,2) + (nelemsX - 1)*nelemsY
            
        ELSEIF (j == nnodesX) THEN
        
            xEDGE(kk,1) = kk - nelemsY
            xEDGE(kk,2) = 0
            xEDGE(kk,3) = xEDGE(kk,1) - (nelemsX - 1)*nelemsY
            
        ELSE
        
            xEDGE(kk,1) = kk - nelemsY
            xEDGE(kk,2) = kk
            
        ENDIF

    ENDDO
    
ENDDO

ALLOCATE(QedgeX(Qelems,2))
ALLOCATE(QedgeY(Qelems,2))

QedgeX = 0.0d0
QedgeY = 0.0d0

k = 1
kk = 0
check = nelemsY

DO i = 1, Qelems

    IF (i > check) THEN
        k  = k + 1
        kk = kk + 1
        check = k*nelemsY
    ENDIF
    
    QedgeX(i,1) = i 
    QedgeX(i,2) = i + nelemsY
    QedgeY(i,1) = i + 1*kk
    QedgeY(i,2) = i + 1*k
    
ENDDO    
 

RETURN
END SUBROUTINE QUAD_MESHy
!-----------------------------------------------------------------------
!   
!                           Subroutine Quad Mesh
!                        Written by Colton J Conroy
!                             @ the C.H.I.L
!                                 1.24.13
!
!-----------------------------------------------------------------------
SUBROUTINE QUAD_MESH(REGION)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: check, nedgesZ, REGION
REAL(real_p) :: b, h
REAL(real_p), DIMENSION(4) :: node

!....Constants

Qelems  = nelemsZ*nelemsXF                       !....# of Quads
Qnnodes = nnodesXF*nnodesZ                       !....# of Quad. nodes 

!....Allocate Matricies

ALLOCATE(Qnode(nnodesZ*nnodesXF,2))
ALLOCATE(CONN(Qelems,4))
ALLOCATE(CONN_jac(Qelems,1))
ALLOCATE(CONNzx(Qelems,1))
ALLOCATE(CONNsigma(Qelems,1))
ALLOCATE(CONNgrid(nelemsZ,nelemsX))
ALLOCATE(CONN_jacZ(Qelems,1))
ALLOCATE(CONN_jacX(Qelems,1))
ALLOCATE(CONN_jacF(Qelems,1))
ALLOCATE(CONN_jacS(Qelems,1))
ALLOCATE(xEDGE(nnodesXF*nelemsZ,3))
ALLOCATE(zEDGE(nelemsXF*nnodesZ,3))
     
!....Initialize Matricies

Qnode     = 0.0d0
CONN      = 0.0d0
CONN_jac  = 0.0d0
CONN_jacZ = 0.0d0
CONN_jacX = 0.0d0
CONN_jacS = 0.0d0
node      = 0.0d0
xEDGE     = 0.0d0
zEDGE     = 0.0d0

!....Nodes of Quads

kk = 0

DO i = 1, nnodesXF

    DO j = 1, nnodesZ
    
        kk = kk + 1
        Qnode(kk,1) = XF(i)
        Qnode(kk,2) = Z(j)
        
    ENDDO
    
ENDDO

!....Node-Element Connectivity

ii = 0
kk = 0

DO i = 1, nelemsXF
    
    DO j = 1, nelemsZ
    
        kk    = kk + 1
        ii    = ii + 1
        check = (i-1)*nnodesZ
         
        IF (kk == check) THEN
        
            CONN(ii,1) = kk + 1
            CONN(ii,2) = kk + 2
            CONN(ii,3) = kk + nnodesZ + 2
            CONN(ii,4) = kk + nnodesZ + 1
            kk = kk + 1
            
        ELSE

            CONN(ii,1) = kk 
            CONN(ii,2) = kk + 1
            CONN(ii,3) = kk + nnodesZ + 1
            CONN(ii,4) = kk + nnodesZ
            
        ENDIF
        
   ENDDO
   
ENDDO

!....xz connectivity

ii = 0

DO i = 1, nelemsX
    
    DO j = 1, VpH
    
        DO kk = 1, nelemsZ
    
            ii         = ii + 1
            CONNzx(ii,1) = i
            
        ENDDO
        
    ENDDO
    
ENDDO    

!....Sigma Layer info

ii = 0

DO i = 1, nelemsX

    DO j = 1, nelemsZ
        
        ii = ii + 1
        CONNsigma(ii,1) = j
        CONNgrid(j,i)   = ii
        
    ENDDO
    
ENDDO    

!....Element Jacobians

kk = 0

DO i = 1, nelemsXF

    DO j = 1, nelemsZ
    
        kk        = kk + 1
        node      = CONN(kk,:)
        b         = Qnode(node(4),1) - Qnode(node(1),1)
        h         = Qnode(node(2),2) - Qnode(node(1),2)

        CONN_jac(kk,1)  = 4.0d0/(b*h)
        CONN_jacZ(kk,1) = 2.0d0/h
        CONN_jacX(kk,1) = 2.0d0/b
        CONN_jacF(kk,1) = CONN_jacZ(kk,1)
        CONN_jacS(kk,1) = h/b

    ENDDO
    
ENDDO

!....Edge Element connectivity

kk = 0

DO j = 1, nelemsXF

    DO i = 1, nnodesZ
    
        kk = kk + 1
        
        IF (i == 1) THEN
        
            zEDGE(kk,1) = -1
            zEDGE(kk,2) = kk - j + 1
            zEDGE(kk,3) = zEDGE(kk,2) + (nelemsZ - 1)
            
        ELSEIF (i == nnodesZ) THEN
        
            zEDGE(kk,1) = kk - j
            zEDGE(kk,2) = 0
            zEDGE(kk,3) = zEDGE(kk,1) - (nelemsZ - 1)
            
        ELSE
        
            zEDGE(kk,1) = kk - j
            zEDGE(kk,2) = kk - j + 1
            
        ENDIF

    ENDDO
    
ENDDO

kk = 0

DO j = 1, nnodesXF

    DO i = 1, nelemsZ
    
        kk = kk + 1
        
        IF (j == 1) THEN
        
            xEDGE(kk,1) = -1
            xEDGE(kk,2) = kk
            xEDGE(kk,3) = xEDGE(kk,2) + (nelemsX - 1)*nelemsZ
            
        ELSEIF (j == nnodesX) THEN
        
            xEDGE(kk,1) = kk - nelemsZ
            xEDGE(kk,2) = 0
            xEDGE(kk,3) = xEDGE(kk,1) - (nelemsX - 1)*nelemsZ
            
        ELSE
        
            xEDGE(kk,1) = kk - nelemsZ
            xEDGE(kk,2) = kk
            
        ENDIF

    ENDDO
    
ENDDO


RETURN
END SUBROUTINE QUAD_MESH
!---------------------------------------------------------------------
!
!                   Subroutine Hexahedral Mesh
!                    Written by Colton J. Conroy
!                         @ the C.H.I.L
!                             7.31.13
!
!----------------------------------------------------------------------
SUBROUTINE HEX_MESH()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: nedgesZ, jj, k, iii, check
INTEGER(int_p), DIMENSION(nelemsXF*nelemsZ,4) :: CONNz
REAL(real_p) :: b, h, w
REAL(real_p), DIMENSION(8) :: node


!....Constants

HEXelems   = nelemsZ*nelemsXF*nelemsY                       !....# of Hex
HEXnnodes  = nnodesXF*nnodesZ*nnodesY                       !....# of Hex nodes 

!....Allocate Matricies

ALLOCATE(HEXnode(nnodesZ*nnodesXF*nnodesY,3))
ALLOCATE(HEXconn(HEXelems,8))
ALLOCATE(HEX_jac(HEXelems,1))
ALLOCATE(HEXzxy(HEXelems,1))
ALLOCATE(HEXsigma(HEXelems,1))
ALLOCATE(HEXgrid(nelemsZ,nelemsX,nelemsY))
ALLOCATE(HEXcolumn(nelemsZ,nelemsX*nelemsY))
ALLOCATE(HEX_jacZ(HEXelems,1))
ALLOCATE(HEX_jacX(HEXelems,1))
ALLOCATE(HEX_jacY(HEXelems,1))
ALLOCATE(HEX_jacZw(HEXelems,1))
ALLOCATE(HEX_jacXw(HEXelems,1))
ALLOCATE(HEX_jacYw(HEXelems,1))
ALLOCATE(HEX_jacF(HEXelems,1))
ALLOCATE(HEX_jacS(HEXelems,1))
ALLOCATE(topLAYER(Qelems))
ALLOCATE(xFACE(nnodesXF*nelemsZ*nelemsY,3))
ALLOCATE(yFACE(nelemsXF*nelemsZ*nnodesY,3))
ALLOCATE(zFACE(nelemsXF*nnodesZ*nelemsY,3))

!....Initialize Matricies

HEXnode   = 0.0d0
HEXconn   = 0.0d0
HEXcolumn = 0.0d0
HEX_jac   = 0.0d0
HEX_jacZ  = 0.0d0
HEX_jacY  = 0.0d0
HEX_jacX  = 0.0d0
HEX_jacZw = 0.0d0
HEX_jacYw = 0.0d0
HEX_jacXw = 0.0d0
topLAYER  = 0.0d0
HEX_jacS  = 0.0d0
xFACE     = 0.0d0
yFACE     = 0.0d0
zFACE     = 0.0d0
node      = 0.0d0
CONNz     = 0.0d0

kk = 0

DO i = 1, nnodesY

    DO j = 1, nnodesXF

        DO ii = 1, nnodesZ
    
            kk = kk + 1
            HEXnode(kk,1) = XF(j)
            HEXnode(kk,2) = Y(i)
            HEXnode(kk,3) = Z(ii)
        
        ENDDO
    
    ENDDO    
    
ENDDO

kk = 0
ii = 0

DO i = 1, nelemsXF

     DO j = 1, nelemsZ
     
        kk = kk + 1
        ii = ii + 1
        check = (i-1)*nnodesZ
        
        IF (kk == check) THEN
        
            CONNz(ii,1) = kk + 1
            CONNz(ii,2) = kk + 2
            CONNz(ii,3) = kk + nnodesZ + 2
            CONNz(ii,4) = kk + nnodesZ + 1
            kk = kk + 1
            
        ELSE
        
            CONNz(ii,1) = kk
            CONNz(ii,2) = kk + 1
            CONNz(ii,3) = kk + nnodesZ + 1
            CONNz(ii,4) = kk + nnodesZ
            
        ENDIF
        
     ENDDO
     
ENDDO        


!....Element-to-node connectivity

kk = 0

DO ii = 1, nelemsY

    DO i = 1, nelemsXF
        
        DO j = 1, nelemsZ
    
            kk    = kk + 1
            
            IF (ii == 1) THEN   !....(0,y,z) Boundary face 

                HEXconn(kk,1) = CONNz(kk,1)
                HEXconn(kk,2) = CONNz(kk,2)
                HEXconn(kk,3) = CONNz(kk,3)
                HEXconn(kk,4) = CONNz(kk,4)
                HEXconn(kk,5) = HEXconn(kk,1) + nnodesXF*nnodesZ
                HEXconn(kk,6) = HEXconn(kk,2) + nnodesXF*nnodesZ
                HEXconn(kk,7) = HEXconn(kk,3) + nnodesXF*nnodesZ
                HEXconn(kk,8) = HEXconn(kk,4) + nnodesXF*nnodesZ
                
           ELSE !....Remaining faces 
           
                jj = nelemsXF*nelemsZ
                HEXconn(kk,1) = HEXconn(kk-jj,5) 
                HEXconn(kk,2) = HEXconn(kk-jj,6) 
                HEXconn(kk,3) = HEXconn(kk-jj,7) 
                HEXconn(kk,4) = HEXconn(kk-jj,8) 
                HEXconn(kk,5) = HEXconn(kk-jj,5) + nnodesXF*nnodesZ
                HEXconn(kk,6) = HEXconn(kk-jj,6) + nnodesXF*nnodesZ
                HEXconn(kk,7) = HEXconn(kk-jj,7) + nnodesXF*nnodesZ
                HEXconn(kk,8) = HEXconn(kk-jj,8) + nnodesXF*nnodesZ
                
           ENDIF     
           
    
        ENDDO
                
   ENDDO
   
ENDDO

!....Element Jacobians

kk = 0

DO ii = 1, nelemsY

    DO i = 1, nelemsXF

        DO j = 1, nelemsZ
    
            kk        = kk + 1
            node      = HEXconn(kk,:)
            b         = HEXnode(node(4),1) - HEXnode(node(1),1)
            h         = HEXnode(node(2),3) - HEXnode(node(1),3)
            w         = HEXnode(node(5),2) - HEXnode(node(1),2)
            
            HEX_jac(kk,1)  = 8.0d0/(b*h*w)
            HEX_jacZ(kk,1) = 2.0d0/h
            HEX_jacX(kk,1) = 2.0d0/b
            HEX_jacY(kk,1) = 2.0d0/w
            
            HEX_jacZw(kk,1) = 1.0d0
            HEX_jacXw(kk,1) = h/b
            HEX_jacYw(kk,1) = h/w
            
        ENDDO    

    ENDDO
    
ENDDO

!....xyz connectivity

ii = 0


DO i = 1, nelemsY

    jj = i

    DO j = 1, nelemsXF
        
 
        DO kk = 1, nelemsZ
    
            ii         = ii + 1
            HEXzxy(ii,1) = jj
            
        ENDDO
        
        jj = jj + nelemsY
        
    ENDDO
    
ENDDO 

!....Column connectivity

ii = 0

DO i = 1, nelemsXF

    DO j = 1, nelemsY
    
        ii = ii + 1
    
        DO k = 1, nelemsZ
        
            HEXcolumn(k,ii) = (i-1)*nelemsZ + (j-1)*nelemsX*nelemsZ + k
            
        ENDDO
        
    ENDDO
        
ENDDO            
        
            


!....Face Element connectivity

!....Z faces

kk = 0
iii = 0

DO ii = 1, nelemsY

    DO j = 1, nelemsXF

        DO i = 1, nnodesZ
    
            kk = kk + 1
            
            IF (ii == 1) THEN
        
                IF (i == 1) THEN
        
                    zFACE(kk,1) = -1
                    zFACE(kk,2) = kk - j + 1
                    zFACE(kk,3) = zFACE(kk,2) + (nelemsZ - 1)
            
                ELSEIF (i == nnodesZ) THEN
        
                    zFACE(kk,1) = kk - j
                    zFACE(kk,2) = 0
                    zFACE(kk,3) = zFACE(kk,1) - (nelemsZ - 1)
                    
                    iii = iii + 1
                    topLAYER(iii) = kk - j
            
                ELSE
        
                    zFACE(kk,1) = kk - j
                    zFACE(kk,2) = kk - j + 1
            
                ENDIF
                
            ELSE
            
                IF (i == 1) THEN
                
                    zFACE(kk,1) = -1
                    zFACE(kk,2) = kk - (ii-1)*nelemsX - j + 1
                    zFACE(kk,3) = zFACE(kk,2) + (nelemsZ - 1)
                    
                ELSEIF (i == nnodesZ) THEN
                
                    zFACE(kk,1) = kk - (ii-1)*nelemsX - j 
                    zFACE(kk,2) = 0
                    zFACE(kk,3) = zFACE(kk,1) - (nelemsZ - 1)
                    
                ELSE
                
                    zFACE(kk,1) = kk - (ii-1)*nelemsX - j 
                    zFACE(kk,2) = kk - (ii-1)*nelemsX - j + 1 
                    
                ENDIF  
                
            ENDIF      
            
        ENDDO    

    ENDDO
    
ENDDO

!....X faces

kk = 0

DO ii = 1, nelemsY

    DO j = 1, nnodesXF

        DO i = 1, nelemsZ
    
            kk = kk + 1
            
            IF (ii == 1) THEN
        
                IF (j == 1) THEN
        
                    xFACE(kk,1) = -1
                    xFACE(kk,2) = kk
                    xFACE(kk,3) = xFACE(kk,2) + (nelemsX - 1)*nelemsZ
            
                ELSEIF (j == nnodesX) THEN
        
                    xFACE(kk,1) = kk - nelemsZ
                    xFACE(kk,2) = 0
                    xFACE(kk,3) = xFACE(kk,1) - (nelemsX - 1)*nelemsZ
            
                ELSE
        
                    xFACE(kk,1) = kk - nelemsZ
                    xFACE(kk,2) = kk
                    
                ENDIF
                
            ELSE
            
                IF (j == 1) THEN
                
                    xFACE(kk,1) = -1
                    xFACE(kk,2) = kk - (ii-1)*nelemsZ
                    xFACE(kk,3) = xFACE(kk,2) + (nelemsX - 1)*nelemsZ
                    
                ELSEIF (j == nnodesX) THEN
                
                    xFACE(kk,1) = kk - nelemsZ  - (ii-1)*nelemsZ
                    xFACE(kk,2) = 0
                    xFACE(kk,3) = xFACE(kk,1) - (nelemsX - 1)*nelemsZ
                    
                ELSE
                
                    xFACE(kk,1) = kk - nelemsZ - (ii-1)*nelemsZ
                    xFACE(kk,2) = kk - (ii-1)*nelemsZ
                    
                ENDIF    
            
            ENDIF
            
        ENDDO

    ENDDO
    
ENDDO

!....Y faces

kk = 0

DO ii = 1, nnodesY

    DO j = 1, nelemsXF

        DO i = 1, nelemsZ
    
            kk = kk + 1
            
            IF (ii == 1) THEN
        
                yFACE(kk,1) = -1
                yFACE(kk,2) = kk
                yFACE(kk,3) = yFACE(kk,2) + (nelemsY - 1)*nelemsX*nelemsY
                
            ELSEIF (ii == nnodesY) THEN
            
                yFACE(kk,1) = kk - nelemsXF*nelemsZ
                yFACE(kk,2) = 0
                yFACE(kk,3) = yFACE(kk,1) - (nelemsY - 1)*nelemsX*nelemsY
            
            ELSE
            
                yFACE(kk,1) = kk - nelemsXF*nelemsZ
                yFACE(kk,2) = kk
            
            ENDIF

        ENDDO

    ENDDO
    
ENDDO

ALLOCATE(HEXfaceX(HEXelems,2))
ALLOCATE(HEXfaceY(HEXelems,2))

HEXfaceX = 0.0d0
HEXfaceY = 0.0d0

k  = 1
kk = 0
check = nelemsZ*nelemsX

DO i = 1, HEXelems

    IF (i > check) THEN
        k  = k + 1
        kk = kk + 1
        check = k*nelemsX*nelemsZ
    ENDIF

    HEXfaceX(i,1) = i + nelemsZ*kk
    HEXfaceX(i,2) = i + nelemsZ*k
    HEXfaceY(i,1) = i 
    HEXfaceY(i,2) = i + nelemsZ*nelemsX
    
ENDDO


IF (debug == 30) THEN
     write(*,*)'XYZconn',HEXzxy(44,1)
ENDIF

END SUBROUTINE HEX_MESH
!----------------------------------------------------------------------
!
!                    Subroutine Get Number of Int. Points
!                          Written by Colton J Conroy
!                                @ the C.H.I.L
!                                   1.28.12
!
!-----------------------------------------------------------------------
SUBROUTINE NUM_INT_PTS(D,N,REGION)

USE precisions
IMPLICIT NONE

INTEGER(int_p) :: I, L, N, D, REGION
REAL(real_p), ALLOCATABLE :: M(:)

!-----------------------------------------------------------------------
!
!    2D Quadrilaterals
!
!-----------------------------------------------------------------------

IF (REGION == 2) THEN
!-----------------------------------------------------------------------
        IF (D.LE.1) THEN                                      ! Ref [1]*
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE( M(L) )
          M(1:L) = (/ 1 /)                                
!-----------------------------------------------------------------------
        ELSEIF (D.LE.2) THEN                                  ! Ref [1]
!-----------------------------------------------------------------------
          L = 3
          ALLOCATE( M(L) )
          M(1:L) = (/ 1,1,1 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.3) THEN                                  ! Ref [1]*
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE( M(L) )
          M(1:L) = (/ 4 /)                                
!-----------------------------------------------------------------------
        ELSEIF (D.LE.4) THEN                                  ! Ref [2]*
!-----------------------------------------------------------------------
          L = 6
          ALLOCATE( M(L) )
          M(1:L) = (/ 1,1,1,1,1,1 /)                          
!-----------------------------------------------------------------------
        ELSEIF (D.LE.5) THEN                                  ! Ref [1]*
!-----------------------------------------------------------------------
          L = 7
          ALLOCATE( M(L) )
          M(1:L) = (/ 1, 1, 1, 1, 1, 1, 1 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.6) THEN                                  ! Ref [2]*
!-----------------------------------------------------------------------
          L = 10
          ALLOCATE( M(L) )
          M(1:L) = (/ 1,1,1,1,1,1,1,1,1,1 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.7) THEN                                  ! Ref [3]*
!-----------------------------------------------------------------------
          L = 4
          ALLOCATE( M(L) )
          M(1:L) = (/ 2,2,2,2 /)                          
!-----------------------------------------------------------------------
        ELSEIF (D.LE.8) THEN                                  ! Ref [2]
!-----------------------------------------------------------------------
          L = 16
          ALLOCATE( M(L) )
          M(1:L) = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.9) THEN                                  ! Ref [5]*
!-----------------------------------------------------------------------
          L = 5
          ALLOCATE( M(L) )
          M(1:L) = (/ 1,4,4,4,4 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.10) THEN                                 ! Ref []
!-----------------------------------------------------------------------
          L = 22
          ALLOCATE( M(L) )
          M(1:L) = 1
!-----------------------------------------------------------------------
        ELSEIF (D.LE.11) THEN                                 ! Ref [5]*
!-----------------------------------------------------------------------
          L = 6
          ALLOCATE( M(L) )
          M(1:L) = (/ 4,4,4,4,4,4 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.12) THEN                                 ! Ref [ ]
!-----------------------------------------------------------------------
          L = 31
          ALLOCATE( M(L) )
          M(1:L) = 1
!-----------------------------------------------------------------------
        ELSEIF (D.LE.13) THEN                                 ! Ref [5]
!-----------------------------------------------------------------------
          L = 9
          ALLOCATE( M(L) )
          M(1:L) = (/ 1,4,4,4,4,4,4,4,4 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.15) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 11
          ALLOCATE( M(L) )
          M(1:L) = (/ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 /)     
!-----------------------------------------------------------------------
        ELSEIF (D.LE.17) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 14
          ALLOCATE( M(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4 /)               
!-----------------------------------------------------------------------
        ELSEIF (D.LE.19) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 17
          ALLOCATE( M(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 /)          
!-----------------------------------------------------------------------
        ELSEIF (D.LE.21) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 21
          ALLOCATE( M(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,1 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.23) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 25
          ALLOCATE( M(L) )
          M(1:L) =(/4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/)
!-----------------------------------------------------------------------
        ELSE
!-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 23 ********  '
          PRINT*,'  ********       for REGION = SQUARE       ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

!.......Count up number of quadrature points

        N = 0
        DO I = 1,L
          N = N + M(I)
        ENDDO
     
ELSEIF (REGION == 3) THEN

!----------------------------------------------------------------------
!
!   3D Hexahedral
!
!----------------------------------------------------------------------
        IF (D.LE.1) THEN                                     ! Ref [1]***
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE( M(L))
          M(1:L) = (/ 1 /)                                ! Multiplicity
!-----------------------------------------------------------------------
        ELSEIF (D.LE.3) THEN                                 ! Ref [1]***
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE(  M(L) )
          M(1:L) = (/ 6 /)                                ! Multiplicity
!-----------------------------------------------------------------------
        ELSEIF (D.LE.5) THEN                                 ! Ref [1]***
!-----------------------------------------------------------------------
          L = 2
          ALLOCATE( M(L) )
          M(1:L) = (/ 6, 8 /)                             ! Multiplicity
!-----------------------------------------------------------------------
        ELSEIF (D.LE.7) THEN                                  ! Ref [2]**
!-----------------------------------------------------------------------
          L = 3
          ALLOCATE( M(L))
          M(1:L) = (/ 6, 8, 24 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.9) THEN                                  ! Ref [3]^
!-----------------------------------------------------------------------
          L = 5
          ALLOCATE(  M(L) )
          M(1:L) = (/ 6, 12, 8, 8, 24 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.11) THEN                                 ! Ref [3]^
!-----------------------------------------------------------------------
          L = 7
          ALLOCATE( M(L) )
          M(1:L) = (/ 6, 12, 8, 8, 8, 24, 24 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.13) THEN                                 ! Ref [1]*
!-----------------------------------------------------------------------
          L = 11
          ALLOCATE( M(L) )
          M(1:L) = (/ 1, 6, 6, 6, 8, 8, 8, 12, 24, 24, 48 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.15) THEN                                  ! Ref [1]
!-----------------------------------------------------------------------
          L = 14
          ALLOCATE( M(L))
          M(1:L) = (/ 1,6,6,6,8,8,8,24,24,24,24,24,24,48 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.17) THEN                                 ! Ref [1]*
!-----------------------------------------------------------------------
          L = 17
          ALLOCATE( M(L))
          M(1:L) = (/ 1,6,6,6,8,8,8,12,12,24,24,24,24,24,24,48,48 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.19) THEN                                 ! Ref [1]*
!-----------------------------------------------------------------------
          L = 24
          ALLOCATE( M(L))
          M(1:L) = (/ 1,6,6,6,6,6,8,8,8,8,12,12,12,24,24,24,24,24,24,24,  &
                                                          24,48,48,48 /)
!-----------------------------------------------------------------------
        ELSE
!-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 19 ********  '
          PRINT*,'  ********       for REGION = CUBE      ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

!.......Count up number of quadrature points

        N = 0
        DO I = 1,L
          N = N + M(I)
        ENDDO        
ENDIF

RETURN
END SUBROUTINE NUM_INT_PTS
!-----------------------------------------------------------------------
! 
!                  3D Tensor Product Quadrature Subroutine
!                        Written by Colton J Conroy
!                              @ the C.H.I.L
!                                  1.3.14
!
!-----------------------------------------------------------------------
SUBROUTINE TENSOR_QUAD_3D(PTS_3D,WTS_3D,PTS_1D,WTS_1D,NPTS_1D)

USE precisions
IMPLICIT NONE

INTEGER :: j, R, Q, P, NPTS_1D
REAL(real_p), DIMENSION(NPTS_1D,1) :: PTS_1D
REAL(real_p), DIMENSION(NPTS_1D) :: WTS_1D
REAL(real_p), DIMENSION(NPTS_1D*NPTS_1D*NPTS_1D,3) :: PTS_3D
REAL(real_p), DIMENSION(NPTS_1D*NPTS_1D*NPTS_1D) :: WTS_3D

PTS_3D = 0.0d0
WTS_3D = 0.0d0

j = 1
DO R = 1, NPTS_1D
    DO Q = 1, NPTS_1D
        DO P = 1, NPTS_1D
            PTS_3D(j,1) = PTS_1D(P,1)
            PTS_3D(j,2) = PTS_1D(Q,1)
            PTS_3D(j,3) = PTS_1D(R,1)
            WTS_3D(j)   = WTS_1D(P)*WTS_1D(Q)*WTS_1D(R)
            j = j+1
        ENDDO
    ENDDO
ENDDO    

RETURN
END SUBROUTINE TENSOR_QUAD_3D
!-----------------------------------------------------------------------
!
!              2D Tensor Product Quadrature Subroutine
!                   Written by Colton J. Conroy
!                          @ the C.H.I.L
!                              1.16.14
!
!-----------------------------------------------------------------------
SUBROUTINE TENSOR_QUAD_2D(PTS_2D,WTS_2D,PTS_1D,WTS_1D,NPTS_1D)

USE precisions
IMPLICIT NONE

INTEGER :: j, Q, P, NPTS_1D
REAL(real_p), DIMENSION(NPTS_1D,1) :: PTS_1D
REAL(real_p), DIMENSION(NPTS_1D) :: WTS_1D
REAL(real_p), DIMENSION(NPTS_1D*NPTS_1D,2) :: PTS_2D
REAL(real_p), DIMENSION(NPTS_1D*NPTS_1D) :: WTS_2D

PTS_2D = 0.0d0
WTS_2D = 0.0d0

j = 1
DO Q = 1, NPTS_1D
    DO P = 1, NPTS_1D
        PTS_2D(j,1) = PTS_1D(P,1)
        PTS_2D(j,2) = PTS_1D(Q,1)
        WTS_2D(j)   = WTS_1D(P)*WTS_1D(Q)
        j = j + 1
    ENDDO
ENDDO    

RETURN
END SUBROUTINE TENSOR_QUAD_2D
!-----------------------------------------------------------------------
!
!                     Subroutine Gauss-Radau Quadrature
!                        Written by Colton J. Conroy
!                               @ the C.H.I.L
!                                   2.5.14
!
!-----------------------------------------------------------------------
SUBROUTINE GR_QUADRATURE(D,PTS,WTS,NPTS)

USE precisions
IMPLICIT NONE

INTEGER(int_p) :: D, NPTS
REAL(real_p), DIMENSION(NPTS,1) :: PTS
REAL(real_p), DIMENSION(NPTS) :: WTS

PTS = 0.0d0
WTS = 0.0d0

IF (D == 1) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) =  0.333333333333333D0
    WTS(1)   =  0.500000000000000D0
    WTS(2)   =  1.500000000000000D0
    
ELSEIF (D == 2) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) = -0.289897948556636D0
    PTS(3,1) =  0.689897948556636D0
    WTS(1)   =  0.222222222222222D0
    WTS(2)   =  1.024971652376843D0
    WTS(3)   =  0.752806125400934D0
    
ELSEIF (D == 3) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) = -0.575318923521694D0
    PTS(3,1) =  0.181066271118531D0
    PTS(4,1) =  0.822824080974592D0
    WTS(1)   =  0.125000000000000D0
    WTS(2)   =  0.657688639960120D0
    WTS(3)   =  0.776386937686344D0        
    WTS(4)   =  0.440924422353536D0

ELSEIF (D == 4) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) = -0.720480271312439D0
    PTS(3,1) = -0.167180864737834D0
    PTS(4,1) =  0.446313972723752D0
    PTS(5,1) =  0.885791607770965D0
    WTS(1)   =  0.080000000000000D0
    WTS(2)   =  0.446207802167142D0
    WTS(3)   =  0.623653045951482D0        
    WTS(4)   =  0.562712030298924D0
    WTS(5)   =  0.287427121582451D0
    
ELSEIF (D == 5) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) = -0.802929828402347D0
    PTS(3,1) = -0.390928546707272D0
    PTS(4,1) =  0.124050379505228D0
    PTS(5,1) =  0.603973164252784D0
    PTS(6,1) =  0.920380285897063D0
    WTS(1)   =  0.055555555555556D0
    WTS(2)   =  0.319640753220511D0
    WTS(3)   =  0.485387188468970D0        
    WTS(4)   =  0.520926783189575D0
    WTS(5)   =  0.416901334311909D0   
    WTS(6)   =  0.201588385253480D0
    
ELSEIF (D == 6) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) = -0.853891342639482D0
    PTS(3,1) = -0.538467724060109D0
    PTS(4,1) = -0.117343037543100D0
    PTS(5,1) =  0.326030619437691D0
    PTS(6,1) =  0.703842800663031D0
    PTS(7,1) =  0.941367145680430D0
    WTS(1)   =  0.040816326530612D0
    WTS(2)   =  0.239227489225312D0
    WTS(3)   =  0.380949873644231D0        
    WTS(4)   =  0.447109829014566D0
    WTS(5)   =  0.424703779005956D0   
    WTS(6)   =  0.318204231467302D0
    WTS(7)   =  0.148988471112020D0
    
ELSEIF (D == 7) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) = -0.887474878926156D0
    PTS(3,1) = -0.639518616526215D0
    PTS(4,1) = -0.294750565773661D0
    PTS(5,1) =  0.094307252661111D0
    PTS(6,1) =  0.468420354430821D0
    PTS(7,1) =  0.770641893678192D0
    PTS(8,1) =  0.955041227122575D0
    WTS(1)   =  0.031250000000000D0
    WTS(2)   =  0.185358154802979D0
    WTS(3)   =  0.304130620646785D0        
    WTS(4)   =  0.376517545389119D0
    WTS(5)   =  0.391572167452494D0   
    WTS(6)   =  0.347014795634501D0
    WTS(7)   =  0.249647901329864D0   
    WTS(8)   =  0.114508814744257D0
    
ELSEIF (D == 8) THEN

    PTS(1,1) = -1.000000000000000D0
    PTS(2,1) = -0.910732089420060D0
    PTS(3,1) = -0.711267485915709D0
    PTS(4,1) = -0.426350485711139D0
    PTS(5,1) = -0.090373369606853D0
    PTS(6,1) =  0.256135670833455D0
    PTS(7,1) =  0.571383041208739D0
    PTS(8,1) =  0.817352784200412D0
    PTS(9,1) =  0.964440169705273D0
    WTS(1)   =  0.024691358024691D0
    WTS(2)   =  0.147654019046315D0
    WTS(3)   =  0.247189378204593D0        
    WTS(4)   =  0.316843775670438D0
    WTS(5)   =  0.348273002772967D0   
    WTS(6)   =  0.337693966975929D0
    WTS(7)   =  0.286386696357231D0   
    WTS(8)   =  0.200553298024551D0
    WTS(9)   =  0.090714504923287D0
    
ELSEIF (D == 9) THEN

    PTS(1,1)  = -1.000000000000000D0
    PTS(2,1)  = -0.927484374233581D0
    PTS(3,1)  = -0.763842042420003D0
    PTS(4,1)  = -0.525646030370079D0
    PTS(5,1)  = -0.236234469390588D0
    PTS(6,1)  =  0.076059197837978D0
    PTS(7,1)  =  0.380664840144724D0
    PTS(8,1)  =  0.647766687674009D0
    PTS(9,1)  =  0.851225220581608D0
    PTS(10,1) =  0.971175180702247D0
    WTS(1)    =  0.020000000000000D0
    WTS(2)    =  0.120296670557482D0
    WTS(3)    =  0.204270131879001D0        
    WTS(4)    =  0.268194837841179D0
    WTS(5)    =  0.305859287724423D0   
    WTS(6)    =  0.313582457226938D0
    WTS(7)    =  0.290610164832918D0   
    WTS(8)    =  0.239193431714380D0
    WTS(9)    =  0.164376012736922D0    
    WTS(10)   =  0.073617005486761D0
  
ENDIF    
    
RETURN
END SUBROUTINE GR_QUADRATURE
!-----------------------------------------------------------------------
!
!                           Subroutine QUADRATURE
!                        Written by Colton J Conroy
!                      @ the Isaac Newton Institute
!                                 12.10.12
!                          
!     This subroutine returns the integration  weights and points 
!     for a given region.  Regions may be one of the
!     following:
!
!      1D:       -1 <= X <= 1
!
!     Edge rules are standard Gauss-Legendre rules. The
!     references for particular rules are noted below.
!
!-----------------------------------------------------------------------
!
!     Input:
!     ------
!       D:      # of points for 1D, Degree of polynomial for 2D
!       REGION:  1D (=1)
!
!     Output:
!     -------
!       PTS: Array of quadrature points of size NPTS by 2
!       WTS: Array of quadrature weights of size NPTS by 1.
!
!-----------------------------------------------------------------------

SUBROUTINE QUADRATURE(D,REGION,PTS,WTS,NPTS,DIM)

USE precisions
IMPLICIT NONE

!.....Declare subroutine input and output

INTEGER(int_p) :: D, REGION, L, I, J, N, NPTS, DIM
REAL(real_p), DIMENSION(NPTS,DIM) :: PTS
REAL(real_p), DIMENSION(NPTS) :: WTS
REAL(real_p), ALLOCATABLE :: A(:), B(:), C(:), M(:), W(:)

!-----------------------------------------------------------------------
!                      1D QUADRATURE RULES
!-----------------------------------------------------------------------
!
!     Notes: (1) These are standard Gauss-Legendre rules.
!            (2) An n point Gauss-Legendre quadrature rule will
!                integrate up to a 2n-1 degree polynomial exactly.
!
!-----------------------------------------------------------------------

IF (REGION == 1) THEN
!-----------------------------------------------------------------------
        IF (D == 1) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = 0.000000000000000D0  ! Points
          WTS(1)   = 1.000000000000000D0  ! Weight
!-----------------------------------------------------------------------
        ELSEIF (D == 2) THEN
!-----------------------------------------------------------------------
          
	      PTS(1,1) = -0.57735026918963D0  ! Points
	      PTS(2,1) = -PTS(1,1)
          WTS(1)   =  1.00000000000000D0  ! Weights
          WTS(2)   =  WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 3) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.77459666924148D0  ! Points
          PTS(2,1) =  0.00000000000000D0
	      PTS(3,1) = -PTS(1,1)
          WTS(1)   =  0.55555555555556D0  ! Weights
          WTS(2)   =  0.88888888888888D0
	      WTS(3)   =  WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 4) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.86113631159405D0  ! Points
          PTS(2,1) = -0.33998104358486D0
	      PTS(3,1) = -PTS(2,1)
	      PTS(4,1) = -PTS(1,1)
          WTS(1)   =  0.34785484513745D0  ! Weights
          WTS(2)   =  0.65214515486255D0
	      WTS(3)   =  WTS(2)
	      WTS(4)   =  WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 5) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.90617984593866D0  ! Points
          PTS(2,1) = -0.53846931010568D0
          PTS(3,1) =  0.00000000000000D0
	      PTS(4,1) = -PTS(2,1)
	      PTS(5,1) = -PTS(1,1)
          WTS(1)   =  0.23692688505619D0  ! Weights
          WTS(2)   =  0.47862867049937D0
          WTS(3)   =  0.56888888888889D0
 	      WTS(4)   =  WTS(2)
	      WTS(5)   =  WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 6) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.93246951420315D0  ! Points
          PTS(2,1) = -0.66120938646626D0
          PTS(3,1) = -0.23861918608320D0
	      PTS(4,1) = -PTS(3,1)
	      PTS(5,1) = -PTS(2,1)
	      PTS(6,1) = -PTS(1,1)
          WTS(1)   =  0.17132449237917D0  ! Weights
          WTS(2)   =  0.36076157304814D0
          WTS(3)   =  0.46791393457269D0
	      WTS(4)   = WTS(3)
	      WTS(5)   = WTS(2)
	      WTS(6)   = WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 7) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.94910791234276D0  ! Points
          PTS(2,1) = -0.74153118559939D0
          PTS(3,1) = -0.40584515137740D0
          PTS(4,1) =  0.D0
	      PTS(5,1) = -PTS(3,1)
	      PTS(6,1) = -PTS(2,1)
	      PTS(7,1) = -PTS(1,1)
          WTS(1)   =  0.12948496616887D0  ! Weights
          WTS(2)   =  0.27970539148928D0
          WTS(3)   =  0.38183005050512D0
          WTS(4)   =  0.41795918367347D0
	      WTS(5)   =  WTS(3)
	      WTS(6)   =  WTS(2)
	      WTS(7)   =  WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 8) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.96028985649754D0  ! Points
          PTS(2,1) = -0.79666647741363D0
          PTS(3,1) = -0.52553240991633D0
          PTS(4,1) = -0.18343464249565D0
	      PTS(5,1) = -PTS(4,1)
	      PTS(6,1) = -PTS(3,1)
	      PTS(7,1) = -PTS(2,1)
	      PTS(8,1) = -PTS(1,1)
          WTS(1)   =  0.10122853629038D0  ! Weights
          WTS(2)   =  0.22238103445337D0
          WTS(3)   =  0.31370664587789D0
          WTS(4)   =  0.36268378337836D0
	      WTS(5)   = WTS(4)
	      WTS(6)   = WTS(3)
	      WTS(7)   = WTS(2)
	      WTS(8)   = WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 9) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.96816023950763D0  ! Points
          PTS(2,1) = -0.83603110732664D0
          PTS(3,1) = -0.61337143270059D0
          PTS(4,1) = -0.32425342340381D0
          PTS(5,1) =  0.00000000000000D0
	      PTS(6,1) = -PTS(4,1)
	      PTS(7,1) = -PTS(3,1)
	      PTS(8,1) = -PTS(2,1)
	      PTS(9,1) = -PTS(1,1)
          WTS(1)   =  0.08127438836163D0  ! Weights
          WTS(2)   =  0.18064816069483D0
          WTS(3)   =  0.26061069640294D0
          WTS(4)   =  0.31234707704000D0
          WTS(5)   =  0.33023935500126D0
	      WTS(6)   =  WTS(4)
	      WTS(7)   =  WTS(3)
	      WTS(8)   =  WTS(2)
	      WTS(9)   =  WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 10) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.97390652851717D0 ! Points
          PTS(2,1) = -0.86506336668898D0
          PTS(3,1) = -0.67940956829902D0
          PTS(4,1) = -0.43339539412925D0
          PTS(5,1) = -0.14887433898163D0
	      PTS(6,1) = -PTS(5,1)
	      PTS(7,1) = -PTS(4,1)
	      PTS(8,1) = -PTS(3,1)
	      PTS(9,1) = -PTS(2,1)
	      PTS(10,1) = -PTS(1,1)
          WTS(1)   =  0.06667134430869D0  ! Weights
          WTS(2)   =  0.14945134915058D0
          WTS(3)   =  0.21908636251598D0
          WTS(4)   =  0.26926671931000D0
          WTS(5)   =  0.29552422471475D0
	      WTS(6)   =  WTS(5)
	      WTS(7)   =  WTS(4)
	      WTS(8)   =  WTS(3)
	      WTS(9)   =  WTS(2)
	      WTS(10)  =  WTS(1)
!-----------------------------------------------------------------------
        ELSEIF (D == 11) THEN
!-----------------------------------------------------------------------
          PTS(1,1) = -0.978228658146072D0  ! Points
          PTS(2,1) = -0.887062599768072D0
          PTS(3,1) = -0.730152005574058D0
          PTS(4,1) = -0.519096129206811D0
          PTS(5,1) = -0.269543155952345D0
          PTS(6,1) =  0.00000000000000D0
	      PTS(7,1) = -PTS(5,1)
	      PTS(8,1) = -PTS(4,1)
	      PTS(9,1) = -PTS(3,1)
	      PTS(10,1) = -PTS(2,1)
	      PTS(11,1) = -PTS(1,1)
          WTS(1)   =  0.055668567116134D0  ! Weights
          WTS(2)   =  0.125580369464930D0
          WTS(3)   =  0.186290210927727D0
          WTS(4)   =  0.233193764591990D0
          WTS(5)   =  0.262804544510247D0
          WTS(6)   =  0.272925086777901D0
	      WTS(7)   =  WTS(5)
	      WTS(8)   =  WTS(4)
	      WTS(9)   =  WTS(3)
	      WTS(10)  =  WTS(2)
	      WTS(11)  =  WTS(1)          
!-----------------------------------------------------------------------
        ELSEIF (D == 12) THEN
!-----------------------------------------------------------------------
          PTS(1,1)  = -0.981560638977557D0 ! Points
          PTS(2,1)  = -0.904117246056125D0
          PTS(3,1)  = -0.769902682704780D0
          PTS(4,1)  = -0.587317950558796D0
          PTS(5,1)  = -0.367831499911997D0
          PTS(6,1)  = -0.125233408390420D0
	      PTS(7,1)  = -PTS(6,1)
	      PTS(8,1)  = -PTS(5,1)
	      PTS(9,1)  = -PTS(4,1)
	      PTS(10,1) = -PTS(3,1)
	      PTS(11,1) = -PTS(2,1)    
	      PTS(12,1) = -PTS(1,1)      
          WTS(1)   =  0.047175326849819D0  ! Weights
          WTS(2)   =  0.106939324308079D0
          WTS(3)   =  0.160078340125383D0
          WTS(4)   =  0.203167417727115D0
          WTS(5)   =  0.233492539742428D0
          WTS(6)   =  0.249147045241822D0
	      WTS(7)   =  WTS(6)
	      WTS(8)   =  WTS(5)
	      WTS(9)   =  WTS(4)
	      WTS(10)  =  WTS(3)
	      WTS(11)  =  WTS(2)    
	      WTS(12)  =  WTS(1)       
!-----------------------------------------------------------------------
        ELSEIF (D == 13) THEN
!-----------------------------------------------------------------------
          PTS(1,1)  = -0.984183135630418D0  ! Points
          PTS(2,1)  = -0.917598204158502D0
          PTS(3,1)  = -0.801578284378796D0
          PTS(4,1)  = -0.642349227505284D0
          PTS(5,1)  = -0.448492791386538D0
          PTS(6,1)  = -0.230458306942889D0
          PTS(7,1)  =  0.00000000000000D0
	      PTS(8,1)  = -PTS(6,1)
	      PTS(9,1)  = -PTS(5,1)
	      PTS(10,1) = -PTS(4,1)
	      PTS(11,1) = -PTS(3,1)
	      PTS(12,1) = -PTS(2,1)    
	      PTS(13,1) = -PTS(1,1)           
          WTS(1)   =  0.040483832138210D0  ! Weights
          WTS(2)   =  0.092121522197771D0
          WTS(3)   =  0.138873686623047D0
          WTS(4)   =  0.178145782612315D0
          WTS(5)   =  0.207816153498783D0
          WTS(6)   =  0.226283148186079D0
          WTS(7)   =  0.232551558746526D0
	      WTS(8)   =  WTS(6)
	      WTS(9)   =  WTS(5)
	      WTS(10)  =  WTS(4)
	      WTS(11)  =  WTS(3)
	      WTS(12)  =  WTS(2)    
	      WTS(13)  =  WTS(1)       
!-----------------------------------------------------------------------
    ELSEIF (D == 14) THEN
!-----------------------------------------------------------------------
          PTS(1,1)  =  -0.986284075744668D0  ! Points
          PTS(2,1)  =  -0.928434167710178D0
          PTS(3,1)  =  -0.827202172187095D0
          PTS(4,1)  =  -0.687292258306539D0
          PTS(5,1)  =  -0.515248966580905D0
          PTS(6,1)  =  -0.319112254056093D0
          PTS(7,1)  =  -0.108054975139041D0
	      PTS(8,1)  = -PTS(7,1)
	      PTS(9,1)  = -PTS(6,1)
	      PTS(10,1) = -PTS(5,1)
	      PTS(11,1) = -PTS(4,1)
	      PTS(12,1) = -PTS(3,1)    
	      PTS(13,1) = -PTS(2,1)  
	      PTS(14,1) = -PTS(1,1)         
          WTS(1)   =  0.035118856655742D0  ! Weights
          WTS(2)   =  0.080158377203965D0
          WTS(3)   =  0.121518968864422D0
          WTS(4)   =  0.157202412024508D0
          WTS(5)   =  0.185539009585654D0
          WTS(6)   =  0.205198166957938D0
          WTS(7)   =  0.215263943525545D0
	      WTS(8)   =  WTS(7)
	      WTS(9)   =  WTS(6)
	      WTS(10)  =  WTS(5)
	      WTS(11)  =  WTS(4)
	      WTS(12)  =  WTS(3)    
	      WTS(13)  =  WTS(2)  
	      WTS(14)  =  WTS(1)
!-----------------------------------------------------------------------
    ELSEIF (D == 15) THEN
!-----------------------------------------------------------------------
          PTS(1,1)  =  -0.987991604145365D0  ! Points
          PTS(2,1)  =  -0.937275532072238D0
          PTS(3,1)  =  -0.848204634377694D0
          PTS(4,1)  =  -0.724418591399470D0
          PTS(5,1)  =  -0.570972125829125D0
          PTS(6,1)  =  -0.394151181172688D0
          PTS(7,1)  =  -0.201194196794467D0
          PTS(8,1)  =   0.000000000000000D0
	      PTS(9,1)  = -PTS(7,1)
	      PTS(10,1) = -PTS(6,1)
	      PTS(11,1) = -PTS(5,1)
	      PTS(12,1) = -PTS(4,1)
	      PTS(13,1) = -PTS(3,1)    
	      PTS(14,1) = -PTS(2,1)  
	      PTS(15,1) = -PTS(1,1)         
          WTS(1)   =  0.030755162349279D0  ! Weights
          WTS(2)   =  0.070366027228910D0
          WTS(3)   =  0.107156874861994D0
          WTS(4)   =  0.139573053648536D0
          WTS(5)   =  0.166268230116555D0
          WTS(6)   =  0.186161004631586D0
          WTS(7)   =  0.198431669034770D0
          WTS(8)   =  0.202578149515106D0
	      WTS(9)   =  WTS(7)
	      WTS(10)  =  WTS(6)
	      WTS(11)  =  WTS(5)
	      WTS(12)  =  WTS(4)
	      WTS(13)  =  WTS(3)    
	      WTS(14)  =  WTS(2)  
	      WTS(15)  =  WTS(1) 	              	           
!-----------------------------------------------------------------------
        ELSEIF (D == 16) THEN
!-----------------------------------------------------------------------
          PTS(1,1)  =  -0.989390690618563D0  ! Points
          PTS(2,1)  =  -0.944601927771250D0
          PTS(3,1)  =  -0.865599829845387D0
          PTS(4,1)  =  -0.755427606961194D0
          PTS(5,1)  =  -0.617864432925138D0
          PTS(6,1)  =  -0.458020961035905D0
          PTS(7,1)  =  -0.281602567183683D0
          PTS(8,1)  =  -0.095012636819780D0
	      PTS(9,1)  = -PTS(8,1)
	      PTS(10,1) = -PTS(7,1)
	      PTS(11,1) = -PTS(6,1)
	      PTS(12,1) = -PTS(5,1)
	      PTS(13,1) = -PTS(4,1)    
	      PTS(14,1) = -PTS(3,1)  
	      PTS(15,1) = -PTS(2,1)   
	      PTS(16,1) = -PTS(1,1)      
          WTS(1)   =  0.027175347944481D0  ! Weights
          WTS(2)   =  0.062244256466172D0
          WTS(3)   =  0.095142131877855D0
          WTS(4)   =  0.124657051969047D0
          WTS(5)   =  0.149574059589351D0
          WTS(6)   =  0.169167214815554D0
          WTS(7)   =  0.182600033491639D0
          WTS(8)   =  0.189451256533539D0
	      WTS(9)   =  WTS(8)
	      WTS(10)  =  WTS(7)
	      WTS(11)  =  WTS(6)
	      WTS(12)  =  WTS(5)
	      WTS(13)  =  WTS(4)    
	      WTS(14)  =  WTS(3)  
	      WTS(15)  =  WTS(2)      
	      WTS(16)  =  WTS(1)  
!-----------------------------------------------------------------------
        ELSE
!-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 25 ********  '
          PRINT*,'  ********        for REGION = 1D          ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

!-----------------------------------------------------------------------        
!                        Shifted Legendre Rules        
!-----------------------------------------------------------------------

    ELSEIF (REGION == 11) THEN

!-----------------------------------------------------------------------    
        IF (D == 1) THEN
!-----------------------------------------------------------------------
        
            PTS(1,1) = 0.5000000000000D0
            WTS(1)   = 1.0000000000000D0
            
!-----------------------------------------------------------------------            
        ELSEIF (D == 2) THEN
!-----------------------------------------------------------------------   
     
            PTS(1,1) = 0.211324865405187D0
            PTS(2,1) = 0.788675134594813D0
            WTS(1)   = 0.500000000000000D0
            WTS(2)   = 0.500000000000000D0
            
!-----------------------------------------------------------------------            
        ELSEIF (D == 3) THEN
!-----------------------------------------------------------------------
        
            PTS(1,1) = 0.112701665379258D0
            PTS(2,1) = 0.500000000000000D0
            PTS(3,1) = 0.887298334620742D0
            WTS(1)   = 0.277777777777778D0
            WTS(2)   = 0.444444444444444D0
            WTS(3)   = 0.277777777777778D0

!----------------------------------------------------------------------            
        ELSEIF (D == 4) THEN
!----------------------------------------------------------------------        
        
            PTS(1,1) = 0.069431844202974D0
            PTS(2,1) = 0.330009478207572D0
            PTS(3,1) = 0.669990521792428D0
            PTS(4,1) = 0.930568155797026D0
            WTS(1)   = 0.173927422568727D0
            WTS(2)   = 0.326072577431273D0
            WTS(3)   = 0.326072577431273D0
            WTS(4)   = 0.173927422568727D0        
            
!----------------------------------------------------------------------            
        ELSEIF(D == 5) THEN
!----------------------------------------------------------------------
            
            PTS(1,1) = 0.046910077030668D0
            PTS(2,1) = 0.230765344947158D0
            PTS(3,1) = 0.500000000000000D0
            PTS(4,1) = 0.769234655052841D0
            PTS(5,1) = 0.953089922969332D0
            WTS(1)   = 0.118463442528095D0
            WTS(2)   = 0.239314335249683D0
            WTS(3)   = 0.284444444444445D0
            WTS(4)   = 0.239314335249683D0   
            WTS(5)   = 0.118463442528095D0         
            
!----------------------------------------------------------------------
        ELSEIF (D == 6) THEN
!----------------------------------------------------------------------

            PTS(1,1) = 0.033765242898424D0
            PTS(2,1) = 0.169395306766868D0
            PTS(3,1) = 0.380690406958402D0
            PTS(4,1) = 0.619309593041598D0
            PTS(5,1) = 0.830604693233132D0
            PTS(6,1) = 0.966234757101576D0
            WTS(1)   = 0.085662246189585D0
            WTS(2)   = 0.180380786524069D0
            WTS(3)   = 0.233956967286345D0
            WTS(4)   = 0.233956967286346D0    
            WTS(5)   = 0.180380786524069D0 
            WTS(6)   = 0.085662246189585D0
            
!---------------------------------------------------------------------
        ELSEIF (D == 7) THEN
!---------------------------------------------------------------------

            PTS(1,1) = 0.025446043828621D0
            PTS(2,1) = 0.129234407200303D0
            PTS(3,1) = 0.297077424311302D0
            PTS(4,1) = 0.500000000000000D0
            PTS(5,1) = 0.702922575688699D0
            PTS(6,1) = 0.870765592799697D0
            PTS(7,1) = 0.974553956171379D0
            WTS(1)   = 0.064742483084435D0
            WTS(2)   = 0.139852695744638D0
            WTS(3)   = 0.190915025252559D0
            WTS(4)   = 0.208979591836735D0   
            WTS(5)   = 0.190915025252559D0 
            WTS(6)   = 0.139852695744639D0      
            WTS(7)   = 0.064742483084435D0        
            
!---------------------------------------------------------------------     
        ELSEIF (D == 8) THEN
!---------------------------------------------------------------------

            PTS(1,1) = 0.019855071751232D0
            PTS(2,1) = 0.101666761293187D0
            PTS(3,1) = 0.237233795041835D0
            PTS(4,1) = 0.408282678752175D0
            PTS(5,1) = 0.591717321247825D0
            PTS(6,1) = 0.762766204958165D0
            PTS(7,1) = 0.898333238706813D0
            PTS(8,1) = 0.980144928248768D0
            WTS(1)   = 0.050614268145188D0
            WTS(2)   = 0.111190517226687D0
            WTS(3)   = 0.156853322938944D0
            WTS(4)   = 0.181341891689181D0   
            WTS(5)   = 0.181341891689181D0
            WTS(6)   = 0.156853322938943D0     
            WTS(7)   = 0.111190517226688D0   
            WTS(8)   = 0.050614268145188D0       
            
!----------------------------------------------------------------------
        ELSEIF (D == 9) THEN
!----------------------------------------------------------------------

            PTS(1,1) = 0.015919880246187D0
            PTS(2,1) = 0.081984446336682D0
            PTS(3,1) = 0.193314283649705D0
            PTS(4,1) = 0.337873288298096D0
            PTS(5,1) = 0.500000000000000D0
            PTS(6,1) = 0.662126711701905D0
            PTS(7,1) = 0.806685716350295D0
            PTS(8,1) = 0.918015553663318D0
            PTS(9,1) = 0.984080119753813D0
            WTS(1)   = 0.040637194180787D0
            WTS(2)   = 0.090324080347429D0
            WTS(3)   = 0.130305348201468D0
            WTS(4)   = 0.156173538520002D0   
            WTS(5)   = 0.165119677500630D0
            WTS(6)   = 0.156173538520001D0   
            WTS(7)   = 0.130305348201468D0    
            WTS(8)   = 0.090324080347428D0    
            WTS(9)   = 0.040637194180787D0   
            
!----------------------------------------------------------------------
        ELSEIF (D == 10) THEN
!----------------------------------------------------------------------

            PTS(1,1)  = 0.013046735741414D0
            PTS(2,1)  = 0.067468316655508D0
            PTS(3,1)  = 0.160295215850488D0
            PTS(4,1)  = 0.283302302935376D0
            PTS(5,1)  = 0.425562830509184D0
            PTS(6,1)  = 0.574437169490816D0
            PTS(7,1)  = 0.716697697064624D0
            PTS(8,1)  = 0.839704784149512D0
            PTS(9,1)  = 0.932531683344492D0
            PTS(10,1) = 0.986953264258586D0
            WTS(1)    = 0.033335672154344D0
            WTS(2)    = 0.074725674575290D0
            WTS(3)    = 0.109543181257991D0
            WTS(4)    = 0.134633359654998D0 
            WTS(5)    = 0.147762112357376D0
            WTS(6)    = 0.147762112357376D0   
            WTS(7)    = 0.134633359654998D0   
            WTS(8)    = 0.109543181257991D0    
            WTS(9)    = 0.074725674575290D0                                             
            WTS(10)   = 0.033335672154344D0
            
!-----------------------------------------------------------------------
        ELSE
!-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 10 ********  '
          PRINT*,'  ********        for REGION = Shifted     ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
          
        ENDIF         
!-----------------------------------------------------------------------
!                        SQUARE QUADRATURE RULES
!-----------------------------------------------------------------------
!
!     All rules have positive weights (with the exception of rule 23,
!     which has 1 negative weight) with all points located inside the
!     triangle (so-called PI rules). Rules marked by * are optimal.  The
!     rest are the best PI rules currently known.
!
!     References:

!     [1] A.H. Stroud, "Approximate calculation of multiple integrals",
!         Prentice-Hall, Englewood Cliffs, N.J., 1971.
!
!     [2] J.W. Wissman and T. Becker, "Partially Symmetric Cubature
!         Formulas for Even Degrees of Exactness", SIAM Journal on
!         Numerical Analysis, 23, 676--685, 1986.
!
!     [3] Construction of Cubature Formulas of Degree Seven and Nine
!         Using Symmetric Planar Regions, Using Orthogonal Polynmials",
!         SIAM Journal on Numerical Analysis, 14, 492--508, 1977.
!
!     [4] I.P. Omelyan and V.B. Solovyan, "Improved cubature formulae of
!         high degrees of exactness for the square", Journal of Compu-
!         tational and Applied Mathematics, 188, 190--204, 2006.
!
!     [5] H.M. Moller, "Minimum-Point Cubature Formula", Numerische
!         Mathematik, 25, 185--200, 1976. (in German)
!
!-----------------------------------------------------------------------

      ELSEIF (REGION == 2) THEN

!-----------------------------------------------------------------------
        IF (D.LE.1) THEN                                      ! Ref [1]*
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1 /)                                ! Multiplicity
          A(1:L) = (/ 0.000000000000000D0 /)              ! Point
          B(1:L) = (/ 0.000000000000000D0 /)
          W(1:L) = (/ 4.000000000000000D0 /)              ! Weight
!-----------------------------------------------------------------------
        ELSEIF (D.LE.2) THEN                                  ! Ref [1]
!-----------------------------------------------------------------------
          L = 3
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1 /)
          A(1:L) = (/ 0.81649658092773D0,          &
                     -0.40824829046386D0,          &      
                    -0.40824829046386D0 /)        
          B(1:L) = (/ 0.00000000000000D0,          &
                      0.70710678118655D0,          &
                     -0.70710678118655D0 /)      
          W(1:L) = (/ 1.33333333333333D0,          &
                      1.33333333333333D0,          &
                      1.33333333333333D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.3) THEN                                  ! Ref [1]*
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4 /)                                ! Multiplicity
          A(1:L) = (/ 0.577350269189626D0 /)              ! Points
          B(1:L) = (/ 0.577350269189626D0 /)
          W(1:L) = (/ 1.000000000000000D0 /)              ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.4) THEN                                  ! Ref [2]*
!-----------------------------------------------------------------------
          L = 6
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1,1,1,1 /)                          ! Multipl. Points
          A(1:L) = (/  0.000000000000000D0,         &
                       0.000000000000000D0,         &
                       0.851914653304601D0,         &
                       0.630912788976754D0,         &
                      -0.851914653304601D0,         &
                     -0.630912788976754D0 /)
          B(1:L) = (/  0.000000000000000D0,         &
                       0.966091783079296D0,         &
                       0.455603727836193D0,         & 
                      -0.731629951573135D0,         & 
                       0.455603727836193D0,         &
                      -0.731629951573135D0 /)
          W(1:L) = (/  1.142857142857143D0,         &                ! Weights
                       0.439560439560440D0,         &
                       0.566072207007532D0,         &
                       0.642719001783677D0,         &
                       0.566072207007532D0,         &
                       0.642719001783677D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.5) THEN                                  ! Ref [1]*
!-----------------------------------------------------------------------
          L = 7
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1, 1, 1, 1, 1, 1, 1 /)
          A(1:L) = (/  0.00000000000000D0,          &
                       0.00000000000000D0,          &
                       0.00000000000000D0,          &
                       0.77459666924148D0,          &
                       0.77459666924148D0,          &
                      -0.77459666924148D0,          &
                      -0.77459666924148D0 /)
          B(1:L) = (/  0.00000000000000D0,          &
                       0.96609178307930D0,          &
                      -0.96609178307930D0,          &
                       0.57735026918963D0,          &
                      -0.57735026918963D0,          &
                       0.57735026918963D0,          &
                      -0.57735026918963D0 /)
          W(1:L) = (/  1.14285714285714D0,          &
                       0.31746031746032D0,          &
                       0.31746031746032D0,          &
                       0.55555555555556D0,          &
                       0.55555555555556D0,          &
                       0.55555555555556D0,          &
                       0.55555555555556D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.6) THEN                                  ! Ref [2]*
!-----------------------------------------------------------------------
          L = 10
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1,1,1,1,1,1,1,1 /)
          A(1:L) = (/    0.000000000000000D0,       &
                         0.000000000000000D0,       &
                         0.888764014654765D0,       &
                         0.604857639464685D0,       & 
                         0.955447506641064D0,       &
                         0.565459993438754D0,       &
                        -0.888764014654765D0,       &
                        -0.604857639464685D0,       &
                        -0.955447506641064D0,       &  
                        -0.565459993438754D0 /)
          B(1:L) =  (/   0.836405633697626D0,       &
                        -0.357460165391307D0,       &
                         0.872101531193131D0,       &
                         0.305985162155427D0,       &
                        -0.410270899466658D0,       &
                        -0.872869311156879D0,       &
                         0.872101531193131D0,       &
                         0.305985162155427D0,       &
                        -0.410270899466658D0,       &
                        -0.872869311156879D0 /)
          W(1:L) =  (/   0.455343245714174D0,       &
                         0.827395973202966D0,       &
                         0.144000884599645D0,       &
                         0.668259104262665D0,       &
                         0.225474004890679D0,       &
                         0.320896396788441D0,       &
                         0.144000884599645D0,       &
                         0.668259104262665D0,       &
                         0.225474004890679D0,       &
                         0.320896396788441D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.7) THEN                                  ! Ref [3]*
!-----------------------------------------------------------------------
          L = 4
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 2,2,2,2 /)                          ! Multiplicity Points
          A(1:L) = (/ 0.91711782231277058626D0,     &      
                      0.61126876646532841440D0,     &
                      0.52942280204265532589D0,     &
                      0.00000000000000000000D0 /)
          B(1:L) = (/ 0.54793120682809232377D0,     &
                      0.93884325665885830459D0,     &
                      0.00000000000000000000D0,     &
                      0.62704137378039531763D0 /)
          W(1:L) = (/ 0.21305721162094912651D0,     &      ! Weights
                      0.17400948894689560610D0,     &
                      0.63585388344327977182D0,     &
                      0.59001271542103076297D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.8) THEN                                  ! Ref [2]
!-----------------------------------------------------------------------
          L = 16
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /)
          A(1:L) = (/   0.000000000000000D0,        &
                        0.000000000000000D0,        &
                        0.000000000000000D0,        &
                        0.000000000000000D0,        &
                        0.639091304900370D0,        & 
                        0.937069076924990D0,        &
                        0.537083530541494D0,        & 
                        0.887188506449625D0,        &
                        0.494698820670197D0,        &
                        0.897495818279768D0,        & 
                       -0.639091304900370D0,        & 
                       -0.937069076924990D0,        &
                       -0.537083530541494D0,        & 
                       -0.887188506449625D0,        &
                       -0.494698820670197D0,        &
                       -0.897495818279768D0 /)
          B(1:L) = (/   0.000000000000000D0,        &
                        0.757629177660505D0,        &
                       -0.236871842255702D0,        &
                       -0.989717929044527D0,        &
                        0.950520955645667D0,        &  
                        0.663882736885633D0,        &
                        0.304210681724104D0,        &  
                       -0.236496718536120D0,        &
                       -0.698953476086564D0,        & 
                       -0.900390774211580D0,        & 
                        0.950520955645667D0,        &  
                        0.663882736885633D0,        &
                        0.304210681724104D0,        &  
                       -0.236496718536120D0,        &
                       -0.698953476086564D0,        & 
                       -0.900390774211580D0 /)
          W(1:L) = (/   0.055364705621440D0,       &     ! Weights
                        0.404389368726076D0,       &
                        0.533546604952635D0,       &  
                        0.117054188786739D0,       & 
                        0.125614417613747D0,       & 
                        0.136544584733588D0,       &
                        0.483408479211257D0,       & 
                        0.252528506429544D0,       &
                        0.361262323882172D0,       &
                        0.085464254086247D0,       &
                        0.125614417613747D0,       & 
                        0.136544584733588D0,       &
                        0.483408479211257D0,       & 
                        0.252528506429544D0,       &
                        0.361262323882172D0,       &
                        0.085464254086247D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.9) THEN                                  ! Ref [5]*
!-----------------------------------------------------------------------
          L = 5
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,4,4,4,4 /)
          A(1:L) = (/ 0.00000000000000000000d0,      &
                      0.96884996636197772072d0,      &  
                      0.75027709997890053354d0,      &
                      0.52373582021442933604d0,      &
                      0.07620832819261717318d0 /)
          B(1:L) = (/ 0.00000000000000000000d0,      &
                      0.63068011973166885417d0,      & 
                      0.92796164595956966740d0,      &
                      0.45333982113564719076d0,      &
                      0.85261572933366230775d0 /)
          W(1:L) = (/ 0.52674897119341563786d0,      &
                      0.08887937817019870697d0,      &
                      0.11209960212959648528d0,      &   
                      0.39828243926207009528d0,      &
                      0.26905133763978080301d0/)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.10) THEN                                 ! Ref []
!-----------------------------------------------------------------------
          L = 22
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = 1
          A(1:L) = (/  4.7324898849276598d-1,     &
                      -3.5072672608918981d-1,     &
                      -4.7113921490701688d-1,     &
                       3.2110023120386515d-2,     &  
                       1.0733227865108741d-1,     & 
                       8.1037492260191812d-1,     &
                      -7.7882541598318511d-1,     & 
                      -6.4763548426267536d-1,     &
                       3.9247487539609610d-1,     &
                       7.6050655071397388d-1,     & 
                      -8.1171510601648733d-1,     &  
                      -1.1429517364223797d-1,     &
                       7.5296563247996007d-1,     & 
                      -6.2402437958984702d-1,     &
                      -2.0650134619887220d-1,     & 
                       5.0891319042960681d-1,     &
                      -9.8171192640479688d-1,     & 
                      -9.4061855719921172d-1,     & 
                      -9.2254816825741193d-1,     &  
                       9.5380192234255112d-1,     &
                       9.6634208368735852d-1,     &
                       9.5774959160007522d-1  /)
          B(1:L) = (/  1.6557852510038315d-1,     &
                       1.8447172062121983d-1,     &
                      -6.6664733059821124d-1,     &
                      -3.1879357593640706d-1,     & 
                       6.1886619139299248d-1,     &
                       6.1159678303492504d-1,     &  
                      -2.1052738914821550d-1,     &
                       6.4749469817525440d-1,     &
                      -7.6311149392438349d-1,     & 
                      -3.6631391678067937d-1,     &
                      -9.2466842429053542d-1,     &
                      -9.4921913140887015d-1,     &
                      -9.7071837396777472d-1,     &  
                       9.8538331193146012d-1,     &
                       9.1195887103573425d-1,     &
                       9.2152907557898267d-1,     &
                      -6.2586619353239670d-1,     & 
                       3.1884535968392919d-1,     &
                       8.7923480439903223d-1,     &
                      -7.5512692061435549d-1,     &
                       1.0431232556636386d-1,     & 
                       9.2621050012583894d-1 /)
          W(1:L) = (/  3.7171764930896173d-1,     &
                       3.8811447402440874d-1,     &
                       2.8395842218278933d-1,     & 
                       4.0824197726154576d-1,     & 
                       3.1498878221231136d-1,     &
                       2.0018320620277513d-1,     &
                       2.6587113477126029d-1,     &
                       2.3057004553370086d-1,     &  
                       2.5827939410342804d-1,     &
                       2.6042397681916851d-1,     &
                       8.8686202216975499d-2,     &
                       1.1861767207465973d-1,     & 
                       5.9110195150351146d-2,     & 
                       3.3387129247073037d-2,     &
                       1.3460977386198075d-1,     &
                       1.3337731192240113d-1,     &  
                       6.3269272761110995d-2,     &
                       1.1984158532391266d-1,     &
                       6.2537941187552140d-2,     & 
                       7.1805898760516657d-2,     &
                       9.7940429484131938d-2,     &
                       3.4467525588983812d-2 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.11) THEN                                 ! Ref [5]*
!-----------------------------------------------------------------------
          L = 6
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4 /)
          A(1:L) = (/ 0.98263922354085547295d0,  &
                      0.82577583590296393730d0,  &
                      0.18858613871864195460d0,  &
                      0.81252054830481310049d0,  & 
                      0.52532025036454776234d0,  &
                      0.41658071912022368274d-1 /)
          B(1:L) = (/ 0.69807610454956756478d0,  &
                      0.93948638281673690721d0,  &
                      0.95353952820153201585d0,  &
                      0.31562343291525419599d0,  &
                      0.71200191307533630655d0,  &
                      0.42484724884866925062d0 /)
          W(1:L) = (/ 0.48020763350723814563d-1,  &
                      0.66071329164550595674d-1,  &
                      0.97386777358668164196d-1,  &
                      0.21173634999894860050d0,  &
                      0.22562606172886338740d0,  &
                      0.35115871839824543766d0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.12) THEN                                 ! Ref [ ]
!-----------------------------------------------------------------------
          L = 31
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = 1
          A(1:L) = (/ -3.8131119459148788d-1,     &
                       6.4573191109453948d-2,     & 
                      -1.0834321329194775d-2,     &
                      -3.1882703020851938d-1,     & 
                       4.9124051440092803d-1,     &  
                      -2.0906799447670715d-1,     &
                       1.8948337189272124d-1,     &
                      -6.9740094114568829d-1,     &
                      -5.8402513589553640d-1,     &
                       4.4701377237977313d-1,     & 
                       3.4867912208378016d-1,     &
                      -8.2562697914640848d-1,     &  
                       8.2839794669248942d-2,     &
                      -4.0351974647584699d-1,     &
                       6.5401042056406233d-1,     &
                      -9.0264377285460362d-1,     & 
                      -8.5295397551927077d-1,     & 
                       8.2216786295901656d-1,     &
                       7.4510555770497011d-1,     &
                      -9.4133986362315170d-1,     & 
                      -5.4461187387612175d-1,     &  
                       7.4407322966908984d-1,     &
                      -7.2254885481492603d-1,     & 
                       6.7524406569352735d-1,     &
                       9.3807676975195520d-1,     &
                      -9.9171408764710778d-1,     &  
                      -9.7391149358748841d-1,     &
                      -9.2231295237113653d-1,     &
                       9.5557708019583965d-1,     &
                       9.8216591295043421d-1,     & 
                       9.2557212335129857d-1  /)
          B(1:L) = (/ -1.5398000784199919d-1,     &
                      -4.9794871568657245d-1,     &
                       1.2224478971657705d-1,     &
                      -7.3294898893865745d-1,     &
                      -7.8409257688727241d-1,     &  
                       6.5611051497298534d-1,     & 
                       8.9911805536738287d-1,     &
                      -4.2248127033807525d-1,     &  
                       3.3100645378716304d-1,     &
                      -2.2659664973677845d-1,     &
                       4.3860404955163790d-1,     & 
                      -7.8339223275580994d-1,     & 
                      -9.4230798695085749d-1,     &  
                       9.6482881096343187d-1,     &
                       7.3857098028675039d-1,     &
                      -1.5881852211379360d-2,     &
                       4.0163251602231481d-1,     & 
                      -5.2660748692256265d-1,     & 
                       1.0988549815212308d-1,     & 
                      -9.5414148142964417d-1,     &
                      -9.5530036379122441d-1,     &
                      -9.7744835460673618d-1,     & 
                       7.9662230649630528d-1,     & 
                       9.9626110485909680d-1,     & 
                       8.8729244026658205d-1,     &
                      -5.2853945477454323d-1,     &
                       5.9955340548238656d-1,     &
                       9.5082592343728289d-1,     &
                      -8.4903353418547223d-1,     &
                      -1.6369711363774472d-1,     &  
                       4.7425880434330858d-1  /)
          W(1:L) = (/  2.3598701292933691d-1,     &
                       2.4110813600847844d-1,     & 
                       2.7048224227639595d-1,     &
                       1.8648166880884143d-1,     &
                       1.7241958393889878d-1,     & 
                       2.4109442092821515d-1,     & 
                       1.4125391453247216d-1,     &
                       1.8292025454973745d-1,     & 
                       2.2222997581864359d-1,     & 
                       2.4904537785770031d-1,     &
                       2.5112362229462498d-1,     &
                       9.2811359998802967d-2,     &
                       9.8238337173087234d-2,     &
                       6.8583817743889233d-2,     &  
                       1.5506938444449642d-1,     & 
                       1.3182320881366066d-1,     & 
                       6.0068518355283269d-2,     &
                       1.5766018519773284d-1,     & 
                       1.9930860602643424d-1,     & 
                       2.1719109379092324d-2,     &
                       6.4596098758824216d-2,     &
                       4.0584519909105145d-2,     &
                       1.2603341007939087d-1,     &
                       3.0006912476756006d-2,     &
                       4.4974753329261395d-2,     &
                       3.8721244577288268d-2,     &
                       4.7718056922379071d-2,     & 
                       3.0056189691189335d-2,     &
                       4.4499282846528827d-2,     &
                       5.6566093223890265d-2,     &
                       9.6814701109561335d-2  /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.13) THEN                                 ! Ref [5]
!-----------------------------------------------------------------------
          L = 9
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,4,4,4,4,4,4,4,4 /)
          A(1:L) = (/ 0.00000000000000000000d0,  &
                      0.77880971155441942252d0,  &
                      0.95729769978630736566d0,  & 
                      0.13818345986246535375d0,  &
                      0.94132722587292523695d0,  &
                      0.47580862521827590507d0,  & 
                      0.75580535657208143627d0,  &
                      0.69625007849174941396d0,  &
                      0.34271655604040678941d0 /)
          B(1:L) = (/ 0.00000000000000000000d0,  &
                      0.98348668243987226379d0,  &
                      0.85955600564163892859d0,  & 
                      0.95892517028753485754d0,  &
                      0.39073621612946100068d0,  &
                      0.85007667369974857597d0,  &
                      0.64782163718701073204d0,  & 
                      0.70741508996444936217d-1,  &
                      0.40930456169403884330d0 /)
          W(1:L) = (/ 0.30038211543122536139d0,  &
                      0.29991838864499131666d-01,  &
                      0.38174421317083669640d-01,  &
                      0.60424923817749980681d-01,  &
                      0.77492738533105339358d-01,  &
                      0.11884466730059560108d0,  &
                      0.12976355037000271129d0,  &
                      0.21334158145718938943d0,  &
                      0.25687074948196783651d0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.15) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 11
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 /)     ! Multipl. Points
          A(1:L) = (/ 0.98798456650771809034922121236542D0,   &
                      0.90815949600657000212015099547736D0,   &
                      0.67928365833453304991325391003325D0,   &
                      0.50911373411758353514778637359508D0,   &
                      0.97675332466910190352385200798077D0,   &
                      0.75619936719149244012005066365912D0,   & 
                      0.89778569328633877480008574782677D0,   &
                      0.20599307074252141729418963873618D0,   &
                      0.45144312511299139017533875886564D0,   & 
                      0.66683824538360873834071399129626D0,   & 
                      0.74295704755765822553432311307323D-1 /)
          B(1:L) = (/ 0.77126821223875533899886933446485D0,   &
                      0.95703183434690690872237176442598D0,   &
                      0.88260197593087253601344445274335D0,   &
                      0.97120312974183699854692313856226D0,   &
                      0.83559862608781647288448365846813D-1,  & 
                      0.75619936719149244012005066365912D0,   &
                      0.46676265923796434848237353795400D0,   &
                      0.84079448454078540426562160968883D0,   &
                      0.56245686233219940637540066298377D0,   &
                      0.19046630243571720761679616635243D0,   &   
                      0.32397702249753019818251854432752D0 /)
          W(1:L) = (/ 0.20881470204497523521771058289754D-1,  &  ! Weights
                      0.25545901574497276542640153395248D-1,  & 
                      0.31203866624933300871149690867662D-1,  &
                      0.38010761595074827467645518285610D-1,  &
                      0.41449061852426148002787373214871D-1,  &
                      0.79320407004083334044710201039891D-1,  &
                      0.88901265758751523303980720079987D-1,  &
                      0.12016982158206027507823569713056D0,   &
                      0.16882043410639799754153511014621D0,   &
                      0.16987162497336185160489786440053D0,   &
                      0.21582538472391594202064661314968D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.17) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 14
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4 /)               ! Multipl.
          A(1:L) = (/ 0.83395914422050762595707520328329D0,    &    ! Points
                      0.96701240760377864958940778776159D0,    &
                      0.98651441086033068583570701960810D0,    &
                      0.17035060808995408160378979789427D-1,   &
                      0.42523896522453030941022314252505D0,    &
                      0.88194109089215356624316729768736D0,    &
                      0.69978907719058600912388980227382D0,    &
                      0.90858412958838344797723547617253D0,    &
                      0.73933759205292015806620475085721D0,    &  
                      0.21246155837885419289305793553062D0,    &
                      0.15644172095846342514322716437912D0,    &
                      0.50398563819427997044830482422504D0,    &
                      0.58243895074467257078462336257058D0,    &
                      0.28970323065541272212557470108155D0 /)
          B(1:L) = (/ 0.99690134998258294114169765276688D0,    &
                      0.92101565015369642619085740576658D0,    &
                      0.56251780244667252352153081833715D0,    &
                      0.98204914256843033449712481990261D0,    &
                      0.95796453236865195741899809390262D0,    &
                      0.74117688569732509081061208741447D0,    &  
                      0.88728024293255774907349992409998D0,    &
                      0.25843122151820770173022364732347D0,    & 
                      0.47295583257297618882242532468520D0,    &
                      0.00000000000000000000000000000000D0,    &
                      0.81259212523912311994061744159953D0,    &
                      0.68201093297792530795370853308359D0,    & 
                      0.11846544560647891209927499369579D0,    &
                      0.40876985953794338411643329876836D0 /)
          W(1:L) = (/ 0.10693483986974526468925667171638D-1,   &     ! Weights
                      0.16771622989325482379964908559201D-1,   &
                      0.21520834803173017585623196363995D-1,   &
                      0.24893201532665059892584476318209D-1,   &
                      0.42463258472030940473501779894230D-1,   &
                      0.53711265037645010830029221647514D-1,   &
                      0.54579479693382460849318673206721D-1,   &
                      0.67375653622461385504403959336007D-1,   &
                      0.98025282885102299426881768287454D-1,   &   
                      0.98325651584666601742691443856667D-1,   &
                      0.10576898319665727200249259468733D0,    &
                      0.10593453574575401283483020818374D0,    &
                      0.14990405484916921950315250394890D0,    &
                      0.15003269160099271050559959853839D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.19) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 17
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 /)          ! Multipl.
          A(1:L) = (/ 0.93790944060174636724164305692148D0,     &    ! Points
                      0.59787966519157168521483636741590D0,     &
                      0.96865781747798834472115427804246D0,     &
                      0.98871713276447330663578670732177D0,     &
                      0.98398534132338386681313796251589D0,     &
                      0.36585775934555586089438861852973D0,     &
                      0.78392149176096602345110934186655D0,     & 
                      0.86937144898957875204792277623522D0,     &
                      0.91889377777573801032232199904313D0,     &
                      0.14078396804456848628363687508206D0,     &
                      0.56379666815446526573028097449771D0,     &
                      0.76200274293070327698230061865243D0,     &
                      0.69745388342191267884442122044750D0,     &
                      0.52310392339404494392757182766006D0,     &
                      0.30566836903929191370033838487568D0,     & 
                      0.44898642628288082765338494532756D0,     &
                      0.98253418759835132805054567507033D-1  /)
          B(1:L) = (/ 0.99998546072852907260205380652647D0,     &
                      0.98732847941781540087026300838966D0,     &
                      0.89837495163532572964949889543009D0,     & 
                      0.62259634389530287767136776393609D0,     &
                      0.81109463487824943619725131193050D-1,    & 
                      0.96462400422891970516096923553526D0,     &
                      0.93756444544378175318295617860873D0,     & 
                      0.74380866034597101765279618918913D0,     & 
                      0.38410875822737883260638757649426D0,     &
                      0.87987957307981017658529064786666D0,     &
                      0.81879210607636015929368719928237D0,     &
                      0.15749920122732269808023356410127D0,     & 
                      0.54657527460177776847334916572487D0,     &   
                      0.00000000000000000000000000000000D0,     &
                      0.65017670270687960549278395428873D0,     &
                      0.34240465380680230945939322837053D0,     &
                      0.23509621629115532325303789785547D0 /)
          W(1:L) = (/ 0.42157189312457273371400997162954D-2,    &       ! Weights
                      0.99237601014741223600089798896376D-2,    &
                      0.15078678879581549295330034135223D-1,    &
                      0.15121496864822956266676863866018D-1,    & 
                      0.23821621047339582724750103846036D-1,    & 
                      0.28746437252189671047629194361029D-1,    &  
                      0.31348715503861464721826917341708D-1,    &
                      0.49762852666717116332578404045926D-1,    &
                      0.55534775159604041101829554407116D-1,    &
                      0.65013710432173970207381031565964D-1,    &
                      0.73819068900731731823980126564345D-1,    &
                      0.90385825968150641564727715959435D-1,    & 
                      0.91881833336425013923433485117405D-1,    &
                      0.94967142638856099564605740671378D-1,    &
                      0.10498174724102843457467466433826D0,     &
                      0.11871336850928058424347422704556D0,     &
                      0.12668324656651729290995285712867D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.21) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 21
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,1 /)
          A(1:L) = (/ 0.99742844318071465788852153329446D0,      &      ! Points
                      0.90092205722857715090631770295789D0,      &
                      0.98137581661152322617866583081353D0,      & 
                      0.47562398846061921636952360049609D0,      & 
                      0.68396217875524370909034162562144D-1,     & 
                      0.96445171170290974019893928634261D0,      &
                      0.96789534067760540982739074113386D0,      &
                      0.71204626183902144323409345553631D0,      &
                      0.86686503563504933746419509606443D0,      & 
                      0.29388895765490255861532312924984D0,      &
                      0.88245289720509646533706783064117D0,      &
                      0.88356465549777630296872390851348D0,      &
                      0.53714330328796591598235477256257D0,      &
                      0.73297737278688849035966507531475D0,      &
                      0.74340815897994389367005251095709D0,      &
                      0.14352379867257862891820869134808D0,      &
                      0.57716374806487034338151001622839D0,      & 
                      0.38010653105519745774291045326124D0,      &
                      0.49380304750704567296572150088764D0,      &
                      0.23300694276964919884562248755586D0,      &
                      0.00000000000000000000000000000000D0 /)
          B(1:L) = (/ 0.52349333540342268677302187698119D0,      &
                      0.99213624611198765984022577926158D0,      &
                      0.93577766500519228442627576803182D0,      &
                      0.98806503981406364167006543205075D0,      &
                      0.98289142121346795894575007552095D0,      &
                      0.73079564192792033198229729840349D0,      & 
                      0.28416943704022564251542547384397D0,      &
                      0.95001737577395264422880843244573D0,      & 
                      0.85754514904426211818642071770359D0,      &
                      0.91478375374420147743639913929870D0,      &
                      0.78310571319347967606817150072422D-1,     &
                      0.50650027819745540446699219481203D0,      &
                      0.82111662321535962574448188332070D0,      &
                      0.67811474157990352881988795983159D0,      & 
                      0.25891212405917942502949658613377D0,      &
                      0.73714284444878233387764988974016D0,      & 
                      0.42605137832252378501649531823938D0,      &
                      0.59190980196005469388929534014808D0,      & 
                      0.46968438389915362845722235027390D-1,     &
                      0.27854054992870057594151995625227D0,      &
                      0.00000000000000000000000000000000D0 /)
          W(1:L) = (/ 0.59245289910274777823163684444026D-2,     &     ! Weights
                      0.63181879993530976749052313712955D-2,     & 
                      0.77654356771885822575525989789117D-2,     &
                      0.13313236524649387539261801099149D-1,     &
                      0.18028023632000303112294832251383D-1,     &    
                      0.20759211084811453614957408545828D-1,     &
                      0.24345594686962939048122829629801D-1,     &
                      0.27291829813157459363923440908365D-1,     &
                      0.31313724489320677669626239379042D-1,     & 
                      0.43028684463625497581941353350492D-1,     & 
                      0.47435107497760585342967184111271D-1,     &   
                      0.49535232791099913534195616930148D-1,     &
                      0.59086609737522776429493659734483D-1,     &
                      0.61598898755573131953283218720393D-1,     &
                      0.68662882737105569078247834904061D-1,     &   
                      0.77511025760449779132628565445376D-1,     &
                      0.77655900861398148704502327066780D-1,     &
                      0.84165159143253389201574044182295D-1,     &
                      0.11664395742356559711226466539621D0,      &
                      0.12587077966701428595400393561138D0,      &
                      0.13498395305263979164774737575571D0 /)
!-----------------------------------------------------------------------
        ELSEIF (D.LE.23) THEN                                 ! Ref [4]
!-----------------------------------------------------------------------
          L = 25
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) =(/4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/)
          A(1:L) = (/  0.22475776435269587875683347524503D0,     & 
                       0.87420793087689680287386901266229D0,     & 
                       0.97602267261364022197431910391898D0,     & 
                       0.98901025758028347071424940975698D0,     &
                       0.99326920924031633712313940841625D0,     &
                       0.56499592316621278816917834165965D0,     &  
                       0.12479933234667809522764885071235D0,     &
                       0.90240439932034384398117518141661D0,     &
                       0.96603416099961004557874518675109D0,     & 
                       0.74327403381257379390317209035852D0,     &
                       0.94397158234337348119479988763974D0,     &   
                       0.36473593323821554671165638117641D0,     &
                       0.90025309287443275420616285586643D0,     & 
                       0.48034665839500246850252725186282D-1,    &
                       0.54919116214309328871373834754894D0,     &
                       0.86487537690899969927500100353866D0,     &
                       0.78769621309348012773202186731599D0,     &
                       0.57671028420945643394316234476328D0,     &
                       0.20013408262762044788586320030931D0,     &
                       0.74602728783032388650354427676695D0,     & 
                       0.64177280458672645848207064078009D0,     &
                       0.40236411309752397363192555539035D0,     &
                       0.46993983120051570507268178090848D0,     &
                       0.22491749438123049571812018957130D0,     &
                       0.22831753386455276245209947736815D0 /)
          B(1:L) = (/  0.29562926123440075779470457695520D0,     &
                       0.99164484170232374329954582252223D0,     &
                       0.96038644389988149213812916422733D0,     & 
                       0.79476170695439095037125111357424D0,     &
                       0.45269906826670511998507618989211D0,     &
                       0.98963443963546680882582623594845D0,     &
                       0.98891634511573107523317221688555D0,     &
                       0.86681747165850322502643568236158D0,     &
                       0.22128837983587450545659391256077D0,     & 
                       0.93739279789140167477396439558006D0,     & 
                       0.62734373811214363563025812471848D0,     &
                       0.93654053553190687038241524198386D0,     &
                       0.25070990341427951813078348104539D-1,    &
                       0.65304687691990057591624038872149D0,     & 
                       0.00000000000000000000000000000000D0,     &
                       0.42021700731140387645536924573814D0,     &
                       0.71486137179637184939863361570348D0,     & 
                       0.82649643106250383189783444260321D0,     & 
                       0.80824267349292209367541999626505D0,     & 
                       0.20581981165646397971556679713472D0,     & 
                       0.51516637687706181290564290685967D0,     &
                       0.64890781819854093908140818020446D0,     & 
                       0.27714799151429808758023024246803D0,     &
                       0.40375503517268762210602535023913D0,     & 
                       0.56116790982355323182358196190726D-1 /)
          W(1:L) = (/ -0.22499144590180737573666435988313D-1,    &
                       0.56757797279709720956600059513997D-2,    &
                       0.61294541632385752564605990154598D-2,    &  
                       0.77307399716921584473393426863691D-2,    & 
                       0.79851030689059360294239183531451D-2,    & 
                       0.10747724909071103372536638574027D-1,    &
                       0.14073722901914630236333688287280D-1,    & 
                       0.22950166656215871772211451659864D-1,    & 
                       0.23104503668552704581666880144042D-1,    &
                       0.24886583460620943481070103801115D-1,    & 
                       0.25495608710781325077549951029551D-1,    &
                       0.35148791457392415926700843453924D-1,    &
                       0.39967612471275640012053163262252D-1,    &
                       0.42210906199290840139296408236284D-1,    & 
                       0.45629308917602135087674549832202D-1,    &
                       0.47184931543983842485535616621711D-1,    &
                       0.48004360511784007777329600464972D-1,    &
                       0.51407132281965816724265717022986D-1,    &
                       0.54210495031621259406088939636629D-1,    & 
                       0.63980226470850872129000560206633D-1,    &
                       0.71999854713967554049756153858067D-1,    &  
                       0.74130368424485001536304396403399D-1,    & 
                       0.87075597857073455631778823034765D-1,    &
                       0.10080008323810791275764368651352D0,     &
                       0.11197008823181576355998539793872D0 /)
!-----------------------------------------------------------------------
        ELSE
!-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 23 ********  '
          PRINT*,'  ********       for REGION = SQUARE       ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

!.......Count up number of quadrature points

        N = 0
        DO I = 1,L
          N = N + M(I)
        ENDDO

!.......Transform quadrature points to master element coordinates

        J = 1
        !DO I = 1, N
        DO I = 1, L
          IF (M(I).EQ.1) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            WTS(J)     =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.2) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) = -A(I)
            PTS(J+1,2) = -B(I)
            WTS(J+1)   =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.4) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) = -B(I)
            PTS(J+1,2) =  A(I)
            WTS(J+1)   =  W(I)
            PTS(J+2,1) = -A(I)
            PTS(J+2,2) = -B(I)
            WTS(J+2)   =  W(I)
            PTS(J+3,1) =  B(I)
            PTS(J+3,2) = -A(I)
            WTS(J+3)   =  W(I)
            J = J + M(I)
            IF (SIZE(M,1) == 1) THEN
                GOTO 10
            ENDIF
          ENDIF
        ENDDO
!-----------------------------------------------------------------------
!                        CUBE QUADRATURE RULES
!-----------------------------------------------------------------------
!
!     A best attempt has been made to find rules which have positive
!     weights with all points located inside the cube, However for some
!     of the higer degrees there have been no such rules developed,
!     consequently some have negative weights and/or points outside the
!     unit cube. Rules marked by ^ are PI, rules marked by * are optimal
!
!     [1] D.A. Dunavant, "Efficient Symmetrical Cubature Rules for
!         Complete Polynomials of High Degree over the Unit Cube",
!         Internatioanl Journal for Numerical Methods in Engineering,
!         23, 397--407, 1986.
!     [2] Kim, Kyoung Joong, & Man Suk Song. "INVARIANT CUBATURE
!         FORMULAS OVER A UNIT CUBE." Comm. Korean Math. Soc 13, No. 4,
!         913--931, 1998.
!     [3] Espelid, Terje O. "On the Construction of Good Fully Symmetric
!         Integration Rules", SIAM Journal on Numerical Analysis,
!         Vol. 24, No. 4, 855-881, 1987.
!
!-----------------------------------------------------------------------

      ELSEIF (REGION == 3) THEN
        
!-----------------------------------------------------------------------
        IF (D.LE.1) THEN                                     ! Ref [1]^*
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1 /)                                ! Multiplicity
          A(1:L) = (/  0.000000000000000D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0 /)
          C(1:L) = (/  0.000000000000000D0 /)
          W(1:L) = (/  8.000000000000000D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.3) THEN                                 ! Ref [1]^*
!-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 6 /)                                ! Multiplicity
          A(1:L) = (/  1.000000000000000D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0 /)
          C(1:L) = (/  0.000000000000000D0 /)
          W(1:L) = (/  1.333333333333333D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.5) THEN                                 ! Ref [1]^*
!-----------------------------------------------------------------------
          L = 2
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 6, 8 /)                             ! Multiplicity
          A(1:L) = (/  0.795822425754221D0,   &
                       0.758786910639328D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0,   &
                       0.758786910639328D0 /)
          C(1:L) = (/  0.000000000000000D0,   &
                       0.758786910639328D0 /)
          W(1:L) = (/  0.886426592797784D0,   &
                       0.335180055401662D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.7) THEN                                  ! Ref [2]^
!-----------------------------------------------------------------------
          L = 3
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 6, 8, 24 /)
          A(1:L) = (/  0.901687807821291D0,    &
                       0.408372221499475D0,    &
                       0.859523090201055D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0,    &
                       0.408372221499475D0,    &
                       0.859523090201055D0 /)
          C(1:L) = (/  0.000000000000000D0,    &
                       0.408372221499475D0,    &
                       0.414735913727988D0 /)
          W(1:L) = (/  0.295189738262623D0,    &
                       0.404055417266202D0,    &
                       0.124850759678944D0 /)             ! Weights
!C-----------------------------------------------------------------------
!C       ELSEIF (D.LE.8) THEN                                 ! Ref [2]^*
!C-----------------------------------------------------------------------
!C          L = 5
!C          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
!C          M(1:L) = (/ 1, 6, 8, 8, 48 /)
!C          A(1:L) = (/  0.000000000000000D0,
!C     &                 0.782460796435947D0,
!C     &                 0.488094669706371D0,
!C     &                 0.862218927661482D0,
!C     &                 0.281113909408340D0 /)             ! Points
!C          B(1:L) = (/  0.000000000000000D0,
!C     &                 0.000000000000000D0,
!C     &                 0.488094669706371D0,
!C    &                 0.862218927661482D0,
!C     &                 0.944196578292009D0 /)
!C          C(1:L) = (/  0.000000000000000D0,
!C     &                 0.000000000000000D0,
!C     &                 0.488094669706371D0,
!C     &                 0.862218927661482D0,
!C     &                 0.697574833707236D0 /)
!C          W(1:L) = (/  0.451903714875209D0,
!C     &                 0.299379177352344D0,
!C     &                 0.300876159371237D0,
!C     &                 0.049484325587704D0,
!C     &                 0.122872389222467D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.9) THEN                                  ! Ref [3]^
!-----------------------------------------------------------------------
          L = 5
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 6, 12, 8, 8, 24 /)
          A(1:L) = (/  0.613681469591708994D0,    &
                       0.877687123257678286D0,    &
                       0.564110807020030054D0,    &
                       0.870099784661975918D0,    &
                       0.432267902630862164D0 /)             ! Points
          B(1:L) = (/  0.000000000000000000D0,    &
                       0.877687123257678286D0,    &
                       0.564110807020030054D0,    &
                       0.870099784661975918D0,    &
                       0.432267902630862164D0 /)
          C(1:L) = (/  0.000000000000000000D0,    &
                       0.000000000000000000D0,    &
                       0.564110807020030054D0,    &
                       0.870099784661975918D0,    &
                       0.938530421864671745D0 /)
          W(1:L) = (/  0.433274995749654543D0,    &
                       0.0917898061361776422D0,   &
                       0.198859838144023500D0,    &
                       0.0501487952993490299D0,   &
                       0.0961168035133733664D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.11) THEN                                 ! Ref [3]^
!-----------------------------------------------------------------------
          L = 7
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 6, 12, 8, 8, 8, 24, 24 /)
          A(1:L) = (/  0.722133038874418502D0,    &
                       0.803933467215284450D0,    &
                       0.280772586651274361D0,    &
                       0.533654008880497079D0,    &
                       0.809488201963098952D0,    &
                       0.980099491009071410D0,    &
                       0.405685980195096387D0 /)             ! Points
          B(1:L) = (/  0.000000000000000000D0,    &
                       0.803933467215284450D0,    &
                       0.280772586651274361D0,    &
                       0.533654008880497079D0,    &
                       0.809488201963098952D0,    &
                       0.980099491009071410D0,    &
                       0.405685980195096387D0 /)
          C(1:L) = (/  0.000000000000000000D0,    &
                       0.000000000000000000D0,    &
                       0.280772586651274361D0,    &
                       0.533654008880497079D0,    &
                       0.809488201963098952D0,    &
                       0.530783831193826387D0,    &
                       0.954583218929566109D0 /)
          W(1:L) = (/  0.236935396192164460D0,    &
                       0.128218254215484766D0,    & 
                       0.135909617156416019D0,    &
                       0.177768064872403869D0,    &
                       0.0620094703001845841D0,   &
                       0.0129505695202585843D0,   &
                       0.0718107368809564271D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.13) THEN                                 ! Ref [1]*
!-----------------------------------------------------------------------
          L = 11
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1, 6, 6, 6, 8, 8, 8, 12, 24, 24, 48 /)
          A(1:L) = (/  0.000000000000000D0,       &
                       0.957697664927095D0,       &
                       0.799168794176284D0,       &
                       0.478190772481903D0,       & 
                       0.863636428036401D0,       &
                       0.625548513457035D0,       &
                       0.367859232104514D0,       &
                       0.869451316547931D0,       &
                       0.767220340241960D0,       &
                       1.012782207207705D0,       &
                       0.944392525641525D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.863636428036401D0,       &
                       0.625548513457035D0,       &
                       0.367859232104514D0,       &
                       0.869451316547931D0,       &
                       0.354778664757558D0,       &
                       1.012782207207705D0,       &
                       0.688378084301308D0 /)
          C(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.863636428036401D0,       &
                       0.625548513457035D0,       &
                       0.367859232104514D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.537923170210472D0,       &
                       0.372427809150424D0 /)
          W(1:L) = (/  0.120671303848159D0,       &
                       0.0827155895731376D0,      &
                      -0.200338143812084D0,       &
                       0.119014821875897D0,       &
                       0.0336821300476000D0,      &
                       0.0991539279234894D0,      &
                       0.0883375530928565D0,      &
                       0.0389329599510302D0,      &
                       0.154603530012075D0,       &
                       0.00380152555202625D0,     &
                       0.0381806114347457D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.15) THEN                                  ! Ref [1]
!-----------------------------------------------------------------------
          L = 14
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,6,6,6,8,8,8,24,24,24,24,24,24,48 /)
          A(1:L) = (/  0.000000000000000D0,       &
                       0.955412599148645D0,       &
                       0.851429576022840D0,       &
                       0.339914995753808D0,       &
                       1.122550051350066D0,       &
                       0.683321227375114D0,       &
                       0.398466735636029D0,       &
                       1.031795356557355D0,       &
                       0.859754056650996D0,       &
                       0.677373166268794D0,       &
                       0.931914101559387D0,       &
                       0.871331764342344D0,       &
                       0.693611775581488D0,       &
                       0.959113457650743D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &   
                       1.122550051350066D0,       &
                       0.683321227375114D0,       &
                       0.398466735636029D0,       &
                       0.870614907692603D0,       &
                       0.417236245841784D0,       &
                       0.234354064403370D0,       &
                       0.931914101559387D0,       &
                       0.871331764342344D0,       &
                       0.693611775581488D0,       &
                       0.622683549553047D0 /)
          C(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       1.122550051350066D0,       &
                       0.683321227375114D0,       & 
                       0.398466735636029D0,       &
                       0.000000000000000D0,       & 
                       0.000000000000000D0,       & 
                       0.000000000000000D0,       &
                       0.710695749151767D0,       &
                       0.277186077584005D0,       &
                       0.258751904925954D0,       &
                       0.406412429958087D0 /)
          W(1:L) = (/  0.0347325185603722D0,      &
                       0.0606208354441015D0,      &
                      -0.0926612667312444D0,      &
                       0.0767301289360290D0,      & 
                       0.0000354673063527669D0,   &
                       0.0613391727954448D0,      &  
                       0.123406122314294D0,       &
                       0.00504204099161268D0,     &
                       0.0727045726671404D0,      &
                       0.0555784362060740D0,      &
                       0.0151521479476415D0,      &
                       0.0270620502572634D0,      & 
                       0.0418922675266216D0,      &
                       0.0208443087896894D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.17) THEN                                 ! Ref [1]*
!-----------------------------------------------------------------------
          L = 17
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,6,6,6,8,8,8,12,12,24,24,24,24,24,24,48,48 /)
          A(1:L) = (/  0.000000000000000D0,       &
                       1.086788874756602D0,       &
                       0.922008663747819D0,       &
                       0.641870173734099D0,       &
                       0.891607017924525D0,       &
                       0.734488336558120D0,       &
                       0.518312376398533D0,       &
                       0.989312892186012D0,       &
                       0.364474560584746D0,       &
                       0.948247130849454D0,       &
                       1.045702906353664D0,       &
                       0.660531051473545D0,       &
                       0.990502583320736D0,       &
                       0.926857250551066D0,       &
                       0.710700783405400D0,       &
                       0.957873006911044D0,       &
                       0.880972079148717D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       & 
                       0.000000000000000D0,       &
                       0.891607017924525D0,       &
                       0.734488336558120D0,       &
                       0.518312376398533D0,       & 
                       0.989312892186012D0,       &
                       0.364474560584746D0,       &
                       0.780127958605308D0,       &
                       0.348818446769418D0,       &
                       0.241437081938889D0,       & 
                       0.990502583320736D0,       &
                       0.926857250551066D0,       & 
                       0.710700783405400D0,       &
                       0.780040197522315D0,       & 
                       0.513930736899189D0 /)
          C(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.891607017924525D0,       &
                       0.734488336558120D0,       &
                       0.518312376398533D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.757226655088057D0,       &  
                       0.315061197287156D0,       &
                       0.208469339137832D0,       &
                       0.527111921367462D0,       & 
                       0.269433491346948D0 /)
          W(1:L) = (/  0.127375023441727D0,       &
                      -0.00641337796713059D0,     &
                       0.0509259994505536D0,      &
                      -0.113461087562205D0,       &
                       0.0111151619355886D0,      &
                       0.0371864272269592D0,      &
                       0.0772186740324679D0,      &
                       0.00213441739735667D0,     &
                       0.0808047120120055D0,      &
                       0.0159035001385071D0,      &
                       0.00765449188677983D0,     &
                       0.0855284062610809D0,      & 
                       0.00197991253834944D0,     &
                       0.00814314291841830D0,     &
                       0.0357632627556923D0,      &
                       0.0155500497623716D0,      &
                       0.0379403443748498D0 /)             ! Weights
!-----------------------------------------------------------------------
        ELSEIF (D.LE.19) THEN                                 ! Ref [1]*
!-----------------------------------------------------------------------
          L = 24
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,6,6,6,6,6,8,8,8,8,12,12,12,24,24,24,24,24,24,24,  &
                                                          24,48,48,48 /)
          A(1:L) = (/  0.000000000000000D0,       &
                       0.950006122416513D0,       &
                       0.914965429121309D0,       &
                       0.822592418605602D0,       &
                       0.652114589201603D0,       &
                       0.157242622070114D0,       &
                       1.001394961101233D0,       &
                       0.949635648752703D0,       &
                       0.594117176975313D0,       & 
                       0.326702572971493D0,       &
                       1.003460985202348D0,       &
                       0.893211074627396D0,       &
                       0.488224948121536D0,       &
                       1.010406331416018D0,       &
                       0.989413303740966D0,       &
                       0.907323294169450D0,       &
                       0.831834985558718D0,       &
                       1.012434400912337D0,       &
                       0.940302135044088D0,       &
                       1.036512025425065D0,       &
                       0.827104008787951D0,       & 
                       0.853891500333261D0,       &
                       0.946965100974417D0,       &
                       0.754378926224437D0 /)             ! Points
          B(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &  
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       1.001394961101233D0,       &
                       0.949635648752703D0,       &
                       0.594117176975313D0,       & 
                       0.326702572971493D0,       &
                       1.003460985202348D0,       &
                       0.893211074627396D0,       &
                       0.488224948121536D0,       &  
                       0.988440552758072D0,       &
                       0.384711778165505D0,       &  
                       0.833136115100468D0,       & 
                       0.220277012103345D0,       &
                       1.012434400912337D0,       &
                       0.940302135044088D0,       &
                       1.036512025425065D0,       &
                       0.827104008787951D0,       &    
                       0.990340675813395D0,       &
                       0.623429573213959D0,       &
                       0.538990774213829D0 /)
          C(1:L) = (/  0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       1.001394961101233D0,       &
                       0.949635648752703D0,       &
                       0.594117176975313D0,       &
                       0.326702572971493D0,       &
                       0.000000000000000D0,       & 
                       0.000000000000000D0,       &
                       0.000000000000000D0,       & 
                       0.000000000000000D0,       & 
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.000000000000000D0,       &
                       0.978322332670673D0,       &
                       0.343148808079936D0,       &
                       0.933501389902834D0,       &
                       0.595273916358421D0,       &
                       0.703401752648627D0,       &
                       0.361644847592954D0,       &
                       0.258107746839302D0 /)
          W(1:L) = (/ -0.296256080038908D0,       &
                      -0.0318213148997343D0,      & 
                       0.0847491973828124D0,      &
                      -0.186007315013108D0,       &
                       0.119234266342962D0,       &
                       0.0927104797995355D0,      &
                       0.0600225327364294D0,      &
                       0.00601997623450650D0,     & 
                       0.0356895101325082D0,      &
                       0.0817900786633523D0,      &
                      -0.0359448438162817D0,      & 
                      -0.0627840019030977D0,      &
                       0.0512421306349611D0,      &
                       0.0209987036073207D0,      & 
                       0.0129343518268887D0,      &
                       0.0508521957924782D0,      &
                       0.0625369819391845D0,      & 
                      -0.0231758216427620D0,      &
                       0.00971367824650409D0,     &
                       0.00266634824440037D0,     & 
                       0.0238739030828427D0,      &
                       0.00449207551990765D0,     &
                       0.0205647395630476D0,      &
                       0.0390081809778401D0/)             ! Weights
!-----------------------------------------------------------------------
        ELSE
!-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 19 ********  '
          PRINT*,'  ********       for REGION = CUBE      ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

!.......Count up number of quadrature points

        N = 0
        DO I = 1,L
          N = N + M(I)
        ENDDO
!.......Transform quadrature points to master element coordinates
        J = 1
        !DO I = 1,N
        DO I = 1, L
          IF (M(I).EQ.1) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            PTS(J,3)   =  C(I)
            WTS(J)     =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.6) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            PTS(J,3)   =  C(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) = -A(I)
            PTS(J+1,2) =  B(I)
            PTS(J+1,3) =  C(I)
            WTS(J+1)   =  W(I)
            PTS(J+2,1) =  B(I)
            PTS(J+2,2) =  A(I)
            PTS(J+2,3) =  C(I)
            WTS(J+2)   =  W(I)
            PTS(J+3,1) =  B(I)
            PTS(J+3,2) = -A(I)
            PTS(J+3,3) =  C(I)
            WTS(J+3)   =  W(I)
            PTS(J+4,1) =  B(I)
            PTS(J+4,2) =  C(I)
            PTS(J+4,3) =  A(I)
            WTS(J+4)   =  W(I)
            PTS(J+5,1) =  B(I)
            PTS(J+5,2) =  C(I)
            PTS(J+5,3) = -A(I)
            WTS(J+5)   =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.8) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            PTS(J,3)   =  C(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) =  A(I)
            PTS(J+1,2) = -B(I)
            PTS(J+1,3) =  C(I)
            WTS(J+1)   =  W(I)
            PTS(J+2,1) =  A(I)
            PTS(J+2,2) =  B(I)
            PTS(J+2,3) = -C(I)
            WTS(J+2)   =  W(I)
            PTS(J+3,1) =  A(I)
            PTS(J+3,2) = -B(I)
            PTS(J+3,3) = -C(I)
            WTS(J+3)   =  W(I)
            PTS(J+4,1) = -A(I)
            PTS(J+4,2) =  B(I)
            PTS(J+4,3) =  C(I)
            WTS(J+4)   =  W(I)
            PTS(J+5,1) = -A(I)
            PTS(J+5,2) = -B(I)
            PTS(J+5,3) =  C(I)
            WTS(J+5)   =  W(I)
            PTS(J+6,1) = -A(I)
            PTS(J+6,2) =  B(I)
            PTS(J+6,3) = -C(I)
            WTS(J+6)   =  W(I)
            PTS(J+7,1) = -A(I)
            PTS(J+7,2) = -B(I)
            PTS(J+7,3) = -C(I)
            WTS(J+7)   =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.12) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            PTS(J,3)   =  C(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) =  A(I)
            PTS(J+1,2) = -B(I)
            PTS(J+1,3) =  C(I)
            WTS(J+1)   =  W(I)
            PTS(J+2,1) = -A(I)
            PTS(J+2,2) =  B(I)
            PTS(J+2,3) =  C(I)
            WTS(J+2)   =  W(I)
            PTS(J+3,1) = -A(I)
            PTS(J+3,2) = -B(I)
            PTS(J+3,3) =  C(I)
            WTS(J+3)   =  W(I)
            PTS(J+4,1) =  C(I)
            PTS(J+4,2) =  A(I)
            PTS(J+4,3) =  B(I)
            WTS(J+4)   =  W(I)
            PTS(J+5,1) =  C(I)
            PTS(J+5,2) =  A(I)
            PTS(J+5,3) = -B(I)
            WTS(J+5)   =  W(I)
            PTS(J+6,1) =  C(I)
            PTS(J+6,2) = -A(I)
            PTS(J+6,3) =  B(I)
            WTS(J+6)   =  W(I)
            PTS(J+7,1) =  C(I)
            PTS(J+7,2) = -A(I)
            PTS(J+7,3) = -B(I)
            WTS(J+7)   =  W(I)
            PTS(J+8,1) =  A(I)
            PTS(J+8,2) =  C(I)
            PTS(J+8,3) =  B(I)
            WTS(J+8)   =  W(I)
            PTS(J+9,1) =  A(I)
            PTS(J+9,2) =  C(I)
            PTS(J+9,3) = -B(I)
            WTS(J+9)   =  W(I)
            PTS(J+10,1) = -A(I)
            PTS(J+10,2) =  C(I)
            PTS(J+10,3) =  B(I)
            WTS(J+10)   =  W(I)
            PTS(J+11,1) = -A(I)
            PTS(J+11,2) =  C(I)
            PTS(J+11,3) = -B(I)
            WTS(J+11)   =  W(I)
            J = J + M(I)
!C            ELSEIF (M(I).EQ.24) THEN !GROUP 4
!C            PTS(J,1)   =  A(I)
!C            PTS(J,2)   =  B(I)
!C            PTS(J,2)   =  C(I)
!C            WTS(J)     =  W(I)
!C            PTS(J+1,1) =  A(I)
!C            PTS(J+1,2) = -B(I)
!C            PTS(J+1,3) =  C(I)
!C            WTS(J+1)   =  W(I)
!C            PTS(J+2,1) = -A(I)
!C            PTS(J+2,2) =  B(I)
!C            PTS(J+2,3) =  C(I)
!C            WTS(J+2)   =  W(I)
!C            PTS(J+3,1) = -A(I)
!C            PTS(J+3,2) = -B(I)
!C            PTS(J+3,3) =  C(I)
!C            WTS(J+3)   =  W(I)
!C            PTS(J+4,1) =  B(I)
!C            PTS(J+4,2) =  A(I)
!C            PTS(J+4,3) =  C(I)
!C            WTS(J+4)   =  W(I)
!C            PTS(J+5,1) =  B(I)
!C            PTS(J+5,2) = -A(I)
!C            PTS(J+5,3) =  C(I)
!C            WTS(J+5)   =  W(I)
!C            PTS(J+6,1) = -B(I)
!C            PTS(J+6,2) =  A(I)
!C            PTS(J+6,3) =  C(I)
!C            WTS(J+6)   =  W(I)
!C            PTS(J+7,1) = -B(I)
!C            PTS(J+7,2) = -A(I)
!C            PTS(J+7,3) =  C(I)
!C            WTS(J+7)   =  W(I)
!C            PTS(J+8,1) =  C(I)
!C            PTS(J+8,2) =  A(I)
!C            PTS(J+8,3) =  B(I)
!C            WTS(J+8)   =  W(I)
!C            PTS(J+9,1) =  C(I)
!C            PTS(J+9,2) =  A(I)
!C            PTS(J+9,3) = -B(I)
!C            WTS(J+9)   =  W(I)
!C            PTS(J+10,1) = C(I)
!C            PTS(J+10,2) = -A(I)
!C            PTS(J+10,3) =  B(I)
!C            WTS(J+10)   =  W(I)
!C            PTS(J+11,1) =  C(I)
!C            PTS(J+11,2) = -A(I)
!C            PTS(J+11,3) = -B(I)
!C            WTS(J+11)   =  W(I)
!C            PTS(J+12,1) =  C(I)
!C            PTS(J+12,2) =  B(I)
!C            PTS(J+12,3) =  A(I)
!C            WTS(J+12)   =  W(I)
!C            PTS(J+13,1) =  C(I)
!C            PTS(J+13,2) =  B(I)
!C            PTS(J+13,3) = -A(I)
!C            WTS(J+13)   =  W(I)
!C            PTS(J+14,1) =  C(I)
!C            PTS(J+14,2) = -B(I)
!C            PTS(J+14,3) =  A(I)
!C            WTS(J+14)   =  W(I)
!C            PTS(J+15,1) =  C(I)
!C            PTS(J+15,2) = -B(I)
!C            PTS(J+15,3) = -A(I)
!C            WTS(J+15)   =  W(I)
!C            PTS(J+16,1) =  A(I)
!C            PTS(J+16,2) =  C(I)
!C            PTS(J+16,3) =  B(I)
!C            WTS(J+16)   =  W(I)
!C            PTS(J+17,1) =  A(I)
!C            PTS(J+17,2) =  C(I)
!C            PTS(J+17,3) = -B(I)
!C            WTS(J+17)   =  W(I)
!C            PTS(J+18,1) = -A(I)
!C            PTS(J+18,2) =  C(I)
!C            PTS(J+18,3) =  B(I)
!C            WTS(J+18)   =  W(I)
!C            PTS(J+19,1) = -A(I)
!C            PTS(J+19,2) =  C(I)
!C            PTS(J+19,3) = -B(I)
!C            WTS(J+19)   =  W(I)
!C            PTS(J+20,1) =  B(I)
!C            PTS(J+20,2) =  C(I)
!C            PTS(J+20,3) =  A(I)
!C            WTS(J+20)   =  W(I)
!C            PTS(J+21,1) =  B(I)
!C            PTS(J+21,2) =  C(I)
!C            PTS(J+21,3) = -A(I)
!C            WTS(J+21)   =  W(I)
!C            PTS(J+22,1) = -B(I)
!C            PTS(J+22,2) =  C(I)
!C            PTS(J+22,3) =  A(I)
!C            WTS(J+22)   =  W(I)
!C            PTS(J+23,1) = -B(I)
!C            PTS(J+23,2) =  C(I)
!C            PTS(J+23,3) = -A(I)
!C            WTS(J+23)   =  W(I)
!C            J = J + M(I)
          ELSEIF (M(I).EQ.24) THEN     !GROUP 5
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            PTS(J,3)   =  C(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) =  A(I)
            PTS(J+1,2) = -B(I)
            PTS(J+1,3) =  C(I)
            WTS(J+1)   =  W(I)
            PTS(J+2,1) = -A(I)
            PTS(J+2,2) =  B(I)
            PTS(J+2,3) =  C(I)
            WTS(J+2)   =  W(I)
            PTS(J+3,1) = -A(I)
            PTS(J+3,2) = -B(I)
            PTS(J+3,3) =  C(I)
            WTS(J+3)   =  W(I)
            PTS(J+4,1) =  B(I)
            PTS(J+4,2) =  A(I)
            PTS(J+4,3) = -C(I)
            WTS(J+4)   =  W(I)
            PTS(J+5,1) =  B(I)
            PTS(J+5,2) = -A(I)
            PTS(J+5,3) = -C(I)
            WTS(J+5)   =  W(I)
            PTS(J+6,1) = -B(I)
            PTS(J+6,2) =  A(I)
            PTS(J+6,3) = -C(I)
            WTS(J+6)   =  W(I)
            PTS(J+7,1) = -B(I)
            PTS(J+7,2) = -A(I)
            PTS(J+7,3) = -C(I)
            WTS(J+7)   =  W(I)
            PTS(J+8,1) =  C(I)
            PTS(J+8,2) =  A(I)
            PTS(J+8,3) =  B(I)
            WTS(J+8)   =  W(I)
            PTS(J+9,1) =  C(I)
            PTS(J+9,2) =  A(I)
            PTS(J+9,3) = -B(I)
            WTS(J+9)   =  W(I)
            PTS(J+10,1) = C(I)
            PTS(J+10,2) = -A(I)
            PTS(J+10,3) =  B(I)
            WTS(J+10)   =  W(I)
            PTS(J+11,1) =  C(I)
            PTS(J+11,2) = -A(I)
            PTS(J+11,3) = -B(I)
            WTS(J+11)   =  W(I)
            PTS(J+12,1) = -C(I)
            PTS(J+12,2) =  B(I)
            PTS(J+12,3) =  A(I)
            WTS(J+12)   =  W(I)
            PTS(J+13,1) = -C(I)
            PTS(J+13,2) =  B(I)
            PTS(J+13,3) = -A(I)
            WTS(J+13)   =  W(I)
            PTS(J+14,1) = -C(I)
            PTS(J+14,2) = -B(I)
            PTS(J+14,3) =  A(I)
            WTS(J+14)   =  W(I)
            PTS(J+15,1) = -C(I)
            PTS(J+15,2) = -B(I)
            PTS(J+15,3) = -A(I)
            WTS(J+15)   =  W(I)
            PTS(J+16,1) =  A(I)
            PTS(J+16,2) =  C(I)
            PTS(J+16,3) =  B(I)
            WTS(J+16)   =  W(I)
            PTS(J+17,1) =  A(I)
            PTS(J+17,2) =  C(I)
            PTS(J+17,3) = -B(I)
            WTS(J+17)   =  W(I)
            PTS(J+18,1) = -A(I)
            PTS(J+18,2) =  C(I)
            PTS(J+18,3) =  B(I)
            WTS(J+18)   =  W(I)
            PTS(J+19,1) = -A(I)
            PTS(J+19,2) =  C(I)
            PTS(J+19,3) = -B(I)
            WTS(J+19)   =  W(I)
            PTS(J+20,1) =  B(I)
            PTS(J+20,2) = -C(I)
            PTS(J+20,3) =  A(I)
            WTS(J+20)   =  W(I)
            PTS(J+21,1) =  B(I)
            PTS(J+21,2) = -C(I)
            PTS(J+21,3) = -A(I)
            WTS(J+21)   =  W(I)
            PTS(J+22,1) = -B(I)
            PTS(J+22,2) = -C(I)
            PTS(J+22,3) =  A(I)
            WTS(J+22)   =  W(I)
            PTS(J+23,1) = -B(I)
            PTS(J+23,2) = -C(I)
            PTS(J+23,3) = -A(I)
            WTS(J+23)   =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.48) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            PTS(J,3)   =  C(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) =  A(I)
            PTS(J+1,2) = -B(I)
            PTS(J+1,3) =  C(I)
            WTS(J+1)   =  W(I)
            PTS(J+2,1) = -A(I)
            PTS(J+2,2) =  B(I)
            PTS(J+2,3) =  C(I)
            WTS(J+2)   =  W(I)
            PTS(J+3,1) = -A(I)
            PTS(J+3,2) = -B(I)
            PTS(J+3,3) =  C(I)
            WTS(J+3)   =  W(I)
            PTS(J+4,1) =  B(I)
            PTS(J+4,2) =  A(I)
            PTS(J+4,3) =  C(I)
            WTS(J+4)   =  W(I)
            PTS(J+5,1) =  B(I)
            PTS(J+5,2) = -A(I)
            PTS(J+5,3) =  C(I)
            WTS(J+5)   =  W(I)
            PTS(J+6,1) = -B(I)
            PTS(J+6,2) =  A(I)
            PTS(J+6,3) =  C(I)
            WTS(J+6)   =  W(I)
            PTS(J+7,1) = -B(I)
            PTS(J+7,2) = -A(I)
            PTS(J+7,3) =  C(I)
            WTS(J+7)   =  W(I)
            PTS(J+8,1) =  A(I)
            PTS(J+8,2) =  B(I)
            PTS(J+8,3) = -C(I)
            WTS(J+8)   =  W(I)
            PTS(J+9,1) =  A(I)
            PTS(J+9,2) = -B(I)
            PTS(J+9,3) = -C(I)
            WTS(J+9)   =  W(I)
            PTS(J+10,1) = -A(I)
            PTS(J+10,2) =  B(I)
            PTS(J+10,3) = -C(I)
            WTS(J+10)   =  W(I)
            PTS(J+11,1) = -A(I)
            PTS(J+11,2) = -B(I)
            PTS(J+11,3) = -C(I)
            WTS(J+11)   =  W(I)
            PTS(J+12,1) =  B(I)
            PTS(J+12,2) =  A(I)
            PTS(J+12,3) = -C(I)
            WTS(J+12)   =  W(I)
            PTS(J+13,1) =  B(I)
            PTS(J+13,2) = -A(I)
            PTS(J+13,3) = -C(I)
            WTS(J+13)   =  W(I)
            PTS(J+14,1) = -B(I)
            PTS(J+14,2) =  A(I)
            PTS(J+14,3) = -C(I)
            WTS(J+14)   =  W(I)
            PTS(J+15,1) = -B(I)
            PTS(J+15,2) = -A(I)
            PTS(J+15,3) = -C(I)
            WTS(J+15)   =  W(I)
            PTS(J+16,1) =  C(I)
            PTS(J+16,2) =  A(I)
            PTS(J+16,3) =  B(I)
            WTS(J+16)   =  W(I)
            PTS(J+17,1) =  C(I)
            PTS(J+17,2) =  A(I)
            PTS(J+17,3) = -B(I)
            WTS(J+17)   =  W(I)
            PTS(J+18,1) =  C(I)
            PTS(J+18,2) = -A(I)
            PTS(J+18,3) =  B(I)
            WTS(J+18)   =  W(I)
            PTS(J+19,1) =  C(I)
            PTS(J+19,2) = -A(I)
            PTS(J+19,3) = -B(I)
            WTS(J+19)   =  W(I)
            PTS(J+20,1) =  C(I)
            PTS(J+20,2) =  B(I)
            PTS(J+20,3) =  A(I)
            WTS(J+20)   =  W(I)
            PTS(J+21,1) =  C(I)
            PTS(J+21,2) =  B(I)
            PTS(J+21,3) = -A(I)
            WTS(J+21)   =  W(I)
            PTS(J+22,1) =  C(I)
            PTS(J+22,2) = -B(I)
            PTS(J+22,3) =  A(I)
            WTS(J+22)   =  W(I)
            PTS(J+23,1) =  C(I)
            PTS(J+23,2) = -B(I)
            PTS(J+23,3) = -A(I)
            WTS(J+23)   =  W(I)
            PTS(J+24,1) = -C(I)
            PTS(J+24,2) =  A(I)
            PTS(J+24,3) =  B(I)
            WTS(J+24)   =  W(I)
            PTS(J+25,1) = -C(I)
            PTS(J+25,2) =  A(I)
            PTS(J+25,3) = -B(I)
            WTS(J+25)   =  W(I)
            PTS(J+26,1) = -C(I)
            PTS(J+26,2) = -A(I)
            PTS(J+26,3) =  B(I)
            WTS(J+26)   =  W(I)
            PTS(J+27,1) = -C(I)
            PTS(J+27,2) = -A(I)
            PTS(J+27,3) = -B(I)
            WTS(J+27)   =  W(I)
            PTS(J+28,1) = -C(I)
            PTS(J+28,2) =  B(I)
            PTS(J+28,3) =  A(I)
            WTS(J+28)   =  W(I)
            PTS(J+29,1) = -C(I)
            PTS(J+29,2) =  B(I)
            PTS(J+29,3) = -A(I)
            WTS(J+29)   =  W(I)
            PTS(J+30,1) = -C(I)
            PTS(J+30,2) = -B(I)
            PTS(J+30,3) =  A(I)
            WTS(J+30)   =  W(I)
            PTS(J+31,1) = -C(I)
            PTS(J+31,2) = -B(I)
            PTS(J+31,3) = -A(I)
            WTS(J+31)   =  W(I)
            PTS(J+32,1) =  A(I)
            PTS(J+32,2) =  C(I)
            PTS(J+32,3) =  B(I)
            WTS(J+32)   =  W(I)
            PTS(J+33,1) =  A(I)
            PTS(J+33,2) =  C(I)
            PTS(J+33,3) = -B(I)
            WTS(J+33)   =  W(I)
            PTS(J+34,1) = -A(I)
            PTS(J+34,2) =  C(I)
            PTS(J+34,3) =  B(I)
            WTS(J+34)   =  W(I)
            PTS(J+35,1) = -A(I)
            PTS(J+35,2) =  C(I)
            PTS(J+35,3) = -B(I)
            WTS(J+35)   =  W(I)
            PTS(J+36,1) =  B(I)
            PTS(J+36,2) =  C(I)
            PTS(J+36,3) =  A(I)
            WTS(J+36)   =  W(I)
            PTS(J+37,1) =  B(I)
            PTS(J+37,2) =  C(I)
            PTS(J+37,3) = -A(I)
            WTS(J+37)   =  W(I)
            PTS(J+38,1) = -B(I)
            PTS(J+38,2) =  C(I)
            PTS(J+38,3) =  A(I)
            WTS(J+38)   =  W(I)
            PTS(J+39,1) = -B(I)
            PTS(J+39,2) =  C(I)
            PTS(J+39,3) = -A(I)
            WTS(J+39)   =  W(I)
            PTS(J+40,1) =  A(I)
            PTS(J+40,2) = -C(I)
            PTS(J+40,3) =  B(I)
            WTS(J+40)   =  W(I)
            PTS(J+41,1) =  A(I)
            PTS(J+41,2) = -C(I)
            PTS(J+41,3) = -B(I)
            WTS(J+41)   =  W(I)
            PTS(J+42,1) = -A(I)
            PTS(J+42,2) = -C(I)
            PTS(J+42,3) =  B(I)
            WTS(J+42)   =  W(I)
            PTS(J+43,1) = -A(I)
            PTS(J+43,2) = -C(I)
            PTS(J+43,3) = -B(I)
            WTS(J+43)   =  W(I)
            PTS(J+44,1) =  B(I)
            PTS(J+44,2) = -C(I)
            PTS(J+44,3) =  A(I)
            WTS(J+44)   =  W(I)
            PTS(J+45,1) =  B(I)
            PTS(J+45,2) = -C(I)
            PTS(J+45,3) = -A(I)
            WTS(J+45)   =  W(I)
            PTS(J+46,1) = -B(I)
            PTS(J+46,2) = -C(I)
            PTS(J+46,3) =  A(I)
            WTS(J+46)   =  W(I)
            PTS(J+47,1) = -B(I)
            PTS(J+47,2) = -C(I)
            PTS(J+47,3) = -A(I)
            WTS(J+47)   =  W(I)
            J = J + M(I)
          ENDIF
        ENDDO
        
!-----------------------------------------------------------------------
      ELSE
!-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 19 ********  '
          PRINT*,'  ********       for REGION = CUBE      ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
      ENDIF
10 RETURN
END SUBROUTINE QUADRATURE
!------------------------------------------------------------------------
!
!                    Subroutine Orthogonal Basis
!                    Written by Colton J. Conroy
!                    @ the Isaac Newton Institute
!                             12.11.12
!
!------------------------------------------------------------------------
SUBROUTINE ORTHOGONAL_BASIS(ELEM, PT, D, DIM, BASIS, DBASIS, BASIS_2D, DBASIS_2D,BASIS_3D, DBASIS_3D)

Use precisions   
USE global, ONLY: FULL3D   
IMPLICIT NONE

!.....Declare subroutine input and output

INTEGER(int_p) :: D, DIM, ELEM
REAL(real_p), DIMENSION(1,DIM) :: PT
REAL(real_p), DIMENSION(D) :: BASIS
REAL(real_p), DIMENSION(D) :: DBASIS
REAL(real_p), DIMENSION(D**2) :: BASIS_2D
REAL(real_p), DIMENSION(D**2,2) :: DBASIS_2D
REAL(real_p), DIMENSION(D**3) :: BASIS_3D
REAL(real_p), DIMENSION(D**3,3) :: DBASIS_3D
      
!.....Explicitly declare remaining variables
      
INTEGER(int_p) :: P, Q, R, RR
REAL(real_p) :: X, Y, Z
REAL(real_p), DIMENSION(0:D-1) :: LEGENDRE_X
REAL(real_p), DIMENSION(0:D-1) :: DLEGENDRE_X
REAL(real_p), DIMENSION(0:D-1) :: LEGENDRE_Y
REAL(real_p), DIMENSION(0:D-1) :: DLEGENDRE_Y
REAL(real_p), DIMENSION(0:D-1) :: LEGENDRE_Z
REAL(real_p), DIMENSION(0:D-1) :: DLEGENDRE_Z
REAL(real_p), DIMENSION(0:D-1,0:D-1) :: PHI
REAL(real_p), DIMENSION(0:D-1,0:D-1,0:D-1) :: PHI3
REAL(real_p), DIMENSION(0:D-1,0:D-1,2) :: DPHI
REAL(real_p), DIMENSION(0:D-1,0:D-1,0:D-1,3) :: DPHI3
REAL(real_p), DIMENSION((D)**2) :: PHI_2D
REAL(real_p), DIMENSION((D)**2,2) :: DPHI_2D
REAL(real_p), DIMENSION((D)**3) :: PHI_3D
REAL(real_p), DIMENSION((D)**3,3) :: DPHI_3D


!.....Check to make sure input is correct

IF (D < 0) THEN
    PRINT*,'  ****** ERROR!!! D must be a positve integer *******'
    PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
    STOP
ELSEIF ((DIM > 3)) THEN
    PRINT*,'  *********** ERROR!!! DIM must be 1, 2, or 3 ***********  '
    PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
    STOP
ELSEIF ((ELEM < 1).OR.(ELEM > 3)) THEN
    IF ((ELEM < 11).OR.(ELEM > 11)) THEN
        PRINT*,'  **** ERROR!!! ELEM must be LINE = 1, TRIA = 2, or QUAD = 3 ****  '
        PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
        STOP
    ENDIF    
ENDIF

!....Assign X, Y, and Z
IF (DIM == 1) THEN
    X = PT(1,1)
    Y = 0.0d0
    Z = 0.0d0
ELSEIF (DIM == 2) THEN
    X = PT(1,1)
    Y = PT(1,2)
    Z = 0.0d0
ELSE
    X = PT(1,1)
    Y = PT(1,2)
    Z = PT(1,3)
ENDIF

!------------------------------------------------------------------------------------
!
! 1D Basis
!
!-----------------------------------------------------------------------------------
IF (ELEM == 1) THEN

    BASIS       = 0.0d0
    DBASIS      = 0.0d0
    LEGENDRE_X  = 0.0d0
    DLEGENDRE_X = 0.0d0
    
    !....D = 0 i.e. constant polynomials
    
    IF (D == 1) THEN
        BASIS = 1.0d0
        DBASIS = 0.0d0
    ENDIF

!....Check to make sure PT is in the master element

    IF (abs(X) > 1) THEN
          PRINT*,' ***** WARNING!! in suroutine orthogonal_basis ***** '
          PRINT*,'       Point is outside of master element!!       '
          PRINT*,'               Execution will continue   '
    ENDIF

!....Legendre polynomials

    LEGENDRE_X(0) = 1.0d0
    DLEGENDRE_X(0) = 0.0d0
    LEGENDRE_X(1) = X
    DLEGENDRE_X(1) = 1.0d0
    IF (D > 2) THEN
        DO P = 1, D-1
           LEGENDRE_X(P+1)  = ((2.0d0*P+1.0d0)*X*LEGENDRE_X(P)-P*LEGENDRE_X(P-1.0d0))/(P+1.0d0)
           DLEGENDRE_X(P+1) = ((2.0d0*P+1.0d0)*(X*DLEGENDRE_X(P)+LEGENDRE_X(P))            &
                            - P*DLEGENDRE_X(P-1.0d0))/(P+1.0d0)
        ENDDO
    ENDIF
    DLEGENDRE_X(0) = 0.d0
    IF (D > 1) THEN
        BASIS = LEGENDRE_X
        DBASIS = DLEGENDRE_X
    ENDIF
        
ELSEIF (ELEM == 11) THEN !....shifted Legendre    
    
    BASIS       = 0.0d0
    DBASIS      = 0.0d0
    LEGENDRE_X  = 0.0d0
    DLEGENDRE_X = 0.0d0    
    
    IF (D == 1) THEN
        BASIS = 1.0d0
        DBASIS = 0.0d0
    ENDIF

!....Check to make sure PT is in the master element

    IF (abs(X) > 1) THEN
          PRINT*,' ***** WARNING!! in suroutine orthogonal_basis ***** '
          PRINT*,'       Point is outside of master element!!       '
          PRINT*,'               Execution will continue   '
    ENDIF

    LEGENDRE_X(0) = 1.0d0
    DLEGENDRE_X(0) = 0.0d0
    LEGENDRE_X(1) = 2.0d0*X - 1.0d0
    DLEGENDRE_X(1) = 2.0d0   
    IF (D == 3) THEN
        LEGENDRE_X(2)  = 6.0d0*X**2-6.0d0*X+1.0d0
        DLEGENDRE_X(2) = 12.0d0*X - 6.0d0
    ELSEIF (D == 4) THEN
        LEGENDRE_X(3)  = 20.0d0*X**3-30.0d0*X**2+12.0d0*X-1.0d0
        DLEGENDRE_X(3) = 60.0d0*X**2-60.0d0*X+12.0d0
    ELSEIF (D == 5) THEN
        LEGENDRE_X(4)  = 70.0d0*X**4 - 140.0d0*X**3 + 90.0d0*X**2 - 20.0d0*X + 1
        DLEGENDRE_X(4) = 280.0d0*X**3 - 420.0d0*X**2 + 180.0d0*X - 20.0d0
    ELSEIF (D > 5) THEN
        PRINT*,'  **** ERROR!!! P shift must be <= 4 ****  '
        PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
        STOP
    ENDIF    
    IF (D > 1) THEN
        BASIS = LEGENDRE_X
        DBASIS = DLEGENDRE_X
    ENDIF    
    
!------------------------------------------------------------------------------------
!
! 2D Quadrilateral Basis
!
!------------------------------------------------------------------------------------
ELSEIF(ELEM == 2) THEN
    
    IF (abs(Y) > 1 .OR. abs(X) > 1) THEN
        PRINT*,' ***** WARNING!! in suroutine orthogonal_basis ***** '
        PRINT*,'       Point is outside of master element!!       '
        PRINT*,'               Execution will continue   '
    ENDIF
    
    BASIS_2D    = 0.0d0
    DBASIS_2D   = 0.0d0
    LEGENDRE_X  = 0.0d0
    DLEGENDRE_X = 0.0d0    
    LEGENDRE_Y  = 0.0d0
    DLEGENDRE_Y = 0.0d0   
    PHI         = 0.0d0
    DPHI        = 0.0d0
    PHI_2D      = 0.0d0
    DPHI_2D     = 0.0d0  
    
!....Legendre Polynomials

    LEGENDRE_X(0) = 1.0d0
    DLEGENDRE_X(0) = 0.0d0
    LEGENDRE_X(1) = X
    DLEGENDRE_X(1) = 1.0d0
    IF (D > 2) THEN
        DO P = 1, D-1
           LEGENDRE_X(P+1)  = ((2.*P+1.)*X*LEGENDRE_X(P)-P*LEGENDRE_X(P-1))/(P+1.0d0)
           DLEGENDRE_X(P+1) = ((2.*P+1.)*(X*DLEGENDRE_X(P)+LEGENDRE_X(P))            &
                            - P*DLEGENDRE_X(P-1))/(P+1.)
        ENDDO
    ENDIF
    DLEGENDRE_X(0) = 0.0d0
    LEGENDRE_Y(0) = 1.0d0
    DLEGENDRE_Y(0) = 0.0d0
    LEGENDRE_Y(1) = Y
    DLEGENDRE_Y(1) = 1.0d0
    IF (D > 2) THEN
        DO P = 1, D-1
           LEGENDRE_Y(P+1)  = ((2.*P+1.)*Y*LEGENDRE_Y(P)-P*LEGENDRE_Y(P-1))/(P+1.0d0)
           DLEGENDRE_Y(P+1) = ((2.*P+1.)*(Y*DLEGENDRE_Y(P)+LEGENDRE_Y(P))            &
                            - P*DLEGENDRE_Y(P-1))/(P+1.)
        ENDDO
    ENDIF
    DLEGENDRE_Y(0) = 0.0d0
!....Basis Functions (Full Tensor Basis)

    DO P = 0, D-1
        DO Q = 0, D-1
            PHI(Q,P)    = LEGENDRE_X(Q)*LEGENDRE_Y(P)
!....DPHI with respect to X
            DPHI(Q,P,1) = DLEGENDRE_X(Q)*LEGENDRE_Y(P)
!....DPHI with respect to Y
            DPHI(Q,P,2) = LEGENDRE_X(Q)*DLEGENDRE_Y(P)  
        ENDDO
    ENDDO
    
!....Re-order basis (1 + eta + zeta + eta*zeta for P = 1)    
    
    R = 1
    DO Q = 0, D-1
        DO P = 0, D-1
            PHI_2D(R)    = PHI(P,Q)
            DPHI_2D(R,1) = DPHI(P,Q,1)
            DPHI_2D(R,2) = DPHI(P,Q,2)
            IF (ABS(PHI_2D(R)) < 1.0D-15) THEN
                PHI_2D(R) = 0.0d0
            ENDIF
            IF (ABS(DPHI_2D(R,1)) < 1.0D-15) THEN
                DPHI_2D(R,1) = 0.0d0
            ENDIF
            IF (ABS(DPHI_2D(R,2)) < 1.0D-15) THEN
                DPHI_2D(R,2) = 0.0d0
            ENDIF
            R = R + 1
        ENDDO
    ENDDO
    
    BASIS_2D  = PHI_2D
    DBASIS_2D = DPHI_2D
    
!----------------------------------------------------------------------------------
!
!   3D Hexahedral Basis 8/1/13 
!
!-----------------------------------------------------------------------------------
ELSEIF (ELEM == 3) THEN

    IF (abs(Y) > 1 .OR. abs(X) > 1 .OR. abs(Z) > 1) THEN
        PRINT*,' ***** WARNING!! in suroutine orthogonal_basis ***** '
        PRINT*,'       Point is outside of master element!!       '
        PRINT*,'               Execution will continue   '
    ENDIF

     BASIS_3D  = 0.0d0
     DBASIS_3D = 0.0d0
     LEGENDRE_X = 0.0d0
     DLEGENDRE_X = 0.0d0
     LEGENDRE_Y = 0.0d0
     DLEGENDRE_Y = 0.0d0
     LEGENDRE_Z = 0.0d0
     DLEGENDRE_Z = 0.0d0   
     PHI_3D = 0.0d0
     DPHI_3D = 0.0d0
     PHI3 = 0.0d0
     DPHI3 = 0.0d0
       
!....Legendre Polynomials

    LEGENDRE_X(0) = 1.0d0
    DLEGENDRE_X(0) = 0.0d0
    LEGENDRE_X(1) = X
    DLEGENDRE_X(1) = 1.0d0
    IF (D > 2) THEN
        DO P = 1, D-1
           LEGENDRE_X(P+1)  = ((2.*P+1.)*X*LEGENDRE_X(P)-P*LEGENDRE_X(P-1))/(P+1.0d0)
           DLEGENDRE_X(P+1) = ((2.*P+1.)*(X*DLEGENDRE_X(P)+LEGENDRE_X(P))            &
                            - P*DLEGENDRE_X(P-1))/(P+1.)
        ENDDO
    ENDIF
    DLEGENDRE_X(0) = 0.0d0
    
    LEGENDRE_Y(0) = 1.0d0
    DLEGENDRE_Y(0) = 0.0d0
    LEGENDRE_Y(1) = Y
    DLEGENDRE_Y(1) = 1.0d0
    IF (D > 2) THEN
        DO P = 1, D-1
           LEGENDRE_Y(P+1)  = ((2.*P+1.)*Y*LEGENDRE_Y(P)-P*LEGENDRE_Y(P-1))/(P+1.0d0)
           DLEGENDRE_Y(P+1) = ((2.*P+1.)*(Y*DLEGENDRE_Y(P)+LEGENDRE_Y(P))            &
                            - P*DLEGENDRE_Y(P-1))/(P+1.)
        ENDDO
    ENDIF
    DLEGENDRE_Y(0) = 0.0d0  
    
    LEGENDRE_Z(0) = 1.0d0
    DLEGENDRE_Z(0) = 0.0d0
    LEGENDRE_Z(1) = Z
    DLEGENDRE_Z(1) = 1.0d0
    IF (D > 2) THEN
        DO P = 1, D-1
           LEGENDRE_Z(P+1)  = ((2.*P+1.)*Z*LEGENDRE_Z(P)-P*LEGENDRE_Z(P-1))/(P+1.0d0)
           DLEGENDRE_Z(P+1) = ((2.*P+1.)*(Z*DLEGENDRE_Z(P)+LEGENDRE_Z(P))            &
                            - P*DLEGENDRE_Z(P-1))/(P+1.)
        ENDDO
    ENDIF
    DLEGENDRE_Z(0) = 0.0d0
    
     IF (FULL3D == 1) THEN
     
!....Basis Functions (Full Tensor Basis)
        DO R = 0, D-1 
            DO P = 0, D-1
                DO Q = 0, D-1
                    PHI3(Q,P,R)    = LEGENDRE_X(Q)*LEGENDRE_Y(P)*LEGENDRE_Z(R)
!....DPHI with respect to X
                    DPHI3(Q,P,R,1) = DLEGENDRE_X(Q)*LEGENDRE_Y(P)*LEGENDRE_Z(R)
!....DPHI with respect to Y
                    DPHI3(Q,P,R,2) = LEGENDRE_X(Q)*DLEGENDRE_Y(P)*LEGENDRE_Z(R)
!....DPHI with respect to Z
                    DPHI3(Q,P,R,3) = LEGENDRE_X(Q)*LEGENDRE_Y(P)*DLEGENDRE_Z(R)  
                ENDDO
            ENDDO    
        ENDDO    
    
!....Re-order basis (1 + eta + xi + eta*xi + zeta + zeta*eta + zeta*xi + zeta*xi*eta for P = 1)    
    
        RR = 1
        DO R = 0, D-1 
            DO Q = 0, D-1
                DO P = 0, D-1
                    PHI_3D(RR)    = PHI3(P,Q,R)
                    DPHI_3D(RR,1) = DPHI3(P,Q,R,1)
                    DPHI_3D(RR,2) = DPHI3(P,Q,R,2)
                    DPHI_3D(RR,3) = DPHI3(P,Q,R,3)
                    IF (ABS(PHI_3D(RR)) < 1.0D-15) THEN
                        PHI_3D(RR) = 0.0d0
                    ENDIF
                    IF (ABS(DPHI_3D(RR,1)) < 1.0D-15) THEN
                        DPHI_3D(RR,1) = 0.0d0
                    ENDIF
                    IF (ABS(DPHI_3D(RR,2)) < 1.0D-15) THEN
                        DPHI_3D(RR,2) = 0.0d0
                    ENDIF
                    IF (ABS(DPHI_3D(RR,3)) < 1.0D-15) THEN
                        DPHI_3D(RR,3) = 0.0d0
                    ENDIF            
                    RR = RR + 1
                ENDDO
            ENDDO    
        ENDDO
    
    ELSEIF (FULL3D == 0) THEN
    
        RR = 1
        DO R = 0, D-1
            DO Q = 0, D-1-R
                DO P = 0, D-1-Q-R
                    PHI_3D(RR)    = LEGENDRE_X(Q)*LEGENDRE_Y(P)*LEGENDRE_Z(R)
                    DPHI_3D(RR,1) = DLEGENDRE_X(Q)*LEGENDRE_Y(P)*LEGENDRE_Z(R)
                    DPHI_3D(RR,2) = LEGENDRE_X(Q)*DLEGENDRE_Y(P)*LEGENDRE_Z(R)
                    DPHI_3D(RR,3) = LEGENDRE_X(Q)*LEGENDRE_Y(P)*DLEGENDRE_Z(R)
                    IF (ABS(PHI_3D(RR)) < 1.0D-15) THEN
                        PHI_3D(RR) = 0.0d0
                    ENDIF
                    IF (ABS(DPHI_3D(RR,1)) < 1.0D-15) THEN
                        DPHI_3D(RR,1) = 0.0d0
                    ENDIF
                    IF (ABS(DPHI_3D(RR,2)) < 1.0D-15) THEN
                        DPHI_3D(RR,2) = 0.0d0
                    ENDIF
                    IF (ABS(DPHI_3D(RR,3)) < 1.0D-15) THEN
                        DPHI_3D(RR,3) = 0.0d0
                    ENDIF            
                    RR = RR + 1
                ENDDO
            ENDDO
        ENDDO
        
    ENDIF    
                        
    BASIS_3D  = PHI_3D
    DBASIS_3D = DPHI_3D    
           
ENDIF

RETURN
END SUBROUTINE ORTHOGONAL_BASIS
!-----------------------------------------------------------------------------------
!
!                 Hot Start Subroutine
!              Written by Colton J. Conroy
!                    @ the C.H.I.L
!
!------------------------------------------------------------------------------------
SUBROUTINE H_START(TIME,REGION)

USE precisions
USE global
IMPLICIT NONE

COMPLEX(16) :: TIME, xCOORD, zCOORD, yCOORD
INTEGER(int_p) :: REGION
REAL(real_p), DIMENSION(L2_2D) :: L2H, Kuse, Nuse
REAL(real_p), DIMENSION(4) :: globalNODES
REAL(real_p), DIMENSION(2) :: Xnode, Znode, Ynode
REAL(real_p) :: grad_SE, Vu, SE, Vavg, xpt, ypt, Huse, k, nn
REAL(real_p) :: zpt, solnT, Vv, Vw, Zb, DSEx, DSEy, Uavg, Fpce
REAL(real_p) :: Fxm, Fym, Fw, FxmA, FymA

solnT = REAL(TIME)
L2H   = 0.0d0
Kuse  = 0.0d0
Nuse  = 0.0d0
Xnode = 0.0d0
Ynode = 0.0d0
Znode = 0.0d0
globalNODES = 0.0d0

!....Read in degrees of freedom

OPEN(10,FILE='ZETA_DOF.2d')
DO i = 1, Qelems
    DO j = 1, NDOF
        READ(10,*)ZETA(j,i,1)
    ENDDO    
ENDDO
CLOSE(10)

IF (REGION == 1) THEN
    OPEN(11,FILE='U_DOF.1d')
    DO i = 1, nelemsX
        DO j = 1, LE
            DO kk = 1, nelemsZ
                DO ii = 1, pu+1
                    READ(11,*)U(ii,kk,i,j,1)
                ENDDO    
            ENDDO
        ENDDO
    ENDDO
    CLOSE(11)
ELSEIF (REGION == 2) THEN
    OPEN(11,FILE='U_DOF.2d')
    DO i = 1, Qelems
        DO j = 1, NDOF
            READ(11,*)U_2D(j,i,1)
        ENDDO
    ENDDO
    CLOSE(11)
ELSEIF (REGION == 3) THEN
    OPEN(11,FILE='U_DOF.3d')
    DO i = 1, HEXelems
        DO j = 1, NDOF_3D
            READ(11,*)U_3D(j,i,1)
        ENDDO    
    ENDDO
    CLOSE(11)
ENDIF

OPEN(12,FILE='Uavg_DOF.2d')
DO i = 1, Qelems
    DO j = 1, NDOF
        READ(12,*)USE(j,i,1)
    ENDDO    
ENDDO
CLOSE(12)

!....L2 projection of H, N, k

DO i = 1, Qelems
    
    globalNODES = CONN(i,:)
    Xnode(1) = Qnode(globalNODES(1),1)
    Xnode(2) = Qnode(globalNODES(3),1)
    Ynode(1) = Qnode(globalNODES(1),2)
    Ynode(2) = Qnode(globalNODES(2),2)
       
    DO j = 1, L2_2D
        
        xpt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,1))*Xnode(1)+(1.0d0+L2PTSz(j,1))*Xnode(2))
        ypt = (1.0d0/2.0d0)*((1.0d0-L2PTSz(j,2))*Ynode(1)+(1.0d0+L2PTSz(j,2))*Ynode(2))
        xCOORD = dcmplx(xpt,0.0d0)
        
        IF (PROB == 1) THEN
        
            Hexact(j,i) = PROB_DATA%Hslope*xpt**2
            zCOORD = dcmplx(Hexact(j,i)*Z(1), 0.0d0)
            Huse = Hexact(j,i)
            CALL LG3D_SOLUTIONS_H(xCOORD,zCOORD,TIME,grad_SE,Vu,SE,Vavg,Huse,nn,k)
            Kuse(j) = k
            Nuse(j) = nn        
            L2H(j)  = Huse    
           
        ELSEIF (PROB == 3) THEN 

            zpt = Z(nelemsZ)
            CALL NL_SMOOTH_3D(xpt,ypt,zpt,solnT,Vu,Vv,Vw,Huse,Zb,SE,DSEx,DSEy,Uavg,Vavg,Fpce,Fxm,Fym,Fw,FxmA,FymA,PROB_DATA)
            Kuse(j) = PROB_DATA%k
            Nuse(j) = PROB_DATA%N0
            L2H(j)  = Zb 

        ENDIF
                          
    ENDDO
        
    Hh(:,i) = MATMUL(L2U,L2H)    
    NX(:,i) = MATMUL(L2U,Nuse)
    KX(:,i) = MATMUL(L2U,Kuse)     
        
ENDDO

RETURN
END SUBROUTINE
!-----------------------------------------------------------------------------------
!
!                 Output Subroutine
!              Written by Colton J. Conroy
!                   @ the C.H.I.L
!
!------------------------------------------------------------------------------------
SUBROUTINE GLOBAL_OUTPUT(REGION)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: REGION
INTEGER(int_p), PARAMETER :: fine = 0
REAL(real_p), DIMENSION(L2_3D) :: W_L2

W_L2 = 0.0d0

!....Write out degrees of freedom

IF (HOT_START == 1) THEN
    OPEN(10,FILE='ZETA_DOF2.2d')
    DO i = 1, Qelems
        WRITE(10,*)ZETA(:,i,1)
    ENDDO
    CLOSE(10)
ELSE
    OPEN(10,FILE='ZETA_DOF.2d')
    DO i = 1, Qelems
        IF (NON_LINEAR == 1) THEN
            WRITE(10,*)ZETA_2D(:,i,1)
        ELSE
            WRITE(10,*)ZETA(:,i,1)
        ENDIF    
    ENDDO
    CLOSE(10)
ENDIF    
!
!IF (REGION == 1) THEN
!    OPEN(11,FILE='U_DOF.1d')
!    DO i = 1, nelemsX
!        DO j = 1, LE
!            DO kk = 1, nelemsZ
!                WRITE(11,*)U(:,kk,i,j,1)
!            ENDDO
!        ENDDO
!    ENDDO
!    CLOSE(11)
!ELSEIF (REGION == 2) THEN

IF (HOT_START == 1) THEN
    OPEN(11,FILE='U_DOF2.3d')
    DO i = 1, HEXelems
        DO j = 1, NDOF_3D
            WRITE(11,*)U_3D(j,i,1)
        ENDDO    
    ENDDO
    CLOSE(11)
ELSE    
    OPEN(11,FILE='U_DOF.3d')
    DO i = 1, HEXelems
        DO j = 1, NDOF_3D
            WRITE(11,*)U_3D(j,i,1)
        ENDDO    
    ENDDO
    CLOSE(11)
ENDIF    

IF (HOT_START == 1) THEN
    OPEN(12,FILE='Uavg_DOF2.2d')
    DO i = 1, Qelems
        WRITE(12,*)USE(:,i,1)
    ENDDO
    CLOSE(12)
ELSE    
    OPEN(12,FILE='Uavg_DOF.2d')
    DO i = 1, Qelems
        WRITE(12,*)USE(:,i,1)
    ENDDO
    CLOSE(12)
ENDIF

IF (HOT_START == 1) THEN
    OPEN(13,FILE='V_DOF2.3d')
    DO i = 1, HEXelems
        WRITE(13,*)V_3D(:,i,irk)
    ENDDO
    CLOSE(13)
ELSE
    OPEN(13,FILE='V_DOF.3d')
    DO i = 1, HEXelems
        WRITE(13,*)V_3D(:,i,irk)
    ENDDO
    CLOSE(13)   
ENDIF 

IF (HOT_START == 1) THEN
    OPEN(12,FILE='Vavg_DOF2.2d')
    DO i = 1, Qelems
        WRITE(12,*)VSE(:,i,1)
    ENDDO
    CLOSE(12)
ELSE    
    OPEN(12,FILE='Vavg_DOF.2d')
    DO i = 1, Qelems
        WRITE(12,*)VSE(:,i,1)
    ENDDO
    CLOSE(12)
ENDIF


IF (HOT_START == 1) THEN
    OPEN(13,FILE='W_DOF2.3d')
    DO i = 1, HEXelems
        WRITE(13,*)W_3D(:,i)
    ENDDO
    CLOSE(13)
ELSE
    OPEN(13,FILE='W_DOF.3d')
    DO i = 1, HEXelems
        WRITE(13,*)W_3D(:,i)
    ENDDO
    CLOSE(13)   
ENDIF      


IF (fine == 1) THEN
    OPEN(14,FILE='W_EX.3d')
    DO i = 1, HEXelems
        W_L2 = MATMUL(L2PHI_3D,W_3D(:,i)) 
        DO j = 1, L2_3D
            WRITE(14,*)W_L2(j)
        ENDDO
    ENDDO
    CLOSE(14)
ENDIF      

IF (HOT_START == 1) THEN
    OPEN(15,FILE='DW_DOF2.3d')
    DO i = 1, HEXelems
        WRITE(15,*)DQW(:,i)
    ENDDO
    CLOSE(15)
ELSE
    OPEN(15,FILE='DW_DOF.3d')
    DO i = 1, HEXelems
        WRITE(15,*)DQW(:,i)
    ENDDO
    CLOSE(15)   
ENDIF  

IF (HOT_START == 1) THEN
    OPEN(16,FILE='RHO_DOF2.3d')
    DO i = 1, HEXelems
        WRITE(16,*)RHO(:,i)
    ENDDO
    CLOSE(16)
ELSE
    OPEN(16,FILE='RHO_DOF.3d')
    DO i = 1, HEXelems
        WRITE(16,*)RHO(:,i)
    ENDDO
    CLOSE(16)   
ENDIF 

IF (HOT_START == 1) THEN
    OPEN(16,FILE='R_DOF2.3d')
    DO i = 1, HEXelems
        WRITE(16,*)Rh(:,i)
    ENDDO
    CLOSE(16)
ELSE
    OPEN(16,FILE='R_DOF.3d')
    DO i = 1, HEXelems
        WRITE(16,*)Rh(:,i)
    ENDDO
    CLOSE(16)   
ENDIF 
                
        
IF (BASIS_OUTPUT == 1) THEN        
    OPEN(17,FILE='2DelemPTS.3d') 
    WRITE(17,*)NPTS_2D
    DO i = 1, NPTS_2D
        DO j = 1, 2
            WRITE(17,*)elemPTSz(i,j)
        ENDDO
    ENDDO        
ENDIF

IF (HOT_START == 1) THEN
    OPEN(18,FILE='N_DOF2.2d')
    DO i = 1, Qelems
        WRITE(18,*)Nx(:,i)
    ENDDO
    CLOSE(18)
ELSE
    OPEN(18,FILE='N_DOF.2d')
    DO i = 1, Qelems
        WRITE(18,*)Nx(:,i)
    ENDDO
    CLOSE(18)   
ENDIF 

RETURN
END SUBROUTINE
!------------------------------------------------------------------------------------
!
!                 Subroutine RK Coefficients
!                 Written by Colton J. Conroy
!                 @ the Isaac Newton Institute
!                           12.12.12
!
!-----------------------------------------------------------------------------------
SUBROUTINE RK_COEFFICIENTS(ORDER,S,NEW)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ORDER, S, PP, SSP1, SSP2, nrk, NEW
REAL(real_p) :: csum
REAL(real_p), DIMENSION(S,S) :: crk
REAL(real_p), ALLOCATABLE :: ATVD(:,:), BTVD(:,:), DTVD(:)

!....Constants

PP = ORDER - 1
SSP1 = S
SSP2 = ORDER

!....Allocate matricies

ALLOCATE(alpha(S,S))
ALLOCATE(beta(S,S))
ALLOCATE(ct(S))

ALLOCATE(ATVD(S,S))
ALLOCATE(BTVD(S,S))
ALLOCATE(DTVD(S))

!....Initialize matricies

alpha = 0.0d0
beta  = 0.0d0
crk   = 0.0d0
ct    = 0.0d0

!....Check input

IF (S < 1 .OR. S > 8) THEN
    PRINT*,'  ****** ERROR!!! Number of stages must be > 1 and < 8 *******'
    PRINT*,'  Execution terminated in subroutine RK_coefficients  '
    STOP
ENDIF

IF (ORDER < 1 .OR. ORDER > 4) THEN
    PRINT*,'  ****** ERROR!!! Order cannot be < 1 or > 4 *******'
    PRINT*,'  Execution terminated in subroutine RK_coefficients  '
    STOP
ENDIF

IF (S < ORDER) THEN
    PRINT*,'  ****** ERROR!!! Stages must be >= order (> for Order 4) *******'
    PRINT*,'  Execution terminated in subroutine RK_coefficients  '
    STOP
ENDIF

IF (ORDER == 4) THEN
    IF (S <= ORDER) THEN
        PRINT*,'  ****** ERROR!!! Stages must be > order for 4th Order *******'
        PRINT*,'  Execution terminated in subroutine RK_coefficients  '
        STOP
    ENDIF
ENDIF

IF (NEW == 1) THEN

    IF (ORDER == 1) THEN
      CFL = 1.0d0/(2.0d0*PP+1.0d0)
      alpha(1,1) = 1.0d0
      beta(1,1)  = 1.0d0
    ENDIF      

    IF ((SSP1.EQ.3).AND.(SSP2.EQ.2)) THEN
        
            NRK = 3
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 0.099308256766310D0
            ATVD(2,2) = 0.900691743233690D0
            ATVD(3,1) = 0.318537560632609D0
            ATVD(3,3) = 0.681462439367391D0
            
            BTVD(1,1) = 0.521093099693580D0
            BTVD(2,2) = 0.469344252350057D0
            BTVD(3,1) = 0.005213666926231D0
            BTVD(3,3) = 0.355105374854702D0
            
    ELSEIF ((SSP1.EQ.4).AND.(SSP2.EQ.2)) THEN
         
            NRK = 4
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.000000000000000D0
            ATVD(2,1) = 0.059418613052739D0
            ATVD(3,1) = 0.162040236890209D0
            ATVD(4,1) = 0.213049467749813D0
            ATVD(2,2) = 0.940581386947261D0
            ATVD(3,2) = 0.000000000000001D0
            ATVD(4,2) = 0.171375745599785D0
            ATVD(3,3) = 0.837959763109790D0
            ATVD(4,4) = 0.615574786650403D0
   
            BTVD(1,1) = 0.396600082220546D0
            BTVD(3,1) = 0.042341962885445D0
            BTVD(4,1) = 0.004443461175963D0
            BTVD(2,2) = 0.373034655398399D0
            BTVD(3,2) = 0.000000000000001D0
            BTVD(4,2) = 0.067967634795482D0
            BTVD(3,3) = 0.332334910946852D0
            BTVD(4,4) = 0.244137010998445D0
            
    ELSEIF ((SSP1.EQ.5).AND.(SSP2.EQ.2)) THEN
         
            NRK = 5
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0

            ATVD(1,1) = 1.000000000000000D0
            ATVD(2,1) = 0.050443378596570D0
            ATVD(3,1) = 0.423050285840291D0
            ATVD(4,1) = 0.053551155948373D0
            ATVD(5,1) = 0.220414081710321D0
            ATVD(2,2) = 0.949556621403430D0
            ATVD(5,2) = 0.000090243705488D0
            ATVD(3,3) = 0.576949714159710D0
            ATVD(4,4) = 0.946448844051627D0
            ATVD(5,5) = 0.779495674584191D0

            BTVD(1,1) = 3.16264691752137D-1
            BTVD(2,2) = 3.00311232169356D-1
            BTVD(3,1) = 1.05127466571023D-1
            BTVD(3,3) = 1.82468823505204D-1
            BTVD(4,1) = 3.6518664289273D-2
            BTVD(4,4) = 2.99328351923154D-1
            BTVD(5,1) = 2.3796996941886D-2
            BTVD(5,2) = 2.8540897699D-5
            BTVD(5,5) = 2.46526959244493D-1
            
    ELSEIF ((SSP1.EQ.6).AND.(SSP2.EQ.2)) THEN

            NRK = 6
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 2.76486751644183D-1
            ATVD(2,2) = 7.23513248355817D-1
            ATVD(3,1) = 9.7984946364171D-2
            ATVD(3,3) = 9.02015053635831D-1
            ATVD(4,1) = 2.59872710659296D-1
            ATVD(4,2) = 7.856020000000001D-9
            ATVD(4,4) = 7.40127281484684D-1
            ATVD(5,1) = 8.7560711616745D-2
            ATVD(5,2) = 6.0713391259179D-2
            ATVD(5,3) = 5.141736418695D-3
            ATVD(5,5) = 8.4658416070538D-1
            ATVD(6,1) = 1.7424655647804D-1
            ATVD(6,2) = 7.0778D-11
            ATVD(6,3) = 6.04526D-9
            ATVD(6,4) = 5.401999866D-6
            ATVD(6,6) = 8.25748035406055D-1
      
            BTVD(1,1) = 2.62807064269699D-1
            BTVD(2,2) = 1.90144392760626D-1
            BTVD(3,1) = 6.9992824632912D-2
            BTVD(3,3) = 2.37055928173108D-1
            BTVD(4,1) = 6.838959780635399D-2
            BTVD(4,2) = 2.064618D-9
            BTVD(4,4) = 1.94510678032903D-1
            BTVD(5,1) = 3.7905872702872D-2
            BTVD(5,2) = 1.5955908118682D-2
            BTVD(5,3) = 1.351284653446D-3
            BTVD(5,5) = 2.22488297932208D-1
            BTVD(6,1) = 1.8749686097568D-2
            BTVD(6,2) = 1.8601D-11
            BTVD(6,3) = 1.588737D-9
            BTVD(6,4) = 1.419683726D-6
            BTVD(6,6) = 2.17012417011537D-1
            
    ELSEIF ((SSP1.EQ.7).AND.(SSP2.EQ.2)) THEN

            NRK = 7
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 8.5997676989D-5
            ATVD(2,2) = 9.99914002323011D-1
            ATVD(3,1) = 3.243699004044D-2
            ATVD(3,3) = 9.6756300995956D-1
            ATVD(4,1) = 1.3583973942284D-2
            ATVD(4,2) = 2.75922050535873D-1
            ATVD(4,4) = 7.10493975521843D-1
            ATVD(5,1) = 2.43475411423557D-1
            ATVD(5,2) = 8.174418516086D-3
            ATVD(5,3) = 1.182963464D-6
            ATVD(5,5) = 7.48348987096893D-1
            ATVD(6,1) = 2.1828984636811D-2
            ATVD(6,2) = 6.7598991294125D-2
            ATVD(6,4) = 1.98272850387789D-1
            ATVD(6,6) = 7.12299173681275D-1
            ATVD(7,1) = 8.455803656946199D-2
            ATVD(7,2) = 1.3632738067068D-2
            ATVD(7,3) = 7.340245076405D-3
            ATVD(7,4) = 1.35057999304014D-1
            ATVD(7,5) = 9.993316295892D-3
            ATVD(7,7) = 7.49417664687159D-1
      
            BTVD(1,1) = 2.26374432502229D-1
            BTVD(2,2) = 2.26354964826904D-1
            BTVD(3,1) = 6.213730480682D-3
            BTVD(3,3) = 2.19031527289744D-1
            BTVD(4,1) = 2.720075254276D-3
            BTVD(4,2) = 6.246169760491D-2
            BTVD(4,4) = 1.6083767050501D-1
            BTVD(5,1) = 2.3562742803953D-2
            BTVD(5,2) = 1.850479352615D-3
            BTVD(5,3) = 2.67792683D-7
            BTVD(5,5) = 1.69407077267677D-1
            BTVD(6,2) = 1.5302683291931D-2
            BTVD(6,4) = 4.4883903987135D-2
            BTVD(6,6) = 1.61246321213905D-1
            BTVD(7,1) = 4.737985600659D-3
            BTVD(7,2) = 3.086103343384D-3
            BTVD(7,3) = 1.661643813598D-3
            BTVD(7,4) = 3.0573677947333D-2
            BTVD(7,5) = 2.262231305298D-3
            BTVD(7,7) = 1.69648998550701D-1
            
    ELSEIF ((SSP1.EQ.8).AND.(SSP2.EQ.2)) THEN

            NRK = 8
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.891067780011D-3
            ATVD(2,2) = 9.95108932219989D-1
            ATVD(3,1) = 8.1093621337591D-2
            ATVD(3,3) = 9.18906378662409D-1
            ATVD(4,1) = 2.067896765374D-3
            ATVD(4,2) = 2.26357903879396D-1
            ATVD(4,4) = 7.7157419935523D-1
            ATVD(5,1) = 1.43804042464509D-1
            ATVD(5,2) = 2.24645670502005D-1
            ATVD(5,3) = 1.90754736D-7
            ATVD(5,5) = 6.3155009627875D-1
            ATVD(6,1) = 8.46689053411D-2
            ATVD(6,2) = 1.83257D-9
            ATVD(6,3) = 3.62842D-10
            ATVD(6,6) = 9.15331092463488D-1
            ATVD(7,1) = 2.4897187565872D-2
            ATVD(7,2) = 2.81616938388197D-1
            ATVD(7,3) = 2.2224486434211D-2
            ATVD(7,4) = 1.291291D-9
            ATVD(7,5) = 3.65980959D-7
            ATVD(7,7) = 6.7126102033947D-1
            ATVD(8,1) = 1.13118501196255D-1
            ATVD(8,2) = 1.423038356041D-3
            ATVD(8,4) = 1.127867D-9
            ATVD(8,5) = 2.8943683461359D-2
            ATVD(8,6) = 1.8450513075117D-2
            ATVD(8,8) = 8.38064262783361D-1
      
            BTVD(1,1) = 1.97619905264336D-1
            BTVD(2,2) = 1.96653332913008D-1
            BTVD(3,1) = 1.38261322304259D-1
            BTVD(3,3) = 1.81594191498059D-1
            BTVD(4,1) = 9.0686509528815D-2
            BTVD(4,2) = 4.473282752048D-2
            BTVD(4,4) = 1.52478420180986D-1
            BTVD(5,1) = 9.7907291032652D-2
            BTVD(5,2) = 4.4394456122649D-2
            BTVD(5,3) = 3.7696933D-8
            BTVD(5,5) = 1.24806870196289D-1
            BTVD(6,1) = 9.178414505525D-3
            BTVD(6,2) = 3.62152D-10
            BTVD(6,3) = 7.1705D-11
            BTVD(6,6) = 1.80887643778135D-1
            BTVD(7,1) = 1.1097719275733D-2
            BTVD(7,2) = 5.5653112685108D-2
            BTVD(7,3) = 4.392000903677D-3
            BTVD(7,4) = 2.55185D-10
            BTVD(7,5) = 7.2325123D-8
            BTVD(7,7) = 1.32654539247127D-1
            BTVD(8,1) = 1.2660624562512D-2
            BTVD(8,2) = 2.81220705108D-4
            BTVD(8,4) = 2.22889D-10
            BTVD(8,5) = 5.719847983635D-3
            BTVD(8,6) = 3.646188645983D-3
            BTVD(8,8) = 1.65618180216673D-1
            
    ELSEIF ((SSP1.EQ.4).AND.(SSP2.EQ.3)) THEN
            
            NRK = 4
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.97246785739694D-1
            ATVD(2,2) = 5.027532142603059D-1
            ATVD(3,1) = 3.86311870571299D-1
            ATVD(3,3) = 6.13688129428701D-1
            ATVD(4,1) = 3.23034811271828D-1
            ATVD(4,2) = 1.6196425547402D-2
            ATVD(4,4) = 6.6076876318077D-1
      
            BTVD(1,1) = 5.89041045727752D-1
            BTVD(2,2) = 2.96142279070879D-1
            BTVD(3,3) = 3.6148749750939D-1
            BTVD(4,1) = 1.12664801028946D-1
            BTVD(4,2) = 9.540359441493999D-3
            BTVD(4,4) = 3.89219923248234D-1
            
    ELSEIF ((SSP1.EQ.5).AND.(SSP2.EQ.3)) THEN

            NRK = 5
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 2.46364980850805D-1
            ATVD(2,2) = 7.53635019149195D-1
            ATVD(3,1) = 3.71117534702807D-1
            ATVD(3,3) = 6.28882465297193D-1
            ATVD(4,1) = 4.30270405975121D-1
            ATVD(4,2) = 5.687592D-9
            ATVD(4,4) = 5.69729588337289D-1
            ATVD(5,1) = 1.40088797078275D-1
            ATVD(5,2) = 2.09584D-9
            ATVD(5,5) = 8.59911200825884D-1
      
            BTVD(1,1) = 4.05845247402079D-1
            BTVD(2,2) = 3.05859190797476D-1
            BTVD(3,1) = 7.9031725047311D-2
            BTVD(3,3) = 2.55228959715369D-1
            BTVD(4,1) = 5.3493423549296D-2
            BTVD(4,2) = 2.308282D-9
            BTVD(4,4) = 2.31222045731032D-1
            BTVD(5,1) = 5.3948649659203D-2
            BTVD(5,2) = 8.50587D-10
            BTVD(5,5) = 3.48990874043D-1
            
    ELSEIF ((SSP1.EQ.6).AND.(SSP2.EQ.3)) THEN

            NRK = 6

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 3.12168767243578D-1
            ATVD(2,2) = 6.87831232756422D-1
            ATVD(3,1) = 7.716265806604999D-2
            ATVD(3,3) = 9.2283734193395D-1
            ATVD(4,1) = 4.12132458656463D-1
            ATVD(4,2) = 2.1717224605D-5
            ATVD(4,4) = 5.87845824118932D-1
            ATVD(5,1) = 3.73425049182364D-1
            ATVD(5,2) = 1.62494172731D-3
            ATVD(5,3) = 6.10182039447D-4
            ATVD(5,5) = 6.2433982705088D-1
            ATVD(6,1) = 1.9475670525709D-2
            ATVD(6,2) = 4.280756369499D-3
            ATVD(6,3) = 1.85954886452D-4
            ATVD(6,4) = 4.30495D-10
            ATVD(6,6) = 9.760576177878439D-1
      
            BTVD(1,1) = 3.13661711287779D-1
            BTVD(2,2) = 2.15746321543562D-1
            BTVD(3,1) = 6.5063517563946D-2
            BTVD(3,3) = 2.89458739911268D-1
            BTVD(4,1) = 9.591954428256901D-2
            BTVD(4,2) = 6.811861834D-6
            BTVD(4,4) = 1.84384727166519D-1
            BTVD(5,1) = 3.5438929973665D-2
            BTVD(5,2) = 5.09682002931D-4
            BTVD(5,3) = 1.9139074269D-4
            BTVD(5,5) = 1.95831498577895D-1
            BTVD(6,1) = 2.3383161968088D-2
            BTVD(6,2) = 1.342709368463D-3
            BTVD(6,3) = 5.8326927907D-5
            BTVD(6,4) = 1.3503D-10
            BTVD(6,6) = 3.06151902710808D-1
            
    ELSEIF ((SSP1.EQ.7).AND.(SSP2.EQ.3)) THEN

            NRK = 7

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 1.47691415319417D-1
            ATVD(2,2) = 8.52308584680583D-1
            ATVD(3,1) = 7.34539305979D-4
            ATVD(3,3) = 9.99265460694021D-1
            ATVD(4,1) = 3.92877396776094D-1
            ATVD(4,2) = 1.47887230450615D-1
            ATVD(4,4) = 4.59235372773291D-1
            ATVD(5,1) = 1.97049996198846D-1
            ATVD(5,2) = 9.368160402672999D-3
            ATVD(5,3) = 8.9935546009991D-2
            ATVD(5,5) = 7.036462973884891D-1
            ATVD(6,1) = 6.1993002028958D-2
            ATVD(6,3) = 6.1395792658554D-2
            ATVD(6,4) = 4.27057986508D-4
            ATVD(6,6) = 8.7618414732598D-1
            ATVD(7,1) = 1.28534993749396D-1
            ATVD(7,3) = 1.053984977328D-2
            ATVD(7,4) = 5.607329238D-6
            ATVD(7,5) = 3.07502271596D-4
            ATVD(7,7) = 8.60612046876491D-1
      
            BTVD(1,1) = 2.580794742818D-1
            BTVD(2,2) = 2.1996335146023D-1
            BTVD(3,1) = 2.319087748D-6
            BTVD(3,3) = 2.57889904763873D-1
            BTVD(4,1) = 4.3608780145257D-2
            BTVD(4,2) = 3.8166658687686D-2
            BTVD(4,4) = 1.18519223576937D-1
            BTVD(5,1) = 7.9389033244D-5
            BTVD(5,2) = 2.41772991171D-3
            BTVD(5,3) = 2.3210518433505D-2
            BTVD(5,5) = 1.81596666510356D-1
            BTVD(6,3) = 1.5844993892434D-2
            BTVD(6,4) = 1.10214900646D-4
            BTVD(6,6) = 2.26125144115936D-1
            BTVD(7,1) = 5.3986013420437D-2
            BTVD(7,3) = 2.720118888497D-3
            BTVD(7,4) = 1.447136582D-6
            BTVD(7,5) = 7.9360024594D-5
            BTVD(7,7) = 2.22106304618468D-1
            
    ELSEIF ((SSP1.EQ.8).AND.(SSP2.EQ.3)) THEN

            NRK = 8

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 1.452835644833D-3
            ATVD(2,2) = 9.98547164355167D-1
            ATVD(3,1) = 1.5824688647D-4
            ATVD(3,3) = 9.9984175311353D-1
            ATVD(4,1) = 1.35300645522171D-1
            ATVD(4,2) = 4.3910257822658D-2
            ATVD(4,4) = 8.20789096655171D-1
            ATVD(5,1) = 2.5418437348026D-1
            ATVD(5,2) = 2.40296919265655D-1
            ATVD(5,3) = 1.73107225811679D-1
            ATVD(5,5) = 3.32411481442406D-1
            ATVD(6,1) = 6.6883284219592D-2
            ATVD(6,2) = 6.766232631260299D-2
            ATVD(6,3) = 6.265014698872D-3
            ATVD(6,4) = 1.099839887838D-3
            ATVD(6,6) = 8.58089534881095D-1
            ATVD(7,1) = 1.7899811147945D-1
            ATVD(7,2) = 6.1236827D-8
            ATVD(7,3) = 2.0509168546D-4
            ATVD(7,4) = 3.5635846554D-5
            ATVD(7,5) = 1.045577373439D-3
            ATVD(7,7) = 8.197155223782711D-1
            ATVD(8,1) = 1.4239397335328D-2
            ATVD(8,2) = 4.979393069413D-3
            ATVD(8,3) = 4.9618753349D-5
            ATVD(8,4) = 4.08154071557D-4
            ATVD(8,5) = 1.95292343D-7
            ATVD(8,6) = 1.40625581752D-4
            ATVD(8,8) = 9.80182615896259D-1
      
            BTVD(1,1) = 2.19133970894904D-1
            BTVD(2,2) = 2.18815605250994D-1
            BTVD(3,1) = 3.7408889951057D-2
            BTVD(3,3) = 2.1909929362629D-1
            BTVD(4,1) = 3.63333123844D-4
            BTVD(4,2) = 9.622229159698001D-3
            BTVD(4,4) = 1.79862774017288D-1
            BTVD(5,1) = 1.4964409637913D-2
            BTVD(5,2) = 5.2657218112495D-2
            BTVD(5,3) = 3.7933673782714D-2
            BTVD(5,5) = 7.2842647899532D-2
            BTVD(6,1) = 7.217055487D-4
            BTVD(6,2) = 1.4827114244867D-2
            BTVD(6,3) = 1.372877548679D-3
            BTVD(6,4) = 2.4101228197D-4
            BTVD(6,6) = 1.88036567161855D-1
            BTVD(7,1) = 1.4581151700251D-2
            BTVD(7,2) = 1.3419069D-8
            BTVD(7,3) = 4.4942555432D-5
            BTVD(7,4) = 7.809024561999999D-6
            BTVD(7,5) = 2.29121521719D-4
            BTVD(7,7) = 1.79627517422941D-1
            BTVD(8,1) = 2.3565431064376D-2
            BTVD(8,2) = 1.091154175947D-3
            BTVD(8,3) = 1.0873154452D-5
            BTVD(8,4) = 8.9440422437D-5
            BTVD(8,5) = 4.2795187D-8
            BTVD(8,6) = 3.0815842139D-5
            BTVD(8,8) = 2.14791308823501D-1
            
    ELSEIF ((SSP1.EQ.5).AND.(SSP2.EQ.4)) THEN
         
            NRK = 5
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 5.78405168674975D-1
            ATVD(2,2) = 4.21594831325025D-1
            ATVD(3,1) = 6.77678199623741D-1
            ATVD(3,2) = 8.822045619988001D-3
            ATVD(3,3) = 3.13499754756272D-1
            ATVD(4,1) = 3.60703695071165D-1
            ATVD(4,2) = 1.5747555930398D-1
            ATVD(4,3) = 3.40817234398D-4
            ATVD(4,4) = 4.81479928390457D-1
            ATVD(5,1) = 4.03281245556668D-1
            ATVD(5,2) = 1.96797187365244D-1
            ATVD(5,3) = 2.9058018283826D-2
            ATVD(5,4) = 1.46343918738107D-1
            ATVD(5,5) = 2.24519630056155D-1
      
            BTVD(1,1) = 3.63990101160008D-1
            BTVD(2,1) = 2.0871083992052D-2
            BTVD(2,2) = 4.17052806266404D-1
            BTVD(3,1) = 1.745020311D-6
            BTVD(3,2) = 8.727001873487D-3
            BTVD(3,3) = 3.10122285119134D-1
            BTVD(4,1) = 5.820402293450001D-4
            BTVD(4,2) = 1.55779006397411D-1
            BTVD(4,3) = 3.37145461634D-4
            BTVD(4,4) = 4.76292735053434D-1
            BTVD(5,1) = 3.6738984755689D-2
            BTVD(5,2) = 1.94677005403644D-1
            BTVD(5,3) = 2.8744963574915D-2
            BTVD(5,4) = 1.44767291852061D-1
            BTVD(5,5) = 2.22100782124214D-1
            
    ELSEIF ((SSP1.EQ.6).AND.(SSP2.EQ.4)) THEN

            NRK = 6
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.67088614334124D-1
            ATVD(2,2) = 5.32911385665876D-1
            ATVD(3,1) = 4.08874458377822D-1
            ATVD(3,3) = 5.91125541622178D-1
            ATVD(4,1) = 4.48601090300132D-1
            ATVD(4,4) = 5.51398909699866D-1
            ATVD(5,1) = 1.12338903697D-4
            ATVD(5,2) = 8.930525958D-6
            ATVD(5,3) = 6.3196058D-8
            ATVD(5,5) = 9.99878667374287D-1
            ATVD(6,1) = 1.32061973305982D-1
            ATVD(6,2) = 1.14804153230642D-1
            ATVD(6,3) = 2.56477290949465D-1
            ATVD(6,4) = 1.16490556977394D-1
            ATVD(6,6) = 3.80166025536516D-1
      
            BTVD(1,1) = 4.42682663664825D-1
            BTVD(2,2) = 2.35910631703883D-1
            BTVD(3,3) = 2.61681029325618D-1
            BTVD(4,1) = 3.46D-13
            BTVD(4,4) = 2.44094738087817D-1
            BTVD(5,1) = 1.0465373390104D-2
            BTVD(5,2) = 3.953389019D-6
            BTVD(5,3) = 2.7975799D-8
            BTVD(5,5) = 4.42628951814885D-1
            BTVD(6,1) = 2.62732028444D-3
            BTVD(6,2) = 5.0821808351925D-2
            BTVD(6,3) = 1.13538050327048D-1
            BTVD(6,4) = 5.1568350054552D-2
            BTVD(6,6) = 1.68292908819375D-1
            
    ELSEIF ((SSP1.EQ.7).AND.(SSP2.EQ.4)) THEN

            NRK = 7
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 5.25693525536155D-1
            ATVD(2,2) = 4.74306474463845D-1
            ATVD(3,1) = 8.3349916526D-5
            ATVD(3,3) = 9.99916650083474D-1
            ATVD(4,1) = 5.24867751177608D-1
            ATVD(4,2) = 3.69729028192D-4
            ATVD(4,4) = 4.747625197942D-1
            ATVD(5,1) = 1.95877892721119D-1
            ATVD(5,2) = 4.30179D-9
            ATVD(5,3) = 3.855559176D-6
            ATVD(5,5) = 8.04118247417915D-1
            ATVD(6,1) = 6.8982360418487D-2
            ATVD(6,2) = 1.3222295413188D-2
            ATVD(6,3) = 7.59402789052D-4
            ATVD(6,4) = 2.82563389428D-4
            ATVD(6,6) = 9.16753377989845D-1
            ATVD(7,1) = 8.509651769937999D-2
            ATVD(7,2) = 7.019462801832D-3
            ATVD(7,3) = 1.51045322680285D-1
            ATVD(7,4) = 3.15632962717013D-1
            ATVD(7,5) = 2.3712178831159D-2
            ATVD(7,7) = 4.17493555270331D-1
      
            BTVD(1,1) = 3.43396366222189D-1
            BTVD(2,2) = 1.62875119806542D-1
            BTVD(3,1) = 4.313430735199D-3
            BTVD(3,3) = 3.43367744163729D-1
            BTVD(4,1) = 2.490553774367D-3
            BTVD(4,2) = 1.26963604768D-4
            BTVD(4,4) = 1.63031724115818D-1
            BTVD(5,1) = 4.766901826487D-3
            BTVD(5,2) = 1.477219D-9
            BTVD(5,3) = 1.323985011D-6
            BTVD(5,5) = 2.76131284176267D-1
            BTVD(6,1) = 4.2233750546714D-2
            BTVD(6,2) = 4.540488198005D-3
            BTVD(6,3) = 2.60776158259D-4
            BTVD(6,4) = 9.7031241157D-5
            BTVD(6,6) = 3.1480977872363D-1
            BTVD(7,1) = 3.7417046429D-5
            BTVD(7,2) = 2.410458018981D-3
            BTVD(7,3) = 5.1868414943268D-2
            BTVD(7,4) = 1.08387212456966D-1
            BTVD(7,5) = 8.142676045831D-3
            BTVD(7,7) = 1.43365769801014D-1
      
    ELSEIF ((SSP1.EQ.8).AND.(SSP2.EQ.4)) THEN

            NRK = 8

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.61181537109233D-1
            ATVD(2,2) = 5.38818462890767D-1
            ATVD(3,1) = 2.0117158505919D-2
            ATVD(3,3) = 9.79882841494082D-1
            ATVD(4,1) = 1.69428274041103D-1
            ATVD(4,2) = 5.822D-12
            ATVD(4,4) = 8.30571725953071D-1
            ATVD(5,1) = 4.31077025374928D-1
            ATVD(5,2) = 9.6737381267313D-2
            ATVD(5,3) = 8.8663023468039D-2
            ATVD(5,5) = 3.83522569889722D-1
            ATVD(6,1) = 8.9850295244391D-2
            ATVD(6,2) = 5.5129D-11
            ATVD(6,3) = 1.2320078306704D-2
            ATVD(6,4) = 2.95511D-10
            ATVD(6,6) = 8.97829626098266D-1
            ATVD(7,1) = 6.244345604206D-3
            ATVD(7,2) = 4.339895694D-4
            ATVD(7,3) = 3.515550716564D-3
            ATVD(7,4) = 9.393008934879D-3
            ATVD(7,5) = 1.339131D-9
            ATVD(7,7) = 9.8041310383582D-1
            ATVD(8,1) = 7.1188337795849D-2
            ATVD(8,2) = 3.4167417741535D-2
            ATVD(8,3) = 1.632796D-9
            ATVD(8,4) = 2.63316641735748D-1
            ATVD(8,5) = 1.51537858712815D-1
            ATVD(8,6) = 1.14976594D-7
            ATVD(8,8) = 4.79789627404664D-1
      
            BTVD(1,1) = 2.78271529217822D-1
            BTVD(2,2) = 1.4993783763941D-1
            BTVD(3,1) = 1.6659D-11
            BTVD(3,3) = 2.72673496756863D-1
            BTVD(4,1) = 2.4301359D-8
            BTVD(4,2) = 1.62D-12
            BTVD(4,4) = 2.31124464306047D-1
            BTVD(5,1) = 1.267923D-9
            BTVD(5,2) = 2.6919259017783D-2
            BTVD(5,3) = 2.4672395125527D-2
            BTVD(5,5) = 1.06723412012762D-1
            BTVD(6,1) = 3.2744D-11
            BTVD(6,2) = 1.5341D-11
            BTVD(6,3) = 3.42832703049D-3
            BTVD(6,4) = 8.2232D-11
            BTVD(6,6) = 2.4984042303143D-1
            BTVD(7,1) = 3.4914608275757D-2
            BTVD(7,2) = 1.20766941141D-4
            BTVD(7,3) = 9.782776739410001D-4
            BTVD(7,4) = 2.613806960265D-3
            BTVD(7,5) = 3.72642D-10
            BTVD(7,7) = 2.72821053669585D-1
            BTVD(8,1) = 4.265D-11
            BTVD(8,2) = 9.507819584361D-3
            BTVD(8,3) = 4.54361D-10
            BTVD(8,4) = 7.3273524564308D-2
            BTVD(8,5) = 4.2168671678409D-2
            BTVD(8,6) = 3.1994713D-8
            BTVD(8,8) = 1.33511793320745D-1
            
    ENDIF            
            
ELSE

!---------------------------------------------------------------------------------
  IF (ORDER == 1) THEN
!---------------------------------------------------------------------------------
      CFL = 1.0d0/(2.0d0*PP+1.0d0)
      alpha(1,1) = 1.0d0
      beta(1,1)  = 1.0d0
!----------------------------------------------------------------------------------
  ELSEIF (ORDER == 2) THEN
!---------------------------------------------------------------------------------
      
!....2 STAGES

      IF (S == 2) THEN

          CFL = 1.0d0/(2.0d0*PP+1.0d0)
          alpha(1,1) = 1.0d0
          alpha(1,2) = 0.0d0
          alpha(2,1) = 1.0d0/2.0d0
          alpha(2,2) = alpha(2,1)
          beta(1,1)  = 1.0d0
          beta(1,2)  = 0.0d0
          beta(2,1)  = 0.0d0
          beta(2,2) = 1.0d0/2.0d0
       
!....3 STAGES

       ELSEIF (S == 3) THEN

          CFL = 0.590449078021362d0
          alpha(1,1) = 1.000000000000000d0
          alpha(1,2) = 0.000000000000000d0
          alpha(1,3) = 0.000000000000000d0
          alpha(2,1) = 0.087353119859156d0
          alpha(2,2) = 0.912646880140844d0
          alpha(2,3) = 0.000000000000000d0
          alpha(3,1) = 0.344956917166841d0
          alpha(3,2) = 0.000000000000000d0
          alpha(3,3) = 0.655043082833159d0
          beta(1,1)  = 0.528005024856522d0
          beta(1,2)  = 0.000000000000000d0
          beta(1,3)  = 0.000000000000000d0
          beta(2,1)  = 0.000000000000000d0
          beta(2,2)  = 0.481882138633993d0
          beta(2,3)  = 0.000000000000000d0
          beta(3,1)  = 0.022826837460491d0
          beta(3,2)  = 0.000000000000000d0
          beta(3,3)  = 0.345866039233415d0

!....8 STAGES

       ELSEIF(S == 8) THEN

          CFL = 1.711397030610134d0
          alpha(1,1) = 1.000000000000000d0
          alpha(1,2) = 0.000000000000000d0
          alpha(1,3) = 0.000000000000000d0
          alpha(1,4) = 0.000000000000000d0
          alpha(1,5) = 0.000000000000000d0
          alpha(1,6) = 0.000000000000000d0
          alpha(1,7) = 0.000000000000000d0
          alpha(1,8) = 0.000000000000000d0
          alpha(2,1) = 0.000022706789062d0
          alpha(2,2) = 0.999977293210938d0
          alpha(2,3) = 0.000000000000000d0
	  alpha(2,4) = 0.000000000000000d0
	  alpha(2,5) = 0.000000000000000d0
      alpha(2,6) = 0.000000000000000d0
	  alpha(2,7) = 0.000000000000000d0
	  alpha(2,8) = 0.000000000000000d0
	  alpha(3,1) = 0.040163352005589d0
	  alpha(3,2) = -0.000000000000000d0
	  alpha(3,3) = 0.959836647994411d0
	  alpha(3,4) = 0.000000000000000d0
 	  alpha(3,5) = 0.000000000000000d0
	  alpha(3,6) = 0.000000000000000d0
	  alpha(3,7) = 0.000000000000000d0
	  alpha(3,8) = 0.000000000000000d0
	  alpha(4,1) = 0.068277481460787d0
      alpha(4,2) = 0.000794112060878d0
	  alpha(4,3) = -0.000000000000000d0
      alpha(4,4) = 0.930928406478335d0
 	  alpha(4,5) = 0.000000000000000d0
	  alpha(4,6) = 0.000000000000000d0
	  alpha(4,7) = 0.000000000000000d0
	  alpha(4,8) = 0.000000000000000d0
	  alpha(5,1) = 0.203201797047416d0
      alpha(5,2) = 0.049250734246472d0
	  alpha(5,3) = 0.028073821192189d0
	  alpha(5,4) = 0.000000000000000d0
	  alpha(5,5) = 0.719473647513923d0
	  alpha(5,6) = 0.000000000000000d0
	  alpha(5,7) = 0.000000000000000d0
	  alpha(5,8) = 0.000000000000000d0
	  alpha(6,1) = 0.135232895616102d0
	  alpha(6,2) = 0.000233842652217d0
	  alpha(6,3) = 0.057212625463321d0
	  alpha(6,4) = 0.005935142172269d0
	  alpha(6,5) = 0.000000000000000d0
	  alpha(6,6) = 0.801385494096090d0
	  alpha(6,7) = 0.000000000000000d0
	  alpha(6,8) = 0.000000000000000d0
	  alpha(7,1) = 0.330663126197853d0
	  alpha(7,2) = -0.000000000000000d0
	  alpha(7,3) = 0.000000000000000d0
	  alpha(7,4) = 0.036269305920717d0
	  alpha(7,5) = 0.087351183650948d0
	  alpha(7,6) = 0.000000000000000d0
	  alpha(7,7) = 0.545716384230481d0
	  alpha(7,8) = 0.000000000000000d0 	  
	  alpha(8,1) = 0.072377654137843d0
	  alpha(8,2) = 0.001047652272129d0
	  alpha(8,3) = -0.000000000000000d0
	  alpha(8,4) = 0.005351824870836d0
	  alpha(8,5) = 0.206003319412558d0
	  alpha(8,6) = 0.025024389094107d0
	  alpha(8,7) = 0.000000000000000d0
	  alpha(8,8) = 0.690195160212526d0
	  beta(1,1)  = 0.203816348874755d0
	  beta(1,2)  = 0.000000000000000d0
	  beta(1,3)  = 0.000000000000000d0
	  beta(1,4)  = 0.000000000000000d0
	  beta(1,5)  = 0.000000000000000d0
	  beta(1,6)  = 0.000000000000000d0
	  beta(1,7)  = 0.000000000000000d0
	  beta(1,8)  = 0.000000000000000d0
	  beta(2,1)  = 0.000000000000000d0
	  beta(2,2)  = 0.203811720859914d0
	  beta(2,3)  = 0.000000000000000d0
	  beta(2,4)  = 0.000000000000000d0
	  beta(2,5)  = 0.000000000000000d0
	  beta(2,6)  = 0.000000000000000d0
	  beta(2,7)  = 0.000000000000000d0
	  beta(2,8)  = 0.000000000000000d0
	  beta(3,1)  = 0.203690674361979d0
	  beta(3,2)  = -0.000000000000000d0
	  beta(3,3)  = 0.195630401110404d0
	  beta(3,4)  = 0.000000000000000d0
	  beta(3,5)  = 0.000000000000000d0
	  beta(3,6)  = 0.000000000000000d0
	  beta(3,7)  = 0.000000000000000d0
	  beta(3,8)  = 0.000000000000000d0
	  beta(4,1)  = 0.036121835699129d0
	  beta(4,2)  = 0.000161853020846d0
	  beta(4,3)  = -0.000000000000000d0
	  beta(4,4)  = 0.189738428872208d0
	  beta(4,5)  = 0.000000000000000d0
	  beta(4,6)  = 0.000000000000000d0
	  beta(4,7)  = 0.000000000000000d0
	  beta(4,8)  = 0.000000000000000d0
	  beta(5,1)  = 0.084801239646524d0
	  beta(5,2)  = 0.010038104833517d0
	  beta(5,3)  = 0.005721903734355d0
	  beta(5,4)  = 0.000000000000000d0
	  beta(5,5)  = 0.146640491947890d0
	  beta(5,6)  = 0.000000000000000d0
	  beta(5,7)  = 0.000000000000000d0
	  beta(5,8)  = 0.000000000000000d0
	  beta(6,1)  = -0.000000000000000d0
	  beta(6,2)  =  0.000047660955586d0
	  beta(6,3)  = 0.011660868431473d0
	  beta(6,4)  = 0.001209679007605d0
	  beta(6,5)  = 0.000000000000000d0
	  beta(6,6)  = 0.163335465447857d0
	  beta(6,7)  = 0.000000000000000d0
	  beta(6,8)  = 0.000000000000000d0
	  beta(7,1)  = 0.058391312304979d0
	  beta(7,2)  = -0.000000000000000d0
	  beta(7,3)  = 0.000000000000000d0
	  beta(7,4)  = 0.007392277508982d0
	  beta(7,5)  = 0.017803599321625d0
	  beta(7,6)  = 0.000000000000000d0
	  beta(7,7)  = 0.111225920954990d0
	  beta(7,8)  = 0.000000000000000d0
	  beta(8,1)  = 0.004977697629959d0
	  beta(8,2)  = 0.000213528660996d0
	  beta(8,3)  = -0.000000000000000d0
	  beta(8,4)  = 0.001090789404991d0
	  beta(8,5)  = 0.041986844418748d0
	  beta(8,6)  = 0.005100379617982d0
	  beta(8,7)  = 0.000000000000000d0
	  beta(8,8)  = 0.140673057565544d0
     
      ENDIF

!----------------------------------------------------------------------------------
  ELSEIF (ORDER == 3) THEN
!---------------------------------------------------------------------------------- 
 
!....3 STAGES

      IF (S == 3) THEN
	  CFL = 1.0d0/(2.0d0*PP+1.0d0)	
      alpha(1,1) = 1.000000000000000d0
	  alpha(1,2) = 0.000000000000000d0
	  alpha(1,3) = 0.000000000000000d0
	  alpha(2,1) = 0.750000000000000d0
	  alpha(2,2) = 0.250000000000000d0
	  alpha(2,3) = 0.000000000000000d0
	  alpha(3,1) = 0.333333333333333d0
	  alpha(3,2) = 0.000000000000000d0
	  alpha(3,3) = 0.666666666666667d0
	  beta(1,1)  = 1.000000000000000d0
	  beta(1,2)  = 0.000000000000000d0
	  beta(1,3)  = 0.000000000000000d0
	  beta(2,1)  = 0.000000000000000d0
	  beta(2,2)  = 0.250000000000000d0
	  beta(2,3)  = 0.000000000000000d0
	  beta(3,1)  = 0.000000000000000d0
	  beta(3,2)  = 0.000000000000000d0
	  beta(3,3)  = 0.666666666666667d0

      ELSEIF (S == 8) THEN


      ENDIF

!------------------------------------------------------------------------------
  ELSEIF (ORDER == 4) THEN
!-------------------------------------------------------------------------------

!....5 Stages (there is no 4th order 4 stage method)

      IF (S == 5) THEN
          
	  CFL = 0.220062927988164d0
      alpha(1,1) = 1.000000000000000d0
	  alpha(1,2) = 0.000000000000000d0
	  alpha(1,3) = 0.000000000000000d0
	  alpha(1,4) = 0.000000000000000d0
	  alpha(1,5) = 0.000000000000000d0
	  alpha(2,1) = 0.261216512493821d0
	  alpha(2,2) = 0.738783487506179d0
 	  alpha(2,3) = 0.000000000000000d0  
	  alpha(2,4) = 0.000000000000000d0
	  alpha(2,5) = 0.000000000000000d0
	  alpha(3,1) = 0.623613752757655d0
	  alpha(3,2) = 0.000000000000000d0
	  alpha(3,3) = 0.376386247242345d0
	  alpha(3,4) = 0.000000000000000d0
	  alpha(3,5) = 0.000000000000000d0
	  alpha(4,1) = 0.444745181201454d0
	  alpha(4,2) = 0.120932584902288d0
	  alpha(4,3) = 0.000000000000000d0
	  alpha(4,4) = 0.434322233896258d0
	  alpha(4,5) = 0.000000000000000d0
	  alpha(5,1) = 0.213357715199957d0
	  alpha(5,2) = 0.209928473023448d0
	  alpha(5,3) = 0.063353148180384d0 
	  alpha(5,4) = -0.000000000000000d0
	  alpha(5,5) = 0.513360663596212d0
	  beta(1,1)  = 0.605491839566400d0
	  beta(1,2)  = 0.000000000000000d0
	  beta(1,3)  = 0.000000000000000d0
	  beta(1,4)  = 0.000000000000000d0
	  beta(1,5)  = 0.000000000000000d0
	  beta(2,1)  = 0.000000000000000d0
	  beta(2,2)  = 0.447327372891397d0
	  beta(2,3)  = 0.000000000000000d0
	  beta(2,4)  = 0.000000000000000d0
	  beta(2,5)  = 0.000000000000000d0
	  beta(3,1)  = 0.000000844149769d0
	  beta(3,2)  = 0.000000000000000d0
	  beta(3,3)  = 0.227898801230261d0
	  beta(3,4)  = 0.000000000000000d0
	  beta(3,5)  = 0.000000000000000d0
	  beta(4,1)  = 0.002856233144485d0
	  beta(4,2)  = 0.073223693296006d0
	  beta(4,3)  = 0.000000000000000d0
	  beta(4,4)  = 0.262978568366434d0
	  beta(4,5)  = 0.000000000000000d0
	  beta(5,1)  = 0.002362549760441d0
	  beta(5,2)  = 0.127109977308333d0 
	  beta(5,3)  = 0.038359814234063d0
	  beta(5,4)  = -0.000000000000000d0
	  beta(5,5)  = 0.310835692561898d0 

      ENDIF
 
   ENDIF
   
ENDIF   

IF (order > 1) THEN
    alpha = ATVD
    beta  = BTVD
ENDIF

IF (debug == 4) THEN
    DO j = 1, S
        write(*,*)beta(:,j)
	write(*,*)'-------'
    ENDDO
ENDIF

!....ct coefficients for exact solution in time

DO i = 0, S-1
    DO j = 1, S
	csum = 0.0d0
	DO ii = i+1, j-1
	    csum = csum + crk(ii,i+1)*alpha(j,ii+1)
	ENDDO
	crk(j,i+1) = beta(j,i+1) + csum
    ENDDO
ENDDO

DO i = 1, S-1
    ct(i+1) = 0.0d0
    DO j = 0, i-1
        ct(i+1) = ct(i+1) + crk(i,j+1)
    ENDDO
ENDDO

IF (debug == 4) THEN
	write(*,*)ct
ENDIF
 
RETURN
END SUBROUTINE RK_COEFFICIENTS
!------------------------------------------------------------------------------------
!
!                 Subroutine WRK Coefficients
!                 Written by Colton J. Conroy
!                 @ the Isaac Newton Institute
!                           12.12.12
!
!-----------------------------------------------------------------------------------
SUBROUTINE WRK_COEFFICIENTS(ORDER,S)

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p) :: ORDER, S, PP, SSP1, SSP2, nrk
REAL(real_p) :: csum
REAL(real_p), DIMENSION(S,S) :: crk
REAL(real_p), ALLOCATABLE :: ATVD(:,:), BTVD(:,:), DTVD(:)

!....Constants

PP = ORDER - 1
SSP1 = S
SSP2 = ORDER

!....Allocate matricies

ALLOCATE(alphaW(S,S))
ALLOCATE(betaW(S,S))
ALLOCATE(ctW(S))

ALLOCATE(ATVD(S,S))
ALLOCATE(BTVD(S,S))
ALLOCATE(DTVD(S))

!....Initialize matricies

alphaW = 0.0d0
betaW  = 0.0d0
crk    = 0.0d0
ctW    = 0.0d0

!....Check input

IF (S < 1 .OR. S > 8) THEN
    PRINT*,'  ****** ERROR!!! Number of stages must be > 1 and < 8 *******'
    PRINT*,'  Execution terminated in subroutine RK_coefficients  '
    STOP
ENDIF

IF (ORDER < 1 .OR. ORDER > 4) THEN
    PRINT*,'  ****** ERROR!!! Order cannot be < 1 or > 4 *******'
    PRINT*,'  Execution terminated in subroutine RK_coefficients  '
    STOP
ENDIF

IF (S < ORDER) THEN
    PRINT*,'  ****** ERROR!!! Stages must be >= order (> for Order 4) *******'
    PRINT*,'  Execution terminated in subroutine RK_coefficients  '
    STOP
ENDIF

IF (ORDER == 4) THEN
    IF (S <= ORDER) THEN
        PRINT*,'  ****** ERROR!!! Stages must be > order for 4th Order *******'
        PRINT*,'  Execution terminated in subroutine RK_coefficients  '
        STOP
    ENDIF
ENDIF


IF (ORDER == 1) THEN
      CFL = 1.0d0/(2.0d0*PP+1.0d0)
      alphaW(1,1) = 1.0d0
      betaW(1,1)  = 1.0d0
ENDIF      

IF ((SSP1.EQ.3).AND.(SSP2.EQ.2)) THEN
        
            NRK = 3
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 0.099308256766310D0
            ATVD(2,2) = 0.900691743233690D0
            ATVD(3,1) = 0.318537560632609D0
            ATVD(3,3) = 0.681462439367391D0
            
            BTVD(1,1) = 0.521093099693580D0
            BTVD(2,2) = 0.469344252350057D0
            BTVD(3,1) = 0.005213666926231D0
            BTVD(3,3) = 0.355105374854702D0
            
ELSEIF ((SSP1.EQ.4).AND.(SSP2.EQ.2)) THEN
         
            NRK = 4
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.000000000000000D0
            ATVD(2,1) = 0.059418613052739D0
            ATVD(3,1) = 0.162040236890209D0
            ATVD(4,1) = 0.213049467749813D0
            ATVD(2,2) = 0.940581386947261D0
            ATVD(3,2) = 0.000000000000001D0
            ATVD(4,2) = 0.171375745599785D0
            ATVD(3,3) = 0.837959763109790D0
            ATVD(4,4) = 0.615574786650403D0
   
            BTVD(1,1) = 0.396600082220546D0
            BTVD(3,1) = 0.042341962885445D0
            BTVD(4,1) = 0.004443461175963D0
            BTVD(2,2) = 0.373034655398399D0
            BTVD(3,2) = 0.000000000000001D0
            BTVD(4,2) = 0.067967634795482D0
            BTVD(3,3) = 0.332334910946852D0
            BTVD(4,4) = 0.244137010998445D0
            
ELSEIF ((SSP1.EQ.5).AND.(SSP2.EQ.2)) THEN
         
            NRK = 5
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0

            ATVD(1,1) = 1.000000000000000D0
            ATVD(2,1) = 0.050443378596570D0
            ATVD(3,1) = 0.423050285840291D0
            ATVD(4,1) = 0.053551155948373D0
            ATVD(5,1) = 0.220414081710321D0
            ATVD(2,2) = 0.949556621403430D0
            ATVD(5,2) = 0.000090243705488D0
            ATVD(3,3) = 0.576949714159710D0
            ATVD(4,4) = 0.946448844051627D0
            ATVD(5,5) = 0.779495674584191D0

            BTVD(1,1) = 3.16264691752137D-1
            BTVD(2,2) = 3.00311232169356D-1
            BTVD(3,1) = 1.05127466571023D-1
            BTVD(3,3) = 1.82468823505204D-1
            BTVD(4,1) = 3.6518664289273D-2
            BTVD(4,4) = 2.99328351923154D-1
            BTVD(5,1) = 2.3796996941886D-2
            BTVD(5,2) = 2.8540897699D-5
            BTVD(5,5) = 2.46526959244493D-1
            
ELSEIF ((SSP1.EQ.6).AND.(SSP2.EQ.2)) THEN

            NRK = 6
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 2.76486751644183D-1
            ATVD(2,2) = 7.23513248355817D-1
            ATVD(3,1) = 9.7984946364171D-2
            ATVD(3,3) = 9.02015053635831D-1
            ATVD(4,1) = 2.59872710659296D-1
            ATVD(4,2) = 7.856020000000001D-9
            ATVD(4,4) = 7.40127281484684D-1
            ATVD(5,1) = 8.7560711616745D-2
            ATVD(5,2) = 6.0713391259179D-2
            ATVD(5,3) = 5.141736418695D-3
            ATVD(5,5) = 8.4658416070538D-1
            ATVD(6,1) = 1.7424655647804D-1
            ATVD(6,2) = 7.0778D-11
            ATVD(6,3) = 6.04526D-9
            ATVD(6,4) = 5.401999866D-6
            ATVD(6,6) = 8.25748035406055D-1
      
            BTVD(1,1) = 2.62807064269699D-1
            BTVD(2,2) = 1.90144392760626D-1
            BTVD(3,1) = 6.9992824632912D-2
            BTVD(3,3) = 2.37055928173108D-1
            BTVD(4,1) = 6.838959780635399D-2
            BTVD(4,2) = 2.064618D-9
            BTVD(4,4) = 1.94510678032903D-1
            BTVD(5,1) = 3.7905872702872D-2
            BTVD(5,2) = 1.5955908118682D-2
            BTVD(5,3) = 1.351284653446D-3
            BTVD(5,5) = 2.22488297932208D-1
            BTVD(6,1) = 1.8749686097568D-2
            BTVD(6,2) = 1.8601D-11
            BTVD(6,3) = 1.588737D-9
            BTVD(6,4) = 1.419683726D-6
            BTVD(6,6) = 2.17012417011537D-1
            
ELSEIF ((SSP1.EQ.7).AND.(SSP2.EQ.2)) THEN

            NRK = 7
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 8.5997676989D-5
            ATVD(2,2) = 9.99914002323011D-1
            ATVD(3,1) = 3.243699004044D-2
            ATVD(3,3) = 9.6756300995956D-1
            ATVD(4,1) = 1.3583973942284D-2
            ATVD(4,2) = 2.75922050535873D-1
            ATVD(4,4) = 7.10493975521843D-1
            ATVD(5,1) = 2.43475411423557D-1
            ATVD(5,2) = 8.174418516086D-3
            ATVD(5,3) = 1.182963464D-6
            ATVD(5,5) = 7.48348987096893D-1
            ATVD(6,1) = 2.1828984636811D-2
            ATVD(6,2) = 6.7598991294125D-2
            ATVD(6,4) = 1.98272850387789D-1
            ATVD(6,6) = 7.12299173681275D-1
            ATVD(7,1) = 8.455803656946199D-2
            ATVD(7,2) = 1.3632738067068D-2
            ATVD(7,3) = 7.340245076405D-3
            ATVD(7,4) = 1.35057999304014D-1
            ATVD(7,5) = 9.993316295892D-3
            ATVD(7,7) = 7.49417664687159D-1
      
            BTVD(1,1) = 2.26374432502229D-1
            BTVD(2,2) = 2.26354964826904D-1
            BTVD(3,1) = 6.213730480682D-3
            BTVD(3,3) = 2.19031527289744D-1
            BTVD(4,1) = 2.720075254276D-3
            BTVD(4,2) = 6.246169760491D-2
            BTVD(4,4) = 1.6083767050501D-1
            BTVD(5,1) = 2.3562742803953D-2
            BTVD(5,2) = 1.850479352615D-3
            BTVD(5,3) = 2.67792683D-7
            BTVD(5,5) = 1.69407077267677D-1
            BTVD(6,2) = 1.5302683291931D-2
            BTVD(6,4) = 4.4883903987135D-2
            BTVD(6,6) = 1.61246321213905D-1
            BTVD(7,1) = 4.737985600659D-3
            BTVD(7,2) = 3.086103343384D-3
            BTVD(7,3) = 1.661643813598D-3
            BTVD(7,4) = 3.0573677947333D-2
            BTVD(7,5) = 2.262231305298D-3
            BTVD(7,7) = 1.69648998550701D-1
            
ELSEIF ((SSP1.EQ.8).AND.(SSP2.EQ.2)) THEN

            NRK = 8
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.891067780011D-3
            ATVD(2,2) = 9.95108932219989D-1
            ATVD(3,1) = 8.1093621337591D-2
            ATVD(3,3) = 9.18906378662409D-1
            ATVD(4,1) = 2.067896765374D-3
            ATVD(4,2) = 2.26357903879396D-1
            ATVD(4,4) = 7.7157419935523D-1
            ATVD(5,1) = 1.43804042464509D-1
            ATVD(5,2) = 2.24645670502005D-1
            ATVD(5,3) = 1.90754736D-7
            ATVD(5,5) = 6.3155009627875D-1
            ATVD(6,1) = 8.46689053411D-2
            ATVD(6,2) = 1.83257D-9
            ATVD(6,3) = 3.62842D-10
            ATVD(6,6) = 9.15331092463488D-1
            ATVD(7,1) = 2.4897187565872D-2
            ATVD(7,2) = 2.81616938388197D-1
            ATVD(7,3) = 2.2224486434211D-2
            ATVD(7,4) = 1.291291D-9
            ATVD(7,5) = 3.65980959D-7
            ATVD(7,7) = 6.7126102033947D-1
            ATVD(8,1) = 1.13118501196255D-1
            ATVD(8,2) = 1.423038356041D-3
            ATVD(8,4) = 1.127867D-9
            ATVD(8,5) = 2.8943683461359D-2
            ATVD(8,6) = 1.8450513075117D-2
            ATVD(8,8) = 8.38064262783361D-1
      
            BTVD(1,1) = 1.97619905264336D-1
            BTVD(2,2) = 1.96653332913008D-1
            BTVD(3,1) = 1.38261322304259D-1
            BTVD(3,3) = 1.81594191498059D-1
            BTVD(4,1) = 9.0686509528815D-2
            BTVD(4,2) = 4.473282752048D-2
            BTVD(4,4) = 1.52478420180986D-1
            BTVD(5,1) = 9.7907291032652D-2
            BTVD(5,2) = 4.4394456122649D-2
            BTVD(5,3) = 3.7696933D-8
            BTVD(5,5) = 1.24806870196289D-1
            BTVD(6,1) = 9.178414505525D-3
            BTVD(6,2) = 3.62152D-10
            BTVD(6,3) = 7.1705D-11
            BTVD(6,6) = 1.80887643778135D-1
            BTVD(7,1) = 1.1097719275733D-2
            BTVD(7,2) = 5.5653112685108D-2
            BTVD(7,3) = 4.392000903677D-3
            BTVD(7,4) = 2.55185D-10
            BTVD(7,5) = 7.2325123D-8
            BTVD(7,7) = 1.32654539247127D-1
            BTVD(8,1) = 1.2660624562512D-2
            BTVD(8,2) = 2.81220705108D-4
            BTVD(8,4) = 2.22889D-10
            BTVD(8,5) = 5.719847983635D-3
            BTVD(8,6) = 3.646188645983D-3
            BTVD(8,8) = 1.65618180216673D-1
            
ELSEIF ((SSP1.EQ.4).AND.(SSP2.EQ.3)) THEN
            
            NRK = 4
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.97246785739694D-1
            ATVD(2,2) = 5.027532142603059D-1
            ATVD(3,1) = 3.86311870571299D-1
            ATVD(3,3) = 6.13688129428701D-1
            ATVD(4,1) = 3.23034811271828D-1
            ATVD(4,2) = 1.6196425547402D-2
            ATVD(4,4) = 6.6076876318077D-1
      
            BTVD(1,1) = 5.89041045727752D-1
            BTVD(2,2) = 2.96142279070879D-1
            BTVD(3,3) = 3.6148749750939D-1
            BTVD(4,1) = 1.12664801028946D-1
            BTVD(4,2) = 9.540359441493999D-3
            BTVD(4,4) = 3.89219923248234D-1
            
ELSEIF ((SSP1.EQ.5).AND.(SSP2.EQ.3)) THEN

            NRK = 5
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 2.46364980850805D-1
            ATVD(2,2) = 7.53635019149195D-1
            ATVD(3,1) = 3.71117534702807D-1
            ATVD(3,3) = 6.28882465297193D-1
            ATVD(4,1) = 4.30270405975121D-1
            ATVD(4,2) = 5.687592D-9
            ATVD(4,4) = 5.69729588337289D-1
            ATVD(5,1) = 1.40088797078275D-1
            ATVD(5,2) = 2.09584D-9
            ATVD(5,5) = 8.59911200825884D-1
      
            BTVD(1,1) = 4.05845247402079D-1
            BTVD(2,2) = 3.05859190797476D-1
            BTVD(3,1) = 7.9031725047311D-2
            BTVD(3,3) = 2.55228959715369D-1
            BTVD(4,1) = 5.3493423549296D-2
            BTVD(4,2) = 2.308282D-9
            BTVD(4,4) = 2.31222045731032D-1
            BTVD(5,1) = 5.3948649659203D-2
            BTVD(5,2) = 8.50587D-10
            BTVD(5,5) = 3.48990874043D-1
            
ELSEIF ((SSP1.EQ.6).AND.(SSP2.EQ.3)) THEN

            NRK = 6

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 3.12168767243578D-1
            ATVD(2,2) = 6.87831232756422D-1
            ATVD(3,1) = 7.716265806604999D-2
            ATVD(3,3) = 9.2283734193395D-1
            ATVD(4,1) = 4.12132458656463D-1
            ATVD(4,2) = 2.1717224605D-5
            ATVD(4,4) = 5.87845824118932D-1
            ATVD(5,1) = 3.73425049182364D-1
            ATVD(5,2) = 1.62494172731D-3
            ATVD(5,3) = 6.10182039447D-4
            ATVD(5,5) = 6.2433982705088D-1
            ATVD(6,1) = 1.9475670525709D-2
            ATVD(6,2) = 4.280756369499D-3
            ATVD(6,3) = 1.85954886452D-4
            ATVD(6,4) = 4.30495D-10
            ATVD(6,6) = 9.760576177878439D-1
      
            BTVD(1,1) = 3.13661711287779D-1
            BTVD(2,2) = 2.15746321543562D-1
            BTVD(3,1) = 6.5063517563946D-2
            BTVD(3,3) = 2.89458739911268D-1
            BTVD(4,1) = 9.591954428256901D-2
            BTVD(4,2) = 6.811861834D-6
            BTVD(4,4) = 1.84384727166519D-1
            BTVD(5,1) = 3.5438929973665D-2
            BTVD(5,2) = 5.09682002931D-4
            BTVD(5,3) = 1.9139074269D-4
            BTVD(5,5) = 1.95831498577895D-1
            BTVD(6,1) = 2.3383161968088D-2
            BTVD(6,2) = 1.342709368463D-3
            BTVD(6,3) = 5.8326927907D-5
            BTVD(6,4) = 1.3503D-10
            BTVD(6,6) = 3.06151902710808D-1
            
ELSEIF ((SSP1.EQ.7).AND.(SSP2.EQ.3)) THEN

            NRK = 7

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 1.47691415319417D-1
            ATVD(2,2) = 8.52308584680583D-1
            ATVD(3,1) = 7.34539305979D-4
            ATVD(3,3) = 9.99265460694021D-1
            ATVD(4,1) = 3.92877396776094D-1
            ATVD(4,2) = 1.47887230450615D-1
            ATVD(4,4) = 4.59235372773291D-1
            ATVD(5,1) = 1.97049996198846D-1
            ATVD(5,2) = 9.368160402672999D-3
            ATVD(5,3) = 8.9935546009991D-2
            ATVD(5,5) = 7.036462973884891D-1
            ATVD(6,1) = 6.1993002028958D-2
            ATVD(6,3) = 6.1395792658554D-2
            ATVD(6,4) = 4.27057986508D-4
            ATVD(6,6) = 8.7618414732598D-1
            ATVD(7,1) = 1.28534993749396D-1
            ATVD(7,3) = 1.053984977328D-2
            ATVD(7,4) = 5.607329238D-6
            ATVD(7,5) = 3.07502271596D-4
            ATVD(7,7) = 8.60612046876491D-1
      
            BTVD(1,1) = 2.580794742818D-1
            BTVD(2,2) = 2.1996335146023D-1
            BTVD(3,1) = 2.319087748D-6
            BTVD(3,3) = 2.57889904763873D-1
            BTVD(4,1) = 4.3608780145257D-2
            BTVD(4,2) = 3.8166658687686D-2
            BTVD(4,4) = 1.18519223576937D-1
            BTVD(5,1) = 7.9389033244D-5
            BTVD(5,2) = 2.41772991171D-3
            BTVD(5,3) = 2.3210518433505D-2
            BTVD(5,5) = 1.81596666510356D-1
            BTVD(6,3) = 1.5844993892434D-2
            BTVD(6,4) = 1.10214900646D-4
            BTVD(6,6) = 2.26125144115936D-1
            BTVD(7,1) = 5.3986013420437D-2
            BTVD(7,3) = 2.720118888497D-3
            BTVD(7,4) = 1.447136582D-6
            BTVD(7,5) = 7.9360024594D-5
            BTVD(7,7) = 2.22106304618468D-1
            
ELSEIF ((SSP1.EQ.8).AND.(SSP2.EQ.3)) THEN

            NRK = 8

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 1.452835644833D-3
            ATVD(2,2) = 9.98547164355167D-1
            ATVD(3,1) = 1.5824688647D-4
            ATVD(3,3) = 9.9984175311353D-1
            ATVD(4,1) = 1.35300645522171D-1
            ATVD(4,2) = 4.3910257822658D-2
            ATVD(4,4) = 8.20789096655171D-1
            ATVD(5,1) = 2.5418437348026D-1
            ATVD(5,2) = 2.40296919265655D-1
            ATVD(5,3) = 1.73107225811679D-1
            ATVD(5,5) = 3.32411481442406D-1
            ATVD(6,1) = 6.6883284219592D-2
            ATVD(6,2) = 6.766232631260299D-2
            ATVD(6,3) = 6.265014698872D-3
            ATVD(6,4) = 1.099839887838D-3
            ATVD(6,6) = 8.58089534881095D-1
            ATVD(7,1) = 1.7899811147945D-1
            ATVD(7,2) = 6.1236827D-8
            ATVD(7,3) = 2.0509168546D-4
            ATVD(7,4) = 3.5635846554D-5
            ATVD(7,5) = 1.045577373439D-3
            ATVD(7,7) = 8.197155223782711D-1
            ATVD(8,1) = 1.4239397335328D-2
            ATVD(8,2) = 4.979393069413D-3
            ATVD(8,3) = 4.9618753349D-5
            ATVD(8,4) = 4.08154071557D-4
            ATVD(8,5) = 1.95292343D-7
            ATVD(8,6) = 1.40625581752D-4
            ATVD(8,8) = 9.80182615896259D-1
      
            BTVD(1,1) = 2.19133970894904D-1
            BTVD(2,2) = 2.18815605250994D-1
            BTVD(3,1) = 3.7408889951057D-2
            BTVD(3,3) = 2.1909929362629D-1
            BTVD(4,1) = 3.63333123844D-4
            BTVD(4,2) = 9.622229159698001D-3
            BTVD(4,4) = 1.79862774017288D-1
            BTVD(5,1) = 1.4964409637913D-2
            BTVD(5,2) = 5.2657218112495D-2
            BTVD(5,3) = 3.7933673782714D-2
            BTVD(5,5) = 7.2842647899532D-2
            BTVD(6,1) = 7.217055487D-4
            BTVD(6,2) = 1.4827114244867D-2
            BTVD(6,3) = 1.372877548679D-3
            BTVD(6,4) = 2.4101228197D-4
            BTVD(6,6) = 1.88036567161855D-1
            BTVD(7,1) = 1.4581151700251D-2
            BTVD(7,2) = 1.3419069D-8
            BTVD(7,3) = 4.4942555432D-5
            BTVD(7,4) = 7.809024561999999D-6
            BTVD(7,5) = 2.29121521719D-4
            BTVD(7,7) = 1.79627517422941D-1
            BTVD(8,1) = 2.3565431064376D-2
            BTVD(8,2) = 1.091154175947D-3
            BTVD(8,3) = 1.0873154452D-5
            BTVD(8,4) = 8.9440422437D-5
            BTVD(8,5) = 4.2795187D-8
            BTVD(8,6) = 3.0815842139D-5
            BTVD(8,8) = 2.14791308823501D-1
            
ELSEIF ((SSP1.EQ.5).AND.(SSP2.EQ.4)) THEN
         
            NRK = 5
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 5.78405168674975D-1
            ATVD(2,2) = 4.21594831325025D-1
            ATVD(3,1) = 6.77678199623741D-1
            ATVD(3,2) = 8.822045619988001D-3
            ATVD(3,3) = 3.13499754756272D-1
            ATVD(4,1) = 3.60703695071165D-1
            ATVD(4,2) = 1.5747555930398D-1
            ATVD(4,3) = 3.40817234398D-4
            ATVD(4,4) = 4.81479928390457D-1
            ATVD(5,1) = 4.03281245556668D-1
            ATVD(5,2) = 1.96797187365244D-1
            ATVD(5,3) = 2.9058018283826D-2
            ATVD(5,4) = 1.46343918738107D-1
            ATVD(5,5) = 2.24519630056155D-1
      
            BTVD(1,1) = 3.63990101160008D-1
            BTVD(2,1) = 2.0871083992052D-2
            BTVD(2,2) = 4.17052806266404D-1
            BTVD(3,1) = 1.745020311D-6
            BTVD(3,2) = 8.727001873487D-3
            BTVD(3,3) = 3.10122285119134D-1
            BTVD(4,1) = 5.820402293450001D-4
            BTVD(4,2) = 1.55779006397411D-1
            BTVD(4,3) = 3.37145461634D-4
            BTVD(4,4) = 4.76292735053434D-1
            BTVD(5,1) = 3.6738984755689D-2
            BTVD(5,2) = 1.94677005403644D-1
            BTVD(5,3) = 2.8744963574915D-2
            BTVD(5,4) = 1.44767291852061D-1
            BTVD(5,5) = 2.22100782124214D-1
            
ELSEIF ((SSP1.EQ.6).AND.(SSP2.EQ.4)) THEN

            NRK = 6
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.67088614334124D-1
            ATVD(2,2) = 5.32911385665876D-1
            ATVD(3,1) = 4.08874458377822D-1
            ATVD(3,3) = 5.91125541622178D-1
            ATVD(4,1) = 4.48601090300132D-1
            ATVD(4,4) = 5.51398909699866D-1
            ATVD(5,1) = 1.12338903697D-4
            ATVD(5,2) = 8.930525958D-6
            ATVD(5,3) = 6.3196058D-8
            ATVD(5,5) = 9.99878667374287D-1
            ATVD(6,1) = 1.32061973305982D-1
            ATVD(6,2) = 1.14804153230642D-1
            ATVD(6,3) = 2.56477290949465D-1
            ATVD(6,4) = 1.16490556977394D-1
            ATVD(6,6) = 3.80166025536516D-1
      
            BTVD(1,1) = 4.42682663664825D-1
            BTVD(2,2) = 2.35910631703883D-1
            BTVD(3,3) = 2.61681029325618D-1
            BTVD(4,1) = 3.46D-13
            BTVD(4,4) = 2.44094738087817D-1
            BTVD(5,1) = 1.0465373390104D-2
            BTVD(5,2) = 3.953389019D-6
            BTVD(5,3) = 2.7975799D-8
            BTVD(5,5) = 4.42628951814885D-1
            BTVD(6,1) = 2.62732028444D-3
            BTVD(6,2) = 5.0821808351925D-2
            BTVD(6,3) = 1.13538050327048D-1
            BTVD(6,4) = 5.1568350054552D-2
            BTVD(6,6) = 1.68292908819375D-1
            
ELSEIF ((SSP1.EQ.7).AND.(SSP2.EQ.4)) THEN

            NRK = 7
            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 5.25693525536155D-1
            ATVD(2,2) = 4.74306474463845D-1
            ATVD(3,1) = 8.3349916526D-5
            ATVD(3,3) = 9.99916650083474D-1
            ATVD(4,1) = 5.24867751177608D-1
            ATVD(4,2) = 3.69729028192D-4
            ATVD(4,4) = 4.747625197942D-1
            ATVD(5,1) = 1.95877892721119D-1
            ATVD(5,2) = 4.30179D-9
            ATVD(5,3) = 3.855559176D-6
            ATVD(5,5) = 8.04118247417915D-1
            ATVD(6,1) = 6.8982360418487D-2
            ATVD(6,2) = 1.3222295413188D-2
            ATVD(6,3) = 7.59402789052D-4
            ATVD(6,4) = 2.82563389428D-4
            ATVD(6,6) = 9.16753377989845D-1
            ATVD(7,1) = 8.509651769937999D-2
            ATVD(7,2) = 7.019462801832D-3
            ATVD(7,3) = 1.51045322680285D-1
            ATVD(7,4) = 3.15632962717013D-1
            ATVD(7,5) = 2.3712178831159D-2
            ATVD(7,7) = 4.17493555270331D-1
      
            BTVD(1,1) = 3.43396366222189D-1
            BTVD(2,2) = 1.62875119806542D-1
            BTVD(3,1) = 4.313430735199D-3
            BTVD(3,3) = 3.43367744163729D-1
            BTVD(4,1) = 2.490553774367D-3
            BTVD(4,2) = 1.26963604768D-4
            BTVD(4,4) = 1.63031724115818D-1
            BTVD(5,1) = 4.766901826487D-3
            BTVD(5,2) = 1.477219D-9
            BTVD(5,3) = 1.323985011D-6
            BTVD(5,5) = 2.76131284176267D-1
            BTVD(6,1) = 4.2233750546714D-2
            BTVD(6,2) = 4.540488198005D-3
            BTVD(6,3) = 2.60776158259D-4
            BTVD(6,4) = 9.7031241157D-5
            BTVD(6,6) = 3.1480977872363D-1
            BTVD(7,1) = 3.7417046429D-5
            BTVD(7,2) = 2.410458018981D-3
            BTVD(7,3) = 5.1868414943268D-2
            BTVD(7,4) = 1.08387212456966D-1
            BTVD(7,5) = 8.142676045831D-3
            BTVD(7,7) = 1.43365769801014D-1
      
ELSEIF ((SSP1.EQ.8).AND.(SSP2.EQ.4)) THEN

            NRK = 8

            ATVD(:,:) = 0.D0
            BTVD(:,:) = 0.D0
            DTVD(:) = 0.D0
            
            ATVD(1,1) = 1.0D0
            ATVD(2,1) = 4.61181537109233D-1
            ATVD(2,2) = 5.38818462890767D-1
            ATVD(3,1) = 2.0117158505919D-2
            ATVD(3,3) = 9.79882841494082D-1
            ATVD(4,1) = 1.69428274041103D-1
            ATVD(4,2) = 5.822D-12
            ATVD(4,4) = 8.30571725953071D-1
            ATVD(5,1) = 4.31077025374928D-1
            ATVD(5,2) = 9.6737381267313D-2
            ATVD(5,3) = 8.8663023468039D-2
            ATVD(5,5) = 3.83522569889722D-1
            ATVD(6,1) = 8.9850295244391D-2
            ATVD(6,2) = 5.5129D-11
            ATVD(6,3) = 1.2320078306704D-2
            ATVD(6,4) = 2.95511D-10
            ATVD(6,6) = 8.97829626098266D-1
            ATVD(7,1) = 6.244345604206D-3
            ATVD(7,2) = 4.339895694D-4
            ATVD(7,3) = 3.515550716564D-3
            ATVD(7,4) = 9.393008934879D-3
            ATVD(7,5) = 1.339131D-9
            ATVD(7,7) = 9.8041310383582D-1
            ATVD(8,1) = 7.1188337795849D-2
            ATVD(8,2) = 3.4167417741535D-2
            ATVD(8,3) = 1.632796D-9
            ATVD(8,4) = 2.63316641735748D-1
            ATVD(8,5) = 1.51537858712815D-1
            ATVD(8,6) = 1.14976594D-7
            ATVD(8,8) = 4.79789627404664D-1
      
            BTVD(1,1) = 2.78271529217822D-1
            BTVD(2,2) = 1.4993783763941D-1
            BTVD(3,1) = 1.6659D-11
            BTVD(3,3) = 2.72673496756863D-1
            BTVD(4,1) = 2.4301359D-8
            BTVD(4,2) = 1.62D-12
            BTVD(4,4) = 2.31124464306047D-1
            BTVD(5,1) = 1.267923D-9
            BTVD(5,2) = 2.6919259017783D-2
            BTVD(5,3) = 2.4672395125527D-2
            BTVD(5,5) = 1.06723412012762D-1
            BTVD(6,1) = 3.2744D-11
            BTVD(6,2) = 1.5341D-11
            BTVD(6,3) = 3.42832703049D-3
            BTVD(6,4) = 8.2232D-11
            BTVD(6,6) = 2.4984042303143D-1
            BTVD(7,1) = 3.4914608275757D-2
            BTVD(7,2) = 1.20766941141D-4
            BTVD(7,3) = 9.782776739410001D-4
            BTVD(7,4) = 2.613806960265D-3
            BTVD(7,5) = 3.72642D-10
            BTVD(7,7) = 2.72821053669585D-1
            BTVD(8,1) = 4.265D-11
            BTVD(8,2) = 9.507819584361D-3
            BTVD(8,3) = 4.54361D-10
            BTVD(8,4) = 7.3273524564308D-2
            BTVD(8,5) = 4.2168671678409D-2
            BTVD(8,6) = 3.1994713D-8
            BTVD(8,8) = 1.33511793320745D-1
            
ENDIF

IF (order > 1) THEN
    alphaW = ATVD
    betaW  = BTVD
ENDIF    

IF (debug == 4) THEN
    DO j = 1, S
        write(*,*)betaW(:,j)
	write(*,*)'-------'
    ENDDO
ENDIF

!....ct coefficients for exact solution in time

DO i = 0, S-1
    DO j = 1, S
	csum = 0.0d0
	DO ii = i+1, j-1
	    csum = csum + crk(ii,i+1)*alphaW(j,ii+1)
	ENDDO
	crk(j,i+1) = betaW(j,i+1) + csum
    ENDDO
ENDDO

DO i = 1, S-1
    ctW(i+1) = 0.0d0
    DO j = 0, i-1
        ctW(i+1) = ctW(i+1) + crk(i,j+1)
    ENDDO
ENDDO

IF (debug == 4) THEN
	write(*,*)ctW
ENDIF
 
RETURN
END SUBROUTINE WRK_COEFFICIENTS
!------------------------------------------------------------------------------------
!
!                 Reads in 1D mesh Subroutine
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                           5.20.14
!
!-----------------------------------------------------------------------------------
SUBROUTINE READ_FORT_1D_MESH()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p), PARAMETER :: debug_check = 0
INTEGER(int_p) :: mesh1d

OPEN(10,FILE='fort1D.14')
! Read in Mesh file (# of Elements and Nodes)
READ(10,'(A24)') mesh1d
READ(10,*) nnodesX, nnodesY, nnodesZ

ALLOCATE(X(nnodesX))
ALLOCATE(Y(nnodesY))
ALLOCATE(Z(nnodesZ))
ALLOCATE(XF(nnodesX))

X  = 0.0d0
Y  = 0.0d0
Z  = 0.0d0
XF = X

DO i = 1, nnodesX
    READ(10,*) ii, X(i)
ENDDO

DO i = 1, nnodesY
    READ(10,*) ii, Y(i)
ENDDO

DO i = 1, nnodesZ
    READ(10,*) ii, Z(i)
ENDDO

CLOSE(10)

PROB_DATA%amp    = 0.500000000000000d0
PROB_DATA%freq   = 0.000140750000000d0
PROB_DATA%phase  = 0.000000000000000d0
PROB_DATA%k      = 0.010000000000000d0
PROB_DATA%g      = 9.806650000000000d0
PROB_DATA%N0     = 0.001500000000000d0
PROB_DATA%h0     = 1.000000000000000d0
PROB_DATA%x0     = X(1)
PROB_DATA%xN     = X(nnodesX)
PROB_DATA%y0     = Y(1)
PROB_DATA%yN     = Y(nnodesY) 
PROB_DATA%Hslope = 0.000000000000000d0 
PROB_DATA%ampU   = 0.1000000000000000d0

nelemsX  = nnodesX - 1
nelemsY  = nnodesY - 1
nelemsZ  = nnodesZ - 1
nelemsXF = nelemsX
nnodesXF = nnodesX

ALLOCATE(xELEM_nodes(nelemsX,2))
ALLOCATE(yELEM_nodes(nelemsY,2))
ALLOCATE(zELEM_nodes(nelemsZ,2))
ALLOCATE(xNODE_elems(nnodesX,2))
ALLOCATE(yNODE_elems(nnodesY,2))
ALLOCATE(zNODE_elems(nnodesZ,2))
ALLOCATE(xELEM_jac(nelemsX))
ALLOCATE(yELEM_jac(nelemsY))
ALLOCATE(zELEM_jac(nelemsZ))
ALLOCATE(xfELEM_nodes(nelemsXF,2))
ALLOCATE(xfNODE_elems(nnodesXF,2))
ALLOCATE(xfELEM_jac(nelemsXF))

!....Initialize Matricies

xELEM_nodes  = 0.0d0
yELEM_nodes  = 0.0d0
zELEM_nodes  = 0.0d0
xfELEM_nodes = 0.0d0
xNODE_elems  = 0.0d0
yNODE_elems  = 0.0d0
zNODE_elems  = 0.0d0
xfNODE_elems = 0.0d0
xELEM_jac    = 0.0d0
yELEM_jac    = 0.0d0
zELEM_jac    = 0.0d0
xfELEM_jac   = 0.0d0

!....Connectivity and Jacobians

DO i = 1, nelemsX 
    xELEM_nodes(i,1) = i
    xELEM_nodes(i,2) = i+1
    xELEM_jac(i) = 2.0d0/(X(i+1)-X(i)) 
ENDDO

DO i = 1, nelemsXF 
    xfELEM_nodes(i,1) = i
    xfELEM_nodes(i,2) = i+1
    xfELEM_jac(i) = 2.0d0/(XF(i+1)-XF(i))
ENDDO

DO i = 1, nelemsY
    yELEM_nodes(i,1) = i
    yELEM_nodes(i,2) = i+1
    yELEM_jac(i) = 2.0d0/(Y(i+1)-Y(i))
ENDDO    

DO i = 1, nelemsZ
    zELEM_nodes(i,1) = i
    zELEM_nodes(i,2) = i+1
    zELEM_jac(i) = 2.0d0/(Z(i+1)-Z(i))
ENDDO

!....Set up node element connectivity for periodic BCs

xNODE_elems(1,1) = nelemsX
xNODE_elems(1,2) = 1
DO i = 2, nnodesX-1
    xNODE_elems(i,1) = i-1
    xNODE_elems(i,2) = i
ENDDO
xNODE_elems(nnodesX,1) = nelemsX
xNODE_elems(nnodesX,2) = 1

xfNODE_elems(1,1) = nelemsXF 
xfNODE_elems(1,2) = 1
DO i = 2, nnodesXF-1
    xfNODE_elems(i,1) = i-1
    xfNODE_elems(i,2) = i
ENDDO
xfNODE_elems(nnodesXF,1) = nelemsXF
xfNODE_elems(nnodesXF,2) = 1

yNODE_elems(1,1) = nelemsY
yNODE_elems(1,2) = 1
DO i = 2, nnodesY-1
    yNODE_elems(i,1) = i-1
    yNODE_elems(i,2) = i
ENDDO
yNODE_elems(nnodesY,1) = nelemsY
yNODE_elems(nnodesY,2) = 1

zNODE_elems(1,1) = nelemsZ
zNODE_elems(1,2) = 1
DO i = 2, nnodesZ-1
    zNODE_elems(i,1) = i-1
    zNODE_elems(i,2) = i
ENDDO
zNODE_elems(nnodesZ,1) = nelemsZ
zNODE_elems(nnodesZ,2) = 1

IF (debug_check == 1) THEN
    OPEN(10,FILE='FORT1Dcheck.14')
    ! Read in Mesh file (# of Elements and Nodes)
    WRITE(10,'(A24)') mesh1d
    WRITE(10,*) nnodesX, nnodesY, nnodesZ

    DO i = 1, nnodesX
        WRITE(10,*) i, X(i)
    ENDDO

    DO i = 1, nnodesY
        WRITE(10,*) i, Y(i)
    ENDDO

    DO i = 1, nnodesZ
        WRITE(10,*) i, Z(i)
    ENDDO

    CLOSE(10)
ENDIF

END SUBROUTINE READ_FORT_1D_MESH
!------------------------------------------------------------------------------------
!
!                 Reads in 2D mesh Subroutine
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                           4.21.14
!
!-----------------------------------------------------------------------------------
SUBROUTINE READ_FORT_2D_MESH()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p), PARAMETER :: debug_check = 1
INTEGER(int_p) :: mesh2d, nedgesX, nedgesY
REAL(real_p) :: dumb, b, h
REAL(real_p), DIMENSION(4) :: node

OPEN(10,FILE='fort2D.14')
! Read in Mesh file (# of Elements and Nodes)
READ(10,'(A24)') mesh2d
READ(10,*) Qelems, Qnnodes

ALLOCATE(CONN(Qelems,4))
ALLOCATE(Qnode(Qnnodes,2))
ALLOCATE(CONNielem(Qelems,2))
ALLOCATE(QedgeX(Qelems,2))
ALLOCATE(QedgeY(Qelems,2))

CONN      = 0.0d0
Qnode     = 0.0d0
CONNielem = 0.0d0
QedgeX    = 0.0d0 
QedgeY    = 0.0d0

DO i = 1, Qnnodes
    READ(10,*) ii, Qnode(i,1), Qnode(i,2), dumb
ENDDO

DO i = 1, Qelems
    READ(10,*) ii, CONN(i,1), CONN(i,2), CONN(i,3), CONN(i,4)
ENDDO

DO i = 1, Qelems
    READ(10,*) ii, CONNielem(i,1), CONNielem(i,2)
ENDDO

DO i = 1, Qelems
    READ(10,*) ii, QedgeX(i,1), QedgeX(i,2)
ENDDO

DO i = 1, Qelems
    READ(10,*) ii, QedgeY(i,1), QedgeY(i,2)
ENDDO

READ(10,*)nelemsZ
ALLOCATE(HEXcolumn(nelemsZ,Qelems))
HEXcolumn = 0
DO i = 1, Qelems
    DO j = 1, nelemsZ
        READ(10,*) ii, HEXcolumn(j,i)
    ENDDO    
ENDDO

READ(10,*)nedgesX
ALLOCATE(xEDGE(nedgesX,2))
xEDGE = 0
DO i = 1, nedgesX
    READ(10,*) ii, xEDGE(i,1), xEDGE(i,2)
ENDDO

READ(10,*)nedgesY
ALLOCATE(yEDGE(nedgesY,2))
yEDGE = 0
DO i = 1, nedgesY
    READ(10,*) ii, yEDGE(i,1), yEDGE(i,2)
ENDDO    

READ(10,*)NEedgesX
ALLOCATE(ExEDGE(NEedgesX))
ExEDGE = 0
DO i = 1, NEedgesX
    READ(10,*) ii, ExEDGE(i)
ENDDO  

READ(10,*)NEedgesY
ALLOCATE(EyEDGE(NEedgesY))
EyEDGE = 0
DO i = 1, NEedgesY
    READ(10,*) ii, EyEDGE(i)
ENDDO  

CLOSE(10)

ALLOCATE(CONN_jac(Qelems,1))
ALLOCATE(CONN_jacY(Qelems,1))
ALLOCATE(CONN_jacX(Qelems,1))
ALLOCATE(CONN_jacF(Qelems,1))
ALLOCATE(CONN_jacS(Qelems,1))
ALLOCATE(CONNsigma(Qelems,1))
ALLOCATE(CONNgrid(Qelems,1))
ALLOCATE(CONNzx(Qelems,1))

CONN_jac  = 0.0d0
CONN_jacY = 0.0d0
CONN_jacX = 0.0d0
CONN_jacF = 0.0d0
CONN_jacS = 0.0d0
CONNsigma = 0.0d0
CONNgrid  = 0.0d0
CONNzx    = 0.0d0
 
DO kk = 1, Qelems   
    node = CONN(kk,:)
    b    = Qnode(node(4),1) - Qnode(node(1),1)  !....Make sure this is grabbing the correct values
    h    = Qnode(node(2),2) - Qnode(node(1),2)
    CONN_jac(kk,1)  = 4.0d0/(b*h)
    CONN_jacY(kk,1) = 2.0d0/h
    CONN_jacX(kk,1) = 2.0d0/b
    CONN_jacF(kk,1) = CONN_jacY(kk,1)
    CONN_jacS(kk,1) = h/b
ENDDO

IF (debug_check == 1) THEN
    OPEN(10,FILE='FORT2Dcheck.14')
    WRITE(10,'(A24)') mesh2d
    WRITE(10,*) Qelems, Qnnodes

    DO i = 1, Qnnodes
        WRITE(10,*) i, Qnode(i,1), Qnode(i,2), 0.0d0
    ENDDO

    DO i = 1, Qelems
        WRITE(10,*) i, CONN(i,1), CONN(i,2), CONN(i,3), CONN(i,4)
    ENDDO

    DO i = 1, Qelems
        WRITE(10,*) i, CONNielem(i,1), CONNielem(i,2)
    ENDDO

    DO i = 1, Qelems
        WRITE(10,*) i, QedgeX(i,1), QedgeX(i,2) 
    ENDDO

    DO i = 1, Qelems
        WRITE(10,*) i, QedgeY(i,1), QedgeY(i,2)
    ENDDO

    WRITE(10,*)nelemsZ
    ii = 0
    DO i = 1, Qelems
        DO j = 1, nelemsZ
            ii = ii + 1
            WRITE(10,*) ii, HEXcolumn(j,i)
        ENDDO    
    ENDDO

    WRITE(10,*)nedgesX
    DO i = 1, nedgesX
        WRITE(10,*) i, xEDGE(i,1), xEDGE(i,2)
    ENDDO

    WRITE(10,*)nedgesY
    DO i = 1, nedgesY
        WRITE(10,*) i, yEDGE(i,1), yEDGE(i,2)
    ENDDO    
    
    WRITE(10,*)NEedgesX
    DO i = 1, NEedgesX
        WRITE(10,*) i, ExEDGE(i)
    ENDDO  

    WRITE(10,*)NEedgesY
    DO i = 1, NEedgesY
        WRITE(10,*) i, EyEDGE(i)
    ENDDO     

    CLOSE(10)
ENDIF

RETURN
END SUBROUTINE READ_FORT_2D_MESH
!------------------------------------------------------------------------------------
!
!                 Reads in 3D mesh Subroutine
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                           4.21.14
!
!-----------------------------------------------------------------------------------
SUBROUTINE READ_FORT_3D_MESH()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p), PARAMETER :: debug_check = 1
INTEGER(int_p) :: mesh3d, nfacesX, nfacesY, nfacesZ  !....Look at vertical mesh for W and see if need to do anything for that
REAL(real_p) :: dumb, b, h, w
REAL(real_p), DIMENSION(8) :: node

OPEN(10,FILE='fort3D.14')
! Read in Mesh file (# of Elements and Nodes)
READ(10,'(A24)') mesh3d
READ(10,*) HEXelems, HEXnnodes

ALLOCATE(HEXnode(HEXnnodes,3))
ALLOCATE(HEXconn(HEXelems,8))
ALLOCATE(HEXzxy(HEXelems,1))
ALLOCATE(HEXfaceX(HEXelems,2))
ALLOCATE(HEXfaceY(HEXelems,2))
ALLOCATE(HEXfaceZ(HEXelems,2))

HEXnode  = 0
HEXconn  = 0
HEXzxy   = 0
HEXfaceX = 0
HEXfaceY = 0
HEXfaceZ = 0

DO i = 1, HEXnnodes
    READ(10,*) ii, HEXnode(i,1), HEXnode(i,2), HEXnode(i,3)
ENDDO

DO i = 1, HEXelems
    READ(10,*) ii, HEXconn(i,1), HEXconn(i,2), HEXconn(i,3), HEXconn(i,4), &
                   HEXconn(i,5), HEXconn(i,6), HEXconn(i,7), HEXconn(i,8)
ENDDO                   

DO i = 1, HEXelems
    READ(10,*) ii, HEXzxy(i,1)
ENDDO

DO i = 1, HEXelems
    READ(10,*) ii, HEXfaceX(i,1), HEXfaceX(i,2)
ENDDO

DO i = 1, HEXelems
    READ(10,*) ii, HEXfaceY(i,1), HEXfaceY(i,2)
ENDDO

DO i = 1, HEXelems
    READ(10,*) ii, HEXfaceZ(i,1), HEXfaceZ(i,2)
ENDDO    

READ(10,*)nfacesZ
ALLOCATE(zFACE(nfacesZ,2))
zFACE = 0.0d0
DO i = 1, nfacesZ
    READ(10,*) ii, zFACE(i,1), zFACE(i,2)
ENDDO

READ(10,*)nfacesX
ALLOCATE(xFACE(nfacesX,2))
xFACE = 0.0d0
DO i = 1, nfacesX
    READ(10,*) ii, xFACE(i,1), xFACE(i,2)
ENDDO

READ(10,*)nfacesY
ALLOCATE(yFACE(nfacesY,2))
yFACE = 0.0d0
DO i = 1, nfacesY
    READ(10,*) ii, yFACE(i,1), yFACE(i,2)
ENDDO

READ(10,*)NEfacesX
ALLOCATE(ExFACE(NEfacesX))
ExFACE = 0.0d0
DO i = 1, NEfacesX
    READ(10,*) ii, ExFACE(i)
ENDDO

READ(10,*)NEfacesY
ALLOCATE(EyFACE(NEfacesY))
EyFACE = 0.0d0
DO i = 1, NEfacesY
    READ(10,*) ii, EyFACE(i)
ENDDO

CLOSE(10)

ALLOCATE(HEX_jac(HEXelems,1))
ALLOCATE(HEXsigma(HEXelems,1))
ALLOCATE(HEXgrid(nelemsZ,nelemsX,nelemsY))   !....need to fix        
ALLOCATE(HEX_jacZ(HEXelems,1))
ALLOCATE(HEX_jacX(HEXelems,1))
ALLOCATE(HEX_jacY(HEXelems,1))
ALLOCATE(HEX_jacZw(HEXelems,1))
ALLOCATE(HEX_jacXw(HEXelems,1))
ALLOCATE(HEX_jacYw(HEXelems,1))
ALLOCATE(HEX_jacF(HEXelems,1))
ALLOCATE(HEX_jacS(HEXelems,1))

HEXgrid   = 0
HEXsigma  = 0
HEX_jac   = 0.0d0
HEX_jacZ  = 0.0d0
HEX_jacY  = 0.0d0
HEX_jacX  = 0.0d0
HEX_jacZw = 0.0d0
HEX_jacYw = 0.0d0
HEX_jacXw = 0.0d0
HEX_jacS  = 0.0d0

DO kk = 1, HEXelems
    
    node      = HEXconn(kk,:)
    b         = HEXnode(node(4),1) - HEXnode(node(1),1) !....check these reference correct nodes
    h         = HEXnode(node(2),3) - HEXnode(node(1),3)
    w         = HEXnode(node(5),2) - HEXnode(node(1),2)
            
    HEX_jac(kk,1)  = 8.0d0/(b*h*w)
    HEX_jacZ(kk,1) = 2.0d0/h
    HEX_jacX(kk,1) = 2.0d0/b
    HEX_jacY(kk,1) = 2.0d0/w
            
    HEX_jacZw(kk,1) = 1.0d0
    HEX_jacXw(kk,1) = h/b
    HEX_jacYw(kk,1) = h/w        
    
ENDDO

IF (debug_check == 1) THEN
    
    OPEN(10,FILE='FORT3Dcheck.14')
    WRITE(10,'(A24)') mesh3d
    WRITE(10,*) HEXelems, HEXnnodes

    DO i = 1, HEXnnodes
        WRITE(10,*) ii, HEXnode(i,1), HEXnode(i,2), HEXnode(i,3)
    ENDDO

    DO i = 1, HEXelems
        WRITE(10,*) i, HEXconn(i,1), HEXconn(i,2), HEXconn(i,3), HEXconn(i,4), &
                   HEXconn(i,5), HEXconn(i,6), HEXconn(i,7), HEXconn(i,8)
    ENDDO                   

    DO i = 1, HEXelems
        WRITE(10,*) i, HEXzxy(i,1)
    ENDDO

    DO i = 1, HEXelems
        WRITE(10,*) i, HEXfaceX(i,1), HEXfaceX(i,2)
    ENDDO

    DO i = 1, HEXelems
        WRITE(10,*) i, HEXfaceY(i,1), HEXfaceY(i,2)
    ENDDO

    DO i = 1, HEXelems
        WRITE(10,*) i, HEXfaceZ(i,1), HEXfaceZ(i,2)
    ENDDO    

    WRITE(10,*)nfacesZ
    DO i = 1, nfacesZ
        WRITE(10,*) i, zFACE(i,1), zFACE(i,2)
    ENDDO

    WRITE(10,*)nfacesX
    DO i = 1, nfacesX
        WRITE(10,*) i, xFACE(i,1), xFACE(i,2)
    ENDDO

    WRITE(10,*)nfacesY
    DO i = 1, nfacesY
        WRITE(10,*) i, yFACE(i,1), yFACE(i,2)
    ENDDO
    
    WRITE(10,*)NEfacesX
    DO i = 1, NEfacesX
        WRITE(10,*) ii, ExFACE(i)
    ENDDO

    WRITE(10,*)NEfacesY
    DO i = 1, NEfacesY
        WRITE(10,*) ii, EyFACE(i)
    ENDDO    

    CLOSE(10)
ENDIF


RETURN
END SUBROUTINE READ_FORT_3D_MESH
!------------------------------------------------------------------------------------
!
!                 Reads in fort 15 Subroutine
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                           4.21.14
!
!-----------------------------------------------------------------------------------
SUBROUTINE READ_FORT_15()

USE precisions
USE global
IMPLICIT NONE

INTEGER(int_p), PARAMETER :: debug_check = 1
INTEGER(int_p) :: Porder, k, PTS_15_2D

!....Read in input data 

IF (pu == 0) THEN
    OPEN(10,FILE='fort3D15_P0.15')
ELSEIF (pu == 1) THEN 
    OPEN(10,FILE='fort3D15_P1.15')
ELSEIF (pu == 2) THEN    
    OPEN(10,FILE='fort3D15_P2.15')
ELSEIF (pu == 3) THEN
    OPEN(10,FILE='fort3D15_P3.15')
ELSE
    WRITE(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*)'~ No 3D.15 file exists for this polynomial approx. ~'
    WRITE(*,*)'~ Create the proper file with ADmesh Quad.         ~'   
    WRITE(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    STOP     
ENDIF

!....Read in open ocean data

READ(10,*) Porder, NCONST

READ(10,*)NEfacesX
ALLOCATE(EAx(Porder,NEfacesX,NCONST))
ALLOCATE(EPx(Porder,NEfacesX,NCONST))
EAx = 0.0d0
EPx = 0.0d0
DO i = 1, NCONST
    DO j = 1, NEfacesX
        DO k = 1, Porder
            READ(10,*) ii, EAx(k,j,i), EPx(k,j,i)
        ENDDO
    ENDDO        
ENDDO

READ(10,*)NEfacesY
ALLOCATE(EAy(Porder,NEfacesY,NCONST))
ALLOCATE(EPy(Porder,NEfacesY,NCONST))
EAy = 0.0d0
EPy = 0.0d0
DO i = 1, NCONST
    DO j = 1, NEfacesY
        DO k = 1, Porder
            READ(10,*) ii, EAy(k,j,i), EPy(k,j,i)
        ENDDO
    ENDDO        
ENDDO

CLOSE(10)

IF (p == 0) THEN
    OPEN(11,FILE='fort2D15_P0.15')
ELSEIF (p == 1) THEN 
    OPEN(11,FILE='fort2D15_P1.15')
ELSEIF (p == 2) THEN    
    OPEN(11,FILE='fort2D15_P2.15')
ELSEIF (p == 3) THEN
    OPEN(11,FILE='fort2D15_P3.15')
ELSE
    WRITE(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*)'~ No 2D.15 file exists for this polynomial approx. ~'
    WRITE(*,*)'~ Create the proper file with ADmesh Quad.         ~'   
    WRITE(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    STOP     
ENDIF

READ(11,*)PTS_15_2D

READ(11,*)NEedgesX
ALLOCATE(EAx2D(PTS_15_2D,NEedgesX,NCONST))
ALLOCATE(EPx2D(PTS_15_2D,NEedgesX,NCONST))
EAx2D = 0.0d0
EPx2D = 0.0d0
DO i = 1, NCONST
    DO j = 1, NEedgesX
        DO k = 1, PTS_15_2D
            READ(11,*) ii, EAx2D(k,j,i), EPx2D(k,j,i)
        ENDDO
    ENDDO        
ENDDO

READ(11,*)NEedgesY
ALLOCATE(EAy2D(Porder,NEedgesY,NCONST))
ALLOCATE(EPy2D(Porder,NEedgesY,NCONST))
EAy2D = 0.0d0
EPy2D = 0.0d0
DO i = 1, NCONST
    DO j = 1, NEedgesY
        DO k = 1, PTS_15_2D
            READ(11,*) ii, EAy2D(k,j,i), EPy2D(k,j,i)
        ENDDO
    ENDDO        
ENDDO

CLOSE(11)

IF (debug_check == 1) THEN

    OPEN(10,FILE='FORT3D15check.15')

    WRITE(10,*) Porder, NCONST

    WRITE(10,*)NEfacesX

    DO i = 1, NCONST
        DO j = 1, NEfacesX
            DO k = 1, Porder
                WRITE(10,*) i, EAx(k,j,i), EPx(k,j,i)
            ENDDO
        ENDDO        
    ENDDO

    WRITE(10,*)NEfacesY

    DO i = 1, NCONST
        DO j = 1, NEfacesY
            DO k = 1, Porder
                WRITE(10,*) i, EAy(k,j,i), EPy(k,j,i)
            ENDDO
        ENDDO        
    ENDDO
    
    WRITE(10,*)PTS_15_2D

    WRITE(10,*)NEedgesX

    DO i = 1, NCONST
        DO j = 1, NEedgesX
            DO k = 1, PTS_15_2D
                WRITE(10,*) ii, EAx2D(k,j,i), EPx2D(k,j,i)
            ENDDO
        ENDDO        
    ENDDO

    WRITE(10,*)NEedgesY

    DO i = 1, NCONST
        DO j = 1, NEedgesY
            DO k = 1, PTS_15_2D
                WRITE(10,*) ii, EAy2D(k,j,i), EPy2D(k,j,i)
            ENDDO
        ENDDO        
    ENDDO

    CLOSE(10)    
    
    
ENDIF

RETURN
END SUBROUTINE READ_FORT_15
!------------------------------------------------------------------------------------
!
!            Subroutine for Calculation of Chebysev Roots
!                      Written by Colton J. Conroy
!                           @ the C.H.I.L
!                                9.6.13
!
!------------------------------------------------------------------------------------
SUBROUTINE CHEBY_ROOTS(n,a,b,xk)

!   n = number of elements
!   a = x0
!   b = xN
!   xk = vector of node locations

USE precisions
IMPLICIT NONE

INTEGER(int_p) :: j, k, n
REAL(real_p) :: top, bot, tk, a, b, C1, C2
REAL(real_p), DIMENSION(n+1) :: xk
REAL(real_p), PARAMETER :: PI = 3.1415926535897932d0

xk = 0.0d0
C1 = (b-a)/2.0d0
C2 = (b+a)/2.0d0

DO k = 0, n
    top   = (2.0d0*n + 1.0d0 - 2.0d0*k)*PI
    bot   = 2.0d0*n + 2.0d0
    tk    = COS(top/bot)
    xk(k+1) = C1*tk + C2
ENDDO
  
RETURN
END SUBROUTINE CHEBY_ROOTS
!------------------------------------------------------------------------------------
!
!                 Gaussian Elimination Subroutine
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                            5.22.13
!
!-----------------------------------------------------------------------------------
SUBROUTINE GAUSS_ELIMINATION(A,x,b,n)

USE precisions
IMPLICIT NONE


INTEGER(int_p) :: i, j, k, n, m
COMPLEX(16), DIMENSION(n,n) :: A, LS
COMPLEX(16), DIMENSION(n) :: x, b 
COMPLEX(16) :: akj, ajj, sum

A  = 0.0d0
LS = 0.0d0
x  = 0.0d0
b  = 0.0d0

DO j = 1, n-1
    DO k = j, n
        akj     = A(k,j)
        ajj     = A(j,j)
        LS(k,j) = -akj/ajj
    ENDDO
    DO k = 1, n
        IF (j < k .and. j <= n) THEN
            b(k) = b(k) + LS(k,j)*b(j)
        ENDIF
        DO m = 1, n
            IF (j < k .and. j <= m) THEN
                A(k,m) = A(k,m) + LS(k,j)*A(j,m)
            ENDIF
        ENDDO
    ENDDO
ENDDO

x(n) = b(n)/A(n,n)
m = n-1
DO i = m,1,-1     
    sum = dcmplx(0.0d0,0.0d0)
    k = i+1
    DO j = k, n
        sum = sum + A(i,j)*x(j)
    ENDDO
    x(i) = (b(i)-sum)/A(i,i)
ENDDO     

END SUBROUTINE GAUSS_ELIMINATION
!------------------------------------------------------------------------------------
!
!                 Tidal Constituent Subroutine
!                 Written by Colton J. Conroy
!                        @ the C.H.I.L
!                            9.17.14
!
!-----------------------------------------------------------------------------------
SUBROUTINE TIDAL_CONST()

USE precisions
USE global
IMPLICIT NONE

REAL(real_p), DIMENSION(10) :: AMIG_FULL

ALLOCATE(AMIG(NCONST))
AMIG = 0.0d0

!....AMIG_FULL 

AMIG_FULL(1)  = 7.2921158358d-5  !....K1
AMIG_FULL(2)  = 6.7597744151d-5  !....O1
AMIG_FULL(3)  = 6.495787d-5      !....Q1
AMIG_FULL(4)  = 1.40518902509d-4 !....M2
AMIG_FULL(5)  = 1.45444104333d-4 !....S2
AMIG_FULL(6)  = 1.37879699487d-4 !....N2
AMIG_FULL(7)  = 1.458408d-4      !....K2
AMIG_FULL(8)  = 2.81037805434d-4 !....M4
AMIG_FULL(9)  = 4.21556708118d-4 !....M6
AMIG_FULL(10) = 7.252219d-5      !....P1

DO i = 1, NCONST
    AMIG(i) = AMIG_FULL(i)
ENDDO    

END SUBROUTINE TIDAL_CONST
