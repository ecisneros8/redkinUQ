MODULE ceqct

  USE ceq_system
  IMPLICIT NONE
  integer             :: ncons  ! number of constraints
  integer             :: nspec  ! number of species
  real*8              :: patm   ! mixture pressure [atm]
  real*8              :: hor    ! mixture enthalpy over gas constant [K]
  real*8, allocatable :: wts(:) ! molecular weights 
  type(sys_type), pointer, save :: sys

CONTAINS

  SUBROUTINE ceqinit(ncs,ne,ng,ns,csi,T0,p0,h0,mw,x0,ther,Ein,Bg,Tceq,zceq)

    ! input/output variables
    integer, intent(in)    :: ncs
    integer, intent(in)    :: ne
    integer, intent(in)    :: ng
    integer, intent(in)    :: ns
    integer, intent(in)    :: csi(ncs)   
    real*8,  intent(in)    :: T0
    real*8,  intent(in)    :: p0
    real*8,  intent(inout) :: h0   
    real*8,  intent(in)    :: mw(ns)
    real*8,  intent(in)    :: x0(ns)
    real*8,  intent(in)    :: ther(ns,15)
    real*8,  intent(in)    :: Ein(ns,ne)
    real*8,  intent(in)    :: Bg(ns,ng)
    real*8,  intent(out)   :: Tceq
    real*8,  intent(out)   :: zceq(ns)
    ! local variables
    integer                :: i, s, e
    integer                :: astat
    integer                :: ierr
    integer                :: lu_op
    integer                :: diag
    integer                :: cs(ncs)
    real*8                 :: z0(ns)
    real*8,  allocatable   :: cg(:,:)
    real*8,  allocatable   :: thermo(:,:)
    real*8                 :: stats(20)

    !=======================================================!
    ! set number of constraints, species
    ncons = ncs + ne + ng
    nspec = ns

    !=======================================================!
    ! allocate
    if(allocated(wts))    deallocate(wts)
    if(allocated(cg))     deallocate(cg)
    if(allocated(thermo)) deallocate(thermo)

    allocate(&
         wts(ns),       &
         cg(ns,ng),     &
         thermo(ns,15), &
         stat=astat)
    if(astat .ne. 0) stop

    !=======================================================!
    ! save p, H/R and molecular weights  
    patm  = p0
    hor   = h0
    wts   = mw

    write(*,*)
    write(*,'(a,1p,2e13.4)') 'H/R, p0 = ', hor, patm
    write(*,*)

    !=======================================================!
    ! set some ceq parameters
    lu_op = 0
    diag  = 0

    !=======================================================!
    ! initial state  
    cs    = csi + 1
    z0    = x0  + 1.0d-15    
    z0    = z0 / sum(z0)
    z0    = z0 / sum(z0 * wts)
    
    !=======================================================!
    ! get thermo data  
    thermo(:,1)    = ther(:,1)
    thermo(:,2:8)  = ther(:,9:15)
    thermo(:,9:15) = ther(:,2:8)
    
    !=======================================================!
    ! get initial condition
    if(ng .gt. 0) then
       cg = Bg
    end if
    
    call ceq_sys_init(         &
         ns, ne, ncs, ng,      &
         Ein, cs,              &
         cg, thermo,           &
         lu_op, diag, sys,     &
         ierr)
    
    call ceq_state(                         &
         sys, N=z0, p_atm=patm, N_h=z0,     &
         T_h=T0, HoR_eq=hor, N_eq=zceq,     &
         T_eq=Tceq, stats=stats, info=ierr)

  END SUBROUTINE ceqinit
  
  SUBROUTINE ceqRecon(r,T,x,flag)

    !% input/output variables
    real*8,  intent(inout) :: r(ncons)
    real*8,  intent(out)   :: T
    real*8,  intent(out)   :: x(nspec)
    integer, intent(out)   :: flag
    !% local variables
    integer                :: i
    integer                :: ierr
    real*8                 :: mmw
    real*8                 :: z(nspec)
    
    !=======================================================!
    ! use ceq to reconstruct full state
    call ceq_param_set(sys, T_low = 100.d0)
    call ceq_state(               &
         sys, c=r, p_atm=patm,    &
         HOR=hor, N_eq=z, T_eq=T, &
         info=ierr)

    !=======================================================!
    ! convert to mass fractions
    mmw = 1.d0 / sum(z)
    x = z * mmw

    !=======================================================!
    ! if ceq failed, set flag to 1
    if(ierr .lt. 0) then
        flag = 1
    else      
        flag = 0
    end if    

  END SUBROUTINE ceqRecon

END MODULE ceqct
