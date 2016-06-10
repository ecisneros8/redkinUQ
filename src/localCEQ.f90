MODULE localCEQ

  USE ceq_system
  IMPLICIT NONE
  type(sys_type), pointer, save :: sys

CONTAINS

  SUBROUTINE ceqinit(ncs,ne,ng,ns,nc,csi,Ti,hi,mw,xi,Ein,Tad,Teq,zeq)

    ! input/output variables
    integer, intent(in)    :: ncs
    integer, intent(in)    :: ne
    integer, intent(in)    :: ng
    integer, intent(in)    :: ns
    integer, intent(in)    :: nc
    integer, intent(in)    :: csi(ncs)
    real*8,  intent(in)    :: Ti
    real*8,  intent(inout) :: hi 
    real*8,  intent(in)    :: mw(ns)
    real*8,  intent(in)    :: xi(ns)
    real*8,  intent(in)    :: Ein(ns,ne)
    real*8,  intent(out)   :: Tad
    real*8,  intent(out)   :: Teq
    real*8,  intent(out)   :: zeq(ns)
    ! local variables
    integer                :: i
    integer                :: astat
    integer                :: ierr
    integer                :: lu_op
    integer                :: diag
    integer                :: CS(ncs)
    real*8                 :: zi(ns)
    real*8,  allocatable   :: Bg(:,:)
    real*8,  allocatable   :: thermo(:,:)
    real*8                 :: stats(20)
    
    lu_op = 0
    diag  = 0
    CS    = csi + 1
    zi    = xi  + 1.0d-15
    zi    = zi/sum(zi)
    zi    = zi/sum(zi*mw)
    
    if(allocated(Bg))     deallocate(Bg)
    if(allocated(thermo)) deallocate(thermo)
    
    allocate(&
         Bg(ns,ng),     &
         thermo(ns,15), &
         stat=astat)
    if(astat .ne. 0) stop

    if(ng .gt. 0) then
       Bg = 1
    end if
    
10  format((1p,5e25.15))
    open(1, file='ceqthermo/thermo.dat', form='formatted')
    read(1,10) thermo   ! thermo data
    close(1)

    ! get equilibrium temperature
    call ceq_sys_init(         &
         ns, ne, 0, 0,         &
         Ein, CS,              &
         BG, thermo,           &
         lu_op, diag, sys,     &
         ierr)
    
    call ceq_state(                       &
         sys, N=zi, p_atm=1.0d0, N_h=zi,  &
         T_h=Ti, HoR_eq=hi,               &
         T_eq=Tad, stats=stats, info=ierr)
    
    ! get initial condition
    call ceq_sys_init(         &
         ns, ne, ncs, ng,      &
         Ein, CS,              &
         BG, thermo,           &
         lu_op, diag, sys,     &
         ierr)    
    
    call ceq_state(                       &
         sys, N=zi, p_atm=1.0d0, N_h=zi,  &
         T_h=Ti, HoR_eq=hi, N_eq=zeq,     &
         T_eq=Teq, stats=stats, info=ierr)
    
  END SUBROUTINE ceqinit
  
  SUBROUTINE ceqcalc(nc,ns,hi,c,T,z)

    !% input/output variables
    integer, intent(in)    :: nc
    integer, intent(in)    :: ns
    real*8,  intent(inout) :: hi
    real*8,  intent(in)    :: c(nc)
    real*8,  intent(out)   :: T
    real*8,  intent(out)   :: z(ns)
    !% local variables
    integer                :: i
    integer                :: ierr
    
    call ceq_state(               &
         sys, c=c, p_atm=1.0d0,   &
         HOR=hi, N_eq=z, T_eq=T,  &
         info=ierr)
    
  END SUBROUTINE ceqcalc

END MODULE localCEQ
