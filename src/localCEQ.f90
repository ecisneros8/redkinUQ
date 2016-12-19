MODULE localCEQ

  USE ceq_system
  IMPLICIT NONE
  type(sys_type), pointer, save :: sys

CONTAINS

  SUBROUTINE ceqinit(ncs,ne,ng,ns,nc,csi,Ti,pi,hi,mw,xi,Ein,Bg,Tad,Teq,zeq)

    ! input/output variables
    integer, intent(in)    :: ncs
    integer, intent(in)    :: ne
    integer, intent(in)    :: ng
    integer, intent(in)    :: ns
    integer, intent(in)    :: nc
    integer, intent(in)    :: csi(ncs)
    real*8,  intent(in)    :: Ti
    real*8,  intent(in)    :: pi
    real*8,  intent(inout) :: hi 
    real*8,  intent(in)    :: mw(ns)
    real*8,  intent(in)    :: xi(ns)
    real*8,  intent(in)    :: Ein(ns,ne)
    real*8,  intent(in)    :: Bg(ns,ng)
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
    real*8,  allocatable   :: Cg(:,:)
    real*8,  allocatable   :: thermo(:,:)
    real*8                 :: stats(20)
    
    lu_op = 0
    diag  = 0
    CS    = csi + 1
    zi    = xi  + 1.0d-15    
    zi    = zi/sum(zi)
    zi    = zi/sum(zi*mw)
    
    if(allocated(Cg))     deallocate(Cg)
    if(allocated(thermo)) deallocate(thermo)
    
    allocate(&
         Cg(ns,ng),     &
         thermo(ns,15), &
         stat=astat)
    if(astat .ne. 0) stop

    Cg = Bg
    
10  format((1p,5e25.15))
    open(1, file='ceqthermo/thermo.dat', form='formatted')
    read(1,10) thermo   ! thermo data
    close(1)

    ! get equilibrium temperature
    call ceq_sys_init(         &
         ns, ne, 0, 0,         &
         Ein, CS,              &
         Cg, thermo,           &
         lu_op, diag, sys,     &
         ierr)
    
    call ceq_state(                       &
         sys, N=zi, p_atm=pi, N_h=zi,     &
         T_h=Ti, HoR_eq=hi,               &
         T_eq=Tad, stats=stats, info=ierr)
    
    ! get initial condition
    call ceq_sys_init(         &
         ns, ne, ncs, ng,      &
         Ein, CS,              &
         Cg, thermo,           &
         lu_op, diag, sys,     &
         ierr)    
    
    call ceq_state(                       &
         sys, N=zi, p_atm=pi, N_h=zi,  &
         T_h=Ti, HoR_eq=hi, N_eq=zeq,     &
         T_eq=Teq, stats=stats, info=ierr)

    if(ierr .eq. -11) then
       Teq = -1.0d0
    end if

  END SUBROUTINE ceqinit
  
  SUBROUTINE ceqRecon(nc,ns,pi,hi,c,T,z,flag)

    !% input/output variables
    integer, intent(in)    :: nc
    integer, intent(in)    :: ns
    real*8,  intent(in)    :: pi
    real*8,  intent(inout) :: hi
    real*8,  intent(inout) :: c(nc)
    real*8,  intent(out)   :: T
    real*8,  intent(out)   :: z(ns)
    integer, intent(out)   :: flag
    !% local variables
    integer                :: i
    integer                :: ierr
    
    call ceq_param_set( sys, T_low = 100.d0 )
    call ceq_state(               &
         sys, c=c, p_atm=pi,      &
         HOR=hi, N_eq=z, T_eq=T,  &
         info=ierr)
    
    if(ierr .lt. 0) then
        flag = 1
    else      
        flag = 0
    end if    

  END SUBROUTINE ceqRecon

END MODULE localCEQ
