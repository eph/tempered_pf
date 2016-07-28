module model_nkmp_t

  use gensys, only: do_gensys
  use fortress_bayesian_model_t, only: fortress_lgss_model
  
  implicit none

  type, public, extends(fortress_lgss_model) :: model 

     integer :: neta  = 4
     double precision, allocatable  :: GAM0(:, :), GAM1(:, :), C(:), PSI(:, :), PPI(:, :), CC(:)     
   contains
     procedure :: system_matrices
  end type model

  interface model
    module procedure new_model
  end interface


contains

  type(model) function new_model() result(m)

    character(len=144) :: name, datafile
    integer :: nobs, T, ns, npara, neps

    name = 'nkmp'
    datafile = '/home/eherbst/Dropbox/SMCExtensions/DSGEModelExtensions/code/src/us.txt'
    
    nobs = 3
    T = 80
    ns = 11
    npara = 13
    neps = 3
    
    call m%construct_model(name, datafile, npara, nobs, T, ns, neps)

    m%p0 = [2.0d0,0.5d0,1.5d0,0.25d0,0.7d0,0.5d0,0.3d0,4.0d0,2.0d0,0.5d0,0.1d0,0.1d0,0.1d0]

   allocate(m%GAM0(m%ns, m%ns), m%GAM1(m%ns, m%ns), m%C(m%ns), m%PSI(m%ns, m%neps), m%PPI(m%ns, m%neta), m%CC(m%ns))

  end function new_Model


  subroutine system_matrices(self, para, error) 
    class(model), intent(inout) :: self

    double precision, intent(in) :: para(self%npara)
    integer, intent(out) :: error

    

    double precision :: tau,kap,psi1,psi2,rhor,rhogg,rhozz,rA,piA,gamQ,sigr,sigg,sigz
    double precision :: gss,piss,gam,phi,bet,nu

    ! gensys
    double precision :: fmat, fwt, ywt, gev, loose, DIV
    integer :: eu(2)


    ! variables 
    integer :: v_c, v_ppi, v_R, v_z, v_y, v_g, v_ylag, v_Eppi, v_Ec, v_Ey, v_Ez 

    ! shocks 
    integer :: e_epsg, e_epsz, e_epsr

    tau = para(1);
    kap = para(2);
    psi1 = para(3);
    psi2 = para(4);
    rhor = para(5);
    rhogg = para(6);
    rhozz = para(7);
    rA = para(8);
    piA = para(9);
    gamQ = para(10);
    sigr = para(11);
    sigg = para(12);
    sigz = para(13);

    ! Helper parameters
    nu = 2.0d0;
    bet = 1.0d0/(1.0d0 + rA/400.0d0);
    gam = 1.0d0 + gamQ/100.0d0;
    piss = 1.0d0 + piA/400.0d0;
    phi = 1.0d0;
    gss = 0.2d0;


    v_c = 1
    v_ppi = 2
    v_R = 3
    v_z = 4
    v_y = 5
    v_g = 6
    v_ylag = 7
    v_Eppi = 8
    v_Ec = 9
    v_Ey = 10
    v_Ez = 11
    e_epsg = 1
    e_epsz = 2
    e_epsr = 3



    self%GAM0 = 0.0d0
    self%GAM1 = 0.0d0
    self%PSI = 0.0d0
    self%PPI = 0.0d0
    self%C = 0.0d0

    self%GAM0(1, v_R) = bet
    self%GAM0(1, v_Ec) = -1.0*bet*tau
    self%GAM0(1, v_c) = bet*tau
    self%GAM0(1, v_Eppi) = -1.0*bet
    self%GAM0(1, v_Ez) = -1.0*bet



    self%GAM0(2, v_Eppi) = -1.0*bet
    self%GAM0(2, v_ppi) = 1.00000000000000
    self%GAM0(2, v_c) = -1.0*kap



    self%GAM0(3, v_g) = -1.00000000000000
    self%GAM0(3, v_c) = -1.00000000000000
    self%GAM0(3, v_y) = 1.00000000000000



    self%GAM0(4, v_R) = -1.00000000000000
    self%GAM0(4, v_ppi) = psi1*(-1.0*rhor + 1.0)
    self%GAM0(4, v_g) = -1.0*psi2*(-1.0*rhor + 1.0)
    self%GAM0(4, v_y) = psi2*(-1.0*rhor + 1.0)
    self%GAM1(4, v_R) = -1.0*rhor
    self%PSI(4, e_epsr) = -1.00000000000000



    self%GAM0(5, v_g) = -1.00000000000000
    self%GAM1(5, v_g) = -1.0*rhogg
    self%PSI(5, e_epsg) = -1.00000000000000



    self%GAM0(6, v_z) = -1.00000000000000
    self%GAM1(6, v_z) = -1.0*rhozz
    self%PSI(6, e_epsz) = -1.00000000000000



    self%GAM0(7, v_ylag) = -1.00000000000000
    self%GAM1(7, v_y) = -1.00000000000000



    self%GAM0(8, v_ppi) = 1.0d0
    self%GAM1(8, v_Eppi) = 1.0d0
    self%PPI(8, 1) = 1.0d0


    self%GAM0(9, v_c) = 1.0d0
    self%GAM1(9, v_Ec) = 1.0d0
    self%PPI(9, 2) = 1.0d0


    self%GAM0(10, v_y) = 1.0d0
    self%GAM1(10, v_Ey) = 1.0d0
    self%PPI(10, 3) = 1.0d0


    self%GAM0(11, v_z) = 1.0d0
    self%GAM1(11, v_Ez) = 1.0d0
    self%PPI(11, 4) = 1.0d0




    call do_gensys(self%TT, self%CC, self%RR, fmat, fwt, ywt, gev, eu, loose, self%GAM0, self%GAM1, self%C, self%PSI, self%PPI, DIV)
    !
    if (eu(1)*eu(2) == 1) error = 0

    self%QQ = 0.0d0
    self%ZZ = 0.0d0
    self%DD = 0.0d0
    self%HH = 0.0d0
    self%QQ(e_epsg, e_epsg) = sigg**2
    self%QQ(e_epsz, e_epsz) = sigz**2
    self%QQ(e_epsr, e_epsr) = sigr**2
    self%ZZ(1, 4) = 1.00000000000000d0;
    self%ZZ(1, 5) = 1.00000000000000d0;
    self%ZZ(1, 7) = -1.00000000000000d0;
    self%ZZ(2, 2) = 4.00000000000000d0;
    self%ZZ(3, 3) = 4.00000000000000d0;
    self%DD(1) = gamQ;
    self%DD(2) = piA;
    self%DD(3) = 4.0d0*gamQ + piA + rA;

    self%HH(1, 1) = 0.0134524274371600d0
    self%HH(2, 2) = 0.0865338708889600d0
    self%HH(3, 3) = 0.200334480638760d0

  end subroutine system_matrices

end module model_nkmp_t
