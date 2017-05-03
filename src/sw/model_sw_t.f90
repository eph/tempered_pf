module model_sw_t
  use, intrinsic :: iso_fortran_env, only: wp => real64


  use gensys, only: do_gensys
  use fortress_bayesian_model_t, only: fortress_lgss_model

  implicit none

  type, public, extends(fortress_lgss_model) :: sw_model 

     integer :: neta  = 12
     real(wp), allocatable  :: GAM0(:, :), GAM1(:, :), C(:), PSI(:, :), PPI(:, :), CC(:)     
   contains
     procedure :: system_matrices
  end type sw_model

  interface sw_model
    module procedure new_sw_model
  end interface


contains

  type(sw_model) function new_sw_model() result(m)

    character(len=144) :: name, datafile
    integer :: nobs, T, ns, npara, neps

    character(len=144) :: prefix
    name = 'sw'

    datafile = './include/tempered_pf/sw/YY.txt'

    nobs = 7
    T = 156
    ns = 53
    npara = 35
    neps = 7
    
    call m%construct_lgss_model_noprior(name, datafile, npara, nobs, T, ns, neps)

    m%p0 = [0.1657d0,0.7869d0,0.5509d0,0.1901d0,1.3333d0,1.6064d0,5.7606d0, &
         0.72d0,0.7d0,1.9d0,0.65d0,0.57d0,0.3d0,0.5462d0,2.0443d0,0.8103d0, &
         0.0882d0,0.2247d0,0.9577d0,0.2194d0,0.9767d0,0.7113d0,0.1479d0, &
         0.8895d0,0.9688d0,0.5d0,0.72d0,0.85d0,0.4582d0,0.24d0,0.5291d0, &
         0.4526d0,0.2449d0,0.141d0,0.2446d0]


   allocate(m%GAM0(m%ns, m%ns), m%GAM1(m%ns, m%ns), m%C(m%ns), m%PSI(m%ns, m%neps), m%PPI(m%ns, m%neta), m%CC(m%ns))

  end function new_sw_model

  subroutine system_matrices(self, para, error) 
    class(sw_model), intent(inout) :: self

    real(wp), intent(in) :: para(self%npara)
    integer, intent(out) :: error

    real(wp) :: constebeta,constepinf,constelab,calfa,csigma,cfc,csadjcost,chabb, &
         cprobw,csigl,cprobp,cindw,cindp,czcap,crpi,crr,cry,crdy,crhoa,crhob, &
         crhog,crhoqs,crhoms,crhopinf,crhow,cgy,cmap,cmaw,siga,sigb,sigg,sigqs,&
         sigm,sigpinf,sigw
    real(wp) :: crk,cpie,cwhlc,cbeta,cg,cky,cr,ctou,clandaw,clandap,cwly,cbetabar,&
         clk,cik,cw,conster,cikbar,ccy,ctrend,crkky,curvp,ciy,cgamma,curvw

    ! gensys 
    real(wp) :: fmat, fwt, ywt, gev, loose, DIV
    integer :: eu(2)


    ! variables 
    integer :: v_r, v_labobs, v_robs, v_pinfobs, v_dy, v_dc, v_dinve, v_dw, v_zcapf, v_rkf, &
         v_kf, v_pkf, v_cf, v_invef, v_yf, v_labf, v_w, v_wf, v_rrf, v_mc, v_zcap, v_rk, &
         v_k, v_pk, v_c, v_inve, v_y, v_lab, v_pinf, v_kpf, v_kp, v_flexgap, v_a, v_b, &
         v_g, v_spinf, v_sw, v_ms, v_qs, v_epinf_VAR, v_ew_VAR, v_Ecf, v_Ew, v_Einvef, &
         v_Ec, v_Erk, v_Epk, v_Elab, v_Epkf, v_Einve, v_Epinf, v_Erkf, v_Elabf 

    ! shocks 
    integer :: e_ea, e_eb, e_eg, e_epinf, e_ew, e_em, e_eqs

    constebeta = para(1);
    constepinf = para(2);
    constelab = para(3);
    calfa = para(4);
    csigma = para(5);
    cfc = para(6);
    csadjcost = para(7);
    chabb = para(8);
    cprobw = para(9);
    csigl = para(10);
    cprobp = para(11);
    cindw = para(12);
    cindp = para(13);
    czcap = para(14);
    crpi = para(15);
    crr = para(16);
    cry = para(17);
    crdy = para(18);
    crhoa = para(19);
    crhob = para(20);
    crhog = para(21);
    crhoqs = para(22);
    crhoms = para(23);
    crhopinf = para(24);
    crhow = para(25);
    cgy = para(26);
    cmap = para(27);
    cmaw = para(28);
    siga = para(29);
    sigb = para(30);
    sigg = para(31);
    sigqs = para(32);
    sigm = para(33);
    sigpinf = para(34);
    sigw = para(35);
    ! Helper parameters
    ctou = 0.025_wp;
    clandaw = 1.5_wp;
    cg = 0.18_wp;
    curvp = 10_wp;
    curvw = 10_wp;
    ctrend = 0.4312_wp;
    cgamma = ctrend/100_wp+1_wp;
    clandap = cfc;
    cbeta = 100_wp/(constebeta+100_wp);
    cpie = constepinf/100_wp+1_wp;
    cbetabar = cbeta*cgamma**(-csigma);
    cr = cpie/(cbeta*cgamma**(-csigma));
    crk = (cbeta**(-1_wp))*(cgamma**csigma) - (1_wp-ctou);
    cw = (calfa**calfa*(1_wp-calfa)**(1_wp-calfa)/(clandap*crk**calfa))**(1_wp/(1_wp-calfa));
    cikbar = (1_wp-(1_wp-ctou)/cgamma);
    cik = (1_wp-(1_wp-ctou)/cgamma)*cgamma;
    clk = ((1_wp-calfa)/calfa)*(crk/cw);
    cky = cfc*(clk)**(calfa-1_wp);
    ciy = cik*cky;
    ccy = 1_wp-cg-cik*cky;
    crkky = crk*cky;
    cwhlc = (1_wp/clandaw)*(1_wp-calfa)/calfa*crk*cky/ccy;
    cwly = 1_wp-crk*cky;
    conster = (cr-1_wp)*100_wp;

    
    v_r = 1
    v_labobs = 2
    v_robs = 3
    v_pinfobs = 4
    v_dy = 5
    v_dc = 6
    v_dinve = 7
    v_dw = 8
    v_zcapf = 9
    v_rkf = 10
    v_kf = 11
    v_pkf = 12
    v_cf = 13
    v_invef = 14
    v_yf = 15
    v_labf = 16
    v_w = 17
    v_wf = 18
    v_rrf = 19
    v_mc = 20
    v_zcap = 21
    v_rk = 22
    v_k = 23
    v_pk = 24
    v_c = 25
    v_inve = 26
    v_y = 27
    v_lab = 28
    v_pinf = 29
    v_kpf = 30
    v_kp = 31
    v_flexgap = 32
    v_a = 33
    v_b = 34
    v_g = 35
    v_spinf = 36
    v_sw = 37
    v_ms = 38
    v_qs = 39
    v_epinf_VAR = 40
    v_ew_VAR = 41
    v_Ecf = 42
    v_Ew = 43
    v_Einvef = 44
    v_Ec = 45
    v_Erk = 46
    v_Epk = 47
    v_Elab = 48
    v_Epkf = 49
    v_Einve = 50
    v_Epinf = 51
    v_Erkf = 52
    v_Elabf = 53
    e_ea = 1
    e_eb = 2
    e_eg = 3
    e_epinf = 4
    e_ew = 5
    e_em = 6
    e_eqs = 7


    self%GAM0 = 0.0_wp
    self%GAM1 = 0.0_wp
    self%PSI = 0.0_wp
    self%PPI = 0.0_wp
    self%C = 0.0_wp

    self%GAM0(1, v_y) = crdy + cry*(-1.0*crr + 1.0)
    self%GAM0(1, v_yf) = -1.0*crdy - 1.0*cry*(-1.0*crr + 1.0)
    self%GAM0(1, v_ms) = 1.00000000000000d0
    self%GAM0(1, v_pinf) = crpi*(-1.0*crr + 1.0)
    self%GAM0(1, v_r) = -1.00000000000000d0
    self%GAM1(1, v_y) = crdy
    self%GAM1(1, v_r) = -1.0*crr
    self%GAM1(1, v_yf) = -1.0*crdy



    self%GAM0(2, v_wf) = -1.0*calfa + 1.0
    self%GAM0(2, v_rkf) = calfa
    self%GAM0(2, v_a) = -1.00000000000000d0



    self%GAM0(3, v_rkf) = 1.0/czcap*(-1.0*czcap + 1.0)
    self%GAM0(3, v_zcapf) = -1.00000000000000d0



    self%GAM0(4, v_wf) = 1.00000000000000d0
    self%GAM0(4, v_labf) = 1.00000000000000d0
    self%GAM0(4, v_rkf) = -1.00000000000000d0
    self%GAM0(4, v_kf) = -1.00000000000000d0



    self%GAM0(5, v_zcapf) = 1.00000000000000d0
    self%GAM0(5, v_kf) = -1.00000000000000d0
    self%GAM1(5, v_kpf) = -1.00000000000000d0



    self%GAM0(6, v_invef) = -1.00000000000000d0
    self%GAM0(6, v_pkf) = cgamma**(-2.0)*1.0/csadjcost*1.0/(cbetabar*cgamma + 1.0)
    self%GAM0(6, v_qs) = 1.00000000000000d0
    self%GAM0(6, v_Einvef) = cbetabar*cgamma*1.0/(cbetabar*cgamma + 1.0)
    self%GAM1(6, v_invef) = -1.0*1.0/(cbetabar*cgamma + 1.0)



    self%GAM0(7, v_Erkf) = crk*1.0/(crk - 1.0*ctou + 1.0)
    self%GAM0(7, v_b) = csigma*1.0/(-1.0*1.0/cgamma*chabb + 1.0)*(1.0/cgamma*chabb + 1.0)
    self%GAM0(7, v_rrf) = -1.00000000000000d0
    self%GAM0(7, v_pkf) = -1.00000000000000d0
    self%GAM0(7, v_Epkf) = (-1.0*ctou + 1.0)*1.0/(crk - 1.0*ctou + 1.0)



    self%GAM0(8, v_Ecf) = 1.0/(1.0/cgamma*chabb + 1.0)
    self%GAM0(8, v_b) = 1.00000000000000d0
    self%GAM0(8, v_rrf) = -1.0*1.0/csigma*(-1.0*1.0/cgamma*chabb + 1.0)*1.0/(1.0/cgamma*chabb + &
         1.0)
    self%GAM0(8, v_Elabf) = -1.0*1.0/csigma*cwhlc*(csigma - 1.0)*1.0/(1.0/cgamma*chabb + 1.0)
    self%GAM0(8, v_cf) = -1.00000000000000d0
    self%GAM0(8, v_labf) = 1.0/csigma*cwhlc*(csigma - 1.0)*1.0/(1.0/cgamma*chabb + 1.0)
    self%GAM1(8, v_cf) = -1.0*1.0/cgamma*chabb*1.0/(1.0/cgamma*chabb + 1.0)



    self%GAM0(9, v_g) = 1.00000000000000d0
    self%GAM0(9, v_invef) = ciy
    self%GAM0(9, v_yf) = -1.00000000000000d0
    self%GAM0(9, v_zcapf) = crkky
    self%GAM0(9, v_cf) = ccy



    self%GAM0(10, v_yf) = -1.00000000000000d0
    self%GAM0(10, v_labf) = cfc*(-1.0*calfa + 1.0)
    self%GAM0(10, v_kf) = calfa*cfc
    self%GAM0(10, v_a) = cfc



    self%GAM0(11, v_wf) = -1.00000000000000d0
    self%GAM0(11, v_cf) = 1.0/(-1.0*1.0/cgamma*chabb + 1.0)
    self%GAM0(11, v_labf) = csigl
    self%GAM1(11, v_cf) = 1.0/cgamma*chabb*1.0/(-1.0*1.0/cgamma*chabb + 1.0)



    self%GAM0(12, v_kpf) = -1.00000000000000d0
    self%GAM0(12, v_invef) = cikbar
    self%GAM0(12, v_qs) = cgamma**2.0*cikbar*csadjcost*(cbetabar*cgamma + 1.0)
    self%GAM1(12, v_kpf) = cikbar - 1.0



    self%GAM0(13, v_mc) = -1.00000000000000d0
    self%GAM0(13, v_rk) = calfa
    self%GAM0(13, v_w) = -1.0*calfa + 1.0
    self%GAM0(13, v_a) = -1.00000000000000d0



    self%GAM0(14, v_zcap) = -1.00000000000000d0
    self%GAM0(14, v_rk) = 1.0/czcap*(-1.0*czcap + 1.0)



    self%GAM0(15, v_w) = 1.00000000000000d0
    self%GAM0(15, v_rk) = -1.00000000000000d0
    self%GAM0(15, v_k) = -1.00000000000000d0
    self%GAM0(15, v_lab) = 1.00000000000000d0



    self%GAM0(16, v_k) = -1.00000000000000d0
    self%GAM0(16, v_zcap) = 1.00000000000000d0
    self%GAM1(16, v_kp) = -1.00000000000000d0



    self%GAM0(17, v_Einve) = cbetabar*cgamma*1.0/(cbetabar*cgamma + 1.0)
    self%GAM0(17, v_pk) = cgamma**(-2.0)*1.0/csadjcost*1.0/(cbetabar*cgamma + 1.0)
    self%GAM0(17, v_inve) = -1.00000000000000d0
    self%GAM0(17, v_qs) = 1.00000000000000d0
    self%GAM1(17, v_inve) = -1.0*1.0/(cbetabar*cgamma + 1.0)



    self%GAM0(18, v_Epinf) = 1.00000000000000d0
    self%GAM0(18, v_Epk) = (-1.0*ctou + 1.0)*1.0/(crk - 1.0*ctou + 1.0)
    self%GAM0(18, v_b) = csigma*1.0/(-1.0*1.0/cgamma*chabb + 1.0)*(1.0/cgamma*chabb + 1.0)
    self%GAM0(18, v_Erk) = crk*1.0/(crk - 1.0*ctou + 1.0)
    self%GAM0(18, v_r) = -1.00000000000000d0
    self%GAM0(18, v_pk) = -1.00000000000000d0



    self%GAM0(19, v_Epinf) = 1.0/csigma*(-1.0*1.0/cgamma*chabb + 1.0)*1.0/(1.0/cgamma*chabb + 1.0)
    self%GAM0(19, v_b) = 1.00000000000000d0
    self%GAM0(19, v_Elab) = -1.0*1.0/csigma*cwhlc*(csigma - 1.0)*1.0/(1.0/cgamma*chabb + 1.0)
    self%GAM0(19, v_Ec) = 1.0/(1.0/cgamma*chabb + 1.0)
    self%GAM0(19, v_c) = -1.00000000000000d0
    self%GAM0(19, v_r) = -1.0*1.0/csigma*(-1.0*1.0/cgamma*chabb + 1.0)*1.0/(1.0/cgamma*chabb + &
         1.0)
    self%GAM0(19, v_lab) = 1.0/csigma*cwhlc*(csigma - 1.0)*1.0/(1.0/cgamma*chabb + 1.0)
    self%GAM1(19, v_c) = -1.0*1.0/cgamma*chabb*1.0/(1.0/cgamma*chabb + 1.0)



    self%GAM0(20, v_g) = 1.00000000000000d0
    self%GAM0(20, v_zcap) = crkky
    self%GAM0(20, v_inve) = ciy
    self%GAM0(20, v_c) = ccy
    self%GAM0(20, v_y) = -1.00000000000000d0



    self%GAM0(21, v_lab) = cfc*(-1.0d0*calfa + 1.0d0)
    self%GAM0(21, v_a) = cfc
    self%GAM0(21, v_k) = calfa*cfc
    self%GAM0(21, v_y) = -1.00000000000000d0



    self%GAM0(22, v_spinf) = 1.000000000000000d0
    self%GAM0(22, v_mc) = 1.0d0/cprobp*(-1.0d0*cprobp + 1.0d0)*1.0d0/(curvp*(cfc - 1.0d0) + 1.0d0)*1.0d0/( &
         cbetabar*cgamma*cindp + 1.0d0)*(-1.0d0*cbetabar*cgamma*cprobp + 1.0d0)
    self%GAM0(22, v_Epinf) = cbetabar*cgamma*1.0d0/(cbetabar*cgamma*cindp + 1.0d0)
    self%GAM0(22, v_pinf) = -1.000000000000000d0
    self%GAM1(22, v_pinf) = -1.0d0*cindp*1.0d0/(cbetabar*cgamma*cindp + 1.0d0)



    self%GAM0(23, v_sw) = 1.000000000000000d0
    self%GAM0(23, v_Ew) = cbetabar*cgamma*1.0d0/(cbetabar*cgamma + 1.0d0)
    self%GAM0(23, v_w) = -1.0d0*1.0d0/cprobw*(-1.0d0*cprobw + 1.0d0)*1.0d0/(cbetabar*cgamma + 1.0d0)*1.0d0/( &
         curvw*(clandaw - 1.0d0) + 1.0d0)*(-1.0d0*cbetabar*cgamma*cprobw + 1.0d0) &
         - 1.0d0
    self%GAM0(23, v_Epinf) = cbetabar*cgamma*1.0d0/(cbetabar*cgamma + 1.0d0)
    self%GAM0(23, v_c) = 1.0d0/cprobw*(-1.0d0*cprobw + 1.0d0)*1.0d0/(cbetabar*cgamma + 1.0d0)*1.0d0/(-1.0d0*1.0d0 &
         /cgamma*chabb + 1.0d0)*1.0d0/(curvw*(clandaw - 1.0d0) + 1.0d0)*(-1.0d0* &
         cbetabar*cgamma*cprobw + 1.0d0)
    self%GAM0(23, v_pinf) = -1.0d0*1.0d0/(cbetabar*cgamma + 1.0d0)*(cbetabar*cgamma*cindw + 1.0d0)
    self%GAM0(23, v_lab) = 1.0d0/cprobw*csigl*(-1.0d0*cprobw + 1.0d0)*1.0d0/(cbetabar*cgamma + 1.0d0)*1.0d0/( &
         curvw*(clandaw - 1.0d0) + 1.0d0)*(-1.0d0*cbetabar*cgamma*cprobw + 1.0d0)
    self%GAM1(23, v_c) = 1.0d0/cgamma*chabb*1.0d0/cprobw*(-1.0d0*cprobw + 1.0d0)*1.0d0/(cbetabar*cgamma + &
         1.0d0)*1.0d0/(-1.0d0*1.0d0/cgamma*chabb + 1.0d0)*1.0d0/(curvw*(clandaw - 1.0d0 &
         ) +  1.0d0)*(-1.0d0*cbetabar*cgamma*cprobw + 1.0d0)
    self%GAM1(23, v_pinf) = -1.0d0*cindw*1.0d0/(cbetabar*cgamma + 1.0d0)
    self%GAM1(23, v_w) = -1.0d0*1.0d0/(cbetabar*cgamma + 1.0d0)



    self%GAM0(24, v_kp) = -1.00000000000000d0
    self%GAM0(24, v_qs) = cgamma**2.0*cikbar*csadjcost*(cbetabar*cgamma + 1.0)
    self%GAM0(24, v_inve) = cikbar
    self%GAM1(24, v_kp) = cikbar - 1.0



    self%GAM0(25, v_flexgap) = -1.00000000000000d0
    self%GAM0(25, v_yf) = -1.00000000000000d0
    self%GAM0(25, v_y) = 1.00000000000000d0



    self%GAM0(26, v_dy) = -1.00000000000000d0
    self%GAM0(26, v_y) = 1.00000000000000d0
    self%GAM1(26, v_y) = 1.00000000000000d0



    self%GAM0(27, v_c) = 1.00000000000000d0
    self%GAM0(27, v_dc) = -1.00000000000000d0
    self%GAM1(27, v_c) = 1.00000000000000d0



    self%GAM0(28, v_dinve) = -1.00000000000000d0
    self%GAM0(28, v_inve) = 1.00000000000000d0
    self%GAM1(28, v_inve) = 1.00000000000000d0



    self%GAM0(29, v_w) = 1.00000000000000d0
    self%GAM0(29, v_dw) = -1.00000000000000d0
    self%GAM1(29, v_w) = 1.00000000000000d0



    self%GAM0(30, v_pinfobs) = -1.00000000000000d0
    self%GAM0(30, v_pinf) = 1.00000000000000d0



    self%GAM0(31, v_r) = 1.00000000000000d0
    self%GAM0(31, v_robs) = -1.00000000000000d0



    self%GAM0(32, v_labobs) = -1.00000000000000d0
    self%GAM0(32, v_lab) = 1.00000000000000d0



    self%GAM0(33, v_a) = -1.00000000000000d0
    self%GAM1(33, v_a) = -1.0d0*crhoa
    self%PSI(33, e_ea) = -1.00000000000000d0



    self%GAM0(34, v_b) = -1.00000000000000d0
    self%GAM1(34, v_b) = -1.0d0*crhob
    self%PSI(34, e_eb) = -1.00000000000000d0



    self%GAM0(35, v_g) = -1.00000000000000d0
    self%GAM1(35, v_g) = -1.0d0*crhog
    self%PSI(35, e_eg) = -1.00000000000000d0
    self%PSI(35, e_ea) = -1.0d0*cgy



    self%GAM0(36, v_spinf) = -1.00000000000000d0
    self%GAM1(36, v_epinf_VAR) = cmap
    self%GAM1(36, v_spinf) = -1.0d0*crhopinf
    self%PSI(36, e_epinf) = -1.00000000000000d0



    self%GAM0(37, v_sw) = -1.00000000000000d0
    self%GAM1(37, v_sw) = -1.0d0*crhow
    self%GAM1(37, v_ew_VAR) = cmaw
    self%PSI(37, e_ew) = -1.00000000000000d0



    self%GAM0(38, v_ms) = -1.00000000000000d0
    self%GAM1(38, v_ms) = -1.0d0*crhoms
    self%PSI(38, e_em) = -1.00000000000000d0



    self%GAM0(39, v_qs) = -1.00000000000000d0
    self%GAM1(39, v_qs) = -1.0d0*crhoqs
    self%PSI(39, e_eqs) = -1.00000000000000d0



    self%GAM0(40, v_epinf_VAR) = -1.00000000000000d0
    self%PSI(40, e_epinf) = -1.00000000000000d0



    self%GAM0(41, v_ew_VAR) = -1.00000000000000d0
    self%PSI(41, e_ew) = -1.00000000000000d0



    self%GAM0(42, v_cf) = 1.0_wp
    self%GAM1(42, v_Ecf) = 1.0_wp
    self%PPI(42, 1) = 1.0_wp


    self%GAM0(43, v_w) = 1.0_wp
    self%GAM1(43, v_Ew) = 1.0_wp
    self%PPI(43, 2) = 1.0_wp


    self%GAM0(44, v_invef) = 1.0_wp
    self%GAM1(44, v_Einvef) = 1.0_wp
    self%PPI(44, 3) = 1.0_wp
      
      
    self%GAM0(45, v_c) = 1.0_wp
    self%GAM1(45, v_Ec) = 1.0_wp
    self%PPI(45, 4) = 1.0_wp
      
      
    self%GAM0(46, v_rk) = 1.0_wp
    self%GAM1(46, v_Erk) = 1.0_wp
    self%PPI(46, 5) = 1.0_wp
      
      
    self%GAM0(47, v_pk) = 1.0_wp
    self%GAM1(47, v_Epk) = 1.0_wp
    self%PPI(47, 6) = 1.0_wp
      
      
    self%GAM0(48, v_lab) = 1.0_wp
    self%GAM1(48, v_Elab) = 1.0_wp
    self%PPI(48, 7) = 1.0_wp
      
      
    self%GAM0(49, v_pkf) = 1.0_wp
    self%GAM1(49, v_Epkf) = 1.0_wp
    self%PPI(49, 8) = 1.0_wp
      
      
    self%GAM0(50, v_inve) = 1.0_wp
    self%GAM1(50, v_Einve) = 1.0_wp
    self%PPI(50, 9) = 1.0_wp


    self%GAM0(51, v_pinf) = 1.0_wp
    self%GAM1(51, v_Epinf) = 1.0_wp
    self%PPI(51, 10) = 1.0_wp

    self%GAM0(52, v_rkf) = 1.0_wp
    self%GAM1(52, v_Erkf) = 1.0_wp
    self%PPI(52, 11) = 1.0_wp


    self%GAM0(53, v_labf) = 1.0_wp
    self%GAM1(53, v_Elabf) = 1.0_wp
    self%PPI(53, 12) = 1.0_wp



    
    call do_gensys(self%TT, self%CC, self%RR, fmat, fwt, ywt, gev, eu, loose, self%GAM0, self%GAM1, self%C, self%PSI, self%PPI, DIV)
! 


    self%QQ = 0.0_wp
    self%ZZ = 0.0_wp
    self%DD = 0.0_wp
    self%HH = 0.0_wp
    self%QQ(e_ea, e_ea) = siga**2
    self%QQ(e_eb, e_eb) = sigb**2
    self%QQ(e_eg, e_eg) = sigg**2
    self%QQ(e_epinf, e_epinf) = sigpinf**2
    self%QQ(e_ew, e_ew) = sigw**2
    self%QQ(e_em, e_em) = sigm**2
    self%QQ(e_eqs, e_eqs) = sigqs**2
    self%ZZ(1, 5) = 1.00000000000000d0
    self%ZZ(2, 6) = 1.00000000000000d0
    self%ZZ(3, 7) = 1.00000000000000d0
    self%ZZ(4, 8) = 1.00000000000000d0
    self%ZZ(5, 2) = 1.00000000000000d0
    self%ZZ(6, 4) = 1.00000000000000d0
    self%ZZ(7, 3) = 1.00000000000000d0
    self%DD(1) = ctrend;
    self%DD(2) = ctrend;
    self%DD(3) = ctrend;
    self%DD(4) = ctrend;
    self%DD(5) = constelab;
    self%DD(6) = constepinf;
    self%DD(7) = 100.0d0*1.0d0/cbeta*cgamma**csigma*cpie - 100.0d0;
    !
    self%HH(1,1) = (0.2_wp*0.865678_wp)**2    
    self%HH(2,2) = (0.2_wp*0.697164_wp)**2    
    self%HH(3,3) = (0.2_wp*2.257345_wp)**2    
    self%HH(4,4) = (0.2_wp*0.563764_wp)**2    
    self%HH(5,5) = (0.2_wp*2.918971_wp)**2    
    self%HH(6,6) = (0.2_wp*0.614915_wp)**2    
    self%HH(7,7) = (0.2_wp*0.826422_wp)**2    

    error = 1
    if (eu(1)*eu(2) == 1) error = 0
    

  end subroutine system_matrices

end module model_sw_t
