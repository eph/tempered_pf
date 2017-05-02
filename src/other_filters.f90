module particle_filter
  use, intrinsic :: iso_fortran_env, only: output_unit, wp => real64

  use fortress_random_t, only: fortress_random
  use fortress_linalg, only: dlyap
  implicit none

  real(wp), parameter :: ONE = 1.0d0, ZERO = 0.0d0, NEG_ONE = -1.0d0, M_PI = 3.141592653589793d0
  real(wp), parameter :: really_small = 1e-10


 contains

   subroutine part_filter(y, TT, RR, QQ, DD, ZZ, HH, t0, npart, filter_type, resampling,ny,nobs,neps,ns, & 
        incloglh,filtered_states,random_seed,resamp_int)


     ! Evaluating the likelihood of LGSS via particle filter
     ! filter_type
     !   0 -- bootstrap
     !   1 -- conditionally optimal
     !   2 -- auxiliary particle filter
     !--------------------------------------------------------------------------------
     integer, intent(in) :: ny, nobs, neps, ns, t0, npart, filter_type, resampling, random_seed, resamp_int
     !f2py threadsafe

     double precision, intent(in) :: y(ny,nobs), TT(ns,ns), RR(ns,neps), QQ(neps,neps)
     double precision, intent(in) :: DD(ny), ZZ(ny,ns), HH(ny,ny)

     double precision, intent(out) :: incloglh(nobs)
     double precision, intent(out) :: filtered_states(nobs,ns)

     double precision :: At(ns), Pt(ns,ns), RQR(ns,ns), Kt(ns,ny), QQRRp(neps,ns)
     double precision :: yhat(ny), nut(ny), Ft(ny,ny), iFt(ny,ny), detFt, TTPt(ns,ns)

     integer :: t, info, i, j, resamp_crit

     double precision :: St(ns, npart), Stold(ns, npart)
     double precision :: wt(npart), wtold(npart), incwt(npart), loglh

     ! BLAS functions
     double precision :: ddot

     ! for singular value decomposition of state variance
     double precision :: U(ns,ns), s(ns), Vt(ns,ns), work(6*ns**2)
     double precision :: Stt(ns), ee(neps), Sto(ns), ZZSt(ny,npart),eta(ny,npart)
     integer :: iwork(8*ns)

     ! for g(s_t|s_{t-1}
     double precision :: iRQR(ns,ns), logdetRQR,RRcQQ(ns,neps)
     double precision :: rnd

     ! for lnpy
     double precision :: iHH(ny,ny), logdetHH, iHHnut(ny)
     double precision :: logdetFt, ZZP0(ny,ns), logdetPt, iPt(ns,ns), cPt(ns,ns), gain(ns)
     double precision :: iPteta(ns), Stbar(ns,npart), Stdiff(ns,npart)
     double precision :: lnp, lng, iFtZZP0(ny,ns)

     



     double precision :: ESS, Zt, uu(npart),cdf(npart)
     integer :: paraind(npart)

     type(fortress_random) :: rng
     ! fro random deviates
     !type (VSL_STREAM_STATE) :: stream
     !integer :: brng, seed, methodu, method, mstorage, time_array(8), errcode
     double precision :: eps(ns,npart), x(npart), x1(npart,1)


     if (resamp_int==0) then
        resamp_crit = npart/2
     else
        resamp_crit = npart+1
     end if

     !------------------------------------------------------------
     ! draw random deviates
     !------------------------------------------------------------
     rng = fortress_random(seed=random_seed)
     eps = rng%norm_rvs(ns,npart)
     x1 = rng%uniform_rvs(npart,1)
     x = x1(:,1)

     !------------------------------------------------------------
     ! get initial variance of states
     !------------------------------------------------------------
     call dgemm('n','t', neps, ns, neps, 1.0d0, QQ, neps, RR, ns, 0.0d0, QQRRp, neps)
     call dgemm('n','n', ns, ns, neps, 1.0d0, RR, ns, QQRRp, neps, 0.0d0, RQR, ns)

     call dlyap(TT, RQR, Pt, ns, info)

     ! Pt = TT*Pt*TT' + RQR
     call dgemm('n','n', ns, ns, ns, 1.0d0, TT, ns, Pt, ns, 0.0d0, TTPt, ns)
     Pt = RQR
     call dgemm('n','t', ns, ns, ns, 1.0d0, TTPt, ns, TT, ns, 1.0d0, Pt, ns)
     !------------------------------------------------------------


     !------------------------------------------------------------
     ! draw initial particles s_0 ~ N(0,P0) L = sqrt(P0)
     !------------------------------------------------------------
     call dgesdd('A', ns, ns, Pt, ns, s, U, ns, Vt, ns, work, 6*ns**2, iwork, info)
     do i = 1,ns
        U(:,i) = U(:,i)*sqrt(s(i))
     end do
     call dgemm('n','n',ns,npart,ns,1.0d0,U,ns,eps(:,:),ns,0.0d0,Stold,ns)

     !------------------------------------------------------------


     !------------------------------------------------------------
     ! set up stuff for g(s_t|s_{t-1})
     !------------------------------------------------------------
     ! get sqrt(RRQQRR') (ASSUMES diagonal QQ)
     do i = 1,neps
        RRcQQ(:,i) = RR(:,i)*sqrt(QQ(i,i))
     end do

     ! get pseudo inverse (RR*QQ*RR')
     logdetRQR = 0.0d0
     call dgesdd('A', ns, ns, RQR, ns, s, U, ns, Vt, ns, work, 6*ns**2, iwork, info)
     do i = 1,ns
        if (s(i) > really_small) then
           U(:,i) = U(:,i)/s(i)
           logdetRQR = log(s(i)) + logdetRQR
        else
           U(:,i) = 0.0d0
        end if
     end do
     call dgemm('t','t',ns,ns,ns,1.0d0,Vt,ns,U,ns,0.0d0,iRQR,ns)

     !------------------------------------------------------------


     !------------------------------------------------------------
     ! set up stuff for conditionally optimal pf
     !------------------------------------------------------------
     if (filter_type == 1) then
        ! SVD destroys (!?!?!) RQR on entry...
        call dgemm('n','t', neps, ns, neps, 1.0d0, QQ, neps, RR, ns, 0.0d0, QQRRp, neps)
        call dgemm('n','n', ns, ns, neps, 1.0d0, RR, ns, QQRRp, neps, 0.0d0, RQR, ns)

        Pt = RQR
        call dcopy(ny*ny, HH, 1, Ft, 1)
        call dsymm('r', 'l', ny, ns, ONE, Pt, ns, ZZ, ny, ZERO, ZZP0, ny)
        call dgemm('n', 't', ny, ny, ns, ONE, ZZP0, ny, ZZ, ny, ONE, Ft, ny)

        ! iFt = inv(Ft)
        call dcopy(ny*ny, Ft, 1, iFt, 1)
        call dpotrf('u', ny, iFt, ny, info)
        call dpotri('u', ny, iFt, ny, info)

        ! det(Ft)
        !detFt = determinant(Ft, ny);
        call determinant(Ft, ny, detFt)
        logdetFt = log(detFt)

        ! Get the conditional variance for s_t|s_{t-1},y_t
        call dsymm('l','u',ny,ns,1.0d0,iFt,ny,ZZP0,ny,0.0d0,iFtZZP0,ny)
        call dgemm('t','n',ns,ns,ny,-1.0d0,ZZP0,ny,iFtZZP0,ny,1.0d0,Pt,ns)

        logdetPt = 0.0d0
        call dgesdd('A', ns, ns, Pt, ns, s, U, ns, Vt, ns, work, 6*ns**2, iwork, info)
        do i = 1,ns
           cPt(:,i) = U(:,i)*sqrt(s(i))
           if (s(i) > really_small) then
              U(:,i) = U(:,i)/s(i)
              logdetPt = log(s(i)) + logdetPt
           else
              U(:,i) = 0.0d0
           end if

        end do
        call dgemm('t','t',ns,ns,ns,1.0d0,Vt,ns,U,ns,0.0d0,iPt,ns)

     end if


     !------------------------------------------------------------
     ! set up stuff for lnpy
     !------------------------------------------------------------
     ! ASSUMES diagonal HH
     iHH = 0.0d0
     logdetHH = 0.0d0
     do i = 1,ny
        iHH(i,i) = 1.0d0 / HH(i,i)
        logdetHH = log(HH(i,i)) + logdetHH
     end do
     !------------------------------------------------------------

     wt = 1.0d0
     incloglh = 0.0d0

     filtered_states = 0.0d0

     do t = 1,nobs
        eps = rng%norm_rvs(ns,npart)
        x1 = rng%uniform_rvs(npart,1)
        x = x1(:,1)

        yhat = y(:,t) - DD

        !------------------------------------------------------------
        ! draw St| St[,y_t]
        !------------------------------------------------------------
        if (filter_type == 0) then
           call dgemm('n','n',ns,npart,neps,1.0d0,RRcQQ,ns,eps(1:neps,:),neps,0.0d0,St,ns)
           call dgemm('n','n',ns,npart,ns,1.0d0,TT,ns,Stold,ns,1.0d0,St,ns)
           call dgemm('n','n',ny,npart,ns,1.0d0,ZZ,ny,St,ns,0.0d0,ZZSt,ny)
        elseif (filter_type == 1) then
           Stold(:,1) = 1.0d0
           call dgemm('n','n',ns,npart,ns,1.0d0,TT,ns,Stold,ns,0.0d0,St,ns)
           call dgemm('n','n',ny,npart,ns,1.0d0,ZZ,ny,St,ns,0.0d0,eta,ny)

           Stbar = St
           do i = 1,npart
              eta(:,i) = yhat - eta(:,i)
           end do
           call dgemm('t','n', ns,npart,ny,1.0d0,iFtZZP0,ny,eta,ny,1.0d0,Stbar,ns)
           St = Stbar
           call dgemm('n','n',ns,npart,ns,1.0d0,cPt,ns,eps,ns,1.0d0,St,ns)
           call dgemm('n','n',ny,npart,ns,1.0d0,ZZ,ny,St,ns,0.0d0,ZZSt,ny)

           Stdiff = St
           call dgemm('n','n',ns,npart,ns,-1.0d0,TT,ns,Stold,ns,1.0d0,Stdiff,ns)
        end if


        !ffdsf$OMP PARALLEL DO PRIVATE(i,ee,Sto,nut,iHHnut,Stt) SHARED(incwt,logdetHH,ny,St,ZZ,TT,RRcQQ,eps,yhat,t,ns,neps,npart,iHH,Stold,ZZSt)
        !------------------------------------------------------------
        ! calculate lnpyi
        !------------------------------------------------------------
        do i = 1,npart
           ! Stt = 0.0d0
           ! Sto = Stold(:,i)

           ! ee = eps(1:neps,i)
           ! call dgemv('n',ns,neps,1.0d0,RRcQQ,ns,ee,1,0.0d0,Stt,1)
           ! call dgemv('n',ns,ns,1.0d0,TT,ns,Sto,1,1.0d0,Stt,1)
           ! Stt = Stt + matmul(TT,Sto)

           ! call dgemv('n',ny,ns,-1.0d0,ZZ,ny,Stt,1,1.0d0,nut,1)
           ! nut = yhat - matmul(ZZ,Stt)
           nut = yhat - ZZSt(:,i)

           call dsymv('u', ny, ONE, iHH, ny, nut, 1, ZERO, iHHnut, 1)
           incwt(i) = -0.5d0*ny*log(2*M_PI) - 0.5d0*logdetHH &
                -0.5d0*ddot(ny,nut,1,iHHnut,1)

           if (filter_type == 1) then

              gain = St(:,i) - Stbar(:,i)

              !call dsymv('u', ns, ONE, iPt, ns, gain, 1, ZERO, iPteta, 1)
              call dgemv('n', ns, ns, ONE, iPt, ns, gain, 1, ZERO, iPteta, 1)

              ! lng
              lng = -0.5d0*neps*log(2*M_PI) - 0.5d0*logdetPt &
                   -0.5d0*ddot(ns,gain,1,iPteta,1)

              gain = Stdiff(:,i)

              !call dsymv('u', ns, ONE, iRQR, ns, gain, 1, ZERO, iPteta, 1)
              call dgemv('n', ns, ns, ONE, iRQR, ns, gain, 1, ZERO, iPteta, 1)

              ! lng
              lnp = -0.5d0*neps*log(2*M_PI) - 0.5d0*logdetRQR &
                   -0.5d0*ddot(ns,gain,1,iPteta,1)
              incwt(i) = incwt(i) + lnp - lng

           end if

        end do


        !------------------------------------------------------------
        ! calculate new weights
        !------------------------------------------------------------
        wt = exp(incwt) * wt
        Zt = sum(wt)/(1.0d0*npart)

        incloglh(t) = log(Zt)

        wt = wt / Zt

        !------------------------------------------------------------
        ! save states
        !------------------------------------------------------------
        do i = 1,npart
           filtered_states(t,:) = filtered_states(t,:) + St(:,i)*wt(i)/(1.0d0*npart)
        end do

        !------------------------------------------------------------
        ! resample, if necessary
        !------------------------------------------------------------
        ESS = npart**2 / sum(wt**2)


        if (ESS < resamp_crit) then
           if (resampling==0) then
              call sys_resampling(npart, wt/sum(wt), x(t), paraind)
              Stold = St(:,paraind)
           elseif (resampling==1) then
              call mult_resampling(npart, wt/sum(wt), x, paraind)
              Stold = St(:,paraind)
           elseif (resampling == 2) then
              !wt = 1.0d0
              wtold = wt/sum(wt)
              cdf(1) = wtold(1)
              do i=2,npart
                 cdf(i) = cdf(i-1) + wtold(i)
              end do

              ! draw a starting point

              rnd = x(1)
              uu = ( rnd -1.0d0 + real( (/ (i, i=1,npart) /) ,8) ) / real(npart,8)
              ! start at the bottom of the CDF
              j=1
              do i=1,npart
                 ! move along the CDF
                 do while (uu(i)>cdf(j))
                    j=j+1
                 end do
                 ! shuffling
                 Stold(:,i) = St(:,j)
              end do
           end if

           wt = 1.0d0
        else
           Stold = St
        end if



     end do

     loglh = sum(incloglh)


   end subroutine part_filter


  subroutine determinant(matrix, r, det)


    ! Computes the determinant of symmetric square matrix, matrix (rank r).
    integer, intent(in) :: r
    double precision, intent(in) :: matrix(r,r)

    double precision, intent(out) :: det

    integer :: info, i, piv(r)
    double precision :: matrix_copy(r, r)

    call dcopy(r*r, matrix, 1, matrix_copy, 1)


    call dpotrf('u', r, matrix_copy, r, info)

    if (info .ne. 0) then
       !write(*,'(a,i4)') 'In determinant(), dgetrf returned error code ', info
       det = -10000.0d0
       return
    end if

    det = 1.0d0

    do i = 1, r

       if (.true.) then !(piv(i) .ne. i) then
          det = det * matrix_copy(i, i) * matrix_copy(i,i)
       else
          det = det * matrix_copy(i, i)
       end if

    end do
  end subroutine determinant


  subroutine mult_resampling(npart, wtsim, randu, paraind)

    integer, intent(in) :: npart

    real(8), intent(in) :: wtsim(npart), randu(npart)
    integer, intent(out) :: paraind(npart)

    real(8) :: cumsum(npart), u
    integer :: i, j


    do i = 1, npart
       cumsum(i) = sum(wtsim(1:i))
    end do

    do i = 1, npart

       u = randu(i)

       j = 1
       do
          if (u < cumsum(j)) exit

          j = j + 1
       end do

       paraind(i) = j

    end do

  end subroutine mult_resampling

  subroutine sys_resampling(npart, wtsim, randu, paraind)

    integer, intent(in) :: npart

    real(8), intent(in) :: wtsim(npart), randu
    integer, intent(out) :: paraind(npart)

    integer :: cweights(npart), m(npart), k, cs

    cweights(1) = floor(npart*wtsim(1) + randu)
    m(1) = cweights(1)

    paraind = 0
    paraind(1:m(1)) = 1

    cs = m(1)
    do k = 2, npart
       cweights(k) = floor(sum(npart*wtsim(1:k))  + randu)

       if (k == npart) then
          cweights(npart) = npart
       end if
       m(k) = cweights(k) - cweights(k-1)
       paraind(cs+1:cs+m(k)) = k

       cs = cs + m(k)


    end do

  end subroutine sys_resampling


end module particle_filter
