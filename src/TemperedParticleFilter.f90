module class_TemperedParticleFilter
  use, intrinsic :: iso_fortran_env, only: wp => real64
 
  use fortress_bayesian_model_t, only: fortress_ss_model
  use fortress_particles_t, only: fortress_particles
  use fortress_random_t, only: fortress_random
 
  use omp_lib
  implicit none
 
  type TemperedParticleFilter
 
     integer :: npart = 1000
     integer :: nproc = 1
     integer :: nlocalpart, nforeignpart, naltpart
     integer :: nintmh = 1
 
     integer :: initialization_T = 1
 
     integer :: JMAX = 20
 
     double precision :: rstar = 2.0d0
 
     logical :: bootstrap = .false.
     class(fortress_ss_model), allocatable :: m
 
     type(fortress_random) :: rng
 
   contains
 
     !     procedure :: filter
     !     procedure :: smoother
     procedure :: lik
 
  end type TemperedParticleFilter
 
  interface TemperedParticleFilter
     module procedure new_TemperedParticleFilter
  end interface TemperedParticleFilter
 
contains
 
  type(TemperedParticleFilter) function new_TemperedParticleFilter &
       (m, npart, seed, nproc, rank, rstar) result(ppf)
 
    class(fortress_ss_model) :: m
 
    integer, optional, intent(in) :: npart, nproc, seed, rank
    double precision, optional :: rstar
    integer :: rng_seed


    if (present(npart)) ppf%npart = npart
    if (present(npart)) ppf%nproc = nproc
    if (present(rstar)) ppf%rstar = rstar
 
     ppf%nlocalpart = ppf%npart / ppf%nproc
    ppf%nforeignpart = 0!ppf%nlocalpart / 1.5d0
    ppf%naltpart = ppf%nlocalpart / 2.0d0
 
    if (present(seed)) then
       rng_seed = seed
    else
       rng_seed = 0
    end if
 
    ppf%rng = fortress_random(seed=rng_seed)
    allocate(ppf%m, source=m)

 
  end function new_TemperedParticleFilter
 
 
  double precision function lik(ppf, para, rank, nproc, resolve, avg_iterations)
    class(TemperedParticleFilter) :: ppf
 
    double precision, intent(in) :: para(ppf%m%npara)
 
    integer, intent(in) :: rank
    integer, intent(in) :: nproc
 
    logical, optional, intent(in) :: resolve
    real(wp), optional, intent(out) :: avg_iterations
    
    logical :: converged, resolve_opt
 
    type(fortress_particles) :: old_local, old_copy, old_foreign, new_local, alt_foreign, alt_copy
   double precision :: py, py_across_procs, relative_weight, py_vector(nproc)
    double precision :: incwt(ppf%nlocalpart)
 
    double precision :: ess, min_ess, randu(ppf%m%T*ppf%JMAX,1), effective_procs, randu2(ppf%nlocalpart,1)
    double precision :: endogsteady(ppf%m%ns), shocks(ppf%m%neps, ppf%nlocalpart)
 
    double precision :: phi_new, phi_old, phi_max, delta_phi, detHH, old_weight(ppf%nlocalpart), wt_kernel(ppf%nlocalpart)
    double precision :: tol, lb, ub, x1, ess0, ess1, essx, alp, pyn1, pyn0, pr_ratio, x1_trans
    logical :: brent_converged, in_bounds
    integer :: i, j, t, phi_j, brent_loops, ii, jj, acpt(ppf%npart), k, nintmh
 
    double precision :: cdf(ppf%nlocalpart), uu(ppf%nlocalpart)
    double precision :: part_rep(ppf%m%ns, ppf%nlocalpart)
    double precision :: shocks_rep(ppf%m%neps, ppf%nlocalpart), shocks_new(ppf%m%neps,ppf%nlocalpart)
    double precision :: part_old_rep(ppf%m%ns, ppf%nlocalpart), st_new(ppf%m%ns)
   
    double precision :: psi_new, psi_min, psi_old, dessx, scale
 
    integer :: total_stages, explode_count

 
 
    character(len=10) :: charphi
 
    resolve_opt = .true.
 
    if (present(resolve)) resolve_opt = resolve
 
    if (resolve_opt) then
       converged = ppf%m%solve(para, nproc, rank)
       if (converged .eqv. .false.) then
          lik = -1000000000000.0d0
          return
       end if
    end if
 
    lik = 0.0d0
    detHH = 1.0d0
    do j = 1, ppf%m%nobs
       detHH = detHH * ppf%m%HH(j,j)
    end do
   
 
    randu = ppf%rng%uniform_rvs(ppf%m%T*ppf%JMAX, 1)
 
 
    ! initialize part
    old_local = fortress_particles(nvars=ppf%m%ns, npart=ppf%nlocalpart)
    new_local = fortress_particles(nvars=ppf%m%ns, npart=ppf%nlocalpart)
 
    endogsteady = ppf%m%steadystate(para)
   
    !old_local%particles = 0.0d0
    do j = 1, ppf%m%ns
       old_local%particles(j,:) = endogsteady(j)
    end do
 

 
    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)
 
       !new_local%particles(:,1) = ppf%m%policy_function(old_local%particles(:,ppf%nlocalpart), shocks(:,1), [1,1])
       j = 2
       explode_count = 1
       do while (j <= ppf%nlocalpart)
 
          old_local%particles(:,j) = ppf%m%policy_function(old_local%particles(:,j-1), shocks(:,j))
          ! check for explosions
          if ( isnan(sum(old_local%particles(:,j))) ) then
             j = max(j-3, 2)
             shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)
             explode_count = explode_count + 1
 
             if (explode_count > 1) then
                print*,'particle',j,'explodes on rank',rank
 
                old_local%particles(:,j-1) = 0.0d0
                do k = 1, ppf%m%ns
                   old_local%particles(k,j-1) = endogsteady(k)
                end do
 
 
             end if
          else
             j = j + 1
             explode_count = 1
          end if
 
       end do
    end do
 
 
    ! second initialization
    !call mpi_barrier(MPI_COMM_WORLD, mpierror)
    do t = 1, ppf%initialization_T
       shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)
 
       !$OMP PARALLEL DO PRIVATE(j)
       do j = 1, ppf%nlocalpart
          new_local%particles(:,j) = ppf%m%policy_function(old_local%particles(:,j), shocks(:,j))
       end do
       !$OMP END PARALLEL DO
       old_local = new_local
    end do
 
    !call mpi_barrier(MPI_COMM_WORLD, mpierror)
 
    old_local%weights = 1.0d0 / (ppf%nlocalpart )!* nproc)
 
    ! ! the particle filter
    total_stages = 0
    !print*,'starting filter'
    do t = 1, ppf%m%T
 
       converged = .false.
 
       psi_old = 0.00d0
       psi_min = 0.001d0
 
       ! draw shocks and regimes
       shocks = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)
 
       !$OMP PARALLEL DO PRIVATE(j)
       do j = 1, ppf%nlocalpart - ppf%nforeignpart
          new_local%particles(:,j) = ppf%m%policy_function(old_local%particles(:,j), &
               shocks(:,j))
       end do
       !$OMP END PARALLEL DO
 
       ! if (t==24) then
       !    open(123,file='gr_forecast.txt',action='write')
       !    do j = 1,ppf%nlocalpart
       !       write(123,'(100f)') new_local%particles(:,j)
       !    end do
       !    close(123)
       ! end if
 
 
       phi_j = 1
       scale = 0.3d0
       do while (.not. converged)
 
          !calculate unnormalized weights
          !$OMP PARALLEL DO PRIVATE(j)
          do j = 1, ppf%nlocalpart
             incwt(j) = ppf%m%logpdfy_kernel(t, new_local%particles(:,j), old_local%particles(:,j), para)   
             old_weight(j) = old_local%weights(j)
!             new_local%weights(j) = ppf%m%pdfy(t
          end do
          !$OMP END PARALLEL DO
 
          ! ! now do the normalization
 
 
          ! ! turn incremental weights into unnormalized kernal
          wt_kernel = incwt
          ess0 = psi_ess(psi_min, psi_old, wt_kernel, ppf%nlocalpart, ppf%rstar)
          ess1 =  psi_ess(1.0d0, psi_old, wt_kernel, ppf%nlocalpart, ppf%rstar)
          if ((ess1 < 0.0d0) .or. (ppf%bootstrap .eqv. .true.)) then
             psi_new = 1.0d0
             converged = .true.
         elseif (ess0 > 0.0d0) then
             print*,'temperature is too high, aborting'
             stop
          else
             tol = 0.1d0
             psi_new = bisection(psi_min, 1.0d0, tol, psi_old, wt_kernel, ppf%npart, ppf%rstar)
          end if
 
          k = ppf%m%nobs
          !$OMP PARALLEL DO PRIVATE(j)
          do j = 1, ppf%nlocalpart
             if (psi_old > 0.0d0) then
                new_local%weights(j) = old_weight(j)  &
                     * (psi_new/psi_old)**(k/2.0d0) * exp(wt_kernel(j)*(psi_new - psi_old))
             else
                new_local%weights(j) = (2.0d0*3.14159265358979323846d0)**(-k/2.0d0) * 1.0d0/sqrt(detHH) &
                     * psi_new**(k/2.0d0) * exp(psi_new*wt_kernel(j)) * old_weight(j)
 
             end if
          end do
          !$OMP END PARALLEL DO
         
          psi_old = psi_new
          psi_min = psi_old
 
 
          ! check the particles
          call new_local%normalize_weights(py)
          ess = new_local%ess()
 
          lik = lik + log(py)
          !print*,t,psi_new,lik
 
          !------------------------------------------------------------
          ! resample
          !------------------------------------------------------------
          !print*,t,psi_new,ess
          if (ess < ppf%nlocalpart) then
 
             cdf(1) = new_local%weights(1)
             do ii=2,ppf%nlocalpart
                cdf(ii) = cdf(ii-1) + new_local%weights(ii)
             end do
 
             uu = ( randu((ppf%JMAX-1)*t+phi_j,1) -1.0d0 + real( (/ (i, i=1,ppf%nlocalpart) /) ,8) ) / real(ppf%nlocalpart,8)
 
             jj=1
             do ii=1,ppf%nlocalpart
                ! move along the CDF
                do while (uu(ii)>cdf(jj))
                   jj=jj+1
                end do
                ! shuffling
                part_rep(:,ii) = new_local%particles(:,jj)
                part_old_rep(:,ii) = old_local%particles(:,jj)
                shocks_rep(:,ii) = shocks(:,jj)
             end do
 
             new_local%particles = part_rep
             new_local%weights = 1.0d0 / ppf%nlocalpart
             old_local%particles = part_old_rep
             old_local%weights = 1.0d0 / ppf%nlocalpart
             shocks = shocks_rep
 
          end if
 
          !------------------------------------------------------------
          ! mutation
          !------------------------------------------------------------
          if (ppf%bootstrap .eqv. .false.) then
             nintmh = ppf%nintmh
 
             !scale = 0.3d0
             do k = 1, nintmh
                shocks_new = ppf%rng%norm_rvs(ppf%m%neps, ppf%nlocalpart)
                shocks_new = shocks_new * scale + shocks
 
                randu2 = ppf%rng%uniform_rvs(ppf%nlocalpart, 1)
 
                acpt = 0
 
                !$OMP PARALLEL DO PRIVATE(j, pyn0, st_new, pyn1, pr_ratio, alp)
                do j = 1, ppf%nlocalpart
                  
                   pyn0 = ppf%m%logpdfy_kernel(t, new_local%particles(:,j), old_local%particles(:,j), para)
                   st_new = ppf%m%policy_function(old_local%particles(:,j), shocks_new(:,j))
                   pyn1 = ppf%m%logpdfy_kernel(t, st_new, old_local%particles(:,j), para)
 
                   pr_ratio = -0.5d0*(sum(shocks_new(:,j)**2) - sum(shocks(:,j)**2))
                   alp = exp(psi_new * (pyn1 - pyn0) + pr_ratio)
 
                   if (randu2(j,1) < alp) then
                      new_local%particles(:,j) = st_new
                      shocks(:,j) = shocks_new(:,j)
                      acpt(j) = 1
                   end if
 
                end do
                !$OMP END PARALLEL DO
 
                !if (sum(acpt)/(1.0d0*ppf%nlocalpart)>0.40d0) scale = scale*1.1d0
                !if (sum(acpt)/(1.0d0*ppf%nlocalpart)<0.40d0) scale = scale*0.9d0
                acpt = 0.95d0 + 0.10d0 * exp(20.0d0(acpt-0.40d0)) / (1.0d0 + exp(20.0d0(acpt-0.40d0)))
                !print*,scale
             end do
             ! if (t==24) then
 
                 !write(*,'(A,I3,A,F8.6,A,F5.3)'),'t = ', t, ' psi_new = ', psi_new, ' acpt = ', sum(acpt)/(
 
             !    write(charphi,'(f8.6)') psi_new
             !    open(123,file='gr_'//trim(adjustl(charphi))//'.txt',action='write')
             !    do j = 1,ppf%nlocalpart
             !       write(123,'(100f)') new_local%particles(:,j)
             !    end do
            !    close(123)
             ! end if
 
          end if
          ! if ((t==24).and. (ppf%bootstrap==.true.)) then
          !    print*,scale,phi_j
          !write(*,'(A,I3,A,F8.6,A,F5.3)'),'t = ', t, ' psi_new = ', psi_new, ' acpt = ', sum(acpt)/(1.0d0*p
 
          !    write(charphi,'(f8.6)') psi_new
          !    open(123,file='grbs.txt',action='write')
          !    do j = 1,ppf%nlocalpart
          !       write(123,'(100f)') new_local%particles(:,j)
          !    end do
          !    close(123)
          ! end if
 
      
          phi_j = phi_j + 1
 
 
          total_stages = total_stages + 1
       end do
       !print*,t,sum(400*log(new_local%particles(9,:)))/ppf%nlocalpart
       old_local = new_local
 
      end do
      if (present(avg_iterations)) avg_iterations = total_stages*1.0d0 / ppf%m%T
      !print*,'number of loops = ', total_stages*1.0d0 / ppf%m%T
    end function lik
 
 
    function psi_ess(psi, psi_old, incwt, npart, r) result(fx)!, fx, dfx)
 
      integer, intent(in) :: npart
      double precision, intent(in) :: psi, psi_old, incwt(npart), r
     
      double precision :: new_weight(npart), new_weight2(npart)
      double precision :: dfx, sw, sw2, a1, a2, da1, da2, fx
 
      if (psi_old == 0.0d0) then
         new_weight = exp(psi * incwt) !* old_weight
         new_weight2 = exp(2.0d0*psi * incwt) !* old_weight
      else
         new_weight = exp((psi - psi_old)*incwt) !* old_weight
         new_weight2 = exp(2.0d0*(psi - psi_old)*incwt) !* old_weight
      end if
      !new_weight = new_weight / ( 1.0d0/npart * sum(new_weight))
      !psi_ess = npart / (1.0 / npart * sum(new_weight**2))
      a1 = sum(new_weight2) / npart
      a2 = (sum(new_weight) / npart)**2
      fx = (a1/a2) - r
 
 
      !da1 = 2.0d0*sum(incwt*new_weight2) / npart
      !da2 = ( 2.0d0 * sum(new_weight) / npart ) * (sum(incwt*new_weight) / npart)
      !print *, psi, psi_old, fx
      !stop
 
 
      !dfx = (da1*a2 - a1*da2) / a2**2
      !dfx = dfx / (a1/a2)
     
 
 
      !if (isnan(phi_ess)) stop
 
    end function psi_ess
 
 
    double precision function newton(lb, ub, tol, psi_old, incwt, npart, rstar, transform)
 
 
      integer, intent(in) :: npart
      double precision, intent(in) :: lb, ub, tol, psi_old, incwt(npart), rstar
 
      logical, intent(in), optional :: transform
 
      logical :: bisection_converged
      integer :: bisection_loops
 
      double precision :: x1, essx
 
             ! call psi_ess(x1, psi_old, wt_kernel, ppf%nlocalpart, ppf%rstar, essx, dessx)
             ! !dessx = dessx * ( exp(-x1_trans/scale)*(ub-lb)/(scale*(1.0d0+exp(-x1_trans/scale)))**2 )
             ! !print*,dessx, exp(-x1_trans), exp(-x1_trans)*(ub-lb)/(1.0d0-exp(-x1_trans))**2
             
             ! brent_converged = abs(essx) < tol
             ! brent_loops = 1
             ! do while (.not. brent_converged)
             !    !x1_trans = x1_trans - essx / dessx
             !    !x1 = lb + (ub - lb) / (1.0d0 + exp(-x1_trans/scale))
             !    x1 = x1 - essx/dessx
             !    ! if (isnan(x1_trans) .and. (essx / dessx < 0)) then
             !    !    x1 = ub - 0.0001d0
             !    !    x1_trans = scale*log( (x1 - lb) / (ub - x1) )
             !    ! end if
             !    in_bounds = (x1 <= 1.0d0) .and. (x1 >= psi_min)
             !    !print*,in_bounds, x1, essx
             !    do while (.not. in_bounds )
             !       ! print*,brent_loops, x1_trans, x1, t
             !       ! open(69,file='TESTWEIGHTS.txt',action='write')
             !       ! do i = 1, ppf%npart
             !       !    write(69,*) wt_kernel(i)
             !       ! end do
             !       ! print*,'phi_old', phi_old
             !       ! x1 = (ub+phi_old)/2.0d0
             !       ! call psi_ess(x1, psi_old, wt_kernel, ppf%nlocalpart, ppf%rstar, essx, dessx)
             !       ! print*,essx
             !       ! print*,dessx
             !       ! close(69)
             !       !stop
 
             !       ! !x1_trans = x1_trans - 0.2 * essx / dessx
             !       ! x1 = lb + (ub - lb) / (1.0d0 + exp(-x1_trans/scale))
             !       ! call psi_ess(1.0d0, psi_old, wt_kernel, ppf%nlocalpart, ppf%rstar, ess1, dessx)
             !       ! print*,'bad guess', x1, psi_min, ess1
             !       ! dessx = dessx * ( exp(-x1_trans/scale)*(ub-lb)/(scale*(1.0d0+exp(-x1_trans/scale)))**
 
             !    end do
             !    print*, 'guessing ', x1, essx
             !    !print*,essx, dessx
             !    ! if (essx > 0.0) then
             !    !    ub = x1
             !    !    x1 = (x1 + lb) / 2.0d0
             !    ! else
             !    !    lb = x1
             !    !    x1 = (x1 + ub) / 2.0d0
             !    ! endif
             !    call psi_ess(x1, psi_old, wt_kernel, ppf%nlocalpart, ppf%rstar, essx, dessx)
             !    !dessx = dessx * ( exp(-x1_trans/scale)*(ub-lb)/(scale*(1.0d0+exp(-x1_trans/scale)))**2 )
             !    !dessx = dessx * ( exp(-x1_trans)*(ub-lb)/(1.0d0-exp(-x1_trans))**2 )
             !    !print*,brent_loops,x1,essx, lb, ub
             !    brent_converged = abs(essx) < tol
 
    end function newton
    double precision function bisection(lb1, ub1, tol, psi_old, incwt, npart, rstar)
 
      integer, intent(in) :: npart
      double precision, intent(in) :: lb1, ub1, tol, psi_old, incwt(npart), rstar
 
      logical :: bisection_converged
      integer :: bisection_loops
     
      double precision :: x1, essx, lb, ub
 
      lb = lb1
      ub = ub1
      x1 = (lb+ub)/2
      essx =  psi_ess(x1, psi_old, incwt, npart, rstar)
      bisection_converged = abs(essx) < tol
      bisection_loops = 1
      do while (.not. bisection_converged)
 
         if (essx > 0.0) then       
            ub = x1                 
            x1 = (x1 + lb) / 2.0d0  
         else                       
            lb = x1                 
            x1 = (x1 + ub) / 2.0d0  
         endif
 
         essx =  psi_ess(x1, psi_old, incwt, npart, rstar)
 
         bisection_converged = abs(essx) < tol
         bisection_loops = bisection_loops + 1
 
 
 
      end do
      !print*,bisection_loops
      bisection = x1
 
    end function bisection
end module class_TemperedParticleFilter
