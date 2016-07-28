program tpf_driver
  
  use model_nkmp_t, only: model
  use class_TemperedParticleFilter, only: TemperedParticleFilter

  use omp_lib

  implicit none


  type(model) :: m 
  type(TemperedParticleFilter) :: ppf

  logical :: convergence, use_bootstrap
  double precision :: l, truth, rstar

  integer :: rank, nproc, mpierror, i, nthreads, npart, nsim, nintmh

  integer :: time0, time1, rate

  double precision, allocatable :: para(:)
  
  character(len=144) :: output_dir, arg, crstar, cnpart, para_file
  character(len=:), allocatable :: output_file

  use_bootstrap = .false.
  npart = 4000
  nsim = 100
  rstar = 2.0d0
  para_file = 'p0.txt'
  crstar = '2.0'
  cnpart = '40000'
  output_dir = '.'
  nintmh = 1
  do i = 1, command_argument_count()

     call get_command_argument(i, arg)

     select case(arg)

     case ('--output_dir', '-odir')
        call get_command_argument(i+1, output_dir)
     case ('--bootstrap')
        use_bootstrap = .true.
     case ('--parameter-file','-p0')
        call get_command_argument(i+1, para_file)
     case ('--npart', '-m')
        call get_command_argument(i+1, arg)
        cnpart = arg
        read(arg, '(i10)') npart
     case ('--nsim', '-n')
        call get_command_argument(i+1, arg)
        read(arg, '(i10)') nsim
     case ('--rstar', '-r')
        call get_command_argument(i+1, arg)
        crstar = arg
        read(arg, '(f16.8)') rstar
     case ('--nintmh', '-mh')
        call get_command_argument(i+1, arg)
        read(arg, '(i10)') nintmh
     end select
  end do

  if (use_bootstrap) then
     output_file = trim(adjustl(output_dir))//'/bootstrap_npart_'//trim(adjustl(cnpart))//'.txt'
  else
     output_file = trim(adjustl(output_dir))//'/tempered_npart_'//trim(adjustl(cnpart))//'_rstar_'//trim(adjustl(crstar))//'.txt'
  end if
  print*,output_file

  print*,'omp_get_num_procs(nthreads) = ', omp_get_num_procs()
  print*,'omp_get_max_threads(nthreads) = ', omp_get_max_threads()
  call omp_set_num_threads(4)
  print*,'omp_get_num_threads(nthreads) = ', omp_get_num_threads()

  !$OMP PARALLEL DO PRIVATE(i)
  do i = 1, 5
     print*, i, omp_get_thread_num()
  end do
  !$END OMP PARALLEL DO
  m = model()
  rank=0; nproc=1

  ppf = TemperedParticleFilter(m, npart=npart, nproc=nproc, rank=rank, seed=rank, rstar=rstar)

  allocate(para(m%npara))
  open(2, file=para_file, action='read')
  do i = 1, m%npara
     read(2,*) para(i)
     print*, 'para(',i,') = ', para(i)
  end do


  truth = m%lik(para)
  if (rank==0) print*,'exact likelihood:', truth
  call system_clock(count_rate=rate)


  call system_clock(time0)
  open(126,file=output_file,action='write')
  do i = 1, nsim
     l = ppf%lik(para, rank=rank, nproc=nproc)
     if (rank==0) print*, 'loglik = ', l
     write(126,*) l - truth
  end do
  close(126)
  call system_clock(time1)
  print*,'average time = ', ( (time1-time0)/real(rate) / i)
  deallocate(para)


end program
