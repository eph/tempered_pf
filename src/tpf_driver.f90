program tpf_driver
  use, intrinsic :: iso_fortran_env, only: wp => real64


  use flap, only: command_line_interface
  use json_module, only: json_core, json_value

  use fortress_util, only: read_array_from_file
  use fortress_bayesian_model_t, only: fortress_lgss_model

  use class_TemperedParticleFilter, only: TemperedParticleFilter

  use model_nkmp_t, only: nkmp_model
  use model_sw_t, only: sw_model

  use omp_lib

  implicit none


  class(fortress_lgss_model), allocatable :: m 
  type(TemperedParticleFilter) :: tpf

  logical :: convergence, use_bootstrap
  double precision :: truth, rstar
  real(wp), allocatable :: pf_liks(:), avg_iterations(:)
  integer :: rank, nproc, mpierror, i, nthreads, npart, nsim, nintmh, seed

  integer :: time0, time1, rate

  double precision, allocatable :: para(:)
  
  character(len=144) :: output_dir, arg, crstar, cnpart, para_file

  character(len=144) :: output_file, prefix
  character(len=4) :: model_choice
  character(len=20) :: sample_choice
  character(len=:), allocatable :: datafile


  type(command_line_interface) :: cli
  integer :: error

  type(json_core) :: json
  type(json_value), pointer :: p, inp, op

  integer :: omp_nthreads, omp_nprocs


  use_bootstrap = .false.

  call cli%init(progname='Tempered Particle Filtering Example', &
       authors='Ed Herbst')
  call cli%add(switch='--model', switch_ab='-m', help='Model', &
       required=.false.,act='store',def='nkmp', choices='nkmp,sw', error=error)
  call cli%add(switch='--sample', switch_ab='-s', help='Sample', &
       required=.false.,act='store',def='great_moderation', &
       choices='great_moderation,great_recession',error=error)
  call cli%add(switch='--npart', switch_ab='-n', help='Number of particles', &
       required=.false.,act='store',def='4000',error=error)
  call cli%add(switch='--pmsv', switch_ab='-p0', help='Parameter File', &
       required=.false.,act='store',def='p0.txt',error=error)
  call cli%add(switch='--nintmh',switch_ab='-mh', help='Number of intermediate MH steps (for TPF)', &
       required=.false.,act='store',def='1',error=error)
  call cli%add(switch='--rstar', switch_ab='-r', help='Inefficiency Ration (for TPF)', &
       required=.false.,act='store',def='2.0',error=error)
  call cli%add(switch='--nsim', help='Inefficiency Ration (for TPF)', &
       required=.false.,act='store',def='100',error=error)
  call cli%add(switch='--seed', help='random seed to use', &
       required=.false.,act='store',def='1848',error=error)

  call cli%add(switch='--output-file', switch_ab='-o', help='Output File', &
       required=.false.,act='store',def='output.json',error=error)


  call cli%parse()

  call cli%get(switch='--model',val=model_choice,error=error) ; if (error/=0) stop
  call cli%get(switch='--sample',val=sample_choice,error=error) ; if (error/=0) stop
  call cli%get(switch='--npart',val=npart,error=error) ; if (error/=0) stop
  call cli%get(switch='--pmsv',val=para_file,error=error) ; if (error/=0) stop
  call cli%get(switch='--nintmh',val=nintmh,error=error) ; if (error/=0) stop
  call cli%get(switch='--rstar',val=rstar,error=error) ; if (error /= 0) stop
  call cli%get(switch='--nsim',val=nsim,error=error) ; if (error /= 0) stop
  call cli%get(switch='--output-file', val=output_file, error=error) ; if (error /= 0) stop
  call cli%get(switch='--seed', val=seed, error=error) ; if (error /= 0) stop
  omp_nthreads = omp_get_num_procs()
  omp_nprocs = omp_get_max_threads()

  
  if (model_choice=='nkmp') then
     allocate(m, source=nkmp_model())

     if (sample_choice=='great_recession') then
#:if CONDA_BUILD
        call get_environment_variable('CONDA', prefix)
        m%datafile = trim(adjustl(prefix))//'/include/tempered_pf/nkmp/yy2003Q1_2013Q4.txt'
#:else
        m%datafile = 'src/nkmp/yy2003Q1_2013Q4.txt'
#:endif
        deallocate(m%yy)
        m%T = 44
        call m%read_data
     end if
  else
     allocate(m, source=sw_model())
  end if

  rank=0; nproc=1


  call json%create_object(p, '')
  call json%create_object(inp, 'inputs')
  call json%add(p, inp)
  
  call json%add(inp, 'model', model_choice)
  call json%add(inp, 'sample', sample_choice)
  call json%add(inp, 'nsim', nsim)
  call json%add(inp, 'npart', npart)
  call json%add(inp, 'omp_nthreads', omp_nthreads)
  call json%add(inp, 'omp_nprocs', omp_nprocs)

  if (use_bootstrap) then
     call json%add(inp, 'filter', 'bootstrap')
     call json%add(inp, 'nintmh', 1)
  else
     call json%add(inp, 'filter', 'tpf')
     call json%add(inp, 'rstar', rstar)
     call json%add(inp, 'nintmh', nintmh)
  end if


  tpf = TemperedParticleFilter(m, npart=npart, nproc=nproc, rank=rank, seed=seed, rstar=rstar)
  
  allocate(para(m%npara))
  call read_array_from_file(para_file, para)
  call json%add(inp, 'para', para)

  truth = m%lik(para)
  if (rank==0) print*,'exact likelihood:', truth


  allocate(pf_liks(nsim), avg_iterations(nsim))
  call system_clock(count_rate=rate) ;
  call system_clock(time0)
  do i = 1, nsim
     pf_liks(i) = tpf%lik(para, rank=rank, nproc=nproc, avg_iterations=avg_iterations(i))
     if (rank==0) print*, 'loglik = ', pf_liks(i)
  end do
  call system_clock(time1)
  print*,'average time = ', ( (time1-time0)/real(rate) / i)


  call json%create_object(op, 'output')
  call json%add(p, op)

  call json%add(op, 'truth', truth)
  call json%add(op, 'average_time', ( (time1-time0)/dble(rate) / i))
  call json%add(op, 'likhat', pf_liks)
  call json%add(op, 'avg_iterations', avg_iterations)

  nullify(inp, op)
  call json%print(p, output_file)

  call json%destroy(p)
  deallocate(para, pf_liks, avg_iterations)


end program
