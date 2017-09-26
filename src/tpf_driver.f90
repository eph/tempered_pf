program tpf_driver
  use, intrinsic :: iso_fortran_env, only: output_unit, wp => real64

  use flap, only: command_line_interface
  use json_module, only: json_core, json_value

  use fortress_util, only: read_array_from_file
  use fortress_bayesian_model_t, only: fortress_lgss_model

  use class_TemperedParticleFilter, only: TemperedParticleFilter
  use particle_filter, only: part_filter

  use model_nkmp_t, only: nkmp_model
  use model_sw_t, only: sw_model

  use omp_lib

  implicit none

  character(len=400), parameter :: placeholder = '${os.getenv('CONDA_PREFIX', './')}$'

  class(fortress_lgss_model), allocatable :: m 
  type(TemperedParticleFilter) :: tpf

  logical :: convergence, save_states
  double precision :: truth, rstar
  real(wp), allocatable :: pf_liks(:), avg_iterations(:)
  integer :: rank, nproc, i, nthreads, npart, nsim, nintmh, seed, i0

  integer :: time0, time1, rate, filt

  double precision, allocatable :: para(:), incloglh(:), filtered_states(:,:)
  

  character(len=144) :: output_dir, arg, crstar, cnpart
  character(len=144) :: output_file, prefix
  character(len=255) :: para_file, base_dir

  character(len=4) :: model_choice
  character(len=20) :: sample_choice, type
  character(len=:), allocatable :: datafile


  type(command_line_interface) :: cli
  integer :: error
  logical :: found

  type(json_core) :: json
  type(json_value), pointer :: p, inp, op, simi_p, node

  integer :: omp_nthreads, omp_nprocs

  character(len=4) :: simi


  call cli%init(progname='tpf_driver', &
       authors='Ed Herbst', &
       description='Program to illustrate the tempered particle filter.')
  call cli%add(switch='--type',help='Use the bootstrap particle filter instead of TPF', &
       required=.false.,act='store', def='tempered', choices='tempered,bootstrap,resample,opt,aux', error=error)
  call cli%add(switch='--model', switch_ab='-m', help='Model', &
       required=.false.,act='store',def='nkmp', choices='nkmp,sw', error=error)
  call cli%add(switch='--sample', switch_ab='-s', help='Sample', &
       required=.false.,act='store',def='great_moderation', &
       choices='great_moderation,great_recession',error=error)
  call cli%add(switch='--npart', switch_ab='-n', help='Number of particles', &
       required=.false.,act='store',def='4000',error=error)
  call cli%add(switch='--pmsv', switch_ab='-p0', help='Parameter File', &
       required=.false.,act='store', def='0' , error=error)
  call cli%add(switch='--nintmh',switch_ab='-mh', help='Number of intermediate MH steps (for TPF)', &
       required=.false.,act='store',def='1',error=error)
  call cli%add(switch='--rstar', switch_ab='-r', help='Inefficiency Ratio (for TPF)', &
       required=.false.,act='store',def='2.0',error=error)
  call cli%add(switch='--nsim', help='Number of repetitions', &
       required=.false.,act='store',def='100',error=error)
  call cli%add(switch='--seed', help='random seed to use', &
       required=.false.,act='store',def='1848',error=error)
  call cli%add(switch='--output-file', switch_ab='-o', help='Output File', &
       required=.false.,act='store',def='output.json',error=error)
  call cli%add(switch='--save-states', help='Output File', &
       required=.false.,act='store_true',def='.false.',error=error)

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
  call cli%get(switch='--type', val=type, error=error) ; if (error /= 0) stop
  call cli%get(switch='--save-states', val=save_states, error=error) ; if (error /= 0) stop


  base_dir = trim(placeholder)//'/include/tempered_pf/                                          '
  if (model_choice=='nkmp') then
     allocate(m, source=nkmp_model())

     if (sample_choice=='great_recession') then
        m%datafile = 'include/tempered_pf/nkmp/yy2003Q1_2013Q4.txt'
        deallocate(m%yy)
        m%T = 44
        call m%read_data
     end if
  else
     allocate(m, source=sw_model())
  end if

  allocate(para(m%npara))

  if (para_file=='0') then
     para_file = 'include/tempered_pf/'//trim(adjustl(model_choice))//'/p0.txt'
  elseif (para_file=='1') then
     para_file = 'include/tempered_pf/'//trim(adjustl(model_choice))//'/p1.txt'
  end if

  call read_array_from_file(trim(para_file), para)


  call json%create_object(p, '')
  call json%create_object(inp, 'inputs')
  call json%add(p, inp)

  call json%add(inp, 'model', model_choice)
  call json%add(inp, 'sample', sample_choice)
  call json%add(inp, 'nsim', nsim)
  call json%add(inp, 'npart', npart)
  call json%add(inp, 'omp_nthreads', omp_nthreads)
  call json%add(inp, 'omp_nprocs', omp_nprocs)

  call json%add(inp, 'filter', type)
  call json%add(inp, 'rstar', rstar)
  call json%add(inp, 'nintmh', nintmh)

  call json%add(inp, 'para', para)

  call json%create_object(op, 'output')
  call json%add(p, op)

  truth = m%lik(para)
  print*,'exact likelihood:', truth




  tpf = TemperedParticleFilter(m, npart=npart, seed=seed, rstar=rstar)
  tpf%type = type
  tpf%nintmh = nintmh

  allocate(pf_liks(nsim), avg_iterations(nsim))
  call system_clock(count_rate=rate)
  call system_clock(time0)
  do i = 1, nsim

     if (save_states) then
        write(simi, '(I3.3)') i
        call json%create_object(simi_p, simi)
        call json%add(op, simi_p)

        pf_liks(i) = tpf%lik(para, avg_iterations=avg_iterations(i), json_states=simi_p)

        ! only save the whole particle swarm for 1 run
        if (i > 1) then
           call json%get_child(simi_p, 'fcst_states', node)
           call json%remove(node, destroy=.true.)

           call json%get_child(simi_p, 'update_states', node)
           call json%remove(node, destroy=.true.)

           call json%get_child(simi_p, 'tempered_states', node)
           call json%remove(node, destroy=.true.)
           end if
        nullify(simi_p)
     else
        if ((type=='tempered').or.(type=='bootstrap').or.(type=='resample')) then 
           pf_liks(i) = tpf%lik(para, avg_iterations=avg_iterations(i))
        else
           call m%system_matrices(para, error)
           allocate(incloglh(m%T), filtered_states(m%T, m%ns))
           associate(yy => m%yy, TT => m%TT, RR => m%RR, QQ => m%QQ, &
                DD => m%DD, ZZ => m%ZZ, HH => m%HH, ny => m%nobs, &
                ns => m%ns, neps => m%neps, t0 => m%t0)
             if (type=='opt') then
                filt = 1
             else
                filt = 0
             end if
             call part_filter(yy, TT, RR, QQ, DD, ZZ, HH, 0, npart, filt, 2, &
                  ny, m%T, neps, ns, incloglh, filtered_states, seed+i-1, 0)

             pf_liks(i) = sum(incloglh)
             deallocate(incloglh, filtered_states)
           end associate
        end if
     end if
     print*, 'loglik = ', pf_liks(i)
  end do
  call system_clock(time1)
  print*,'average time = ', ( (time1-time0)/real(rate) / i)


  call json%add(op, 'truth', truth)
  call json%add(op, 'average_time', ( (time1-time0)/dble(rate) / i))
  call json%add(op, 'likhat', pf_liks)
  call json%add(op, 'avg_iterations', avg_iterations)

  nullify(inp, op)
  call json%print(p, output_file)

  call json%destroy(p)
  deallocate(para, pf_liks, avg_iterations)


end program
