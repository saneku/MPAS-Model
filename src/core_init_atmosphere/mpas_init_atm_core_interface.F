! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module init_atm_core_interface

   use mpas_attlist
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units
   use mpas_log, only : mpas_log_write

   contains


   !***********************************************************************
   !
   !  routine init_atm_setup_core
   !
   !> \brief   Init atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_core(core)

      use mpas_derived_types, only : core_type
      use init_atm_core, only : init_atm_core_init, init_atm_core_run, init_atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => init_atm_core_init
      core % core_run => init_atm_core_run
      core % core_finalize => init_atm_core_finalize
      core % define_packages => init_atm_define_packages
      core % setup_packages => init_atm_setup_packages
      core % setup_decompositions => init_atm_setup_decompositions
      core % setup_clock => init_atm_setup_clock
      core % setup_log => init_atm_setup_log
      core % get_mesh_stream => init_atm_get_mesh_stream
      core % setup_immutable_streams => init_atm_setup_immutable_streams
      core % setup_derived_dimensions => init_atm_setup_derived_dimensions
      core % setup_decomposed_dimensions => init_atm_setup_decomposed_dimensions
      core % setup_block => init_atm_setup_block
      core % setup_namelist => init_atm_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'

#include "core_variables.inc"

   end subroutine init_atm_setup_core


   !***********************************************************************
   !
   !  routine init_atm_setup_domain
   !
   !> \brief   Init atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_domain(domain)

      use mpas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

#include "domain_variables.inc"

   end subroutine init_atm_setup_domain


   !***********************************************************************
   !
   !  function init_atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  init_atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_packages(configs, streamInfo, packages, iocontext) result(ierr)

      use mpas_derived_types, only : mpas_pool_type, mpas_io_context_type, MPAS_streamInfo_type
      use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_package

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (MPAS_streamInfo_type), intent(inout) :: streamInfo
      type (mpas_pool_type), intent(inout) :: packages
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr
      logical :: lexist

      logical, pointer :: initial_conds, sfc_update, lbcs
      logical, pointer :: gwd_stage_in, gwd_stage_out, vertical_stage_in, vertical_stage_out, met_stage_in, met_stage_out
      logical, pointer :: gwd_gsl_stage_out
      logical, pointer :: config_native_gwd_static, config_static_interp, config_vertical_grid, config_met_interp
      logical, pointer :: config_native_gwd_gsl_static
      logical, pointer :: first_guess_field
      logical, pointer :: mp_thompson_aers_in
      integer, pointer :: config_init_case
      logical, pointer :: noahmp, config_noahmp_static


      ierr = init_atm_setup_packages_when(configs, packages)
      if (ierr /= 0) then
         return
      end if

      call mpas_pool_get_config(configs, 'config_init_case', config_init_case)
      call mpas_pool_get_config(configs, 'config_static_interp', config_static_interp)
      call mpas_pool_get_config(configs, 'config_native_gwd_static', config_native_gwd_static)
      call mpas_pool_get_config(configs, 'config_native_gwd_gsl_static', config_native_gwd_gsl_static)
      call mpas_pool_get_config(configs, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_get_config(configs, 'config_met_interp', config_met_interp)

      nullify(initial_conds)
      call mpas_pool_get_package(packages, 'initial_condsActive', initial_conds)

      nullify(sfc_update)
      call mpas_pool_get_package(packages, 'sfc_updateActive', sfc_update)

      nullify(lbcs)
      call mpas_pool_get_package(packages, 'lbcsActive', lbcs)

      nullify(gwd_stage_in)
      call mpas_pool_get_package(packages, 'gwd_stage_inActive', gwd_stage_in)

      nullify(gwd_stage_out)
      call mpas_pool_get_package(packages, 'gwd_stage_outActive', gwd_stage_out)

      nullify(gwd_gsl_stage_out)
      call mpas_pool_get_package(packages, 'gwd_gsl_stage_outActive', gwd_gsl_stage_out)

      nullify(vertical_stage_in)
      call mpas_pool_get_package(packages, 'vertical_stage_inActive', vertical_stage_in)

      nullify(vertical_stage_out)
      call mpas_pool_get_package(packages, 'vertical_stage_outActive', vertical_stage_out)

      nullify(met_stage_in)
      call mpas_pool_get_package(packages, 'met_stage_inActive', met_stage_in)

      nullify(met_stage_out)
      call mpas_pool_get_package(packages, 'met_stage_outActive', met_stage_out)

      nullify(mp_thompson_aers_in)
      call mpas_pool_get_package(packages, 'mp_thompson_aers_inActive', mp_thompson_aers_in)

      if (.not. associated(initial_conds) .or. &
          .not. associated(sfc_update) .or. &
          .not. associated(gwd_stage_in) .or. &
          .not. associated(gwd_stage_out) .or. &
          .not. associated(gwd_gsl_stage_out) .or. &
          .not. associated(vertical_stage_in) .or. &
          .not. associated(vertical_stage_out) .or. &
          .not. associated(met_stage_in) .or. &
          .not. associated(met_stage_out) .or. &
          .not. associated(mp_thompson_aers_in)) then
         call mpas_log_write('********************************************************************************', messageType=MPAS_LOG_ERR)
         call mpas_log_write('* Error while setting up packages for init_atmosphere core.',                      messageType=MPAS_LOG_ERR)
         call mpas_log_write('********************************************************************************', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      if (config_init_case == 8) then
         initial_conds = .false.
         sfc_update = .true.
      else
         initial_conds = .true.
         sfc_update = .false.
      end if

      if (config_init_case == 9) then
         lbcs = .true.
         mp_thompson_aers_in = .false.
         inquire(file="QNWFA_QNIFA_SIGMA_MONTHLY.dat",exist=lexist)
         if(lexist) mp_thompson_aers_in = .true.
      else
         lbcs = .false.
         mp_thompson_aers_in = .false.
      end if

      if (config_init_case == 7) then

         !
         ! The logic here is a little convoluted
         ! For input, we want to read in fields from all earlier stages, except if those earlier stages are being run now
         ! For output, we want to output the fields that were computed in a stage and all of those from earlier stages
         !
         gwd_stage_in = config_native_gwd_static .and. &
                        (.not. config_static_interp)
         gwd_stage_out = config_native_gwd_static
         gwd_gsl_stage_out = config_native_gwd_gsl_static
         vertical_stage_in = config_vertical_grid .and. &
                             (.not. config_native_gwd_static) .and. &
                             (.not. config_static_interp)
         vertical_stage_out = config_vertical_grid
         met_stage_in = config_met_interp .and. &
                        (.not. config_native_gwd_static) .and. &
                        (.not. config_static_interp) .and. &
                        (.not. config_vertical_grid)
         met_stage_out = config_met_interp

         mp_thompson_aers_in = .false.
         inquire(file="QNWFA_QNIFA_SIGMA_MONTHLY.dat",exist=lexist)
         if((lexist .and. met_stage_out) .or. (lexist .and. met_stage_in)) mp_thompson_aers_in = .true.

      else if (config_init_case == 8) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .true.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .false.

      !
      ! When interpolating LBC fields, we need all inputs that would be needed for the interpolation
      ! of ICs, so met_stage_in = .true.
      !
      else if (config_init_case == 9) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .true.
         met_stage_out = .true.

         mp_thompson_aers_in = .false.
         inquire(file="QNWFA_QNIFA_SIGMA_MONTHLY.dat",exist=lexist)
         if((lexist .and. met_stage_out) .or. (lexist .and. met_stage_in)) mp_thompson_aers_in = .true.

         initial_conds = .false.   ! Also, turn off the initial_conds package to avoid writing the IC "output" stream

      else if (config_init_case == 13) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .true.
         met_stage_in = .false.
         met_stage_out = .false.

      else
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .true.
      end if

      !
      ! Package for 3-d first-guess atmospheric and land-surface fields is only active if
      ! we are interpolating real-data ICs or LBCs
      !
      nullify(first_guess_field)
      call mpas_pool_get_package(packages, 'first_guess_fieldActive', first_guess_field)

      first_guess_field = .false.
      if ((config_init_case == 7 .and. config_met_interp) .or. config_init_case == 9) then
         first_guess_field = .true.
      end if

      !
      ! Noah-MP
      !
      nullify(config_noahmp_static)
      call mpas_pool_get_config(configs, 'config_noahmp_static', config_noahmp_static)

      nullify(noahmp)
      call mpas_pool_get_package(packages, 'noahmpActive', noahmp)

      if (associated(config_noahmp_static) .and. associated(noahmp)) then
         noahmp = config_noahmp_static
      else
         call mpas_log_write('********************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('* Error while setting up packages for init_atmosphere core:',                      &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('*    Either the ''noahmp'' package or the ''config_noahmp_static'' namelist',      &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('*    option is not defined.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('********************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

   end function init_atm_setup_packages


   !***********************************************************************
   !
   !  function init_atm_setup_clock
   !
   !> \brief   Simulation clock setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_clock(core_clock, configs) result(ierr)

      use mpas_timekeeping, only : mpas_set_time, mpas_set_timeInterval, mpas_create_clock
      use mpas_derived_types, only : MPAS_Clock_type, MPAS_Time_type, MPAS_Timeinterval_type, mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config
      use mpas_kind_types, only : StrKIND

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      character(len=StrKIND), pointer :: config_start_time, config_stop_time
      integer, pointer :: config_fg_interval

      type (MPAS_Time_type)         :: start_time, stop_time
      type (MPAS_TimeInterval_type) :: dt


      ierr = 0

      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_fg_interval', config_fg_interval)

      call mpas_set_time(start_time, dateTimeString=trim(config_start_time))
      call mpas_set_time(stop_time, dateTimeString=trim(config_stop_time))
      call mpas_set_timeInterval(dt, S=config_fg_interval)
      call mpas_create_clock(core_clock, start_time, dt, stopTime=stop_time)

   end function init_atm_setup_clock


   !***********************************************************************
   !
   !  function init_atm_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_log(logInfo, domain, unitNumbers) result(iErr)!{{{

      use mpas_derived_types, only : mpas_log_type, domain_type
      use mpas_log, only : mpas_log_init, mpas_log_open
      use mpas_framework, only : mpas_framework_report_settings

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo  !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain  !< domain object to provide info for setting up log manager
      integer, dimension(2), intent(in), optional :: unitNumbers !< Fortran unit numbers to use for output and error logs
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, unitNumbers=unitNumbers, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

      call mpas_log_write('')
      call mpas_log_write('MPAS Init-Atmosphere Version '//trim(domain % core % modelVersion))
      call mpas_log_write('')

      call mpas_framework_report_settings(domain)

   end function init_atm_setup_log!}}}


   !***********************************************************************
   !
   !  function init_atm_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function init_atm_get_mesh_stream(configs, streamInfo, stream) result(ierr)

      use mpas_kind_types, only : StrKIND
      use mpas_derived_types, only : mpas_pool_type, MPAS_streamInfo_type
      use mpas_pool_routines, only : mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (MPAS_streamInfo_type), intent(inout) :: streamInfo
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      ierr = 0

      write(stream,'(a)') 'input'

   end function init_atm_get_mesh_stream


   !***********************************************************************
   !
   !  function init_atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_decompositions(decompList) result(ierr)

      use mpas_derived_types, only : mpas_decomp_list, mpas_decomp_function, MPAS_DECOMP_NOERR
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method, &
                              mpas_uniform_decomp 

      implicit none

      type (mpas_decomp_list), pointer :: decompList
      integer :: ierr

      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      end if

   end function init_atm_setup_decompositions


   !***********************************************************************
   !
   !  function init_atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be a
   !>  init atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_block(block) result(ierr)

      use mpas_derived_types, only : block_type

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      ierr = 0

      call init_atm_generate_structs(block, block % structs, block % dimensions, block % packages)
#ifdef MPAS_USE_MUSICA
      ierr = init_atm_extend_scalars_for_chemistry(block)
#endif

   end function init_atm_setup_block

#ifdef MPAS_USE_MUSICA
   !-----------------------------------------------------------------------
   !
   !  function init_atm_extend_scalars_for_chemistry
   !
   !> \brief   Build init-atmosphere chemistry scalar lists from MICM+emissions
   !> \details
   !>  Runtime chemistry species are discovered from config_micm_file,
   !>  config_chem_emission_species, and config_chem_additional_species.
   !>  Their union (without duplicates) defines
   !>  chemistry constituents in state/scalars and lbc_state/lbc_scalars.
   !
   !-----------------------------------------------------------------------
   function init_atm_extend_scalars_for_chemistry(block) result(ierr)

      use mpas_kind_types,  only : StrKIND
      use mpas_derived_types, only : MPAS_LOG_CRIT
      use mpas_log,         only : mpas_log_write

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      type (mpas_pool_type), pointer :: statePool, lbcPool
      type (field3dReal), pointer :: scalarsField
      character(len=StrKIND), pointer :: config_micm_file
      character(len=StrKIND), pointer :: config_chem_emission_species
      character(len=StrKIND), pointer :: config_chem_additional_species
      character(len=StrKIND) :: micm_file
      character(len=StrKIND) :: emission_species
      character(len=StrKIND) :: additional_species
      character(len=64), allocatable :: merged_names(:)
      integer :: n_merged
      integer :: ierr_query
      character(len=:), allocatable :: query_message
      integer :: err_level

      ierr = 0
      micm_file = ''
      emission_species = ''
      additional_species = ''

      nullify(config_micm_file)
      nullify(config_chem_emission_species)
      nullify(config_chem_additional_species)
      err_level = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)
      call mpas_pool_get_config(block % domain % configs, 'config_micm_file', config_micm_file)
      call mpas_pool_get_config(block % domain % configs, 'config_chem_emission_species', config_chem_emission_species)
      call mpas_pool_get_config(block % domain % configs, 'config_chem_additional_species', config_chem_additional_species)
      call mpas_pool_set_error_level(err_level)

      if (associated(config_micm_file)) micm_file = trim(config_micm_file)
      if (associated(config_chem_emission_species)) emission_species = trim(config_chem_emission_species)
      if (associated(config_chem_additional_species)) additional_species = trim(config_chem_additional_species)

      if (len_trim(micm_file) == 0 .and. len_trim(emission_species) == 0 .and. &
          len_trim(additional_species) == 0) return

      call init_atm_query_dynamic_species(micm_file, emission_species, additional_species, &
                                          n_merged, merged_names, ierr_query, query_message)
      if (ierr_query /= 0) then
         call mpas_log_write('[InitChem] Could not build chemistry species list: ' // trim(query_message), &
                             messageType=MPAS_LOG_CRIT)
         ierr = ierr_query
         return
      end if

      nullify(statePool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      if (associated(statePool)) then
         nullify(scalarsField)
         call mpas_pool_get_field(statePool, 'scalars', scalarsField, 1)
         if (associated(scalarsField)) then
            call init_atm_rebuild_var_array(statePool, scalarsField, 'num_scalars', .false., n_merged, merged_names)
         end if
      end if

      nullify(lbcPool)
      call mpas_pool_get_subpool(block % structs, 'lbc_state', lbcPool)
      if (associated(lbcPool)) then
         nullify(scalarsField)
         call mpas_pool_get_field(lbcPool, 'lbc_scalars', scalarsField, 1)
         if (associated(scalarsField)) then
            call init_atm_rebuild_var_array(lbcPool, scalarsField, 'num_lbc_scalars', .true., n_merged, merged_names)
         end if
      end if

      if (allocated(merged_names)) deallocate(merged_names)

      contains

      subroutine init_atm_query_dynamic_species(micm_config, emission_list, additional_list, &
                                                n_species, species, error_code, error_message)

         character(len=*), intent(in) :: micm_config
         character(len=*), intent(in) :: emission_list
         character(len=*), intent(in) :: additional_list
         integer, intent(out) :: n_species
         character(len=64), allocatable, intent(out) :: species(:)
         integer, intent(out) :: error_code
         character(len=:), allocatable, intent(out) :: error_message

         character(len=64), allocatable :: micm_names(:)
         character(len=64), allocatable :: emission_names(:)
         character(len=64), allocatable :: additional_names(:)
         character(len=64), allocatable :: work_names(:)
         character(len=64) :: target_name
         integer :: n_micm, n_emission, n_additional, i

         error_code = 0
         error_message = ''
         n_species = 0

         call init_atm_query_micm_species(micm_config, n_micm, micm_names, error_code, error_message)
         if (error_code /= 0) return

         call init_atm_parse_species_list(emission_list, n_emission, emission_names)
         call init_atm_parse_species_list(additional_list, n_additional, additional_names)

         allocate(work_names(max(1, n_micm + n_emission + n_additional)))
         do i = 1, n_micm
            call init_atm_append_unique(work_names, n_species, micm_names(i))
         end do
         do i = 1, n_emission
            target_name = init_atm_target_tracer_for_emission(emission_names(i), micm_names, n_micm)
            call init_atm_append_unique(work_names, n_species, target_name)
         end do
         do i = 1, n_additional
            call init_atm_append_unique(work_names, n_species, additional_names(i))
         end do

         allocate(species(max(1, n_species)))
         if (n_species > 0) species(1:n_species) = work_names(1:n_species)

         if (allocated(micm_names)) deallocate(micm_names)
         if (allocated(emission_names)) deallocate(emission_names)
         if (allocated(additional_names)) deallocate(additional_names)
         if (allocated(work_names)) deallocate(work_names)

      end subroutine init_atm_query_dynamic_species

      subroutine init_atm_query_micm_species(config_path, n_species, species_names, error_code, error_message)

         use iso_fortran_env, only : real64
         use musica_micm, only : micm_t, RosenbrockStandardOrder
         use musica_state, only : state_t
         use musica_util, only : error_t

         character(len=*), intent(in) :: config_path
         integer, intent(out) :: n_species
         character(len=64), allocatable, intent(out) :: species_names(:)
         integer, intent(out) :: error_code
         character(len=:), allocatable, intent(out) :: error_message

         character(len=*), parameter :: molar_mass_property = "molecular weight [kg mol-1]"
         character(len=*), parameter :: molar_mass_fallback = "__molar mass"

         type(micm_t), pointer :: tmp_micm => null()
         type(state_t), pointer :: tmp_state => null()
         type(error_t) :: error
         type(error_t) :: fallback_error

         character(len=64), allocatable :: local_names(:)
         character(len=:), allocatable :: micm_name
         character(len=64) :: lower_name
         integer :: n_total, i
         real(kind=real64) :: molar_mass_val
         logical :: is_background

         error_code = 0
         error_message = ''
         n_species = 0

         if (len_trim(config_path) == 0) then
            allocate(species_names(0))
            return
         end if

         tmp_micm => micm_t(trim(config_path), RosenbrockStandardOrder, error)
         if (.not. error % is_success()) then
            error_code = 1
            error_message = '[InitChem] Failed to open MICM configuration: ' // trim(config_path)
            allocate(species_names(0))
            goto 100
         end if

         tmp_state => tmp_micm % get_state(1, error)
         if (.not. error % is_success()) then
            error_code = 1
            error_message = '[InitChem] Failed to instantiate MICM state for species discovery'
            allocate(species_names(0))
            goto 100
         end if

         n_total = tmp_state % species_ordering % size()
         allocate(local_names(max(1, n_total)))

         do i = 1, n_total
            micm_name = trim(tmp_state % species_ordering % name(i))
            lower_name = init_atm_to_lower(micm_name)
            is_background = trim(lower_name) == 'o2' .or. trim(lower_name) == 'n2' .or. trim(lower_name) == 'h2o'

            molar_mass_val = tmp_micm % get_species_property_double(micm_name, molar_mass_property, error)
            if (.not. error % is_success()) then
               molar_mass_val = tmp_micm % get_species_property_double(micm_name, molar_mass_fallback, fallback_error)
               error = fallback_error
               if (.not. error % is_success()) cycle
            end if

            if (molar_mass_val <= 0.0_real64 .and. .not. is_background) cycle

            n_species = n_species + 1
            local_names(n_species) = trim(init_atm_to_lower(micm_name))
         end do

         allocate(species_names(max(1, n_species)))
         if (n_species > 0) species_names(1:n_species) = local_names(1:n_species)

         if (allocated(local_names)) deallocate(local_names)

100      continue
         if (associated(tmp_state)) deallocate(tmp_state)
         if (associated(tmp_micm)) deallocate(tmp_micm)

      end subroutine init_atm_query_micm_species

      subroutine init_atm_rebuild_var_array(pool, fld, count_dim_name, is_lbc, n_chem, chem_names)

         use mpas_attlist, only : mpas_add_att

         type(mpas_pool_type), pointer, intent(inout) :: pool
         type(field3dReal), pointer, intent(inout) :: fld
         character(len=*), intent(in) :: count_dim_name
         logical, intent(in) :: is_lbc
         integer, intent(in) :: n_chem
         character(len=64), intent(in) :: chem_names(:)

         character(len=64), allocatable :: merged(:)
         character(len=64), allocatable :: old_names(:)
         type(att_lists_type), dimension(:), pointer :: old_atts
         integer :: n_final
         integer :: i
         integer :: base_count
         integer :: chem_start, chem_end

         if (.not. associated(fld)) return

         allocate(old_names(size(fld % constituentNames)))
         old_names(:) = fld % constituentNames(:)

         if (is_lbc) then
            base_count = 5
         else
            base_count = 6
         end if

         allocate(merged(max(1, base_count + n_chem)))
         n_final = 0
         if (is_lbc) then
            call init_atm_append_unique(merged, n_final, 'lbc_qv')
            call init_atm_append_unique(merged, n_final, 'lbc_qc')
            call init_atm_append_unique(merged, n_final, 'lbc_qr')
            call init_atm_append_unique(merged, n_final, 'lbc_nifa')
            call init_atm_append_unique(merged, n_final, 'lbc_nwfa')
         else
            call init_atm_append_unique(merged, n_final, 'qv')
            call init_atm_append_unique(merged, n_final, 'qc')
            call init_atm_append_unique(merged, n_final, 'qr')
            call init_atm_append_unique(merged, n_final, 'nifa')
            call init_atm_append_unique(merged, n_final, 'nwfa')
            call init_atm_append_unique(merged, n_final, 'tke')
         end if

         do i = 1, n_chem
            if (is_lbc) then
               call init_atm_append_unique(merged, n_final, 'lbc_' // trim(chem_names(i)))
            else
               call init_atm_append_unique(merged, n_final, trim(chem_names(i)))
            end if
         end do

         do i = 1, size(old_names)
            call init_atm_set_dimension(pool, 'index_' // trim(old_names(i)), -1)
         end do

         call init_atm_set_dimension(pool, count_dim_name, n_final)

         call init_atm_set_dimension(pool, 'moist_start', 1)
         call init_atm_set_dimension(pool, 'moist_end', 3)
         if (is_lbc) then
            call init_atm_set_dimension(pool, 'number_start', 4)
            call init_atm_set_dimension(pool, 'number_end', 5)
            call init_atm_set_dimension(pool, 'turbulence_start', -1)
            call init_atm_set_dimension(pool, 'turbulence_end', -1)
         else
            call init_atm_set_dimension(pool, 'number_start', 4)
            call init_atm_set_dimension(pool, 'number_end', 5)
            call init_atm_set_dimension(pool, 'turbulence_start', 6)
            call init_atm_set_dimension(pool, 'turbulence_end', 6)
         end if

         chem_start = base_count + 1
         chem_end = n_final
         if (chem_start > chem_end) then
            chem_start = -1
            chem_end = -1
         end if
         call init_atm_set_dimension(pool, 'chem_start', chem_start)
         call init_atm_set_dimension(pool, 'chem_end', chem_end)

         if (associated(fld % constituentNames)) deallocate(fld % constituentNames)
         allocate(fld % constituentNames(n_final))
         do i = 1, n_final
            fld % constituentNames(i) = trim(merged(i))
            call init_atm_set_dimension(pool, 'index_' // trim(merged(i)), i)
         end do

         if (associated(fld % attLists)) then
            old_atts => fld % attLists
            do i = 1, size(old_atts)
               if (associated(old_atts(i) % attList)) deallocate(old_atts(i) % attList)
            end do
            deallocate(old_atts)
         end if

         allocate(fld % attLists(n_final))
         do i = 1, n_final
            allocate(fld % attLists(i) % attList)
            call init_atm_set_species_attributes(fld % attLists(i) % attList, trim(fld % constituentNames(i)), is_lbc)
         end do

         if (allocated(merged)) deallocate(merged)
         if (allocated(old_names)) deallocate(old_names)

      end subroutine init_atm_rebuild_var_array

      subroutine init_atm_set_species_attributes(attList, name, is_lbc)

         use mpas_attlist, only : att_list_type, mpas_add_att

         type(att_list_type), pointer, intent(inout) :: attList
         character(len=*), intent(in) :: name
         logical, intent(in) :: is_lbc

         character(len=64) :: base_name

         base_name = trim(name)
         if (is_lbc) then
            if (len_trim(base_name) > 4 .and. base_name(1:4) == 'lbc_') base_name = base_name(5:)
         end if

         select case (trim(base_name))
         case ('qv')
            call mpas_add_att(attList, 'units', 'kg kg^{-1}')
            call mpas_add_att(attList, 'long_name', init_atm_lbc_text(is_lbc, 'Water vapor mixing ratio'))
         case ('qc')
            call mpas_add_att(attList, 'units', 'kg kg^{-1}')
            call mpas_add_att(attList, 'long_name', init_atm_lbc_text(is_lbc, 'Cloud water mixing ratio'))
         case ('qr')
            call mpas_add_att(attList, 'units', 'kg kg^{-1}')
            call mpas_add_att(attList, 'long_name', init_atm_lbc_text(is_lbc, 'Rain water mixing ratio'))
         case ('nifa')
            call mpas_add_att(attList, 'units', 'nb kg^{-1}')
            call mpas_add_att(attList, 'long_name', init_atm_lbc_text(is_lbc, 'Gocart ice-friendly aerosol number concentration'))
         case ('nwfa')
            call mpas_add_att(attList, 'units', 'nb kg^{-1}')
            call mpas_add_att(attList, 'long_name', init_atm_lbc_text(is_lbc, 'Gocart water-friendly aerosol number concentration'))
         case ('tke')
            call mpas_add_att(attList, 'units', 'm^2 s^{-2}')
            call mpas_add_att(attList, 'long_name', 'Turbulent kinetic energy for the prognostic tke LES scheme')
         case default
            call mpas_add_att(attList, 'units', 'kg kg^{-1}')
            call mpas_add_att(attList, 'long_name', init_atm_lbc_text(is_lbc, trim(base_name) // ' chemistry mixing ratio'))
         end select

      end subroutine init_atm_set_species_attributes

      function init_atm_lbc_text(is_lbc, base_text) result(out_text)

         logical, intent(in) :: is_lbc
         character(len=*), intent(in) :: base_text
         character(len=:), allocatable :: out_text

         if (is_lbc) then
            out_text = trim(base_text) // ' on lateral boundary cells'
         else
            out_text = trim(base_text)
         end if

      end function init_atm_lbc_text

      subroutine init_atm_set_dimension(pool, dim_name, dim_val)

         type(mpas_pool_type), pointer, intent(inout) :: pool
         character(len=*), intent(in) :: dim_name
         integer, intent(in) :: dim_val
         integer, pointer :: dim_ptr

         nullify(dim_ptr)
         call mpas_pool_get_dimension(pool, trim(dim_name), dim_ptr)
         if (associated(dim_ptr)) then
            dim_ptr = dim_val
         else
            call mpas_pool_add_dimension(pool, trim(dim_name), dim_val)
         end if

      end subroutine init_atm_set_dimension

      subroutine init_atm_parse_species_list(list_text, n_names, names)

         character(len=*), intent(in) :: list_text
         integer, intent(out) :: n_names
         character(len=64), allocatable, intent(out) :: names(:)

         character(len=:), allocatable :: clean
         integer :: i, start_pos

         clean = trim(list_text)
         do i = 1, len(clean)
            if (clean(i:i) == ',' .or. clean(i:i) == ';' .or. clean(i:i) == char(9)) clean(i:i) = ' '
         end do

         allocate(names(max(1, len_trim(clean))))
         n_names = 0
         i = 1
         do while (i <= len_trim(clean))
            do while (i <= len_trim(clean) .and. clean(i:i) == ' ')
               i = i + 1
            end do
            if (i > len_trim(clean)) exit
            start_pos = i
            do while (i <= len_trim(clean) .and. clean(i:i) /= ' ')
               i = i + 1
            end do
            n_names = n_names + 1
            names(n_names) = clean(start_pos:i-1)
            if (len_trim(names(n_names)) >= 4) then
               if (init_atm_to_lower(names(n_names)(1:4)) == 'emi_') names(n_names) = names(n_names)(5:)
            end if
         end do

      end subroutine init_atm_parse_species_list

      subroutine init_atm_append_unique(names, n_names, new_name)

         character(len=64), intent(inout) :: names(:)
         integer, intent(inout) :: n_names
         character(len=*), intent(in) :: new_name
         integer :: i

         if (len_trim(new_name) == 0) return
         do i = 1, n_names
            if (trim(init_atm_to_lower(names(i))) == trim(init_atm_to_lower(new_name))) return
         end do
         n_names = n_names + 1
         names(n_names) = trim(new_name)

      end subroutine init_atm_append_unique

      function init_atm_target_tracer_for_emission(emission_name, micm_names, n_micm) result(target_name)

         character(len=*), intent(in) :: emission_name
         character(len=64), intent(in) :: micm_names(:)
         integer, intent(in) :: n_micm
         character(len=64) :: target_name

         character(len=64) :: lower_emit, lower_tracer
         integer :: i

         lower_emit = trim(init_atm_to_lower(emission_name))
         if (init_atm_is_passive_emission(lower_emit)) then
            target_name = trim(emission_name)
         else
            target_name = trim(lower_emit)
         end if

         do i = 1, n_micm
            lower_tracer = trim(init_atm_to_lower(micm_names(i)))
            if (trim(lower_tracer) == trim(lower_emit)) then
               target_name = trim(micm_names(i))
               return
            end if
         end do

      end function init_atm_target_tracer_for_emission

      pure logical function init_atm_is_passive_emission(spec_name)

         character(len=*), intent(in) :: spec_name
         character(len=64) :: lower_name

         lower_name = trim(init_atm_to_lower(spec_name))
         init_atm_is_passive_emission = trim(lower_name) == 'dust_1' .or. &
                                        trim(lower_name) == 'dust1' .or. &
                                        trim(lower_name) == 'bc' .or. &
                                        trim(lower_name) == 'oc' .or. &
                                        trim(lower_name) == 'pm25' .or. &
                                        trim(lower_name) == 'pm2.5' .or. &
                                        trim(lower_name) == 'pm10'

      end function init_atm_is_passive_emission

      pure function init_atm_upper_name(text) result(out_text)

         character(len=*), intent(in) :: text
         character(len=len(text)) :: out_text
         integer :: i, ia

         do i = 1, len(text)
            ia = iachar(text(i:i))
            if (ia >= iachar('a') .and. ia <= iachar('z')) then
               out_text(i:i) = achar(ia - 32)
            else
               out_text(i:i) = text(i:i)
            end if
         end do

      end function init_atm_upper_name

      pure function init_atm_to_lower(text) result(out_text)

         character(len=*), intent(in) :: text
         character(len=len(text)) :: out_text
         integer :: i, ia

         do i = 1, len(text)
            ia = iachar(text(i:i))
            if (ia >= iachar('A') .and. ia <= iachar('Z')) then
               out_text(i:i) = achar(ia + 32)
            else
               out_text(i:i) = text(i:i)
            end if
         end do

      end function init_atm_to_lower

   end function init_atm_extend_scalars_for_chemistry
#endif


#include "setup_immutable_streams.inc"

#include "block_dimension_routines.inc"

#include "define_packages.inc"

#include "setup_packages.inc"

#include "structs_and_variables.inc"

#include "namelist_call.inc"

#include "namelist_defines.inc"

end module init_atm_core_interface
