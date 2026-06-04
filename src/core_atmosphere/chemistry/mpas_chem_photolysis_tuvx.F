! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_chem_photolysis_tuvx
!
!> \brief MPAS-Chem backend for MUSICA/TUV-x photolysis rates
!> \details
!>  Uses the public MUSICA/TUV-x Fortran API to construct a host-driven
!>  photolysis solver, update MPAS atmospheric columns, run TUV-x, and
!>  return mechanism-specific rates to the MPAS-Chem photolysis driver.
!>  MPAS supplies height, air density, temperature, O3, O2, N2, H2O, and
!>  hydrometeor fields; TUV-x configuration files supply spectral physics.
!>
!>  When config_tuvx_config_file is empty, this module is a no-op and
!>  no photolysis rates are supplied by this path.
!
!-----------------------------------------------------------------------
module mpas_chem_photolysis_tuvx

#ifdef MPAS_USE_MUSICA
    use musica_util, only: dk => musica_dk, error_t, mappings_t
    use musica_tuvx, only: tuvx_t
    use musica_tuvx_grid, only: grid_t
    use musica_tuvx_grid_map, only: grid_map_t
    use musica_tuvx_profile, only: profile_t
    use musica_tuvx_profile_map, only: profile_map_t
    use musica_tuvx_radiator, only: radiator_t
    use musica_tuvx_radiator_map, only: radiator_map_t
#endif

    use mpas_kind_types, only: RKIND

    implicit none

    private

    public :: tuvx_init, tuvx_compute_photolysis, tuvx_is_enabled, tuvx_finalize
    public :: tuvx_query_rate_names
    public :: tuvx_n_rates, tuvx_rate_names

#ifdef MPAS_USE_MUSICA
    type :: tuvx_object_set_t
        type(tuvx_t),     pointer :: solver => null()
        type(grid_t),     pointer :: height_grid => null()
        type(grid_t),     pointer :: cloud_wavelength_grid => null()
        type(profile_t),  pointer :: air_profile => null()
        type(profile_t),  pointer :: temp_profile => null()
        type(profile_t),  pointer :: o3_profile => null()
        type(profile_t),  pointer :: o2_profile => null()
        type(profile_t),  pointer :: n2_profile => null()
        type(profile_t),  pointer :: h2o_profile => null()
        type(radiator_t), pointer :: cloud_radiator => null()
    end type tuvx_object_set_t

    type :: tuvx_column_layout_t
        integer :: n_mpas_layers = 0
        integer :: n_ext_layers = 0
        integer :: n_total_layers = 0
    end type tuvx_column_layout_t

    type :: upper_column_profile_t
        real(dk), allocatable :: z_km_edges(:)
        real(dk), allocatable :: temperature_edges(:)
        real(dk), allocatable :: air_edges(:)
        real(dk), allocatable :: o3_edges(:)
    end type upper_column_profile_t

    type :: tuvx_workspace_t
        integer :: n_wavelength_bins = 0
        integer :: n_photo_rates = 0
        integer :: n_heating_rates = 0
        character(len=64), allocatable :: photo_names(:)
        integer, allocatable :: tuvx_indices(:)
        real(dk), allocatable :: cloud_od(:,:)
        real(dk), allocatable :: cloud_ssa(:,:)
        real(dk), allocatable :: cloud_g(:,:,:)
        real(dk), allocatable :: photo_rates(:,:)
        ! TUV-x returns heating rates through the public solver API. MPAS-Chem
        ! keeps this as solver workspace only; it is not a diagnostic or forcing.
        real(dk), allocatable :: heating_rates(:,:)
    end type tuvx_workspace_t

    type(tuvx_object_set_t), save :: tuvx_objects
    type(tuvx_column_layout_t), save :: column_layout
    type(upper_column_profile_t), save :: extension_profile
    type(tuvx_workspace_t), save :: tuvx_work

    ! Cloud optical property constants (warm cloud, wavelength-independent)
    real(dk), parameter :: CLOUD_SSA_VAL = 0.999999_dk  ! pure scattering
    real(dk), parameter :: CLOUD_G_VAL   = 0.85_dk      ! strong forward scattering

    ! Cloud microphysics constants
    real(dk), parameter :: R_EFF_CLOUD = 10.0e-6_dk   ! cloud droplet r_eff [m]
    real(dk), parameter :: R_EFF_RAIN  = 500.0e-6_dk  ! rain drop r_eff [m]
    real(dk), parameter :: RHO_WATER   = 1000.0_dk    ! water density [kg/m^3]

    ! Physical constants
    real(dk), parameter :: NA    = 6.02214076e23_dk  ! Avogadro [mol^-1]
    real(dk), parameter :: M_AIR = 0.02897_dk        ! dry air molar mass [kg/mol]
    real(dk), parameter :: M_O3  = 0.04800_dk        ! ozone molar mass [kg/mol]
    real(dk), parameter :: M_O2  = 0.031998_dk       ! oxygen molar mass [kg/mol]
    real(dk), parameter :: M_N2  = 0.0280134_dk      ! nitrogen molar mass [kg/mol]
    real(dk), parameter :: M_H2O = 0.01801528_dk     ! water molar mass [kg/mol]
    real(dk), parameter :: VMR_O2 = 0.2095_dk        ! O2 volume mixing ratio
    real(dk), parameter :: VMR_N2 = 0.78084_dk       ! N2 volume mixing ratio
#endif

contains

    !------------------------------------------------------------------------
    !  subroutine tuvx_init
    !
    !> \brief Initialize TUV-x with from-host grids and profiles
    !> \details
    !>  Creates height grid and atmospheric profiles (air, temperature, O3, O2, N2, H2O)
    !>  as from-host objects, registers them, constructs the TUV-x solver,
    !>  then retrieves updatable handles. Caches all configured photolysis
    !>  reaction names and allocates persistent work arrays.
    !------------------------------------------------------------------------
    subroutine tuvx_init(config_file, nVertLevels_in, enable_extension, &
                         extension_file, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
#endif

        character(len=*), intent(in)               :: config_file
        integer, intent(in)                        :: nVertLevels_in
        logical, intent(in)                        :: enable_extension
        character(len=*), intent(in)               :: extension_file
        integer, intent(out)                       :: error_code
        character(len=:), allocatable, intent(out)  :: error_message

#ifdef MPAS_USE_MUSICA
        character(len=256) :: msg

        error_code = 0
        error_message = ''

        if (len_trim(config_file) == 0) then
            call mpas_log_write('[TUV-x] No config file specified; TUV-x disabled.')
            return
        end if

        column_layout%n_mpas_layers = nVertLevels_in
        column_layout%n_ext_layers  = 0

        call mpas_log_write('[TUV-x] Initializing with config: ' // trim(config_file))
        write(msg, '(A,I0)') '[TUV-x] MPAS nVertLevels = ', column_layout%n_mpas_layers
        call mpas_log_write(trim(msg))

        ! Load upper-atmosphere climatology if the extension is enabled
        if (enable_extension .and. len_trim(extension_file) > 0) then
            call read_upper_column_extension(trim(extension_file), error_code, error_message)
            if (error_code /= 0) return
            write(msg, '(A,I0,A,A)') '[TUV-x] Extension enabled: column_layout%n_ext_layers = ', column_layout%n_ext_layers, &
                                     ' from ', trim(extension_file)
            call mpas_log_write(trim(msg))
        else if (enable_extension) then
            error_code = 1
            error_message = '[TUV-x] config_tuvx_top_extension is .true. but config_tuvx_extension_file is empty'
            return
        end if

        column_layout%n_total_layers = column_layout%n_mpas_layers + column_layout%n_ext_layers

        write(msg, '(A,I0)') '[TUV-x] Total column layers (grid sections) = ', column_layout%n_total_layers
        call mpas_log_write(trim(msg))

        call create_host_column_objects(error_code, error_message)
        if (error_code /= 0) return

        call create_tuvx_solver(trim(config_file), error_code, error_message)
        if (error_code /= 0) return

        call cache_tuvx_outputs(error_code, error_message)
        if (error_code /= 0) return

        write(msg, '(A,I0,A,I0)') '[TUV-x] Initialized: n_photo=', tuvx_work%n_photo_rates, &
            ' n_heating=', tuvx_work%n_heating_rates
        call mpas_log_write(trim(msg))
#else
        error_code = 0
        error_message = ''
#endif

    end subroutine tuvx_init


    !------------------------------------------------------------------------
    !  subroutine tuvx_query_rate_names
    !
    !> \brief Discover TUV-x photolysis rate names from a config file
    !> \details
    !>  Initializes a temporary TUV-x solver, copies the discovered photolysis
    !>  names, then releases the solver. This is used during MPAS block setup
    !>  so output streams can request mechanism-specific rate names.
    !------------------------------------------------------------------------
    subroutine tuvx_query_rate_names(config_file, nVertLevels_in, enable_extension, &
                                    extension_file, n_rates, names, error_code, error_message)

        character(len=*), intent(in)                :: config_file
        integer, intent(in)                         :: nVertLevels_in
        logical, intent(in)                         :: enable_extension
        character(len=*), intent(in)                :: extension_file
        integer, intent(out)                        :: n_rates
        character(len=64), allocatable, intent(out) :: names(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)  :: error_message

        error_code = 0
        error_message = ''
        n_rates = 0

#ifdef MPAS_USE_MUSICA
        if (len_trim(config_file) == 0) then
            allocate(names(0))
            return
        end if

        call tuvx_init(config_file, nVertLevels_in, enable_extension, extension_file, &
                       error_code, error_message)
        if (error_code /= 0) then
            allocate(names(0))
            call tuvx_finalize()
            return
        end if

        n_rates = tuvx_n_rates()
        allocate(names(max(1, n_rates)))
        if (n_rates > 0) call tuvx_rate_names(names)

        call tuvx_finalize()
#else
        allocate(names(0))
#endif

    end subroutine tuvx_query_rate_names


    !------------------------------------------------------------------------
    !  subroutine tuvx_compute_photolysis
    !
    !> \brief Compute photolysis rates for a single atmospheric column
    !> \details
    !>  Converts MPAS fields to TUV-x units, updates from-host profiles,
    !>  calls TUV-x run, and extracts every cached photolysis rate
    !>  (edge values averaged to midpoints) for the MPAS slice.
    !------------------------------------------------------------------------
    subroutine tuvx_compute_photolysis(cos_sza, z_edges_m, rho_air, temperature, &
                                        qo3, qo2, qn2, qh2o, qc, qr, nVertLevels_in, &
                                        photo_rates_out, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write
        use mpas_timer, only : mpas_timer_start, mpas_timer_stop
#endif

        real(kind=RKIND), intent(in)  :: cos_sza
        real(kind=RKIND), intent(in)  :: z_edges_m(:)       ! nVertLevels+1 edges [m]
        real(kind=RKIND), intent(in)  :: rho_air(:)          ! nVertLevels midpoints [kg/m^3]
        real(kind=RKIND), intent(in)  :: temperature(:)      ! nVertLevels midpoints [K]
        real(kind=RKIND), intent(in)  :: qo3(:)              ! nVertLevels [kg/kg]
        real(kind=RKIND), intent(in)  :: qo2(:)              ! nVertLevels [kg/kg]
        real(kind=RKIND), intent(in)  :: qn2(:)              ! nVertLevels [kg/kg]
        real(kind=RKIND), intent(in)  :: qh2o(:)             ! nVertLevels [kg/kg]
        real(kind=RKIND), intent(in)  :: qc(:)               ! nVertLevels cloud water [kg/kg]
        real(kind=RKIND), intent(in)  :: qr(:)               ! nVertLevels rain water [kg/kg]
        integer, intent(in)           :: nVertLevels_in
        real(kind=RKIND), intent(out) :: photo_rates_out(:,:) ! (n_rates, nVertLevels)
        integer, intent(out)          :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        real(dk) :: z_km(column_layout%n_total_layers + 1), z_mid_km(column_layout%n_total_layers)
        real(dk) :: air_nd(column_layout%n_total_layers), air_edges(column_layout%n_total_layers + 1)
        real(dk) :: temp_mid(column_layout%n_total_layers), temp_edges(column_layout%n_total_layers + 1)
        real(dk) :: o3_nd(column_layout%n_total_layers), o2_nd(column_layout%n_total_layers), n2_nd(column_layout%n_total_layers), h2o_nd(column_layout%n_total_layers)
        real(dk) :: o3_edges(column_layout%n_total_layers + 1), o2_edges(column_layout%n_total_layers + 1)
        real(dk) :: n2_edges(column_layout%n_total_layers + 1), h2o_edges(column_layout%n_total_layers + 1)
        real(dk) :: air_ld(column_layout%n_total_layers), o3_ld(column_layout%n_total_layers), o2_ld(column_layout%n_total_layers)
        real(dk) :: n2_ld(column_layout%n_total_layers), h2o_ld(column_layout%n_total_layers)
        real(dk) :: dz_cm(column_layout%n_total_layers), dz_m_mpas(column_layout%n_mpas_layers)
        real(dk) :: cloud_tau_1d(column_layout%n_mpas_layers)
        real(dk) :: sza_rad, esd, cos_sza_clamped
        type(error_t) :: err

        error_code = 0
        error_message = ''

        cos_sza_clamped = max(-1.0_dk, min(1.0_dk, real(cos_sza, dk)))

        ! Nighttime: skip TUV-x call, return zero for every cached rate
        if (cos_sza_clamped <= 0.0_dk) then
            photo_rates_out(1:tuvx_work%n_photo_rates, 1:column_layout%n_mpas_layers) = 0.0_RKIND
            return
        end if


        call build_column_profiles(z_edges_m, rho_air, temperature, qo3, qo2, qn2, qh2o, &
                                   z_km, z_mid_km, dz_cm, dz_m_mpas, &
                                   air_nd, temp_mid, o3_nd, o2_nd, n2_nd, h2o_nd, &
                                   air_edges, temp_edges, o3_edges, o2_edges, n2_edges, h2o_edges, &
                                   air_ld, o3_ld, o2_ld, n2_ld, h2o_ld)

        call update_tuvx_column_state(z_km, z_mid_km, &
                                      air_nd, air_edges, air_ld, &
                                      temp_mid, temp_edges, &
                                      o3_nd, o3_edges, o3_ld, &
                                      o2_nd, o2_edges, o2_ld, &
                                      n2_nd, n2_edges, n2_ld, &
                                      h2o_nd, h2o_edges, h2o_ld, &
                                      error_code, error_message)
        if (error_code /= 0) return

        call update_cloud_radiator(qc, qr, rho_air, dz_m_mpas, cloud_tau_1d, error_code, error_message)
        if (error_code /= 0) return

        ! SZA in radians; earth-sun distance ~1 AU
        sza_rad = acos(cos_sza_clamped)
        esd = 1.0_dk

        call mpas_timer_start('chem_tuvx_solver_run')
        call tuvx_objects%solver%run(sza_rad, esd, tuvx_work%photo_rates, tuvx_work%heating_rates, err)
        if (capture_tuvx_error(err, error_code, error_message)) then
            call mpas_timer_stop('chem_tuvx_solver_run')
            return
        end if
        call mpas_timer_stop('chem_tuvx_solver_run')

        call copy_tuvx_rates_to_mpas(photo_rates_out)
#else
        error_code = 0
        error_message = ''
        photo_rates_out = 0.0_RKIND
#endif

    end subroutine tuvx_compute_photolysis


    !------------------------------------------------------------------------
    !  function tuvx_n_rates
    !
    !> \brief Returns the number of photolysis reactions TUV-x reported.
    !> \details Zero when TUV-x is disabled; otherwise equals tuvx_work%n_photo_rates.
    !------------------------------------------------------------------------
    integer function tuvx_n_rates()

#ifdef MPAS_USE_MUSICA
        tuvx_n_rates = tuvx_work%n_photo_rates
#else
        tuvx_n_rates = 0
#endif

    end function tuvx_n_rates


    !------------------------------------------------------------------------
    !  subroutine tuvx_rate_names
    !
    !> \brief Return the cached list of TUV-x photolysis reaction names
    !> \details Caller-allocated character array must be at least
    !>  tuvx_n_rates() elements long. Fills entries 1..tuvx_work%n_photo_rates.
    !------------------------------------------------------------------------
    subroutine tuvx_rate_names(names_out)

        character(len=*), intent(out) :: names_out(:)

#ifdef MPAS_USE_MUSICA
        integer :: i
        do i = 1, tuvx_work%n_photo_rates
            names_out(i) = tuvx_work%photo_names(i)
        end do
#endif

    end subroutine tuvx_rate_names


    !------------------------------------------------------------------------
    !  function tuvx_is_enabled
    !
    !> \brief Returns .true. if TUV-x solver is initialized
    !------------------------------------------------------------------------
    logical function tuvx_is_enabled()

#ifdef MPAS_USE_MUSICA
        tuvx_is_enabled = associated(tuvx_objects%solver)
#else
        tuvx_is_enabled = .false.
#endif

    end function tuvx_is_enabled


    !------------------------------------------------------------------------
    !  subroutine tuvx_finalize
    !
    !> \brief Clean up TUV-x resources
    !------------------------------------------------------------------------
    subroutine tuvx_finalize()

#ifdef MPAS_USE_MUSICA
        use mpas_log, only: mpas_log_write

        if (allocated(tuvx_work%photo_rates))     deallocate(tuvx_work%photo_rates)
        if (allocated(tuvx_work%heating_rates))   deallocate(tuvx_work%heating_rates)
        if (allocated(tuvx_work%cloud_od))        deallocate(tuvx_work%cloud_od)
        if (allocated(tuvx_work%cloud_ssa))       deallocate(tuvx_work%cloud_ssa)
        if (allocated(tuvx_work%cloud_g))         deallocate(tuvx_work%cloud_g)
        if (allocated(extension_profile%z_km_edges))  deallocate(extension_profile%z_km_edges)
        if (allocated(extension_profile%temperature_edges))   deallocate(extension_profile%temperature_edges)
        if (allocated(extension_profile%air_edges)) deallocate(extension_profile%air_edges)
        if (allocated(extension_profile%o3_edges))  deallocate(extension_profile%o3_edges)
        if (allocated(tuvx_work%photo_names))     deallocate(tuvx_work%photo_names)
        if (allocated(tuvx_work%tuvx_indices))    deallocate(tuvx_work%tuvx_indices)

        if (associated(tuvx_objects%solver)) then
            deallocate(tuvx_objects%solver)
            nullify(tuvx_objects%solver)
        end if

        if (associated(tuvx_objects%cloud_wavelength_grid)) then
            deallocate(tuvx_objects%cloud_wavelength_grid)
            nullify(tuvx_objects%cloud_wavelength_grid)
        end if

        nullify(tuvx_objects%height_grid)
        nullify(tuvx_objects%air_profile)
        nullify(tuvx_objects%temp_profile)
        nullify(tuvx_objects%o3_profile)
        nullify(tuvx_objects%o2_profile)
        nullify(tuvx_objects%n2_profile)
        nullify(tuvx_objects%h2o_profile)
        nullify(tuvx_objects%cloud_radiator)

        tuvx_work%n_photo_rates = 0
        tuvx_work%n_heating_rates = 0
        tuvx_work%n_wavelength_bins = 0
        column_layout%n_mpas_layers = 0
        column_layout%n_ext_layers = 0
        column_layout%n_total_layers = 0

        call mpas_log_write('[TUV-x] Finalized.')
#endif

    end subroutine tuvx_finalize

#ifdef MPAS_USE_MUSICA
    !------------------------------------------------------------------------
    !  subroutine create_host_column_objects
    !
    !> \brief Build host-controlled TUV-x grids, profiles, and radiators.
    !------------------------------------------------------------------------
    subroutine create_host_column_objects(error_code, error_message)

        use mpas_log, only: mpas_log_write

        integer, intent(out)                       :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err
        character(len=256) :: msg

        error_code = 0
        error_message = ''

        tuvx_objects%height_grid => grid_t("height", "km", column_layout%n_total_layers, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        tuvx_objects%air_profile  => profile_t("air",         "molecule cm-3", tuvx_objects%height_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%temp_profile => profile_t("temperature", "K",             tuvx_objects%height_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%o3_profile   => profile_t("O3",          "molecule cm-3", tuvx_objects%height_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%o2_profile   => profile_t("O2",          "molecule cm-3", tuvx_objects%height_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%n2_profile   => profile_t("N2",          "molecule cm-3", tuvx_objects%height_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%h2o_profile  => profile_t("H2O",         "molecule cm-3", tuvx_objects%height_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        tuvx_work%n_wavelength_bins = 102
        tuvx_objects%cloud_wavelength_grid => grid_t("wavelength", "nm", tuvx_work%n_wavelength_bins, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        tuvx_objects%cloud_radiator => radiator_t("clouds", tuvx_objects%height_grid, tuvx_objects%cloud_wavelength_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        write(msg, '(A,I0)') '[TUV-x] Cloud radiator created with tuvx_work%n_wavelength_bins = ', tuvx_work%n_wavelength_bins
        call mpas_log_write(trim(msg))

    end subroutine create_host_column_objects

    !------------------------------------------------------------------------
    !  subroutine create_tuvx_solver
    !
    !> \brief Construct the TUV-x solver and retrieve solver-owned handles.
    !------------------------------------------------------------------------
    subroutine create_tuvx_solver(config_file, error_code, error_message)

        use mpas_log, only: mpas_log_write

        character(len=*), intent(in)                :: config_file
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)  :: error_message

        type(error_t) :: err
        type(grid_map_t),     pointer :: grids     => null()
        type(profile_map_t),  pointer :: profiles  => null()
        type(radiator_map_t), pointer :: radiators => null()

        error_code = 0
        error_message = ''

        grids => grid_map_t(err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call grids%add(tuvx_objects%height_grid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        profiles => profile_map_t(err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profiles%add(tuvx_objects%air_profile, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profiles%add(tuvx_objects%temp_profile, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profiles%add(tuvx_objects%o3_profile, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profiles%add(tuvx_objects%o2_profile, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profiles%add(tuvx_objects%n2_profile, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profiles%add(tuvx_objects%h2o_profile, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        radiators => radiator_map_t(err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call radiators%add(tuvx_objects%cloud_radiator, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        tuvx_objects%solver => tuvx_t(trim(config_file), grids, profiles, radiators, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        grids => tuvx_objects%solver%get_grids(err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        profiles => tuvx_objects%solver%get_profiles(err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        radiators => tuvx_objects%solver%get_radiators(err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        tuvx_objects%height_grid  => grids%get("height", "km", err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%air_profile  => profiles%get("air", "molecule cm-3", err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%temp_profile => profiles%get("temperature", "K", err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%o3_profile   => profiles%get("O3", "molecule cm-3", err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%o2_profile   => profiles%get("O2", "molecule cm-3", err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%n2_profile   => profiles%get("N2", "molecule cm-3", err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_objects%h2o_profile  => profiles%get("H2O", "molecule cm-3", err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        tuvx_objects%cloud_radiator => radiators%get("clouds", err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call mpas_log_write('[TUV-x] Cloud radiator registered and retrieved from solver.')

    end subroutine create_tuvx_solver

    !------------------------------------------------------------------------
    !  subroutine cache_tuvx_outputs
    !
    !> \brief Cache photolysis/heating output orderings and work arrays.
    !------------------------------------------------------------------------
    subroutine cache_tuvx_outputs(error_code, error_message)

        use mpas_log, only: mpas_log_write

        integer, intent(out)                       :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err
        type(mappings_t), pointer :: photo_ordering => null()
        type(mappings_t), pointer :: heat_ordering  => null()
        integer :: i
        character(len=256) :: msg

        error_code = 0
        error_message = ''

        allocate(tuvx_work%cloud_od(column_layout%n_total_layers, tuvx_work%n_wavelength_bins))
        allocate(tuvx_work%cloud_ssa(column_layout%n_total_layers, tuvx_work%n_wavelength_bins))
        allocate(tuvx_work%cloud_g(column_layout%n_total_layers, tuvx_work%n_wavelength_bins, 1))
        tuvx_work%cloud_od = 0.0_dk
        tuvx_work%cloud_ssa = CLOUD_SSA_VAL
        tuvx_work%cloud_g(:,:,1) = CLOUD_G_VAL

        photo_ordering => tuvx_objects%solver%get_photolysis_rate_constants_ordering(err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_work%n_photo_rates = photo_ordering%size()
        if (tuvx_work%n_photo_rates < 1) then
            error_code = 1
            error_message = '[TUV-x] TUV-x config registered zero photolysis reactions'
            deallocate(photo_ordering)
            return
        end if

        if (allocated(tuvx_work%photo_names))  deallocate(tuvx_work%photo_names)
        if (allocated(tuvx_work%tuvx_indices)) deallocate(tuvx_work%tuvx_indices)
        allocate(tuvx_work%photo_names(tuvx_work%n_photo_rates))
        allocate(tuvx_work%tuvx_indices(tuvx_work%n_photo_rates))
        do i = 1, tuvx_work%n_photo_rates
            tuvx_work%photo_names(i)  = trim(photo_ordering%name(i))
            tuvx_work%tuvx_indices(i) = photo_ordering%index(i)
            write(msg, '(A,A,A,I0)') '[TUV-x] Photolysis reaction ', &
                trim(tuvx_work%photo_names(i)), ' index = ', tuvx_work%tuvx_indices(i)
            call mpas_log_write(trim(msg))
        end do
        deallocate(photo_ordering)

        heat_ordering => tuvx_objects%solver%get_heating_rates_ordering(err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        tuvx_work%n_heating_rates = heat_ordering%size()
        deallocate(heat_ordering)

        allocate(tuvx_work%photo_rates(column_layout%n_total_layers + 1, max(tuvx_work%n_photo_rates, 1)))
        allocate(tuvx_work%heating_rates(column_layout%n_total_layers + 1, max(tuvx_work%n_heating_rates, 1)))
        tuvx_work%photo_rates = 0.0_dk
        tuvx_work%heating_rates = 0.0_dk

    end subroutine cache_tuvx_outputs

    !------------------------------------------------------------------------
    !  subroutine build_column_profiles
    !
    !> \brief Convert one MPAS column into TUV-x grid/profile arrays.
    !------------------------------------------------------------------------
    subroutine build_column_profiles(z_edges_m, rho_air, temperature, qo3, qo2, qn2, qh2o, &
                                     z_km, z_mid_km, dz_cm, dz_m_mpas, &
                                     air_nd, temp_mid, o3_nd, o2_nd, n2_nd, h2o_nd, &
                                     air_edges, temp_edges, o3_edges, o2_edges, n2_edges, h2o_edges, &
                                     air_ld, o3_ld, o2_ld, n2_ld, h2o_ld)

        real(kind=RKIND), intent(in) :: z_edges_m(:)
        real(kind=RKIND), intent(in) :: rho_air(:), temperature(:)
        real(kind=RKIND), intent(in) :: qo3(:), qo2(:), qn2(:), qh2o(:)
        real(dk), intent(out) :: z_km(:), z_mid_km(:), dz_cm(:), dz_m_mpas(:)
        real(dk), intent(out) :: air_nd(:), temp_mid(:), o3_nd(:), o2_nd(:), n2_nd(:), h2o_nd(:)
        real(dk), intent(out) :: air_edges(:), temp_edges(:), o3_edges(:), o2_edges(:), n2_edges(:), h2o_edges(:)
        real(dk), intent(out) :: air_ld(:), o3_ld(:), o2_ld(:), n2_ld(:), h2o_ld(:)

        integer :: k, j

        z_km(1:column_layout%n_mpas_layers + 1) = real(z_edges_m(1:column_layout%n_mpas_layers + 1), dk) * 1.0e-3_dk
        do j = 1, column_layout%n_ext_layers
            z_km(column_layout%n_mpas_layers + 1 + j) = extension_profile%z_km_edges(j + 1)
        end do

        z_mid_km = 0.5_dk * (z_km(1:column_layout%n_total_layers) + z_km(2:column_layout%n_total_layers + 1))
        dz_cm = (z_km(2:column_layout%n_total_layers + 1) - z_km(1:column_layout%n_total_layers)) * 1.0e5_dk
        dz_m_mpas = (z_km(2:column_layout%n_mpas_layers + 1) - z_km(1:column_layout%n_mpas_layers)) * 1.0e3_dk

        air_nd(1:column_layout%n_mpas_layers) = real(rho_air(1:column_layout%n_mpas_layers), dk) * (NA / M_AIR) * 1.0e-6_dk
        temp_mid(1:column_layout%n_mpas_layers) = real(temperature(1:column_layout%n_mpas_layers), dk)
        o3_nd(1:column_layout%n_mpas_layers) = real(qo3(1:column_layout%n_mpas_layers), dk) * &
                                 real(rho_air(1:column_layout%n_mpas_layers), dk) * (NA / M_O3) * 1.0e-6_dk

        if (column_layout%n_ext_layers > 0) then
            block
                real(dk) :: T_anchor_offset, ext_first_mid_T
                ext_first_mid_T = 0.5_dk * (extension_profile%temperature_edges(1) + extension_profile%temperature_edges(2))
                T_anchor_offset = temp_mid(column_layout%n_mpas_layers) - ext_first_mid_T
                do j = 1, column_layout%n_ext_layers
                    air_nd(column_layout%n_mpas_layers + j) = 0.5_dk * (extension_profile%air_edges(j) + extension_profile%air_edges(j + 1))
                    temp_mid(column_layout%n_mpas_layers + j) = 0.5_dk * (extension_profile%temperature_edges(j) + extension_profile%temperature_edges(j + 1)) + &
                                                  T_anchor_offset
                    o3_nd(column_layout%n_mpas_layers + j) = 0.5_dk * (extension_profile%o3_edges(j) + extension_profile%o3_edges(j + 1))
                end do
            end block
        end if

        o2_nd = VMR_O2 * air_nd
        n2_nd = VMR_N2 * air_nd
        h2o_nd = 0.0_dk
        do k = 1, column_layout%n_mpas_layers
            if (real(qo2(k), dk) > 0.0_dk) then
                o2_nd(k) = real(qo2(k), dk) * real(rho_air(k), dk) * (NA / M_O2) * 1.0e-6_dk
            end if
            if (real(qn2(k), dk) > 0.0_dk) then
                n2_nd(k) = real(qn2(k), dk) * real(rho_air(k), dk) * (NA / M_N2) * 1.0e-6_dk
            end if
            if (real(qh2o(k), dk) > 0.0_dk) then
                h2o_nd(k) = real(qh2o(k), dk) * real(rho_air(k), dk) * (NA / M_H2O) * 1.0e-6_dk
            end if
        end do

        call midpoint_to_edge_values(air_nd, air_edges)
        call midpoint_to_edge_values(temp_mid, temp_edges)
        call midpoint_to_edge_values(o3_nd, o3_edges)
        call midpoint_to_edge_values(o2_nd, o2_edges)
        call midpoint_to_edge_values(n2_nd, n2_edges)
        call midpoint_to_edge_values(h2o_nd, h2o_edges)

        air_ld = air_nd * dz_cm
        o3_ld  = o3_nd  * dz_cm
        o2_ld  = o2_nd  * dz_cm
        n2_ld  = n2_nd  * dz_cm
        h2o_ld = h2o_nd * dz_cm

    end subroutine build_column_profiles

    !------------------------------------------------------------------------
    !  subroutine update_tuvx_column_state
    !
    !> \brief Push one converted column into the public TUV-x objects.
    !------------------------------------------------------------------------
    subroutine update_tuvx_column_state(z_km, z_mid_km, &
                                        air_nd, air_edges, air_ld, &
                                        temp_mid, temp_edges, &
                                        o3_nd, o3_edges, o3_ld, &
                                        o2_nd, o2_edges, o2_ld, &
                                        n2_nd, n2_edges, n2_ld, &
                                        h2o_nd, h2o_edges, h2o_ld, &
                                        error_code, error_message)

        real(dk), intent(in) :: z_km(:), z_mid_km(:)
        real(dk), intent(in) :: air_nd(:), air_edges(:), air_ld(:)
        real(dk), intent(in) :: temp_mid(:), temp_edges(:)
        real(dk), intent(in) :: o3_nd(:), o3_edges(:), o3_ld(:)
        real(dk), intent(in) :: o2_nd(:), o2_edges(:), o2_ld(:)
        real(dk), intent(in) :: n2_nd(:), n2_edges(:), n2_ld(:)
        real(dk), intent(in) :: h2o_nd(:), h2o_edges(:), h2o_ld(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err

        error_code = 0
        error_message = ''

        call tuvx_objects%height_grid%set_edges(z_km, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call tuvx_objects%height_grid%set_midpoints(z_mid_km, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        call set_density_profile(tuvx_objects%air_profile, air_nd, air_edges, air_ld, error_code, error_message)
        if (error_code /= 0) return

        call tuvx_objects%air_profile%calculate_exo_layer_density(7.0_dk, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        call tuvx_objects%temp_profile%set_midpoint_values(temp_mid, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call tuvx_objects%temp_profile%set_edge_values(temp_edges, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

        call set_density_profile(tuvx_objects%o3_profile, o3_nd, o3_edges, o3_ld, error_code, error_message)
        if (error_code /= 0) return
        call set_density_profile(tuvx_objects%o2_profile, o2_nd, o2_edges, o2_ld, error_code, error_message)
        if (error_code /= 0) return
        call set_density_profile(tuvx_objects%n2_profile, n2_nd, n2_edges, n2_ld, error_code, error_message)
        if (error_code /= 0) return
        call set_density_profile(tuvx_objects%h2o_profile, h2o_nd, h2o_edges, h2o_ld, error_code, error_message)

    end subroutine update_tuvx_column_state

    !------------------------------------------------------------------------
    !  subroutine set_density_profile
    !
    !> \brief Set midpoint, edge, and layer-density values on a TUV-x profile.
    !------------------------------------------------------------------------
    subroutine set_density_profile(profile, midpoint_values, edge_values, layer_densities, &
                                   error_code, error_message)

        type(profile_t), pointer, intent(inout) :: profile
        real(dk), intent(in) :: midpoint_values(:), edge_values(:), layer_densities(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err

        error_code = 0
        error_message = ''

        call profile%set_midpoint_values(midpoint_values, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profile%set_edge_values(edge_values, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call profile%set_layer_densities(layer_densities, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

    end subroutine set_density_profile

    !------------------------------------------------------------------------
    !  subroutine update_cloud_radiator
    !
    !> \brief Convert MPAS cloud/rain water to TUV-x cloud optical properties.
    !------------------------------------------------------------------------
    subroutine update_cloud_radiator(qc, qr, rho_air, dz_m_mpas, cloud_tau_1d, error_code, error_message)

        real(kind=RKIND), intent(in) :: qc(:), qr(:), rho_air(:)
        real(dk), intent(in) :: dz_m_mpas(:)
        real(dk), intent(out) :: cloud_tau_1d(:)
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        type(error_t) :: err
        integer :: iWave

        error_code = 0
        error_message = ''

        call estimate_layer_cloud_tau( &
            real(qc(1:column_layout%n_mpas_layers), dk), real(qr(1:column_layout%n_mpas_layers), dk), &
            real(rho_air(1:column_layout%n_mpas_layers), dk), dz_m_mpas, cloud_tau_1d)

        do iWave = 1, tuvx_work%n_wavelength_bins
            tuvx_work%cloud_od(1:column_layout%n_mpas_layers, iWave) = cloud_tau_1d(:)
        end do

        call tuvx_objects%cloud_radiator%set_optical_depths(tuvx_work%cloud_od, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call tuvx_objects%cloud_radiator%set_single_scattering_albedos(tuvx_work%cloud_ssa, err)
        if (capture_tuvx_error(err, error_code, error_message)) return
        call tuvx_objects%cloud_radiator%set_asymmetry_factors(tuvx_work%cloud_g, err)
        if (capture_tuvx_error(err, error_code, error_message)) return

    end subroutine update_cloud_radiator

    !------------------------------------------------------------------------
    !  subroutine copy_tuvx_rates_to_mpas
    !
    !> \brief Average TUV-x edge rates to MPAS layers.
    !------------------------------------------------------------------------
    subroutine copy_tuvx_rates_to_mpas(photo_rates_out)

        real(kind=RKIND), intent(out) :: photo_rates_out(:,:)

        integer :: ir, k

        do ir = 1, tuvx_work%n_photo_rates
            do k = 1, column_layout%n_mpas_layers
                photo_rates_out(ir, k) = real( &
                    0.5_dk * (tuvx_work%photo_rates(k, tuvx_work%tuvx_indices(ir)) + &
                              tuvx_work%photo_rates(k + 1, tuvx_work%tuvx_indices(ir))), RKIND)
            end do
        end do

    end subroutine copy_tuvx_rates_to_mpas

    !------------------------------------------------------------------------
    !  subroutine midpoint_to_edge_values
    !
    !> \brief Create simple edge values from midpoint values.
    !------------------------------------------------------------------------
    subroutine midpoint_to_edge_values(midpoints, edges)

        real(dk), intent(in) :: midpoints(:)
        real(dk), intent(out) :: edges(:)

        integer :: k, n

        n = size(midpoints)
        edges(1) = midpoints(1)
        do k = 2, n
            edges(k) = 0.5_dk * (midpoints(k - 1) + midpoints(k))
        end do
        edges(n + 1) = midpoints(n)

    end subroutine midpoint_to_edge_values

    !------------------------------------------------------------------------
    !  subroutine read_upper_column_extension
    !
    !> \brief Read upper-atmosphere climatology CSV for the TUV-x extension
    !> \details
    !>  Expects a header line followed by edge values, one row per edge,
    !>  ordered bottom-up. Columns:
    !>    z_km, T_K, n_air_molec_cm3, n_O3_molec_cm3
    !>  With N rows after the header, column_layout%n_ext_layers = N - 1.
    !------------------------------------------------------------------------
    subroutine read_upper_column_extension(filename, error_code, error_message)

        use mpas_io_units, only: mpas_new_unit, mpas_release_unit
        use mpas_log,      only: mpas_log_write

        character(len=*), intent(in)               :: filename
        integer, intent(out)                       :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: iunit, ios, nrows
        character(len=256) :: msg

        error_code = 0
        error_message = ''

        call mpas_new_unit(iunit)
        open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            error_code = 1
            error_message = '[TUV-x] Failed to open extension file: ' // trim(filename)
            call mpas_release_unit(iunit)
            return
        end if

        call count_extension_rows(iunit, nrows, ios)
        if (ios /= 0) then
            error_code = 1
            error_message = '[TUV-x] Failed reading extension file header: ' // trim(filename)
            close(iunit)
            call mpas_release_unit(iunit)
            return
        end if

        if (nrows < 2) then
            error_code = 1
            error_message = '[TUV-x] Extension file must contain at least 2 data rows (1 layer)'
            close(iunit)
            call mpas_release_unit(iunit)
            return
        end if

        column_layout%n_ext_layers = nrows - 1

        call clear_extension_profile()
        allocate(extension_profile%z_km_edges (nrows))
        allocate(extension_profile%temperature_edges  (nrows))
        allocate(extension_profile%air_edges(nrows))
        allocate(extension_profile%o3_edges (nrows))

        rewind(iunit)
        call read_extension_profile_rows(iunit, filename, nrows, error_code, error_message)

        close(iunit)
        call mpas_release_unit(iunit)

        if (error_code /= 0) return

        call validate_extension_profile(nrows, error_code, error_message)
        if (error_code /= 0) return

        write(msg, '(A,F6.1,A,F6.1,A,I0,A)') '[TUV-x] Extension: z = ', &
             extension_profile%z_km_edges(1), ' -> ', extension_profile%z_km_edges(nrows), ' km (', &
             column_layout%n_ext_layers, ' layers)'
        call mpas_log_write(trim(msg))

    end subroutine read_upper_column_extension

    !------------------------------------------------------------------------
    !  subroutine clear_extension_profile
    !------------------------------------------------------------------------
    subroutine clear_extension_profile()

        if (allocated(extension_profile%z_km_edges)) deallocate(extension_profile%z_km_edges)
        if (allocated(extension_profile%temperature_edges)) deallocate(extension_profile%temperature_edges)
        if (allocated(extension_profile%air_edges)) deallocate(extension_profile%air_edges)
        if (allocated(extension_profile%o3_edges)) deallocate(extension_profile%o3_edges)

    end subroutine clear_extension_profile

    !------------------------------------------------------------------------
    !  subroutine count_extension_rows
    !------------------------------------------------------------------------
    subroutine count_extension_rows(iunit, nrows, header_status)

        integer, intent(in) :: iunit
        integer, intent(out) :: nrows
        integer, intent(out) :: header_status

        integer :: ios
        character(len=512) :: line

        nrows = 0
        read(iunit, '(A)', iostat=header_status) line
        if (header_status /= 0) return

        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) nrows = nrows + 1
        end do

    end subroutine count_extension_rows

    !------------------------------------------------------------------------
    !  subroutine read_extension_profile_rows
    !------------------------------------------------------------------------
    subroutine read_extension_profile_rows(iunit, filename, nrows, error_code, error_message)

        integer, intent(in) :: iunit, nrows
        character(len=*), intent(in) :: filename
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: ios, row
        character(len=512) :: line
        character(len=256) :: msg

        error_code = 0
        error_message = ''

        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) then
            error_code = 1
            error_message = '[TUV-x] Failed reading extension file header: ' // trim(filename)
            return
        end if

        row = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            row = row + 1
            if (row > nrows) exit
            read(line, *, iostat=ios) extension_profile%z_km_edges(row), &
                                     extension_profile%temperature_edges(row), &
                                     extension_profile%air_edges(row), &
                                     extension_profile%o3_edges(row)
            if (ios /= 0) then
                error_code = 1
                write(msg, '(A,I0,A,A)') '[TUV-x] Failed parsing data row ', row, ' of ', trim(filename)
                error_message = trim(msg)
                return
            end if
        end do

        if (row /= nrows) then
            error_code = 1
            write(msg, '(A,I0,A,I0,A,A)') '[TUV-x] Expected ', nrows, ' extension rows; read ', row, &
                                          ' from ', trim(filename)
            error_message = trim(msg)
        end if

    end subroutine read_extension_profile_rows

    !------------------------------------------------------------------------
    !  subroutine validate_extension_profile
    !------------------------------------------------------------------------
    subroutine validate_extension_profile(nrows, error_code, error_message)

        integer, intent(in) :: nrows
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        integer :: row

        error_code = 0
        error_message = ''

        do row = 2, nrows
            if (extension_profile%z_km_edges(row) <= extension_profile%z_km_edges(row - 1)) then
                error_code = 1
                error_message = '[TUV-x] Extension z_km values must be strictly increasing'
                return
            end if
        end do

    end subroutine validate_extension_profile


    !------------------------------------------------------------------------
    !  subroutine estimate_layer_cloud_tau
    !
    !> \brief Estimate bulk cloud optical depth per MPAS layer.
    !> \details
    !>  tau = 3 * LWC * dz / (2 * r_eff * rho_water)
    !>  Cloud droplets (r_eff ~10 um) dominate; rain drops (r_eff ~500 um)
    !>  contribute ~50x less per unit mass.
    !------------------------------------------------------------------------
    subroutine estimate_layer_cloud_tau(qc_mmr, qr_mmr, rho_air, dz_m, tau_layer)

        real(dk), intent(in)  :: qc_mmr(:), qr_mmr(:), rho_air(:), dz_m(:)
        real(dk), intent(out) :: tau_layer(:)

        integer :: layer
        real(dk) :: condensate_column
        real(dk) :: rain_column

        tau_layer = 0.0_dk
        do layer = 1, size(tau_layer)
            condensate_column = max(qc_mmr(layer), 0.0_dk) * max(rho_air(layer), 0.0_dk) * max(dz_m(layer), 0.0_dk)
            rain_column = max(qr_mmr(layer), 0.0_dk) * max(rho_air(layer), 0.0_dk) * max(dz_m(layer), 0.0_dk)

            tau_layer(layer) = 1.5_dk * (condensate_column / (R_EFF_CLOUD * RHO_WATER) + &
                                         rain_column / (R_EFF_RAIN * RHO_WATER))
        end do

    end subroutine estimate_layer_cloud_tau

    !------------------------------------------------------------------------
    !  function capture_tuvx_error
    !------------------------------------------------------------------------
    logical function capture_tuvx_error(err, error_code, error_message)

        type(error_t), intent(in) :: err
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

        character(len=32) :: code_text

        capture_tuvx_error = .not. err%is_success()
        error_code = 0
        error_message = ''

        if (capture_tuvx_error) then
            error_code = err%code()
            write(code_text, '(I0)') error_code
            error_message = '[TUV-x Error] ' // trim(err%category()) // ' code=' // &
                            trim(code_text) // ': ' // trim(err%message())
        end if

    end function capture_tuvx_error
#endif

end module mpas_chem_photolysis_tuvx
