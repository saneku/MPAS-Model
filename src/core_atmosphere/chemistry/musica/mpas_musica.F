! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_musica
!
!> \brief Manages interactions with the MUSICA chemistry package
!> \author CheMPAS-A Developers
!> \date   20 August 2025
!> \details
!>  MPAS wrapper around MUSICA/MICM. It initializes MICM from the selected
!>  mechanism file, discovers mechanism species at runtime, maps MPAS
!>  q<SPECIES> tracers to MICM species, converts between MPAS mass mixing
!>  ratios and MICM concentrations, injects photolysis rate parameters,
!>  advances the MICM solver, and copies updated chemistry back to MPAS
!>  scalars.
!--------------------------------------------------------------------------
module mpas_musica

    use iso_fortran_env, only: real64
    use musica_micm,  only: micm_t, solver_stats_t, get_micm_version
    use musica_state, only: state_t

    implicit none

    private

    public :: musica_init, musica_step, musica_finalize
    public :: musica_cache_photo_indices, musica_set_photolysis_rates
    public :: musica_query_species, resolve_mpas_indices, musica_mpas_index_for_species
    public :: n_chem_species, n_photo_rp, photo_rp_names

    type :: species_entry_t
        character(len=64) :: micm_name
        character(len=64) :: mpas_name
        integer           :: micm_index
        integer           :: mpas_index
        real(kind=real64) :: molar_mass
        logical           :: maps_to_mpas
        logical           :: is_o2
        logical           :: is_n2
        logical           :: is_h2o
    end type species_entry_t

    type(micm_t),  pointer :: micm  => null()
    type(state_t), pointer :: state => null()
    logical :: musica_is_initialized = .false.
    integer :: number_of_micm_grid_cells = 0

    integer, save :: n_chem_species = 0
    type(species_entry_t), allocatable :: chem_species(:)

    integer, save :: n_photo_rp = 0
    integer, allocatable :: photo_rp_indices(:)
    character(len=64), allocatable :: photo_rp_names(:)
    integer, save :: n_micm_stall_retries_total = 0

    real(kind=real64), parameter :: mw_air_kg_per_mol = 28.9647e-3_real64
    real(kind=real64), parameter :: mw_o2_kg_per_mol = 31.998e-3_real64
    real(kind=real64), parameter :: mw_n2_kg_per_mol = 28.0134e-3_real64
    real(kind=real64), parameter :: mw_h2o_kg_per_mol = 18.01528e-3_real64
    real(kind=real64), parameter :: vmr_o2 = 0.20946_real64
    real(kind=real64), parameter :: vmr_n2 = 0.78084_real64

    contains

    !------------------------------------------------------------------------
    !  subroutine musica_init
    !------------------------------------------------------------------------
    subroutine musica_init(filename_of_micm_configuration, number_of_grid_cells, &
            error_code, error_message)

        use musica_micm, only : RosenbrockStandardOrder
        use musica_util, only : error_t, string_t
        use mpas_kind_types, only : RKIND
        use mpas_log, only : mpas_log_write

        character(len=*),              intent(in)    :: filename_of_micm_configuration
        integer,                       intent(in)    :: number_of_grid_cells
        integer,                       intent(out)   :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        character(len=*), parameter :: MOLAR_MASS_PROPERTY = "molecular weight [kg mol-1]"
        character(len=*), parameter :: MOLAR_MASS_FALLBACK_PROPERTY = "__molar mass"
        type(error_t)  :: error
        type(error_t)  :: fallback_error
        type(string_t) :: micm_version
        integer :: solver_type = RosenbrockStandardOrder
        integer :: i_species, spec_idx
        real(kind=real64) :: molar_mass_val
        character(len=:), allocatable :: sp_name
        character(len=64) :: sp_lower
        character(len=4096) :: skipped_zero_weight_names
        integer :: n_skipped_zero_weight
        logical :: is_background_species

        error_code = 0
        error_message = ''

        if (len_trim(filename_of_micm_configuration) == 0) then
            call mpas_log_write('MUSICA chemistry disabled: no MICM configuration file specified')
            return
        end if

        micm_version = get_micm_version()

        call mpas_log_write('Initializing MUSICA chemistry package...')
        call mpas_log_write('MICM configuration file: ' // trim(filename_of_micm_configuration))
        call mpas_log_write('MICM version: ' // micm_version%value_)
        call mpas_log_write('MICM number of grid cells: $i', intArgs=[number_of_grid_cells])

        micm => micm_t(trim(filename_of_micm_configuration), solver_type, error)
        if (has_error_occurred(error, error_message, error_code)) return

        state => micm%get_state(number_of_grid_cells, error)
        if (has_error_occurred(error, error_message, error_code)) return

        number_of_micm_grid_cells = number_of_grid_cells
        call initialize_rate_parameters_defaults(state)

        n_chem_species = state%species_ordering%size()
        if (allocated(chem_species)) deallocate(chem_species)
        allocate(chem_species(n_chem_species))
        skipped_zero_weight_names = ''
        n_skipped_zero_weight = 0

        do i_species = 1, n_chem_species
            sp_name = trim(state%species_ordering%name(i_species))
            sp_lower = to_lower(trim(sp_name))
            is_background_species = is_background_micm_species(sp_lower)

            spec_idx = state%species_ordering%index(sp_name, error)
            if (has_error_occurred(error, error_message, error_code)) return

            molar_mass_val = micm%get_species_property_double(sp_name, MOLAR_MASS_PROPERTY, error)
            if (.not. error%is_success()) then
                molar_mass_val = micm%get_species_property_double(sp_name, MOLAR_MASS_FALLBACK_PROPERTY, fallback_error)
                error = fallback_error
                if (has_error_occurred(error, error_message, error_code)) return
            end if
            if (molar_mass_val <= 0.0_real64) then
                if (is_background_species) then
                    molar_mass_val = background_molar_mass(sp_lower)
                else
                    chem_species(i_species)%micm_name = sp_name
                    chem_species(i_species)%mpas_name = trim(to_lower(sp_name))
                    chem_species(i_species)%micm_index = spec_idx
                    chem_species(i_species)%mpas_index = 0
                    chem_species(i_species)%molar_mass = 0.0_real64
                    chem_species(i_species)%maps_to_mpas = .false.
                    chem_species(i_species)%is_o2 = .false.
                    chem_species(i_species)%is_n2 = .false.
                    chem_species(i_species)%is_h2o = .false.
                    n_skipped_zero_weight = n_skipped_zero_weight + 1
                    call append_name_to_string(skipped_zero_weight_names, sp_name)
                    cycle
                end if
            end if

            chem_species(i_species)%micm_name = sp_name
            chem_species(i_species)%mpas_name = trim(to_lower(sp_name))
            chem_species(i_species)%micm_index = spec_idx
            chem_species(i_species)%mpas_index = 0
            chem_species(i_species)%molar_mass = molar_mass_val
            chem_species(i_species)%maps_to_mpas = .true.
            chem_species(i_species)%is_o2 = trim(sp_lower) == 'o2'
            chem_species(i_species)%is_n2 = trim(sp_lower) == 'n2'
            chem_species(i_species)%is_h2o = trim(sp_lower) == 'h2o'

            call mpas_log_write('[MUSICA] MICM species: ' // trim(sp_name) // &
                ' -> ' // trim(chem_species(i_species)%mpas_name) // &
                ' molar_mass=$r kg/mol', realArgs=[real(molar_mass_val, kind=RKIND)])
        end do

        if (n_skipped_zero_weight > 0) then
            call mpas_log_write('[MUSICA] Skipped MPAS mapping for zero-weight MICM species: $i', &
                                intArgs=[n_skipped_zero_weight])
            call mpas_log_write('[MUSICA] Zero-weight MICM species not mapped to MPAS: ' // &
                                trim(skipped_zero_weight_names))
        end if

        call mpas_log_write('MICM rate parameters: $i (storage size=$i)', &
                            intArgs=[state%number_of_rate_parameters, size(state%rate_parameters)])

        musica_is_initialized = .true.

    end subroutine musica_init

    !------------------------------------------------------------------------
    !  subroutine musica_query_species
    !------------------------------------------------------------------------
    subroutine musica_query_species(config_path, num_species, species_names, error_code, error_message)

        use musica_micm,  only: RosenbrockStandardOrder
        use musica_state, only: state_t
        use musica_util,  only: error_t
        use mpas_log,     only: mpas_log_write

        character(len=*), intent(in)                 :: config_path
        integer, intent(out)                         :: num_species
        character(len=64), allocatable, intent(out)  :: species_names(:)
        integer, intent(out)                         :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        type(micm_t), pointer  :: tmp_micm => null()
        type(state_t), pointer :: tmp_state => null()
        type(error_t)          :: error
        type(error_t)          :: fallback_error
        integer                :: i
        integer                :: n_total
        character(len=*), parameter :: MOLAR_MASS_PROPERTY = "molecular weight [kg mol-1]"
        character(len=*), parameter :: MOLAR_MASS_FALLBACK_PROPERTY = "__molar mass"
        real(kind=real64)      :: molar_mass_val
        character(len=:), allocatable :: sp_name
        character(len=64)      :: sp_lower
        character(len=4096)    :: skipped_zero_weight_names
        integer                :: n_skipped_zero_weight
        logical                :: is_background_species

        error_code = 0
        error_message = ''
        num_species = 0
        skipped_zero_weight_names = ''
        n_skipped_zero_weight = 0

        if (len_trim(config_path) == 0) then
            call mpas_log_write('[MUSICA] No MICM config path; skipping chemistry tracer discovery')
            return
        end if

        tmp_micm => micm_t(trim(config_path), RosenbrockStandardOrder, error)
        if (has_error_occurred(error, error_message, error_code)) goto 100

        tmp_state => tmp_micm%get_state(1, error)
        if (has_error_occurred(error, error_message, error_code)) goto 100

        n_total = tmp_state%species_ordering%size()
        allocate(species_names(n_total))

        do i = 1, n_total
            sp_name = trim(tmp_state%species_ordering%name(i))
            sp_lower = to_lower(trim(sp_name))
            is_background_species = is_background_micm_species(sp_lower)
            molar_mass_val = tmp_micm%get_species_property_double(sp_name, MOLAR_MASS_PROPERTY, error)
            if (.not. error%is_success()) then
                molar_mass_val = tmp_micm%get_species_property_double(sp_name, MOLAR_MASS_FALLBACK_PROPERTY, fallback_error)
                error = fallback_error
                if (has_error_occurred(error, error_message, error_code)) goto 100
            end if

            if (molar_mass_val <= 0.0_real64) then
                if (is_background_species) then
                    molar_mass_val = background_molar_mass(sp_lower)
                else
                    n_skipped_zero_weight = n_skipped_zero_weight + 1
                    call append_name_to_string(skipped_zero_weight_names, sp_name)
                    cycle
                end if
            end if

            num_species = num_species + 1
            species_names(num_species) = trim(to_lower(sp_name))
            call mpas_log_write('[MUSICA] Discovered chemistry tracer: ' // trim(species_names(num_species)))
        end do

        if (n_skipped_zero_weight > 0) then
            call mpas_log_write('[MUSICA] Skipped dynamic MPAS tracers for zero-weight MICM species: $i', &
                                intArgs=[n_skipped_zero_weight])
            call mpas_log_write('[MUSICA] Zero-weight MICM species not allocated as MPAS tracers: ' // &
                                trim(skipped_zero_weight_names))
        end if

100     continue
        if (associated(tmp_state)) deallocate(tmp_state)
        if (associated(tmp_micm)) deallocate(tmp_micm)

    end subroutine musica_query_species

    !------------------------------------------------------------------------
    !  subroutine resolve_mpas_indices
    !------------------------------------------------------------------------
    subroutine resolve_mpas_indices(statePool, error_code, error_message)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_dimension
        use mpas_log, only: mpas_log_write

        type(mpas_pool_type), intent(in)             :: statePool
        integer, intent(out)                         :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        integer :: i_sp
        integer, pointer :: idx_ptr

        error_code = 0
        error_message = ''

        do i_sp = 1, n_chem_species
            if (.not. chem_species(i_sp)%maps_to_mpas) cycle
            nullify(idx_ptr)
            call mpas_pool_get_dimension(statePool, 'index_' // trim(chem_species(i_sp)%mpas_name), idx_ptr)

            if (.not. associated(idx_ptr)) then
                error_code = 1
                error_message = '[MUSICA] No MPAS tracer index for ' // trim(chem_species(i_sp)%mpas_name)
                call mpas_log_write(error_message)
                return
            end if

            if (idx_ptr < 1) then
                error_code = 1
                error_message = '[MUSICA] Invalid MPAS tracer index for ' // trim(chem_species(i_sp)%mpas_name)
                call mpas_log_write(error_message)
                return
            end if

            chem_species(i_sp)%mpas_index = idx_ptr
            call mpas_log_write('[MUSICA] Resolved ' // trim(chem_species(i_sp)%micm_name) // &
                ' -> MPAS index $i', intArgs=[chem_species(i_sp)%mpas_index])
        end do

    end subroutine resolve_mpas_indices

    !------------------------------------------------------------------------
    !  function musica_mpas_index_for_species
    !------------------------------------------------------------------------
    integer function musica_mpas_index_for_species(species_name) result(idx)

        character(len=*), intent(in) :: species_name
        character(len=64) :: target
        integer :: i_sp

        idx = -1
        target = to_lower(trim(species_name))

        do i_sp = 1, n_chem_species
            if (.not. chem_species(i_sp)%maps_to_mpas) cycle
            if (chem_species(i_sp)%mpas_index <= 0) cycle
            if (trim(to_lower(chem_species(i_sp)%micm_name)) == trim(target)) then
                idx = chem_species(i_sp)%mpas_index
                return
            end if
        end do

    end function musica_mpas_index_for_species

    !------------------------------------------------------------------------
    !  subroutine musica_step
    !------------------------------------------------------------------------
    subroutine musica_step(scalars, qv_mmr, rho_zz, zz, temperature, pressure, time_step, &
            error_code, error_message)

        use mpas_kind_types, only : RKIND
        use mpas_log, only : mpas_log_write
        use mpas_timer, only : mpas_timer_start, mpas_timer_stop
        use musica_util, only : error_t, string_t
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        real(kind=RKIND), dimension(:,:), intent(in)      :: qv_mmr
        real(kind=RKIND), dimension(:,:), intent(in)      :: rho_zz
        real(kind=RKIND), dimension(:,:), intent(in)      :: zz
        real(kind=RKIND), dimension(:,:), intent(in)      :: temperature
        real(kind=RKIND), dimension(:,:), intent(in)      :: pressure
        real(kind=RKIND),                intent(in)       :: time_step
        integer,                         intent(out)      :: error_code
        character(len=:), allocatable,   intent(out)      :: error_message

        type(error_t)        :: error
        type(string_t)       :: solver_state
        type(solver_stats_t) :: solver_stats
        integer              :: nVertLevels, nCells
        integer              :: iCell, k, iGrid, iConc, i_sp
        real(kind=real64)    :: rho_air, t_air, p_air, air_mol, q_mmr, qv_air
        real(kind=real64)    :: elapsed, remaining, advanced
        real(kind=real64)    :: attempt_dt, min_attempt_dt
        integer              :: sub_call
        integer, parameter   :: MAX_SUB_CALLS = 100
        integer, parameter   :: MAX_STALL_RETRIES = 20
        integer              :: stall_retry
        integer              :: n_bad_conc, n_bad_rate
        integer              :: cell_stride, var_stride
        real(kind=real64), parameter :: rho_air_floor = 1.0e-12_real64
        real(kind=real64), parameter :: t_air_min = 150.0_real64
        real(kind=real64), parameter :: t_air_max = 350.0_real64
        real(kind=real64), parameter :: p_air_min = 100.0_real64
        real(kind=real64), parameter :: p_air_max = 200000.0_real64

        error_code = 0
        error_message = ''

        if (.not. musica_is_initialized) return

        nVertLevels = size(qv_mmr, 1)
        nCells = size(qv_mmr, 2)

        if (size(scalars, 2) < nVertLevels .or. size(scalars, 3) < nCells) then
            error_code = 1
            error_message = 'MUSICA step failed: scalars dimensions do not match chemistry grid.'
            return
        end if

        if (size(rho_zz, 1) /= nVertLevels .or. size(rho_zz, 2) /= nCells .or. &
            size(zz, 1) /= nVertLevels .or. size(zz, 2) /= nCells .or. &
            size(temperature, 1) /= nVertLevels .or. size(temperature, 2) /= nCells .or. &
            size(pressure, 1) /= nVertLevels .or. size(pressure, 2) /= nCells) then
            error_code = 1
            error_message = 'MUSICA step failed: environment dimensions do not match chemistry grid.'
            return
        end if

        if (nVertLevels * nCells /= number_of_micm_grid_cells) then
            error_code = 1
            error_message = 'MUSICA step failed: MPAS (nVertLevels*nCells) does not match MICM grid size.'
            return
        end if

        call mpas_timer_start('chem_mpas_to_micm')

        call sanitize_mmr_tracers(scalars, nVertLevels, nCells)

        cell_stride = state%species_strides%grid_cell
        var_stride  = state%species_strides%variable

        do iCell = 1, nCells
            do k = 1, nVertLevels
                iGrid = (iCell - 1) * nVertLevels + k

                rho_air = real(rho_zz(k, iCell), kind=real64) * real(zz(k, iCell), kind=real64)
                if (.not. ieee_is_finite(rho_air) .or. rho_air <= tiny(1.0_real64)) rho_air = rho_air_floor
                rho_air = max(rho_air, rho_air_floor)

                t_air = real(temperature(k, iCell), kind=real64)
                if (.not. ieee_is_finite(t_air) .or. t_air <= 0.0_real64) t_air = 280.0_real64
                t_air = min(max(t_air, t_air_min), t_air_max)

                p_air = real(pressure(k, iCell), kind=real64)
                if (.not. ieee_is_finite(p_air) .or. p_air <= 0.0_real64) p_air = 101325.0_real64
                p_air = min(max(p_air, p_air_min), p_air_max)

                state%conditions(iGrid)%temperature = t_air
                state%conditions(iGrid)%pressure = p_air
                state%conditions(iGrid)%air_density = rho_air / mw_air_kg_per_mol
                air_mol = rho_air / mw_air_kg_per_mol

                do i_sp = 1, n_chem_species
                    if (.not. chem_species(i_sp)%maps_to_mpas .and. .not. chem_species(i_sp)%is_o2 .and. &
                        .not. chem_species(i_sp)%is_n2 .and. .not. chem_species(i_sp)%is_h2o) cycle
                    iConc = (iGrid - 1) * cell_stride + &
                            (chem_species(i_sp)%micm_index - 1) * var_stride + 1

                    if (chem_species(i_sp)%is_o2) then
                        q_mmr = 0.0_real64
                        if (chem_species(i_sp)%mpas_index >= 1) then
                            q_mmr = real(scalars(chem_species(i_sp)%mpas_index, k, iCell), kind=real64)
                        end if
                        if (.not. ieee_is_finite(q_mmr) .or. q_mmr <= 0.0_real64) then
                            q_mmr = vmr_o2 * mw_o2_kg_per_mol / mw_air_kg_per_mol
                        end if
                        state%concentrations(iConc) = q_mmr * rho_air / chem_species(i_sp)%molar_mass
                    else if (chem_species(i_sp)%is_n2) then
                        q_mmr = 0.0_real64
                        if (chem_species(i_sp)%mpas_index >= 1) then
                            q_mmr = real(scalars(chem_species(i_sp)%mpas_index, k, iCell), kind=real64)
                        end if
                        if (.not. ieee_is_finite(q_mmr) .or. q_mmr <= 0.0_real64) then
                            q_mmr = vmr_n2 * mw_n2_kg_per_mol / mw_air_kg_per_mol
                        end if
                        state%concentrations(iConc) = q_mmr * rho_air / chem_species(i_sp)%molar_mass
                    else if (chem_species(i_sp)%is_h2o) then
                        qv_air = real(qv_mmr(k, iCell), kind=real64)
                        if (.not. ieee_is_finite(qv_air) .or. qv_air < 0.0_real64) qv_air = 0.0_real64
                        q_mmr = qv_air
                        if (q_mmr <= 0.0_real64 .and. chem_species(i_sp)%mpas_index >= 1) then
                            q_mmr = real(scalars(chem_species(i_sp)%mpas_index, k, iCell), kind=real64)
                            if (.not. ieee_is_finite(q_mmr) .or. q_mmr < 0.0_real64) q_mmr = 0.0_real64
                        end if
                        state%concentrations(iConc) = q_mmr * rho_air / chem_species(i_sp)%molar_mass
                    else
                        if (chem_species(i_sp)%mpas_index < 1) cycle
                        q_mmr = max(0.0_real64, real(scalars(chem_species(i_sp)%mpas_index, k, iCell), kind=real64))
                        state%concentrations(iConc) = q_mmr * rho_air / chem_species(i_sp)%molar_mass
                    end if
                end do
            end do
        end do

        call sanitize_micm_state(state, n_bad_conc, n_bad_rate)
        if (n_bad_conc > 0 .or. n_bad_rate > 0) then
            call mpas_log_write('[MUSICA] Sanitized MICM inputs (NaN/Inf): concentrations=$i rate_parameters=$i', &
                                intArgs=[n_bad_conc, n_bad_rate])
        end if

        call mpas_timer_stop('chem_mpas_to_micm')

        elapsed = 0.0_real64
        sub_call = 0
        min_attempt_dt = max(1.0e-12_real64, 1.0e-6_real64 * real(time_step, kind=real64))

        call mpas_timer_start('chem_micm_solve')

        do while (elapsed < real(time_step, kind=real64) .and. sub_call < MAX_SUB_CALLS)
            remaining = real(time_step, kind=real64) - elapsed
            attempt_dt = remaining
            stall_retry = 0

            do
                sub_call = sub_call + 1
                call micm%solve(attempt_dt, state, solver_state, solver_stats, error)
                if (has_error_occurred(error, error_message, error_code)) then
                    call mpas_timer_stop('chem_micm_solve')
                    return
                end if

                advanced = solver_stats%final_time()
                if (ieee_is_finite(advanced) .and. advanced > 0.0_real64) exit

                call sanitize_micm_state(state, n_bad_conc, n_bad_rate)
                stall_retry = stall_retry + 1
                n_micm_stall_retries_total = n_micm_stall_retries_total + 1

                if (stall_retry > MAX_STALL_RETRIES .or. attempt_dt <= min_attempt_dt) then
                    error_code = 1
                    error_message = '[MUSICA] MICM stalled at sub-step; solver state: ' // &
                                    solver_state%get_char_array()
                    call mpas_timer_stop('chem_micm_solve')
                    return
                end if

                attempt_dt = max(0.5_real64 * attempt_dt, min_attempt_dt)
            end do

            if (stall_retry > 0) then
                call mpas_log_write('[MUSICA] MICM recovered after stall retries: $i', &
                                    intArgs=[stall_retry])
            end if

            elapsed = elapsed + advanced
        end do

        call mpas_timer_stop('chem_micm_solve')

        if (elapsed < real(time_step, kind=real64)) then
            error_code = 1
            error_message = '[MUSICA] MICM sub-step budget exhausted.'
            return
        end if

        call mpas_timer_start('chem_micm_to_mpas')

        do iCell = 1, nCells
            do k = 1, nVertLevels
                iGrid = (iCell - 1) * nVertLevels + k
                rho_air = real(rho_zz(k, iCell), kind=real64) * real(zz(k, iCell), kind=real64)
                if (.not. ieee_is_finite(rho_air) .or. rho_air <= tiny(1.0_real64)) rho_air = rho_air_floor
                rho_air = max(rho_air, rho_air_floor)

                do i_sp = 1, n_chem_species
                    if (.not. chem_species(i_sp)%maps_to_mpas .and. .not. chem_species(i_sp)%is_o2 .and. &
                        .not. chem_species(i_sp)%is_n2 .and. .not. chem_species(i_sp)%is_h2o) cycle
                    iConc = (iGrid - 1) * cell_stride + &
                            (chem_species(i_sp)%micm_index - 1) * var_stride + 1
                    if (chem_species(i_sp)%mpas_index < 1) cycle
                    q_mmr = max(0.0_real64, safe_finite(state%concentrations(iConc))) * &
                            chem_species(i_sp)%molar_mass / rho_air
                    scalars(chem_species(i_sp)%mpas_index, k, iCell) = real(q_mmr, kind=RKIND)
                end do
            end do
        end do

        call mpas_timer_stop('chem_micm_to_mpas')

    end subroutine musica_step

    !------------------------------------------------------------------------
    !  subroutine musica_cache_photo_indices
    !------------------------------------------------------------------------
    subroutine musica_cache_photo_indices(names, error_code, error_message)

        use mpas_log, only : mpas_log_write
        use musica_util, only : error_t

        character(len=*), intent(in)                :: names(:)
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)  :: error_message

        integer :: i
        character(len=:), allocatable :: target_name
        type(error_t) :: error

        error_code = 0
        error_message = ''

        if (.not. associated(state)) then
            error_code = 1
            error_message = '[MUSICA] Cannot cache photolysis rates before MICM state is initialized.'
            return
        end if

        if (state%number_of_rate_parameters == 0) then
            error_code = 1
            error_message = '[MUSICA] MICM mechanism has no rate parameters for photolysis.'
            return
        end if

        n_photo_rp = size(names)
        if (allocated(photo_rp_names)) deallocate(photo_rp_names)
        if (allocated(photo_rp_indices)) deallocate(photo_rp_indices)
        allocate(photo_rp_names(n_photo_rp))
        allocate(photo_rp_indices(n_photo_rp))

        do i = 1, n_photo_rp
            photo_rp_names(i) = trim(names(i))
            target_name = 'PHOTO.' // trim(names(i))
            photo_rp_indices(i) = state%rate_parameters_ordering%index(target_name, error)
            if (.not. error%is_success()) then
                error_code = 1
                error_message = '[MUSICA] Missing MICM rate parameter ' // target_name
                return
            end if
            call mpas_log_write('[MUSICA] Cached ' // target_name // &
                ' rate parameter index = $i', intArgs=[photo_rp_indices(i)])
        end do

    end subroutine musica_cache_photo_indices

    !------------------------------------------------------------------------
    !  subroutine musica_set_photolysis_rates
    !------------------------------------------------------------------------
    subroutine musica_set_photolysis_rates(rates, nCells, nVertLevels, n_rates, &
                                           error_code, error_message)

        use mpas_kind_types, only : RKIND

        real(kind=RKIND), intent(in)                :: rates(:,:,:) ! (rate, level, cell)
        integer, intent(in)                         :: nCells
        integer, intent(in)                         :: nVertLevels
        integer, intent(in)                         :: n_rates
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)  :: error_message

        integer :: cell_stride, var_stride
        integer :: iCell, k, ir, iGrid, iRate

        error_code = 0
        error_message = ''

        if (n_photo_rp < 1) return
        if (.not. associated(state)) return

        if (n_rates /= n_photo_rp) then
            error_code = 1
            error_message = '[MUSICA] Photolysis rate-count mismatch.'
            return
        end if

        if (nCells * nVertLevels /= number_of_micm_grid_cells) then
            error_code = 1
            error_message = '[MUSICA] Photolysis grid size does not match MICM grid size.'
            return
        end if

        cell_stride = state%rate_parameters_strides%grid_cell
        var_stride = state%rate_parameters_strides%variable

        do ir = 1, n_photo_rp
            iRate = photo_rp_indices(ir)
            do iCell = 1, nCells
                do k = 1, nVertLevels
                    iGrid = (iCell - 1) * nVertLevels + k
                    state%rate_parameters((iGrid - 1) * cell_stride + &
                        (iRate - 1) * var_stride + 1) = rates(ir, k, iCell)
                end do
            end do
        end do

    end subroutine musica_set_photolysis_rates

    !------------------------------------------------------------------------
    !  routine initialize_rate_parameters_defaults
    !------------------------------------------------------------------------
    subroutine initialize_rate_parameters_defaults(state_in)

        use musica_util, only : error_t

        type(state_t), pointer, intent(inout) :: state_in

        integer :: cell_stride, var_stride, rp_index, cell, storage_index, n_rp
        character(len=:), allocatable :: rp_name
        real(kind=real64) :: rate_val
        type(error_t) :: error

        if (.not. associated(state_in)) return
        if (state_in%number_of_rate_parameters <= 0) return

        cell_stride = state_in%rate_parameters_strides%grid_cell
        var_stride  = state_in%rate_parameters_strides%variable
        n_rp = state_in%rate_parameters_ordering%size()

        do rp_index = 1, n_rp
            rp_name = trim(state_in%rate_parameters_ordering%name(rp_index))
            storage_index = state_in%rate_parameters_ordering%index(rp_name, error)
            if (.not. error%is_success()) cycle

            if (rp_name(1:min(5, len(rp_name))) == 'USER.') then
                rate_val = 1.0_real64
            else
                rate_val = 0.0_real64
            end if

            do cell = 1, number_of_micm_grid_cells
                state_in%rate_parameters((cell - 1) * cell_stride + &
                                         (storage_index - 1) * var_stride + 1) = rate_val
            end do
        end do

    end subroutine initialize_rate_parameters_defaults

    !------------------------------------------------------------------------
    !  routine sanitize_micm_state
    !------------------------------------------------------------------------
    subroutine sanitize_micm_state(state_in, n_bad_conc, n_bad_rate)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        type(state_t), pointer, intent(inout) :: state_in
        integer, intent(out) :: n_bad_conc, n_bad_rate
        integer :: i

        n_bad_conc = 0
        n_bad_rate = 0

        if (.not. associated(state_in)) return

        do i = 1, size(state_in%concentrations)
            if (.not. ieee_is_finite(state_in%concentrations(i)) .or. state_in%concentrations(i) < 0.0_real64) then
                state_in%concentrations(i) = 0.0_real64
                n_bad_conc = n_bad_conc + 1
            end if
        end do

        do i = 1, size(state_in%rate_parameters)
            if (.not. ieee_is_finite(state_in%rate_parameters(i))) then
                state_in%rate_parameters(i) = 0.0_real64
                n_bad_rate = n_bad_rate + 1
            end if
        end do

    end subroutine sanitize_micm_state

    !------------------------------------------------------------------------
    !  routine sanitize_mmr_tracers
    !------------------------------------------------------------------------
    subroutine sanitize_mmr_tracers(scalars, nVertLevels, nCells)

        use mpas_kind_types, only : RKIND
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        integer, intent(in) :: nVertLevels, nCells
        integer :: i_sp, iCell, k, idx

        do i_sp = 1, n_chem_species
            if (.not. chem_species(i_sp)%maps_to_mpas) cycle
            idx = chem_species(i_sp)%mpas_index
            if (idx < 1 .or. idx > size(scalars, 1)) cycle
            do iCell = 1, nCells
                do k = 1, nVertLevels
                    if (.not. ieee_is_finite(scalars(idx, k, iCell)) .or. scalars(idx, k, iCell) < 0.0_RKIND) then
                        scalars(idx, k, iCell) = 0.0_RKIND
                    else if (scalars(idx, k, iCell) > 1.0_RKIND) then
                        scalars(idx, k, iCell) = 1.0_RKIND
                    end if
                end do
            end do
        end do

    end subroutine sanitize_mmr_tracers

    !------------------------------------------------------------------------
    !  function safe_finite
    !------------------------------------------------------------------------
    pure real(kind=real64) function safe_finite(value)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(kind=real64), intent(in) :: value

        if (ieee_is_finite(value)) then
            safe_finite = value
        else
            safe_finite = 0.0_real64
        end if
    end function safe_finite

    !------------------------------------------------------------------------
    !  function to_lower
    !------------------------------------------------------------------------
    pure function to_lower(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i, ia

        do i = 1, len(text)
            ia = iachar(text(i:i))
            if (ia >= iachar('A') .and. ia <= iachar('Z')) then
                out(i:i) = achar(ia + 32)
            else
                out(i:i) = text(i:i)
            end if
        end do
    end function to_lower

    !------------------------------------------------------------------------
    !  function is_background_micm_species
    !------------------------------------------------------------------------
    pure logical function is_background_micm_species(species_name)

        character(len=*), intent(in) :: species_name
        character(len=len(species_name)) :: lower_name

        lower_name = to_lower(trim(species_name))
        is_background_micm_species = trim(lower_name) == 'o2' .or. &
                                     trim(lower_name) == 'n2' .or. &
                                     trim(lower_name) == 'h2o'

    end function is_background_micm_species

    !------------------------------------------------------------------------
    !  function background_molar_mass
    !------------------------------------------------------------------------
    pure real(kind=real64) function background_molar_mass(species_name)

        character(len=*), intent(in) :: species_name
        character(len=len(species_name)) :: lower_name

        lower_name = to_lower(trim(species_name))
        select case (trim(lower_name))
        case ('o2')
            background_molar_mass = mw_o2_kg_per_mol
        case ('n2')
            background_molar_mass = mw_n2_kg_per_mol
        case ('h2o')
            background_molar_mass = mw_h2o_kg_per_mol
        case default
            background_molar_mass = 0.0_real64
        end select

    end function background_molar_mass

    !------------------------------------------------------------------------
    !  subroutine append_name_to_string
    !------------------------------------------------------------------------
    subroutine append_name_to_string(list, name)

        character(len=*), intent(inout) :: list
        character(len=*), intent(in)    :: name

        integer :: used, ncopy, room

        if (len_trim(name) == 0) return

        used = len_trim(list)
        if (used == 0) then
            room = len(list)
            ncopy = min(len_trim(name), room)
            if (ncopy > 0) list(1:ncopy) = name(1:ncopy)
            return
        end if

        if (used + 2 > len(list)) return
        list(used + 1:used + 2) = ', '
        used = used + 2

        room = len(list) - used
        ncopy = min(len_trim(name), room)
        if (ncopy > 0) list(used + 1:used + ncopy) = name(1:ncopy)

    end subroutine append_name_to_string

    !------------------------------------------------------------------------
    !  subroutine musica_finalize
    !------------------------------------------------------------------------
    subroutine musica_finalize()

        use mpas_log, only : mpas_log_write

        if (musica_is_initialized) then
            call mpas_log_write('Finalizing MUSICA chemistry package...')
        end if
        if (n_micm_stall_retries_total > 0) then
            call mpas_log_write('[MUSICA] MICM stall retries total: $i', &
                                intArgs=[n_micm_stall_retries_total])
        end if

        number_of_micm_grid_cells = 0
        n_photo_rp = 0
        n_chem_species = 0
        n_micm_stall_retries_total = 0
        if (allocated(photo_rp_names)) deallocate(photo_rp_names)
        if (allocated(photo_rp_indices)) deallocate(photo_rp_indices)
        if (allocated(chem_species)) deallocate(chem_species)
        if (associated(state)) deallocate(state)
        if (associated(micm)) deallocate(micm)
        musica_is_initialized = .false.

    end subroutine musica_finalize

    !------------------------------------------------------------------------
    !  function has_error_occurred
    !------------------------------------------------------------------------
    logical function has_error_occurred(error, error_message, error_code)
        use musica_util, only: error_t

        type(error_t),                 intent(in)  :: error
        character(len=:), allocatable, intent(out) :: error_message
        integer,                       intent(out) :: error_code

        character(len=30) :: error_code_str

        if ( error%is_success( ) ) then
          error_code = 0
          error_message = ''
          has_error_occurred = .false.
          return
        end if
        error_code = error%code( )
        write(error_code_str, '(I30)') error%code( )
        error_message = '[MUSICA Error]: ' // error%category( ) // '[' // &
                        trim( adjustl( error_code_str ) ) // ']: ' // error%message( )
        has_error_occurred = .true.

    end function has_error_occurred

end module mpas_musica
