! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_chem_photolysis_driver
!
!> \brief Drives solar geometry and TUV-x photolysis for MPAS chemistry
!> \details
!>  Owns the high-level photolysis workflow: initializes optional TUV-x,
!>  updates cos_sza diagnostics, prepares vertical column inputs for TUV-x,
!>  sanitizes returned photolysis rates, writes photolysis diagnostics, and
!>  passes the rates to MUSICA/MICM before chemistry is stepped.
!-----------------------------------------------------------------------
module mpas_chem_photolysis_driver

    use mpas_kind_types, only: RKIND, StrKIND

    implicit none

    private

    public :: photolysis_init
    public :: photolysis_step
    public :: photolysis_update_due
    public :: photolysis_advance_time
    public :: photolysis_enabled
    public :: photolysis_finalize

    real(kind=RKIND), save :: tuvx_update_interval = 0.0_RKIND
    real(kind=RKIND), save :: tuvx_time_since_last = huge(1.0_RKIND)
    logical, save :: use_tuvx = .false.

contains

    !------------------------------------------------------------------------
    !  routine photolysis_init
    !------------------------------------------------------------------------
    subroutine photolysis_init(configs, structs, clock, nVertLevels, micm_mode, error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_cache_photo_indices
        use mpas_tuvx, only: tuvx_init, tuvx_is_enabled, tuvx_n_rates, tuvx_rate_names
#endif
        use mpas_chem_utils, only: chemistry_days_in_year
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT, MPAS_NOW, MPAS_Time_type
        use mpas_log, only : mpas_log_write
        use mpas_pool_routines, only: mpas_pool_get_config, mpas_pool_get_subpool, mpas_pool_get_dimension
        use mpas_timekeeping, only: MPAS_Clock_type, mpas_get_clock_time, mpas_get_time

        type(mpas_pool_type), intent(in) :: configs
        type(mpas_pool_type), intent(in) :: structs
        type(MPAS_Clock_type), intent(in) :: clock
        integer, intent(in) :: nVertLevels
        logical, intent(in) :: micm_mode
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        type(mpas_pool_type), pointer :: mesh
        type(mpas_pool_type), pointer :: diag
        type(mpas_pool_type), pointer :: state
        integer, pointer :: nCells_ptr
        integer :: nCells
        real(kind=RKIND), pointer :: real_ptr
        character(len=StrKIND), pointer :: tuvx_filepath_ptr
        logical, pointer :: tuvx_top_extension_ptr
        character(len=StrKIND), pointer :: tuvx_extension_file_ptr
        logical :: tuvx_top_extension
        character(len=StrKIND) :: tuvx_extension_file
        type(MPAS_Time_type) :: currTime
        integer :: ierr, DoY, hour, minute, second, days_in_year
        real(kind=RKIND) :: hour_utc
#endif

        error_code = 0
        error_message = ''
        use_tuvx = .false.
        tuvx_time_since_last = huge(1.0_RKIND)

#ifdef MPAS_USE_MUSICA
        nullify(real_ptr)
        call mpas_pool_get_config(configs, 'config_tuvx_update_interval', real_ptr)
        if (associated(real_ptr)) tuvx_update_interval = max(0.0_RKIND, real_ptr)

        nullify(mesh)
        nullify(diag)
        nullify(state)
        call mpas_pool_get_subpool(structs, 'mesh', mesh)
        call mpas_pool_get_subpool(structs, 'diag', diag)
        call mpas_pool_get_subpool(structs, 'state', state)

        nCells = 0
        if (associated(state)) then
            nullify(nCells_ptr)
            call mpas_pool_get_dimension(state, 'nCells', nCells_ptr)
            if (associated(nCells_ptr)) nCells = nCells_ptr
        end if

        if (associated(diag)) then
            call zero_photolysis_diag(diag, ['j_NO2    ', 'j_O2     ', 'j_O3_O   ', 'j_O3_O1D '])
        end if

        currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
        call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
        days_in_year = chemistry_days_in_year(configs, currTime)
        hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND + &
                   real(second, RKIND) / 3600.0_RKIND
        call update_cos_sza_diag(mesh, diag, DoY, hour_utc, days_in_year, nCells)

        nullify(tuvx_filepath_ptr)
        call mpas_pool_get_config(configs, 'config_tuvx_config_file', tuvx_filepath_ptr)
        if (associated(tuvx_filepath_ptr)) then
        if (len_trim(tuvx_filepath_ptr) > 0) then
            tuvx_top_extension = .false.
            tuvx_extension_file = ''
            nullify(tuvx_top_extension_ptr)
            nullify(tuvx_extension_file_ptr)
            call mpas_pool_get_config(configs, 'config_tuvx_top_extension', tuvx_top_extension_ptr)
            call mpas_pool_get_config(configs, 'config_tuvx_extension_file', tuvx_extension_file_ptr)
            if (associated(tuvx_top_extension_ptr)) tuvx_top_extension = tuvx_top_extension_ptr
            if (associated(tuvx_extension_file_ptr)) tuvx_extension_file = tuvx_extension_file_ptr

            call tuvx_init(tuvx_filepath_ptr, nVertLevels, tuvx_top_extension, &
                           tuvx_extension_file, error_code, error_message)
            if (error_code /= 0) return
            use_tuvx = tuvx_is_enabled()
        end if
        end if

        if (use_tuvx .and. micm_mode) then
            block
                character(len=64), allocatable :: photo_names(:)
                integer :: n_rates
                n_rates = tuvx_n_rates()
                allocate(photo_names(n_rates))
                call tuvx_rate_names(photo_names)
                call musica_cache_photo_indices(photo_names, error_code, error_message)
                deallocate(photo_names)
            end block
            if (error_code /= 0) return
        else if (use_tuvx .and. .not. micm_mode) then
            call mpas_log_write('[Chemistry] Photolysis-only mode active: TUV-x enabled, MICM disabled.')
        end if
#else
        error_code = 0
        error_message = ''
#endif

    end subroutine photolysis_init


    !------------------------------------------------------------------------
    !  function photolysis_enabled
    !------------------------------------------------------------------------
    logical function photolysis_enabled() result(enabled)

        enabled = use_tuvx

    end function photolysis_enabled


    !------------------------------------------------------------------------
    !  function photolysis_update_due
    !------------------------------------------------------------------------
    logical function photolysis_update_due() result(update_now)

        update_now = use_tuvx .and. (tuvx_time_since_last >= tuvx_update_interval)

    end function photolysis_update_due


    !------------------------------------------------------------------------
    !  routine photolysis_advance_time
    !------------------------------------------------------------------------
    subroutine photolysis_advance_time(dt, updated_now)

        real(kind=RKIND), intent(in) :: dt
        logical, intent(in) :: updated_now

        if (.not. use_tuvx) return

        if (updated_now) then
            tuvx_time_since_last = 0.0_RKIND
        else
            tuvx_time_since_last = tuvx_time_since_last + dt
        end if

    end subroutine photolysis_advance_time


    !------------------------------------------------------------------------
    !  routine photolysis_step
    !------------------------------------------------------------------------
    subroutine photolysis_step(state, mesh, diag, scalars, rho, qv_chem, temperature_chem, &
                               theta_m, exner, zgrid, nCells, nVertLevels, idx_qv_val, &
                               idx_o3_val, idx_o2_val, idx_n2_val, idx_h2o_val, DoY, &
                               hour_utc, days_in_year, update_photolysis, micm_mode, &
                               error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_set_photolysis_rates, n_photo_rp
        use mpas_tuvx, only: tuvx_compute_photolysis, tuvx_n_rates
#endif
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_log, only : mpas_log_write
        use mpas_pool_routines, only: mpas_pool_get_array, mpas_pool_get_dimension
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        type(mpas_pool_type), intent(in) :: state
        type(mpas_pool_type), pointer, intent(in) :: mesh
        type(mpas_pool_type), pointer, intent(inout) :: diag
        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        real(kind=RKIND), dimension(:,:), intent(in) :: rho
        real(kind=RKIND), dimension(:,:), intent(in) :: qv_chem
        real(kind=RKIND), dimension(:,:), intent(in) :: temperature_chem
        real(kind=RKIND), dimension(:,:), pointer, intent(in) :: theta_m
        real(kind=RKIND), dimension(:,:), pointer, intent(in) :: exner
        real(kind=RKIND), dimension(:,:), pointer, intent(in) :: zgrid
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(in) :: idx_qv_val
        integer, intent(in) :: idx_o3_val
        integer, intent(in) :: idx_o2_val
        integer, intent(in) :: idx_n2_val
        integer, intent(in) :: idx_h2o_val
        integer, intent(in) :: DoY
        real(kind=RKIND), intent(in) :: hour_utc
        integer, intent(in) :: days_in_year
        logical, intent(in) :: update_photolysis
        logical, intent(in) :: micm_mode
        integer, intent(out) :: error_code
        character(len=:), allocatable, intent(out) :: error_message

#ifdef MPAS_USE_MUSICA
        real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
        integer, pointer :: index_qc, index_qr
        integer :: n_photo_rates_local
        integer :: iCell, k
        integer :: idx_qc_val, idx_qr_val
        integer :: n_bad_photo, n_clamped_photo
        integer :: n_bad_photo_cell, n_clamped_photo_cell
        real(kind=RKIND) :: cos_sza, qv_loc
        real(kind=RKIND), parameter :: RAD2DEG = 57.29577951308232_RKIND
        real(kind=RKIND), parameter :: rvord = 461.5_RKIND / 287.0_RKIND
        real(kind=RKIND), allocatable :: photo_rates(:,:,:)
        real(kind=RKIND), allocatable :: rho_air_col(:), temperature_col(:)
        real(kind=RKIND), allocatable :: qc_col(:), qr_col(:), o3_col(:), o2_col(:), n2_col(:), h2o_col(:)
#endif

        error_code = 0
        error_message = ''

#ifdef MPAS_USE_MUSICA
        call update_cos_sza_diag(mesh, diag, DoY, hour_utc, days_in_year, nCells)

        if (.not. update_photolysis) return

        n_photo_rates_local = 0
        if (use_tuvx) n_photo_rates_local = tuvx_n_rates()
        if (n_photo_rates_local <= 0) return

        allocate(photo_rates(n_photo_rates_local, nVertLevels, nCells))
        photo_rates = 0.0_RKIND
        n_bad_photo = 0
        n_clamped_photo = 0

        nullify(index_qc)
        nullify(index_qr)
        call mpas_pool_get_dimension(state, 'index_qc', index_qc)
        call mpas_pool_get_dimension(state, 'index_qr', index_qr)
        idx_qc_val = -1
        idx_qr_val = -1
        if (associated(index_qc)) idx_qc_val = index_qc
        if (associated(index_qr)) idx_qr_val = index_qr

        allocate(rho_air_col(nVertLevels))
        allocate(temperature_col(nVertLevels))
        allocate(qc_col(nVertLevels))
        allocate(qr_col(nVertLevels))
        allocate(o3_col(nVertLevels))
        allocate(o2_col(nVertLevels))
        allocate(n2_col(nVertLevels))
        allocate(h2o_col(nVertLevels))

        nullify(latCell)
        nullify(lonCell)
        call mpas_pool_get_array(mesh, 'latCell', latCell)
        call mpas_pool_get_array(mesh, 'lonCell', lonCell)
        if (.not. associated(latCell) .or. .not. associated(lonCell)) then
            error_code = 1
            error_message = '[TUV-x] Could not resolve latCell/lonCell for photolysis.'
            deallocate(photo_rates, rho_air_col, temperature_col, qc_col, qr_col, o3_col, o2_col, n2_col, h2o_col)
            return
        end if

        do iCell = 1, nCells
            cos_sza = solar_cos_sza_safe(DoY, hour_utc, latCell(iCell) * RAD2DEG, lonCell(iCell) * RAD2DEG, days_in_year)

            do k = 1, nVertLevels
                qv_loc = 0.0_RKIND
                if (idx_qv_val > 0) then
                    qv_loc = scalars(idx_qv_val, k, iCell)
                    if (.not. ieee_is_finite(qv_loc) .or. qv_loc < 0.0_RKIND) qv_loc = 0.0_RKIND
                end if

                if (idx_qv_val > 0) then
                    rho_air_col(k) = rho(k, iCell) * (1.0_RKIND + qv_loc)
                    temperature_col(k) = (theta_m(k, iCell) / (1.0_RKIND + rvord * qv_loc)) * exner(k, iCell)
                else
                    rho_air_col(k) = rho(k, iCell)
                    temperature_col(k) = theta_m(k, iCell) * exner(k, iCell)
                end if

                if (.not. ieee_is_finite(rho_air_col(k)) .or. rho_air_col(k) <= 0.0_RKIND) then
                    rho_air_col(k) = max(1.0e-12_RKIND, rho(k, iCell))
                end if
                if (.not. ieee_is_finite(temperature_col(k)) .or. temperature_col(k) <= 0.0_RKIND) then
                    temperature_col(k) = temperature_chem(k, iCell)
                end if
            end do

            if (idx_qc_val > 0) then
                qc_col(:) = scalars(idx_qc_val, :, iCell)
            else
                qc_col(:) = 0.0_RKIND
            end if
            if (idx_qr_val > 0) then
                qr_col(:) = scalars(idx_qr_val, :, iCell)
            else
                qr_col(:) = 0.0_RKIND
            end if
            if (idx_o3_val > 0) then
                o3_col(:) = scalars(idx_o3_val, :, iCell)
            else
                o3_col(:) = 0.0_RKIND
            end if
            if (idx_o2_val > 0) then
                o2_col(:) = scalars(idx_o2_val, :, iCell)
            else
                o2_col(:) = 0.0_RKIND
            end if
            if (idx_n2_val > 0) then
                n2_col(:) = scalars(idx_n2_val, :, iCell)
            else
                n2_col(:) = 0.0_RKIND
            end if
            h2o_col(:) = qv_chem(:, iCell)
            if (idx_h2o_val > 0) then
                do k = 1, nVertLevels
                    if (scalars(idx_h2o_val, k, iCell) > 0.0_RKIND) h2o_col(k) = scalars(idx_h2o_val, k, iCell)
                end do
            end if

            call tuvx_compute_photolysis(cos_sza, zgrid(1:nVertLevels+1, iCell), rho_air_col, &
                                         temperature_col, o3_col, o2_col, n2_col, h2o_col, &
                                         qc_col, qr_col, nVertLevels, &
                                         photo_rates(:, :, iCell), error_code, error_message)
            if (error_code /= 0) then
                deallocate(photo_rates, rho_air_col, temperature_col, qc_col, qr_col, o3_col, o2_col, &
                           n2_col, h2o_col)
                return
            end if
            call sanitize_photo_rates(photo_rates(:, :, iCell), n_bad_photo_cell, n_clamped_photo_cell)
            n_bad_photo = n_bad_photo + n_bad_photo_cell
            n_clamped_photo = n_clamped_photo + n_clamped_photo_cell
        end do

        deallocate(rho_air_col, temperature_col, qc_col, qr_col, o3_col, o2_col, n2_col, h2o_col)
        if (n_bad_photo > 0 .or. n_clamped_photo > 0) then
            call mpas_log_write('[TUV-x] Sanitized photolysis rates: replaced_nonfinite_or_negative=$i clamped_high=$i', &
                                intArgs=[n_bad_photo, n_clamped_photo])
        end if

        if (micm_mode) then
            call musica_set_photolysis_rates(photo_rates, nCells, nVertLevels, n_photo_rp, error_code, error_message)
            if (error_code /= 0) then
                deallocate(photo_rates)
                return
            end if
        end if

        call set_photolysis_diag(diag, photo_rates, nCells, nVertLevels, micm_mode)
        deallocate(photo_rates)
#endif

    end subroutine photolysis_step


    !------------------------------------------------------------------------
    !  routine update_cos_sza_diag
    !------------------------------------------------------------------------
    subroutine update_cos_sza_diag(mesh, diag, DoY, hour_utc, days_in_year, nCells)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_array

        type(mpas_pool_type), pointer, intent(in) :: mesh
        type(mpas_pool_type), pointer, intent(inout) :: diag
        integer, intent(in) :: DoY
        real(kind=RKIND), intent(in) :: hour_utc
        integer, intent(in) :: days_in_year
        integer, intent(in) :: nCells

        real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
        real(kind=RKIND), dimension(:), pointer :: diag_cos_sza
        integer :: iCell
        real(kind=RKIND), parameter :: RAD2DEG = 57.29577951308232_RKIND

        if (.not. associated(mesh)) return
        if (.not. associated(diag)) return

        nullify(latCell)
        nullify(lonCell)
        nullify(diag_cos_sza)
        call mpas_pool_get_array(mesh, 'latCell', latCell)
        call mpas_pool_get_array(mesh, 'lonCell', lonCell)
        call mpas_pool_get_array(diag, 'cos_sza', diag_cos_sza)

        if (.not. associated(latCell)) return
        if (.not. associated(lonCell)) return
        if (.not. associated(diag_cos_sza)) return

        do iCell = 1, min(nCells, size(diag_cos_sza), size(latCell), size(lonCell))
            diag_cos_sza(iCell) = solar_cos_sza_safe(DoY, hour_utc, latCell(iCell) * RAD2DEG, &
                                                     lonCell(iCell) * RAD2DEG, days_in_year)
        end do

    end subroutine update_cos_sza_diag


    !------------------------------------------------------------------------
    !  function solar_cos_sza_safe
    !------------------------------------------------------------------------
    real(kind=RKIND) function solar_cos_sza_safe(DoY, hour_utc, lat_deg, lon_deg, days_in_year) result(cos_sza)

        use mpas_solar_geometry, only: solar_cos_sza
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        integer, intent(in) :: DoY
        real(kind=RKIND), intent(in) :: hour_utc
        real(kind=RKIND), intent(in) :: lat_deg
        real(kind=RKIND), intent(in) :: lon_deg
        integer, intent(in) :: days_in_year

        cos_sza = solar_cos_sza(DoY, hour_utc, lat_deg, lon_deg, days_in_year)
        if (.not. ieee_is_finite(cos_sza)) then
            cos_sza = 0.0_RKIND
        else
            cos_sza = max(-1.0_RKIND, min(1.0_RKIND, cos_sza))
        end if

    end function solar_cos_sza_safe


    !------------------------------------------------------------------------
    !  routine set_photolysis_diag
    !------------------------------------------------------------------------
    subroutine set_photolysis_diag(diag, photo_rates, nCells, nVertLevels, micm_mode)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: photo_rp_names, n_photo_rp
        use mpas_tuvx, only: tuvx_n_rates, tuvx_rate_names
#endif
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_log, only : mpas_log_write
        use mpas_pool_routines, only: mpas_pool_get_array

        type(mpas_pool_type), intent(inout), pointer :: diag
        real(kind=RKIND), intent(in) :: photo_rates(:,:,:)
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        logical, intent(in) :: micm_mode

#ifdef MPAS_USE_MUSICA
        real(kind=RKIND), dimension(:,:), pointer :: diag_field
        integer :: ir, iCell, k, n_rates_diag
        character(len=64), allocatable :: rate_names_diag(:)
        character(len=64) :: rp_name
        integer :: idx_jno2, idx_jo2_a, idx_jo2_b, idx_jo3_a, idx_jo3_b
        character(len=256) :: missing_rates

        if (.not. associated(diag)) return

        if (micm_mode) then
            n_rates_diag = n_photo_rp
            if (n_rates_diag <= 0) return
            allocate(rate_names_diag(n_rates_diag))
            rate_names_diag(:) = photo_rp_names(1:n_rates_diag)
        else
            n_rates_diag = tuvx_n_rates()
            if (n_rates_diag <= 0) return
            allocate(rate_names_diag(n_rates_diag))
            call tuvx_rate_names(rate_names_diag)
        end if

        do ir = 1, min(n_rates_diag, size(photo_rates, 1))
            nullify(diag_field)
            call mpas_pool_get_array(diag, 'j_' // trim(rate_names_diag(ir)), diag_field)
            if (.not. associated(diag_field)) cycle

            do iCell = 1, min(nCells, size(diag_field, 2))
                do k = 1, min(nVertLevels, size(diag_field, 1))
                    diag_field(k, iCell) = photo_rates(ir, k, iCell)
                end do
            end do
        end do

        idx_jno2 = 0
        idx_jo2_a = 0
        idx_jo2_b = 0
        idx_jo3_a = 0
        idx_jo3_b = 0

        do ir = 1, min(n_rates_diag, size(photo_rates, 1))
            rp_name = trim(rate_names_diag(ir))
            select case (rp_name)
            case ('jno2', 'JNO2')
                idx_jno2 = ir
            case ('jo2_a', 'JO2_A')
                idx_jo2_a = ir
            case ('jo2_b', 'JO2_B')
                idx_jo2_b = ir
            case ('jo3_a', 'JO3_A')
                idx_jo3_a = ir
            case ('jo3_b', 'JO3_B')
                idx_jo3_b = ir
            end select
        end do

        missing_rates = ''
        if (idx_jno2 == 0) missing_rates = trim(missing_rates) // ' jno2'
        if (idx_jo2_a == 0) missing_rates = trim(missing_rates) // ' jo2_a'
        if (idx_jo2_b == 0) missing_rates = trim(missing_rates) // ' jo2_b'
        if (idx_jo3_a == 0) missing_rates = trim(missing_rates) // ' jo3_a'
        if (idx_jo3_b == 0) missing_rates = trim(missing_rates) // ' jo3_b'
        if (len_trim(missing_rates) > 0) then
            call mpas_log_write('[TUV-x] Missing required photolysis rate(s) for diagnostics:' // trim(missing_rates), &
                                messageType=MPAS_LOG_CRIT)
            deallocate(rate_names_diag)
            return
        end if

        nullify(diag_field)
        call mpas_pool_get_array(diag, 'j_NO2', diag_field)
        if (associated(diag_field)) then
            do iCell = 1, min(nCells, size(diag_field, 2))
                do k = 1, min(nVertLevels, size(diag_field, 1))
                    diag_field(k, iCell) = photo_rates(idx_jno2, k, iCell)
                end do
            end do
        end if

        nullify(diag_field)
        call mpas_pool_get_array(diag, 'j_O2', diag_field)
        if (associated(diag_field)) then
            do iCell = 1, min(nCells, size(diag_field, 2))
                do k = 1, min(nVertLevels, size(diag_field, 1))
                    diag_field(k, iCell) = photo_rates(idx_jo2_a, k, iCell) + photo_rates(idx_jo2_b, k, iCell)
                end do
            end do
        end if

        nullify(diag_field)
        call mpas_pool_get_array(diag, 'j_O3_O', diag_field)
        if (associated(diag_field)) then
            do iCell = 1, min(nCells, size(diag_field, 2))
                do k = 1, min(nVertLevels, size(diag_field, 1))
                    diag_field(k, iCell) = photo_rates(idx_jo3_a, k, iCell)
                end do
            end do
        end if

        nullify(diag_field)
        call mpas_pool_get_array(diag, 'j_O3_O1D', diag_field)
        if (associated(diag_field)) then
            do iCell = 1, min(nCells, size(diag_field, 2))
                do k = 1, min(nVertLevels, size(diag_field, 1))
                    diag_field(k, iCell) = photo_rates(idx_jo3_b, k, iCell)
                end do
            end do
        end if

        deallocate(rate_names_diag)
#endif

    end subroutine set_photolysis_diag


    !------------------------------------------------------------------------
    !  routine zero_photolysis_diag
    !------------------------------------------------------------------------
    subroutine zero_photolysis_diag(diag, names)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_array

        type(mpas_pool_type), intent(inout), pointer :: diag
        character(len=*), intent(in) :: names(:)

        real(kind=RKIND), dimension(:,:), pointer :: diag_field
        integer :: i

        if (.not. associated(diag)) return

        do i = 1, size(names)
            nullify(diag_field)
            call mpas_pool_get_array(diag, trim(names(i)), diag_field)
            if (associated(diag_field)) diag_field = 0.0_RKIND
        end do

    end subroutine zero_photolysis_diag


    !------------------------------------------------------------------------
    !  routine sanitize_photo_rates
    !------------------------------------------------------------------------
    subroutine sanitize_photo_rates(rates, n_bad, n_clamped)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(kind=RKIND), dimension(:,:), intent(inout) :: rates
        integer, intent(out) :: n_bad, n_clamped
        integer :: ir, k
        real(kind=RKIND), parameter :: photo_rate_max = 1.0_RKIND

        n_bad = 0
        n_clamped = 0
        do ir = 1, size(rates, 1)
            do k = 1, size(rates, 2)
                if (.not. ieee_is_finite(rates(ir, k)) .or. rates(ir, k) < 0.0_RKIND) then
                    rates(ir, k) = 0.0_RKIND
                    n_bad = n_bad + 1
                else if (rates(ir, k) > photo_rate_max) then
                    rates(ir, k) = photo_rate_max
                    n_clamped = n_clamped + 1
                end if
            end do
        end do

    end subroutine sanitize_photo_rates


    !------------------------------------------------------------------------
    !  routine photolysis_finalize
    !------------------------------------------------------------------------
    subroutine photolysis_finalize()

#ifdef MPAS_USE_MUSICA
        use mpas_tuvx, only: tuvx_finalize

        call tuvx_finalize()
#endif
        use_tuvx = .false.
        tuvx_time_since_last = huge(1.0_RKIND)

    end subroutine photolysis_finalize

end module mpas_chem_photolysis_driver
