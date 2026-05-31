! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_atm_chemistry
!
!> \brief Top-level MPAS atmosphere chemistry orchestrator
!> \details
!>  Coordinates MPAS-side chemistry work without owning low-level details:
!>  initializes MUSICA/MICM, dynamic emission tracers, and optional TUV-x;
!>  prepares MPAS thermodynamic/environment fields; applies chemistry
!>  emissions; advances photolysis diagnostics/rates; and steps MICM.
!-----------------------------------------------------------------------
module mpas_atm_chemistry

    use mpas_kind_types, only: RKIND, StrKIND

    implicit none

    private

    public :: chemistry_init, chemistry_step, chemistry_finalize
    public :: chemistry_query_species

    logical, save :: micm_enabled = .false.
    integer, save :: chem_substeps = 1

contains

    !------------------------------------------------------------------------
    !  routine chemistry_init
    !------------------------------------------------------------------------
    subroutine chemistry_init(configs, structs, clock)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_init, resolve_mpas_indices
        use mpas_chem_emissions_driver, only: emissions_init
        use mpas_chem_photolysis_driver, only: photolysis_init
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: mpas_pool_type, MPAS_LOG_CRIT
        use mpas_pool_routines, only: mpas_pool_get_config, mpas_pool_get_subpool, mpas_pool_get_dimension
        use mpas_timekeeping, only: MPAS_Clock_type

        type (mpas_pool_type), intent(in) :: configs
        type (mpas_pool_type), intent(in) :: structs
        type (MPAS_Clock_type), intent(in) :: clock

#ifdef MPAS_USE_MUSICA
        character(len=StrKIND), pointer :: filepath_ptr
        character(len=StrKIND), pointer :: emission_list_ptr
        character(len=StrKIND)          :: micm_file
        character(len=StrKIND)          :: emission_list
        integer                         :: error_code
        character(len=:), allocatable   :: error_message
        integer                         :: nCells
        integer                         :: nVertLevels
        integer                         :: nChemGridCells
        type(mpas_pool_type), pointer   :: state
        integer, pointer                :: nCells_ptr, nVertLevels_ptr
        integer, pointer                :: int_ptr
#endif

        call mpas_log_write('Initializing chemistry packages...')

#ifdef MPAS_USE_MUSICA
        call mpas_pool_get_subpool(structs, 'state', state)
        nullify(nCells_ptr)
        nullify(nVertLevels_ptr)
        call mpas_pool_get_dimension(state, 'nCells', nCells_ptr)
        call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels_ptr)
        if (.not. associated(nCells_ptr) .or. .not. associated(nVertLevels_ptr)) then
            call mpas_log_write('[Chemistry] Could not resolve nCells/nVertLevels in state pool.', &
                                messageType=MPAS_LOG_CRIT)
            return
        end if

        nCells = nCells_ptr
        nVertLevels = nVertLevels_ptr
        nChemGridCells = nCells * nVertLevels

        nullify(filepath_ptr)
        call mpas_pool_get_config(configs, 'config_micm_file', filepath_ptr)
        micm_file = ''
        if (associated(filepath_ptr)) micm_file = filepath_ptr
        micm_enabled = (len_trim(micm_file) /= 0)

        nullify(emission_list_ptr)
        call mpas_pool_get_config(configs, 'config_chem_emission_species', emission_list_ptr)
        emission_list = ''
        if (associated(emission_list_ptr)) emission_list = emission_list_ptr

        if (micm_enabled) then
            call musica_init(micm_file, nChemGridCells, error_code, error_message)
            if (error_code /= 0) then
                call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                return
            end if

            call resolve_mpas_indices(state, error_code, error_message)
            if (error_code /= 0) then
                call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                return
            end if
        else
            call mpas_log_write('[Chemistry] MICM disabled (config_micm_file is empty).')
        end if

        call emissions_init(micm_file, emission_list, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            return
        end if

        nullify(int_ptr)
        call mpas_pool_get_config(configs, 'config_chem_substeps', int_ptr)
        if (associated(int_ptr)) chem_substeps = max(1, int_ptr)

        call photolysis_init(configs, structs, clock, nVertLevels, micm_enabled, error_code, error_message)
        if (error_code /= 0) then
            call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
            return
        end if
#endif

    end subroutine chemistry_init


    !------------------------------------------------------------------------
    !  routine chemistry_step
    !------------------------------------------------------------------------
    subroutine chemistry_step(domain, dt)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_step, musica_mpas_index_for_species
        use mpas_chem_emissions_driver, only: emissions_count, emissions_prepare, emissions_apply
        use mpas_chem_photolysis_driver, only: photolysis_update_due, photolysis_step, photolysis_advance_time
        use mpas_chem_utils, only: resolve_scalar_index, chemistry_days_in_year
#endif
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only: domain_type, block_type, mpas_pool_type, MPAS_LOG_CRIT, MPAS_NOW, MPAS_Time_type
        use mpas_pool_routines, only: mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension, mpas_pool_get_config
        use mpas_timekeeping, only: mpas_get_clock_time, mpas_get_time
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        type(domain_type), intent(inout) :: domain
        real(kind=RKIND), intent(in)     :: dt

        type(block_type), pointer                   :: block_ptr
        type(mpas_pool_type), pointer               :: state, mesh, diag
        real(kind=RKIND), dimension(:,:,:), pointer :: scalars
        real(kind=RKIND), dimension(:,:,:), pointer :: chem_emissions
        real(kind=RKIND), dimension(:,:), pointer   :: rho_zz, zz

#ifdef MPAS_USE_MUSICA
        integer                         :: error_code
        character(len=:), allocatable   :: error_message
        character(len=StrKIND), pointer :: filepath_ptr
        integer                         :: iCell, k
        integer                         :: nCells_local, nVertLevels_local
        real(kind=RKIND), parameter     :: aircol_eps = 1.0e-30_RKIND
        real(kind=RKIND), parameter     :: p0_ref = 100000.0_RKIND
        real(kind=RKIND), parameter     :: cp_over_rd = 3.5_RKIND
        real(kind=RKIND), parameter     :: rvord = 461.5_RKIND / 287.0_RKIND
        real(kind=RKIND), dimension(:,:), allocatable :: inv_aircol
        real(kind=RKIND), dimension(:,:), allocatable :: rho
        real(kind=RKIND), dimension(:,:), allocatable :: temperature_chem, pressure_chem, qv_chem
        real(kind=RKIND), dimension(:,:), pointer :: zgrid
        real(kind=RKIND)                :: dz, aircol, qv_loc
        real(kind=RKIND), dimension(:,:), pointer :: theta_m, exner
        integer, pointer :: index_qv
        integer, pointer :: nCells_ptr, nVertLevels_ptr
        type(MPAS_Time_type) :: currTime
        integer :: ierr, DoY, hour, minute, second, days_in_year
        integer :: i_sub
        integer :: idx_o3_val, idx_o2_val, idx_n2_val, idx_h2o_val
        integer :: idx_qv_val
        real(kind=RKIND) :: hour_utc, sub_dt
        logical :: update_photolysis, micm_enabled_local, apply_emissions
        integer, allocatable :: emission_indices(:)
#endif

#ifdef MPAS_USE_MUSICA
        call mpas_log_write('Stepping chemistry packages...')

        block_ptr => domain % blocklist
        if (.not. associated(block_ptr)) return

        nullify(filepath_ptr)
        call mpas_pool_get_config(domain % blocklist % configs, 'config_micm_file', filepath_ptr)
        micm_enabled_local = micm_enabled
        if (associated(filepath_ptr)) then
            micm_enabled_local = micm_enabled_local .and. (len_trim(filepath_ptr) /= 0)
        else
            micm_enabled_local = .false.
        end if
        update_photolysis = photolysis_update_due()

        currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
        call mpas_get_time(currTime, DoY=DoY, H=hour, M=minute, S=second)
        days_in_year = chemistry_days_in_year(domain % blocklist % configs, currTime)
        hour_utc = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND + &
                   real(second, RKIND) / 3600.0_RKIND

        do while (associated(block_ptr))
            nullify(state)
            nullify(mesh)
            nullify(diag)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)

            call mpas_pool_get_array(state, 'scalars', scalars, 1)
            call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
            call mpas_pool_get_array(mesh, 'zz', zz)
            call mpas_pool_get_array(mesh, 'zgrid', zgrid)
            call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
            call mpas_pool_get_array(diag, 'exner', exner)

            nullify(chem_emissions)
            if (emissions_count() > 0) then
                call mpas_pool_get_array(state, 'chem_emissions', chem_emissions, 1)
            end if

            nullify(nCells_ptr)
            nullify(nVertLevels_ptr)
            call mpas_pool_get_dimension(mesh, 'nCells', nCells_ptr)
            call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels_ptr)
            if (.not. associated(nCells_ptr) .or. .not. associated(nVertLevels_ptr)) then
                call mpas_log_write('[Chemistry] Could not resolve nCells/nVertLevels in mesh pool.', &
                                    messageType=MPAS_LOG_CRIT)
                return
            end if
            nCells_local = nCells_ptr
            nVertLevels_local = nVertLevels_ptr

            nullify(index_qv)
            call mpas_pool_get_dimension(state, 'index_qv', index_qv)
            idx_qv_val = -1
            if (associated(index_qv)) idx_qv_val = index_qv

            idx_o3_val = -1
            if (micm_enabled_local) idx_o3_val = musica_mpas_index_for_species('O3')
            if (idx_o3_val <= 0) idx_o3_val = resolve_scalar_index(state, 'qO3')
            if (idx_o3_val <= 0) idx_o3_val = resolve_scalar_index(state, 'qo3')
            if (idx_o3_val <= 0) idx_o3_val = resolve_scalar_index(state, 'o3')

            idx_o2_val = -1
            if (micm_enabled_local) idx_o2_val = musica_mpas_index_for_species('O2')
            if (idx_o2_val <= 0) idx_o2_val = resolve_scalar_index(state, 'qO2')
            if (idx_o2_val <= 0) idx_o2_val = resolve_scalar_index(state, 'qo2')

            idx_n2_val = -1
            if (micm_enabled_local) idx_n2_val = musica_mpas_index_for_species('N2')
            if (idx_n2_val <= 0) idx_n2_val = resolve_scalar_index(state, 'qN2')
            if (idx_n2_val <= 0) idx_n2_val = resolve_scalar_index(state, 'qn2')

            idx_h2o_val = -1
            if (micm_enabled_local) idx_h2o_val = musica_mpas_index_for_species('H2O')
            if (idx_h2o_val <= 0) idx_h2o_val = resolve_scalar_index(state, 'qH2O')
            if (idx_h2o_val <= 0) idx_h2o_val = resolve_scalar_index(state, 'qh2o')

            call emissions_prepare(state, scalars, chem_emissions, emission_indices, apply_emissions)

            if (micm_enabled_local .or. apply_emissions) then
                allocate(inv_aircol(nVertLevels_local, nCells_local))
            end if

            allocate(rho(nVertLevels_local, nCells_local))
            allocate(temperature_chem(nVertLevels_local, nCells_local))
            allocate(pressure_chem(nVertLevels_local, nCells_local))
            allocate(qv_chem(nVertLevels_local, nCells_local))

            do iCell = 1, nCells_local
                do k = 1, nVertLevels_local
                    rho(k, iCell) = rho_zz(k, iCell) * zz(k, iCell)
                    if (micm_enabled_local .or. apply_emissions) then
                        dz = zgrid(k+1, iCell) - zgrid(k, iCell)
                        if (dz < 0.0_RKIND) dz = -dz
                        aircol = rho(k, iCell) * dz
                        if (aircol > aircol_eps) then
                            inv_aircol(k, iCell) = 1.0_RKIND / aircol
                        else
                            inv_aircol(k, iCell) = 0.0_RKIND
                        end if
                    end if

                    qv_loc = 0.0_RKIND
                    if (idx_qv_val > 0) then
                        qv_loc = scalars(idx_qv_val, k, iCell)
                        if (.not. ieee_is_finite(qv_loc) .or. qv_loc < 0.0_RKIND) qv_loc = 0.0_RKIND
                    end if
                    qv_chem(k, iCell) = qv_loc

                    if (idx_qv_val > 0) then
                        temperature_chem(k, iCell) = (theta_m(k, iCell) / (1.0_RKIND + rvord * qv_loc)) * exner(k, iCell)
                    else
                        temperature_chem(k, iCell) = theta_m(k, iCell) * exner(k, iCell)
                    end if
                    if (.not. ieee_is_finite(temperature_chem(k, iCell)) .or. temperature_chem(k, iCell) <= 0.0_RKIND) then
                        temperature_chem(k, iCell) = 280.0_RKIND
                    end if

                    if (ieee_is_finite(exner(k, iCell)) .and. exner(k, iCell) > 0.0_RKIND) then
                        pressure_chem(k, iCell) = p0_ref * exner(k, iCell) ** cp_over_rd
                    else
                        pressure_chem(k, iCell) = p0_ref
                    end if
                end do
            end do

            if (apply_emissions) then
                call emissions_apply(scalars, chem_emissions, inv_aircol, emission_indices, dt, nCells_local)
            end if

            call photolysis_step(state, mesh, diag, scalars, rho, qv_chem, temperature_chem, &
                                 theta_m, exner, zgrid, nCells_local, nVertLevels_local, idx_qv_val, &
                                 idx_o3_val, idx_o2_val, idx_n2_val, idx_h2o_val, DoY, &
                                 hour_utc, days_in_year, update_photolysis, micm_enabled_local, &
                                 error_code, error_message)
            if (error_code /= 0) then
                if (allocated(inv_aircol)) deallocate(inv_aircol)
                if (allocated(emission_indices)) deallocate(emission_indices)
                deallocate(rho, temperature_chem, pressure_chem, qv_chem)
                call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                return
            end if

            if (micm_enabled_local) then
                sub_dt = dt / real(chem_substeps, RKIND)
                do i_sub = 1, chem_substeps
                    call musica_step( &
                        scalars(:,1:nVertLevels_local,1:nCells_local), &
                        qv_chem, &
                        rho_zz(1:nVertLevels_local,1:nCells_local), &
                        zz(1:nVertLevels_local,1:nCells_local), &
                        temperature_chem, pressure_chem, sub_dt, error_code, error_message &
                    )

                    if (error_code /= 0) then
                        if (allocated(inv_aircol)) deallocate(inv_aircol)
                        if (allocated(emission_indices)) deallocate(emission_indices)
                        deallocate(rho, temperature_chem, pressure_chem, qv_chem)
                        call mpas_log_write(error_message, messageType=MPAS_LOG_CRIT)
                        return
                    end if
                end do
            end if

            if (allocated(inv_aircol)) deallocate(inv_aircol)
            if (allocated(emission_indices)) deallocate(emission_indices)
            deallocate(rho, temperature_chem, pressure_chem, qv_chem)

            block_ptr => block_ptr % next
        end do

        call photolysis_advance_time(dt, update_photolysis)
#endif

    end subroutine chemistry_step


    !------------------------------------------------------------------------
    !  routine chemistry_query_species
    !------------------------------------------------------------------------
    subroutine chemistry_query_species(config_micm_file, config_emission_species, &
                                       num_species, species_names, &
                                       num_emissions, emission_names, &
                                       error_code, error_message)

        use mpas_chem_emissions_driver, only: emissions_query_species

        character(len=*), intent(in)                 :: config_micm_file
        character(len=*), intent(in)                 :: config_emission_species
        integer, intent(out)                         :: num_species
        character(len=64), allocatable, intent(out)  :: species_names(:)
        integer, intent(out)                         :: num_emissions
        character(len=64), allocatable, intent(out)  :: emission_names(:)
        integer, intent(out)                         :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        call emissions_query_species(config_micm_file, config_emission_species, &
                                     num_species, species_names, num_emissions, emission_names, &
                                     error_code, error_message)

    end subroutine chemistry_query_species


    !------------------------------------------------------------------------
    !  routine chemistry_finalize
    !------------------------------------------------------------------------
    subroutine chemistry_finalize()

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_finalize
        use mpas_chem_emissions_driver, only: emissions_finalize
        use mpas_chem_photolysis_driver, only: photolysis_finalize
#endif
        use mpas_log, only : mpas_log_write

#ifdef MPAS_USE_MUSICA
        call mpas_log_write('Finalizing chemistry packages...')
        call photolysis_finalize()
        call musica_finalize()
        call emissions_finalize()
#endif

    end subroutine chemistry_finalize

end module mpas_atm_chemistry
