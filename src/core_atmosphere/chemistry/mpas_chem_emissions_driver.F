! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_chem_emissions_driver
!
!> \brief Applies chemistry emissions to MPAS chemistry tracers
!> \details
!>  Owns the runtime emission species table, maps input fields named
!>  emi_<species> to their MPAS scalar tracers, sanitizes emission fluxes,
!>  and converts fluxes from kg m-2 s-1 into tracer increments in kg kg-1
!>  using the model layer air column mass.
!-----------------------------------------------------------------------
module mpas_chem_emissions_driver

    use mpas_kind_types, only: RKIND

    implicit none

    private

    public :: emissions_query_species
    public :: emissions_init
    public :: emissions_count
    public :: emissions_prepare
    public :: emissions_apply
    public :: emissions_finalize

    integer, save :: n_emission_species = 0
    character(len=64), allocatable, save :: emission_species(:)
    character(len=64), allocatable, save :: emission_tracers(:)

contains

    !------------------------------------------------------------------------
    !  routine emissions_query_species
    !------------------------------------------------------------------------
    subroutine emissions_query_species(config_micm_file, config_emission_species, config_additional_species, &
                                       num_species, species_names, &
                                       num_emissions, emission_names, &
                                       error_code, error_message)

#ifdef MPAS_USE_MUSICA
        use mpas_musica, only: musica_query_species
#endif
        use mpas_chem_utils, only: parse_species_list, append_unique, target_tracer_for_emission
        use mpas_log, only : mpas_log_write

        character(len=*), intent(in)                 :: config_micm_file
        character(len=*), intent(in)                 :: config_emission_species
        character(len=*), intent(in)                 :: config_additional_species
        integer, intent(out)                         :: num_species
        character(len=64), allocatable, intent(out)  :: species_names(:)
        integer, intent(out)                         :: num_emissions
        character(len=64), allocatable, intent(out)  :: emission_names(:)
        integer, intent(out)                         :: error_code
        character(len=:), allocatable, intent(out)   :: error_message

        character(len=64), allocatable :: micm_names(:)
        character(len=64), allocatable :: additional_names(:)
        character(len=64), allocatable :: work_names(:)
        character(len=64) :: target_name
        integer :: n_micm, n_additional, i, max_names

        error_code = 0
        error_message = ''
        num_species = 0
        num_emissions = 0

#ifdef MPAS_USE_MUSICA
        n_micm = 0
        allocate(micm_names(0))
        if (len_trim(config_micm_file) > 0) then
            call musica_query_species(config_micm_file, n_micm, micm_names, error_code, error_message)
            if (error_code /= 0) return
        end if

        call parse_species_list(config_emission_species, num_emissions, emission_names)
        call parse_species_list(config_additional_species, n_additional, additional_names)

        max_names = n_micm + num_emissions + n_additional
        allocate(work_names(max(1, max_names)))

        do i = 1, n_micm
            call append_unique(work_names, num_species, trim(micm_names(i)))
        end do

        do i = 1, num_emissions
            target_name = target_tracer_for_emission(emission_names(i), micm_names, n_micm)
            call append_unique(work_names, num_species, trim(target_name))
        end do
        do i = 1, n_additional
            call append_unique(work_names, num_species, trim(additional_names(i)))
        end do

        allocate(species_names(num_species))
        if (num_species > 0) species_names(:) = work_names(1:num_species)

        if (num_species > 0) then
            call mpas_log_write('[Chemistry] Dynamic tracer count: $i', intArgs=[num_species])
        end if

        if (allocated(micm_names)) deallocate(micm_names)
        if (allocated(additional_names)) deallocate(additional_names)
        if (allocated(work_names)) deallocate(work_names)
#else
        allocate(species_names(0))
        allocate(emission_names(0))
#endif

    end subroutine emissions_query_species


    !------------------------------------------------------------------------
    !  routine emissions_init
    !------------------------------------------------------------------------
    subroutine emissions_init(config_micm_file, config_emission_species, config_additional_species, error_code, error_message)

        use mpas_chem_utils, only: target_tracer_for_emission
        use mpas_log, only : mpas_log_write

        character(len=*), intent(in)                :: config_micm_file
        character(len=*), intent(in)                :: config_emission_species
        character(len=*), intent(in)                :: config_additional_species
        integer, intent(out)                        :: error_code
        character(len=:), allocatable, intent(out)  :: error_message

        character(len=64), allocatable :: tracer_names(:)
        character(len=64), allocatable :: local_emission_names(:)
        integer :: n_tracers, i

        error_code = 0
        error_message = ''

        call emissions_finalize()

        call emissions_query_species(config_micm_file, config_emission_species, config_additional_species, &
                                     n_tracers, tracer_names, n_emission_species, local_emission_names, &
                                     error_code, error_message)
        if (error_code /= 0) return

        allocate(emission_species(n_emission_species))
        allocate(emission_tracers(n_emission_species))

        do i = 1, n_emission_species
            emission_species(i) = local_emission_names(i)
            emission_tracers(i) = target_tracer_for_emission(local_emission_names(i), tracer_names, n_tracers)
            call mpas_log_write('[Chemistry] Emission field emi_' // trim(emission_species(i)) // &
                                ' updates tracer ' // trim(emission_tracers(i)))
        end do

        if (allocated(tracer_names)) deallocate(tracer_names)
        if (allocated(local_emission_names)) deallocate(local_emission_names)

    end subroutine emissions_init


    !------------------------------------------------------------------------
    !  function emissions_count
    !------------------------------------------------------------------------
    integer function emissions_count() result(n)

        n = n_emission_species

    end function emissions_count


    !------------------------------------------------------------------------
    !  routine emissions_prepare
    !------------------------------------------------------------------------
    subroutine emissions_prepare(state, scalars, chem_emissions, emission_indices, apply_emissions)

        use mpas_chem_utils, only: resolve_scalar_index, sanitize_species_by_index
        use mpas_derived_types, only: mpas_pool_type

        type(mpas_pool_type), intent(in) :: state
        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        real(kind=RKIND), dimension(:,:,:), pointer, intent(inout) :: chem_emissions
        integer, allocatable, intent(out) :: emission_indices(:)
        logical, intent(out) :: apply_emissions

        integer :: i_emit

        apply_emissions = .false.

        if (n_emission_species <= 0) return
        if (.not. associated(chem_emissions)) return

        allocate(emission_indices(n_emission_species))
        emission_indices(:) = -1

        do i_emit = 1, n_emission_species
            emission_indices(i_emit) = resolve_scalar_index(state, trim(emission_tracers(i_emit)))
            if (emission_indices(i_emit) > 0 .and. i_emit <= size(chem_emissions, 1)) apply_emissions = .true.
        end do

        if (.not. apply_emissions) return

        do i_emit = 1, n_emission_species
            if (emission_indices(i_emit) > 0 .and. i_emit <= size(chem_emissions, 1)) then
                call sanitize_species_by_index(scalars, emission_indices(i_emit))
                call sanitize_emission_field(chem_emissions(i_emit,:,:))
            end if
        end do

    end subroutine emissions_prepare


    !------------------------------------------------------------------------
    !  routine emissions_apply
    !------------------------------------------------------------------------
    subroutine emissions_apply(scalars, chem_emissions, inv_aircol, emission_indices, dt, nCells)

        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        real(kind=RKIND), dimension(:,:,:), pointer, intent(in) :: chem_emissions
        real(kind=RKIND), dimension(:,:), intent(in) :: inv_aircol
        integer, intent(in) :: emission_indices(:)
        real(kind=RKIND), intent(in) :: dt
        integer, intent(in) :: nCells

        integer :: i_emit

        if (.not. associated(chem_emissions)) return

        do i_emit = 1, n_emission_species
            if (emission_indices(i_emit) > 0 .and. i_emit <= size(chem_emissions, 1)) then
                scalars(emission_indices(i_emit),:,1:nCells) = &
                    scalars(emission_indices(i_emit),:,1:nCells) + &
                    dt * chem_emissions(i_emit,:,1:nCells) * inv_aircol(:,:)
            end if
        end do

    end subroutine emissions_apply


    !------------------------------------------------------------------------
    !  routine sanitize_emission_field
    !------------------------------------------------------------------------
    subroutine sanitize_emission_field(field)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(kind=RKIND), dimension(:,:), intent(inout) :: field
        integer :: iCell, k

        do iCell = 1, size(field, 2)
            do k = 1, size(field, 1)
                if (.not. ieee_is_finite(field(k, iCell)) .or. field(k, iCell) < 0.0_RKIND) then
                    field(k, iCell) = 0.0_RKIND
                end if
            end do
        end do

    end subroutine sanitize_emission_field


    !------------------------------------------------------------------------
    !  routine emissions_finalize
    !------------------------------------------------------------------------
    subroutine emissions_finalize()

        if (allocated(emission_species)) deallocate(emission_species)
        if (allocated(emission_tracers)) deallocate(emission_tracers)
        n_emission_species = 0

    end subroutine emissions_finalize

end module mpas_chem_emissions_driver
