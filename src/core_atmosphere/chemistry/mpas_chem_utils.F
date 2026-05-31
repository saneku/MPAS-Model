! Copyright (c) 2025 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_chem_utils
!
!> \brief Shared helper routines for MPAS chemistry drivers
!> \details
!>  Collects small stateless utilities used by the chemistry, emissions,
!>  and photolysis drivers: scalar-index lookup, species-list parsing,
!>  string normalization, simple tracer sanitization, and active-calendar
!>  year-length lookup.
!-----------------------------------------------------------------------
module mpas_chem_utils

    use mpas_kind_types, only: RKIND, StrKIND

    implicit none

    private

    public :: resolve_scalar_index
    public :: sanitize_species_by_index
    public :: parse_species_list
    public :: append_unique
    public :: target_tracer_for_emission
    public :: to_lower
    public :: chemistry_days_in_year

contains

    !------------------------------------------------------------------------
    !  function resolve_scalar_index
    !------------------------------------------------------------------------
    integer function resolve_scalar_index(state, tracer_name) result(idx)

        use mpas_derived_types, only: mpas_pool_type
        use mpas_pool_routines, only: mpas_pool_get_dimension

        type(mpas_pool_type), intent(in) :: state
        character(len=*), intent(in) :: tracer_name
        integer, pointer :: idx_ptr

        idx = -1
        nullify(idx_ptr)
        call mpas_pool_get_dimension(state, 'index_' // trim(tracer_name), idx_ptr)
        if (associated(idx_ptr)) idx = idx_ptr

    end function resolve_scalar_index


    !------------------------------------------------------------------------
    !  routine sanitize_species_by_index
    !------------------------------------------------------------------------
    subroutine sanitize_species_by_index(scalars, idx)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(kind=RKIND), dimension(:,:,:), intent(inout) :: scalars
        integer, intent(in) :: idx
        integer :: iCell, k

        if (idx < 1 .or. idx > size(scalars, 1)) return

        do iCell = 1, size(scalars, 3)
            do k = 1, size(scalars, 2)
                if (.not. ieee_is_finite(scalars(idx, k, iCell)) .or. scalars(idx, k, iCell) < 0.0_RKIND) then
                    scalars(idx, k, iCell) = 0.0_RKIND
                end if
            end do
        end do

    end subroutine sanitize_species_by_index


    !------------------------------------------------------------------------
    !  routine parse_species_list
    !------------------------------------------------------------------------
    subroutine parse_species_list(list_text, n_names, names)

        character(len=*), intent(in)                 :: list_text
        integer, intent(out)                         :: n_names
        character(len=64), allocatable, intent(out)  :: names(:)

        character(len=:), allocatable :: clean
        integer :: i, start_pos, nmax

        clean = trim(list_text)
        do i = 1, len(clean)
            if (clean(i:i) == ',' .or. clean(i:i) == ';' .or. clean(i:i) == char(9)) clean(i:i) = ' '
        end do

        nmax = max(1, len_trim(clean))
        allocate(names(nmax))
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
        end do

        call shrink_name_array(names, n_names)

    end subroutine parse_species_list


    !------------------------------------------------------------------------
    !  routine shrink_name_array
    !------------------------------------------------------------------------
    subroutine shrink_name_array(names, n_names)

        character(len=64), allocatable, intent(inout) :: names(:)
        integer, intent(in) :: n_names
        character(len=64), allocatable :: tmp(:)

        allocate(tmp(n_names))
        if (n_names > 0) tmp(:) = names(1:n_names)
        deallocate(names)
        call move_alloc(tmp, names)

    end subroutine shrink_name_array


    !------------------------------------------------------------------------
    !  routine append_unique
    !------------------------------------------------------------------------
    subroutine append_unique(names, n_names, new_name)

        character(len=64), intent(inout) :: names(:)
        integer, intent(inout) :: n_names
        character(len=*), intent(in) :: new_name
        integer :: i

        if (len_trim(new_name) == 0) return
        do i = 1, n_names
            if (trim(to_lower(names(i))) == trim(to_lower(new_name))) return
        end do
        n_names = n_names + 1
        names(n_names) = trim(new_name)

    end subroutine append_unique


    !------------------------------------------------------------------------
    !  function target_tracer_for_emission
    !------------------------------------------------------------------------
    function target_tracer_for_emission(emission_name, tracer_names, n_tracers) result(target)

        character(len=*), intent(in) :: emission_name
        character(len=64), intent(in) :: tracer_names(:)
        integer, intent(in) :: n_tracers
        character(len=64) :: target
        character(len=64) :: lower_emit, lower_tracer
        integer :: i

        lower_emit = to_lower(trim(emission_name))
        if (is_passive_emission(lower_emit)) then
            target = trim(emission_name)
        else
            target = 'q' // upper_chem_name(trim(emission_name))
        end if

        do i = 1, n_tracers
            lower_tracer = to_lower(trim(tracer_names(i)))
            if (lower_tracer(1:1) == 'q') then
                if (trim(lower_tracer(2:)) == trim(lower_emit)) then
                    target = tracer_names(i)
                    return
                end if
            end if
            if (trim(lower_tracer) == trim(lower_emit)) then
                target = tracer_names(i)
                return
            end if
        end do

    end function target_tracer_for_emission


    !------------------------------------------------------------------------
    !  function is_passive_emission
    !------------------------------------------------------------------------
    pure logical function is_passive_emission(emission_name)

        character(len=*), intent(in) :: emission_name
        character(len=len(emission_name)) :: lower_name

        lower_name = to_lower(trim(emission_name))

        is_passive_emission = trim(lower_name) == 'dust_1' .or. &
                              trim(lower_name) == 'dust1' .or. &
                              trim(lower_name) == 'bc' .or. &
                              trim(lower_name) == 'oc' .or. &
                              trim(lower_name) == 'pm25' .or. &
                              trim(lower_name) == 'pm2.5' .or. &
                              trim(lower_name) == 'pm10'

    end function is_passive_emission


    !------------------------------------------------------------------------
    !  function upper_chem_name
    !------------------------------------------------------------------------
    pure function upper_chem_name(text) result(out)

        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i, ia

        do i = 1, len(text)
            ia = iachar(text(i:i))
            if (ia >= iachar('a') .and. ia <= iachar('z')) then
                out(i:i) = achar(ia - 32)
            else
                out(i:i) = text(i:i)
            end if
        end do

    end function upper_chem_name


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
    !  function chemistry_days_in_year
    !------------------------------------------------------------------------
    integer function chemistry_days_in_year(configs, currTime) result(days_in_year)

        use mpas_derived_types, only: mpas_pool_type, MPAS_Time_type
        use mpas_pool_routines, only: mpas_pool_get_config
        use mpas_timekeeping, only: mpas_get_time, isLeapYear

        type(mpas_pool_type), intent(in) :: configs
        type(MPAS_Time_type), intent(in) :: currTime

        character(len=StrKIND), pointer :: calendar_type
        integer :: ierr, year

        days_in_year = 365
        year = 0
        call mpas_get_time(currTime, YYYY=year, ierr=ierr)

        nullify(calendar_type)
        call mpas_pool_get_config(configs, 'config_calendar_type', calendar_type)

        if (ierr == 0 .and. associated(calendar_type)) then
            if (trim(calendar_type) == 'gregorian' .and. isLeapYear(year)) days_in_year = 366
        end if

    end function chemistry_days_in_year

end module mpas_chem_utils
