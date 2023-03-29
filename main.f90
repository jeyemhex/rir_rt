!
! main.f90
! Copyright (C) 2023 Edward Higgins <ed.higgins@york.ac.uk>
!
! Distributed under terms of the MIT license.
!

program rir_rt
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)

  integer, parameter :: num_rays = 5000
  integer, parameter :: num_frequency_bands = 6
  integer, parameter :: num_surfaces = 6

  real(dp), parameter :: pi = 3.1415926535_dp
  real(dp), parameter :: c = 343.0_dp ! speed of sound

  real(dp), parameter :: histogram_timestep = 0.001_dp
  real(dp), parameter :: impulse_response_time = 1.0_dp
  integer, parameter :: num_time_bins = ceiling(impulse_response_time / histogram_timestep)

  type :: ray_type
    real(dp) :: pos(3)
    real(dp) :: dir(3)
    real(dp) :: time
    real(dp) :: energy
  end type ray_type

  ! Description of the Room, the source and reciever
  real(dp) :: source_coord(3) = [2.0_dp, 2.0_dp, 2.0_dp]
  real(dp) :: reciever_coord(3) = [5.0_dp, 5.0_dp, 1.8_dp]
  real(dp) :: reciever_radius = 0.0875
  real(dp) :: room_dimensions(3) = [10.0_dp, 8.0_dp, 4.0_dp]

  real(dp), allocatable :: frequencies(:) ! Frequency bands to be sampled
  real(dp), allocatable :: A(:,:)         ! Absorption coefficients of the walls
  real(dp), allocatable :: R(:,:)         ! Reflection coefficients "  "   "
  real(dp), allocatable :: D(:,:)         ! Scattering coefficients "  "   "

  real(dp), allocatable :: energy_histogram(:,:)   ! Time/frequency energy histogram

  type(ray_type) :: ray
  real(dp)  :: displacement(3)
  real(dp)  :: distance
  real(dp)  :: ray_recv_energy
  real(dp)  :: ray_recv_vector(3)
  real(dp)  :: recv_time_of_arrival

  real(dp)  :: normal(3)
  real(dp)  :: cos_theta
  real(dp)  :: cos_alpha
  real(dp)  :: E
  real(dp)  :: random_dir(3)
  real(dp)  :: reflected_dir(3)
  integer   :: t_bin

  integer :: iband, iray, isurface

  call initialise_arrays()

  ! Calculate the energy histogram using ray tracing
  do iband = 1, num_frequency_bands
    do iray = 1, num_rays
      ray = ray_type(pos = source_coord,            &
                     dir = random_sample_sphere(),  &
                     time = 0.0_dp,                 &
                     energy = 1.0_dp)

       do while (ray%time < impulse_response_time)
         call get_impact_surface(ray, isurface, displacement)
         distance = norm2(displacement)

         ray%pos = ray%pos + displacement
         ray%time = ray%time + distance/c
         ray%energy = ray%energy * R(isurface, iband)

         ray_recv_energy = ray%energy*D(isurface, iband)
         ray_recv_vector = reciever_coord - ray%pos

         distance = norm2(ray_recv_vector)
         recv_time_of_arrival = ray%time + distance/c

         if (recv_time_of_arrival > impulse_response_time) cycle

         normal = get_surface_normal(isurface)
         cos_theta = dot_product(ray_recv_vector, normal) / norm2(ray_recv_vector)
         cos_alpha = sqrt(norm2(ray_recv_vector)**2 - reciever_radius**2) / norm2(ray_recv_vector)
         E = (1-cos_alpha) * 2*cos_theta * ray_recv_energy
         t_bin = nint(recv_time_of_arrival/histogram_timestep)
         energy_histogram(t_bin, iband) = E

         random_dir = random_sample_sphere()
         if (dot_product(random_dir, normal) < 0) random_dir = -random_dir

         reflected_dir = ray%dir - 2*dot_product(ray%dir, normal) * normal
         reflected_dir = reflected_dir / norm2(reflected_dir)

         ray%dir = D(isurface,iband)*random_dir + (1-D(isurface,iband))*reflected_dir
         ray%dir = ray%dir / norm2(ray%dir)

       end do
    end do
  end do



contains
  subroutine initialise_arrays()
    allocate(frequencies(num_frequency_bands))
    frequencies = [125, 250, 500, 1000, 2000, 4000]

    allocate(A(num_surfaces, num_frequency_bands))
    A = reshape([0.08_dp, 0.09_dp, 0.12_dp,  0.16_dp, 0.22_dp, 0.24_dp, &
                 0.08_dp, 0.09_dp, 0.12_dp,  0.16_dp, 0.22_dp, 0.24_dp, &
                 0.08_dp, 0.09_dp, 0.12_dp,  0.16_dp, 0.22_dp, 0.24_dp, &
                 0.08_dp, 0.09_dp, 0.12_dp,  0.16_dp, 0.22_dp, 0.24_dp, &
                 0.08_dp, 0.09_dp, 0.12_dp,  0.16_dp, 0.22_dp, 0.24_dp, &
                 0.08_dp, 0.09_dp, 0.12_dp,  0.16_dp, 0.22_dp, 0.24_dp], shape(A))

    allocate(R(num_surfaces, num_frequency_bands))
    R = sqrt(1-A)

    allocate(D(num_surfaces, num_frequency_bands))
    D = reshape([0.05_dp, 0.3_dp, 0.7_dp, 0.9_dp, 0.92_dp, 0.94_dp, &
                 0.05_dp, 0.3_dp, 0.7_dp, 0.9_dp, 0.92_dp, 0.94_dp, & 
                 0.05_dp, 0.3_dp, 0.7_dp, 0.9_dp, 0.92_dp, 0.94_dp, & 
                 0.05_dp, 0.3_dp, 0.7_dp, 0.9_dp, 0.92_dp, 0.94_dp, & 
                 0.01_dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.5_dp, &
                 0.01_dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.5_dp], shape(D))


    allocate(energy_histogram(num_time_bins, num_frequency_bands))
    energy_histogram = 0.0_dp

  end subroutine initialise_arrays

  function random_sample_sphere() result(r)
    real(dp) :: r(3)
    real(dp)  :: z
    real(dp)  :: lon
    real(dp)  :: lat
    real(dp)  :: s, x, y

    call random_number(z)
    z = 2*z - 1
    z = min(z, 1.0_dp)
    z = max(z, -1.0_dp)
    lat = acos(z)

    call random_number(lon)
    lon = 2*pi*lon

    s = sin(lat)
    x = cos(lon)*s
    y = sin(lon)*s
    r = [x,y,z]

  end function random_sample_sphere

  subroutine get_impact_surface(ray, isurface, displacement)
    type(ray_type), intent(in)  :: ray
    integer,        intent(out) :: isurface
    real(dp),       intent(out) :: displacement(3)

    real(dp) :: d
    real(dp) :: tmp

    isurface = 0
    d = huge(1.0_dp)

    ! Check for collisions with the x-surfaces
    if (ray%dir(1) < 0) then
      d = -ray%pos(1) / ray%dir(1)
      if (d == 0) d = huge(1.0_dp)
      isurface = 1
    else if (ray%dir(1) > 0) then
      d = (room_dimensions(1)-ray%pos(1)) / ray%dir(1)
      if (d == 0) d = huge(1.0_dp)
      isurface = 2
    end if

    ! Check for collisions with the y-surfaces
    if (ray%dir(2) < 0) then
      tmp = -ray%pos(2) / ray%dir(2)
      if (tmp < d .and. tmp > 0) then
        d = tmp
        isurface = 3
      end if
    else if (ray%dir(2) > 0) then
      tmp = (room_dimensions(2)-ray%pos(2)) / ray%dir(2)
      if (tmp < d .and. tmp > 0) then
        d = tmp
        isurface = 4
      end if
    end if

    ! Check for collisions with the z-surfaces
    if (ray%dir(3) < 0) then
      tmp = -ray%pos(3) / ray%dir(3)
      if (tmp < d .and. tmp > 0) then
        d = tmp
        isurface = 5
      end if
    else if (ray%dir(3) > 0) then
      tmp = (room_dimensions(3)-ray%pos(3)) / ray%dir(3)
      if (tmp < d .and. tmp > 0) then
        d = tmp
        isurface = 6
      end if
    end if

    displacement = d*ray%dir

  end subroutine get_impact_surface

  function get_surface_normal(i) result(normal)
    real(dp) :: normal(3)
    integer, intent(in) :: i

    select case(i)
    case(1)
      normal = [1,0,0]
    case(2)
      normal = [-1,0,0]
    case(3)
      normal = [0,1,0]
    case(4)
      normal = [0,-1,0]
    case(5)
      normal = [0,0,1]
    case(6)
      normal = [0,0,-1]
    case default
      error stop "Invalid surface"
    end select

  end function get_surface_normal
end program rir_rt
