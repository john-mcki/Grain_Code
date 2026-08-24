  !! ----- Voronoi tessellation mictrostucture generation ----------------------------------------------------------------------- !!
  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  subroutine InitialState( Phi )
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    implicit none

    !! ----- Input/output arrays ----- !!
    real(real64), intent(inout), dimension(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,0:np) :: Phi

    !! ----- Subroutine-specific variables and arrays ----- !!
    integer(int32) :: ip, ix, iy, iz                              !! Various counters
    integer(int32) :: zMin
    real(real64) :: rand, dist2, dist2min                         !! Variables for constructing the Voronoi tessellation
    real(real64), dimension(np,3) :: centers                      !! Locations of grain centers
    real(real64), dimension(np) :: angles                       !! Angles for calculating distances
    real(real64), allocatable :: theta(:,:,:)        !! Local angle
    real(real64), allocatable :: pert(:,:), d_pert(:,:)        !! Array for the sinusoidal perturbations to the thickness
    real(real64) :: area_local, area_total
    real(real64) :: x_frac, x_left, x_right
    real(real64) :: Lx, Lz, x_loc, z_loc

    !! ----- Useful part of subroutine ----- !!
    !! allocate necessary arrays
    allocate( theta(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )

    Phi = 0.d0
    Lx = dble(nx) * dh

    !! Generate locations of grain centers for the Voronoi tessellation
    call InitializeRandomSeed()
    do ip = 1, np, 1
      call random_number(rand)
      centers(ip,1) = rand*real(nx,real64)*dh
      call random_number(rand)
      centers(ip,2) = rand*real(ny,real64)*dh
      call random_number(rand)
      centers(ip,3) = rand*real(nz,real64)*dh
      call random_number(rand)
      angles(ip) = pi*(2.0d0*rand - 1.0d0)
    end do

    !! Calculate distances of each point in the microstructure relative to the grain centers
    !! Note: is there a better way to do this that will be faster with large numbers of grains?
    do iz = z1-1,z2+1
      do iy = y1-1,y2+1
        do ix = x1-1,x2+1
          dist2min = huge(1.d0)
          do ip = 1, np, 1
            dist2 = (centers(ip,1)-(real(ix,real64)*dh))**2.d0 &
              + (centers(ip,2)-(real(iy,real64)*dh))**2.d0 &
              + (centers(ip,3)-(real(iz,real64)*dh))**2.d0
            if ( dist2 .lt. dist2min ) then
              dist2min = dist2
              theta(ix,iy,iz) = angles(ip)
            endif
          enddo
        enddo
      enddo
    enddo

    !! Set the initial order parameters
    do ip = 1,np
      do iz = z1-1,z2+1
        do iy = y1-1,y2+1
          do ix = x1-1,x2+1
            if ( theta(ix,iy,iz) .eq. angles(ip) ) then
              Phi(ix,iy,iz,ip) = 1.0d0
            else
              Phi(ix,iy,iz,ip) = 0.0d0
            endif
          enddo
        enddo
      enddo
    enddo

    deallocate( theta )
    CALL update_BCs_par_PF( Phi, 'N', 'N', 'N', 'N', 'N', 'N' )
  end subroutine InitialState

  !! ----- Allen-Cahn smoothing ------------------------------------------------------------------------------------------------- !!
  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  subroutine SmoothAC( Phi )
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none

    !! ----- Input/output arrays ----- !!
    real(real64), intent(inout), dimension(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,0:np) :: Phi

    !! ----- Subroutine-specific arrays and variables ----- !!
    integer(int32) :: ip, ix, iy, iz, it                          !! Various counter variables
    real(real64) :: dfdphi, lapl                                  !! Variational derivative of free energy and the laplacian
    real(real64), allocatable :: Phi_new(:,:,:,:) !! The updated order parameter
    real(real64), allocatable :: smsq(:,:,:)                      !! Sum of squared order parameters
    integer(int32) :: is

    !! ----- Useful part of subroutine ----- !!

    !! Allocate the arrays
    allocate( Phi_new(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,0:np) )
    allocate( smsq(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )

    is = 1
    if (smooth_Li) is = 0

    !if (rank .eq. 0) print *, is

    !! Evolve order parameters through time
    do it = 0,ssteps
      if ( any(ieee_is_nan(Phi)) ) then
        print *, 'ERROR: NaN value detected during smoothing on task, iteration: ', rank, it
        call MPI_Abort( cartcomm, 1, ierror )
      endif
      smsq = sum(Phi*Phi,4)
      do iz = z1,z2
        do iy = y1,y2
          do ix = x1,x2
            do ip = is,np
              !! Update the derivative of the free energy
              dfdphi = Phi(ix,iy,iz,ip)**3.d0 - Phi(ix,iy,iz,ip) &
                + 2.0*gamma*Phi(ix,iy,iz,ip)*(smsq(ix,iy,iz) - Phi(ix,iy,iz,ip)*Phi(ix,iy,iz,ip))

              !! Calculate the interfacial energy laplacian
              lapl = ( Phi(ix+1,iy,iz,ip) + Phi(ix-1,iy,iz,ip) &
                + Phi(ix,iy+1,iz,ip) + Phi(ix,iy-1,iz,ip) &
                + Phi(ix,iy,iz+1,ip) + Phi(ix,iy,iz-1,ip) &
                - 6.0d0*Phi(ix,iy,iz,ip) )/dhsq

              !! Evolve order parameters
              Phi_new(ix,iy,iz,ip) = Phi(ix,iy,iz,ip) - sLdt*( W*dfdphi - eps_sq*lapl )

              ! Keep Evolved Solid & Liquid Order Parameters Within Model Bounds
              ! if ( Phi_new(ix,iy,iz,ip) .lt. PhiCutoff ) then
              !   Phi_new(ix,iy,iz,ip) = PhiCutoff
              ! else if (Phi_new(ix,iy,iz,ip) .gt. 1.0d0 ) then
              !   Phi_new(ix,iy,iz,ip) = 1.0d0
              ! end if
            enddo
          enddo
        enddo
      enddo

      ! Apply Order Parameter Boundary Conditions in Buffer Layers
      do ip = 0,np
        CALL update_BCs_par( Phi_new(:,:,:,ip), Sol_W, Sol_E, Sol_N, Sol_S, Sol_U, Sol_D )
      end do
      
      !if (mod(it, 100) .eq. 0) print *, it, rank, maxval(abs(phi_new(:,:,:,1:) - phi(:,:,:,1:)))

      ! Update Order Parameter
      Phi(:,:,:,is:np) = Phi_new(:,:,:,is:np)

    end do

    where (Phi .lt. 0.d0) Phi = 0.d0
    where (Phi .gt. 1.d0) Phi = 1.d0

    !! Normalize all of the order parameters to sum to one
    do concurrent( ix=x1-1:x2+1, iy=y1-1:y2+1, iz=z1-1:z2+1 )
      Phi(ix,iy,iz,:) = Phi(ix,iy,iz,:) / sum(Phi(ix,iy,iz,:))
    enddo

    deallocate( Phi_new, smsq )

  end subroutine SmoothAC
  !------------------------------------------------------------------------------------------------------------------------------!
  ! End Order Parameter Smoothing Using Allen-Cahn Dynamics ---------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------------!
