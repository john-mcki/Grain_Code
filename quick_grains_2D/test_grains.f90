program main
IMPLICIT NONE

    INTEGER,PARAMETER :: x1 = 1, x2 = 256, y1 = 1, y2 = 256, z1 = 1, z2 = 1
    INTEGER,PARAMETER :: nx = 256, ny = 256, nz = 1, np = 200, ssteps =  4600
    REAL(8),PARAMETER :: pi = 3.14159265358979323846
    REAL(8),PARAMETER :: dh = 1.0d0
    REAL(8),PARAMETER :: gamma = 1.5d0, W = 1.0d0, eps_sq = 1d0, sLdt = 0.1d0
    
    real(8) :: dfdphi, lapl                                  !! Variational derivative of free energy and the laplacian
    real(8), allocatable :: Phi_new(:,:,:,:) !! The updated order parameter
    real(8), allocatable :: smsq(:,:,:)                      !! Sum of squared order parameters
    integer(4) :: is
    !! ----- Input/output arrays ----- !!
 !   real(8), intent(inout), dimension(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,0:np) :: Phi
    real(8), allocatable :: Phi(:,:,:,:), Phi2(:,:,:,:)
    integer(4), allocatable :: featureID(:,:,:)
    logical, allocatable :: mask(:,:,:)
    integer(4), allocatable :: index_to_featureID(:), index_to_index(:)

    !! ----- Subroutine-specific variables and arrays ----- !!
    integer(4) :: ip, jp, ix, iy, iz, di, it, npi                          !! Various counters
    integer(4) :: npeff, npeff_fin
    real(8) :: rand, dist2, dist2min   !! Variables for constructing the Voronoi tessellation
    real(8), dimension(np,3) :: centers                      !! Locations of grain centers
    logical, dimension(np) :: to_keep                      !! throw away grains based on center location
    real(8), dimension(3) :: dist, length, coord
    real(8), allocatable :: theta(:,:,:)        !! Local angle

    !! ----- Useful part of subroutine ----- !!
    !! allocate necessary arrays
    allocate( featureID(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )

    featureID = 0

    !! Generate locations of grain centers for the Voronoi tessellation
    length(1) = real(nx,8)*dh
    length(2) = real(ny,8)*dh
    
    call InitializeRandomSeed()
    do ip = 1, np, 1
      do di = 1, 2
        call random_number(rand)
        centers(ip,di) = rand*length(di)
      end do
    end do

    !! Calculate distances of each point in the microstructure relative to the grain centers
    !! Note: is there a better way to do this that will be faster with large numbers of grains?
    do iz = z1-1,z2+1
      do iy = y1-1,y2+1
        do ix = x1-1,x2+1
          dist2min = huge(1.d0)
          do ip = 1, np, 1
            coord(1)  = real(ix,8)*dh
            coord(2)  = real(iy,8)*dh
          ! periodic
            dist2 = 0d0
            do di = 1,2
              dist(di) = coord(di) - centers(ip,di)
              ! comment out line below for non periodic
              !dist(di) = modulo(modulo(dist(di) - length(di)/2d0, length(di)) + length(di), length(di))
              dist2 = dist2 + dist(di)*dist(di)
            end do
            if ( dist2 .lt. dist2min ) then
              dist2min = dist2
              featureID(ix,iy,iz) = ip
            endif
          enddo
        enddo
      enddo
    enddo
    
    ! knock out grains
    
    allocate( mask(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )
    jp = 1
    allocate( index_to_featureID(1:np) )
    do ip = 1,np
      mask = featureID == ip
      dist2 = 0d0
      do di = 1, 2
        dist2 = dist2 + (centers(ip,di) - dble(nx/2))**2
      end do
      if ((count(mask) < 10*10*3) .or. (dist2 > (dble(nx)/2*0.85)**2)) then
        where (mask)
          featureID = np+1
        end where
      else
        index_to_featureID(jp) = ip
        jp = jp+1
      endif
    enddo
    deallocate( mask)
    index_to_featureID(jp+1) = np+1
    npeff = jp+1
    print*, np, npeff
    
    allocate( Phi(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,1:npeff) )
    !! Set the initial order parameters
    do ip = 1,npeff
      do iz = z1-1,z2+1
        do iy = y1-1,y2+1
          do ix = x1-1,x2+1
            if ( featureID(ix,iy,iz) .eq. index_to_featureID(ip)) then
              Phi(ix,iy,iz,ip) = 1.0d0
            else
              Phi(ix,iy,iz,ip) = 0.0d0
            endif
          enddo
        enddo
      enddo
    enddo
    deallocate(featureID)
    deallocate(index_to_featureID)
    
!deallocate( theta )

   !! Allocate the arrays
    allocate( Phi_new(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,1:npeff) )
    allocate( smsq(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )

    !! Evolve order parameters through time
    do it = 0,ssteps
      smsq = sum(Phi*Phi,4)
      do iz = z1,z2
        do iy = y1,y2
          do ix = x1,x2
            do ip = 1,npeff
              !! Update the derivative of the free energy
              dfdphi = Phi(ix,iy,iz,ip)**3.d0 - Phi(ix,iy,iz,ip) &
                + 2.0*gamma*Phi(ix,iy,iz,ip)*(smsq(ix,iy,iz) - Phi(ix,iy,iz,ip)*Phi(ix,iy,iz,ip))

              !! Calculate the interfacial energy laplacian
              lapl = ( Phi(ix+1,iy,iz,ip) + Phi(ix-1,iy,iz,ip) &
                + Phi(ix,iy+1,iz,ip) + Phi(ix,iy-1,iz,ip) &
                + Phi(ix,iy,iz+1,ip) + Phi(ix,iy,iz-1,ip) &
                - 6.0d0*Phi(ix,iy,iz,ip) )/(dh*dh)

              !! Evolve order parameters
              Phi_new(ix,iy,iz,ip) = Phi(ix,iy,iz,ip) - sLdt*( W*dfdphi - eps_sq*lapl )

            enddo
          enddo
        enddo
      enddo

      ! Apply Order Parameter Boundary Conditions in Buffer Layers
      do ip = 1,npeff
        Phi_new(x1-1,:,:,ip) = Phi_new(x2,:,:,ip)
        Phi_new(x2+1,:,:,ip) = Phi_new(x1,:,:,ip)
        Phi_new(:,y1-1,:,ip) = Phi_new(:,y2,:,ip)
        Phi_new(:,y2+1,:,ip) = Phi_new(:,y1,:,ip)
        Phi_new(:,:,z1-1,ip) = Phi_new(:,:,z2,ip)
        Phi_new(:,:,z2+1,ip) = Phi_new(:,:,z1,ip)
      end do

      ! Update Order Parameter
      Phi(:,:,:,1:npeff) = Phi_new(:,:,:,1:npeff)
      print*, sum(phi(:,:,:,1))
      
      
    if (mod(it, 100) == 0) then
      deallocate( Phi_new )
      allocate( index_to_index(1:npeff) )
      jp = 0
      do ip = 1,npeff
         if (sum(Phi(:,:,:,ip)) > 20.d0) then
           jp = jp + 1
           index_to_index(jp) = ip
         endif
      end do
      npeff_fin = jp
      print*, npeff, npeff_fin
      allocate(Phi_new(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,1:npeff_fin))
      do ip = 1,npeff_fin
        Phi_new(:,:,:,ip) = Phi(:,:,:,index_to_index(ip))
      end do
      deallocate(Phi)
      allocate(Phi(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,1:npeff_fin))
      Phi = Phi_new
      deallocate(index_to_index)
      npeff = npeff_fin
    end if
    
    end do
    
    deallocate( Phi_new )
    ! more grain reduction
    allocate( index_to_index(1:npeff) )
    jp = 0
    do ip = 1,npeff
       if (sum(Phi(:,:,:,ip)) > 20.d0) then
         jp = jp + 1
         index_to_index(jp) = ip
       endif
    end do
    npeff_fin = jp
    print*, npeff, npeff_fin
    allocate(Phi_new(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,1:npeff_fin))
    do ip = 1,npeff_fin
      Phi_new(:,:,:,ip) = Phi(:,:,:,index_to_index(ip))
    end do
    deallocate(Phi)
    allocate(Phi(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,1:npeff_fin))
    Phi = Phi_new
    deallocate(Phi_new)

    where (Phi .lt. 0.d0) Phi = 0.d0
    where (Phi .gt. 1.d0) Phi = 1.d0

    !! Normalize all of the order parameters to sum to one
    do concurrent( ix=x1-1:x2+1, iy=y1-1:y2+1, iz=z1-1:z2+1 )
      Phi(ix,iy,iz,:) = Phi(ix,iy,iz,:) / sum(Phi(ix,iy,iz,:))
    enddo

    smsq(:,:,:) = 0d0
    do iz = z1,z2
      do iy = y1,y2
        do ix = x1,x2
          do ip = 1,npeff_fin
            do jp = 1,npeff_fin
              if (ip /= jp) smsq(ix,iy,iz) = smsq(ix,iy,iz) + Phi(ix,iy,iz,ip)*Phi(ix,iy,iz,jp)
            enddo
          enddo
        enddo
      enddo
    enddo
    OPEN(UNIT=20,FILE='grain_vis.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) smsq(x1:x2,y1:y2,z1:z2)
    CLOSE(20)
    
    OPEN(UNIT=20,FILE='grains.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) Phi(x1:x2,y1:y2,z1:z2,1:npeff_fin)
    CLOSE(20)    
    
  
end program main

subroutine InitializeRandomSeed()

    ! Define Variables & Arrays
    implicit none
    integer(4) :: sze
    integer(4), dimension(:), ALLOCATABLE :: seed

    ! Create Random Seed for Random Number Generation
    call random_seed(size=sze)
    allocate(seed(sze))

    seed = 44392
    
    call random_seed(put=seed)
    DEallocate(seed)

  end subroutine InitializeRandomSeed
  
  
