!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!! FILE: tensor_mod.f90
!! LANG: Fortran 2008
!! AUTHOR: Alexander Chadwick
!! EMAIL: afchadwi@umich.edu
!! AFFILIATION: University of Michigan
!! DATE: January 12, 2017
!!
!! DESCRIPTION: This file includes a module and subroutines that are necessary for rotating the various tensors
!!      in the mechanical behavior model. It is straightforward enough to add additional conventions for the Euler
!!      angle rotations, but because DREAM.3D uses the 'zxz' convention that is the default here.
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module tensor_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! ----- module variables and paramters, if needed ---------------------------

  private :: Cel_tnsr, Cel_mtrx, pi

  real(real64), dimension(6,6) :: Cel_mtrx                                  !! Unrotated matrix of elastic constants
  real(real64), dimension(3,3,3,3) :: Cel_tnsr                              !! Unrotated tensor of elastic constants

!  real(real64), dimension(3,3) :: Bel_tnsr                                  !! Unrotated tensor of volumetric expansion coefficients

  real(real64), parameter :: pi = acos(-1._real64)                          !! Pi


  !! ----- module interfaces, as needed ----------------------------------------

  !! ----- module subroutines --------------------------------------------------

contains
  subroutine load_Cel_tnsr( Cel_vec_in )
    implicit none

    !! ----- input variables ---------------------------------------------------
    
    real(real64), dimension(1:21), intent(in) :: Cel_vec_in

    !! ----- subroutine variables ----------------------------------------------

    integer :: i, j, k, l, a, b                                             !! Pretty much just indexing variables

    !! ----- useful part of subroutine -----------------------------------------

    !! Populate the elastic constant matrix from the vector
    Cel_mtrx(1,1) = Cel_vec_in(1)
    Cel_mtrx(2,2) = Cel_vec_in(2)
    Cel_mtrx(3,3) = Cel_vec_in(3)
    Cel_mtrx(4,4) = Cel_vec_in(4)
    Cel_mtrx(5,5) = Cel_vec_in(5)
    Cel_mtrx(6,6) = Cel_vec_in(6)
    Cel_mtrx(1,2) = Cel_vec_in(7)
    Cel_mtrx(1,3) = Cel_vec_in(8)
    Cel_mtrx(1,4) = Cel_vec_in(9)
    Cel_mtrx(1,5) = Cel_vec_in(10)
    Cel_mtrx(1,6) = Cel_vec_in(11)
    Cel_mtrx(2,3) = Cel_vec_in(12)
    Cel_mtrx(2,4) = Cel_vec_in(13)
    Cel_mtrx(2,5) = Cel_vec_in(14)
    Cel_mtrx(2,6) = Cel_vec_in(15)
    Cel_mtrx(3,4) = Cel_vec_in(16)
    Cel_mtrx(3,5) = Cel_vec_in(17)
    Cel_mtrx(3,6) = Cel_vec_in(18)
    Cel_mtrx(4,5) = Cel_vec_in(19)
    Cel_mtrx(4,6) = Cel_vec_in(20)
    Cel_mtrx(5,6) = Cel_vec_in(21)
    
    !! Generate the rest of the matrix by symmetry
    do i = 2,6; do j = 1,i-1
      Cel_mtrx(i,j) = Cel_mtrx(j,i)
    enddo; enddo

    !! Convert the 6x6 matrix to the 3x3x3x3 tensor form
    do i = 1,3; do j = 1,3; do k = 1,3; do l = 1,3
      if ( i .eq. j ) then
        a = i
      elseif ( ( ( i .eq. 1 ) .and. ( j .eq. 2 ) ) .or. ( ( i .eq. 2 ) .and. ( j .eq. 1 ) ) ) then
        a = 6
      elseif ( ( ( i .eq. 1 ) .and. ( j .eq. 3 ) ) .or. ( ( i .eq. 3 ) .and. ( j .eq. 1 ) ) ) then
        a = 5
      else
        a = 4
      endif

      if ( k .eq. l ) then
        b = k
      elseif ( ( ( k .eq. 1 ) .and. ( l .eq. 2 ) ) .or. ( ( k .eq. 2 ) .and. ( l .eq. 1 ) ) ) then
        b = 6
      elseif ( ( ( k .eq. 1 ) .and. ( l .eq. 3 ) ) .or. ( ( k .eq. 3 ) .and. ( l .eq. 1 ) ) ) then
        b = 5
      else
        b = 4
      endif

      Cel_tnsr(i,j,k,l) = Cel_mtrx(a,b)
    enddo; enddo; enddo; enddo

  end subroutine load_Cel_tnsr

  function load_Bel_tnsr( Bel_vec_in ) result( Bel_tnsr1 )
    !! This subroutine loads the unrotated tensor for the expansion coefficients
    implicit none

    !! ----- input variables ---------------------------------------------------
    
    real(real64), dimension(1:6), intent(in) :: Bel_vec_in
    real(real64), dimension(3,3) :: Bel_tnsr1


    !! ----- subroutine variables ----------------------------------------------

    integer :: i, j                                                         !! Pretty much just indexing variables

    !! ----- useful part of subroutine -----------------------------------------

    !! Populate the elastic constant matrix from the vector
    Bel_tnsr1(1,1) = Bel_vec_in(1)
    Bel_tnsr1(2,2) = Bel_vec_in(2)
    Bel_tnsr1(3,3) = Bel_vec_in(3)
    Bel_tnsr1(2,3) = Bel_vec_in(4)
    Bel_tnsr1(1,3) = Bel_vec_in(5)
    Bel_tnsr1(1,2) = Bel_vec_in(6)
    
    !! Generate the rest of the matrix by symmetry
    do i = 2,3; do j = 1,i-1
      Bel_tnsr1(i,j) = Bel_tnsr1(j,i)
    enddo; enddo

  end function load_Bel_tnsr

  function rot_Cel_tnsr( ang1, ang2, ang3 ) result( Cel_vec_rot )
    !! This function actually rotates the elastic tensor for a given set of zxz 
    !! Euler angles and then returns a vector of the rotated constants.
    implicit none

    !! ----- input variables ---------------------------------------------------
    real(real64), intent(in) :: ang1, ang2, ang3                            !! Euler angles (in radians)

    !! ----- function variables ------------------------------------------------
    real(real64), dimension(21) :: Cel_vec_rot                              !! Rotated elastic constant vector (for output)
    real(real64), dimension(3,3,3,3) :: Cel_tnsr_rot                        !! Rotated elastic constant tensor
    real(real64), dimension(3,3) :: Q                                       !! Rotation matrix
    integer :: i, j, k, l, m, n, o, p                                       !! Just indexing variables

    !! ----- useful part of function -------------------------------------------

    !! ~~~~~~~~~~~~~~~~~~~~
    !! Build the rotation matrix
    !! ~~~~~~~~~~~~~~~~~~~~

    Q = matmul( RotZ(ang1), matmul( RotX(ang2), RotZ(ang3) ) )              !! Multiply the individual matrices together
    Q = transpose(Q)                                                        !! Matrix must be transposed for this to work

    !! ~~~~~~~~~~~~~~~~~~~~
    !! Rotate the tensor (do math!)
    !! ~~~~~~~~~~~~~~~~~~~~

    Cel_tnsr_rot = 0._real64
    do i = 1,3; do j = 1,3; do k = 1,3; do l = 1,3
      do m = 1,3; do n = 1,3; do o = 1,3; do p = 1,3
        Cel_tnsr_rot(i,j,k,l) = Cel_tnsr_rot(i,j,k,l) + Q(i,m)*Q(j,n)*Q(k,o)*Q(l,p)*Cel_tnsr(m,n,o,p)
      enddo; enddo; enddo; enddo
    enddo; enddo; enddo; enddo

    where ( abs(Cel_tnsr_rot) .le. 1.e10_real64*epsilon(1._real64)*maxval(abs(Cel_tnsr_rot)) ) Cel_tnsr_rot = 0._real64

    Cel_vec_rot(1) = Cel_tnsr_rot(1,1,1,1)
    Cel_vec_rot(2) = Cel_tnsr_rot(2,2,2,2)
    Cel_vec_rot(3) = Cel_tnsr_rot(3,3,3,3)
    Cel_vec_rot(4) = Cel_tnsr_rot(2,3,2,3)
    Cel_vec_rot(5) = Cel_tnsr_rot(1,3,1,3)
    Cel_vec_rot(6) = Cel_tnsr_rot(1,2,1,2)
    Cel_vec_rot(7) = Cel_tnsr_rot(1,1,2,2)
    Cel_vec_rot(8) = Cel_tnsr_rot(1,1,3,3)
    Cel_vec_rot(9) = Cel_tnsr_rot(1,1,2,3)
    Cel_vec_rot(10) = Cel_tnsr_rot(1,1,1,3)
    Cel_vec_rot(11) = Cel_tnsr_rot(1,1,1,2)
    Cel_vec_rot(12) = Cel_tnsr_rot(2,2,3,3)
    Cel_vec_rot(13) = Cel_tnsr_rot(2,2,2,3)
    Cel_vec_rot(14) = Cel_tnsr_rot(2,2,1,3)
    Cel_vec_rot(15) = Cel_tnsr_rot(2,2,1,2)
    Cel_vec_rot(16) = Cel_tnsr_rot(3,3,2,3)
    Cel_vec_rot(17) = Cel_tnsr_rot(3,3,1,3)
    Cel_vec_rot(18) = Cel_tnsr_rot(3,3,1,2)
    Cel_vec_rot(19) = Cel_tnsr_rot(2,3,1,3)
    Cel_vec_rot(20) = Cel_tnsr_rot(2,3,1,2)
    Cel_vec_rot(21) = Cel_tnsr_rot(1,3,1,2)

  end function rot_Cel_tnsr

  function rot_Bel_tnsr( ang1, ang2, ang3 , Bel_tnsr) result( Bel_vec_rot )
    !! This function actually rotates the expansion coefficient tensor for a given set of zxz 
    !! Euler angles and then returns a vector of the rotated constants.
    implicit none

    !! ----- input variables ---------------------------------------------------
    real(real64), intent(in) :: ang1, ang2, ang3                            !! Euler angles (in radians)
    real(real64), dimension(3,3), intent(in) :: Bel_tnsr
    !! ----- function variables ------------------------------------------------
    real(real64), dimension(6) :: Bel_vec_rot                              !! Rotated expansion constant vector (for output)
    real(real64), dimension(3,3) :: Bel_tnsr_rot                            !! Rotated expansion constant tensor
    real(real64), dimension(3,3) :: Q                                       !! Rotation matrix
    integer :: i, j, m, n                                       !! Just indexing variables

    !! ----- useful part of function -------------------------------------------

    !! ~~~~~~~~~~~~~~~~~~~~
    !! Build the rotation matrix
    !! ~~~~~~~~~~~~~~~~~~~~

    Q = matmul( RotZ(ang1), matmul( RotX(ang2), RotZ(ang3) ) )              !! Multiply the individual matrices together
    Q = transpose(Q)                                                        !! Matrix must be transposed for this to work

    !! ~~~~~~~~~~~~~~~~~~~~
    !! Rotate the tensor (do math!)
    !! ~~~~~~~~~~~~~~~~~~~~

    Bel_tnsr_rot = 0._real64
    do i = 1,3; do j = 1,3
      do m = 1,3; do n = 1,3
        Bel_tnsr_rot(i,j) = Bel_tnsr_rot(i,j) + Q(i,m)*Q(j,n)*Bel_tnsr(m,n)
      enddo; enddo
    enddo; enddo

    where ( abs(Bel_tnsr_rot) .le. 1.e10_real64*epsilon(1._real64)*maxval(abs(Bel_tnsr_rot)) ) Bel_tnsr_rot = 0._real64

    Bel_vec_rot(1) = Bel_tnsr_rot(1,1)
    Bel_vec_rot(2) = Bel_tnsr_rot(2,2)
    Bel_vec_rot(3) = Bel_tnsr_rot(3,3)
    Bel_vec_rot(4) = Bel_tnsr_rot(2,3)
    Bel_vec_rot(5) = Bel_tnsr_rot(1,3)
    Bel_vec_rot(6) = Bel_tnsr_rot(1,2)

  end function rot_Bel_tnsr

  function RotX( ang ) result( RotMtrx )
    !! This function generates a matrix around the x-axis
    implicit none
    
    !! ----- input variables ---------------------------------------------------
    real(real64), intent(in) :: ang                                         !! Input angle in radians

    !! ----- function variables ------------------------------------------------
    real(real64), dimension(3,3) :: RotMtrx                                 !! Output rotation matrix

    !! ----- useful part of function -------------------------------------------
    RotMtrx = 0._real64
    RotMtrx(1,1) = 1._real64
    RotMtrx(2,2) = cos(ang)
    RotMtrx(2,3) = -sin(ang)
    RotMtrx(3,2) = sin(ang)
    RotMtrx(3,3) = cos(ang)
  end function RotX

  function RotZ( ang ) result( RotMtrx )
    !! This function generates a matrix around the z-axis
    implicit none
    
    !! ----- input variables ---------------------------------------------------
    real(real64), intent(in) :: ang                                         !! Input angle in radians

    !! ----- function variables ------------------------------------------------
    real(real64), dimension(3,3) :: RotMtrx                                 !! Output rotation matrix

    !! ----- useful part of function -------------------------------------------
    RotMtrx = 0._real64
    RotMtrx(1,1) = cos(ang)
    RotMtrx(1,2) = -sin(ang)
    RotMtrx(2,1) = sin(ang)
    RotMtrx(2,2) = cos(ang)
    RotMtrx(3,3) = 1._real64
  end function RotZ

end module tensor_mod
