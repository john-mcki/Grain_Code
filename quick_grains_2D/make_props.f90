program main
use tensor_mod
IMPLICIT NONE
  INTEGER,PARAMETER :: x1 = 1, x2 = 256, y1 = 1, y2 = 256, z1 = 1, z2 = 1
  INTEGER,PARAMETER :: nx = 256, ny = 256, nz = 1, np = 33
  INTEGER :: i, j, k, l, outgrain
  real(8) :: v1(3), v2(3), v3(3)
  real(8) :: Cel_ref(1:21), Bel_ref(1:6), D_ref(1:6)
  real(8), allocatable :: Cel(:,:,:,:), Bel(:,:,:,:), D(:,:,:,:), CelRot(:,:), BelRot(:,:), DRot(:,:)
  real(8), allocatable :: Cel2D(:,:,:,:), Bel2D(:,:,:,:), D2D(:,:,:,:)
  real(8), allocatable :: euler1(:), euler2(:), euler3(:)    
  real(8), allocatable :: PhiSolid(:,:,:,:)
  real(8), allocatable :: PsiOut(:,:)
  real(real64), dimension(3,3) :: Bel_tnsr, D_tnsr
  real(8) :: phisum
   !! Generate the elastic constant arrays

  allocate( euler1(1:np), euler2(1:np), euler3(1:np) )
  allocate( CelRot(1:np,1:21), BelRot(1:np,1:6), DRot(1:np,1:6) )
  call InitializeRandomSeed()

! ref: https://math.stackexchange.com/questions/442418/random-generation-of-rotation-matrices/1288873
    do i = 1,np
      call random_number( v1 )
      v1 = v1/norm2(v1)
      do while (.true.)
        call random_number( v2 )
        v2 = v2/norm2(v2)
        if (abs(DOT_PRODUCT(v1,v2)) < 0.99d0) EXIT
      end do
      v2 = v2 - DOT_PRODUCT(v1,v2)*v1
      v2 = v2/norm2(v2)
      v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
      v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
      v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
! ref: https://math.stackexchange.com/questions/3328656/convert-rotation-matrix-to-euler-angles-zyz-y-convention-analytically
      euler1(i) = atan2(v3(2), v3(1))
      euler2(i) = atan2(v3(1)*cos(euler1(i)) + v3(2)*sin(euler1(i)), v3(3))
      euler3(i) = atan2(v2(3),-v1(3))
      print*, v3
    end do
    Bel_ref = 0.d0
    Cel_ref = 0.d0
    D_ref = 0.d0
    
! mapping key
! Cel_mtrx(1,1) = Cel_vec_in(1)
! Cel_mtrx(2,2) = Cel_vec_in(2)
! Cel_mtrx(3,3) = Cel_vec_in(3)
! Cel_mtrx(4,4) = Cel_vec_in(4)
! Cel_mtrx(5,5) = Cel_vec_in(5)
! Cel_mtrx(6,6) = Cel_vec_in(6)
! Cel_mtrx(1,2) = Cel_vec_in(7)
! Cel_mtrx(1,3) = Cel_vec_in(8)
! Cel_mtrx(1,4) = Cel_vec_in(9)
! Cel_mtrx(1,5) = Cel_vec_in(10)
! Cel_mtrx(1,6) = Cel_vec_in(11)
! Cel_mtrx(2,3) = Cel_vec_in(12)
! Cel_mtrx(2,4) = Cel_vec_in(13)
! Cel_mtrx(2,5) = Cel_vec_in(14)
! Cel_mtrx(2,6) = Cel_vec_in(15)
! Cel_mtrx(3,4) = Cel_vec_in(16)
! Cel_mtrx(3,5) = Cel_vec_in(17)
! Cel_mtrx(3,6) = Cel_vec_in(18)
! Cel_mtrx(4,5) = Cel_vec_in(19)
! Cel_mtrx(4,6) = Cel_vec_in(20)
! Cel_mtrx(5,6) = Cel_vec_in(21)
    print*, "ch1"
    Cel_ref(1) = 260.d0
    Cel_ref(2) = 260.d0
    Cel_ref(3) = 200.d0
    Cel_ref(4) = 45.d0
    Cel_ref(5) = 45.d0
    Cel_ref(6) = 94.d0
    Cel_ref(7) = 84.d0
    Cel_ref(8) = 46.d0
    Cel_ref(9) = -16.d0
    Cel_ref(10:11) = 0.d0
    Cel_ref(12) = 46.d0
    Cel_ref(13) = 16.d0
    Cel_ref(14:15) = 0.d0
    Cel_ref(16:18) = 0.d0
    Cel_ref(19:20) = 0.d0
    Cel_ref(21) = -14.d0
! mapping key for rank 2 tensor
!    Bel_tnsr(1,1) = Bel_vec_in(1)
!    Bel_tnsr(2,2) = Bel_vec_in(2)
!    Bel_tnsr(3,3) = Bel_vec_in(3)
!    Bel_tnsr(2,3) = Bel_vec_in(4)
!    Bel_tnsr(1,3) = Bel_vec_in(5)
!    Bel_tnsr(1,2) = Bel_vec_in(6)
    Bel_ref(1) = 0.028
    Bel_ref(2) = 0.028
    Bel_ref(3) = -0.030
    Bel_ref(4) = 0.d0
    Bel_ref(5) = 0.d0
    Bel_ref(6) = 0.d0
    
    D_ref(1) = 1d-3
    D_ref(2) = 1d-3
    D_ref(3) = 1d-4
    D_ref(4) = 0.d0
    D_ref(5) = 0.d0
    D_ref(6) = 0.d0
    print*, "ch2"
    call load_Cel_tnsr( Cel_ref )
    Bel_tnsr =  load_Bel_tnsr( Bel_ref )
    D_tnsr = load_Bel_tnsr( D_ref )
    print*, "ch3"
    do i = 1,np
      CelRot(i,1:21) = rot_Cel_tnsr( euler1(i), euler2(i), euler3(i) )
      BelRot(i,1:6) = rot_Bel_tnsr( euler1(i), euler2(i), euler3(i), Bel_tnsr)
      DRot(i,1:6)  = rot_Bel_tnsr( euler1(i), euler2(i), euler3(i), D_tnsr)
    enddo

    print*, "ch4"
    
    print*, "reading in structure"
    allocate( PhiSolid(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,1:np) )
    
    OPEN(UNIT=20,FILE='grains.dat',FORM='UNFORMATTED',STATUS='OLD',ACTION='READ', access='stream')
    READ(20) PhiSolid(x1:x2,y1:y2,z1:z2,1:np)
    CLOSE(20)
    
    outgrain = maxloc(PhiSolid(1,1,1,1:np), 1)
    print*, "outer grainID: ", outgrain

    allocate( PsiOut(x1:x2+1,y1:y2+1))
    PsiOut(x1:x2+1,y1:y2+1) = 0.d0
    PsiOut(x1:x2,y1:y2) = 1.d0 - PhiSolid(x1:x2,y1:y2,1, outgrain)
!    PsiOut(x1:x2+1,y1:y2+1) = 1d0
    OPEN(UNIT=20,FILE='psi.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) PsiOut
    CLOSE(20)   

    PhiSolid(x1-1:x2+1,y1-1:y2+1,z1-1:z2+1,outgrain) = 0d0
    print*, "renormalizing order parameters"
    do k = z1-1,z2+1
      do j = y1-1,y2+1
        do i = x1-1,x2+1
          phisum = sum(PhiSolid(i,j,k,:))
          phisum = max(phisum, 1d-6)
          PhiSolid(i,j,k,:) = PhiSolid(i,j,k,:)/phisum
        enddo
      enddo
    enddo
    
    allocate( Cel(1:21,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )
    allocate( Bel(1:6,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )
    allocate( D(1:6,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )

    Cel = 0.d0
    print*, "starting elastic constants calculation"
    
    deallocate(PsiOut)
    
    do l = 1,np
      do k = z1-1,z2+1
        do j = y1-1,y2+1
          do i = x1-1,x2+1
            Cel(1:21,i,j,k) = Cel(1:21,i,j,k) + PhiSolid(i,j,k,l)*CelRot(l,1:21)
            Bel(1:6, i,j,k) = Bel(1:6, i,j,k) + PhiSolid(i,j,k,l)*BelRot(l,1:6)
            D(1:6, i,j,k) = D(1:6, i,j,k) + PhiSolid(i,j,k,l)*DRot(l,1:6)
          enddo
        enddo
      enddo
    enddo
    
    ! periodic boundary conditions
    Cel(1:21,x2+1,y1:y2+1,z1:z2+1) = Cel(1:21,x1,y1:y2+1,z1:z2+1)
    Cel(1:21,x1:x2+1,y2+1,z1:z2+1) = Cel(1:21,x1:x2+1,y1,z1:z2+1)
    Cel(1:21,x1:x2+1,y1:y2+1,z2+1) = Cel(1:21,x1:x2+1,y1:y2+1,z1)
    
    allocate( Cel2D(1:6,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )
    allocate( Bel2D(1:3,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )
    allocate( D2D(1:3,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) )   

    do k = z1-1,z2+1
      do j = y1-1,y2+1
        do i = x1-1,x2+1
          call reduce2D_4thrank(Cel(1:21,i,j,k) , Cel2D(1:6,i,j,k))
          call reduce2D_2ndrank(Bel(1:6,i,j,k) , Bel2D(1:3, i,j,k))
          call reduce2D_2ndrank(D(1:6,i,j,k) , D2D(1:3, i,j,k))
        enddo
      enddo
    enddo
    
! constant test cases
!    Cel2D(1:6,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 0.0
!    Cel2D(1:2,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 0.7
!    Cel2D(3,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 0.3
!    Cel2D(4:5,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 0.2
!    Bel2D(1:3,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 0.0
!    Bel2D(1:2,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 0.01
!    D2D(1:3,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 0.0
!    D2D(1:2,x1-1:x2+1,y1-1:y2+1,z1-1:z2+1) = 1d-4
    
    
    do i = 1,6
      print*, sum(Cel2D(i,x1:x2+1,y1:y2+1,z1:z2+1))/((nx+1)*(ny+1)*(nz+1))
    enddo
    do i = 1,3
      print*, sum(Bel2D(i,x1:x2+1,y1:y2+1,z1:z2+1))/((nx+1)*(ny+1)*(nz+1))
    enddo
    do i = 1,3
      print*, sum(D2D(i,x1:x2+1,y1:y2+1,z1:z2+1))/((nx+1)*(ny+1)*(nz+1))
    enddo
    
    Cel2D(1:3,x1:x2+1,y1:y2+1,1) = MAX(Cel2D(1:3,x1:x2+1,y1:y2+1,1), 1.d0)
    
    OPEN(UNIT=20,FILE='Cv2D_1.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) Cel2D(1:2,x1:x2+1,y1:y2+1,1)
    CLOSE(20)
    
    OPEN(UNIT=20,FILE='Cv2D_2.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) Cel2D(3:4,x1:x2+1,y1:y2+1,1)
    CLOSE(20)
    
    OPEN(UNIT=20,FILE='Cv2D_3.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) Cel2D(5:6,x1:x2+1,y1:y2+1,1)
    CLOSE(20)

    OPEN(UNIT=20,FILE='eigv2D_1.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) Bel2D(1:2,x1:x2+1,y1:y2+1,1)
    CLOSE(20)
    
    OPEN(UNIT=20,FILE='eigv2D_2.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) Bel2D(3,x1:x2+1,y1:y2+1,1)
    CLOSE(20)
    
    OPEN(UNIT=20,FILE='Dv2D_1.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) D2D(1:2,x1:x2+1,y1:y2+1,1)
    CLOSE(20)
    
    OPEN(UNIT=20,FILE='Dv2D_2.dat',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE', access='stream')
    WRITE(20) D2D(3,x1:x2+1,y1:y2+1,1)
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

    seed = 2000
    
    call random_seed(put=seed)
    DEallocate(seed)

end subroutine InitializeRandomSeed

subroutine reduce2D_2ndrank(Bel_vec_in, Bel_vec_out)
    ! Define Variables & Arrays
    implicit none
    real(8), dimension(1:6), intent(in) :: Bel_vec_in
    real(8), dimension(1:3), intent(out) :: Bel_vec_out
    
    Bel_vec_out(1) = Bel_vec_in(1)
    Bel_vec_out(2) = Bel_vec_in(2)
    Bel_vec_out(3) = Bel_vec_in(6)
    
end subroutine reduce2D_2ndrank

subroutine reduce2D_4thrank(Cel_vec_in, Cel_vec_out)
    ! Define Variables & Arrays
    implicit none
    real(8), dimension(1:21), intent(in) :: Cel_vec_in
    real(8), dimension(1:6), intent(out) :: Cel_vec_out
    
    Cel_vec_out(1) = Cel_vec_in(1) ! C11
    Cel_vec_out(2) = Cel_vec_in(2) ! C22
    Cel_vec_out(3) = Cel_vec_in(6) ! C66
    Cel_vec_out(4) = Cel_vec_in(7) ! C12
    Cel_vec_out(5) = Cel_vec_in(11) ! C16
    Cel_vec_out(6) = Cel_vec_in(15) ! C26

end subroutine reduce2D_4thrank

