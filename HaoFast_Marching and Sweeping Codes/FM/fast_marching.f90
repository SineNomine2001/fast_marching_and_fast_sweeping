!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><
    module fmdata
      ! Range of coordinates will be -M to M
      integer(4) :: M = 201
      integer(4) :: N

      ! NUmber of boundary points
      integer(4) :: nbp 

      integer(4), allocatable :: s(:,:)
      real(8),    allocatable :: T(:,:)
      real(8),    allocatable :: b(:,:,:)
    end module
!---------------------------------------------------------<


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><
    program fm_ab
      use fmdata
      implicit none

      integer(4) :: i, j, k, ij(2)
      real(8)    :: deg
      real(8)    :: cmp, Tw, Te, Ts, Tn, Tij

      real(8), external :: bsch

      character  :: fname*256, c*5
      real(8)    :: L, dr(2)

      N = (M-1)*2+1
    
      ! global initialisation
      allocate(s(N,N), T(N,N))
      T = 1.0D10
      s = 2

      ! define the BCs & scaling
      read (*,'(A)') fname
      open (1, file=trim(fname))

      read (1,*) c, nbp
      allocate(b(2,nbp,1))

      read (1,*) b(:,:,1)
      close(1)

!#if defined DEBUG
!print*, nbp, b(:,nbp,1)
!#endif
       L = max(maxval(b(1,:,1))-minval(b(1,:,1)), maxval(b(2,:,1))-minval(b(2,:,1)))
!!     dr(1) = sum(b(1,:,1)) / dble(nbp)
!1     dr(2) = sum(b(2,:,1)) / dble(nbp)
       dr(1) = (maxval(b(1,:,1)) + minval(b(1,:,1))) / 2
       dr(2) = (maxval(b(2,:,1)) + minval(b(2,:,1))) / 2
!#if defined DEBUG
!print*, 'dr and L are:', dr, L
!#endif
      b(1,:,1) = (b(1,:,1) - dr(1)) / L * 1.00D0*dble(M-1) + dble(M-1)
      b(2,:,1) = (b(2,:,1) - dr(2)) / L * 1.00D0*dble(M-1) + dble(M-1)

      open (2,file='bcurve.dat')
      write(2,*) 'TITLE="boundary"'
      write(2,*) 'VARIABLES="X", "Y"'
      write(2,*) 'ZONE T="boundary"'
      write(2,'(X,3(A,I0),A)') 'I=', nbp, ', J=', 1, ', K=', 1, ', ZONETYPE=Ordered'
      write(2,*) 'DATAPACKING=POINT'
      write(2,'(2E18.10)') b(1:2,:,1) - dble(M)
      close(2)

      do k=1, nbp
        i = int(b(1,k,1))
        j = int(b(2,K,1))
        T(i  ,j  ) = bsch(i  ,j  ) !sqrt((dble(i  )-b(1,k,1))**2 + (dble(j  )-b(2,k,1))**2)
        T(i+1,j  ) = bsch(i+1,j  ) !sqrt((dble(i+1)-b(1,k,1))**2 + (dble(j  )-b(2,k,1))**2)
        T(i  ,j+1) = bsch(i  ,j+1) !sqrt((dble(i  )-b(1,k,1))**2 + (dble(j+1)-b(2,k,1))**2)
        T(i+1,j+1) = bsch(i+1,j+1) !sqrt((dble(i+1)-b(1,k,1))**2 + (dble(j+1)-b(2,k,1))**2)
        s(i:i+1,j:j+1) = 0
      end do

      do k=1, nbp
        ij = int(b(1:2,k,1))

        i = ij(1) - 1
        j = ij(2)
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if

        i = ij(1) - 1
        j = ij(2) + 1
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if

        i = ij(1)
        j = ij(2) - 1
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if

        i = ij(1) 
        j = ij(2) + 2
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if

        i = ij(1) + 1
        j = ij(2) - 1
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if

        i = ij(1) + 1
        j = ij(2) + 2
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if

        i = ij(1) + 2
        j = ij(2) 
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if

        i = ij(1) + 2
        j = ij(2) + 1
        Tw = T(i-1,j  )
        Te = T(i+1,j  )
        Ts = T(i  ,j-1)
        Tn = T(i  ,j+1)
        call quad(Tw, Te, Ts, Tn, Tij)
        if(s(i,j) /= 0) then
          T(i,j) = Tij
          s(i,j) = 1
        end if
      end do

      ! call the FM subroutine
      call fm

      ! output to plot3d


    end program
!---------------------------------------------------------<


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><
    subroutine fm
      use fmdata
      implicit none

      integer(4) :: ij(2), i, j
      real(8)    :: cmp, Tw, Te, Ts, Tn, Tij

!     integer(4), allocatable :: s(:,:)
!     real(8),    allocatable :: T(:,:)

!     allocate(s(N,N), T(N,N))

      ! initialisation
!     T = 1.0D10
!     s = 2

!     T(M-1,M) = 1.0D0
!     T(M+1,M) = 1.0D0
!     T(M,M-1) = 1.0D0
!     T(M,M+1) = 1.0D0
!     s(M-1,M) = 1
!     s(M+1,M) = 1
!     s(M,M-1) = 1
!     s(M,M+1) = 1

!     T(M,M) = 0.0D0
!     s(M,M) = 0

!      T(1+1,:) = 1.0D0
!      T(N-1,:) = 1.0D0
!      T(:,1+1) = 1.0D0
!      T(:,N-1) = 1.0D0
!      s(1+1,:) = 1
!      s(N-1,:) = 1
!      s(:,1+1) = 1
!      s(:,N-1) = 1
!
!      T(1,:) = 0.0D0
!      T(N,:) = 0.0D0
!      T(:,1) = 0.0D0
!      T(:,N) = 0.0D0
!      s(1,:) = 0
!      s(N,:) = 0
!      s(:,1) = 0
!      s(:,N) = 0

      ! marching
      do while(sum(s) > 0)
        if(mod(sum(s),500 ) == 0) write(*,'(1A,$)') '.'
        if(mod(sum(s),5000) == 0) write(*,'(I0)') sum(s)

        ! find the trial point wiht smallest T
        cmp = 2.0D10
        do i=1, N
          do j=1, N
            if(s(i,j) /= 1) cycle
            if(cmp > T(i,j)) then 
              cmp  = T(i,j)
              ij   =  [i,j]
            end if
          end do
        end do

        ! update knowns, trials, and unknowns
        i = ij(1)
        j = ij(2)
        s(i,j) = 0
        if(i > 1) then; if(s(i-1,j) /= 0) s(i-1,j) = 1; end if
        if(i < N) then; if(s(i+1,j) /= 0) s(i+1,j) = 1; end if
        if(j > 1) then; if(s(i,j-1) /= 0) s(i,j-1) = 1; end if
        if(j < N) then; if(s(i,j+1) /= 0) s(i,j+1) = 1; end if

        ! solve the quadratic
        ! i-1,j
        i = ij(1) - 1
        j = ij(2)
        if(i >= 1) then
          if(i > 1) then
            Tw = T(i-1,j)
          else
            Tw = T(i+1,j)
          end if

          if(i < N) then
            Te = T(i+1,j)
          else
            Tw = T(i-1,j)
          end if

          if(j > 1) then
            Ts = T(i,j-1)
          else
            Ts = T(i,j+1)
          end if

          if(j < N) then
            Tn = T(i,j+1)
          else
            Tn = T(i,j-1)
          end if
          call quad(Tw, Te, Ts, Tn, Tij)
          if(s(i,j) == 1) T(i,j) = Tij
        end if

        ! i+1,j
        i = ij(1) + 1
        j = ij(2)
        if(i <= N) then
          if(i > 1) then
            Tw = T(i-1,j)
          else
            Tw = T(i+1,j)
          end if

          if(i < N) then
            Te = T(i+1,j)
          else
            Tw = T(i-1,j)
          end if

          if(j > 1) then
            Ts = T(i,j-1)
          else
            Ts = T(i,j+1)
          end if

          if(j < N) then
            Tn = T(i,j+1)
          else
            Tn = T(i,j-1)
          end if
          call quad(Tw, Te, Ts, Tn, Tij)
          if(s(i,j) == 1) T(i,j) = Tij
        end if

        ! i,j-1
        i = ij(1)
        j = ij(2) - 1
        if(j >= 1) then
          if(i > 1) then
            Tw = T(i-1,j)
          else
            Tw = T(i+1,j)
          end if

          if(i < N) then
            Te = T(i+1,j)
          else
            Tw = T(i-1,j)
          end if

          if(j > 1) then
            Ts = T(i,j-1)
          else
            Ts = T(i,j+1)
          end if

          if(j < N) then
            Tn = T(i,j+1)
          else
            Tn = T(i,j-1)
          end if
          call quad(Tw, Te, Ts, Tn, Tij)
          if(s(i,j) == 1) T(i,j) = Tij
        end if

        ! i,j+1
        i = ij(1)
        j = ij(2) + 1
        if(j <= N) then
          if(i > 1) then
            Tw = T(i-1,j)
          else
            Tw = T(i+1,j)
          end if

          if(i < N) then
            Te = T(i+1,j)
          else
            Tw = T(i-1,j)
          end if

          if(j > 1) then
            Ts = T(i,j-1)
          else
            Ts = T(i,j+1)
          end if

          if(j < N) then
            Tn = T(i,j+1)
          else
            Tn = T(i,j-1)
          end if
          call quad(Tw, Te, Ts, Tn, Tij)
          if(s(i,j) == 1) T(i,j) = Tij
        end if
      end do

      print*
!!!   write(*,'(21F6.2)') ((T(i,j), i=M, N), j=N, M, -1)

      open (1,file='fm.g')
      write(1,*) N, N
      write(1,*) ((DBLE(i-M), i=1, N), j=1, N), ((DBLE(j-M), i=1, N), j=1, N)
      close(1)

      open (2,file='fm.q')
      write(2,*) N, N
      write(2,*) 1.0, 1.0, 1.0, 1.0
      write(2,*) ((T(i,j), i=1, N), j=1, N)
      write(2,*) ((T(i,j), i=1, N), j=1, N)
      write(2,*) ((T(i,j), i=1, N), j=1, N)
      write(2,*) ((T(i,j), i=1, N), j=1, N)
      close(2)

      open (3,file='fm.dat')
      write(3,*) 'VARIABLES="X", "Y", "FM"'
      do i = 1, N
        do j = 1, N
          if (DBLE(i-M) > -101 .and. DBLE(i-M) < 101 .and. DBLE(j-M) > -51 .and. DBLE(j-M) < 51) then
            write(3,*) DBLE(i-M)/10, DBLE(j-M)/10, T(i,j)/10
          end if
        end do
      end do
    end subroutine


    subroutine quad(Tw,Te,Ts,Tn,Tij)
      implicit none
      real(8) :: Tw, Te, Ts, Tn, Tij
      real(8) :: a, b, x

      a = min(Tw, Te)
      b = min(Ts, Tn)

      if(dabs(a-b) >= 1.0D0) then
        x = min(a,b) + 1.0D0
      else
        x = (a + b + dsqrt(2.0D0*1.0D0 - (a-b)**2)) / 2.0D0
      end if

      Tij = x

      return
    end subroutine


    real(8) function xp(x)
      implicit none
      real(8) :: x

      if(x > 0.0D0) then
        xp = x
      else
        xp = 0.0D0
      end if

      return
    end function


    real(8) function bsch(i,j)
      use fmdata
      implicit none
      integer(4) :: i
      integer(4) :: j
      integer(4) :: k

      bsch = 1.0D10

      do k=1, nbp
        bsch = min(bsch, sqrt((dble(i)-b(1,k,1))**2 + (dble(j)-b(2,k,1))**2))
      end do

      return
    end function
