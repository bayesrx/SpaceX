subroutine bigtdsub(n,m,o,x,z,b,s,L)


  use omp_lib

  implicit none
  integer :: i
  integer :: j
  integer :: k
  integer :: n
  integer :: m
  integer :: p
  integer :: o
  integer :: ii
  integer :: L
  
  real(kind = 4)   :: x(n*m*o)
  real(kind = 4)  :: z(n*m*o*L) 
  real  :: b(n*n)
  real  :: s(n*n*L)
  real , allocatable :: a(:,:,:)  
  real , allocatable :: y(:,:,:)
  real , allocatable :: f(:,:,:,:)
  real , allocatable :: u(:,:,:,:)
  real , allocatable :: d(:,:,:,:)      
  real , allocatable :: c(:,:,:)
  real , allocatable :: g(:,:,:)
  real , allocatable :: v(:,:,:)
  real , allocatable :: h(:,:)
  real , allocatable :: e(:,:)
  real , allocatable :: r(:,:)
  real , allocatable :: t(:,:,:)
  real , allocatable :: sq(:)
  real , allocatable :: bsq(:,:)
  real , allocatable :: w(:,:)
  real , allocatable :: q(:,:,:)

  real , parameter :: pi = 3.141592653589793D+00
  integer :: thread_num
  real  wtime


  thread_num = omp_get_max_threads ( )
  p=L


  allocate(y(1:n,1:m,1:o))
  allocate(t(1:m,1:n,1:o))
  allocate(d(1:n,1:m,1:o,1:p))
  allocate(c(1:n,1:n,1:o))
  allocate(sq(1:n))
  allocate(bsq(1:n,1:n))
  allocate(w(1:n,1:n))
  allocate(v(1:n,1:n,1:o))
  allocate(h(1:n,1:n))
  allocate(e(1:n,1:n))
  allocate(r(1:n,1:n))
  allocate(f(1:m,1:n,1:o,1:p))
  allocate(u(1:n,1:n,1:o,1:p))
  allocate(q(1:n,1:n,1:p))  
  wtime = omp_get_wtime ( )

      y = reshape(x,(/n,m,o/))
      d = reshape(z,(/n,m,o,p/))


      
  
  
  do ii=1,o
     t(1:m,1:n,ii) = TRANSPOSE(y(1:n,1:m,ii))
     do j=1,p
        f(1:m,1:n,ii,j) = TRANSPOSE(d(1:n,1:m,ii,j))
        enddo
  enddo
  

  
  !$omp parallel shared ( y,t,c, n, m,v,p,h,e,r,w,u) default(none) 
  do ii=1,o

     call dot(y(1:n,1:m,ii),t(1:m,1:n,ii),n,m,c(1:n,1:n,ii))
     call cov2acor(c(1:n,1:n,ii),n,v(1:n,1:n,ii))
     do j=1,3
        call dot(y(1:n,1:m,ii),t(1:m,1:n,ii),n,m,h)
        call dot(d(1:n,1:m,ii,j),f(1:m,1:n,ii,j),n,m,e)

        r = h + e

        call cov2acor(r,n,w)

        u(1:n,1:n,ii,j)=w
enddo
enddo


call mean1sub(n,n,o,v,w)
        call mean2sub(n,n,o,L,u,q)
    

!$omp end parallel 

  wtime = omp_get_wtime ( ) - wtime
  b = pack(w,.true.)
  s = pack(q,.true.)
  
!
!  Free memory.
!
  deallocate ( y )
  deallocate ( t )

!
!  Terminate.
!


Contains
  Subroutine dot( a, b,n,m,c)
    Real, Dimension(:,:), Intent( In    ) :: a
    Real, Dimension(:,:), Intent(In  ):: b
    Real, Dimension(:,:),Intent(   Out ) :: c
    Integer :: i,p,j,n,m
  !$omp do
  do i = 1, n
    do j = 1, n
      c(i,j) = 0.0
      do k = 1, m
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
enddo

 !$omp end do
return
End Subroutine dot


  subroutine cov2acor( a,n,b)
    Real, Dimension(:,:), Intent( In    ) :: a
    Real, Dimension(:,:),Intent(   Out ) :: b
    Real :: sq(n)
    Real :: bsq(n,n)
    Integer :: i,p,j,n,m


 !$omp do
 do i=1,n
        sq(i)=sqrt((1/a(i,i)))
    bsq(1:n,i) = sq(i)
 enddo
 do i=1,n
    do j=1,n

       b(i,j) = sq(i) * a(i,j) * bsq(i,j)
    enddo
 enddo
 
 !$omp end do
 

end subroutine cov2acor

subroutine mean1sub(n,m,o,a,b)
  implicit none
  integer :: n,m,i,j,k,o
  real(kind = 4) :: f(n*m*o)
  
  real :: b(n,m),sum1
  real :: d(n*m)
  real :: a(n,m,o)
!  a = reshape(f,(/n,m,o/))

  
!$omp do  
  do i=1,n
     do j=1,m
        sum1 = 0.0
        do k=1,o
           sum1 = sum1 + a(i,j,k)
        enddo
        b(i,j) = sum1/o
     enddo
  enddo

  !$omp end do
  
end subroutine mean1sub

subroutine mean2sub(n,m,o,p,a,b)
  implicit none
  integer :: n,m,i,j,k,o,l,p
  real(kind = 4) :: f(n*m*o*p)
  real :: b(n,m,o),sum1
  real :: d(n*m*o)
  real :: a(n,m,o,p)
!  a = reshape(f,(/n,m,o,p/))
  
!$omp do  
  do i=1,n
     do j=1,m
        do k=1,p
           sum1 = 0.0
           do l=1,o
           sum1 = sum1 + a(i,j,l,k)
        enddo
        b(i,j,k) = sum1/o
     enddo
  enddo
enddo

!$omp end do
end subroutine mean2sub


end subroutine bigtdsub



