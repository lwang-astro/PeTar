      Program test

      use iso_c_binding
      use omp_lib

      implicit none
      Real*8 count
      integer ip
      integer*8 seed
*      real*8 rand_f64
*      external rand_f64, srand_parallel
      include 'rand.h'

      seed = 1234
      call srand_parallel(seed)

!$omp parallel private(ip,count)
      ip = omp_get_thread_num()
      count = rand_f64()
      write(*,*) ip, count
!$omp end parallel 

      END
      
