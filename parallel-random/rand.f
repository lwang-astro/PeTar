      real*8 function rand_f64() 
      use iso_c_binding
      
      interface 
      function rand_f64_bind() bind(c, name='rand_f64')
      import 
      implicit none
      real(c_double) rand_f64_bind
      end function
      end interface

      rand_f64 = rand_f64_bind()
      
      return
      end

      integer*8 function rand_uint64() 
      use iso_c_binding
      
      interface 
      function rand_uint64_bind() bind(c, name='rand_uint64')
      import 
      implicit none
      integer(c_long_long) rand_uint64_bind
      end function
      end interface

      rand_uint64 = rand_uint64_bind()
      
      return
      end

      subroutine srand_parallel(seed, rank) 
      use iso_c_binding
      integer*8 seed
      integer*4 rank
      
      interface 
      subroutine srand_parallel_bind(seed, rank) bind(c, 
     &        name='srand_parallel')
      import 
      integer(c_long_long) seed
      integer(c_int) rank
      end subroutine
      end interface

      call srand_parallel_bind(seed, rank)

      end
