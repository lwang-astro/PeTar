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

      subroutine srand_parallel(seed) 
      use iso_c_binding
      integer*8 seed
      
      interface 
      subroutine srand_parallel_bind(seed) bind(c, 
     &        name='srand_parallel')
      import 
      implicit none
      integer(c_long_long) seed
      end subroutine
      end interface

      call srand_parallel_bind(seed)

      end
