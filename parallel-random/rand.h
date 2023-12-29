      interface 

      function rand_f64() bind(c, name='rand_f64')
      import 
      implicit none
      real(c_double) rand_f64
      end function

      function rand_uint64() bind(c, name='rand_uint64')
      import 
      implicit none
      integer(c_long_long) rand_uint64
      end function

      subroutine srand_parallel(seed) bind(c, 
     &        name='srand_parallel')
      import 
      implicit none
      integer(c_long_long) seed
      end subroutine

      end interface
