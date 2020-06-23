subroutine krome_initab(gama, mu)
   use krome_main
   use krome_user
   implicit none
   real*8::gama, mu

   call krome_init()
   call krome_set_user_Av(1.d2)
   call krome_set_user_dgomega(5.d-1)
   call krome_set_user_zeta(7.6923d0)
   call krome_set_user_rad(1.d0)
   call krome_set_user_gRad(1.d-5)
   call krome_set_user_gArea(2.4d-22)
   call krome_set_user_fr(1.d0)
   call krome_set_user_desorb(1.d0)
   call krome_set_user_therm(1.d0)
   call krome_set_user_h2desorb(1.d0)
   call krome_set_user_crdesorb(0.d0)
   call krome_set_user_uvcr(1.d0)
   call krome_set_gamma(gama)
   call krome_set_mu(mu)

end subroutine krome_initab
