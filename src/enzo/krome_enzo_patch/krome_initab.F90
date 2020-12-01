subroutine krome_initab(gama, mu)
  use krome_main
  use krome_user
  implicit none
  real*8::gama, mu

  call krome_init()
  ! call krome_set_user_Av(1.d2)
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


! C++ interface wrapper of conserveLinGetRef_x in subs
! subroutine krome_initref(x, ref)
!   use krome_commons
!   use krome_subs
!   implicit none
!   real*8 :: x(nmols), ref(natoms)

!   ref(:) = conserveLinGetRef_x(x(:))
! end subroutine krome_initref


! subroutine krome_renormref(x, ref)
!   use krome_user
!   use krome_commons
!   use krome_subs
!   implicit none
!   real*8 :: x(nmols), ref(natoms)

!   call krome_conserveLin_x(x(:),ref(:))
! end subroutine krome_renormref


subroutine krome_initref(d, &
   De, CHI, OI, HNCI, &
   HCNI, H2I, CI, HI, &
   H2OI, OHI, O2I, CH2I, &
   H2COI, HCOI, MGI, NH3I, &
   NOI, CNI, COI, N2I, &
   NH2I, CH3I, CH4I, NI, &
   NHI, HeI, HNOI, CH3OHI, &
   CO2I, H2CNI, HNCOI, NO2I, &
   O2HI, OCNI, CH3OH_DUSTI, HNCO_DUSTI, &
   H2CO_DUSTI, CH4_DUSTI, CO_DUSTI, H2O_DUSTI, &
   NO_DUSTI, CO2_DUSTI, N2_DUSTI, HCN_DUSTI, &
   NH3_DUSTI, O2_DUSTI, NO2_DUSTI, HNO_DUSTI, &
   O2H_DUSTI, H2CN_DUSTI, MG_DUSTI, HNC_DUSTI, &
   E_DUSTI, HCOII, HII, HOCII, &
   CII, CH2II, CHII, H2COII, &
   MGII, NH3II, NOII, CNII, &
   COII, N2II, O2II, H2OII, &
   NH2II, OII, OHII, CH3II, &
   CH4II, NII, HCNII, NHII, &
   H2II, HeII, HNOII, H2NOII, &
   H3II, H3COII, H3OII, HCNHII, &
   HCO2II, HeHII, N2HII, O2HII, &
   in, jn, kn, idim, &
   is, js, ks, ie, je, ke, &
   ATOM_C, ATOM_E, ATOM_H, ATOM_MG, &
   ATOM_O, ATOM_N, ATOM_HE)

!     SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
!
!     2014, KROME DEVELOPERS to interface the package with ENZO
!
!     PURPOSE:
!     Solve the multi-species rate and cool equations via KROME.
!
!     INPUTS:
!     in,jn,kn - dimensions of 3D fields
!
!     d        - total density field
!
!     is,ie    - start and end indices of active region (zero based)
!     idim     - dimensionality (rank) of problem
!
!     OUTPUTS:
!     update chemical abundances densities (HI, HII, etc) and energy
!
!     PARAMETERS:
!     mh      - H mass in cgs units
!
!-----------------------------------------------------------------------
!     USE KROME
use krome_main
use krome_user
use krome_constants
use krome_commons

implicit none

real*8,parameter::mh=p_mass !mass_h
real*8::dom,krome_x(krome_nmols),ref(krome_natoms)
real*8::d(in,jn,kn)
integer::in,jn,kn,is,js,ks,ie,je,ke,idim
integer::i,j,k

real*8::De(in,jn,kn)
real*8::CHI(in,jn,kn)
real*8::OI(in,jn,kn)
real*8::HNCI(in,jn,kn)
real*8::HCNI(in,jn,kn)
real*8::H2I(in,jn,kn)
real*8::CI(in,jn,kn)
real*8::HI(in,jn,kn)
real*8::H2OI(in,jn,kn)
real*8::OHI(in,jn,kn)
real*8::O2I(in,jn,kn)
real*8::CH2I(in,jn,kn)
real*8::H2COI(in,jn,kn)
real*8::HCOI(in,jn,kn)
real*8::MGI(in,jn,kn)
real*8::NH3I(in,jn,kn)
real*8::NOI(in,jn,kn)
real*8::CNI(in,jn,kn)
real*8::COI(in,jn,kn)
real*8::N2I(in,jn,kn)
real*8::NH2I(in,jn,kn)
real*8::CH3I(in,jn,kn)
real*8::CH4I(in,jn,kn)
real*8::NI(in,jn,kn)
real*8::NHI(in,jn,kn)
real*8::HeI(in,jn,kn)
real*8::HNOI(in,jn,kn)
real*8::CH3OHI(in,jn,kn)
real*8::CO2I(in,jn,kn)
real*8::H2CNI(in,jn,kn)
real*8::HNCOI(in,jn,kn)
real*8::NO2I(in,jn,kn)
real*8::O2HI(in,jn,kn)
real*8::OCNI(in,jn,kn)
real*8::CH3OH_DUSTI(in,jn,kn)
real*8::HNCO_DUSTI(in,jn,kn)
real*8::H2CO_DUSTI(in,jn,kn)
real*8::CH4_DUSTI(in,jn,kn)
real*8::CO_DUSTI(in,jn,kn)
real*8::H2O_DUSTI(in,jn,kn)
real*8::NO_DUSTI(in,jn,kn)
real*8::CO2_DUSTI(in,jn,kn)
real*8::N2_DUSTI(in,jn,kn)
real*8::HCN_DUSTI(in,jn,kn)
real*8::NH3_DUSTI(in,jn,kn)
real*8::O2_DUSTI(in,jn,kn)
real*8::NO2_DUSTI(in,jn,kn)
real*8::HNO_DUSTI(in,jn,kn)
real*8::O2H_DUSTI(in,jn,kn)
real*8::H2CN_DUSTI(in,jn,kn)
real*8::MG_DUSTI(in,jn,kn)
real*8::HNC_DUSTI(in,jn,kn)
real*8::E_DUSTI(in,jn,kn)
real*8::HCOII(in,jn,kn)
real*8::HII(in,jn,kn)
real*8::HOCII(in,jn,kn)
real*8::CII(in,jn,kn)
real*8::CH2II(in,jn,kn)
real*8::CHII(in,jn,kn)
real*8::H2COII(in,jn,kn)
real*8::MGII(in,jn,kn)
real*8::NH3II(in,jn,kn)
real*8::NOII(in,jn,kn)
real*8::CNII(in,jn,kn)
real*8::COII(in,jn,kn)
real*8::N2II(in,jn,kn)
real*8::O2II(in,jn,kn)
real*8::H2OII(in,jn,kn)
real*8::NH2II(in,jn,kn)
real*8::OII(in,jn,kn)
real*8::OHII(in,jn,kn)
real*8::CH3II(in,jn,kn)
real*8::CH4II(in,jn,kn)
real*8::NII(in,jn,kn)
real*8::HCNII(in,jn,kn)
real*8::NHII(in,jn,kn)
real*8::H2II(in,jn,kn)
real*8::HeII(in,jn,kn)
real*8::HNOII(in,jn,kn)
real*8::H2NOII(in,jn,kn)
real*8::H3II(in,jn,kn)
real*8::H3COII(in,jn,kn)
real*8::H3OII(in,jn,kn)
real*8::HCNHII(in,jn,kn)
real*8::HCO2II(in,jn,kn)
real*8::HeHII(in,jn,kn)
real*8::N2HII(in,jn,kn)
real*8::O2HII(in,jn,kn)

real*8::ATOM_H(in,jn,kn)
real*8::ATOM_E(in,jn,kn)
real*8::ATOM_HE(in,jn,kn)
real*8::ATOM_C(in,jn,kn)
real*8::ATOM_N(in,jn,kn)
real*8::ATOM_O(in,jn,kn)
real*8::ATOM_MG(in,jn,kn)

!******************************

!set units
! dom = urho/mh

!loop over zones
do k = ks+1, ke+1
  do j = js+1, je+1
    do i = is+1, ie+1

      krome_x(krome_idx_E) = De(i,j,k)
      krome_x(krome_idx_CH) = CHI(i,j,k)
      krome_x(krome_idx_O) = OI(i,j,k)
      krome_x(krome_idx_HNC) = HNCI(i,j,k)
      krome_x(krome_idx_HCN) = HCNI(i,j,k)
      krome_x(krome_idx_H2) = H2I(i,j,k)
      krome_x(krome_idx_C) = CI(i,j,k)
      krome_x(krome_idx_H) = HI(i,j,k)
      krome_x(krome_idx_H2O) = H2OI(i,j,k)
      krome_x(krome_idx_OH) = OHI(i,j,k)
      krome_x(krome_idx_O2) = O2I(i,j,k)
      krome_x(krome_idx_CH2) = CH2I(i,j,k)
      krome_x(krome_idx_H2CO) = H2COI(i,j,k)
      krome_x(krome_idx_HCO) = HCOI(i,j,k)
      krome_x(krome_idx_MG) = MGI(i,j,k)
      krome_x(krome_idx_NH3) = NH3I(i,j,k)
      krome_x(krome_idx_NO) = NOI(i,j,k)
      krome_x(krome_idx_CN) = CNI(i,j,k)
      krome_x(krome_idx_CO) = COI(i,j,k)
      krome_x(krome_idx_N2) = N2I(i,j,k)
      krome_x(krome_idx_NH2) = NH2I(i,j,k)
      krome_x(krome_idx_CH3) = CH3I(i,j,k)
      krome_x(krome_idx_CH4) = CH4I(i,j,k)
      krome_x(krome_idx_N) = NI(i,j,k)
      krome_x(krome_idx_NH) = NHI(i,j,k)
      krome_x(krome_idx_HE) = HeI(i,j,k)
      krome_x(krome_idx_HNO) = HNOI(i,j,k)
      krome_x(krome_idx_CH3OH) = CH3OHI(i,j,k)
      krome_x(krome_idx_CO2) = CO2I(i,j,k)
      krome_x(krome_idx_H2CN) = H2CNI(i,j,k)
      krome_x(krome_idx_HNCO) = HNCOI(i,j,k)
      krome_x(krome_idx_NO2) = NO2I(i,j,k)
      krome_x(krome_idx_O2H) = O2HI(i,j,k)
      krome_x(krome_idx_OCN) = OCNI(i,j,k)
      krome_x(krome_idx_CH3OH_DUST) = CH3OH_DUSTI(i,j,k)
      krome_x(krome_idx_HNCO_DUST) = HNCO_DUSTI(i,j,k)
      krome_x(krome_idx_H2CO_DUST) = H2CO_DUSTI(i,j,k)
      krome_x(krome_idx_CH4_DUST) = CH4_DUSTI(i,j,k)
      krome_x(krome_idx_CO_DUST) = CO_DUSTI(i,j,k)
      krome_x(krome_idx_H2O_DUST) = H2O_DUSTI(i,j,k)
      krome_x(krome_idx_NO_DUST) = NO_DUSTI(i,j,k)
      krome_x(krome_idx_CO2_DUST) = CO2_DUSTI(i,j,k)
      krome_x(krome_idx_N2_DUST) = N2_DUSTI(i,j,k)
      krome_x(krome_idx_HCN_DUST) = HCN_DUSTI(i,j,k)
      krome_x(krome_idx_NH3_DUST) = NH3_DUSTI(i,j,k)
      krome_x(krome_idx_O2_DUST) = O2_DUSTI(i,j,k)
      krome_x(krome_idx_NO2_DUST) = NO2_DUSTI(i,j,k)
      krome_x(krome_idx_HNO_DUST) = HNO_DUSTI(i,j,k)
      krome_x(krome_idx_O2H_DUST) = O2H_DUSTI(i,j,k)
      krome_x(krome_idx_H2CN_DUST) = H2CN_DUSTI(i,j,k)
      krome_x(krome_idx_MG_DUST) = MG_DUSTI(i,j,k)
      krome_x(krome_idx_HNC_DUST) = HNC_DUSTI(i,j,k)
      krome_x(krome_idx_E_DUST) = E_DUSTI(i,j,k)
      krome_x(krome_idx_HCOj) = HCOII(i,j,k)
      krome_x(krome_idx_Hj) = HII(i,j,k)
      krome_x(krome_idx_HOCj) = HOCII(i,j,k)
      krome_x(krome_idx_Cj) = CII(i,j,k)
      krome_x(krome_idx_CH2j) = CH2II(i,j,k)
      krome_x(krome_idx_CHj) = CHII(i,j,k)
      krome_x(krome_idx_H2COj) = H2COII(i,j,k)
      krome_x(krome_idx_MGj) = MGII(i,j,k)
      krome_x(krome_idx_NH3j) = NH3II(i,j,k)
      krome_x(krome_idx_NOj) = NOII(i,j,k)
      krome_x(krome_idx_CNj) = CNII(i,j,k)
      krome_x(krome_idx_COj) = COII(i,j,k)
      krome_x(krome_idx_N2j) = N2II(i,j,k)
      krome_x(krome_idx_O2j) = O2II(i,j,k)
      krome_x(krome_idx_H2Oj) = H2OII(i,j,k)
      krome_x(krome_idx_NH2j) = NH2II(i,j,k)
      krome_x(krome_idx_Oj) = OII(i,j,k)
      krome_x(krome_idx_OHj) = OHII(i,j,k)
      krome_x(krome_idx_CH3j) = CH3II(i,j,k)
      krome_x(krome_idx_CH4j) = CH4II(i,j,k)
      krome_x(krome_idx_Nj) = NII(i,j,k)
      krome_x(krome_idx_HCNj) = HCNII(i,j,k)
      krome_x(krome_idx_NHj) = NHII(i,j,k)
      krome_x(krome_idx_H2j) = H2II(i,j,k)
      krome_x(krome_idx_HEj) = HeII(i,j,k)
      krome_x(krome_idx_HNOj) = HNOII(i,j,k)
      krome_x(krome_idx_H2NOj) = H2NOII(i,j,k)
      krome_x(krome_idx_H3j) = H3II(i,j,k)
      krome_x(krome_idx_H3COj) = H3COII(i,j,k)
      krome_x(krome_idx_H3Oj) = H3OII(i,j,k)
      krome_x(krome_idx_HCNHj) = HCNHII(i,j,k)
      krome_x(krome_idx_HCO2j) = HCO2II(i,j,k)
      krome_x(krome_idx_HEHj) = HeHII(i,j,k)
      krome_x(krome_idx_N2Hj) = N2HII(i,j,k)
      krome_x(krome_idx_O2Hj) = O2HII(i,j,k)

      ref(:) = krome_conserveLinGetRef_x(krome_x(:))
      
      ATOM_H(i,j,k) = ref(idx_atom_H)
      ATOM_HE(i,j,k) = ref(idx_atom_HE)
      ATOM_C(i,j,k) = ref(idx_atom_C)
      ATOM_N(i,j,k) = ref(idx_atom_N)
      ATOM_O(i,j,k) = ref(idx_atom_O)
      ATOM_MG(i,j,k) = ref(idx_atom_MG)
      ATOM_E(i,j,k) = ref(idx_atom_E)

    end do
  end do
end do

end subroutine krome_initref

subroutine krome_renormref(d, &
   De, CHI, OI, HNCI, &
   HCNI, H2I, CI, HI, &
   H2OI, OHI, O2I, CH2I, &
   H2COI, HCOI, MGI, NH3I, &
   NOI, CNI, COI, N2I, &
   NH2I, CH3I, CH4I, NI, &
   NHI, HeI, HNOI, CH3OHI, &
   CO2I, H2CNI, HNCOI, NO2I, &
   O2HI, OCNI, CH3OH_DUSTI, HNCO_DUSTI, &
   H2CO_DUSTI, CH4_DUSTI, CO_DUSTI, H2O_DUSTI, &
   NO_DUSTI, CO2_DUSTI, N2_DUSTI, HCN_DUSTI, &
   NH3_DUSTI, O2_DUSTI, NO2_DUSTI, HNO_DUSTI, &
   O2H_DUSTI, H2CN_DUSTI, MG_DUSTI, HNC_DUSTI, &
   E_DUSTI, HCOII, HII, HOCII, &
   CII, CH2II, CHII, H2COII, &
   MGII, NH3II, NOII, CNII, &
   COII, N2II, O2II, H2OII, &
   NH2II, OII, OHII, CH3II, &
   CH4II, NII, HCNII, NHII, &
   H2II, HeII, HNOII, H2NOII, &
   H3II, H3COII, H3OII, HCNHII, &
   HCO2II, HeHII, N2HII, O2HII, &
   in, jn, kn, idim, &
   is, js, ks, ie, je, ke, &
   ATOM_C, ATOM_E, ATOM_H, ATOM_MG, &
   ATOM_O, ATOM_N, ATOM_HE)

!     SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
!
!     2014, KROME DEVELOPERS to interface the package with ENZO
!
!     PURPOSE:
!     Solve the multi-species rate and cool equations via KROME.
!
!     INPUTS:
!     in,jn,kn - dimensions of 3D fields
!
!     d        - total density field
!
!     is,ie    - start and end indices of active region (zero based)
!     idim     - dimensionality (rank) of problem
!
!     OUTPUTS:
!     update chemical abundances densities (HI, HII, etc) and energy
!
!     PARAMETERS:
!     mh      - H mass in cgs units
!
!-----------------------------------------------------------------------
!     USE KROME
use krome_main
use krome_user
use krome_constants
use krome_commons

implicit none

real*8,parameter::mh=p_mass !mass_h
real*8::dom,krome_x(krome_nmols),ref(krome_natoms)
real*8::d(in,jn,kn)
integer::in,jn,kn,is,js,ks,ie,je,ke,idim
integer::i,j,k

real*8::De(in,jn,kn)
real*8::CHI(in,jn,kn)
real*8::OI(in,jn,kn)
real*8::HNCI(in,jn,kn)
real*8::HCNI(in,jn,kn)
real*8::H2I(in,jn,kn)
real*8::CI(in,jn,kn)
real*8::HI(in,jn,kn)
real*8::H2OI(in,jn,kn)
real*8::OHI(in,jn,kn)
real*8::O2I(in,jn,kn)
real*8::CH2I(in,jn,kn)
real*8::H2COI(in,jn,kn)
real*8::HCOI(in,jn,kn)
real*8::MGI(in,jn,kn)
real*8::NH3I(in,jn,kn)
real*8::NOI(in,jn,kn)
real*8::CNI(in,jn,kn)
real*8::COI(in,jn,kn)
real*8::N2I(in,jn,kn)
real*8::NH2I(in,jn,kn)
real*8::CH3I(in,jn,kn)
real*8::CH4I(in,jn,kn)
real*8::NI(in,jn,kn)
real*8::NHI(in,jn,kn)
real*8::HeI(in,jn,kn)
real*8::HNOI(in,jn,kn)
real*8::CH3OHI(in,jn,kn)
real*8::CO2I(in,jn,kn)
real*8::H2CNI(in,jn,kn)
real*8::HNCOI(in,jn,kn)
real*8::NO2I(in,jn,kn)
real*8::O2HI(in,jn,kn)
real*8::OCNI(in,jn,kn)
real*8::CH3OH_DUSTI(in,jn,kn)
real*8::HNCO_DUSTI(in,jn,kn)
real*8::H2CO_DUSTI(in,jn,kn)
real*8::CH4_DUSTI(in,jn,kn)
real*8::CO_DUSTI(in,jn,kn)
real*8::H2O_DUSTI(in,jn,kn)
real*8::NO_DUSTI(in,jn,kn)
real*8::CO2_DUSTI(in,jn,kn)
real*8::N2_DUSTI(in,jn,kn)
real*8::HCN_DUSTI(in,jn,kn)
real*8::NH3_DUSTI(in,jn,kn)
real*8::O2_DUSTI(in,jn,kn)
real*8::NO2_DUSTI(in,jn,kn)
real*8::HNO_DUSTI(in,jn,kn)
real*8::O2H_DUSTI(in,jn,kn)
real*8::H2CN_DUSTI(in,jn,kn)
real*8::MG_DUSTI(in,jn,kn)
real*8::HNC_DUSTI(in,jn,kn)
real*8::E_DUSTI(in,jn,kn)
real*8::HCOII(in,jn,kn)
real*8::HII(in,jn,kn)
real*8::HOCII(in,jn,kn)
real*8::CII(in,jn,kn)
real*8::CH2II(in,jn,kn)
real*8::CHII(in,jn,kn)
real*8::H2COII(in,jn,kn)
real*8::MGII(in,jn,kn)
real*8::NH3II(in,jn,kn)
real*8::NOII(in,jn,kn)
real*8::CNII(in,jn,kn)
real*8::COII(in,jn,kn)
real*8::N2II(in,jn,kn)
real*8::O2II(in,jn,kn)
real*8::H2OII(in,jn,kn)
real*8::NH2II(in,jn,kn)
real*8::OII(in,jn,kn)
real*8::OHII(in,jn,kn)
real*8::CH3II(in,jn,kn)
real*8::CH4II(in,jn,kn)
real*8::NII(in,jn,kn)
real*8::HCNII(in,jn,kn)
real*8::NHII(in,jn,kn)
real*8::H2II(in,jn,kn)
real*8::HeII(in,jn,kn)
real*8::HNOII(in,jn,kn)
real*8::H2NOII(in,jn,kn)
real*8::H3II(in,jn,kn)
real*8::H3COII(in,jn,kn)
real*8::H3OII(in,jn,kn)
real*8::HCNHII(in,jn,kn)
real*8::HCO2II(in,jn,kn)
real*8::HeHII(in,jn,kn)
real*8::N2HII(in,jn,kn)
real*8::O2HII(in,jn,kn)

real*8::ATOM_H(in,jn,kn)
real*8::ATOM_E(in,jn,kn)
real*8::ATOM_HE(in,jn,kn)
real*8::ATOM_C(in,jn,kn)
real*8::ATOM_N(in,jn,kn)
real*8::ATOM_O(in,jn,kn)
real*8::ATOM_MG(in,jn,kn)

!******************************

!set units
! dom = urho/mh

!loop over zones
do k = ks+1, ke+1
  do j = js+1, je+1
    do i = is+1, ie+1

      ref(idx_atom_H) = ATOM_H(i,j,k)
      ref(idx_atom_HE) = ATOM_HE(i,j,k)
      ref(idx_atom_C) = ATOM_C(i,j,k)
      ref(idx_atom_N) = ATOM_N(i,j,k)
      ref(idx_atom_O) = ATOM_O(i,j,k)
      ref(idx_atom_MG) = ATOM_MG(i,j,k)
      ref(idx_atom_E) = ATOM_E(i,j,k)

      krome_x(krome_idx_E) = De(i,j,k)
      krome_x(krome_idx_CH) = CHI(i,j,k)
      krome_x(krome_idx_O) = OI(i,j,k)
      krome_x(krome_idx_HNC) = HNCI(i,j,k)
      krome_x(krome_idx_HCN) = HCNI(i,j,k)
      krome_x(krome_idx_H2) = H2I(i,j,k)
      krome_x(krome_idx_C) = CI(i,j,k)
      krome_x(krome_idx_H) = HI(i,j,k)
      krome_x(krome_idx_H2O) = H2OI(i,j,k)
      krome_x(krome_idx_OH) = OHI(i,j,k)
      krome_x(krome_idx_O2) = O2I(i,j,k)
      krome_x(krome_idx_CH2) = CH2I(i,j,k)
      krome_x(krome_idx_H2CO) = H2COI(i,j,k)
      krome_x(krome_idx_HCO) = HCOI(i,j,k)
      krome_x(krome_idx_MG) = MGI(i,j,k)
      krome_x(krome_idx_NH3) = NH3I(i,j,k)
      krome_x(krome_idx_NO) = NOI(i,j,k)
      krome_x(krome_idx_CN) = CNI(i,j,k)
      krome_x(krome_idx_CO) = COI(i,j,k)
      krome_x(krome_idx_N2) = N2I(i,j,k)
      krome_x(krome_idx_NH2) = NH2I(i,j,k)
      krome_x(krome_idx_CH3) = CH3I(i,j,k)
      krome_x(krome_idx_CH4) = CH4I(i,j,k)
      krome_x(krome_idx_N) = NI(i,j,k)
      krome_x(krome_idx_NH) = NHI(i,j,k)
      krome_x(krome_idx_HE) = HeI(i,j,k)
      krome_x(krome_idx_HNO) = HNOI(i,j,k)
      krome_x(krome_idx_CH3OH) = CH3OHI(i,j,k)
      krome_x(krome_idx_CO2) = CO2I(i,j,k)
      krome_x(krome_idx_H2CN) = H2CNI(i,j,k)
      krome_x(krome_idx_HNCO) = HNCOI(i,j,k)
      krome_x(krome_idx_NO2) = NO2I(i,j,k)
      krome_x(krome_idx_O2H) = O2HI(i,j,k)
      krome_x(krome_idx_OCN) = OCNI(i,j,k)
      krome_x(krome_idx_CH3OH_DUST) = CH3OH_DUSTI(i,j,k)
      krome_x(krome_idx_HNCO_DUST) = HNCO_DUSTI(i,j,k)
      krome_x(krome_idx_H2CO_DUST) = H2CO_DUSTI(i,j,k)
      krome_x(krome_idx_CH4_DUST) = CH4_DUSTI(i,j,k)
      krome_x(krome_idx_CO_DUST) = CO_DUSTI(i,j,k)
      krome_x(krome_idx_H2O_DUST) = H2O_DUSTI(i,j,k)
      krome_x(krome_idx_NO_DUST) = NO_DUSTI(i,j,k)
      krome_x(krome_idx_CO2_DUST) = CO2_DUSTI(i,j,k)
      krome_x(krome_idx_N2_DUST) = N2_DUSTI(i,j,k)
      krome_x(krome_idx_HCN_DUST) = HCN_DUSTI(i,j,k)
      krome_x(krome_idx_NH3_DUST) = NH3_DUSTI(i,j,k)
      krome_x(krome_idx_O2_DUST) = O2_DUSTI(i,j,k)
      krome_x(krome_idx_NO2_DUST) = NO2_DUSTI(i,j,k)
      krome_x(krome_idx_HNO_DUST) = HNO_DUSTI(i,j,k)
      krome_x(krome_idx_O2H_DUST) = O2H_DUSTI(i,j,k)
      krome_x(krome_idx_H2CN_DUST) = H2CN_DUSTI(i,j,k)
      krome_x(krome_idx_MG_DUST) = MG_DUSTI(i,j,k)
      krome_x(krome_idx_HNC_DUST) = HNC_DUSTI(i,j,k)
      krome_x(krome_idx_E_DUST) = E_DUSTI(i,j,k)
      krome_x(krome_idx_HCOj) = HCOII(i,j,k)
      krome_x(krome_idx_Hj) = HII(i,j,k)
      krome_x(krome_idx_HOCj) = HOCII(i,j,k)
      krome_x(krome_idx_Cj) = CII(i,j,k)
      krome_x(krome_idx_CH2j) = CH2II(i,j,k)
      krome_x(krome_idx_CHj) = CHII(i,j,k)
      krome_x(krome_idx_H2COj) = H2COII(i,j,k)
      krome_x(krome_idx_MGj) = MGII(i,j,k)
      krome_x(krome_idx_NH3j) = NH3II(i,j,k)
      krome_x(krome_idx_NOj) = NOII(i,j,k)
      krome_x(krome_idx_CNj) = CNII(i,j,k)
      krome_x(krome_idx_COj) = COII(i,j,k)
      krome_x(krome_idx_N2j) = N2II(i,j,k)
      krome_x(krome_idx_O2j) = O2II(i,j,k)
      krome_x(krome_idx_H2Oj) = H2OII(i,j,k)
      krome_x(krome_idx_NH2j) = NH2II(i,j,k)
      krome_x(krome_idx_Oj) = OII(i,j,k)
      krome_x(krome_idx_OHj) = OHII(i,j,k)
      krome_x(krome_idx_CH3j) = CH3II(i,j,k)
      krome_x(krome_idx_CH4j) = CH4II(i,j,k)
      krome_x(krome_idx_Nj) = NII(i,j,k)
      krome_x(krome_idx_HCNj) = HCNII(i,j,k)
      krome_x(krome_idx_NHj) = NHII(i,j,k)
      krome_x(krome_idx_H2j) = H2II(i,j,k)
      krome_x(krome_idx_HEj) = HeII(i,j,k)
      krome_x(krome_idx_HNOj) = HNOII(i,j,k)
      krome_x(krome_idx_H2NOj) = H2NOII(i,j,k)
      krome_x(krome_idx_H3j) = H3II(i,j,k)
      krome_x(krome_idx_H3COj) = H3COII(i,j,k)
      krome_x(krome_idx_H3Oj) = H3OII(i,j,k)
      krome_x(krome_idx_HCNHj) = HCNHII(i,j,k)
      krome_x(krome_idx_HCO2j) = HCO2II(i,j,k)
      krome_x(krome_idx_HEHj) = HeHII(i,j,k)
      krome_x(krome_idx_N2Hj) = N2HII(i,j,k)
      krome_x(krome_idx_O2Hj) = O2HII(i,j,k)

      call krome_conserveLin_x(krome_x(:),ref(:))

      De(i,j,k) = krome_x(krome_idx_E)
      CHI(i,j,k) = krome_x(krome_idx_CH)
      OI(i,j,k) = krome_x(krome_idx_O)
      HNCI(i,j,k) = krome_x(krome_idx_HNC)
      HCNI(i,j,k) = krome_x(krome_idx_HCN)
      H2I(i,j,k) = krome_x(krome_idx_H2)
      CI(i,j,k) = krome_x(krome_idx_C)
      HI(i,j,k) = krome_x(krome_idx_H)
      H2OI(i,j,k) = krome_x(krome_idx_H2O)
      OHI(i,j,k) = krome_x(krome_idx_OH)
      O2I(i,j,k) = krome_x(krome_idx_O2)
      CH2I(i,j,k) = krome_x(krome_idx_CH2)
      H2COI(i,j,k) = krome_x(krome_idx_H2CO)
      HCOI(i,j,k) = krome_x(krome_idx_HCO)
      MGI(i,j,k) = krome_x(krome_idx_MG)
      NH3I(i,j,k) = krome_x(krome_idx_NH3)
      NOI(i,j,k) = krome_x(krome_idx_NO)
      CNI(i,j,k) = krome_x(krome_idx_CN)
      COI(i,j,k) = krome_x(krome_idx_CO)
      N2I(i,j,k) = krome_x(krome_idx_N2)
      NH2I(i,j,k) = krome_x(krome_idx_NH2)
      CH3I(i,j,k) = krome_x(krome_idx_CH3)
      CH4I(i,j,k) = krome_x(krome_idx_CH4)
      NI(i,j,k) = krome_x(krome_idx_N)
      NHI(i,j,k) = krome_x(krome_idx_NH)
      HeI(i,j,k) = krome_x(krome_idx_HE)
      HNOI(i,j,k) = krome_x(krome_idx_HNO)
      CH3OHI(i,j,k) = krome_x(krome_idx_CH3OH)
      CO2I(i,j,k) = krome_x(krome_idx_CO2)
      H2CNI(i,j,k) = krome_x(krome_idx_H2CN)
      HNCOI(i,j,k) = krome_x(krome_idx_HNCO)
      NO2I(i,j,k) = krome_x(krome_idx_NO2)
      O2HI(i,j,k) = krome_x(krome_idx_O2H)
      OCNI(i,j,k) = krome_x(krome_idx_OCN)
      CH3OH_DUSTI(i,j,k) = krome_x(krome_idx_CH3OH_DUST)
      HNCO_DUSTI(i,j,k) = krome_x(krome_idx_HNCO_DUST)
      H2CO_DUSTI(i,j,k) = krome_x(krome_idx_H2CO_DUST)
      CH4_DUSTI(i,j,k) = krome_x(krome_idx_CH4_DUST)
      CO_DUSTI(i,j,k) = krome_x(krome_idx_CO_DUST)
      H2O_DUSTI(i,j,k) = krome_x(krome_idx_H2O_DUST)
      NO_DUSTI(i,j,k) = krome_x(krome_idx_NO_DUST)
      CO2_DUSTI(i,j,k) = krome_x(krome_idx_CO2_DUST)
      N2_DUSTI(i,j,k) = krome_x(krome_idx_N2_DUST)
      HCN_DUSTI(i,j,k) = krome_x(krome_idx_HCN_DUST)
      NH3_DUSTI(i,j,k) = krome_x(krome_idx_NH3_DUST)
      O2_DUSTI(i,j,k) = krome_x(krome_idx_O2_DUST)
      NO2_DUSTI(i,j,k) = krome_x(krome_idx_NO2_DUST)
      HNO_DUSTI(i,j,k) = krome_x(krome_idx_HNO_DUST)
      O2H_DUSTI(i,j,k) = krome_x(krome_idx_O2H_DUST)
      H2CN_DUSTI(i,j,k) = krome_x(krome_idx_H2CN_DUST)
      MG_DUSTI(i,j,k) = krome_x(krome_idx_MG_DUST)
      HNC_DUSTI(i,j,k) = krome_x(krome_idx_HNC_DUST)
      E_DUSTI(i,j,k) = krome_x(krome_idx_E_DUST)
      HCOII(i,j,k) = krome_x(krome_idx_HCOj)
      HII(i,j,k) = krome_x(krome_idx_Hj)
      HOCII(i,j,k) = krome_x(krome_idx_HOCj)
      CII(i,j,k) = krome_x(krome_idx_Cj)
      CH2II(i,j,k) = krome_x(krome_idx_CH2j)
      CHII(i,j,k) = krome_x(krome_idx_CHj)
      H2COII(i,j,k) = krome_x(krome_idx_H2COj)
      MGII(i,j,k) = krome_x(krome_idx_MGj)
      NH3II(i,j,k) = krome_x(krome_idx_NH3j)
      NOII(i,j,k) = krome_x(krome_idx_NOj)
      CNII(i,j,k) = krome_x(krome_idx_CNj)
      COII(i,j,k) = krome_x(krome_idx_COj)
      N2II(i,j,k) = krome_x(krome_idx_N2j)
      O2II(i,j,k) = krome_x(krome_idx_O2j)
      H2OII(i,j,k) = krome_x(krome_idx_H2Oj)
      NH2II(i,j,k) = krome_x(krome_idx_NH2j)
      OII(i,j,k) = krome_x(krome_idx_Oj)
      OHII(i,j,k) = krome_x(krome_idx_OHj)
      CH3II(i,j,k) = krome_x(krome_idx_CH3j)
      CH4II(i,j,k) = krome_x(krome_idx_CH4j)
      NII(i,j,k) = krome_x(krome_idx_Nj)
      HCNII(i,j,k) = krome_x(krome_idx_HCNj)
      NHII(i,j,k) = krome_x(krome_idx_NHj)
      H2II(i,j,k) = krome_x(krome_idx_H2j)
      HeII(i,j,k) = krome_x(krome_idx_HEj)
      HNOII(i,j,k) = krome_x(krome_idx_HNOj)
      H2NOII(i,j,k) = krome_x(krome_idx_H2NOj)
      H3II(i,j,k) = krome_x(krome_idx_H3j)
      H3COII(i,j,k) = krome_x(krome_idx_H3COj)
      H3OII(i,j,k) = krome_x(krome_idx_H3Oj)
      HCNHII(i,j,k) = krome_x(krome_idx_HCNHj)
      HCO2II(i,j,k) = krome_x(krome_idx_HCO2j)
      HeHII(i,j,k) = krome_x(krome_idx_HEHj)
      N2HII(i,j,k) = krome_x(krome_idx_N2Hj)
      O2HII(i,j,k) = krome_x(krome_idx_O2Hj)

    end do
  end do
end do

end subroutine krome_renormref
