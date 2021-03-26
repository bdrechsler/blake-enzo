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








subroutine krome_initref(d, urho, &
   De, CHI, OI, HNCI, &
   HCNI, H2I, CI, HI, &
   H2OI, OHI, O2I, CH2I, &
   H2COI, HCOI, MGI, NH3I, &
   NOI, SII, SIC2I, SIC3I, &
   SICI, SIH2I, SIH3I, CNI, &
   COI, N2I, NH2I, CH3I, &
   CH4I, NI, NHI, SIH4I, &
   SIHI, SIOI, HeI, HNOI, &
   CH3OHI, CO2I, H2CNI, H2SIOI, &
   HNCOI, NO2I, O2HI, OCNI, &
   CH3OH_DUSTI, HNCO_DUSTI, H2CO_DUSTI, SIH4_DUSTI, &
   H2SIO_DUSTI, SIC_DUSTI, SIC2_DUSTI, SIC3_DUSTI, &
   CH4_DUSTI, CO_DUSTI, H2O_DUSTI, NO_DUSTI, &
   CO2_DUSTI, N2_DUSTI, HCN_DUSTI, NH3_DUSTI, &
   O2_DUSTI, NO2_DUSTI, HNO_DUSTI, O2H_DUSTI, &
   H2CN_DUSTI, MG_DUSTI, HNC_DUSTI, E_DUSTI, &
   SIO_DUSTI, HCOII, HII, HOCII, &
   CII, CH2II, CHII, H2COII, &
   MGII, NH3II, NOII, SIII, &
   SIC2II, SIC3II, SICII, SIH2II, &
   SIH3II, CNII, COII, N2II, &
   O2II, H2OII, NH2II, OII, &
   OHII, CH3II, CH4II, NII, &
   HCNII, NHII, SIH4II, SIHII, &
   SIOII, H2II, HeII, HNOII, &
   H2NOII, H3II, H3COII, H3OII, &
   HCNHII, HCO2II, HeHII, N2HII, &
   O2HII, SIH5II, SIOHII, &
   in, jn, kn, idim, &
   is, js, ks, ie, je, ke, &
   ATOM_C, ATOM_H, ATOM_MG, &
   ATOM_O, ATOM_N, ATOM_SI, ATOM_HE)

   use krome_main
   use krome_user
   use krome_constants
   use krome_commons
 
   implicit none
 
   real*8,parameter::mh=p_mass !mass_h
   real*8::urho,dom,krome_x(krome_nmols),ref(krome_natoms)
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
   real*8::SII(in,jn,kn)
   real*8::SIC2I(in,jn,kn)
   real*8::SIC3I(in,jn,kn)
   real*8::SICI(in,jn,kn)
   real*8::SIH2I(in,jn,kn)
   real*8::SIH3I(in,jn,kn)
   real*8::CNI(in,jn,kn)
   real*8::COI(in,jn,kn)
   real*8::N2I(in,jn,kn)
   real*8::NH2I(in,jn,kn)
   real*8::CH3I(in,jn,kn)
   real*8::CH4I(in,jn,kn)
   real*8::NI(in,jn,kn)
   real*8::NHI(in,jn,kn)
   real*8::SIH4I(in,jn,kn)
   real*8::SIHI(in,jn,kn)
   real*8::SIOI(in,jn,kn)
   real*8::HeI(in,jn,kn)
   real*8::HNOI(in,jn,kn)
   real*8::CH3OHI(in,jn,kn)
   real*8::CO2I(in,jn,kn)
   real*8::H2CNI(in,jn,kn)
   real*8::H2SIOI(in,jn,kn)
   real*8::HNCOI(in,jn,kn)
   real*8::NO2I(in,jn,kn)
   real*8::O2HI(in,jn,kn)
   real*8::OCNI(in,jn,kn)
   real*8::CH3OH_DUSTI(in,jn,kn)
   real*8::HNCO_DUSTI(in,jn,kn)
   real*8::H2CO_DUSTI(in,jn,kn)
   real*8::SIH4_DUSTI(in,jn,kn)
   real*8::H2SIO_DUSTI(in,jn,kn)
   real*8::SIC_DUSTI(in,jn,kn)
   real*8::SIC2_DUSTI(in,jn,kn)
   real*8::SIC3_DUSTI(in,jn,kn)
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
   real*8::SIO_DUSTI(in,jn,kn)
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
   real*8::SIII(in,jn,kn)
   real*8::SIC2II(in,jn,kn)
   real*8::SIC3II(in,jn,kn)
   real*8::SICII(in,jn,kn)
   real*8::SIH2II(in,jn,kn)
   real*8::SIH3II(in,jn,kn)
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
   real*8::SIH4II(in,jn,kn)
   real*8::SIHII(in,jn,kn)
   real*8::SIOII(in,jn,kn)
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
   real*8::SIH5II(in,jn,kn)
   real*8::SIOHII(in,jn,kn)
 
   real*8::ATOM_H(in,jn,kn)
   real*8::ATOM_HE(in,jn,kn)
   real*8::ATOM_C(in,jn,kn)
   real*8::ATOM_N(in,jn,kn)
   real*8::ATOM_O(in,jn,kn)
   real*8::ATOM_MG(in,jn,kn)
   real*8::ATOM_SI(in,jn,kn)

   do k = ks+1, ke+1
      do j = js+1, je+1
        do i = is+1, ie+1

          !convert to number densities
          krome_x(krome_idx_E) = De(i,j,k) * dom
          krome_x(krome_idx_CH) = CHI(i,j,k) * dom * 0.0769230769231d0
          krome_x(krome_idx_O) = OI(i,j,k) * dom * 0.0625d0
          krome_x(krome_idx_HNC) = HNCI(i,j,k) * dom * 0.037037037037d0
          krome_x(krome_idx_HCN) = HCNI(i,j,k) * dom * 0.037037037037d0
          krome_x(krome_idx_H2) = H2I(i,j,k) * dom * 0.5d0
          krome_x(krome_idx_C) = CI(i,j,k) * dom * 0.0833333333333d0
          krome_x(krome_idx_H) = HI(i,j,k) * dom
          krome_x(krome_idx_H2O) = H2OI(i,j,k) * dom * 0.0555555555556d0
          krome_x(krome_idx_OH) = OHI(i,j,k) * dom * 0.0588235294118d0
          krome_x(krome_idx_O2) = O2I(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_CH2) = CH2I(i,j,k) * dom * 0.0714285714286d0
          krome_x(krome_idx_H2CO) = H2COI(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_HCO) = HCOI(i,j,k) * dom * 0.0344827586207d0
          krome_x(krome_idx_MG) = MGI(i,j,k) * dom * 0.0416666666667d0
          krome_x(krome_idx_NH3) = NH3I(i,j,k) * dom * 0.0588235294118d0
          krome_x(krome_idx_NO) = NOI(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_SI) = SII(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_SIC2) = SIC2I(i,j,k) * dom * 0.0192307692308d0
          krome_x(krome_idx_SIC3) = SIC3I(i,j,k) * dom * 0.015625d0
          krome_x(krome_idx_SIC) = SICI(i,j,k) * dom * 0.025d0
          krome_x(krome_idx_SIH2) = SIH2I(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_SIH3) = SIH3I(i,j,k) * dom * 0.0322580645161d0
          krome_x(krome_idx_CN) = CNI(i,j,k) * dom * 0.0384615384615d0
          krome_x(krome_idx_CO) = COI(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_N2) = N2I(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_NH2) = NH2I(i,j,k) * dom * 0.0625d0
          krome_x(krome_idx_CH3) = CH3I(i,j,k) * dom * 0.0666666666667d0
          krome_x(krome_idx_CH4) = CH4I(i,j,k) * dom * 0.0625d0
          krome_x(krome_idx_N) = NI(i,j,k) * dom * 0.0714285714286d0
          krome_x(krome_idx_NH) = NHI(i,j,k) * dom * 0.0666666666667d0
          krome_x(krome_idx_SIH4) = SIH4I(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_SIH) = SIHI(i,j,k) * dom * 0.0344827586207d0
          krome_x(krome_idx_SIO) = SIOI(i,j,k) * dom * 0.0227272727273d0
          krome_x(krome_idx_HE) = HeI(i,j,k) * dom * 0.25d0
          krome_x(krome_idx_HNO) = HNOI(i,j,k) * dom * 0.0322580645161d0
          krome_x(krome_idx_CH3OH) = CH3OHI(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_CO2) = CO2I(i,j,k) * dom * 0.0227272727273d0
          krome_x(krome_idx_H2CN) = H2CNI(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_H2SIO) = H2SIOI(i,j,k) * dom * 0.0217391304348d0
          krome_x(krome_idx_HNCO) = HNCOI(i,j,k) * dom * 0.0232558139535d0
          krome_x(krome_idx_NO2) = NO2I(i,j,k) * dom * 0.0217391304348d0
          krome_x(krome_idx_O2H) = O2HI(i,j,k) * dom * 0.030303030303d0
          krome_x(krome_idx_OCN) = OCNI(i,j,k) * dom * 0.0238095238095d0
          krome_x(krome_idx_CH3OH_DUST) = CH3OH_DUSTI(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_HNCO_DUST) = HNCO_DUSTI(i,j,k) * dom * 0.0232558139535d0
          krome_x(krome_idx_H2CO_DUST) = H2CO_DUSTI(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_SIH4_DUST) = SIH4_DUSTI(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_H2SIO_DUST) = H2SIO_DUSTI(i,j,k) * dom * 0.0217391304348d0
          krome_x(krome_idx_SIC_DUST) = SIC_DUSTI(i,j,k) * dom * 0.025d0
          krome_x(krome_idx_SIC2_DUST) = SIC2_DUSTI(i,j,k) * dom * 0.0192307692308d0
          krome_x(krome_idx_SIC3_DUST) = SIC3_DUSTI(i,j,k) * dom * 0.015625d0
          krome_x(krome_idx_CH4_DUST) = CH4_DUSTI(i,j,k) * dom * 0.0625d0
          krome_x(krome_idx_CO_DUST) = CO_DUSTI(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_H2O_DUST) = H2O_DUSTI(i,j,k) * dom * 0.0555555555556d0
          krome_x(krome_idx_NO_DUST) = NO_DUSTI(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_CO2_DUST) = CO2_DUSTI(i,j,k) * dom * 0.0227272727273d0
          krome_x(krome_idx_N2_DUST) = N2_DUSTI(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_HCN_DUST) = HCN_DUSTI(i,j,k) * dom * 0.037037037037d0
          krome_x(krome_idx_NH3_DUST) = NH3_DUSTI(i,j,k) * dom * 0.0588235294118d0
          krome_x(krome_idx_O2_DUST) = O2_DUSTI(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_NO2_DUST) = NO2_DUSTI(i,j,k) * dom * 0.0217391304348d0
          krome_x(krome_idx_HNO_DUST) = HNO_DUSTI(i,j,k) * dom * 0.0322580645161d0
          krome_x(krome_idx_O2H_DUST) = O2H_DUSTI(i,j,k) * dom * 0.030303030303d0
          krome_x(krome_idx_H2CN_DUST) = H2CN_DUSTI(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_MG_DUST) = MG_DUSTI(i,j,k) * dom * 0.0416666666667d0
          krome_x(krome_idx_HNC_DUST) = HNC_DUSTI(i,j,k) * dom * 0.037037037037d0
          krome_x(krome_idx_E_DUST) = E_DUSTI(i,j,k) * dom
          krome_x(krome_idx_SIO_DUST) = SIO_DUSTI(i,j,k) * dom * 0.0227272727273d0
          krome_x(krome_idx_HCOj) = HCOII(i,j,k) * dom * 0.0344827586207d0
          krome_x(krome_idx_Hj) = HII(i,j,k) * dom
          krome_x(krome_idx_HOCj) = HOCII(i,j,k) * dom * 0.0344827586207d0
          krome_x(krome_idx_Cj) = CII(i,j,k) * dom * 0.0833333333333d0
          krome_x(krome_idx_CH2j) = CH2II(i,j,k) * dom * 0.0714285714286d0
          krome_x(krome_idx_CHj) = CHII(i,j,k) * dom * 0.0769230769231d0
          krome_x(krome_idx_H2COj) = H2COII(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_MGj) = MGII(i,j,k) * dom * 0.0416666666667d0
          krome_x(krome_idx_NH3j) = NH3II(i,j,k) * dom * 0.0588235294118d0
          krome_x(krome_idx_NOj) = NOII(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_SIj) = SIII(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_SIC2j) = SIC2II(i,j,k) * dom * 0.0192307692308d0
          krome_x(krome_idx_SIC3j) = SIC3II(i,j,k) * dom * 0.015625d0
          krome_x(krome_idx_SICj) = SICII(i,j,k) * dom * 0.025d0
          krome_x(krome_idx_SIH2j) = SIH2II(i,j,k) * dom * 0.0333333333333d0
          krome_x(krome_idx_SIH3j) = SIH3II(i,j,k) * dom * 0.0322580645161d0
          krome_x(krome_idx_CNj) = CNII(i,j,k) * dom * 0.0384615384615d0
          krome_x(krome_idx_COj) = COII(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_N2j) = N2II(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_O2j) = O2II(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_H2Oj) = H2OII(i,j,k) * dom * 0.0555555555556d0
          krome_x(krome_idx_NH2j) = NH2II(i,j,k) * dom * 0.0625d0
          krome_x(krome_idx_Oj) = OII(i,j,k) * dom * 0.0625d0
          krome_x(krome_idx_OHj) = OHII(i,j,k) * dom * 0.0588235294118d0
          krome_x(krome_idx_CH3j) = CH3II(i,j,k) * dom * 0.0666666666667d0
          krome_x(krome_idx_CH4j) = CH4II(i,j,k) * dom * 0.0625d0
          krome_x(krome_idx_Nj) = NII(i,j,k) * dom * 0.0714285714286d0
          krome_x(krome_idx_HCNj) = HCNII(i,j,k) * dom * 0.037037037037d0
          krome_x(krome_idx_NHj) = NHII(i,j,k) * dom * 0.0666666666667d0
          krome_x(krome_idx_SIH4j) = SIH4II(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_SIHj) = SIHII(i,j,k) * dom * 0.0344827586207d0
          krome_x(krome_idx_SIOj) = SIOII(i,j,k) * dom * 0.0227272727273d0
          krome_x(krome_idx_H2j) = H2II(i,j,k) * dom * 0.5d0
          krome_x(krome_idx_HEj) = HeII(i,j,k) * dom * 0.25d0
          krome_x(krome_idx_HNOj) = HNOII(i,j,k) * dom * 0.0322580645161d0
          krome_x(krome_idx_H2NOj) = H2NOII(i,j,k) * dom * 0.03125d0
          krome_x(krome_idx_H3j) = H3II(i,j,k) * dom * 0.333333333333d0
          krome_x(krome_idx_H3COj) = H3COII(i,j,k) * dom * 0.0322580645161d0
          krome_x(krome_idx_H3Oj) = H3OII(i,j,k) * dom * 0.0526315789474d0
          krome_x(krome_idx_HCNHj) = HCNHII(i,j,k) * dom * 0.0357142857143d0
          krome_x(krome_idx_HCO2j) = HCO2II(i,j,k) * dom * 0.0222222222222d0
          krome_x(krome_idx_HEHj) = HeHII(i,j,k) * dom * 0.2d0
          krome_x(krome_idx_N2Hj) = N2HII(i,j,k) * dom * 0.0344827586207d0
          krome_x(krome_idx_O2Hj) = O2HII(i,j,k) * dom * 0.030303030303d0
          krome_x(krome_idx_SIH5j) = SIH5II(i,j,k) * dom * 0.030303030303d0
          krome_x(krome_idx_SIOHj) = SIOHII(i,j,k) * dom * 0.0222222222222d0

          ref(:) = krome_conserveLinGetRef_x(krome_x(:))
      
          ATOM_H(i,j,k) = ref(idx_atom_H)
          ATOM_HE(i,j,k) = ref(idx_atom_HE)
          ATOM_C(i,j,k) = ref(idx_atom_C)
          ATOM_N(i,j,k) = ref(idx_atom_N)
          ATOM_O(i,j,k) = ref(idx_atom_O)
          ATOM_MG(i,j,k) = ref(idx_atom_MG)
          ATOM_SI(i,j,k) = ref(idx_atom_SI)
        end do
      end do
   end do

end subroutine krome_initref

subroutine krome_renormref(d, urho, &
   De, CHI, OI, HNCI, &
   HCNI, H2I, CI, HI, &
   H2OI, OHI, O2I, CH2I, &
   H2COI, HCOI, MGI, NH3I, &
   NOI, SII, SIC2I, SIC3I, &
   SICI, SIH2I, SIH3I, CNI, &
   COI, N2I, NH2I, CH3I, &
   CH4I, NI, NHI, SIH4I, &
   SIHI, SIOI, HeI, HNOI, &
   CH3OHI, CO2I, H2CNI, H2SIOI, &
   HNCOI, NO2I, O2HI, OCNI, &
   CH3OH_DUSTI, HNCO_DUSTI, H2CO_DUSTI, SIH4_DUSTI, &
   H2SIO_DUSTI, SIC_DUSTI, SIC2_DUSTI, SIC3_DUSTI, &
   CH4_DUSTI, CO_DUSTI, H2O_DUSTI, NO_DUSTI, &
   CO2_DUSTI, N2_DUSTI, HCN_DUSTI, NH3_DUSTI, &
   O2_DUSTI, NO2_DUSTI, HNO_DUSTI, O2H_DUSTI, &
   H2CN_DUSTI, MG_DUSTI, HNC_DUSTI, E_DUSTI, &
   SIO_DUSTI, HCOII, HII, HOCII, &
   CII, CH2II, CHII, H2COII, &
   MGII, NH3II, NOII, SIII, &
   SIC2II, SIC3II, SICII, SIH2II, &
   SIH3II, CNII, COII, N2II, &
   O2II, H2OII, NH2II, OII, &
   OHII, CH3II, CH4II, NII, &
   HCNII, NHII, SIH4II, SIHII, &
   SIOII, H2II, HeII, HNOII, &
   H2NOII, H3II, H3COII, H3OII, &
   HCNHII, HCO2II, HeHII, N2HII, &
   O2HII, SIH5II, SIOHII, &
   in, jn, kn, idim, &
   is, js, ks, ie, je, ke, &
   ATOM_C, ATOM_H, ATOM_MG, &
   ATOM_O, ATOM_N, ATOM_SI, ATOM_HE)


!     USE KROME
use krome_main
use krome_user
use krome_constants
use krome_commons

implicit none

real*8,parameter::mh=p_mass !mass_h
real*8::urho,dom,idom,krome_x(krome_nmols),ref(krome_natoms)
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
real*8::SII(in,jn,kn)
real*8::SIC2I(in,jn,kn)
real*8::SIC3I(in,jn,kn)
real*8::SICI(in,jn,kn)
real*8::SIH2I(in,jn,kn)
real*8::SIH3I(in,jn,kn)
real*8::CNI(in,jn,kn)
real*8::COI(in,jn,kn)
real*8::N2I(in,jn,kn)
real*8::NH2I(in,jn,kn)
real*8::CH3I(in,jn,kn)
real*8::CH4I(in,jn,kn)
real*8::NI(in,jn,kn)
real*8::NHI(in,jn,kn)
real*8::SIH4I(in,jn,kn)
real*8::SIHI(in,jn,kn)
real*8::SIOI(in,jn,kn)
real*8::HeI(in,jn,kn)
real*8::HNOI(in,jn,kn)
real*8::CH3OHI(in,jn,kn)
real*8::CO2I(in,jn,kn)
real*8::H2CNI(in,jn,kn)
real*8::H2SIOI(in,jn,kn)
real*8::HNCOI(in,jn,kn)
real*8::NO2I(in,jn,kn)
real*8::O2HI(in,jn,kn)
real*8::OCNI(in,jn,kn)
real*8::CH3OH_DUSTI(in,jn,kn)
real*8::HNCO_DUSTI(in,jn,kn)
real*8::H2CO_DUSTI(in,jn,kn)
real*8::SIH4_DUSTI(in,jn,kn)
real*8::H2SIO_DUSTI(in,jn,kn)
real*8::SIC_DUSTI(in,jn,kn)
real*8::SIC2_DUSTI(in,jn,kn)
real*8::SIC3_DUSTI(in,jn,kn)
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
real*8::SIO_DUSTI(in,jn,kn)
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
real*8::SIII(in,jn,kn)
real*8::SIC2II(in,jn,kn)
real*8::SIC3II(in,jn,kn)
real*8::SICII(in,jn,kn)
real*8::SIH2II(in,jn,kn)
real*8::SIH3II(in,jn,kn)
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
real*8::SIH4II(in,jn,kn)
real*8::SIHII(in,jn,kn)
real*8::SIOII(in,jn,kn)
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
real*8::SIH5II(in,jn,kn)
real*8::SIOHII(in,jn,kn)

real*8::ATOM_H(in,jn,kn)
real*8::ATOM_HE(in,jn,kn)
real*8::ATOM_C(in,jn,kn)
real*8::ATOM_N(in,jn,kn)
real*8::ATOM_O(in,jn,kn)
real*8::ATOM_MG(in,jn,kn)
real*8::ATOM_SI(in,jn,kn)

!******************************

!set units
dom = urho/mh

!loop over zones
do k = ks+1, ke+1
  do j = js+1, je+1
    do i = is+1, ie+1

      ref(idx_atom_H) = ATOM_H(i,j,k)
      ref(idx_atom_HE) = ATOM_HE(i,j,k)
      ref(idx_atom_C) = ATOM_C(i,j,k)
      ref(idx_atom_N) = ATOM_N(i,j,k)
      ref(idx_atom_O) = ATOM_O(i,j,k)
      ref(idx_atom_SI) = ATOM_SI(i,j,k)
      ref(idx_atom_MG) = ATOM_MG(i,j,k)

      krome_x(krome_idx_E) = De(i,j,k) * dom
      krome_x(krome_idx_CH) = CHI(i,j,k) * dom * 0.0769230769231d0
      krome_x(krome_idx_O) = OI(i,j,k) * dom * 0.0625d0
      krome_x(krome_idx_HNC) = HNCI(i,j,k) * dom * 0.037037037037d0
      krome_x(krome_idx_HCN) = HCNI(i,j,k) * dom * 0.037037037037d0
      krome_x(krome_idx_H2) = H2I(i,j,k) * dom * 0.5d0
      krome_x(krome_idx_C) = CI(i,j,k) * dom * 0.0833333333333d0
      krome_x(krome_idx_H) = HI(i,j,k) * dom
      krome_x(krome_idx_H2O) = H2OI(i,j,k) * dom * 0.0555555555556d0
      krome_x(krome_idx_OH) = OHI(i,j,k) * dom * 0.0588235294118d0
      krome_x(krome_idx_O2) = O2I(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_CH2) = CH2I(i,j,k) * dom * 0.0714285714286d0
      krome_x(krome_idx_H2CO) = H2COI(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_HCO) = HCOI(i,j,k) * dom * 0.0344827586207d0
      krome_x(krome_idx_MG) = MGI(i,j,k) * dom * 0.0416666666667d0
      krome_x(krome_idx_NH3) = NH3I(i,j,k) * dom * 0.0588235294118d0
      krome_x(krome_idx_NO) = NOI(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_SI) = SII(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_SIC2) = SIC2I(i,j,k) * dom * 0.0192307692308d0
      krome_x(krome_idx_SIC3) = SIC3I(i,j,k) * dom * 0.015625d0
      krome_x(krome_idx_SIC) = SICI(i,j,k) * dom * 0.025d0
      krome_x(krome_idx_SIH2) = SIH2I(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_SIH3) = SIH3I(i,j,k) * dom * 0.0322580645161d0
      krome_x(krome_idx_CN) = CNI(i,j,k) * dom * 0.0384615384615d0
      krome_x(krome_idx_CO) = COI(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_N2) = N2I(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_NH2) = NH2I(i,j,k) * dom * 0.0625d0
      krome_x(krome_idx_CH3) = CH3I(i,j,k) * dom * 0.0666666666667d0
      krome_x(krome_idx_CH4) = CH4I(i,j,k) * dom * 0.0625d0
      krome_x(krome_idx_N) = NI(i,j,k) * dom * 0.0714285714286d0
      krome_x(krome_idx_NH) = NHI(i,j,k) * dom * 0.0666666666667d0
      krome_x(krome_idx_SIH4) = SIH4I(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_SIH) = SIHI(i,j,k) * dom * 0.0344827586207d0
      krome_x(krome_idx_SIO) = SIOI(i,j,k) * dom * 0.0227272727273d0
      krome_x(krome_idx_HE) = HeI(i,j,k) * dom * 0.25d0
      krome_x(krome_idx_HNO) = HNOI(i,j,k) * dom * 0.0322580645161d0
      krome_x(krome_idx_CH3OH) = CH3OHI(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_CO2) = CO2I(i,j,k) * dom * 0.0227272727273d0
      krome_x(krome_idx_H2CN) = H2CNI(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_H2SIO) = H2SIOI(i,j,k) * dom * 0.0217391304348d0
      krome_x(krome_idx_HNCO) = HNCOI(i,j,k) * dom * 0.0232558139535d0
      krome_x(krome_idx_NO2) = NO2I(i,j,k) * dom * 0.0217391304348d0
      krome_x(krome_idx_O2H) = O2HI(i,j,k) * dom * 0.030303030303d0
      krome_x(krome_idx_OCN) = OCNI(i,j,k) * dom * 0.0238095238095d0
      krome_x(krome_idx_CH3OH_DUST) = CH3OH_DUSTI(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_HNCO_DUST) = HNCO_DUSTI(i,j,k) * dom * 0.0232558139535d0
      krome_x(krome_idx_H2CO_DUST) = H2CO_DUSTI(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_SIH4_DUST) = SIH4_DUSTI(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_H2SIO_DUST) = H2SIO_DUSTI(i,j,k) * dom * 0.0217391304348d0
      krome_x(krome_idx_SIC_DUST) = SIC_DUSTI(i,j,k) * dom * 0.025d0
      krome_x(krome_idx_SIC2_DUST) = SIC2_DUSTI(i,j,k) * dom * 0.0192307692308d0
      krome_x(krome_idx_SIC3_DUST) = SIC3_DUSTI(i,j,k) * dom * 0.015625d0
      krome_x(krome_idx_CH4_DUST) = CH4_DUSTI(i,j,k) * dom * 0.0625d0
      krome_x(krome_idx_CO_DUST) = CO_DUSTI(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_H2O_DUST) = H2O_DUSTI(i,j,k) * dom * 0.0555555555556d0
      krome_x(krome_idx_NO_DUST) = NO_DUSTI(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_CO2_DUST) = CO2_DUSTI(i,j,k) * dom * 0.0227272727273d0
      krome_x(krome_idx_N2_DUST) = N2_DUSTI(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_HCN_DUST) = HCN_DUSTI(i,j,k) * dom * 0.037037037037d0
      krome_x(krome_idx_NH3_DUST) = NH3_DUSTI(i,j,k) * dom * 0.0588235294118d0
      krome_x(krome_idx_O2_DUST) = O2_DUSTI(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_NO2_DUST) = NO2_DUSTI(i,j,k) * dom * 0.0217391304348d0
      krome_x(krome_idx_HNO_DUST) = HNO_DUSTI(i,j,k) * dom * 0.0322580645161d0
      krome_x(krome_idx_O2H_DUST) = O2H_DUSTI(i,j,k) * dom * 0.030303030303d0
      krome_x(krome_idx_H2CN_DUST) = H2CN_DUSTI(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_MG_DUST) = MG_DUSTI(i,j,k) * dom * 0.0416666666667d0
      krome_x(krome_idx_HNC_DUST) = HNC_DUSTI(i,j,k) * dom * 0.037037037037d0
      krome_x(krome_idx_E_DUST) = E_DUSTI(i,j,k) * dom
      krome_x(krome_idx_SIO_DUST) = SIO_DUSTI(i,j,k) * dom * 0.0227272727273d0
      krome_x(krome_idx_HCOj) = HCOII(i,j,k) * dom * 0.0344827586207d0
      krome_x(krome_idx_Hj) = HII(i,j,k) * dom
      krome_x(krome_idx_HOCj) = HOCII(i,j,k) * dom * 0.0344827586207d0
      krome_x(krome_idx_Cj) = CII(i,j,k) * dom * 0.0833333333333d0
      krome_x(krome_idx_CH2j) = CH2II(i,j,k) * dom * 0.0714285714286d0
      krome_x(krome_idx_CHj) = CHII(i,j,k) * dom * 0.0769230769231d0
      krome_x(krome_idx_H2COj) = H2COII(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_MGj) = MGII(i,j,k) * dom * 0.0416666666667d0
      krome_x(krome_idx_NH3j) = NH3II(i,j,k) * dom * 0.0588235294118d0
      krome_x(krome_idx_NOj) = NOII(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_SIj) = SIII(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_SIC2j) = SIC2II(i,j,k) * dom * 0.0192307692308d0
      krome_x(krome_idx_SIC3j) = SIC3II(i,j,k) * dom * 0.015625d0
      krome_x(krome_idx_SICj) = SICII(i,j,k) * dom * 0.025d0
      krome_x(krome_idx_SIH2j) = SIH2II(i,j,k) * dom * 0.0333333333333d0
      krome_x(krome_idx_SIH3j) = SIH3II(i,j,k) * dom * 0.0322580645161d0
      krome_x(krome_idx_CNj) = CNII(i,j,k) * dom * 0.0384615384615d0
      krome_x(krome_idx_COj) = COII(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_N2j) = N2II(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_O2j) = O2II(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_H2Oj) = H2OII(i,j,k) * dom * 0.0555555555556d0
      krome_x(krome_idx_NH2j) = NH2II(i,j,k) * dom * 0.0625d0
      krome_x(krome_idx_Oj) = OII(i,j,k) * dom * 0.0625d0
      krome_x(krome_idx_OHj) = OHII(i,j,k) * dom * 0.0588235294118d0
      krome_x(krome_idx_CH3j) = CH3II(i,j,k) * dom * 0.0666666666667d0
      krome_x(krome_idx_CH4j) = CH4II(i,j,k) * dom * 0.0625d0
      krome_x(krome_idx_Nj) = NII(i,j,k) * dom * 0.0714285714286d0
      krome_x(krome_idx_HCNj) = HCNII(i,j,k) * dom * 0.037037037037d0
      krome_x(krome_idx_NHj) = NHII(i,j,k) * dom * 0.0666666666667d0
      krome_x(krome_idx_SIH4j) = SIH4II(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_SIHj) = SIHII(i,j,k) * dom * 0.0344827586207d0
      krome_x(krome_idx_SIOj) = SIOII(i,j,k) * dom * 0.0227272727273d0
      krome_x(krome_idx_H2j) = H2II(i,j,k) * dom * 0.5d0
      krome_x(krome_idx_HEj) = HeII(i,j,k) * dom * 0.25d0
      krome_x(krome_idx_HNOj) = HNOII(i,j,k) * dom * 0.0322580645161d0
      krome_x(krome_idx_H2NOj) = H2NOII(i,j,k) * dom * 0.03125d0
      krome_x(krome_idx_H3j) = H3II(i,j,k) * dom * 0.333333333333d0
      krome_x(krome_idx_H3COj) = H3COII(i,j,k) * dom * 0.0322580645161d0
      krome_x(krome_idx_H3Oj) = H3OII(i,j,k) * dom * 0.0526315789474d0
      krome_x(krome_idx_HCNHj) = HCNHII(i,j,k) * dom * 0.0357142857143d0
      krome_x(krome_idx_HCO2j) = HCO2II(i,j,k) * dom * 0.0222222222222d0
      krome_x(krome_idx_HEHj) = HeHII(i,j,k) * dom * 0.2d0
      krome_x(krome_idx_N2Hj) = N2HII(i,j,k) * dom * 0.0344827586207d0
      krome_x(krome_idx_O2Hj) = O2HII(i,j,k) * dom * 0.030303030303d0
      krome_x(krome_idx_SIH5j) = SIH5II(i,j,k) * dom * 0.030303030303d0
      krome_x(krome_idx_SIOHj) = SIOHII(i,j,k) * dom * 0.0222222222222d0


      call krome_conserveLin_x(krome_x(:),ref(:))

      idom = 1.d0/dom
      !convert back to code units
      De(i,j,k) = krome_x(krome_idx_E) * idom
      CHI(i,j,k) = krome_x(krome_idx_CH) * idom * 13d0
      OI(i,j,k) = krome_x(krome_idx_O) * idom * 16d0
      HNCI(i,j,k) = krome_x(krome_idx_HNC) * idom * 27d0
      HCNI(i,j,k) = krome_x(krome_idx_HCN) * idom * 27d0
      H2I(i,j,k) = krome_x(krome_idx_H2) * idom * 2d0
      CI(i,j,k) = krome_x(krome_idx_C) * idom * 12d0
      HI(i,j,k) = krome_x(krome_idx_H) * idom
      H2OI(i,j,k) = krome_x(krome_idx_H2O) * idom * 18d0
      OHI(i,j,k) = krome_x(krome_idx_OH) * idom * 17d0
      O2I(i,j,k) = krome_x(krome_idx_O2) * idom * 32d0
      CH2I(i,j,k) = krome_x(krome_idx_CH2) * idom * 14d0
      H2COI(i,j,k) = krome_x(krome_idx_H2CO) * idom * 30d0
      HCOI(i,j,k) = krome_x(krome_idx_HCO) * idom * 29d0
      MGI(i,j,k) = krome_x(krome_idx_MG) * idom * 24d0
      NH3I(i,j,k) = krome_x(krome_idx_NH3) * idom * 17d0
      NOI(i,j,k) = krome_x(krome_idx_NO) * idom * 30d0
      SII(i,j,k) = krome_x(krome_idx_SI) * idom * 28d0
      SIC2I(i,j,k) = krome_x(krome_idx_SIC2) * idom * 52d0
      SIC3I(i,j,k) = krome_x(krome_idx_SIC3) * idom * 64d0
      SICI(i,j,k) = krome_x(krome_idx_SIC) * idom * 40d0
      SIH2I(i,j,k) = krome_x(krome_idx_SIH2) * idom * 30d0
      SIH3I(i,j,k) = krome_x(krome_idx_SIH3) * idom * 31d0
      CNI(i,j,k) = krome_x(krome_idx_CN) * idom * 26d0
      COI(i,j,k) = krome_x(krome_idx_CO) * idom * 28d0
      N2I(i,j,k) = krome_x(krome_idx_N2) * idom * 28d0
      NH2I(i,j,k) = krome_x(krome_idx_NH2) * idom * 16d0
      CH3I(i,j,k) = krome_x(krome_idx_CH3) * idom * 15d0
      CH4I(i,j,k) = krome_x(krome_idx_CH4) * idom * 16d0
      NI(i,j,k) = krome_x(krome_idx_N) * idom * 14d0
      NHI(i,j,k) = krome_x(krome_idx_NH) * idom * 15d0
      SIH4I(i,j,k) = krome_x(krome_idx_SIH4) * idom * 32d0
      SIHI(i,j,k) = krome_x(krome_idx_SIH) * idom * 29d0
      SIOI(i,j,k) = krome_x(krome_idx_SIO) * idom * 44d0
      HeI(i,j,k) = krome_x(krome_idx_HE) * idom * 4d0
      HNOI(i,j,k) = krome_x(krome_idx_HNO) * idom * 31d0
      CH3OHI(i,j,k) = krome_x(krome_idx_CH3OH) * idom * 32d0
      CO2I(i,j,k) = krome_x(krome_idx_CO2) * idom * 44d0
      H2CNI(i,j,k) = krome_x(krome_idx_H2CN) * idom * 28d0
      H2SIOI(i,j,k) = krome_x(krome_idx_H2SIO) * idom * 46d0
      HNCOI(i,j,k) = krome_x(krome_idx_HNCO) * idom * 43d0
      NO2I(i,j,k) = krome_x(krome_idx_NO2) * idom * 46d0
      O2HI(i,j,k) = krome_x(krome_idx_O2H) * idom * 33d0
      OCNI(i,j,k) = krome_x(krome_idx_OCN) * idom * 42d0
      CH3OH_DUSTI(i,j,k) = krome_x(krome_idx_CH3OH_DUST) * idom * 32d0
      HNCO_DUSTI(i,j,k) = krome_x(krome_idx_HNCO_DUST) * idom * 43d0
      H2CO_DUSTI(i,j,k) = krome_x(krome_idx_H2CO_DUST) * idom * 30d0
      SIH4_DUSTI(i,j,k) = krome_x(krome_idx_SIH4_DUST) * idom * 32d0
      H2SIO_DUSTI(i,j,k) = krome_x(krome_idx_H2SIO_DUST) * idom * 46d0
      SIC_DUSTI(i,j,k) = krome_x(krome_idx_SIC_DUST) * idom * 40d0
      SIC2_DUSTI(i,j,k) = krome_x(krome_idx_SIC2_DUST) * idom * 52d0
      SIC3_DUSTI(i,j,k) = krome_x(krome_idx_SIC3_DUST) * idom * 64d0
      CH4_DUSTI(i,j,k) = krome_x(krome_idx_CH4_DUST) * idom * 16d0
      CO_DUSTI(i,j,k) = krome_x(krome_idx_CO_DUST) * idom * 28d0
      H2O_DUSTI(i,j,k) = krome_x(krome_idx_H2O_DUST) * idom * 18d0
      NO_DUSTI(i,j,k) = krome_x(krome_idx_NO_DUST) * idom * 30d0
      CO2_DUSTI(i,j,k) = krome_x(krome_idx_CO2_DUST) * idom * 44d0
      N2_DUSTI(i,j,k) = krome_x(krome_idx_N2_DUST) * idom * 28d0
      HCN_DUSTI(i,j,k) = krome_x(krome_idx_HCN_DUST) * idom * 27d0
      NH3_DUSTI(i,j,k) = krome_x(krome_idx_NH3_DUST) * idom * 17d0
      O2_DUSTI(i,j,k) = krome_x(krome_idx_O2_DUST) * idom * 32d0
      NO2_DUSTI(i,j,k) = krome_x(krome_idx_NO2_DUST) * idom * 46d0
      HNO_DUSTI(i,j,k) = krome_x(krome_idx_HNO_DUST) * idom * 31d0
      O2H_DUSTI(i,j,k) = krome_x(krome_idx_O2H_DUST) * idom * 33d0
      H2CN_DUSTI(i,j,k) = krome_x(krome_idx_H2CN_DUST) * idom * 28d0
      MG_DUSTI(i,j,k) = krome_x(krome_idx_MG_DUST) * idom * 24d0
      HNC_DUSTI(i,j,k) = krome_x(krome_idx_HNC_DUST) * idom * 27d0
      E_DUSTI(i,j,k) = krome_x(krome_idx_E_DUST) * idom
      SIO_DUSTI(i,j,k) = krome_x(krome_idx_SIO_DUST) * idom * 44d0
      HCOII(i,j,k) = krome_x(krome_idx_HCOj) * idom * 29d0
      HII(i,j,k) = krome_x(krome_idx_Hj) * idom
      HOCII(i,j,k) = krome_x(krome_idx_HOCj) * idom * 29d0
      CII(i,j,k) = krome_x(krome_idx_Cj) * idom * 12d0
      CH2II(i,j,k) = krome_x(krome_idx_CH2j) * idom * 14d0
      CHII(i,j,k) = krome_x(krome_idx_CHj) * idom * 13d0
      H2COII(i,j,k) = krome_x(krome_idx_H2COj) * idom * 30d0
      MGII(i,j,k) = krome_x(krome_idx_MGj) * idom * 24d0
      NH3II(i,j,k) = krome_x(krome_idx_NH3j) * idom * 17d0
      NOII(i,j,k) = krome_x(krome_idx_NOj) * idom * 30d0
      SIII(i,j,k) = krome_x(krome_idx_SIj) * idom * 28d0
      SIC2II(i,j,k) = krome_x(krome_idx_SIC2j) * idom * 52d0
      SIC3II(i,j,k) = krome_x(krome_idx_SIC3j) * idom * 64d0
      SICII(i,j,k) = krome_x(krome_idx_SICj) * idom * 40d0
      SIH2II(i,j,k) = krome_x(krome_idx_SIH2j) * idom * 30d0
      SIH3II(i,j,k) = krome_x(krome_idx_SIH3j) * idom * 31d0
      CNII(i,j,k) = krome_x(krome_idx_CNj) * idom * 26d0
      COII(i,j,k) = krome_x(krome_idx_COj) * idom * 28d0
      N2II(i,j,k) = krome_x(krome_idx_N2j) * idom * 28d0
      O2II(i,j,k) = krome_x(krome_idx_O2j) * idom * 32d0
      H2OII(i,j,k) = krome_x(krome_idx_H2Oj) * idom * 18d0
      NH2II(i,j,k) = krome_x(krome_idx_NH2j) * idom * 16d0
      OII(i,j,k) = krome_x(krome_idx_Oj) * idom * 16d0
      OHII(i,j,k) = krome_x(krome_idx_OHj) * idom * 17d0
      CH3II(i,j,k) = krome_x(krome_idx_CH3j) * idom * 15d0
      CH4II(i,j,k) = krome_x(krome_idx_CH4j) * idom * 16d0
      NII(i,j,k) = krome_x(krome_idx_Nj) * idom * 14d0
      HCNII(i,j,k) = krome_x(krome_idx_HCNj) * idom * 27d0
      NHII(i,j,k) = krome_x(krome_idx_NHj) * idom * 15d0
      SIH4II(i,j,k) = krome_x(krome_idx_SIH4j) * idom * 32d0
      SIHII(i,j,k) = krome_x(krome_idx_SIHj) * idom * 29d0
      SIOII(i,j,k) = krome_x(krome_idx_SIOj) * idom * 44d0
      H2II(i,j,k) = krome_x(krome_idx_H2j) * idom * 2d0
      HeII(i,j,k) = krome_x(krome_idx_HEj) * idom * 4d0
      HNOII(i,j,k) = krome_x(krome_idx_HNOj) * idom * 31d0
      H2NOII(i,j,k) = krome_x(krome_idx_H2NOj) * idom * 32d0
      H3II(i,j,k) = krome_x(krome_idx_H3j) * idom * 3d0
      H3COII(i,j,k) = krome_x(krome_idx_H3COj) * idom * 31d0
      H3OII(i,j,k) = krome_x(krome_idx_H3Oj) * idom * 19d0
      HCNHII(i,j,k) = krome_x(krome_idx_HCNHj) * idom * 28d0
      HCO2II(i,j,k) = krome_x(krome_idx_HCO2j) * idom * 45d0
      HeHII(i,j,k) = krome_x(krome_idx_HEHj) * idom * 5d0
      N2HII(i,j,k) = krome_x(krome_idx_N2Hj) * idom * 29d0
      O2HII(i,j,k) = krome_x(krome_idx_O2Hj) * idom * 33d0
      SIH5II(i,j,k) = krome_x(krome_idx_SIH5j) * idom * 33d0
      SIOHII(i,j,k) = krome_x(krome_idx_SIOHj) * idom * 45d0
    end do
  end do
end do

end subroutine krome_renormref
