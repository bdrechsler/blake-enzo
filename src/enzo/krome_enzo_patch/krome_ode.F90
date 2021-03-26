
!############### MODULE ##############
module krome_ode
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2021-03-24 02:55:08
  !  Changeset xxxxxxx
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    use krome_fit
    use krome_user_commons
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
    real*8::rr
    integer::i,r1,r2,p1,p2,p3,p4

    n(:) = nin(:)

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),2.73d0)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

    k(:) = coe_tab(n(:)) !compute coefficients

    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      p1 = arr_p1(i)
      p2 = arr_p2(i)
      p3 = arr_p3(i)
      p4 = arr_p4(i)
      rr = k(i)*n(r1)*n(r2)
      dn(r1) = dn(r1) - rr
      dn(r2) = dn(r2) - rr
      dn(p1) = dn(p1) + rr
      dn(p2) = dn(p2) + rr
      dn(p3) = dn(p3) + rr
      dn(p4) = dn(p4) + rr
    end do

    dn(idx_H) = dn(idx_H) - 2.0*( 1.d-17*sqrt(Tgas)*n(idx_H)*get_Hnuclei(n(:)) - h2d(n(:))*n(idx_H2) )
    dn(idx_H2) = dn(idx_H2) + 1.d-17*sqrt(Tgas)*n(idx_H)*get_Hnuclei(n(:)) - h2d(n(:))*n(idx_H2)

    last_coe(:) = k(:)

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_gadiab
    use krome_getphys
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

    k(:) = last_coe(:) !get rate coefficients

    if(j==1) then
    elseif(j==1) then
      pdj(1) =  &
          -k(327)*n(idx_HCNj)  &
          -k(1218)*n(idx_H2COj)  &
          -k(328)*n(idx_HCNHj)  &
          -k(304)*n(idx_CNj)  &
          -k(357)*n(idx_SIH3j)  &
          -k(323)*n(idx_H3Oj)  &
          -k(345)*n(idx_NH3j)  &
          -k(295)*n(idx_CHj)  &
          -k(356)*n(idx_SIH2j)  &
          -k(313)*n(idx_H2Oj)  &
          -k(347)*n(idx_O2j)  &
          -k(1220)*n(idx_MGj)  &
          -k(1216)*n(idx_CH3j)  &
          -k(311)*n(idx_H2NOj)  &
          -k(1223)*n(idx_SIj)  &
          -k(300)*n(idx_CH3j)  &
          -k(335)*n(idx_HNOj)  &
          -k(341)*n(idx_NHj)  &
          -k(358)*n(idx_SIH3j)  &
          -k(352)*n(idx_SIC3j)  &
          -k(326)*n(idx_H3Oj)  &
          -k(354)*n(idx_SIH2j)  &
          -k(305)*n(idx_COj)  &
          -k(325)*n(idx_H3Oj)  &
          -k(339)*n(idx_N2Hj)  &
          -k(309)*n(idx_H2COj)  &
          -k(343)*n(idx_NH2j)  &
          -k(310)*n(idx_H2COj)  &
          -k(359)*n(idx_SIH4j)  &
          -k(355)*n(idx_SIH2j)  &
          -k(1307)  &
          -k(296)*n(idx_CH2j)  &
          -k(349)*n(idx_OHj)  &
          -k(308)*n(idx_H2COj)  &
          -k(9)*n(idx_H2)  &
          -k(1219)*n(idx_HEj)  &
          +k(9)*n(idx_H2)  &
          -k(329)*n(idx_HCNHj)  &
          -k(1221)*n(idx_Nj)  &
          -k(307)*n(idx_H2COj)  &
          -k(333)*n(idx_HCO2j)  &
          -k(1215)*n(idx_Cj)  &
          -k(337)*n(idx_HEHj)  &
          -k(346)*n(idx_NOj)  &
          -k(353)*n(idx_SIHj)  &
          -k(314)*n(idx_H2Oj)  &
          -k(351)*n(idx_SIC2j)  &
          -k(330)*n(idx_HCNHj)  &
          -k(312)*n(idx_H2NOj)  &
          -k(297)*n(idx_CH2j)  &
          -k(320)*n(idx_H3COj)  &
          -k(321)*n(idx_H3COj)  &
          -k(342)*n(idx_NH2j)  &
          -k(298)*n(idx_CH2j)  &
          -k(361)*n(idx_SIH5j)  &
          -k(322)*n(idx_H3COj)  &
          -k(301)*n(idx_CH3j)  &
          -k(348)*n(idx_O2Hj)  &
          -k(350)*n(idx_SICj)  &
          -k(344)*n(idx_NH3j)  &
          -k(316)*n(idx_H3j)  &
          -k(319)*n(idx_H3COj)  &
          -k(318)*n(idx_H3COj)  &
          -k(334)*n(idx_HCO2j)  &
          -k(1217)*n(idx_Hj)  &
          -k(365)*n(idx_SIOHj)  &
          -k(331)*n(idx_HCOj)  &
          -k(302)*n(idx_CH4j)  &
          -k(363)*n(idx_SIOj)  &
          -k(324)*n(idx_H3Oj)  &
          -k(360)*n(idx_SIH4j)  &
          -k(317)*n(idx_H3j)  &
          -k(1222)*n(idx_Oj)  &
          -k(340)*n(idx_N2Hj)  &
          -k(338)*n(idx_N2j)  &
          -k(303)*n(idx_CH4j)  &
          -k(332)*n(idx_HCO2j)  &
          -k(306)*n(idx_H2j)  &
          -k(362)*n(idx_SIH5j)  &
          -k(336)*n(idx_HOCj)  &
          -k(364)*n(idx_SIOHj)  &
          -k(299)*n(idx_CH3j)  &
          -k(315)*n(idx_H2Oj)
      pdj(2) =  &
          +k(298)*n(idx_CH2j)  &
          +k(300)*n(idx_CH3j)  &
          +k(319)*n(idx_H3COj)  &
          +k(301)*n(idx_CH3j)
      pdj(3) =  &
          +k(313)*n(idx_H2Oj)  &
          +k(314)*n(idx_H2Oj)  &
          +k(349)*n(idx_OHj)  &
          +k(346)*n(idx_NOj)  &
          +k(307)*n(idx_H2COj)  &
          +2.d0*k(347)*n(idx_O2j)  &
          +k(305)*n(idx_COj)  &
          +k(333)*n(idx_HCO2j)  &
          +k(1222)*n(idx_Oj)  &
          +k(363)*n(idx_SIOj)  &
          +k(324)*n(idx_H3Oj)
      pdj(4) =  &
          +k(330)*n(idx_HCNHj)
      pdj(5) =  &
          +k(329)*n(idx_HCNHj)
      pdj(6) =  &
          +k(313)*n(idx_H2Oj)  &
          +k(325)*n(idx_H3Oj)  &
          +k(300)*n(idx_CH3j)  &
          +k(354)*n(idx_SIH2j)  &
          +k(361)*n(idx_SIH5j)  &
          +k(316)*n(idx_H3j)  &
          -k(9)*n(idx_H2)  &
          +k(312)*n(idx_H2NOj)  &
          +k(358)*n(idx_SIH3j)  &
          +k(320)*n(idx_H3COj)  &
          +k(359)*n(idx_SIH4j)  &
          +k(296)*n(idx_CH2j)  &
          +k(308)*n(idx_H2COj)  &
          +k(324)*n(idx_H3Oj)
      pdj(7) =  &
          +k(296)*n(idx_CH2j)  &
          +k(297)*n(idx_CH2j)  &
          +k(352)*n(idx_SIC3j)  &
          +k(304)*n(idx_CNj)  &
          +k(305)*n(idx_COj)  &
          +k(1215)*n(idx_Cj)  &
          +k(295)*n(idx_CHj)  &
          +k(350)*n(idx_SICj)  &
          +k(351)*n(idx_SIC2j)
      pdj(8) =  &
          +2.d0*k(328)*n(idx_HCNHj)  &
          +k(299)*n(idx_CH3j)  &
          +2.d0*k(301)*n(idx_CH3j)  &
          +k(357)*n(idx_SIH3j)  &
          +k(329)*n(idx_HCNHj)  &
          +k(331)*n(idx_HCOj)  &
          +k(333)*n(idx_HCO2j)  &
          +2.d0*k(326)*n(idx_H3Oj)  &
          +k(353)*n(idx_SIHj)  &
          +2.d0*k(302)*n(idx_CH4j)  &
          +2.d0*k(342)*n(idx_NH2j)  &
          +k(344)*n(idx_NH3j)  &
          +2.d0*k(297)*n(idx_CH2j)  &
          +k(349)*n(idx_OHj)  &
          +k(356)*n(idx_SIH2j)  &
          +k(365)*n(idx_SIOHj)  &
          +k(360)*n(idx_SIH4j)  &
          +2.d0*k(345)*n(idx_NH3j)  &
          +k(330)*n(idx_HCNHj)  &
          +k(310)*n(idx_H2COj)  &
          +k(348)*n(idx_O2Hj)  &
          +2.d0*k(9)*n(idx_H2)  &
          +k(336)*n(idx_HOCj)  &
          +k(303)*n(idx_CH4j)  &
          +k(335)*n(idx_HNOj)  &
          +k(324)*n(idx_H3Oj)  &
          +k(343)*n(idx_NH2j)  &
          +k(298)*n(idx_CH2j)  &
          +k(339)*n(idx_N2Hj)  &
          +k(362)*n(idx_SIH5j)  &
          +k(316)*n(idx_H3j)  &
          +k(332)*n(idx_HCO2j)  &
          +k(321)*n(idx_H3COj)  &
          +3.d0*k(317)*n(idx_H3j)  &
          +2.d0*k(355)*n(idx_SIH2j)  &
          +2.d0*k(306)*n(idx_H2j)  &
          +k(341)*n(idx_NHj)  &
          +2.d0*k(314)*n(idx_H2Oj)  &
          +k(315)*n(idx_H2Oj)  &
          +2.d0*k(322)*n(idx_H3COj)  &
          +k(327)*n(idx_HCNj)  &
          +k(337)*n(idx_HEHj)  &
          +2.d0*k(309)*n(idx_H2COj)  &
          +k(295)*n(idx_CHj)  &
          +k(320)*n(idx_H3COj)  &
          +k(1217)*n(idx_Hj)  &
          +k(311)*n(idx_H2NOj)  &
          +k(323)*n(idx_H3Oj)
      pdj(9) =  &
          +k(319)*n(idx_H3COj)  &
          +k(323)*n(idx_H3Oj)
      pdj(10) =  &
          +k(334)*n(idx_HCO2j)  &
          +k(364)*n(idx_SIOHj)  &
          +k(315)*n(idx_H2Oj)  &
          +k(326)*n(idx_H3Oj)  &
          +k(325)*n(idx_H3Oj)  &
          +k(318)*n(idx_H3COj)
      pdj(11) =  &
          +k(348)*n(idx_O2Hj)
      pdj(12) =  &
          +k(318)*n(idx_H3COj)  &
          +k(302)*n(idx_CH4j)  &
          +k(307)*n(idx_H2COj)  &
          +k(299)*n(idx_CH3j)
      pdj(13) =  &
          +k(1218)*n(idx_H2COj)  &
          +k(321)*n(idx_H3COj)
      pdj(14) =  &
          +k(322)*n(idx_H3COj)  &
          +k(310)*n(idx_H2COj)
      pdj(15) =  &
          +k(1220)*n(idx_MGj)
      pdj(17) =  &
          +k(312)*n(idx_H2NOj)  &
          +k(335)*n(idx_HNOj)
      pdj(18) =  &
          +k(364)*n(idx_SIOHj)  &
          +k(353)*n(idx_SIHj)  &
          +k(354)*n(idx_SIH2j)  &
          +k(1223)*n(idx_SIj)  &
          +k(363)*n(idx_SIOj)  &
          +k(350)*n(idx_SICj)  &
          +k(355)*n(idx_SIH2j)
      pdj(19) =  &
          +k(352)*n(idx_SIC3j)
      pdj(21) =  &
          +k(351)*n(idx_SIC2j)
      pdj(22) =  &
          +k(357)*n(idx_SIH3j)  &
          +k(359)*n(idx_SIH4j)
      pdj(23) =  &
          +k(360)*n(idx_SIH4j)  &
          +k(361)*n(idx_SIH5j)
      pdj(24) =  &
          +k(327)*n(idx_HCNj)  &
          +k(328)*n(idx_HCNHj)
      pdj(25) =  &
          +k(334)*n(idx_HCO2j)  &
          +k(333)*n(idx_HCO2j)  &
          +k(309)*n(idx_H2COj)  &
          +k(320)*n(idx_H3COj)  &
          +k(331)*n(idx_HCOj)  &
          +k(336)*n(idx_HOCj)  &
          +k(308)*n(idx_H2COj)
      pdj(26) =  &
          +k(339)*n(idx_N2Hj)
      pdj(27) =  &
          +k(344)*n(idx_NH3j)
      pdj(28) =  &
          +k(1216)*n(idx_CH3j)  &
          +k(303)*n(idx_CH4j)
      pdj(30) =  &
          +2.d0*k(338)*n(idx_N2j)  &
          +k(341)*n(idx_NHj)  &
          +k(346)*n(idx_NOj)  &
          +k(1221)*n(idx_Nj)  &
          +k(304)*n(idx_CNj)  &
          +k(342)*n(idx_NH2j)  &
          +k(340)*n(idx_N2Hj)
      pdj(31) =  &
          +k(345)*n(idx_NH3j)  &
          +k(340)*n(idx_N2Hj)  &
          +k(343)*n(idx_NH2j)
      pdj(32) =  &
          +k(362)*n(idx_SIH5j)
      pdj(33) =  &
          +k(356)*n(idx_SIH2j)  &
          +k(358)*n(idx_SIH3j)
      pdj(34) =  &
          +k(365)*n(idx_SIOHj)
      pdj(35) =  &
          +k(337)*n(idx_HEHj)  &
          +k(1219)*n(idx_HEj)
      pdj(36) =  &
          +k(311)*n(idx_H2NOj)
      pdj(38) =  &
          +k(332)*n(idx_HCO2j)
      pdj(68) =  &
          +k(1307)
      pdj(70) =  &
          -k(331)*n(idx_HCOj)
      pdj(71) =  &
          -k(1217)*n(idx_Hj)
      pdj(72) =  &
          -k(336)*n(idx_HOCj)
      pdj(73) =  &
          -k(1215)*n(idx_Cj)
      pdj(74) =  &
          -k(296)*n(idx_CH2j)  &
          -k(297)*n(idx_CH2j)  &
          -k(298)*n(idx_CH2j)
      pdj(75) =  &
          -k(295)*n(idx_CHj)
      pdj(76) =  &
          -k(310)*n(idx_H2COj)  &
          -k(1218)*n(idx_H2COj)  &
          -k(307)*n(idx_H2COj)  &
          -k(309)*n(idx_H2COj)  &
          -k(308)*n(idx_H2COj)
      pdj(77) =  &
          -k(1220)*n(idx_MGj)
      pdj(78) =  &
          -k(344)*n(idx_NH3j)  &
          -k(345)*n(idx_NH3j)
      pdj(79) =  &
          -k(346)*n(idx_NOj)
      pdj(80) =  &
          -k(1223)*n(idx_SIj)
      pdj(81) =  &
          -k(351)*n(idx_SIC2j)
      pdj(82) =  &
          -k(352)*n(idx_SIC3j)
      pdj(83) =  &
          -k(350)*n(idx_SICj)
      pdj(84) =  &
          -k(354)*n(idx_SIH2j)  &
          -k(355)*n(idx_SIH2j)  &
          -k(356)*n(idx_SIH2j)
      pdj(85) =  &
          -k(357)*n(idx_SIH3j)  &
          -k(358)*n(idx_SIH3j)
      pdj(86) =  &
          -k(304)*n(idx_CNj)
      pdj(87) =  &
          -k(305)*n(idx_COj)
      pdj(88) =  &
          -k(338)*n(idx_N2j)
      pdj(89) =  &
          -k(347)*n(idx_O2j)
      pdj(90) =  &
          -k(314)*n(idx_H2Oj)  &
          -k(313)*n(idx_H2Oj)  &
          -k(315)*n(idx_H2Oj)
      pdj(91) =  &
          -k(343)*n(idx_NH2j)  &
          -k(342)*n(idx_NH2j)
      pdj(92) =  &
          -k(1222)*n(idx_Oj)
      pdj(93) =  &
          -k(349)*n(idx_OHj)
      pdj(94) =  &
          -k(300)*n(idx_CH3j)  &
          -k(301)*n(idx_CH3j)  &
          -k(299)*n(idx_CH3j)  &
          -k(1216)*n(idx_CH3j)
      pdj(95) =  &
          -k(303)*n(idx_CH4j)  &
          -k(302)*n(idx_CH4j)
      pdj(96) =  &
          -k(1221)*n(idx_Nj)
      pdj(97) =  &
          -k(327)*n(idx_HCNj)
      pdj(98) =  &
          -k(341)*n(idx_NHj)
      pdj(99) =  &
          -k(359)*n(idx_SIH4j)  &
          -k(360)*n(idx_SIH4j)
      pdj(100) =  &
          -k(353)*n(idx_SIHj)
      pdj(101) =  &
          -k(363)*n(idx_SIOj)
      pdj(102) =  &
          -k(306)*n(idx_H2j)
      pdj(103) =  &
          -k(1219)*n(idx_HEj)
      pdj(104) =  &
          -k(335)*n(idx_HNOj)
      pdj(105) =  &
          -k(312)*n(idx_H2NOj)  &
          -k(311)*n(idx_H2NOj)
      pdj(106) =  &
          -k(317)*n(idx_H3j)  &
          -k(316)*n(idx_H3j)
      pdj(107) =  &
          -k(320)*n(idx_H3COj)  &
          -k(322)*n(idx_H3COj)  &
          -k(319)*n(idx_H3COj)  &
          -k(318)*n(idx_H3COj)  &
          -k(321)*n(idx_H3COj)
      pdj(108) =  &
          -k(324)*n(idx_H3Oj)  &
          -k(323)*n(idx_H3Oj)  &
          -k(326)*n(idx_H3Oj)  &
          -k(325)*n(idx_H3Oj)
      pdj(109) =  &
          -k(329)*n(idx_HCNHj)  &
          -k(328)*n(idx_HCNHj)  &
          -k(330)*n(idx_HCNHj)
      pdj(110) =  &
          -k(332)*n(idx_HCO2j)  &
          -k(333)*n(idx_HCO2j)  &
          -k(334)*n(idx_HCO2j)
      pdj(111) =  &
          -k(337)*n(idx_HEHj)
      pdj(112) =  &
          -k(340)*n(idx_N2Hj)  &
          -k(339)*n(idx_N2Hj)
      pdj(113) =  &
          -k(348)*n(idx_O2Hj)
      pdj(114) =  &
          -k(362)*n(idx_SIH5j)  &
          -k(361)*n(idx_SIH5j)
      pdj(115) =  &
          -k(364)*n(idx_SIOHj)  &
          -k(365)*n(idx_SIOHj)
    elseif(j==2) then
      pdj(1) =  &
          +k(1)*n(idx_O)  &
          +k(1130)
      pdj(2) =  &
          -k(937)*n(idx_O2)  &
          -k(931)*n(idx_NO)  &
          -k(926)*n(idx_HCO)  &
          -k(477)*n(idx_SIj)  &
          -k(54)*n(idx_CNj)  &
          -k(1129)  &
          -k(467)*n(idx_HCOj)  &
          -k(79)*n(idx_Hj)  &
          -k(463)*n(idx_H3Oj)  &
          -k(103)*n(idx_H2j)  &
          -k(932)*n(idx_NO)  &
          -k(56)*n(idx_H2COj)  &
          -k(959)*n(idx_H2)  &
          -k(471)*n(idx_NHj)  &
          -k(925)*n(idx_H2CO)  &
          -k(465)*n(idx_HCNHj)  &
          -k(941)*n(idx_O)  &
          -k(469)*n(idx_Nj)  &
          -k(935)*n(idx_O2)  &
          -k(1130)  &
          -k(251)  &
          -k(464)*n(idx_HCNj)  &
          -k(60)*n(idx_NH2j)  &
          -k(942)*n(idx_OH)  &
          -k(461)*n(idx_H2Oj)  &
          -k(460)*n(idx_H2COj)  &
          -k(466)*n(idx_HCNHj)  &
          -k(142)*n(idx_HEj)  &
          -k(939)*n(idx_O2H)  &
          -k(1254)  &
          -k(58)*n(idx_Nj)  &
          -k(933)*n(idx_NO)  &
          -k(3)*n(idx_H2)  &
          -k(55)*n(idx_COj)  &
          -k(62)*n(idx_O2j)  &
          -k(473)*n(idx_Oj)  &
          -k(16)*n(idx_Cj)  &
          -k(513)*n(idx_H2j)  &
          -k(475)*n(idx_O2Hj)  &
          -k(663)*n(idx_HEj)  &
          -k(1202)*n(idx_H2)  &
          -k(59)*n(idx_N2j)  &
          -k(472)*n(idx_NH2j)  &
          -k(940)*n(idx_O)  &
          -k(459)*n(idx_COj)  &
          -k(930)*n(idx_N)  &
          -k(462)*n(idx_H3COj)  &
          -k(468)*n(idx_HNOj)  &
          -k(938)*n(idx_O2H)  &
          -k(936)*n(idx_O2)  &
          -k(479)*n(idx_SIOj)  &
          -k(927)*n(idx_HNO)  &
          -k(1)*n(idx_O)  &
          -k(57)*n(idx_H2Oj)  &
          -k(476)*n(idx_OHj)  &
          -k(478)*n(idx_SIHj)  &
          -k(929)*n(idx_N)  &
          -k(61)*n(idx_Oj)  &
          -k(581)*n(idx_H3j)  &
          -k(10)*n(idx_H)  &
          -k(971)*n(idx_H)  &
          -k(928)*n(idx_N2)  &
          -k(474)*n(idx_O2j)  &
          -k(924)*n(idx_CO2)  &
          -k(934)*n(idx_O2)  &
          -k(470)*n(idx_N2Hj)  &
          -k(63)*n(idx_OHj)
      pdj(3) =  &
          +k(937)*n(idx_O2)  &
          -k(941)*n(idx_O)  &
          -k(1)*n(idx_O)  &
          +k(476)*n(idx_OHj)  &
          +k(474)*n(idx_O2j)  &
          +k(931)*n(idx_NO)  &
          -k(940)*n(idx_O)  &
          +k(935)*n(idx_O2)  &
          +k(61)*n(idx_Oj)
      pdj(4) =  &
          +k(466)*n(idx_HCNHj)
      pdj(5) =  &
          +k(465)*n(idx_HCNHj)  &
          +k(931)*n(idx_NO)  &
          +k(928)*n(idx_N2)
      pdj(6) =  &
          -k(3)*n(idx_H2)  &
          +k(3)*n(idx_H2)  &
          +k(581)*n(idx_H3j)  &
          -k(1202)*n(idx_H2)  &
          +k(103)*n(idx_H2j)  &
          +k(971)*n(idx_H)  &
          -k(959)*n(idx_H2)
      pdj(7) =  &
          +k(10)*n(idx_H)  &
          +k(941)*n(idx_O)  &
          +k(930)*n(idx_N)  &
          +k(459)*n(idx_COj)  &
          +k(16)*n(idx_Cj)  &
          +k(3)*n(idx_H2)  &
          +k(971)*n(idx_H)  &
          +k(251)  &
          +k(1129)
      pdj(8) =  &
          +k(473)*n(idx_Oj)  &
          +k(929)*n(idx_N)  &
          +2.d0*k(10)*n(idx_H)  &
          +k(942)*n(idx_OH)  &
          +k(1129)  &
          +k(477)*n(idx_SIj)  &
          +k(3)*n(idx_H2)  &
          +k(935)*n(idx_O2)  &
          +k(933)*n(idx_NO)  &
          -k(10)*n(idx_H)  &
          +k(513)*n(idx_H2j)  &
          +k(251)  &
          -k(971)*n(idx_H)  &
          +k(79)*n(idx_Hj)  &
          +k(959)*n(idx_H2)  &
          +k(469)*n(idx_Nj)  &
          +k(940)*n(idx_O)  &
          +k(663)*n(idx_HEj)  &
          +k(934)*n(idx_O2)
      pdj(9) =  &
          +k(463)*n(idx_H3Oj)  &
          +k(57)*n(idx_H2Oj)
      pdj(10) =  &
          +k(461)*n(idx_H2Oj)  &
          +k(941)*n(idx_O)  &
          +k(936)*n(idx_O2)  &
          -k(942)*n(idx_OH)  &
          +k(63)*n(idx_OHj)  &
          +k(938)*n(idx_O2H)
      pdj(11) =  &
          -k(937)*n(idx_O2)  &
          -k(936)*n(idx_O2)  &
          +k(62)*n(idx_O2j)  &
          -k(935)*n(idx_O2)  &
          +k(939)*n(idx_O2H)  &
          +k(475)*n(idx_O2Hj)  &
          -k(934)*n(idx_O2)
      pdj(12) =  &
          +k(939)*n(idx_O2H)  &
          +k(959)*n(idx_H2)  &
          +k(926)*n(idx_HCO)  &
          +k(927)*n(idx_HNO)  &
          +k(925)*n(idx_H2CO)
      pdj(13) =  &
          +k(462)*n(idx_H3COj)  &
          +k(56)*n(idx_H2COj)  &
          -k(925)*n(idx_H2CO)
      pdj(14) =  &
          +k(937)*n(idx_O2)  &
          +k(942)*n(idx_OH)  &
          +k(925)*n(idx_H2CO)  &
          +k(460)*n(idx_H2COj)  &
          +k(924)*n(idx_CO2)  &
          -k(926)*n(idx_HCO)  &
          +k(938)*n(idx_O2H)  &
          +k(932)*n(idx_NO)
      pdj(17) =  &
          -k(932)*n(idx_NO)  &
          -k(931)*n(idx_NO)  &
          -k(933)*n(idx_NO)  &
          +k(468)*n(idx_HNOj)  &
          +k(927)*n(idx_HNO)
      pdj(18) =  &
          +k(478)*n(idx_SIHj)  &
          +k(479)*n(idx_SIOj)
      pdj(24) =  &
          +k(929)*n(idx_N)  &
          +k(54)*n(idx_CNj)  &
          +k(464)*n(idx_HCNj)
      pdj(25) =  &
          +k(55)*n(idx_COj)  &
          +k(924)*n(idx_CO2)  &
          +k(936)*n(idx_O2)  &
          +k(926)*n(idx_HCO)  &
          +k(940)*n(idx_O)  &
          +k(467)*n(idx_HCOj)  &
          +k(935)*n(idx_O2)
      pdj(26) =  &
          +k(470)*n(idx_N2Hj)  &
          -k(928)*n(idx_N2)  &
          +k(59)*n(idx_N2j)
      pdj(27) =  &
          +k(60)*n(idx_NH2j)
      pdj(28) =  &
          +k(1202)*n(idx_H2)
      pdj(30) =  &
          +k(58)*n(idx_Nj)  &
          +k(471)*n(idx_NHj)  &
          +k(928)*n(idx_N2)  &
          -k(930)*n(idx_N)  &
          -k(929)*n(idx_N)  &
          +k(932)*n(idx_NO)
      pdj(31) =  &
          +k(472)*n(idx_NH2j)  &
          +k(930)*n(idx_N)
      pdj(35) =  &
          +k(142)*n(idx_HEj)  &
          +k(663)*n(idx_HEj)
      pdj(36) =  &
          -k(927)*n(idx_HNO)
      pdj(38) =  &
          -k(924)*n(idx_CO2)  &
          +k(934)*n(idx_O2)
      pdj(43) =  &
          -k(938)*n(idx_O2H)  &
          -k(939)*n(idx_O2H)
      pdj(44) =  &
          +k(933)*n(idx_NO)
      pdj(53) =  &
          +k(1254)
      pdj(70) =  &
          +k(459)*n(idx_COj)  &
          +k(479)*n(idx_SIOj)  &
          +k(1)*n(idx_O)  &
          -k(467)*n(idx_HCOj)  &
          +k(474)*n(idx_O2j)
      pdj(71) =  &
          -k(79)*n(idx_Hj)
      pdj(73) =  &
          +k(663)*n(idx_HEj)  &
          -k(16)*n(idx_Cj)
      pdj(74) =  &
          +k(470)*n(idx_N2Hj)  &
          +k(472)*n(idx_NH2j)  &
          +k(468)*n(idx_HNOj)  &
          +k(471)*n(idx_NHj)  &
          +k(460)*n(idx_H2COj)  &
          +k(476)*n(idx_OHj)  &
          +k(462)*n(idx_H3COj)  &
          +k(464)*n(idx_HCNj)  &
          +k(581)*n(idx_H3j)  &
          +k(465)*n(idx_HCNHj)  &
          +k(475)*n(idx_O2Hj)  &
          +k(513)*n(idx_H2j)  &
          +k(463)*n(idx_H3Oj)  &
          +k(478)*n(idx_SIHj)  &
          +k(467)*n(idx_HCOj)  &
          +k(461)*n(idx_H2Oj)  &
          +k(466)*n(idx_HCNHj)
      pdj(75) =  &
          +k(58)*n(idx_Nj)  &
          +k(56)*n(idx_H2COj)  &
          +k(59)*n(idx_N2j)  &
          +k(57)*n(idx_H2Oj)  &
          +k(62)*n(idx_O2j)  &
          +k(55)*n(idx_COj)  &
          +k(1130)  &
          +k(103)*n(idx_H2j)  &
          +k(142)*n(idx_HEj)  &
          +k(54)*n(idx_CNj)  &
          +k(16)*n(idx_Cj)  &
          +k(79)*n(idx_Hj)  &
          +k(60)*n(idx_NH2j)  &
          +k(63)*n(idx_OHj)  &
          +k(61)*n(idx_Oj)
      pdj(76) =  &
          -k(460)*n(idx_H2COj)  &
          -k(56)*n(idx_H2COj)
      pdj(80) =  &
          -k(477)*n(idx_SIj)
      pdj(83) =  &
          +k(477)*n(idx_SIj)
      pdj(86) =  &
          -k(54)*n(idx_CNj)  &
          +k(469)*n(idx_Nj)
      pdj(87) =  &
          +k(473)*n(idx_Oj)  &
          -k(55)*n(idx_COj)  &
          -k(459)*n(idx_COj)
      pdj(88) =  &
          -k(59)*n(idx_N2j)
      pdj(89) =  &
          -k(474)*n(idx_O2j)  &
          -k(62)*n(idx_O2j)
      pdj(90) =  &
          -k(57)*n(idx_H2Oj)  &
          -k(461)*n(idx_H2Oj)
      pdj(91) =  &
          -k(472)*n(idx_NH2j)  &
          -k(60)*n(idx_NH2j)
      pdj(92) =  &
          -k(61)*n(idx_Oj)  &
          -k(473)*n(idx_Oj)
      pdj(93) =  &
          -k(476)*n(idx_OHj)  &
          -k(63)*n(idx_OHj)
      pdj(96) =  &
          -k(469)*n(idx_Nj)  &
          -k(58)*n(idx_Nj)
      pdj(97) =  &
          -k(464)*n(idx_HCNj)
      pdj(98) =  &
          -k(471)*n(idx_NHj)
      pdj(100) =  &
          -k(478)*n(idx_SIHj)
      pdj(101) =  &
          -k(479)*n(idx_SIOj)
      pdj(102) =  &
          -k(513)*n(idx_H2j)  &
          -k(103)*n(idx_H2j)
      pdj(103) =  &
          -k(663)*n(idx_HEj)  &
          -k(142)*n(idx_HEj)
      pdj(104) =  &
          -k(468)*n(idx_HNOj)
      pdj(106) =  &
          -k(581)*n(idx_H3j)
      pdj(107) =  &
          -k(462)*n(idx_H3COj)
      pdj(108) =  &
          -k(463)*n(idx_H3Oj)
      pdj(109) =  &
          -k(466)*n(idx_HCNHj)  &
          -k(465)*n(idx_HCNHj)
      pdj(112) =  &
          -k(470)*n(idx_N2Hj)
      pdj(113) =  &
          -k(475)*n(idx_O2Hj)
    elseif(j==3) then
      pdj(1) =  &
          +k(240)  &
          +k(1)*n(idx_CH)  &
          +k(283)
      pdj(2) =  &
          +k(898)*n(idx_CH2)  &
          -k(940)*n(idx_CH)  &
          -k(941)*n(idx_CH)  &
          -k(1)*n(idx_CH)
      pdj(3) =  &
          -k(1214)*n(idx_SI)  &
          -k(1197)*n(idx_C)  &
          -k(824)*n(idx_H2Oj)  &
          -k(1075)*n(idx_NH3)  &
          -k(834)*n(idx_SIH2j)  &
          -k(1057)*n(idx_CH4)  &
          -k(1071)*n(idx_HNO)  &
          -k(90)*n(idx_Hj)  &
          -k(1213)*n(idx_SIj)  &
          -k(444)*n(idx_CH3j)  &
          -k(826)*n(idx_N2j)  &
          -k(830)*n(idx_O2Hj)  &
          -k(1088)*n(idx_SIH3)  &
          -k(1080)*n(idx_OCN)  &
          -k(1063)*n(idx_H2O)  &
          -k(218)*n(idx_COj)  &
          -k(1087)*n(idx_SIH2)  &
          -k(1081)*n(idx_OH)  &
          -k(1065)*n(idx_HCN)  &
          -k(283)  &
          -k(1084)*n(idx_SIC)  &
          -k(1059)*n(idx_CN)  &
          -k(1060)*n(idx_CO2)  &
          -k(415)*n(idx_CHj)  &
          -k(825)*n(idx_HCO2j)  &
          -k(1086)*n(idx_SIH2)  &
          -k(966)*n(idx_H2)  &
          -k(1048)*n(idx_NH)  &
          -k(895)*n(idx_CH2)  &
          -k(445)*n(idx_CH3j)  &
          -k(836)*n(idx_SIOj)  &
          -k(1073)*n(idx_NH2)  &
          -k(1076)*n(idx_NO2)  &
          -k(600)*n(idx_H3j)  &
          -k(1078)*n(idx_O2H)  &
          -k(829)*n(idx_NH3j)  &
          -k(941)*n(idx_CH)  &
          -k(897)*n(idx_CH2)  &
          -k(1292)  &
          -k(1079)*n(idx_OCN)  &
          -k(832)*n(idx_SICj)  &
          -k(1089)*n(idx_SIH4)  &
          -k(1068)*n(idx_HCO)  &
          -k(940)*n(idx_CH)  &
          -k(1194)*n(idx_Cj)  &
          -k(1067)*n(idx_HCO)  &
          -k(217)*n(idx_CNj)  &
          -k(1208)*n(idx_H)  &
          -k(1058)*n(idx_CN)  &
          -k(828)*n(idx_NH2j)  &
          -k(768)*n(idx_NHj)  &
          -k(219)*n(idx_N2j)  &
          -k(1062)*n(idx_H2CO)  &
          -k(1072)*n(idx_N2)  &
          -k(823)*n(idx_CH4j)  &
          -k(833)*n(idx_SIHj)  &
          -k(896)*n(idx_CH2)  &
          -k(1069)*n(idx_HNO)  &
          -k(1070)*n(idx_HNO)  &
          -k(1074)*n(idx_NH2)  &
          -k(827)*n(idx_N2Hj)  &
          -k(1)*n(idx_CH)  &
          -k(916)*n(idx_CH3)  &
          -k(1082)*n(idx_SIC2)  &
          -k(835)*n(idx_SIH3j)  &
          -4.d0*k(1212)*n(idx_O)  &
          -k(527)*n(idx_H2j)  &
          -k(1061)*n(idx_H2CN)  &
          -k(1047)*n(idx_NH)  &
          -k(831)*n(idx_OHj)  &
          -k(1066)*n(idx_HCN)  &
          -k(898)*n(idx_CH2)  &
          -k(1083)*n(idx_SIC3)  &
          -k(1077)*n(idx_NO)  &
          -k(1064)*n(idx_HCN)  &
          -k(240)  &
          -k(1085)*n(idx_SIC)  &
          -k(422)*n(idx_CH2j)  &
          -k(1090)*n(idx_SIH)  &
          -k(599)*n(idx_H3j)  &
          -k(917)*n(idx_CH3)
      pdj(5) =  &
          -k(1065)*n(idx_HCN)  &
          -k(1064)*n(idx_HCN)  &
          -k(1066)*n(idx_HCN)
      pdj(6) =  &
          -k(966)*n(idx_H2)  &
          +k(916)*n(idx_CH3)  &
          +k(829)*n(idx_NH3j)  &
          +k(600)*n(idx_H3j)  &
          +k(1061)*n(idx_H2CN)  &
          +k(445)*n(idx_CH3j)  &
          +k(895)*n(idx_CH2)  &
          +k(824)*n(idx_H2Oj)  &
          +k(835)*n(idx_SIH3j)  &
          +k(1086)*n(idx_SIH2)
      pdj(7) =  &
          +k(941)*n(idx_CH)  &
          +k(1085)*n(idx_SIC)  &
          -k(1197)*n(idx_C)  &
          +k(832)*n(idx_SICj)  &
          +k(1059)*n(idx_CN)
      pdj(8) =  &
          +k(1047)*n(idx_NH)  &
          +k(1090)*n(idx_SIH)  &
          +k(1073)*n(idx_NH2)  &
          +k(831)*n(idx_OHj)  &
          +k(415)*n(idx_CHj)  &
          +k(917)*n(idx_CH3)  &
          +k(1069)*n(idx_HNO)  &
          +k(940)*n(idx_CH)  &
          +2.d0*k(896)*n(idx_CH2)  &
          +k(834)*n(idx_SIH2j)  &
          +k(1081)*n(idx_OH)  &
          +k(599)*n(idx_H3j)  &
          +k(444)*n(idx_CH3j)  &
          +k(897)*n(idx_CH2)  &
          +k(1067)*n(idx_HCO)  &
          +k(527)*n(idx_H2j)  &
          +k(828)*n(idx_NH2j)  &
          +k(422)*n(idx_CH2j)  &
          +k(1066)*n(idx_HCN)  &
          +k(90)*n(idx_Hj)  &
          +k(916)*n(idx_CH3)  &
          +k(966)*n(idx_H2)  &
          +2.d0*k(1087)*n(idx_SIH2)  &
          +k(833)*n(idx_SIHj)  &
          -k(1208)*n(idx_H)  &
          +k(1088)*n(idx_SIH3)
      pdj(9) =  &
          -k(1063)*n(idx_H2O)
      pdj(10) =  &
          +k(1078)*n(idx_O2H)  &
          +k(1062)*n(idx_H2CO)  &
          +2.d0*k(1063)*n(idx_H2O)  &
          +k(1070)*n(idx_HNO)  &
          +k(823)*n(idx_CH4j)  &
          -k(1081)*n(idx_OH)  &
          +k(966)*n(idx_H2)  &
          +k(1208)*n(idx_H)  &
          +k(1089)*n(idx_SIH4)  &
          +k(898)*n(idx_CH2)  &
          +k(1075)*n(idx_NH3)  &
          +k(1068)*n(idx_HCO)  &
          +k(1057)*n(idx_CH4)  &
          +k(1048)*n(idx_NH)  &
          +k(941)*n(idx_CH)  &
          +k(1074)*n(idx_NH2)  &
          +k(1064)*n(idx_HCN)
      pdj(11) =  &
          +k(1078)*n(idx_O2H)  &
          +k(1076)*n(idx_NO2)  &
          +k(825)*n(idx_HCO2j)  &
          +k(836)*n(idx_SIOj)  &
          +k(1080)*n(idx_OCN)  &
          +k(1060)*n(idx_CO2)  &
          +k(1081)*n(idx_OH)  &
          +2.d0*k(1212)*n(idx_O)  &
          +k(830)*n(idx_O2Hj)  &
          +k(1077)*n(idx_NO)  &
          +k(1071)*n(idx_HNO)
      pdj(12) =  &
          -k(895)*n(idx_CH2)  &
          -k(898)*n(idx_CH2)  &
          -k(896)*n(idx_CH2)  &
          -k(897)*n(idx_CH2)
      pdj(13) =  &
          -k(1062)*n(idx_H2CO)  &
          +k(917)*n(idx_CH3)
      pdj(14) =  &
          +k(897)*n(idx_CH2)  &
          +k(1062)*n(idx_H2CO)  &
          -k(1068)*n(idx_HCO)  &
          -k(1067)*n(idx_HCO)
      pdj(16) =  &
          -k(1075)*n(idx_NH3)
      pdj(17) =  &
          +k(1047)*n(idx_NH)  &
          +k(1059)*n(idx_CN)  &
          +k(1076)*n(idx_NO2)  &
          +k(1070)*n(idx_HNO)  &
          +k(1079)*n(idx_OCN)  &
          +k(1072)*n(idx_N2)  &
          -k(1077)*n(idx_NO)
      pdj(18) =  &
          -k(1214)*n(idx_SI)  &
          +k(1084)*n(idx_SIC)
      pdj(19) =  &
          +k(1083)*n(idx_SIC3)  &
          -k(1082)*n(idx_SIC2)
      pdj(20) =  &
          -k(1083)*n(idx_SIC3)
      pdj(21) =  &
          -k(1084)*n(idx_SIC)  &
          +k(1082)*n(idx_SIC2)  &
          -k(1085)*n(idx_SIC)
      pdj(22) =  &
          -k(1086)*n(idx_SIH2)  &
          -k(1087)*n(idx_SIH2)
      pdj(23) =  &
          -k(1088)*n(idx_SIH3)  &
          +k(1089)*n(idx_SIH4)
      pdj(24) =  &
          -k(1059)*n(idx_CN)  &
          +k(1080)*n(idx_OCN)  &
          +k(1064)*n(idx_HCN)  &
          +k(217)*n(idx_CNj)  &
          -k(1058)*n(idx_CN)
      pdj(25) =  &
          +k(1083)*n(idx_SIC3)  &
          +k(916)*n(idx_CH3)  &
          +k(1079)*n(idx_OCN)  &
          +k(1084)*n(idx_SIC)  &
          +k(1065)*n(idx_HCN)  &
          +k(940)*n(idx_CH)  &
          +k(896)*n(idx_CH2)  &
          +k(218)*n(idx_COj)  &
          +k(1068)*n(idx_HCO)  &
          +k(1058)*n(idx_CN)  &
          +k(1060)*n(idx_CO2)  &
          +k(895)*n(idx_CH2)  &
          +k(1082)*n(idx_SIC2)  &
          +k(1197)*n(idx_C)
      pdj(26) =  &
          +k(219)*n(idx_N2j)  &
          +k(827)*n(idx_N2Hj)  &
          -k(1072)*n(idx_N2)
      pdj(27) =  &
          -k(1073)*n(idx_NH2)  &
          -k(1074)*n(idx_NH2)  &
          +k(1075)*n(idx_NH3)
      pdj(28) =  &
          -k(916)*n(idx_CH3)  &
          -k(917)*n(idx_CH3)  &
          +k(1057)*n(idx_CH4)
      pdj(29) =  &
          -k(1057)*n(idx_CH4)
      pdj(30) =  &
          +k(768)*n(idx_NHj)  &
          +k(826)*n(idx_N2j)  &
          +k(1072)*n(idx_N2)  &
          +k(1058)*n(idx_CN)  &
          +k(1048)*n(idx_NH)  &
          +k(1077)*n(idx_NO)
      pdj(31) =  &
          +k(1074)*n(idx_NH2)  &
          -k(1047)*n(idx_NH)  &
          -k(1048)*n(idx_NH)  &
          +k(1065)*n(idx_HCN)  &
          +k(1071)*n(idx_HNO)
      pdj(32) =  &
          -k(1089)*n(idx_SIH4)
      pdj(33) =  &
          -k(1090)*n(idx_SIH)
      pdj(34) =  &
          +k(1085)*n(idx_SIC)  &
          +k(1090)*n(idx_SIH)  &
          +k(1087)*n(idx_SIH2)  &
          +k(1214)*n(idx_SI)  &
          +k(1086)*n(idx_SIH2)
      pdj(36) =  &
          -k(1071)*n(idx_HNO)  &
          +k(1073)*n(idx_NH2)  &
          -k(1070)*n(idx_HNO)  &
          -k(1069)*n(idx_HNO)
      pdj(38) =  &
          +k(1067)*n(idx_HCO)  &
          -k(1060)*n(idx_CO2)
      pdj(39) =  &
          -k(1061)*n(idx_H2CN)
      pdj(40) =  &
          +k(1088)*n(idx_SIH3)
      pdj(42) =  &
          -k(1076)*n(idx_NO2)  &
          +k(1069)*n(idx_HNO)
      pdj(43) =  &
          -k(1078)*n(idx_O2H)
      pdj(44) =  &
          -k(1080)*n(idx_OCN)  &
          -k(1079)*n(idx_OCN)  &
          +k(1061)*n(idx_H2CN)  &
          +k(1066)*n(idx_HCN)
      pdj(55) =  &
          +k(1292)
      pdj(70) =  &
          +k(445)*n(idx_CH3j)  &
          +k(422)*n(idx_CH2j)  &
          +k(1)*n(idx_CH)  &
          +k(825)*n(idx_HCO2j)
      pdj(71) =  &
          -k(90)*n(idx_Hj)
      pdj(73) =  &
          -k(1194)*n(idx_Cj)
      pdj(74) =  &
          -k(422)*n(idx_CH2j)
      pdj(75) =  &
          -k(415)*n(idx_CHj)
      pdj(76) =  &
          +k(444)*n(idx_CH3j)
      pdj(78) =  &
          -k(829)*n(idx_NH3j)
      pdj(79) =  &
          +k(826)*n(idx_N2j)
      pdj(80) =  &
          -k(1213)*n(idx_SIj)  &
          +k(836)*n(idx_SIOj)
      pdj(83) =  &
          -k(832)*n(idx_SICj)
      pdj(84) =  &
          -k(834)*n(idx_SIH2j)
      pdj(85) =  &
          -k(835)*n(idx_SIH3j)
      pdj(86) =  &
          -k(217)*n(idx_CNj)
      pdj(87) =  &
          +k(415)*n(idx_CHj)  &
          -k(218)*n(idx_COj)  &
          +k(1194)*n(idx_Cj)
      pdj(88) =  &
          -k(219)*n(idx_N2j)  &
          -k(826)*n(idx_N2j)
      pdj(89) =  &
          +k(831)*n(idx_OHj)  &
          +k(824)*n(idx_H2Oj)
      pdj(90) =  &
          -k(824)*n(idx_H2Oj)  &
          +k(599)*n(idx_H3j)
      pdj(91) =  &
          -k(828)*n(idx_NH2j)
      pdj(92) =  &
          +k(219)*n(idx_N2j)  &
          +k(90)*n(idx_Hj)  &
          +k(217)*n(idx_CNj)  &
          +k(283)  &
          +k(218)*n(idx_COj)  &
          +k(240)
      pdj(93) =  &
          +k(827)*n(idx_N2Hj)  &
          +k(600)*n(idx_H3j)  &
          -k(831)*n(idx_OHj)  &
          +k(527)*n(idx_H2j)  &
          +k(768)*n(idx_NHj)  &
          +k(830)*n(idx_O2Hj)
      pdj(94) =  &
          -k(445)*n(idx_CH3j)  &
          +k(823)*n(idx_CH4j)  &
          -k(444)*n(idx_CH3j)
      pdj(95) =  &
          -k(823)*n(idx_CH4j)
      pdj(98) =  &
          -k(768)*n(idx_NHj)
      pdj(100) =  &
          -k(833)*n(idx_SIHj)
      pdj(101) =  &
          -k(836)*n(idx_SIOj)  &
          +k(1213)*n(idx_SIj)  &
          +k(832)*n(idx_SICj)  &
          +k(833)*n(idx_SIHj)
      pdj(102) =  &
          -k(527)*n(idx_H2j)
      pdj(104) =  &
          +k(828)*n(idx_NH2j)  &
          +k(829)*n(idx_NH3j)
      pdj(106) =  &
          -k(600)*n(idx_H3j)  &
          -k(599)*n(idx_H3j)
      pdj(110) =  &
          -k(825)*n(idx_HCO2j)
      pdj(112) =  &
          -k(827)*n(idx_N2Hj)
      pdj(113) =  &
          -k(830)*n(idx_O2Hj)
      pdj(115) =  &
          +k(834)*n(idx_SIH2j)  &
          +k(835)*n(idx_SIH3j)
    elseif(j==4) then
      pdj(3) =  &
          +k(845)*n(idx_OHj)
      pdj(4) =  &
          -k(652)*n(idx_O2Hj)  &
          -k(649)*n(idx_HCOj)  &
          -k(1306)  &
          -k(684)*n(idx_HEj)  &
          -k(610)*n(idx_H3Oj)  &
          -k(651)*n(idx_N2Hj)  &
          -k(686)*n(idx_HEj)  &
          -k(980)*n(idx_H)  &
          -k(648)*n(idx_H3COj)  &
          -k(560)*n(idx_H2Oj)  &
          -k(263)  &
          -k(845)*n(idx_OHj)  &
          -k(590)*n(idx_H3j)  &
          -k(1152)  &
          -k(776)*n(idx_NH2j)  &
          -k(2)*n(idx_Hj)  &
          -k(647)*n(idx_H2COj)  &
          -k(685)*n(idx_HEj)  &
          -k(650)*n(idx_HNOj)  &
          -k(761)*n(idx_NHj)  &
          -k(627)*n(idx_HCNj)  &
          -k(408)*n(idx_CHj)
      pdj(5) =  &
          +k(2)*n(idx_Hj)  &
          +k(980)*n(idx_H)
      pdj(6) =  &
          +k(590)*n(idx_H3j)
      pdj(7) =  &
          +k(408)*n(idx_CHj)  &
          +k(686)*n(idx_HEj)
      pdj(8) =  &
          +k(263)  &
          +k(980)*n(idx_H)  &
          +k(685)*n(idx_HEj)  &
          -k(980)*n(idx_H)  &
          +k(1152)  &
          +k(684)*n(idx_HEj)
      pdj(9) =  &
          +k(610)*n(idx_H3Oj)
      pdj(10) =  &
          +k(560)*n(idx_H2Oj)
      pdj(11) =  &
          +k(652)*n(idx_O2Hj)
      pdj(13) =  &
          +k(648)*n(idx_H3COj)
      pdj(14) =  &
          +k(647)*n(idx_H2COj)
      pdj(17) =  &
          +k(650)*n(idx_HNOj)
      pdj(24) =  &
          +k(627)*n(idx_HCNj)  &
          +k(263)  &
          +k(1152)
      pdj(25) =  &
          +k(649)*n(idx_HCOj)
      pdj(26) =  &
          +k(651)*n(idx_N2Hj)
      pdj(30) =  &
          +k(761)*n(idx_NHj)  &
          +k(685)*n(idx_HEj)
      pdj(31) =  &
          +k(776)*n(idx_NH2j)
      pdj(35) =  &
          +k(686)*n(idx_HEj)  &
          +k(685)*n(idx_HEj)  &
          +k(684)*n(idx_HEj)
      pdj(67) =  &
          +k(1306)
      pdj(70) =  &
          -k(649)*n(idx_HCOj)
      pdj(71) =  &
          +k(2)*n(idx_Hj)  &
          -k(2)*n(idx_Hj)
      pdj(73) =  &
          +k(685)*n(idx_HEj)
      pdj(75) =  &
          -k(408)*n(idx_CHj)
      pdj(76) =  &
          -k(647)*n(idx_H2COj)
      pdj(86) =  &
          +k(684)*n(idx_HEj)
      pdj(90) =  &
          -k(560)*n(idx_H2Oj)
      pdj(91) =  &
          -k(776)*n(idx_NH2j)
      pdj(93) =  &
          -k(845)*n(idx_OHj)
      pdj(97) =  &
          -k(627)*n(idx_HCNj)
      pdj(98) =  &
          -k(761)*n(idx_NHj)  &
          +k(686)*n(idx_HEj)
      pdj(103) =  &
          -k(685)*n(idx_HEj)  &
          -k(684)*n(idx_HEj)  &
          -k(686)*n(idx_HEj)
      pdj(104) =  &
          -k(650)*n(idx_HNOj)
      pdj(106) =  &
          -k(590)*n(idx_H3j)
      pdj(107) =  &
          -k(648)*n(idx_H3COj)
      pdj(108) =  &
          -k(610)*n(idx_H3Oj)
      pdj(109) =  &
          +k(648)*n(idx_H3COj)  &
          +k(647)*n(idx_H2COj)  &
          +k(649)*n(idx_HCOj)  &
          +k(627)*n(idx_HCNj)  &
          +k(776)*n(idx_NH2j)  &
          +k(650)*n(idx_HNOj)  &
          +k(761)*n(idx_NHj)  &
          +k(651)*n(idx_N2Hj)  &
          +k(845)*n(idx_OHj)  &
          +k(408)*n(idx_CHj)  &
          +k(652)*n(idx_O2Hj)  &
          +k(590)*n(idx_H3j)  &
          +k(610)*n(idx_H3Oj)  &
          +k(560)*n(idx_H2Oj)
      pdj(112) =  &
          -k(651)*n(idx_N2Hj)
      pdj(113) =  &
          -k(652)*n(idx_O2Hj)
    elseif(j==5) then
      pdj(2) =  &
          +k(816)*n(idx_Oj)  &
          +k(678)*n(idx_HEj)
      pdj(3) =  &
          -k(1066)*n(idx_O)  &
          -k(1065)*n(idx_O)  &
          -k(1064)*n(idx_O)  &
          +k(842)*n(idx_OHj)
      pdj(5) =  &
          -k(628)*n(idx_H2COj)  &
          -k(1065)*n(idx_O)  &
          -k(816)*n(idx_Oj)  &
          -k(162)*n(idx_Nj)  &
          -k(609)*n(idx_H3Oj)  &
          -k(633)*n(idx_O2Hj)  &
          -k(1267)  &
          -k(815)*n(idx_Oj)  &
          -k(1096)*n(idx_OH)  &
          -k(977)*n(idx_H)  &
          -k(135)*n(idx_COj)  &
          -k(66)*n(idx_CNj)  &
          -k(629)*n(idx_H3COj)  &
          -k(680)*n(idx_HEj)  &
          -k(774)*n(idx_NH2j)  &
          -k(260)  &
          -k(624)*n(idx_HCNj)  &
          -k(136)*n(idx_N2j)  &
          -k(678)*n(idx_HEj)  &
          -k(1095)*n(idx_OH)  &
          -k(677)*n(idx_HEj)  &
          -k(679)*n(idx_HEj)  &
          -k(630)*n(idx_HCOj)  &
          -k(82)*n(idx_Hj)  &
          -k(108)*n(idx_H2j)  &
          -k(406)*n(idx_CHj)  &
          -k(588)*n(idx_H3j)  &
          -k(759)*n(idx_NHj)  &
          -k(631)*n(idx_HNOj)  &
          -k(1064)*n(idx_O)  &
          -k(1148)  &
          -k(1066)*n(idx_O)  &
          -k(557)*n(idx_H2Oj)  &
          -k(632)*n(idx_N2Hj)  &
          -k(842)*n(idx_OHj)
      pdj(6) =  &
          +k(977)*n(idx_H)  &
          +k(588)*n(idx_H3j)  &
          +k(108)*n(idx_H2j)
      pdj(7) =  &
          +k(406)*n(idx_CHj)
      pdj(8) =  &
          +k(82)*n(idx_Hj)  &
          +k(1066)*n(idx_O)  &
          -k(977)*n(idx_H)  &
          +k(260)  &
          +k(1148)  &
          +k(677)*n(idx_HEj)  &
          +k(679)*n(idx_HEj)
      pdj(9) =  &
          +k(1095)*n(idx_OH)  &
          +k(609)*n(idx_H3Oj)
      pdj(10) =  &
          -k(1095)*n(idx_OH)  &
          -k(1096)*n(idx_OH)  &
          +k(557)*n(idx_H2Oj)  &
          +k(1064)*n(idx_O)
      pdj(11) =  &
          +k(633)*n(idx_O2Hj)
      pdj(13) =  &
          +k(629)*n(idx_H3COj)
      pdj(14) =  &
          +k(628)*n(idx_H2COj)
      pdj(17) =  &
          +k(631)*n(idx_HNOj)
      pdj(24) =  &
          +k(624)*n(idx_HCNj)  &
          +k(260)  &
          +k(66)*n(idx_CNj)  &
          +k(977)*n(idx_H)  &
          +k(1148)  &
          +k(1095)*n(idx_OH)  &
          +k(1064)*n(idx_O)
      pdj(25) =  &
          +k(1065)*n(idx_O)  &
          +k(1096)*n(idx_OH)  &
          +k(135)*n(idx_COj)  &
          +k(630)*n(idx_HCOj)
      pdj(26) =  &
          +k(632)*n(idx_N2Hj)  &
          +k(136)*n(idx_N2j)
      pdj(27) =  &
          +k(1096)*n(idx_OH)
      pdj(30) =  &
          +k(815)*n(idx_Oj)  &
          +k(759)*n(idx_NHj)  &
          +k(680)*n(idx_HEj)  &
          +k(162)*n(idx_Nj)  &
          +k(679)*n(idx_HEj)
      pdj(31) =  &
          +k(1065)*n(idx_O)  &
          +k(774)*n(idx_NH2j)
      pdj(35) =  &
          +k(677)*n(idx_HEj)  &
          +k(680)*n(idx_HEj)  &
          +k(678)*n(idx_HEj)  &
          +k(679)*n(idx_HEj)
      pdj(44) =  &
          +k(1066)*n(idx_O)
      pdj(59) =  &
          +k(1267)
      pdj(70) =  &
          +k(815)*n(idx_Oj)  &
          -k(630)*n(idx_HCOj)
      pdj(71) =  &
          -k(82)*n(idx_Hj)
      pdj(73) =  &
          +k(679)*n(idx_HEj)
      pdj(75) =  &
          +k(680)*n(idx_HEj)  &
          -k(406)*n(idx_CHj)
      pdj(76) =  &
          -k(628)*n(idx_H2COj)
      pdj(79) =  &
          +k(816)*n(idx_Oj)
      pdj(86) =  &
          +k(677)*n(idx_HEj)  &
          -k(66)*n(idx_CNj)
      pdj(87) =  &
          -k(135)*n(idx_COj)
      pdj(88) =  &
          -k(136)*n(idx_N2j)
      pdj(90) =  &
          -k(557)*n(idx_H2Oj)
      pdj(91) =  &
          -k(774)*n(idx_NH2j)
      pdj(92) =  &
          -k(816)*n(idx_Oj)  &
          -k(815)*n(idx_Oj)
      pdj(93) =  &
          -k(842)*n(idx_OHj)
      pdj(96) =  &
          +k(678)*n(idx_HEj)  &
          -k(162)*n(idx_Nj)
      pdj(97) =  &
          +k(82)*n(idx_Hj)  &
          +k(135)*n(idx_COj)  &
          +k(162)*n(idx_Nj)  &
          +k(66)*n(idx_CNj)  &
          +k(108)*n(idx_H2j)  &
          +k(136)*n(idx_N2j)  &
          -k(624)*n(idx_HCNj)
      pdj(98) =  &
          -k(759)*n(idx_NHj)
      pdj(102) =  &
          -k(108)*n(idx_H2j)
      pdj(103) =  &
          -k(677)*n(idx_HEj)  &
          -k(678)*n(idx_HEj)  &
          -k(679)*n(idx_HEj)  &
          -k(680)*n(idx_HEj)
      pdj(104) =  &
          -k(631)*n(idx_HNOj)
      pdj(106) =  &
          -k(588)*n(idx_H3j)
      pdj(107) =  &
          -k(629)*n(idx_H3COj)
      pdj(108) =  &
          -k(609)*n(idx_H3Oj)
      pdj(109) =  &
          +k(624)*n(idx_HCNj)  &
          +k(557)*n(idx_H2Oj)  &
          +k(588)*n(idx_H3j)  &
          +k(631)*n(idx_HNOj)  &
          +k(628)*n(idx_H2COj)  &
          +k(630)*n(idx_HCOj)  &
          +k(632)*n(idx_N2Hj)  &
          +k(629)*n(idx_H3COj)  &
          +k(842)*n(idx_OHj)  &
          +k(774)*n(idx_NH2j)  &
          +k(759)*n(idx_NHj)  &
          +k(609)*n(idx_H3Oj)  &
          +k(633)*n(idx_O2Hj)  &
          +k(406)*n(idx_CHj)
      pdj(112) =  &
          -k(632)*n(idx_N2Hj)
      pdj(113) =  &
          -k(633)*n(idx_O2Hj)
    elseif(j==6) then
      pdj(1) =  &
          +k(9)*n(idx_E)  &
          +k(234)  &
          -k(9)*n(idx_E)  &
          +k(235)
      pdj(2) =  &
          -k(3)*n(idx_CH)  &
          -k(1202)*n(idx_CH)  &
          +k(956)*n(idx_C)  &
          -k(959)*n(idx_CH)
      pdj(3) =  &
          +2.d0*k(7)*n(idx_O2)  &
          -k(966)*n(idx_O)  &
          +k(8)*n(idx_OH)
      pdj(5) =  &
          +k(960)*n(idx_CN)
      pdj(6) =  &
          -k(1204)*n(idx_SIHj)  &
          -k(9)*n(idx_E)  &
          -k(956)*n(idx_C)  &
          -k(960)*n(idx_CN)  &
          -k(8)*n(idx_OH)  &
          -k(529)*n(idx_Cj)  &
          -k(963)*n(idx_NH)  &
          +k(7)*n(idx_O2)  &
          +k(5)*n(idx_H2O)  &
          -k(964)*n(idx_O2)  &
          -k(537)*n(idx_HEj)  &
          -k(540)*n(idx_N2j)  &
          -k(536)*n(idx_HCNj)  &
          -k(967)*n(idx_OH)  &
          -k(236)  &
          -k(958)*n(idx_CH3)  &
          -k(7)*n(idx_O2)  &
          -k(1202)*n(idx_CH)  &
          -k(539)*n(idx_Nj)  &
          -k(531)*n(idx_CH2j)  &
          -k(533)*n(idx_COj)  &
          -k(546)*n(idx_OHj)  &
          -k(535)*n(idx_H2Oj)  &
          -k(962)*n(idx_NH2)  &
          -k(1200)*n(idx_Cj)  &
          -k(1205)*n(idx_SIH3j)  &
          -k(544)*n(idx_Oj)  &
          -k(532)*n(idx_CNj)  &
          -k(234)  &
          -k(961)*n(idx_N)  &
          -k(543)*n(idx_NH2j)  &
          -k(548)*n(idx_SIOj)  &
          -k(966)*n(idx_O)  &
          -k(534)*n(idx_COj)  &
          -k(1203)*n(idx_SIj)  &
          -k(965)*n(idx_O2)  &
          -k(959)*n(idx_CH)  &
          -k(5)*n(idx_H2O)  &
          -k(547)*n(idx_SIH4j)  &
          +2.d0*k(4)*n(idx_H2)  &
          +k(3)*n(idx_CH)  &
          -k(3)*n(idx_CH)  &
          -k(6)*n(idx_HOCj)  &
          -k(517)*n(idx_H2j)  &
          -k(1201)*n(idx_C)  &
          -4.d0*k(4)*n(idx_H2)  &
          -k(545)*n(idx_O2Hj)  &
          -k(542)*n(idx_NHj)  &
          -k(541)*n(idx_NHj)  &
          -k(957)*n(idx_CH2)  &
          -k(530)*n(idx_CHj)  &
          -k(116)*n(idx_HEj)  &
          +k(8)*n(idx_OH)  &
          -k(538)*n(idx_HEHj)  &
          -k(235)  &
          +k(6)*n(idx_HOCj)  &
          -k(11)*n(idx_H)
      pdj(7) =  &
          -k(956)*n(idx_C)  &
          +k(3)*n(idx_CH)  &
          -k(1201)*n(idx_C)
      pdj(8) =  &
          +k(967)*n(idx_OH)  &
          +k(959)*n(idx_CH)  &
          +k(542)*n(idx_NHj)  &
          +k(957)*n(idx_CH2)  &
          +k(966)*n(idx_O)  &
          +k(964)*n(idx_O2)  &
          +k(5)*n(idx_H2O)  &
          +k(531)*n(idx_CH2j)  &
          +k(961)*n(idx_N)  &
          +k(543)*n(idx_NH2j)  &
          +k(547)*n(idx_SIH4j)  &
          +k(958)*n(idx_CH3)  &
          +2.d0*k(236)  &
          +k(535)*n(idx_H2Oj)  &
          +k(517)*n(idx_H2j)  &
          +k(234)  &
          +k(534)*n(idx_COj)  &
          +k(530)*n(idx_CHj)  &
          +k(529)*n(idx_Cj)  &
          +k(540)*n(idx_N2j)  &
          +2.d0*k(9)*n(idx_E)  &
          +k(960)*n(idx_CN)  &
          -k(11)*n(idx_H)  &
          +k(544)*n(idx_Oj)  &
          +k(537)*n(idx_HEj)  &
          +k(3)*n(idx_CH)  &
          +k(956)*n(idx_C)  &
          +k(539)*n(idx_Nj)  &
          +k(532)*n(idx_CNj)  &
          +4.d0*k(4)*n(idx_H2)  &
          +k(533)*n(idx_COj)  &
          +k(536)*n(idx_HCNj)  &
          +k(548)*n(idx_SIOj)  &
          +3.d0*k(11)*n(idx_H)  &
          +k(962)*n(idx_NH2)  &
          +k(8)*n(idx_OH)  &
          +k(546)*n(idx_OHj)  &
          +k(963)*n(idx_NH)
      pdj(9) =  &
          +k(967)*n(idx_OH)  &
          -k(5)*n(idx_H2O)
      pdj(10) =  &
          -k(8)*n(idx_OH)  &
          +k(966)*n(idx_O)  &
          +2.d0*k(965)*n(idx_O2)  &
          -k(967)*n(idx_OH)  &
          +k(5)*n(idx_H2O)
      pdj(11) =  &
          -k(964)*n(idx_O2)  &
          -k(965)*n(idx_O2)  &
          -k(7)*n(idx_O2)  &
          +k(545)*n(idx_O2Hj)
      pdj(12) =  &
          +k(959)*n(idx_CH)  &
          +k(1201)*n(idx_C)  &
          -k(957)*n(idx_CH2)
      pdj(16) =  &
          +k(962)*n(idx_NH2)
      pdj(24) =  &
          -k(960)*n(idx_CN)
      pdj(27) =  &
          -k(962)*n(idx_NH2)  &
          +k(963)*n(idx_NH)
      pdj(28) =  &
          +k(957)*n(idx_CH2)  &
          -k(958)*n(idx_CH3)  &
          +k(1202)*n(idx_CH)
      pdj(29) =  &
          +k(958)*n(idx_CH3)
      pdj(30) =  &
          +k(541)*n(idx_NHj)  &
          -k(961)*n(idx_N)
      pdj(31) =  &
          -k(963)*n(idx_NH)  &
          +k(961)*n(idx_N)
      pdj(35) =  &
          +k(116)*n(idx_HEj)  &
          +k(537)*n(idx_HEj)  &
          +k(538)*n(idx_HEHj)
      pdj(43) =  &
          +k(964)*n(idx_O2)
      pdj(70) =  &
          +k(6)*n(idx_HOCj)  &
          +k(533)*n(idx_COj)
      pdj(71) =  &
          +k(234)  &
          +k(537)*n(idx_HEj)
      pdj(72) =  &
          +k(534)*n(idx_COj)  &
          -k(6)*n(idx_HOCj)
      pdj(73) =  &
          -k(1200)*n(idx_Cj)  &
          -k(529)*n(idx_Cj)
      pdj(74) =  &
          +k(530)*n(idx_CHj)  &
          +k(1200)*n(idx_Cj)  &
          -k(531)*n(idx_CH2j)
      pdj(75) =  &
          +k(529)*n(idx_Cj)  &
          -k(530)*n(idx_CHj)
      pdj(78) =  &
          +k(543)*n(idx_NH2j)
      pdj(80) =  &
          -k(1203)*n(idx_SIj)
      pdj(84) =  &
          +k(1203)*n(idx_SIj)
      pdj(85) =  &
          -k(1205)*n(idx_SIH3j)  &
          +k(1204)*n(idx_SIHj)
      pdj(86) =  &
          -k(532)*n(idx_CNj)
      pdj(87) =  &
          -k(534)*n(idx_COj)  &
          -k(533)*n(idx_COj)
      pdj(88) =  &
          -k(540)*n(idx_N2j)
      pdj(90) =  &
          -k(535)*n(idx_H2Oj)  &
          +k(546)*n(idx_OHj)
      pdj(91) =  &
          +k(542)*n(idx_NHj)  &
          -k(543)*n(idx_NH2j)
      pdj(92) =  &
          -k(544)*n(idx_Oj)
      pdj(93) =  &
          +k(544)*n(idx_Oj)  &
          -k(546)*n(idx_OHj)
      pdj(94) =  &
          +k(531)*n(idx_CH2j)
      pdj(96) =  &
          -k(539)*n(idx_Nj)
      pdj(97) =  &
          +k(532)*n(idx_CNj)  &
          -k(536)*n(idx_HCNj)
      pdj(98) =  &
          +k(539)*n(idx_Nj)  &
          -k(542)*n(idx_NHj)  &
          -k(541)*n(idx_NHj)
      pdj(99) =  &
          -k(547)*n(idx_SIH4j)
      pdj(100) =  &
          -k(1204)*n(idx_SIHj)
      pdj(101) =  &
          -k(548)*n(idx_SIOj)
      pdj(102) =  &
          +k(116)*n(idx_HEj)  &
          +k(235)  &
          -k(517)*n(idx_H2j)
      pdj(103) =  &
          -k(537)*n(idx_HEj)  &
          -k(116)*n(idx_HEj)
      pdj(106) =  &
          +k(541)*n(idx_NHj)  &
          +k(538)*n(idx_HEHj)  &
          +k(545)*n(idx_O2Hj)  &
          +k(517)*n(idx_H2j)
      pdj(108) =  &
          +k(535)*n(idx_H2Oj)
      pdj(109) =  &
          +k(536)*n(idx_HCNj)
      pdj(111) =  &
          -k(538)*n(idx_HEHj)
      pdj(112) =  &
          +k(540)*n(idx_N2j)
      pdj(113) =  &
          -k(545)*n(idx_O2Hj)
      pdj(114) =  &
          +k(1205)*n(idx_SIH3j)  &
          +k(547)*n(idx_SIH4j)
      pdj(115) =  &
          +k(548)*n(idx_SIOj)
    elseif(j==7) then
      pdj(1) =  &
          +k(241)  &
          +k(1108)  &
          +k(232)
      pdj(2) =  &
          +k(871)*n(idx_NH)  &
          +k(865)*n(idx_HCO)  &
          +k(956)*n(idx_H2)  &
          +2.d0*k(864)*n(idx_CH2)  &
          +k(1207)*n(idx_H)  &
          +k(877)*n(idx_OH)  &
          +k(869)*n(idx_NH2)
      pdj(3) =  &
          +k(874)*n(idx_O2)  &
          +k(394)*n(idx_OHj)  &
          -k(1197)*n(idx_O)  &
          +k(392)*n(idx_O2j)  &
          +k(877)*n(idx_OH)  &
          +k(872)*n(idx_NO)
      pdj(4) =  &
          +k(868)*n(idx_NH2)  &
          +k(1005)*n(idx_HNCO)
      pdj(5) =  &
          +k(867)*n(idx_NH2)
      pdj(6) =  &
          +k(577)*n(idx_H3j)  &
          -k(1201)*n(idx_H2)  &
          -k(956)*n(idx_H2)  &
          +k(385)*n(idx_H3Oj)
      pdj(7) =  &
          -k(874)*n(idx_O2)  &
          -k(29)*n(idx_COj)  &
          -k(872)*n(idx_NO)  &
          -k(876)*n(idx_OH)  &
          -k(1251)  &
          -k(390)*n(idx_N2Hj)  &
          -k(384)*n(idx_H2Oj)  &
          -k(391)*n(idx_NHj)  &
          -k(1207)*n(idx_H)  &
          -k(1108)  &
          -k(877)*n(idx_OH)  &
          -k(869)*n(idx_NH2)  &
          -k(232)  &
          -k(140)*n(idx_HEj)  &
          -k(1195)*n(idx_N)  &
          -k(386)*n(idx_HCNj)  &
          -k(395)*n(idx_SIHj)  &
          -k(388)*n(idx_HCO2j)  &
          -k(577)*n(idx_H3j)  &
          -k(387)*n(idx_HCOj)  &
          -k(956)*n(idx_H2)  &
          -k(1201)*n(idx_H2)  &
          -k(871)*n(idx_NH)  &
          -k(389)*n(idx_HNOj)  &
          -k(864)*n(idx_CH2)  &
          -k(1196)*n(idx_Oj)  &
          -k(241)  &
          -k(385)*n(idx_H3Oj)  &
          -k(510)*n(idx_H2j)  &
          -k(1197)*n(idx_O)  &
          -k(1005)*n(idx_HNCO)  &
          -k(393)*n(idx_O2Hj)  &
          -k(867)*n(idx_NH2)  &
          -k(870)*n(idx_NH)  &
          -k(865)*n(idx_HCO)  &
          -k(392)*n(idx_O2j)  &
          -k(873)*n(idx_NO)  &
          -k(878)*n(idx_SIH)  &
          -k(30)*n(idx_N2j)  &
          -k(28)*n(idx_CNj)  &
          -k(875)*n(idx_OCN)  &
          -k(868)*n(idx_NH2)  &
          -k(866)*n(idx_N2)  &
          -k(394)*n(idx_OHj)  &
          -k(396)*n(idx_SIOj)  &
          -k(31)*n(idx_O2j)
      pdj(8) =  &
          +k(510)*n(idx_H2j)  &
          +k(878)*n(idx_SIH)  &
          +k(876)*n(idx_OH)  &
          +k(867)*n(idx_NH2)  &
          +k(395)*n(idx_SIHj)  &
          +k(868)*n(idx_NH2)  &
          -k(1207)*n(idx_H)  &
          +k(870)*n(idx_NH)  &
          +k(956)*n(idx_H2)
      pdj(10) =  &
          +k(384)*n(idx_H2Oj)  &
          -k(877)*n(idx_OH)  &
          -k(876)*n(idx_OH)
      pdj(11) =  &
          -k(874)*n(idx_O2)  &
          +k(31)*n(idx_O2j)  &
          +k(393)*n(idx_O2Hj)
      pdj(12) =  &
          -k(864)*n(idx_CH2)  &
          +k(1201)*n(idx_H2)
      pdj(14) =  &
          -k(865)*n(idx_HCO)
      pdj(17) =  &
          +k(389)*n(idx_HNOj)  &
          -k(873)*n(idx_NO)  &
          -k(872)*n(idx_NO)
      pdj(21) =  &
          +k(878)*n(idx_SIH)
      pdj(24) =  &
          +k(875)*n(idx_OCN)  &
          +k(386)*n(idx_HCNj)  &
          +k(1195)*n(idx_N)  &
          +k(28)*n(idx_CNj)  &
          +k(872)*n(idx_NO)  &
          +k(870)*n(idx_NH)  &
          +k(866)*n(idx_N2)
      pdj(25) =  &
          +k(874)*n(idx_O2)  &
          +k(875)*n(idx_OCN)  &
          +k(865)*n(idx_HCO)  &
          +k(876)*n(idx_OH)  &
          +k(396)*n(idx_SIOj)  &
          +k(1005)*n(idx_HNCO)  &
          +k(1197)*n(idx_O)  &
          +k(873)*n(idx_NO)  &
          +k(387)*n(idx_HCOj)  &
          +k(29)*n(idx_COj)
      pdj(26) =  &
          -k(866)*n(idx_N2)  &
          +k(30)*n(idx_N2j)  &
          +k(390)*n(idx_N2Hj)
      pdj(27) =  &
          -k(869)*n(idx_NH2)  &
          -k(867)*n(idx_NH2)  &
          -k(868)*n(idx_NH2)
      pdj(30) =  &
          +k(871)*n(idx_NH)  &
          +k(873)*n(idx_NO)  &
          +k(391)*n(idx_NHj)  &
          -k(1195)*n(idx_N)  &
          +k(866)*n(idx_N2)
      pdj(31) =  &
          -k(870)*n(idx_NH)  &
          -k(871)*n(idx_NH)  &
          +k(869)*n(idx_NH2)
      pdj(33) =  &
          -k(878)*n(idx_SIH)
      pdj(35) =  &
          +k(140)*n(idx_HEj)
      pdj(38) =  &
          +k(388)*n(idx_HCO2j)
      pdj(41) =  &
          -k(1005)*n(idx_HNCO)
      pdj(44) =  &
          -k(875)*n(idx_OCN)
      pdj(53) =  &
          +k(1251)
      pdj(70) =  &
          -k(387)*n(idx_HCOj)  &
          +k(385)*n(idx_H3Oj)
      pdj(73) =  &
          +k(1108)  &
          +k(30)*n(idx_N2j)  &
          +k(31)*n(idx_O2j)  &
          +k(28)*n(idx_CNj)  &
          +k(140)*n(idx_HEj)  &
          +k(241)  &
          +k(29)*n(idx_COj)  &
          +k(232)
      pdj(75) =  &
          +k(510)*n(idx_H2j)  &
          +k(386)*n(idx_HCNj)  &
          +k(394)*n(idx_OHj)  &
          +k(577)*n(idx_H3j)  &
          +k(391)*n(idx_NHj)  &
          +k(384)*n(idx_H2Oj)  &
          +k(388)*n(idx_HCO2j)  &
          +k(393)*n(idx_O2Hj)  &
          +k(389)*n(idx_HNOj)  &
          +k(387)*n(idx_HCOj)  &
          +k(390)*n(idx_N2Hj)
      pdj(80) =  &
          +k(396)*n(idx_SIOj)
      pdj(83) =  &
          +k(395)*n(idx_SIHj)
      pdj(86) =  &
          -k(28)*n(idx_CNj)
      pdj(87) =  &
          -k(29)*n(idx_COj)  &
          +k(1196)*n(idx_Oj)  &
          +k(392)*n(idx_O2j)
      pdj(88) =  &
          -k(30)*n(idx_N2j)
      pdj(89) =  &
          -k(392)*n(idx_O2j)  &
          -k(31)*n(idx_O2j)
      pdj(90) =  &
          -k(384)*n(idx_H2Oj)
      pdj(92) =  &
          -k(1196)*n(idx_Oj)
      pdj(93) =  &
          -k(394)*n(idx_OHj)
      pdj(97) =  &
          -k(386)*n(idx_HCNj)
      pdj(98) =  &
          -k(391)*n(idx_NHj)
      pdj(100) =  &
          -k(395)*n(idx_SIHj)
      pdj(101) =  &
          -k(396)*n(idx_SIOj)
      pdj(102) =  &
          -k(510)*n(idx_H2j)
      pdj(103) =  &
          -k(140)*n(idx_HEj)
      pdj(104) =  &
          -k(389)*n(idx_HNOj)
      pdj(106) =  &
          -k(577)*n(idx_H3j)
      pdj(108) =  &
          -k(385)*n(idx_H3Oj)
      pdj(110) =  &
          -k(388)*n(idx_HCO2j)
      pdj(112) =  &
          -k(390)*n(idx_N2Hj)
      pdj(113) =  &
          -k(393)*n(idx_O2Hj)
    elseif(j==8) then
      pdj(1) =  &
          +k(259)  &
          +k(237)
      pdj(2) =  &
          +k(968)*n(idx_CH2)  &
          +k(1207)*n(idx_C)  &
          -k(971)*n(idx_CH)  &
          -k(10)*n(idx_CH)
      pdj(3) =  &
          +k(132)*n(idx_Oj)  &
          +k(990)*n(idx_O2)  &
          +k(988)*n(idx_NO)  &
          +k(979)*n(idx_HCO)  &
          +k(997)*n(idx_OH)  &
          +k(14)*n(idx_OH)  &
          +2.d0*k(13)*n(idx_O2)  &
          +k(981)*n(idx_HNO)  &
          +k(991)*n(idx_O2H)  &
          -k(1208)*n(idx_O)  &
          +k(994)*n(idx_OCN)
      pdj(4) =  &
          -k(980)*n(idx_HNC)
      pdj(5) =  &
          +k(980)*n(idx_HNC)  &
          +k(994)*n(idx_OCN)  &
          +k(974)*n(idx_H2CN)  &
          +k(130)*n(idx_HCNj)  &
          -k(977)*n(idx_HCN)
      pdj(6) =  &
          +k(615)*n(idx_CHj)  &
          +k(971)*n(idx_CH)  &
          +k(129)*n(idx_H2j)  &
          +k(976)*n(idx_H2O)  &
          +k(984)*n(idx_NH2)  &
          +k(997)*n(idx_OH)  &
          +k(975)*n(idx_H2CO)  &
          +k(974)*n(idx_H2CN)  &
          +k(986)*n(idx_NH)  &
          +k(968)*n(idx_CH2)  &
          +k(616)*n(idx_CH2j)  &
          +k(977)*n(idx_HCN)  &
          +k(970)*n(idx_CH4)  &
          +k(620)*n(idx_SIHj)  &
          +k(982)*n(idx_HNO)  &
          +k(978)*n(idx_HCO)  &
          -k(11)*n(idx_H2)  &
          +k(985)*n(idx_NH3)  &
          +k(618)*n(idx_CH4j)  &
          +k(617)*n(idx_CH3j)  &
          +k(992)*n(idx_O2H)  &
          +k(969)*n(idx_CH3)
      pdj(7) =  &
          +k(973)*n(idx_CO)  &
          -k(1207)*n(idx_C)  &
          +k(971)*n(idx_CH)  &
          +k(10)*n(idx_CH)
      pdj(8) =  &
          -k(132)*n(idx_Oj)  &
          -k(982)*n(idx_HNO)  &
          -k(976)*n(idx_H2O)  &
          -k(993)*n(idx_O2H)  &
          -k(617)*n(idx_CH3j)  &
          -k(971)*n(idx_CH)  &
          -k(131)*n(idx_HEj)  &
          -k(985)*n(idx_NH3)  &
          -k(12)*n(idx_H2O)  &
          -k(129)*n(idx_H2j)  &
          -k(1208)*n(idx_O)  &
          -k(259)  &
          -k(10)*n(idx_CH)  &
          -k(977)*n(idx_HCN)  &
          +2.d0*k(10)*n(idx_CH)  &
          -k(983)*n(idx_HNO)  &
          -k(969)*n(idx_CH3)  &
          -k(970)*n(idx_CH4)  &
          -k(996)*n(idx_OCN)  &
          -k(974)*n(idx_H2CN)  &
          -k(1207)*n(idx_C)  &
          -k(972)*n(idx_CO2)  &
          -k(13)*n(idx_O2)  &
          -k(1198)*n(idx_Hj)  &
          -k(1206)*n(idx_Cj)  &
          -k(615)*n(idx_CHj)  &
          -k(994)*n(idx_OCN)  &
          -k(981)*n(idx_HNO)  &
          -k(991)*n(idx_O2H)  &
          -k(992)*n(idx_O2H)  &
          -k(619)*n(idx_HEHj)  &
          -k(968)*n(idx_CH2)  &
          -k(616)*n(idx_CH2j)  &
          +k(13)*n(idx_O2)  &
          -k(975)*n(idx_H2CO)  &
          -k(620)*n(idx_SIHj)  &
          -k(979)*n(idx_HCO)  &
          -k(987)*n(idx_NO2)  &
          -k(978)*n(idx_HCO)  &
          +3.d0*k(11)*n(idx_H2)  &
          -k(128)*n(idx_COj)  &
          -k(127)*n(idx_CNj)  &
          -k(1210)*n(idx_SIj)  &
          -k(986)*n(idx_NH)  &
          -k(984)*n(idx_NH2)  &
          -k(237)  &
          -k(995)*n(idx_OCN)  &
          -k(618)*n(idx_CH4j)  &
          -k(14)*n(idx_OH)  &
          -k(988)*n(idx_NO)  &
          -k(11)*n(idx_H2)  &
          -k(973)*n(idx_CO)  &
          -k(130)*n(idx_HCNj)  &
          +k(980)*n(idx_HNC)  &
          -k(980)*n(idx_HNC)  &
          +2.d0*k(14)*n(idx_OH)  &
          -k(1209)*n(idx_OH)  &
          -k(997)*n(idx_OH)  &
          -k(990)*n(idx_O2)  &
          +2.d0*k(12)*n(idx_H2O)  &
          -k(989)*n(idx_NO)
      pdj(9) =  &
          -k(12)*n(idx_H2O)  &
          +k(991)*n(idx_O2H)  &
          +k(1209)*n(idx_OH)  &
          -k(976)*n(idx_H2O)
      pdj(10) =  &
          -k(997)*n(idx_OH)  &
          -k(14)*n(idx_OH)  &
          +k(973)*n(idx_CO)  &
          +k(990)*n(idx_O2)  &
          +k(976)*n(idx_H2O)  &
          +k(972)*n(idx_CO2)  &
          -k(1209)*n(idx_OH)  &
          +k(996)*n(idx_OCN)  &
          +2.d0*k(993)*n(idx_O2H)  &
          +k(1208)*n(idx_O)  &
          +k(983)*n(idx_HNO)  &
          +k(12)*n(idx_H2O)  &
          +k(989)*n(idx_NO)  &
          +k(987)*n(idx_NO2)
      pdj(11) =  &
          -k(13)*n(idx_O2)  &
          -k(990)*n(idx_O2)  &
          +k(992)*n(idx_O2H)
      pdj(12) =  &
          -k(968)*n(idx_CH2)  &
          +k(969)*n(idx_CH3)  &
          +k(979)*n(idx_HCO)
      pdj(13) =  &
          -k(975)*n(idx_H2CO)
      pdj(14) =  &
          -k(978)*n(idx_HCO)  &
          +k(975)*n(idx_H2CO)  &
          -k(979)*n(idx_HCO)
      pdj(16) =  &
          -k(985)*n(idx_NH3)
      pdj(17) =  &
          -k(988)*n(idx_NO)  &
          +k(982)*n(idx_HNO)  &
          -k(989)*n(idx_NO)  &
          +k(987)*n(idx_NO2)
      pdj(24) =  &
          +k(977)*n(idx_HCN)  &
          +k(996)*n(idx_OCN)  &
          +k(127)*n(idx_CNj)
      pdj(25) =  &
          +k(978)*n(idx_HCO)  &
          -k(973)*n(idx_CO)  &
          +k(972)*n(idx_CO2)  &
          +k(128)*n(idx_COj)  &
          +k(995)*n(idx_OCN)
      pdj(27) =  &
          +k(981)*n(idx_HNO)  &
          -k(984)*n(idx_NH2)  &
          +k(985)*n(idx_NH3)
      pdj(28) =  &
          +k(970)*n(idx_CH4)  &
          -k(969)*n(idx_CH3)
      pdj(29) =  &
          -k(970)*n(idx_CH4)
      pdj(30) =  &
          +k(989)*n(idx_NO)  &
          +k(986)*n(idx_NH)
      pdj(31) =  &
          +k(988)*n(idx_NO)  &
          +k(984)*n(idx_NH2)  &
          +k(995)*n(idx_OCN)  &
          -k(986)*n(idx_NH)  &
          +k(983)*n(idx_HNO)
      pdj(35) =  &
          +k(131)*n(idx_HEj)  &
          +k(619)*n(idx_HEHj)
      pdj(36) =  &
          -k(982)*n(idx_HNO)  &
          -k(981)*n(idx_HNO)  &
          -k(983)*n(idx_HNO)
      pdj(38) =  &
          -k(972)*n(idx_CO2)
      pdj(39) =  &
          -k(974)*n(idx_H2CN)
      pdj(42) =  &
          -k(987)*n(idx_NO2)
      pdj(43) =  &
          -k(993)*n(idx_O2H)  &
          -k(991)*n(idx_O2H)  &
          -k(992)*n(idx_O2H)
      pdj(44) =  &
          -k(994)*n(idx_OCN)  &
          -k(995)*n(idx_OCN)  &
          -k(996)*n(idx_OCN)
      pdj(71) =  &
          +k(132)*n(idx_Oj)  &
          -k(1198)*n(idx_Hj)  &
          +k(129)*n(idx_H2j)  &
          +k(127)*n(idx_CNj)  &
          +k(259)  &
          +k(128)*n(idx_COj)  &
          +k(131)*n(idx_HEj)  &
          +k(130)*n(idx_HCNj)  &
          +k(237)
      pdj(73) =  &
          +k(615)*n(idx_CHj)  &
          -k(1206)*n(idx_Cj)
      pdj(74) =  &
          -k(616)*n(idx_CH2j)  &
          +k(617)*n(idx_CH3j)
      pdj(75) =  &
          +k(1206)*n(idx_Cj)  &
          +k(616)*n(idx_CH2j)  &
          -k(615)*n(idx_CHj)
      pdj(80) =  &
          -k(1210)*n(idx_SIj)  &
          +k(620)*n(idx_SIHj)
      pdj(86) =  &
          -k(127)*n(idx_CNj)
      pdj(87) =  &
          -k(128)*n(idx_COj)
      pdj(92) =  &
          -k(132)*n(idx_Oj)
      pdj(94) =  &
          +k(618)*n(idx_CH4j)  &
          -k(617)*n(idx_CH3j)
      pdj(95) =  &
          -k(618)*n(idx_CH4j)
      pdj(97) =  &
          -k(130)*n(idx_HCNj)
      pdj(100) =  &
          +k(1210)*n(idx_SIj)  &
          -k(620)*n(idx_SIHj)
      pdj(102) =  &
          +k(619)*n(idx_HEHj)  &
          -k(129)*n(idx_H2j)  &
          +k(1198)*n(idx_Hj)
      pdj(103) =  &
          -k(131)*n(idx_HEj)
      pdj(111) =  &
          -k(619)*n(idx_HEHj)
    elseif(j==9) then
      pdj(1) =  &
          +k(1142)
      pdj(3) =  &
          +k(841)*n(idx_OHj)  &
          -k(1063)*n(idx_O)  &
          +k(757)*n(idx_NHj)  &
          +k(211)*n(idx_Oj)
      pdj(5) =  &
          +k(125)*n(idx_HCNj)
      pdj(6) =  &
          +k(107)*n(idx_H2j)  &
          +k(587)*n(idx_H3j)  &
          +k(405)*n(idx_CHj)  &
          +k(5)*n(idx_H2)  &
          -k(5)*n(idx_H2)  &
          +k(976)*n(idx_H)  &
          +k(756)*n(idx_NHj)
      pdj(7) =  &
          +k(404)*n(idx_CHj)
      pdj(8) =  &
          +k(1143)  &
          +k(257)  &
          +k(371)*n(idx_Cj)  &
          -k(976)*n(idx_H)  &
          +k(519)*n(idx_H2j)  &
          +k(5)*n(idx_H2)  &
          +k(674)*n(idx_HEj)  &
          +2.d0*k(12)*n(idx_H)  &
          +k(573)*n(idx_SIj)  &
          +k(419)*n(idx_CH2j)  &
          +k(403)*n(idx_CHj)  &
          +k(81)*n(idx_Hj)  &
          -k(12)*n(idx_H)  &
          +k(372)*n(idx_Cj)
      pdj(9) =  &
          -k(674)*n(idx_HEj)  &
          -k(905)*n(idx_CH3)  &
          -k(758)*n(idx_NHj)  &
          -k(124)*n(idx_COj)  &
          -k(81)*n(idx_Hj)  &
          -k(5)*n(idx_H2)  &
          -k(565)*n(idx_H3COj)  &
          -k(561)*n(idx_CNj)  &
          -k(755)*n(idx_NHj)  &
          -k(675)*n(idx_HEj)  &
          -k(403)*n(idx_CHj)  &
          -k(404)*n(idx_CHj)  &
          -k(405)*n(idx_CHj)  &
          -k(519)*n(idx_H2j)  &
          -k(570)*n(idx_N2j)  &
          -k(126)*n(idx_N2j)  &
          -k(976)*n(idx_H)  &
          -k(772)*n(idx_NH2j)  &
          -k(756)*n(idx_NHj)  &
          -k(1258)  &
          -k(562)*n(idx_CNj)  &
          -k(563)*n(idx_COj)  &
          -k(573)*n(idx_SIj)  &
          -k(177)*n(idx_NHj)  &
          -k(257)  &
          -k(1143)  &
          -k(1037)*n(idx_NH)  &
          -k(107)*n(idx_H2j)  &
          -k(841)*n(idx_OHj)  &
          -k(125)*n(idx_HCNj)  &
          -k(773)*n(idx_NH2j)  &
          -k(757)*n(idx_NHj)  &
          -k(556)*n(idx_H2Oj)  &
          -k(419)*n(idx_CH2j)  &
          -k(1063)*n(idx_O)  &
          -k(372)*n(idx_Cj)  &
          -k(566)*n(idx_HCNj)  &
          -k(568)*n(idx_HCO2j)  &
          -k(12)*n(idx_H)  &
          -k(221)*n(idx_OHj)  &
          -k(567)*n(idx_HCOj)  &
          -k(569)*n(idx_HNOj)  &
          -k(575)*n(idx_SIH4j)  &
          -k(587)*n(idx_H3j)  &
          -k(161)*n(idx_Nj)  &
          -k(144)*n(idx_HEj)  &
          -k(574)*n(idx_SIHj)  &
          -k(371)*n(idx_Cj)  &
          -k(576)*n(idx_SIH5j)  &
          -k(211)*n(idx_Oj)  &
          -k(571)*n(idx_N2Hj)  &
          -k(451)*n(idx_CH4j)  &
          -k(1142)  &
          -k(564)*n(idx_H2COj)  &
          -k(572)*n(idx_O2Hj)
      pdj(10) =  &
          +k(556)*n(idx_H2Oj)  &
          +k(570)*n(idx_N2j)  &
          +k(257)  &
          +k(1143)  &
          +k(773)*n(idx_NH2j)  &
          +k(1037)*n(idx_NH)  &
          +k(905)*n(idx_CH3)  &
          +k(675)*n(idx_HEj)  &
          +k(5)*n(idx_H2)  &
          +k(563)*n(idx_COj)  &
          +k(12)*n(idx_H)  &
          +k(221)*n(idx_OHj)  &
          +k(561)*n(idx_CNj)  &
          +k(976)*n(idx_H)  &
          +k(758)*n(idx_NHj)  &
          +2.d0*k(1063)*n(idx_O)
      pdj(11) =  &
          +k(572)*n(idx_O2Hj)
      pdj(13) =  &
          +k(565)*n(idx_H3COj)
      pdj(14) =  &
          +k(564)*n(idx_H2COj)
      pdj(17) =  &
          +k(569)*n(idx_HNOj)
      pdj(18) =  &
          +k(574)*n(idx_SIHj)
      pdj(23) =  &
          +k(575)*n(idx_SIH4j)
      pdj(24) =  &
          +k(566)*n(idx_HCNj)
      pdj(25) =  &
          +k(124)*n(idx_COj)  &
          +k(567)*n(idx_HCOj)
      pdj(26) =  &
          +k(571)*n(idx_N2Hj)  &
          +k(126)*n(idx_N2j)
      pdj(27) =  &
          +k(1037)*n(idx_NH)
      pdj(28) =  &
          +k(451)*n(idx_CH4j)  &
          -k(905)*n(idx_CH3)
      pdj(29) =  &
          +k(905)*n(idx_CH3)
      pdj(30) =  &
          +k(161)*n(idx_Nj)  &
          +k(755)*n(idx_NHj)
      pdj(31) =  &
          +k(177)*n(idx_NHj)  &
          +k(772)*n(idx_NH2j)  &
          +k(562)*n(idx_CNj)  &
          -k(1037)*n(idx_NH)
      pdj(32) =  &
          +k(576)*n(idx_SIH5j)
      pdj(35) =  &
          +k(674)*n(idx_HEj)  &
          +k(675)*n(idx_HEj)  &
          +k(144)*n(idx_HEj)
      pdj(38) =  &
          +k(568)*n(idx_HCO2j)
      pdj(55) =  &
          +k(1258)
      pdj(70) =  &
          +k(405)*n(idx_CHj)  &
          -k(567)*n(idx_HCOj)  &
          +k(563)*n(idx_COj)  &
          +k(371)*n(idx_Cj)  &
          +k(562)*n(idx_CNj)
      pdj(71) =  &
          +k(675)*n(idx_HEj)  &
          -k(81)*n(idx_Hj)
      pdj(72) =  &
          +k(372)*n(idx_Cj)
      pdj(73) =  &
          -k(372)*n(idx_Cj)  &
          -k(371)*n(idx_Cj)
      pdj(74) =  &
          -k(419)*n(idx_CH2j)
      pdj(75) =  &
          -k(403)*n(idx_CHj)  &
          -k(405)*n(idx_CHj)  &
          -k(404)*n(idx_CHj)
      pdj(76) =  &
          +k(403)*n(idx_CHj)  &
          -k(564)*n(idx_H2COj)
      pdj(78) =  &
          +k(773)*n(idx_NH2j)  &
          +k(757)*n(idx_NHj)
      pdj(80) =  &
          -k(573)*n(idx_SIj)
      pdj(86) =  &
          -k(561)*n(idx_CNj)  &
          -k(562)*n(idx_CNj)
      pdj(87) =  &
          -k(563)*n(idx_COj)  &
          -k(124)*n(idx_COj)
      pdj(88) =  &
          -k(570)*n(idx_N2j)  &
          -k(126)*n(idx_N2j)
      pdj(90) =  &
          +k(107)*n(idx_H2j)  &
          +k(177)*n(idx_NHj)  &
          +k(1142)  &
          +k(144)*n(idx_HEj)  &
          +k(161)*n(idx_Nj)  &
          -k(556)*n(idx_H2Oj)  &
          +k(126)*n(idx_N2j)  &
          +k(221)*n(idx_OHj)  &
          +k(125)*n(idx_HCNj)  &
          +k(81)*n(idx_Hj)  &
          +k(211)*n(idx_Oj)  &
          +k(124)*n(idx_COj)
      pdj(91) =  &
          -k(772)*n(idx_NH2j)  &
          -k(773)*n(idx_NH2j)  &
          +k(758)*n(idx_NHj)
      pdj(92) =  &
          -k(211)*n(idx_Oj)
      pdj(93) =  &
          -k(221)*n(idx_OHj)  &
          +k(674)*n(idx_HEj)  &
          -k(841)*n(idx_OHj)
      pdj(95) =  &
          -k(451)*n(idx_CH4j)
      pdj(96) =  &
          -k(161)*n(idx_Nj)
      pdj(97) =  &
          -k(125)*n(idx_HCNj)  &
          +k(561)*n(idx_CNj)  &
          -k(566)*n(idx_HCNj)
      pdj(98) =  &
          -k(756)*n(idx_NHj)  &
          -k(177)*n(idx_NHj)  &
          -k(758)*n(idx_NHj)  &
          -k(755)*n(idx_NHj)  &
          -k(757)*n(idx_NHj)
      pdj(99) =  &
          -k(575)*n(idx_SIH4j)
      pdj(100) =  &
          -k(574)*n(idx_SIHj)
      pdj(102) =  &
          -k(519)*n(idx_H2j)  &
          -k(107)*n(idx_H2j)
      pdj(103) =  &
          -k(674)*n(idx_HEj)  &
          -k(675)*n(idx_HEj)  &
          -k(144)*n(idx_HEj)
      pdj(104) =  &
          +k(756)*n(idx_NHj)  &
          -k(569)*n(idx_HNOj)
      pdj(106) =  &
          -k(587)*n(idx_H3j)
      pdj(107) =  &
          +k(419)*n(idx_CH2j)  &
          -k(565)*n(idx_H3COj)
      pdj(108) =  &
          +k(587)*n(idx_H3j)  &
          +k(556)*n(idx_H2Oj)  &
          +k(575)*n(idx_SIH4j)  &
          +k(567)*n(idx_HCOj)  &
          +k(564)*n(idx_H2COj)  &
          +k(841)*n(idx_OHj)  &
          +k(519)*n(idx_H2j)  &
          +k(451)*n(idx_CH4j)  &
          +k(772)*n(idx_NH2j)  &
          +k(568)*n(idx_HCO2j)  &
          +k(755)*n(idx_NHj)  &
          +k(569)*n(idx_HNOj)  &
          +k(404)*n(idx_CHj)  &
          +k(572)*n(idx_O2Hj)  &
          +k(565)*n(idx_H3COj)  &
          +k(571)*n(idx_N2Hj)  &
          +k(566)*n(idx_HCNj)  &
          +k(576)*n(idx_SIH5j)  &
          +k(574)*n(idx_SIHj)
      pdj(110) =  &
          -k(568)*n(idx_HCO2j)
      pdj(112) =  &
          +k(570)*n(idx_N2j)  &
          -k(571)*n(idx_N2Hj)
      pdj(113) =  &
          -k(572)*n(idx_O2Hj)
      pdj(114) =  &
          -k(576)*n(idx_SIH5j)
      pdj(115) =  &
          +k(573)*n(idx_SIj)
    elseif(j==10) then
      pdj(1) =  &
          +k(1176)
      pdj(2) =  &
          +k(900)*n(idx_CH2)  &
          -k(942)*n(idx_CH)  &
          +k(877)*n(idx_C)
      pdj(3) =  &
          +k(848)*n(idx_OHj)  &
          +k(1091)*n(idx_CN)  &
          +k(852)*n(idx_COj)  &
          +2.d0*k(1102)*n(idx_OH)  &
          +k(1033)*n(idx_NH2)  &
          +k(877)*n(idx_C)  &
          +k(285)  &
          +k(918)*n(idx_CH3)  &
          +k(997)*n(idx_H)  &
          +k(901)*n(idx_CH2)  &
          +k(1051)*n(idx_NH)  &
          -k(1081)*n(idx_O)  &
          +k(1175)  &
          +k(853)*n(idx_H2Oj)  &
          +k(216)*n(idx_Oj)  &
          +k(14)*n(idx_H)  &
          +k(8)*n(idx_H2)  &
          +k(1027)*n(idx_N)
      pdj(5) =  &
          +k(1091)*n(idx_CN)  &
          -k(1095)*n(idx_HCN)  &
          -k(1096)*n(idx_HCN)
      pdj(6) =  &
          +k(919)*n(idx_CH3)  &
          -k(967)*n(idx_H2)  &
          +k(601)*n(idx_H3j)  &
          +k(997)*n(idx_H)  &
          +k(115)*n(idx_H2j)  &
          -k(8)*n(idx_H2)  &
          +k(446)*n(idx_CH3j)  &
          +k(8)*n(idx_H2)  &
          +k(416)*n(idx_CHj)
      pdj(7) =  &
          -k(877)*n(idx_C)  &
          -k(876)*n(idx_C)
      pdj(8) =  &
          +k(528)*n(idx_H2j)  &
          +k(285)  &
          +k(899)*n(idx_CH2)  &
          -k(14)*n(idx_H)  &
          +k(860)*n(idx_SIj)  &
          +k(1100)*n(idx_NO)  &
          +k(1092)*n(idx_CN)  &
          +k(1081)*n(idx_O)  &
          +k(1175)  &
          +k(876)*n(idx_C)  &
          +k(967)*n(idx_H2)  &
          +k(1103)*n(idx_SI)  &
          +k(380)*n(idx_Cj)  &
          +2.d0*k(14)*n(idx_H)  &
          -k(1209)*n(idx_H)  &
          +k(856)*n(idx_HCOj)  &
          +k(1050)*n(idx_NH)  &
          +k(1093)*n(idx_CO)  &
          +k(8)*n(idx_H2)  &
          +k(1026)*n(idx_N)  &
          +k(820)*n(idx_Oj)  &
          +k(942)*n(idx_CH)  &
          +k(91)*n(idx_Hj)  &
          -k(997)*n(idx_H)  &
          +k(700)*n(idx_HEj)
      pdj(9) =  &
          +k(920)*n(idx_CH3)  &
          +2.d0*k(1102)*n(idx_OH)  &
          +k(1098)*n(idx_HNO)  &
          +k(1209)*n(idx_H)  &
          +k(1097)*n(idx_HCO)  &
          +k(1099)*n(idx_NH3)  &
          +k(1094)*n(idx_H2CO)  &
          +k(1032)*n(idx_NH2)  &
          +k(967)*n(idx_H2)  &
          +k(900)*n(idx_CH2)  &
          +k(923)*n(idx_CH4)  &
          +k(1049)*n(idx_NH)  &
          +k(1101)*n(idx_O2H)  &
          +k(1095)*n(idx_HCN)
      pdj(10) =  &
          -k(1103)*n(idx_SI)  &
          -k(115)*n(idx_H2j)  &
          -k(856)*n(idx_HCOj)  &
          -k(1093)*n(idx_CO)  &
          -k(528)*n(idx_H2j)  &
          -k(1255)  &
          -k(876)*n(idx_C)  &
          -k(1081)*n(idx_O)  &
          -k(380)*n(idx_Cj)  &
          -k(14)*n(idx_H)  &
          -k(1097)*n(idx_HCO)  &
          -k(1209)*n(idx_H)  &
          -k(228)*n(idx_N2j)  &
          -k(848)*n(idx_OHj)  &
          -k(857)*n(idx_HNOj)  &
          -k(91)*n(idx_Hj)  &
          -k(1092)*n(idx_CN)  &
          -k(967)*n(idx_H2)  &
          -k(942)*n(idx_CH)  &
          -k(1091)*n(idx_CN)  &
          -k(285)  &
          -k(1050)*n(idx_NH)  &
          -k(854)*n(idx_HCNj)  &
          -k(1101)*n(idx_O2H)  &
          -k(852)*n(idx_COj)  &
          -k(1032)*n(idx_NH2)  &
          -k(1094)*n(idx_H2CO)  &
          -k(853)*n(idx_H2Oj)  &
          -k(859)*n(idx_O2Hj)  &
          -k(877)*n(idx_C)  &
          -k(1100)*n(idx_NO)  &
          -k(769)*n(idx_NHj)  &
          -k(1099)*n(idx_NH3)  &
          -k(1176)  &
          -k(226)*n(idx_CNj)  &
          -k(601)*n(idx_H3j)  &
          -k(923)*n(idx_CH4)  &
          -k(900)*n(idx_CH2)  &
          -k(997)*n(idx_H)  &
          -k(1026)*n(idx_N)  &
          -k(899)*n(idx_CH2)  &
          -k(1033)*n(idx_NH2)  &
          -k(8)*n(idx_H2)  &
          -k(918)*n(idx_CH3)  &
          -k(920)*n(idx_CH3)  &
          -k(1051)*n(idx_NH)  &
          -k(700)*n(idx_HEj)  &
          -k(1095)*n(idx_HCN)  &
          -k(227)*n(idx_COj)  &
          -k(416)*n(idx_CHj)  &
          -k(901)*n(idx_CH2)  &
          -k(858)*n(idx_N2Hj)  &
          -k(1096)*n(idx_HCN)  &
          -k(820)*n(idx_Oj)  &
          -k(860)*n(idx_SIj)  &
          -k(170)*n(idx_Nj)  &
          -k(1049)*n(idx_NH)  &
          -k(1175)  &
          -k(1027)*n(idx_N)  &
          -4.d0*k(1102)*n(idx_OH)  &
          -k(216)*n(idx_Oj)  &
          -k(855)*n(idx_HCOj)  &
          -k(1098)*n(idx_HNO)  &
          -k(919)*n(idx_CH3)  &
          -k(446)*n(idx_CH3j)
      pdj(11) =  &
          +k(859)*n(idx_O2Hj)  &
          +k(1101)*n(idx_O2H)  &
          +k(1081)*n(idx_O)
      pdj(12) =  &
          +k(920)*n(idx_CH3)  &
          -k(900)*n(idx_CH2)  &
          -k(901)*n(idx_CH2)  &
          -k(899)*n(idx_CH2)
      pdj(13) =  &
          +k(919)*n(idx_CH3)  &
          +k(899)*n(idx_CH2)  &
          -k(1094)*n(idx_H2CO)
      pdj(14) =  &
          -k(1097)*n(idx_HCO)  &
          +k(942)*n(idx_CH)  &
          +k(1094)*n(idx_H2CO)
      pdj(16) =  &
          +k(1033)*n(idx_NH2)  &
          -k(1099)*n(idx_NH3)
      pdj(17) =  &
          +k(1098)*n(idx_HNO)  &
          -k(1100)*n(idx_NO)  &
          +k(1026)*n(idx_N)  &
          +k(857)*n(idx_HNOj)
      pdj(18) =  &
          -k(1103)*n(idx_SI)
      pdj(24) =  &
          -k(1092)*n(idx_CN)  &
          +k(226)*n(idx_CNj)  &
          +k(1095)*n(idx_HCN)  &
          -k(1091)*n(idx_CN)  &
          +k(854)*n(idx_HCNj)
      pdj(25) =  &
          -k(1093)*n(idx_CO)  &
          +k(1097)*n(idx_HCO)  &
          +k(876)*n(idx_C)  &
          +k(1096)*n(idx_HCN)  &
          +k(855)*n(idx_HCOj)  &
          +k(227)*n(idx_COj)
      pdj(26) =  &
          +k(858)*n(idx_N2Hj)  &
          +k(228)*n(idx_N2j)
      pdj(27) =  &
          +k(1096)*n(idx_HCN)  &
          -k(1032)*n(idx_NH2)  &
          -k(1033)*n(idx_NH2)  &
          +k(1051)*n(idx_NH)  &
          +k(1099)*n(idx_NH3)
      pdj(28) =  &
          +k(901)*n(idx_CH2)  &
          -k(919)*n(idx_CH3)  &
          +k(923)*n(idx_CH4)  &
          -k(918)*n(idx_CH3)  &
          -k(920)*n(idx_CH3)
      pdj(29) =  &
          -k(923)*n(idx_CH4)  &
          +k(918)*n(idx_CH3)
      pdj(30) =  &
          +k(769)*n(idx_NHj)  &
          +k(170)*n(idx_Nj)  &
          +k(1049)*n(idx_NH)  &
          -k(1026)*n(idx_N)  &
          -k(1027)*n(idx_N)
      pdj(31) =  &
          -k(1050)*n(idx_NH)  &
          -k(1049)*n(idx_NH)  &
          +k(1032)*n(idx_NH2)  &
          +k(1027)*n(idx_N)  &
          -k(1051)*n(idx_NH)
      pdj(34) =  &
          +k(1103)*n(idx_SI)
      pdj(35) =  &
          +k(700)*n(idx_HEj)
      pdj(36) =  &
          +k(1050)*n(idx_NH)  &
          -k(1098)*n(idx_HNO)
      pdj(38) =  &
          +k(1093)*n(idx_CO)
      pdj(42) =  &
          +k(1100)*n(idx_NO)
      pdj(43) =  &
          -k(1101)*n(idx_O2H)
      pdj(44) =  &
          +k(1092)*n(idx_CN)
      pdj(55) =  &
          +k(1255)
      pdj(70) =  &
          -k(855)*n(idx_HCOj)  &
          +k(852)*n(idx_COj)  &
          -k(856)*n(idx_HCOj)
      pdj(71) =  &
          -k(91)*n(idx_Hj)
      pdj(73) =  &
          -k(380)*n(idx_Cj)
      pdj(75) =  &
          -k(416)*n(idx_CHj)
      pdj(76) =  &
          +k(446)*n(idx_CH3j)
      pdj(80) =  &
          -k(860)*n(idx_SIj)
      pdj(86) =  &
          -k(226)*n(idx_CNj)
      pdj(87) =  &
          +k(380)*n(idx_Cj)  &
          +k(416)*n(idx_CHj)  &
          -k(852)*n(idx_COj)  &
          -k(227)*n(idx_COj)
      pdj(88) =  &
          -k(228)*n(idx_N2j)
      pdj(89) =  &
          +k(820)*n(idx_Oj)
      pdj(90) =  &
          +k(859)*n(idx_O2Hj)  &
          +k(528)*n(idx_H2j)  &
          +k(601)*n(idx_H3j)  &
          +k(848)*n(idx_OHj)  &
          -k(853)*n(idx_H2Oj)  &
          +k(854)*n(idx_HCNj)  &
          +k(855)*n(idx_HCOj)  &
          +k(769)*n(idx_NHj)  &
          +k(858)*n(idx_N2Hj)  &
          +k(857)*n(idx_HNOj)
      pdj(92) =  &
          -k(216)*n(idx_Oj)  &
          +k(700)*n(idx_HEj)  &
          -k(820)*n(idx_Oj)
      pdj(93) =  &
          -k(848)*n(idx_OHj)  &
          +k(227)*n(idx_COj)  &
          +k(170)*n(idx_Nj)  &
          +k(1176)  &
          +k(115)*n(idx_H2j)  &
          +k(91)*n(idx_Hj)  &
          +k(216)*n(idx_Oj)  &
          +k(226)*n(idx_CNj)  &
          +k(228)*n(idx_N2j)
      pdj(94) =  &
          -k(446)*n(idx_CH3j)
      pdj(96) =  &
          -k(170)*n(idx_Nj)
      pdj(97) =  &
          -k(854)*n(idx_HCNj)
      pdj(98) =  &
          -k(769)*n(idx_NHj)
      pdj(101) =  &
          +k(860)*n(idx_SIj)
      pdj(102) =  &
          -k(115)*n(idx_H2j)  &
          -k(528)*n(idx_H2j)
      pdj(103) =  &
          -k(700)*n(idx_HEj)
      pdj(104) =  &
          -k(857)*n(idx_HNOj)
      pdj(106) =  &
          -k(601)*n(idx_H3j)
      pdj(108) =  &
          +k(853)*n(idx_H2Oj)
      pdj(110) =  &
          +k(856)*n(idx_HCOj)
      pdj(112) =  &
          -k(858)*n(idx_N2Hj)
      pdj(113) =  &
          -k(859)*n(idx_O2Hj)
    elseif(j==11) then
      pdj(1) =  &
          +k(1169)  &
          +k(280)
      pdj(2) =  &
          -k(937)*n(idx_CH)  &
          -k(936)*n(idx_CH)  &
          -k(935)*n(idx_CH)  &
          -k(934)*n(idx_CH)
      pdj(3) =  &
          +k(937)*n(idx_CH)  &
          +k(1107)*n(idx_SI)  &
          +k(1045)*n(idx_NH)  &
          +k(215)*n(idx_Oj)  &
          +2.d0*k(13)*n(idx_H)  &
          +k(954)*n(idx_CO)  &
          +2.d0*k(281)  &
          +k(1053)*n(idx_NO)  &
          +2.d0*k(7)*n(idx_H2)  &
          +k(697)*n(idx_HEj)  &
          +k(1024)*n(idx_N)  &
          +k(729)*n(idx_Nj)  &
          +k(377)*n(idx_Cj)  &
          +k(443)*n(idx_CH3j)  &
          +k(935)*n(idx_CH)  &
          +k(413)*n(idx_CHj)  &
          +k(950)*n(idx_CN)  &
          +k(874)*n(idx_C)  &
          +k(990)*n(idx_H)  &
          +k(778)*n(idx_NH2j)  &
          +k(893)*n(idx_CH2)  &
          +2.d0*k(1170)
      pdj(5) =  &
          +k(134)*n(idx_HCNj)
      pdj(6) =  &
          -k(7)*n(idx_H2)  &
          +k(114)*n(idx_H2j)  &
          +k(598)*n(idx_H3j)  &
          -k(965)*n(idx_H2)  &
          -k(964)*n(idx_H2)  &
          +k(7)*n(idx_H2)  &
          +k(890)*n(idx_CH2)
      pdj(7) =  &
          -k(874)*n(idx_C)
      pdj(8) =  &
          +k(964)*n(idx_H2)  &
          +k(89)*n(idx_Hj)  &
          +k(934)*n(idx_CH)  &
          +k(935)*n(idx_CH)  &
          +2.d0*k(891)*n(idx_CH2)  &
          -k(990)*n(idx_H)  &
          -k(13)*n(idx_H)  &
          +k(13)*n(idx_H)  &
          +k(526)*n(idx_H2j)
      pdj(9) =  &
          +k(122)*n(idx_H2Oj)  &
          +k(913)*n(idx_CH3)  &
          +k(892)*n(idx_CH2)
      pdj(10) =  &
          +k(990)*n(idx_H)  &
          +k(863)*n(idx_SIH2j)  &
          +k(779)*n(idx_NH2j)  &
          +2.d0*k(965)*n(idx_H2)  &
          +k(412)*n(idx_CHj)  &
          +k(766)*n(idx_NHj)  &
          +k(912)*n(idx_CH3)  &
          +k(894)*n(idx_CH2)  &
          +k(421)*n(idx_CH2j)  &
          +k(936)*n(idx_CH)  &
          +k(1002)*n(idx_HCO)  &
          +k(225)*n(idx_OHj)  &
          +k(1046)*n(idx_NH)
      pdj(11) =  &
          -k(7)*n(idx_H2)  &
          -k(215)*n(idx_Oj)  &
          -k(122)*n(idx_H2Oj)  &
          -k(89)*n(idx_Hj)  &
          -k(1169)  &
          -k(377)*n(idx_Cj)  &
          -k(443)*n(idx_CH3j)  &
          -k(1002)*n(idx_HCO)  &
          -k(730)*n(idx_Nj)  &
          -k(990)*n(idx_H)  &
          -k(52)*n(idx_CH4j)  &
          -k(281)  &
          -k(1053)*n(idx_NO)  &
          -k(729)*n(idx_Nj)  &
          -k(598)*n(idx_H3j)  &
          -k(914)*n(idx_CH3)  &
          -k(937)*n(idx_CH)  &
          -k(949)*n(idx_CN)  &
          -k(964)*n(idx_H2)  &
          -k(954)*n(idx_CO)  &
          -k(1293)  &
          -k(378)*n(idx_Cj)  &
          -k(922)*n(idx_CH4)  &
          -k(13)*n(idx_H)  &
          -k(147)*n(idx_HEj)  &
          -k(913)*n(idx_CH3)  &
          -k(1107)*n(idx_SI)  &
          -k(950)*n(idx_CN)  &
          -k(1170)  &
          -k(912)*n(idx_CH3)  &
          -k(934)*n(idx_CH)  &
          -k(412)*n(idx_CHj)  &
          -k(482)*n(idx_CNj)  &
          -k(779)*n(idx_NH2j)  &
          -k(965)*n(idx_H2)  &
          -k(74)*n(idx_COj)  &
          -k(778)*n(idx_NH2j)  &
          -k(280)  &
          -k(891)*n(idx_CH2)  &
          -k(134)*n(idx_HCNj)  &
          -k(526)*n(idx_H2j)  &
          -k(69)*n(idx_CNj)  &
          -k(114)*n(idx_H2j)  &
          -k(863)*n(idx_SIH2j)  &
          -k(1045)*n(idx_NH)  &
          -k(892)*n(idx_CH2)  &
          -k(1003)*n(idx_HCO)  &
          -k(169)*n(idx_Nj)  &
          -k(894)*n(idx_CH2)  &
          -k(767)*n(idx_NHj)  &
          -k(414)*n(idx_CHj)  &
          -k(1056)*n(idx_OCN)  &
          -k(874)*n(idx_C)  &
          -k(1055)*n(idx_OCN)  &
          -k(697)*n(idx_HEj)  &
          -k(413)*n(idx_CHj)  &
          -k(1024)*n(idx_N)  &
          -k(893)*n(idx_CH2)  &
          -k(936)*n(idx_CH)  &
          -k(935)*n(idx_CH)  &
          -k(225)*n(idx_OHj)  &
          -k(1046)*n(idx_NH)  &
          -k(174)*n(idx_N2j)  &
          -k(180)*n(idx_NHj)  &
          -k(550)*n(idx_H2COj)  &
          -k(890)*n(idx_CH2)  &
          -k(421)*n(idx_CH2j)  &
          -k(766)*n(idx_NHj)
      pdj(12) =  &
          +k(914)*n(idx_CH3)  &
          -k(894)*n(idx_CH2)  &
          -k(893)*n(idx_CH2)  &
          -k(891)*n(idx_CH2)  &
          -k(892)*n(idx_CH2)  &
          -k(890)*n(idx_CH2)
      pdj(13) =  &
          +k(893)*n(idx_CH2)  &
          +k(912)*n(idx_CH3)
      pdj(14) =  &
          -k(1002)*n(idx_HCO)  &
          +k(937)*n(idx_CH)  &
          +k(894)*n(idx_CH2)  &
          +k(913)*n(idx_CH3)  &
          +k(414)*n(idx_CHj)  &
          -k(1003)*n(idx_HCO)
      pdj(17) =  &
          +k(1055)*n(idx_OCN)  &
          +k(949)*n(idx_CN)  &
          +k(730)*n(idx_Nj)  &
          -k(1053)*n(idx_NO)  &
          +k(1024)*n(idx_N)  &
          +k(1046)*n(idx_NH)
      pdj(18) =  &
          -k(1107)*n(idx_SI)
      pdj(24) =  &
          -k(950)*n(idx_CN)  &
          -k(949)*n(idx_CN)  &
          +k(69)*n(idx_CNj)
      pdj(25) =  &
          +k(74)*n(idx_COj)  &
          +k(949)*n(idx_CN)  &
          +k(378)*n(idx_Cj)  &
          +k(874)*n(idx_C)  &
          +k(1003)*n(idx_HCO)  &
          +k(936)*n(idx_CH)  &
          +k(1056)*n(idx_OCN)  &
          +k(935)*n(idx_CH)  &
          +k(482)*n(idx_CNj)  &
          -k(954)*n(idx_CO)  &
          +k(892)*n(idx_CH2)
      pdj(26) =  &
          +k(174)*n(idx_N2j)
      pdj(28) =  &
          -k(913)*n(idx_CH3)  &
          -k(914)*n(idx_CH3)  &
          -k(912)*n(idx_CH3)  &
          +k(922)*n(idx_CH4)
      pdj(29) =  &
          +k(52)*n(idx_CH4j)  &
          -k(922)*n(idx_CH4)
      pdj(30) =  &
          -k(1024)*n(idx_N)  &
          +k(767)*n(idx_NHj)  &
          +k(169)*n(idx_Nj)
      pdj(31) =  &
          -k(1046)*n(idx_NH)  &
          -k(1045)*n(idx_NH)  &
          +k(180)*n(idx_NHj)
      pdj(34) =  &
          +k(1107)*n(idx_SI)
      pdj(35) =  &
          +k(697)*n(idx_HEj)  &
          +k(147)*n(idx_HEj)
      pdj(36) =  &
          +k(1045)*n(idx_NH)
      pdj(38) =  &
          +k(1055)*n(idx_OCN)  &
          +k(934)*n(idx_CH)  &
          +k(954)*n(idx_CO)  &
          +k(891)*n(idx_CH2)  &
          +k(1002)*n(idx_HCO)  &
          +k(890)*n(idx_CH2)
      pdj(42) =  &
          +k(1056)*n(idx_OCN)  &
          +k(1053)*n(idx_NO)
      pdj(43) =  &
          +k(914)*n(idx_CH3)  &
          +k(964)*n(idx_H2)  &
          +k(550)*n(idx_H2COj)  &
          +k(1003)*n(idx_HCO)  &
          +k(922)*n(idx_CH4)
      pdj(44) =  &
          +k(950)*n(idx_CN)  &
          -k(1056)*n(idx_OCN)  &
          -k(1055)*n(idx_OCN)
      pdj(61) =  &
          +k(1293)
      pdj(70) =  &
          +k(413)*n(idx_CHj)  &
          +k(550)*n(idx_H2COj)  &
          +k(421)*n(idx_CH2j)
      pdj(71) =  &
          -k(89)*n(idx_Hj)
      pdj(73) =  &
          -k(378)*n(idx_Cj)  &
          -k(377)*n(idx_Cj)
      pdj(74) =  &
          -k(421)*n(idx_CH2j)
      pdj(75) =  &
          -k(414)*n(idx_CHj)  &
          -k(412)*n(idx_CHj)  &
          -k(413)*n(idx_CHj)
      pdj(76) =  &
          -k(550)*n(idx_H2COj)
      pdj(79) =  &
          +k(766)*n(idx_NHj)  &
          +k(482)*n(idx_CNj)  &
          +k(729)*n(idx_Nj)
      pdj(84) =  &
          -k(863)*n(idx_SIH2j)
      pdj(86) =  &
          -k(69)*n(idx_CNj)  &
          -k(482)*n(idx_CNj)
      pdj(87) =  &
          +k(412)*n(idx_CHj)  &
          +k(377)*n(idx_Cj)  &
          -k(74)*n(idx_COj)
      pdj(88) =  &
          -k(174)*n(idx_N2j)
      pdj(89) =  &
          +k(74)*n(idx_COj)  &
          +k(114)*n(idx_H2j)  &
          +k(215)*n(idx_Oj)  &
          +k(180)*n(idx_NHj)  &
          +k(122)*n(idx_H2Oj)  &
          +k(225)*n(idx_OHj)  &
          +k(147)*n(idx_HEj)  &
          +k(169)*n(idx_Nj)  &
          +k(1169)  &
          +k(280)  &
          +k(69)*n(idx_CNj)  &
          +k(134)*n(idx_HCNj)  &
          +k(52)*n(idx_CH4j)  &
          +k(89)*n(idx_Hj)  &
          +k(174)*n(idx_N2j)
      pdj(90) =  &
          -k(122)*n(idx_H2Oj)
      pdj(91) =  &
          -k(779)*n(idx_NH2j)  &
          -k(778)*n(idx_NH2j)
      pdj(92) =  &
          -k(215)*n(idx_Oj)  &
          +k(414)*n(idx_CHj)  &
          +k(378)*n(idx_Cj)  &
          +k(730)*n(idx_Nj)  &
          +k(697)*n(idx_HEj)
      pdj(93) =  &
          -k(225)*n(idx_OHj)
      pdj(94) =  &
          -k(443)*n(idx_CH3j)
      pdj(95) =  &
          -k(52)*n(idx_CH4j)
      pdj(96) =  &
          -k(729)*n(idx_Nj)  &
          -k(169)*n(idx_Nj)  &
          -k(730)*n(idx_Nj)
      pdj(97) =  &
          -k(134)*n(idx_HCNj)
      pdj(98) =  &
          -k(180)*n(idx_NHj)  &
          -k(767)*n(idx_NHj)  &
          -k(766)*n(idx_NHj)
      pdj(102) =  &
          -k(526)*n(idx_H2j)  &
          -k(114)*n(idx_H2j)
      pdj(103) =  &
          -k(147)*n(idx_HEj)  &
          -k(697)*n(idx_HEj)
      pdj(104) =  &
          +k(779)*n(idx_NH2j)
      pdj(105) =  &
          +k(778)*n(idx_NH2j)
      pdj(106) =  &
          -k(598)*n(idx_H3j)
      pdj(107) =  &
          +k(443)*n(idx_CH3j)
      pdj(113) =  &
          +k(767)*n(idx_NHj)  &
          +k(598)*n(idx_H3j)  &
          +k(526)*n(idx_H2j)
      pdj(115) =  &
          +k(863)*n(idx_SIH2j)
    elseif(j==12) then
      pdj(1) =  &
          +k(1113)  &
          +k(243)
      pdj(2) =  &
          +k(881)*n(idx_CN)  &
          +k(244)  &
          +2.d0*k(864)*n(idx_C)  &
          +2.d0*k(879)*n(idx_CH2)  &
          +k(968)*n(idx_H)  &
          +k(1114)  &
          +k(898)*n(idx_O)  &
          +k(1008)*n(idx_N)  &
          +k(423)*n(idx_COj)  &
          +k(900)*n(idx_OH)
      pdj(3) =  &
          -k(896)*n(idx_O)  &
          +k(44)*n(idx_Oj)  &
          +k(438)*n(idx_OHj)  &
          -k(895)*n(idx_O)  &
          +k(893)*n(idx_O2)  &
          -k(897)*n(idx_O)  &
          +k(901)*n(idx_OH)  &
          -k(898)*n(idx_O)  &
          +k(436)*n(idx_O2j)
      pdj(4) =  &
          +k(429)*n(idx_HCNHj)  &
          +k(1007)*n(idx_N)
      pdj(5) =  &
          +k(888)*n(idx_NO)  &
          +k(885)*n(idx_N2)  &
          +k(428)*n(idx_HCNHj)  &
          +k(1006)*n(idx_N)  &
          +k(881)*n(idx_CN)
      pdj(6) =  &
          -k(957)*n(idx_H2)  &
          +k(492)*n(idx_Hj)  &
          +k(968)*n(idx_H)  &
          +k(101)*n(idx_H2j)  &
          +k(890)*n(idx_O2)  &
          +k(654)*n(idx_HEj)  &
          +k(578)*n(idx_H3j)  &
          +k(895)*n(idx_O)
      pdj(7) =  &
          +k(15)*n(idx_Cj)  &
          -k(864)*n(idx_C)
      pdj(8) =  &
          -k(968)*n(idx_H)  &
          +k(76)*n(idx_Hj)  &
          +k(957)*n(idx_H2)  &
          +k(244)  &
          +k(511)*n(idx_H2j)  &
          +k(897)*n(idx_O)  &
          +k(1114)  &
          +k(1006)*n(idx_N)  &
          +k(1007)*n(idx_N)  &
          +2.d0*k(891)*n(idx_O2)  &
          +k(655)*n(idx_HEj)  &
          +k(899)*n(idx_OH)  &
          +k(889)*n(idx_NO)  &
          +2.d0*k(896)*n(idx_O)
      pdj(9) =  &
          +k(41)*n(idx_H2Oj)  &
          +k(426)*n(idx_H3Oj)  &
          +k(900)*n(idx_OH)  &
          +k(892)*n(idx_O2)
      pdj(10) =  &
          +k(894)*n(idx_O2)  &
          -k(899)*n(idx_OH)  &
          +k(898)*n(idx_O)  &
          +k(425)*n(idx_H2Oj)  &
          -k(901)*n(idx_OH)  &
          +k(888)*n(idx_NO)  &
          +k(46)*n(idx_OHj)  &
          -k(900)*n(idx_OH)
      pdj(11) =  &
          -k(891)*n(idx_O2)  &
          -k(892)*n(idx_O2)  &
          -k(893)*n(idx_O2)  &
          +k(437)*n(idx_O2Hj)  &
          -k(890)*n(idx_O2)  &
          +k(45)*n(idx_O2j)  &
          -k(894)*n(idx_O2)
      pdj(12) =  &
          -k(430)*n(idx_HCOj)  &
          -k(968)*n(idx_H)  &
          -k(433)*n(idx_NHj)  &
          -k(40)*n(idx_H2COj)  &
          -k(897)*n(idx_O)  &
          -k(895)*n(idx_O)  &
          -k(899)*n(idx_OH)  &
          -k(886)*n(idx_NO2)  &
          -k(424)*n(idx_H2COj)  &
          -k(883)*n(idx_HCO)  &
          -k(901)*n(idx_OH)  &
          -k(864)*n(idx_C)  &
          -k(893)*n(idx_O2)  &
          -k(46)*n(idx_OHj)  &
          -k(76)*n(idx_Hj)  &
          -4.d0*k(879)*n(idx_CH2)  &
          -k(423)*n(idx_COj)  &
          -k(654)*n(idx_HEj)  &
          -k(15)*n(idx_Cj)  &
          -k(45)*n(idx_O2j)  &
          -k(492)*n(idx_Hj)  &
          -k(1006)*n(idx_N)  &
          -k(896)*n(idx_O)  &
          -k(511)*n(idx_H2j)  &
          -k(427)*n(idx_HCNj)  &
          -k(1113)  &
          -k(435)*n(idx_NH3j)  &
          -k(425)*n(idx_H2Oj)  &
          -k(428)*n(idx_HCNHj)  &
          -k(578)*n(idx_H3j)  &
          -k(892)*n(idx_O2)  &
          -k(38)*n(idx_CNj)  &
          -k(101)*n(idx_H2j)  &
          -k(439)*n(idx_SIOj)  &
          -k(39)*n(idx_COj)  &
          -k(887)*n(idx_NO)  &
          -k(243)  &
          -k(880)*n(idx_CH4)  &
          -k(438)*n(idx_OHj)  &
          -k(432)*n(idx_N2Hj)  &
          -k(888)*n(idx_NO)  &
          -k(885)*n(idx_N2)  &
          -k(882)*n(idx_H2CO)  &
          -k(244)  &
          -k(1257)  &
          -k(957)*n(idx_H2)  &
          -k(884)*n(idx_HNO)  &
          -k(655)*n(idx_HEj)  &
          -k(1008)*n(idx_N)  &
          -k(889)*n(idx_NO)  &
          -k(894)*n(idx_O2)  &
          -k(429)*n(idx_HCNHj)  &
          -k(43)*n(idx_NH2j)  &
          -k(436)*n(idx_O2j)  &
          -k(881)*n(idx_CN)  &
          -k(156)*n(idx_Nj)  &
          -k(891)*n(idx_O2)  &
          -k(898)*n(idx_O)  &
          -k(1114)  &
          -k(434)*n(idx_NH2j)  &
          -k(42)*n(idx_N2j)  &
          -k(41)*n(idx_H2Oj)  &
          -k(44)*n(idx_Oj)  &
          -k(431)*n(idx_HNOj)  &
          -k(890)*n(idx_O2)  &
          -k(900)*n(idx_OH)  &
          -k(426)*n(idx_H3Oj)  &
          -k(1007)*n(idx_N)  &
          -k(437)*n(idx_O2Hj)
      pdj(13) =  &
          +k(439)*n(idx_SIOj)  &
          +k(886)*n(idx_NO2)  &
          -k(882)*n(idx_H2CO)  &
          +k(893)*n(idx_O2)  &
          +k(40)*n(idx_H2COj)  &
          +k(887)*n(idx_NO)  &
          +k(899)*n(idx_OH)
      pdj(14) =  &
          +k(424)*n(idx_H2COj)  &
          +k(897)*n(idx_O)  &
          +k(882)*n(idx_H2CO)  &
          -k(883)*n(idx_HCO)  &
          +k(894)*n(idx_O2)
      pdj(17) =  &
          +k(886)*n(idx_NO2)  &
          -k(889)*n(idx_NO)  &
          +k(431)*n(idx_HNOj)  &
          +k(884)*n(idx_HNO)  &
          -k(888)*n(idx_NO)  &
          -k(887)*n(idx_NO)
      pdj(24) =  &
          +k(427)*n(idx_HCNj)  &
          +k(38)*n(idx_CNj)  &
          -k(881)*n(idx_CN)
      pdj(25) =  &
          +k(39)*n(idx_COj)  &
          +k(883)*n(idx_HCO)  &
          +k(895)*n(idx_O)  &
          +k(430)*n(idx_HCOj)  &
          +k(892)*n(idx_O2)  &
          +k(896)*n(idx_O)
      pdj(26) =  &
          +k(42)*n(idx_N2j)  &
          +k(432)*n(idx_N2Hj)  &
          -k(885)*n(idx_N2)
      pdj(27) =  &
          +k(435)*n(idx_NH3j)  &
          +k(43)*n(idx_NH2j)
      pdj(28) =  &
          +k(957)*n(idx_H2)  &
          +k(883)*n(idx_HCO)  &
          +2.d0*k(879)*n(idx_CH2)  &
          +k(884)*n(idx_HNO)  &
          +k(901)*n(idx_OH)  &
          +k(882)*n(idx_H2CO)  &
          +2.d0*k(880)*n(idx_CH4)
      pdj(29) =  &
          -k(880)*n(idx_CH4)
      pdj(30) =  &
          +k(433)*n(idx_NHj)  &
          -k(1007)*n(idx_N)  &
          -k(1006)*n(idx_N)  &
          +k(887)*n(idx_NO)  &
          +k(156)*n(idx_Nj)  &
          -k(1008)*n(idx_N)
      pdj(31) =  &
          +k(1008)*n(idx_N)  &
          +k(885)*n(idx_N2)  &
          +k(434)*n(idx_NH2j)
      pdj(35) =  &
          +k(654)*n(idx_HEj)  &
          +k(655)*n(idx_HEj)
      pdj(36) =  &
          -k(884)*n(idx_HNO)
      pdj(38) =  &
          +k(891)*n(idx_O2)  &
          +k(890)*n(idx_O2)
      pdj(41) =  &
          +k(889)*n(idx_NO)
      pdj(42) =  &
          -k(886)*n(idx_NO2)
      pdj(53) =  &
          +k(1257)
      pdj(70) =  &
          -k(430)*n(idx_HCOj)  &
          +k(423)*n(idx_COj)
      pdj(71) =  &
          -k(492)*n(idx_Hj)  &
          -k(76)*n(idx_Hj)
      pdj(73) =  &
          +k(654)*n(idx_HEj)  &
          -k(15)*n(idx_Cj)
      pdj(74) =  &
          +k(39)*n(idx_COj)  &
          +k(44)*n(idx_Oj)  &
          +k(76)*n(idx_Hj)  &
          +k(41)*n(idx_H2Oj)  &
          +k(40)*n(idx_H2COj)  &
          +k(101)*n(idx_H2j)  &
          +k(243)  &
          +k(1113)  &
          +k(156)*n(idx_Nj)  &
          +k(15)*n(idx_Cj)  &
          +k(45)*n(idx_O2j)  &
          +k(42)*n(idx_N2j)  &
          +k(43)*n(idx_NH2j)  &
          +k(46)*n(idx_OHj)  &
          +k(38)*n(idx_CNj)
      pdj(75) =  &
          +k(655)*n(idx_HEj)  &
          +k(492)*n(idx_Hj)
      pdj(76) =  &
          -k(424)*n(idx_H2COj)  &
          -k(40)*n(idx_H2COj)  &
          +k(436)*n(idx_O2j)
      pdj(78) =  &
          -k(435)*n(idx_NH3j)
      pdj(80) =  &
          +k(439)*n(idx_SIOj)
      pdj(86) =  &
          -k(38)*n(idx_CNj)
      pdj(87) =  &
          -k(423)*n(idx_COj)  &
          -k(39)*n(idx_COj)
      pdj(88) =  &
          -k(42)*n(idx_N2j)
      pdj(89) =  &
          -k(436)*n(idx_O2j)  &
          -k(45)*n(idx_O2j)
      pdj(90) =  &
          -k(425)*n(idx_H2Oj)  &
          -k(41)*n(idx_H2Oj)
      pdj(91) =  &
          -k(43)*n(idx_NH2j)  &
          -k(434)*n(idx_NH2j)
      pdj(92) =  &
          -k(44)*n(idx_Oj)
      pdj(93) =  &
          -k(46)*n(idx_OHj)  &
          -k(438)*n(idx_OHj)
      pdj(94) =  &
          +k(433)*n(idx_NHj)  &
          +k(438)*n(idx_OHj)  &
          +k(429)*n(idx_HCNHj)  &
          +k(431)*n(idx_HNOj)  &
          +k(432)*n(idx_N2Hj)  &
          +k(428)*n(idx_HCNHj)  &
          +k(424)*n(idx_H2COj)  &
          +k(511)*n(idx_H2j)  &
          +k(430)*n(idx_HCOj)  &
          +k(425)*n(idx_H2Oj)  &
          +k(437)*n(idx_O2Hj)  &
          +k(427)*n(idx_HCNj)  &
          +k(578)*n(idx_H3j)  &
          +k(426)*n(idx_H3Oj)  &
          +k(434)*n(idx_NH2j)  &
          +k(435)*n(idx_NH3j)
      pdj(96) =  &
          -k(156)*n(idx_Nj)
      pdj(97) =  &
          -k(427)*n(idx_HCNj)
      pdj(98) =  &
          -k(433)*n(idx_NHj)
      pdj(101) =  &
          -k(439)*n(idx_SIOj)
      pdj(102) =  &
          -k(101)*n(idx_H2j)  &
          -k(511)*n(idx_H2j)
      pdj(103) =  &
          -k(654)*n(idx_HEj)  &
          -k(655)*n(idx_HEj)
      pdj(104) =  &
          -k(431)*n(idx_HNOj)
      pdj(106) =  &
          -k(578)*n(idx_H3j)
      pdj(108) =  &
          -k(426)*n(idx_H3Oj)
      pdj(109) =  &
          -k(428)*n(idx_HCNHj)  &
          -k(429)*n(idx_HCNHj)
      pdj(112) =  &
          -k(432)*n(idx_N2Hj)
      pdj(113) =  &
          -k(437)*n(idx_O2Hj)
    elseif(j==13) then
      pdj(1) =  &
          +k(1139)  &
          +k(1140)
      pdj(2) =  &
          -k(925)*n(idx_CH)  &
          +k(370)*n(idx_Cj)
      pdj(3) =  &
          +k(840)*n(idx_OHj)  &
          -k(1062)*n(idx_O)  &
          +k(210)*n(idx_Oj)  &
          +k(673)*n(idx_HEj)
      pdj(4) =  &
          +k(635)*n(idx_HCNHj)
      pdj(5) =  &
          +k(634)*n(idx_HCNHj)  &
          +k(943)*n(idx_CN)  &
          +k(480)*n(idx_CNj)
      pdj(6) =  &
          +k(498)*n(idx_Hj)  &
          +k(1137)  &
          +k(256)  &
          +k(586)*n(idx_H3j)  &
          +k(499)*n(idx_Hj)  &
          +k(106)*n(idx_H2j)  &
          +k(975)*n(idx_H)  &
          +k(671)*n(idx_HEj)  &
          +k(518)*n(idx_H2j)
      pdj(7) =  &
          +k(401)*n(idx_CHj)  &
          +k(17)*n(idx_Cj)
      pdj(8) =  &
          +k(80)*n(idx_Hj)  &
          +k(498)*n(idx_Hj)  &
          +k(672)*n(idx_HEj)  &
          +2.d0*k(1138)  &
          +k(731)*n(idx_N2j)  &
          +k(552)*n(idx_O2j)  &
          +k(518)*n(idx_H2j)  &
          +k(1140)  &
          -k(975)*n(idx_H)
      pdj(9) =  &
          +k(1094)*n(idx_OH)  &
          +k(118)*n(idx_H2Oj)  &
          +k(608)*n(idx_H3Oj)
      pdj(10) =  &
          -k(1094)*n(idx_OH)  &
          +k(220)*n(idx_OHj)  &
          +k(555)*n(idx_H2Oj)  &
          +k(814)*n(idx_Oj)  &
          +k(1062)*n(idx_O)
      pdj(11) =  &
          +k(117)*n(idx_O2j)  &
          +k(553)*n(idx_O2Hj)  &
          +k(552)*n(idx_O2j)
      pdj(12) =  &
          +k(925)*n(idx_CH)  &
          +k(402)*n(idx_CHj)  &
          +k(723)*n(idx_Nj)  &
          -k(882)*n(idx_CH2)
      pdj(13) =  &
          -k(636)*n(idx_HCOj)  &
          -k(17)*n(idx_Cj)  &
          -k(143)*n(idx_HEj)  &
          -k(401)*n(idx_CHj)  &
          -k(731)*n(idx_N2j)  &
          -k(975)*n(idx_H)  &
          -k(753)*n(idx_NHj)  &
          -k(814)*n(idx_Oj)  &
          -k(1253)  &
          -k(771)*n(idx_NH2j)  &
          -k(723)*n(idx_Nj)  &
          -k(925)*n(idx_CH)  &
          -k(608)*n(idx_H3Oj)  &
          -k(1139)  &
          -k(369)*n(idx_Cj)  &
          -k(904)*n(idx_CH3)  &
          -k(220)*n(idx_OHj)  &
          -k(943)*n(idx_CN)  &
          -k(71)*n(idx_COj)  &
          -k(722)*n(idx_Nj)  &
          -k(402)*n(idx_CHj)  &
          -k(118)*n(idx_H2Oj)  &
          -k(498)*n(idx_Hj)  &
          -k(160)*n(idx_Nj)  &
          -k(1094)*n(idx_OH)  &
          -k(840)*n(idx_OHj)  &
          -k(586)*n(idx_H3j)  &
          -k(400)*n(idx_CHj)  &
          -k(754)*n(idx_NHj)  &
          -k(480)*n(idx_CNj)  &
          -k(553)*n(idx_O2Hj)  &
          -k(671)*n(idx_HEj)  &
          -k(1140)  &
          -k(673)*n(idx_HEj)  &
          -k(80)*n(idx_Hj)  &
          -k(549)*n(idx_H2COj)  &
          -k(441)*n(idx_CH3j)  &
          -k(634)*n(idx_HCNHj)  &
          -k(176)*n(idx_NHj)  &
          -k(623)*n(idx_HCNj)  &
          -k(65)*n(idx_CNj)  &
          -k(672)*n(idx_HEj)  &
          -k(635)*n(idx_HCNHj)  &
          -k(418)*n(idx_CH2j)  &
          -k(1062)*n(idx_O)  &
          -k(210)*n(idx_Oj)  &
          -k(117)*n(idx_O2j)  &
          -k(551)*n(idx_HNOj)  &
          -k(1138)  &
          -k(171)*n(idx_N2j)  &
          -k(555)*n(idx_H2Oj)  &
          -k(106)*n(idx_H2j)  &
          -k(736)*n(idx_N2Hj)  &
          -k(50)*n(idx_CH4j)  &
          -k(370)*n(idx_Cj)  &
          -k(485)*n(idx_COj)  &
          -k(450)*n(idx_CH4j)  &
          -k(1137)  &
          -k(256)  &
          -k(499)*n(idx_Hj)  &
          -k(552)*n(idx_O2j)  &
          -k(518)*n(idx_H2j)  &
          -k(770)*n(idx_NH2j)  &
          -k(882)*n(idx_CH2)
      pdj(14) =  &
          +k(1094)*n(idx_OH)  &
          +k(943)*n(idx_CN)  &
          +k(1062)*n(idx_O)  &
          +k(549)*n(idx_H2COj)  &
          +k(485)*n(idx_COj)  &
          +k(975)*n(idx_H)  &
          +k(771)*n(idx_NH2j)  &
          +k(882)*n(idx_CH2)  &
          +k(925)*n(idx_CH)  &
          +k(904)*n(idx_CH3)
      pdj(17) =  &
          +k(551)*n(idx_HNOj)
      pdj(24) =  &
          +k(65)*n(idx_CNj)  &
          +k(623)*n(idx_HCNj)  &
          -k(943)*n(idx_CN)
      pdj(25) =  &
          +k(1137)  &
          +k(256)  &
          +k(636)*n(idx_HCOj)  &
          +k(400)*n(idx_CHj)  &
          +k(71)*n(idx_COj)  &
          +k(369)*n(idx_Cj)  &
          +k(1138)
      pdj(26) =  &
          +k(736)*n(idx_N2Hj)  &
          +k(731)*n(idx_N2j)  &
          +k(171)*n(idx_N2j)
      pdj(27) =  &
          +k(754)*n(idx_NHj)
      pdj(28) =  &
          -k(904)*n(idx_CH3)  &
          +k(882)*n(idx_CH2)  &
          +k(450)*n(idx_CH4j)  &
          +k(418)*n(idx_CH2j)
      pdj(29) =  &
          +k(441)*n(idx_CH3j)  &
          +k(50)*n(idx_CH4j)  &
          +k(904)*n(idx_CH3)
      pdj(30) =  &
          +k(753)*n(idx_NHj)  &
          +k(160)*n(idx_Nj)
      pdj(31) =  &
          +k(770)*n(idx_NH2j)  &
          +k(722)*n(idx_Nj)  &
          +k(176)*n(idx_NHj)
      pdj(35) =  &
          +k(143)*n(idx_HEj)  &
          +k(672)*n(idx_HEj)  &
          +k(673)*n(idx_HEj)  &
          +k(671)*n(idx_HEj)
      pdj(47) =  &
          +k(1253)
      pdj(70) =  &
          -k(636)*n(idx_HCOj)  &
          +k(731)*n(idx_N2j)  &
          +k(370)*n(idx_Cj)  &
          +k(485)*n(idx_COj)  &
          +k(722)*n(idx_Nj)  &
          +k(402)*n(idx_CHj)  &
          +k(499)*n(idx_Hj)  &
          +k(754)*n(idx_NHj)  &
          +k(480)*n(idx_CNj)  &
          +k(814)*n(idx_Oj)  &
          +k(418)*n(idx_CH2j)  &
          +k(518)*n(idx_H2j)  &
          +k(672)*n(idx_HEj)  &
          +k(552)*n(idx_O2j)  &
          +k(441)*n(idx_CH3j)  &
          +k(1140)
      pdj(71) =  &
          -k(80)*n(idx_Hj)  &
          -k(499)*n(idx_Hj)  &
          -k(498)*n(idx_Hj)
      pdj(73) =  &
          -k(369)*n(idx_Cj)  &
          -k(17)*n(idx_Cj)  &
          -k(370)*n(idx_Cj)
      pdj(74) =  &
          -k(418)*n(idx_CH2j)  &
          +k(673)*n(idx_HEj)  &
          +k(369)*n(idx_Cj)
      pdj(75) =  &
          -k(402)*n(idx_CHj)  &
          -k(401)*n(idx_CHj)  &
          -k(400)*n(idx_CHj)
      pdj(76) =  &
          +k(80)*n(idx_Hj)  &
          +k(220)*n(idx_OHj)  &
          +k(160)*n(idx_Nj)  &
          -k(549)*n(idx_H2COj)  &
          +k(118)*n(idx_H2Oj)  &
          +k(143)*n(idx_HEj)  &
          +k(176)*n(idx_NHj)  &
          +k(17)*n(idx_Cj)  &
          +k(106)*n(idx_H2j)  &
          +k(117)*n(idx_O2j)  &
          +k(71)*n(idx_COj)  &
          +k(1139)  &
          +k(65)*n(idx_CNj)  &
          +k(50)*n(idx_CH4j)  &
          +k(171)*n(idx_N2j)  &
          +k(210)*n(idx_Oj)
      pdj(78) =  &
          +k(771)*n(idx_NH2j)
      pdj(79) =  &
          +k(723)*n(idx_Nj)
      pdj(86) =  &
          -k(480)*n(idx_CNj)  &
          -k(65)*n(idx_CNj)
      pdj(87) =  &
          +k(498)*n(idx_Hj)  &
          -k(485)*n(idx_COj)  &
          +k(671)*n(idx_HEj)  &
          -k(71)*n(idx_COj)
      pdj(88) =  &
          -k(731)*n(idx_N2j)  &
          -k(171)*n(idx_N2j)
      pdj(89) =  &
          -k(552)*n(idx_O2j)  &
          -k(117)*n(idx_O2j)
      pdj(90) =  &
          -k(555)*n(idx_H2Oj)  &
          -k(118)*n(idx_H2Oj)
      pdj(91) =  &
          -k(771)*n(idx_NH2j)  &
          -k(770)*n(idx_NH2j)
      pdj(92) =  &
          -k(210)*n(idx_Oj)  &
          -k(814)*n(idx_Oj)
      pdj(93) =  &
          -k(220)*n(idx_OHj)  &
          -k(840)*n(idx_OHj)
      pdj(94) =  &
          +k(400)*n(idx_CHj)  &
          -k(441)*n(idx_CH3j)
      pdj(95) =  &
          -k(450)*n(idx_CH4j)  &
          -k(50)*n(idx_CH4j)
      pdj(96) =  &
          -k(160)*n(idx_Nj)  &
          -k(722)*n(idx_Nj)  &
          -k(723)*n(idx_Nj)
      pdj(97) =  &
          -k(623)*n(idx_HCNj)
      pdj(98) =  &
          -k(754)*n(idx_NHj)  &
          -k(753)*n(idx_NHj)  &
          -k(176)*n(idx_NHj)
      pdj(102) =  &
          -k(518)*n(idx_H2j)  &
          -k(106)*n(idx_H2j)
      pdj(103) =  &
          -k(673)*n(idx_HEj)  &
          -k(143)*n(idx_HEj)  &
          -k(672)*n(idx_HEj)  &
          -k(671)*n(idx_HEj)
      pdj(104) =  &
          -k(551)*n(idx_HNOj)
      pdj(106) =  &
          -k(586)*n(idx_H3j)
      pdj(107) =  &
          +k(551)*n(idx_HNOj)  &
          +k(586)*n(idx_H3j)  &
          +k(636)*n(idx_HCOj)  &
          +k(553)*n(idx_O2Hj)  &
          +k(623)*n(idx_HCNj)  &
          +k(753)*n(idx_NHj)  &
          +k(549)*n(idx_H2COj)  &
          +k(555)*n(idx_H2Oj)  &
          +k(401)*n(idx_CHj)  &
          +k(840)*n(idx_OHj)  &
          +k(634)*n(idx_HCNHj)  &
          +k(770)*n(idx_NH2j)  &
          +k(635)*n(idx_HCNHj)  &
          +k(736)*n(idx_N2Hj)  &
          +k(450)*n(idx_CH4j)  &
          +k(608)*n(idx_H3Oj)
      pdj(108) =  &
          -k(608)*n(idx_H3Oj)
      pdj(109) =  &
          -k(634)*n(idx_HCNHj)  &
          -k(635)*n(idx_HCNHj)
      pdj(112) =  &
          -k(736)*n(idx_N2Hj)
      pdj(113) =  &
          -k(553)*n(idx_O2Hj)
    elseif(j==14) then
      pdj(1) =  &
          +k(1151)  &
          +k(262)
      pdj(2) =  &
          +k(865)*n(idx_C)  &
          +k(32)*n(idx_CHj)  &
          -k(926)*n(idx_CH)
      pdj(3) =  &
          +k(683)*n(idx_HEj)  &
          -k(1068)*n(idx_O)  &
          +k(1016)*n(idx_N)  &
          +k(212)*n(idx_Oj)  &
          +k(979)*n(idx_H)  &
          +k(844)*n(idx_OHj)  &
          -k(1067)*n(idx_O)
      pdj(5) =  &
          +k(1016)*n(idx_N)  &
          +k(944)*n(idx_CN)
      pdj(6) =  &
          +k(501)*n(idx_Hj)  &
          +2.d0*k(998)*n(idx_HCO)  &
          +k(978)*n(idx_H)  &
          +k(109)*n(idx_H2j)  &
          +k(589)*n(idx_H3j)
      pdj(7) =  &
          -k(865)*n(idx_C)  &
          +k(18)*n(idx_Cj)
      pdj(8) =  &
          +k(1150)  &
          +k(261)  &
          +k(1017)*n(idx_N)  &
          -k(979)*n(idx_H)  &
          +k(1067)*n(idx_O)  &
          +k(83)*n(idx_Hj)  &
          +k(681)*n(idx_HEj)  &
          -k(978)*n(idx_H)
      pdj(9) =  &
          +k(1097)*n(idx_OH)  &
          +k(119)*n(idx_H2Oj)
      pdj(10) =  &
          +k(1068)*n(idx_O)  &
          +k(559)*n(idx_H2Oj)  &
          -k(1097)*n(idx_OH)  &
          +k(1002)*n(idx_O2)  &
          +k(222)*n(idx_OHj)
      pdj(11) =  &
          +k(138)*n(idx_O2j)  &
          -k(1002)*n(idx_O2)  &
          -k(1003)*n(idx_O2)  &
          +k(646)*n(idx_O2Hj)  &
          +k(1004)*n(idx_O2H)
      pdj(12) =  &
          +k(979)*n(idx_H)  &
          -k(883)*n(idx_CH2)  &
          +k(926)*n(idx_CH)
      pdj(13) =  &
          +k(1000)*n(idx_HNO)  &
          +k(1004)*n(idx_O2H)  &
          +k(137)*n(idx_H2COj)  &
          +2.d0*k(999)*n(idx_HCO)
      pdj(14) =  &
          -k(642)*n(idx_H2COj)  &
          -k(1001)*n(idx_NO)  &
          -k(262)  &
          -k(138)*n(idx_O2j)  &
          -k(865)*n(idx_C)  &
          -k(119)*n(idx_H2Oj)  &
          -k(724)*n(idx_Nj)  &
          -k(18)*n(idx_Cj)  &
          -k(1004)*n(idx_O2H)  &
          -k(1151)  &
          -k(139)*n(idx_SIOj)  &
          -4.d0*k(998)*n(idx_HCO)  &
          -k(637)*n(idx_HCOj)  &
          -k(775)*n(idx_NH2j)  &
          -k(212)*n(idx_Oj)  &
          -k(978)*n(idx_H)  &
          -k(682)*n(idx_HEj)  &
          -k(47)*n(idx_CH3j)  &
          -k(1003)*n(idx_O2)  &
          -k(407)*n(idx_CHj)  &
          -k(1016)*n(idx_N)  &
          -k(32)*n(idx_CHj)  &
          -k(646)*n(idx_O2Hj)  &
          -k(163)*n(idx_Nj)  &
          -k(520)*n(idx_H2j)  &
          -k(979)*n(idx_H)  &
          -k(190)*n(idx_NH3j)  &
          -k(732)*n(idx_N2j)  &
          -k(626)*n(idx_HCNj)  &
          -k(67)*n(idx_CNj)  &
          -k(926)*n(idx_CH)  &
          -k(420)*n(idx_CH2j)  &
          -4.d0*k(999)*n(idx_HCO)  &
          -k(644)*n(idx_N2Hj)  &
          -k(843)*n(idx_OHj)  &
          -k(181)*n(idx_NH2j)  &
          -k(1097)*n(idx_OH)  &
          -k(72)*n(idx_COj)  &
          -k(1262)  &
          -k(442)*n(idx_CH3j)  &
          -k(683)*n(idx_HEj)  &
          -k(1150)  &
          -k(760)*n(idx_NHj)  &
          -k(109)*n(idx_H2j)  &
          -k(261)  &
          -k(1017)*n(idx_N)  &
          -k(83)*n(idx_Hj)  &
          -k(625)*n(idx_HCNj)  &
          -k(481)*n(idx_CNj)  &
          -k(906)*n(idx_CH3)  &
          -k(645)*n(idx_O2j)  &
          -k(643)*n(idx_HNOj)  &
          -k(817)*n(idx_Oj)  &
          -k(1067)*n(idx_O)  &
          -k(844)*n(idx_OHj)  &
          -k(559)*n(idx_H2Oj)  &
          -k(137)*n(idx_H2COj)  &
          -k(1000)*n(idx_HNO)  &
          -k(883)*n(idx_CH2)  &
          -k(1068)*n(idx_O)  &
          -k(681)*n(idx_HEj)  &
          -k(222)*n(idx_OHj)  &
          -k(502)*n(idx_Hj)  &
          -k(944)*n(idx_CN)  &
          -k(1015)*n(idx_N)  &
          -k(1002)*n(idx_O2)  &
          -k(172)*n(idx_N2j)  &
          -k(589)*n(idx_H3j)  &
          -k(373)*n(idx_Cj)  &
          -k(558)*n(idx_H2Oj)  &
          -k(501)*n(idx_Hj)
      pdj(16) =  &
          +k(190)*n(idx_NH3j)
      pdj(17) =  &
          +k(1000)*n(idx_HNO)  &
          -k(1001)*n(idx_NO)  &
          +k(643)*n(idx_HNOj)
      pdj(24) =  &
          +k(625)*n(idx_HCNj)  &
          +k(67)*n(idx_CNj)  &
          -k(944)*n(idx_CN)
      pdj(25) =  &
          +k(682)*n(idx_HEj)  &
          +k(1015)*n(idx_N)  &
          +k(502)*n(idx_Hj)  &
          +k(944)*n(idx_CN)  &
          +k(978)*n(idx_H)  &
          +k(420)*n(idx_CH2j)  &
          +k(520)*n(idx_H2j)  &
          +k(724)*n(idx_Nj)  &
          +k(1150)  &
          +k(373)*n(idx_Cj)  &
          +k(261)  &
          +k(926)*n(idx_CH)  &
          +k(642)*n(idx_H2COj)  &
          +k(732)*n(idx_N2j)  &
          +k(1068)*n(idx_O)  &
          +4.d0*k(998)*n(idx_HCO)  &
          +k(843)*n(idx_OHj)  &
          +k(637)*n(idx_HCOj)  &
          +k(558)*n(idx_H2Oj)  &
          +k(481)*n(idx_CNj)  &
          +k(72)*n(idx_COj)  &
          +k(407)*n(idx_CHj)  &
          +2.d0*k(999)*n(idx_HCO)  &
          +k(1001)*n(idx_NO)  &
          +k(626)*n(idx_HCNj)  &
          +k(1003)*n(idx_O2)  &
          +k(1097)*n(idx_OH)  &
          +k(645)*n(idx_O2j)  &
          +k(906)*n(idx_CH3)  &
          +k(865)*n(idx_C)  &
          +k(883)*n(idx_CH2)  &
          +k(442)*n(idx_CH3j)  &
          +k(817)*n(idx_Oj)
      pdj(26) =  &
          +k(172)*n(idx_N2j)  &
          +k(644)*n(idx_N2Hj)
      pdj(27) =  &
          +k(181)*n(idx_NH2j)
      pdj(28) =  &
          -k(906)*n(idx_CH3)  &
          +k(883)*n(idx_CH2)  &
          +k(47)*n(idx_CH3j)
      pdj(29) =  &
          +k(906)*n(idx_CH3)
      pdj(30) =  &
          -k(1015)*n(idx_N)  &
          -k(1016)*n(idx_N)  &
          -k(1017)*n(idx_N)  &
          +k(760)*n(idx_NHj)  &
          +k(163)*n(idx_Nj)
      pdj(31) =  &
          +k(775)*n(idx_NH2j)  &
          +k(1015)*n(idx_N)
      pdj(34) =  &
          +k(139)*n(idx_SIOj)
      pdj(35) =  &
          +k(683)*n(idx_HEj)  &
          +k(681)*n(idx_HEj)
      pdj(36) =  &
          -k(1000)*n(idx_HNO)  &
          +k(1001)*n(idx_NO)
      pdj(38) =  &
          +k(1067)*n(idx_O)  &
          +k(1002)*n(idx_O2)
      pdj(43) =  &
          +k(1003)*n(idx_O2)  &
          -k(1004)*n(idx_O2H)
      pdj(44) =  &
          +k(1017)*n(idx_N)
      pdj(47) =  &
          +k(1262)
      pdj(70) =  &
          +k(137)*n(idx_H2COj)  &
          +k(138)*n(idx_O2j)  &
          +k(190)*n(idx_NH3j)  &
          +k(109)*n(idx_H2j)  &
          +k(18)*n(idx_Cj)  &
          +k(47)*n(idx_CH3j)  &
          +k(181)*n(idx_NH2j)  &
          +k(32)*n(idx_CHj)  &
          +k(163)*n(idx_Nj)  &
          +k(172)*n(idx_N2j)  &
          +k(262)  &
          -k(637)*n(idx_HCOj)  &
          +k(1151)  &
          +k(72)*n(idx_COj)  &
          +k(139)*n(idx_SIOj)  &
          +k(83)*n(idx_Hj)  &
          +k(119)*n(idx_H2Oj)  &
          +k(212)*n(idx_Oj)  &
          +k(222)*n(idx_OHj)  &
          +k(67)*n(idx_CNj)
      pdj(71) =  &
          -k(83)*n(idx_Hj)  &
          -k(502)*n(idx_Hj)  &
          -k(501)*n(idx_Hj)
      pdj(73) =  &
          -k(373)*n(idx_Cj)  &
          -k(18)*n(idx_Cj)
      pdj(74) =  &
          -k(420)*n(idx_CH2j)  &
          +k(407)*n(idx_CHj)
      pdj(75) =  &
          +k(683)*n(idx_HEj)  &
          -k(407)*n(idx_CHj)  &
          +k(373)*n(idx_Cj)  &
          -k(32)*n(idx_CHj)
      pdj(76) =  &
          -k(642)*n(idx_H2COj)  &
          +k(775)*n(idx_NH2j)  &
          -k(137)*n(idx_H2COj)  &
          +k(625)*n(idx_HCNj)  &
          +k(643)*n(idx_HNOj)  &
          +k(646)*n(idx_O2Hj)  &
          +k(760)*n(idx_NHj)  &
          +k(844)*n(idx_OHj)  &
          +k(637)*n(idx_HCOj)  &
          +k(589)*n(idx_H3j)  &
          +k(644)*n(idx_N2Hj)  &
          +k(559)*n(idx_H2Oj)
      pdj(78) =  &
          -k(190)*n(idx_NH3j)
      pdj(86) =  &
          -k(481)*n(idx_CNj)  &
          -k(67)*n(idx_CNj)
      pdj(87) =  &
          -k(72)*n(idx_COj)  &
          +k(501)*n(idx_Hj)  &
          +k(681)*n(idx_HEj)
      pdj(88) =  &
          -k(172)*n(idx_N2j)  &
          -k(732)*n(idx_N2j)
      pdj(89) =  &
          -k(138)*n(idx_O2j)  &
          -k(645)*n(idx_O2j)
      pdj(90) =  &
          -k(119)*n(idx_H2Oj)  &
          +k(843)*n(idx_OHj)  &
          -k(558)*n(idx_H2Oj)  &
          -k(559)*n(idx_H2Oj)
      pdj(91) =  &
          -k(775)*n(idx_NH2j)  &
          -k(181)*n(idx_NH2j)
      pdj(92) =  &
          -k(817)*n(idx_Oj)  &
          -k(212)*n(idx_Oj)
      pdj(93) =  &
          -k(843)*n(idx_OHj)  &
          -k(222)*n(idx_OHj)  &
          +k(817)*n(idx_Oj)  &
          -k(844)*n(idx_OHj)
      pdj(94) =  &
          +k(420)*n(idx_CH2j)  &
          -k(47)*n(idx_CH3j)  &
          -k(442)*n(idx_CH3j)
      pdj(95) =  &
          +k(442)*n(idx_CH3j)
      pdj(96) =  &
          -k(163)*n(idx_Nj)  &
          -k(724)*n(idx_Nj)
      pdj(97) =  &
          -k(626)*n(idx_HCNj)  &
          +k(481)*n(idx_CNj)  &
          -k(625)*n(idx_HCNj)
      pdj(98) =  &
          +k(724)*n(idx_Nj)  &
          -k(760)*n(idx_NHj)
      pdj(101) =  &
          -k(139)*n(idx_SIOj)
      pdj(102) =  &
          -k(109)*n(idx_H2j)  &
          +k(502)*n(idx_Hj)  &
          -k(520)*n(idx_H2j)
      pdj(103) =  &
          -k(681)*n(idx_HEj)  &
          -k(683)*n(idx_HEj)  &
          -k(682)*n(idx_HEj)
      pdj(104) =  &
          -k(643)*n(idx_HNOj)
      pdj(106) =  &
          +k(520)*n(idx_H2j)  &
          -k(589)*n(idx_H3j)
      pdj(107) =  &
          +k(642)*n(idx_H2COj)
      pdj(108) =  &
          +k(558)*n(idx_H2Oj)
      pdj(109) =  &
          +k(626)*n(idx_HCNj)
      pdj(111) =  &
          +k(682)*n(idx_HEj)
      pdj(112) =  &
          -k(644)*n(idx_N2Hj)  &
          +k(732)*n(idx_N2j)
      pdj(113) =  &
          -k(646)*n(idx_O2Hj)  &
          +k(645)*n(idx_O2j)
    elseif(j==15) then
      pdj(1) =  &
          +k(267)  &
          +k(1155)
      pdj(2) =  &
          +k(33)*n(idx_CHj)
      pdj(6) =  &
          +k(592)*n(idx_H3j)
      pdj(7) =  &
          +k(19)*n(idx_Cj)
      pdj(8) =  &
          +k(84)*n(idx_Hj)  &
          +k(592)*n(idx_H3j)
      pdj(9) =  &
          +k(120)*n(idx_H2Oj)
      pdj(11) =  &
          +k(153)*n(idx_O2j)
      pdj(13) =  &
          +k(149)*n(idx_H2COj)
      pdj(14) =  &
          +k(150)*n(idx_HCOj)
      pdj(15) =  &
          -k(151)*n(idx_N2j)  &
          -k(153)*n(idx_O2j)  &
          -k(592)*n(idx_H3j)  &
          -k(1301)  &
          -k(19)*n(idx_Cj)  &
          -k(152)*n(idx_NOj)  &
          -k(149)*n(idx_H2COj)  &
          -k(191)*n(idx_NH3j)  &
          -k(155)*n(idx_SIOj)  &
          -k(150)*n(idx_HCOj)  &
          -k(154)*n(idx_SIj)  &
          -k(120)*n(idx_H2Oj)  &
          -k(84)*n(idx_Hj)  &
          -k(1155)  &
          -k(164)*n(idx_Nj)  &
          -k(267)  &
          -k(48)*n(idx_CH3j)  &
          -k(33)*n(idx_CHj)
      pdj(16) =  &
          +k(191)*n(idx_NH3j)
      pdj(17) =  &
          +k(152)*n(idx_NOj)
      pdj(18) =  &
          +k(154)*n(idx_SIj)
      pdj(26) =  &
          +k(151)*n(idx_N2j)
      pdj(28) =  &
          +k(48)*n(idx_CH3j)
      pdj(30) =  &
          +k(164)*n(idx_Nj)
      pdj(34) =  &
          +k(155)*n(idx_SIOj)
      pdj(66) =  &
          +k(1301)
      pdj(70) =  &
          -k(150)*n(idx_HCOj)
      pdj(71) =  &
          -k(84)*n(idx_Hj)
      pdj(73) =  &
          -k(19)*n(idx_Cj)
      pdj(75) =  &
          -k(33)*n(idx_CHj)
      pdj(76) =  &
          -k(149)*n(idx_H2COj)
      pdj(77) =  &
          +k(267)  &
          +k(151)*n(idx_N2j)  &
          +k(164)*n(idx_Nj)  &
          +k(120)*n(idx_H2Oj)  &
          +k(150)*n(idx_HCOj)  &
          +k(84)*n(idx_Hj)  &
          +k(149)*n(idx_H2COj)  &
          +k(155)*n(idx_SIOj)  &
          +k(154)*n(idx_SIj)  &
          +k(191)*n(idx_NH3j)  &
          +k(1155)  &
          +k(19)*n(idx_Cj)  &
          +k(153)*n(idx_O2j)  &
          +k(592)*n(idx_H3j)  &
          +k(33)*n(idx_CHj)  &
          +k(152)*n(idx_NOj)  &
          +k(48)*n(idx_CH3j)
      pdj(78) =  &
          -k(191)*n(idx_NH3j)
      pdj(79) =  &
          -k(152)*n(idx_NOj)
      pdj(80) =  &
          -k(154)*n(idx_SIj)
      pdj(88) =  &
          -k(151)*n(idx_N2j)
      pdj(89) =  &
          -k(153)*n(idx_O2j)
      pdj(90) =  &
          -k(120)*n(idx_H2Oj)
      pdj(94) =  &
          -k(48)*n(idx_CH3j)
      pdj(96) =  &
          -k(164)*n(idx_Nj)
      pdj(101) =  &
          -k(155)*n(idx_SIOj)
      pdj(106) =  &
          -k(592)*n(idx_H3j)
    elseif(j==16) then
      pdj(1) =  &
          +k(1161)  &
          +k(273)
      pdj(2) =  &
          +k(34)*n(idx_CHj)
      pdj(3) =  &
          +k(214)*n(idx_Oj)  &
          -k(1075)*n(idx_O)
      pdj(5) =  &
          +k(197)*n(idx_HCNj)  &
          +k(1034)*n(idx_CN)
      pdj(6) =  &
          +k(1162)  &
          +k(274)  &
          +k(985)*n(idx_H)  &
          +k(725)*n(idx_Nj)  &
          +k(111)*n(idx_H2j)  &
          +k(692)*n(idx_HEj)  &
          +k(375)*n(idx_Cj)
      pdj(7) =  &
          +k(20)*n(idx_Cj)
      pdj(8) =  &
          +k(272)  &
          +k(693)*n(idx_HEj)  &
          -k(985)*n(idx_H)  &
          +k(86)*n(idx_Hj)  &
          +k(1160)
      pdj(9) =  &
          +k(1099)*n(idx_OH)  &
          +k(196)*n(idx_H2Oj)
      pdj(10) =  &
          +k(1075)*n(idx_O)  &
          -k(1099)*n(idx_OH)  &
          +k(223)*n(idx_OHj)
      pdj(11) =  &
          +k(199)*n(idx_O2j)
      pdj(13) =  &
          +k(195)*n(idx_H2COj)
      pdj(16) =  &
          -k(1160)  &
          -k(86)*n(idx_Hj)  &
          -k(34)*n(idx_CHj)  &
          -k(272)  &
          -k(274)  &
          -k(909)*n(idx_CH3)  &
          -k(1075)*n(idx_O)  &
          -k(1038)*n(idx_NH)  &
          -k(182)*n(idx_NH2j)  &
          -k(692)*n(idx_HEj)  &
          -k(20)*n(idx_Cj)  &
          -k(726)*n(idx_Nj)  &
          -k(111)*n(idx_H2j)  &
          -k(793)*n(idx_COj)  &
          -k(794)*n(idx_HCNj)  &
          -k(194)*n(idx_COj)  &
          -k(195)*n(idx_H2COj)  &
          -k(166)*n(idx_Nj)  &
          -k(1161)  &
          -k(375)*n(idx_Cj)  &
          -k(985)*n(idx_H)  &
          -k(214)*n(idx_Oj)  &
          -k(1268)  &
          -k(196)*n(idx_H2Oj)  &
          -k(273)  &
          -k(199)*n(idx_O2j)  &
          -k(197)*n(idx_HCNj)  &
          -k(1034)*n(idx_CN)  &
          -k(1099)*n(idx_OH)  &
          -k(178)*n(idx_NHj)  &
          -k(693)*n(idx_HEj)  &
          -k(223)*n(idx_OHj)  &
          -k(725)*n(idx_Nj)  &
          -k(1162)  &
          -k(51)*n(idx_CH4j)  &
          -k(198)*n(idx_N2j)  &
          -k(146)*n(idx_HEj)
      pdj(24) =  &
          -k(1034)*n(idx_CN)
      pdj(25) =  &
          +k(194)*n(idx_COj)
      pdj(26) =  &
          +k(198)*n(idx_N2j)
      pdj(27) =  &
          +k(1099)*n(idx_OH)  &
          +k(1034)*n(idx_CN)  &
          +k(909)*n(idx_CH3)  &
          +k(794)*n(idx_HCNj)  &
          +k(1160)  &
          +k(985)*n(idx_H)  &
          +k(1075)*n(idx_O)  &
          +k(793)*n(idx_COj)  &
          +2.d0*k(1038)*n(idx_NH)  &
          +k(272)  &
          +k(182)*n(idx_NH2j)
      pdj(28) =  &
          -k(909)*n(idx_CH3)
      pdj(29) =  &
          +k(51)*n(idx_CH4j)  &
          +k(909)*n(idx_CH3)
      pdj(30) =  &
          +k(166)*n(idx_Nj)
      pdj(31) =  &
          -k(1038)*n(idx_NH)  &
          +k(178)*n(idx_NHj)  &
          +k(274)  &
          +k(726)*n(idx_Nj)  &
          +k(1162)
      pdj(35) =  &
          +k(692)*n(idx_HEj)  &
          +k(693)*n(idx_HEj)  &
          +k(146)*n(idx_HEj)
      pdj(60) =  &
          +k(1268)
      pdj(70) =  &
          +k(793)*n(idx_COj)
      pdj(71) =  &
          -k(86)*n(idx_Hj)
      pdj(73) =  &
          -k(375)*n(idx_Cj)  &
          -k(20)*n(idx_Cj)
      pdj(75) =  &
          -k(34)*n(idx_CHj)
      pdj(76) =  &
          -k(195)*n(idx_H2COj)
      pdj(78) =  &
          +k(214)*n(idx_Oj)  &
          +k(198)*n(idx_N2j)  &
          +k(197)*n(idx_HCNj)  &
          +k(223)*n(idx_OHj)  &
          +k(1161)  &
          +k(194)*n(idx_COj)  &
          +k(195)*n(idx_H2COj)  &
          +k(178)*n(idx_NHj)  &
          +k(51)*n(idx_CH4j)  &
          +k(182)*n(idx_NH2j)  &
          +k(146)*n(idx_HEj)  &
          +k(199)*n(idx_O2j)  &
          +k(196)*n(idx_H2Oj)  &
          +k(111)*n(idx_H2j)  &
          +k(34)*n(idx_CHj)  &
          +k(166)*n(idx_Nj)  &
          +k(20)*n(idx_Cj)  &
          +k(273)  &
          +k(86)*n(idx_Hj)
      pdj(87) =  &
          -k(793)*n(idx_COj)  &
          -k(194)*n(idx_COj)
      pdj(88) =  &
          -k(198)*n(idx_N2j)
      pdj(89) =  &
          -k(199)*n(idx_O2j)
      pdj(90) =  &
          -k(196)*n(idx_H2Oj)
      pdj(91) =  &
          +k(693)*n(idx_HEj)  &
          -k(182)*n(idx_NH2j)  &
          +k(726)*n(idx_Nj)
      pdj(92) =  &
          -k(214)*n(idx_Oj)
      pdj(93) =  &
          -k(223)*n(idx_OHj)
      pdj(95) =  &
          -k(51)*n(idx_CH4j)
      pdj(96) =  &
          -k(166)*n(idx_Nj)  &
          -k(726)*n(idx_Nj)  &
          -k(725)*n(idx_Nj)
      pdj(97) =  &
          -k(197)*n(idx_HCNj)  &
          -k(794)*n(idx_HCNj)  &
          +k(375)*n(idx_Cj)
      pdj(98) =  &
          -k(178)*n(idx_NHj)  &
          +k(692)*n(idx_HEj)
      pdj(102) =  &
          -k(111)*n(idx_H2j)
      pdj(103) =  &
          -k(146)*n(idx_HEj)  &
          -k(692)*n(idx_HEj)  &
          -k(693)*n(idx_HEj)
      pdj(109) =  &
          +k(794)*n(idx_HCNj)
      pdj(112) =  &
          +k(725)*n(idx_Nj)
    elseif(j==17) then
      pdj(1) =  &
          +k(278)  &
          +k(1166)
      pdj(2) =  &
          -k(932)*n(idx_CH)  &
          -k(931)*n(idx_CH)  &
          +k(35)*n(idx_CHj)  &
          -k(933)*n(idx_CH)
      pdj(3) =  &
          +k(728)*n(idx_Nj)  &
          +k(1167)  &
          +k(765)*n(idx_NHj)  &
          +k(931)*n(idx_CH)  &
          +k(1053)*n(idx_O2)  &
          +k(1023)*n(idx_N)  &
          +k(279)  &
          +k(1043)*n(idx_NH)  &
          +k(696)*n(idx_HEj)  &
          +k(988)*n(idx_H)  &
          +k(847)*n(idx_OHj)  &
          -k(1077)*n(idx_O)  &
          +k(872)*n(idx_C)
      pdj(5) =  &
          +k(931)*n(idx_CH)  &
          +k(133)*n(idx_HCNj)  &
          +k(911)*n(idx_CH3)  &
          +k(888)*n(idx_CH2)
      pdj(6) =  &
          +k(597)*n(idx_H3j)  &
          +k(113)*n(idx_H2j)
      pdj(7) =  &
          -k(872)*n(idx_C)  &
          -k(873)*n(idx_C)  &
          +k(21)*n(idx_Cj)
      pdj(8) =  &
          +k(933)*n(idx_CH)  &
          -k(989)*n(idx_H)  &
          +k(88)*n(idx_Hj)  &
          +k(1100)*n(idx_OH)  &
          +k(889)*n(idx_CH2)  &
          +k(1031)*n(idx_NH2)  &
          +k(525)*n(idx_H2j)  &
          +k(1043)*n(idx_NH)  &
          -k(988)*n(idx_H)
      pdj(9) =  &
          +k(911)*n(idx_CH3)  &
          +k(121)*n(idx_H2Oj)  &
          +k(1030)*n(idx_NH2)
      pdj(10) =  &
          +k(888)*n(idx_CH2)  &
          +k(1044)*n(idx_NH)  &
          -k(1100)*n(idx_OH)  &
          +k(989)*n(idx_H)  &
          +k(1031)*n(idx_NH2)  &
          +k(224)*n(idx_OHj)
      pdj(11) =  &
          +2.d0*k(1052)*n(idx_NO)  &
          +k(808)*n(idx_O2Hj)  &
          +k(206)*n(idx_O2j)  &
          +k(1077)*n(idx_O)  &
          -k(1053)*n(idx_O2)
      pdj(12) =  &
          -k(887)*n(idx_CH2)  &
          -k(889)*n(idx_CH2)  &
          -k(888)*n(idx_CH2)  &
          +k(37)*n(idx_CH2j)
      pdj(13) =  &
          +k(204)*n(idx_H2COj)  &
          +k(887)*n(idx_CH2)
      pdj(14) =  &
          -k(1001)*n(idx_HCO)  &
          +k(932)*n(idx_CH)
      pdj(16) =  &
          +k(192)*n(idx_NH3j)
      pdj(17) =  &
          -k(808)*n(idx_O2Hj)  &
          -k(179)*n(idx_NHj)  &
          -k(1001)*n(idx_HCO)  &
          -k(765)*n(idx_NHj)  &
          -k(696)*n(idx_HEj)  &
          -k(173)*n(idx_N2j)  &
          -k(847)*n(idx_OHj)  &
          -k(525)*n(idx_H2j)  &
          -k(1256)  &
          -k(1044)*n(idx_NH)  &
          -k(1077)*n(idx_O)  &
          -k(988)*n(idx_H)  &
          -k(133)*n(idx_HCNj)  &
          -k(113)*n(idx_H2j)  &
          -k(1053)*n(idx_O2)  &
          -k(1167)  &
          -k(931)*n(idx_CH)  &
          -k(947)*n(idx_CN)  &
          -k(888)*n(idx_CH2)  &
          -k(35)*n(idx_CHj)  &
          -k(889)*n(idx_CH2)  &
          -k(168)*n(idx_Nj)  &
          -k(88)*n(idx_Hj)  &
          -k(1043)*n(idx_NH)  &
          -k(1023)*n(idx_N)  &
          -k(1100)*n(idx_OH)  &
          -k(1106)*n(idx_SI)  &
          -k(278)  &
          -k(183)*n(idx_NH2j)  &
          -k(1054)*n(idx_OCN)  &
          -k(73)*n(idx_COj)  &
          -k(989)*n(idx_H)  &
          -k(911)*n(idx_CH3)  &
          -k(728)*n(idx_Nj)  &
          -4.d0*k(1052)*n(idx_NO)  &
          -k(21)*n(idx_Cj)  &
          -k(948)*n(idx_CN)  &
          -k(1030)*n(idx_NH2)  &
          -k(207)*n(idx_SIOj)  &
          -k(224)*n(idx_OHj)  &
          -k(887)*n(idx_CH2)  &
          -k(37)*n(idx_CH2j)  &
          -k(205)*n(idx_HNOj)  &
          -k(49)*n(idx_CH3j)  &
          -k(695)*n(idx_HEj)  &
          -k(68)*n(idx_CNj)  &
          -k(872)*n(idx_C)  &
          -k(1166)  &
          -k(933)*n(idx_CH)  &
          -k(597)*n(idx_H3j)  &
          -k(121)*n(idx_H2Oj)  &
          -k(192)*n(idx_NH3j)  &
          -k(1031)*n(idx_NH2)  &
          -k(932)*n(idx_CH)  &
          -k(206)*n(idx_O2j)  &
          -k(204)*n(idx_H2COj)  &
          -k(279)  &
          -k(873)*n(idx_C)
      pdj(18) =  &
          -k(1106)*n(idx_SI)
      pdj(24) =  &
          -k(947)*n(idx_CN)  &
          -k(948)*n(idx_CN)  &
          +k(872)*n(idx_C)  &
          +k(68)*n(idx_CNj)
      pdj(25) =  &
          +k(1001)*n(idx_HCO)  &
          +k(73)*n(idx_COj)  &
          +k(947)*n(idx_CN)  &
          +k(873)*n(idx_C)
      pdj(26) =  &
          +k(1044)*n(idx_NH)  &
          +2.d0*k(1052)*n(idx_NO)  &
          +k(1023)*n(idx_N)  &
          +k(947)*n(idx_CN)  &
          +k(173)*n(idx_N2j)  &
          +k(1030)*n(idx_NH2)  &
          +k(1054)*n(idx_OCN)  &
          +k(1031)*n(idx_NH2)  &
          +k(1043)*n(idx_NH)
      pdj(27) =  &
          +k(183)*n(idx_NH2j)  &
          -k(1030)*n(idx_NH2)  &
          -k(1031)*n(idx_NH2)
      pdj(28) =  &
          -k(911)*n(idx_CH3)  &
          +k(49)*n(idx_CH3j)
      pdj(30) =  &
          +k(948)*n(idx_CN)  &
          +k(168)*n(idx_Nj)  &
          +k(873)*n(idx_C)  &
          +k(989)*n(idx_H)  &
          +k(695)*n(idx_HEj)  &
          +k(1077)*n(idx_O)  &
          +k(932)*n(idx_CH)  &
          +k(1167)  &
          -k(1023)*n(idx_N)  &
          +k(1106)*n(idx_SI)  &
          +k(887)*n(idx_CH2)  &
          +k(279)
      pdj(31) =  &
          -k(1043)*n(idx_NH)  &
          +k(988)*n(idx_H)  &
          -k(1044)*n(idx_NH)  &
          +k(179)*n(idx_NHj)
      pdj(34) =  &
          +k(1106)*n(idx_SI)  &
          +k(207)*n(idx_SIOj)
      pdj(35) =  &
          +k(695)*n(idx_HEj)  &
          +k(696)*n(idx_HEj)
      pdj(36) =  &
          +k(1001)*n(idx_HCO)  &
          +k(205)*n(idx_HNOj)
      pdj(38) =  &
          +k(1054)*n(idx_OCN)
      pdj(41) =  &
          +k(889)*n(idx_CH2)
      pdj(42) =  &
          +k(1053)*n(idx_O2)  &
          +k(1100)*n(idx_OH)
      pdj(44) =  &
          +k(933)*n(idx_CH)  &
          -k(1054)*n(idx_OCN)  &
          +k(948)*n(idx_CN)
      pdj(56) =  &
          +k(1256)
      pdj(71) =  &
          -k(88)*n(idx_Hj)
      pdj(73) =  &
          -k(21)*n(idx_Cj)
      pdj(74) =  &
          -k(37)*n(idx_CH2j)
      pdj(75) =  &
          -k(35)*n(idx_CHj)
      pdj(76) =  &
          -k(204)*n(idx_H2COj)
      pdj(78) =  &
          -k(192)*n(idx_NH3j)
      pdj(79) =  &
          +k(168)*n(idx_Nj)  &
          +k(207)*n(idx_SIOj)  &
          +k(192)*n(idx_NH3j)  &
          +k(73)*n(idx_COj)  &
          +k(88)*n(idx_Hj)  &
          +k(68)*n(idx_CNj)  &
          +k(204)*n(idx_H2COj)  &
          +k(179)*n(idx_NHj)  &
          +k(278)  &
          +k(35)*n(idx_CHj)  &
          +k(173)*n(idx_N2j)  &
          +k(49)*n(idx_CH3j)  &
          +k(205)*n(idx_HNOj)  &
          +k(224)*n(idx_OHj)  &
          +k(183)*n(idx_NH2j)  &
          +k(21)*n(idx_Cj)  &
          +k(121)*n(idx_H2Oj)  &
          +k(37)*n(idx_CH2j)  &
          +k(206)*n(idx_O2j)  &
          +k(1166)  &
          +k(113)*n(idx_H2j)  &
          +k(133)*n(idx_HCNj)
      pdj(86) =  &
          -k(68)*n(idx_CNj)
      pdj(87) =  &
          -k(73)*n(idx_COj)
      pdj(88) =  &
          +k(728)*n(idx_Nj)  &
          -k(173)*n(idx_N2j)
      pdj(89) =  &
          -k(206)*n(idx_O2j)
      pdj(90) =  &
          -k(121)*n(idx_H2Oj)
      pdj(91) =  &
          -k(183)*n(idx_NH2j)
      pdj(92) =  &
          +k(695)*n(idx_HEj)
      pdj(93) =  &
          -k(847)*n(idx_OHj)  &
          -k(224)*n(idx_OHj)
      pdj(94) =  &
          -k(49)*n(idx_CH3j)
      pdj(96) =  &
          -k(168)*n(idx_Nj)  &
          -k(728)*n(idx_Nj)  &
          +k(696)*n(idx_HEj)
      pdj(97) =  &
          -k(133)*n(idx_HCNj)
      pdj(98) =  &
          -k(765)*n(idx_NHj)  &
          -k(179)*n(idx_NHj)
      pdj(101) =  &
          -k(207)*n(idx_SIOj)
      pdj(102) =  &
          -k(525)*n(idx_H2j)  &
          -k(113)*n(idx_H2j)
      pdj(103) =  &
          -k(696)*n(idx_HEj)  &
          -k(695)*n(idx_HEj)
      pdj(104) =  &
          +k(847)*n(idx_OHj)  &
          +k(597)*n(idx_H3j)  &
          +k(525)*n(idx_H2j)  &
          -k(205)*n(idx_HNOj)  &
          +k(808)*n(idx_O2Hj)
      pdj(106) =  &
          -k(597)*n(idx_H3j)
      pdj(112) =  &
          +k(765)*n(idx_NHj)
      pdj(113) =  &
          -k(808)*n(idx_O2Hj)
    elseif(j==18) then
      pdj(1) =  &
          +k(286)  &
          +k(1177)
      pdj(2) =  &
          +k(36)*n(idx_CHj)
      pdj(3) =  &
          +k(849)*n(idx_OHj)  &
          +k(1107)*n(idx_O2)  &
          -k(1214)*n(idx_O)
      pdj(6) =  &
          +k(602)*n(idx_H3j)
      pdj(7) =  &
          +k(22)*n(idx_Cj)  &
          +k(1105)*n(idx_CO)
      pdj(8) =  &
          +k(1103)*n(idx_OH)  &
          +k(92)*n(idx_Hj)
      pdj(9) =  &
          +k(611)*n(idx_H3Oj)  &
          +k(123)*n(idx_H2Oj)
      pdj(10) =  &
          -k(1103)*n(idx_OH)
      pdj(11) =  &
          -k(1107)*n(idx_O2)  &
          +k(231)*n(idx_O2j)
      pdj(13) =  &
          +k(229)*n(idx_H2COj)
      pdj(16) =  &
          +k(193)*n(idx_NH3j)
      pdj(17) =  &
          -k(1106)*n(idx_NO)  &
          +k(230)*n(idx_NOj)
      pdj(18) =  &
          -k(1107)*n(idx_O2)  &
          -k(1214)*n(idx_O)  &
          -k(1106)*n(idx_NO)  &
          -k(1229)  &
          -k(849)*n(idx_OHj)  &
          -k(36)*n(idx_CHj)  &
          -k(1177)  &
          -k(1105)*n(idx_CO)  &
          -k(862)*n(idx_HCOj)  &
          -k(148)*n(idx_HEj)  &
          -k(602)*n(idx_H3j)  &
          -k(229)*n(idx_H2COj)  &
          -k(286)  &
          -k(611)*n(idx_H3Oj)  &
          -k(1104)*n(idx_CO2)  &
          -k(92)*n(idx_Hj)  &
          -k(231)*n(idx_O2j)  &
          -k(123)*n(idx_H2Oj)  &
          -k(193)*n(idx_NH3j)  &
          -k(1103)*n(idx_OH)  &
          -k(22)*n(idx_Cj)  &
          -k(230)*n(idx_NOj)
      pdj(25) =  &
          +k(1104)*n(idx_CO2)  &
          +k(862)*n(idx_HCOj)  &
          -k(1105)*n(idx_CO)
      pdj(30) =  &
          +k(1106)*n(idx_NO)
      pdj(34) =  &
          +k(1214)*n(idx_O)  &
          +k(1103)*n(idx_OH)  &
          +k(1104)*n(idx_CO2)  &
          +k(1107)*n(idx_O2)  &
          +k(1106)*n(idx_NO)  &
          +k(1105)*n(idx_CO)
      pdj(35) =  &
          +k(148)*n(idx_HEj)
      pdj(38) =  &
          -k(1104)*n(idx_CO2)
      pdj(48) =  &
          +k(1229)
      pdj(70) =  &
          -k(862)*n(idx_HCOj)
      pdj(71) =  &
          -k(92)*n(idx_Hj)
      pdj(73) =  &
          -k(22)*n(idx_Cj)
      pdj(75) =  &
          -k(36)*n(idx_CHj)
      pdj(76) =  &
          -k(229)*n(idx_H2COj)
      pdj(78) =  &
          -k(193)*n(idx_NH3j)
      pdj(79) =  &
          -k(230)*n(idx_NOj)
      pdj(80) =  &
          +k(148)*n(idx_HEj)  &
          +k(286)  &
          +k(123)*n(idx_H2Oj)  &
          +k(92)*n(idx_Hj)  &
          +k(231)*n(idx_O2j)  &
          +k(230)*n(idx_NOj)  &
          +k(1177)  &
          +k(36)*n(idx_CHj)  &
          +k(22)*n(idx_Cj)  &
          +k(229)*n(idx_H2COj)  &
          +k(193)*n(idx_NH3j)
      pdj(89) =  &
          -k(231)*n(idx_O2j)
      pdj(90) =  &
          -k(123)*n(idx_H2Oj)
      pdj(93) =  &
          -k(849)*n(idx_OHj)
      pdj(100) =  &
          +k(849)*n(idx_OHj)  &
          +k(611)*n(idx_H3Oj)  &
          +k(862)*n(idx_HCOj)  &
          +k(602)*n(idx_H3j)
      pdj(103) =  &
          -k(148)*n(idx_HEj)
      pdj(106) =  &
          -k(602)*n(idx_H3j)
      pdj(108) =  &
          -k(611)*n(idx_H3Oj)
    elseif(j==19) then
      pdj(3) =  &
          -k(1082)*n(idx_O)
      pdj(7) =  &
          +k(23)*n(idx_Cj)  &
          +k(287)
      pdj(8) =  &
          +k(93)*n(idx_Hj)
      pdj(19) =  &
          -k(1241)  &
          -k(93)*n(idx_Hj)  &
          -k(287)  &
          -k(23)*n(idx_Cj)  &
          -k(1082)*n(idx_O)
      pdj(21) =  &
          +k(287)  &
          +k(1082)*n(idx_O)
      pdj(25) =  &
          +k(1082)*n(idx_O)
      pdj(51) =  &
          +k(1241)
      pdj(71) =  &
          -k(93)*n(idx_Hj)
      pdj(73) =  &
          -k(23)*n(idx_Cj)
      pdj(81) =  &
          +k(23)*n(idx_Cj)  &
          +k(93)*n(idx_Hj)
    elseif(j==20) then
      pdj(3) =  &
          -k(1083)*n(idx_O)
      pdj(7) =  &
          +k(1178)  &
          +k(701)*n(idx_HEj)  &
          +k(24)*n(idx_Cj)  &
          +k(288)
      pdj(8) =  &
          +k(94)*n(idx_Hj)
      pdj(19) =  &
          +k(1083)*n(idx_O)  &
          +k(1178)  &
          +k(288)
      pdj(20) =  &
          -k(1244)  &
          -k(94)*n(idx_Hj)  &
          -k(288)  &
          -k(1178)  &
          -k(24)*n(idx_Cj)  &
          -k(1083)*n(idx_O)  &
          -k(701)*n(idx_HEj)
      pdj(25) =  &
          +k(1083)*n(idx_O)
      pdj(35) =  &
          +k(701)*n(idx_HEj)
      pdj(52) =  &
          +k(1244)
      pdj(71) =  &
          -k(94)*n(idx_Hj)
      pdj(73) =  &
          -k(24)*n(idx_Cj)
      pdj(81) =  &
          +k(701)*n(idx_HEj)
      pdj(82) =  &
          +k(24)*n(idx_Cj)  &
          +k(94)*n(idx_Hj)
      pdj(103) =  &
          -k(701)*n(idx_HEj)
    elseif(j==21) then
      pdj(3) =  &
          -k(1085)*n(idx_O)  &
          -k(1084)*n(idx_O)
      pdj(7) =  &
          +k(702)*n(idx_HEj)  &
          +k(25)*n(idx_Cj)  &
          +k(1179)  &
          +k(289)  &
          +k(1085)*n(idx_O)
      pdj(8) =  &
          +k(95)*n(idx_Hj)
      pdj(18) =  &
          +k(703)*n(idx_HEj)  &
          +k(1028)*n(idx_N)  &
          +k(1179)  &
          +k(289)  &
          +k(1084)*n(idx_O)
      pdj(21) =  &
          -k(1179)  &
          -k(95)*n(idx_Hj)  &
          -k(1084)*n(idx_O)  &
          -k(702)*n(idx_HEj)  &
          -k(289)  &
          -k(703)*n(idx_HEj)  &
          -k(1085)*n(idx_O)  &
          -k(1028)*n(idx_N)  &
          -k(1240)  &
          -k(25)*n(idx_Cj)
      pdj(24) =  &
          +k(1028)*n(idx_N)
      pdj(25) =  &
          +k(1084)*n(idx_O)
      pdj(30) =  &
          -k(1028)*n(idx_N)
      pdj(34) =  &
          +k(1085)*n(idx_O)
      pdj(35) =  &
          +k(702)*n(idx_HEj)  &
          +k(703)*n(idx_HEj)
      pdj(50) =  &
          +k(1240)
      pdj(71) =  &
          -k(95)*n(idx_Hj)
      pdj(73) =  &
          +k(703)*n(idx_HEj)  &
          -k(25)*n(idx_Cj)
      pdj(80) =  &
          +k(702)*n(idx_HEj)
      pdj(83) =  &
          +k(25)*n(idx_Cj)  &
          +k(95)*n(idx_Hj)
      pdj(103) =  &
          -k(702)*n(idx_HEj)  &
          -k(703)*n(idx_HEj)
    elseif(j==22) then
      pdj(1) =  &
          +k(1181)
      pdj(3) =  &
          -k(1086)*n(idx_O)  &
          -k(1087)*n(idx_O)
      pdj(6) =  &
          +k(381)*n(idx_Cj)  &
          +k(506)*n(idx_Hj)  &
          +k(704)*n(idx_HEj)  &
          +k(603)*n(idx_H3j)  &
          +k(1086)*n(idx_O)
      pdj(7) =  &
          +k(26)*n(idx_Cj)
      pdj(8) =  &
          +2.d0*k(1087)*n(idx_O)  &
          +k(96)*n(idx_Hj)  &
          +k(705)*n(idx_HEj)  &
          +k(1182)  &
          +k(290)
      pdj(9) =  &
          +k(612)*n(idx_H3Oj)
      pdj(22) =  &
          -k(381)*n(idx_Cj)  &
          -k(603)*n(idx_H3j)  &
          -k(96)*n(idx_Hj)  &
          -k(1086)*n(idx_O)  &
          -k(290)  &
          -k(612)*n(idx_H3Oj)  &
          -k(1182)  &
          -k(704)*n(idx_HEj)  &
          -k(705)*n(idx_HEj)  &
          -k(506)*n(idx_Hj)  &
          -k(1181)  &
          -k(638)*n(idx_HCOj)  &
          -k(1087)*n(idx_O)  &
          -k(1234)  &
          -k(26)*n(idx_Cj)
      pdj(25) =  &
          +k(638)*n(idx_HCOj)
      pdj(33) =  &
          +k(1182)  &
          +k(290)
      pdj(34) =  &
          +k(1086)*n(idx_O)  &
          +k(1087)*n(idx_O)
      pdj(35) =  &
          +k(704)*n(idx_HEj)  &
          +k(705)*n(idx_HEj)
      pdj(48) =  &
          +k(1234)
      pdj(70) =  &
          -k(638)*n(idx_HCOj)
      pdj(71) =  &
          -k(96)*n(idx_Hj)  &
          -k(506)*n(idx_Hj)
      pdj(73) =  &
          -k(26)*n(idx_Cj)  &
          -k(381)*n(idx_Cj)
      pdj(80) =  &
          +k(704)*n(idx_HEj)
      pdj(83) =  &
          +k(381)*n(idx_Cj)
      pdj(84) =  &
          +k(1181)  &
          +k(96)*n(idx_Hj)  &
          +k(26)*n(idx_Cj)
      pdj(85) =  &
          +k(638)*n(idx_HCOj)  &
          +k(612)*n(idx_H3Oj)  &
          +k(603)*n(idx_H3j)
      pdj(100) =  &
          +k(506)*n(idx_Hj)  &
          +k(705)*n(idx_HEj)
      pdj(103) =  &
          -k(704)*n(idx_HEj)  &
          -k(705)*n(idx_HEj)
      pdj(106) =  &
          -k(603)*n(idx_H3j)
      pdj(108) =  &
          -k(612)*n(idx_H3Oj)
    elseif(j==23) then
      pdj(1) =  &
          +k(1184)
      pdj(3) =  &
          -k(1088)*n(idx_O)
      pdj(6) =  &
          +k(604)*n(idx_H3j)  &
          +k(706)*n(idx_HEj)  &
          +k(1185)  &
          +k(507)*n(idx_Hj)
      pdj(7) =  &
          +k(27)*n(idx_Cj)
      pdj(8) =  &
          +k(97)*n(idx_Hj)  &
          +k(291)  &
          +k(1088)*n(idx_O)  &
          +k(707)*n(idx_HEj)  &
          +k(1183)
      pdj(22) =  &
          +k(291)  &
          +k(1183)
      pdj(23) =  &
          -k(507)*n(idx_Hj)  &
          -k(707)*n(idx_HEj)  &
          -k(1184)  &
          -k(1185)  &
          -k(706)*n(idx_HEj)  &
          -k(291)  &
          -k(604)*n(idx_H3j)  &
          -k(1088)*n(idx_O)  &
          -k(1183)  &
          -k(97)*n(idx_Hj)  &
          -k(27)*n(idx_Cj)  &
          -k(1236)
      pdj(33) =  &
          +k(1185)
      pdj(35) =  &
          +k(706)*n(idx_HEj)  &
          +k(707)*n(idx_HEj)
      pdj(40) =  &
          +k(1088)*n(idx_O)
      pdj(48) =  &
          +k(1236)
      pdj(71) =  &
          -k(507)*n(idx_Hj)  &
          -k(97)*n(idx_Hj)
      pdj(73) =  &
          -k(27)*n(idx_Cj)
      pdj(84) =  &
          +k(707)*n(idx_HEj)  &
          +k(507)*n(idx_Hj)
      pdj(85) =  &
          +k(97)*n(idx_Hj)  &
          +k(27)*n(idx_Cj)  &
          +k(1184)
      pdj(99) =  &
          +k(604)*n(idx_H3j)
      pdj(100) =  &
          +k(706)*n(idx_HEj)
      pdj(103) =  &
          -k(706)*n(idx_HEj)  &
          -k(707)*n(idx_HEj)
      pdj(106) =  &
          -k(604)*n(idx_H3j)
    elseif(j==24) then
      pdj(2) =  &
          +k(881)*n(idx_CH2)
      pdj(3) =  &
          +k(950)*n(idx_O2)  &
          +k(1091)*n(idx_OH)  &
          +k(837)*n(idx_OHj)  &
          -k(1058)*n(idx_O)  &
          -k(1059)*n(idx_O)
      pdj(5) =  &
          +k(944)*n(idx_HCO)  &
          +k(1034)*n(idx_NH3)  &
          +k(943)*n(idx_H2CO)  &
          +k(1091)*n(idx_OH)  &
          +k(945)*n(idx_HNO)  &
          +k(1036)*n(idx_NH)  &
          +k(881)*n(idx_CH2)  &
          +k(903)*n(idx_CH3)  &
          +k(921)*n(idx_CH4)  &
          +k(951)*n(idx_SIH4)  &
          +k(960)*n(idx_H2)
      pdj(6) =  &
          +k(104)*n(idx_H2j)  &
          +k(582)*n(idx_H3j)  &
          -k(960)*n(idx_H2)
      pdj(7) =  &
          +k(812)*n(idx_Oj)  &
          +k(1059)*n(idx_O)  &
          +k(1131)  &
          +k(252)  &
          +k(664)*n(idx_HEj)  &
          +k(1012)*n(idx_N)
      pdj(8) =  &
          +k(960)*n(idx_H2)  &
          +k(514)*n(idx_H2j)  &
          +k(1092)*n(idx_OH)
      pdj(10) =  &
          -k(1091)*n(idx_OH)  &
          -k(1092)*n(idx_OH)
      pdj(11) =  &
          +k(484)*n(idx_O2Hj)  &
          -k(949)*n(idx_O2)  &
          -k(950)*n(idx_O2)
      pdj(12) =  &
          -k(881)*n(idx_CH2)  &
          +k(903)*n(idx_CH3)
      pdj(13) =  &
          -k(943)*n(idx_H2CO)
      pdj(14) =  &
          -k(944)*n(idx_HCO)  &
          +k(943)*n(idx_H2CO)
      pdj(16) =  &
          -k(1034)*n(idx_NH3)
      pdj(17) =  &
          +k(946)*n(idx_NO2)  &
          -k(947)*n(idx_NO)  &
          +k(945)*n(idx_HNO)  &
          +k(949)*n(idx_O2)  &
          +k(483)*n(idx_HNOj)  &
          +k(1059)*n(idx_O)  &
          -k(948)*n(idx_NO)
      pdj(23) =  &
          +k(951)*n(idx_SIH4)
      pdj(24) =  &
          -k(748)*n(idx_NHj)  &
          -k(1034)*n(idx_NH3)  &
          -k(950)*n(idx_O2)  &
          -k(664)*n(idx_HEj)  &
          -k(837)*n(idx_OHj)  &
          -k(949)*n(idx_O2)  &
          -k(903)*n(idx_CH3)  &
          -k(484)*n(idx_O2Hj)  &
          -k(960)*n(idx_H2)  &
          -k(944)*n(idx_HCO)  &
          -k(943)*n(idx_H2CO)  &
          -k(921)*n(idx_CH4)  &
          -k(947)*n(idx_NO)  &
          -k(1036)*n(idx_NH)  &
          -k(951)*n(idx_SIH4)  &
          -k(812)*n(idx_Oj)  &
          -k(665)*n(idx_HEj)  &
          -k(1091)*n(idx_OH)  &
          -k(881)*n(idx_CH2)  &
          -k(1012)*n(idx_N)  &
          -k(70)*n(idx_N2j)  &
          -k(252)  &
          -k(483)*n(idx_HNOj)  &
          -k(1058)*n(idx_O)  &
          -k(1059)*n(idx_O)  &
          -k(1264)  &
          -k(1092)*n(idx_OH)  &
          -k(948)*n(idx_NO)  &
          -k(158)*n(idx_Nj)  &
          -k(945)*n(idx_HNO)  &
          -k(104)*n(idx_H2j)  &
          -k(1131)  &
          -k(582)*n(idx_H3j)  &
          -k(514)*n(idx_H2j)  &
          -k(946)*n(idx_NO2)
      pdj(25) =  &
          +k(1058)*n(idx_O)  &
          +k(944)*n(idx_HCO)  &
          +k(947)*n(idx_NO)  &
          +k(949)*n(idx_O2)
      pdj(26) =  &
          +k(70)*n(idx_N2j)  &
          +k(947)*n(idx_NO)  &
          +k(1012)*n(idx_N)
      pdj(27) =  &
          +k(1034)*n(idx_NH3)
      pdj(28) =  &
          -k(903)*n(idx_CH3)  &
          +k(921)*n(idx_CH4)
      pdj(29) =  &
          -k(921)*n(idx_CH4)
      pdj(30) =  &
          +k(158)*n(idx_Nj)  &
          +k(1058)*n(idx_O)  &
          +k(748)*n(idx_NHj)  &
          +k(1036)*n(idx_NH)  &
          +k(1131)  &
          +k(948)*n(idx_NO)  &
          +k(252)  &
          +k(665)*n(idx_HEj)  &
          -k(1012)*n(idx_N)
      pdj(31) =  &
          -k(1036)*n(idx_NH)
      pdj(32) =  &
          -k(951)*n(idx_SIH4)
      pdj(35) =  &
          +k(665)*n(idx_HEj)  &
          +k(664)*n(idx_HEj)
      pdj(36) =  &
          -k(945)*n(idx_HNO)
      pdj(42) =  &
          -k(946)*n(idx_NO2)
      pdj(44) =  &
          +k(950)*n(idx_O2)  &
          +k(948)*n(idx_NO)  &
          +k(1092)*n(idx_OH)  &
          +k(946)*n(idx_NO2)
      pdj(59) =  &
          +k(1264)
      pdj(73) =  &
          +k(665)*n(idx_HEj)
      pdj(79) =  &
          +k(812)*n(idx_Oj)
      pdj(86) =  &
          +k(70)*n(idx_N2j)  &
          +k(104)*n(idx_H2j)  &
          +k(158)*n(idx_Nj)
      pdj(88) =  &
          -k(70)*n(idx_N2j)
      pdj(92) =  &
          -k(812)*n(idx_Oj)
      pdj(93) =  &
          -k(837)*n(idx_OHj)
      pdj(96) =  &
          +k(664)*n(idx_HEj)  &
          -k(158)*n(idx_Nj)
      pdj(97) =  &
          +k(484)*n(idx_O2Hj)  &
          +k(748)*n(idx_NHj)  &
          +k(483)*n(idx_HNOj)  &
          +k(582)*n(idx_H3j)  &
          +k(514)*n(idx_H2j)  &
          +k(837)*n(idx_OHj)
      pdj(98) =  &
          -k(748)*n(idx_NHj)
      pdj(102) =  &
          -k(104)*n(idx_H2j)  &
          -k(514)*n(idx_H2j)
      pdj(103) =  &
          -k(665)*n(idx_HEj)  &
          -k(664)*n(idx_HEj)
      pdj(104) =  &
          -k(483)*n(idx_HNOj)
      pdj(106) =  &
          -k(582)*n(idx_H3j)
      pdj(113) =  &
          -k(484)*n(idx_O2Hj)
    elseif(j==25) then
      pdj(1) =  &
          +k(233)
      pdj(3) =  &
          +k(209)*n(idx_Oj)  &
          +k(839)*n(idx_OHj)  &
          +k(670)*n(idx_HEj)  &
          +k(954)*n(idx_O2)  &
          +k(1134)  &
          +k(254)
      pdj(6) =  &
          +k(105)*n(idx_H2j)  &
          +k(585)*n(idx_H3j)  &
          +k(584)*n(idx_H3j)
      pdj(7) =  &
          +k(973)*n(idx_H)  &
          +k(721)*n(idx_Nj)  &
          +k(1105)*n(idx_SI)  &
          +k(1134)  &
          +k(254)
      pdj(8) =  &
          -k(973)*n(idx_H)  &
          +k(1093)*n(idx_OH)  &
          +k(516)*n(idx_H2j)
      pdj(10) =  &
          +k(554)*n(idx_H2Oj)  &
          -k(1093)*n(idx_OH)  &
          +k(955)*n(idx_O2H)  &
          +k(973)*n(idx_H)
      pdj(11) =  &
          +k(489)*n(idx_O2Hj)  &
          -k(954)*n(idx_O2)
      pdj(17) =  &
          +k(953)*n(idx_NO2)  &
          +k(487)*n(idx_HNOj)
      pdj(18) =  &
          -k(1105)*n(idx_SI)
      pdj(23) =  &
          +k(490)*n(idx_SIH4j)
      pdj(24) =  &
          +k(64)*n(idx_CNj)  &
          +k(622)*n(idx_HCNj)
      pdj(25) =  &
          -k(64)*n(idx_CNj)  &
          -k(105)*n(idx_H2j)  &
          -k(839)*n(idx_OHj)  &
          -k(953)*n(idx_NO2)  &
          -k(488)*n(idx_N2Hj)  &
          -k(622)*n(idx_HCNj)  &
          -k(721)*n(idx_Nj)  &
          -k(554)*n(idx_H2Oj)  &
          -k(973)*n(idx_H)  &
          -k(209)*n(idx_Oj)  &
          -k(1105)*n(idx_SI)  &
          -k(1093)*n(idx_OH)  &
          -k(75)*n(idx_N2j)  &
          -k(159)*n(idx_Nj)  &
          -k(954)*n(idx_O2)  &
          -k(491)*n(idx_SIOj)  &
          -k(486)*n(idx_HCO2j)  &
          -k(585)*n(idx_H3j)  &
          -k(1252)  &
          -k(1228)  &
          -k(955)*n(idx_O2H)  &
          -k(670)*n(idx_HEj)  &
          -k(1305)  &
          -k(516)*n(idx_H2j)  &
          -k(752)*n(idx_NHj)  &
          -k(487)*n(idx_HNOj)  &
          -k(584)*n(idx_H3j)  &
          -k(1134)  &
          -k(952)*n(idx_HNO)  &
          -k(233)  &
          -k(254)  &
          -k(489)*n(idx_O2Hj)  &
          -k(490)*n(idx_SIH4j)  &
          -k(449)*n(idx_CH4j)
      pdj(26) =  &
          +k(75)*n(idx_N2j)  &
          +k(488)*n(idx_N2Hj)
      pdj(28) =  &
          +k(449)*n(idx_CH4j)
      pdj(30) =  &
          +k(159)*n(idx_Nj)  &
          +k(752)*n(idx_NHj)
      pdj(31) =  &
          +k(952)*n(idx_HNO)
      pdj(34) =  &
          +k(1105)*n(idx_SI)
      pdj(35) =  &
          +k(670)*n(idx_HEj)
      pdj(36) =  &
          -k(952)*n(idx_HNO)
      pdj(38) =  &
          +k(491)*n(idx_SIOj)  &
          +k(952)*n(idx_HNO)  &
          +k(955)*n(idx_O2H)  &
          +k(486)*n(idx_HCO2j)  &
          +k(1093)*n(idx_OH)  &
          +k(953)*n(idx_NO2)  &
          +k(954)*n(idx_O2)
      pdj(42) =  &
          -k(953)*n(idx_NO2)
      pdj(43) =  &
          -k(955)*n(idx_O2H)
      pdj(45) =  &
          +k(1305)
      pdj(47) =  &
          +k(1228)
      pdj(54) =  &
          +k(1252)
      pdj(70) =  &
          +k(584)*n(idx_H3j)  &
          +k(488)*n(idx_N2Hj)  &
          +k(554)*n(idx_H2Oj)  &
          +k(489)*n(idx_O2Hj)  &
          +k(487)*n(idx_HNOj)  &
          +k(486)*n(idx_HCO2j)  &
          +k(516)*n(idx_H2j)  &
          +k(490)*n(idx_SIH4j)  &
          +k(449)*n(idx_CH4j)  &
          +k(839)*n(idx_OHj)  &
          +k(752)*n(idx_NHj)  &
          +k(622)*n(idx_HCNj)
      pdj(72) =  &
          +k(585)*n(idx_H3j)
      pdj(73) =  &
          +k(670)*n(idx_HEj)
      pdj(79) =  &
          +k(721)*n(idx_Nj)
      pdj(80) =  &
          +k(491)*n(idx_SIOj)
      pdj(86) =  &
          -k(64)*n(idx_CNj)
      pdj(87) =  &
          +k(233)  &
          +k(159)*n(idx_Nj)  &
          +k(209)*n(idx_Oj)  &
          +k(105)*n(idx_H2j)  &
          +k(64)*n(idx_CNj)  &
          +k(75)*n(idx_N2j)
      pdj(88) =  &
          -k(75)*n(idx_N2j)
      pdj(90) =  &
          -k(554)*n(idx_H2Oj)
      pdj(92) =  &
          -k(209)*n(idx_Oj)
      pdj(93) =  &
          -k(839)*n(idx_OHj)
      pdj(95) =  &
          -k(449)*n(idx_CH4j)
      pdj(96) =  &
          -k(721)*n(idx_Nj)  &
          -k(159)*n(idx_Nj)
      pdj(97) =  &
          -k(622)*n(idx_HCNj)
      pdj(98) =  &
          -k(752)*n(idx_NHj)
      pdj(99) =  &
          -k(490)*n(idx_SIH4j)
      pdj(101) =  &
          -k(491)*n(idx_SIOj)
      pdj(102) =  &
          -k(516)*n(idx_H2j)  &
          -k(105)*n(idx_H2j)
      pdj(103) =  &
          -k(670)*n(idx_HEj)
      pdj(104) =  &
          -k(487)*n(idx_HNOj)
      pdj(106) =  &
          -k(584)*n(idx_H3j)  &
          -k(585)*n(idx_H3j)
      pdj(110) =  &
          -k(486)*n(idx_HCO2j)
      pdj(112) =  &
          -k(488)*n(idx_N2Hj)
      pdj(113) =  &
          -k(489)*n(idx_O2Hj)
    elseif(j==26) then
      pdj(2) =  &
          -k(928)*n(idx_CH)
      pdj(3) =  &
          -k(1072)*n(idx_O)  &
          +k(846)*n(idx_OHj)
      pdj(5) =  &
          +k(885)*n(idx_CH2)  &
          +k(928)*n(idx_CH)
      pdj(6) =  &
          +k(593)*n(idx_H3j)
      pdj(7) =  &
          -k(866)*n(idx_C)
      pdj(8) =  &
          +k(522)*n(idx_H2j)
      pdj(11) =  &
          +k(734)*n(idx_O2Hj)
      pdj(12) =  &
          -k(885)*n(idx_CH2)
      pdj(17) =  &
          +k(733)*n(idx_HNOj)  &
          +k(1072)*n(idx_O)
      pdj(24) =  &
          +k(866)*n(idx_C)
      pdj(26) =  &
          -k(1156)  &
          -k(818)*n(idx_Oj)  &
          -k(268)  &
          -k(762)*n(idx_NHj)  &
          -k(689)*n(idx_HEj)  &
          -k(928)*n(idx_CH)  &
          -k(846)*n(idx_OHj)  &
          -k(1072)*n(idx_O)  &
          -k(145)*n(idx_HEj)  &
          -k(885)*n(idx_CH2)  &
          -k(734)*n(idx_O2Hj)  &
          -k(866)*n(idx_C)  &
          -k(1263)  &
          -k(593)*n(idx_H3j)  &
          -k(522)*n(idx_H2j)  &
          -k(733)*n(idx_HNOj)
      pdj(30) =  &
          +k(818)*n(idx_Oj)  &
          +k(1072)*n(idx_O)  &
          +k(928)*n(idx_CH)  &
          +2.d0*k(268)  &
          +k(866)*n(idx_C)  &
          +k(689)*n(idx_HEj)  &
          +2.d0*k(1156)  &
          +k(762)*n(idx_NHj)
      pdj(31) =  &
          +k(885)*n(idx_CH2)
      pdj(35) =  &
          +k(689)*n(idx_HEj)  &
          +k(145)*n(idx_HEj)
      pdj(58) =  &
          +k(1263)
      pdj(79) =  &
          +k(818)*n(idx_Oj)
      pdj(88) =  &
          +k(145)*n(idx_HEj)
      pdj(92) =  &
          -k(818)*n(idx_Oj)
      pdj(93) =  &
          -k(846)*n(idx_OHj)
      pdj(96) =  &
          +k(689)*n(idx_HEj)
      pdj(98) =  &
          -k(762)*n(idx_NHj)
      pdj(102) =  &
          -k(522)*n(idx_H2j)
      pdj(103) =  &
          -k(145)*n(idx_HEj)  &
          -k(689)*n(idx_HEj)
      pdj(104) =  &
          -k(733)*n(idx_HNOj)
      pdj(106) =  &
          -k(593)*n(idx_H3j)
      pdj(112) =  &
          +k(846)*n(idx_OHj)  &
          +k(734)*n(idx_O2Hj)  &
          +k(733)*n(idx_HNOj)  &
          +k(522)*n(idx_H2j)  &
          +k(762)*n(idx_NHj)  &
          +k(593)*n(idx_H3j)
      pdj(113) =  &
          -k(734)*n(idx_O2Hj)
    elseif(j==27) then
      pdj(1) =  &
          +k(270)  &
          +k(1158)
      pdj(2) =  &
          +k(869)*n(idx_C)
      pdj(3) =  &
          -k(1074)*n(idx_O)  &
          +k(1033)*n(idx_OH)  &
          +k(792)*n(idx_OHj)  &
          +k(213)*n(idx_Oj)  &
          -k(1073)*n(idx_O)
      pdj(4) =  &
          +k(868)*n(idx_C)  &
          +k(787)*n(idx_HCNHj)
      pdj(5) =  &
          +k(786)*n(idx_HCNHj)  &
          +k(867)*n(idx_C)
      pdj(6) =  &
          +k(690)*n(idx_HEj)  &
          +k(594)*n(idx_H3j)  &
          +k(984)*n(idx_H)  &
          -k(962)*n(idx_H2)  &
          +k(110)*n(idx_H2j)  &
          +k(410)*n(idx_CHj)
      pdj(7) =  &
          -k(868)*n(idx_C)  &
          -k(869)*n(idx_C)  &
          -k(867)*n(idx_C)
      pdj(8) =  &
          +k(691)*n(idx_HEj)  &
          +k(271)  &
          +k(868)*n(idx_C)  &
          +k(1031)*n(idx_NO)  &
          +k(1073)*n(idx_O)  &
          -k(984)*n(idx_H)  &
          +k(867)*n(idx_C)  &
          +k(962)*n(idx_H2)  &
          +k(85)*n(idx_Hj)  &
          +k(374)*n(idx_Cj)  &
          +k(1159)
      pdj(9) =  &
          +k(784)*n(idx_H3Oj)  &
          +k(1030)*n(idx_NO)  &
          +k(1032)*n(idx_OH)  &
          +k(186)*n(idx_H2Oj)
      pdj(10) =  &
          +k(782)*n(idx_H2Oj)  &
          +k(1031)*n(idx_NO)  &
          -k(1032)*n(idx_OH)  &
          -k(1033)*n(idx_OH)  &
          +k(1074)*n(idx_O)  &
          +k(189)*n(idx_OHj)
      pdj(11) =  &
          +k(791)*n(idx_O2Hj)  &
          +k(188)*n(idx_O2j)
      pdj(13) =  &
          +k(783)*n(idx_H3COj)
      pdj(14) =  &
          +k(781)*n(idx_H2COj)
      pdj(16) =  &
          +k(962)*n(idx_H2)  &
          +k(1033)*n(idx_OH)  &
          +k(1029)*n(idx_CH4)
      pdj(17) =  &
          -k(1031)*n(idx_NO)  &
          +k(789)*n(idx_HNOj)  &
          -k(1030)*n(idx_NO)
      pdj(24) =  &
          +k(184)*n(idx_CNj)  &
          +k(785)*n(idx_HCNj)
      pdj(25) =  &
          +k(788)*n(idx_HCOj)  &
          +k(185)*n(idx_COj)
      pdj(26) =  &
          +k(1030)*n(idx_NO)  &
          +k(790)*n(idx_N2Hj)  &
          +k(1031)*n(idx_NO)  &
          +k(187)*n(idx_N2j)
      pdj(27) =  &
          -k(270)  &
          -k(1031)*n(idx_NO)  &
          -k(782)*n(idx_H2Oj)  &
          -k(790)*n(idx_N2Hj)  &
          -k(781)*n(idx_H2COj)  &
          -k(908)*n(idx_CH3)  &
          -k(187)*n(idx_N2j)  &
          -k(777)*n(idx_NH2j)  &
          -k(185)*n(idx_COj)  &
          -k(867)*n(idx_C)  &
          -k(962)*n(idx_H2)  &
          -k(594)*n(idx_H3j)  &
          -k(792)*n(idx_OHj)  &
          -k(213)*n(idx_Oj)  &
          -k(189)*n(idx_OHj)  &
          -k(85)*n(idx_Hj)  &
          -k(788)*n(idx_HCOj)  &
          -k(1159)  &
          -k(1029)*n(idx_CH4)  &
          -k(783)*n(idx_H3COj)  &
          -k(789)*n(idx_HNOj)  &
          -k(186)*n(idx_H2Oj)  &
          -k(690)*n(idx_HEj)  &
          -k(784)*n(idx_H3Oj)  &
          -k(1290)  &
          -k(271)  &
          -k(1032)*n(idx_OH)  &
          -k(984)*n(idx_H)  &
          -k(374)*n(idx_Cj)  &
          -k(763)*n(idx_NHj)  &
          -k(165)*n(idx_Nj)  &
          -k(410)*n(idx_CHj)  &
          -k(868)*n(idx_C)  &
          -k(785)*n(idx_HCNj)  &
          -k(1074)*n(idx_O)  &
          -k(184)*n(idx_CNj)  &
          -k(791)*n(idx_O2Hj)  &
          -k(188)*n(idx_O2j)  &
          -k(1033)*n(idx_OH)  &
          -k(780)*n(idx_COj)  &
          -k(869)*n(idx_C)  &
          -k(1158)  &
          -k(1073)*n(idx_O)  &
          -k(787)*n(idx_HCNHj)  &
          -k(691)*n(idx_HEj)  &
          -k(1030)*n(idx_NO)  &
          -k(786)*n(idx_HCNHj)  &
          -k(110)*n(idx_H2j)
      pdj(28) =  &
          +k(1029)*n(idx_CH4)  &
          -k(908)*n(idx_CH3)
      pdj(29) =  &
          -k(1029)*n(idx_CH4)  &
          +k(908)*n(idx_CH3)
      pdj(30) =  &
          +k(763)*n(idx_NHj)  &
          +k(165)*n(idx_Nj)
      pdj(31) =  &
          +k(271)  &
          +k(1159)  &
          +k(777)*n(idx_NH2j)  &
          +k(780)*n(idx_COj)  &
          +k(908)*n(idx_CH3)  &
          +k(984)*n(idx_H)  &
          +k(1074)*n(idx_O)  &
          +k(869)*n(idx_C)  &
          +k(1032)*n(idx_OH)
      pdj(35) =  &
          +k(690)*n(idx_HEj)  &
          +k(691)*n(idx_HEj)
      pdj(36) =  &
          +k(1073)*n(idx_O)
      pdj(60) =  &
          +k(1290)
      pdj(70) =  &
          +k(780)*n(idx_COj)  &
          -k(788)*n(idx_HCOj)
      pdj(71) =  &
          -k(85)*n(idx_Hj)
      pdj(73) =  &
          -k(374)*n(idx_Cj)
      pdj(75) =  &
          -k(410)*n(idx_CHj)
      pdj(76) =  &
          -k(781)*n(idx_H2COj)
      pdj(78) =  &
          +k(594)*n(idx_H3j)  &
          +k(777)*n(idx_NH2j)  &
          +k(782)*n(idx_H2Oj)  &
          +k(785)*n(idx_HCNj)  &
          +k(789)*n(idx_HNOj)  &
          +k(786)*n(idx_HCNHj)  &
          +k(790)*n(idx_N2Hj)  &
          +k(791)*n(idx_O2Hj)  &
          +k(781)*n(idx_H2COj)  &
          +k(784)*n(idx_H3Oj)  &
          +k(783)*n(idx_H3COj)  &
          +k(787)*n(idx_HCNHj)  &
          +k(763)*n(idx_NHj)  &
          +k(792)*n(idx_OHj)  &
          +k(788)*n(idx_HCOj)
      pdj(86) =  &
          -k(184)*n(idx_CNj)
      pdj(87) =  &
          -k(780)*n(idx_COj)  &
          -k(185)*n(idx_COj)
      pdj(88) =  &
          -k(187)*n(idx_N2j)
      pdj(89) =  &
          -k(188)*n(idx_O2j)
      pdj(90) =  &
          -k(186)*n(idx_H2Oj)  &
          -k(782)*n(idx_H2Oj)
      pdj(91) =  &
          +k(270)  &
          +k(184)*n(idx_CNj)  &
          +k(186)*n(idx_H2Oj)  &
          +k(188)*n(idx_O2j)  &
          +k(213)*n(idx_Oj)  &
          +k(85)*n(idx_Hj)  &
          -k(777)*n(idx_NH2j)  &
          +k(110)*n(idx_H2j)  &
          +k(187)*n(idx_N2j)  &
          +k(165)*n(idx_Nj)  &
          +k(1158)  &
          +k(185)*n(idx_COj)  &
          +k(189)*n(idx_OHj)
      pdj(92) =  &
          -k(213)*n(idx_Oj)
      pdj(93) =  &
          -k(792)*n(idx_OHj)  &
          -k(189)*n(idx_OHj)
      pdj(96) =  &
          +k(690)*n(idx_HEj)  &
          -k(165)*n(idx_Nj)
      pdj(97) =  &
          +k(374)*n(idx_Cj)  &
          -k(785)*n(idx_HCNj)  &
          +k(410)*n(idx_CHj)
      pdj(98) =  &
          +k(691)*n(idx_HEj)  &
          -k(763)*n(idx_NHj)
      pdj(102) =  &
          -k(110)*n(idx_H2j)
      pdj(103) =  &
          -k(690)*n(idx_HEj)  &
          -k(691)*n(idx_HEj)
      pdj(104) =  &
          -k(789)*n(idx_HNOj)
      pdj(106) =  &
          -k(594)*n(idx_H3j)
      pdj(107) =  &
          -k(783)*n(idx_H3COj)
      pdj(108) =  &
          -k(784)*n(idx_H3Oj)
      pdj(109) =  &
          -k(787)*n(idx_HCNHj)  &
          -k(786)*n(idx_HCNHj)
      pdj(112) =  &
          -k(790)*n(idx_N2Hj)
      pdj(113) =  &
          -k(791)*n(idx_O2Hj)
    elseif(j==28) then
      pdj(1) =  &
          +k(246)  &
          +k(1118)
      pdj(2) =  &
          +k(1119)  &
          +k(247)
      pdj(3) =  &
          +k(918)*n(idx_OH)  &
          -k(916)*n(idx_O)  &
          -k(917)*n(idx_O)
      pdj(5) =  &
          +k(903)*n(idx_CN)  &
          +k(911)*n(idx_NO)  &
          +k(1011)*n(idx_N)  &
          +k(1010)*n(idx_N)
      pdj(6) =  &
          +k(579)*n(idx_H3j)  &
          +k(1119)  &
          +k(656)*n(idx_HEj)  &
          +k(1010)*n(idx_N)  &
          +k(916)*n(idx_O)  &
          -k(958)*n(idx_H2)  &
          +k(919)*n(idx_OH)  &
          +k(969)*n(idx_H)  &
          +k(247)
      pdj(8) =  &
          +k(917)*n(idx_O)  &
          +k(958)*n(idx_H2)  &
          +k(245)  &
          +k(916)*n(idx_O)  &
          +k(77)*n(idx_Hj)  &
          -k(969)*n(idx_H)  &
          +2.d0*k(1011)*n(idx_N)  &
          +k(1009)*n(idx_N)  &
          +k(1117)
      pdj(9) =  &
          +k(920)*n(idx_OH)  &
          -k(905)*n(idx_H2O)  &
          +k(911)*n(idx_NO)  &
          +k(913)*n(idx_O2)
      pdj(10) =  &
          -k(918)*n(idx_OH)  &
          +k(905)*n(idx_H2O)  &
          -k(920)*n(idx_OH)  &
          +k(912)*n(idx_O2)  &
          -k(919)*n(idx_OH)
      pdj(11) =  &
          -k(914)*n(idx_O2)  &
          -k(913)*n(idx_O2)  &
          -k(912)*n(idx_O2)  &
          +k(915)*n(idx_O2H)
      pdj(12) =  &
          +k(903)*n(idx_CN)  &
          +k(245)  &
          +2.d0*k(902)*n(idx_CH3)  &
          +k(969)*n(idx_H)  &
          +k(914)*n(idx_O2)  &
          +k(920)*n(idx_OH)  &
          +k(1117)
      pdj(13) =  &
          -k(904)*n(idx_H2CO)  &
          +k(917)*n(idx_O)  &
          +k(910)*n(idx_NO2)  &
          +k(912)*n(idx_O2)  &
          +k(919)*n(idx_OH)
      pdj(14) =  &
          -k(906)*n(idx_HCO)  &
          +k(913)*n(idx_O2)  &
          +k(904)*n(idx_H2CO)
      pdj(16) =  &
          -k(909)*n(idx_NH3)
      pdj(17) =  &
          +k(907)*n(idx_HNO)  &
          -k(911)*n(idx_NO)
      pdj(24) =  &
          -k(903)*n(idx_CN)
      pdj(25) =  &
          +k(906)*n(idx_HCO)  &
          +k(916)*n(idx_O)
      pdj(27) =  &
          -k(908)*n(idx_NH2)  &
          +k(909)*n(idx_NH3)
      pdj(28) =  &
          -k(656)*n(idx_HEj)  &
          -k(907)*n(idx_HNO)  &
          -k(920)*n(idx_OH)  &
          -k(1010)*n(idx_N)  &
          -k(917)*n(idx_O)  &
          -k(245)  &
          -4.d0*k(902)*n(idx_CH3)  &
          -k(905)*n(idx_H2O)  &
          -k(958)*n(idx_H2)  &
          -k(908)*n(idx_NH2)  &
          -k(77)*n(idx_Hj)  &
          -k(916)*n(idx_O)  &
          -k(912)*n(idx_O2)  &
          -k(1117)  &
          -k(906)*n(idx_HCO)  &
          -k(247)  &
          -k(903)*n(idx_CN)  &
          -k(913)*n(idx_O2)  &
          -k(969)*n(idx_H)  &
          -k(1009)*n(idx_N)  &
          -k(1119)  &
          -k(579)*n(idx_H3j)  &
          -k(1260)  &
          -k(918)*n(idx_OH)  &
          -k(910)*n(idx_NO2)  &
          -k(914)*n(idx_O2)  &
          -k(904)*n(idx_H2CO)  &
          -k(1011)*n(idx_N)  &
          -k(909)*n(idx_NH3)  &
          -k(919)*n(idx_OH)  &
          -k(911)*n(idx_NO)  &
          -k(915)*n(idx_O2H)  &
          -k(246)  &
          -k(1118)
      pdj(29) =  &
          +k(958)*n(idx_H2)  &
          +k(909)*n(idx_NH3)  &
          +k(904)*n(idx_H2CO)  &
          +k(908)*n(idx_NH2)  &
          +2.d0*k(902)*n(idx_CH3)  &
          +k(907)*n(idx_HNO)  &
          +k(915)*n(idx_O2H)  &
          +k(905)*n(idx_H2O)  &
          +k(918)*n(idx_OH)  &
          +k(906)*n(idx_HCO)
      pdj(30) =  &
          -k(1011)*n(idx_N)  &
          -k(1010)*n(idx_N)  &
          -k(1009)*n(idx_N)
      pdj(31) =  &
          +k(908)*n(idx_NH2)
      pdj(35) =  &
          +k(656)*n(idx_HEj)
      pdj(36) =  &
          -k(907)*n(idx_HNO)  &
          +k(910)*n(idx_NO2)
      pdj(39) =  &
          +k(1009)*n(idx_N)
      pdj(42) =  &
          -k(910)*n(idx_NO2)
      pdj(43) =  &
          +k(914)*n(idx_O2)  &
          -k(915)*n(idx_O2H)
      pdj(53) =  &
          +k(1260)
      pdj(71) =  &
          -k(77)*n(idx_Hj)
      pdj(75) =  &
          +k(656)*n(idx_HEj)
      pdj(94) =  &
          +k(246)  &
          +k(1118)  &
          +k(77)*n(idx_Hj)
      pdj(95) =  &
          +k(579)*n(idx_H3j)
      pdj(103) =  &
          -k(656)*n(idx_HEj)
      pdj(106) =  &
          -k(579)*n(idx_H3j)
    elseif(j==29) then
      pdj(1) =  &
          +k(1127)
      pdj(2) =  &
          +k(1128)
      pdj(3) =  &
          -k(1057)*n(idx_O)  &
          +k(208)*n(idx_Oj)
      pdj(5) =  &
          +k(921)*n(idx_CN)
      pdj(6) =  &
          +k(1125)  &
          +k(456)*n(idx_N2j)  &
          +k(718)*n(idx_Nj)  &
          +k(496)*n(idx_Hj)  &
          +k(970)*n(idx_H)  &
          +k(1128)  &
          +k(660)*n(idx_HEj)  &
          +k(250)  &
          +k(659)*n(idx_HEj)  &
          +k(512)*n(idx_H2j)  &
          +k(102)*n(idx_H2j)
      pdj(8) =  &
          +k(457)*n(idx_N2j)  &
          +2.d0*k(719)*n(idx_Nj)  &
          +k(717)*n(idx_Nj)  &
          -k(970)*n(idx_H)  &
          +k(512)*n(idx_H2j)  &
          +k(78)*n(idx_Hj)  &
          +k(1128)  &
          +k(661)*n(idx_HEj)  &
          +k(718)*n(idx_Nj)  &
          +k(659)*n(idx_HEj)  &
          +k(1126)
      pdj(9) =  &
          +k(923)*n(idx_OH)
      pdj(10) =  &
          +k(811)*n(idx_Oj)  &
          -k(923)*n(idx_OH)  &
          +k(1057)*n(idx_O)
      pdj(11) =  &
          -k(922)*n(idx_O2)
      pdj(12) =  &
          +k(1125)  &
          -k(880)*n(idx_CH2)  &
          +k(250)  &
          +k(458)*n(idx_OHj)
      pdj(16) =  &
          +k(1029)*n(idx_NH2)
      pdj(24) =  &
          -k(921)*n(idx_CN)
      pdj(25) =  &
          +k(53)*n(idx_COj)
      pdj(26) =  &
          +k(457)*n(idx_N2j)  &
          +k(456)*n(idx_N2j)
      pdj(27) =  &
          -k(1029)*n(idx_NH2)  &
          +k(1035)*n(idx_NH)
      pdj(28) =  &
          +k(662)*n(idx_HEj)  &
          +k(453)*n(idx_H2COj)  &
          +k(454)*n(idx_H2Oj)  &
          +k(923)*n(idx_OH)  &
          +2.d0*k(880)*n(idx_CH2)  &
          +k(921)*n(idx_CN)  &
          +k(452)*n(idx_COj)  &
          +k(1029)*n(idx_NH2)  &
          +k(455)*n(idx_HCNj)  &
          +k(970)*n(idx_H)  &
          +k(1126)  &
          +k(1035)*n(idx_NH)  &
          +k(1057)*n(idx_O)  &
          +k(922)*n(idx_O2)
      pdj(29) =  &
          -k(496)*n(idx_Hj)  &
          -k(923)*n(idx_OH)  &
          -k(1261)  &
          -k(455)*n(idx_HCNj)  &
          -k(78)*n(idx_Hj)  &
          -k(208)*n(idx_Oj)  &
          -k(661)*n(idx_HEj)  &
          -k(454)*n(idx_H2Oj)  &
          -k(157)*n(idx_Nj)  &
          -k(141)*n(idx_HEj)  &
          -k(512)*n(idx_H2j)  &
          -k(921)*n(idx_CN)  &
          -k(458)*n(idx_OHj)  &
          -k(718)*n(idx_Nj)  &
          -k(659)*n(idx_HEj)  &
          -k(453)*n(idx_H2COj)  &
          -k(1057)*n(idx_O)  &
          -k(1126)  &
          -k(1029)*n(idx_NH2)  &
          -k(102)*n(idx_H2j)  &
          -k(660)*n(idx_HEj)  &
          -k(1035)*n(idx_NH)  &
          -k(1128)  &
          -k(880)*n(idx_CH2)  &
          -k(719)*n(idx_Nj)  &
          -k(452)*n(idx_COj)  &
          -k(457)*n(idx_N2j)  &
          -k(662)*n(idx_HEj)  &
          -k(970)*n(idx_H)  &
          -k(456)*n(idx_N2j)  &
          -k(250)  &
          -k(717)*n(idx_Nj)  &
          -k(1125)  &
          -k(53)*n(idx_COj)  &
          -k(811)*n(idx_Oj)  &
          -k(1127)  &
          -k(922)*n(idx_O2)
      pdj(30) =  &
          +k(157)*n(idx_Nj)  &
          +k(717)*n(idx_Nj)
      pdj(31) =  &
          -k(1035)*n(idx_NH)
      pdj(35) =  &
          +k(661)*n(idx_HEj)  &
          +k(141)*n(idx_HEj)  &
          +k(662)*n(idx_HEj)  &
          +k(659)*n(idx_HEj)  &
          +k(660)*n(idx_HEj)
      pdj(43) =  &
          +k(922)*n(idx_O2)
      pdj(53) =  &
          +k(1261)
      pdj(70) =  &
          +k(452)*n(idx_COj)
      pdj(71) =  &
          -k(78)*n(idx_Hj)  &
          +k(662)*n(idx_HEj)  &
          -k(496)*n(idx_Hj)
      pdj(74) =  &
          +k(660)*n(idx_HEj)  &
          +k(456)*n(idx_N2j)
      pdj(75) =  &
          +k(659)*n(idx_HEj)
      pdj(76) =  &
          -k(453)*n(idx_H2COj)
      pdj(87) =  &
          -k(452)*n(idx_COj)  &
          -k(53)*n(idx_COj)
      pdj(88) =  &
          -k(457)*n(idx_N2j)  &
          -k(456)*n(idx_N2j)
      pdj(90) =  &
          -k(454)*n(idx_H2Oj)
      pdj(92) =  &
          -k(208)*n(idx_Oj)  &
          -k(811)*n(idx_Oj)
      pdj(93) =  &
          -k(458)*n(idx_OHj)
      pdj(94) =  &
          +k(457)*n(idx_N2j)  &
          +k(811)*n(idx_Oj)  &
          +k(717)*n(idx_Nj)  &
          +k(496)*n(idx_Hj)  &
          +k(661)*n(idx_HEj)  &
          +k(512)*n(idx_H2j)
      pdj(95) =  &
          +k(157)*n(idx_Nj)  &
          +k(208)*n(idx_Oj)  &
          +k(78)*n(idx_Hj)  &
          +k(102)*n(idx_H2j)  &
          +k(141)*n(idx_HEj)  &
          +k(1127)  &
          +k(53)*n(idx_COj)
      pdj(96) =  &
          -k(718)*n(idx_Nj)  &
          -k(717)*n(idx_Nj)  &
          -k(719)*n(idx_Nj)  &
          -k(157)*n(idx_Nj)
      pdj(97) =  &
          -k(455)*n(idx_HCNj)  &
          +k(718)*n(idx_Nj)
      pdj(102) =  &
          -k(512)*n(idx_H2j)  &
          -k(102)*n(idx_H2j)
      pdj(103) =  &
          -k(660)*n(idx_HEj)  &
          -k(659)*n(idx_HEj)  &
          -k(662)*n(idx_HEj)  &
          -k(661)*n(idx_HEj)  &
          -k(141)*n(idx_HEj)
      pdj(107) =  &
          +k(453)*n(idx_H2COj)
      pdj(108) =  &
          +k(458)*n(idx_OHj)  &
          +k(454)*n(idx_H2Oj)
      pdj(109) =  &
          +k(719)*n(idx_Nj)  &
          +k(455)*n(idx_HCNj)
    elseif(j==30) then
      pdj(1) =  &
          +k(239)  &
          +k(269)
      pdj(2) =  &
          -k(929)*n(idx_CH)  &
          -k(930)*n(idx_CH)  &
          +k(1008)*n(idx_CH2)
      pdj(3) =  &
          +k(1023)*n(idx_NO)  &
          +k(1016)*n(idx_HCO)  &
          +k(1027)*n(idx_OH)  &
          +k(1024)*n(idx_O2)  &
          +2.d0*k(1020)*n(idx_NO2)  &
          +k(743)*n(idx_O2j)
      pdj(4) =  &
          +k(1007)*n(idx_CH2)
      pdj(5) =  &
          +k(1016)*n(idx_HCO)  &
          +k(1014)*n(idx_H2CN)  &
          +k(1010)*n(idx_CH3)  &
          +k(1006)*n(idx_CH2)  &
          +k(1011)*n(idx_CH3)
      pdj(6) =  &
          +k(1010)*n(idx_CH3)  &
          -k(961)*n(idx_H2)  &
          +k(740)*n(idx_H2Oj)
      pdj(7) =  &
          +k(930)*n(idx_CH)  &
          +k(1012)*n(idx_CN)  &
          +k(738)*n(idx_CNj)  &
          -k(1195)*n(idx_C)
      pdj(8) =  &
          +k(742)*n(idx_NH2j)  &
          +2.d0*k(1011)*n(idx_CH3)  &
          +k(744)*n(idx_OHj)  &
          +k(739)*n(idx_H2Oj)  &
          +k(1009)*n(idx_CH3)  &
          +k(929)*n(idx_CH)  &
          +k(741)*n(idx_NHj)  &
          +k(1019)*n(idx_NH)  &
          +k(523)*n(idx_H2j)  &
          +k(961)*n(idx_H2)  &
          +k(1007)*n(idx_CH2)  &
          +k(1026)*n(idx_OH)  &
          +k(737)*n(idx_CH2j)  &
          +k(409)*n(idx_CHj)  &
          +k(1006)*n(idx_CH2)  &
          +k(1017)*n(idx_HCO)
      pdj(10) =  &
          -k(1027)*n(idx_OH)  &
          -k(1026)*n(idx_OH)
      pdj(11) =  &
          -k(1024)*n(idx_O2)  &
          +k(1025)*n(idx_O2H)  &
          +k(1022)*n(idx_NO2)
      pdj(12) =  &
          -k(1008)*n(idx_CH2)  &
          -k(1007)*n(idx_CH2)  &
          -k(1006)*n(idx_CH2)
      pdj(14) =  &
          -k(1017)*n(idx_HCO)  &
          -k(1016)*n(idx_HCO)  &
          -k(1015)*n(idx_HCO)
      pdj(17) =  &
          +2.d0*k(1021)*n(idx_NO2)  &
          +k(1013)*n(idx_CO2)  &
          +k(1024)*n(idx_O2)  &
          +k(1026)*n(idx_OH)  &
          +k(747)*n(idx_SIOj)  &
          -k(1023)*n(idx_NO)  &
          +k(1018)*n(idx_HNO)
      pdj(18) =  &
          +k(746)*n(idx_SIOj)  &
          +k(1028)*n(idx_SIC)
      pdj(21) =  &
          -k(1028)*n(idx_SIC)
      pdj(24) =  &
          +k(1195)*n(idx_C)  &
          -k(1012)*n(idx_CN)  &
          +k(745)*n(idx_SICj)  &
          +k(1028)*n(idx_SIC)  &
          +k(929)*n(idx_CH)
      pdj(25) =  &
          +k(1013)*n(idx_CO2)  &
          +k(1015)*n(idx_HCO)
      pdj(26) =  &
          +k(1023)*n(idx_NO)  &
          +k(1019)*n(idx_NH)  &
          +k(175)*n(idx_N2j)  &
          +k(1022)*n(idx_NO2)  &
          +k(1012)*n(idx_CN)  &
          +k(1020)*n(idx_NO2)
      pdj(28) =  &
          -k(1010)*n(idx_CH3)  &
          -k(1009)*n(idx_CH3)  &
          -k(1011)*n(idx_CH3)
      pdj(30) =  &
          -k(1020)*n(idx_NO2)  &
          -k(1026)*n(idx_OH)  &
          -k(1013)*n(idx_CO2)  &
          -k(745)*n(idx_SICj)  &
          -k(742)*n(idx_NH2j)  &
          -k(929)*n(idx_CH)  &
          -k(1018)*n(idx_HNO)  &
          -k(1211)*n(idx_Nj)  &
          -k(1195)*n(idx_C)  &
          -k(523)*n(idx_H2j)  &
          -k(961)*n(idx_H2)  &
          -k(741)*n(idx_NHj)  &
          -k(1027)*n(idx_OH)  &
          -k(1007)*n(idx_CH2)  &
          -k(1014)*n(idx_H2CN)  &
          -k(409)*n(idx_CHj)  &
          -k(743)*n(idx_O2j)  &
          -k(1017)*n(idx_HCO)  &
          -k(1011)*n(idx_CH3)  &
          -k(740)*n(idx_H2Oj)  &
          -k(1019)*n(idx_NH)  &
          -k(738)*n(idx_CNj)  &
          -k(744)*n(idx_OHj)  &
          -k(1021)*n(idx_NO2)  &
          -k(1006)*n(idx_CH2)  &
          -k(1291)  &
          -k(1010)*n(idx_CH3)  &
          -k(1028)*n(idx_SIC)  &
          -k(1022)*n(idx_NO2)  &
          -k(1012)*n(idx_CN)  &
          -k(747)*n(idx_SIOj)  &
          -k(1024)*n(idx_O2)  &
          -k(1193)*n(idx_Cj)  &
          -k(930)*n(idx_CH)  &
          -k(746)*n(idx_SIOj)  &
          -k(1008)*n(idx_CH2)  &
          -k(269)  &
          -k(175)*n(idx_N2j)  &
          -k(739)*n(idx_H2Oj)  &
          -k(1025)*n(idx_O2H)  &
          -k(737)*n(idx_CH2j)  &
          -k(239)  &
          -k(1009)*n(idx_CH3)  &
          -k(1016)*n(idx_HCO)  &
          -k(1015)*n(idx_HCO)  &
          -k(1023)*n(idx_NO)
      pdj(31) =  &
          +k(1015)*n(idx_HCO)  &
          +k(1027)*n(idx_OH)  &
          +k(1008)*n(idx_CH2)  &
          +k(1025)*n(idx_O2H)  &
          +k(930)*n(idx_CH)  &
          +k(1014)*n(idx_H2CN)  &
          +k(961)*n(idx_H2)  &
          +k(1018)*n(idx_HNO)  &
          -k(1019)*n(idx_NH)
      pdj(36) =  &
          -k(1018)*n(idx_HNO)
      pdj(38) =  &
          -k(1013)*n(idx_CO2)
      pdj(39) =  &
          -k(1014)*n(idx_H2CN)  &
          +k(1009)*n(idx_CH3)
      pdj(42) =  &
          -k(1022)*n(idx_NO2)  &
          -k(1020)*n(idx_NO2)  &
          -k(1021)*n(idx_NO2)
      pdj(43) =  &
          -k(1025)*n(idx_O2H)
      pdj(44) =  &
          +k(1017)*n(idx_HCO)
      pdj(60) =  &
          +k(1291)
      pdj(73) =  &
          -k(1193)*n(idx_Cj)
      pdj(74) =  &
          -k(737)*n(idx_CH2j)
      pdj(75) =  &
          -k(409)*n(idx_CHj)
      pdj(79) =  &
          +k(746)*n(idx_SIOj)  &
          +k(743)*n(idx_O2j)  &
          +k(744)*n(idx_OHj)  &
          +k(740)*n(idx_H2Oj)
      pdj(80) =  &
          +k(747)*n(idx_SIOj)  &
          +k(745)*n(idx_SICj)
      pdj(83) =  &
          -k(745)*n(idx_SICj)
      pdj(86) =  &
          -k(738)*n(idx_CNj)  &
          +k(409)*n(idx_CHj)  &
          +k(1193)*n(idx_Cj)
      pdj(88) =  &
          +k(741)*n(idx_NHj)  &
          -k(175)*n(idx_N2j)  &
          +k(738)*n(idx_CNj)  &
          +k(1211)*n(idx_Nj)
      pdj(89) =  &
          -k(743)*n(idx_O2j)
      pdj(90) =  &
          -k(739)*n(idx_H2Oj)  &
          -k(740)*n(idx_H2Oj)
      pdj(91) =  &
          -k(742)*n(idx_NH2j)
      pdj(93) =  &
          -k(744)*n(idx_OHj)
      pdj(96) =  &
          +k(175)*n(idx_N2j)  &
          -k(1211)*n(idx_Nj)  &
          +k(239)  &
          +k(269)
      pdj(97) =  &
          +k(737)*n(idx_CH2j)
      pdj(98) =  &
          +k(523)*n(idx_H2j)  &
          -k(741)*n(idx_NHj)
      pdj(101) =  &
          -k(747)*n(idx_SIOj)  &
          -k(746)*n(idx_SIOj)
      pdj(102) =  &
          -k(523)*n(idx_H2j)
      pdj(104) =  &
          +k(739)*n(idx_H2Oj)
      pdj(112) =  &
          +k(742)*n(idx_NH2j)
    elseif(j==31) then
      pdj(1) =  &
          +k(1164)  &
          +k(276)
      pdj(2) =  &
          +k(871)*n(idx_C)
      pdj(3) =  &
          -k(1047)*n(idx_O)  &
          +k(1043)*n(idx_NO)  &
          +k(807)*n(idx_OHj)  &
          +k(805)*n(idx_O2j)  &
          +k(1051)*n(idx_OH)  &
          +k(203)*n(idx_Oj)  &
          -k(1048)*n(idx_O)  &
          +k(1045)*n(idx_O2)
      pdj(5) =  &
          +k(1036)*n(idx_CN)
      pdj(6) =  &
          +k(795)*n(idx_CH3j)  &
          +k(986)*n(idx_H)  &
          +k(595)*n(idx_H3j)  &
          +k(411)*n(idx_CHj)  &
          +k(112)*n(idx_H2j)  &
          +2.d0*k(1039)*n(idx_NH)  &
          -k(963)*n(idx_H2)
      pdj(7) =  &
          -k(871)*n(idx_C)  &
          -k(870)*n(idx_C)
      pdj(8) =  &
          +k(1047)*n(idx_O)  &
          +k(727)*n(idx_Nj)  &
          +k(1050)*n(idx_OH)  &
          +k(963)*n(idx_H2)  &
          +k(1163)  &
          +k(376)*n(idx_Cj)  &
          +k(1043)*n(idx_NO)  &
          +k(694)*n(idx_HEj)  &
          +k(275)  &
          +4.d0*k(1040)*n(idx_NH)  &
          +k(870)*n(idx_C)  &
          +k(87)*n(idx_Hj)  &
          -k(986)*n(idx_H)  &
          +k(804)*n(idx_Oj)  &
          +k(524)*n(idx_H2j)  &
          +k(1019)*n(idx_N)
      pdj(9) =  &
          -k(1037)*n(idx_H2O)  &
          +k(1049)*n(idx_OH)
      pdj(10) =  &
          +k(1046)*n(idx_O2)  &
          +k(1037)*n(idx_H2O)  &
          -k(1049)*n(idx_OH)  &
          +k(1044)*n(idx_NO)  &
          +k(1048)*n(idx_O)  &
          -k(1050)*n(idx_OH)  &
          -k(1051)*n(idx_OH)
      pdj(11) =  &
          -k(1045)*n(idx_O2)  &
          -k(1046)*n(idx_O2)  &
          +k(806)*n(idx_O2Hj)
      pdj(16) =  &
          -k(1038)*n(idx_NH3)
      pdj(17) =  &
          +k(1047)*n(idx_O)  &
          +k(1046)*n(idx_O2)  &
          -k(1044)*n(idx_NO)  &
          +k(801)*n(idx_HNOj)  &
          -k(1043)*n(idx_NO)  &
          +k(1042)*n(idx_NO2)
      pdj(24) =  &
          +k(870)*n(idx_C)  &
          +k(200)*n(idx_CNj)  &
          -k(1036)*n(idx_CN)  &
          +k(799)*n(idx_HCNj)
      pdj(25) =  &
          +k(201)*n(idx_COj)  &
          +k(800)*n(idx_HCOj)
      pdj(26) =  &
          +k(202)*n(idx_N2j)  &
          +k(802)*n(idx_N2Hj)  &
          +k(1043)*n(idx_NO)  &
          +2.d0*k(1040)*n(idx_NH)  &
          +k(1044)*n(idx_NO)  &
          +2.d0*k(1039)*n(idx_NH)  &
          +k(1019)*n(idx_N)
      pdj(27) =  &
          +k(963)*n(idx_H2)  &
          +k(1035)*n(idx_CH4)  &
          +2.d0*k(1038)*n(idx_NH3)  &
          +2.d0*k(1041)*n(idx_NH)  &
          +k(1037)*n(idx_H2O)  &
          +k(1051)*n(idx_OH)
      pdj(28) =  &
          +k(1035)*n(idx_CH4)
      pdj(29) =  &
          -k(1035)*n(idx_CH4)
      pdj(30) =  &
          +k(986)*n(idx_H)  &
          +k(1163)  &
          +k(167)*n(idx_Nj)  &
          +2.d0*k(1041)*n(idx_NH)  &
          +k(275)  &
          +k(1048)*n(idx_O)  &
          +k(764)*n(idx_NHj)  &
          +k(871)*n(idx_C)  &
          +k(803)*n(idx_NH2j)  &
          +k(1036)*n(idx_CN)  &
          +k(796)*n(idx_COj)  &
          +k(797)*n(idx_H2COj)  &
          -k(1019)*n(idx_N)  &
          +k(1049)*n(idx_OH)  &
          +k(798)*n(idx_H2Oj)
      pdj(31) =  &
          -k(806)*n(idx_O2Hj)  &
          -k(871)*n(idx_C)  &
          -k(1038)*n(idx_NH3)  &
          -k(799)*n(idx_HCNj)  &
          -k(801)*n(idx_HNOj)  &
          -k(1043)*n(idx_NO)  &
          -k(1046)*n(idx_O2)  &
          -4.d0*k(1039)*n(idx_NH)  &
          -k(800)*n(idx_HCOj)  &
          -k(1037)*n(idx_H2O)  &
          -4.d0*k(1040)*n(idx_NH)  &
          -k(1035)*n(idx_CH4)  &
          -k(87)*n(idx_Hj)  &
          -k(595)*n(idx_H3j)  &
          -k(1042)*n(idx_NO2)  &
          -k(1045)*n(idx_O2)  &
          -k(798)*n(idx_H2Oj)  &
          -k(276)  &
          -k(1163)  &
          -k(1051)*n(idx_OH)  &
          -k(200)*n(idx_CNj)  &
          -k(1019)*n(idx_N)  &
          -k(1048)*n(idx_O)  &
          -4.d0*k(1041)*n(idx_NH)  &
          -k(112)*n(idx_H2j)  &
          -k(795)*n(idx_CH3j)  &
          -k(986)*n(idx_H)  &
          -k(1266)  &
          -k(727)*n(idx_Nj)  &
          -k(805)*n(idx_O2j)  &
          -k(1164)  &
          -k(803)*n(idx_NH2j)  &
          -k(376)*n(idx_Cj)  &
          -k(1050)*n(idx_OH)  &
          -k(524)*n(idx_H2j)  &
          -k(796)*n(idx_COj)  &
          -k(201)*n(idx_COj)  &
          -k(411)*n(idx_CHj)  &
          -k(1036)*n(idx_CN)  &
          -k(804)*n(idx_Oj)  &
          -k(1044)*n(idx_NO)  &
          -k(1047)*n(idx_O)  &
          -k(203)*n(idx_Oj)  &
          -k(764)*n(idx_NHj)  &
          -k(202)*n(idx_N2j)  &
          -k(275)  &
          -k(870)*n(idx_C)  &
          -k(1049)*n(idx_OH)  &
          -k(167)*n(idx_Nj)  &
          -k(694)*n(idx_HEj)  &
          -k(807)*n(idx_OHj)  &
          -k(802)*n(idx_N2Hj)  &
          -k(963)*n(idx_H2)  &
          -k(797)*n(idx_H2COj)
      pdj(35) =  &
          +k(694)*n(idx_HEj)
      pdj(36) =  &
          +k(1042)*n(idx_NO2)  &
          +k(1050)*n(idx_OH)  &
          +k(1045)*n(idx_O2)
      pdj(42) =  &
          -k(1042)*n(idx_NO2)
      pdj(60) =  &
          +k(1266)
      pdj(70) =  &
          +k(796)*n(idx_COj)  &
          -k(800)*n(idx_HCOj)
      pdj(71) =  &
          -k(87)*n(idx_Hj)
      pdj(73) =  &
          -k(376)*n(idx_Cj)
      pdj(75) =  &
          -k(411)*n(idx_CHj)
      pdj(76) =  &
          -k(797)*n(idx_H2COj)
      pdj(78) =  &
          +k(803)*n(idx_NH2j)
      pdj(79) =  &
          +k(804)*n(idx_Oj)
      pdj(86) =  &
          +k(376)*n(idx_Cj)  &
          -k(200)*n(idx_CNj)  &
          +k(411)*n(idx_CHj)
      pdj(87) =  &
          -k(796)*n(idx_COj)  &
          -k(201)*n(idx_COj)
      pdj(88) =  &
          +k(727)*n(idx_Nj)  &
          -k(202)*n(idx_N2j)
      pdj(89) =  &
          -k(805)*n(idx_O2j)
      pdj(90) =  &
          -k(798)*n(idx_H2Oj)
      pdj(91) =  &
          +k(802)*n(idx_N2Hj)  &
          +k(800)*n(idx_HCOj)  &
          +k(595)*n(idx_H3j)  &
          +k(801)*n(idx_HNOj)  &
          +k(807)*n(idx_OHj)  &
          +k(799)*n(idx_HCNj)  &
          +k(764)*n(idx_NHj)  &
          -k(803)*n(idx_NH2j)  &
          +k(806)*n(idx_O2Hj)  &
          +k(524)*n(idx_H2j)
      pdj(92) =  &
          -k(804)*n(idx_Oj)  &
          -k(203)*n(idx_Oj)
      pdj(93) =  &
          -k(807)*n(idx_OHj)
      pdj(94) =  &
          -k(795)*n(idx_CH3j)
      pdj(96) =  &
          +k(694)*n(idx_HEj)  &
          -k(167)*n(idx_Nj)  &
          -k(727)*n(idx_Nj)
      pdj(97) =  &
          -k(799)*n(idx_HCNj)
      pdj(98) =  &
          +k(202)*n(idx_N2j)  &
          +k(201)*n(idx_COj)  &
          +k(276)  &
          +k(167)*n(idx_Nj)  &
          -k(764)*n(idx_NHj)  &
          +k(112)*n(idx_H2j)  &
          +k(87)*n(idx_Hj)  &
          +k(203)*n(idx_Oj)  &
          +k(200)*n(idx_CNj)  &
          +k(1164)
      pdj(102) =  &
          -k(524)*n(idx_H2j)  &
          -k(112)*n(idx_H2j)
      pdj(103) =  &
          -k(694)*n(idx_HEj)
      pdj(104) =  &
          +k(805)*n(idx_O2j)  &
          -k(801)*n(idx_HNOj)
      pdj(106) =  &
          -k(595)*n(idx_H3j)
      pdj(107) =  &
          +k(797)*n(idx_H2COj)
      pdj(108) =  &
          +k(798)*n(idx_H2Oj)
      pdj(109) =  &
          +k(795)*n(idx_CH3j)
      pdj(112) =  &
          -k(802)*n(idx_N2Hj)
      pdj(113) =  &
          -k(806)*n(idx_O2Hj)
    elseif(j==32) then
      pdj(3) =  &
          -k(1089)*n(idx_O)
      pdj(5) =  &
          +k(951)*n(idx_CN)
      pdj(6) =  &
          +k(1186)  &
          +k(1188)  &
          +k(508)*n(idx_Hj)  &
          +k(709)*n(idx_HEj)  &
          +2.d0*k(708)*n(idx_HEj)  &
          +k(292)  &
          +k(605)*n(idx_H3j)
      pdj(8) =  &
          +k(98)*n(idx_Hj)  &
          +k(709)*n(idx_HEj)  &
          +k(1187)  &
          +k(1188)
      pdj(10) =  &
          +k(1089)*n(idx_O)
      pdj(22) =  &
          +k(1186)  &
          +k(292)
      pdj(23) =  &
          +k(1187)  &
          +k(1089)*n(idx_O)  &
          +k(951)*n(idx_CN)
      pdj(24) =  &
          -k(951)*n(idx_CN)
      pdj(25) =  &
          +k(639)*n(idx_HCOj)
      pdj(29) =  &
          +k(447)*n(idx_CH3j)
      pdj(32) =  &
          -k(951)*n(idx_CN)  &
          -k(1187)  &
          -k(98)*n(idx_Hj)  &
          -k(708)*n(idx_HEj)  &
          -k(292)  &
          -k(1188)  &
          -k(1238)  &
          -k(1186)  &
          -k(447)*n(idx_CH3j)  &
          -k(709)*n(idx_HEj)  &
          -k(508)*n(idx_Hj)  &
          -k(639)*n(idx_HCOj)  &
          -k(1089)*n(idx_O)  &
          -k(605)*n(idx_H3j)
      pdj(33) =  &
          +k(1188)
      pdj(35) =  &
          +k(708)*n(idx_HEj)  &
          +k(709)*n(idx_HEj)
      pdj(48) =  &
          +k(1238)
      pdj(70) =  &
          -k(639)*n(idx_HCOj)
      pdj(71) =  &
          -k(508)*n(idx_Hj)  &
          -k(98)*n(idx_Hj)
      pdj(80) =  &
          +k(708)*n(idx_HEj)
      pdj(85) =  &
          +k(508)*n(idx_Hj)  &
          +k(447)*n(idx_CH3j)
      pdj(94) =  &
          -k(447)*n(idx_CH3j)
      pdj(99) =  &
          +k(98)*n(idx_Hj)
      pdj(100) =  &
          +k(709)*n(idx_HEj)
      pdj(103) =  &
          -k(708)*n(idx_HEj)  &
          -k(709)*n(idx_HEj)
      pdj(106) =  &
          -k(605)*n(idx_H3j)
      pdj(114) =  &
          +k(605)*n(idx_H3j)  &
          +k(639)*n(idx_HCOj)
    elseif(j==33) then
      pdj(3) =  &
          +k(850)*n(idx_OHj)  &
          -k(1090)*n(idx_O)
      pdj(6) =  &
          +k(509)*n(idx_Hj)  &
          +k(606)*n(idx_H3j)
      pdj(7) =  &
          -k(878)*n(idx_C)
      pdj(8) =  &
          +k(878)*n(idx_C)  &
          +k(382)*n(idx_Cj)  &
          +k(1189)  &
          +k(1090)*n(idx_O)  &
          +k(710)*n(idx_HEj)  &
          +k(99)*n(idx_Hj)  &
          +k(293)
      pdj(9) =  &
          +k(613)*n(idx_H3Oj)
      pdj(18) =  &
          +k(1189)  &
          +k(293)
      pdj(21) =  &
          +k(878)*n(idx_C)
      pdj(25) =  &
          +k(640)*n(idx_HCOj)
      pdj(33) =  &
          -k(1231)  &
          -k(509)*n(idx_Hj)  &
          -k(293)  &
          -k(1189)  &
          -k(710)*n(idx_HEj)  &
          -k(99)*n(idx_Hj)  &
          -k(878)*n(idx_C)  &
          -k(640)*n(idx_HCOj)  &
          -k(613)*n(idx_H3Oj)  &
          -k(1090)*n(idx_O)  &
          -k(382)*n(idx_Cj)  &
          -k(850)*n(idx_OHj)  &
          -k(606)*n(idx_H3j)
      pdj(34) =  &
          +k(1090)*n(idx_O)
      pdj(35) =  &
          +k(710)*n(idx_HEj)
      pdj(48) =  &
          +k(1231)
      pdj(70) =  &
          -k(640)*n(idx_HCOj)
      pdj(71) =  &
          -k(99)*n(idx_Hj)  &
          -k(509)*n(idx_Hj)
      pdj(73) =  &
          -k(382)*n(idx_Cj)
      pdj(80) =  &
          +k(509)*n(idx_Hj)  &
          +k(710)*n(idx_HEj)
      pdj(83) =  &
          +k(382)*n(idx_Cj)
      pdj(84) =  &
          +k(640)*n(idx_HCOj)  &
          +k(850)*n(idx_OHj)  &
          +k(606)*n(idx_H3j)  &
          +k(613)*n(idx_H3Oj)
      pdj(93) =  &
          -k(850)*n(idx_OHj)
      pdj(100) =  &
          +k(99)*n(idx_Hj)
      pdj(103) =  &
          -k(710)*n(idx_HEj)
      pdj(106) =  &
          -k(606)*n(idx_H3j)
      pdj(108) =  &
          -k(613)*n(idx_H3Oj)
    elseif(j==34) then
      pdj(1) =  &
          +k(1192)
      pdj(3) =  &
          +k(1191)  &
          +k(294)  &
          +k(711)*n(idx_HEj)  &
          +k(851)*n(idx_OHj)
      pdj(6) =  &
          +k(607)*n(idx_H3j)
      pdj(8) =  &
          +k(100)*n(idx_Hj)
      pdj(9) =  &
          +k(614)*n(idx_H3Oj)
      pdj(18) =  &
          +k(1191)  &
          +k(294)  &
          +k(712)*n(idx_HEj)
      pdj(25) =  &
          +k(641)*n(idx_HCOj)  &
          +k(383)*n(idx_Cj)
      pdj(34) =  &
          -k(1192)  &
          -k(294)  &
          -k(1230)  &
          -k(711)*n(idx_HEj)  &
          -k(383)*n(idx_Cj)  &
          -k(641)*n(idx_HCOj)  &
          -k(851)*n(idx_OHj)  &
          -k(607)*n(idx_H3j)  &
          -k(1191)  &
          -k(614)*n(idx_H3Oj)  &
          -k(100)*n(idx_Hj)  &
          -k(712)*n(idx_HEj)
      pdj(35) =  &
          +k(711)*n(idx_HEj)  &
          +k(712)*n(idx_HEj)
      pdj(49) =  &
          +k(1230)
      pdj(70) =  &
          -k(641)*n(idx_HCOj)
      pdj(71) =  &
          -k(100)*n(idx_Hj)
      pdj(73) =  &
          -k(383)*n(idx_Cj)
      pdj(80) =  &
          +k(383)*n(idx_Cj)  &
          +k(711)*n(idx_HEj)
      pdj(92) =  &
          +k(712)*n(idx_HEj)
      pdj(93) =  &
          -k(851)*n(idx_OHj)
      pdj(101) =  &
          +k(1192)  &
          +k(100)*n(idx_Hj)
      pdj(103) =  &
          -k(712)*n(idx_HEj)  &
          -k(711)*n(idx_HEj)
      pdj(106) =  &
          -k(607)*n(idx_H3j)
      pdj(108) =  &
          -k(614)*n(idx_H3Oj)
      pdj(115) =  &
          +k(641)*n(idx_HCOj)  &
          +k(607)*n(idx_H3j)  &
          +k(851)*n(idx_OHj)  &
          +k(614)*n(idx_H3Oj)
    elseif(j==35) then
      pdj(1) =  &
          +k(266)  &
          +k(238)
      pdj(8) =  &
          +k(521)*n(idx_H2j)
      pdj(35) =  &
          -k(266)  &
          -k(1199)*n(idx_Hj)  &
          -k(521)*n(idx_H2j)  &
          -k(238)
      pdj(71) =  &
          -k(1199)*n(idx_Hj)
      pdj(102) =  &
          -k(521)*n(idx_H2j)
      pdj(103) =  &
          +k(266)  &
          +k(238)
      pdj(111) =  &
          +k(521)*n(idx_H2j)  &
          +k(1199)*n(idx_Hj)
    elseif(j==36) then
      pdj(2) =  &
          -k(927)*n(idx_CH)
      pdj(3) =  &
          -k(1070)*n(idx_O)  &
          +k(981)*n(idx_H)  &
          -k(1069)*n(idx_O)  &
          -k(1071)*n(idx_O)
      pdj(5) =  &
          +k(945)*n(idx_CN)
      pdj(6) =  &
          +k(982)*n(idx_H)  &
          +k(591)*n(idx_H3j)  &
          +k(504)*n(idx_Hj)
      pdj(8) =  &
          -k(983)*n(idx_H)  &
          +k(1069)*n(idx_O)  &
          +k(687)*n(idx_HEj)  &
          +k(1154)  &
          -k(981)*n(idx_H)  &
          +k(265)  &
          -k(982)*n(idx_H)
      pdj(9) =  &
          +k(1098)*n(idx_OH)
      pdj(10) =  &
          +k(1070)*n(idx_O)  &
          -k(1098)*n(idx_OH)  &
          +k(983)*n(idx_H)
      pdj(11) =  &
          +k(1071)*n(idx_O)
      pdj(12) =  &
          +k(927)*n(idx_CH)  &
          -k(884)*n(idx_CH2)
      pdj(13) =  &
          +k(1000)*n(idx_HCO)
      pdj(14) =  &
          -k(1000)*n(idx_HCO)
      pdj(17) =  &
          +k(982)*n(idx_H)  &
          +k(1000)*n(idx_HCO)  &
          +k(1098)*n(idx_OH)  &
          +k(688)*n(idx_HEj)  &
          +k(907)*n(idx_CH3)  &
          +k(1154)  &
          +k(1070)*n(idx_O)  &
          +k(884)*n(idx_CH2)  &
          +k(927)*n(idx_CH)  &
          +k(265)  &
          +k(945)*n(idx_CN)  &
          +k(1018)*n(idx_N)
      pdj(24) =  &
          -k(945)*n(idx_CN)
      pdj(25) =  &
          -k(952)*n(idx_CO)
      pdj(27) =  &
          +k(981)*n(idx_H)
      pdj(28) =  &
          +k(884)*n(idx_CH2)  &
          -k(907)*n(idx_CH3)
      pdj(29) =  &
          +k(907)*n(idx_CH3)
      pdj(30) =  &
          -k(1018)*n(idx_N)
      pdj(31) =  &
          +k(1071)*n(idx_O)  &
          +k(952)*n(idx_CO)  &
          +k(983)*n(idx_H)  &
          +k(1018)*n(idx_N)
      pdj(35) =  &
          +k(688)*n(idx_HEj)  &
          +k(687)*n(idx_HEj)
      pdj(36) =  &
          -k(983)*n(idx_H)  &
          -k(952)*n(idx_CO)  &
          -k(1070)*n(idx_O)  &
          -k(688)*n(idx_HEj)  &
          -k(591)*n(idx_H3j)  &
          -k(945)*n(idx_CN)  &
          -k(1295)  &
          -k(1098)*n(idx_OH)  &
          -k(1154)  &
          -k(982)*n(idx_H)  &
          -k(1018)*n(idx_N)  &
          -k(1069)*n(idx_O)  &
          -k(504)*n(idx_Hj)  &
          -k(265)  &
          -k(927)*n(idx_CH)  &
          -k(884)*n(idx_CH2)  &
          -k(1071)*n(idx_O)  &
          -k(687)*n(idx_HEj)  &
          -k(1000)*n(idx_HCO)  &
          -k(981)*n(idx_H)  &
          -k(907)*n(idx_CH3)
      pdj(38) =  &
          +k(952)*n(idx_CO)
      pdj(42) =  &
          +k(1069)*n(idx_O)
      pdj(63) =  &
          +k(1295)
      pdj(71) =  &
          +k(688)*n(idx_HEj)  &
          -k(504)*n(idx_Hj)
      pdj(79) =  &
          +k(687)*n(idx_HEj)  &
          +k(504)*n(idx_Hj)
      pdj(103) =  &
          -k(687)*n(idx_HEj)  &
          -k(688)*n(idx_HEj)
      pdj(105) =  &
          +k(591)*n(idx_H3j)
      pdj(106) =  &
          -k(591)*n(idx_H3j)
    elseif(j==37) then
      pdj(1) =  &
          +k(1121)
      pdj(2) =  &
          +k(366)*n(idx_Cj)
      pdj(6) =  &
          +k(1120)  &
          +k(248)  &
          +2.d0*k(495)*n(idx_Hj)  &
          +k(494)*n(idx_Hj)  &
          +k(580)*n(idx_H3j)
      pdj(8) =  &
          +k(1121)  &
          +k(821)*n(idx_O2j)  &
          +k(715)*n(idx_Nj)  &
          +k(716)*n(idx_Nj)  &
          +k(713)*n(idx_Nj)
      pdj(9) =  &
          +k(493)*n(idx_Hj)  &
          +k(580)*n(idx_H3j)  &
          +k(809)*n(idx_Oj)
      pdj(10) =  &
          +k(249)  &
          +k(658)*n(idx_HEj)  &
          +k(1122)  &
          +k(810)*n(idx_Oj)
      pdj(11) =  &
          +k(821)*n(idx_O2j)
      pdj(12) =  &
          +k(398)*n(idx_CHj)
      pdj(13) =  &
          +k(248)  &
          +k(1120)  &
          +k(397)*n(idx_CHj)
      pdj(14) =  &
          +k(367)*n(idx_Cj)
      pdj(17) =  &
          +k(716)*n(idx_Nj)
      pdj(28) =  &
          +k(249)  &
          +k(657)*n(idx_HEj)  &
          +k(715)*n(idx_Nj)  &
          +k(1122)  &
          +k(861)*n(idx_SIj)
      pdj(29) =  &
          +k(440)*n(idx_CH3j)
      pdj(31) =  &
          +k(713)*n(idx_Nj)  &
          +k(714)*n(idx_Nj)
      pdj(35) =  &
          +k(657)*n(idx_HEj)  &
          +k(658)*n(idx_HEj)
      pdj(37) =  &
          -k(716)*n(idx_Nj)  &
          -k(495)*n(idx_Hj)  &
          -k(809)*n(idx_Oj)  &
          -k(810)*n(idx_Oj)  &
          -k(249)  &
          -k(1224)  &
          -k(494)*n(idx_Hj)  &
          -k(1120)  &
          -k(1122)  &
          -k(398)*n(idx_CHj)  &
          -k(714)*n(idx_Nj)  &
          -k(248)  &
          -k(657)*n(idx_HEj)  &
          -k(397)*n(idx_CHj)  &
          -k(658)*n(idx_HEj)  &
          -k(366)*n(idx_Cj)  &
          -k(821)*n(idx_O2j)  &
          -k(440)*n(idx_CH3j)  &
          -k(713)*n(idx_Nj)  &
          -k(861)*n(idx_SIj)  &
          -k(1121)  &
          -k(367)*n(idx_Cj)  &
          -k(493)*n(idx_Hj)  &
          -k(715)*n(idx_Nj)  &
          -k(580)*n(idx_H3j)
      pdj(45) =  &
          +k(1224)
      pdj(70) =  &
          +k(495)*n(idx_Hj)
      pdj(71) =  &
          -k(495)*n(idx_Hj)  &
          -k(494)*n(idx_Hj)  &
          -k(493)*n(idx_Hj)
      pdj(73) =  &
          -k(366)*n(idx_Cj)  &
          -k(367)*n(idx_Cj)
      pdj(75) =  &
          -k(397)*n(idx_CHj)  &
          -k(398)*n(idx_CHj)
      pdj(76) =  &
          +k(809)*n(idx_Oj)  &
          +k(713)*n(idx_Nj)
      pdj(79) =  &
          +k(715)*n(idx_Nj)
      pdj(80) =  &
          -k(861)*n(idx_SIj)
      pdj(89) =  &
          -k(821)*n(idx_O2j)
      pdj(92) =  &
          -k(810)*n(idx_Oj)  &
          -k(809)*n(idx_Oj)
      pdj(93) =  &
          +k(657)*n(idx_HEj)
      pdj(94) =  &
          +k(493)*n(idx_Hj)  &
          +k(367)*n(idx_Cj)  &
          -k(440)*n(idx_CH3j)  &
          +k(658)*n(idx_HEj)  &
          +k(580)*n(idx_H3j)  &
          +k(397)*n(idx_CHj)  &
          +k(716)*n(idx_Nj)
      pdj(96) =  &
          -k(713)*n(idx_Nj)  &
          -k(715)*n(idx_Nj)  &
          -k(714)*n(idx_Nj)  &
          -k(716)*n(idx_Nj)
      pdj(103) =  &
          -k(657)*n(idx_HEj)  &
          -k(658)*n(idx_HEj)
      pdj(106) =  &
          -k(580)*n(idx_H3j)
      pdj(107) =  &
          +k(714)*n(idx_Nj)  &
          +k(1121)  &
          +k(440)*n(idx_CH3j)  &
          +k(494)*n(idx_Hj)  &
          +k(398)*n(idx_CHj)  &
          +k(366)*n(idx_Cj)  &
          +k(821)*n(idx_O2j)  &
          +k(810)*n(idx_Oj)
      pdj(115) =  &
          +k(861)*n(idx_SIj)
    elseif(j==38) then
      pdj(2) =  &
          -k(924)*n(idx_CH)
      pdj(3) =  &
          +k(666)*n(idx_HEj)  &
          -k(1060)*n(idx_O)  &
          +k(1133)  &
          +k(838)*n(idx_OHj)  &
          +k(497)*n(idx_Hj)  &
          +k(253)
      pdj(6) =  &
          +k(583)*n(idx_H3j)
      pdj(7) =  &
          +k(668)*n(idx_HEj)
      pdj(8) =  &
          -k(972)*n(idx_H)  &
          +k(515)*n(idx_H2j)
      pdj(10) =  &
          +k(972)*n(idx_H)
      pdj(11) =  &
          +k(669)*n(idx_HEj)  &
          +k(1060)*n(idx_O)  &
          +k(822)*n(idx_O2Hj)
      pdj(14) =  &
          +k(924)*n(idx_CH)  &
          +k(751)*n(idx_NHj)
      pdj(17) =  &
          +k(653)*n(idx_HNOj)  &
          +k(1013)*n(idx_N)  &
          +k(720)*n(idx_Nj)
      pdj(18) =  &
          -k(1104)*n(idx_SI)
      pdj(24) =  &
          +k(621)*n(idx_HCNj)
      pdj(25) =  &
          +k(1060)*n(idx_O)  &
          +k(667)*n(idx_HEj)  &
          +k(750)*n(idx_NHj)  &
          +k(924)*n(idx_CH)  &
          +k(813)*n(idx_Oj)  &
          +k(399)*n(idx_CHj)  &
          +k(1104)*n(idx_SI)  &
          +k(972)*n(idx_H)  &
          +k(1133)  &
          +k(253)  &
          +k(1013)*n(idx_N)  &
          +k(368)*n(idx_Cj)  &
          +k(417)*n(idx_CH2j)
      pdj(26) =  &
          +k(735)*n(idx_N2Hj)
      pdj(28) =  &
          +k(448)*n(idx_CH4j)
      pdj(30) =  &
          +k(749)*n(idx_NHj)  &
          -k(1013)*n(idx_N)
      pdj(34) =  &
          +k(1104)*n(idx_SI)
      pdj(35) =  &
          +k(669)*n(idx_HEj)  &
          +k(668)*n(idx_HEj)  &
          +k(667)*n(idx_HEj)  &
          +k(666)*n(idx_HEj)
      pdj(38) =  &
          -k(621)*n(idx_HCNj)  &
          -k(253)  &
          -k(720)*n(idx_Nj)  &
          -k(669)*n(idx_HEj)  &
          -k(666)*n(idx_HEj)  &
          -k(399)*n(idx_CHj)  &
          -k(924)*n(idx_CH)  &
          -k(583)*n(idx_H3j)  &
          -k(1104)*n(idx_SI)  &
          -k(813)*n(idx_Oj)  &
          -k(497)*n(idx_Hj)  &
          -k(751)*n(idx_NHj)  &
          -k(822)*n(idx_O2Hj)  &
          -k(417)*n(idx_CH2j)  &
          -k(515)*n(idx_H2j)  &
          -k(368)*n(idx_Cj)  &
          -k(448)*n(idx_CH4j)  &
          -k(735)*n(idx_N2Hj)  &
          -k(838)*n(idx_OHj)  &
          -k(749)*n(idx_NHj)  &
          -k(667)*n(idx_HEj)  &
          -k(972)*n(idx_H)  &
          -k(1259)  &
          -k(653)*n(idx_HNOj)  &
          -k(1060)*n(idx_O)  &
          -k(750)*n(idx_NHj)  &
          -k(1133)  &
          -k(1013)*n(idx_N)  &
          -k(668)*n(idx_HEj)
      pdj(57) =  &
          +k(1259)
      pdj(70) =  &
          +k(497)*n(idx_Hj)  &
          +k(399)*n(idx_CHj)
      pdj(71) =  &
          -k(497)*n(idx_Hj)
      pdj(73) =  &
          +k(669)*n(idx_HEj)  &
          -k(368)*n(idx_Cj)
      pdj(74) =  &
          -k(417)*n(idx_CH2j)
      pdj(75) =  &
          -k(399)*n(idx_CHj)
      pdj(76) =  &
          +k(417)*n(idx_CH2j)
      pdj(79) =  &
          +k(751)*n(idx_NHj)
      pdj(87) =  &
          +k(720)*n(idx_Nj)  &
          +k(368)*n(idx_Cj)  &
          +k(666)*n(idx_HEj)
      pdj(89) =  &
          +k(668)*n(idx_HEj)  &
          +k(813)*n(idx_Oj)
      pdj(92) =  &
          -k(813)*n(idx_Oj)  &
          +k(667)*n(idx_HEj)
      pdj(93) =  &
          -k(838)*n(idx_OHj)
      pdj(95) =  &
          -k(448)*n(idx_CH4j)
      pdj(96) =  &
          -k(720)*n(idx_Nj)
      pdj(97) =  &
          -k(621)*n(idx_HCNj)
      pdj(98) =  &
          -k(751)*n(idx_NHj)  &
          -k(750)*n(idx_NHj)  &
          -k(749)*n(idx_NHj)
      pdj(102) =  &
          -k(515)*n(idx_H2j)
      pdj(103) =  &
          -k(667)*n(idx_HEj)  &
          -k(666)*n(idx_HEj)  &
          -k(669)*n(idx_HEj)  &
          -k(668)*n(idx_HEj)
      pdj(104) =  &
          -k(653)*n(idx_HNOj)  &
          +k(750)*n(idx_NHj)
      pdj(106) =  &
          -k(583)*n(idx_H3j)
      pdj(110) =  &
          +k(735)*n(idx_N2Hj)  &
          +k(822)*n(idx_O2Hj)  &
          +k(653)*n(idx_HNOj)  &
          +k(749)*n(idx_NHj)  &
          +k(515)*n(idx_H2j)  &
          +k(838)*n(idx_OHj)  &
          +k(448)*n(idx_CH4j)  &
          +k(583)*n(idx_H3j)  &
          +k(621)*n(idx_HCNj)
      pdj(112) =  &
          -k(735)*n(idx_N2Hj)
      pdj(113) =  &
          -k(822)*n(idx_O2Hj)
    elseif(j==39) then
      pdj(3) =  &
          -k(1061)*n(idx_O)
      pdj(5) =  &
          +k(255)  &
          +k(974)*n(idx_H)  &
          +k(1014)*n(idx_N)  &
          +k(1136)
      pdj(6) =  &
          +k(974)*n(idx_H)  &
          +k(1061)*n(idx_O)
      pdj(8) =  &
          +k(255)  &
          -k(974)*n(idx_H)  &
          +k(1136)
      pdj(30) =  &
          -k(1014)*n(idx_N)
      pdj(31) =  &
          +k(1014)*n(idx_N)
      pdj(39) =  &
          -k(974)*n(idx_H)  &
          -k(1299)  &
          -k(1014)*n(idx_N)  &
          -k(255)  &
          -k(1136)  &
          -k(1061)*n(idx_O)
      pdj(44) =  &
          +k(1061)*n(idx_O)
      pdj(65) =  &
          +k(1299)
    elseif(j==40) then
      pdj(6) =  &
          +k(1144)  &
          +k(500)*n(idx_Hj)  &
          +k(258)
      pdj(8) =  &
          +k(676)*n(idx_HEj)  &
          +2.d0*k(1145)
      pdj(34) =  &
          +k(1144)  &
          +k(1145)  &
          +k(258)
      pdj(35) =  &
          +k(676)*n(idx_HEj)
      pdj(40) =  &
          -k(1248)  &
          -k(1144)  &
          -k(1145)  &
          -k(500)*n(idx_Hj)  &
          -k(676)*n(idx_HEj)  &
          -k(258)
      pdj(49) =  &
          +k(1248)
      pdj(71) =  &
          -k(500)*n(idx_Hj)
      pdj(103) =  &
          -k(676)*n(idx_HEj)
      pdj(115) =  &
          +k(676)*n(idx_HEj)  &
          +k(500)*n(idx_Hj)
    elseif(j==41) then
      pdj(4) =  &
          +k(1005)*n(idx_C)
      pdj(7) =  &
          -k(1005)*n(idx_C)
      pdj(25) =  &
          +k(503)*n(idx_Hj)  &
          +k(264)  &
          +k(1153)  &
          +k(1005)*n(idx_C)
      pdj(31) =  &
          +k(264)  &
          +k(1153)
      pdj(41) =  &
          -k(1153)  &
          -k(1005)*n(idx_C)  &
          -k(503)*n(idx_Hj)  &
          -k(264)  &
          -k(1225)
      pdj(46) =  &
          +k(1225)
      pdj(71) =  &
          -k(503)*n(idx_Hj)
      pdj(91) =  &
          +k(503)*n(idx_Hj)
    elseif(j==42) then
      pdj(3) =  &
          +k(277)  &
          +2.d0*k(1020)*n(idx_N)  &
          +k(1165)  &
          -k(1076)*n(idx_O)
      pdj(6) =  &
          +k(596)*n(idx_H3j)
      pdj(8) =  &
          -k(987)*n(idx_H)
      pdj(10) =  &
          +k(987)*n(idx_H)  &
          +k(596)*n(idx_H3j)  &
          +k(505)*n(idx_Hj)
      pdj(11) =  &
          +k(819)*n(idx_Oj)  &
          +k(1022)*n(idx_N)  &
          +k(1076)*n(idx_O)
      pdj(12) =  &
          -k(886)*n(idx_CH2)
      pdj(13) =  &
          +k(886)*n(idx_CH2)  &
          +k(910)*n(idx_CH3)
      pdj(17) =  &
          +k(946)*n(idx_CN)  &
          +k(277)  &
          +k(987)*n(idx_H)  &
          +k(1076)*n(idx_O)  &
          +2.d0*k(1021)*n(idx_N)  &
          +k(953)*n(idx_CO)  &
          +k(886)*n(idx_CH2)  &
          +k(1042)*n(idx_NH)  &
          +k(1165)
      pdj(24) =  &
          -k(946)*n(idx_CN)
      pdj(25) =  &
          -k(953)*n(idx_CO)
      pdj(26) =  &
          +k(1022)*n(idx_N)  &
          +k(1020)*n(idx_N)
      pdj(28) =  &
          -k(910)*n(idx_CH3)
      pdj(30) =  &
          -k(1020)*n(idx_N)  &
          -k(1021)*n(idx_N)  &
          -k(1022)*n(idx_N)
      pdj(31) =  &
          -k(1042)*n(idx_NH)
      pdj(36) =  &
          +k(1042)*n(idx_NH)  &
          +k(910)*n(idx_CH3)
      pdj(38) =  &
          +k(953)*n(idx_CO)
      pdj(42) =  &
          -k(819)*n(idx_Oj)  &
          -k(505)*n(idx_Hj)  &
          -k(953)*n(idx_CO)  &
          -k(1022)*n(idx_N)  &
          -k(1076)*n(idx_O)  &
          -k(596)*n(idx_H3j)  &
          -k(1294)  &
          -k(946)*n(idx_CN)  &
          -k(1165)  &
          -k(1042)*n(idx_NH)  &
          -k(987)*n(idx_H)  &
          -k(277)  &
          -k(886)*n(idx_CH2)  &
          -k(1020)*n(idx_N)  &
          -k(1021)*n(idx_N)  &
          -k(910)*n(idx_CH3)
      pdj(44) =  &
          +k(946)*n(idx_CN)
      pdj(62) =  &
          +k(1294)
      pdj(71) =  &
          -k(505)*n(idx_Hj)
      pdj(79) =  &
          +k(819)*n(idx_Oj)  &
          +k(596)*n(idx_H3j)  &
          +k(505)*n(idx_Hj)
      pdj(92) =  &
          -k(819)*n(idx_Oj)
      pdj(106) =  &
          -k(596)*n(idx_H3j)
    elseif(j==43) then
      pdj(2) =  &
          -k(939)*n(idx_CH)  &
          -k(938)*n(idx_CH)
      pdj(3) =  &
          +k(1172)  &
          +k(991)*n(idx_H)  &
          -k(1078)*n(idx_O)
      pdj(6) =  &
          +k(992)*n(idx_H)
      pdj(8) =  &
          +k(282)  &
          -k(991)*n(idx_H)  &
          -k(993)*n(idx_H)  &
          +k(1171)  &
          -k(992)*n(idx_H)
      pdj(9) =  &
          +k(1101)*n(idx_OH)  &
          +k(991)*n(idx_H)
      pdj(10) =  &
          +k(955)*n(idx_CO)  &
          -k(1101)*n(idx_OH)  &
          +k(1078)*n(idx_O)  &
          +2.d0*k(993)*n(idx_H)  &
          +k(938)*n(idx_CH)  &
          +k(1172)
      pdj(11) =  &
          +k(992)*n(idx_H)  &
          +k(915)*n(idx_CH3)  &
          +k(1078)*n(idx_O)  &
          +k(939)*n(idx_CH)  &
          +k(282)  &
          +k(1171)  &
          +k(1025)*n(idx_N)  &
          +k(1004)*n(idx_HCO)  &
          +k(1101)*n(idx_OH)
      pdj(12) =  &
          +k(939)*n(idx_CH)
      pdj(13) =  &
          +k(1004)*n(idx_HCO)
      pdj(14) =  &
          +k(938)*n(idx_CH)  &
          -k(1004)*n(idx_HCO)
      pdj(25) =  &
          -k(955)*n(idx_CO)
      pdj(28) =  &
          -k(915)*n(idx_CH3)
      pdj(29) =  &
          +k(915)*n(idx_CH3)
      pdj(30) =  &
          -k(1025)*n(idx_N)
      pdj(31) =  &
          +k(1025)*n(idx_N)
      pdj(38) =  &
          +k(955)*n(idx_CO)
      pdj(43) =  &
          -k(1101)*n(idx_OH)  &
          -k(939)*n(idx_CH)  &
          -k(955)*n(idx_CO)  &
          -k(991)*n(idx_H)  &
          -k(1004)*n(idx_HCO)  &
          -k(1025)*n(idx_N)  &
          -k(993)*n(idx_H)  &
          -k(1078)*n(idx_O)  &
          -k(1172)  &
          -k(915)*n(idx_CH3)  &
          -k(938)*n(idx_CH)  &
          -k(1171)  &
          -k(1304)  &
          -k(282)  &
          -k(992)*n(idx_H)
      pdj(64) =  &
          +k(1304)
    elseif(j==44) then
      pdj(3) =  &
          -k(1079)*n(idx_O)  &
          -k(1080)*n(idx_O)  &
          +k(284)  &
          +k(698)*n(idx_HEj)  &
          +k(1173)  &
          +k(994)*n(idx_H)
      pdj(5) =  &
          +k(994)*n(idx_H)
      pdj(7) =  &
          -k(875)*n(idx_C)
      pdj(8) =  &
          -k(994)*n(idx_H)  &
          -k(996)*n(idx_H)  &
          -k(995)*n(idx_H)
      pdj(10) =  &
          +k(996)*n(idx_H)
      pdj(11) =  &
          -k(1056)*n(idx_O2)  &
          +k(1080)*n(idx_O)  &
          -k(1055)*n(idx_O2)
      pdj(17) =  &
          +k(1055)*n(idx_O2)  &
          -k(1054)*n(idx_NO)  &
          +k(1079)*n(idx_O)
      pdj(24) =  &
          +k(1173)  &
          +k(284)  &
          +k(379)*n(idx_Cj)  &
          +k(996)*n(idx_H)  &
          +k(875)*n(idx_C)  &
          +k(699)*n(idx_HEj)  &
          +k(1080)*n(idx_O)
      pdj(25) =  &
          +k(1079)*n(idx_O)  &
          +k(875)*n(idx_C)  &
          +k(1056)*n(idx_O2)  &
          +k(995)*n(idx_H)
      pdj(26) =  &
          +k(1054)*n(idx_NO)
      pdj(31) =  &
          +k(995)*n(idx_H)
      pdj(35) =  &
          +k(698)*n(idx_HEj)  &
          +k(699)*n(idx_HEj)
      pdj(38) =  &
          +k(1055)*n(idx_O2)  &
          +k(1054)*n(idx_NO)
      pdj(42) =  &
          +k(1056)*n(idx_O2)
      pdj(44) =  &
          -k(284)  &
          -k(994)*n(idx_H)  &
          -k(1079)*n(idx_O)  &
          -k(875)*n(idx_C)  &
          -k(1226)  &
          -k(1080)*n(idx_O)  &
          -k(699)*n(idx_HEj)  &
          -k(1173)  &
          -k(1056)*n(idx_O2)  &
          -k(698)*n(idx_HEj)  &
          -k(1055)*n(idx_O2)  &
          -k(379)*n(idx_Cj)  &
          -k(996)*n(idx_H)  &
          -k(995)*n(idx_H)  &
          -k(1054)*n(idx_NO)
      pdj(46) =  &
          +k(1226)
      pdj(73) =  &
          -k(379)*n(idx_Cj)
      pdj(86) =  &
          +k(698)*n(idx_HEj)
      pdj(87) =  &
          +k(379)*n(idx_Cj)
      pdj(92) =  &
          +k(699)*n(idx_HEj)
      pdj(103) =  &
          -k(698)*n(idx_HEj)  &
          -k(699)*n(idx_HEj)
    elseif(j==45) then
      pdj(37) =  &
          +k(1321)
      pdj(45) =  &
          -k(1321)
    elseif(j==46) then
      pdj(41) =  &
          +k(1323)
      pdj(46) =  &
          -k(1323)
    elseif(j==47) then
      pdj(13) =  &
          +k(1318)
      pdj(47) =  &
          -k(1318)
    elseif(j==48) then
      pdj(32) =  &
          +k(1326)
      pdj(48) =  &
          -k(1326)
    elseif(j==49) then
      pdj(40) =  &
          +k(1329)
      pdj(49) =  &
          -k(1329)
    elseif(j==50) then
      pdj(21) =  &
          +k(1327)
      pdj(50) =  &
          -k(1327)
    elseif(j==51) then
      pdj(19) =  &
          +k(1330)
      pdj(51) =  &
          -k(1330)
    elseif(j==52) then
      pdj(20) =  &
          +k(1331)
      pdj(52) =  &
          -k(1331)
    elseif(j==53) then
      pdj(29) =  &
          +k(1308)
      pdj(53) =  &
          -k(1308)
    elseif(j==54) then
      pdj(25) =  &
          +k(1314)
      pdj(54) =  &
          -k(1314)
    elseif(j==55) then
      pdj(9) =  &
          +k(1310)
      pdj(55) =  &
          -k(1310)
    elseif(j==56) then
      pdj(17) =  &
          +k(1317)
      pdj(56) =  &
          -k(1317)
    elseif(j==57) then
      pdj(38) =  &
          +k(1324)
      pdj(57) =  &
          -k(1324)
    elseif(j==58) then
      pdj(26) =  &
          +k(1315)
      pdj(58) =  &
          -k(1315)
    elseif(j==59) then
      pdj(5) =  &
          +k(1312)
      pdj(59) =  &
          -k(1312)
    elseif(j==60) then
      pdj(16) =  &
          +k(1309)
      pdj(60) =  &
          -k(1309)
    elseif(j==61) then
      pdj(11) =  &
          +k(1320)
      pdj(61) =  &
          -k(1320)
    elseif(j==62) then
      pdj(42) =  &
          +k(1325)
      pdj(62) =  &
          -k(1325)
    elseif(j==63) then
      pdj(36) =  &
          +k(1319)
      pdj(63) =  &
          -k(1319)
    elseif(j==64) then
      pdj(43) =  &
          +k(1322)
      pdj(64) =  &
          -k(1322)
    elseif(j==65) then
      pdj(39) =  &
          +k(1316)
      pdj(65) =  &
          -k(1316)
    elseif(j==66) then
      pdj(15) =  &
          +k(1311)
      pdj(66) =  &
          -k(1311)
    elseif(j==67) then
      pdj(4) =  &
          +k(1313)
      pdj(67) =  &
          -k(1313)
    elseif(j==68) then
    elseif(j==69) then
      pdj(34) =  &
          +k(1328)
      pdj(69) =  &
          -k(1328)
    elseif(j==70) then
      pdj(1) =  &
          -k(331)*n(idx_E)
      pdj(2) =  &
          -k(467)*n(idx_CH)
      pdj(4) =  &
          -k(649)*n(idx_HNC)
      pdj(5) =  &
          -k(630)*n(idx_HCN)
      pdj(7) =  &
          -k(387)*n(idx_C)
      pdj(8) =  &
          +k(331)*n(idx_E)  &
          +k(856)*n(idx_OH)  &
          +k(1149)
      pdj(9) =  &
          -k(567)*n(idx_H2O)
      pdj(10) =  &
          -k(855)*n(idx_OH)  &
          -k(856)*n(idx_OH)
      pdj(12) =  &
          -k(430)*n(idx_CH2)
      pdj(13) =  &
          -k(636)*n(idx_H2CO)
      pdj(14) =  &
          +k(150)*n(idx_MG)  &
          -k(637)*n(idx_HCO)
      pdj(15) =  &
          -k(150)*n(idx_MG)
      pdj(18) =  &
          -k(862)*n(idx_SI)
      pdj(22) =  &
          -k(638)*n(idx_SIH2)
      pdj(25) =  &
          +k(430)*n(idx_CH2)  &
          +k(630)*n(idx_HCN)  &
          +k(862)*n(idx_SI)  &
          +k(649)*n(idx_HNC)  &
          +k(567)*n(idx_H2O)  &
          +k(387)*n(idx_C)  &
          +k(467)*n(idx_CH)  &
          +k(638)*n(idx_SIH2)  &
          +k(855)*n(idx_OH)  &
          +k(331)*n(idx_E)  &
          +k(800)*n(idx_NH)  &
          +k(636)*n(idx_H2CO)  &
          +k(637)*n(idx_HCO)  &
          +k(788)*n(idx_NH2)  &
          +k(640)*n(idx_SIH)  &
          +k(641)*n(idx_SIO)  &
          +k(639)*n(idx_SIH4)
      pdj(27) =  &
          -k(788)*n(idx_NH2)
      pdj(31) =  &
          -k(800)*n(idx_NH)
      pdj(32) =  &
          -k(639)*n(idx_SIH4)
      pdj(33) =  &
          -k(640)*n(idx_SIH)
      pdj(34) =  &
          -k(641)*n(idx_SIO)
      pdj(47) =  &
          +k(1282)
      pdj(70) =  &
          -k(638)*n(idx_SIH2)  &
          -k(640)*n(idx_SIH)  &
          -k(862)*n(idx_SI)  &
          -k(430)*n(idx_CH2)  &
          -k(150)*n(idx_MG)  &
          -k(636)*n(idx_H2CO)  &
          -k(855)*n(idx_OH)  &
          -k(630)*n(idx_HCN)  &
          -k(567)*n(idx_H2O)  &
          -k(649)*n(idx_HNC)  &
          -k(641)*n(idx_SIO)  &
          -k(331)*n(idx_E)  &
          -k(800)*n(idx_NH)  &
          -k(856)*n(idx_OH)  &
          -k(467)*n(idx_CH)  &
          -k(637)*n(idx_HCO)  &
          -k(387)*n(idx_C)  &
          -k(1149)  &
          -k(788)*n(idx_NH2)  &
          -k(639)*n(idx_SIH4)  &
          -k(1282)
      pdj(74) =  &
          +k(467)*n(idx_CH)
      pdj(75) =  &
          +k(387)*n(idx_C)
      pdj(76) =  &
          +k(637)*n(idx_HCO)
      pdj(77) =  &
          +k(150)*n(idx_MG)
      pdj(78) =  &
          +k(788)*n(idx_NH2)
      pdj(84) =  &
          +k(640)*n(idx_SIH)
      pdj(85) =  &
          +k(638)*n(idx_SIH2)
      pdj(87) =  &
          +k(1149)
      pdj(90) =  &
          +k(855)*n(idx_OH)
      pdj(91) =  &
          +k(800)*n(idx_NH)
      pdj(94) =  &
          +k(430)*n(idx_CH2)
      pdj(100) =  &
          +k(862)*n(idx_SI)
      pdj(107) =  &
          +k(636)*n(idx_H2CO)
      pdj(108) =  &
          +k(567)*n(idx_H2O)
      pdj(109) =  &
          +k(630)*n(idx_HCN)  &
          +k(649)*n(idx_HNC)
      pdj(110) =  &
          +k(856)*n(idx_OH)
      pdj(114) =  &
          +k(639)*n(idx_SIH4)
      pdj(115) =  &
          +k(641)*n(idx_SIO)
    elseif(j==71) then
      pdj(1) =  &
          -k(1217)*n(idx_E)
      pdj(2) =  &
          -k(79)*n(idx_CH)
      pdj(3) =  &
          +k(497)*n(idx_CO2)  &
          -k(90)*n(idx_O)
      pdj(4) =  &
          -k(2)*n(idx_HNC)
      pdj(5) =  &
          -k(82)*n(idx_HCN)  &
          +k(2)*n(idx_HNC)
      pdj(6) =  &
          +k(492)*n(idx_CH2)  &
          +2.d0*k(495)*n(idx_CH3OH)  &
          +k(508)*n(idx_SIH4)  &
          +k(506)*n(idx_SIH2)  &
          +k(501)*n(idx_HCO)  &
          +k(496)*n(idx_CH4)  &
          +k(504)*n(idx_HNO)  &
          +k(509)*n(idx_SIH)  &
          +k(500)*n(idx_H2SIO)  &
          +k(507)*n(idx_SIH3)  &
          +k(498)*n(idx_H2CO)  &
          +k(499)*n(idx_H2CO)  &
          +k(494)*n(idx_CH3OH)
      pdj(8) =  &
          +k(83)*n(idx_HCO)  &
          +k(79)*n(idx_CH)  &
          -k(1198)*n(idx_H)  &
          +k(78)*n(idx_CH4)  &
          +k(98)*n(idx_SIH4)  &
          +k(88)*n(idx_NO)  &
          +k(95)*n(idx_SIC)  &
          +k(100)*n(idx_SIO)  &
          +k(498)*n(idx_H2CO)  &
          +k(80)*n(idx_H2CO)  &
          +k(96)*n(idx_SIH2)  &
          +k(84)*n(idx_MG)  &
          +k(89)*n(idx_O2)  &
          +k(97)*n(idx_SIH3)  &
          +k(86)*n(idx_NH3)  &
          +k(87)*n(idx_NH)  &
          +k(1217)*n(idx_E)  &
          +k(85)*n(idx_NH2)  &
          +k(93)*n(idx_SIC2)  &
          +k(82)*n(idx_HCN)  &
          +k(99)*n(idx_SIH)  &
          +k(94)*n(idx_SIC3)  &
          +k(92)*n(idx_SI)  &
          +k(90)*n(idx_O)  &
          +k(91)*n(idx_OH)  &
          +k(76)*n(idx_CH2)  &
          +k(77)*n(idx_CH3)  &
          +k(81)*n(idx_H2O)
      pdj(9) =  &
          -k(81)*n(idx_H2O)  &
          +k(493)*n(idx_CH3OH)
      pdj(10) =  &
          +k(505)*n(idx_NO2)  &
          -k(91)*n(idx_OH)
      pdj(11) =  &
          -k(89)*n(idx_O2)
      pdj(12) =  &
          -k(492)*n(idx_CH2)  &
          -k(76)*n(idx_CH2)
      pdj(13) =  &
          -k(498)*n(idx_H2CO)  &
          -k(499)*n(idx_H2CO)  &
          -k(80)*n(idx_H2CO)
      pdj(14) =  &
          -k(502)*n(idx_HCO)  &
          -k(501)*n(idx_HCO)  &
          -k(83)*n(idx_HCO)
      pdj(15) =  &
          -k(84)*n(idx_MG)
      pdj(16) =  &
          -k(86)*n(idx_NH3)
      pdj(17) =  &
          -k(88)*n(idx_NO)
      pdj(18) =  &
          -k(92)*n(idx_SI)
      pdj(19) =  &
          -k(93)*n(idx_SIC2)
      pdj(20) =  &
          -k(94)*n(idx_SIC3)
      pdj(21) =  &
          -k(95)*n(idx_SIC)
      pdj(22) =  &
          -k(506)*n(idx_SIH2)  &
          -k(96)*n(idx_SIH2)
      pdj(23) =  &
          -k(97)*n(idx_SIH3)  &
          -k(507)*n(idx_SIH3)
      pdj(25) =  &
          +k(502)*n(idx_HCO)  &
          +k(503)*n(idx_HNCO)
      pdj(27) =  &
          -k(85)*n(idx_NH2)
      pdj(28) =  &
          -k(77)*n(idx_CH3)
      pdj(29) =  &
          -k(78)*n(idx_CH4)  &
          -k(496)*n(idx_CH4)
      pdj(31) =  &
          -k(87)*n(idx_NH)
      pdj(32) =  &
          -k(98)*n(idx_SIH4)  &
          -k(508)*n(idx_SIH4)
      pdj(33) =  &
          -k(509)*n(idx_SIH)  &
          -k(99)*n(idx_SIH)
      pdj(34) =  &
          -k(100)*n(idx_SIO)
      pdj(35) =  &
          -k(1199)*n(idx_HE)
      pdj(36) =  &
          -k(504)*n(idx_HNO)
      pdj(37) =  &
          -k(495)*n(idx_CH3OH)  &
          -k(493)*n(idx_CH3OH)  &
          -k(494)*n(idx_CH3OH)
      pdj(38) =  &
          -k(497)*n(idx_CO2)
      pdj(40) =  &
          -k(500)*n(idx_H2SIO)
      pdj(41) =  &
          -k(503)*n(idx_HNCO)
      pdj(42) =  &
          -k(505)*n(idx_NO2)
      pdj(70) =  &
          +k(83)*n(idx_HCO)  &
          +k(495)*n(idx_CH3OH)  &
          +k(499)*n(idx_H2CO)  &
          +k(497)*n(idx_CO2)
      pdj(71) =  &
          -k(1198)*n(idx_H)  &
          -k(498)*n(idx_H2CO)  &
          -k(84)*n(idx_MG)  &
          -k(86)*n(idx_NH3)  &
          -k(493)*n(idx_CH3OH)  &
          -k(504)*n(idx_HNO)  &
          -k(87)*n(idx_NH)  &
          -k(509)*n(idx_SIH)  &
          -k(99)*n(idx_SIH)  &
          -k(92)*n(idx_SI)  &
          -k(79)*n(idx_CH)  &
          -k(500)*n(idx_H2SIO)  &
          -k(2)*n(idx_HNC)  &
          -k(89)*n(idx_O2)  &
          -k(497)*n(idx_CO2)  &
          -k(503)*n(idx_HNCO)  &
          -k(78)*n(idx_CH4)  &
          -k(85)*n(idx_NH2)  &
          -k(501)*n(idx_HCO)  &
          -k(495)*n(idx_CH3OH)  &
          -k(1217)*n(idx_E)  &
          -k(506)*n(idx_SIH2)  &
          -k(91)*n(idx_OH)  &
          -k(80)*n(idx_H2CO)  &
          -k(82)*n(idx_HCN)  &
          -k(95)*n(idx_SIC)  &
          -k(76)*n(idx_CH2)  &
          -k(96)*n(idx_SIH2)  &
          -k(81)*n(idx_H2O)  &
          -k(492)*n(idx_CH2)  &
          -k(508)*n(idx_SIH4)  &
          -k(505)*n(idx_NO2)  &
          -k(507)*n(idx_SIH3)  &
          -k(83)*n(idx_HCO)  &
          -k(502)*n(idx_HCO)  &
          -k(98)*n(idx_SIH4)  &
          -k(496)*n(idx_CH4)  &
          -k(1199)*n(idx_HE)  &
          -k(77)*n(idx_CH3)  &
          -k(100)*n(idx_SIO)  &
          +k(2)*n(idx_HNC)  &
          -k(88)*n(idx_NO)  &
          -k(94)*n(idx_SIC3)  &
          -k(494)*n(idx_CH3OH)  &
          -k(90)*n(idx_O)  &
          -k(93)*n(idx_SIC2)  &
          -k(97)*n(idx_SIH3)  &
          -k(499)*n(idx_H2CO)
      pdj(74) =  &
          +k(76)*n(idx_CH2)
      pdj(75) =  &
          +k(492)*n(idx_CH2)  &
          +k(79)*n(idx_CH)
      pdj(76) =  &
          +k(80)*n(idx_H2CO)
      pdj(77) =  &
          +k(84)*n(idx_MG)
      pdj(78) =  &
          +k(86)*n(idx_NH3)
      pdj(79) =  &
          +k(88)*n(idx_NO)  &
          +k(504)*n(idx_HNO)  &
          +k(505)*n(idx_NO2)
      pdj(80) =  &
          +k(92)*n(idx_SI)  &
          +k(509)*n(idx_SIH)
      pdj(81) =  &
          +k(93)*n(idx_SIC2)
      pdj(82) =  &
          +k(94)*n(idx_SIC3)
      pdj(83) =  &
          +k(95)*n(idx_SIC)
      pdj(84) =  &
          +k(96)*n(idx_SIH2)  &
          +k(507)*n(idx_SIH3)
      pdj(85) =  &
          +k(508)*n(idx_SIH4)  &
          +k(97)*n(idx_SIH3)
      pdj(87) =  &
          +k(498)*n(idx_H2CO)  &
          +k(501)*n(idx_HCO)
      pdj(89) =  &
          +k(89)*n(idx_O2)
      pdj(90) =  &
          +k(81)*n(idx_H2O)
      pdj(91) =  &
          +k(503)*n(idx_HNCO)  &
          +k(85)*n(idx_NH2)
      pdj(92) =  &
          +k(90)*n(idx_O)
      pdj(93) =  &
          +k(91)*n(idx_OH)
      pdj(94) =  &
          +k(496)*n(idx_CH4)  &
          +k(493)*n(idx_CH3OH)  &
          +k(77)*n(idx_CH3)
      pdj(95) =  &
          +k(78)*n(idx_CH4)
      pdj(97) =  &
          +k(82)*n(idx_HCN)
      pdj(98) =  &
          +k(87)*n(idx_NH)
      pdj(99) =  &
          +k(98)*n(idx_SIH4)
      pdj(100) =  &
          +k(99)*n(idx_SIH)  &
          +k(506)*n(idx_SIH2)
      pdj(101) =  &
          +k(100)*n(idx_SIO)
      pdj(102) =  &
          +k(1198)*n(idx_H)  &
          +k(502)*n(idx_HCO)
      pdj(107) =  &
          +k(494)*n(idx_CH3OH)
      pdj(111) =  &
          +k(1199)*n(idx_HE)
      pdj(115) =  &
          +k(500)*n(idx_H2SIO)
    elseif(j==72) then
      pdj(1) =  &
          -k(336)*n(idx_E)
      pdj(6) =  &
          -k(6)*n(idx_H2)  &
          +k(6)*n(idx_H2)
      pdj(8) =  &
          +k(336)*n(idx_E)
      pdj(25) =  &
          +k(336)*n(idx_E)
      pdj(47) =  &
          +k(1227)
      pdj(70) =  &
          +k(6)*n(idx_H2)
      pdj(72) =  &
          -k(1227)  &
          -k(6)*n(idx_H2)  &
          -k(336)*n(idx_E)
    elseif(j==73) then
      pdj(1) =  &
          -k(1215)*n(idx_E)
      pdj(2) =  &
          +k(366)*n(idx_CH3OH)  &
          -k(16)*n(idx_CH)  &
          +k(370)*n(idx_H2CO)
      pdj(3) =  &
          +k(377)*n(idx_O2)  &
          -k(1194)*n(idx_O)
      pdj(6) =  &
          +k(381)*n(idx_SIH2)  &
          +k(375)*n(idx_NH3)  &
          -k(529)*n(idx_H2)  &
          -k(1200)*n(idx_H2)
      pdj(7) =  &
          +k(24)*n(idx_SIC3)  &
          +k(23)*n(idx_SIC2)  &
          +k(1215)*n(idx_E)  &
          +k(25)*n(idx_SIC)  &
          +k(20)*n(idx_NH3)  &
          +k(21)*n(idx_NO)  &
          +k(26)*n(idx_SIH2)  &
          +k(19)*n(idx_MG)  &
          +k(22)*n(idx_SI)  &
          +k(17)*n(idx_H2CO)  &
          +k(18)*n(idx_HCO)  &
          +k(27)*n(idx_SIH3)  &
          +k(15)*n(idx_CH2)  &
          +k(16)*n(idx_CH)
      pdj(8) =  &
          +k(374)*n(idx_NH2)  &
          +k(380)*n(idx_OH)  &
          +k(529)*n(idx_H2)  &
          +k(376)*n(idx_NH)  &
          -k(1206)*n(idx_H)  &
          +k(371)*n(idx_H2O)  &
          +k(382)*n(idx_SIH)  &
          +k(372)*n(idx_H2O)
      pdj(9) =  &
          -k(371)*n(idx_H2O)  &
          -k(372)*n(idx_H2O)
      pdj(10) =  &
          -k(380)*n(idx_OH)
      pdj(11) =  &
          -k(378)*n(idx_O2)  &
          -k(377)*n(idx_O2)
      pdj(12) =  &
          -k(15)*n(idx_CH2)
      pdj(13) =  &
          -k(370)*n(idx_H2CO)  &
          -k(17)*n(idx_H2CO)  &
          -k(369)*n(idx_H2CO)
      pdj(14) =  &
          -k(373)*n(idx_HCO)  &
          +k(367)*n(idx_CH3OH)  &
          -k(18)*n(idx_HCO)
      pdj(15) =  &
          -k(19)*n(idx_MG)
      pdj(16) =  &
          -k(375)*n(idx_NH3)  &
          -k(20)*n(idx_NH3)
      pdj(17) =  &
          -k(21)*n(idx_NO)
      pdj(18) =  &
          -k(22)*n(idx_SI)
      pdj(19) =  &
          -k(23)*n(idx_SIC2)
      pdj(20) =  &
          -k(24)*n(idx_SIC3)
      pdj(21) =  &
          -k(25)*n(idx_SIC)
      pdj(22) =  &
          -k(381)*n(idx_SIH2)  &
          -k(26)*n(idx_SIH2)
      pdj(23) =  &
          -k(27)*n(idx_SIH3)
      pdj(24) =  &
          +k(379)*n(idx_OCN)
      pdj(25) =  &
          +k(378)*n(idx_O2)  &
          +k(369)*n(idx_H2CO)  &
          +k(383)*n(idx_SIO)  &
          +k(373)*n(idx_HCO)  &
          +k(368)*n(idx_CO2)
      pdj(27) =  &
          -k(374)*n(idx_NH2)
      pdj(30) =  &
          -k(1193)*n(idx_N)
      pdj(31) =  &
          -k(376)*n(idx_NH)
      pdj(33) =  &
          -k(382)*n(idx_SIH)
      pdj(34) =  &
          -k(383)*n(idx_SIO)
      pdj(37) =  &
          -k(367)*n(idx_CH3OH)  &
          -k(366)*n(idx_CH3OH)
      pdj(38) =  &
          -k(368)*n(idx_CO2)
      pdj(44) =  &
          -k(379)*n(idx_OCN)
      pdj(53) =  &
          +k(1265)
      pdj(70) =  &
          +k(371)*n(idx_H2O)  &
          +k(18)*n(idx_HCO)  &
          +k(370)*n(idx_H2CO)
      pdj(72) =  &
          +k(372)*n(idx_H2O)
      pdj(73) =  &
          -k(382)*n(idx_SIH)  &
          -k(374)*n(idx_NH2)  &
          -k(23)*n(idx_SIC2)  &
          -k(27)*n(idx_SIH3)  &
          -k(15)*n(idx_CH2)  &
          -k(377)*n(idx_O2)  &
          -k(25)*n(idx_SIC)  &
          -k(370)*n(idx_H2CO)  &
          -k(17)*n(idx_H2CO)  &
          -k(376)*n(idx_NH)  &
          -k(21)*n(idx_NO)  &
          -k(1194)*n(idx_O)  &
          -k(16)*n(idx_CH)  &
          -k(383)*n(idx_SIO)  &
          -k(1206)*n(idx_H)  &
          -k(379)*n(idx_OCN)  &
          -k(19)*n(idx_MG)  &
          -k(372)*n(idx_H2O)  &
          -k(18)*n(idx_HCO)  &
          -k(378)*n(idx_O2)  &
          -k(1200)*n(idx_H2)  &
          -k(381)*n(idx_SIH2)  &
          -k(1193)*n(idx_N)  &
          -k(371)*n(idx_H2O)  &
          -k(380)*n(idx_OH)  &
          -k(26)*n(idx_SIH2)  &
          -k(368)*n(idx_CO2)  &
          -k(20)*n(idx_NH3)  &
          -k(22)*n(idx_SI)  &
          -k(369)*n(idx_H2CO)  &
          -k(529)*n(idx_H2)  &
          -k(1265)  &
          -k(24)*n(idx_SIC3)  &
          -k(373)*n(idx_HCO)  &
          -k(375)*n(idx_NH3)  &
          -k(366)*n(idx_CH3OH)  &
          -k(367)*n(idx_CH3OH)  &
          -k(1215)*n(idx_E)
      pdj(74) =  &
          +k(15)*n(idx_CH2)  &
          +k(369)*n(idx_H2CO)  &
          +k(1200)*n(idx_H2)
      pdj(75) =  &
          +k(529)*n(idx_H2)  &
          +k(1206)*n(idx_H)  &
          +k(373)*n(idx_HCO)  &
          +k(16)*n(idx_CH)
      pdj(76) =  &
          +k(17)*n(idx_H2CO)
      pdj(77) =  &
          +k(19)*n(idx_MG)
      pdj(78) =  &
          +k(20)*n(idx_NH3)
      pdj(79) =  &
          +k(21)*n(idx_NO)
      pdj(80) =  &
          +k(383)*n(idx_SIO)  &
          +k(22)*n(idx_SI)
      pdj(81) =  &
          +k(23)*n(idx_SIC2)
      pdj(82) =  &
          +k(24)*n(idx_SIC3)
      pdj(83) =  &
          +k(381)*n(idx_SIH2)  &
          +k(382)*n(idx_SIH)  &
          +k(25)*n(idx_SIC)
      pdj(84) =  &
          +k(26)*n(idx_SIH2)
      pdj(85) =  &
          +k(27)*n(idx_SIH3)
      pdj(86) =  &
          +k(376)*n(idx_NH)  &
          +k(1193)*n(idx_N)
      pdj(87) =  &
          +k(379)*n(idx_OCN)  &
          +k(380)*n(idx_OH)  &
          +k(1194)*n(idx_O)  &
          +k(377)*n(idx_O2)  &
          +k(368)*n(idx_CO2)
      pdj(92) =  &
          +k(378)*n(idx_O2)
      pdj(94) =  &
          +k(367)*n(idx_CH3OH)
      pdj(97) =  &
          +k(374)*n(idx_NH2)  &
          +k(375)*n(idx_NH3)
      pdj(107) =  &
          +k(366)*n(idx_CH3OH)
    elseif(j==74) then
      pdj(1) =  &
          -k(298)*n(idx_E)  &
          -k(296)*n(idx_E)  &
          -k(297)*n(idx_E)
      pdj(2) =  &
          +k(1112)  &
          +k(298)*n(idx_E)
      pdj(3) =  &
          -k(422)*n(idx_O)
      pdj(6) =  &
          -k(531)*n(idx_H2)  &
          +k(1110)  &
          +k(616)*n(idx_H)  &
          +k(296)*n(idx_E)
      pdj(7) =  &
          +k(296)*n(idx_E)  &
          +k(297)*n(idx_E)
      pdj(8) =  &
          +2.d0*k(297)*n(idx_E)  &
          +k(419)*n(idx_H2O)  &
          +k(298)*n(idx_E)  &
          -k(616)*n(idx_H)  &
          +k(422)*n(idx_O)  &
          +k(737)*n(idx_N)  &
          +k(531)*n(idx_H2)  &
          +k(1111)
      pdj(9) =  &
          -k(419)*n(idx_H2O)
      pdj(10) =  &
          +k(421)*n(idx_O2)
      pdj(11) =  &
          -k(421)*n(idx_O2)
      pdj(12) =  &
          +k(37)*n(idx_NO)
      pdj(13) =  &
          -k(418)*n(idx_H2CO)
      pdj(14) =  &
          -k(420)*n(idx_HCO)
      pdj(17) =  &
          -k(37)*n(idx_NO)
      pdj(25) =  &
          +k(420)*n(idx_HCO)  &
          +k(417)*n(idx_CO2)
      pdj(28) =  &
          +k(418)*n(idx_H2CO)
      pdj(30) =  &
          -k(737)*n(idx_N)
      pdj(38) =  &
          -k(417)*n(idx_CO2)
      pdj(53) =  &
          +k(1279)
      pdj(70) =  &
          +k(418)*n(idx_H2CO)  &
          +k(422)*n(idx_O)  &
          +k(421)*n(idx_O2)
      pdj(71) =  &
          +k(1112)
      pdj(73) =  &
          +k(1110)
      pdj(74) =  &
          -k(417)*n(idx_CO2)  &
          -k(531)*n(idx_H2)  &
          -k(296)*n(idx_E)  &
          -k(1111)  &
          -k(298)*n(idx_E)  &
          -k(737)*n(idx_N)  &
          -k(1110)  &
          -k(297)*n(idx_E)  &
          -k(37)*n(idx_NO)  &
          -k(419)*n(idx_H2O)  &
          -k(420)*n(idx_HCO)  &
          -k(1112)  &
          -k(418)*n(idx_H2CO)  &
          -k(422)*n(idx_O)  &
          -k(616)*n(idx_H)  &
          -k(421)*n(idx_O2)  &
          -k(1279)
      pdj(75) =  &
          +k(616)*n(idx_H)  &
          +k(1111)
      pdj(76) =  &
          +k(417)*n(idx_CO2)
      pdj(79) =  &
          +k(37)*n(idx_NO)
      pdj(94) =  &
          +k(420)*n(idx_HCO)  &
          +k(531)*n(idx_H2)
      pdj(97) =  &
          +k(737)*n(idx_N)
      pdj(107) =  &
          +k(419)*n(idx_H2O)
    elseif(j==75) then
      pdj(1) =  &
          -k(295)*n(idx_E)
      pdj(2) =  &
          +k(33)*n(idx_MG)  &
          +k(36)*n(idx_SI)  &
          +k(34)*n(idx_NH3)  &
          +k(32)*n(idx_HCO)  &
          +k(35)*n(idx_NO)
      pdj(3) =  &
          -k(415)*n(idx_O)  &
          +k(413)*n(idx_O2)
      pdj(4) =  &
          -k(408)*n(idx_HNC)
      pdj(5) =  &
          -k(406)*n(idx_HCN)
      pdj(6) =  &
          +k(405)*n(idx_H2O)  &
          +k(411)*n(idx_NH)  &
          +k(410)*n(idx_NH2)  &
          +k(416)*n(idx_OH)  &
          +k(615)*n(idx_H)  &
          -k(530)*n(idx_H2)
      pdj(7) =  &
          +k(1109)  &
          +k(406)*n(idx_HCN)  &
          +k(295)*n(idx_E)  &
          +k(404)*n(idx_H2O)  &
          +k(401)*n(idx_H2CO)  &
          +k(408)*n(idx_HNC)
      pdj(8) =  &
          +k(403)*n(idx_H2O)  &
          +k(295)*n(idx_E)  &
          +k(242)  &
          -k(615)*n(idx_H)  &
          +k(409)*n(idx_N)  &
          +k(415)*n(idx_O)  &
          +k(530)*n(idx_H2)
      pdj(9) =  &
          -k(403)*n(idx_H2O)  &
          -k(405)*n(idx_H2O)  &
          -k(404)*n(idx_H2O)
      pdj(10) =  &
          +k(412)*n(idx_O2)  &
          -k(416)*n(idx_OH)
      pdj(11) =  &
          -k(414)*n(idx_O2)  &
          -k(412)*n(idx_O2)  &
          -k(413)*n(idx_O2)
      pdj(12) =  &
          +k(398)*n(idx_CH3OH)  &
          +k(402)*n(idx_H2CO)
      pdj(13) =  &
          +k(397)*n(idx_CH3OH)  &
          -k(401)*n(idx_H2CO)  &
          -k(400)*n(idx_H2CO)  &
          -k(402)*n(idx_H2CO)
      pdj(14) =  &
          -k(407)*n(idx_HCO)  &
          -k(32)*n(idx_HCO)  &
          +k(414)*n(idx_O2)
      pdj(15) =  &
          -k(33)*n(idx_MG)
      pdj(16) =  &
          -k(34)*n(idx_NH3)
      pdj(17) =  &
          -k(35)*n(idx_NO)
      pdj(18) =  &
          -k(36)*n(idx_SI)
      pdj(25) =  &
          +k(400)*n(idx_H2CO)  &
          +k(407)*n(idx_HCO)  &
          +k(399)*n(idx_CO2)
      pdj(27) =  &
          -k(410)*n(idx_NH2)
      pdj(30) =  &
          -k(409)*n(idx_N)
      pdj(31) =  &
          -k(411)*n(idx_NH)
      pdj(37) =  &
          -k(398)*n(idx_CH3OH)  &
          -k(397)*n(idx_CH3OH)
      pdj(38) =  &
          -k(399)*n(idx_CO2)
      pdj(53) =  &
          +k(1273)
      pdj(70) =  &
          +k(405)*n(idx_H2O)  &
          +k(402)*n(idx_H2CO)  &
          +k(32)*n(idx_HCO)  &
          +k(413)*n(idx_O2)  &
          +k(399)*n(idx_CO2)
      pdj(71) =  &
          +k(1109)
      pdj(73) =  &
          +k(615)*n(idx_H)  &
          +k(242)
      pdj(74) =  &
          +k(530)*n(idx_H2)  &
          +k(407)*n(idx_HCO)
      pdj(75) =  &
          -k(295)*n(idx_E)  &
          -k(404)*n(idx_H2O)  &
          -k(1273)  &
          -k(415)*n(idx_O)  &
          -k(615)*n(idx_H)  &
          -k(409)*n(idx_N)  &
          -k(397)*n(idx_CH3OH)  &
          -k(35)*n(idx_NO)  &
          -k(403)*n(idx_H2O)  &
          -k(411)*n(idx_NH)  &
          -k(402)*n(idx_H2CO)  &
          -k(408)*n(idx_HNC)  &
          -k(399)*n(idx_CO2)  &
          -k(36)*n(idx_SI)  &
          -k(412)*n(idx_O2)  &
          -k(407)*n(idx_HCO)  &
          -k(33)*n(idx_MG)  &
          -k(416)*n(idx_OH)  &
          -k(400)*n(idx_H2CO)  &
          -k(32)*n(idx_HCO)  &
          -k(414)*n(idx_O2)  &
          -k(530)*n(idx_H2)  &
          -k(405)*n(idx_H2O)  &
          -k(1109)  &
          -k(401)*n(idx_H2CO)  &
          -k(34)*n(idx_NH3)  &
          -k(406)*n(idx_HCN)  &
          -k(413)*n(idx_O2)  &
          -k(398)*n(idx_CH3OH)  &
          -k(242)  &
          -k(410)*n(idx_NH2)
      pdj(76) =  &
          +k(403)*n(idx_H2O)
      pdj(77) =  &
          +k(33)*n(idx_MG)
      pdj(78) =  &
          +k(34)*n(idx_NH3)
      pdj(79) =  &
          +k(35)*n(idx_NO)
      pdj(80) =  &
          +k(36)*n(idx_SI)
      pdj(86) =  &
          +k(411)*n(idx_NH)  &
          +k(409)*n(idx_N)
      pdj(87) =  &
          +k(415)*n(idx_O)  &
          +k(416)*n(idx_OH)  &
          +k(412)*n(idx_O2)
      pdj(92) =  &
          +k(414)*n(idx_O2)
      pdj(94) =  &
          +k(397)*n(idx_CH3OH)  &
          +k(400)*n(idx_H2CO)
      pdj(97) =  &
          +k(410)*n(idx_NH2)
      pdj(107) =  &
          +k(398)*n(idx_CH3OH)  &
          +k(401)*n(idx_H2CO)
      pdj(108) =  &
          +k(404)*n(idx_H2O)
      pdj(109) =  &
          +k(408)*n(idx_HNC)  &
          +k(406)*n(idx_HCN)
    elseif(j==76) then
      pdj(1) =  &
          -k(307)*n(idx_E)  &
          -k(308)*n(idx_E)  &
          -k(309)*n(idx_E)  &
          -k(310)*n(idx_E)  &
          -k(1218)*n(idx_E)
      pdj(2) =  &
          -k(56)*n(idx_CH)  &
          -k(460)*n(idx_CH)
      pdj(3) =  &
          +k(307)*n(idx_E)
      pdj(4) =  &
          -k(647)*n(idx_HNC)
      pdj(5) =  &
          -k(628)*n(idx_HCN)
      pdj(6) =  &
          +k(308)*n(idx_E)
      pdj(8) =  &
          +2.d0*k(309)*n(idx_E)  &
          +k(310)*n(idx_E)
      pdj(9) =  &
          -k(564)*n(idx_H2O)
      pdj(11) =  &
          -k(550)*n(idx_O2)
      pdj(12) =  &
          -k(424)*n(idx_CH2)  &
          +k(307)*n(idx_E)  &
          -k(40)*n(idx_CH2)
      pdj(13) =  &
          +k(195)*n(idx_NH3)  &
          +k(1218)*n(idx_E)  &
          +k(40)*n(idx_CH2)  &
          -k(549)*n(idx_H2CO)  &
          +k(56)*n(idx_CH)  &
          +k(149)*n(idx_MG)  &
          +k(229)*n(idx_SI)  &
          +k(204)*n(idx_NO)  &
          +k(137)*n(idx_HCO)
      pdj(14) =  &
          +k(424)*n(idx_CH2)  &
          +k(781)*n(idx_NH2)  &
          +k(460)*n(idx_CH)  &
          -k(642)*n(idx_HCO)  &
          +k(628)*n(idx_HCN)  &
          +k(549)*n(idx_H2CO)  &
          -k(137)*n(idx_HCO)  &
          +k(564)*n(idx_H2O)  &
          +k(310)*n(idx_E)  &
          +k(647)*n(idx_HNC)
      pdj(15) =  &
          -k(149)*n(idx_MG)
      pdj(16) =  &
          -k(195)*n(idx_NH3)
      pdj(17) =  &
          -k(204)*n(idx_NO)
      pdj(18) =  &
          -k(229)*n(idx_SI)
      pdj(25) =  &
          +k(308)*n(idx_E)  &
          +k(309)*n(idx_E)  &
          +k(642)*n(idx_HCO)
      pdj(27) =  &
          -k(781)*n(idx_NH2)
      pdj(28) =  &
          +k(453)*n(idx_CH4)
      pdj(29) =  &
          -k(453)*n(idx_CH4)
      pdj(30) =  &
          +k(797)*n(idx_NH)
      pdj(31) =  &
          -k(797)*n(idx_NH)
      pdj(43) =  &
          +k(550)*n(idx_O2)
      pdj(47) =  &
          +k(1285)
      pdj(70) =  &
          +k(550)*n(idx_O2)  &
          +k(137)*n(idx_HCO)
      pdj(74) =  &
          +k(40)*n(idx_CH2)  &
          +k(460)*n(idx_CH)
      pdj(75) =  &
          +k(56)*n(idx_CH)
      pdj(76) =  &
          -k(642)*n(idx_HCO)  &
          -k(549)*n(idx_H2CO)  &
          -k(453)*n(idx_CH4)  &
          -k(308)*n(idx_E)  &
          -k(424)*n(idx_CH2)  &
          -k(309)*n(idx_E)  &
          -k(310)*n(idx_E)  &
          -k(195)*n(idx_NH3)  &
          -k(40)*n(idx_CH2)  &
          -k(628)*n(idx_HCN)  &
          -k(797)*n(idx_NH)  &
          -k(781)*n(idx_NH2)  &
          -k(137)*n(idx_HCO)  &
          -k(204)*n(idx_NO)  &
          -k(149)*n(idx_MG)  &
          -k(564)*n(idx_H2O)  &
          -k(307)*n(idx_E)  &
          -k(229)*n(idx_SI)  &
          -k(550)*n(idx_O2)  &
          -k(56)*n(idx_CH)  &
          -k(647)*n(idx_HNC)  &
          -k(1285)  &
          -k(460)*n(idx_CH)  &
          -k(1218)*n(idx_E)
      pdj(77) =  &
          +k(149)*n(idx_MG)
      pdj(78) =  &
          +k(781)*n(idx_NH2)  &
          +k(195)*n(idx_NH3)
      pdj(79) =  &
          +k(204)*n(idx_NO)
      pdj(80) =  &
          +k(229)*n(idx_SI)
      pdj(94) =  &
          +k(424)*n(idx_CH2)
      pdj(107) =  &
          +k(797)*n(idx_NH)  &
          +k(642)*n(idx_HCO)  &
          +k(549)*n(idx_H2CO)  &
          +k(453)*n(idx_CH4)
      pdj(108) =  &
          +k(564)*n(idx_H2O)
      pdj(109) =  &
          +k(647)*n(idx_HNC)  &
          +k(628)*n(idx_HCN)
    elseif(j==77) then
      pdj(1) =  &
          -k(1220)*n(idx_E)
      pdj(15) =  &
          +k(1220)*n(idx_E)
      pdj(66) =  &
          +k(1300)
      pdj(77) =  &
          -k(1220)*n(idx_E)  &
          -k(1300)
    elseif(j==78) then
      pdj(1) =  &
          -k(344)*n(idx_E)  &
          -k(345)*n(idx_E)
      pdj(3) =  &
          -k(829)*n(idx_O)
      pdj(6) =  &
          +k(829)*n(idx_O)
      pdj(8) =  &
          +k(344)*n(idx_E)  &
          +2.d0*k(345)*n(idx_E)
      pdj(12) =  &
          -k(435)*n(idx_CH2)
      pdj(14) =  &
          -k(190)*n(idx_HCO)
      pdj(15) =  &
          -k(191)*n(idx_MG)
      pdj(16) =  &
          +k(193)*n(idx_SI)  &
          +k(190)*n(idx_HCO)  &
          +k(192)*n(idx_NO)  &
          +k(191)*n(idx_MG)
      pdj(17) =  &
          -k(192)*n(idx_NO)
      pdj(18) =  &
          -k(193)*n(idx_SI)
      pdj(27) =  &
          +k(344)*n(idx_E)  &
          +k(435)*n(idx_CH2)
      pdj(31) =  &
          +k(345)*n(idx_E)
      pdj(60) =  &
          +k(1284)
      pdj(70) =  &
          +k(190)*n(idx_HCO)
      pdj(77) =  &
          +k(191)*n(idx_MG)
      pdj(78) =  &
          -k(191)*n(idx_MG)  &
          -k(344)*n(idx_E)  &
          -k(435)*n(idx_CH2)  &
          -k(829)*n(idx_O)  &
          -k(192)*n(idx_NO)  &
          -k(345)*n(idx_E)  &
          -k(1284)  &
          -k(190)*n(idx_HCO)  &
          -k(193)*n(idx_SI)
      pdj(79) =  &
          +k(192)*n(idx_NO)
      pdj(80) =  &
          +k(193)*n(idx_SI)
      pdj(94) =  &
          +k(435)*n(idx_CH2)
      pdj(104) =  &
          +k(829)*n(idx_O)
    elseif(j==79) then
      pdj(1) =  &
          -k(346)*n(idx_E)
      pdj(3) =  &
          +k(346)*n(idx_E)
      pdj(15) =  &
          -k(152)*n(idx_MG)
      pdj(17) =  &
          +k(230)*n(idx_SI)  &
          +k(152)*n(idx_MG)
      pdj(18) =  &
          -k(230)*n(idx_SI)
      pdj(30) =  &
          +k(346)*n(idx_E)
      pdj(56) =  &
          +k(1278)
      pdj(77) =  &
          +k(152)*n(idx_MG)
      pdj(79) =  &
          -k(230)*n(idx_SI)  &
          -k(152)*n(idx_MG)  &
          -k(346)*n(idx_E)  &
          -k(1278)
      pdj(80) =  &
          +k(230)*n(idx_SI)
    elseif(j==80) then
      pdj(1) =  &
          -k(1223)*n(idx_E)
      pdj(2) =  &
          -k(477)*n(idx_CH)
      pdj(3) =  &
          -k(1213)*n(idx_O)
      pdj(6) =  &
          -k(1203)*n(idx_H2)
      pdj(8) =  &
          +k(573)*n(idx_H2O)  &
          +k(860)*n(idx_OH)  &
          +k(477)*n(idx_CH)  &
          -k(1210)*n(idx_H)
      pdj(9) =  &
          -k(573)*n(idx_H2O)
      pdj(10) =  &
          -k(860)*n(idx_OH)
      pdj(15) =  &
          -k(154)*n(idx_MG)
      pdj(18) =  &
          +k(154)*n(idx_MG)  &
          +k(1223)*n(idx_E)
      pdj(28) =  &
          +k(861)*n(idx_CH3OH)
      pdj(37) =  &
          -k(861)*n(idx_CH3OH)
      pdj(48) =  &
          +k(1232)
      pdj(77) =  &
          +k(154)*n(idx_MG)
      pdj(80) =  &
          -k(154)*n(idx_MG)  &
          -k(861)*n(idx_CH3OH)  &
          -k(477)*n(idx_CH)  &
          -k(1232)  &
          -k(860)*n(idx_OH)  &
          -k(1223)*n(idx_E)  &
          -k(1203)*n(idx_H2)  &
          -k(1213)*n(idx_O)  &
          -k(573)*n(idx_H2O)  &
          -k(1210)*n(idx_H)
      pdj(83) =  &
          +k(477)*n(idx_CH)
      pdj(84) =  &
          +k(1203)*n(idx_H2)
      pdj(100) =  &
          +k(1210)*n(idx_H)
      pdj(101) =  &
          +k(860)*n(idx_OH)  &
          +k(1213)*n(idx_O)
      pdj(115) =  &
          +k(861)*n(idx_CH3OH)  &
          +k(573)*n(idx_H2O)
    elseif(j==81) then
      pdj(1) =  &
          -k(351)*n(idx_E)
      pdj(7) =  &
          +k(351)*n(idx_E)
      pdj(21) =  &
          +k(351)*n(idx_E)
      pdj(51) =  &
          +k(1243)
      pdj(81) =  &
          -k(351)*n(idx_E)  &
          -k(1243)
    elseif(j==82) then
      pdj(1) =  &
          -k(352)*n(idx_E)
      pdj(7) =  &
          +k(352)*n(idx_E)
      pdj(19) =  &
          +k(352)*n(idx_E)
      pdj(52) =  &
          +k(1245)
      pdj(82) =  &
          -k(1245)  &
          -k(352)*n(idx_E)
    elseif(j==83) then
      pdj(1) =  &
          -k(350)*n(idx_E)
      pdj(3) =  &
          -k(832)*n(idx_O)
      pdj(7) =  &
          +k(832)*n(idx_O)  &
          +k(350)*n(idx_E)
      pdj(18) =  &
          +k(350)*n(idx_E)
      pdj(24) =  &
          +k(745)*n(idx_N)
      pdj(30) =  &
          -k(745)*n(idx_N)
      pdj(50) =  &
          +k(1242)
      pdj(80) =  &
          +k(745)*n(idx_N)
      pdj(83) =  &
          -k(832)*n(idx_O)  &
          -k(350)*n(idx_E)  &
          -k(745)*n(idx_N)  &
          -k(1242)
      pdj(101) =  &
          +k(832)*n(idx_O)
    elseif(j==84) then
      pdj(1) =  &
          -k(356)*n(idx_E)  &
          -k(355)*n(idx_E)  &
          -k(354)*n(idx_E)
      pdj(3) =  &
          -k(834)*n(idx_O)
      pdj(6) =  &
          +k(354)*n(idx_E)
      pdj(8) =  &
          +2.d0*k(355)*n(idx_E)  &
          +k(834)*n(idx_O)  &
          +k(356)*n(idx_E)
      pdj(10) =  &
          +k(863)*n(idx_O2)
      pdj(11) =  &
          -k(863)*n(idx_O2)
      pdj(18) =  &
          +k(355)*n(idx_E)  &
          +k(354)*n(idx_E)
      pdj(33) =  &
          +k(356)*n(idx_E)
      pdj(48) =  &
          +k(1235)
      pdj(84) =  &
          -k(355)*n(idx_E)  &
          -k(354)*n(idx_E)  &
          -k(356)*n(idx_E)  &
          -k(863)*n(idx_O2)  &
          -k(834)*n(idx_O)  &
          -k(1235)
      pdj(115) =  &
          +k(834)*n(idx_O)  &
          +k(863)*n(idx_O2)
    elseif(j==85) then
      pdj(1) =  &
          -k(358)*n(idx_E)  &
          -k(357)*n(idx_E)
      pdj(3) =  &
          -k(835)*n(idx_O)
      pdj(6) =  &
          +k(835)*n(idx_O)  &
          -k(1205)*n(idx_H2)  &
          +k(358)*n(idx_E)
      pdj(8) =  &
          +k(357)*n(idx_E)
      pdj(22) =  &
          +k(357)*n(idx_E)
      pdj(33) =  &
          +k(358)*n(idx_E)
      pdj(48) =  &
          +k(1237)
      pdj(85) =  &
          -k(1237)  &
          -k(358)*n(idx_E)  &
          -k(357)*n(idx_E)  &
          -k(1205)*n(idx_H2)  &
          -k(835)*n(idx_O)
      pdj(114) =  &
          +k(1205)*n(idx_H2)
      pdj(115) =  &
          +k(835)*n(idx_O)
    elseif(j==86) then
      pdj(1) =  &
          -k(304)*n(idx_E)
      pdj(2) =  &
          -k(54)*n(idx_CH)
      pdj(3) =  &
          -k(217)*n(idx_O)
      pdj(5) =  &
          +k(480)*n(idx_H2CO)  &
          -k(66)*n(idx_HCN)
      pdj(6) =  &
          -k(532)*n(idx_H2)
      pdj(7) =  &
          +k(738)*n(idx_N)  &
          -k(28)*n(idx_C)  &
          +k(304)*n(idx_E)
      pdj(8) =  &
          -k(127)*n(idx_H)  &
          +k(532)*n(idx_H2)
      pdj(9) =  &
          -k(561)*n(idx_H2O)  &
          -k(562)*n(idx_H2O)
      pdj(10) =  &
          +k(561)*n(idx_H2O)  &
          -k(226)*n(idx_OH)
      pdj(11) =  &
          -k(69)*n(idx_O2)  &
          -k(482)*n(idx_O2)
      pdj(12) =  &
          -k(38)*n(idx_CH2)
      pdj(13) =  &
          -k(65)*n(idx_H2CO)  &
          -k(480)*n(idx_H2CO)
      pdj(14) =  &
          -k(67)*n(idx_HCO)  &
          -k(481)*n(idx_HCO)
      pdj(17) =  &
          -k(68)*n(idx_NO)
      pdj(24) =  &
          +k(226)*n(idx_OH)  &
          +k(184)*n(idx_NH2)  &
          +k(28)*n(idx_C)  &
          +k(127)*n(idx_H)  &
          +k(200)*n(idx_NH)  &
          +k(54)*n(idx_CH)  &
          +k(68)*n(idx_NO)  &
          +k(67)*n(idx_HCO)  &
          +k(217)*n(idx_O)  &
          +k(69)*n(idx_O2)  &
          +k(65)*n(idx_H2CO)  &
          +k(66)*n(idx_HCN)  &
          +k(38)*n(idx_CH2)  &
          +k(64)*n(idx_CO)
      pdj(25) =  &
          -k(64)*n(idx_CO)  &
          +k(482)*n(idx_O2)  &
          +k(481)*n(idx_HCO)
      pdj(27) =  &
          -k(184)*n(idx_NH2)
      pdj(30) =  &
          +k(304)*n(idx_E)  &
          -k(738)*n(idx_N)
      pdj(31) =  &
          +k(562)*n(idx_H2O)  &
          -k(200)*n(idx_NH)
      pdj(59) =  &
          +k(1277)
      pdj(70) =  &
          +k(562)*n(idx_H2O)  &
          +k(480)*n(idx_H2CO)  &
          +k(67)*n(idx_HCO)
      pdj(71) =  &
          +k(127)*n(idx_H)
      pdj(73) =  &
          +k(28)*n(idx_C)
      pdj(74) =  &
          +k(38)*n(idx_CH2)
      pdj(75) =  &
          +k(54)*n(idx_CH)
      pdj(76) =  &
          +k(65)*n(idx_H2CO)
      pdj(79) =  &
          +k(482)*n(idx_O2)  &
          +k(68)*n(idx_NO)
      pdj(86) =  &
          -k(304)*n(idx_E)  &
          -k(67)*n(idx_HCO)  &
          -k(482)*n(idx_O2)  &
          -k(68)*n(idx_NO)  &
          -k(127)*n(idx_H)  &
          -k(65)*n(idx_H2CO)  &
          -k(38)*n(idx_CH2)  &
          -k(481)*n(idx_HCO)  &
          -k(66)*n(idx_HCN)  &
          -k(562)*n(idx_H2O)  &
          -k(69)*n(idx_O2)  &
          -k(561)*n(idx_H2O)  &
          -k(217)*n(idx_O)  &
          -k(64)*n(idx_CO)  &
          -k(738)*n(idx_N)  &
          -k(28)*n(idx_C)  &
          -k(184)*n(idx_NH2)  &
          -k(200)*n(idx_NH)  &
          -k(54)*n(idx_CH)  &
          -k(1277)  &
          -k(480)*n(idx_H2CO)  &
          -k(226)*n(idx_OH)  &
          -k(532)*n(idx_H2)
      pdj(87) =  &
          +k(64)*n(idx_CO)
      pdj(88) =  &
          +k(738)*n(idx_N)
      pdj(89) =  &
          +k(69)*n(idx_O2)
      pdj(91) =  &
          +k(184)*n(idx_NH2)
      pdj(92) =  &
          +k(217)*n(idx_O)
      pdj(93) =  &
          +k(226)*n(idx_OH)
      pdj(97) =  &
          +k(561)*n(idx_H2O)  &
          +k(66)*n(idx_HCN)  &
          +k(481)*n(idx_HCO)  &
          +k(532)*n(idx_H2)
      pdj(98) =  &
          +k(200)*n(idx_NH)
    elseif(j==87) then
      pdj(1) =  &
          -k(305)*n(idx_E)
      pdj(2) =  &
          +k(423)*n(idx_CH2)  &
          -k(55)*n(idx_CH)  &
          -k(459)*n(idx_CH)
      pdj(3) =  &
          +k(1132)  &
          -k(218)*n(idx_O)  &
          +k(305)*n(idx_E)  &
          +k(852)*n(idx_OH)
      pdj(5) =  &
          -k(135)*n(idx_HCN)
      pdj(6) =  &
          -k(534)*n(idx_H2)  &
          -k(533)*n(idx_H2)
      pdj(7) =  &
          +k(459)*n(idx_CH)  &
          +k(305)*n(idx_E)  &
          -k(29)*n(idx_C)
      pdj(8) =  &
          +k(533)*n(idx_H2)  &
          +k(534)*n(idx_H2)  &
          -k(128)*n(idx_H)
      pdj(9) =  &
          -k(124)*n(idx_H2O)  &
          -k(563)*n(idx_H2O)
      pdj(10) =  &
          -k(852)*n(idx_OH)  &
          -k(227)*n(idx_OH)  &
          +k(563)*n(idx_H2O)
      pdj(11) =  &
          -k(74)*n(idx_O2)
      pdj(12) =  &
          -k(39)*n(idx_CH2)  &
          -k(423)*n(idx_CH2)
      pdj(13) =  &
          -k(71)*n(idx_H2CO)  &
          -k(485)*n(idx_H2CO)
      pdj(14) =  &
          -k(72)*n(idx_HCO)  &
          +k(485)*n(idx_H2CO)
      pdj(16) =  &
          -k(194)*n(idx_NH3)  &
          -k(793)*n(idx_NH3)
      pdj(17) =  &
          -k(73)*n(idx_NO)
      pdj(25) =  &
          +k(72)*n(idx_HCO)  &
          +k(124)*n(idx_H2O)  &
          +k(185)*n(idx_NH2)  &
          +k(73)*n(idx_NO)  &
          +k(201)*n(idx_NH)  &
          +k(128)*n(idx_H)  &
          +k(53)*n(idx_CH4)  &
          +k(227)*n(idx_OH)  &
          +k(55)*n(idx_CH)  &
          +k(71)*n(idx_H2CO)  &
          +k(29)*n(idx_C)  &
          +k(39)*n(idx_CH2)  &
          +k(74)*n(idx_O2)  &
          +k(135)*n(idx_HCN)  &
          +k(218)*n(idx_O)  &
          +k(194)*n(idx_NH3)
      pdj(27) =  &
          -k(185)*n(idx_NH2)  &
          +k(793)*n(idx_NH3)  &
          -k(780)*n(idx_NH2)
      pdj(28) =  &
          +k(452)*n(idx_CH4)
      pdj(29) =  &
          -k(452)*n(idx_CH4)  &
          -k(53)*n(idx_CH4)
      pdj(30) =  &
          +k(796)*n(idx_NH)
      pdj(31) =  &
          +k(780)*n(idx_NH2)  &
          -k(201)*n(idx_NH)  &
          -k(796)*n(idx_NH)
      pdj(54) =  &
          +k(1276)
      pdj(70) =  &
          +k(423)*n(idx_CH2)  &
          +k(533)*n(idx_H2)  &
          +k(793)*n(idx_NH3)  &
          +k(796)*n(idx_NH)  &
          +k(563)*n(idx_H2O)  &
          +k(72)*n(idx_HCO)  &
          +k(780)*n(idx_NH2)  &
          +k(485)*n(idx_H2CO)  &
          +k(452)*n(idx_CH4)  &
          +k(852)*n(idx_OH)  &
          +k(459)*n(idx_CH)
      pdj(71) =  &
          +k(128)*n(idx_H)
      pdj(72) =  &
          +k(534)*n(idx_H2)
      pdj(73) =  &
          +k(1132)  &
          +k(29)*n(idx_C)
      pdj(74) =  &
          +k(39)*n(idx_CH2)
      pdj(75) =  &
          +k(55)*n(idx_CH)
      pdj(76) =  &
          +k(71)*n(idx_H2CO)
      pdj(78) =  &
          +k(194)*n(idx_NH3)
      pdj(79) =  &
          +k(73)*n(idx_NO)
      pdj(87) =  &
          -k(55)*n(idx_CH)  &
          -k(74)*n(idx_O2)  &
          -k(53)*n(idx_CH4)  &
          -k(305)*n(idx_E)  &
          -k(201)*n(idx_NH)  &
          -k(485)*n(idx_H2CO)  &
          -k(73)*n(idx_NO)  &
          -k(71)*n(idx_H2CO)  &
          -k(128)*n(idx_H)  &
          -k(534)*n(idx_H2)  &
          -k(1132)  &
          -k(29)*n(idx_C)  &
          -k(459)*n(idx_CH)  &
          -k(227)*n(idx_OH)  &
          -k(793)*n(idx_NH3)  &
          -k(780)*n(idx_NH2)  &
          -k(124)*n(idx_H2O)  &
          -k(194)*n(idx_NH3)  &
          -k(135)*n(idx_HCN)  &
          -k(796)*n(idx_NH)  &
          -k(185)*n(idx_NH2)  &
          -k(72)*n(idx_HCO)  &
          -k(1276)  &
          -k(852)*n(idx_OH)  &
          -k(39)*n(idx_CH2)  &
          -k(563)*n(idx_H2O)  &
          -k(218)*n(idx_O)  &
          -k(452)*n(idx_CH4)  &
          -k(423)*n(idx_CH2)  &
          -k(533)*n(idx_H2)
      pdj(89) =  &
          +k(74)*n(idx_O2)
      pdj(90) =  &
          +k(124)*n(idx_H2O)
      pdj(91) =  &
          +k(185)*n(idx_NH2)
      pdj(92) =  &
          +k(218)*n(idx_O)
      pdj(93) =  &
          +k(227)*n(idx_OH)
      pdj(95) =  &
          +k(53)*n(idx_CH4)
      pdj(97) =  &
          +k(135)*n(idx_HCN)
      pdj(98) =  &
          +k(201)*n(idx_NH)
    elseif(j==88) then
      pdj(1) =  &
          -k(338)*n(idx_E)
      pdj(2) =  &
          -k(59)*n(idx_CH)
      pdj(3) =  &
          -k(826)*n(idx_O)  &
          -k(219)*n(idx_O)
      pdj(5) =  &
          -k(136)*n(idx_HCN)
      pdj(6) =  &
          -k(540)*n(idx_H2)  &
          +k(456)*n(idx_CH4)
      pdj(7) =  &
          -k(30)*n(idx_C)
      pdj(8) =  &
          +k(457)*n(idx_CH4)  &
          +k(731)*n(idx_H2CO)  &
          +k(540)*n(idx_H2)
      pdj(9) =  &
          -k(126)*n(idx_H2O)  &
          -k(570)*n(idx_H2O)
      pdj(10) =  &
          +k(570)*n(idx_H2O)  &
          -k(228)*n(idx_OH)
      pdj(11) =  &
          -k(174)*n(idx_O2)
      pdj(12) =  &
          -k(42)*n(idx_CH2)
      pdj(13) =  &
          -k(171)*n(idx_H2CO)  &
          -k(731)*n(idx_H2CO)
      pdj(14) =  &
          -k(172)*n(idx_HCO)  &
          -k(732)*n(idx_HCO)
      pdj(15) =  &
          -k(151)*n(idx_MG)
      pdj(16) =  &
          -k(198)*n(idx_NH3)
      pdj(17) =  &
          -k(173)*n(idx_NO)
      pdj(24) =  &
          -k(70)*n(idx_CN)
      pdj(25) =  &
          +k(732)*n(idx_HCO)  &
          -k(75)*n(idx_CO)
      pdj(26) =  &
          +k(75)*n(idx_CO)  &
          +k(219)*n(idx_O)  &
          +k(174)*n(idx_O2)  &
          +k(228)*n(idx_OH)  &
          +k(456)*n(idx_CH4)  &
          +k(30)*n(idx_C)  &
          +k(187)*n(idx_NH2)  &
          +k(175)*n(idx_N)  &
          +k(172)*n(idx_HCO)  &
          +k(136)*n(idx_HCN)  &
          +k(151)*n(idx_MG)  &
          +k(171)*n(idx_H2CO)  &
          +k(198)*n(idx_NH3)  &
          +k(42)*n(idx_CH2)  &
          +k(173)*n(idx_NO)  &
          +k(202)*n(idx_NH)  &
          +k(70)*n(idx_CN)  &
          +k(126)*n(idx_H2O)  &
          +k(731)*n(idx_H2CO)  &
          +k(59)*n(idx_CH)  &
          +k(457)*n(idx_CH4)
      pdj(27) =  &
          -k(187)*n(idx_NH2)
      pdj(29) =  &
          -k(456)*n(idx_CH4)  &
          -k(457)*n(idx_CH4)
      pdj(30) =  &
          +2.d0*k(338)*n(idx_E)  &
          -k(175)*n(idx_N)  &
          +k(826)*n(idx_O)
      pdj(31) =  &
          -k(202)*n(idx_NH)
      pdj(58) =  &
          +k(1272)
      pdj(70) =  &
          +k(172)*n(idx_HCO)  &
          +k(731)*n(idx_H2CO)
      pdj(73) =  &
          +k(30)*n(idx_C)
      pdj(74) =  &
          +k(42)*n(idx_CH2)  &
          +k(456)*n(idx_CH4)
      pdj(75) =  &
          +k(59)*n(idx_CH)
      pdj(76) =  &
          +k(171)*n(idx_H2CO)
      pdj(77) =  &
          +k(151)*n(idx_MG)
      pdj(78) =  &
          +k(198)*n(idx_NH3)
      pdj(79) =  &
          +k(826)*n(idx_O)  &
          +k(173)*n(idx_NO)
      pdj(86) =  &
          +k(70)*n(idx_CN)
      pdj(87) =  &
          +k(75)*n(idx_CO)
      pdj(88) =  &
          -k(187)*n(idx_NH2)  &
          -k(172)*n(idx_HCO)  &
          -k(59)*n(idx_CH)  &
          -k(173)*n(idx_NO)  &
          -k(198)*n(idx_NH3)  &
          -k(75)*n(idx_CO)  &
          -k(30)*n(idx_C)  &
          -k(338)*n(idx_E)  &
          -k(202)*n(idx_NH)  &
          -k(570)*n(idx_H2O)  &
          -k(70)*n(idx_CN)  &
          -k(457)*n(idx_CH4)  &
          -k(219)*n(idx_O)  &
          -k(151)*n(idx_MG)  &
          -k(175)*n(idx_N)  &
          -k(1272)  &
          -k(228)*n(idx_OH)  &
          -k(126)*n(idx_H2O)  &
          -k(540)*n(idx_H2)  &
          -k(731)*n(idx_H2CO)  &
          -k(826)*n(idx_O)  &
          -k(732)*n(idx_HCO)  &
          -k(174)*n(idx_O2)  &
          -k(136)*n(idx_HCN)  &
          -k(456)*n(idx_CH4)  &
          -k(171)*n(idx_H2CO)  &
          -k(42)*n(idx_CH2)
      pdj(89) =  &
          +k(174)*n(idx_O2)
      pdj(90) =  &
          +k(126)*n(idx_H2O)
      pdj(91) =  &
          +k(187)*n(idx_NH2)
      pdj(92) =  &
          +k(219)*n(idx_O)
      pdj(93) =  &
          +k(228)*n(idx_OH)
      pdj(94) =  &
          +k(457)*n(idx_CH4)
      pdj(96) =  &
          +k(175)*n(idx_N)
      pdj(97) =  &
          +k(136)*n(idx_HCN)
      pdj(98) =  &
          +k(202)*n(idx_NH)
      pdj(112) =  &
          +k(732)*n(idx_HCO)  &
          +k(570)*n(idx_H2O)  &
          +k(540)*n(idx_H2)
    elseif(j==89) then
      pdj(1) =  &
          -k(347)*n(idx_E)
      pdj(2) =  &
          -k(62)*n(idx_CH)  &
          -k(474)*n(idx_CH)
      pdj(3) =  &
          +2.d0*k(347)*n(idx_E)  &
          +k(436)*n(idx_CH2)  &
          +k(474)*n(idx_CH)  &
          +k(1168)  &
          +k(392)*n(idx_C)  &
          +k(743)*n(idx_N)  &
          +k(805)*n(idx_NH)
      pdj(7) =  &
          -k(31)*n(idx_C)  &
          -k(392)*n(idx_C)
      pdj(8) =  &
          +k(821)*n(idx_CH3OH)  &
          +k(552)*n(idx_H2CO)
      pdj(11) =  &
          +k(138)*n(idx_HCO)  &
          +k(31)*n(idx_C)  &
          +k(153)*n(idx_MG)  &
          +k(188)*n(idx_NH2)  &
          +k(62)*n(idx_CH)  &
          +k(117)*n(idx_H2CO)  &
          +k(206)*n(idx_NO)  &
          +k(231)*n(idx_SI)  &
          +k(199)*n(idx_NH3)  &
          +k(821)*n(idx_CH3OH)  &
          +k(552)*n(idx_H2CO)  &
          +k(45)*n(idx_CH2)
      pdj(12) =  &
          -k(45)*n(idx_CH2)  &
          -k(436)*n(idx_CH2)
      pdj(13) =  &
          -k(552)*n(idx_H2CO)  &
          -k(117)*n(idx_H2CO)
      pdj(14) =  &
          -k(645)*n(idx_HCO)  &
          -k(138)*n(idx_HCO)
      pdj(15) =  &
          -k(153)*n(idx_MG)
      pdj(16) =  &
          -k(199)*n(idx_NH3)
      pdj(17) =  &
          -k(206)*n(idx_NO)
      pdj(18) =  &
          -k(231)*n(idx_SI)
      pdj(25) =  &
          +k(645)*n(idx_HCO)
      pdj(27) =  &
          -k(188)*n(idx_NH2)
      pdj(30) =  &
          -k(743)*n(idx_N)
      pdj(31) =  &
          -k(805)*n(idx_NH)
      pdj(37) =  &
          -k(821)*n(idx_CH3OH)
      pdj(61) =  &
          +k(1271)
      pdj(70) =  &
          +k(138)*n(idx_HCO)  &
          +k(474)*n(idx_CH)  &
          +k(552)*n(idx_H2CO)
      pdj(73) =  &
          +k(31)*n(idx_C)
      pdj(74) =  &
          +k(45)*n(idx_CH2)
      pdj(75) =  &
          +k(62)*n(idx_CH)
      pdj(76) =  &
          +k(436)*n(idx_CH2)  &
          +k(117)*n(idx_H2CO)
      pdj(77) =  &
          +k(153)*n(idx_MG)
      pdj(78) =  &
          +k(199)*n(idx_NH3)
      pdj(79) =  &
          +k(206)*n(idx_NO)  &
          +k(743)*n(idx_N)
      pdj(80) =  &
          +k(231)*n(idx_SI)
      pdj(87) =  &
          +k(392)*n(idx_C)
      pdj(89) =  &
          -k(31)*n(idx_C)  &
          -k(805)*n(idx_NH)  &
          -k(1271)  &
          -k(62)*n(idx_CH)  &
          -k(117)*n(idx_H2CO)  &
          -k(347)*n(idx_E)  &
          -k(153)*n(idx_MG)  &
          -k(436)*n(idx_CH2)  &
          -k(45)*n(idx_CH2)  &
          -k(206)*n(idx_NO)  &
          -k(645)*n(idx_HCO)  &
          -k(743)*n(idx_N)  &
          -k(231)*n(idx_SI)  &
          -k(552)*n(idx_H2CO)  &
          -k(199)*n(idx_NH3)  &
          -k(138)*n(idx_HCO)  &
          -k(392)*n(idx_C)  &
          -k(188)*n(idx_NH2)  &
          -k(821)*n(idx_CH3OH)  &
          -k(474)*n(idx_CH)  &
          -k(1168)
      pdj(91) =  &
          +k(188)*n(idx_NH2)
      pdj(92) =  &
          +k(1168)
      pdj(104) =  &
          +k(805)*n(idx_NH)
      pdj(107) =  &
          +k(821)*n(idx_CH3OH)
      pdj(113) =  &
          +k(645)*n(idx_HCO)
    elseif(j==90) then
      pdj(1) =  &
          -k(315)*n(idx_E)  &
          -k(314)*n(idx_E)  &
          -k(313)*n(idx_E)
      pdj(2) =  &
          -k(57)*n(idx_CH)  &
          -k(461)*n(idx_CH)
      pdj(3) =  &
          +k(853)*n(idx_OH)  &
          +k(313)*n(idx_E)  &
          -k(824)*n(idx_O)  &
          +k(314)*n(idx_E)
      pdj(4) =  &
          -k(560)*n(idx_HNC)
      pdj(5) =  &
          -k(557)*n(idx_HCN)
      pdj(6) =  &
          -k(535)*n(idx_H2)  &
          +k(824)*n(idx_O)  &
          +k(313)*n(idx_E)  &
          +k(740)*n(idx_N)
      pdj(7) =  &
          -k(384)*n(idx_C)
      pdj(8) =  &
          +k(1141)  &
          +k(535)*n(idx_H2)  &
          +k(739)*n(idx_N)  &
          +k(315)*n(idx_E)  &
          +2.d0*k(314)*n(idx_E)
      pdj(9) =  &
          +k(186)*n(idx_NH2)  &
          +k(118)*n(idx_H2CO)  &
          +k(120)*n(idx_MG)  &
          -k(556)*n(idx_H2O)  &
          +k(57)*n(idx_CH)  &
          +k(196)*n(idx_NH3)  &
          +k(41)*n(idx_CH2)  &
          +k(119)*n(idx_HCO)  &
          +k(121)*n(idx_NO)  &
          +k(122)*n(idx_O2)  &
          +k(123)*n(idx_SI)
      pdj(10) =  &
          +k(384)*n(idx_C)  &
          +k(782)*n(idx_NH2)  &
          +k(425)*n(idx_CH2)  &
          -k(853)*n(idx_OH)  &
          +k(315)*n(idx_E)  &
          +k(557)*n(idx_HCN)  &
          +k(554)*n(idx_CO)  &
          +k(560)*n(idx_HNC)  &
          +k(556)*n(idx_H2O)  &
          +k(461)*n(idx_CH)  &
          +k(555)*n(idx_H2CO)  &
          +k(559)*n(idx_HCO)
      pdj(11) =  &
          -k(122)*n(idx_O2)
      pdj(12) =  &
          -k(425)*n(idx_CH2)  &
          -k(41)*n(idx_CH2)
      pdj(13) =  &
          -k(118)*n(idx_H2CO)  &
          -k(555)*n(idx_H2CO)
      pdj(14) =  &
          -k(119)*n(idx_HCO)  &
          -k(558)*n(idx_HCO)  &
          -k(559)*n(idx_HCO)
      pdj(15) =  &
          -k(120)*n(idx_MG)
      pdj(16) =  &
          -k(196)*n(idx_NH3)
      pdj(17) =  &
          -k(121)*n(idx_NO)
      pdj(18) =  &
          -k(123)*n(idx_SI)
      pdj(25) =  &
          +k(558)*n(idx_HCO)  &
          -k(554)*n(idx_CO)
      pdj(27) =  &
          -k(782)*n(idx_NH2)  &
          -k(186)*n(idx_NH2)
      pdj(28) =  &
          +k(454)*n(idx_CH4)
      pdj(29) =  &
          -k(454)*n(idx_CH4)
      pdj(30) =  &
          +k(798)*n(idx_NH)  &
          -k(740)*n(idx_N)  &
          -k(739)*n(idx_N)
      pdj(31) =  &
          -k(798)*n(idx_NH)
      pdj(55) =  &
          +k(1281)
      pdj(70) =  &
          +k(554)*n(idx_CO)  &
          +k(119)*n(idx_HCO)
      pdj(74) =  &
          +k(461)*n(idx_CH)  &
          +k(41)*n(idx_CH2)
      pdj(75) =  &
          +k(384)*n(idx_C)  &
          +k(57)*n(idx_CH)
      pdj(76) =  &
          +k(118)*n(idx_H2CO)  &
          +k(559)*n(idx_HCO)
      pdj(77) =  &
          +k(120)*n(idx_MG)
      pdj(78) =  &
          +k(196)*n(idx_NH3)  &
          +k(782)*n(idx_NH2)
      pdj(79) =  &
          +k(740)*n(idx_N)  &
          +k(121)*n(idx_NO)
      pdj(80) =  &
          +k(123)*n(idx_SI)
      pdj(89) =  &
          +k(122)*n(idx_O2)  &
          +k(824)*n(idx_O)
      pdj(90) =  &
          -k(1141)  &
          -k(57)*n(idx_CH)  &
          -k(535)*n(idx_H2)  &
          -k(123)*n(idx_SI)  &
          -k(122)*n(idx_O2)  &
          -k(425)*n(idx_CH2)  &
          -k(554)*n(idx_CO)  &
          -k(558)*n(idx_HCO)  &
          -k(853)*n(idx_OH)  &
          -k(461)*n(idx_CH)  &
          -k(314)*n(idx_E)  &
          -k(313)*n(idx_E)  &
          -k(560)*n(idx_HNC)  &
          -k(186)*n(idx_NH2)  &
          -k(782)*n(idx_NH2)  &
          -k(196)*n(idx_NH3)  &
          -k(1281)  &
          -k(315)*n(idx_E)  &
          -k(118)*n(idx_H2CO)  &
          -k(740)*n(idx_N)  &
          -k(798)*n(idx_NH)  &
          -k(555)*n(idx_H2CO)  &
          -k(384)*n(idx_C)  &
          -k(119)*n(idx_HCO)  &
          -k(41)*n(idx_CH2)  &
          -k(454)*n(idx_CH4)  &
          -k(120)*n(idx_MG)  &
          -k(557)*n(idx_HCN)  &
          -k(824)*n(idx_O)  &
          -k(121)*n(idx_NO)  &
          -k(559)*n(idx_HCO)  &
          -k(556)*n(idx_H2O)  &
          -k(739)*n(idx_N)
      pdj(91) =  &
          +k(186)*n(idx_NH2)
      pdj(93) =  &
          +k(1141)
      pdj(94) =  &
          +k(425)*n(idx_CH2)
      pdj(104) =  &
          +k(739)*n(idx_N)
      pdj(107) =  &
          +k(555)*n(idx_H2CO)
      pdj(108) =  &
          +k(853)*n(idx_OH)  &
          +k(558)*n(idx_HCO)  &
          +k(798)*n(idx_NH)  &
          +k(535)*n(idx_H2)  &
          +k(454)*n(idx_CH4)  &
          +k(556)*n(idx_H2O)
      pdj(109) =  &
          +k(557)*n(idx_HCN)  &
          +k(560)*n(idx_HNC)
    elseif(j==91) then
      pdj(1) =  &
          -k(342)*n(idx_E)  &
          -k(343)*n(idx_E)
      pdj(2) =  &
          -k(60)*n(idx_CH)  &
          -k(472)*n(idx_CH)
      pdj(3) =  &
          -k(828)*n(idx_O)  &
          +k(778)*n(idx_O2)
      pdj(4) =  &
          -k(776)*n(idx_HNC)
      pdj(5) =  &
          -k(774)*n(idx_HCN)
      pdj(6) =  &
          -k(543)*n(idx_H2)
      pdj(8) =  &
          +2.d0*k(342)*n(idx_E)  &
          +k(343)*n(idx_E)  &
          +k(828)*n(idx_O)  &
          +k(742)*n(idx_N)  &
          +k(543)*n(idx_H2)
      pdj(9) =  &
          -k(772)*n(idx_H2O)  &
          -k(773)*n(idx_H2O)
      pdj(10) =  &
          +k(773)*n(idx_H2O)  &
          +k(779)*n(idx_O2)
      pdj(11) =  &
          -k(779)*n(idx_O2)  &
          -k(778)*n(idx_O2)
      pdj(12) =  &
          -k(434)*n(idx_CH2)  &
          -k(43)*n(idx_CH2)
      pdj(13) =  &
          -k(771)*n(idx_H2CO)  &
          -k(770)*n(idx_H2CO)
      pdj(14) =  &
          +k(771)*n(idx_H2CO)  &
          -k(775)*n(idx_HCO)  &
          -k(181)*n(idx_HCO)
      pdj(16) =  &
          -k(182)*n(idx_NH3)
      pdj(17) =  &
          -k(183)*n(idx_NO)
      pdj(27) =  &
          -k(777)*n(idx_NH2)  &
          +k(60)*n(idx_CH)  &
          +k(181)*n(idx_HCO)  &
          +k(43)*n(idx_CH2)  &
          +k(182)*n(idx_NH3)  &
          +k(183)*n(idx_NO)
      pdj(30) =  &
          +k(803)*n(idx_NH)  &
          +k(342)*n(idx_E)  &
          -k(742)*n(idx_N)
      pdj(31) =  &
          +k(774)*n(idx_HCN)  &
          +k(777)*n(idx_NH2)  &
          +k(343)*n(idx_E)  &
          +k(434)*n(idx_CH2)  &
          +k(775)*n(idx_HCO)  &
          +k(776)*n(idx_HNC)  &
          +k(472)*n(idx_CH)  &
          -k(803)*n(idx_NH)  &
          +k(770)*n(idx_H2CO)  &
          +k(772)*n(idx_H2O)
      pdj(60) =  &
          +k(1280)
      pdj(70) =  &
          +k(181)*n(idx_HCO)
      pdj(74) =  &
          +k(472)*n(idx_CH)  &
          +k(43)*n(idx_CH2)
      pdj(75) =  &
          +k(60)*n(idx_CH)
      pdj(76) =  &
          +k(775)*n(idx_HCO)
      pdj(78) =  &
          +k(543)*n(idx_H2)  &
          +k(803)*n(idx_NH)  &
          +k(773)*n(idx_H2O)  &
          +k(182)*n(idx_NH3)  &
          +k(771)*n(idx_H2CO)  &
          +k(777)*n(idx_NH2)
      pdj(79) =  &
          +k(183)*n(idx_NO)
      pdj(91) =  &
          -k(183)*n(idx_NO)  &
          -k(1280)  &
          -k(343)*n(idx_E)  &
          -k(828)*n(idx_O)  &
          -k(43)*n(idx_CH2)  &
          -k(772)*n(idx_H2O)  &
          -k(771)*n(idx_H2CO)  &
          -k(742)*n(idx_N)  &
          -k(776)*n(idx_HNC)  &
          -k(803)*n(idx_NH)  &
          -k(778)*n(idx_O2)  &
          -k(342)*n(idx_E)  &
          -k(60)*n(idx_CH)  &
          -k(775)*n(idx_HCO)  &
          -k(472)*n(idx_CH)  &
          -k(774)*n(idx_HCN)  &
          -k(181)*n(idx_HCO)  &
          -k(770)*n(idx_H2CO)  &
          -k(777)*n(idx_NH2)  &
          -k(182)*n(idx_NH3)  &
          -k(434)*n(idx_CH2)  &
          -k(773)*n(idx_H2O)  &
          -k(543)*n(idx_H2)  &
          -k(779)*n(idx_O2)
      pdj(94) =  &
          +k(434)*n(idx_CH2)
      pdj(104) =  &
          +k(828)*n(idx_O)  &
          +k(779)*n(idx_O2)
      pdj(105) =  &
          +k(778)*n(idx_O2)
      pdj(107) =  &
          +k(770)*n(idx_H2CO)
      pdj(108) =  &
          +k(772)*n(idx_H2O)
      pdj(109) =  &
          +k(774)*n(idx_HCN)  &
          +k(776)*n(idx_HNC)
      pdj(112) =  &
          +k(742)*n(idx_N)
    elseif(j==92) then
      pdj(1) =  &
          -k(1222)*n(idx_E)
      pdj(2) =  &
          -k(473)*n(idx_CH)  &
          +k(816)*n(idx_HCN)  &
          -k(61)*n(idx_CH)
      pdj(3) =  &
          +k(209)*n(idx_CO)  &
          +k(214)*n(idx_NH3)  &
          +k(211)*n(idx_H2O)  &
          +k(203)*n(idx_NH)  &
          +k(212)*n(idx_HCO)  &
          +k(132)*n(idx_H)  &
          +k(1222)*n(idx_E)  &
          +k(216)*n(idx_OH)  &
          +k(213)*n(idx_NH2)  &
          +k(208)*n(idx_CH4)  &
          +k(44)*n(idx_CH2)  &
          +k(61)*n(idx_CH)  &
          +k(210)*n(idx_H2CO)  &
          +k(215)*n(idx_O2)
      pdj(5) =  &
          -k(816)*n(idx_HCN)  &
          -k(815)*n(idx_HCN)
      pdj(6) =  &
          -k(544)*n(idx_H2)
      pdj(7) =  &
          +k(812)*n(idx_CN)  &
          -k(1196)*n(idx_C)
      pdj(8) =  &
          -k(132)*n(idx_H)  &
          +k(473)*n(idx_CH)  &
          +k(544)*n(idx_H2)  &
          +k(820)*n(idx_OH)  &
          +k(804)*n(idx_NH)
      pdj(9) =  &
          +k(809)*n(idx_CH3OH)  &
          -k(211)*n(idx_H2O)
      pdj(10) =  &
          +k(814)*n(idx_H2CO)  &
          +k(811)*n(idx_CH4)  &
          -k(820)*n(idx_OH)  &
          -k(216)*n(idx_OH)  &
          +k(810)*n(idx_CH3OH)
      pdj(11) =  &
          -k(215)*n(idx_O2)  &
          +k(819)*n(idx_NO2)
      pdj(12) =  &
          -k(44)*n(idx_CH2)
      pdj(13) =  &
          -k(210)*n(idx_H2CO)  &
          -k(814)*n(idx_H2CO)
      pdj(14) =  &
          -k(212)*n(idx_HCO)  &
          -k(817)*n(idx_HCO)
      pdj(16) =  &
          -k(214)*n(idx_NH3)
      pdj(24) =  &
          -k(812)*n(idx_CN)
      pdj(25) =  &
          +k(817)*n(idx_HCO)  &
          +k(813)*n(idx_CO2)  &
          -k(209)*n(idx_CO)
      pdj(26) =  &
          -k(818)*n(idx_N2)
      pdj(27) =  &
          -k(213)*n(idx_NH2)
      pdj(29) =  &
          -k(811)*n(idx_CH4)  &
          -k(208)*n(idx_CH4)
      pdj(30) =  &
          +k(818)*n(idx_N2)  &
          +k(815)*n(idx_HCN)
      pdj(31) =  &
          -k(804)*n(idx_NH)  &
          -k(203)*n(idx_NH)
      pdj(37) =  &
          -k(809)*n(idx_CH3OH)  &
          -k(810)*n(idx_CH3OH)
      pdj(38) =  &
          -k(813)*n(idx_CO2)
      pdj(42) =  &
          -k(819)*n(idx_NO2)
      pdj(55) =  &
          +k(1270)
      pdj(70) =  &
          +k(814)*n(idx_H2CO)  &
          +k(815)*n(idx_HCN)  &
          +k(212)*n(idx_HCO)
      pdj(71) =  &
          +k(132)*n(idx_H)
      pdj(74) =  &
          +k(44)*n(idx_CH2)
      pdj(75) =  &
          +k(61)*n(idx_CH)
      pdj(76) =  &
          +k(210)*n(idx_H2CO)  &
          +k(809)*n(idx_CH3OH)
      pdj(78) =  &
          +k(214)*n(idx_NH3)
      pdj(79) =  &
          +k(812)*n(idx_CN)  &
          +k(818)*n(idx_N2)  &
          +k(804)*n(idx_NH)  &
          +k(816)*n(idx_HCN)  &
          +k(819)*n(idx_NO2)
      pdj(87) =  &
          +k(209)*n(idx_CO)  &
          +k(473)*n(idx_CH)  &
          +k(1196)*n(idx_C)
      pdj(89) =  &
          +k(215)*n(idx_O2)  &
          +k(820)*n(idx_OH)  &
          +k(813)*n(idx_CO2)
      pdj(90) =  &
          +k(211)*n(idx_H2O)
      pdj(91) =  &
          +k(213)*n(idx_NH2)
      pdj(92) =  &
          -k(820)*n(idx_OH)  &
          -k(804)*n(idx_NH)  &
          -k(203)*n(idx_NH)  &
          -k(1270)  &
          -k(213)*n(idx_NH2)  &
          -k(816)*n(idx_HCN)  &
          -k(216)*n(idx_OH)  &
          -k(819)*n(idx_NO2)  &
          -k(815)*n(idx_HCN)  &
          -k(813)*n(idx_CO2)  &
          -k(810)*n(idx_CH3OH)  &
          -k(132)*n(idx_H)  &
          -k(812)*n(idx_CN)  &
          -k(212)*n(idx_HCO)  &
          -k(814)*n(idx_H2CO)  &
          -k(818)*n(idx_N2)  &
          -k(61)*n(idx_CH)  &
          -k(210)*n(idx_H2CO)  &
          -k(214)*n(idx_NH3)  &
          -k(544)*n(idx_H2)  &
          -k(473)*n(idx_CH)  &
          -k(209)*n(idx_CO)  &
          -k(809)*n(idx_CH3OH)  &
          -k(817)*n(idx_HCO)  &
          -k(44)*n(idx_CH2)  &
          -k(811)*n(idx_CH4)  &
          -k(1222)*n(idx_E)  &
          -k(215)*n(idx_O2)  &
          -k(208)*n(idx_CH4)  &
          -k(211)*n(idx_H2O)  &
          -k(1196)*n(idx_C)
      pdj(93) =  &
          +k(544)*n(idx_H2)  &
          +k(817)*n(idx_HCO)  &
          +k(216)*n(idx_OH)
      pdj(94) =  &
          +k(811)*n(idx_CH4)
      pdj(95) =  &
          +k(208)*n(idx_CH4)
      pdj(98) =  &
          +k(203)*n(idx_NH)
      pdj(107) =  &
          +k(810)*n(idx_CH3OH)
    elseif(j==93) then
      pdj(1) =  &
          -k(349)*n(idx_E)
      pdj(2) =  &
          -k(63)*n(idx_CH)  &
          -k(476)*n(idx_CH)
      pdj(3) =  &
          +k(845)*n(idx_HNC)  &
          +k(840)*n(idx_H2CO)  &
          -k(831)*n(idx_O)  &
          +k(476)*n(idx_CH)  &
          +k(846)*n(idx_N2)  &
          +k(792)*n(idx_NH2)  &
          +k(841)*n(idx_H2O)  &
          +k(848)*n(idx_OH)  &
          +k(349)*n(idx_E)  &
          +k(849)*n(idx_SI)  &
          +k(842)*n(idx_HCN)  &
          +k(837)*n(idx_CN)  &
          +k(847)*n(idx_NO)  &
          +k(851)*n(idx_SIO)  &
          +k(438)*n(idx_CH2)  &
          +k(844)*n(idx_HCO)  &
          +k(850)*n(idx_SIH)  &
          +k(838)*n(idx_CO2)  &
          +k(807)*n(idx_NH)  &
          +k(839)*n(idx_CO)  &
          +k(394)*n(idx_C)
      pdj(4) =  &
          -k(845)*n(idx_HNC)
      pdj(5) =  &
          -k(842)*n(idx_HCN)
      pdj(6) =  &
          -k(546)*n(idx_H2)
      pdj(7) =  &
          -k(394)*n(idx_C)
      pdj(8) =  &
          +k(831)*n(idx_O)  &
          +k(1174)  &
          +k(546)*n(idx_H2)  &
          +k(349)*n(idx_E)  &
          +k(744)*n(idx_N)
      pdj(9) =  &
          -k(841)*n(idx_H2O)  &
          -k(221)*n(idx_H2O)
      pdj(10) =  &
          +k(223)*n(idx_NH3)  &
          +k(220)*n(idx_H2CO)  &
          +k(224)*n(idx_NO)  &
          +k(222)*n(idx_HCO)  &
          -k(848)*n(idx_OH)  &
          +k(46)*n(idx_CH2)  &
          +k(221)*n(idx_H2O)  &
          +k(189)*n(idx_NH2)  &
          +k(225)*n(idx_O2)  &
          +k(63)*n(idx_CH)
      pdj(11) =  &
          -k(225)*n(idx_O2)
      pdj(12) =  &
          +k(458)*n(idx_CH4)  &
          -k(438)*n(idx_CH2)  &
          -k(46)*n(idx_CH2)
      pdj(13) =  &
          -k(220)*n(idx_H2CO)  &
          -k(840)*n(idx_H2CO)
      pdj(14) =  &
          -k(843)*n(idx_HCO)  &
          -k(222)*n(idx_HCO)  &
          -k(844)*n(idx_HCO)
      pdj(16) =  &
          -k(223)*n(idx_NH3)
      pdj(17) =  &
          -k(847)*n(idx_NO)  &
          -k(224)*n(idx_NO)
      pdj(18) =  &
          -k(849)*n(idx_SI)
      pdj(24) =  &
          -k(837)*n(idx_CN)
      pdj(25) =  &
          -k(839)*n(idx_CO)  &
          +k(843)*n(idx_HCO)
      pdj(26) =  &
          -k(846)*n(idx_N2)
      pdj(27) =  &
          -k(792)*n(idx_NH2)  &
          -k(189)*n(idx_NH2)
      pdj(29) =  &
          -k(458)*n(idx_CH4)
      pdj(30) =  &
          -k(744)*n(idx_N)
      pdj(31) =  &
          -k(807)*n(idx_NH)
      pdj(33) =  &
          -k(850)*n(idx_SIH)
      pdj(34) =  &
          -k(851)*n(idx_SIO)
      pdj(38) =  &
          -k(838)*n(idx_CO2)
      pdj(55) =  &
          +k(1275)
      pdj(70) =  &
          +k(839)*n(idx_CO)  &
          +k(222)*n(idx_HCO)
      pdj(74) =  &
          +k(46)*n(idx_CH2)  &
          +k(476)*n(idx_CH)
      pdj(75) =  &
          +k(394)*n(idx_C)  &
          +k(63)*n(idx_CH)
      pdj(76) =  &
          +k(220)*n(idx_H2CO)  &
          +k(844)*n(idx_HCO)
      pdj(78) =  &
          +k(792)*n(idx_NH2)  &
          +k(223)*n(idx_NH3)
      pdj(79) =  &
          +k(224)*n(idx_NO)  &
          +k(744)*n(idx_N)
      pdj(84) =  &
          +k(850)*n(idx_SIH)
      pdj(89) =  &
          +k(225)*n(idx_O2)  &
          +k(831)*n(idx_O)
      pdj(90) =  &
          +k(546)*n(idx_H2)  &
          +k(843)*n(idx_HCO)  &
          +k(221)*n(idx_H2O)  &
          +k(848)*n(idx_OH)
      pdj(91) =  &
          +k(189)*n(idx_NH2)  &
          +k(807)*n(idx_NH)
      pdj(92) =  &
          +k(1174)
      pdj(93) =  &
          -k(189)*n(idx_NH2)  &
          -k(851)*n(idx_SIO)  &
          -k(546)*n(idx_H2)  &
          -k(349)*n(idx_E)  &
          -k(792)*n(idx_NH2)  &
          -k(844)*n(idx_HCO)  &
          -k(476)*n(idx_CH)  &
          -k(843)*n(idx_HCO)  &
          -k(839)*n(idx_CO)  &
          -k(222)*n(idx_HCO)  &
          -k(221)*n(idx_H2O)  &
          -k(831)*n(idx_O)  &
          -k(845)*n(idx_HNC)  &
          -k(223)*n(idx_NH3)  &
          -k(842)*n(idx_HCN)  &
          -k(63)*n(idx_CH)  &
          -k(840)*n(idx_H2CO)  &
          -k(458)*n(idx_CH4)  &
          -k(807)*n(idx_NH)  &
          -k(850)*n(idx_SIH)  &
          -k(438)*n(idx_CH2)  &
          -k(394)*n(idx_C)  &
          -k(841)*n(idx_H2O)  &
          -k(846)*n(idx_N2)  &
          -k(744)*n(idx_N)  &
          -k(848)*n(idx_OH)  &
          -k(849)*n(idx_SI)  &
          -k(220)*n(idx_H2CO)  &
          -k(46)*n(idx_CH2)  &
          -k(837)*n(idx_CN)  &
          -k(1174)  &
          -k(224)*n(idx_NO)  &
          -k(847)*n(idx_NO)  &
          -k(838)*n(idx_CO2)  &
          -k(1275)  &
          -k(225)*n(idx_O2)
      pdj(94) =  &
          +k(438)*n(idx_CH2)
      pdj(97) =  &
          +k(837)*n(idx_CN)
      pdj(100) =  &
          +k(849)*n(idx_SI)
      pdj(104) =  &
          +k(847)*n(idx_NO)
      pdj(107) =  &
          +k(840)*n(idx_H2CO)
      pdj(108) =  &
          +k(841)*n(idx_H2O)  &
          +k(458)*n(idx_CH4)
      pdj(109) =  &
          +k(845)*n(idx_HNC)  &
          +k(842)*n(idx_HCN)
      pdj(110) =  &
          +k(838)*n(idx_CO2)
      pdj(112) =  &
          +k(846)*n(idx_N2)
      pdj(115) =  &
          +k(851)*n(idx_SIO)
    elseif(j==94) then
      pdj(1) =  &
          -k(299)*n(idx_E)  &
          -k(300)*n(idx_E)  &
          -k(301)*n(idx_E)  &
          -k(1216)*n(idx_E)
      pdj(2) =  &
          +k(300)*n(idx_E)  &
          +k(301)*n(idx_E)
      pdj(3) =  &
          -k(445)*n(idx_O)  &
          +k(443)*n(idx_O2)  &
          -k(444)*n(idx_O)
      pdj(6) =  &
          +k(445)*n(idx_O)  &
          +k(617)*n(idx_H)  &
          +k(1115)  &
          +k(795)*n(idx_NH)  &
          +k(300)*n(idx_E)  &
          +k(446)*n(idx_OH)
      pdj(8) =  &
          +2.d0*k(301)*n(idx_E)  &
          -k(617)*n(idx_H)  &
          +k(444)*n(idx_O)  &
          +k(1116)  &
          +k(299)*n(idx_E)
      pdj(10) =  &
          -k(446)*n(idx_OH)
      pdj(11) =  &
          -k(443)*n(idx_O2)
      pdj(12) =  &
          +k(299)*n(idx_E)
      pdj(13) =  &
          -k(441)*n(idx_H2CO)
      pdj(14) =  &
          -k(442)*n(idx_HCO)  &
          -k(47)*n(idx_HCO)
      pdj(15) =  &
          -k(48)*n(idx_MG)
      pdj(17) =  &
          -k(49)*n(idx_NO)
      pdj(25) =  &
          +k(442)*n(idx_HCO)
      pdj(28) =  &
          +k(47)*n(idx_HCO)  &
          +k(48)*n(idx_MG)  &
          +k(1216)*n(idx_E)  &
          +k(49)*n(idx_NO)
      pdj(29) =  &
          +k(447)*n(idx_SIH4)  &
          +k(441)*n(idx_H2CO)  &
          +k(440)*n(idx_CH3OH)
      pdj(31) =  &
          -k(795)*n(idx_NH)
      pdj(32) =  &
          -k(447)*n(idx_SIH4)
      pdj(37) =  &
          -k(440)*n(idx_CH3OH)
      pdj(53) =  &
          +k(1286)
      pdj(70) =  &
          +k(47)*n(idx_HCO)  &
          +k(445)*n(idx_O)  &
          +k(441)*n(idx_H2CO)
      pdj(74) =  &
          +k(1116)  &
          +k(617)*n(idx_H)
      pdj(75) =  &
          +k(1115)
      pdj(76) =  &
          +k(444)*n(idx_O)  &
          +k(446)*n(idx_OH)
      pdj(77) =  &
          +k(48)*n(idx_MG)
      pdj(79) =  &
          +k(49)*n(idx_NO)
      pdj(85) =  &
          +k(447)*n(idx_SIH4)
      pdj(94) =  &
          -k(1216)*n(idx_E)  &
          -k(443)*n(idx_O2)  &
          -k(795)*n(idx_NH)  &
          -k(49)*n(idx_NO)  &
          -k(444)*n(idx_O)  &
          -k(446)*n(idx_OH)  &
          -k(299)*n(idx_E)  &
          -k(300)*n(idx_E)  &
          -k(617)*n(idx_H)  &
          -k(1115)  &
          -k(447)*n(idx_SIH4)  &
          -k(301)*n(idx_E)  &
          -k(48)*n(idx_MG)  &
          -k(1116)  &
          -k(440)*n(idx_CH3OH)  &
          -k(442)*n(idx_HCO)  &
          -k(441)*n(idx_H2CO)  &
          -k(445)*n(idx_O)  &
          -k(1286)  &
          -k(47)*n(idx_HCO)
      pdj(95) =  &
          +k(442)*n(idx_HCO)
      pdj(107) =  &
          +k(443)*n(idx_O2)  &
          +k(440)*n(idx_CH3OH)
      pdj(109) =  &
          +k(795)*n(idx_NH)
    elseif(j==95) then
      pdj(1) =  &
          -k(302)*n(idx_E)  &
          -k(303)*n(idx_E)
      pdj(3) =  &
          -k(823)*n(idx_O)
      pdj(6) =  &
          +k(618)*n(idx_H)  &
          +k(1123)
      pdj(8) =  &
          -k(618)*n(idx_H)  &
          +k(1124)  &
          +2.d0*k(302)*n(idx_E)  &
          +k(303)*n(idx_E)
      pdj(9) =  &
          -k(451)*n(idx_H2O)
      pdj(10) =  &
          +k(823)*n(idx_O)
      pdj(11) =  &
          -k(52)*n(idx_O2)
      pdj(12) =  &
          +k(302)*n(idx_E)
      pdj(13) =  &
          -k(50)*n(idx_H2CO)  &
          -k(450)*n(idx_H2CO)
      pdj(16) =  &
          -k(51)*n(idx_NH3)
      pdj(25) =  &
          -k(449)*n(idx_CO)
      pdj(28) =  &
          +k(448)*n(idx_CO2)  &
          +k(451)*n(idx_H2O)  &
          +k(450)*n(idx_H2CO)  &
          +k(449)*n(idx_CO)  &
          +k(303)*n(idx_E)
      pdj(29) =  &
          +k(52)*n(idx_O2)  &
          +k(51)*n(idx_NH3)  &
          +k(50)*n(idx_H2CO)
      pdj(38) =  &
          -k(448)*n(idx_CO2)
      pdj(53) =  &
          +k(1289)
      pdj(70) =  &
          +k(449)*n(idx_CO)
      pdj(74) =  &
          +k(1123)
      pdj(76) =  &
          +k(50)*n(idx_H2CO)
      pdj(78) =  &
          +k(51)*n(idx_NH3)
      pdj(89) =  &
          +k(52)*n(idx_O2)
      pdj(94) =  &
          +k(823)*n(idx_O)  &
          +k(618)*n(idx_H)  &
          +k(1124)
      pdj(95) =  &
          -k(448)*n(idx_CO2)  &
          -k(823)*n(idx_O)  &
          -k(449)*n(idx_CO)  &
          -k(618)*n(idx_H)  &
          -k(50)*n(idx_H2CO)  &
          -k(302)*n(idx_E)  &
          -k(303)*n(idx_E)  &
          -k(51)*n(idx_NH3)  &
          -k(1289)  &
          -k(1123)  &
          -k(451)*n(idx_H2O)  &
          -k(52)*n(idx_O2)  &
          -k(1124)  &
          -k(450)*n(idx_H2CO)
      pdj(107) =  &
          +k(450)*n(idx_H2CO)
      pdj(108) =  &
          +k(451)*n(idx_H2O)
      pdj(110) =  &
          +k(448)*n(idx_CO2)
    elseif(j==96) then
      pdj(1) =  &
          -k(1221)*n(idx_E)
      pdj(2) =  &
          -k(58)*n(idx_CH)  &
          -k(469)*n(idx_CH)
      pdj(3) =  &
          +k(729)*n(idx_O2)  &
          +k(728)*n(idx_NO)
      pdj(5) =  &
          -k(162)*n(idx_HCN)
      pdj(6) =  &
          -k(539)*n(idx_H2)  &
          +k(718)*n(idx_CH4)  &
          +k(725)*n(idx_NH3)
      pdj(7) =  &
          +k(721)*n(idx_CO)
      pdj(8) =  &
          +2.d0*k(719)*n(idx_CH4)  &
          +k(539)*n(idx_H2)  &
          +k(716)*n(idx_CH3OH)  &
          +k(727)*n(idx_NH)  &
          +k(469)*n(idx_CH)  &
          +k(717)*n(idx_CH4)  &
          +k(713)*n(idx_CH3OH)  &
          +k(718)*n(idx_CH4)  &
          +k(715)*n(idx_CH3OH)
      pdj(9) =  &
          -k(161)*n(idx_H2O)
      pdj(10) =  &
          -k(170)*n(idx_OH)
      pdj(11) =  &
          -k(730)*n(idx_O2)  &
          -k(729)*n(idx_O2)  &
          -k(169)*n(idx_O2)
      pdj(12) =  &
          +k(723)*n(idx_H2CO)  &
          -k(156)*n(idx_CH2)
      pdj(13) =  &
          -k(160)*n(idx_H2CO)  &
          -k(722)*n(idx_H2CO)  &
          -k(723)*n(idx_H2CO)
      pdj(14) =  &
          -k(163)*n(idx_HCO)  &
          -k(724)*n(idx_HCO)
      pdj(15) =  &
          -k(164)*n(idx_MG)
      pdj(16) =  &
          -k(726)*n(idx_NH3)  &
          -k(166)*n(idx_NH3)  &
          -k(725)*n(idx_NH3)
      pdj(17) =  &
          +k(730)*n(idx_O2)  &
          +k(716)*n(idx_CH3OH)  &
          +k(720)*n(idx_CO2)  &
          -k(728)*n(idx_NO)  &
          -k(168)*n(idx_NO)
      pdj(24) =  &
          -k(158)*n(idx_CN)
      pdj(25) =  &
          +k(724)*n(idx_HCO)  &
          -k(159)*n(idx_CO)  &
          -k(721)*n(idx_CO)
      pdj(27) =  &
          -k(165)*n(idx_NH2)
      pdj(28) =  &
          +k(715)*n(idx_CH3OH)
      pdj(29) =  &
          -k(157)*n(idx_CH4)  &
          -k(719)*n(idx_CH4)  &
          -k(717)*n(idx_CH4)  &
          -k(718)*n(idx_CH4)
      pdj(30) =  &
          +k(165)*n(idx_NH2)  &
          +k(167)*n(idx_NH)  &
          +k(170)*n(idx_OH)  &
          +k(160)*n(idx_H2CO)  &
          +k(157)*n(idx_CH4)  &
          +k(168)*n(idx_NO)  &
          +k(58)*n(idx_CH)  &
          +k(156)*n(idx_CH2)  &
          +k(162)*n(idx_HCN)  &
          +k(159)*n(idx_CO)  &
          +k(1221)*n(idx_E)  &
          +k(161)*n(idx_H2O)  &
          +k(169)*n(idx_O2)  &
          +k(158)*n(idx_CN)  &
          +k(164)*n(idx_MG)  &
          +k(717)*n(idx_CH4)  &
          +k(163)*n(idx_HCO)  &
          +k(166)*n(idx_NH3)  &
          -k(1211)*n(idx_N)
      pdj(31) =  &
          +k(714)*n(idx_CH3OH)  &
          -k(727)*n(idx_NH)  &
          +k(722)*n(idx_H2CO)  &
          +k(726)*n(idx_NH3)  &
          +k(713)*n(idx_CH3OH)  &
          -k(167)*n(idx_NH)
      pdj(37) =  &
          -k(714)*n(idx_CH3OH)  &
          -k(713)*n(idx_CH3OH)  &
          -k(716)*n(idx_CH3OH)  &
          -k(715)*n(idx_CH3OH)
      pdj(38) =  &
          -k(720)*n(idx_CO2)
      pdj(60) =  &
          +k(1269)
      pdj(70) =  &
          +k(163)*n(idx_HCO)  &
          +k(722)*n(idx_H2CO)
      pdj(74) =  &
          +k(156)*n(idx_CH2)
      pdj(75) =  &
          +k(58)*n(idx_CH)
      pdj(76) =  &
          +k(713)*n(idx_CH3OH)  &
          +k(160)*n(idx_H2CO)
      pdj(77) =  &
          +k(164)*n(idx_MG)
      pdj(78) =  &
          +k(166)*n(idx_NH3)
      pdj(79) =  &
          +k(168)*n(idx_NO)  &
          +k(729)*n(idx_O2)  &
          +k(721)*n(idx_CO)  &
          +k(715)*n(idx_CH3OH)  &
          +k(723)*n(idx_H2CO)
      pdj(86) =  &
          +k(158)*n(idx_CN)  &
          +k(469)*n(idx_CH)
      pdj(87) =  &
          +k(159)*n(idx_CO)  &
          +k(720)*n(idx_CO2)
      pdj(88) =  &
          +k(1211)*n(idx_N)  &
          +k(727)*n(idx_NH)  &
          +k(728)*n(idx_NO)
      pdj(89) =  &
          +k(169)*n(idx_O2)
      pdj(90) =  &
          +k(161)*n(idx_H2O)
      pdj(91) =  &
          +k(165)*n(idx_NH2)  &
          +k(726)*n(idx_NH3)
      pdj(92) =  &
          +k(730)*n(idx_O2)
      pdj(93) =  &
          +k(170)*n(idx_OH)
      pdj(94) =  &
          +k(716)*n(idx_CH3OH)  &
          +k(717)*n(idx_CH4)
      pdj(95) =  &
          +k(157)*n(idx_CH4)
      pdj(96) =  &
          -k(163)*n(idx_HCO)  &
          -k(719)*n(idx_CH4)  &
          -k(720)*n(idx_CO2)  &
          -k(165)*n(idx_NH2)  &
          -k(168)*n(idx_NO)  &
          -k(730)*n(idx_O2)  &
          -k(714)*n(idx_CH3OH)  &
          -k(159)*n(idx_CO)  &
          -k(1269)  &
          -k(160)*n(idx_H2CO)  &
          -k(727)*n(idx_NH)  &
          -k(1221)*n(idx_E)  &
          -k(169)*n(idx_O2)  &
          -k(1211)*n(idx_N)  &
          -k(156)*n(idx_CH2)  &
          -k(729)*n(idx_O2)  &
          -k(170)*n(idx_OH)  &
          -k(166)*n(idx_NH3)  &
          -k(726)*n(idx_NH3)  &
          -k(58)*n(idx_CH)  &
          -k(717)*n(idx_CH4)  &
          -k(723)*n(idx_H2CO)  &
          -k(158)*n(idx_CN)  &
          -k(728)*n(idx_NO)  &
          -k(722)*n(idx_H2CO)  &
          -k(167)*n(idx_NH)  &
          -k(157)*n(idx_CH4)  &
          -k(721)*n(idx_CO)  &
          -k(162)*n(idx_HCN)  &
          -k(469)*n(idx_CH)  &
          -k(718)*n(idx_CH4)  &
          -k(725)*n(idx_NH3)  &
          -k(713)*n(idx_CH3OH)  &
          -k(716)*n(idx_CH3OH)  &
          -k(715)*n(idx_CH3OH)  &
          -k(539)*n(idx_H2)  &
          -k(724)*n(idx_HCO)  &
          -k(164)*n(idx_MG)  &
          -k(161)*n(idx_H2O)
      pdj(97) =  &
          +k(162)*n(idx_HCN)  &
          +k(718)*n(idx_CH4)
      pdj(98) =  &
          +k(167)*n(idx_NH)  &
          +k(539)*n(idx_H2)  &
          +k(724)*n(idx_HCO)
      pdj(107) =  &
          +k(714)*n(idx_CH3OH)
      pdj(109) =  &
          +k(719)*n(idx_CH4)
      pdj(112) =  &
          +k(725)*n(idx_NH3)
    elseif(j==97) then
      pdj(1) =  &
          -k(327)*n(idx_E)
      pdj(2) =  &
          -k(464)*n(idx_CH)
      pdj(4) =  &
          -k(627)*n(idx_HNC)
      pdj(5) =  &
          +k(197)*n(idx_NH3)  &
          +k(125)*n(idx_H2O)  &
          +k(130)*n(idx_H)  &
          +k(134)*n(idx_O2)  &
          +k(133)*n(idx_NO)  &
          -k(624)*n(idx_HCN)
      pdj(6) =  &
          -k(536)*n(idx_H2)
      pdj(7) =  &
          -k(386)*n(idx_C)
      pdj(8) =  &
          -k(130)*n(idx_H)  &
          +k(536)*n(idx_H2)  &
          +k(327)*n(idx_E)
      pdj(9) =  &
          -k(125)*n(idx_H2O)  &
          -k(566)*n(idx_H2O)
      pdj(10) =  &
          -k(854)*n(idx_OH)
      pdj(11) =  &
          -k(134)*n(idx_O2)
      pdj(12) =  &
          -k(427)*n(idx_CH2)
      pdj(13) =  &
          -k(623)*n(idx_H2CO)
      pdj(14) =  &
          -k(626)*n(idx_HCO)  &
          -k(625)*n(idx_HCO)
      pdj(16) =  &
          -k(197)*n(idx_NH3)  &
          -k(794)*n(idx_NH3)
      pdj(17) =  &
          -k(133)*n(idx_NO)
      pdj(24) =  &
          +k(624)*n(idx_HCN)  &
          +k(386)*n(idx_C)  &
          +k(464)*n(idx_CH)  &
          +k(623)*n(idx_H2CO)  &
          +k(566)*n(idx_H2O)  &
          +k(427)*n(idx_CH2)  &
          +k(799)*n(idx_NH)  &
          +k(622)*n(idx_CO)  &
          +k(785)*n(idx_NH2)  &
          +k(627)*n(idx_HNC)  &
          +k(327)*n(idx_E)  &
          +k(621)*n(idx_CO2)  &
          +k(854)*n(idx_OH)  &
          +k(625)*n(idx_HCO)
      pdj(25) =  &
          -k(622)*n(idx_CO)  &
          +k(626)*n(idx_HCO)
      pdj(27) =  &
          +k(794)*n(idx_NH3)  &
          -k(785)*n(idx_NH2)
      pdj(28) =  &
          +k(455)*n(idx_CH4)
      pdj(29) =  &
          -k(455)*n(idx_CH4)
      pdj(31) =  &
          -k(799)*n(idx_NH)
      pdj(38) =  &
          -k(621)*n(idx_CO2)
      pdj(59) =  &
          +k(1283)
      pdj(70) =  &
          +k(622)*n(idx_CO)
      pdj(71) =  &
          +k(130)*n(idx_H)
      pdj(74) =  &
          +k(464)*n(idx_CH)
      pdj(75) =  &
          +k(386)*n(idx_C)
      pdj(76) =  &
          +k(625)*n(idx_HCO)
      pdj(78) =  &
          +k(197)*n(idx_NH3)  &
          +k(785)*n(idx_NH2)
      pdj(79) =  &
          +k(133)*n(idx_NO)
      pdj(89) =  &
          +k(134)*n(idx_O2)
      pdj(90) =  &
          +k(854)*n(idx_OH)  &
          +k(125)*n(idx_H2O)
      pdj(91) =  &
          +k(799)*n(idx_NH)
      pdj(94) =  &
          +k(427)*n(idx_CH2)
      pdj(97) =  &
          -k(125)*n(idx_H2O)  &
          -k(197)*n(idx_NH3)  &
          -k(427)*n(idx_CH2)  &
          -k(794)*n(idx_NH3)  &
          -k(455)*n(idx_CH4)  &
          -k(566)*n(idx_H2O)  &
          -k(464)*n(idx_CH)  &
          -k(130)*n(idx_H)  &
          -k(386)*n(idx_C)  &
          -k(621)*n(idx_CO2)  &
          -k(625)*n(idx_HCO)  &
          -k(536)*n(idx_H2)  &
          -k(626)*n(idx_HCO)  &
          -k(627)*n(idx_HNC)  &
          -k(623)*n(idx_H2CO)  &
          -k(799)*n(idx_NH)  &
          -k(1283)  &
          -k(622)*n(idx_CO)  &
          -k(134)*n(idx_O2)  &
          -k(785)*n(idx_NH2)  &
          -k(854)*n(idx_OH)  &
          -k(133)*n(idx_NO)  &
          -k(327)*n(idx_E)  &
          -k(624)*n(idx_HCN)
      pdj(107) =  &
          +k(623)*n(idx_H2CO)
      pdj(108) =  &
          +k(566)*n(idx_H2O)
      pdj(109) =  &
          +k(794)*n(idx_NH3)  &
          +k(536)*n(idx_H2)  &
          +k(627)*n(idx_HNC)  &
          +k(626)*n(idx_HCO)  &
          +k(624)*n(idx_HCN)  &
          +k(455)*n(idx_CH4)
      pdj(110) =  &
          +k(621)*n(idx_CO2)
    elseif(j==98) then
      pdj(1) =  &
          -k(341)*n(idx_E)
      pdj(2) =  &
          -k(471)*n(idx_CH)
      pdj(3) =  &
          +k(765)*n(idx_NO)  &
          -k(768)*n(idx_O)  &
          +k(757)*n(idx_H2O)
      pdj(4) =  &
          -k(761)*n(idx_HNC)
      pdj(5) =  &
          -k(759)*n(idx_HCN)
      pdj(6) =  &
          +k(756)*n(idx_H2O)  &
          -k(542)*n(idx_H2)  &
          -k(541)*n(idx_H2)
      pdj(7) =  &
          -k(391)*n(idx_C)
      pdj(8) =  &
          +k(542)*n(idx_H2)  &
          +k(741)*n(idx_N)  &
          +k(341)*n(idx_E)
      pdj(9) =  &
          -k(757)*n(idx_H2O)  &
          -k(177)*n(idx_H2O)  &
          -k(756)*n(idx_H2O)  &
          -k(755)*n(idx_H2O)  &
          -k(758)*n(idx_H2O)
      pdj(10) =  &
          +k(758)*n(idx_H2O)  &
          -k(769)*n(idx_OH)  &
          +k(766)*n(idx_O2)
      pdj(11) =  &
          -k(766)*n(idx_O2)  &
          -k(767)*n(idx_O2)  &
          -k(180)*n(idx_O2)
      pdj(12) =  &
          -k(433)*n(idx_CH2)
      pdj(13) =  &
          -k(176)*n(idx_H2CO)  &
          -k(754)*n(idx_H2CO)  &
          -k(753)*n(idx_H2CO)
      pdj(14) =  &
          -k(760)*n(idx_HCO)  &
          +k(751)*n(idx_CO2)
      pdj(16) =  &
          -k(178)*n(idx_NH3)
      pdj(17) =  &
          -k(765)*n(idx_NO)  &
          -k(179)*n(idx_NO)
      pdj(24) =  &
          -k(748)*n(idx_CN)
      pdj(25) =  &
          -k(752)*n(idx_CO)  &
          +k(750)*n(idx_CO2)
      pdj(26) =  &
          -k(762)*n(idx_N2)
      pdj(27) =  &
          +k(754)*n(idx_H2CO)  &
          -k(763)*n(idx_NH2)
      pdj(30) =  &
          +k(769)*n(idx_OH)  &
          +k(748)*n(idx_CN)  &
          +k(341)*n(idx_E)  &
          +k(1157)  &
          +k(433)*n(idx_CH2)  &
          +k(471)*n(idx_CH)  &
          +k(764)*n(idx_NH)  &
          +k(768)*n(idx_O)  &
          +k(541)*n(idx_H2)  &
          -k(741)*n(idx_N)  &
          +k(763)*n(idx_NH2)  &
          +k(760)*n(idx_HCO)  &
          +k(759)*n(idx_HCN)  &
          +k(749)*n(idx_CO2)  &
          +k(761)*n(idx_HNC)  &
          +k(755)*n(idx_H2O)  &
          +k(753)*n(idx_H2CO)  &
          +k(752)*n(idx_CO)  &
          +k(391)*n(idx_C)  &
          +k(762)*n(idx_N2)  &
          +k(767)*n(idx_O2)
      pdj(31) =  &
          -k(764)*n(idx_NH)  &
          +k(176)*n(idx_H2CO)  &
          +k(177)*n(idx_H2O)  &
          +k(180)*n(idx_O2)  &
          +k(178)*n(idx_NH3)  &
          +k(179)*n(idx_NO)
      pdj(38) =  &
          -k(750)*n(idx_CO2)  &
          -k(751)*n(idx_CO2)  &
          -k(749)*n(idx_CO2)
      pdj(60) =  &
          +k(1274)
      pdj(70) =  &
          +k(752)*n(idx_CO)  &
          +k(754)*n(idx_H2CO)
      pdj(71) =  &
          +k(1157)
      pdj(74) =  &
          +k(471)*n(idx_CH)
      pdj(75) =  &
          +k(391)*n(idx_C)
      pdj(76) =  &
          +k(760)*n(idx_HCO)  &
          +k(176)*n(idx_H2CO)
      pdj(78) =  &
          +k(763)*n(idx_NH2)  &
          +k(178)*n(idx_NH3)  &
          +k(757)*n(idx_H2O)
      pdj(79) =  &
          +k(766)*n(idx_O2)  &
          +k(751)*n(idx_CO2)  &
          +k(179)*n(idx_NO)
      pdj(88) =  &
          +k(741)*n(idx_N)
      pdj(89) =  &
          +k(180)*n(idx_O2)
      pdj(90) =  &
          +k(769)*n(idx_OH)  &
          +k(177)*n(idx_H2O)
      pdj(91) =  &
          +k(542)*n(idx_H2)  &
          +k(758)*n(idx_H2O)  &
          +k(764)*n(idx_NH)
      pdj(93) =  &
          +k(768)*n(idx_O)
      pdj(94) =  &
          +k(433)*n(idx_CH2)
      pdj(97) =  &
          +k(748)*n(idx_CN)
      pdj(98) =  &
          -k(757)*n(idx_H2O)  &
          -k(433)*n(idx_CH2)  &
          -k(542)*n(idx_H2)  &
          -k(752)*n(idx_CO)  &
          -k(755)*n(idx_H2O)  &
          -k(769)*n(idx_OH)  &
          -k(766)*n(idx_O2)  &
          -k(759)*n(idx_HCN)  &
          -k(760)*n(idx_HCO)  &
          -k(753)*n(idx_H2CO)  &
          -k(748)*n(idx_CN)  &
          -k(765)*n(idx_NO)  &
          -k(176)*n(idx_H2CO)  &
          -k(391)*n(idx_C)  &
          -k(754)*n(idx_H2CO)  &
          -k(762)*n(idx_N2)  &
          -k(180)*n(idx_O2)  &
          -k(750)*n(idx_CO2)  &
          -k(764)*n(idx_NH)  &
          -k(756)*n(idx_H2O)  &
          -k(178)*n(idx_NH3)  &
          -k(1274)  &
          -k(749)*n(idx_CO2)  &
          -k(758)*n(idx_H2O)  &
          -k(751)*n(idx_CO2)  &
          -k(471)*n(idx_CH)  &
          -k(767)*n(idx_O2)  &
          -k(341)*n(idx_E)  &
          -k(177)*n(idx_H2O)  &
          -k(179)*n(idx_NO)  &
          -k(741)*n(idx_N)  &
          -k(541)*n(idx_H2)  &
          -k(768)*n(idx_O)  &
          -k(761)*n(idx_HNC)  &
          -k(1157)  &
          -k(763)*n(idx_NH2)
      pdj(104) =  &
          +k(756)*n(idx_H2O)  &
          +k(750)*n(idx_CO2)
      pdj(106) =  &
          +k(541)*n(idx_H2)
      pdj(107) =  &
          +k(753)*n(idx_H2CO)
      pdj(108) =  &
          +k(755)*n(idx_H2O)
      pdj(109) =  &
          +k(759)*n(idx_HCN)  &
          +k(761)*n(idx_HNC)
      pdj(110) =  &
          +k(749)*n(idx_CO2)
      pdj(112) =  &
          +k(765)*n(idx_NO)  &
          +k(762)*n(idx_N2)
      pdj(113) =  &
          +k(767)*n(idx_O2)
    elseif(j==99) then
      pdj(1) =  &
          -k(359)*n(idx_E)  &
          -k(360)*n(idx_E)
      pdj(6) =  &
          -k(547)*n(idx_H2)  &
          +k(359)*n(idx_E)
      pdj(8) =  &
          +k(360)*n(idx_E)  &
          +k(547)*n(idx_H2)
      pdj(9) =  &
          -k(575)*n(idx_H2O)
      pdj(22) =  &
          +k(359)*n(idx_E)
      pdj(23) =  &
          +k(490)*n(idx_CO)  &
          +k(575)*n(idx_H2O)  &
          +k(360)*n(idx_E)
      pdj(25) =  &
          -k(490)*n(idx_CO)
      pdj(48) =  &
          +k(1239)
      pdj(70) =  &
          +k(490)*n(idx_CO)
      pdj(99) =  &
          -k(359)*n(idx_E)  &
          -k(547)*n(idx_H2)  &
          -k(1239)  &
          -k(360)*n(idx_E)  &
          -k(575)*n(idx_H2O)  &
          -k(490)*n(idx_CO)
      pdj(108) =  &
          +k(575)*n(idx_H2O)
      pdj(114) =  &
          +k(547)*n(idx_H2)
    elseif(j==100) then
      pdj(1) =  &
          -k(353)*n(idx_E)
      pdj(2) =  &
          -k(478)*n(idx_CH)
      pdj(3) =  &
          -k(833)*n(idx_O)
      pdj(6) =  &
          -k(1204)*n(idx_H2)  &
          +k(620)*n(idx_H)
      pdj(7) =  &
          -k(395)*n(idx_C)
      pdj(8) =  &
          +k(353)*n(idx_E)  &
          +k(395)*n(idx_C)  &
          -k(620)*n(idx_H)  &
          +k(833)*n(idx_O)  &
          +k(1180)
      pdj(9) =  &
          -k(574)*n(idx_H2O)
      pdj(18) =  &
          +k(353)*n(idx_E)  &
          +k(574)*n(idx_H2O)  &
          +k(478)*n(idx_CH)
      pdj(48) =  &
          +k(1233)
      pdj(74) =  &
          +k(478)*n(idx_CH)
      pdj(80) =  &
          +k(1180)  &
          +k(620)*n(idx_H)
      pdj(83) =  &
          +k(395)*n(idx_C)
      pdj(85) =  &
          +k(1204)*n(idx_H2)
      pdj(100) =  &
          -k(395)*n(idx_C)  &
          -k(1204)*n(idx_H2)  &
          -k(574)*n(idx_H2O)  &
          -k(1233)  &
          -k(1180)  &
          -k(478)*n(idx_CH)  &
          -k(620)*n(idx_H)  &
          -k(833)*n(idx_O)  &
          -k(353)*n(idx_E)
      pdj(101) =  &
          +k(833)*n(idx_O)
      pdj(108) =  &
          +k(574)*n(idx_H2O)
    elseif(j==101) then
      pdj(1) =  &
          -k(363)*n(idx_E)
      pdj(2) =  &
          -k(479)*n(idx_CH)
      pdj(3) =  &
          +k(1190)  &
          -k(836)*n(idx_O)  &
          +k(363)*n(idx_E)
      pdj(6) =  &
          -k(548)*n(idx_H2)
      pdj(7) =  &
          -k(396)*n(idx_C)
      pdj(8) =  &
          +k(548)*n(idx_H2)
      pdj(11) =  &
          +k(836)*n(idx_O)
      pdj(12) =  &
          -k(439)*n(idx_CH2)
      pdj(13) =  &
          +k(439)*n(idx_CH2)
      pdj(14) =  &
          -k(139)*n(idx_HCO)
      pdj(15) =  &
          -k(155)*n(idx_MG)
      pdj(17) =  &
          -k(207)*n(idx_NO)  &
          +k(747)*n(idx_N)
      pdj(18) =  &
          +k(479)*n(idx_CH)  &
          +k(363)*n(idx_E)  &
          +k(746)*n(idx_N)
      pdj(25) =  &
          -k(491)*n(idx_CO)  &
          +k(396)*n(idx_C)
      pdj(30) =  &
          -k(747)*n(idx_N)  &
          -k(746)*n(idx_N)
      pdj(34) =  &
          +k(155)*n(idx_MG)  &
          +k(207)*n(idx_NO)  &
          +k(139)*n(idx_HCO)
      pdj(38) =  &
          +k(491)*n(idx_CO)
      pdj(49) =  &
          +k(1246)
      pdj(70) =  &
          +k(479)*n(idx_CH)  &
          +k(139)*n(idx_HCO)
      pdj(77) =  &
          +k(155)*n(idx_MG)
      pdj(79) =  &
          +k(207)*n(idx_NO)  &
          +k(746)*n(idx_N)
      pdj(80) =  &
          +k(747)*n(idx_N)  &
          +k(439)*n(idx_CH2)  &
          +k(836)*n(idx_O)  &
          +k(491)*n(idx_CO)  &
          +k(1190)  &
          +k(396)*n(idx_C)
      pdj(101) =  &
          -k(747)*n(idx_N)  &
          -k(491)*n(idx_CO)  &
          -k(139)*n(idx_HCO)  &
          -k(746)*n(idx_N)  &
          -k(155)*n(idx_MG)  &
          -k(1246)  &
          -k(479)*n(idx_CH)  &
          -k(836)*n(idx_O)  &
          -k(207)*n(idx_NO)  &
          -k(363)*n(idx_E)  &
          -k(439)*n(idx_CH2)  &
          -k(1190)  &
          -k(548)*n(idx_H2)  &
          -k(396)*n(idx_C)
      pdj(115) =  &
          +k(548)*n(idx_H2)
    elseif(j==102) then
      pdj(1) =  &
          -k(306)*n(idx_E)
      pdj(2) =  &
          -k(103)*n(idx_CH)  &
          -k(513)*n(idx_CH)
      pdj(3) =  &
          -k(527)*n(idx_O)
      pdj(5) =  &
          -k(108)*n(idx_HCN)
      pdj(6) =  &
          +k(102)*n(idx_CH4)  &
          -k(517)*n(idx_H2)  &
          +k(103)*n(idx_CH)  &
          +k(112)*n(idx_NH)  &
          +k(114)*n(idx_O2)  &
          +k(111)*n(idx_NH3)  &
          +k(113)*n(idx_NO)  &
          +k(108)*n(idx_HCN)  &
          +k(110)*n(idx_NH2)  &
          +k(101)*n(idx_CH2)  &
          +k(106)*n(idx_H2CO)  &
          +k(115)*n(idx_OH)  &
          +k(512)*n(idx_CH4)  &
          +k(104)*n(idx_CN)  &
          +k(518)*n(idx_H2CO)  &
          +k(109)*n(idx_HCO)  &
          +k(105)*n(idx_CO)  &
          +k(129)*n(idx_H)  &
          +k(107)*n(idx_H2O)
      pdj(7) =  &
          -k(510)*n(idx_C)
      pdj(8) =  &
          -k(129)*n(idx_H)  &
          +k(516)*n(idx_CO)  &
          +k(511)*n(idx_CH2)  &
          +k(528)*n(idx_OH)  &
          +k(524)*n(idx_NH)  &
          +k(527)*n(idx_O)  &
          +k(1135)  &
          +k(510)*n(idx_C)  &
          +k(521)*n(idx_HE)  &
          +k(523)*n(idx_N)  &
          +k(526)*n(idx_O2)  &
          +k(519)*n(idx_H2O)  &
          +k(522)*n(idx_N2)  &
          +k(515)*n(idx_CO2)  &
          +k(517)*n(idx_H2)  &
          +k(512)*n(idx_CH4)  &
          +2.d0*k(306)*n(idx_E)  &
          +k(518)*n(idx_H2CO)  &
          +k(514)*n(idx_CN)  &
          +k(525)*n(idx_NO)  &
          +k(513)*n(idx_CH)
      pdj(9) =  &
          -k(519)*n(idx_H2O)  &
          -k(107)*n(idx_H2O)
      pdj(10) =  &
          -k(115)*n(idx_OH)  &
          -k(528)*n(idx_OH)
      pdj(11) =  &
          -k(114)*n(idx_O2)  &
          -k(526)*n(idx_O2)
      pdj(12) =  &
          -k(101)*n(idx_CH2)  &
          -k(511)*n(idx_CH2)
      pdj(13) =  &
          -k(518)*n(idx_H2CO)  &
          -k(106)*n(idx_H2CO)
      pdj(14) =  &
          -k(520)*n(idx_HCO)  &
          -k(109)*n(idx_HCO)
      pdj(16) =  &
          -k(111)*n(idx_NH3)
      pdj(17) =  &
          -k(525)*n(idx_NO)  &
          -k(113)*n(idx_NO)
      pdj(24) =  &
          -k(104)*n(idx_CN)  &
          -k(514)*n(idx_CN)
      pdj(25) =  &
          +k(520)*n(idx_HCO)  &
          -k(516)*n(idx_CO)  &
          -k(105)*n(idx_CO)
      pdj(26) =  &
          -k(522)*n(idx_N2)
      pdj(27) =  &
          -k(110)*n(idx_NH2)
      pdj(29) =  &
          -k(512)*n(idx_CH4)  &
          -k(102)*n(idx_CH4)
      pdj(30) =  &
          -k(523)*n(idx_N)
      pdj(31) =  &
          -k(524)*n(idx_NH)  &
          -k(112)*n(idx_NH)
      pdj(35) =  &
          -k(521)*n(idx_HE)
      pdj(38) =  &
          -k(515)*n(idx_CO2)
      pdj(70) =  &
          +k(516)*n(idx_CO)  &
          +k(109)*n(idx_HCO)  &
          +k(518)*n(idx_H2CO)
      pdj(71) =  &
          +k(1135)  &
          +k(129)*n(idx_H)
      pdj(74) =  &
          +k(101)*n(idx_CH2)  &
          +k(513)*n(idx_CH)
      pdj(75) =  &
          +k(103)*n(idx_CH)  &
          +k(510)*n(idx_C)
      pdj(76) =  &
          +k(106)*n(idx_H2CO)
      pdj(78) =  &
          +k(111)*n(idx_NH3)
      pdj(79) =  &
          +k(113)*n(idx_NO)
      pdj(86) =  &
          +k(104)*n(idx_CN)
      pdj(87) =  &
          +k(105)*n(idx_CO)
      pdj(89) =  &
          +k(114)*n(idx_O2)
      pdj(90) =  &
          +k(528)*n(idx_OH)  &
          +k(107)*n(idx_H2O)
      pdj(91) =  &
          +k(110)*n(idx_NH2)  &
          +k(524)*n(idx_NH)
      pdj(93) =  &
          +k(527)*n(idx_O)  &
          +k(115)*n(idx_OH)
      pdj(94) =  &
          +k(511)*n(idx_CH2)  &
          +k(512)*n(idx_CH4)
      pdj(95) =  &
          +k(102)*n(idx_CH4)
      pdj(97) =  &
          +k(514)*n(idx_CN)  &
          +k(108)*n(idx_HCN)
      pdj(98) =  &
          +k(112)*n(idx_NH)  &
          +k(523)*n(idx_N)
      pdj(102) =  &
          -k(523)*n(idx_N)  &
          -k(515)*n(idx_CO2)  &
          -k(522)*n(idx_N2)  &
          -k(516)*n(idx_CO)  &
          -k(306)*n(idx_E)  &
          -k(101)*n(idx_CH2)  &
          -k(102)*n(idx_CH4)  &
          -k(112)*n(idx_NH)  &
          -k(113)*n(idx_NO)  &
          -k(528)*n(idx_OH)  &
          -k(521)*n(idx_HE)  &
          -k(518)*n(idx_H2CO)  &
          -k(111)*n(idx_NH3)  &
          -k(510)*n(idx_C)  &
          -k(107)*n(idx_H2O)  &
          -k(526)*n(idx_O2)  &
          -k(527)*n(idx_O)  &
          -k(115)*n(idx_OH)  &
          -k(108)*n(idx_HCN)  &
          -k(520)*n(idx_HCO)  &
          -k(511)*n(idx_CH2)  &
          -k(129)*n(idx_H)  &
          -k(512)*n(idx_CH4)  &
          -k(104)*n(idx_CN)  &
          -k(525)*n(idx_NO)  &
          -k(103)*n(idx_CH)  &
          -k(114)*n(idx_O2)  &
          -k(514)*n(idx_CN)  &
          -k(513)*n(idx_CH)  &
          -k(524)*n(idx_NH)  &
          -k(109)*n(idx_HCO)  &
          -k(105)*n(idx_CO)  &
          -k(517)*n(idx_H2)  &
          -k(106)*n(idx_H2CO)  &
          -k(110)*n(idx_NH2)  &
          -k(1135)  &
          -k(519)*n(idx_H2O)
      pdj(104) =  &
          +k(525)*n(idx_NO)
      pdj(106) =  &
          +k(520)*n(idx_HCO)  &
          +k(517)*n(idx_H2)
      pdj(108) =  &
          +k(519)*n(idx_H2O)
      pdj(110) =  &
          +k(515)*n(idx_CO2)
      pdj(111) =  &
          +k(521)*n(idx_HE)
      pdj(112) =  &
          +k(522)*n(idx_N2)
      pdj(113) =  &
          +k(526)*n(idx_O2)
    elseif(j==103) then
      pdj(1) =  &
          -k(1219)*n(idx_E)
      pdj(2) =  &
          -k(663)*n(idx_CH)  &
          -k(142)*n(idx_CH)  &
          +k(678)*n(idx_HCN)
      pdj(3) =  &
          +k(666)*n(idx_CO2)  &
          +k(711)*n(idx_SIO)  &
          +k(673)*n(idx_H2CO)  &
          +k(696)*n(idx_NO)  &
          +k(670)*n(idx_CO)  &
          +k(698)*n(idx_OCN)  &
          +k(697)*n(idx_O2)  &
          +k(683)*n(idx_HCO)
      pdj(4) =  &
          -k(684)*n(idx_HNC)  &
          -k(685)*n(idx_HNC)  &
          -k(686)*n(idx_HNC)
      pdj(5) =  &
          -k(677)*n(idx_HCN)  &
          -k(680)*n(idx_HCN)  &
          -k(678)*n(idx_HCN)  &
          -k(679)*n(idx_HCN)
      pdj(6) =  &
          +2.d0*k(708)*n(idx_SIH4)  &
          -k(116)*n(idx_H2)  &
          +k(654)*n(idx_CH2)  &
          +k(659)*n(idx_CH4)  &
          +k(656)*n(idx_CH3)  &
          +k(690)*n(idx_NH2)  &
          +k(706)*n(idx_SIH3)  &
          -k(537)*n(idx_H2)  &
          +k(671)*n(idx_H2CO)  &
          +k(692)*n(idx_NH3)  &
          +k(704)*n(idx_SIH2)  &
          +k(660)*n(idx_CH4)  &
          +k(709)*n(idx_SIH4)
      pdj(7) =  &
          +k(686)*n(idx_HNC)  &
          +k(668)*n(idx_CO2)  &
          +k(701)*n(idx_SIC3)  &
          +k(702)*n(idx_SIC)  &
          +k(664)*n(idx_CN)  &
          -k(140)*n(idx_C)
      pdj(8) =  &
          +k(700)*n(idx_OH)  &
          +k(677)*n(idx_HCN)  &
          +k(672)*n(idx_H2CO)  &
          +k(687)*n(idx_HNO)  &
          +k(679)*n(idx_HCN)  &
          +k(691)*n(idx_NH2)  &
          +k(537)*n(idx_H2)  &
          +k(659)*n(idx_CH4)  &
          +k(685)*n(idx_HNC)  &
          +k(710)*n(idx_SIH)  &
          +k(707)*n(idx_SIH3)  &
          +k(674)*n(idx_H2O)  &
          +k(681)*n(idx_HCO)  &
          +k(676)*n(idx_H2SIO)  &
          +k(693)*n(idx_NH3)  &
          +k(694)*n(idx_NH)  &
          +k(705)*n(idx_SIH2)  &
          +k(661)*n(idx_CH4)  &
          +k(684)*n(idx_HNC)  &
          +k(663)*n(idx_CH)  &
          +k(709)*n(idx_SIH4)  &
          +k(655)*n(idx_CH2)  &
          -k(131)*n(idx_H)
      pdj(9) =  &
          -k(675)*n(idx_H2O)  &
          -k(144)*n(idx_H2O)  &
          -k(674)*n(idx_H2O)
      pdj(10) =  &
          +k(658)*n(idx_CH3OH)  &
          -k(700)*n(idx_OH)  &
          +k(675)*n(idx_H2O)
      pdj(11) =  &
          -k(147)*n(idx_O2)  &
          -k(697)*n(idx_O2)  &
          +k(669)*n(idx_CO2)
      pdj(12) =  &
          -k(654)*n(idx_CH2)  &
          -k(655)*n(idx_CH2)
      pdj(13) =  &
          -k(672)*n(idx_H2CO)  &
          -k(673)*n(idx_H2CO)  &
          -k(143)*n(idx_H2CO)  &
          -k(671)*n(idx_H2CO)
      pdj(14) =  &
          -k(681)*n(idx_HCO)  &
          -k(682)*n(idx_HCO)  &
          -k(683)*n(idx_HCO)
      pdj(16) =  &
          -k(146)*n(idx_NH3)  &
          -k(692)*n(idx_NH3)  &
          -k(693)*n(idx_NH3)
      pdj(17) =  &
          -k(695)*n(idx_NO)  &
          +k(688)*n(idx_HNO)  &
          -k(696)*n(idx_NO)
      pdj(18) =  &
          +k(712)*n(idx_SIO)  &
          -k(148)*n(idx_SI)  &
          +k(703)*n(idx_SIC)
      pdj(20) =  &
          -k(701)*n(idx_SIC3)
      pdj(21) =  &
          -k(702)*n(idx_SIC)  &
          -k(703)*n(idx_SIC)
      pdj(22) =  &
          -k(704)*n(idx_SIH2)  &
          -k(705)*n(idx_SIH2)
      pdj(23) =  &
          -k(706)*n(idx_SIH3)  &
          -k(707)*n(idx_SIH3)
      pdj(24) =  &
          -k(664)*n(idx_CN)  &
          +k(699)*n(idx_OCN)  &
          -k(665)*n(idx_CN)
      pdj(25) =  &
          +k(682)*n(idx_HCO)  &
          -k(670)*n(idx_CO)  &
          +k(667)*n(idx_CO2)
      pdj(26) =  &
          -k(689)*n(idx_N2)  &
          -k(145)*n(idx_N2)
      pdj(27) =  &
          -k(690)*n(idx_NH2)  &
          -k(691)*n(idx_NH2)
      pdj(28) =  &
          +k(662)*n(idx_CH4)  &
          +k(657)*n(idx_CH3OH)  &
          -k(656)*n(idx_CH3)
      pdj(29) =  &
          -k(661)*n(idx_CH4)  &
          -k(660)*n(idx_CH4)  &
          -k(662)*n(idx_CH4)  &
          -k(141)*n(idx_CH4)  &
          -k(659)*n(idx_CH4)
      pdj(30) =  &
          +k(689)*n(idx_N2)  &
          +k(695)*n(idx_NO)  &
          +k(679)*n(idx_HCN)  &
          +k(680)*n(idx_HCN)  &
          +k(665)*n(idx_CN)  &
          +k(685)*n(idx_HNC)
      pdj(31) =  &
          -k(694)*n(idx_NH)
      pdj(32) =  &
          -k(709)*n(idx_SIH4)  &
          -k(708)*n(idx_SIH4)
      pdj(33) =  &
          -k(710)*n(idx_SIH)
      pdj(34) =  &
          -k(712)*n(idx_SIO)  &
          -k(711)*n(idx_SIO)
      pdj(35) =  &
          +k(710)*n(idx_SIH)  &
          +k(688)*n(idx_HNO)  &
          +k(706)*n(idx_SIH3)  &
          +k(686)*n(idx_HNC)  &
          +k(711)*n(idx_SIO)  &
          +k(694)*n(idx_NH)  &
          +k(662)*n(idx_CH4)  &
          +k(678)*n(idx_HCN)  &
          +k(131)*n(idx_H)  &
          +k(692)*n(idx_NH3)  &
          +k(680)*n(idx_HCN)  &
          +k(146)*n(idx_NH3)  &
          +k(677)*n(idx_HCN)  &
          +k(147)*n(idx_O2)  &
          +k(657)*n(idx_CH3OH)  &
          +k(672)*n(idx_H2CO)  &
          +k(687)*n(idx_HNO)  &
          +k(679)*n(idx_HCN)  &
          +k(691)*n(idx_NH2)  &
          +k(700)*n(idx_OH)  &
          +k(709)*n(idx_SIH4)  &
          +k(689)*n(idx_N2)  &
          +k(537)*n(idx_H2)  &
          +k(697)*n(idx_O2)  &
          +k(659)*n(idx_CH4)  &
          +k(656)*n(idx_CH3)  &
          +k(664)*n(idx_CN)  &
          +k(684)*n(idx_HNC)  &
          +k(683)*n(idx_HCO)  &
          +k(704)*n(idx_SIH2)  &
          +k(116)*n(idx_H2)  &
          +k(707)*n(idx_SIH3)  &
          +k(658)*n(idx_CH3OH)  &
          +k(698)*n(idx_OCN)  &
          +k(696)*n(idx_NO)  &
          +k(144)*n(idx_H2O)  &
          +k(142)*n(idx_CH)  &
          +k(674)*n(idx_H2O)  &
          +k(681)*n(idx_HCO)  &
          +k(695)*n(idx_NO)  &
          +k(676)*n(idx_H2SIO)  &
          +k(669)*n(idx_CO2)  &
          +k(701)*n(idx_SIC3)  &
          +k(660)*n(idx_CH4)  &
          +k(673)*n(idx_H2CO)  &
          +k(690)*n(idx_NH2)  &
          +k(143)*n(idx_H2CO)  &
          +k(666)*n(idx_CO2)  &
          +k(693)*n(idx_NH3)  &
          +k(670)*n(idx_CO)  &
          +k(654)*n(idx_CH2)  &
          +k(702)*n(idx_SIC)  &
          +k(148)*n(idx_SI)  &
          +k(665)*n(idx_CN)  &
          +k(685)*n(idx_HNC)  &
          +k(708)*n(idx_SIH4)  &
          +k(705)*n(idx_SIH2)  &
          +k(661)*n(idx_CH4)  &
          +k(699)*n(idx_OCN)  &
          +k(140)*n(idx_C)  &
          +k(668)*n(idx_CO2)  &
          +k(675)*n(idx_H2O)  &
          +k(703)*n(idx_SIC)  &
          +k(141)*n(idx_CH4)  &
          +k(712)*n(idx_SIO)  &
          +k(671)*n(idx_H2CO)  &
          +k(667)*n(idx_CO2)  &
          +k(1219)*n(idx_E)  &
          +k(663)*n(idx_CH)  &
          +k(145)*n(idx_N2)  &
          +k(655)*n(idx_CH2)
      pdj(36) =  &
          -k(688)*n(idx_HNO)  &
          -k(687)*n(idx_HNO)
      pdj(37) =  &
          -k(658)*n(idx_CH3OH)  &
          -k(657)*n(idx_CH3OH)
      pdj(38) =  &
          -k(668)*n(idx_CO2)  &
          -k(669)*n(idx_CO2)  &
          -k(666)*n(idx_CO2)  &
          -k(667)*n(idx_CO2)
      pdj(40) =  &
          -k(676)*n(idx_H2SIO)
      pdj(44) =  &
          -k(699)*n(idx_OCN)  &
          -k(698)*n(idx_OCN)
      pdj(70) =  &
          +k(672)*n(idx_H2CO)
      pdj(71) =  &
          +k(662)*n(idx_CH4)  &
          +k(688)*n(idx_HNO)  &
          +k(537)*n(idx_H2)  &
          +k(131)*n(idx_H)  &
          +k(675)*n(idx_H2O)
      pdj(73) =  &
          +k(654)*n(idx_CH2)  &
          +k(669)*n(idx_CO2)  &
          +k(703)*n(idx_SIC)  &
          +k(679)*n(idx_HCN)  &
          +k(663)*n(idx_CH)  &
          +k(670)*n(idx_CO)  &
          +k(140)*n(idx_C)  &
          +k(665)*n(idx_CN)  &
          +k(685)*n(idx_HNC)
      pdj(74) =  &
          +k(673)*n(idx_H2CO)  &
          +k(660)*n(idx_CH4)
      pdj(75) =  &
          +k(142)*n(idx_CH)  &
          +k(659)*n(idx_CH4)  &
          +k(656)*n(idx_CH3)  &
          +k(680)*n(idx_HCN)  &
          +k(683)*n(idx_HCO)  &
          +k(655)*n(idx_CH2)
      pdj(76) =  &
          +k(143)*n(idx_H2CO)
      pdj(78) =  &
          +k(146)*n(idx_NH3)
      pdj(79) =  &
          +k(687)*n(idx_HNO)
      pdj(80) =  &
          +k(708)*n(idx_SIH4)  &
          +k(711)*n(idx_SIO)  &
          +k(702)*n(idx_SIC)  &
          +k(704)*n(idx_SIH2)  &
          +k(148)*n(idx_SI)  &
          +k(710)*n(idx_SIH)
      pdj(81) =  &
          +k(701)*n(idx_SIC3)
      pdj(84) =  &
          +k(707)*n(idx_SIH3)
      pdj(86) =  &
          +k(698)*n(idx_OCN)  &
          +k(677)*n(idx_HCN)  &
          +k(684)*n(idx_HNC)
      pdj(87) =  &
          +k(681)*n(idx_HCO)  &
          +k(671)*n(idx_H2CO)  &
          +k(666)*n(idx_CO2)
      pdj(88) =  &
          +k(145)*n(idx_N2)
      pdj(89) =  &
          +k(147)*n(idx_O2)  &
          +k(668)*n(idx_CO2)
      pdj(90) =  &
          +k(144)*n(idx_H2O)
      pdj(91) =  &
          +k(693)*n(idx_NH3)
      pdj(92) =  &
          +k(695)*n(idx_NO)  &
          +k(700)*n(idx_OH)  &
          +k(712)*n(idx_SIO)  &
          +k(667)*n(idx_CO2)  &
          +k(699)*n(idx_OCN)  &
          +k(697)*n(idx_O2)
      pdj(93) =  &
          +k(657)*n(idx_CH3OH)  &
          +k(674)*n(idx_H2O)
      pdj(94) =  &
          +k(658)*n(idx_CH3OH)  &
          +k(661)*n(idx_CH4)
      pdj(95) =  &
          +k(141)*n(idx_CH4)
      pdj(96) =  &
          +k(689)*n(idx_N2)  &
          +k(690)*n(idx_NH2)  &
          +k(678)*n(idx_HCN)  &
          +k(664)*n(idx_CN)  &
          +k(696)*n(idx_NO)  &
          +k(694)*n(idx_NH)
      pdj(98) =  &
          +k(692)*n(idx_NH3)  &
          +k(686)*n(idx_HNC)  &
          +k(691)*n(idx_NH2)
      pdj(100) =  &
          +k(705)*n(idx_SIH2)  &
          +k(706)*n(idx_SIH3)  &
          +k(709)*n(idx_SIH4)
      pdj(102) =  &
          +k(116)*n(idx_H2)
      pdj(103) =  &
          -k(658)*n(idx_CH3OH)  &
          -k(697)*n(idx_O2)  &
          -k(144)*n(idx_H2O)  &
          -k(669)*n(idx_CO2)  &
          -k(707)*n(idx_SIH3)  &
          -k(700)*n(idx_OH)  &
          -k(677)*n(idx_HCN)  &
          -k(116)*n(idx_H2)  &
          -k(685)*n(idx_HNC)  &
          -k(143)*n(idx_H2CO)  &
          -k(683)*n(idx_HCO)  &
          -k(710)*n(idx_SIH)  &
          -k(696)*n(idx_NO)  &
          -k(682)*n(idx_HCO)  &
          -k(142)*n(idx_CH)  &
          -k(673)*n(idx_H2CO)  &
          -k(140)*n(idx_C)  &
          -k(667)*n(idx_CO2)  &
          -k(661)*n(idx_CH4)  &
          -k(694)*n(idx_NH)  &
          -k(698)*n(idx_OCN)  &
          -k(703)*n(idx_SIC)  &
          -k(699)*n(idx_OCN)  &
          -k(665)*n(idx_CN)  &
          -k(692)*n(idx_NH3)  &
          -k(655)*n(idx_CH2)  &
          -k(695)*n(idx_NO)  &
          -k(705)*n(idx_SIH2)  &
          -k(672)*n(idx_H2CO)  &
          -k(688)*n(idx_HNO)  &
          -k(537)*n(idx_H2)  &
          -k(668)*n(idx_CO2)  &
          -k(679)*n(idx_HCN)  &
          -k(674)*n(idx_H2O)  &
          -k(662)*n(idx_CH4)  &
          -k(660)*n(idx_CH4)  &
          -k(701)*n(idx_SIC3)  &
          -k(664)*n(idx_CN)  &
          -k(681)*n(idx_HCO)  &
          -k(657)*n(idx_CH3OH)  &
          -k(687)*n(idx_HNO)  &
          -k(708)*n(idx_SIH4)  &
          -k(689)*n(idx_N2)  &
          -k(684)*n(idx_HNC)  &
          -k(147)*n(idx_O2)  &
          -k(686)*n(idx_HNC)  &
          -k(670)*n(idx_CO)  &
          -k(711)*n(idx_SIO)  &
          -k(680)*n(idx_HCN)  &
          -k(706)*n(idx_SIH3)  &
          -k(141)*n(idx_CH4)  &
          -k(146)*n(idx_NH3)  &
          -k(691)*n(idx_NH2)  &
          -k(148)*n(idx_SI)  &
          -k(145)*n(idx_N2)  &
          -k(678)*n(idx_HCN)  &
          -k(704)*n(idx_SIH2)  &
          -k(675)*n(idx_H2O)  &
          -k(671)*n(idx_H2CO)  &
          -k(690)*n(idx_NH2)  &
          -k(654)*n(idx_CH2)  &
          -k(663)*n(idx_CH)  &
          -k(676)*n(idx_H2SIO)  &
          -k(1219)*n(idx_E)  &
          -k(666)*n(idx_CO2)  &
          -k(656)*n(idx_CH3)  &
          -k(702)*n(idx_SIC)  &
          -k(693)*n(idx_NH3)  &
          -k(712)*n(idx_SIO)  &
          -k(659)*n(idx_CH4)  &
          -k(709)*n(idx_SIH4)  &
          -k(131)*n(idx_H)
      pdj(111) =  &
          +k(682)*n(idx_HCO)
      pdj(115) =  &
          +k(676)*n(idx_H2SIO)
    elseif(j==104) then
      pdj(1) =  &
          -k(335)*n(idx_E)
      pdj(2) =  &
          -k(468)*n(idx_CH)
      pdj(4) =  &
          -k(650)*n(idx_HNC)
      pdj(5) =  &
          -k(631)*n(idx_HCN)
      pdj(7) =  &
          -k(389)*n(idx_C)
      pdj(8) =  &
          +k(335)*n(idx_E)
      pdj(9) =  &
          -k(569)*n(idx_H2O)
      pdj(10) =  &
          -k(857)*n(idx_OH)
      pdj(12) =  &
          -k(431)*n(idx_CH2)
      pdj(13) =  &
          -k(551)*n(idx_H2CO)
      pdj(14) =  &
          -k(643)*n(idx_HCO)
      pdj(17) =  &
          +k(487)*n(idx_CO)  &
          +k(551)*n(idx_H2CO)  &
          +k(468)*n(idx_CH)  &
          +k(653)*n(idx_CO2)  &
          +k(643)*n(idx_HCO)  &
          +k(650)*n(idx_HNC)  &
          +k(389)*n(idx_C)  &
          +k(789)*n(idx_NH2)  &
          +k(483)*n(idx_CN)  &
          +k(631)*n(idx_HCN)  &
          +k(569)*n(idx_H2O)  &
          +k(335)*n(idx_E)  &
          +k(801)*n(idx_NH)  &
          +k(431)*n(idx_CH2)  &
          -k(205)*n(idx_NO)  &
          +k(857)*n(idx_OH)  &
          +k(733)*n(idx_N2)
      pdj(24) =  &
          -k(483)*n(idx_CN)
      pdj(25) =  &
          -k(487)*n(idx_CO)
      pdj(26) =  &
          -k(733)*n(idx_N2)
      pdj(27) =  &
          -k(789)*n(idx_NH2)
      pdj(31) =  &
          -k(801)*n(idx_NH)
      pdj(36) =  &
          +k(205)*n(idx_NO)
      pdj(38) =  &
          -k(653)*n(idx_CO2)
      pdj(63) =  &
          +k(1296)
      pdj(70) =  &
          +k(487)*n(idx_CO)
      pdj(74) =  &
          +k(468)*n(idx_CH)
      pdj(75) =  &
          +k(389)*n(idx_C)
      pdj(76) =  &
          +k(643)*n(idx_HCO)
      pdj(78) =  &
          +k(789)*n(idx_NH2)
      pdj(79) =  &
          +k(205)*n(idx_NO)
      pdj(90) =  &
          +k(857)*n(idx_OH)
      pdj(91) =  &
          +k(801)*n(idx_NH)
      pdj(94) =  &
          +k(431)*n(idx_CH2)
      pdj(97) =  &
          +k(483)*n(idx_CN)
      pdj(104) =  &
          -k(389)*n(idx_C)  &
          -k(733)*n(idx_N2)  &
          -k(631)*n(idx_HCN)  &
          -k(468)*n(idx_CH)  &
          -k(1296)  &
          -k(487)*n(idx_CO)  &
          -k(569)*n(idx_H2O)  &
          -k(551)*n(idx_H2CO)  &
          -k(801)*n(idx_NH)  &
          -k(431)*n(idx_CH2)  &
          -k(643)*n(idx_HCO)  &
          -k(335)*n(idx_E)  &
          -k(789)*n(idx_NH2)  &
          -k(857)*n(idx_OH)  &
          -k(653)*n(idx_CO2)  &
          -k(205)*n(idx_NO)  &
          -k(483)*n(idx_CN)  &
          -k(650)*n(idx_HNC)
      pdj(107) =  &
          +k(551)*n(idx_H2CO)
      pdj(108) =  &
          +k(569)*n(idx_H2O)
      pdj(109) =  &
          +k(650)*n(idx_HNC)  &
          +k(631)*n(idx_HCN)
      pdj(110) =  &
          +k(653)*n(idx_CO2)
      pdj(112) =  &
          +k(733)*n(idx_N2)
    elseif(j==105) then
      pdj(1) =  &
          -k(312)*n(idx_E)  &
          -k(311)*n(idx_E)
      pdj(6) =  &
          +k(312)*n(idx_E)
      pdj(8) =  &
          +k(311)*n(idx_E)  &
          +k(1297)
      pdj(17) =  &
          +k(312)*n(idx_E)
      pdj(36) =  &
          +k(311)*n(idx_E)
      pdj(63) =  &
          +k(1297)
      pdj(105) =  &
          -k(312)*n(idx_E)  &
          -k(311)*n(idx_E)  &
          -k(1297)
    elseif(j==106) then
      pdj(1) =  &
          -k(317)*n(idx_E)  &
          -k(316)*n(idx_E)
      pdj(2) =  &
          -k(581)*n(idx_CH)
      pdj(3) =  &
          -k(599)*n(idx_O)  &
          -k(600)*n(idx_O)
      pdj(4) =  &
          -k(590)*n(idx_HNC)
      pdj(5) =  &
          -k(588)*n(idx_HCN)
      pdj(6) =  &
          +k(581)*n(idx_CH)  &
          +k(593)*n(idx_N2)  &
          +k(606)*n(idx_SIH)  &
          +k(577)*n(idx_C)  &
          +k(598)*n(idx_O2)  &
          +k(591)*n(idx_HNO)  &
          +k(600)*n(idx_O)  &
          +k(585)*n(idx_CO)  &
          +k(595)*n(idx_NH)  &
          +k(596)*n(idx_NO2)  &
          +k(588)*n(idx_HCN)  &
          +k(592)*n(idx_MG)  &
          +k(586)*n(idx_H2CO)  &
          +k(579)*n(idx_CH3)  &
          +k(597)*n(idx_NO)  &
          +k(316)*n(idx_E)  &
          +k(584)*n(idx_CO)  &
          +k(578)*n(idx_CH2)  &
          +k(607)*n(idx_SIO)  &
          +k(605)*n(idx_SIH4)  &
          +k(587)*n(idx_H2O)  &
          +k(580)*n(idx_CH3OH)  &
          +k(1147)  &
          +k(590)*n(idx_HNC)  &
          +k(603)*n(idx_SIH2)  &
          +k(602)*n(idx_SI)  &
          +k(589)*n(idx_HCO)  &
          +k(604)*n(idx_SIH3)  &
          +k(601)*n(idx_OH)  &
          +k(583)*n(idx_CO2)  &
          +k(594)*n(idx_NH2)  &
          +k(582)*n(idx_CN)
      pdj(7) =  &
          -k(577)*n(idx_C)
      pdj(8) =  &
          +3.d0*k(317)*n(idx_E)  &
          +k(316)*n(idx_E)  &
          +k(1146)  &
          +k(599)*n(idx_O)  &
          +k(592)*n(idx_MG)
      pdj(9) =  &
          -k(587)*n(idx_H2O)  &
          +k(580)*n(idx_CH3OH)
      pdj(10) =  &
          +k(596)*n(idx_NO2)  &
          -k(601)*n(idx_OH)
      pdj(11) =  &
          -k(598)*n(idx_O2)
      pdj(12) =  &
          -k(578)*n(idx_CH2)
      pdj(13) =  &
          -k(586)*n(idx_H2CO)
      pdj(14) =  &
          -k(589)*n(idx_HCO)
      pdj(15) =  &
          -k(592)*n(idx_MG)
      pdj(17) =  &
          -k(597)*n(idx_NO)
      pdj(18) =  &
          -k(602)*n(idx_SI)
      pdj(22) =  &
          -k(603)*n(idx_SIH2)
      pdj(23) =  &
          -k(604)*n(idx_SIH3)
      pdj(24) =  &
          -k(582)*n(idx_CN)
      pdj(25) =  &
          -k(585)*n(idx_CO)  &
          -k(584)*n(idx_CO)
      pdj(26) =  &
          -k(593)*n(idx_N2)
      pdj(27) =  &
          -k(594)*n(idx_NH2)
      pdj(28) =  &
          -k(579)*n(idx_CH3)
      pdj(31) =  &
          -k(595)*n(idx_NH)
      pdj(32) =  &
          -k(605)*n(idx_SIH4)
      pdj(33) =  &
          -k(606)*n(idx_SIH)
      pdj(34) =  &
          -k(607)*n(idx_SIO)
      pdj(36) =  &
          -k(591)*n(idx_HNO)
      pdj(37) =  &
          -k(580)*n(idx_CH3OH)
      pdj(38) =  &
          -k(583)*n(idx_CO2)
      pdj(42) =  &
          -k(596)*n(idx_NO2)
      pdj(70) =  &
          +k(584)*n(idx_CO)
      pdj(71) =  &
          +k(1147)
      pdj(72) =  &
          +k(585)*n(idx_CO)
      pdj(74) =  &
          +k(581)*n(idx_CH)
      pdj(75) =  &
          +k(577)*n(idx_C)
      pdj(76) =  &
          +k(589)*n(idx_HCO)
      pdj(77) =  &
          +k(592)*n(idx_MG)
      pdj(78) =  &
          +k(594)*n(idx_NH2)
      pdj(79) =  &
          +k(596)*n(idx_NO2)
      pdj(84) =  &
          +k(606)*n(idx_SIH)
      pdj(85) =  &
          +k(603)*n(idx_SIH2)
      pdj(90) =  &
          +k(601)*n(idx_OH)  &
          +k(599)*n(idx_O)
      pdj(91) =  &
          +k(595)*n(idx_NH)
      pdj(93) =  &
          +k(600)*n(idx_O)
      pdj(94) =  &
          +k(578)*n(idx_CH2)  &
          +k(580)*n(idx_CH3OH)
      pdj(95) =  &
          +k(579)*n(idx_CH3)
      pdj(97) =  &
          +k(582)*n(idx_CN)
      pdj(99) =  &
          +k(604)*n(idx_SIH3)
      pdj(100) =  &
          +k(602)*n(idx_SI)
      pdj(102) =  &
          +k(1146)
      pdj(104) =  &
          +k(597)*n(idx_NO)
      pdj(105) =  &
          +k(591)*n(idx_HNO)
      pdj(106) =  &
          -k(584)*n(idx_CO)  &
          -k(579)*n(idx_CH3)  &
          -k(1147)  &
          -k(592)*n(idx_MG)  &
          -k(589)*n(idx_HCO)  &
          -k(591)*n(idx_HNO)  &
          -k(593)*n(idx_N2)  &
          -k(582)*n(idx_CN)  &
          -k(599)*n(idx_O)  &
          -k(585)*n(idx_CO)  &
          -k(597)*n(idx_NO)  &
          -k(607)*n(idx_SIO)  &
          -k(581)*n(idx_CH)  &
          -k(577)*n(idx_C)  &
          -k(600)*n(idx_O)  &
          -k(605)*n(idx_SIH4)  &
          -k(1146)  &
          -k(317)*n(idx_E)  &
          -k(587)*n(idx_H2O)  &
          -k(594)*n(idx_NH2)  &
          -k(588)*n(idx_HCN)  &
          -k(595)*n(idx_NH)  &
          -k(316)*n(idx_E)  &
          -k(583)*n(idx_CO2)  &
          -k(578)*n(idx_CH2)  &
          -k(603)*n(idx_SIH2)  &
          -k(601)*n(idx_OH)  &
          -k(580)*n(idx_CH3OH)  &
          -k(602)*n(idx_SI)  &
          -k(596)*n(idx_NO2)  &
          -k(586)*n(idx_H2CO)  &
          -k(604)*n(idx_SIH3)  &
          -k(590)*n(idx_HNC)  &
          -k(598)*n(idx_O2)  &
          -k(606)*n(idx_SIH)
      pdj(107) =  &
          +k(586)*n(idx_H2CO)
      pdj(108) =  &
          +k(587)*n(idx_H2O)
      pdj(109) =  &
          +k(588)*n(idx_HCN)  &
          +k(590)*n(idx_HNC)
      pdj(110) =  &
          +k(583)*n(idx_CO2)
      pdj(112) =  &
          +k(593)*n(idx_N2)
      pdj(113) =  &
          +k(598)*n(idx_O2)
      pdj(114) =  &
          +k(605)*n(idx_SIH4)
      pdj(115) =  &
          +k(607)*n(idx_SIO)
    elseif(j==107) then
      pdj(1) =  &
          -k(320)*n(idx_E)  &
          -k(322)*n(idx_E)  &
          -k(319)*n(idx_E)  &
          -k(318)*n(idx_E)  &
          -k(321)*n(idx_E)
      pdj(2) =  &
          +k(319)*n(idx_E)  &
          -k(462)*n(idx_CH)
      pdj(4) =  &
          -k(648)*n(idx_HNC)
      pdj(5) =  &
          -k(629)*n(idx_HCN)
      pdj(6) =  &
          +k(320)*n(idx_E)
      pdj(8) =  &
          +2.d0*k(322)*n(idx_E)  &
          +k(320)*n(idx_E)  &
          +k(321)*n(idx_E)
      pdj(9) =  &
          +k(319)*n(idx_E)  &
          -k(565)*n(idx_H2O)
      pdj(10) =  &
          +k(318)*n(idx_E)
      pdj(12) =  &
          +k(318)*n(idx_E)
      pdj(13) =  &
          +k(783)*n(idx_NH2)  &
          +k(462)*n(idx_CH)  &
          +k(321)*n(idx_E)  &
          +k(648)*n(idx_HNC)  &
          +k(629)*n(idx_HCN)  &
          +k(565)*n(idx_H2O)
      pdj(14) =  &
          +k(322)*n(idx_E)
      pdj(25) =  &
          +k(320)*n(idx_E)
      pdj(27) =  &
          -k(783)*n(idx_NH2)
      pdj(45) =  &
          +k(1250)
      pdj(74) =  &
          +k(462)*n(idx_CH)
      pdj(78) =  &
          +k(783)*n(idx_NH2)
      pdj(107) =  &
          -k(319)*n(idx_E)  &
          -k(318)*n(idx_E)  &
          -k(648)*n(idx_HNC)  &
          -k(565)*n(idx_H2O)  &
          -k(783)*n(idx_NH2)  &
          -k(322)*n(idx_E)  &
          -k(1250)  &
          -k(320)*n(idx_E)  &
          -k(321)*n(idx_E)  &
          -k(629)*n(idx_HCN)  &
          -k(462)*n(idx_CH)
      pdj(108) =  &
          +k(565)*n(idx_H2O)
      pdj(109) =  &
          +k(648)*n(idx_HNC)  &
          +k(629)*n(idx_HCN)
    elseif(j==108) then
      pdj(1) =  &
          -k(324)*n(idx_E)  &
          -k(325)*n(idx_E)  &
          -k(323)*n(idx_E)  &
          -k(326)*n(idx_E)
      pdj(2) =  &
          -k(463)*n(idx_CH)
      pdj(3) =  &
          +k(324)*n(idx_E)
      pdj(4) =  &
          -k(610)*n(idx_HNC)
      pdj(5) =  &
          -k(609)*n(idx_HCN)
      pdj(6) =  &
          +k(324)*n(idx_E)  &
          +k(325)*n(idx_E)  &
          +k(385)*n(idx_C)
      pdj(7) =  &
          -k(385)*n(idx_C)
      pdj(8) =  &
          +k(323)*n(idx_E)  &
          +k(324)*n(idx_E)  &
          +2.d0*k(326)*n(idx_E)  &
          +k(1287)
      pdj(9) =  &
          +k(609)*n(idx_HCN)  &
          +k(611)*n(idx_SI)  &
          +k(610)*n(idx_HNC)  &
          +k(426)*n(idx_CH2)  &
          +k(608)*n(idx_H2CO)  &
          +k(784)*n(idx_NH2)  &
          +k(323)*n(idx_E)  &
          +k(612)*n(idx_SIH2)  &
          +k(613)*n(idx_SIH)  &
          +k(614)*n(idx_SIO)  &
          +k(463)*n(idx_CH)
      pdj(10) =  &
          +k(325)*n(idx_E)  &
          +k(326)*n(idx_E)
      pdj(12) =  &
          -k(426)*n(idx_CH2)
      pdj(13) =  &
          -k(608)*n(idx_H2CO)
      pdj(18) =  &
          -k(611)*n(idx_SI)
      pdj(22) =  &
          -k(612)*n(idx_SIH2)
      pdj(27) =  &
          -k(784)*n(idx_NH2)
      pdj(33) =  &
          -k(613)*n(idx_SIH)
      pdj(34) =  &
          -k(614)*n(idx_SIO)
      pdj(55) =  &
          +k(1287)
      pdj(70) =  &
          +k(385)*n(idx_C)
      pdj(74) =  &
          +k(463)*n(idx_CH)
      pdj(78) =  &
          +k(784)*n(idx_NH2)
      pdj(84) =  &
          +k(613)*n(idx_SIH)
      pdj(85) =  &
          +k(612)*n(idx_SIH2)
      pdj(94) =  &
          +k(426)*n(idx_CH2)
      pdj(100) =  &
          +k(611)*n(idx_SI)
      pdj(107) =  &
          +k(608)*n(idx_H2CO)
      pdj(108) =  &
          -k(612)*n(idx_SIH2)  &
          -k(784)*n(idx_NH2)  &
          -k(1287)  &
          -k(614)*n(idx_SIO)  &
          -k(326)*n(idx_E)  &
          -k(463)*n(idx_CH)  &
          -k(324)*n(idx_E)  &
          -k(325)*n(idx_E)  &
          -k(613)*n(idx_SIH)  &
          -k(608)*n(idx_H2CO)  &
          -k(609)*n(idx_HCN)  &
          -k(426)*n(idx_CH2)  &
          -k(610)*n(idx_HNC)  &
          -k(323)*n(idx_E)  &
          -k(385)*n(idx_C)  &
          -k(611)*n(idx_SI)
      pdj(109) =  &
          +k(609)*n(idx_HCN)  &
          +k(610)*n(idx_HNC)
      pdj(115) =  &
          +k(614)*n(idx_SIO)
    elseif(j==109) then
      pdj(1) =  &
          -k(328)*n(idx_E)  &
          -k(329)*n(idx_E)  &
          -k(330)*n(idx_E)
      pdj(2) =  &
          -k(466)*n(idx_CH)  &
          -k(465)*n(idx_CH)
      pdj(4) =  &
          +k(330)*n(idx_E)  &
          +k(787)*n(idx_NH2)  &
          +k(635)*n(idx_H2CO)  &
          +k(466)*n(idx_CH)  &
          +k(429)*n(idx_CH2)
      pdj(5) =  &
          +k(634)*n(idx_H2CO)  &
          +k(786)*n(idx_NH2)  &
          +k(428)*n(idx_CH2)  &
          +k(465)*n(idx_CH)  &
          +k(329)*n(idx_E)
      pdj(8) =  &
          +k(330)*n(idx_E)  &
          +k(1302)  &
          +2.d0*k(328)*n(idx_E)  &
          +k(329)*n(idx_E)
      pdj(12) =  &
          -k(429)*n(idx_CH2)  &
          -k(428)*n(idx_CH2)
      pdj(13) =  &
          -k(634)*n(idx_H2CO)  &
          -k(635)*n(idx_H2CO)
      pdj(24) =  &
          +k(328)*n(idx_E)
      pdj(27) =  &
          -k(786)*n(idx_NH2)  &
          -k(787)*n(idx_NH2)
      pdj(59) =  &
          +k(1302)
      pdj(74) =  &
          +k(466)*n(idx_CH)  &
          +k(465)*n(idx_CH)
      pdj(78) =  &
          +k(787)*n(idx_NH2)  &
          +k(786)*n(idx_NH2)
      pdj(94) =  &
          +k(429)*n(idx_CH2)  &
          +k(428)*n(idx_CH2)
      pdj(107) =  &
          +k(634)*n(idx_H2CO)  &
          +k(635)*n(idx_H2CO)
      pdj(109) =  &
          -k(786)*n(idx_NH2)  &
          -k(466)*n(idx_CH)  &
          -k(1302)  &
          -k(635)*n(idx_H2CO)  &
          -k(330)*n(idx_E)  &
          -k(787)*n(idx_NH2)  &
          -k(428)*n(idx_CH2)  &
          -k(634)*n(idx_H2CO)  &
          -k(429)*n(idx_CH2)  &
          -k(465)*n(idx_CH)  &
          -k(328)*n(idx_E)  &
          -k(329)*n(idx_E)
    elseif(j==110) then
      pdj(1) =  &
          -k(333)*n(idx_E)  &
          -k(332)*n(idx_E)  &
          -k(334)*n(idx_E)
      pdj(3) =  &
          -k(825)*n(idx_O)  &
          +k(333)*n(idx_E)
      pdj(7) =  &
          -k(388)*n(idx_C)
      pdj(8) =  &
          +k(1288)  &
          +k(333)*n(idx_E)  &
          +k(332)*n(idx_E)
      pdj(9) =  &
          -k(568)*n(idx_H2O)
      pdj(10) =  &
          +k(334)*n(idx_E)
      pdj(11) =  &
          +k(825)*n(idx_O)
      pdj(25) =  &
          +k(334)*n(idx_E)  &
          +k(333)*n(idx_E)  &
          -k(486)*n(idx_CO)
      pdj(38) =  &
          +k(486)*n(idx_CO)  &
          +k(568)*n(idx_H2O)  &
          +k(388)*n(idx_C)  &
          +k(332)*n(idx_E)
      pdj(57) =  &
          +k(1288)
      pdj(70) =  &
          +k(825)*n(idx_O)  &
          +k(486)*n(idx_CO)
      pdj(75) =  &
          +k(388)*n(idx_C)
      pdj(108) =  &
          +k(568)*n(idx_H2O)
      pdj(110) =  &
          -k(825)*n(idx_O)  &
          -k(1288)  &
          -k(333)*n(idx_E)  &
          -k(332)*n(idx_E)  &
          -k(568)*n(idx_H2O)  &
          -k(334)*n(idx_E)  &
          -k(388)*n(idx_C)  &
          -k(486)*n(idx_CO)
    elseif(j==111) then
      pdj(1) =  &
          -k(337)*n(idx_E)
      pdj(6) =  &
          -k(538)*n(idx_H2)
      pdj(8) =  &
          +k(337)*n(idx_E)  &
          -k(619)*n(idx_H)
      pdj(35) =  &
          +k(619)*n(idx_H)  &
          +k(337)*n(idx_E)  &
          +k(538)*n(idx_H2)
      pdj(102) =  &
          +k(619)*n(idx_H)
      pdj(106) =  &
          +k(538)*n(idx_H2)
      pdj(111) =  &
          -k(619)*n(idx_H)  &
          -k(337)*n(idx_E)  &
          -k(538)*n(idx_H2)
    elseif(j==112) then
      pdj(1) =  &
          -k(340)*n(idx_E)  &
          -k(339)*n(idx_E)
      pdj(2) =  &
          -k(470)*n(idx_CH)
      pdj(3) =  &
          -k(827)*n(idx_O)
      pdj(4) =  &
          -k(651)*n(idx_HNC)
      pdj(5) =  &
          -k(632)*n(idx_HCN)
      pdj(7) =  &
          -k(390)*n(idx_C)
      pdj(8) =  &
          +k(339)*n(idx_E)  &
          +k(1303)
      pdj(9) =  &
          -k(571)*n(idx_H2O)
      pdj(10) =  &
          -k(858)*n(idx_OH)
      pdj(12) =  &
          -k(432)*n(idx_CH2)
      pdj(13) =  &
          -k(736)*n(idx_H2CO)
      pdj(14) =  &
          -k(644)*n(idx_HCO)
      pdj(25) =  &
          -k(488)*n(idx_CO)
      pdj(26) =  &
          +k(390)*n(idx_C)  &
          +k(790)*n(idx_NH2)  &
          +k(571)*n(idx_H2O)  &
          +k(488)*n(idx_CO)  &
          +k(644)*n(idx_HCO)  &
          +k(858)*n(idx_OH)  &
          +k(632)*n(idx_HCN)  &
          +k(827)*n(idx_O)  &
          +k(735)*n(idx_CO2)  &
          +k(736)*n(idx_H2CO)  &
          +k(802)*n(idx_NH)  &
          +k(651)*n(idx_HNC)  &
          +k(339)*n(idx_E)  &
          +k(432)*n(idx_CH2)  &
          +k(470)*n(idx_CH)
      pdj(27) =  &
          -k(790)*n(idx_NH2)
      pdj(30) =  &
          +k(340)*n(idx_E)
      pdj(31) =  &
          +k(340)*n(idx_E)  &
          -k(802)*n(idx_NH)
      pdj(38) =  &
          -k(735)*n(idx_CO2)
      pdj(58) =  &
          +k(1303)
      pdj(70) =  &
          +k(488)*n(idx_CO)
      pdj(74) =  &
          +k(470)*n(idx_CH)
      pdj(75) =  &
          +k(390)*n(idx_C)
      pdj(76) =  &
          +k(644)*n(idx_HCO)
      pdj(78) =  &
          +k(790)*n(idx_NH2)
      pdj(90) =  &
          +k(858)*n(idx_OH)
      pdj(91) =  &
          +k(802)*n(idx_NH)
      pdj(93) =  &
          +k(827)*n(idx_O)
      pdj(94) =  &
          +k(432)*n(idx_CH2)
      pdj(107) =  &
          +k(736)*n(idx_H2CO)
      pdj(108) =  &
          +k(571)*n(idx_H2O)
      pdj(109) =  &
          +k(632)*n(idx_HCN)  &
          +k(651)*n(idx_HNC)
      pdj(110) =  &
          +k(735)*n(idx_CO2)
      pdj(112) =  &
          -k(632)*n(idx_HCN)  &
          -k(858)*n(idx_OH)  &
          -k(827)*n(idx_O)  &
          -k(432)*n(idx_CH2)  &
          -k(1303)  &
          -k(571)*n(idx_H2O)  &
          -k(735)*n(idx_CO2)  &
          -k(644)*n(idx_HCO)  &
          -k(390)*n(idx_C)  &
          -k(470)*n(idx_CH)  &
          -k(651)*n(idx_HNC)  &
          -k(339)*n(idx_E)  &
          -k(488)*n(idx_CO)  &
          -k(790)*n(idx_NH2)  &
          -k(736)*n(idx_H2CO)  &
          -k(802)*n(idx_NH)  &
          -k(340)*n(idx_E)
    elseif(j==113) then
      pdj(1) =  &
          -k(348)*n(idx_E)
      pdj(2) =  &
          -k(475)*n(idx_CH)
      pdj(3) =  &
          -k(830)*n(idx_O)
      pdj(4) =  &
          -k(652)*n(idx_HNC)
      pdj(5) =  &
          -k(633)*n(idx_HCN)
      pdj(6) =  &
          -k(545)*n(idx_H2)
      pdj(7) =  &
          -k(393)*n(idx_C)
      pdj(8) =  &
          +k(348)*n(idx_E)
      pdj(9) =  &
          -k(572)*n(idx_H2O)
      pdj(10) =  &
          -k(859)*n(idx_OH)
      pdj(11) =  &
          +k(545)*n(idx_H2)  &
          +k(437)*n(idx_CH2)  &
          +k(475)*n(idx_CH)  &
          +k(822)*n(idx_CO2)  &
          +k(646)*n(idx_HCO)  &
          +k(553)*n(idx_H2CO)  &
          +k(348)*n(idx_E)  &
          +k(489)*n(idx_CO)  &
          +k(791)*n(idx_NH2)  &
          +k(830)*n(idx_O)  &
          +k(652)*n(idx_HNC)  &
          +k(572)*n(idx_H2O)  &
          +k(734)*n(idx_N2)  &
          +k(393)*n(idx_C)  &
          +k(633)*n(idx_HCN)  &
          +k(808)*n(idx_NO)  &
          +k(859)*n(idx_OH)  &
          +k(484)*n(idx_CN)  &
          +k(806)*n(idx_NH)
      pdj(12) =  &
          -k(437)*n(idx_CH2)
      pdj(13) =  &
          -k(553)*n(idx_H2CO)
      pdj(14) =  &
          -k(646)*n(idx_HCO)
      pdj(17) =  &
          -k(808)*n(idx_NO)
      pdj(24) =  &
          -k(484)*n(idx_CN)
      pdj(25) =  &
          -k(489)*n(idx_CO)
      pdj(26) =  &
          -k(734)*n(idx_N2)
      pdj(27) =  &
          -k(791)*n(idx_NH2)
      pdj(31) =  &
          -k(806)*n(idx_NH)
      pdj(38) =  &
          -k(822)*n(idx_CO2)
      pdj(64) =  &
          +k(1298)
      pdj(70) =  &
          +k(489)*n(idx_CO)
      pdj(74) =  &
          +k(475)*n(idx_CH)
      pdj(75) =  &
          +k(393)*n(idx_C)
      pdj(76) =  &
          +k(646)*n(idx_HCO)
      pdj(78) =  &
          +k(791)*n(idx_NH2)
      pdj(90) =  &
          +k(859)*n(idx_OH)
      pdj(91) =  &
          +k(806)*n(idx_NH)
      pdj(93) =  &
          +k(830)*n(idx_O)
      pdj(94) =  &
          +k(437)*n(idx_CH2)
      pdj(97) =  &
          +k(484)*n(idx_CN)
      pdj(104) =  &
          +k(808)*n(idx_NO)
      pdj(106) =  &
          +k(545)*n(idx_H2)
      pdj(107) =  &
          +k(553)*n(idx_H2CO)
      pdj(108) =  &
          +k(572)*n(idx_H2O)
      pdj(109) =  &
          +k(633)*n(idx_HCN)  &
          +k(652)*n(idx_HNC)
      pdj(110) =  &
          +k(822)*n(idx_CO2)
      pdj(112) =  &
          +k(734)*n(idx_N2)
      pdj(113) =  &
          -k(553)*n(idx_H2CO)  &
          -k(806)*n(idx_NH)  &
          -k(652)*n(idx_HNC)  &
          -k(348)*n(idx_E)  &
          -k(633)*n(idx_HCN)  &
          -k(489)*n(idx_CO)  &
          -k(646)*n(idx_HCO)  &
          -k(1298)  &
          -k(822)*n(idx_CO2)  &
          -k(437)*n(idx_CH2)  &
          -k(545)*n(idx_H2)  &
          -k(572)*n(idx_H2O)  &
          -k(734)*n(idx_N2)  &
          -k(393)*n(idx_C)  &
          -k(808)*n(idx_NO)  &
          -k(484)*n(idx_CN)  &
          -k(830)*n(idx_O)  &
          -k(859)*n(idx_OH)  &
          -k(791)*n(idx_NH2)  &
          -k(475)*n(idx_CH)
    elseif(j==114) then
      pdj(1) =  &
          -k(361)*n(idx_E)  &
          -k(362)*n(idx_E)
      pdj(6) =  &
          +k(361)*n(idx_E)
      pdj(8) =  &
          +k(1249)  &
          +k(362)*n(idx_E)
      pdj(9) =  &
          -k(576)*n(idx_H2O)
      pdj(23) =  &
          +k(361)*n(idx_E)
      pdj(32) =  &
          +k(576)*n(idx_H2O)  &
          +k(362)*n(idx_E)
      pdj(48) =  &
          +k(1249)
      pdj(108) =  &
          +k(576)*n(idx_H2O)
      pdj(114) =  &
          -k(361)*n(idx_E)  &
          -k(576)*n(idx_H2O)  &
          -k(362)*n(idx_E)  &
          -k(1249)
    elseif(j==115) then
      pdj(1) =  &
          -k(364)*n(idx_E)  &
          -k(365)*n(idx_E)
      pdj(8) =  &
          +k(365)*n(idx_E)
      pdj(10) =  &
          +k(364)*n(idx_E)
      pdj(18) =  &
          +k(364)*n(idx_E)
      pdj(34) =  &
          +k(365)*n(idx_E)
      pdj(49) =  &
          +k(1247)
      pdj(115) =  &
          -k(1247)  &
          -k(364)*n(idx_E)  &
          -k(365)*n(idx_E)
    elseif(j==116) then
    elseif(j==117) then
    elseif(j==118) then

    elseif(j==119) then
    end if

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    use krome_gadiab
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

    !d[E_dot]/d[E]
    pd(1,1) =  &
        -k(327)*n(idx_HCNj)  &
        -k(1218)*n(idx_H2COj)  &
        -k(328)*n(idx_HCNHj)  &
        -k(304)*n(idx_CNj)  &
        -k(357)*n(idx_SIH3j)  &
        -k(323)*n(idx_H3Oj)  &
        -k(345)*n(idx_NH3j)  &
        -k(295)*n(idx_CHj)  &
        -k(356)*n(idx_SIH2j)  &
        -k(313)*n(idx_H2Oj)  &
        -k(347)*n(idx_O2j)  &
        -k(1220)*n(idx_MGj)  &
        -k(1216)*n(idx_CH3j)  &
        -k(311)*n(idx_H2NOj)  &
        -k(1223)*n(idx_SIj)  &
        -k(300)*n(idx_CH3j)  &
        -k(335)*n(idx_HNOj)  &
        -k(341)*n(idx_NHj)  &
        -k(358)*n(idx_SIH3j)  &
        -k(352)*n(idx_SIC3j)  &
        -k(326)*n(idx_H3Oj)  &
        -k(354)*n(idx_SIH2j)  &
        -k(305)*n(idx_COj)  &
        -k(325)*n(idx_H3Oj)  &
        -k(339)*n(idx_N2Hj)  &
        -k(309)*n(idx_H2COj)  &
        -k(343)*n(idx_NH2j)  &
        -k(310)*n(idx_H2COj)  &
        -k(359)*n(idx_SIH4j)  &
        -k(355)*n(idx_SIH2j)  &
        -k(1307)  &
        -k(296)*n(idx_CH2j)  &
        -k(349)*n(idx_OHj)  &
        -k(308)*n(idx_H2COj)  &
        -k(9)*n(idx_H2)  &
        -k(1219)*n(idx_HEj)  &
        +k(9)*n(idx_H2)  &
        -k(329)*n(idx_HCNHj)  &
        -k(1221)*n(idx_Nj)  &
        -k(307)*n(idx_H2COj)  &
        -k(333)*n(idx_HCO2j)  &
        -k(1215)*n(idx_Cj)  &
        -k(337)*n(idx_HEHj)  &
        -k(346)*n(idx_NOj)  &
        -k(353)*n(idx_SIHj)  &
        -k(314)*n(idx_H2Oj)  &
        -k(351)*n(idx_SIC2j)  &
        -k(330)*n(idx_HCNHj)  &
        -k(312)*n(idx_H2NOj)  &
        -k(297)*n(idx_CH2j)  &
        -k(320)*n(idx_H3COj)  &
        -k(321)*n(idx_H3COj)  &
        -k(342)*n(idx_NH2j)  &
        -k(298)*n(idx_CH2j)  &
        -k(361)*n(idx_SIH5j)  &
        -k(322)*n(idx_H3COj)  &
        -k(301)*n(idx_CH3j)  &
        -k(348)*n(idx_O2Hj)  &
        -k(350)*n(idx_SICj)  &
        -k(344)*n(idx_NH3j)  &
        -k(316)*n(idx_H3j)  &
        -k(319)*n(idx_H3COj)  &
        -k(318)*n(idx_H3COj)  &
        -k(334)*n(idx_HCO2j)  &
        -k(1217)*n(idx_Hj)  &
        -k(365)*n(idx_SIOHj)  &
        -k(331)*n(idx_HCOj)  &
        -k(302)*n(idx_CH4j)  &
        -k(363)*n(idx_SIOj)  &
        -k(324)*n(idx_H3Oj)  &
        -k(360)*n(idx_SIH4j)  &
        -k(317)*n(idx_H3j)  &
        -k(1222)*n(idx_Oj)  &
        -k(340)*n(idx_N2Hj)  &
        -k(338)*n(idx_N2j)  &
        -k(303)*n(idx_CH4j)  &
        -k(332)*n(idx_HCO2j)  &
        -k(306)*n(idx_H2j)  &
        -k(362)*n(idx_SIH5j)  &
        -k(336)*n(idx_HOCj)  &
        -k(364)*n(idx_SIOHj)  &
        -k(299)*n(idx_CH3j)  &
        -k(315)*n(idx_H2Oj)

    !d[CH_dot]/d[E]
    pd(2,1) =  &
        +k(298)*n(idx_CH2j)  &
        +k(300)*n(idx_CH3j)  &
        +k(319)*n(idx_H3COj)  &
        +k(301)*n(idx_CH3j)

    !d[O_dot]/d[E]
    pd(3,1) =  &
        +k(313)*n(idx_H2Oj)  &
        +k(314)*n(idx_H2Oj)  &
        +k(349)*n(idx_OHj)  &
        +k(346)*n(idx_NOj)  &
        +k(307)*n(idx_H2COj)  &
        +2.d0*k(347)*n(idx_O2j)  &
        +k(305)*n(idx_COj)  &
        +k(333)*n(idx_HCO2j)  &
        +k(1222)*n(idx_Oj)  &
        +k(363)*n(idx_SIOj)  &
        +k(324)*n(idx_H3Oj)

    !d[HNC_dot]/d[E]
    pd(4,1) =  &
        +k(330)*n(idx_HCNHj)

    !d[HCN_dot]/d[E]
    pd(5,1) =  &
        +k(329)*n(idx_HCNHj)

    !d[H2_dot]/d[E]
    pd(6,1) =  &
        +k(313)*n(idx_H2Oj)  &
        +k(325)*n(idx_H3Oj)  &
        +k(300)*n(idx_CH3j)  &
        +k(354)*n(idx_SIH2j)  &
        +k(361)*n(idx_SIH5j)  &
        +k(316)*n(idx_H3j)  &
        -k(9)*n(idx_H2)  &
        +k(312)*n(idx_H2NOj)  &
        +k(358)*n(idx_SIH3j)  &
        +k(320)*n(idx_H3COj)  &
        +k(359)*n(idx_SIH4j)  &
        +k(296)*n(idx_CH2j)  &
        +k(308)*n(idx_H2COj)  &
        +k(324)*n(idx_H3Oj)

    !d[C_dot]/d[E]
    pd(7,1) =  &
        +k(296)*n(idx_CH2j)  &
        +k(297)*n(idx_CH2j)  &
        +k(352)*n(idx_SIC3j)  &
        +k(304)*n(idx_CNj)  &
        +k(305)*n(idx_COj)  &
        +k(1215)*n(idx_Cj)  &
        +k(295)*n(idx_CHj)  &
        +k(350)*n(idx_SICj)  &
        +k(351)*n(idx_SIC2j)

    !d[H_dot]/d[E]
    pd(8,1) =  &
        +2.d0*k(328)*n(idx_HCNHj)  &
        +k(299)*n(idx_CH3j)  &
        +2.d0*k(301)*n(idx_CH3j)  &
        +k(357)*n(idx_SIH3j)  &
        +k(329)*n(idx_HCNHj)  &
        +k(331)*n(idx_HCOj)  &
        +k(333)*n(idx_HCO2j)  &
        +2.d0*k(326)*n(idx_H3Oj)  &
        +k(353)*n(idx_SIHj)  &
        +2.d0*k(302)*n(idx_CH4j)  &
        +2.d0*k(342)*n(idx_NH2j)  &
        +k(344)*n(idx_NH3j)  &
        +2.d0*k(297)*n(idx_CH2j)  &
        +k(349)*n(idx_OHj)  &
        +k(356)*n(idx_SIH2j)  &
        +k(365)*n(idx_SIOHj)  &
        +k(360)*n(idx_SIH4j)  &
        +2.d0*k(345)*n(idx_NH3j)  &
        +k(330)*n(idx_HCNHj)  &
        +k(310)*n(idx_H2COj)  &
        +k(348)*n(idx_O2Hj)  &
        +2.d0*k(9)*n(idx_H2)  &
        +k(336)*n(idx_HOCj)  &
        +k(303)*n(idx_CH4j)  &
        +k(335)*n(idx_HNOj)  &
        +k(324)*n(idx_H3Oj)  &
        +k(343)*n(idx_NH2j)  &
        +k(298)*n(idx_CH2j)  &
        +k(339)*n(idx_N2Hj)  &
        +k(362)*n(idx_SIH5j)  &
        +k(316)*n(idx_H3j)  &
        +k(332)*n(idx_HCO2j)  &
        +k(321)*n(idx_H3COj)  &
        +3.d0*k(317)*n(idx_H3j)  &
        +2.d0*k(355)*n(idx_SIH2j)  &
        +2.d0*k(306)*n(idx_H2j)  &
        +k(341)*n(idx_NHj)  &
        +2.d0*k(314)*n(idx_H2Oj)  &
        +k(315)*n(idx_H2Oj)  &
        +2.d0*k(322)*n(idx_H3COj)  &
        +k(327)*n(idx_HCNj)  &
        +k(337)*n(idx_HEHj)  &
        +2.d0*k(309)*n(idx_H2COj)  &
        +k(295)*n(idx_CHj)  &
        +k(320)*n(idx_H3COj)  &
        +k(1217)*n(idx_Hj)  &
        +k(311)*n(idx_H2NOj)  &
        +k(323)*n(idx_H3Oj)

    !d[H2O_dot]/d[E]
    pd(9,1) =  &
        +k(319)*n(idx_H3COj)  &
        +k(323)*n(idx_H3Oj)

    !d[OH_dot]/d[E]
    pd(10,1) =  &
        +k(334)*n(idx_HCO2j)  &
        +k(364)*n(idx_SIOHj)  &
        +k(315)*n(idx_H2Oj)  &
        +k(326)*n(idx_H3Oj)  &
        +k(325)*n(idx_H3Oj)  &
        +k(318)*n(idx_H3COj)

    !d[O2_dot]/d[E]
    pd(11,1) =  &
        +k(348)*n(idx_O2Hj)

    !d[CH2_dot]/d[E]
    pd(12,1) =  &
        +k(318)*n(idx_H3COj)  &
        +k(302)*n(idx_CH4j)  &
        +k(307)*n(idx_H2COj)  &
        +k(299)*n(idx_CH3j)

    !d[H2CO_dot]/d[E]
    pd(13,1) =  &
        +k(1218)*n(idx_H2COj)  &
        +k(321)*n(idx_H3COj)

    !d[HCO_dot]/d[E]
    pd(14,1) =  &
        +k(322)*n(idx_H3COj)  &
        +k(310)*n(idx_H2COj)

    !d[MG_dot]/d[E]
    pd(15,1) =  &
        +k(1220)*n(idx_MGj)

    !d[NO_dot]/d[E]
    pd(17,1) =  &
        +k(312)*n(idx_H2NOj)  &
        +k(335)*n(idx_HNOj)

    !d[SI_dot]/d[E]
    pd(18,1) =  &
        +k(364)*n(idx_SIOHj)  &
        +k(353)*n(idx_SIHj)  &
        +k(354)*n(idx_SIH2j)  &
        +k(1223)*n(idx_SIj)  &
        +k(363)*n(idx_SIOj)  &
        +k(350)*n(idx_SICj)  &
        +k(355)*n(idx_SIH2j)

    !d[SIC2_dot]/d[E]
    pd(19,1) =  &
        +k(352)*n(idx_SIC3j)

    !d[SIC_dot]/d[E]
    pd(21,1) =  &
        +k(351)*n(idx_SIC2j)

    !d[SIH2_dot]/d[E]
    pd(22,1) =  &
        +k(357)*n(idx_SIH3j)  &
        +k(359)*n(idx_SIH4j)

    !d[SIH3_dot]/d[E]
    pd(23,1) =  &
        +k(360)*n(idx_SIH4j)  &
        +k(361)*n(idx_SIH5j)

    !d[CN_dot]/d[E]
    pd(24,1) =  &
        +k(327)*n(idx_HCNj)  &
        +k(328)*n(idx_HCNHj)

    !d[CO_dot]/d[E]
    pd(25,1) =  &
        +k(334)*n(idx_HCO2j)  &
        +k(333)*n(idx_HCO2j)  &
        +k(309)*n(idx_H2COj)  &
        +k(320)*n(idx_H3COj)  &
        +k(331)*n(idx_HCOj)  &
        +k(336)*n(idx_HOCj)  &
        +k(308)*n(idx_H2COj)

    !d[N2_dot]/d[E]
    pd(26,1) =  &
        +k(339)*n(idx_N2Hj)

    !d[NH2_dot]/d[E]
    pd(27,1) =  &
        +k(344)*n(idx_NH3j)

    !d[CH3_dot]/d[E]
    pd(28,1) =  &
        +k(1216)*n(idx_CH3j)  &
        +k(303)*n(idx_CH4j)

    !d[N_dot]/d[E]
    pd(30,1) =  &
        +2.d0*k(338)*n(idx_N2j)  &
        +k(341)*n(idx_NHj)  &
        +k(346)*n(idx_NOj)  &
        +k(1221)*n(idx_Nj)  &
        +k(304)*n(idx_CNj)  &
        +k(342)*n(idx_NH2j)  &
        +k(340)*n(idx_N2Hj)

    !d[NH_dot]/d[E]
    pd(31,1) =  &
        +k(345)*n(idx_NH3j)  &
        +k(340)*n(idx_N2Hj)  &
        +k(343)*n(idx_NH2j)

    !d[SIH4_dot]/d[E]
    pd(32,1) =  &
        +k(362)*n(idx_SIH5j)

    !d[SIH_dot]/d[E]
    pd(33,1) =  &
        +k(356)*n(idx_SIH2j)  &
        +k(358)*n(idx_SIH3j)

    !d[SIO_dot]/d[E]
    pd(34,1) =  &
        +k(365)*n(idx_SIOHj)

    !d[HE_dot]/d[E]
    pd(35,1) =  &
        +k(337)*n(idx_HEHj)  &
        +k(1219)*n(idx_HEj)

    !d[HNO_dot]/d[E]
    pd(36,1) =  &
        +k(311)*n(idx_H2NOj)

    !d[CO2_dot]/d[E]
    pd(38,1) =  &
        +k(332)*n(idx_HCO2j)

    !d[E_DUST_dot]/d[E]
    pd(68,1) =  &
        +k(1307)

    !d[HCO+_dot]/d[E]
    pd(70,1) =  &
        -k(331)*n(idx_HCOj)

    !d[H+_dot]/d[E]
    pd(71,1) =  &
        -k(1217)*n(idx_Hj)

    !d[HOC+_dot]/d[E]
    pd(72,1) =  &
        -k(336)*n(idx_HOCj)

    !d[C+_dot]/d[E]
    pd(73,1) =  &
        -k(1215)*n(idx_Cj)

    !d[CH2+_dot]/d[E]
    pd(74,1) =  &
        -k(296)*n(idx_CH2j)  &
        -k(297)*n(idx_CH2j)  &
        -k(298)*n(idx_CH2j)

    !d[CH+_dot]/d[E]
    pd(75,1) =  &
        -k(295)*n(idx_CHj)

    !d[H2CO+_dot]/d[E]
    pd(76,1) =  &
        -k(310)*n(idx_H2COj)  &
        -k(1218)*n(idx_H2COj)  &
        -k(307)*n(idx_H2COj)  &
        -k(309)*n(idx_H2COj)  &
        -k(308)*n(idx_H2COj)

    !d[MG+_dot]/d[E]
    pd(77,1) =  &
        -k(1220)*n(idx_MGj)

    !d[NH3+_dot]/d[E]
    pd(78,1) =  &
        -k(344)*n(idx_NH3j)  &
        -k(345)*n(idx_NH3j)

    !d[NO+_dot]/d[E]
    pd(79,1) =  &
        -k(346)*n(idx_NOj)

    !d[SI+_dot]/d[E]
    pd(80,1) =  &
        -k(1223)*n(idx_SIj)

    !d[SIC2+_dot]/d[E]
    pd(81,1) =  &
        -k(351)*n(idx_SIC2j)

    !d[SIC3+_dot]/d[E]
    pd(82,1) =  &
        -k(352)*n(idx_SIC3j)

    !d[SIC+_dot]/d[E]
    pd(83,1) =  &
        -k(350)*n(idx_SICj)

    !d[SIH2+_dot]/d[E]
    pd(84,1) =  &
        -k(354)*n(idx_SIH2j)  &
        -k(355)*n(idx_SIH2j)  &
        -k(356)*n(idx_SIH2j)

    !d[SIH3+_dot]/d[E]
    pd(85,1) =  &
        -k(357)*n(idx_SIH3j)  &
        -k(358)*n(idx_SIH3j)

    !d[CN+_dot]/d[E]
    pd(86,1) =  &
        -k(304)*n(idx_CNj)

    !d[CO+_dot]/d[E]
    pd(87,1) =  &
        -k(305)*n(idx_COj)

    !d[N2+_dot]/d[E]
    pd(88,1) =  &
        -k(338)*n(idx_N2j)

    !d[O2+_dot]/d[E]
    pd(89,1) =  &
        -k(347)*n(idx_O2j)

    !d[H2O+_dot]/d[E]
    pd(90,1) =  &
        -k(314)*n(idx_H2Oj)  &
        -k(313)*n(idx_H2Oj)  &
        -k(315)*n(idx_H2Oj)

    !d[NH2+_dot]/d[E]
    pd(91,1) =  &
        -k(343)*n(idx_NH2j)  &
        -k(342)*n(idx_NH2j)

    !d[O+_dot]/d[E]
    pd(92,1) =  &
        -k(1222)*n(idx_Oj)

    !d[OH+_dot]/d[E]
    pd(93,1) =  &
        -k(349)*n(idx_OHj)

    !d[CH3+_dot]/d[E]
    pd(94,1) =  &
        -k(300)*n(idx_CH3j)  &
        -k(301)*n(idx_CH3j)  &
        -k(299)*n(idx_CH3j)  &
        -k(1216)*n(idx_CH3j)

    !d[CH4+_dot]/d[E]
    pd(95,1) =  &
        -k(303)*n(idx_CH4j)  &
        -k(302)*n(idx_CH4j)

    !d[N+_dot]/d[E]
    pd(96,1) =  &
        -k(1221)*n(idx_Nj)

    !d[HCN+_dot]/d[E]
    pd(97,1) =  &
        -k(327)*n(idx_HCNj)

    !d[NH+_dot]/d[E]
    pd(98,1) =  &
        -k(341)*n(idx_NHj)

    !d[SIH4+_dot]/d[E]
    pd(99,1) =  &
        -k(359)*n(idx_SIH4j)  &
        -k(360)*n(idx_SIH4j)

    !d[SIH+_dot]/d[E]
    pd(100,1) =  &
        -k(353)*n(idx_SIHj)

    !d[SIO+_dot]/d[E]
    pd(101,1) =  &
        -k(363)*n(idx_SIOj)

    !d[H2+_dot]/d[E]
    pd(102,1) =  &
        -k(306)*n(idx_H2j)

    !d[HE+_dot]/d[E]
    pd(103,1) =  &
        -k(1219)*n(idx_HEj)

    !d[HNO+_dot]/d[E]
    pd(104,1) =  &
        -k(335)*n(idx_HNOj)

    !d[H2NO+_dot]/d[E]
    pd(105,1) =  &
        -k(312)*n(idx_H2NOj)  &
        -k(311)*n(idx_H2NOj)

    !d[H3+_dot]/d[E]
    pd(106,1) =  &
        -k(317)*n(idx_H3j)  &
        -k(316)*n(idx_H3j)

    !d[H3CO+_dot]/d[E]
    pd(107,1) =  &
        -k(320)*n(idx_H3COj)  &
        -k(322)*n(idx_H3COj)  &
        -k(319)*n(idx_H3COj)  &
        -k(318)*n(idx_H3COj)  &
        -k(321)*n(idx_H3COj)

    !d[H3O+_dot]/d[E]
    pd(108,1) =  &
        -k(324)*n(idx_H3Oj)  &
        -k(323)*n(idx_H3Oj)  &
        -k(326)*n(idx_H3Oj)  &
        -k(325)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[E]
    pd(109,1) =  &
        -k(329)*n(idx_HCNHj)  &
        -k(328)*n(idx_HCNHj)  &
        -k(330)*n(idx_HCNHj)

    !d[HCO2+_dot]/d[E]
    pd(110,1) =  &
        -k(332)*n(idx_HCO2j)  &
        -k(333)*n(idx_HCO2j)  &
        -k(334)*n(idx_HCO2j)

    !d[HEH+_dot]/d[E]
    pd(111,1) =  &
        -k(337)*n(idx_HEHj)

    !d[N2H+_dot]/d[E]
    pd(112,1) =  &
        -k(340)*n(idx_N2Hj)  &
        -k(339)*n(idx_N2Hj)

    !d[O2H+_dot]/d[E]
    pd(113,1) =  &
        -k(348)*n(idx_O2Hj)

    !d[SIH5+_dot]/d[E]
    pd(114,1) =  &
        -k(362)*n(idx_SIH5j)  &
        -k(361)*n(idx_SIH5j)

    !d[SIOH+_dot]/d[E]
    pd(115,1) =  &
        -k(364)*n(idx_SIOHj)  &
        -k(365)*n(idx_SIOHj)

    !d[E_dot]/d[CH]
    pd(1,2) =  &
        +k(1)*n(idx_O)  &
        +k(1130)

    !d[CH_dot]/d[CH]
    pd(2,2) =  &
        -k(937)*n(idx_O2)  &
        -k(931)*n(idx_NO)  &
        -k(926)*n(idx_HCO)  &
        -k(477)*n(idx_SIj)  &
        -k(54)*n(idx_CNj)  &
        -k(1129)  &
        -k(467)*n(idx_HCOj)  &
        -k(79)*n(idx_Hj)  &
        -k(463)*n(idx_H3Oj)  &
        -k(103)*n(idx_H2j)  &
        -k(932)*n(idx_NO)  &
        -k(56)*n(idx_H2COj)  &
        -k(959)*n(idx_H2)  &
        -k(471)*n(idx_NHj)  &
        -k(925)*n(idx_H2CO)  &
        -k(465)*n(idx_HCNHj)  &
        -k(941)*n(idx_O)  &
        -k(469)*n(idx_Nj)  &
        -k(935)*n(idx_O2)  &
        -k(1130)  &
        -k(251)  &
        -k(464)*n(idx_HCNj)  &
        -k(60)*n(idx_NH2j)  &
        -k(942)*n(idx_OH)  &
        -k(461)*n(idx_H2Oj)  &
        -k(460)*n(idx_H2COj)  &
        -k(466)*n(idx_HCNHj)  &
        -k(142)*n(idx_HEj)  &
        -k(939)*n(idx_O2H)  &
        -k(1254)  &
        -k(58)*n(idx_Nj)  &
        -k(933)*n(idx_NO)  &
        -k(3)*n(idx_H2)  &
        -k(55)*n(idx_COj)  &
        -k(62)*n(idx_O2j)  &
        -k(473)*n(idx_Oj)  &
        -k(16)*n(idx_Cj)  &
        -k(513)*n(idx_H2j)  &
        -k(475)*n(idx_O2Hj)  &
        -k(663)*n(idx_HEj)  &
        -k(1202)*n(idx_H2)  &
        -k(59)*n(idx_N2j)  &
        -k(472)*n(idx_NH2j)  &
        -k(940)*n(idx_O)  &
        -k(459)*n(idx_COj)  &
        -k(930)*n(idx_N)  &
        -k(462)*n(idx_H3COj)  &
        -k(468)*n(idx_HNOj)  &
        -k(938)*n(idx_O2H)  &
        -k(936)*n(idx_O2)  &
        -k(479)*n(idx_SIOj)  &
        -k(927)*n(idx_HNO)  &
        -k(1)*n(idx_O)  &
        -k(57)*n(idx_H2Oj)  &
        -k(476)*n(idx_OHj)  &
        -k(478)*n(idx_SIHj)  &
        -k(929)*n(idx_N)  &
        -k(61)*n(idx_Oj)  &
        -k(581)*n(idx_H3j)  &
        -k(10)*n(idx_H)  &
        -k(971)*n(idx_H)  &
        -k(928)*n(idx_N2)  &
        -k(474)*n(idx_O2j)  &
        -k(924)*n(idx_CO2)  &
        -k(934)*n(idx_O2)  &
        -k(470)*n(idx_N2Hj)  &
        -k(63)*n(idx_OHj)

    !d[O_dot]/d[CH]
    pd(3,2) =  &
        +k(937)*n(idx_O2)  &
        -k(941)*n(idx_O)  &
        -k(1)*n(idx_O)  &
        +k(476)*n(idx_OHj)  &
        +k(474)*n(idx_O2j)  &
        +k(931)*n(idx_NO)  &
        -k(940)*n(idx_O)  &
        +k(935)*n(idx_O2)  &
        +k(61)*n(idx_Oj)

    !d[HNC_dot]/d[CH]
    pd(4,2) =  &
        +k(466)*n(idx_HCNHj)

    !d[HCN_dot]/d[CH]
    pd(5,2) =  &
        +k(465)*n(idx_HCNHj)  &
        +k(931)*n(idx_NO)  &
        +k(928)*n(idx_N2)

    !d[H2_dot]/d[CH]
    pd(6,2) =  &
        -k(3)*n(idx_H2)  &
        +k(3)*n(idx_H2)  &
        +k(581)*n(idx_H3j)  &
        -k(1202)*n(idx_H2)  &
        +k(103)*n(idx_H2j)  &
        +k(971)*n(idx_H)  &
        -k(959)*n(idx_H2)

    !d[C_dot]/d[CH]
    pd(7,2) =  &
        +k(10)*n(idx_H)  &
        +k(941)*n(idx_O)  &
        +k(930)*n(idx_N)  &
        +k(459)*n(idx_COj)  &
        +k(16)*n(idx_Cj)  &
        +k(3)*n(idx_H2)  &
        +k(971)*n(idx_H)  &
        +k(251)  &
        +k(1129)

    !d[H_dot]/d[CH]
    pd(8,2) =  &
        +k(473)*n(idx_Oj)  &
        +k(929)*n(idx_N)  &
        +2.d0*k(10)*n(idx_H)  &
        +k(942)*n(idx_OH)  &
        +k(1129)  &
        +k(477)*n(idx_SIj)  &
        +k(3)*n(idx_H2)  &
        +k(935)*n(idx_O2)  &
        +k(933)*n(idx_NO)  &
        -k(10)*n(idx_H)  &
        +k(513)*n(idx_H2j)  &
        +k(251)  &
        -k(971)*n(idx_H)  &
        +k(79)*n(idx_Hj)  &
        +k(959)*n(idx_H2)  &
        +k(469)*n(idx_Nj)  &
        +k(940)*n(idx_O)  &
        +k(663)*n(idx_HEj)  &
        +k(934)*n(idx_O2)

    !d[H2O_dot]/d[CH]
    pd(9,2) =  &
        +k(463)*n(idx_H3Oj)  &
        +k(57)*n(idx_H2Oj)

    !d[OH_dot]/d[CH]
    pd(10,2) =  &
        +k(461)*n(idx_H2Oj)  &
        +k(941)*n(idx_O)  &
        +k(936)*n(idx_O2)  &
        -k(942)*n(idx_OH)  &
        +k(63)*n(idx_OHj)  &
        +k(938)*n(idx_O2H)

    !d[O2_dot]/d[CH]
    pd(11,2) =  &
        -k(937)*n(idx_O2)  &
        -k(936)*n(idx_O2)  &
        +k(62)*n(idx_O2j)  &
        -k(935)*n(idx_O2)  &
        +k(939)*n(idx_O2H)  &
        +k(475)*n(idx_O2Hj)  &
        -k(934)*n(idx_O2)

    !d[CH2_dot]/d[CH]
    pd(12,2) =  &
        +k(939)*n(idx_O2H)  &
        +k(959)*n(idx_H2)  &
        +k(926)*n(idx_HCO)  &
        +k(927)*n(idx_HNO)  &
        +k(925)*n(idx_H2CO)

    !d[H2CO_dot]/d[CH]
    pd(13,2) =  &
        +k(462)*n(idx_H3COj)  &
        +k(56)*n(idx_H2COj)  &
        -k(925)*n(idx_H2CO)

    !d[HCO_dot]/d[CH]
    pd(14,2) =  &
        +k(937)*n(idx_O2)  &
        +k(942)*n(idx_OH)  &
        +k(925)*n(idx_H2CO)  &
        +k(460)*n(idx_H2COj)  &
        +k(924)*n(idx_CO2)  &
        -k(926)*n(idx_HCO)  &
        +k(938)*n(idx_O2H)  &
        +k(932)*n(idx_NO)

    !d[NO_dot]/d[CH]
    pd(17,2) =  &
        -k(932)*n(idx_NO)  &
        -k(931)*n(idx_NO)  &
        -k(933)*n(idx_NO)  &
        +k(468)*n(idx_HNOj)  &
        +k(927)*n(idx_HNO)

    !d[SI_dot]/d[CH]
    pd(18,2) =  &
        +k(478)*n(idx_SIHj)  &
        +k(479)*n(idx_SIOj)

    !d[CN_dot]/d[CH]
    pd(24,2) =  &
        +k(929)*n(idx_N)  &
        +k(54)*n(idx_CNj)  &
        +k(464)*n(idx_HCNj)

    !d[CO_dot]/d[CH]
    pd(25,2) =  &
        +k(55)*n(idx_COj)  &
        +k(924)*n(idx_CO2)  &
        +k(936)*n(idx_O2)  &
        +k(926)*n(idx_HCO)  &
        +k(940)*n(idx_O)  &
        +k(467)*n(idx_HCOj)  &
        +k(935)*n(idx_O2)

    !d[N2_dot]/d[CH]
    pd(26,2) =  &
        +k(470)*n(idx_N2Hj)  &
        -k(928)*n(idx_N2)  &
        +k(59)*n(idx_N2j)

    !d[NH2_dot]/d[CH]
    pd(27,2) =  &
        +k(60)*n(idx_NH2j)

    !d[CH3_dot]/d[CH]
    pd(28,2) =  &
        +k(1202)*n(idx_H2)

    !d[N_dot]/d[CH]
    pd(30,2) =  &
        +k(58)*n(idx_Nj)  &
        +k(471)*n(idx_NHj)  &
        +k(928)*n(idx_N2)  &
        -k(930)*n(idx_N)  &
        -k(929)*n(idx_N)  &
        +k(932)*n(idx_NO)

    !d[NH_dot]/d[CH]
    pd(31,2) =  &
        +k(472)*n(idx_NH2j)  &
        +k(930)*n(idx_N)

    !d[HE_dot]/d[CH]
    pd(35,2) =  &
        +k(142)*n(idx_HEj)  &
        +k(663)*n(idx_HEj)

    !d[HNO_dot]/d[CH]
    pd(36,2) =  &
        -k(927)*n(idx_HNO)

    !d[CO2_dot]/d[CH]
    pd(38,2) =  &
        -k(924)*n(idx_CO2)  &
        +k(934)*n(idx_O2)

    !d[O2H_dot]/d[CH]
    pd(43,2) =  &
        -k(938)*n(idx_O2H)  &
        -k(939)*n(idx_O2H)

    !d[OCN_dot]/d[CH]
    pd(44,2) =  &
        +k(933)*n(idx_NO)

    !d[CH4_DUST_dot]/d[CH]
    pd(53,2) =  &
        +k(1254)

    !d[HCO+_dot]/d[CH]
    pd(70,2) =  &
        +k(459)*n(idx_COj)  &
        +k(479)*n(idx_SIOj)  &
        +k(1)*n(idx_O)  &
        -k(467)*n(idx_HCOj)  &
        +k(474)*n(idx_O2j)

    !d[H+_dot]/d[CH]
    pd(71,2) =  &
        -k(79)*n(idx_Hj)

    !d[C+_dot]/d[CH]
    pd(73,2) =  &
        +k(663)*n(idx_HEj)  &
        -k(16)*n(idx_Cj)

    !d[CH2+_dot]/d[CH]
    pd(74,2) =  &
        +k(470)*n(idx_N2Hj)  &
        +k(472)*n(idx_NH2j)  &
        +k(468)*n(idx_HNOj)  &
        +k(471)*n(idx_NHj)  &
        +k(460)*n(idx_H2COj)  &
        +k(476)*n(idx_OHj)  &
        +k(462)*n(idx_H3COj)  &
        +k(464)*n(idx_HCNj)  &
        +k(581)*n(idx_H3j)  &
        +k(465)*n(idx_HCNHj)  &
        +k(475)*n(idx_O2Hj)  &
        +k(513)*n(idx_H2j)  &
        +k(463)*n(idx_H3Oj)  &
        +k(478)*n(idx_SIHj)  &
        +k(467)*n(idx_HCOj)  &
        +k(461)*n(idx_H2Oj)  &
        +k(466)*n(idx_HCNHj)

    !d[CH+_dot]/d[CH]
    pd(75,2) =  &
        +k(58)*n(idx_Nj)  &
        +k(56)*n(idx_H2COj)  &
        +k(59)*n(idx_N2j)  &
        +k(57)*n(idx_H2Oj)  &
        +k(62)*n(idx_O2j)  &
        +k(55)*n(idx_COj)  &
        +k(1130)  &
        +k(103)*n(idx_H2j)  &
        +k(142)*n(idx_HEj)  &
        +k(54)*n(idx_CNj)  &
        +k(16)*n(idx_Cj)  &
        +k(79)*n(idx_Hj)  &
        +k(60)*n(idx_NH2j)  &
        +k(63)*n(idx_OHj)  &
        +k(61)*n(idx_Oj)

    !d[H2CO+_dot]/d[CH]
    pd(76,2) =  &
        -k(460)*n(idx_H2COj)  &
        -k(56)*n(idx_H2COj)

    !d[SI+_dot]/d[CH]
    pd(80,2) =  &
        -k(477)*n(idx_SIj)

    !d[SIC+_dot]/d[CH]
    pd(83,2) =  &
        +k(477)*n(idx_SIj)

    !d[CN+_dot]/d[CH]
    pd(86,2) =  &
        -k(54)*n(idx_CNj)  &
        +k(469)*n(idx_Nj)

    !d[CO+_dot]/d[CH]
    pd(87,2) =  &
        +k(473)*n(idx_Oj)  &
        -k(55)*n(idx_COj)  &
        -k(459)*n(idx_COj)

    !d[N2+_dot]/d[CH]
    pd(88,2) =  &
        -k(59)*n(idx_N2j)

    !d[O2+_dot]/d[CH]
    pd(89,2) =  &
        -k(474)*n(idx_O2j)  &
        -k(62)*n(idx_O2j)

    !d[H2O+_dot]/d[CH]
    pd(90,2) =  &
        -k(57)*n(idx_H2Oj)  &
        -k(461)*n(idx_H2Oj)

    !d[NH2+_dot]/d[CH]
    pd(91,2) =  &
        -k(472)*n(idx_NH2j)  &
        -k(60)*n(idx_NH2j)

    !d[O+_dot]/d[CH]
    pd(92,2) =  &
        -k(61)*n(idx_Oj)  &
        -k(473)*n(idx_Oj)

    !d[OH+_dot]/d[CH]
    pd(93,2) =  &
        -k(476)*n(idx_OHj)  &
        -k(63)*n(idx_OHj)

    !d[N+_dot]/d[CH]
    pd(96,2) =  &
        -k(469)*n(idx_Nj)  &
        -k(58)*n(idx_Nj)

    !d[HCN+_dot]/d[CH]
    pd(97,2) =  &
        -k(464)*n(idx_HCNj)

    !d[NH+_dot]/d[CH]
    pd(98,2) =  &
        -k(471)*n(idx_NHj)

    !d[SIH+_dot]/d[CH]
    pd(100,2) =  &
        -k(478)*n(idx_SIHj)

    !d[SIO+_dot]/d[CH]
    pd(101,2) =  &
        -k(479)*n(idx_SIOj)

    !d[H2+_dot]/d[CH]
    pd(102,2) =  &
        -k(513)*n(idx_H2j)  &
        -k(103)*n(idx_H2j)

    !d[HE+_dot]/d[CH]
    pd(103,2) =  &
        -k(663)*n(idx_HEj)  &
        -k(142)*n(idx_HEj)

    !d[HNO+_dot]/d[CH]
    pd(104,2) =  &
        -k(468)*n(idx_HNOj)

    !d[H3+_dot]/d[CH]
    pd(106,2) =  &
        -k(581)*n(idx_H3j)

    !d[H3CO+_dot]/d[CH]
    pd(107,2) =  &
        -k(462)*n(idx_H3COj)

    !d[H3O+_dot]/d[CH]
    pd(108,2) =  &
        -k(463)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[CH]
    pd(109,2) =  &
        -k(466)*n(idx_HCNHj)  &
        -k(465)*n(idx_HCNHj)

    !d[N2H+_dot]/d[CH]
    pd(112,2) =  &
        -k(470)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CH]
    pd(113,2) =  &
        -k(475)*n(idx_O2Hj)

    !d[E_dot]/d[O]
    pd(1,3) =  &
        +k(240)  &
        +k(1)*n(idx_CH)  &
        +k(283)

    !d[CH_dot]/d[O]
    pd(2,3) =  &
        +k(898)*n(idx_CH2)  &
        -k(940)*n(idx_CH)  &
        -k(941)*n(idx_CH)  &
        -k(1)*n(idx_CH)

    !d[O_dot]/d[O]
    pd(3,3) =  &
        -k(1214)*n(idx_SI)  &
        -k(1197)*n(idx_C)  &
        -k(824)*n(idx_H2Oj)  &
        -k(1075)*n(idx_NH3)  &
        -k(834)*n(idx_SIH2j)  &
        -k(1057)*n(idx_CH4)  &
        -k(1071)*n(idx_HNO)  &
        -k(90)*n(idx_Hj)  &
        -k(1213)*n(idx_SIj)  &
        -k(444)*n(idx_CH3j)  &
        -k(826)*n(idx_N2j)  &
        -k(830)*n(idx_O2Hj)  &
        -k(1088)*n(idx_SIH3)  &
        -k(1080)*n(idx_OCN)  &
        -k(1063)*n(idx_H2O)  &
        -k(218)*n(idx_COj)  &
        -k(1087)*n(idx_SIH2)  &
        -k(1081)*n(idx_OH)  &
        -k(1065)*n(idx_HCN)  &
        -k(283)  &
        -k(1084)*n(idx_SIC)  &
        -k(1059)*n(idx_CN)  &
        -k(1060)*n(idx_CO2)  &
        -k(415)*n(idx_CHj)  &
        -k(825)*n(idx_HCO2j)  &
        -k(1086)*n(idx_SIH2)  &
        -k(966)*n(idx_H2)  &
        -k(1048)*n(idx_NH)  &
        -k(895)*n(idx_CH2)  &
        -k(445)*n(idx_CH3j)  &
        -k(836)*n(idx_SIOj)  &
        -k(1073)*n(idx_NH2)  &
        -k(1076)*n(idx_NO2)  &
        -k(600)*n(idx_H3j)  &
        -k(1078)*n(idx_O2H)  &
        -k(829)*n(idx_NH3j)  &
        -k(941)*n(idx_CH)  &
        -k(897)*n(idx_CH2)  &
        -k(1292)  &
        -k(1079)*n(idx_OCN)  &
        -k(832)*n(idx_SICj)  &
        -k(1089)*n(idx_SIH4)  &
        -k(1068)*n(idx_HCO)  &
        -k(940)*n(idx_CH)  &
        -k(1194)*n(idx_Cj)  &
        -k(1067)*n(idx_HCO)  &
        -k(217)*n(idx_CNj)  &
        -k(1208)*n(idx_H)  &
        -k(1058)*n(idx_CN)  &
        -k(828)*n(idx_NH2j)  &
        -k(768)*n(idx_NHj)  &
        -k(219)*n(idx_N2j)  &
        -k(1062)*n(idx_H2CO)  &
        -k(1072)*n(idx_N2)  &
        -k(823)*n(idx_CH4j)  &
        -k(833)*n(idx_SIHj)  &
        -k(896)*n(idx_CH2)  &
        -k(1069)*n(idx_HNO)  &
        -k(1070)*n(idx_HNO)  &
        -k(1074)*n(idx_NH2)  &
        -k(827)*n(idx_N2Hj)  &
        -k(1)*n(idx_CH)  &
        -k(916)*n(idx_CH3)  &
        -k(1082)*n(idx_SIC2)  &
        -k(835)*n(idx_SIH3j)  &
        -4.d0*k(1212)*n(idx_O)  &
        -k(527)*n(idx_H2j)  &
        -k(1061)*n(idx_H2CN)  &
        -k(1047)*n(idx_NH)  &
        -k(831)*n(idx_OHj)  &
        -k(1066)*n(idx_HCN)  &
        -k(898)*n(idx_CH2)  &
        -k(1083)*n(idx_SIC3)  &
        -k(1077)*n(idx_NO)  &
        -k(1064)*n(idx_HCN)  &
        -k(240)  &
        -k(1085)*n(idx_SIC)  &
        -k(422)*n(idx_CH2j)  &
        -k(1090)*n(idx_SIH)  &
        -k(599)*n(idx_H3j)  &
        -k(917)*n(idx_CH3)

    !d[HCN_dot]/d[O]
    pd(5,3) =  &
        -k(1065)*n(idx_HCN)  &
        -k(1064)*n(idx_HCN)  &
        -k(1066)*n(idx_HCN)

    !d[H2_dot]/d[O]
    pd(6,3) =  &
        -k(966)*n(idx_H2)  &
        +k(916)*n(idx_CH3)  &
        +k(829)*n(idx_NH3j)  &
        +k(600)*n(idx_H3j)  &
        +k(1061)*n(idx_H2CN)  &
        +k(445)*n(idx_CH3j)  &
        +k(895)*n(idx_CH2)  &
        +k(824)*n(idx_H2Oj)  &
        +k(835)*n(idx_SIH3j)  &
        +k(1086)*n(idx_SIH2)

    !d[C_dot]/d[O]
    pd(7,3) =  &
        +k(941)*n(idx_CH)  &
        +k(1085)*n(idx_SIC)  &
        -k(1197)*n(idx_C)  &
        +k(832)*n(idx_SICj)  &
        +k(1059)*n(idx_CN)

    !d[H_dot]/d[O]
    pd(8,3) =  &
        +k(1047)*n(idx_NH)  &
        +k(1090)*n(idx_SIH)  &
        +k(1073)*n(idx_NH2)  &
        +k(831)*n(idx_OHj)  &
        +k(415)*n(idx_CHj)  &
        +k(917)*n(idx_CH3)  &
        +k(1069)*n(idx_HNO)  &
        +k(940)*n(idx_CH)  &
        +2.d0*k(896)*n(idx_CH2)  &
        +k(834)*n(idx_SIH2j)  &
        +k(1081)*n(idx_OH)  &
        +k(599)*n(idx_H3j)  &
        +k(444)*n(idx_CH3j)  &
        +k(897)*n(idx_CH2)  &
        +k(1067)*n(idx_HCO)  &
        +k(527)*n(idx_H2j)  &
        +k(828)*n(idx_NH2j)  &
        +k(422)*n(idx_CH2j)  &
        +k(1066)*n(idx_HCN)  &
        +k(90)*n(idx_Hj)  &
        +k(916)*n(idx_CH3)  &
        +k(966)*n(idx_H2)  &
        +2.d0*k(1087)*n(idx_SIH2)  &
        +k(833)*n(idx_SIHj)  &
        -k(1208)*n(idx_H)  &
        +k(1088)*n(idx_SIH3)

    !d[H2O_dot]/d[O]
    pd(9,3) =  &
        -k(1063)*n(idx_H2O)

    !d[OH_dot]/d[O]
    pd(10,3) =  &
        +k(1078)*n(idx_O2H)  &
        +k(1062)*n(idx_H2CO)  &
        +2.d0*k(1063)*n(idx_H2O)  &
        +k(1070)*n(idx_HNO)  &
        +k(823)*n(idx_CH4j)  &
        -k(1081)*n(idx_OH)  &
        +k(966)*n(idx_H2)  &
        +k(1208)*n(idx_H)  &
        +k(1089)*n(idx_SIH4)  &
        +k(898)*n(idx_CH2)  &
        +k(1075)*n(idx_NH3)  &
        +k(1068)*n(idx_HCO)  &
        +k(1057)*n(idx_CH4)  &
        +k(1048)*n(idx_NH)  &
        +k(941)*n(idx_CH)  &
        +k(1074)*n(idx_NH2)  &
        +k(1064)*n(idx_HCN)

    !d[O2_dot]/d[O]
    pd(11,3) =  &
        +k(1078)*n(idx_O2H)  &
        +k(1076)*n(idx_NO2)  &
        +k(825)*n(idx_HCO2j)  &
        +k(836)*n(idx_SIOj)  &
        +k(1080)*n(idx_OCN)  &
        +k(1060)*n(idx_CO2)  &
        +k(1081)*n(idx_OH)  &
        +2.d0*k(1212)*n(idx_O)  &
        +k(830)*n(idx_O2Hj)  &
        +k(1077)*n(idx_NO)  &
        +k(1071)*n(idx_HNO)

    !d[CH2_dot]/d[O]
    pd(12,3) =  &
        -k(895)*n(idx_CH2)  &
        -k(898)*n(idx_CH2)  &
        -k(896)*n(idx_CH2)  &
        -k(897)*n(idx_CH2)

    !d[H2CO_dot]/d[O]
    pd(13,3) =  &
        -k(1062)*n(idx_H2CO)  &
        +k(917)*n(idx_CH3)

    !d[HCO_dot]/d[O]
    pd(14,3) =  &
        +k(897)*n(idx_CH2)  &
        +k(1062)*n(idx_H2CO)  &
        -k(1068)*n(idx_HCO)  &
        -k(1067)*n(idx_HCO)

    !d[NH3_dot]/d[O]
    pd(16,3) =  &
        -k(1075)*n(idx_NH3)

    !d[NO_dot]/d[O]
    pd(17,3) =  &
        +k(1047)*n(idx_NH)  &
        +k(1059)*n(idx_CN)  &
        +k(1076)*n(idx_NO2)  &
        +k(1070)*n(idx_HNO)  &
        +k(1079)*n(idx_OCN)  &
        +k(1072)*n(idx_N2)  &
        -k(1077)*n(idx_NO)

    !d[SI_dot]/d[O]
    pd(18,3) =  &
        -k(1214)*n(idx_SI)  &
        +k(1084)*n(idx_SIC)

    !d[SIC2_dot]/d[O]
    pd(19,3) =  &
        +k(1083)*n(idx_SIC3)  &
        -k(1082)*n(idx_SIC2)

    !d[SIC3_dot]/d[O]
    pd(20,3) =  &
        -k(1083)*n(idx_SIC3)

    !d[SIC_dot]/d[O]
    pd(21,3) =  &
        -k(1084)*n(idx_SIC)  &
        +k(1082)*n(idx_SIC2)  &
        -k(1085)*n(idx_SIC)

    !d[SIH2_dot]/d[O]
    pd(22,3) =  &
        -k(1086)*n(idx_SIH2)  &
        -k(1087)*n(idx_SIH2)

    !d[SIH3_dot]/d[O]
    pd(23,3) =  &
        -k(1088)*n(idx_SIH3)  &
        +k(1089)*n(idx_SIH4)

    !d[CN_dot]/d[O]
    pd(24,3) =  &
        -k(1059)*n(idx_CN)  &
        +k(1080)*n(idx_OCN)  &
        +k(1064)*n(idx_HCN)  &
        +k(217)*n(idx_CNj)  &
        -k(1058)*n(idx_CN)

    !d[CO_dot]/d[O]
    pd(25,3) =  &
        +k(1083)*n(idx_SIC3)  &
        +k(916)*n(idx_CH3)  &
        +k(1079)*n(idx_OCN)  &
        +k(1084)*n(idx_SIC)  &
        +k(1065)*n(idx_HCN)  &
        +k(940)*n(idx_CH)  &
        +k(896)*n(idx_CH2)  &
        +k(218)*n(idx_COj)  &
        +k(1068)*n(idx_HCO)  &
        +k(1058)*n(idx_CN)  &
        +k(1060)*n(idx_CO2)  &
        +k(895)*n(idx_CH2)  &
        +k(1082)*n(idx_SIC2)  &
        +k(1197)*n(idx_C)

    !d[N2_dot]/d[O]
    pd(26,3) =  &
        +k(219)*n(idx_N2j)  &
        +k(827)*n(idx_N2Hj)  &
        -k(1072)*n(idx_N2)

    !d[NH2_dot]/d[O]
    pd(27,3) =  &
        -k(1073)*n(idx_NH2)  &
        -k(1074)*n(idx_NH2)  &
        +k(1075)*n(idx_NH3)

    !d[CH3_dot]/d[O]
    pd(28,3) =  &
        -k(916)*n(idx_CH3)  &
        -k(917)*n(idx_CH3)  &
        +k(1057)*n(idx_CH4)

    !d[CH4_dot]/d[O]
    pd(29,3) =  &
        -k(1057)*n(idx_CH4)

    !d[N_dot]/d[O]
    pd(30,3) =  &
        +k(768)*n(idx_NHj)  &
        +k(826)*n(idx_N2j)  &
        +k(1072)*n(idx_N2)  &
        +k(1058)*n(idx_CN)  &
        +k(1048)*n(idx_NH)  &
        +k(1077)*n(idx_NO)

    !d[NH_dot]/d[O]
    pd(31,3) =  &
        +k(1074)*n(idx_NH2)  &
        -k(1047)*n(idx_NH)  &
        -k(1048)*n(idx_NH)  &
        +k(1065)*n(idx_HCN)  &
        +k(1071)*n(idx_HNO)

    !d[SIH4_dot]/d[O]
    pd(32,3) =  &
        -k(1089)*n(idx_SIH4)

    !d[SIH_dot]/d[O]
    pd(33,3) =  &
        -k(1090)*n(idx_SIH)

    !d[SIO_dot]/d[O]
    pd(34,3) =  &
        +k(1085)*n(idx_SIC)  &
        +k(1090)*n(idx_SIH)  &
        +k(1087)*n(idx_SIH2)  &
        +k(1214)*n(idx_SI)  &
        +k(1086)*n(idx_SIH2)

    !d[HNO_dot]/d[O]
    pd(36,3) =  &
        -k(1071)*n(idx_HNO)  &
        +k(1073)*n(idx_NH2)  &
        -k(1070)*n(idx_HNO)  &
        -k(1069)*n(idx_HNO)

    !d[CO2_dot]/d[O]
    pd(38,3) =  &
        +k(1067)*n(idx_HCO)  &
        -k(1060)*n(idx_CO2)

    !d[H2CN_dot]/d[O]
    pd(39,3) =  &
        -k(1061)*n(idx_H2CN)

    !d[H2SIO_dot]/d[O]
    pd(40,3) =  &
        +k(1088)*n(idx_SIH3)

    !d[NO2_dot]/d[O]
    pd(42,3) =  &
        -k(1076)*n(idx_NO2)  &
        +k(1069)*n(idx_HNO)

    !d[O2H_dot]/d[O]
    pd(43,3) =  &
        -k(1078)*n(idx_O2H)

    !d[OCN_dot]/d[O]
    pd(44,3) =  &
        -k(1080)*n(idx_OCN)  &
        -k(1079)*n(idx_OCN)  &
        +k(1061)*n(idx_H2CN)  &
        +k(1066)*n(idx_HCN)

    !d[H2O_DUST_dot]/d[O]
    pd(55,3) =  &
        +k(1292)

    !d[HCO+_dot]/d[O]
    pd(70,3) =  &
        +k(445)*n(idx_CH3j)  &
        +k(422)*n(idx_CH2j)  &
        +k(1)*n(idx_CH)  &
        +k(825)*n(idx_HCO2j)

    !d[H+_dot]/d[O]
    pd(71,3) =  &
        -k(90)*n(idx_Hj)

    !d[C+_dot]/d[O]
    pd(73,3) =  &
        -k(1194)*n(idx_Cj)

    !d[CH2+_dot]/d[O]
    pd(74,3) =  &
        -k(422)*n(idx_CH2j)

    !d[CH+_dot]/d[O]
    pd(75,3) =  &
        -k(415)*n(idx_CHj)

    !d[H2CO+_dot]/d[O]
    pd(76,3) =  &
        +k(444)*n(idx_CH3j)

    !d[NH3+_dot]/d[O]
    pd(78,3) =  &
        -k(829)*n(idx_NH3j)

    !d[NO+_dot]/d[O]
    pd(79,3) =  &
        +k(826)*n(idx_N2j)

    !d[SI+_dot]/d[O]
    pd(80,3) =  &
        -k(1213)*n(idx_SIj)  &
        +k(836)*n(idx_SIOj)

    !d[SIC+_dot]/d[O]
    pd(83,3) =  &
        -k(832)*n(idx_SICj)

    !d[SIH2+_dot]/d[O]
    pd(84,3) =  &
        -k(834)*n(idx_SIH2j)

    !d[SIH3+_dot]/d[O]
    pd(85,3) =  &
        -k(835)*n(idx_SIH3j)

    !d[CN+_dot]/d[O]
    pd(86,3) =  &
        -k(217)*n(idx_CNj)

    !d[CO+_dot]/d[O]
    pd(87,3) =  &
        +k(415)*n(idx_CHj)  &
        -k(218)*n(idx_COj)  &
        +k(1194)*n(idx_Cj)

    !d[N2+_dot]/d[O]
    pd(88,3) =  &
        -k(219)*n(idx_N2j)  &
        -k(826)*n(idx_N2j)

    !d[O2+_dot]/d[O]
    pd(89,3) =  &
        +k(831)*n(idx_OHj)  &
        +k(824)*n(idx_H2Oj)

    !d[H2O+_dot]/d[O]
    pd(90,3) =  &
        -k(824)*n(idx_H2Oj)  &
        +k(599)*n(idx_H3j)

    !d[NH2+_dot]/d[O]
    pd(91,3) =  &
        -k(828)*n(idx_NH2j)

    !d[O+_dot]/d[O]
    pd(92,3) =  &
        +k(219)*n(idx_N2j)  &
        +k(90)*n(idx_Hj)  &
        +k(217)*n(idx_CNj)  &
        +k(283)  &
        +k(218)*n(idx_COj)  &
        +k(240)

    !d[OH+_dot]/d[O]
    pd(93,3) =  &
        +k(827)*n(idx_N2Hj)  &
        +k(600)*n(idx_H3j)  &
        -k(831)*n(idx_OHj)  &
        +k(527)*n(idx_H2j)  &
        +k(768)*n(idx_NHj)  &
        +k(830)*n(idx_O2Hj)

    !d[CH3+_dot]/d[O]
    pd(94,3) =  &
        -k(445)*n(idx_CH3j)  &
        +k(823)*n(idx_CH4j)  &
        -k(444)*n(idx_CH3j)

    !d[CH4+_dot]/d[O]
    pd(95,3) =  &
        -k(823)*n(idx_CH4j)

    !d[NH+_dot]/d[O]
    pd(98,3) =  &
        -k(768)*n(idx_NHj)

    !d[SIH+_dot]/d[O]
    pd(100,3) =  &
        -k(833)*n(idx_SIHj)

    !d[SIO+_dot]/d[O]
    pd(101,3) =  &
        -k(836)*n(idx_SIOj)  &
        +k(1213)*n(idx_SIj)  &
        +k(832)*n(idx_SICj)  &
        +k(833)*n(idx_SIHj)

    !d[H2+_dot]/d[O]
    pd(102,3) =  &
        -k(527)*n(idx_H2j)

    !d[HNO+_dot]/d[O]
    pd(104,3) =  &
        +k(828)*n(idx_NH2j)  &
        +k(829)*n(idx_NH3j)

    !d[H3+_dot]/d[O]
    pd(106,3) =  &
        -k(600)*n(idx_H3j)  &
        -k(599)*n(idx_H3j)

    !d[HCO2+_dot]/d[O]
    pd(110,3) =  &
        -k(825)*n(idx_HCO2j)

    !d[N2H+_dot]/d[O]
    pd(112,3) =  &
        -k(827)*n(idx_N2Hj)

    !d[O2H+_dot]/d[O]
    pd(113,3) =  &
        -k(830)*n(idx_O2Hj)

    !d[SIOH+_dot]/d[O]
    pd(115,3) =  &
        +k(834)*n(idx_SIH2j)  &
        +k(835)*n(idx_SIH3j)

    !d[O_dot]/d[HNC]
    pd(3,4) =  &
        +k(845)*n(idx_OHj)

    !d[HNC_dot]/d[HNC]
    pd(4,4) =  &
        -k(652)*n(idx_O2Hj)  &
        -k(649)*n(idx_HCOj)  &
        -k(1306)  &
        -k(684)*n(idx_HEj)  &
        -k(610)*n(idx_H3Oj)  &
        -k(651)*n(idx_N2Hj)  &
        -k(686)*n(idx_HEj)  &
        -k(980)*n(idx_H)  &
        -k(648)*n(idx_H3COj)  &
        -k(560)*n(idx_H2Oj)  &
        -k(263)  &
        -k(845)*n(idx_OHj)  &
        -k(590)*n(idx_H3j)  &
        -k(1152)  &
        -k(776)*n(idx_NH2j)  &
        -k(2)*n(idx_Hj)  &
        -k(647)*n(idx_H2COj)  &
        -k(685)*n(idx_HEj)  &
        -k(650)*n(idx_HNOj)  &
        -k(761)*n(idx_NHj)  &
        -k(627)*n(idx_HCNj)  &
        -k(408)*n(idx_CHj)

    !d[HCN_dot]/d[HNC]
    pd(5,4) =  &
        +k(2)*n(idx_Hj)  &
        +k(980)*n(idx_H)

    !d[H2_dot]/d[HNC]
    pd(6,4) =  &
        +k(590)*n(idx_H3j)

    !d[C_dot]/d[HNC]
    pd(7,4) =  &
        +k(408)*n(idx_CHj)  &
        +k(686)*n(idx_HEj)

    !d[H_dot]/d[HNC]
    pd(8,4) =  &
        +k(263)  &
        +k(980)*n(idx_H)  &
        +k(685)*n(idx_HEj)  &
        -k(980)*n(idx_H)  &
        +k(1152)  &
        +k(684)*n(idx_HEj)

    !d[H2O_dot]/d[HNC]
    pd(9,4) =  &
        +k(610)*n(idx_H3Oj)

    !d[OH_dot]/d[HNC]
    pd(10,4) =  &
        +k(560)*n(idx_H2Oj)

    !d[O2_dot]/d[HNC]
    pd(11,4) =  &
        +k(652)*n(idx_O2Hj)

    !d[H2CO_dot]/d[HNC]
    pd(13,4) =  &
        +k(648)*n(idx_H3COj)

    !d[HCO_dot]/d[HNC]
    pd(14,4) =  &
        +k(647)*n(idx_H2COj)

    !d[NO_dot]/d[HNC]
    pd(17,4) =  &
        +k(650)*n(idx_HNOj)

    !d[CN_dot]/d[HNC]
    pd(24,4) =  &
        +k(627)*n(idx_HCNj)  &
        +k(263)  &
        +k(1152)

    !d[CO_dot]/d[HNC]
    pd(25,4) =  &
        +k(649)*n(idx_HCOj)

    !d[N2_dot]/d[HNC]
    pd(26,4) =  &
        +k(651)*n(idx_N2Hj)

    !d[N_dot]/d[HNC]
    pd(30,4) =  &
        +k(761)*n(idx_NHj)  &
        +k(685)*n(idx_HEj)

    !d[NH_dot]/d[HNC]
    pd(31,4) =  &
        +k(776)*n(idx_NH2j)

    !d[HE_dot]/d[HNC]
    pd(35,4) =  &
        +k(686)*n(idx_HEj)  &
        +k(685)*n(idx_HEj)  &
        +k(684)*n(idx_HEj)

    !d[HNC_DUST_dot]/d[HNC]
    pd(67,4) =  &
        +k(1306)

    !d[HCO+_dot]/d[HNC]
    pd(70,4) =  &
        -k(649)*n(idx_HCOj)

    !d[H+_dot]/d[HNC]
    pd(71,4) =  &
        +k(2)*n(idx_Hj)  &
        -k(2)*n(idx_Hj)

    !d[C+_dot]/d[HNC]
    pd(73,4) =  &
        +k(685)*n(idx_HEj)

    !d[CH+_dot]/d[HNC]
    pd(75,4) =  &
        -k(408)*n(idx_CHj)

    !d[H2CO+_dot]/d[HNC]
    pd(76,4) =  &
        -k(647)*n(idx_H2COj)

    !d[CN+_dot]/d[HNC]
    pd(86,4) =  &
        +k(684)*n(idx_HEj)

    !d[H2O+_dot]/d[HNC]
    pd(90,4) =  &
        -k(560)*n(idx_H2Oj)

    !d[NH2+_dot]/d[HNC]
    pd(91,4) =  &
        -k(776)*n(idx_NH2j)

    !d[OH+_dot]/d[HNC]
    pd(93,4) =  &
        -k(845)*n(idx_OHj)

    !d[HCN+_dot]/d[HNC]
    pd(97,4) =  &
        -k(627)*n(idx_HCNj)

    !d[NH+_dot]/d[HNC]
    pd(98,4) =  &
        -k(761)*n(idx_NHj)  &
        +k(686)*n(idx_HEj)

    !d[HE+_dot]/d[HNC]
    pd(103,4) =  &
        -k(685)*n(idx_HEj)  &
        -k(684)*n(idx_HEj)  &
        -k(686)*n(idx_HEj)

    !d[HNO+_dot]/d[HNC]
    pd(104,4) =  &
        -k(650)*n(idx_HNOj)

    !d[H3+_dot]/d[HNC]
    pd(106,4) =  &
        -k(590)*n(idx_H3j)

    !d[H3CO+_dot]/d[HNC]
    pd(107,4) =  &
        -k(648)*n(idx_H3COj)

    !d[H3O+_dot]/d[HNC]
    pd(108,4) =  &
        -k(610)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[HNC]
    pd(109,4) =  &
        +k(648)*n(idx_H3COj)  &
        +k(647)*n(idx_H2COj)  &
        +k(649)*n(idx_HCOj)  &
        +k(627)*n(idx_HCNj)  &
        +k(776)*n(idx_NH2j)  &
        +k(650)*n(idx_HNOj)  &
        +k(761)*n(idx_NHj)  &
        +k(651)*n(idx_N2Hj)  &
        +k(845)*n(idx_OHj)  &
        +k(408)*n(idx_CHj)  &
        +k(652)*n(idx_O2Hj)  &
        +k(590)*n(idx_H3j)  &
        +k(610)*n(idx_H3Oj)  &
        +k(560)*n(idx_H2Oj)

    !d[N2H+_dot]/d[HNC]
    pd(112,4) =  &
        -k(651)*n(idx_N2Hj)

    !d[O2H+_dot]/d[HNC]
    pd(113,4) =  &
        -k(652)*n(idx_O2Hj)

    !d[CH_dot]/d[HCN]
    pd(2,5) =  &
        +k(816)*n(idx_Oj)  &
        +k(678)*n(idx_HEj)

    !d[O_dot]/d[HCN]
    pd(3,5) =  &
        -k(1066)*n(idx_O)  &
        -k(1065)*n(idx_O)  &
        -k(1064)*n(idx_O)  &
        +k(842)*n(idx_OHj)

    !d[HCN_dot]/d[HCN]
    pd(5,5) =  &
        -k(628)*n(idx_H2COj)  &
        -k(1065)*n(idx_O)  &
        -k(816)*n(idx_Oj)  &
        -k(162)*n(idx_Nj)  &
        -k(609)*n(idx_H3Oj)  &
        -k(633)*n(idx_O2Hj)  &
        -k(1267)  &
        -k(815)*n(idx_Oj)  &
        -k(1096)*n(idx_OH)  &
        -k(977)*n(idx_H)  &
        -k(135)*n(idx_COj)  &
        -k(66)*n(idx_CNj)  &
        -k(629)*n(idx_H3COj)  &
        -k(680)*n(idx_HEj)  &
        -k(774)*n(idx_NH2j)  &
        -k(260)  &
        -k(624)*n(idx_HCNj)  &
        -k(136)*n(idx_N2j)  &
        -k(678)*n(idx_HEj)  &
        -k(1095)*n(idx_OH)  &
        -k(677)*n(idx_HEj)  &
        -k(679)*n(idx_HEj)  &
        -k(630)*n(idx_HCOj)  &
        -k(82)*n(idx_Hj)  &
        -k(108)*n(idx_H2j)  &
        -k(406)*n(idx_CHj)  &
        -k(588)*n(idx_H3j)  &
        -k(759)*n(idx_NHj)  &
        -k(631)*n(idx_HNOj)  &
        -k(1064)*n(idx_O)  &
        -k(1148)  &
        -k(1066)*n(idx_O)  &
        -k(557)*n(idx_H2Oj)  &
        -k(632)*n(idx_N2Hj)  &
        -k(842)*n(idx_OHj)

    !d[H2_dot]/d[HCN]
    pd(6,5) =  &
        +k(977)*n(idx_H)  &
        +k(588)*n(idx_H3j)  &
        +k(108)*n(idx_H2j)

    !d[C_dot]/d[HCN]
    pd(7,5) =  &
        +k(406)*n(idx_CHj)

    !d[H_dot]/d[HCN]
    pd(8,5) =  &
        +k(82)*n(idx_Hj)  &
        +k(1066)*n(idx_O)  &
        -k(977)*n(idx_H)  &
        +k(260)  &
        +k(1148)  &
        +k(677)*n(idx_HEj)  &
        +k(679)*n(idx_HEj)

    !d[H2O_dot]/d[HCN]
    pd(9,5) =  &
        +k(1095)*n(idx_OH)  &
        +k(609)*n(idx_H3Oj)

    !d[OH_dot]/d[HCN]
    pd(10,5) =  &
        -k(1095)*n(idx_OH)  &
        -k(1096)*n(idx_OH)  &
        +k(557)*n(idx_H2Oj)  &
        +k(1064)*n(idx_O)

    !d[O2_dot]/d[HCN]
    pd(11,5) =  &
        +k(633)*n(idx_O2Hj)

    !d[H2CO_dot]/d[HCN]
    pd(13,5) =  &
        +k(629)*n(idx_H3COj)

    !d[HCO_dot]/d[HCN]
    pd(14,5) =  &
        +k(628)*n(idx_H2COj)

    !d[NO_dot]/d[HCN]
    pd(17,5) =  &
        +k(631)*n(idx_HNOj)

    !d[CN_dot]/d[HCN]
    pd(24,5) =  &
        +k(624)*n(idx_HCNj)  &
        +k(260)  &
        +k(66)*n(idx_CNj)  &
        +k(977)*n(idx_H)  &
        +k(1148)  &
        +k(1095)*n(idx_OH)  &
        +k(1064)*n(idx_O)

    !d[CO_dot]/d[HCN]
    pd(25,5) =  &
        +k(1065)*n(idx_O)  &
        +k(1096)*n(idx_OH)  &
        +k(135)*n(idx_COj)  &
        +k(630)*n(idx_HCOj)

    !d[N2_dot]/d[HCN]
    pd(26,5) =  &
        +k(632)*n(idx_N2Hj)  &
        +k(136)*n(idx_N2j)

    !d[NH2_dot]/d[HCN]
    pd(27,5) =  &
        +k(1096)*n(idx_OH)

    !d[N_dot]/d[HCN]
    pd(30,5) =  &
        +k(815)*n(idx_Oj)  &
        +k(759)*n(idx_NHj)  &
        +k(680)*n(idx_HEj)  &
        +k(162)*n(idx_Nj)  &
        +k(679)*n(idx_HEj)

    !d[NH_dot]/d[HCN]
    pd(31,5) =  &
        +k(1065)*n(idx_O)  &
        +k(774)*n(idx_NH2j)

    !d[HE_dot]/d[HCN]
    pd(35,5) =  &
        +k(677)*n(idx_HEj)  &
        +k(680)*n(idx_HEj)  &
        +k(678)*n(idx_HEj)  &
        +k(679)*n(idx_HEj)

    !d[OCN_dot]/d[HCN]
    pd(44,5) =  &
        +k(1066)*n(idx_O)

    !d[HCN_DUST_dot]/d[HCN]
    pd(59,5) =  &
        +k(1267)

    !d[HCO+_dot]/d[HCN]
    pd(70,5) =  &
        +k(815)*n(idx_Oj)  &
        -k(630)*n(idx_HCOj)

    !d[H+_dot]/d[HCN]
    pd(71,5) =  &
        -k(82)*n(idx_Hj)

    !d[C+_dot]/d[HCN]
    pd(73,5) =  &
        +k(679)*n(idx_HEj)

    !d[CH+_dot]/d[HCN]
    pd(75,5) =  &
        +k(680)*n(idx_HEj)  &
        -k(406)*n(idx_CHj)

    !d[H2CO+_dot]/d[HCN]
    pd(76,5) =  &
        -k(628)*n(idx_H2COj)

    !d[NO+_dot]/d[HCN]
    pd(79,5) =  &
        +k(816)*n(idx_Oj)

    !d[CN+_dot]/d[HCN]
    pd(86,5) =  &
        +k(677)*n(idx_HEj)  &
        -k(66)*n(idx_CNj)

    !d[CO+_dot]/d[HCN]
    pd(87,5) =  &
        -k(135)*n(idx_COj)

    !d[N2+_dot]/d[HCN]
    pd(88,5) =  &
        -k(136)*n(idx_N2j)

    !d[H2O+_dot]/d[HCN]
    pd(90,5) =  &
        -k(557)*n(idx_H2Oj)

    !d[NH2+_dot]/d[HCN]
    pd(91,5) =  &
        -k(774)*n(idx_NH2j)

    !d[O+_dot]/d[HCN]
    pd(92,5) =  &
        -k(816)*n(idx_Oj)  &
        -k(815)*n(idx_Oj)

    !d[OH+_dot]/d[HCN]
    pd(93,5) =  &
        -k(842)*n(idx_OHj)

    !d[N+_dot]/d[HCN]
    pd(96,5) =  &
        +k(678)*n(idx_HEj)  &
        -k(162)*n(idx_Nj)

    !d[HCN+_dot]/d[HCN]
    pd(97,5) =  &
        +k(82)*n(idx_Hj)  &
        +k(135)*n(idx_COj)  &
        +k(162)*n(idx_Nj)  &
        +k(66)*n(idx_CNj)  &
        +k(108)*n(idx_H2j)  &
        +k(136)*n(idx_N2j)  &
        -k(624)*n(idx_HCNj)

    !d[NH+_dot]/d[HCN]
    pd(98,5) =  &
        -k(759)*n(idx_NHj)

    !d[H2+_dot]/d[HCN]
    pd(102,5) =  &
        -k(108)*n(idx_H2j)

    !d[HE+_dot]/d[HCN]
    pd(103,5) =  &
        -k(677)*n(idx_HEj)  &
        -k(678)*n(idx_HEj)  &
        -k(679)*n(idx_HEj)  &
        -k(680)*n(idx_HEj)

    !d[HNO+_dot]/d[HCN]
    pd(104,5) =  &
        -k(631)*n(idx_HNOj)

    !d[H3+_dot]/d[HCN]
    pd(106,5) =  &
        -k(588)*n(idx_H3j)

    !d[H3CO+_dot]/d[HCN]
    pd(107,5) =  &
        -k(629)*n(idx_H3COj)

    !d[H3O+_dot]/d[HCN]
    pd(108,5) =  &
        -k(609)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[HCN]
    pd(109,5) =  &
        +k(624)*n(idx_HCNj)  &
        +k(557)*n(idx_H2Oj)  &
        +k(588)*n(idx_H3j)  &
        +k(631)*n(idx_HNOj)  &
        +k(628)*n(idx_H2COj)  &
        +k(630)*n(idx_HCOj)  &
        +k(632)*n(idx_N2Hj)  &
        +k(629)*n(idx_H3COj)  &
        +k(842)*n(idx_OHj)  &
        +k(774)*n(idx_NH2j)  &
        +k(759)*n(idx_NHj)  &
        +k(609)*n(idx_H3Oj)  &
        +k(633)*n(idx_O2Hj)  &
        +k(406)*n(idx_CHj)

    !d[N2H+_dot]/d[HCN]
    pd(112,5) =  &
        -k(632)*n(idx_N2Hj)

    !d[O2H+_dot]/d[HCN]
    pd(113,5) =  &
        -k(633)*n(idx_O2Hj)

    !d[E_dot]/d[H2]
    pd(1,6) =  &
        +k(9)*n(idx_E)  &
        +k(234)  &
        -k(9)*n(idx_E)  &
        +k(235)

    !d[CH_dot]/d[H2]
    pd(2,6) =  &
        -k(3)*n(idx_CH)  &
        -k(1202)*n(idx_CH)  &
        +k(956)*n(idx_C)  &
        -k(959)*n(idx_CH)

    !d[O_dot]/d[H2]
    pd(3,6) =  &
        +2.d0*k(7)*n(idx_O2)  &
        -k(966)*n(idx_O)  &
        +k(8)*n(idx_OH)

    !d[HCN_dot]/d[H2]
    pd(5,6) =  &
        +k(960)*n(idx_CN)

    !d[H2_dot]/d[H2]
    pd(6,6) =  &
        -k(1204)*n(idx_SIHj)  &
        -k(9)*n(idx_E)  &
        -k(956)*n(idx_C)  &
        -k(960)*n(idx_CN)  &
        -k(8)*n(idx_OH)  &
        -k(529)*n(idx_Cj)  &
        -k(963)*n(idx_NH)  &
        +k(7)*n(idx_O2)  &
        +k(5)*n(idx_H2O)  &
        -k(964)*n(idx_O2)  &
        -k(537)*n(idx_HEj)  &
        -k(540)*n(idx_N2j)  &
        -k(536)*n(idx_HCNj)  &
        -k(967)*n(idx_OH)  &
        -k(236)  &
        -k(958)*n(idx_CH3)  &
        -k(7)*n(idx_O2)  &
        -k(1202)*n(idx_CH)  &
        -k(539)*n(idx_Nj)  &
        -k(531)*n(idx_CH2j)  &
        -k(533)*n(idx_COj)  &
        -k(546)*n(idx_OHj)  &
        -k(535)*n(idx_H2Oj)  &
        -k(962)*n(idx_NH2)  &
        -k(1200)*n(idx_Cj)  &
        -k(1205)*n(idx_SIH3j)  &
        -k(544)*n(idx_Oj)  &
        -k(532)*n(idx_CNj)  &
        -k(234)  &
        -k(961)*n(idx_N)  &
        -k(543)*n(idx_NH2j)  &
        -k(548)*n(idx_SIOj)  &
        -k(966)*n(idx_O)  &
        -k(534)*n(idx_COj)  &
        -k(1203)*n(idx_SIj)  &
        -k(965)*n(idx_O2)  &
        -k(959)*n(idx_CH)  &
        -k(5)*n(idx_H2O)  &
        -k(547)*n(idx_SIH4j)  &
        +2.d0*k(4)*n(idx_H2)  &
        +k(3)*n(idx_CH)  &
        -k(3)*n(idx_CH)  &
        -k(6)*n(idx_HOCj)  &
        -k(517)*n(idx_H2j)  &
        -k(1201)*n(idx_C)  &
        -4.d0*k(4)*n(idx_H2)  &
        -k(545)*n(idx_O2Hj)  &
        -k(542)*n(idx_NHj)  &
        -k(541)*n(idx_NHj)  &
        -k(957)*n(idx_CH2)  &
        -k(530)*n(idx_CHj)  &
        -k(116)*n(idx_HEj)  &
        +k(8)*n(idx_OH)  &
        -k(538)*n(idx_HEHj)  &
        -k(235)  &
        +k(6)*n(idx_HOCj)  &
        -k(11)*n(idx_H)

    !d[C_dot]/d[H2]
    pd(7,6) =  &
        -k(956)*n(idx_C)  &
        +k(3)*n(idx_CH)  &
        -k(1201)*n(idx_C)

    !d[H_dot]/d[H2]
    pd(8,6) =  &
        +k(967)*n(idx_OH)  &
        +k(959)*n(idx_CH)  &
        +k(542)*n(idx_NHj)  &
        +k(957)*n(idx_CH2)  &
        +k(966)*n(idx_O)  &
        +k(964)*n(idx_O2)  &
        +k(5)*n(idx_H2O)  &
        +k(531)*n(idx_CH2j)  &
        +k(961)*n(idx_N)  &
        +k(543)*n(idx_NH2j)  &
        +k(547)*n(idx_SIH4j)  &
        +k(958)*n(idx_CH3)  &
        +2.d0*k(236)  &
        +k(535)*n(idx_H2Oj)  &
        +k(517)*n(idx_H2j)  &
        +k(234)  &
        +k(534)*n(idx_COj)  &
        +k(530)*n(idx_CHj)  &
        +k(529)*n(idx_Cj)  &
        +k(540)*n(idx_N2j)  &
        +2.d0*k(9)*n(idx_E)  &
        +k(960)*n(idx_CN)  &
        -k(11)*n(idx_H)  &
        +k(544)*n(idx_Oj)  &
        +k(537)*n(idx_HEj)  &
        +k(3)*n(idx_CH)  &
        +k(956)*n(idx_C)  &
        +k(539)*n(idx_Nj)  &
        +k(532)*n(idx_CNj)  &
        +4.d0*k(4)*n(idx_H2)  &
        +k(533)*n(idx_COj)  &
        +k(536)*n(idx_HCNj)  &
        +k(548)*n(idx_SIOj)  &
        +3.d0*k(11)*n(idx_H)  &
        +k(962)*n(idx_NH2)  &
        +k(8)*n(idx_OH)  &
        +k(546)*n(idx_OHj)  &
        +k(963)*n(idx_NH)

    !d[H2O_dot]/d[H2]
    pd(9,6) =  &
        +k(967)*n(idx_OH)  &
        -k(5)*n(idx_H2O)

    !d[OH_dot]/d[H2]
    pd(10,6) =  &
        -k(8)*n(idx_OH)  &
        +k(966)*n(idx_O)  &
        +2.d0*k(965)*n(idx_O2)  &
        -k(967)*n(idx_OH)  &
        +k(5)*n(idx_H2O)

    !d[O2_dot]/d[H2]
    pd(11,6) =  &
        -k(964)*n(idx_O2)  &
        -k(965)*n(idx_O2)  &
        -k(7)*n(idx_O2)  &
        +k(545)*n(idx_O2Hj)

    !d[CH2_dot]/d[H2]
    pd(12,6) =  &
        +k(959)*n(idx_CH)  &
        +k(1201)*n(idx_C)  &
        -k(957)*n(idx_CH2)

    !d[NH3_dot]/d[H2]
    pd(16,6) =  &
        +k(962)*n(idx_NH2)

    !d[CN_dot]/d[H2]
    pd(24,6) =  &
        -k(960)*n(idx_CN)

    !d[NH2_dot]/d[H2]
    pd(27,6) =  &
        -k(962)*n(idx_NH2)  &
        +k(963)*n(idx_NH)

    !d[CH3_dot]/d[H2]
    pd(28,6) =  &
        +k(957)*n(idx_CH2)  &
        -k(958)*n(idx_CH3)  &
        +k(1202)*n(idx_CH)

    !d[CH4_dot]/d[H2]
    pd(29,6) =  &
        +k(958)*n(idx_CH3)

    !d[N_dot]/d[H2]
    pd(30,6) =  &
        +k(541)*n(idx_NHj)  &
        -k(961)*n(idx_N)

    !d[NH_dot]/d[H2]
    pd(31,6) =  &
        -k(963)*n(idx_NH)  &
        +k(961)*n(idx_N)

    !d[HE_dot]/d[H2]
    pd(35,6) =  &
        +k(116)*n(idx_HEj)  &
        +k(537)*n(idx_HEj)  &
        +k(538)*n(idx_HEHj)

    !d[O2H_dot]/d[H2]
    pd(43,6) =  &
        +k(964)*n(idx_O2)

    !d[HCO+_dot]/d[H2]
    pd(70,6) =  &
        +k(6)*n(idx_HOCj)  &
        +k(533)*n(idx_COj)

    !d[H+_dot]/d[H2]
    pd(71,6) =  &
        +k(234)  &
        +k(537)*n(idx_HEj)

    !d[HOC+_dot]/d[H2]
    pd(72,6) =  &
        +k(534)*n(idx_COj)  &
        -k(6)*n(idx_HOCj)

    !d[C+_dot]/d[H2]
    pd(73,6) =  &
        -k(1200)*n(idx_Cj)  &
        -k(529)*n(idx_Cj)

    !d[CH2+_dot]/d[H2]
    pd(74,6) =  &
        +k(530)*n(idx_CHj)  &
        +k(1200)*n(idx_Cj)  &
        -k(531)*n(idx_CH2j)

    !d[CH+_dot]/d[H2]
    pd(75,6) =  &
        +k(529)*n(idx_Cj)  &
        -k(530)*n(idx_CHj)

    !d[NH3+_dot]/d[H2]
    pd(78,6) =  &
        +k(543)*n(idx_NH2j)

    !d[SI+_dot]/d[H2]
    pd(80,6) =  &
        -k(1203)*n(idx_SIj)

    !d[SIH2+_dot]/d[H2]
    pd(84,6) =  &
        +k(1203)*n(idx_SIj)

    !d[SIH3+_dot]/d[H2]
    pd(85,6) =  &
        -k(1205)*n(idx_SIH3j)  &
        +k(1204)*n(idx_SIHj)

    !d[CN+_dot]/d[H2]
    pd(86,6) =  &
        -k(532)*n(idx_CNj)

    !d[CO+_dot]/d[H2]
    pd(87,6) =  &
        -k(534)*n(idx_COj)  &
        -k(533)*n(idx_COj)

    !d[N2+_dot]/d[H2]
    pd(88,6) =  &
        -k(540)*n(idx_N2j)

    !d[H2O+_dot]/d[H2]
    pd(90,6) =  &
        -k(535)*n(idx_H2Oj)  &
        +k(546)*n(idx_OHj)

    !d[NH2+_dot]/d[H2]
    pd(91,6) =  &
        +k(542)*n(idx_NHj)  &
        -k(543)*n(idx_NH2j)

    !d[O+_dot]/d[H2]
    pd(92,6) =  &
        -k(544)*n(idx_Oj)

    !d[OH+_dot]/d[H2]
    pd(93,6) =  &
        +k(544)*n(idx_Oj)  &
        -k(546)*n(idx_OHj)

    !d[CH3+_dot]/d[H2]
    pd(94,6) =  &
        +k(531)*n(idx_CH2j)

    !d[N+_dot]/d[H2]
    pd(96,6) =  &
        -k(539)*n(idx_Nj)

    !d[HCN+_dot]/d[H2]
    pd(97,6) =  &
        +k(532)*n(idx_CNj)  &
        -k(536)*n(idx_HCNj)

    !d[NH+_dot]/d[H2]
    pd(98,6) =  &
        +k(539)*n(idx_Nj)  &
        -k(542)*n(idx_NHj)  &
        -k(541)*n(idx_NHj)

    !d[SIH4+_dot]/d[H2]
    pd(99,6) =  &
        -k(547)*n(idx_SIH4j)

    !d[SIH+_dot]/d[H2]
    pd(100,6) =  &
        -k(1204)*n(idx_SIHj)

    !d[SIO+_dot]/d[H2]
    pd(101,6) =  &
        -k(548)*n(idx_SIOj)

    !d[H2+_dot]/d[H2]
    pd(102,6) =  &
        +k(116)*n(idx_HEj)  &
        +k(235)  &
        -k(517)*n(idx_H2j)

    !d[HE+_dot]/d[H2]
    pd(103,6) =  &
        -k(537)*n(idx_HEj)  &
        -k(116)*n(idx_HEj)

    !d[H3+_dot]/d[H2]
    pd(106,6) =  &
        +k(541)*n(idx_NHj)  &
        +k(538)*n(idx_HEHj)  &
        +k(545)*n(idx_O2Hj)  &
        +k(517)*n(idx_H2j)

    !d[H3O+_dot]/d[H2]
    pd(108,6) =  &
        +k(535)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[H2]
    pd(109,6) =  &
        +k(536)*n(idx_HCNj)

    !d[HEH+_dot]/d[H2]
    pd(111,6) =  &
        -k(538)*n(idx_HEHj)

    !d[N2H+_dot]/d[H2]
    pd(112,6) =  &
        +k(540)*n(idx_N2j)

    !d[O2H+_dot]/d[H2]
    pd(113,6) =  &
        -k(545)*n(idx_O2Hj)

    !d[SIH5+_dot]/d[H2]
    pd(114,6) =  &
        +k(1205)*n(idx_SIH3j)  &
        +k(547)*n(idx_SIH4j)

    !d[SIOH+_dot]/d[H2]
    pd(115,6) =  &
        +k(548)*n(idx_SIOj)

    !d[E_dot]/d[C]
    pd(1,7) =  &
        +k(241)  &
        +k(1108)  &
        +k(232)

    !d[CH_dot]/d[C]
    pd(2,7) =  &
        +k(871)*n(idx_NH)  &
        +k(865)*n(idx_HCO)  &
        +k(956)*n(idx_H2)  &
        +2.d0*k(864)*n(idx_CH2)  &
        +k(1207)*n(idx_H)  &
        +k(877)*n(idx_OH)  &
        +k(869)*n(idx_NH2)

    !d[O_dot]/d[C]
    pd(3,7) =  &
        +k(874)*n(idx_O2)  &
        +k(394)*n(idx_OHj)  &
        -k(1197)*n(idx_O)  &
        +k(392)*n(idx_O2j)  &
        +k(877)*n(idx_OH)  &
        +k(872)*n(idx_NO)

    !d[HNC_dot]/d[C]
    pd(4,7) =  &
        +k(868)*n(idx_NH2)  &
        +k(1005)*n(idx_HNCO)

    !d[HCN_dot]/d[C]
    pd(5,7) =  &
        +k(867)*n(idx_NH2)

    !d[H2_dot]/d[C]
    pd(6,7) =  &
        +k(577)*n(idx_H3j)  &
        -k(1201)*n(idx_H2)  &
        -k(956)*n(idx_H2)  &
        +k(385)*n(idx_H3Oj)

    !d[C_dot]/d[C]
    pd(7,7) =  &
        -k(874)*n(idx_O2)  &
        -k(29)*n(idx_COj)  &
        -k(872)*n(idx_NO)  &
        -k(876)*n(idx_OH)  &
        -k(1251)  &
        -k(390)*n(idx_N2Hj)  &
        -k(384)*n(idx_H2Oj)  &
        -k(391)*n(idx_NHj)  &
        -k(1207)*n(idx_H)  &
        -k(1108)  &
        -k(877)*n(idx_OH)  &
        -k(869)*n(idx_NH2)  &
        -k(232)  &
        -k(140)*n(idx_HEj)  &
        -k(1195)*n(idx_N)  &
        -k(386)*n(idx_HCNj)  &
        -k(395)*n(idx_SIHj)  &
        -k(388)*n(idx_HCO2j)  &
        -k(577)*n(idx_H3j)  &
        -k(387)*n(idx_HCOj)  &
        -k(956)*n(idx_H2)  &
        -k(1201)*n(idx_H2)  &
        -k(871)*n(idx_NH)  &
        -k(389)*n(idx_HNOj)  &
        -k(864)*n(idx_CH2)  &
        -k(1196)*n(idx_Oj)  &
        -k(241)  &
        -k(385)*n(idx_H3Oj)  &
        -k(510)*n(idx_H2j)  &
        -k(1197)*n(idx_O)  &
        -k(1005)*n(idx_HNCO)  &
        -k(393)*n(idx_O2Hj)  &
        -k(867)*n(idx_NH2)  &
        -k(870)*n(idx_NH)  &
        -k(865)*n(idx_HCO)  &
        -k(392)*n(idx_O2j)  &
        -k(873)*n(idx_NO)  &
        -k(878)*n(idx_SIH)  &
        -k(30)*n(idx_N2j)  &
        -k(28)*n(idx_CNj)  &
        -k(875)*n(idx_OCN)  &
        -k(868)*n(idx_NH2)  &
        -k(866)*n(idx_N2)  &
        -k(394)*n(idx_OHj)  &
        -k(396)*n(idx_SIOj)  &
        -k(31)*n(idx_O2j)

    !d[H_dot]/d[C]
    pd(8,7) =  &
        +k(510)*n(idx_H2j)  &
        +k(878)*n(idx_SIH)  &
        +k(876)*n(idx_OH)  &
        +k(867)*n(idx_NH2)  &
        +k(395)*n(idx_SIHj)  &
        +k(868)*n(idx_NH2)  &
        -k(1207)*n(idx_H)  &
        +k(870)*n(idx_NH)  &
        +k(956)*n(idx_H2)

    !d[OH_dot]/d[C]
    pd(10,7) =  &
        +k(384)*n(idx_H2Oj)  &
        -k(877)*n(idx_OH)  &
        -k(876)*n(idx_OH)

    !d[O2_dot]/d[C]
    pd(11,7) =  &
        -k(874)*n(idx_O2)  &
        +k(31)*n(idx_O2j)  &
        +k(393)*n(idx_O2Hj)

    !d[CH2_dot]/d[C]
    pd(12,7) =  &
        -k(864)*n(idx_CH2)  &
        +k(1201)*n(idx_H2)

    !d[HCO_dot]/d[C]
    pd(14,7) =  &
        -k(865)*n(idx_HCO)

    !d[NO_dot]/d[C]
    pd(17,7) =  &
        +k(389)*n(idx_HNOj)  &
        -k(873)*n(idx_NO)  &
        -k(872)*n(idx_NO)

    !d[SIC_dot]/d[C]
    pd(21,7) =  &
        +k(878)*n(idx_SIH)

    !d[CN_dot]/d[C]
    pd(24,7) =  &
        +k(875)*n(idx_OCN)  &
        +k(386)*n(idx_HCNj)  &
        +k(1195)*n(idx_N)  &
        +k(28)*n(idx_CNj)  &
        +k(872)*n(idx_NO)  &
        +k(870)*n(idx_NH)  &
        +k(866)*n(idx_N2)

    !d[CO_dot]/d[C]
    pd(25,7) =  &
        +k(874)*n(idx_O2)  &
        +k(875)*n(idx_OCN)  &
        +k(865)*n(idx_HCO)  &
        +k(876)*n(idx_OH)  &
        +k(396)*n(idx_SIOj)  &
        +k(1005)*n(idx_HNCO)  &
        +k(1197)*n(idx_O)  &
        +k(873)*n(idx_NO)  &
        +k(387)*n(idx_HCOj)  &
        +k(29)*n(idx_COj)

    !d[N2_dot]/d[C]
    pd(26,7) =  &
        -k(866)*n(idx_N2)  &
        +k(30)*n(idx_N2j)  &
        +k(390)*n(idx_N2Hj)

    !d[NH2_dot]/d[C]
    pd(27,7) =  &
        -k(869)*n(idx_NH2)  &
        -k(867)*n(idx_NH2)  &
        -k(868)*n(idx_NH2)

    !d[N_dot]/d[C]
    pd(30,7) =  &
        +k(871)*n(idx_NH)  &
        +k(873)*n(idx_NO)  &
        +k(391)*n(idx_NHj)  &
        -k(1195)*n(idx_N)  &
        +k(866)*n(idx_N2)

    !d[NH_dot]/d[C]
    pd(31,7) =  &
        -k(870)*n(idx_NH)  &
        -k(871)*n(idx_NH)  &
        +k(869)*n(idx_NH2)

    !d[SIH_dot]/d[C]
    pd(33,7) =  &
        -k(878)*n(idx_SIH)

    !d[HE_dot]/d[C]
    pd(35,7) =  &
        +k(140)*n(idx_HEj)

    !d[CO2_dot]/d[C]
    pd(38,7) =  &
        +k(388)*n(idx_HCO2j)

    !d[HNCO_dot]/d[C]
    pd(41,7) =  &
        -k(1005)*n(idx_HNCO)

    !d[OCN_dot]/d[C]
    pd(44,7) =  &
        -k(875)*n(idx_OCN)

    !d[CH4_DUST_dot]/d[C]
    pd(53,7) =  &
        +k(1251)

    !d[HCO+_dot]/d[C]
    pd(70,7) =  &
        -k(387)*n(idx_HCOj)  &
        +k(385)*n(idx_H3Oj)

    !d[C+_dot]/d[C]
    pd(73,7) =  &
        +k(1108)  &
        +k(30)*n(idx_N2j)  &
        +k(31)*n(idx_O2j)  &
        +k(28)*n(idx_CNj)  &
        +k(140)*n(idx_HEj)  &
        +k(241)  &
        +k(29)*n(idx_COj)  &
        +k(232)

    !d[CH+_dot]/d[C]
    pd(75,7) =  &
        +k(510)*n(idx_H2j)  &
        +k(386)*n(idx_HCNj)  &
        +k(394)*n(idx_OHj)  &
        +k(577)*n(idx_H3j)  &
        +k(391)*n(idx_NHj)  &
        +k(384)*n(idx_H2Oj)  &
        +k(388)*n(idx_HCO2j)  &
        +k(393)*n(idx_O2Hj)  &
        +k(389)*n(idx_HNOj)  &
        +k(387)*n(idx_HCOj)  &
        +k(390)*n(idx_N2Hj)

    !d[SI+_dot]/d[C]
    pd(80,7) =  &
        +k(396)*n(idx_SIOj)

    !d[SIC+_dot]/d[C]
    pd(83,7) =  &
        +k(395)*n(idx_SIHj)

    !d[CN+_dot]/d[C]
    pd(86,7) =  &
        -k(28)*n(idx_CNj)

    !d[CO+_dot]/d[C]
    pd(87,7) =  &
        -k(29)*n(idx_COj)  &
        +k(1196)*n(idx_Oj)  &
        +k(392)*n(idx_O2j)

    !d[N2+_dot]/d[C]
    pd(88,7) =  &
        -k(30)*n(idx_N2j)

    !d[O2+_dot]/d[C]
    pd(89,7) =  &
        -k(392)*n(idx_O2j)  &
        -k(31)*n(idx_O2j)

    !d[H2O+_dot]/d[C]
    pd(90,7) =  &
        -k(384)*n(idx_H2Oj)

    !d[O+_dot]/d[C]
    pd(92,7) =  &
        -k(1196)*n(idx_Oj)

    !d[OH+_dot]/d[C]
    pd(93,7) =  &
        -k(394)*n(idx_OHj)

    !d[HCN+_dot]/d[C]
    pd(97,7) =  &
        -k(386)*n(idx_HCNj)

    !d[NH+_dot]/d[C]
    pd(98,7) =  &
        -k(391)*n(idx_NHj)

    !d[SIH+_dot]/d[C]
    pd(100,7) =  &
        -k(395)*n(idx_SIHj)

    !d[SIO+_dot]/d[C]
    pd(101,7) =  &
        -k(396)*n(idx_SIOj)

    !d[H2+_dot]/d[C]
    pd(102,7) =  &
        -k(510)*n(idx_H2j)

    !d[HE+_dot]/d[C]
    pd(103,7) =  &
        -k(140)*n(idx_HEj)

    !d[HNO+_dot]/d[C]
    pd(104,7) =  &
        -k(389)*n(idx_HNOj)

    !d[H3+_dot]/d[C]
    pd(106,7) =  &
        -k(577)*n(idx_H3j)

    !d[H3O+_dot]/d[C]
    pd(108,7) =  &
        -k(385)*n(idx_H3Oj)

    !d[HCO2+_dot]/d[C]
    pd(110,7) =  &
        -k(388)*n(idx_HCO2j)

    !d[N2H+_dot]/d[C]
    pd(112,7) =  &
        -k(390)*n(idx_N2Hj)

    !d[O2H+_dot]/d[C]
    pd(113,7) =  &
        -k(393)*n(idx_O2Hj)

    !d[E_dot]/d[H]
    pd(1,8) =  &
        +k(259)  &
        +k(237)

    !d[CH_dot]/d[H]
    pd(2,8) =  &
        +k(968)*n(idx_CH2)  &
        +k(1207)*n(idx_C)  &
        -k(971)*n(idx_CH)  &
        -k(10)*n(idx_CH)

    !d[O_dot]/d[H]
    pd(3,8) =  &
        +k(132)*n(idx_Oj)  &
        +k(990)*n(idx_O2)  &
        +k(988)*n(idx_NO)  &
        +k(979)*n(idx_HCO)  &
        +k(997)*n(idx_OH)  &
        +k(14)*n(idx_OH)  &
        +2.d0*k(13)*n(idx_O2)  &
        +k(981)*n(idx_HNO)  &
        +k(991)*n(idx_O2H)  &
        -k(1208)*n(idx_O)  &
        +k(994)*n(idx_OCN)

    !d[HNC_dot]/d[H]
    pd(4,8) =  &
        -k(980)*n(idx_HNC)

    !d[HCN_dot]/d[H]
    pd(5,8) =  &
        +k(980)*n(idx_HNC)  &
        +k(994)*n(idx_OCN)  &
        +k(974)*n(idx_H2CN)  &
        +k(130)*n(idx_HCNj)  &
        -k(977)*n(idx_HCN)

    !d[H2_dot]/d[H]
    pd(6,8) =  &
        +k(615)*n(idx_CHj)  &
        +k(971)*n(idx_CH)  &
        +k(129)*n(idx_H2j)  &
        +k(976)*n(idx_H2O)  &
        +k(984)*n(idx_NH2)  &
        +k(997)*n(idx_OH)  &
        +k(975)*n(idx_H2CO)  &
        +k(974)*n(idx_H2CN)  &
        +k(986)*n(idx_NH)  &
        +k(968)*n(idx_CH2)  &
        +k(616)*n(idx_CH2j)  &
        +k(977)*n(idx_HCN)  &
        +k(970)*n(idx_CH4)  &
        +k(620)*n(idx_SIHj)  &
        +k(982)*n(idx_HNO)  &
        +k(978)*n(idx_HCO)  &
        -k(11)*n(idx_H2)  &
        +k(985)*n(idx_NH3)  &
        +k(618)*n(idx_CH4j)  &
        +k(617)*n(idx_CH3j)  &
        +k(992)*n(idx_O2H)  &
        +k(969)*n(idx_CH3)

    !d[C_dot]/d[H]
    pd(7,8) =  &
        +k(973)*n(idx_CO)  &
        -k(1207)*n(idx_C)  &
        +k(971)*n(idx_CH)  &
        +k(10)*n(idx_CH)

    !d[H_dot]/d[H]
    pd(8,8) =  &
        -k(132)*n(idx_Oj)  &
        -k(982)*n(idx_HNO)  &
        -k(976)*n(idx_H2O)  &
        -k(993)*n(idx_O2H)  &
        -k(617)*n(idx_CH3j)  &
        -k(971)*n(idx_CH)  &
        -k(131)*n(idx_HEj)  &
        -k(985)*n(idx_NH3)  &
        -k(12)*n(idx_H2O)  &
        -k(129)*n(idx_H2j)  &
        -k(1208)*n(idx_O)  &
        -k(259)  &
        -k(10)*n(idx_CH)  &
        -k(977)*n(idx_HCN)  &
        +2.d0*k(10)*n(idx_CH)  &
        -k(983)*n(idx_HNO)  &
        -k(969)*n(idx_CH3)  &
        -k(970)*n(idx_CH4)  &
        -k(996)*n(idx_OCN)  &
        -k(974)*n(idx_H2CN)  &
        -k(1207)*n(idx_C)  &
        -k(972)*n(idx_CO2)  &
        -k(13)*n(idx_O2)  &
        -k(1198)*n(idx_Hj)  &
        -k(1206)*n(idx_Cj)  &
        -k(615)*n(idx_CHj)  &
        -k(994)*n(idx_OCN)  &
        -k(981)*n(idx_HNO)  &
        -k(991)*n(idx_O2H)  &
        -k(992)*n(idx_O2H)  &
        -k(619)*n(idx_HEHj)  &
        -k(968)*n(idx_CH2)  &
        -k(616)*n(idx_CH2j)  &
        +k(13)*n(idx_O2)  &
        -k(975)*n(idx_H2CO)  &
        -k(620)*n(idx_SIHj)  &
        -k(979)*n(idx_HCO)  &
        -k(987)*n(idx_NO2)  &
        -k(978)*n(idx_HCO)  &
        +3.d0*k(11)*n(idx_H2)  &
        -k(128)*n(idx_COj)  &
        -k(127)*n(idx_CNj)  &
        -k(1210)*n(idx_SIj)  &
        -k(986)*n(idx_NH)  &
        -k(984)*n(idx_NH2)  &
        -k(237)  &
        -k(995)*n(idx_OCN)  &
        -k(618)*n(idx_CH4j)  &
        -k(14)*n(idx_OH)  &
        -k(988)*n(idx_NO)  &
        -k(11)*n(idx_H2)  &
        -k(973)*n(idx_CO)  &
        -k(130)*n(idx_HCNj)  &
        +k(980)*n(idx_HNC)  &
        -k(980)*n(idx_HNC)  &
        +2.d0*k(14)*n(idx_OH)  &
        -k(1209)*n(idx_OH)  &
        -k(997)*n(idx_OH)  &
        -k(990)*n(idx_O2)  &
        +2.d0*k(12)*n(idx_H2O)  &
        -k(989)*n(idx_NO)

    !d[H2O_dot]/d[H]
    pd(9,8) =  &
        -k(12)*n(idx_H2O)  &
        +k(991)*n(idx_O2H)  &
        +k(1209)*n(idx_OH)  &
        -k(976)*n(idx_H2O)

    !d[OH_dot]/d[H]
    pd(10,8) =  &
        -k(997)*n(idx_OH)  &
        -k(14)*n(idx_OH)  &
        +k(973)*n(idx_CO)  &
        +k(990)*n(idx_O2)  &
        +k(976)*n(idx_H2O)  &
        +k(972)*n(idx_CO2)  &
        -k(1209)*n(idx_OH)  &
        +k(996)*n(idx_OCN)  &
        +2.d0*k(993)*n(idx_O2H)  &
        +k(1208)*n(idx_O)  &
        +k(983)*n(idx_HNO)  &
        +k(12)*n(idx_H2O)  &
        +k(989)*n(idx_NO)  &
        +k(987)*n(idx_NO2)

    !d[O2_dot]/d[H]
    pd(11,8) =  &
        -k(13)*n(idx_O2)  &
        -k(990)*n(idx_O2)  &
        +k(992)*n(idx_O2H)

    !d[CH2_dot]/d[H]
    pd(12,8) =  &
        -k(968)*n(idx_CH2)  &
        +k(969)*n(idx_CH3)  &
        +k(979)*n(idx_HCO)

    !d[H2CO_dot]/d[H]
    pd(13,8) =  &
        -k(975)*n(idx_H2CO)

    !d[HCO_dot]/d[H]
    pd(14,8) =  &
        -k(978)*n(idx_HCO)  &
        +k(975)*n(idx_H2CO)  &
        -k(979)*n(idx_HCO)

    !d[NH3_dot]/d[H]
    pd(16,8) =  &
        -k(985)*n(idx_NH3)

    !d[NO_dot]/d[H]
    pd(17,8) =  &
        -k(988)*n(idx_NO)  &
        +k(982)*n(idx_HNO)  &
        -k(989)*n(idx_NO)  &
        +k(987)*n(idx_NO2)

    !d[CN_dot]/d[H]
    pd(24,8) =  &
        +k(977)*n(idx_HCN)  &
        +k(996)*n(idx_OCN)  &
        +k(127)*n(idx_CNj)

    !d[CO_dot]/d[H]
    pd(25,8) =  &
        +k(978)*n(idx_HCO)  &
        -k(973)*n(idx_CO)  &
        +k(972)*n(idx_CO2)  &
        +k(128)*n(idx_COj)  &
        +k(995)*n(idx_OCN)

    !d[NH2_dot]/d[H]
    pd(27,8) =  &
        +k(981)*n(idx_HNO)  &
        -k(984)*n(idx_NH2)  &
        +k(985)*n(idx_NH3)

    !d[CH3_dot]/d[H]
    pd(28,8) =  &
        +k(970)*n(idx_CH4)  &
        -k(969)*n(idx_CH3)

    !d[CH4_dot]/d[H]
    pd(29,8) =  &
        -k(970)*n(idx_CH4)

    !d[N_dot]/d[H]
    pd(30,8) =  &
        +k(989)*n(idx_NO)  &
        +k(986)*n(idx_NH)

    !d[NH_dot]/d[H]
    pd(31,8) =  &
        +k(988)*n(idx_NO)  &
        +k(984)*n(idx_NH2)  &
        +k(995)*n(idx_OCN)  &
        -k(986)*n(idx_NH)  &
        +k(983)*n(idx_HNO)

    !d[HE_dot]/d[H]
    pd(35,8) =  &
        +k(131)*n(idx_HEj)  &
        +k(619)*n(idx_HEHj)

    !d[HNO_dot]/d[H]
    pd(36,8) =  &
        -k(982)*n(idx_HNO)  &
        -k(981)*n(idx_HNO)  &
        -k(983)*n(idx_HNO)

    !d[CO2_dot]/d[H]
    pd(38,8) =  &
        -k(972)*n(idx_CO2)

    !d[H2CN_dot]/d[H]
    pd(39,8) =  &
        -k(974)*n(idx_H2CN)

    !d[NO2_dot]/d[H]
    pd(42,8) =  &
        -k(987)*n(idx_NO2)

    !d[O2H_dot]/d[H]
    pd(43,8) =  &
        -k(993)*n(idx_O2H)  &
        -k(991)*n(idx_O2H)  &
        -k(992)*n(idx_O2H)

    !d[OCN_dot]/d[H]
    pd(44,8) =  &
        -k(994)*n(idx_OCN)  &
        -k(995)*n(idx_OCN)  &
        -k(996)*n(idx_OCN)

    !d[H+_dot]/d[H]
    pd(71,8) =  &
        +k(132)*n(idx_Oj)  &
        -k(1198)*n(idx_Hj)  &
        +k(129)*n(idx_H2j)  &
        +k(127)*n(idx_CNj)  &
        +k(259)  &
        +k(128)*n(idx_COj)  &
        +k(131)*n(idx_HEj)  &
        +k(130)*n(idx_HCNj)  &
        +k(237)

    !d[C+_dot]/d[H]
    pd(73,8) =  &
        +k(615)*n(idx_CHj)  &
        -k(1206)*n(idx_Cj)

    !d[CH2+_dot]/d[H]
    pd(74,8) =  &
        -k(616)*n(idx_CH2j)  &
        +k(617)*n(idx_CH3j)

    !d[CH+_dot]/d[H]
    pd(75,8) =  &
        +k(1206)*n(idx_Cj)  &
        +k(616)*n(idx_CH2j)  &
        -k(615)*n(idx_CHj)

    !d[SI+_dot]/d[H]
    pd(80,8) =  &
        -k(1210)*n(idx_SIj)  &
        +k(620)*n(idx_SIHj)

    !d[CN+_dot]/d[H]
    pd(86,8) =  &
        -k(127)*n(idx_CNj)

    !d[CO+_dot]/d[H]
    pd(87,8) =  &
        -k(128)*n(idx_COj)

    !d[O+_dot]/d[H]
    pd(92,8) =  &
        -k(132)*n(idx_Oj)

    !d[CH3+_dot]/d[H]
    pd(94,8) =  &
        +k(618)*n(idx_CH4j)  &
        -k(617)*n(idx_CH3j)

    !d[CH4+_dot]/d[H]
    pd(95,8) =  &
        -k(618)*n(idx_CH4j)

    !d[HCN+_dot]/d[H]
    pd(97,8) =  &
        -k(130)*n(idx_HCNj)

    !d[SIH+_dot]/d[H]
    pd(100,8) =  &
        +k(1210)*n(idx_SIj)  &
        -k(620)*n(idx_SIHj)

    !d[H2+_dot]/d[H]
    pd(102,8) =  &
        +k(619)*n(idx_HEHj)  &
        -k(129)*n(idx_H2j)  &
        +k(1198)*n(idx_Hj)

    !d[HE+_dot]/d[H]
    pd(103,8) =  &
        -k(131)*n(idx_HEj)

    !d[HEH+_dot]/d[H]
    pd(111,8) =  &
        -k(619)*n(idx_HEHj)

    !d[E_dot]/d[H2O]
    pd(1,9) =  &
        +k(1142)

    !d[O_dot]/d[H2O]
    pd(3,9) =  &
        +k(841)*n(idx_OHj)  &
        -k(1063)*n(idx_O)  &
        +k(757)*n(idx_NHj)  &
        +k(211)*n(idx_Oj)

    !d[HCN_dot]/d[H2O]
    pd(5,9) =  &
        +k(125)*n(idx_HCNj)

    !d[H2_dot]/d[H2O]
    pd(6,9) =  &
        +k(107)*n(idx_H2j)  &
        +k(587)*n(idx_H3j)  &
        +k(405)*n(idx_CHj)  &
        +k(5)*n(idx_H2)  &
        -k(5)*n(idx_H2)  &
        +k(976)*n(idx_H)  &
        +k(756)*n(idx_NHj)

    !d[C_dot]/d[H2O]
    pd(7,9) =  &
        +k(404)*n(idx_CHj)

    !d[H_dot]/d[H2O]
    pd(8,9) =  &
        +k(1143)  &
        +k(257)  &
        +k(371)*n(idx_Cj)  &
        -k(976)*n(idx_H)  &
        +k(519)*n(idx_H2j)  &
        +k(5)*n(idx_H2)  &
        +k(674)*n(idx_HEj)  &
        +2.d0*k(12)*n(idx_H)  &
        +k(573)*n(idx_SIj)  &
        +k(419)*n(idx_CH2j)  &
        +k(403)*n(idx_CHj)  &
        +k(81)*n(idx_Hj)  &
        -k(12)*n(idx_H)  &
        +k(372)*n(idx_Cj)

    !d[H2O_dot]/d[H2O]
    pd(9,9) =  &
        -k(674)*n(idx_HEj)  &
        -k(905)*n(idx_CH3)  &
        -k(758)*n(idx_NHj)  &
        -k(124)*n(idx_COj)  &
        -k(81)*n(idx_Hj)  &
        -k(5)*n(idx_H2)  &
        -k(565)*n(idx_H3COj)  &
        -k(561)*n(idx_CNj)  &
        -k(755)*n(idx_NHj)  &
        -k(675)*n(idx_HEj)  &
        -k(403)*n(idx_CHj)  &
        -k(404)*n(idx_CHj)  &
        -k(405)*n(idx_CHj)  &
        -k(519)*n(idx_H2j)  &
        -k(570)*n(idx_N2j)  &
        -k(126)*n(idx_N2j)  &
        -k(976)*n(idx_H)  &
        -k(772)*n(idx_NH2j)  &
        -k(756)*n(idx_NHj)  &
        -k(1258)  &
        -k(562)*n(idx_CNj)  &
        -k(563)*n(idx_COj)  &
        -k(573)*n(idx_SIj)  &
        -k(177)*n(idx_NHj)  &
        -k(257)  &
        -k(1143)  &
        -k(1037)*n(idx_NH)  &
        -k(107)*n(idx_H2j)  &
        -k(841)*n(idx_OHj)  &
        -k(125)*n(idx_HCNj)  &
        -k(773)*n(idx_NH2j)  &
        -k(757)*n(idx_NHj)  &
        -k(556)*n(idx_H2Oj)  &
        -k(419)*n(idx_CH2j)  &
        -k(1063)*n(idx_O)  &
        -k(372)*n(idx_Cj)  &
        -k(566)*n(idx_HCNj)  &
        -k(568)*n(idx_HCO2j)  &
        -k(12)*n(idx_H)  &
        -k(221)*n(idx_OHj)  &
        -k(567)*n(idx_HCOj)  &
        -k(569)*n(idx_HNOj)  &
        -k(575)*n(idx_SIH4j)  &
        -k(587)*n(idx_H3j)  &
        -k(161)*n(idx_Nj)  &
        -k(144)*n(idx_HEj)  &
        -k(574)*n(idx_SIHj)  &
        -k(371)*n(idx_Cj)  &
        -k(576)*n(idx_SIH5j)  &
        -k(211)*n(idx_Oj)  &
        -k(571)*n(idx_N2Hj)  &
        -k(451)*n(idx_CH4j)  &
        -k(1142)  &
        -k(564)*n(idx_H2COj)  &
        -k(572)*n(idx_O2Hj)

    !d[OH_dot]/d[H2O]
    pd(10,9) =  &
        +k(556)*n(idx_H2Oj)  &
        +k(570)*n(idx_N2j)  &
        +k(257)  &
        +k(1143)  &
        +k(773)*n(idx_NH2j)  &
        +k(1037)*n(idx_NH)  &
        +k(905)*n(idx_CH3)  &
        +k(675)*n(idx_HEj)  &
        +k(5)*n(idx_H2)  &
        +k(563)*n(idx_COj)  &
        +k(12)*n(idx_H)  &
        +k(221)*n(idx_OHj)  &
        +k(561)*n(idx_CNj)  &
        +k(976)*n(idx_H)  &
        +k(758)*n(idx_NHj)  &
        +2.d0*k(1063)*n(idx_O)

    !d[O2_dot]/d[H2O]
    pd(11,9) =  &
        +k(572)*n(idx_O2Hj)

    !d[H2CO_dot]/d[H2O]
    pd(13,9) =  &
        +k(565)*n(idx_H3COj)

    !d[HCO_dot]/d[H2O]
    pd(14,9) =  &
        +k(564)*n(idx_H2COj)

    !d[NO_dot]/d[H2O]
    pd(17,9) =  &
        +k(569)*n(idx_HNOj)

    !d[SI_dot]/d[H2O]
    pd(18,9) =  &
        +k(574)*n(idx_SIHj)

    !d[SIH3_dot]/d[H2O]
    pd(23,9) =  &
        +k(575)*n(idx_SIH4j)

    !d[CN_dot]/d[H2O]
    pd(24,9) =  &
        +k(566)*n(idx_HCNj)

    !d[CO_dot]/d[H2O]
    pd(25,9) =  &
        +k(124)*n(idx_COj)  &
        +k(567)*n(idx_HCOj)

    !d[N2_dot]/d[H2O]
    pd(26,9) =  &
        +k(571)*n(idx_N2Hj)  &
        +k(126)*n(idx_N2j)

    !d[NH2_dot]/d[H2O]
    pd(27,9) =  &
        +k(1037)*n(idx_NH)

    !d[CH3_dot]/d[H2O]
    pd(28,9) =  &
        +k(451)*n(idx_CH4j)  &
        -k(905)*n(idx_CH3)

    !d[CH4_dot]/d[H2O]
    pd(29,9) =  &
        +k(905)*n(idx_CH3)

    !d[N_dot]/d[H2O]
    pd(30,9) =  &
        +k(161)*n(idx_Nj)  &
        +k(755)*n(idx_NHj)

    !d[NH_dot]/d[H2O]
    pd(31,9) =  &
        +k(177)*n(idx_NHj)  &
        +k(772)*n(idx_NH2j)  &
        +k(562)*n(idx_CNj)  &
        -k(1037)*n(idx_NH)

    !d[SIH4_dot]/d[H2O]
    pd(32,9) =  &
        +k(576)*n(idx_SIH5j)

    !d[HE_dot]/d[H2O]
    pd(35,9) =  &
        +k(674)*n(idx_HEj)  &
        +k(675)*n(idx_HEj)  &
        +k(144)*n(idx_HEj)

    !d[CO2_dot]/d[H2O]
    pd(38,9) =  &
        +k(568)*n(idx_HCO2j)

    !d[H2O_DUST_dot]/d[H2O]
    pd(55,9) =  &
        +k(1258)

    !d[HCO+_dot]/d[H2O]
    pd(70,9) =  &
        +k(405)*n(idx_CHj)  &
        -k(567)*n(idx_HCOj)  &
        +k(563)*n(idx_COj)  &
        +k(371)*n(idx_Cj)  &
        +k(562)*n(idx_CNj)

    !d[H+_dot]/d[H2O]
    pd(71,9) =  &
        +k(675)*n(idx_HEj)  &
        -k(81)*n(idx_Hj)

    !d[HOC+_dot]/d[H2O]
    pd(72,9) =  &
        +k(372)*n(idx_Cj)

    !d[C+_dot]/d[H2O]
    pd(73,9) =  &
        -k(372)*n(idx_Cj)  &
        -k(371)*n(idx_Cj)

    !d[CH2+_dot]/d[H2O]
    pd(74,9) =  &
        -k(419)*n(idx_CH2j)

    !d[CH+_dot]/d[H2O]
    pd(75,9) =  &
        -k(403)*n(idx_CHj)  &
        -k(405)*n(idx_CHj)  &
        -k(404)*n(idx_CHj)

    !d[H2CO+_dot]/d[H2O]
    pd(76,9) =  &
        +k(403)*n(idx_CHj)  &
        -k(564)*n(idx_H2COj)

    !d[NH3+_dot]/d[H2O]
    pd(78,9) =  &
        +k(773)*n(idx_NH2j)  &
        +k(757)*n(idx_NHj)

    !d[SI+_dot]/d[H2O]
    pd(80,9) =  &
        -k(573)*n(idx_SIj)

    !d[CN+_dot]/d[H2O]
    pd(86,9) =  &
        -k(561)*n(idx_CNj)  &
        -k(562)*n(idx_CNj)

    !d[CO+_dot]/d[H2O]
    pd(87,9) =  &
        -k(563)*n(idx_COj)  &
        -k(124)*n(idx_COj)

    !d[N2+_dot]/d[H2O]
    pd(88,9) =  &
        -k(570)*n(idx_N2j)  &
        -k(126)*n(idx_N2j)

    !d[H2O+_dot]/d[H2O]
    pd(90,9) =  &
        +k(107)*n(idx_H2j)  &
        +k(177)*n(idx_NHj)  &
        +k(1142)  &
        +k(144)*n(idx_HEj)  &
        +k(161)*n(idx_Nj)  &
        -k(556)*n(idx_H2Oj)  &
        +k(126)*n(idx_N2j)  &
        +k(221)*n(idx_OHj)  &
        +k(125)*n(idx_HCNj)  &
        +k(81)*n(idx_Hj)  &
        +k(211)*n(idx_Oj)  &
        +k(124)*n(idx_COj)

    !d[NH2+_dot]/d[H2O]
    pd(91,9) =  &
        -k(772)*n(idx_NH2j)  &
        -k(773)*n(idx_NH2j)  &
        +k(758)*n(idx_NHj)

    !d[O+_dot]/d[H2O]
    pd(92,9) =  &
        -k(211)*n(idx_Oj)

    !d[OH+_dot]/d[H2O]
    pd(93,9) =  &
        -k(221)*n(idx_OHj)  &
        +k(674)*n(idx_HEj)  &
        -k(841)*n(idx_OHj)

    !d[CH4+_dot]/d[H2O]
    pd(95,9) =  &
        -k(451)*n(idx_CH4j)

    !d[N+_dot]/d[H2O]
    pd(96,9) =  &
        -k(161)*n(idx_Nj)

    !d[HCN+_dot]/d[H2O]
    pd(97,9) =  &
        -k(125)*n(idx_HCNj)  &
        +k(561)*n(idx_CNj)  &
        -k(566)*n(idx_HCNj)

    !d[NH+_dot]/d[H2O]
    pd(98,9) =  &
        -k(756)*n(idx_NHj)  &
        -k(177)*n(idx_NHj)  &
        -k(758)*n(idx_NHj)  &
        -k(755)*n(idx_NHj)  &
        -k(757)*n(idx_NHj)

    !d[SIH4+_dot]/d[H2O]
    pd(99,9) =  &
        -k(575)*n(idx_SIH4j)

    !d[SIH+_dot]/d[H2O]
    pd(100,9) =  &
        -k(574)*n(idx_SIHj)

    !d[H2+_dot]/d[H2O]
    pd(102,9) =  &
        -k(519)*n(idx_H2j)  &
        -k(107)*n(idx_H2j)

    !d[HE+_dot]/d[H2O]
    pd(103,9) =  &
        -k(674)*n(idx_HEj)  &
        -k(675)*n(idx_HEj)  &
        -k(144)*n(idx_HEj)

    !d[HNO+_dot]/d[H2O]
    pd(104,9) =  &
        +k(756)*n(idx_NHj)  &
        -k(569)*n(idx_HNOj)

    !d[H3+_dot]/d[H2O]
    pd(106,9) =  &
        -k(587)*n(idx_H3j)

    !d[H3CO+_dot]/d[H2O]
    pd(107,9) =  &
        +k(419)*n(idx_CH2j)  &
        -k(565)*n(idx_H3COj)

    !d[H3O+_dot]/d[H2O]
    pd(108,9) =  &
        +k(587)*n(idx_H3j)  &
        +k(556)*n(idx_H2Oj)  &
        +k(575)*n(idx_SIH4j)  &
        +k(567)*n(idx_HCOj)  &
        +k(564)*n(idx_H2COj)  &
        +k(841)*n(idx_OHj)  &
        +k(519)*n(idx_H2j)  &
        +k(451)*n(idx_CH4j)  &
        +k(772)*n(idx_NH2j)  &
        +k(568)*n(idx_HCO2j)  &
        +k(755)*n(idx_NHj)  &
        +k(569)*n(idx_HNOj)  &
        +k(404)*n(idx_CHj)  &
        +k(572)*n(idx_O2Hj)  &
        +k(565)*n(idx_H3COj)  &
        +k(571)*n(idx_N2Hj)  &
        +k(566)*n(idx_HCNj)  &
        +k(576)*n(idx_SIH5j)  &
        +k(574)*n(idx_SIHj)

    !d[HCO2+_dot]/d[H2O]
    pd(110,9) =  &
        -k(568)*n(idx_HCO2j)

    !d[N2H+_dot]/d[H2O]
    pd(112,9) =  &
        +k(570)*n(idx_N2j)  &
        -k(571)*n(idx_N2Hj)

    !d[O2H+_dot]/d[H2O]
    pd(113,9) =  &
        -k(572)*n(idx_O2Hj)

    !d[SIH5+_dot]/d[H2O]
    pd(114,9) =  &
        -k(576)*n(idx_SIH5j)

    !d[SIOH+_dot]/d[H2O]
    pd(115,9) =  &
        +k(573)*n(idx_SIj)

    !d[E_dot]/d[OH]
    pd(1,10) =  &
        +k(1176)

    !d[CH_dot]/d[OH]
    pd(2,10) =  &
        +k(900)*n(idx_CH2)  &
        -k(942)*n(idx_CH)  &
        +k(877)*n(idx_C)

    !d[O_dot]/d[OH]
    pd(3,10) =  &
        +k(848)*n(idx_OHj)  &
        +k(1091)*n(idx_CN)  &
        +k(852)*n(idx_COj)  &
        +2.d0*k(1102)*n(idx_OH)  &
        +k(1033)*n(idx_NH2)  &
        +k(877)*n(idx_C)  &
        +k(285)  &
        +k(918)*n(idx_CH3)  &
        +k(997)*n(idx_H)  &
        +k(901)*n(idx_CH2)  &
        +k(1051)*n(idx_NH)  &
        -k(1081)*n(idx_O)  &
        +k(1175)  &
        +k(853)*n(idx_H2Oj)  &
        +k(216)*n(idx_Oj)  &
        +k(14)*n(idx_H)  &
        +k(8)*n(idx_H2)  &
        +k(1027)*n(idx_N)

    !d[HCN_dot]/d[OH]
    pd(5,10) =  &
        +k(1091)*n(idx_CN)  &
        -k(1095)*n(idx_HCN)  &
        -k(1096)*n(idx_HCN)

    !d[H2_dot]/d[OH]
    pd(6,10) =  &
        +k(919)*n(idx_CH3)  &
        -k(967)*n(idx_H2)  &
        +k(601)*n(idx_H3j)  &
        +k(997)*n(idx_H)  &
        +k(115)*n(idx_H2j)  &
        -k(8)*n(idx_H2)  &
        +k(446)*n(idx_CH3j)  &
        +k(8)*n(idx_H2)  &
        +k(416)*n(idx_CHj)

    !d[C_dot]/d[OH]
    pd(7,10) =  &
        -k(877)*n(idx_C)  &
        -k(876)*n(idx_C)

    !d[H_dot]/d[OH]
    pd(8,10) =  &
        +k(528)*n(idx_H2j)  &
        +k(285)  &
        +k(899)*n(idx_CH2)  &
        -k(14)*n(idx_H)  &
        +k(860)*n(idx_SIj)  &
        +k(1100)*n(idx_NO)  &
        +k(1092)*n(idx_CN)  &
        +k(1081)*n(idx_O)  &
        +k(1175)  &
        +k(876)*n(idx_C)  &
        +k(967)*n(idx_H2)  &
        +k(1103)*n(idx_SI)  &
        +k(380)*n(idx_Cj)  &
        +2.d0*k(14)*n(idx_H)  &
        -k(1209)*n(idx_H)  &
        +k(856)*n(idx_HCOj)  &
        +k(1050)*n(idx_NH)  &
        +k(1093)*n(idx_CO)  &
        +k(8)*n(idx_H2)  &
        +k(1026)*n(idx_N)  &
        +k(820)*n(idx_Oj)  &
        +k(942)*n(idx_CH)  &
        +k(91)*n(idx_Hj)  &
        -k(997)*n(idx_H)  &
        +k(700)*n(idx_HEj)

    !d[H2O_dot]/d[OH]
    pd(9,10) =  &
        +k(920)*n(idx_CH3)  &
        +2.d0*k(1102)*n(idx_OH)  &
        +k(1098)*n(idx_HNO)  &
        +k(1209)*n(idx_H)  &
        +k(1097)*n(idx_HCO)  &
        +k(1099)*n(idx_NH3)  &
        +k(1094)*n(idx_H2CO)  &
        +k(1032)*n(idx_NH2)  &
        +k(967)*n(idx_H2)  &
        +k(900)*n(idx_CH2)  &
        +k(923)*n(idx_CH4)  &
        +k(1049)*n(idx_NH)  &
        +k(1101)*n(idx_O2H)  &
        +k(1095)*n(idx_HCN)

    !d[OH_dot]/d[OH]
    pd(10,10) =  &
        -k(1103)*n(idx_SI)  &
        -k(115)*n(idx_H2j)  &
        -k(856)*n(idx_HCOj)  &
        -k(1093)*n(idx_CO)  &
        -k(528)*n(idx_H2j)  &
        -k(1255)  &
        -k(876)*n(idx_C)  &
        -k(1081)*n(idx_O)  &
        -k(380)*n(idx_Cj)  &
        -k(14)*n(idx_H)  &
        -k(1097)*n(idx_HCO)  &
        -k(1209)*n(idx_H)  &
        -k(228)*n(idx_N2j)  &
        -k(848)*n(idx_OHj)  &
        -k(857)*n(idx_HNOj)  &
        -k(91)*n(idx_Hj)  &
        -k(1092)*n(idx_CN)  &
        -k(967)*n(idx_H2)  &
        -k(942)*n(idx_CH)  &
        -k(1091)*n(idx_CN)  &
        -k(285)  &
        -k(1050)*n(idx_NH)  &
        -k(854)*n(idx_HCNj)  &
        -k(1101)*n(idx_O2H)  &
        -k(852)*n(idx_COj)  &
        -k(1032)*n(idx_NH2)  &
        -k(1094)*n(idx_H2CO)  &
        -k(853)*n(idx_H2Oj)  &
        -k(859)*n(idx_O2Hj)  &
        -k(877)*n(idx_C)  &
        -k(1100)*n(idx_NO)  &
        -k(769)*n(idx_NHj)  &
        -k(1099)*n(idx_NH3)  &
        -k(1176)  &
        -k(226)*n(idx_CNj)  &
        -k(601)*n(idx_H3j)  &
        -k(923)*n(idx_CH4)  &
        -k(900)*n(idx_CH2)  &
        -k(997)*n(idx_H)  &
        -k(1026)*n(idx_N)  &
        -k(899)*n(idx_CH2)  &
        -k(1033)*n(idx_NH2)  &
        -k(8)*n(idx_H2)  &
        -k(918)*n(idx_CH3)  &
        -k(920)*n(idx_CH3)  &
        -k(1051)*n(idx_NH)  &
        -k(700)*n(idx_HEj)  &
        -k(1095)*n(idx_HCN)  &
        -k(227)*n(idx_COj)  &
        -k(416)*n(idx_CHj)  &
        -k(901)*n(idx_CH2)  &
        -k(858)*n(idx_N2Hj)  &
        -k(1096)*n(idx_HCN)  &
        -k(820)*n(idx_Oj)  &
        -k(860)*n(idx_SIj)  &
        -k(170)*n(idx_Nj)  &
        -k(1049)*n(idx_NH)  &
        -k(1175)  &
        -k(1027)*n(idx_N)  &
        -4.d0*k(1102)*n(idx_OH)  &
        -k(216)*n(idx_Oj)  &
        -k(855)*n(idx_HCOj)  &
        -k(1098)*n(idx_HNO)  &
        -k(919)*n(idx_CH3)  &
        -k(446)*n(idx_CH3j)

    !d[O2_dot]/d[OH]
    pd(11,10) =  &
        +k(859)*n(idx_O2Hj)  &
        +k(1101)*n(idx_O2H)  &
        +k(1081)*n(idx_O)

    !d[CH2_dot]/d[OH]
    pd(12,10) =  &
        +k(920)*n(idx_CH3)  &
        -k(900)*n(idx_CH2)  &
        -k(901)*n(idx_CH2)  &
        -k(899)*n(idx_CH2)

    !d[H2CO_dot]/d[OH]
    pd(13,10) =  &
        +k(919)*n(idx_CH3)  &
        +k(899)*n(idx_CH2)  &
        -k(1094)*n(idx_H2CO)

    !d[HCO_dot]/d[OH]
    pd(14,10) =  &
        -k(1097)*n(idx_HCO)  &
        +k(942)*n(idx_CH)  &
        +k(1094)*n(idx_H2CO)

    !d[NH3_dot]/d[OH]
    pd(16,10) =  &
        +k(1033)*n(idx_NH2)  &
        -k(1099)*n(idx_NH3)

    !d[NO_dot]/d[OH]
    pd(17,10) =  &
        +k(1098)*n(idx_HNO)  &
        -k(1100)*n(idx_NO)  &
        +k(1026)*n(idx_N)  &
        +k(857)*n(idx_HNOj)

    !d[SI_dot]/d[OH]
    pd(18,10) =  &
        -k(1103)*n(idx_SI)

    !d[CN_dot]/d[OH]
    pd(24,10) =  &
        -k(1092)*n(idx_CN)  &
        +k(226)*n(idx_CNj)  &
        +k(1095)*n(idx_HCN)  &
        -k(1091)*n(idx_CN)  &
        +k(854)*n(idx_HCNj)

    !d[CO_dot]/d[OH]
    pd(25,10) =  &
        -k(1093)*n(idx_CO)  &
        +k(1097)*n(idx_HCO)  &
        +k(876)*n(idx_C)  &
        +k(1096)*n(idx_HCN)  &
        +k(855)*n(idx_HCOj)  &
        +k(227)*n(idx_COj)

    !d[N2_dot]/d[OH]
    pd(26,10) =  &
        +k(858)*n(idx_N2Hj)  &
        +k(228)*n(idx_N2j)

    !d[NH2_dot]/d[OH]
    pd(27,10) =  &
        +k(1096)*n(idx_HCN)  &
        -k(1032)*n(idx_NH2)  &
        -k(1033)*n(idx_NH2)  &
        +k(1051)*n(idx_NH)  &
        +k(1099)*n(idx_NH3)

    !d[CH3_dot]/d[OH]
    pd(28,10) =  &
        +k(901)*n(idx_CH2)  &
        -k(919)*n(idx_CH3)  &
        +k(923)*n(idx_CH4)  &
        -k(918)*n(idx_CH3)  &
        -k(920)*n(idx_CH3)

    !d[CH4_dot]/d[OH]
    pd(29,10) =  &
        -k(923)*n(idx_CH4)  &
        +k(918)*n(idx_CH3)

    !d[N_dot]/d[OH]
    pd(30,10) =  &
        +k(769)*n(idx_NHj)  &
        +k(170)*n(idx_Nj)  &
        +k(1049)*n(idx_NH)  &
        -k(1026)*n(idx_N)  &
        -k(1027)*n(idx_N)

    !d[NH_dot]/d[OH]
    pd(31,10) =  &
        -k(1050)*n(idx_NH)  &
        -k(1049)*n(idx_NH)  &
        +k(1032)*n(idx_NH2)  &
        +k(1027)*n(idx_N)  &
        -k(1051)*n(idx_NH)

    !d[SIO_dot]/d[OH]
    pd(34,10) =  &
        +k(1103)*n(idx_SI)

    !d[HE_dot]/d[OH]
    pd(35,10) =  &
        +k(700)*n(idx_HEj)

    !d[HNO_dot]/d[OH]
    pd(36,10) =  &
        +k(1050)*n(idx_NH)  &
        -k(1098)*n(idx_HNO)

    !d[CO2_dot]/d[OH]
    pd(38,10) =  &
        +k(1093)*n(idx_CO)

    !d[NO2_dot]/d[OH]
    pd(42,10) =  &
        +k(1100)*n(idx_NO)

    !d[O2H_dot]/d[OH]
    pd(43,10) =  &
        -k(1101)*n(idx_O2H)

    !d[OCN_dot]/d[OH]
    pd(44,10) =  &
        +k(1092)*n(idx_CN)

    !d[H2O_DUST_dot]/d[OH]
    pd(55,10) =  &
        +k(1255)

    !d[HCO+_dot]/d[OH]
    pd(70,10) =  &
        -k(855)*n(idx_HCOj)  &
        +k(852)*n(idx_COj)  &
        -k(856)*n(idx_HCOj)

    !d[H+_dot]/d[OH]
    pd(71,10) =  &
        -k(91)*n(idx_Hj)

    !d[C+_dot]/d[OH]
    pd(73,10) =  &
        -k(380)*n(idx_Cj)

    !d[CH+_dot]/d[OH]
    pd(75,10) =  &
        -k(416)*n(idx_CHj)

    !d[H2CO+_dot]/d[OH]
    pd(76,10) =  &
        +k(446)*n(idx_CH3j)

    !d[SI+_dot]/d[OH]
    pd(80,10) =  &
        -k(860)*n(idx_SIj)

    !d[CN+_dot]/d[OH]
    pd(86,10) =  &
        -k(226)*n(idx_CNj)

    !d[CO+_dot]/d[OH]
    pd(87,10) =  &
        +k(380)*n(idx_Cj)  &
        +k(416)*n(idx_CHj)  &
        -k(852)*n(idx_COj)  &
        -k(227)*n(idx_COj)

    !d[N2+_dot]/d[OH]
    pd(88,10) =  &
        -k(228)*n(idx_N2j)

    !d[O2+_dot]/d[OH]
    pd(89,10) =  &
        +k(820)*n(idx_Oj)

    !d[H2O+_dot]/d[OH]
    pd(90,10) =  &
        +k(859)*n(idx_O2Hj)  &
        +k(528)*n(idx_H2j)  &
        +k(601)*n(idx_H3j)  &
        +k(848)*n(idx_OHj)  &
        -k(853)*n(idx_H2Oj)  &
        +k(854)*n(idx_HCNj)  &
        +k(855)*n(idx_HCOj)  &
        +k(769)*n(idx_NHj)  &
        +k(858)*n(idx_N2Hj)  &
        +k(857)*n(idx_HNOj)

    !d[O+_dot]/d[OH]
    pd(92,10) =  &
        -k(216)*n(idx_Oj)  &
        +k(700)*n(idx_HEj)  &
        -k(820)*n(idx_Oj)

    !d[OH+_dot]/d[OH]
    pd(93,10) =  &
        -k(848)*n(idx_OHj)  &
        +k(227)*n(idx_COj)  &
        +k(170)*n(idx_Nj)  &
        +k(1176)  &
        +k(115)*n(idx_H2j)  &
        +k(91)*n(idx_Hj)  &
        +k(216)*n(idx_Oj)  &
        +k(226)*n(idx_CNj)  &
        +k(228)*n(idx_N2j)

    !d[CH3+_dot]/d[OH]
    pd(94,10) =  &
        -k(446)*n(idx_CH3j)

    !d[N+_dot]/d[OH]
    pd(96,10) =  &
        -k(170)*n(idx_Nj)

    !d[HCN+_dot]/d[OH]
    pd(97,10) =  &
        -k(854)*n(idx_HCNj)

    !d[NH+_dot]/d[OH]
    pd(98,10) =  &
        -k(769)*n(idx_NHj)

    !d[SIO+_dot]/d[OH]
    pd(101,10) =  &
        +k(860)*n(idx_SIj)

    !d[H2+_dot]/d[OH]
    pd(102,10) =  &
        -k(115)*n(idx_H2j)  &
        -k(528)*n(idx_H2j)

    !d[HE+_dot]/d[OH]
    pd(103,10) =  &
        -k(700)*n(idx_HEj)

    !d[HNO+_dot]/d[OH]
    pd(104,10) =  &
        -k(857)*n(idx_HNOj)

    !d[H3+_dot]/d[OH]
    pd(106,10) =  &
        -k(601)*n(idx_H3j)

    !d[H3O+_dot]/d[OH]
    pd(108,10) =  &
        +k(853)*n(idx_H2Oj)

    !d[HCO2+_dot]/d[OH]
    pd(110,10) =  &
        +k(856)*n(idx_HCOj)

    !d[N2H+_dot]/d[OH]
    pd(112,10) =  &
        -k(858)*n(idx_N2Hj)

    !d[O2H+_dot]/d[OH]
    pd(113,10) =  &
        -k(859)*n(idx_O2Hj)

    !d[E_dot]/d[O2]
    pd(1,11) =  &
        +k(1169)  &
        +k(280)

    !d[CH_dot]/d[O2]
    pd(2,11) =  &
        -k(937)*n(idx_CH)  &
        -k(936)*n(idx_CH)  &
        -k(935)*n(idx_CH)  &
        -k(934)*n(idx_CH)

    !d[O_dot]/d[O2]
    pd(3,11) =  &
        +k(937)*n(idx_CH)  &
        +k(1107)*n(idx_SI)  &
        +k(1045)*n(idx_NH)  &
        +k(215)*n(idx_Oj)  &
        +2.d0*k(13)*n(idx_H)  &
        +k(954)*n(idx_CO)  &
        +2.d0*k(281)  &
        +k(1053)*n(idx_NO)  &
        +2.d0*k(7)*n(idx_H2)  &
        +k(697)*n(idx_HEj)  &
        +k(1024)*n(idx_N)  &
        +k(729)*n(idx_Nj)  &
        +k(377)*n(idx_Cj)  &
        +k(443)*n(idx_CH3j)  &
        +k(935)*n(idx_CH)  &
        +k(413)*n(idx_CHj)  &
        +k(950)*n(idx_CN)  &
        +k(874)*n(idx_C)  &
        +k(990)*n(idx_H)  &
        +k(778)*n(idx_NH2j)  &
        +k(893)*n(idx_CH2)  &
        +2.d0*k(1170)

    !d[HCN_dot]/d[O2]
    pd(5,11) =  &
        +k(134)*n(idx_HCNj)

    !d[H2_dot]/d[O2]
    pd(6,11) =  &
        -k(7)*n(idx_H2)  &
        +k(114)*n(idx_H2j)  &
        +k(598)*n(idx_H3j)  &
        -k(965)*n(idx_H2)  &
        -k(964)*n(idx_H2)  &
        +k(7)*n(idx_H2)  &
        +k(890)*n(idx_CH2)

    !d[C_dot]/d[O2]
    pd(7,11) =  &
        -k(874)*n(idx_C)

    !d[H_dot]/d[O2]
    pd(8,11) =  &
        +k(964)*n(idx_H2)  &
        +k(89)*n(idx_Hj)  &
        +k(934)*n(idx_CH)  &
        +k(935)*n(idx_CH)  &
        +2.d0*k(891)*n(idx_CH2)  &
        -k(990)*n(idx_H)  &
        -k(13)*n(idx_H)  &
        +k(13)*n(idx_H)  &
        +k(526)*n(idx_H2j)

    !d[H2O_dot]/d[O2]
    pd(9,11) =  &
        +k(122)*n(idx_H2Oj)  &
        +k(913)*n(idx_CH3)  &
        +k(892)*n(idx_CH2)

    !d[OH_dot]/d[O2]
    pd(10,11) =  &
        +k(990)*n(idx_H)  &
        +k(863)*n(idx_SIH2j)  &
        +k(779)*n(idx_NH2j)  &
        +2.d0*k(965)*n(idx_H2)  &
        +k(412)*n(idx_CHj)  &
        +k(766)*n(idx_NHj)  &
        +k(912)*n(idx_CH3)  &
        +k(894)*n(idx_CH2)  &
        +k(421)*n(idx_CH2j)  &
        +k(936)*n(idx_CH)  &
        +k(1002)*n(idx_HCO)  &
        +k(225)*n(idx_OHj)  &
        +k(1046)*n(idx_NH)

    !d[O2_dot]/d[O2]
    pd(11,11) =  &
        -k(7)*n(idx_H2)  &
        -k(215)*n(idx_Oj)  &
        -k(122)*n(idx_H2Oj)  &
        -k(89)*n(idx_Hj)  &
        -k(1169)  &
        -k(377)*n(idx_Cj)  &
        -k(443)*n(idx_CH3j)  &
        -k(1002)*n(idx_HCO)  &
        -k(730)*n(idx_Nj)  &
        -k(990)*n(idx_H)  &
        -k(52)*n(idx_CH4j)  &
        -k(281)  &
        -k(1053)*n(idx_NO)  &
        -k(729)*n(idx_Nj)  &
        -k(598)*n(idx_H3j)  &
        -k(914)*n(idx_CH3)  &
        -k(937)*n(idx_CH)  &
        -k(949)*n(idx_CN)  &
        -k(964)*n(idx_H2)  &
        -k(954)*n(idx_CO)  &
        -k(1293)  &
        -k(378)*n(idx_Cj)  &
        -k(922)*n(idx_CH4)  &
        -k(13)*n(idx_H)  &
        -k(147)*n(idx_HEj)  &
        -k(913)*n(idx_CH3)  &
        -k(1107)*n(idx_SI)  &
        -k(950)*n(idx_CN)  &
        -k(1170)  &
        -k(912)*n(idx_CH3)  &
        -k(934)*n(idx_CH)  &
        -k(412)*n(idx_CHj)  &
        -k(482)*n(idx_CNj)  &
        -k(779)*n(idx_NH2j)  &
        -k(965)*n(idx_H2)  &
        -k(74)*n(idx_COj)  &
        -k(778)*n(idx_NH2j)  &
        -k(280)  &
        -k(891)*n(idx_CH2)  &
        -k(134)*n(idx_HCNj)  &
        -k(526)*n(idx_H2j)  &
        -k(69)*n(idx_CNj)  &
        -k(114)*n(idx_H2j)  &
        -k(863)*n(idx_SIH2j)  &
        -k(1045)*n(idx_NH)  &
        -k(892)*n(idx_CH2)  &
        -k(1003)*n(idx_HCO)  &
        -k(169)*n(idx_Nj)  &
        -k(894)*n(idx_CH2)  &
        -k(767)*n(idx_NHj)  &
        -k(414)*n(idx_CHj)  &
        -k(1056)*n(idx_OCN)  &
        -k(874)*n(idx_C)  &
        -k(1055)*n(idx_OCN)  &
        -k(697)*n(idx_HEj)  &
        -k(413)*n(idx_CHj)  &
        -k(1024)*n(idx_N)  &
        -k(893)*n(idx_CH2)  &
        -k(936)*n(idx_CH)  &
        -k(935)*n(idx_CH)  &
        -k(225)*n(idx_OHj)  &
        -k(1046)*n(idx_NH)  &
        -k(174)*n(idx_N2j)  &
        -k(180)*n(idx_NHj)  &
        -k(550)*n(idx_H2COj)  &
        -k(890)*n(idx_CH2)  &
        -k(421)*n(idx_CH2j)  &
        -k(766)*n(idx_NHj)

    !d[CH2_dot]/d[O2]
    pd(12,11) =  &
        +k(914)*n(idx_CH3)  &
        -k(894)*n(idx_CH2)  &
        -k(893)*n(idx_CH2)  &
        -k(891)*n(idx_CH2)  &
        -k(892)*n(idx_CH2)  &
        -k(890)*n(idx_CH2)

    !d[H2CO_dot]/d[O2]
    pd(13,11) =  &
        +k(893)*n(idx_CH2)  &
        +k(912)*n(idx_CH3)

    !d[HCO_dot]/d[O2]
    pd(14,11) =  &
        -k(1002)*n(idx_HCO)  &
        +k(937)*n(idx_CH)  &
        +k(894)*n(idx_CH2)  &
        +k(913)*n(idx_CH3)  &
        +k(414)*n(idx_CHj)  &
        -k(1003)*n(idx_HCO)

    !d[NO_dot]/d[O2]
    pd(17,11) =  &
        +k(1055)*n(idx_OCN)  &
        +k(949)*n(idx_CN)  &
        +k(730)*n(idx_Nj)  &
        -k(1053)*n(idx_NO)  &
        +k(1024)*n(idx_N)  &
        +k(1046)*n(idx_NH)

    !d[SI_dot]/d[O2]
    pd(18,11) =  &
        -k(1107)*n(idx_SI)

    !d[CN_dot]/d[O2]
    pd(24,11) =  &
        -k(950)*n(idx_CN)  &
        -k(949)*n(idx_CN)  &
        +k(69)*n(idx_CNj)

    !d[CO_dot]/d[O2]
    pd(25,11) =  &
        +k(74)*n(idx_COj)  &
        +k(949)*n(idx_CN)  &
        +k(378)*n(idx_Cj)  &
        +k(874)*n(idx_C)  &
        +k(1003)*n(idx_HCO)  &
        +k(936)*n(idx_CH)  &
        +k(1056)*n(idx_OCN)  &
        +k(935)*n(idx_CH)  &
        +k(482)*n(idx_CNj)  &
        -k(954)*n(idx_CO)  &
        +k(892)*n(idx_CH2)

    !d[N2_dot]/d[O2]
    pd(26,11) =  &
        +k(174)*n(idx_N2j)

    !d[CH3_dot]/d[O2]
    pd(28,11) =  &
        -k(913)*n(idx_CH3)  &
        -k(914)*n(idx_CH3)  &
        -k(912)*n(idx_CH3)  &
        +k(922)*n(idx_CH4)

    !d[CH4_dot]/d[O2]
    pd(29,11) =  &
        +k(52)*n(idx_CH4j)  &
        -k(922)*n(idx_CH4)

    !d[N_dot]/d[O2]
    pd(30,11) =  &
        -k(1024)*n(idx_N)  &
        +k(767)*n(idx_NHj)  &
        +k(169)*n(idx_Nj)

    !d[NH_dot]/d[O2]
    pd(31,11) =  &
        -k(1046)*n(idx_NH)  &
        -k(1045)*n(idx_NH)  &
        +k(180)*n(idx_NHj)

    !d[SIO_dot]/d[O2]
    pd(34,11) =  &
        +k(1107)*n(idx_SI)

    !d[HE_dot]/d[O2]
    pd(35,11) =  &
        +k(697)*n(idx_HEj)  &
        +k(147)*n(idx_HEj)

    !d[HNO_dot]/d[O2]
    pd(36,11) =  &
        +k(1045)*n(idx_NH)

    !d[CO2_dot]/d[O2]
    pd(38,11) =  &
        +k(1055)*n(idx_OCN)  &
        +k(934)*n(idx_CH)  &
        +k(954)*n(idx_CO)  &
        +k(891)*n(idx_CH2)  &
        +k(1002)*n(idx_HCO)  &
        +k(890)*n(idx_CH2)

    !d[NO2_dot]/d[O2]
    pd(42,11) =  &
        +k(1056)*n(idx_OCN)  &
        +k(1053)*n(idx_NO)

    !d[O2H_dot]/d[O2]
    pd(43,11) =  &
        +k(914)*n(idx_CH3)  &
        +k(964)*n(idx_H2)  &
        +k(550)*n(idx_H2COj)  &
        +k(1003)*n(idx_HCO)  &
        +k(922)*n(idx_CH4)

    !d[OCN_dot]/d[O2]
    pd(44,11) =  &
        +k(950)*n(idx_CN)  &
        -k(1056)*n(idx_OCN)  &
        -k(1055)*n(idx_OCN)

    !d[O2_DUST_dot]/d[O2]
    pd(61,11) =  &
        +k(1293)

    !d[HCO+_dot]/d[O2]
    pd(70,11) =  &
        +k(413)*n(idx_CHj)  &
        +k(550)*n(idx_H2COj)  &
        +k(421)*n(idx_CH2j)

    !d[H+_dot]/d[O2]
    pd(71,11) =  &
        -k(89)*n(idx_Hj)

    !d[C+_dot]/d[O2]
    pd(73,11) =  &
        -k(378)*n(idx_Cj)  &
        -k(377)*n(idx_Cj)

    !d[CH2+_dot]/d[O2]
    pd(74,11) =  &
        -k(421)*n(idx_CH2j)

    !d[CH+_dot]/d[O2]
    pd(75,11) =  &
        -k(414)*n(idx_CHj)  &
        -k(412)*n(idx_CHj)  &
        -k(413)*n(idx_CHj)

    !d[H2CO+_dot]/d[O2]
    pd(76,11) =  &
        -k(550)*n(idx_H2COj)

    !d[NO+_dot]/d[O2]
    pd(79,11) =  &
        +k(766)*n(idx_NHj)  &
        +k(482)*n(idx_CNj)  &
        +k(729)*n(idx_Nj)

    !d[SIH2+_dot]/d[O2]
    pd(84,11) =  &
        -k(863)*n(idx_SIH2j)

    !d[CN+_dot]/d[O2]
    pd(86,11) =  &
        -k(69)*n(idx_CNj)  &
        -k(482)*n(idx_CNj)

    !d[CO+_dot]/d[O2]
    pd(87,11) =  &
        +k(412)*n(idx_CHj)  &
        +k(377)*n(idx_Cj)  &
        -k(74)*n(idx_COj)

    !d[N2+_dot]/d[O2]
    pd(88,11) =  &
        -k(174)*n(idx_N2j)

    !d[O2+_dot]/d[O2]
    pd(89,11) =  &
        +k(74)*n(idx_COj)  &
        +k(114)*n(idx_H2j)  &
        +k(215)*n(idx_Oj)  &
        +k(180)*n(idx_NHj)  &
        +k(122)*n(idx_H2Oj)  &
        +k(225)*n(idx_OHj)  &
        +k(147)*n(idx_HEj)  &
        +k(169)*n(idx_Nj)  &
        +k(1169)  &
        +k(280)  &
        +k(69)*n(idx_CNj)  &
        +k(134)*n(idx_HCNj)  &
        +k(52)*n(idx_CH4j)  &
        +k(89)*n(idx_Hj)  &
        +k(174)*n(idx_N2j)

    !d[H2O+_dot]/d[O2]
    pd(90,11) =  &
        -k(122)*n(idx_H2Oj)

    !d[NH2+_dot]/d[O2]
    pd(91,11) =  &
        -k(779)*n(idx_NH2j)  &
        -k(778)*n(idx_NH2j)

    !d[O+_dot]/d[O2]
    pd(92,11) =  &
        -k(215)*n(idx_Oj)  &
        +k(414)*n(idx_CHj)  &
        +k(378)*n(idx_Cj)  &
        +k(730)*n(idx_Nj)  &
        +k(697)*n(idx_HEj)

    !d[OH+_dot]/d[O2]
    pd(93,11) =  &
        -k(225)*n(idx_OHj)

    !d[CH3+_dot]/d[O2]
    pd(94,11) =  &
        -k(443)*n(idx_CH3j)

    !d[CH4+_dot]/d[O2]
    pd(95,11) =  &
        -k(52)*n(idx_CH4j)

    !d[N+_dot]/d[O2]
    pd(96,11) =  &
        -k(729)*n(idx_Nj)  &
        -k(169)*n(idx_Nj)  &
        -k(730)*n(idx_Nj)

    !d[HCN+_dot]/d[O2]
    pd(97,11) =  &
        -k(134)*n(idx_HCNj)

    !d[NH+_dot]/d[O2]
    pd(98,11) =  &
        -k(180)*n(idx_NHj)  &
        -k(767)*n(idx_NHj)  &
        -k(766)*n(idx_NHj)

    !d[H2+_dot]/d[O2]
    pd(102,11) =  &
        -k(526)*n(idx_H2j)  &
        -k(114)*n(idx_H2j)

    !d[HE+_dot]/d[O2]
    pd(103,11) =  &
        -k(147)*n(idx_HEj)  &
        -k(697)*n(idx_HEj)

    !d[HNO+_dot]/d[O2]
    pd(104,11) =  &
        +k(779)*n(idx_NH2j)

    !d[H2NO+_dot]/d[O2]
    pd(105,11) =  &
        +k(778)*n(idx_NH2j)

    !d[H3+_dot]/d[O2]
    pd(106,11) =  &
        -k(598)*n(idx_H3j)

    !d[H3CO+_dot]/d[O2]
    pd(107,11) =  &
        +k(443)*n(idx_CH3j)

    !d[O2H+_dot]/d[O2]
    pd(113,11) =  &
        +k(767)*n(idx_NHj)  &
        +k(598)*n(idx_H3j)  &
        +k(526)*n(idx_H2j)

    !d[SIOH+_dot]/d[O2]
    pd(115,11) =  &
        +k(863)*n(idx_SIH2j)

    !d[E_dot]/d[CH2]
    pd(1,12) =  &
        +k(1113)  &
        +k(243)

    !d[CH_dot]/d[CH2]
    pd(2,12) =  &
        +k(881)*n(idx_CN)  &
        +k(244)  &
        +2.d0*k(864)*n(idx_C)  &
        +2.d0*k(879)*n(idx_CH2)  &
        +k(968)*n(idx_H)  &
        +k(1114)  &
        +k(898)*n(idx_O)  &
        +k(1008)*n(idx_N)  &
        +k(423)*n(idx_COj)  &
        +k(900)*n(idx_OH)

    !d[O_dot]/d[CH2]
    pd(3,12) =  &
        -k(896)*n(idx_O)  &
        +k(44)*n(idx_Oj)  &
        +k(438)*n(idx_OHj)  &
        -k(895)*n(idx_O)  &
        +k(893)*n(idx_O2)  &
        -k(897)*n(idx_O)  &
        +k(901)*n(idx_OH)  &
        -k(898)*n(idx_O)  &
        +k(436)*n(idx_O2j)

    !d[HNC_dot]/d[CH2]
    pd(4,12) =  &
        +k(429)*n(idx_HCNHj)  &
        +k(1007)*n(idx_N)

    !d[HCN_dot]/d[CH2]
    pd(5,12) =  &
        +k(888)*n(idx_NO)  &
        +k(885)*n(idx_N2)  &
        +k(428)*n(idx_HCNHj)  &
        +k(1006)*n(idx_N)  &
        +k(881)*n(idx_CN)

    !d[H2_dot]/d[CH2]
    pd(6,12) =  &
        -k(957)*n(idx_H2)  &
        +k(492)*n(idx_Hj)  &
        +k(968)*n(idx_H)  &
        +k(101)*n(idx_H2j)  &
        +k(890)*n(idx_O2)  &
        +k(654)*n(idx_HEj)  &
        +k(578)*n(idx_H3j)  &
        +k(895)*n(idx_O)

    !d[C_dot]/d[CH2]
    pd(7,12) =  &
        +k(15)*n(idx_Cj)  &
        -k(864)*n(idx_C)

    !d[H_dot]/d[CH2]
    pd(8,12) =  &
        -k(968)*n(idx_H)  &
        +k(76)*n(idx_Hj)  &
        +k(957)*n(idx_H2)  &
        +k(244)  &
        +k(511)*n(idx_H2j)  &
        +k(897)*n(idx_O)  &
        +k(1114)  &
        +k(1006)*n(idx_N)  &
        +k(1007)*n(idx_N)  &
        +2.d0*k(891)*n(idx_O2)  &
        +k(655)*n(idx_HEj)  &
        +k(899)*n(idx_OH)  &
        +k(889)*n(idx_NO)  &
        +2.d0*k(896)*n(idx_O)

    !d[H2O_dot]/d[CH2]
    pd(9,12) =  &
        +k(41)*n(idx_H2Oj)  &
        +k(426)*n(idx_H3Oj)  &
        +k(900)*n(idx_OH)  &
        +k(892)*n(idx_O2)

    !d[OH_dot]/d[CH2]
    pd(10,12) =  &
        +k(894)*n(idx_O2)  &
        -k(899)*n(idx_OH)  &
        +k(898)*n(idx_O)  &
        +k(425)*n(idx_H2Oj)  &
        -k(901)*n(idx_OH)  &
        +k(888)*n(idx_NO)  &
        +k(46)*n(idx_OHj)  &
        -k(900)*n(idx_OH)

    !d[O2_dot]/d[CH2]
    pd(11,12) =  &
        -k(891)*n(idx_O2)  &
        -k(892)*n(idx_O2)  &
        -k(893)*n(idx_O2)  &
        +k(437)*n(idx_O2Hj)  &
        -k(890)*n(idx_O2)  &
        +k(45)*n(idx_O2j)  &
        -k(894)*n(idx_O2)

    !d[CH2_dot]/d[CH2]
    pd(12,12) =  &
        -k(430)*n(idx_HCOj)  &
        -k(968)*n(idx_H)  &
        -k(433)*n(idx_NHj)  &
        -k(40)*n(idx_H2COj)  &
        -k(897)*n(idx_O)  &
        -k(895)*n(idx_O)  &
        -k(899)*n(idx_OH)  &
        -k(886)*n(idx_NO2)  &
        -k(424)*n(idx_H2COj)  &
        -k(883)*n(idx_HCO)  &
        -k(901)*n(idx_OH)  &
        -k(864)*n(idx_C)  &
        -k(893)*n(idx_O2)  &
        -k(46)*n(idx_OHj)  &
        -k(76)*n(idx_Hj)  &
        -4.d0*k(879)*n(idx_CH2)  &
        -k(423)*n(idx_COj)  &
        -k(654)*n(idx_HEj)  &
        -k(15)*n(idx_Cj)  &
        -k(45)*n(idx_O2j)  &
        -k(492)*n(idx_Hj)  &
        -k(1006)*n(idx_N)  &
        -k(896)*n(idx_O)  &
        -k(511)*n(idx_H2j)  &
        -k(427)*n(idx_HCNj)  &
        -k(1113)  &
        -k(435)*n(idx_NH3j)  &
        -k(425)*n(idx_H2Oj)  &
        -k(428)*n(idx_HCNHj)  &
        -k(578)*n(idx_H3j)  &
        -k(892)*n(idx_O2)  &
        -k(38)*n(idx_CNj)  &
        -k(101)*n(idx_H2j)  &
        -k(439)*n(idx_SIOj)  &
        -k(39)*n(idx_COj)  &
        -k(887)*n(idx_NO)  &
        -k(243)  &
        -k(880)*n(idx_CH4)  &
        -k(438)*n(idx_OHj)  &
        -k(432)*n(idx_N2Hj)  &
        -k(888)*n(idx_NO)  &
        -k(885)*n(idx_N2)  &
        -k(882)*n(idx_H2CO)  &
        -k(244)  &
        -k(1257)  &
        -k(957)*n(idx_H2)  &
        -k(884)*n(idx_HNO)  &
        -k(655)*n(idx_HEj)  &
        -k(1008)*n(idx_N)  &
        -k(889)*n(idx_NO)  &
        -k(894)*n(idx_O2)  &
        -k(429)*n(idx_HCNHj)  &
        -k(43)*n(idx_NH2j)  &
        -k(436)*n(idx_O2j)  &
        -k(881)*n(idx_CN)  &
        -k(156)*n(idx_Nj)  &
        -k(891)*n(idx_O2)  &
        -k(898)*n(idx_O)  &
        -k(1114)  &
        -k(434)*n(idx_NH2j)  &
        -k(42)*n(idx_N2j)  &
        -k(41)*n(idx_H2Oj)  &
        -k(44)*n(idx_Oj)  &
        -k(431)*n(idx_HNOj)  &
        -k(890)*n(idx_O2)  &
        -k(900)*n(idx_OH)  &
        -k(426)*n(idx_H3Oj)  &
        -k(1007)*n(idx_N)  &
        -k(437)*n(idx_O2Hj)

    !d[H2CO_dot]/d[CH2]
    pd(13,12) =  &
        +k(439)*n(idx_SIOj)  &
        +k(886)*n(idx_NO2)  &
        -k(882)*n(idx_H2CO)  &
        +k(893)*n(idx_O2)  &
        +k(40)*n(idx_H2COj)  &
        +k(887)*n(idx_NO)  &
        +k(899)*n(idx_OH)

    !d[HCO_dot]/d[CH2]
    pd(14,12) =  &
        +k(424)*n(idx_H2COj)  &
        +k(897)*n(idx_O)  &
        +k(882)*n(idx_H2CO)  &
        -k(883)*n(idx_HCO)  &
        +k(894)*n(idx_O2)

    !d[NO_dot]/d[CH2]
    pd(17,12) =  &
        +k(886)*n(idx_NO2)  &
        -k(889)*n(idx_NO)  &
        +k(431)*n(idx_HNOj)  &
        +k(884)*n(idx_HNO)  &
        -k(888)*n(idx_NO)  &
        -k(887)*n(idx_NO)

    !d[CN_dot]/d[CH2]
    pd(24,12) =  &
        +k(427)*n(idx_HCNj)  &
        +k(38)*n(idx_CNj)  &
        -k(881)*n(idx_CN)

    !d[CO_dot]/d[CH2]
    pd(25,12) =  &
        +k(39)*n(idx_COj)  &
        +k(883)*n(idx_HCO)  &
        +k(895)*n(idx_O)  &
        +k(430)*n(idx_HCOj)  &
        +k(892)*n(idx_O2)  &
        +k(896)*n(idx_O)

    !d[N2_dot]/d[CH2]
    pd(26,12) =  &
        +k(42)*n(idx_N2j)  &
        +k(432)*n(idx_N2Hj)  &
        -k(885)*n(idx_N2)

    !d[NH2_dot]/d[CH2]
    pd(27,12) =  &
        +k(435)*n(idx_NH3j)  &
        +k(43)*n(idx_NH2j)

    !d[CH3_dot]/d[CH2]
    pd(28,12) =  &
        +k(957)*n(idx_H2)  &
        +k(883)*n(idx_HCO)  &
        +2.d0*k(879)*n(idx_CH2)  &
        +k(884)*n(idx_HNO)  &
        +k(901)*n(idx_OH)  &
        +k(882)*n(idx_H2CO)  &
        +2.d0*k(880)*n(idx_CH4)

    !d[CH4_dot]/d[CH2]
    pd(29,12) =  &
        -k(880)*n(idx_CH4)

    !d[N_dot]/d[CH2]
    pd(30,12) =  &
        +k(433)*n(idx_NHj)  &
        -k(1007)*n(idx_N)  &
        -k(1006)*n(idx_N)  &
        +k(887)*n(idx_NO)  &
        +k(156)*n(idx_Nj)  &
        -k(1008)*n(idx_N)

    !d[NH_dot]/d[CH2]
    pd(31,12) =  &
        +k(1008)*n(idx_N)  &
        +k(885)*n(idx_N2)  &
        +k(434)*n(idx_NH2j)

    !d[HE_dot]/d[CH2]
    pd(35,12) =  &
        +k(654)*n(idx_HEj)  &
        +k(655)*n(idx_HEj)

    !d[HNO_dot]/d[CH2]
    pd(36,12) =  &
        -k(884)*n(idx_HNO)

    !d[CO2_dot]/d[CH2]
    pd(38,12) =  &
        +k(891)*n(idx_O2)  &
        +k(890)*n(idx_O2)

    !d[HNCO_dot]/d[CH2]
    pd(41,12) =  &
        +k(889)*n(idx_NO)

    !d[NO2_dot]/d[CH2]
    pd(42,12) =  &
        -k(886)*n(idx_NO2)

    !d[CH4_DUST_dot]/d[CH2]
    pd(53,12) =  &
        +k(1257)

    !d[HCO+_dot]/d[CH2]
    pd(70,12) =  &
        -k(430)*n(idx_HCOj)  &
        +k(423)*n(idx_COj)

    !d[H+_dot]/d[CH2]
    pd(71,12) =  &
        -k(492)*n(idx_Hj)  &
        -k(76)*n(idx_Hj)

    !d[C+_dot]/d[CH2]
    pd(73,12) =  &
        +k(654)*n(idx_HEj)  &
        -k(15)*n(idx_Cj)

    !d[CH2+_dot]/d[CH2]
    pd(74,12) =  &
        +k(39)*n(idx_COj)  &
        +k(44)*n(idx_Oj)  &
        +k(76)*n(idx_Hj)  &
        +k(41)*n(idx_H2Oj)  &
        +k(40)*n(idx_H2COj)  &
        +k(101)*n(idx_H2j)  &
        +k(243)  &
        +k(1113)  &
        +k(156)*n(idx_Nj)  &
        +k(15)*n(idx_Cj)  &
        +k(45)*n(idx_O2j)  &
        +k(42)*n(idx_N2j)  &
        +k(43)*n(idx_NH2j)  &
        +k(46)*n(idx_OHj)  &
        +k(38)*n(idx_CNj)

    !d[CH+_dot]/d[CH2]
    pd(75,12) =  &
        +k(655)*n(idx_HEj)  &
        +k(492)*n(idx_Hj)

    !d[H2CO+_dot]/d[CH2]
    pd(76,12) =  &
        -k(424)*n(idx_H2COj)  &
        -k(40)*n(idx_H2COj)  &
        +k(436)*n(idx_O2j)

    !d[NH3+_dot]/d[CH2]
    pd(78,12) =  &
        -k(435)*n(idx_NH3j)

    !d[SI+_dot]/d[CH2]
    pd(80,12) =  &
        +k(439)*n(idx_SIOj)

    !d[CN+_dot]/d[CH2]
    pd(86,12) =  &
        -k(38)*n(idx_CNj)

    !d[CO+_dot]/d[CH2]
    pd(87,12) =  &
        -k(423)*n(idx_COj)  &
        -k(39)*n(idx_COj)

    !d[N2+_dot]/d[CH2]
    pd(88,12) =  &
        -k(42)*n(idx_N2j)

    !d[O2+_dot]/d[CH2]
    pd(89,12) =  &
        -k(436)*n(idx_O2j)  &
        -k(45)*n(idx_O2j)

    !d[H2O+_dot]/d[CH2]
    pd(90,12) =  &
        -k(425)*n(idx_H2Oj)  &
        -k(41)*n(idx_H2Oj)

    !d[NH2+_dot]/d[CH2]
    pd(91,12) =  &
        -k(43)*n(idx_NH2j)  &
        -k(434)*n(idx_NH2j)

    !d[O+_dot]/d[CH2]
    pd(92,12) =  &
        -k(44)*n(idx_Oj)

    !d[OH+_dot]/d[CH2]
    pd(93,12) =  &
        -k(46)*n(idx_OHj)  &
        -k(438)*n(idx_OHj)

    !d[CH3+_dot]/d[CH2]
    pd(94,12) =  &
        +k(433)*n(idx_NHj)  &
        +k(438)*n(idx_OHj)  &
        +k(429)*n(idx_HCNHj)  &
        +k(431)*n(idx_HNOj)  &
        +k(432)*n(idx_N2Hj)  &
        +k(428)*n(idx_HCNHj)  &
        +k(424)*n(idx_H2COj)  &
        +k(511)*n(idx_H2j)  &
        +k(430)*n(idx_HCOj)  &
        +k(425)*n(idx_H2Oj)  &
        +k(437)*n(idx_O2Hj)  &
        +k(427)*n(idx_HCNj)  &
        +k(578)*n(idx_H3j)  &
        +k(426)*n(idx_H3Oj)  &
        +k(434)*n(idx_NH2j)  &
        +k(435)*n(idx_NH3j)

    !d[N+_dot]/d[CH2]
    pd(96,12) =  &
        -k(156)*n(idx_Nj)

    !d[HCN+_dot]/d[CH2]
    pd(97,12) =  &
        -k(427)*n(idx_HCNj)

    !d[NH+_dot]/d[CH2]
    pd(98,12) =  &
        -k(433)*n(idx_NHj)

    !d[SIO+_dot]/d[CH2]
    pd(101,12) =  &
        -k(439)*n(idx_SIOj)

    !d[H2+_dot]/d[CH2]
    pd(102,12) =  &
        -k(101)*n(idx_H2j)  &
        -k(511)*n(idx_H2j)

    !d[HE+_dot]/d[CH2]
    pd(103,12) =  &
        -k(654)*n(idx_HEj)  &
        -k(655)*n(idx_HEj)

    !d[HNO+_dot]/d[CH2]
    pd(104,12) =  &
        -k(431)*n(idx_HNOj)

    !d[H3+_dot]/d[CH2]
    pd(106,12) =  &
        -k(578)*n(idx_H3j)

    !d[H3O+_dot]/d[CH2]
    pd(108,12) =  &
        -k(426)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[CH2]
    pd(109,12) =  &
        -k(428)*n(idx_HCNHj)  &
        -k(429)*n(idx_HCNHj)

    !d[N2H+_dot]/d[CH2]
    pd(112,12) =  &
        -k(432)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CH2]
    pd(113,12) =  &
        -k(437)*n(idx_O2Hj)

    !d[E_dot]/d[H2CO]
    pd(1,13) =  &
        +k(1139)  &
        +k(1140)

    !d[CH_dot]/d[H2CO]
    pd(2,13) =  &
        -k(925)*n(idx_CH)  &
        +k(370)*n(idx_Cj)

    !d[O_dot]/d[H2CO]
    pd(3,13) =  &
        +k(840)*n(idx_OHj)  &
        -k(1062)*n(idx_O)  &
        +k(210)*n(idx_Oj)  &
        +k(673)*n(idx_HEj)

    !d[HNC_dot]/d[H2CO]
    pd(4,13) =  &
        +k(635)*n(idx_HCNHj)

    !d[HCN_dot]/d[H2CO]
    pd(5,13) =  &
        +k(634)*n(idx_HCNHj)  &
        +k(943)*n(idx_CN)  &
        +k(480)*n(idx_CNj)

    !d[H2_dot]/d[H2CO]
    pd(6,13) =  &
        +k(498)*n(idx_Hj)  &
        +k(1137)  &
        +k(256)  &
        +k(586)*n(idx_H3j)  &
        +k(499)*n(idx_Hj)  &
        +k(106)*n(idx_H2j)  &
        +k(975)*n(idx_H)  &
        +k(671)*n(idx_HEj)  &
        +k(518)*n(idx_H2j)

    !d[C_dot]/d[H2CO]
    pd(7,13) =  &
        +k(401)*n(idx_CHj)  &
        +k(17)*n(idx_Cj)

    !d[H_dot]/d[H2CO]
    pd(8,13) =  &
        +k(80)*n(idx_Hj)  &
        +k(498)*n(idx_Hj)  &
        +k(672)*n(idx_HEj)  &
        +2.d0*k(1138)  &
        +k(731)*n(idx_N2j)  &
        +k(552)*n(idx_O2j)  &
        +k(518)*n(idx_H2j)  &
        +k(1140)  &
        -k(975)*n(idx_H)

    !d[H2O_dot]/d[H2CO]
    pd(9,13) =  &
        +k(1094)*n(idx_OH)  &
        +k(118)*n(idx_H2Oj)  &
        +k(608)*n(idx_H3Oj)

    !d[OH_dot]/d[H2CO]
    pd(10,13) =  &
        -k(1094)*n(idx_OH)  &
        +k(220)*n(idx_OHj)  &
        +k(555)*n(idx_H2Oj)  &
        +k(814)*n(idx_Oj)  &
        +k(1062)*n(idx_O)

    !d[O2_dot]/d[H2CO]
    pd(11,13) =  &
        +k(117)*n(idx_O2j)  &
        +k(553)*n(idx_O2Hj)  &
        +k(552)*n(idx_O2j)

    !d[CH2_dot]/d[H2CO]
    pd(12,13) =  &
        +k(925)*n(idx_CH)  &
        +k(402)*n(idx_CHj)  &
        +k(723)*n(idx_Nj)  &
        -k(882)*n(idx_CH2)

    !d[H2CO_dot]/d[H2CO]
    pd(13,13) =  &
        -k(636)*n(idx_HCOj)  &
        -k(17)*n(idx_Cj)  &
        -k(143)*n(idx_HEj)  &
        -k(401)*n(idx_CHj)  &
        -k(731)*n(idx_N2j)  &
        -k(975)*n(idx_H)  &
        -k(753)*n(idx_NHj)  &
        -k(814)*n(idx_Oj)  &
        -k(1253)  &
        -k(771)*n(idx_NH2j)  &
        -k(723)*n(idx_Nj)  &
        -k(925)*n(idx_CH)  &
        -k(608)*n(idx_H3Oj)  &
        -k(1139)  &
        -k(369)*n(idx_Cj)  &
        -k(904)*n(idx_CH3)  &
        -k(220)*n(idx_OHj)  &
        -k(943)*n(idx_CN)  &
        -k(71)*n(idx_COj)  &
        -k(722)*n(idx_Nj)  &
        -k(402)*n(idx_CHj)  &
        -k(118)*n(idx_H2Oj)  &
        -k(498)*n(idx_Hj)  &
        -k(160)*n(idx_Nj)  &
        -k(1094)*n(idx_OH)  &
        -k(840)*n(idx_OHj)  &
        -k(586)*n(idx_H3j)  &
        -k(400)*n(idx_CHj)  &
        -k(754)*n(idx_NHj)  &
        -k(480)*n(idx_CNj)  &
        -k(553)*n(idx_O2Hj)  &
        -k(671)*n(idx_HEj)  &
        -k(1140)  &
        -k(673)*n(idx_HEj)  &
        -k(80)*n(idx_Hj)  &
        -k(549)*n(idx_H2COj)  &
        -k(441)*n(idx_CH3j)  &
        -k(634)*n(idx_HCNHj)  &
        -k(176)*n(idx_NHj)  &
        -k(623)*n(idx_HCNj)  &
        -k(65)*n(idx_CNj)  &
        -k(672)*n(idx_HEj)  &
        -k(635)*n(idx_HCNHj)  &
        -k(418)*n(idx_CH2j)  &
        -k(1062)*n(idx_O)  &
        -k(210)*n(idx_Oj)  &
        -k(117)*n(idx_O2j)  &
        -k(551)*n(idx_HNOj)  &
        -k(1138)  &
        -k(171)*n(idx_N2j)  &
        -k(555)*n(idx_H2Oj)  &
        -k(106)*n(idx_H2j)  &
        -k(736)*n(idx_N2Hj)  &
        -k(50)*n(idx_CH4j)  &
        -k(370)*n(idx_Cj)  &
        -k(485)*n(idx_COj)  &
        -k(450)*n(idx_CH4j)  &
        -k(1137)  &
        -k(256)  &
        -k(499)*n(idx_Hj)  &
        -k(552)*n(idx_O2j)  &
        -k(518)*n(idx_H2j)  &
        -k(770)*n(idx_NH2j)  &
        -k(882)*n(idx_CH2)

    !d[HCO_dot]/d[H2CO]
    pd(14,13) =  &
        +k(1094)*n(idx_OH)  &
        +k(943)*n(idx_CN)  &
        +k(1062)*n(idx_O)  &
        +k(549)*n(idx_H2COj)  &
        +k(485)*n(idx_COj)  &
        +k(975)*n(idx_H)  &
        +k(771)*n(idx_NH2j)  &
        +k(882)*n(idx_CH2)  &
        +k(925)*n(idx_CH)  &
        +k(904)*n(idx_CH3)

    !d[NO_dot]/d[H2CO]
    pd(17,13) =  &
        +k(551)*n(idx_HNOj)

    !d[CN_dot]/d[H2CO]
    pd(24,13) =  &
        +k(65)*n(idx_CNj)  &
        +k(623)*n(idx_HCNj)  &
        -k(943)*n(idx_CN)

    !d[CO_dot]/d[H2CO]
    pd(25,13) =  &
        +k(1137)  &
        +k(256)  &
        +k(636)*n(idx_HCOj)  &
        +k(400)*n(idx_CHj)  &
        +k(71)*n(idx_COj)  &
        +k(369)*n(idx_Cj)  &
        +k(1138)

    !d[N2_dot]/d[H2CO]
    pd(26,13) =  &
        +k(736)*n(idx_N2Hj)  &
        +k(731)*n(idx_N2j)  &
        +k(171)*n(idx_N2j)

    !d[NH2_dot]/d[H2CO]
    pd(27,13) =  &
        +k(754)*n(idx_NHj)

    !d[CH3_dot]/d[H2CO]
    pd(28,13) =  &
        -k(904)*n(idx_CH3)  &
        +k(882)*n(idx_CH2)  &
        +k(450)*n(idx_CH4j)  &
        +k(418)*n(idx_CH2j)

    !d[CH4_dot]/d[H2CO]
    pd(29,13) =  &
        +k(441)*n(idx_CH3j)  &
        +k(50)*n(idx_CH4j)  &
        +k(904)*n(idx_CH3)

    !d[N_dot]/d[H2CO]
    pd(30,13) =  &
        +k(753)*n(idx_NHj)  &
        +k(160)*n(idx_Nj)

    !d[NH_dot]/d[H2CO]
    pd(31,13) =  &
        +k(770)*n(idx_NH2j)  &
        +k(722)*n(idx_Nj)  &
        +k(176)*n(idx_NHj)

    !d[HE_dot]/d[H2CO]
    pd(35,13) =  &
        +k(143)*n(idx_HEj)  &
        +k(672)*n(idx_HEj)  &
        +k(673)*n(idx_HEj)  &
        +k(671)*n(idx_HEj)

    !d[H2CO_DUST_dot]/d[H2CO]
    pd(47,13) =  &
        +k(1253)

    !d[HCO+_dot]/d[H2CO]
    pd(70,13) =  &
        -k(636)*n(idx_HCOj)  &
        +k(731)*n(idx_N2j)  &
        +k(370)*n(idx_Cj)  &
        +k(485)*n(idx_COj)  &
        +k(722)*n(idx_Nj)  &
        +k(402)*n(idx_CHj)  &
        +k(499)*n(idx_Hj)  &
        +k(754)*n(idx_NHj)  &
        +k(480)*n(idx_CNj)  &
        +k(814)*n(idx_Oj)  &
        +k(418)*n(idx_CH2j)  &
        +k(518)*n(idx_H2j)  &
        +k(672)*n(idx_HEj)  &
        +k(552)*n(idx_O2j)  &
        +k(441)*n(idx_CH3j)  &
        +k(1140)

    !d[H+_dot]/d[H2CO]
    pd(71,13) =  &
        -k(80)*n(idx_Hj)  &
        -k(499)*n(idx_Hj)  &
        -k(498)*n(idx_Hj)

    !d[C+_dot]/d[H2CO]
    pd(73,13) =  &
        -k(369)*n(idx_Cj)  &
        -k(17)*n(idx_Cj)  &
        -k(370)*n(idx_Cj)

    !d[CH2+_dot]/d[H2CO]
    pd(74,13) =  &
        -k(418)*n(idx_CH2j)  &
        +k(673)*n(idx_HEj)  &
        +k(369)*n(idx_Cj)

    !d[CH+_dot]/d[H2CO]
    pd(75,13) =  &
        -k(402)*n(idx_CHj)  &
        -k(401)*n(idx_CHj)  &
        -k(400)*n(idx_CHj)

    !d[H2CO+_dot]/d[H2CO]
    pd(76,13) =  &
        +k(80)*n(idx_Hj)  &
        +k(220)*n(idx_OHj)  &
        +k(160)*n(idx_Nj)  &
        -k(549)*n(idx_H2COj)  &
        +k(118)*n(idx_H2Oj)  &
        +k(143)*n(idx_HEj)  &
        +k(176)*n(idx_NHj)  &
        +k(17)*n(idx_Cj)  &
        +k(106)*n(idx_H2j)  &
        +k(117)*n(idx_O2j)  &
        +k(71)*n(idx_COj)  &
        +k(1139)  &
        +k(65)*n(idx_CNj)  &
        +k(50)*n(idx_CH4j)  &
        +k(171)*n(idx_N2j)  &
        +k(210)*n(idx_Oj)

    !d[NH3+_dot]/d[H2CO]
    pd(78,13) =  &
        +k(771)*n(idx_NH2j)

    !d[NO+_dot]/d[H2CO]
    pd(79,13) =  &
        +k(723)*n(idx_Nj)

    !d[CN+_dot]/d[H2CO]
    pd(86,13) =  &
        -k(480)*n(idx_CNj)  &
        -k(65)*n(idx_CNj)

    !d[CO+_dot]/d[H2CO]
    pd(87,13) =  &
        +k(498)*n(idx_Hj)  &
        -k(485)*n(idx_COj)  &
        +k(671)*n(idx_HEj)  &
        -k(71)*n(idx_COj)

    !d[N2+_dot]/d[H2CO]
    pd(88,13) =  &
        -k(731)*n(idx_N2j)  &
        -k(171)*n(idx_N2j)

    !d[O2+_dot]/d[H2CO]
    pd(89,13) =  &
        -k(552)*n(idx_O2j)  &
        -k(117)*n(idx_O2j)

    !d[H2O+_dot]/d[H2CO]
    pd(90,13) =  &
        -k(555)*n(idx_H2Oj)  &
        -k(118)*n(idx_H2Oj)

    !d[NH2+_dot]/d[H2CO]
    pd(91,13) =  &
        -k(771)*n(idx_NH2j)  &
        -k(770)*n(idx_NH2j)

    !d[O+_dot]/d[H2CO]
    pd(92,13) =  &
        -k(210)*n(idx_Oj)  &
        -k(814)*n(idx_Oj)

    !d[OH+_dot]/d[H2CO]
    pd(93,13) =  &
        -k(220)*n(idx_OHj)  &
        -k(840)*n(idx_OHj)

    !d[CH3+_dot]/d[H2CO]
    pd(94,13) =  &
        +k(400)*n(idx_CHj)  &
        -k(441)*n(idx_CH3j)

    !d[CH4+_dot]/d[H2CO]
    pd(95,13) =  &
        -k(450)*n(idx_CH4j)  &
        -k(50)*n(idx_CH4j)

    !d[N+_dot]/d[H2CO]
    pd(96,13) =  &
        -k(160)*n(idx_Nj)  &
        -k(722)*n(idx_Nj)  &
        -k(723)*n(idx_Nj)

    !d[HCN+_dot]/d[H2CO]
    pd(97,13) =  &
        -k(623)*n(idx_HCNj)

    !d[NH+_dot]/d[H2CO]
    pd(98,13) =  &
        -k(754)*n(idx_NHj)  &
        -k(753)*n(idx_NHj)  &
        -k(176)*n(idx_NHj)

    !d[H2+_dot]/d[H2CO]
    pd(102,13) =  &
        -k(518)*n(idx_H2j)  &
        -k(106)*n(idx_H2j)

    !d[HE+_dot]/d[H2CO]
    pd(103,13) =  &
        -k(673)*n(idx_HEj)  &
        -k(143)*n(idx_HEj)  &
        -k(672)*n(idx_HEj)  &
        -k(671)*n(idx_HEj)

    !d[HNO+_dot]/d[H2CO]
    pd(104,13) =  &
        -k(551)*n(idx_HNOj)

    !d[H3+_dot]/d[H2CO]
    pd(106,13) =  &
        -k(586)*n(idx_H3j)

    !d[H3CO+_dot]/d[H2CO]
    pd(107,13) =  &
        +k(551)*n(idx_HNOj)  &
        +k(586)*n(idx_H3j)  &
        +k(636)*n(idx_HCOj)  &
        +k(553)*n(idx_O2Hj)  &
        +k(623)*n(idx_HCNj)  &
        +k(753)*n(idx_NHj)  &
        +k(549)*n(idx_H2COj)  &
        +k(555)*n(idx_H2Oj)  &
        +k(401)*n(idx_CHj)  &
        +k(840)*n(idx_OHj)  &
        +k(634)*n(idx_HCNHj)  &
        +k(770)*n(idx_NH2j)  &
        +k(635)*n(idx_HCNHj)  &
        +k(736)*n(idx_N2Hj)  &
        +k(450)*n(idx_CH4j)  &
        +k(608)*n(idx_H3Oj)

    !d[H3O+_dot]/d[H2CO]
    pd(108,13) =  &
        -k(608)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[H2CO]
    pd(109,13) =  &
        -k(634)*n(idx_HCNHj)  &
        -k(635)*n(idx_HCNHj)

    !d[N2H+_dot]/d[H2CO]
    pd(112,13) =  &
        -k(736)*n(idx_N2Hj)

    !d[O2H+_dot]/d[H2CO]
    pd(113,13) =  &
        -k(553)*n(idx_O2Hj)

    !d[E_dot]/d[HCO]
    pd(1,14) =  &
        +k(1151)  &
        +k(262)

    !d[CH_dot]/d[HCO]
    pd(2,14) =  &
        +k(865)*n(idx_C)  &
        +k(32)*n(idx_CHj)  &
        -k(926)*n(idx_CH)

    !d[O_dot]/d[HCO]
    pd(3,14) =  &
        +k(683)*n(idx_HEj)  &
        -k(1068)*n(idx_O)  &
        +k(1016)*n(idx_N)  &
        +k(212)*n(idx_Oj)  &
        +k(979)*n(idx_H)  &
        +k(844)*n(idx_OHj)  &
        -k(1067)*n(idx_O)

    !d[HCN_dot]/d[HCO]
    pd(5,14) =  &
        +k(1016)*n(idx_N)  &
        +k(944)*n(idx_CN)

    !d[H2_dot]/d[HCO]
    pd(6,14) =  &
        +k(501)*n(idx_Hj)  &
        +2.d0*k(998)*n(idx_HCO)  &
        +k(978)*n(idx_H)  &
        +k(109)*n(idx_H2j)  &
        +k(589)*n(idx_H3j)

    !d[C_dot]/d[HCO]
    pd(7,14) =  &
        -k(865)*n(idx_C)  &
        +k(18)*n(idx_Cj)

    !d[H_dot]/d[HCO]
    pd(8,14) =  &
        +k(1150)  &
        +k(261)  &
        +k(1017)*n(idx_N)  &
        -k(979)*n(idx_H)  &
        +k(1067)*n(idx_O)  &
        +k(83)*n(idx_Hj)  &
        +k(681)*n(idx_HEj)  &
        -k(978)*n(idx_H)

    !d[H2O_dot]/d[HCO]
    pd(9,14) =  &
        +k(1097)*n(idx_OH)  &
        +k(119)*n(idx_H2Oj)

    !d[OH_dot]/d[HCO]
    pd(10,14) =  &
        +k(1068)*n(idx_O)  &
        +k(559)*n(idx_H2Oj)  &
        -k(1097)*n(idx_OH)  &
        +k(1002)*n(idx_O2)  &
        +k(222)*n(idx_OHj)

    !d[O2_dot]/d[HCO]
    pd(11,14) =  &
        +k(138)*n(idx_O2j)  &
        -k(1002)*n(idx_O2)  &
        -k(1003)*n(idx_O2)  &
        +k(646)*n(idx_O2Hj)  &
        +k(1004)*n(idx_O2H)

    !d[CH2_dot]/d[HCO]
    pd(12,14) =  &
        +k(979)*n(idx_H)  &
        -k(883)*n(idx_CH2)  &
        +k(926)*n(idx_CH)

    !d[H2CO_dot]/d[HCO]
    pd(13,14) =  &
        +k(1000)*n(idx_HNO)  &
        +k(1004)*n(idx_O2H)  &
        +k(137)*n(idx_H2COj)  &
        +2.d0*k(999)*n(idx_HCO)

    !d[HCO_dot]/d[HCO]
    pd(14,14) =  &
        -k(642)*n(idx_H2COj)  &
        -k(1001)*n(idx_NO)  &
        -k(262)  &
        -k(138)*n(idx_O2j)  &
        -k(865)*n(idx_C)  &
        -k(119)*n(idx_H2Oj)  &
        -k(724)*n(idx_Nj)  &
        -k(18)*n(idx_Cj)  &
        -k(1004)*n(idx_O2H)  &
        -k(1151)  &
        -k(139)*n(idx_SIOj)  &
        -4.d0*k(998)*n(idx_HCO)  &
        -k(637)*n(idx_HCOj)  &
        -k(775)*n(idx_NH2j)  &
        -k(212)*n(idx_Oj)  &
        -k(978)*n(idx_H)  &
        -k(682)*n(idx_HEj)  &
        -k(47)*n(idx_CH3j)  &
        -k(1003)*n(idx_O2)  &
        -k(407)*n(idx_CHj)  &
        -k(1016)*n(idx_N)  &
        -k(32)*n(idx_CHj)  &
        -k(646)*n(idx_O2Hj)  &
        -k(163)*n(idx_Nj)  &
        -k(520)*n(idx_H2j)  &
        -k(979)*n(idx_H)  &
        -k(190)*n(idx_NH3j)  &
        -k(732)*n(idx_N2j)  &
        -k(626)*n(idx_HCNj)  &
        -k(67)*n(idx_CNj)  &
        -k(926)*n(idx_CH)  &
        -k(420)*n(idx_CH2j)  &
        -4.d0*k(999)*n(idx_HCO)  &
        -k(644)*n(idx_N2Hj)  &
        -k(843)*n(idx_OHj)  &
        -k(181)*n(idx_NH2j)  &
        -k(1097)*n(idx_OH)  &
        -k(72)*n(idx_COj)  &
        -k(1262)  &
        -k(442)*n(idx_CH3j)  &
        -k(683)*n(idx_HEj)  &
        -k(1150)  &
        -k(760)*n(idx_NHj)  &
        -k(109)*n(idx_H2j)  &
        -k(261)  &
        -k(1017)*n(idx_N)  &
        -k(83)*n(idx_Hj)  &
        -k(625)*n(idx_HCNj)  &
        -k(481)*n(idx_CNj)  &
        -k(906)*n(idx_CH3)  &
        -k(645)*n(idx_O2j)  &
        -k(643)*n(idx_HNOj)  &
        -k(817)*n(idx_Oj)  &
        -k(1067)*n(idx_O)  &
        -k(844)*n(idx_OHj)  &
        -k(559)*n(idx_H2Oj)  &
        -k(137)*n(idx_H2COj)  &
        -k(1000)*n(idx_HNO)  &
        -k(883)*n(idx_CH2)  &
        -k(1068)*n(idx_O)  &
        -k(681)*n(idx_HEj)  &
        -k(222)*n(idx_OHj)  &
        -k(502)*n(idx_Hj)  &
        -k(944)*n(idx_CN)  &
        -k(1015)*n(idx_N)  &
        -k(1002)*n(idx_O2)  &
        -k(172)*n(idx_N2j)  &
        -k(589)*n(idx_H3j)  &
        -k(373)*n(idx_Cj)  &
        -k(558)*n(idx_H2Oj)  &
        -k(501)*n(idx_Hj)

    !d[NH3_dot]/d[HCO]
    pd(16,14) =  &
        +k(190)*n(idx_NH3j)

    !d[NO_dot]/d[HCO]
    pd(17,14) =  &
        +k(1000)*n(idx_HNO)  &
        -k(1001)*n(idx_NO)  &
        +k(643)*n(idx_HNOj)

    !d[CN_dot]/d[HCO]
    pd(24,14) =  &
        +k(625)*n(idx_HCNj)  &
        +k(67)*n(idx_CNj)  &
        -k(944)*n(idx_CN)

    !d[CO_dot]/d[HCO]
    pd(25,14) =  &
        +k(682)*n(idx_HEj)  &
        +k(1015)*n(idx_N)  &
        +k(502)*n(idx_Hj)  &
        +k(944)*n(idx_CN)  &
        +k(978)*n(idx_H)  &
        +k(420)*n(idx_CH2j)  &
        +k(520)*n(idx_H2j)  &
        +k(724)*n(idx_Nj)  &
        +k(1150)  &
        +k(373)*n(idx_Cj)  &
        +k(261)  &
        +k(926)*n(idx_CH)  &
        +k(642)*n(idx_H2COj)  &
        +k(732)*n(idx_N2j)  &
        +k(1068)*n(idx_O)  &
        +4.d0*k(998)*n(idx_HCO)  &
        +k(843)*n(idx_OHj)  &
        +k(637)*n(idx_HCOj)  &
        +k(558)*n(idx_H2Oj)  &
        +k(481)*n(idx_CNj)  &
        +k(72)*n(idx_COj)  &
        +k(407)*n(idx_CHj)  &
        +2.d0*k(999)*n(idx_HCO)  &
        +k(1001)*n(idx_NO)  &
        +k(626)*n(idx_HCNj)  &
        +k(1003)*n(idx_O2)  &
        +k(1097)*n(idx_OH)  &
        +k(645)*n(idx_O2j)  &
        +k(906)*n(idx_CH3)  &
        +k(865)*n(idx_C)  &
        +k(883)*n(idx_CH2)  &
        +k(442)*n(idx_CH3j)  &
        +k(817)*n(idx_Oj)

    !d[N2_dot]/d[HCO]
    pd(26,14) =  &
        +k(172)*n(idx_N2j)  &
        +k(644)*n(idx_N2Hj)

    !d[NH2_dot]/d[HCO]
    pd(27,14) =  &
        +k(181)*n(idx_NH2j)

    !d[CH3_dot]/d[HCO]
    pd(28,14) =  &
        -k(906)*n(idx_CH3)  &
        +k(883)*n(idx_CH2)  &
        +k(47)*n(idx_CH3j)

    !d[CH4_dot]/d[HCO]
    pd(29,14) =  &
        +k(906)*n(idx_CH3)

    !d[N_dot]/d[HCO]
    pd(30,14) =  &
        -k(1015)*n(idx_N)  &
        -k(1016)*n(idx_N)  &
        -k(1017)*n(idx_N)  &
        +k(760)*n(idx_NHj)  &
        +k(163)*n(idx_Nj)

    !d[NH_dot]/d[HCO]
    pd(31,14) =  &
        +k(775)*n(idx_NH2j)  &
        +k(1015)*n(idx_N)

    !d[SIO_dot]/d[HCO]
    pd(34,14) =  &
        +k(139)*n(idx_SIOj)

    !d[HE_dot]/d[HCO]
    pd(35,14) =  &
        +k(683)*n(idx_HEj)  &
        +k(681)*n(idx_HEj)

    !d[HNO_dot]/d[HCO]
    pd(36,14) =  &
        -k(1000)*n(idx_HNO)  &
        +k(1001)*n(idx_NO)

    !d[CO2_dot]/d[HCO]
    pd(38,14) =  &
        +k(1067)*n(idx_O)  &
        +k(1002)*n(idx_O2)

    !d[O2H_dot]/d[HCO]
    pd(43,14) =  &
        +k(1003)*n(idx_O2)  &
        -k(1004)*n(idx_O2H)

    !d[OCN_dot]/d[HCO]
    pd(44,14) =  &
        +k(1017)*n(idx_N)

    !d[H2CO_DUST_dot]/d[HCO]
    pd(47,14) =  &
        +k(1262)

    !d[HCO+_dot]/d[HCO]
    pd(70,14) =  &
        +k(137)*n(idx_H2COj)  &
        +k(138)*n(idx_O2j)  &
        +k(190)*n(idx_NH3j)  &
        +k(109)*n(idx_H2j)  &
        +k(18)*n(idx_Cj)  &
        +k(47)*n(idx_CH3j)  &
        +k(181)*n(idx_NH2j)  &
        +k(32)*n(idx_CHj)  &
        +k(163)*n(idx_Nj)  &
        +k(172)*n(idx_N2j)  &
        +k(262)  &
        -k(637)*n(idx_HCOj)  &
        +k(1151)  &
        +k(72)*n(idx_COj)  &
        +k(139)*n(idx_SIOj)  &
        +k(83)*n(idx_Hj)  &
        +k(119)*n(idx_H2Oj)  &
        +k(212)*n(idx_Oj)  &
        +k(222)*n(idx_OHj)  &
        +k(67)*n(idx_CNj)

    !d[H+_dot]/d[HCO]
    pd(71,14) =  &
        -k(83)*n(idx_Hj)  &
        -k(502)*n(idx_Hj)  &
        -k(501)*n(idx_Hj)

    !d[C+_dot]/d[HCO]
    pd(73,14) =  &
        -k(373)*n(idx_Cj)  &
        -k(18)*n(idx_Cj)

    !d[CH2+_dot]/d[HCO]
    pd(74,14) =  &
        -k(420)*n(idx_CH2j)  &
        +k(407)*n(idx_CHj)

    !d[CH+_dot]/d[HCO]
    pd(75,14) =  &
        +k(683)*n(idx_HEj)  &
        -k(407)*n(idx_CHj)  &
        +k(373)*n(idx_Cj)  &
        -k(32)*n(idx_CHj)

    !d[H2CO+_dot]/d[HCO]
    pd(76,14) =  &
        -k(642)*n(idx_H2COj)  &
        +k(775)*n(idx_NH2j)  &
        -k(137)*n(idx_H2COj)  &
        +k(625)*n(idx_HCNj)  &
        +k(643)*n(idx_HNOj)  &
        +k(646)*n(idx_O2Hj)  &
        +k(760)*n(idx_NHj)  &
        +k(844)*n(idx_OHj)  &
        +k(637)*n(idx_HCOj)  &
        +k(589)*n(idx_H3j)  &
        +k(644)*n(idx_N2Hj)  &
        +k(559)*n(idx_H2Oj)

    !d[NH3+_dot]/d[HCO]
    pd(78,14) =  &
        -k(190)*n(idx_NH3j)

    !d[CN+_dot]/d[HCO]
    pd(86,14) =  &
        -k(481)*n(idx_CNj)  &
        -k(67)*n(idx_CNj)

    !d[CO+_dot]/d[HCO]
    pd(87,14) =  &
        -k(72)*n(idx_COj)  &
        +k(501)*n(idx_Hj)  &
        +k(681)*n(idx_HEj)

    !d[N2+_dot]/d[HCO]
    pd(88,14) =  &
        -k(172)*n(idx_N2j)  &
        -k(732)*n(idx_N2j)

    !d[O2+_dot]/d[HCO]
    pd(89,14) =  &
        -k(138)*n(idx_O2j)  &
        -k(645)*n(idx_O2j)

    !d[H2O+_dot]/d[HCO]
    pd(90,14) =  &
        -k(119)*n(idx_H2Oj)  &
        +k(843)*n(idx_OHj)  &
        -k(558)*n(idx_H2Oj)  &
        -k(559)*n(idx_H2Oj)

    !d[NH2+_dot]/d[HCO]
    pd(91,14) =  &
        -k(775)*n(idx_NH2j)  &
        -k(181)*n(idx_NH2j)

    !d[O+_dot]/d[HCO]
    pd(92,14) =  &
        -k(817)*n(idx_Oj)  &
        -k(212)*n(idx_Oj)

    !d[OH+_dot]/d[HCO]
    pd(93,14) =  &
        -k(843)*n(idx_OHj)  &
        -k(222)*n(idx_OHj)  &
        +k(817)*n(idx_Oj)  &
        -k(844)*n(idx_OHj)

    !d[CH3+_dot]/d[HCO]
    pd(94,14) =  &
        +k(420)*n(idx_CH2j)  &
        -k(47)*n(idx_CH3j)  &
        -k(442)*n(idx_CH3j)

    !d[CH4+_dot]/d[HCO]
    pd(95,14) =  &
        +k(442)*n(idx_CH3j)

    !d[N+_dot]/d[HCO]
    pd(96,14) =  &
        -k(163)*n(idx_Nj)  &
        -k(724)*n(idx_Nj)

    !d[HCN+_dot]/d[HCO]
    pd(97,14) =  &
        -k(626)*n(idx_HCNj)  &
        +k(481)*n(idx_CNj)  &
        -k(625)*n(idx_HCNj)

    !d[NH+_dot]/d[HCO]
    pd(98,14) =  &
        +k(724)*n(idx_Nj)  &
        -k(760)*n(idx_NHj)

    !d[SIO+_dot]/d[HCO]
    pd(101,14) =  &
        -k(139)*n(idx_SIOj)

    !d[H2+_dot]/d[HCO]
    pd(102,14) =  &
        -k(109)*n(idx_H2j)  &
        +k(502)*n(idx_Hj)  &
        -k(520)*n(idx_H2j)

    !d[HE+_dot]/d[HCO]
    pd(103,14) =  &
        -k(681)*n(idx_HEj)  &
        -k(683)*n(idx_HEj)  &
        -k(682)*n(idx_HEj)

    !d[HNO+_dot]/d[HCO]
    pd(104,14) =  &
        -k(643)*n(idx_HNOj)

    !d[H3+_dot]/d[HCO]
    pd(106,14) =  &
        +k(520)*n(idx_H2j)  &
        -k(589)*n(idx_H3j)

    !d[H3CO+_dot]/d[HCO]
    pd(107,14) =  &
        +k(642)*n(idx_H2COj)

    !d[H3O+_dot]/d[HCO]
    pd(108,14) =  &
        +k(558)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[HCO]
    pd(109,14) =  &
        +k(626)*n(idx_HCNj)

    !d[HEH+_dot]/d[HCO]
    pd(111,14) =  &
        +k(682)*n(idx_HEj)

    !d[N2H+_dot]/d[HCO]
    pd(112,14) =  &
        -k(644)*n(idx_N2Hj)  &
        +k(732)*n(idx_N2j)

    !d[O2H+_dot]/d[HCO]
    pd(113,14) =  &
        -k(646)*n(idx_O2Hj)  &
        +k(645)*n(idx_O2j)

    !d[E_dot]/d[MG]
    pd(1,15) =  &
        +k(267)  &
        +k(1155)

    !d[CH_dot]/d[MG]
    pd(2,15) =  &
        +k(33)*n(idx_CHj)

    !d[H2_dot]/d[MG]
    pd(6,15) =  &
        +k(592)*n(idx_H3j)

    !d[C_dot]/d[MG]
    pd(7,15) =  &
        +k(19)*n(idx_Cj)

    !d[H_dot]/d[MG]
    pd(8,15) =  &
        +k(84)*n(idx_Hj)  &
        +k(592)*n(idx_H3j)

    !d[H2O_dot]/d[MG]
    pd(9,15) =  &
        +k(120)*n(idx_H2Oj)

    !d[O2_dot]/d[MG]
    pd(11,15) =  &
        +k(153)*n(idx_O2j)

    !d[H2CO_dot]/d[MG]
    pd(13,15) =  &
        +k(149)*n(idx_H2COj)

    !d[HCO_dot]/d[MG]
    pd(14,15) =  &
        +k(150)*n(idx_HCOj)

    !d[MG_dot]/d[MG]
    pd(15,15) =  &
        -k(151)*n(idx_N2j)  &
        -k(153)*n(idx_O2j)  &
        -k(592)*n(idx_H3j)  &
        -k(1301)  &
        -k(19)*n(idx_Cj)  &
        -k(152)*n(idx_NOj)  &
        -k(149)*n(idx_H2COj)  &
        -k(191)*n(idx_NH3j)  &
        -k(155)*n(idx_SIOj)  &
        -k(150)*n(idx_HCOj)  &
        -k(154)*n(idx_SIj)  &
        -k(120)*n(idx_H2Oj)  &
        -k(84)*n(idx_Hj)  &
        -k(1155)  &
        -k(164)*n(idx_Nj)  &
        -k(267)  &
        -k(48)*n(idx_CH3j)  &
        -k(33)*n(idx_CHj)

    !d[NH3_dot]/d[MG]
    pd(16,15) =  &
        +k(191)*n(idx_NH3j)

    !d[NO_dot]/d[MG]
    pd(17,15) =  &
        +k(152)*n(idx_NOj)

    !d[SI_dot]/d[MG]
    pd(18,15) =  &
        +k(154)*n(idx_SIj)

    !d[N2_dot]/d[MG]
    pd(26,15) =  &
        +k(151)*n(idx_N2j)

    !d[CH3_dot]/d[MG]
    pd(28,15) =  &
        +k(48)*n(idx_CH3j)

    !d[N_dot]/d[MG]
    pd(30,15) =  &
        +k(164)*n(idx_Nj)

    !d[SIO_dot]/d[MG]
    pd(34,15) =  &
        +k(155)*n(idx_SIOj)

    !d[MG_DUST_dot]/d[MG]
    pd(66,15) =  &
        +k(1301)

    !d[HCO+_dot]/d[MG]
    pd(70,15) =  &
        -k(150)*n(idx_HCOj)

    !d[H+_dot]/d[MG]
    pd(71,15) =  &
        -k(84)*n(idx_Hj)

    !d[C+_dot]/d[MG]
    pd(73,15) =  &
        -k(19)*n(idx_Cj)

    !d[CH+_dot]/d[MG]
    pd(75,15) =  &
        -k(33)*n(idx_CHj)

    !d[H2CO+_dot]/d[MG]
    pd(76,15) =  &
        -k(149)*n(idx_H2COj)

    !d[MG+_dot]/d[MG]
    pd(77,15) =  &
        +k(267)  &
        +k(151)*n(idx_N2j)  &
        +k(164)*n(idx_Nj)  &
        +k(120)*n(idx_H2Oj)  &
        +k(150)*n(idx_HCOj)  &
        +k(84)*n(idx_Hj)  &
        +k(149)*n(idx_H2COj)  &
        +k(155)*n(idx_SIOj)  &
        +k(154)*n(idx_SIj)  &
        +k(191)*n(idx_NH3j)  &
        +k(1155)  &
        +k(19)*n(idx_Cj)  &
        +k(153)*n(idx_O2j)  &
        +k(592)*n(idx_H3j)  &
        +k(33)*n(idx_CHj)  &
        +k(152)*n(idx_NOj)  &
        +k(48)*n(idx_CH3j)

    !d[NH3+_dot]/d[MG]
    pd(78,15) =  &
        -k(191)*n(idx_NH3j)

    !d[NO+_dot]/d[MG]
    pd(79,15) =  &
        -k(152)*n(idx_NOj)

    !d[SI+_dot]/d[MG]
    pd(80,15) =  &
        -k(154)*n(idx_SIj)

    !d[N2+_dot]/d[MG]
    pd(88,15) =  &
        -k(151)*n(idx_N2j)

    !d[O2+_dot]/d[MG]
    pd(89,15) =  &
        -k(153)*n(idx_O2j)

    !d[H2O+_dot]/d[MG]
    pd(90,15) =  &
        -k(120)*n(idx_H2Oj)

    !d[CH3+_dot]/d[MG]
    pd(94,15) =  &
        -k(48)*n(idx_CH3j)

    !d[N+_dot]/d[MG]
    pd(96,15) =  &
        -k(164)*n(idx_Nj)

    !d[SIO+_dot]/d[MG]
    pd(101,15) =  &
        -k(155)*n(idx_SIOj)

    !d[H3+_dot]/d[MG]
    pd(106,15) =  &
        -k(592)*n(idx_H3j)

    !d[E_dot]/d[NH3]
    pd(1,16) =  &
        +k(1161)  &
        +k(273)

    !d[CH_dot]/d[NH3]
    pd(2,16) =  &
        +k(34)*n(idx_CHj)

    !d[O_dot]/d[NH3]
    pd(3,16) =  &
        +k(214)*n(idx_Oj)  &
        -k(1075)*n(idx_O)

    !d[HCN_dot]/d[NH3]
    pd(5,16) =  &
        +k(197)*n(idx_HCNj)  &
        +k(1034)*n(idx_CN)

    !d[H2_dot]/d[NH3]
    pd(6,16) =  &
        +k(1162)  &
        +k(274)  &
        +k(985)*n(idx_H)  &
        +k(725)*n(idx_Nj)  &
        +k(111)*n(idx_H2j)  &
        +k(692)*n(idx_HEj)  &
        +k(375)*n(idx_Cj)

    !d[C_dot]/d[NH3]
    pd(7,16) =  &
        +k(20)*n(idx_Cj)

    !d[H_dot]/d[NH3]
    pd(8,16) =  &
        +k(272)  &
        +k(693)*n(idx_HEj)  &
        -k(985)*n(idx_H)  &
        +k(86)*n(idx_Hj)  &
        +k(1160)

    !d[H2O_dot]/d[NH3]
    pd(9,16) =  &
        +k(1099)*n(idx_OH)  &
        +k(196)*n(idx_H2Oj)

    !d[OH_dot]/d[NH3]
    pd(10,16) =  &
        +k(1075)*n(idx_O)  &
        -k(1099)*n(idx_OH)  &
        +k(223)*n(idx_OHj)

    !d[O2_dot]/d[NH3]
    pd(11,16) =  &
        +k(199)*n(idx_O2j)

    !d[H2CO_dot]/d[NH3]
    pd(13,16) =  &
        +k(195)*n(idx_H2COj)

    !d[NH3_dot]/d[NH3]
    pd(16,16) =  &
        -k(1160)  &
        -k(86)*n(idx_Hj)  &
        -k(34)*n(idx_CHj)  &
        -k(272)  &
        -k(274)  &
        -k(909)*n(idx_CH3)  &
        -k(1075)*n(idx_O)  &
        -k(1038)*n(idx_NH)  &
        -k(182)*n(idx_NH2j)  &
        -k(692)*n(idx_HEj)  &
        -k(20)*n(idx_Cj)  &
        -k(726)*n(idx_Nj)  &
        -k(111)*n(idx_H2j)  &
        -k(793)*n(idx_COj)  &
        -k(794)*n(idx_HCNj)  &
        -k(194)*n(idx_COj)  &
        -k(195)*n(idx_H2COj)  &
        -k(166)*n(idx_Nj)  &
        -k(1161)  &
        -k(375)*n(idx_Cj)  &
        -k(985)*n(idx_H)  &
        -k(214)*n(idx_Oj)  &
        -k(1268)  &
        -k(196)*n(idx_H2Oj)  &
        -k(273)  &
        -k(199)*n(idx_O2j)  &
        -k(197)*n(idx_HCNj)  &
        -k(1034)*n(idx_CN)  &
        -k(1099)*n(idx_OH)  &
        -k(178)*n(idx_NHj)  &
        -k(693)*n(idx_HEj)  &
        -k(223)*n(idx_OHj)  &
        -k(725)*n(idx_Nj)  &
        -k(1162)  &
        -k(51)*n(idx_CH4j)  &
        -k(198)*n(idx_N2j)  &
        -k(146)*n(idx_HEj)

    !d[CN_dot]/d[NH3]
    pd(24,16) =  &
        -k(1034)*n(idx_CN)

    !d[CO_dot]/d[NH3]
    pd(25,16) =  &
        +k(194)*n(idx_COj)

    !d[N2_dot]/d[NH3]
    pd(26,16) =  &
        +k(198)*n(idx_N2j)

    !d[NH2_dot]/d[NH3]
    pd(27,16) =  &
        +k(1099)*n(idx_OH)  &
        +k(1034)*n(idx_CN)  &
        +k(909)*n(idx_CH3)  &
        +k(794)*n(idx_HCNj)  &
        +k(1160)  &
        +k(985)*n(idx_H)  &
        +k(1075)*n(idx_O)  &
        +k(793)*n(idx_COj)  &
        +2.d0*k(1038)*n(idx_NH)  &
        +k(272)  &
        +k(182)*n(idx_NH2j)

    !d[CH3_dot]/d[NH3]
    pd(28,16) =  &
        -k(909)*n(idx_CH3)

    !d[CH4_dot]/d[NH3]
    pd(29,16) =  &
        +k(51)*n(idx_CH4j)  &
        +k(909)*n(idx_CH3)

    !d[N_dot]/d[NH3]
    pd(30,16) =  &
        +k(166)*n(idx_Nj)

    !d[NH_dot]/d[NH3]
    pd(31,16) =  &
        -k(1038)*n(idx_NH)  &
        +k(178)*n(idx_NHj)  &
        +k(274)  &
        +k(726)*n(idx_Nj)  &
        +k(1162)

    !d[HE_dot]/d[NH3]
    pd(35,16) =  &
        +k(692)*n(idx_HEj)  &
        +k(693)*n(idx_HEj)  &
        +k(146)*n(idx_HEj)

    !d[NH3_DUST_dot]/d[NH3]
    pd(60,16) =  &
        +k(1268)

    !d[HCO+_dot]/d[NH3]
    pd(70,16) =  &
        +k(793)*n(idx_COj)

    !d[H+_dot]/d[NH3]
    pd(71,16) =  &
        -k(86)*n(idx_Hj)

    !d[C+_dot]/d[NH3]
    pd(73,16) =  &
        -k(375)*n(idx_Cj)  &
        -k(20)*n(idx_Cj)

    !d[CH+_dot]/d[NH3]
    pd(75,16) =  &
        -k(34)*n(idx_CHj)

    !d[H2CO+_dot]/d[NH3]
    pd(76,16) =  &
        -k(195)*n(idx_H2COj)

    !d[NH3+_dot]/d[NH3]
    pd(78,16) =  &
        +k(214)*n(idx_Oj)  &
        +k(198)*n(idx_N2j)  &
        +k(197)*n(idx_HCNj)  &
        +k(223)*n(idx_OHj)  &
        +k(1161)  &
        +k(194)*n(idx_COj)  &
        +k(195)*n(idx_H2COj)  &
        +k(178)*n(idx_NHj)  &
        +k(51)*n(idx_CH4j)  &
        +k(182)*n(idx_NH2j)  &
        +k(146)*n(idx_HEj)  &
        +k(199)*n(idx_O2j)  &
        +k(196)*n(idx_H2Oj)  &
        +k(111)*n(idx_H2j)  &
        +k(34)*n(idx_CHj)  &
        +k(166)*n(idx_Nj)  &
        +k(20)*n(idx_Cj)  &
        +k(273)  &
        +k(86)*n(idx_Hj)

    !d[CO+_dot]/d[NH3]
    pd(87,16) =  &
        -k(793)*n(idx_COj)  &
        -k(194)*n(idx_COj)

    !d[N2+_dot]/d[NH3]
    pd(88,16) =  &
        -k(198)*n(idx_N2j)

    !d[O2+_dot]/d[NH3]
    pd(89,16) =  &
        -k(199)*n(idx_O2j)

    !d[H2O+_dot]/d[NH3]
    pd(90,16) =  &
        -k(196)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NH3]
    pd(91,16) =  &
        +k(693)*n(idx_HEj)  &
        -k(182)*n(idx_NH2j)  &
        +k(726)*n(idx_Nj)

    !d[O+_dot]/d[NH3]
    pd(92,16) =  &
        -k(214)*n(idx_Oj)

    !d[OH+_dot]/d[NH3]
    pd(93,16) =  &
        -k(223)*n(idx_OHj)

    !d[CH4+_dot]/d[NH3]
    pd(95,16) =  &
        -k(51)*n(idx_CH4j)

    !d[N+_dot]/d[NH3]
    pd(96,16) =  &
        -k(166)*n(idx_Nj)  &
        -k(726)*n(idx_Nj)  &
        -k(725)*n(idx_Nj)

    !d[HCN+_dot]/d[NH3]
    pd(97,16) =  &
        -k(197)*n(idx_HCNj)  &
        -k(794)*n(idx_HCNj)  &
        +k(375)*n(idx_Cj)

    !d[NH+_dot]/d[NH3]
    pd(98,16) =  &
        -k(178)*n(idx_NHj)  &
        +k(692)*n(idx_HEj)

    !d[H2+_dot]/d[NH3]
    pd(102,16) =  &
        -k(111)*n(idx_H2j)

    !d[HE+_dot]/d[NH3]
    pd(103,16) =  &
        -k(146)*n(idx_HEj)  &
        -k(692)*n(idx_HEj)  &
        -k(693)*n(idx_HEj)

    !d[HCNH+_dot]/d[NH3]
    pd(109,16) =  &
        +k(794)*n(idx_HCNj)

    !d[N2H+_dot]/d[NH3]
    pd(112,16) =  &
        +k(725)*n(idx_Nj)

    !d[E_dot]/d[NO]
    pd(1,17) =  &
        +k(278)  &
        +k(1166)

    !d[CH_dot]/d[NO]
    pd(2,17) =  &
        -k(932)*n(idx_CH)  &
        -k(931)*n(idx_CH)  &
        +k(35)*n(idx_CHj)  &
        -k(933)*n(idx_CH)

    !d[O_dot]/d[NO]
    pd(3,17) =  &
        +k(728)*n(idx_Nj)  &
        +k(1167)  &
        +k(765)*n(idx_NHj)  &
        +k(931)*n(idx_CH)  &
        +k(1053)*n(idx_O2)  &
        +k(1023)*n(idx_N)  &
        +k(279)  &
        +k(1043)*n(idx_NH)  &
        +k(696)*n(idx_HEj)  &
        +k(988)*n(idx_H)  &
        +k(847)*n(idx_OHj)  &
        -k(1077)*n(idx_O)  &
        +k(872)*n(idx_C)

    !d[HCN_dot]/d[NO]
    pd(5,17) =  &
        +k(931)*n(idx_CH)  &
        +k(133)*n(idx_HCNj)  &
        +k(911)*n(idx_CH3)  &
        +k(888)*n(idx_CH2)

    !d[H2_dot]/d[NO]
    pd(6,17) =  &
        +k(597)*n(idx_H3j)  &
        +k(113)*n(idx_H2j)

    !d[C_dot]/d[NO]
    pd(7,17) =  &
        -k(872)*n(idx_C)  &
        -k(873)*n(idx_C)  &
        +k(21)*n(idx_Cj)

    !d[H_dot]/d[NO]
    pd(8,17) =  &
        +k(933)*n(idx_CH)  &
        -k(989)*n(idx_H)  &
        +k(88)*n(idx_Hj)  &
        +k(1100)*n(idx_OH)  &
        +k(889)*n(idx_CH2)  &
        +k(1031)*n(idx_NH2)  &
        +k(525)*n(idx_H2j)  &
        +k(1043)*n(idx_NH)  &
        -k(988)*n(idx_H)

    !d[H2O_dot]/d[NO]
    pd(9,17) =  &
        +k(911)*n(idx_CH3)  &
        +k(121)*n(idx_H2Oj)  &
        +k(1030)*n(idx_NH2)

    !d[OH_dot]/d[NO]
    pd(10,17) =  &
        +k(888)*n(idx_CH2)  &
        +k(1044)*n(idx_NH)  &
        -k(1100)*n(idx_OH)  &
        +k(989)*n(idx_H)  &
        +k(1031)*n(idx_NH2)  &
        +k(224)*n(idx_OHj)

    !d[O2_dot]/d[NO]
    pd(11,17) =  &
        +2.d0*k(1052)*n(idx_NO)  &
        +k(808)*n(idx_O2Hj)  &
        +k(206)*n(idx_O2j)  &
        +k(1077)*n(idx_O)  &
        -k(1053)*n(idx_O2)

    !d[CH2_dot]/d[NO]
    pd(12,17) =  &
        -k(887)*n(idx_CH2)  &
        -k(889)*n(idx_CH2)  &
        -k(888)*n(idx_CH2)  &
        +k(37)*n(idx_CH2j)

    !d[H2CO_dot]/d[NO]
    pd(13,17) =  &
        +k(204)*n(idx_H2COj)  &
        +k(887)*n(idx_CH2)

    !d[HCO_dot]/d[NO]
    pd(14,17) =  &
        -k(1001)*n(idx_HCO)  &
        +k(932)*n(idx_CH)

    !d[NH3_dot]/d[NO]
    pd(16,17) =  &
        +k(192)*n(idx_NH3j)

    !d[NO_dot]/d[NO]
    pd(17,17) =  &
        -k(808)*n(idx_O2Hj)  &
        -k(179)*n(idx_NHj)  &
        -k(1001)*n(idx_HCO)  &
        -k(765)*n(idx_NHj)  &
        -k(696)*n(idx_HEj)  &
        -k(173)*n(idx_N2j)  &
        -k(847)*n(idx_OHj)  &
        -k(525)*n(idx_H2j)  &
        -k(1256)  &
        -k(1044)*n(idx_NH)  &
        -k(1077)*n(idx_O)  &
        -k(988)*n(idx_H)  &
        -k(133)*n(idx_HCNj)  &
        -k(113)*n(idx_H2j)  &
        -k(1053)*n(idx_O2)  &
        -k(1167)  &
        -k(931)*n(idx_CH)  &
        -k(947)*n(idx_CN)  &
        -k(888)*n(idx_CH2)  &
        -k(35)*n(idx_CHj)  &
        -k(889)*n(idx_CH2)  &
        -k(168)*n(idx_Nj)  &
        -k(88)*n(idx_Hj)  &
        -k(1043)*n(idx_NH)  &
        -k(1023)*n(idx_N)  &
        -k(1100)*n(idx_OH)  &
        -k(1106)*n(idx_SI)  &
        -k(278)  &
        -k(183)*n(idx_NH2j)  &
        -k(1054)*n(idx_OCN)  &
        -k(73)*n(idx_COj)  &
        -k(989)*n(idx_H)  &
        -k(911)*n(idx_CH3)  &
        -k(728)*n(idx_Nj)  &
        -4.d0*k(1052)*n(idx_NO)  &
        -k(21)*n(idx_Cj)  &
        -k(948)*n(idx_CN)  &
        -k(1030)*n(idx_NH2)  &
        -k(207)*n(idx_SIOj)  &
        -k(224)*n(idx_OHj)  &
        -k(887)*n(idx_CH2)  &
        -k(37)*n(idx_CH2j)  &
        -k(205)*n(idx_HNOj)  &
        -k(49)*n(idx_CH3j)  &
        -k(695)*n(idx_HEj)  &
        -k(68)*n(idx_CNj)  &
        -k(872)*n(idx_C)  &
        -k(1166)  &
        -k(933)*n(idx_CH)  &
        -k(597)*n(idx_H3j)  &
        -k(121)*n(idx_H2Oj)  &
        -k(192)*n(idx_NH3j)  &
        -k(1031)*n(idx_NH2)  &
        -k(932)*n(idx_CH)  &
        -k(206)*n(idx_O2j)  &
        -k(204)*n(idx_H2COj)  &
        -k(279)  &
        -k(873)*n(idx_C)

    !d[SI_dot]/d[NO]
    pd(18,17) =  &
        -k(1106)*n(idx_SI)

    !d[CN_dot]/d[NO]
    pd(24,17) =  &
        -k(947)*n(idx_CN)  &
        -k(948)*n(idx_CN)  &
        +k(872)*n(idx_C)  &
        +k(68)*n(idx_CNj)

    !d[CO_dot]/d[NO]
    pd(25,17) =  &
        +k(1001)*n(idx_HCO)  &
        +k(73)*n(idx_COj)  &
        +k(947)*n(idx_CN)  &
        +k(873)*n(idx_C)

    !d[N2_dot]/d[NO]
    pd(26,17) =  &
        +k(1044)*n(idx_NH)  &
        +2.d0*k(1052)*n(idx_NO)  &
        +k(1023)*n(idx_N)  &
        +k(947)*n(idx_CN)  &
        +k(173)*n(idx_N2j)  &
        +k(1030)*n(idx_NH2)  &
        +k(1054)*n(idx_OCN)  &
        +k(1031)*n(idx_NH2)  &
        +k(1043)*n(idx_NH)

    !d[NH2_dot]/d[NO]
    pd(27,17) =  &
        +k(183)*n(idx_NH2j)  &
        -k(1030)*n(idx_NH2)  &
        -k(1031)*n(idx_NH2)

    !d[CH3_dot]/d[NO]
    pd(28,17) =  &
        -k(911)*n(idx_CH3)  &
        +k(49)*n(idx_CH3j)

    !d[N_dot]/d[NO]
    pd(30,17) =  &
        +k(948)*n(idx_CN)  &
        +k(168)*n(idx_Nj)  &
        +k(873)*n(idx_C)  &
        +k(989)*n(idx_H)  &
        +k(695)*n(idx_HEj)  &
        +k(1077)*n(idx_O)  &
        +k(932)*n(idx_CH)  &
        +k(1167)  &
        -k(1023)*n(idx_N)  &
        +k(1106)*n(idx_SI)  &
        +k(887)*n(idx_CH2)  &
        +k(279)

    !d[NH_dot]/d[NO]
    pd(31,17) =  &
        -k(1043)*n(idx_NH)  &
        +k(988)*n(idx_H)  &
        -k(1044)*n(idx_NH)  &
        +k(179)*n(idx_NHj)

    !d[SIO_dot]/d[NO]
    pd(34,17) =  &
        +k(1106)*n(idx_SI)  &
        +k(207)*n(idx_SIOj)

    !d[HE_dot]/d[NO]
    pd(35,17) =  &
        +k(695)*n(idx_HEj)  &
        +k(696)*n(idx_HEj)

    !d[HNO_dot]/d[NO]
    pd(36,17) =  &
        +k(1001)*n(idx_HCO)  &
        +k(205)*n(idx_HNOj)

    !d[CO2_dot]/d[NO]
    pd(38,17) =  &
        +k(1054)*n(idx_OCN)

    !d[HNCO_dot]/d[NO]
    pd(41,17) =  &
        +k(889)*n(idx_CH2)

    !d[NO2_dot]/d[NO]
    pd(42,17) =  &
        +k(1053)*n(idx_O2)  &
        +k(1100)*n(idx_OH)

    !d[OCN_dot]/d[NO]
    pd(44,17) =  &
        +k(933)*n(idx_CH)  &
        -k(1054)*n(idx_OCN)  &
        +k(948)*n(idx_CN)

    !d[NO_DUST_dot]/d[NO]
    pd(56,17) =  &
        +k(1256)

    !d[H+_dot]/d[NO]
    pd(71,17) =  &
        -k(88)*n(idx_Hj)

    !d[C+_dot]/d[NO]
    pd(73,17) =  &
        -k(21)*n(idx_Cj)

    !d[CH2+_dot]/d[NO]
    pd(74,17) =  &
        -k(37)*n(idx_CH2j)

    !d[CH+_dot]/d[NO]
    pd(75,17) =  &
        -k(35)*n(idx_CHj)

    !d[H2CO+_dot]/d[NO]
    pd(76,17) =  &
        -k(204)*n(idx_H2COj)

    !d[NH3+_dot]/d[NO]
    pd(78,17) =  &
        -k(192)*n(idx_NH3j)

    !d[NO+_dot]/d[NO]
    pd(79,17) =  &
        +k(168)*n(idx_Nj)  &
        +k(207)*n(idx_SIOj)  &
        +k(192)*n(idx_NH3j)  &
        +k(73)*n(idx_COj)  &
        +k(88)*n(idx_Hj)  &
        +k(68)*n(idx_CNj)  &
        +k(204)*n(idx_H2COj)  &
        +k(179)*n(idx_NHj)  &
        +k(278)  &
        +k(35)*n(idx_CHj)  &
        +k(173)*n(idx_N2j)  &
        +k(49)*n(idx_CH3j)  &
        +k(205)*n(idx_HNOj)  &
        +k(224)*n(idx_OHj)  &
        +k(183)*n(idx_NH2j)  &
        +k(21)*n(idx_Cj)  &
        +k(121)*n(idx_H2Oj)  &
        +k(37)*n(idx_CH2j)  &
        +k(206)*n(idx_O2j)  &
        +k(1166)  &
        +k(113)*n(idx_H2j)  &
        +k(133)*n(idx_HCNj)

    !d[CN+_dot]/d[NO]
    pd(86,17) =  &
        -k(68)*n(idx_CNj)

    !d[CO+_dot]/d[NO]
    pd(87,17) =  &
        -k(73)*n(idx_COj)

    !d[N2+_dot]/d[NO]
    pd(88,17) =  &
        +k(728)*n(idx_Nj)  &
        -k(173)*n(idx_N2j)

    !d[O2+_dot]/d[NO]
    pd(89,17) =  &
        -k(206)*n(idx_O2j)

    !d[H2O+_dot]/d[NO]
    pd(90,17) =  &
        -k(121)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NO]
    pd(91,17) =  &
        -k(183)*n(idx_NH2j)

    !d[O+_dot]/d[NO]
    pd(92,17) =  &
        +k(695)*n(idx_HEj)

    !d[OH+_dot]/d[NO]
    pd(93,17) =  &
        -k(847)*n(idx_OHj)  &
        -k(224)*n(idx_OHj)

    !d[CH3+_dot]/d[NO]
    pd(94,17) =  &
        -k(49)*n(idx_CH3j)

    !d[N+_dot]/d[NO]
    pd(96,17) =  &
        -k(168)*n(idx_Nj)  &
        -k(728)*n(idx_Nj)  &
        +k(696)*n(idx_HEj)

    !d[HCN+_dot]/d[NO]
    pd(97,17) =  &
        -k(133)*n(idx_HCNj)

    !d[NH+_dot]/d[NO]
    pd(98,17) =  &
        -k(765)*n(idx_NHj)  &
        -k(179)*n(idx_NHj)

    !d[SIO+_dot]/d[NO]
    pd(101,17) =  &
        -k(207)*n(idx_SIOj)

    !d[H2+_dot]/d[NO]
    pd(102,17) =  &
        -k(525)*n(idx_H2j)  &
        -k(113)*n(idx_H2j)

    !d[HE+_dot]/d[NO]
    pd(103,17) =  &
        -k(696)*n(idx_HEj)  &
        -k(695)*n(idx_HEj)

    !d[HNO+_dot]/d[NO]
    pd(104,17) =  &
        +k(847)*n(idx_OHj)  &
        +k(597)*n(idx_H3j)  &
        +k(525)*n(idx_H2j)  &
        -k(205)*n(idx_HNOj)  &
        +k(808)*n(idx_O2Hj)

    !d[H3+_dot]/d[NO]
    pd(106,17) =  &
        -k(597)*n(idx_H3j)

    !d[N2H+_dot]/d[NO]
    pd(112,17) =  &
        +k(765)*n(idx_NHj)

    !d[O2H+_dot]/d[NO]
    pd(113,17) =  &
        -k(808)*n(idx_O2Hj)

    !d[E_dot]/d[SI]
    pd(1,18) =  &
        +k(286)  &
        +k(1177)

    !d[CH_dot]/d[SI]
    pd(2,18) =  &
        +k(36)*n(idx_CHj)

    !d[O_dot]/d[SI]
    pd(3,18) =  &
        +k(849)*n(idx_OHj)  &
        +k(1107)*n(idx_O2)  &
        -k(1214)*n(idx_O)

    !d[H2_dot]/d[SI]
    pd(6,18) =  &
        +k(602)*n(idx_H3j)

    !d[C_dot]/d[SI]
    pd(7,18) =  &
        +k(22)*n(idx_Cj)  &
        +k(1105)*n(idx_CO)

    !d[H_dot]/d[SI]
    pd(8,18) =  &
        +k(1103)*n(idx_OH)  &
        +k(92)*n(idx_Hj)

    !d[H2O_dot]/d[SI]
    pd(9,18) =  &
        +k(611)*n(idx_H3Oj)  &
        +k(123)*n(idx_H2Oj)

    !d[OH_dot]/d[SI]
    pd(10,18) =  &
        -k(1103)*n(idx_OH)

    !d[O2_dot]/d[SI]
    pd(11,18) =  &
        -k(1107)*n(idx_O2)  &
        +k(231)*n(idx_O2j)

    !d[H2CO_dot]/d[SI]
    pd(13,18) =  &
        +k(229)*n(idx_H2COj)

    !d[NH3_dot]/d[SI]
    pd(16,18) =  &
        +k(193)*n(idx_NH3j)

    !d[NO_dot]/d[SI]
    pd(17,18) =  &
        -k(1106)*n(idx_NO)  &
        +k(230)*n(idx_NOj)

    !d[SI_dot]/d[SI]
    pd(18,18) =  &
        -k(1107)*n(idx_O2)  &
        -k(1214)*n(idx_O)  &
        -k(1106)*n(idx_NO)  &
        -k(1229)  &
        -k(849)*n(idx_OHj)  &
        -k(36)*n(idx_CHj)  &
        -k(1177)  &
        -k(1105)*n(idx_CO)  &
        -k(862)*n(idx_HCOj)  &
        -k(148)*n(idx_HEj)  &
        -k(602)*n(idx_H3j)  &
        -k(229)*n(idx_H2COj)  &
        -k(286)  &
        -k(611)*n(idx_H3Oj)  &
        -k(1104)*n(idx_CO2)  &
        -k(92)*n(idx_Hj)  &
        -k(231)*n(idx_O2j)  &
        -k(123)*n(idx_H2Oj)  &
        -k(193)*n(idx_NH3j)  &
        -k(1103)*n(idx_OH)  &
        -k(22)*n(idx_Cj)  &
        -k(230)*n(idx_NOj)

    !d[CO_dot]/d[SI]
    pd(25,18) =  &
        +k(1104)*n(idx_CO2)  &
        +k(862)*n(idx_HCOj)  &
        -k(1105)*n(idx_CO)

    !d[N_dot]/d[SI]
    pd(30,18) =  &
        +k(1106)*n(idx_NO)

    !d[SIO_dot]/d[SI]
    pd(34,18) =  &
        +k(1214)*n(idx_O)  &
        +k(1103)*n(idx_OH)  &
        +k(1104)*n(idx_CO2)  &
        +k(1107)*n(idx_O2)  &
        +k(1106)*n(idx_NO)  &
        +k(1105)*n(idx_CO)

    !d[HE_dot]/d[SI]
    pd(35,18) =  &
        +k(148)*n(idx_HEj)

    !d[CO2_dot]/d[SI]
    pd(38,18) =  &
        -k(1104)*n(idx_CO2)

    !d[SIH4_DUST_dot]/d[SI]
    pd(48,18) =  &
        +k(1229)

    !d[HCO+_dot]/d[SI]
    pd(70,18) =  &
        -k(862)*n(idx_HCOj)

    !d[H+_dot]/d[SI]
    pd(71,18) =  &
        -k(92)*n(idx_Hj)

    !d[C+_dot]/d[SI]
    pd(73,18) =  &
        -k(22)*n(idx_Cj)

    !d[CH+_dot]/d[SI]
    pd(75,18) =  &
        -k(36)*n(idx_CHj)

    !d[H2CO+_dot]/d[SI]
    pd(76,18) =  &
        -k(229)*n(idx_H2COj)

    !d[NH3+_dot]/d[SI]
    pd(78,18) =  &
        -k(193)*n(idx_NH3j)

    !d[NO+_dot]/d[SI]
    pd(79,18) =  &
        -k(230)*n(idx_NOj)

    !d[SI+_dot]/d[SI]
    pd(80,18) =  &
        +k(148)*n(idx_HEj)  &
        +k(286)  &
        +k(123)*n(idx_H2Oj)  &
        +k(92)*n(idx_Hj)  &
        +k(231)*n(idx_O2j)  &
        +k(230)*n(idx_NOj)  &
        +k(1177)  &
        +k(36)*n(idx_CHj)  &
        +k(22)*n(idx_Cj)  &
        +k(229)*n(idx_H2COj)  &
        +k(193)*n(idx_NH3j)

    !d[O2+_dot]/d[SI]
    pd(89,18) =  &
        -k(231)*n(idx_O2j)

    !d[H2O+_dot]/d[SI]
    pd(90,18) =  &
        -k(123)*n(idx_H2Oj)

    !d[OH+_dot]/d[SI]
    pd(93,18) =  &
        -k(849)*n(idx_OHj)

    !d[SIH+_dot]/d[SI]
    pd(100,18) =  &
        +k(849)*n(idx_OHj)  &
        +k(611)*n(idx_H3Oj)  &
        +k(862)*n(idx_HCOj)  &
        +k(602)*n(idx_H3j)

    !d[HE+_dot]/d[SI]
    pd(103,18) =  &
        -k(148)*n(idx_HEj)

    !d[H3+_dot]/d[SI]
    pd(106,18) =  &
        -k(602)*n(idx_H3j)

    !d[H3O+_dot]/d[SI]
    pd(108,18) =  &
        -k(611)*n(idx_H3Oj)

    !d[O_dot]/d[SIC2]
    pd(3,19) =  &
        -k(1082)*n(idx_O)

    !d[C_dot]/d[SIC2]
    pd(7,19) =  &
        +k(23)*n(idx_Cj)  &
        +k(287)

    !d[H_dot]/d[SIC2]
    pd(8,19) =  &
        +k(93)*n(idx_Hj)

    !d[SIC2_dot]/d[SIC2]
    pd(19,19) =  &
        -k(1241)  &
        -k(93)*n(idx_Hj)  &
        -k(287)  &
        -k(23)*n(idx_Cj)  &
        -k(1082)*n(idx_O)

    !d[SIC_dot]/d[SIC2]
    pd(21,19) =  &
        +k(287)  &
        +k(1082)*n(idx_O)

    !d[CO_dot]/d[SIC2]
    pd(25,19) =  &
        +k(1082)*n(idx_O)

    !d[SIC2_DUST_dot]/d[SIC2]
    pd(51,19) =  &
        +k(1241)

    !d[H+_dot]/d[SIC2]
    pd(71,19) =  &
        -k(93)*n(idx_Hj)

    !d[C+_dot]/d[SIC2]
    pd(73,19) =  &
        -k(23)*n(idx_Cj)

    !d[SIC2+_dot]/d[SIC2]
    pd(81,19) =  &
        +k(23)*n(idx_Cj)  &
        +k(93)*n(idx_Hj)

    !d[O_dot]/d[SIC3]
    pd(3,20) =  &
        -k(1083)*n(idx_O)

    !d[C_dot]/d[SIC3]
    pd(7,20) =  &
        +k(1178)  &
        +k(701)*n(idx_HEj)  &
        +k(24)*n(idx_Cj)  &
        +k(288)

    !d[H_dot]/d[SIC3]
    pd(8,20) =  &
        +k(94)*n(idx_Hj)

    !d[SIC2_dot]/d[SIC3]
    pd(19,20) =  &
        +k(1083)*n(idx_O)  &
        +k(1178)  &
        +k(288)

    !d[SIC3_dot]/d[SIC3]
    pd(20,20) =  &
        -k(1244)  &
        -k(94)*n(idx_Hj)  &
        -k(288)  &
        -k(1178)  &
        -k(24)*n(idx_Cj)  &
        -k(1083)*n(idx_O)  &
        -k(701)*n(idx_HEj)

    !d[CO_dot]/d[SIC3]
    pd(25,20) =  &
        +k(1083)*n(idx_O)

    !d[HE_dot]/d[SIC3]
    pd(35,20) =  &
        +k(701)*n(idx_HEj)

    !d[SIC3_DUST_dot]/d[SIC3]
    pd(52,20) =  &
        +k(1244)

    !d[H+_dot]/d[SIC3]
    pd(71,20) =  &
        -k(94)*n(idx_Hj)

    !d[C+_dot]/d[SIC3]
    pd(73,20) =  &
        -k(24)*n(idx_Cj)

    !d[SIC2+_dot]/d[SIC3]
    pd(81,20) =  &
        +k(701)*n(idx_HEj)

    !d[SIC3+_dot]/d[SIC3]
    pd(82,20) =  &
        +k(24)*n(idx_Cj)  &
        +k(94)*n(idx_Hj)

    !d[HE+_dot]/d[SIC3]
    pd(103,20) =  &
        -k(701)*n(idx_HEj)

    !d[O_dot]/d[SIC]
    pd(3,21) =  &
        -k(1085)*n(idx_O)  &
        -k(1084)*n(idx_O)

    !d[C_dot]/d[SIC]
    pd(7,21) =  &
        +k(702)*n(idx_HEj)  &
        +k(25)*n(idx_Cj)  &
        +k(1179)  &
        +k(289)  &
        +k(1085)*n(idx_O)

    !d[H_dot]/d[SIC]
    pd(8,21) =  &
        +k(95)*n(idx_Hj)

    !d[SI_dot]/d[SIC]
    pd(18,21) =  &
        +k(703)*n(idx_HEj)  &
        +k(1028)*n(idx_N)  &
        +k(1179)  &
        +k(289)  &
        +k(1084)*n(idx_O)

    !d[SIC_dot]/d[SIC]
    pd(21,21) =  &
        -k(1179)  &
        -k(95)*n(idx_Hj)  &
        -k(1084)*n(idx_O)  &
        -k(702)*n(idx_HEj)  &
        -k(289)  &
        -k(703)*n(idx_HEj)  &
        -k(1085)*n(idx_O)  &
        -k(1028)*n(idx_N)  &
        -k(1240)  &
        -k(25)*n(idx_Cj)

    !d[CN_dot]/d[SIC]
    pd(24,21) =  &
        +k(1028)*n(idx_N)

    !d[CO_dot]/d[SIC]
    pd(25,21) =  &
        +k(1084)*n(idx_O)

    !d[N_dot]/d[SIC]
    pd(30,21) =  &
        -k(1028)*n(idx_N)

    !d[SIO_dot]/d[SIC]
    pd(34,21) =  &
        +k(1085)*n(idx_O)

    !d[HE_dot]/d[SIC]
    pd(35,21) =  &
        +k(702)*n(idx_HEj)  &
        +k(703)*n(idx_HEj)

    !d[SIC_DUST_dot]/d[SIC]
    pd(50,21) =  &
        +k(1240)

    !d[H+_dot]/d[SIC]
    pd(71,21) =  &
        -k(95)*n(idx_Hj)

    !d[C+_dot]/d[SIC]
    pd(73,21) =  &
        +k(703)*n(idx_HEj)  &
        -k(25)*n(idx_Cj)

    !d[SI+_dot]/d[SIC]
    pd(80,21) =  &
        +k(702)*n(idx_HEj)

    !d[SIC+_dot]/d[SIC]
    pd(83,21) =  &
        +k(25)*n(idx_Cj)  &
        +k(95)*n(idx_Hj)

    !d[HE+_dot]/d[SIC]
    pd(103,21) =  &
        -k(702)*n(idx_HEj)  &
        -k(703)*n(idx_HEj)

    !d[E_dot]/d[SIH2]
    pd(1,22) =  &
        +k(1181)

    !d[O_dot]/d[SIH2]
    pd(3,22) =  &
        -k(1086)*n(idx_O)  &
        -k(1087)*n(idx_O)

    !d[H2_dot]/d[SIH2]
    pd(6,22) =  &
        +k(381)*n(idx_Cj)  &
        +k(506)*n(idx_Hj)  &
        +k(704)*n(idx_HEj)  &
        +k(603)*n(idx_H3j)  &
        +k(1086)*n(idx_O)

    !d[C_dot]/d[SIH2]
    pd(7,22) =  &
        +k(26)*n(idx_Cj)

    !d[H_dot]/d[SIH2]
    pd(8,22) =  &
        +2.d0*k(1087)*n(idx_O)  &
        +k(96)*n(idx_Hj)  &
        +k(705)*n(idx_HEj)  &
        +k(1182)  &
        +k(290)

    !d[H2O_dot]/d[SIH2]
    pd(9,22) =  &
        +k(612)*n(idx_H3Oj)

    !d[SIH2_dot]/d[SIH2]
    pd(22,22) =  &
        -k(381)*n(idx_Cj)  &
        -k(603)*n(idx_H3j)  &
        -k(96)*n(idx_Hj)  &
        -k(1086)*n(idx_O)  &
        -k(290)  &
        -k(612)*n(idx_H3Oj)  &
        -k(1182)  &
        -k(704)*n(idx_HEj)  &
        -k(705)*n(idx_HEj)  &
        -k(506)*n(idx_Hj)  &
        -k(1181)  &
        -k(638)*n(idx_HCOj)  &
        -k(1087)*n(idx_O)  &
        -k(1234)  &
        -k(26)*n(idx_Cj)

    !d[CO_dot]/d[SIH2]
    pd(25,22) =  &
        +k(638)*n(idx_HCOj)

    !d[SIH_dot]/d[SIH2]
    pd(33,22) =  &
        +k(1182)  &
        +k(290)

    !d[SIO_dot]/d[SIH2]
    pd(34,22) =  &
        +k(1086)*n(idx_O)  &
        +k(1087)*n(idx_O)

    !d[HE_dot]/d[SIH2]
    pd(35,22) =  &
        +k(704)*n(idx_HEj)  &
        +k(705)*n(idx_HEj)

    !d[SIH4_DUST_dot]/d[SIH2]
    pd(48,22) =  &
        +k(1234)

    !d[HCO+_dot]/d[SIH2]
    pd(70,22) =  &
        -k(638)*n(idx_HCOj)

    !d[H+_dot]/d[SIH2]
    pd(71,22) =  &
        -k(96)*n(idx_Hj)  &
        -k(506)*n(idx_Hj)

    !d[C+_dot]/d[SIH2]
    pd(73,22) =  &
        -k(26)*n(idx_Cj)  &
        -k(381)*n(idx_Cj)

    !d[SI+_dot]/d[SIH2]
    pd(80,22) =  &
        +k(704)*n(idx_HEj)

    !d[SIC+_dot]/d[SIH2]
    pd(83,22) =  &
        +k(381)*n(idx_Cj)

    !d[SIH2+_dot]/d[SIH2]
    pd(84,22) =  &
        +k(1181)  &
        +k(96)*n(idx_Hj)  &
        +k(26)*n(idx_Cj)

    !d[SIH3+_dot]/d[SIH2]
    pd(85,22) =  &
        +k(638)*n(idx_HCOj)  &
        +k(612)*n(idx_H3Oj)  &
        +k(603)*n(idx_H3j)

    !d[SIH+_dot]/d[SIH2]
    pd(100,22) =  &
        +k(506)*n(idx_Hj)  &
        +k(705)*n(idx_HEj)

    !d[HE+_dot]/d[SIH2]
    pd(103,22) =  &
        -k(704)*n(idx_HEj)  &
        -k(705)*n(idx_HEj)

    !d[H3+_dot]/d[SIH2]
    pd(106,22) =  &
        -k(603)*n(idx_H3j)

    !d[H3O+_dot]/d[SIH2]
    pd(108,22) =  &
        -k(612)*n(idx_H3Oj)

    !d[E_dot]/d[SIH3]
    pd(1,23) =  &
        +k(1184)

    !d[O_dot]/d[SIH3]
    pd(3,23) =  &
        -k(1088)*n(idx_O)

    !d[H2_dot]/d[SIH3]
    pd(6,23) =  &
        +k(604)*n(idx_H3j)  &
        +k(706)*n(idx_HEj)  &
        +k(1185)  &
        +k(507)*n(idx_Hj)

    !d[C_dot]/d[SIH3]
    pd(7,23) =  &
        +k(27)*n(idx_Cj)

    !d[H_dot]/d[SIH3]
    pd(8,23) =  &
        +k(97)*n(idx_Hj)  &
        +k(291)  &
        +k(1088)*n(idx_O)  &
        +k(707)*n(idx_HEj)  &
        +k(1183)

    !d[SIH2_dot]/d[SIH3]
    pd(22,23) =  &
        +k(291)  &
        +k(1183)

    !d[SIH3_dot]/d[SIH3]
    pd(23,23) =  &
        -k(507)*n(idx_Hj)  &
        -k(707)*n(idx_HEj)  &
        -k(1184)  &
        -k(1185)  &
        -k(706)*n(idx_HEj)  &
        -k(291)  &
        -k(604)*n(idx_H3j)  &
        -k(1088)*n(idx_O)  &
        -k(1183)  &
        -k(97)*n(idx_Hj)  &
        -k(27)*n(idx_Cj)  &
        -k(1236)

    !d[SIH_dot]/d[SIH3]
    pd(33,23) =  &
        +k(1185)

    !d[HE_dot]/d[SIH3]
    pd(35,23) =  &
        +k(706)*n(idx_HEj)  &
        +k(707)*n(idx_HEj)

    !d[H2SIO_dot]/d[SIH3]
    pd(40,23) =  &
        +k(1088)*n(idx_O)

    !d[SIH4_DUST_dot]/d[SIH3]
    pd(48,23) =  &
        +k(1236)

    !d[H+_dot]/d[SIH3]
    pd(71,23) =  &
        -k(507)*n(idx_Hj)  &
        -k(97)*n(idx_Hj)

    !d[C+_dot]/d[SIH3]
    pd(73,23) =  &
        -k(27)*n(idx_Cj)

    !d[SIH2+_dot]/d[SIH3]
    pd(84,23) =  &
        +k(707)*n(idx_HEj)  &
        +k(507)*n(idx_Hj)

    !d[SIH3+_dot]/d[SIH3]
    pd(85,23) =  &
        +k(97)*n(idx_Hj)  &
        +k(27)*n(idx_Cj)  &
        +k(1184)

    !d[SIH4+_dot]/d[SIH3]
    pd(99,23) =  &
        +k(604)*n(idx_H3j)

    !d[SIH+_dot]/d[SIH3]
    pd(100,23) =  &
        +k(706)*n(idx_HEj)

    !d[HE+_dot]/d[SIH3]
    pd(103,23) =  &
        -k(706)*n(idx_HEj)  &
        -k(707)*n(idx_HEj)

    !d[H3+_dot]/d[SIH3]
    pd(106,23) =  &
        -k(604)*n(idx_H3j)

    !d[CH_dot]/d[CN]
    pd(2,24) =  &
        +k(881)*n(idx_CH2)

    !d[O_dot]/d[CN]
    pd(3,24) =  &
        +k(950)*n(idx_O2)  &
        +k(1091)*n(idx_OH)  &
        +k(837)*n(idx_OHj)  &
        -k(1058)*n(idx_O)  &
        -k(1059)*n(idx_O)

    !d[HCN_dot]/d[CN]
    pd(5,24) =  &
        +k(944)*n(idx_HCO)  &
        +k(1034)*n(idx_NH3)  &
        +k(943)*n(idx_H2CO)  &
        +k(1091)*n(idx_OH)  &
        +k(945)*n(idx_HNO)  &
        +k(1036)*n(idx_NH)  &
        +k(881)*n(idx_CH2)  &
        +k(903)*n(idx_CH3)  &
        +k(921)*n(idx_CH4)  &
        +k(951)*n(idx_SIH4)  &
        +k(960)*n(idx_H2)

    !d[H2_dot]/d[CN]
    pd(6,24) =  &
        +k(104)*n(idx_H2j)  &
        +k(582)*n(idx_H3j)  &
        -k(960)*n(idx_H2)

    !d[C_dot]/d[CN]
    pd(7,24) =  &
        +k(812)*n(idx_Oj)  &
        +k(1059)*n(idx_O)  &
        +k(1131)  &
        +k(252)  &
        +k(664)*n(idx_HEj)  &
        +k(1012)*n(idx_N)

    !d[H_dot]/d[CN]
    pd(8,24) =  &
        +k(960)*n(idx_H2)  &
        +k(514)*n(idx_H2j)  &
        +k(1092)*n(idx_OH)

    !d[OH_dot]/d[CN]
    pd(10,24) =  &
        -k(1091)*n(idx_OH)  &
        -k(1092)*n(idx_OH)

    !d[O2_dot]/d[CN]
    pd(11,24) =  &
        +k(484)*n(idx_O2Hj)  &
        -k(949)*n(idx_O2)  &
        -k(950)*n(idx_O2)

    !d[CH2_dot]/d[CN]
    pd(12,24) =  &
        -k(881)*n(idx_CH2)  &
        +k(903)*n(idx_CH3)

    !d[H2CO_dot]/d[CN]
    pd(13,24) =  &
        -k(943)*n(idx_H2CO)

    !d[HCO_dot]/d[CN]
    pd(14,24) =  &
        -k(944)*n(idx_HCO)  &
        +k(943)*n(idx_H2CO)

    !d[NH3_dot]/d[CN]
    pd(16,24) =  &
        -k(1034)*n(idx_NH3)

    !d[NO_dot]/d[CN]
    pd(17,24) =  &
        +k(946)*n(idx_NO2)  &
        -k(947)*n(idx_NO)  &
        +k(945)*n(idx_HNO)  &
        +k(949)*n(idx_O2)  &
        +k(483)*n(idx_HNOj)  &
        +k(1059)*n(idx_O)  &
        -k(948)*n(idx_NO)

    !d[SIH3_dot]/d[CN]
    pd(23,24) =  &
        +k(951)*n(idx_SIH4)

    !d[CN_dot]/d[CN]
    pd(24,24) =  &
        -k(748)*n(idx_NHj)  &
        -k(1034)*n(idx_NH3)  &
        -k(950)*n(idx_O2)  &
        -k(664)*n(idx_HEj)  &
        -k(837)*n(idx_OHj)  &
        -k(949)*n(idx_O2)  &
        -k(903)*n(idx_CH3)  &
        -k(484)*n(idx_O2Hj)  &
        -k(960)*n(idx_H2)  &
        -k(944)*n(idx_HCO)  &
        -k(943)*n(idx_H2CO)  &
        -k(921)*n(idx_CH4)  &
        -k(947)*n(idx_NO)  &
        -k(1036)*n(idx_NH)  &
        -k(951)*n(idx_SIH4)  &
        -k(812)*n(idx_Oj)  &
        -k(665)*n(idx_HEj)  &
        -k(1091)*n(idx_OH)  &
        -k(881)*n(idx_CH2)  &
        -k(1012)*n(idx_N)  &
        -k(70)*n(idx_N2j)  &
        -k(252)  &
        -k(483)*n(idx_HNOj)  &
        -k(1058)*n(idx_O)  &
        -k(1059)*n(idx_O)  &
        -k(1264)  &
        -k(1092)*n(idx_OH)  &
        -k(948)*n(idx_NO)  &
        -k(158)*n(idx_Nj)  &
        -k(945)*n(idx_HNO)  &
        -k(104)*n(idx_H2j)  &
        -k(1131)  &
        -k(582)*n(idx_H3j)  &
        -k(514)*n(idx_H2j)  &
        -k(946)*n(idx_NO2)

    !d[CO_dot]/d[CN]
    pd(25,24) =  &
        +k(1058)*n(idx_O)  &
        +k(944)*n(idx_HCO)  &
        +k(947)*n(idx_NO)  &
        +k(949)*n(idx_O2)

    !d[N2_dot]/d[CN]
    pd(26,24) =  &
        +k(70)*n(idx_N2j)  &
        +k(947)*n(idx_NO)  &
        +k(1012)*n(idx_N)

    !d[NH2_dot]/d[CN]
    pd(27,24) =  &
        +k(1034)*n(idx_NH3)

    !d[CH3_dot]/d[CN]
    pd(28,24) =  &
        -k(903)*n(idx_CH3)  &
        +k(921)*n(idx_CH4)

    !d[CH4_dot]/d[CN]
    pd(29,24) =  &
        -k(921)*n(idx_CH4)

    !d[N_dot]/d[CN]
    pd(30,24) =  &
        +k(158)*n(idx_Nj)  &
        +k(1058)*n(idx_O)  &
        +k(748)*n(idx_NHj)  &
        +k(1036)*n(idx_NH)  &
        +k(1131)  &
        +k(948)*n(idx_NO)  &
        +k(252)  &
        +k(665)*n(idx_HEj)  &
        -k(1012)*n(idx_N)

    !d[NH_dot]/d[CN]
    pd(31,24) =  &
        -k(1036)*n(idx_NH)

    !d[SIH4_dot]/d[CN]
    pd(32,24) =  &
        -k(951)*n(idx_SIH4)

    !d[HE_dot]/d[CN]
    pd(35,24) =  &
        +k(665)*n(idx_HEj)  &
        +k(664)*n(idx_HEj)

    !d[HNO_dot]/d[CN]
    pd(36,24) =  &
        -k(945)*n(idx_HNO)

    !d[NO2_dot]/d[CN]
    pd(42,24) =  &
        -k(946)*n(idx_NO2)

    !d[OCN_dot]/d[CN]
    pd(44,24) =  &
        +k(950)*n(idx_O2)  &
        +k(948)*n(idx_NO)  &
        +k(1092)*n(idx_OH)  &
        +k(946)*n(idx_NO2)

    !d[HCN_DUST_dot]/d[CN]
    pd(59,24) =  &
        +k(1264)

    !d[C+_dot]/d[CN]
    pd(73,24) =  &
        +k(665)*n(idx_HEj)

    !d[NO+_dot]/d[CN]
    pd(79,24) =  &
        +k(812)*n(idx_Oj)

    !d[CN+_dot]/d[CN]
    pd(86,24) =  &
        +k(70)*n(idx_N2j)  &
        +k(104)*n(idx_H2j)  &
        +k(158)*n(idx_Nj)

    !d[N2+_dot]/d[CN]
    pd(88,24) =  &
        -k(70)*n(idx_N2j)

    !d[O+_dot]/d[CN]
    pd(92,24) =  &
        -k(812)*n(idx_Oj)

    !d[OH+_dot]/d[CN]
    pd(93,24) =  &
        -k(837)*n(idx_OHj)

    !d[N+_dot]/d[CN]
    pd(96,24) =  &
        +k(664)*n(idx_HEj)  &
        -k(158)*n(idx_Nj)

    !d[HCN+_dot]/d[CN]
    pd(97,24) =  &
        +k(484)*n(idx_O2Hj)  &
        +k(748)*n(idx_NHj)  &
        +k(483)*n(idx_HNOj)  &
        +k(582)*n(idx_H3j)  &
        +k(514)*n(idx_H2j)  &
        +k(837)*n(idx_OHj)

    !d[NH+_dot]/d[CN]
    pd(98,24) =  &
        -k(748)*n(idx_NHj)

    !d[H2+_dot]/d[CN]
    pd(102,24) =  &
        -k(104)*n(idx_H2j)  &
        -k(514)*n(idx_H2j)

    !d[HE+_dot]/d[CN]
    pd(103,24) =  &
        -k(665)*n(idx_HEj)  &
        -k(664)*n(idx_HEj)

    !d[HNO+_dot]/d[CN]
    pd(104,24) =  &
        -k(483)*n(idx_HNOj)

    !d[H3+_dot]/d[CN]
    pd(106,24) =  &
        -k(582)*n(idx_H3j)

    !d[O2H+_dot]/d[CN]
    pd(113,24) =  &
        -k(484)*n(idx_O2Hj)

    !d[E_dot]/d[CO]
    pd(1,25) =  &
        +k(233)

    !d[O_dot]/d[CO]
    pd(3,25) =  &
        +k(209)*n(idx_Oj)  &
        +k(839)*n(idx_OHj)  &
        +k(670)*n(idx_HEj)  &
        +k(954)*n(idx_O2)  &
        +k(1134)  &
        +k(254)

    !d[H2_dot]/d[CO]
    pd(6,25) =  &
        +k(105)*n(idx_H2j)  &
        +k(585)*n(idx_H3j)  &
        +k(584)*n(idx_H3j)

    !d[C_dot]/d[CO]
    pd(7,25) =  &
        +k(973)*n(idx_H)  &
        +k(721)*n(idx_Nj)  &
        +k(1105)*n(idx_SI)  &
        +k(1134)  &
        +k(254)

    !d[H_dot]/d[CO]
    pd(8,25) =  &
        -k(973)*n(idx_H)  &
        +k(1093)*n(idx_OH)  &
        +k(516)*n(idx_H2j)

    !d[OH_dot]/d[CO]
    pd(10,25) =  &
        +k(554)*n(idx_H2Oj)  &
        -k(1093)*n(idx_OH)  &
        +k(955)*n(idx_O2H)  &
        +k(973)*n(idx_H)

    !d[O2_dot]/d[CO]
    pd(11,25) =  &
        +k(489)*n(idx_O2Hj)  &
        -k(954)*n(idx_O2)

    !d[NO_dot]/d[CO]
    pd(17,25) =  &
        +k(953)*n(idx_NO2)  &
        +k(487)*n(idx_HNOj)

    !d[SI_dot]/d[CO]
    pd(18,25) =  &
        -k(1105)*n(idx_SI)

    !d[SIH3_dot]/d[CO]
    pd(23,25) =  &
        +k(490)*n(idx_SIH4j)

    !d[CN_dot]/d[CO]
    pd(24,25) =  &
        +k(64)*n(idx_CNj)  &
        +k(622)*n(idx_HCNj)

    !d[CO_dot]/d[CO]
    pd(25,25) =  &
        -k(64)*n(idx_CNj)  &
        -k(105)*n(idx_H2j)  &
        -k(839)*n(idx_OHj)  &
        -k(953)*n(idx_NO2)  &
        -k(488)*n(idx_N2Hj)  &
        -k(622)*n(idx_HCNj)  &
        -k(721)*n(idx_Nj)  &
        -k(554)*n(idx_H2Oj)  &
        -k(973)*n(idx_H)  &
        -k(209)*n(idx_Oj)  &
        -k(1105)*n(idx_SI)  &
        -k(1093)*n(idx_OH)  &
        -k(75)*n(idx_N2j)  &
        -k(159)*n(idx_Nj)  &
        -k(954)*n(idx_O2)  &
        -k(491)*n(idx_SIOj)  &
        -k(486)*n(idx_HCO2j)  &
        -k(585)*n(idx_H3j)  &
        -k(1252)  &
        -k(1228)  &
        -k(955)*n(idx_O2H)  &
        -k(670)*n(idx_HEj)  &
        -k(1305)  &
        -k(516)*n(idx_H2j)  &
        -k(752)*n(idx_NHj)  &
        -k(487)*n(idx_HNOj)  &
        -k(584)*n(idx_H3j)  &
        -k(1134)  &
        -k(952)*n(idx_HNO)  &
        -k(233)  &
        -k(254)  &
        -k(489)*n(idx_O2Hj)  &
        -k(490)*n(idx_SIH4j)  &
        -k(449)*n(idx_CH4j)

    !d[N2_dot]/d[CO]
    pd(26,25) =  &
        +k(75)*n(idx_N2j)  &
        +k(488)*n(idx_N2Hj)

    !d[CH3_dot]/d[CO]
    pd(28,25) =  &
        +k(449)*n(idx_CH4j)

    !d[N_dot]/d[CO]
    pd(30,25) =  &
        +k(159)*n(idx_Nj)  &
        +k(752)*n(idx_NHj)

    !d[NH_dot]/d[CO]
    pd(31,25) =  &
        +k(952)*n(idx_HNO)

    !d[SIO_dot]/d[CO]
    pd(34,25) =  &
        +k(1105)*n(idx_SI)

    !d[HE_dot]/d[CO]
    pd(35,25) =  &
        +k(670)*n(idx_HEj)

    !d[HNO_dot]/d[CO]
    pd(36,25) =  &
        -k(952)*n(idx_HNO)

    !d[CO2_dot]/d[CO]
    pd(38,25) =  &
        +k(491)*n(idx_SIOj)  &
        +k(952)*n(idx_HNO)  &
        +k(955)*n(idx_O2H)  &
        +k(486)*n(idx_HCO2j)  &
        +k(1093)*n(idx_OH)  &
        +k(953)*n(idx_NO2)  &
        +k(954)*n(idx_O2)

    !d[NO2_dot]/d[CO]
    pd(42,25) =  &
        -k(953)*n(idx_NO2)

    !d[O2H_dot]/d[CO]
    pd(43,25) =  &
        -k(955)*n(idx_O2H)

    !d[CH3OH_DUST_dot]/d[CO]
    pd(45,25) =  &
        +k(1305)

    !d[H2CO_DUST_dot]/d[CO]
    pd(47,25) =  &
        +k(1228)

    !d[CO_DUST_dot]/d[CO]
    pd(54,25) =  &
        +k(1252)

    !d[HCO+_dot]/d[CO]
    pd(70,25) =  &
        +k(584)*n(idx_H3j)  &
        +k(488)*n(idx_N2Hj)  &
        +k(554)*n(idx_H2Oj)  &
        +k(489)*n(idx_O2Hj)  &
        +k(487)*n(idx_HNOj)  &
        +k(486)*n(idx_HCO2j)  &
        +k(516)*n(idx_H2j)  &
        +k(490)*n(idx_SIH4j)  &
        +k(449)*n(idx_CH4j)  &
        +k(839)*n(idx_OHj)  &
        +k(752)*n(idx_NHj)  &
        +k(622)*n(idx_HCNj)

    !d[HOC+_dot]/d[CO]
    pd(72,25) =  &
        +k(585)*n(idx_H3j)

    !d[C+_dot]/d[CO]
    pd(73,25) =  &
        +k(670)*n(idx_HEj)

    !d[NO+_dot]/d[CO]
    pd(79,25) =  &
        +k(721)*n(idx_Nj)

    !d[SI+_dot]/d[CO]
    pd(80,25) =  &
        +k(491)*n(idx_SIOj)

    !d[CN+_dot]/d[CO]
    pd(86,25) =  &
        -k(64)*n(idx_CNj)

    !d[CO+_dot]/d[CO]
    pd(87,25) =  &
        +k(233)  &
        +k(159)*n(idx_Nj)  &
        +k(209)*n(idx_Oj)  &
        +k(105)*n(idx_H2j)  &
        +k(64)*n(idx_CNj)  &
        +k(75)*n(idx_N2j)

    !d[N2+_dot]/d[CO]
    pd(88,25) =  &
        -k(75)*n(idx_N2j)

    !d[H2O+_dot]/d[CO]
    pd(90,25) =  &
        -k(554)*n(idx_H2Oj)

    !d[O+_dot]/d[CO]
    pd(92,25) =  &
        -k(209)*n(idx_Oj)

    !d[OH+_dot]/d[CO]
    pd(93,25) =  &
        -k(839)*n(idx_OHj)

    !d[CH4+_dot]/d[CO]
    pd(95,25) =  &
        -k(449)*n(idx_CH4j)

    !d[N+_dot]/d[CO]
    pd(96,25) =  &
        -k(721)*n(idx_Nj)  &
        -k(159)*n(idx_Nj)

    !d[HCN+_dot]/d[CO]
    pd(97,25) =  &
        -k(622)*n(idx_HCNj)

    !d[NH+_dot]/d[CO]
    pd(98,25) =  &
        -k(752)*n(idx_NHj)

    !d[SIH4+_dot]/d[CO]
    pd(99,25) =  &
        -k(490)*n(idx_SIH4j)

    !d[SIO+_dot]/d[CO]
    pd(101,25) =  &
        -k(491)*n(idx_SIOj)

    !d[H2+_dot]/d[CO]
    pd(102,25) =  &
        -k(516)*n(idx_H2j)  &
        -k(105)*n(idx_H2j)

    !d[HE+_dot]/d[CO]
    pd(103,25) =  &
        -k(670)*n(idx_HEj)

    !d[HNO+_dot]/d[CO]
    pd(104,25) =  &
        -k(487)*n(idx_HNOj)

    !d[H3+_dot]/d[CO]
    pd(106,25) =  &
        -k(584)*n(idx_H3j)  &
        -k(585)*n(idx_H3j)

    !d[HCO2+_dot]/d[CO]
    pd(110,25) =  &
        -k(486)*n(idx_HCO2j)

    !d[N2H+_dot]/d[CO]
    pd(112,25) =  &
        -k(488)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CO]
    pd(113,25) =  &
        -k(489)*n(idx_O2Hj)

    !d[CH_dot]/d[N2]
    pd(2,26) =  &
        -k(928)*n(idx_CH)

    !d[O_dot]/d[N2]
    pd(3,26) =  &
        -k(1072)*n(idx_O)  &
        +k(846)*n(idx_OHj)

    !d[HCN_dot]/d[N2]
    pd(5,26) =  &
        +k(885)*n(idx_CH2)  &
        +k(928)*n(idx_CH)

    !d[H2_dot]/d[N2]
    pd(6,26) =  &
        +k(593)*n(idx_H3j)

    !d[C_dot]/d[N2]
    pd(7,26) =  &
        -k(866)*n(idx_C)

    !d[H_dot]/d[N2]
    pd(8,26) =  &
        +k(522)*n(idx_H2j)

    !d[O2_dot]/d[N2]
    pd(11,26) =  &
        +k(734)*n(idx_O2Hj)

    !d[CH2_dot]/d[N2]
    pd(12,26) =  &
        -k(885)*n(idx_CH2)

    !d[NO_dot]/d[N2]
    pd(17,26) =  &
        +k(733)*n(idx_HNOj)  &
        +k(1072)*n(idx_O)

    !d[CN_dot]/d[N2]
    pd(24,26) =  &
        +k(866)*n(idx_C)

    !d[N2_dot]/d[N2]
    pd(26,26) =  &
        -k(1156)  &
        -k(818)*n(idx_Oj)  &
        -k(268)  &
        -k(762)*n(idx_NHj)  &
        -k(689)*n(idx_HEj)  &
        -k(928)*n(idx_CH)  &
        -k(846)*n(idx_OHj)  &
        -k(1072)*n(idx_O)  &
        -k(145)*n(idx_HEj)  &
        -k(885)*n(idx_CH2)  &
        -k(734)*n(idx_O2Hj)  &
        -k(866)*n(idx_C)  &
        -k(1263)  &
        -k(593)*n(idx_H3j)  &
        -k(522)*n(idx_H2j)  &
        -k(733)*n(idx_HNOj)

    !d[N_dot]/d[N2]
    pd(30,26) =  &
        +k(818)*n(idx_Oj)  &
        +k(1072)*n(idx_O)  &
        +k(928)*n(idx_CH)  &
        +2.d0*k(268)  &
        +k(866)*n(idx_C)  &
        +k(689)*n(idx_HEj)  &
        +2.d0*k(1156)  &
        +k(762)*n(idx_NHj)

    !d[NH_dot]/d[N2]
    pd(31,26) =  &
        +k(885)*n(idx_CH2)

    !d[HE_dot]/d[N2]
    pd(35,26) =  &
        +k(689)*n(idx_HEj)  &
        +k(145)*n(idx_HEj)

    !d[N2_DUST_dot]/d[N2]
    pd(58,26) =  &
        +k(1263)

    !d[NO+_dot]/d[N2]
    pd(79,26) =  &
        +k(818)*n(idx_Oj)

    !d[N2+_dot]/d[N2]
    pd(88,26) =  &
        +k(145)*n(idx_HEj)

    !d[O+_dot]/d[N2]
    pd(92,26) =  &
        -k(818)*n(idx_Oj)

    !d[OH+_dot]/d[N2]
    pd(93,26) =  &
        -k(846)*n(idx_OHj)

    !d[N+_dot]/d[N2]
    pd(96,26) =  &
        +k(689)*n(idx_HEj)

    !d[NH+_dot]/d[N2]
    pd(98,26) =  &
        -k(762)*n(idx_NHj)

    !d[H2+_dot]/d[N2]
    pd(102,26) =  &
        -k(522)*n(idx_H2j)

    !d[HE+_dot]/d[N2]
    pd(103,26) =  &
        -k(145)*n(idx_HEj)  &
        -k(689)*n(idx_HEj)

    !d[HNO+_dot]/d[N2]
    pd(104,26) =  &
        -k(733)*n(idx_HNOj)

    !d[H3+_dot]/d[N2]
    pd(106,26) =  &
        -k(593)*n(idx_H3j)

    !d[N2H+_dot]/d[N2]
    pd(112,26) =  &
        +k(846)*n(idx_OHj)  &
        +k(734)*n(idx_O2Hj)  &
        +k(733)*n(idx_HNOj)  &
        +k(522)*n(idx_H2j)  &
        +k(762)*n(idx_NHj)  &
        +k(593)*n(idx_H3j)

    !d[O2H+_dot]/d[N2]
    pd(113,26) =  &
        -k(734)*n(idx_O2Hj)

    !d[E_dot]/d[NH2]
    pd(1,27) =  &
        +k(270)  &
        +k(1158)

    !d[CH_dot]/d[NH2]
    pd(2,27) =  &
        +k(869)*n(idx_C)

    !d[O_dot]/d[NH2]
    pd(3,27) =  &
        -k(1074)*n(idx_O)  &
        +k(1033)*n(idx_OH)  &
        +k(792)*n(idx_OHj)  &
        +k(213)*n(idx_Oj)  &
        -k(1073)*n(idx_O)

    !d[HNC_dot]/d[NH2]
    pd(4,27) =  &
        +k(868)*n(idx_C)  &
        +k(787)*n(idx_HCNHj)

    !d[HCN_dot]/d[NH2]
    pd(5,27) =  &
        +k(786)*n(idx_HCNHj)  &
        +k(867)*n(idx_C)

    !d[H2_dot]/d[NH2]
    pd(6,27) =  &
        +k(690)*n(idx_HEj)  &
        +k(594)*n(idx_H3j)  &
        +k(984)*n(idx_H)  &
        -k(962)*n(idx_H2)  &
        +k(110)*n(idx_H2j)  &
        +k(410)*n(idx_CHj)

    !d[C_dot]/d[NH2]
    pd(7,27) =  &
        -k(868)*n(idx_C)  &
        -k(869)*n(idx_C)  &
        -k(867)*n(idx_C)

    !d[H_dot]/d[NH2]
    pd(8,27) =  &
        +k(691)*n(idx_HEj)  &
        +k(271)  &
        +k(868)*n(idx_C)  &
        +k(1031)*n(idx_NO)  &
        +k(1073)*n(idx_O)  &
        -k(984)*n(idx_H)  &
        +k(867)*n(idx_C)  &
        +k(962)*n(idx_H2)  &
        +k(85)*n(idx_Hj)  &
        +k(374)*n(idx_Cj)  &
        +k(1159)

    !d[H2O_dot]/d[NH2]
    pd(9,27) =  &
        +k(784)*n(idx_H3Oj)  &
        +k(1030)*n(idx_NO)  &
        +k(1032)*n(idx_OH)  &
        +k(186)*n(idx_H2Oj)

    !d[OH_dot]/d[NH2]
    pd(10,27) =  &
        +k(782)*n(idx_H2Oj)  &
        +k(1031)*n(idx_NO)  &
        -k(1032)*n(idx_OH)  &
        -k(1033)*n(idx_OH)  &
        +k(1074)*n(idx_O)  &
        +k(189)*n(idx_OHj)

    !d[O2_dot]/d[NH2]
    pd(11,27) =  &
        +k(791)*n(idx_O2Hj)  &
        +k(188)*n(idx_O2j)

    !d[H2CO_dot]/d[NH2]
    pd(13,27) =  &
        +k(783)*n(idx_H3COj)

    !d[HCO_dot]/d[NH2]
    pd(14,27) =  &
        +k(781)*n(idx_H2COj)

    !d[NH3_dot]/d[NH2]
    pd(16,27) =  &
        +k(962)*n(idx_H2)  &
        +k(1033)*n(idx_OH)  &
        +k(1029)*n(idx_CH4)

    !d[NO_dot]/d[NH2]
    pd(17,27) =  &
        -k(1031)*n(idx_NO)  &
        +k(789)*n(idx_HNOj)  &
        -k(1030)*n(idx_NO)

    !d[CN_dot]/d[NH2]
    pd(24,27) =  &
        +k(184)*n(idx_CNj)  &
        +k(785)*n(idx_HCNj)

    !d[CO_dot]/d[NH2]
    pd(25,27) =  &
        +k(788)*n(idx_HCOj)  &
        +k(185)*n(idx_COj)

    !d[N2_dot]/d[NH2]
    pd(26,27) =  &
        +k(1030)*n(idx_NO)  &
        +k(790)*n(idx_N2Hj)  &
        +k(1031)*n(idx_NO)  &
        +k(187)*n(idx_N2j)

    !d[NH2_dot]/d[NH2]
    pd(27,27) =  &
        -k(270)  &
        -k(1031)*n(idx_NO)  &
        -k(782)*n(idx_H2Oj)  &
        -k(790)*n(idx_N2Hj)  &
        -k(781)*n(idx_H2COj)  &
        -k(908)*n(idx_CH3)  &
        -k(187)*n(idx_N2j)  &
        -k(777)*n(idx_NH2j)  &
        -k(185)*n(idx_COj)  &
        -k(867)*n(idx_C)  &
        -k(962)*n(idx_H2)  &
        -k(594)*n(idx_H3j)  &
        -k(792)*n(idx_OHj)  &
        -k(213)*n(idx_Oj)  &
        -k(189)*n(idx_OHj)  &
        -k(85)*n(idx_Hj)  &
        -k(788)*n(idx_HCOj)  &
        -k(1159)  &
        -k(1029)*n(idx_CH4)  &
        -k(783)*n(idx_H3COj)  &
        -k(789)*n(idx_HNOj)  &
        -k(186)*n(idx_H2Oj)  &
        -k(690)*n(idx_HEj)  &
        -k(784)*n(idx_H3Oj)  &
        -k(1290)  &
        -k(271)  &
        -k(1032)*n(idx_OH)  &
        -k(984)*n(idx_H)  &
        -k(374)*n(idx_Cj)  &
        -k(763)*n(idx_NHj)  &
        -k(165)*n(idx_Nj)  &
        -k(410)*n(idx_CHj)  &
        -k(868)*n(idx_C)  &
        -k(785)*n(idx_HCNj)  &
        -k(1074)*n(idx_O)  &
        -k(184)*n(idx_CNj)  &
        -k(791)*n(idx_O2Hj)  &
        -k(188)*n(idx_O2j)  &
        -k(1033)*n(idx_OH)  &
        -k(780)*n(idx_COj)  &
        -k(869)*n(idx_C)  &
        -k(1158)  &
        -k(1073)*n(idx_O)  &
        -k(787)*n(idx_HCNHj)  &
        -k(691)*n(idx_HEj)  &
        -k(1030)*n(idx_NO)  &
        -k(786)*n(idx_HCNHj)  &
        -k(110)*n(idx_H2j)

    !d[CH3_dot]/d[NH2]
    pd(28,27) =  &
        +k(1029)*n(idx_CH4)  &
        -k(908)*n(idx_CH3)

    !d[CH4_dot]/d[NH2]
    pd(29,27) =  &
        -k(1029)*n(idx_CH4)  &
        +k(908)*n(idx_CH3)

    !d[N_dot]/d[NH2]
    pd(30,27) =  &
        +k(763)*n(idx_NHj)  &
        +k(165)*n(idx_Nj)

    !d[NH_dot]/d[NH2]
    pd(31,27) =  &
        +k(271)  &
        +k(1159)  &
        +k(777)*n(idx_NH2j)  &
        +k(780)*n(idx_COj)  &
        +k(908)*n(idx_CH3)  &
        +k(984)*n(idx_H)  &
        +k(1074)*n(idx_O)  &
        +k(869)*n(idx_C)  &
        +k(1032)*n(idx_OH)

    !d[HE_dot]/d[NH2]
    pd(35,27) =  &
        +k(690)*n(idx_HEj)  &
        +k(691)*n(idx_HEj)

    !d[HNO_dot]/d[NH2]
    pd(36,27) =  &
        +k(1073)*n(idx_O)

    !d[NH3_DUST_dot]/d[NH2]
    pd(60,27) =  &
        +k(1290)

    !d[HCO+_dot]/d[NH2]
    pd(70,27) =  &
        +k(780)*n(idx_COj)  &
        -k(788)*n(idx_HCOj)

    !d[H+_dot]/d[NH2]
    pd(71,27) =  &
        -k(85)*n(idx_Hj)

    !d[C+_dot]/d[NH2]
    pd(73,27) =  &
        -k(374)*n(idx_Cj)

    !d[CH+_dot]/d[NH2]
    pd(75,27) =  &
        -k(410)*n(idx_CHj)

    !d[H2CO+_dot]/d[NH2]
    pd(76,27) =  &
        -k(781)*n(idx_H2COj)

    !d[NH3+_dot]/d[NH2]
    pd(78,27) =  &
        +k(594)*n(idx_H3j)  &
        +k(777)*n(idx_NH2j)  &
        +k(782)*n(idx_H2Oj)  &
        +k(785)*n(idx_HCNj)  &
        +k(789)*n(idx_HNOj)  &
        +k(786)*n(idx_HCNHj)  &
        +k(790)*n(idx_N2Hj)  &
        +k(791)*n(idx_O2Hj)  &
        +k(781)*n(idx_H2COj)  &
        +k(784)*n(idx_H3Oj)  &
        +k(783)*n(idx_H3COj)  &
        +k(787)*n(idx_HCNHj)  &
        +k(763)*n(idx_NHj)  &
        +k(792)*n(idx_OHj)  &
        +k(788)*n(idx_HCOj)

    !d[CN+_dot]/d[NH2]
    pd(86,27) =  &
        -k(184)*n(idx_CNj)

    !d[CO+_dot]/d[NH2]
    pd(87,27) =  &
        -k(780)*n(idx_COj)  &
        -k(185)*n(idx_COj)

    !d[N2+_dot]/d[NH2]
    pd(88,27) =  &
        -k(187)*n(idx_N2j)

    !d[O2+_dot]/d[NH2]
    pd(89,27) =  &
        -k(188)*n(idx_O2j)

    !d[H2O+_dot]/d[NH2]
    pd(90,27) =  &
        -k(186)*n(idx_H2Oj)  &
        -k(782)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NH2]
    pd(91,27) =  &
        +k(270)  &
        +k(184)*n(idx_CNj)  &
        +k(186)*n(idx_H2Oj)  &
        +k(188)*n(idx_O2j)  &
        +k(213)*n(idx_Oj)  &
        +k(85)*n(idx_Hj)  &
        -k(777)*n(idx_NH2j)  &
        +k(110)*n(idx_H2j)  &
        +k(187)*n(idx_N2j)  &
        +k(165)*n(idx_Nj)  &
        +k(1158)  &
        +k(185)*n(idx_COj)  &
        +k(189)*n(idx_OHj)

    !d[O+_dot]/d[NH2]
    pd(92,27) =  &
        -k(213)*n(idx_Oj)

    !d[OH+_dot]/d[NH2]
    pd(93,27) =  &
        -k(792)*n(idx_OHj)  &
        -k(189)*n(idx_OHj)

    !d[N+_dot]/d[NH2]
    pd(96,27) =  &
        +k(690)*n(idx_HEj)  &
        -k(165)*n(idx_Nj)

    !d[HCN+_dot]/d[NH2]
    pd(97,27) =  &
        +k(374)*n(idx_Cj)  &
        -k(785)*n(idx_HCNj)  &
        +k(410)*n(idx_CHj)

    !d[NH+_dot]/d[NH2]
    pd(98,27) =  &
        +k(691)*n(idx_HEj)  &
        -k(763)*n(idx_NHj)

    !d[H2+_dot]/d[NH2]
    pd(102,27) =  &
        -k(110)*n(idx_H2j)

    !d[HE+_dot]/d[NH2]
    pd(103,27) =  &
        -k(690)*n(idx_HEj)  &
        -k(691)*n(idx_HEj)

    !d[HNO+_dot]/d[NH2]
    pd(104,27) =  &
        -k(789)*n(idx_HNOj)

    !d[H3+_dot]/d[NH2]
    pd(106,27) =  &
        -k(594)*n(idx_H3j)

    !d[H3CO+_dot]/d[NH2]
    pd(107,27) =  &
        -k(783)*n(idx_H3COj)

    !d[H3O+_dot]/d[NH2]
    pd(108,27) =  &
        -k(784)*n(idx_H3Oj)

    !d[HCNH+_dot]/d[NH2]
    pd(109,27) =  &
        -k(787)*n(idx_HCNHj)  &
        -k(786)*n(idx_HCNHj)

    !d[N2H+_dot]/d[NH2]
    pd(112,27) =  &
        -k(790)*n(idx_N2Hj)

    !d[O2H+_dot]/d[NH2]
    pd(113,27) =  &
        -k(791)*n(idx_O2Hj)

    !d[E_dot]/d[CH3]
    pd(1,28) =  &
        +k(246)  &
        +k(1118)

    !d[CH_dot]/d[CH3]
    pd(2,28) =  &
        +k(1119)  &
        +k(247)

    !d[O_dot]/d[CH3]
    pd(3,28) =  &
        +k(918)*n(idx_OH)  &
        -k(916)*n(idx_O)  &
        -k(917)*n(idx_O)

    !d[HCN_dot]/d[CH3]
    pd(5,28) =  &
        +k(903)*n(idx_CN)  &
        +k(911)*n(idx_NO)  &
        +k(1011)*n(idx_N)  &
        +k(1010)*n(idx_N)

    !d[H2_dot]/d[CH3]
    pd(6,28) =  &
        +k(579)*n(idx_H3j)  &
        +k(1119)  &
        +k(656)*n(idx_HEj)  &
        +k(1010)*n(idx_N)  &
        +k(916)*n(idx_O)  &
        -k(958)*n(idx_H2)  &
        +k(919)*n(idx_OH)  &
        +k(969)*n(idx_H)  &
        +k(247)

    !d[H_dot]/d[CH3]
    pd(8,28) =  &
        +k(917)*n(idx_O)  &
        +k(958)*n(idx_H2)  &
        +k(245)  &
        +k(916)*n(idx_O)  &
        +k(77)*n(idx_Hj)  &
        -k(969)*n(idx_H)  &
        +2.d0*k(1011)*n(idx_N)  &
        +k(1009)*n(idx_N)  &
        +k(1117)

    !d[H2O_dot]/d[CH3]
    pd(9,28) =  &
        +k(920)*n(idx_OH)  &
        -k(905)*n(idx_H2O)  &
        +k(911)*n(idx_NO)  &
        +k(913)*n(idx_O2)

    !d[OH_dot]/d[CH3]
    pd(10,28) =  &
        -k(918)*n(idx_OH)  &
        +k(905)*n(idx_H2O)  &
        -k(920)*n(idx_OH)  &
        +k(912)*n(idx_O2)  &
        -k(919)*n(idx_OH)

    !d[O2_dot]/d[CH3]
    pd(11,28) =  &
        -k(914)*n(idx_O2)  &
        -k(913)*n(idx_O2)  &
        -k(912)*n(idx_O2)  &
        +k(915)*n(idx_O2H)

    !d[CH2_dot]/d[CH3]
    pd(12,28) =  &
        +k(903)*n(idx_CN)  &
        +k(245)  &
        +2.d0*k(902)*n(idx_CH3)  &
        +k(969)*n(idx_H)  &
        +k(914)*n(idx_O2)  &
        +k(920)*n(idx_OH)  &
        +k(1117)

    !d[H2CO_dot]/d[CH3]
    pd(13,28) =  &
        -k(904)*n(idx_H2CO)  &
        +k(917)*n(idx_O)  &
        +k(910)*n(idx_NO2)  &
        +k(912)*n(idx_O2)  &
        +k(919)*n(idx_OH)

    !d[HCO_dot]/d[CH3]
    pd(14,28) =  &
        -k(906)*n(idx_HCO)  &
        +k(913)*n(idx_O2)  &
        +k(904)*n(idx_H2CO)

    !d[NH3_dot]/d[CH3]
    pd(16,28) =  &
        -k(909)*n(idx_NH3)

    !d[NO_dot]/d[CH3]
    pd(17,28) =  &
        +k(907)*n(idx_HNO)  &
        -k(911)*n(idx_NO)

    !d[CN_dot]/d[CH3]
    pd(24,28) =  &
        -k(903)*n(idx_CN)

    !d[CO_dot]/d[CH3]
    pd(25,28) =  &
        +k(906)*n(idx_HCO)  &
        +k(916)*n(idx_O)

    !d[NH2_dot]/d[CH3]
    pd(27,28) =  &
        -k(908)*n(idx_NH2)  &
        +k(909)*n(idx_NH3)

    !d[CH3_dot]/d[CH3]
    pd(28,28) =  &
        -k(656)*n(idx_HEj)  &
        -k(907)*n(idx_HNO)  &
        -k(920)*n(idx_OH)  &
        -k(1010)*n(idx_N)  &
        -k(917)*n(idx_O)  &
        -k(245)  &
        -4.d0*k(902)*n(idx_CH3)  &
        -k(905)*n(idx_H2O)  &
        -k(958)*n(idx_H2)  &
        -k(908)*n(idx_NH2)  &
        -k(77)*n(idx_Hj)  &
        -k(916)*n(idx_O)  &
        -k(912)*n(idx_O2)  &
        -k(1117)  &
        -k(906)*n(idx_HCO)  &
        -k(247)  &
        -k(903)*n(idx_CN)  &
        -k(913)*n(idx_O2)  &
        -k(969)*n(idx_H)  &
        -k(1009)*n(idx_N)  &
        -k(1119)  &
        -k(579)*n(idx_H3j)  &
        -k(1260)  &
        -k(918)*n(idx_OH)  &
        -k(910)*n(idx_NO2)  &
        -k(914)*n(idx_O2)  &
        -k(904)*n(idx_H2CO)  &
        -k(1011)*n(idx_N)  &
        -k(909)*n(idx_NH3)  &
        -k(919)*n(idx_OH)  &
        -k(911)*n(idx_NO)  &
        -k(915)*n(idx_O2H)  &
        -k(246)  &
        -k(1118)

    !d[CH4_dot]/d[CH3]
    pd(29,28) =  &
        +k(958)*n(idx_H2)  &
        +k(909)*n(idx_NH3)  &
        +k(904)*n(idx_H2CO)  &
        +k(908)*n(idx_NH2)  &
        +2.d0*k(902)*n(idx_CH3)  &
        +k(907)*n(idx_HNO)  &
        +k(915)*n(idx_O2H)  &
        +k(905)*n(idx_H2O)  &
        +k(918)*n(idx_OH)  &
        +k(906)*n(idx_HCO)

    !d[N_dot]/d[CH3]
    pd(30,28) =  &
        -k(1011)*n(idx_N)  &
        -k(1010)*n(idx_N)  &
        -k(1009)*n(idx_N)

    !d[NH_dot]/d[CH3]
    pd(31,28) =  &
        +k(908)*n(idx_NH2)

    !d[HE_dot]/d[CH3]
    pd(35,28) =  &
        +k(656)*n(idx_HEj)

    !d[HNO_dot]/d[CH3]
    pd(36,28) =  &
        -k(907)*n(idx_HNO)  &
        +k(910)*n(idx_NO2)

    !d[H2CN_dot]/d[CH3]
    pd(39,28) =  &
        +k(1009)*n(idx_N)

    !d[NO2_dot]/d[CH3]
    pd(42,28) =  &
        -k(910)*n(idx_NO2)

    !d[O2H_dot]/d[CH3]
    pd(43,28) =  &
        +k(914)*n(idx_O2)  &
        -k(915)*n(idx_O2H)

    !d[CH4_DUST_dot]/d[CH3]
    pd(53,28) =  &
        +k(1260)

    !d[H+_dot]/d[CH3]
    pd(71,28) =  &
        -k(77)*n(idx_Hj)

    !d[CH+_dot]/d[CH3]
    pd(75,28) =  &
        +k(656)*n(idx_HEj)

    !d[CH3+_dot]/d[CH3]
    pd(94,28) =  &
        +k(246)  &
        +k(1118)  &
        +k(77)*n(idx_Hj)

    !d[CH4+_dot]/d[CH3]
    pd(95,28) =  &
        +k(579)*n(idx_H3j)

    !d[HE+_dot]/d[CH3]
    pd(103,28) =  &
        -k(656)*n(idx_HEj)

    !d[H3+_dot]/d[CH3]
    pd(106,28) =  &
        -k(579)*n(idx_H3j)

    !d[E_dot]/d[CH4]
    pd(1,29) =  &
        +k(1127)

    !d[CH_dot]/d[CH4]
    pd(2,29) =  &
        +k(1128)

    !d[O_dot]/d[CH4]
    pd(3,29) =  &
        -k(1057)*n(idx_O)  &
        +k(208)*n(idx_Oj)

    !d[HCN_dot]/d[CH4]
    pd(5,29) =  &
        +k(921)*n(idx_CN)

    !d[H2_dot]/d[CH4]
    pd(6,29) =  &
        +k(1125)  &
        +k(456)*n(idx_N2j)  &
        +k(718)*n(idx_Nj)  &
        +k(496)*n(idx_Hj)  &
        +k(970)*n(idx_H)  &
        +k(1128)  &
        +k(660)*n(idx_HEj)  &
        +k(250)  &
        +k(659)*n(idx_HEj)  &
        +k(512)*n(idx_H2j)  &
        +k(102)*n(idx_H2j)

    !d[H_dot]/d[CH4]
    pd(8,29) =  &
        +k(457)*n(idx_N2j)  &
        +2.d0*k(719)*n(idx_Nj)  &
        +k(717)*n(idx_Nj)  &
        -k(970)*n(idx_H)  &
        +k(512)*n(idx_H2j)  &
        +k(78)*n(idx_Hj)  &
        +k(1128)  &
        +k(661)*n(idx_HEj)  &
        +k(718)*n(idx_Nj)  &
        +k(659)*n(idx_HEj)  &
        +k(1126)

    !d[H2O_dot]/d[CH4]
    pd(9,29) =  &
        +k(923)*n(idx_OH)

    !d[OH_dot]/d[CH4]
    pd(10,29) =  &
        +k(811)*n(idx_Oj)  &
        -k(923)*n(idx_OH)  &
        +k(1057)*n(idx_O)

    !d[O2_dot]/d[CH4]
    pd(11,29) =  &
        -k(922)*n(idx_O2)

    !d[CH2_dot]/d[CH4]
    pd(12,29) =  &
        +k(1125)  &
        -k(880)*n(idx_CH2)  &
        +k(250)  &
        +k(458)*n(idx_OHj)

    !d[NH3_dot]/d[CH4]
    pd(16,29) =  &
        +k(1029)*n(idx_NH2)

    !d[CN_dot]/d[CH4]
    pd(24,29) =  &
        -k(921)*n(idx_CN)

    !d[CO_dot]/d[CH4]
    pd(25,29) =  &
        +k(53)*n(idx_COj)

    !d[N2_dot]/d[CH4]
    pd(26,29) =  &
        +k(457)*n(idx_N2j)  &
        +k(456)*n(idx_N2j)

    !d[NH2_dot]/d[CH4]
    pd(27,29) =  &
        -k(1029)*n(idx_NH2)  &
        +k(1035)*n(idx_NH)

    !d[CH3_dot]/d[CH4]
    pd(28,29) =  &
        +k(662)*n(idx_HEj)  &
        +k(453)*n(idx_H2COj)  &
        +k(454)*n(idx_H2Oj)  &
        +k(923)*n(idx_OH)  &
        +2.d0*k(880)*n(idx_CH2)  &
        +k(921)*n(idx_CN)  &
        +k(452)*n(idx_COj)  &
        +k(1029)*n(idx_NH2)  &
        +k(455)*n(idx_HCNj)  &
        +k(970)*n(idx_H)  &
        +k(1126)  &
        +k(1035)*n(idx_NH)  &
        +k(1057)*n(idx_O)  &
        +k(922)*n(idx_O2)

    !d[CH4_dot]/d[CH4]
    pd(29,29) =  &
        -k(496)*n(idx_Hj)  &
        -k(923)*n(idx_OH)  &
        -k(1261)  &
        -k(455)*n(idx_HCNj)  &
        -k(78)*n(idx_Hj)  &
        -k(208)*n(idx_Oj)  &
        -k(661)*n(idx_HEj)  &
        -k(454)*n(idx_H2Oj)  &
        -k(157)*n(idx_Nj)  &
        -k(141)*n(idx_HEj)  &
        -k(512)*n(idx_H2j)  &
        -k(921)*n(idx_CN)  &
        -k(458)*n(idx_OHj)  &
        -k(718)*n(idx_Nj)  &
        -k(659)*n(idx_HEj)  &
        -k(453)*n(idx_H2COj)  &
        -k(1057)*n(idx_O)  &
        -k(1126)  &
        -k(1029)*n(idx_NH2)  &
        -k(102)*n(idx_H2j)  &
        -k(660)*n(idx_HEj)  &
        -k(1035)*n(idx_NH)  &
        -k(1128)  &
        -k(880)*n(idx_CH2)  &
        -k(719)*n(idx_Nj)  &
        -k(452)*n(idx_COj)  &
        -k(457)*n(idx_N2j)  &
        -k(662)*n(idx_HEj)  &
        -k(970)*n(idx_H)  &
        -k(456)*n(idx_N2j)  &
        -k(250)  &
        -k(717)*n(idx_Nj)  &
        -k(1125)  &
        -k(53)*n(idx_COj)  &
        -k(811)*n(idx_Oj)  &
        -k(1127)  &
        -k(922)*n(idx_O2)

    !d[N_dot]/d[CH4]
    pd(30,29) =  &
        +k(157)*n(idx_Nj)  &
        +k(717)*n(idx_Nj)

    !d[NH_dot]/d[CH4]
    pd(31,29) =  &
        -k(1035)*n(idx_NH)

    !d[HE_dot]/d[CH4]
    pd(35,29) =  &
        +k(661)*n(idx_HEj)  &
        +k(141)*n(idx_HEj)  &
        +k(662)*n(idx_HEj)  &
        +k(659)*n(idx_HEj)  &
        +k(660)*n(idx_HEj)

    !d[O2H_dot]/d[CH4]
    pd(43,29) =  &
        +k(922)*n(idx_O2)

    !d[CH4_DUST_dot]/d[CH4]
    pd(53,29) =  &
        +k(1261)

    !d[HCO+_dot]/d[CH4]
    pd(70,29) =  &
        +k(452)*n(idx_COj)

    !d[H+_dot]/d[CH4]
    pd(71,29) =  &
        -k(78)*n(idx_Hj)  &
        +k(662)*n(idx_HEj)  &
        -k(496)*n(idx_Hj)

    !d[CH2+_dot]/d[CH4]
    pd(74,29) =  &
        +k(660)*n(idx_HEj)  &
        +k(456)*n(idx_N2j)

    !d[CH+_dot]/d[CH4]
    pd(75,29) =  &
        +k(659)*n(idx_HEj)

    !d[H2CO+_dot]/d[CH4]
    pd(76,29) =  &
        -k(453)*n(idx_H2COj)

    !d[CO+_dot]/d[CH4]
    pd(87,29) =  &
        -k(452)*n(idx_COj)  &
        -k(53)*n(idx_COj)

    !d[N2+_dot]/d[CH4]
    pd(88,29) =  &
        -k(457)*n(idx_N2j)  &
        -k(456)*n(idx_N2j)

    !d[H2O+_dot]/d[CH4]
    pd(90,29) =  &
        -k(454)*n(idx_H2Oj)

    !d[O+_dot]/d[CH4]
    pd(92,29) =  &
        -k(208)*n(idx_Oj)  &
        -k(811)*n(idx_Oj)

    !d[OH+_dot]/d[CH4]
    pd(93,29) =  &
        -k(458)*n(idx_OHj)

    !d[CH3+_dot]/d[CH4]
    pd(94,29) =  &
        +k(457)*n(idx_N2j)  &
        +k(811)*n(idx_Oj)  &
        +k(717)*n(idx_Nj)  &
        +k(496)*n(idx_Hj)  &
        +k(661)*n(idx_HEj)  &
        +k(512)*n(idx_H2j)

    !d[CH4+_dot]/d[CH4]
    pd(95,29) =  &
        +k(157)*n(idx_Nj)  &
        +k(208)*n(idx_Oj)  &
        +k(78)*n(idx_Hj)  &
        +k(102)*n(idx_H2j)  &
        +k(141)*n(idx_HEj)  &
        +k(1127)  &
        +k(53)*n(idx_COj)

    !d[N+_dot]/d[CH4]
    pd(96,29) =  &
        -k(718)*n(idx_Nj)  &
        -k(717)*n(idx_Nj)  &
        -k(719)*n(idx_Nj)  &
        -k(157)*n(idx_Nj)

    !d[HCN+_dot]/d[CH4]
    pd(97,29) =  &
        -k(455)*n(idx_HCNj)  &
        +k(718)*n(idx_Nj)

    !d[H2+_dot]/d[CH4]
    pd(102,29) =  &
        -k(512)*n(idx_H2j)  &
        -k(102)*n(idx_H2j)

    !d[HE+_dot]/d[CH4]
    pd(103,29) =  &
        -k(660)*n(idx_HEj)  &
        -k(659)*n(idx_HEj)  &
        -k(662)*n(idx_HEj)  &
        -k(661)*n(idx_HEj)  &
        -k(141)*n(idx_HEj)

    !d[H3CO+_dot]/d[CH4]
    pd(107,29) =  &
        +k(453)*n(idx_H2COj)

    !d[H3O+_dot]/d[CH4]
    pd(108,29) =  &
        +k(458)*n(idx_OHj)  &
        +k(454)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[CH4]
    pd(109,29) =  &
        +k(719)*n(idx_Nj)  &
        +k(455)*n(idx_HCNj)

    !d[E_dot]/d[N]
    pd(1,30) =  &
        +k(239)  &
        +k(269)

    !d[CH_dot]/d[N]
    pd(2,30) =  &
        -k(929)*n(idx_CH)  &
        -k(930)*n(idx_CH)  &
        +k(1008)*n(idx_CH2)

    !d[O_dot]/d[N]
    pd(3,30) =  &
        +k(1023)*n(idx_NO)  &
        +k(1016)*n(idx_HCO)  &
        +k(1027)*n(idx_OH)  &
        +k(1024)*n(idx_O2)  &
        +2.d0*k(1020)*n(idx_NO2)  &
        +k(743)*n(idx_O2j)

    !d[HNC_dot]/d[N]
    pd(4,30) =  &
        +k(1007)*n(idx_CH2)

    !d[HCN_dot]/d[N]
    pd(5,30) =  &
        +k(1016)*n(idx_HCO)  &
        +k(1014)*n(idx_H2CN)  &
        +k(1010)*n(idx_CH3)  &
        +k(1006)*n(idx_CH2)  &
        +k(1011)*n(idx_CH3)

    !d[H2_dot]/d[N]
    pd(6,30) =  &
        +k(1010)*n(idx_CH3)  &
        -k(961)*n(idx_H2)  &
        +k(740)*n(idx_H2Oj)

    !d[C_dot]/d[N]
    pd(7,30) =  &
        +k(930)*n(idx_CH)  &
        +k(1012)*n(idx_CN)  &
        +k(738)*n(idx_CNj)  &
        -k(1195)*n(idx_C)

    !d[H_dot]/d[N]
    pd(8,30) =  &
        +k(742)*n(idx_NH2j)  &
        +2.d0*k(1011)*n(idx_CH3)  &
        +k(744)*n(idx_OHj)  &
        +k(739)*n(idx_H2Oj)  &
        +k(1009)*n(idx_CH3)  &
        +k(929)*n(idx_CH)  &
        +k(741)*n(idx_NHj)  &
        +k(1019)*n(idx_NH)  &
        +k(523)*n(idx_H2j)  &
        +k(961)*n(idx_H2)  &
        +k(1007)*n(idx_CH2)  &
        +k(1026)*n(idx_OH)  &
        +k(737)*n(idx_CH2j)  &
        +k(409)*n(idx_CHj)  &
        +k(1006)*n(idx_CH2)  &
        +k(1017)*n(idx_HCO)

    !d[OH_dot]/d[N]
    pd(10,30) =  &
        -k(1027)*n(idx_OH)  &
        -k(1026)*n(idx_OH)

    !d[O2_dot]/d[N]
    pd(11,30) =  &
        -k(1024)*n(idx_O2)  &
        +k(1025)*n(idx_O2H)  &
        +k(1022)*n(idx_NO2)

    !d[CH2_dot]/d[N]
    pd(12,30) =  &
        -k(1008)*n(idx_CH2)  &
        -k(1007)*n(idx_CH2)  &
        -k(1006)*n(idx_CH2)

    !d[HCO_dot]/d[N]
    pd(14,30) =  &
        -k(1017)*n(idx_HCO)  &
        -k(1016)*n(idx_HCO)  &
        -k(1015)*n(idx_HCO)

    !d[NO_dot]/d[N]
    pd(17,30) =  &
        +2.d0*k(1021)*n(idx_NO2)  &
        +k(1013)*n(idx_CO2)  &
        +k(1024)*n(idx_O2)  &
        +k(1026)*n(idx_OH)  &
        +k(747)*n(idx_SIOj)  &
        -k(1023)*n(idx_NO)  &
        +k(1018)*n(idx_HNO)

    !d[SI_dot]/d[N]
    pd(18,30) =  &
        +k(746)*n(idx_SIOj)  &
        +k(1028)*n(idx_SIC)

    !d[SIC_dot]/d[N]
    pd(21,30) =  &
        -k(1028)*n(idx_SIC)

    !d[CN_dot]/d[N]
    pd(24,30) =  &
        +k(1195)*n(idx_C)  &
        -k(1012)*n(idx_CN)  &
        +k(745)*n(idx_SICj)  &
        +k(1028)*n(idx_SIC)  &
        +k(929)*n(idx_CH)

    !d[CO_dot]/d[N]
    pd(25,30) =  &
        +k(1013)*n(idx_CO2)  &
        +k(1015)*n(idx_HCO)

    !d[N2_dot]/d[N]
    pd(26,30) =  &
        +k(1023)*n(idx_NO)  &
        +k(1019)*n(idx_NH)  &
        +k(175)*n(idx_N2j)  &
        +k(1022)*n(idx_NO2)  &
        +k(1012)*n(idx_CN)  &
        +k(1020)*n(idx_NO2)

    !d[CH3_dot]/d[N]
    pd(28,30) =  &
        -k(1010)*n(idx_CH3)  &
        -k(1009)*n(idx_CH3)  &
        -k(1011)*n(idx_CH3)

    !d[N_dot]/d[N]
    pd(30,30) =  &
        -k(1020)*n(idx_NO2)  &
        -k(1026)*n(idx_OH)  &
        -k(1013)*n(idx_CO2)  &
        -k(745)*n(idx_SICj)  &
        -k(742)*n(idx_NH2j)  &
        -k(929)*n(idx_CH)  &
        -k(1018)*n(idx_HNO)  &
        -k(1211)*n(idx_Nj)  &
        -k(1195)*n(idx_C)  &
        -k(523)*n(idx_H2j)  &
        -k(961)*n(idx_H2)  &
        -k(741)*n(idx_NHj)  &
        -k(1027)*n(idx_OH)  &
        -k(1007)*n(idx_CH2)  &
        -k(1014)*n(idx_H2CN)  &
        -k(409)*n(idx_CHj)  &
        -k(743)*n(idx_O2j)  &
        -k(1017)*n(idx_HCO)  &
        -k(1011)*n(idx_CH3)  &
        -k(740)*n(idx_H2Oj)  &
        -k(1019)*n(idx_NH)  &
        -k(738)*n(idx_CNj)  &
        -k(744)*n(idx_OHj)  &
        -k(1021)*n(idx_NO2)  &
        -k(1006)*n(idx_CH2)  &
        -k(1291)  &
        -k(1010)*n(idx_CH3)  &
        -k(1028)*n(idx_SIC)  &
        -k(1022)*n(idx_NO2)  &
        -k(1012)*n(idx_CN)  &
        -k(747)*n(idx_SIOj)  &
        -k(1024)*n(idx_O2)  &
        -k(1193)*n(idx_Cj)  &
        -k(930)*n(idx_CH)  &
        -k(746)*n(idx_SIOj)  &
        -k(1008)*n(idx_CH2)  &
        -k(269)  &
        -k(175)*n(idx_N2j)  &
        -k(739)*n(idx_H2Oj)  &
        -k(1025)*n(idx_O2H)  &
        -k(737)*n(idx_CH2j)  &
        -k(239)  &
        -k(1009)*n(idx_CH3)  &
        -k(1016)*n(idx_HCO)  &
        -k(1015)*n(idx_HCO)  &
        -k(1023)*n(idx_NO)

    !d[NH_dot]/d[N]
    pd(31,30) =  &
        +k(1015)*n(idx_HCO)  &
        +k(1027)*n(idx_OH)  &
        +k(1008)*n(idx_CH2)  &
        +k(1025)*n(idx_O2H)  &
        +k(930)*n(idx_CH)  &
        +k(1014)*n(idx_H2CN)  &
        +k(961)*n(idx_H2)  &
        +k(1018)*n(idx_HNO)  &
        -k(1019)*n(idx_NH)

    !d[HNO_dot]/d[N]
    pd(36,30) =  &
        -k(1018)*n(idx_HNO)

    !d[CO2_dot]/d[N]
    pd(38,30) =  &
        -k(1013)*n(idx_CO2)

    !d[H2CN_dot]/d[N]
    pd(39,30) =  &
        -k(1014)*n(idx_H2CN)  &
        +k(1009)*n(idx_CH3)

    !d[NO2_dot]/d[N]
    pd(42,30) =  &
        -k(1022)*n(idx_NO2)  &
        -k(1020)*n(idx_NO2)  &
        -k(1021)*n(idx_NO2)

    !d[O2H_dot]/d[N]
    pd(43,30) =  &
        -k(1025)*n(idx_O2H)

    !d[OCN_dot]/d[N]
    pd(44,30) =  &
        +k(1017)*n(idx_HCO)

    !d[NH3_DUST_dot]/d[N]
    pd(60,30) =  &
        +k(1291)

    !d[C+_dot]/d[N]
    pd(73,30) =  &
        -k(1193)*n(idx_Cj)

    !d[CH2+_dot]/d[N]
    pd(74,30) =  &
        -k(737)*n(idx_CH2j)

    !d[CH+_dot]/d[N]
    pd(75,30) =  &
        -k(409)*n(idx_CHj)

    !d[NO+_dot]/d[N]
    pd(79,30) =  &
        +k(746)*n(idx_SIOj)  &
        +k(743)*n(idx_O2j)  &
        +k(744)*n(idx_OHj)  &
        +k(740)*n(idx_H2Oj)

    !d[SI+_dot]/d[N]
    pd(80,30) =  &
        +k(747)*n(idx_SIOj)  &
        +k(745)*n(idx_SICj)

    !d[SIC+_dot]/d[N]
    pd(83,30) =  &
        -k(745)*n(idx_SICj)

    !d[CN+_dot]/d[N]
    pd(86,30) =  &
        -k(738)*n(idx_CNj)  &
        +k(409)*n(idx_CHj)  &
        +k(1193)*n(idx_Cj)

    !d[N2+_dot]/d[N]
    pd(88,30) =  &
        +k(741)*n(idx_NHj)  &
        -k(175)*n(idx_N2j)  &
        +k(738)*n(idx_CNj)  &
        +k(1211)*n(idx_Nj)

    !d[O2+_dot]/d[N]
    pd(89,30) =  &
        -k(743)*n(idx_O2j)

    !d[H2O+_dot]/d[N]
    pd(90,30) =  &
        -k(739)*n(idx_H2Oj)  &
        -k(740)*n(idx_H2Oj)

    !d[NH2+_dot]/d[N]
    pd(91,30) =  &
        -k(742)*n(idx_NH2j)

    !d[OH+_dot]/d[N]
    pd(93,30) =  &
        -k(744)*n(idx_OHj)

    !d[N+_dot]/d[N]
    pd(96,30) =  &
        +k(175)*n(idx_N2j)  &
        -k(1211)*n(idx_Nj)  &
        +k(239)  &
        +k(269)

    !d[HCN+_dot]/d[N]
    pd(97,30) =  &
        +k(737)*n(idx_CH2j)

    !d[NH+_dot]/d[N]
    pd(98,30) =  &
        +k(523)*n(idx_H2j)  &
        -k(741)*n(idx_NHj)

    !d[SIO+_dot]/d[N]
    pd(101,30) =  &
        -k(747)*n(idx_SIOj)  &
        -k(746)*n(idx_SIOj)

    !d[H2+_dot]/d[N]
    pd(102,30) =  &
        -k(523)*n(idx_H2j)

    !d[HNO+_dot]/d[N]
    pd(104,30) =  &
        +k(739)*n(idx_H2Oj)

    !d[N2H+_dot]/d[N]
    pd(112,30) =  &
        +k(742)*n(idx_NH2j)

    !d[E_dot]/d[NH]
    pd(1,31) =  &
        +k(1164)  &
        +k(276)

    !d[CH_dot]/d[NH]
    pd(2,31) =  &
        +k(871)*n(idx_C)

    !d[O_dot]/d[NH]
    pd(3,31) =  &
        -k(1047)*n(idx_O)  &
        +k(1043)*n(idx_NO)  &
        +k(807)*n(idx_OHj)  &
        +k(805)*n(idx_O2j)  &
        +k(1051)*n(idx_OH)  &
        +k(203)*n(idx_Oj)  &
        -k(1048)*n(idx_O)  &
        +k(1045)*n(idx_O2)

    !d[HCN_dot]/d[NH]
    pd(5,31) =  &
        +k(1036)*n(idx_CN)

    !d[H2_dot]/d[NH]
    pd(6,31) =  &
        +k(795)*n(idx_CH3j)  &
        +k(986)*n(idx_H)  &
        +k(595)*n(idx_H3j)  &
        +k(411)*n(idx_CHj)  &
        +k(112)*n(idx_H2j)  &
        +2.d0*k(1039)*n(idx_NH)  &
        -k(963)*n(idx_H2)

    !d[C_dot]/d[NH]
    pd(7,31) =  &
        -k(871)*n(idx_C)  &
        -k(870)*n(idx_C)

    !d[H_dot]/d[NH]
    pd(8,31) =  &
        +k(1047)*n(idx_O)  &
        +k(727)*n(idx_Nj)  &
        +k(1050)*n(idx_OH)  &
        +k(963)*n(idx_H2)  &
        +k(1163)  &
        +k(376)*n(idx_Cj)  &
        +k(1043)*n(idx_NO)  &
        +k(694)*n(idx_HEj)  &
        +k(275)  &
        +4.d0*k(1040)*n(idx_NH)  &
        +k(870)*n(idx_C)  &
        +k(87)*n(idx_Hj)  &
        -k(986)*n(idx_H)  &
        +k(804)*n(idx_Oj)  &
        +k(524)*n(idx_H2j)  &
        +k(1019)*n(idx_N)

    !d[H2O_dot]/d[NH]
    pd(9,31) =  &
        -k(1037)*n(idx_H2O)  &
        +k(1049)*n(idx_OH)

    !d[OH_dot]/d[NH]
    pd(10,31) =  &
        +k(1046)*n(idx_O2)  &
        +k(1037)*n(idx_H2O)  &
        -k(1049)*n(idx_OH)  &
        +k(1044)*n(idx_NO)  &
        +k(1048)*n(idx_O)  &
        -k(1050)*n(idx_OH)  &
        -k(1051)*n(idx_OH)

    !d[O2_dot]/d[NH]
    pd(11,31) =  &
        -k(1045)*n(idx_O2)  &
        -k(1046)*n(idx_O2)  &
        +k(806)*n(idx_O2Hj)

    !d[NH3_dot]/d[NH]
    pd(16,31) =  &
        -k(1038)*n(idx_NH3)

    !d[NO_dot]/d[NH]
    pd(17,31) =  &
        +k(1047)*n(idx_O)  &
        +k(1046)*n(idx_O2)  &
        -k(1044)*n(idx_NO)  &
        +k(801)*n(idx_HNOj)  &
        -k(1043)*n(idx_NO)  &
        +k(1042)*n(idx_NO2)

    !d[CN_dot]/d[NH]
    pd(24,31) =  &
        +k(870)*n(idx_C)  &
        +k(200)*n(idx_CNj)  &
        -k(1036)*n(idx_CN)  &
        +k(799)*n(idx_HCNj)

    !d[CO_dot]/d[NH]
    pd(25,31) =  &
        +k(201)*n(idx_COj)  &
        +k(800)*n(idx_HCOj)

    !d[N2_dot]/d[NH]
    pd(26,31) =  &
        +k(202)*n(idx_N2j)  &
        +k(802)*n(idx_N2Hj)  &
        +k(1043)*n(idx_NO)  &
        +2.d0*k(1040)*n(idx_NH)  &
        +k(1044)*n(idx_NO)  &
        +2.d0*k(1039)*n(idx_NH)  &
        +k(1019)*n(idx_N)

    !d[NH2_dot]/d[NH]
    pd(27,31) =  &
        +k(963)*n(idx_H2)  &
        +k(1035)*n(idx_CH4)  &
        +2.d0*k(1038)*n(idx_NH3)  &
        +2.d0*k(1041)*n(idx_NH)  &
        +k(1037)*n(idx_H2O)  &
        +k(1051)*n(idx_OH)

    !d[CH3_dot]/d[NH]
    pd(28,31) =  &
        +k(1035)*n(idx_CH4)

    !d[CH4_dot]/d[NH]
    pd(29,31) =  &
        -k(1035)*n(idx_CH4)

    !d[N_dot]/d[NH]
    pd(30,31) =  &
        +k(986)*n(idx_H)  &
        +k(1163)  &
        +k(167)*n(idx_Nj)  &
        +2.d0*k(1041)*n(idx_NH)  &
        +k(275)  &
        +k(1048)*n(idx_O)  &
        +k(764)*n(idx_NHj)  &
        +k(871)*n(idx_C)  &
        +k(803)*n(idx_NH2j)  &
        +k(1036)*n(idx_CN)  &
        +k(796)*n(idx_COj)  &
        +k(797)*n(idx_H2COj)  &
        -k(1019)*n(idx_N)  &
        +k(1049)*n(idx_OH)  &
        +k(798)*n(idx_H2Oj)

    !d[NH_dot]/d[NH]
    pd(31,31) =  &
        -k(806)*n(idx_O2Hj)  &
        -k(871)*n(idx_C)  &
        -k(1038)*n(idx_NH3)  &
        -k(799)*n(idx_HCNj)  &
        -k(801)*n(idx_HNOj)  &
        -k(1043)*n(idx_NO)  &
        -k(1046)*n(idx_O2)  &
        -4.d0*k(1039)*n(idx_NH)  &
        -k(800)*n(idx_HCOj)  &
        -k(1037)*n(idx_H2O)  &
        -4.d0*k(1040)*n(idx_NH)  &
        -k(1035)*n(idx_CH4)  &
        -k(87)*n(idx_Hj)  &
        -k(595)*n(idx_H3j)  &
        -k(1042)*n(idx_NO2)  &
        -k(1045)*n(idx_O2)  &
        -k(798)*n(idx_H2Oj)  &
        -k(276)  &
        -k(1163)  &
        -k(1051)*n(idx_OH)  &
        -k(200)*n(idx_CNj)  &
        -k(1019)*n(idx_N)  &
        -k(1048)*n(idx_O)  &
        -4.d0*k(1041)*n(idx_NH)  &
        -k(112)*n(idx_H2j)  &
        -k(795)*n(idx_CH3j)  &
        -k(986)*n(idx_H)  &
        -k(1266)  &
        -k(727)*n(idx_Nj)  &
        -k(805)*n(idx_O2j)  &
        -k(1164)  &
        -k(803)*n(idx_NH2j)  &
        -k(376)*n(idx_Cj)  &
        -k(1050)*n(idx_OH)  &
        -k(524)*n(idx_H2j)  &
        -k(796)*n(idx_COj)  &
        -k(201)*n(idx_COj)  &
        -k(411)*n(idx_CHj)  &
        -k(1036)*n(idx_CN)  &
        -k(804)*n(idx_Oj)  &
        -k(1044)*n(idx_NO)  &
        -k(1047)*n(idx_O)  &
        -k(203)*n(idx_Oj)  &
        -k(764)*n(idx_NHj)  &
        -k(202)*n(idx_N2j)  &
        -k(275)  &
        -k(870)*n(idx_C)  &
        -k(1049)*n(idx_OH)  &
        -k(167)*n(idx_Nj)  &
        -k(694)*n(idx_HEj)  &
        -k(807)*n(idx_OHj)  &
        -k(802)*n(idx_N2Hj)  &
        -k(963)*n(idx_H2)  &
        -k(797)*n(idx_H2COj)

    !d[HE_dot]/d[NH]
    pd(35,31) =  &
        +k(694)*n(idx_HEj)

    !d[HNO_dot]/d[NH]
    pd(36,31) =  &
        +k(1042)*n(idx_NO2)  &
        +k(1050)*n(idx_OH)  &
        +k(1045)*n(idx_O2)

    !d[NO2_dot]/d[NH]
    pd(42,31) =  &
        -k(1042)*n(idx_NO2)

    !d[NH3_DUST_dot]/d[NH]
    pd(60,31) =  &
        +k(1266)

    !d[HCO+_dot]/d[NH]
    pd(70,31) =  &
        +k(796)*n(idx_COj)  &
        -k(800)*n(idx_HCOj)

    !d[H+_dot]/d[NH]
    pd(71,31) =  &
        -k(87)*n(idx_Hj)

    !d[C+_dot]/d[NH]
    pd(73,31) =  &
        -k(376)*n(idx_Cj)

    !d[CH+_dot]/d[NH]
    pd(75,31) =  &
        -k(411)*n(idx_CHj)

    !d[H2CO+_dot]/d[NH]
    pd(76,31) =  &
        -k(797)*n(idx_H2COj)

    !d[NH3+_dot]/d[NH]
    pd(78,31) =  &
        +k(803)*n(idx_NH2j)

    !d[NO+_dot]/d[NH]
    pd(79,31) =  &
        +k(804)*n(idx_Oj)

    !d[CN+_dot]/d[NH]
    pd(86,31) =  &
        +k(376)*n(idx_Cj)  &
        -k(200)*n(idx_CNj)  &
        +k(411)*n(idx_CHj)

    !d[CO+_dot]/d[NH]
    pd(87,31) =  &
        -k(796)*n(idx_COj)  &
        -k(201)*n(idx_COj)

    !d[N2+_dot]/d[NH]
    pd(88,31) =  &
        +k(727)*n(idx_Nj)  &
        -k(202)*n(idx_N2j)

    !d[O2+_dot]/d[NH]
    pd(89,31) =  &
        -k(805)*n(idx_O2j)

    !d[H2O+_dot]/d[NH]
    pd(90,31) =  &
        -k(798)*n(idx_H2Oj)

    !d[NH2+_dot]/d[NH]
    pd(91,31) =  &
        +k(802)*n(idx_N2Hj)  &
        +k(800)*n(idx_HCOj)  &
        +k(595)*n(idx_H3j)  &
        +k(801)*n(idx_HNOj)  &
        +k(807)*n(idx_OHj)  &
        +k(799)*n(idx_HCNj)  &
        +k(764)*n(idx_NHj)  &
        -k(803)*n(idx_NH2j)  &
        +k(806)*n(idx_O2Hj)  &
        +k(524)*n(idx_H2j)

    !d[O+_dot]/d[NH]
    pd(92,31) =  &
        -k(804)*n(idx_Oj)  &
        -k(203)*n(idx_Oj)

    !d[OH+_dot]/d[NH]
    pd(93,31) =  &
        -k(807)*n(idx_OHj)

    !d[CH3+_dot]/d[NH]
    pd(94,31) =  &
        -k(795)*n(idx_CH3j)

    !d[N+_dot]/d[NH]
    pd(96,31) =  &
        +k(694)*n(idx_HEj)  &
        -k(167)*n(idx_Nj)  &
        -k(727)*n(idx_Nj)

    !d[HCN+_dot]/d[NH]
    pd(97,31) =  &
        -k(799)*n(idx_HCNj)

    !d[NH+_dot]/d[NH]
    pd(98,31) =  &
        +k(202)*n(idx_N2j)  &
        +k(201)*n(idx_COj)  &
        +k(276)  &
        +k(167)*n(idx_Nj)  &
        -k(764)*n(idx_NHj)  &
        +k(112)*n(idx_H2j)  &
        +k(87)*n(idx_Hj)  &
        +k(203)*n(idx_Oj)  &
        +k(200)*n(idx_CNj)  &
        +k(1164)

    !d[H2+_dot]/d[NH]
    pd(102,31) =  &
        -k(524)*n(idx_H2j)  &
        -k(112)*n(idx_H2j)

    !d[HE+_dot]/d[NH]
    pd(103,31) =  &
        -k(694)*n(idx_HEj)

    !d[HNO+_dot]/d[NH]
    pd(104,31) =  &
        +k(805)*n(idx_O2j)  &
        -k(801)*n(idx_HNOj)

    !d[H3+_dot]/d[NH]
    pd(106,31) =  &
        -k(595)*n(idx_H3j)

    !d[H3CO+_dot]/d[NH]
    pd(107,31) =  &
        +k(797)*n(idx_H2COj)

    !d[H3O+_dot]/d[NH]
    pd(108,31) =  &
        +k(798)*n(idx_H2Oj)

    !d[HCNH+_dot]/d[NH]
    pd(109,31) =  &
        +k(795)*n(idx_CH3j)

    !d[N2H+_dot]/d[NH]
    pd(112,31) =  &
        -k(802)*n(idx_N2Hj)

    !d[O2H+_dot]/d[NH]
    pd(113,31) =  &
        -k(806)*n(idx_O2Hj)

    !d[O_dot]/d[SIH4]
    pd(3,32) =  &
        -k(1089)*n(idx_O)

    !d[HCN_dot]/d[SIH4]
    pd(5,32) =  &
        +k(951)*n(idx_CN)

    !d[H2_dot]/d[SIH4]
    pd(6,32) =  &
        +k(1186)  &
        +k(1188)  &
        +k(508)*n(idx_Hj)  &
        +k(709)*n(idx_HEj)  &
        +2.d0*k(708)*n(idx_HEj)  &
        +k(292)  &
        +k(605)*n(idx_H3j)

    !d[H_dot]/d[SIH4]
    pd(8,32) =  &
        +k(98)*n(idx_Hj)  &
        +k(709)*n(idx_HEj)  &
        +k(1187)  &
        +k(1188)

    !d[OH_dot]/d[SIH4]
    pd(10,32) =  &
        +k(1089)*n(idx_O)

    !d[SIH2_dot]/d[SIH4]
    pd(22,32) =  &
        +k(1186)  &
        +k(292)

    !d[SIH3_dot]/d[SIH4]
    pd(23,32) =  &
        +k(1187)  &
        +k(1089)*n(idx_O)  &
        +k(951)*n(idx_CN)

    !d[CN_dot]/d[SIH4]
    pd(24,32) =  &
        -k(951)*n(idx_CN)

    !d[CO_dot]/d[SIH4]
    pd(25,32) =  &
        +k(639)*n(idx_HCOj)

    !d[CH4_dot]/d[SIH4]
    pd(29,32) =  &
        +k(447)*n(idx_CH3j)

    !d[SIH4_dot]/d[SIH4]
    pd(32,32) =  &
        -k(951)*n(idx_CN)  &
        -k(1187)  &
        -k(98)*n(idx_Hj)  &
        -k(708)*n(idx_HEj)  &
        -k(292)  &
        -k(1188)  &
        -k(1238)  &
        -k(1186)  &
        -k(447)*n(idx_CH3j)  &
        -k(709)*n(idx_HEj)  &
        -k(508)*n(idx_Hj)  &
        -k(639)*n(idx_HCOj)  &
        -k(1089)*n(idx_O)  &
        -k(605)*n(idx_H3j)

    !d[SIH_dot]/d[SIH4]
    pd(33,32) =  &
        +k(1188)

    !d[HE_dot]/d[SIH4]
    pd(35,32) =  &
        +k(708)*n(idx_HEj)  &
        +k(709)*n(idx_HEj)

    !d[SIH4_DUST_dot]/d[SIH4]
    pd(48,32) =  &
        +k(1238)

    !d[HCO+_dot]/d[SIH4]
    pd(70,32) =  &
        -k(639)*n(idx_HCOj)

    !d[H+_dot]/d[SIH4]
    pd(71,32) =  &
        -k(508)*n(idx_Hj)  &
        -k(98)*n(idx_Hj)

    !d[SI+_dot]/d[SIH4]
    pd(80,32) =  &
        +k(708)*n(idx_HEj)

    !d[SIH3+_dot]/d[SIH4]
    pd(85,32) =  &
        +k(508)*n(idx_Hj)  &
        +k(447)*n(idx_CH3j)

    !d[CH3+_dot]/d[SIH4]
    pd(94,32) =  &
        -k(447)*n(idx_CH3j)

    !d[SIH4+_dot]/d[SIH4]
    pd(99,32) =  &
        +k(98)*n(idx_Hj)

    !d[SIH+_dot]/d[SIH4]
    pd(100,32) =  &
        +k(709)*n(idx_HEj)

    !d[HE+_dot]/d[SIH4]
    pd(103,32) =  &
        -k(708)*n(idx_HEj)  &
        -k(709)*n(idx_HEj)

    !d[H3+_dot]/d[SIH4]
    pd(106,32) =  &
        -k(605)*n(idx_H3j)

    !d[SIH5+_dot]/d[SIH4]
    pd(114,32) =  &
        +k(605)*n(idx_H3j)  &
        +k(639)*n(idx_HCOj)

    !d[O_dot]/d[SIH]
    pd(3,33) =  &
        +k(850)*n(idx_OHj)  &
        -k(1090)*n(idx_O)

    !d[H2_dot]/d[SIH]
    pd(6,33) =  &
        +k(509)*n(idx_Hj)  &
        +k(606)*n(idx_H3j)

    !d[C_dot]/d[SIH]
    pd(7,33) =  &
        -k(878)*n(idx_C)

    !d[H_dot]/d[SIH]
    pd(8,33) =  &
        +k(878)*n(idx_C)  &
        +k(382)*n(idx_Cj)  &
        +k(1189)  &
        +k(1090)*n(idx_O)  &
        +k(710)*n(idx_HEj)  &
        +k(99)*n(idx_Hj)  &
        +k(293)

    !d[H2O_dot]/d[SIH]
    pd(9,33) =  &
        +k(613)*n(idx_H3Oj)

    !d[SI_dot]/d[SIH]
    pd(18,33) =  &
        +k(1189)  &
        +k(293)

    !d[SIC_dot]/d[SIH]
    pd(21,33) =  &
        +k(878)*n(idx_C)

    !d[CO_dot]/d[SIH]
    pd(25,33) =  &
        +k(640)*n(idx_HCOj)

    !d[SIH_dot]/d[SIH]
    pd(33,33) =  &
        -k(1231)  &
        -k(509)*n(idx_Hj)  &
        -k(293)  &
        -k(1189)  &
        -k(710)*n(idx_HEj)  &
        -k(99)*n(idx_Hj)  &
        -k(878)*n(idx_C)  &
        -k(640)*n(idx_HCOj)  &
        -k(613)*n(idx_H3Oj)  &
        -k(1090)*n(idx_O)  &
        -k(382)*n(idx_Cj)  &
        -k(850)*n(idx_OHj)  &
        -k(606)*n(idx_H3j)

    !d[SIO_dot]/d[SIH]
    pd(34,33) =  &
        +k(1090)*n(idx_O)

    !d[HE_dot]/d[SIH]
    pd(35,33) =  &
        +k(710)*n(idx_HEj)

    !d[SIH4_DUST_dot]/d[SIH]
    pd(48,33) =  &
        +k(1231)

    !d[HCO+_dot]/d[SIH]
    pd(70,33) =  &
        -k(640)*n(idx_HCOj)

    !d[H+_dot]/d[SIH]
    pd(71,33) =  &
        -k(99)*n(idx_Hj)  &
        -k(509)*n(idx_Hj)

    !d[C+_dot]/d[SIH]
    pd(73,33) =  &
        -k(382)*n(idx_Cj)

    !d[SI+_dot]/d[SIH]
    pd(80,33) =  &
        +k(509)*n(idx_Hj)  &
        +k(710)*n(idx_HEj)

    !d[SIC+_dot]/d[SIH]
    pd(83,33) =  &
        +k(382)*n(idx_Cj)

    !d[SIH2+_dot]/d[SIH]
    pd(84,33) =  &
        +k(640)*n(idx_HCOj)  &
        +k(850)*n(idx_OHj)  &
        +k(606)*n(idx_H3j)  &
        +k(613)*n(idx_H3Oj)

    !d[OH+_dot]/d[SIH]
    pd(93,33) =  &
        -k(850)*n(idx_OHj)

    !d[SIH+_dot]/d[SIH]
    pd(100,33) =  &
        +k(99)*n(idx_Hj)

    !d[HE+_dot]/d[SIH]
    pd(103,33) =  &
        -k(710)*n(idx_HEj)

    !d[H3+_dot]/d[SIH]
    pd(106,33) =  &
        -k(606)*n(idx_H3j)

    !d[H3O+_dot]/d[SIH]
    pd(108,33) =  &
        -k(613)*n(idx_H3Oj)

    !d[E_dot]/d[SIO]
    pd(1,34) =  &
        +k(1192)

    !d[O_dot]/d[SIO]
    pd(3,34) =  &
        +k(1191)  &
        +k(294)  &
        +k(711)*n(idx_HEj)  &
        +k(851)*n(idx_OHj)

    !d[H2_dot]/d[SIO]
    pd(6,34) =  &
        +k(607)*n(idx_H3j)

    !d[H_dot]/d[SIO]
    pd(8,34) =  &
        +k(100)*n(idx_Hj)

    !d[H2O_dot]/d[SIO]
    pd(9,34) =  &
        +k(614)*n(idx_H3Oj)

    !d[SI_dot]/d[SIO]
    pd(18,34) =  &
        +k(1191)  &
        +k(294)  &
        +k(712)*n(idx_HEj)

    !d[CO_dot]/d[SIO]
    pd(25,34) =  &
        +k(641)*n(idx_HCOj)  &
        +k(383)*n(idx_Cj)

    !d[SIO_dot]/d[SIO]
    pd(34,34) =  &
        -k(1192)  &
        -k(294)  &
        -k(1230)  &
        -k(711)*n(idx_HEj)  &
        -k(383)*n(idx_Cj)  &
        -k(641)*n(idx_HCOj)  &
        -k(851)*n(idx_OHj)  &
        -k(607)*n(idx_H3j)  &
        -k(1191)  &
        -k(614)*n(idx_H3Oj)  &
        -k(100)*n(idx_Hj)  &
        -k(712)*n(idx_HEj)

    !d[HE_dot]/d[SIO]
    pd(35,34) =  &
        +k(711)*n(idx_HEj)  &
        +k(712)*n(idx_HEj)

    !d[H2SIO_DUST_dot]/d[SIO]
    pd(49,34) =  &
        +k(1230)

    !d[HCO+_dot]/d[SIO]
    pd(70,34) =  &
        -k(641)*n(idx_HCOj)

    !d[H+_dot]/d[SIO]
    pd(71,34) =  &
        -k(100)*n(idx_Hj)

    !d[C+_dot]/d[SIO]
    pd(73,34) =  &
        -k(383)*n(idx_Cj)

    !d[SI+_dot]/d[SIO]
    pd(80,34) =  &
        +k(383)*n(idx_Cj)  &
        +k(711)*n(idx_HEj)

    !d[O+_dot]/d[SIO]
    pd(92,34) =  &
        +k(712)*n(idx_HEj)

    !d[OH+_dot]/d[SIO]
    pd(93,34) =  &
        -k(851)*n(idx_OHj)

    !d[SIO+_dot]/d[SIO]
    pd(101,34) =  &
        +k(1192)  &
        +k(100)*n(idx_Hj)

    !d[HE+_dot]/d[SIO]
    pd(103,34) =  &
        -k(712)*n(idx_HEj)  &
        -k(711)*n(idx_HEj)

    !d[H3+_dot]/d[SIO]
    pd(106,34) =  &
        -k(607)*n(idx_H3j)

    !d[H3O+_dot]/d[SIO]
    pd(108,34) =  &
        -k(614)*n(idx_H3Oj)

    !d[SIOH+_dot]/d[SIO]
    pd(115,34) =  &
        +k(641)*n(idx_HCOj)  &
        +k(607)*n(idx_H3j)  &
        +k(851)*n(idx_OHj)  &
        +k(614)*n(idx_H3Oj)

    !d[E_dot]/d[HE]
    pd(1,35) =  &
        +k(266)  &
        +k(238)

    !d[H_dot]/d[HE]
    pd(8,35) =  &
        +k(521)*n(idx_H2j)

    !d[HE_dot]/d[HE]
    pd(35,35) =  &
        -k(266)  &
        -k(1199)*n(idx_Hj)  &
        -k(521)*n(idx_H2j)  &
        -k(238)

    !d[H+_dot]/d[HE]
    pd(71,35) =  &
        -k(1199)*n(idx_Hj)

    !d[H2+_dot]/d[HE]
    pd(102,35) =  &
        -k(521)*n(idx_H2j)

    !d[HE+_dot]/d[HE]
    pd(103,35) =  &
        +k(266)  &
        +k(238)

    !d[HEH+_dot]/d[HE]
    pd(111,35) =  &
        +k(521)*n(idx_H2j)  &
        +k(1199)*n(idx_Hj)

    !d[CH_dot]/d[HNO]
    pd(2,36) =  &
        -k(927)*n(idx_CH)

    !d[O_dot]/d[HNO]
    pd(3,36) =  &
        -k(1070)*n(idx_O)  &
        +k(981)*n(idx_H)  &
        -k(1069)*n(idx_O)  &
        -k(1071)*n(idx_O)

    !d[HCN_dot]/d[HNO]
    pd(5,36) =  &
        +k(945)*n(idx_CN)

    !d[H2_dot]/d[HNO]
    pd(6,36) =  &
        +k(982)*n(idx_H)  &
        +k(591)*n(idx_H3j)  &
        +k(504)*n(idx_Hj)

    !d[H_dot]/d[HNO]
    pd(8,36) =  &
        -k(983)*n(idx_H)  &
        +k(1069)*n(idx_O)  &
        +k(687)*n(idx_HEj)  &
        +k(1154)  &
        -k(981)*n(idx_H)  &
        +k(265)  &
        -k(982)*n(idx_H)

    !d[H2O_dot]/d[HNO]
    pd(9,36) =  &
        +k(1098)*n(idx_OH)

    !d[OH_dot]/d[HNO]
    pd(10,36) =  &
        +k(1070)*n(idx_O)  &
        -k(1098)*n(idx_OH)  &
        +k(983)*n(idx_H)

    !d[O2_dot]/d[HNO]
    pd(11,36) =  &
        +k(1071)*n(idx_O)

    !d[CH2_dot]/d[HNO]
    pd(12,36) =  &
        +k(927)*n(idx_CH)  &
        -k(884)*n(idx_CH2)

    !d[H2CO_dot]/d[HNO]
    pd(13,36) =  &
        +k(1000)*n(idx_HCO)

    !d[HCO_dot]/d[HNO]
    pd(14,36) =  &
        -k(1000)*n(idx_HCO)

    !d[NO_dot]/d[HNO]
    pd(17,36) =  &
        +k(982)*n(idx_H)  &
        +k(1000)*n(idx_HCO)  &
        +k(1098)*n(idx_OH)  &
        +k(688)*n(idx_HEj)  &
        +k(907)*n(idx_CH3)  &
        +k(1154)  &
        +k(1070)*n(idx_O)  &
        +k(884)*n(idx_CH2)  &
        +k(927)*n(idx_CH)  &
        +k(265)  &
        +k(945)*n(idx_CN)  &
        +k(1018)*n(idx_N)

    !d[CN_dot]/d[HNO]
    pd(24,36) =  &
        -k(945)*n(idx_CN)

    !d[CO_dot]/d[HNO]
    pd(25,36) =  &
        -k(952)*n(idx_CO)

    !d[NH2_dot]/d[HNO]
    pd(27,36) =  &
        +k(981)*n(idx_H)

    !d[CH3_dot]/d[HNO]
    pd(28,36) =  &
        +k(884)*n(idx_CH2)  &
        -k(907)*n(idx_CH3)

    !d[CH4_dot]/d[HNO]
    pd(29,36) =  &
        +k(907)*n(idx_CH3)

    !d[N_dot]/d[HNO]
    pd(30,36) =  &
        -k(1018)*n(idx_N)

    !d[NH_dot]/d[HNO]
    pd(31,36) =  &
        +k(1071)*n(idx_O)  &
        +k(952)*n(idx_CO)  &
        +k(983)*n(idx_H)  &
        +k(1018)*n(idx_N)

    !d[HE_dot]/d[HNO]
    pd(35,36) =  &
        +k(688)*n(idx_HEj)  &
        +k(687)*n(idx_HEj)

    !d[HNO_dot]/d[HNO]
    pd(36,36) =  &
        -k(983)*n(idx_H)  &
        -k(952)*n(idx_CO)  &
        -k(1070)*n(idx_O)  &
        -k(688)*n(idx_HEj)  &
        -k(591)*n(idx_H3j)  &
        -k(945)*n(idx_CN)  &
        -k(1295)  &
        -k(1098)*n(idx_OH)  &
        -k(1154)  &
        -k(982)*n(idx_H)  &
        -k(1018)*n(idx_N)  &
        -k(1069)*n(idx_O)  &
        -k(504)*n(idx_Hj)  &
        -k(265)  &
        -k(927)*n(idx_CH)  &
        -k(884)*n(idx_CH2)  &
        -k(1071)*n(idx_O)  &
        -k(687)*n(idx_HEj)  &
        -k(1000)*n(idx_HCO)  &
        -k(981)*n(idx_H)  &
        -k(907)*n(idx_CH3)

    !d[CO2_dot]/d[HNO]
    pd(38,36) =  &
        +k(952)*n(idx_CO)

    !d[NO2_dot]/d[HNO]
    pd(42,36) =  &
        +k(1069)*n(idx_O)

    !d[HNO_DUST_dot]/d[HNO]
    pd(63,36) =  &
        +k(1295)

    !d[H+_dot]/d[HNO]
    pd(71,36) =  &
        +k(688)*n(idx_HEj)  &
        -k(504)*n(idx_Hj)

    !d[NO+_dot]/d[HNO]
    pd(79,36) =  &
        +k(687)*n(idx_HEj)  &
        +k(504)*n(idx_Hj)

    !d[HE+_dot]/d[HNO]
    pd(103,36) =  &
        -k(687)*n(idx_HEj)  &
        -k(688)*n(idx_HEj)

    !d[H2NO+_dot]/d[HNO]
    pd(105,36) =  &
        +k(591)*n(idx_H3j)

    !d[H3+_dot]/d[HNO]
    pd(106,36) =  &
        -k(591)*n(idx_H3j)

    !d[E_dot]/d[CH3OH]
    pd(1,37) =  &
        +k(1121)

    !d[CH_dot]/d[CH3OH]
    pd(2,37) =  &
        +k(366)*n(idx_Cj)

    !d[H2_dot]/d[CH3OH]
    pd(6,37) =  &
        +k(1120)  &
        +k(248)  &
        +2.d0*k(495)*n(idx_Hj)  &
        +k(494)*n(idx_Hj)  &
        +k(580)*n(idx_H3j)

    !d[H_dot]/d[CH3OH]
    pd(8,37) =  &
        +k(1121)  &
        +k(821)*n(idx_O2j)  &
        +k(715)*n(idx_Nj)  &
        +k(716)*n(idx_Nj)  &
        +k(713)*n(idx_Nj)

    !d[H2O_dot]/d[CH3OH]
    pd(9,37) =  &
        +k(493)*n(idx_Hj)  &
        +k(580)*n(idx_H3j)  &
        +k(809)*n(idx_Oj)

    !d[OH_dot]/d[CH3OH]
    pd(10,37) =  &
        +k(249)  &
        +k(658)*n(idx_HEj)  &
        +k(1122)  &
        +k(810)*n(idx_Oj)

    !d[O2_dot]/d[CH3OH]
    pd(11,37) =  &
        +k(821)*n(idx_O2j)

    !d[CH2_dot]/d[CH3OH]
    pd(12,37) =  &
        +k(398)*n(idx_CHj)

    !d[H2CO_dot]/d[CH3OH]
    pd(13,37) =  &
        +k(248)  &
        +k(1120)  &
        +k(397)*n(idx_CHj)

    !d[HCO_dot]/d[CH3OH]
    pd(14,37) =  &
        +k(367)*n(idx_Cj)

    !d[NO_dot]/d[CH3OH]
    pd(17,37) =  &
        +k(716)*n(idx_Nj)

    !d[CH3_dot]/d[CH3OH]
    pd(28,37) =  &
        +k(249)  &
        +k(657)*n(idx_HEj)  &
        +k(715)*n(idx_Nj)  &
        +k(1122)  &
        +k(861)*n(idx_SIj)

    !d[CH4_dot]/d[CH3OH]
    pd(29,37) =  &
        +k(440)*n(idx_CH3j)

    !d[NH_dot]/d[CH3OH]
    pd(31,37) =  &
        +k(713)*n(idx_Nj)  &
        +k(714)*n(idx_Nj)

    !d[HE_dot]/d[CH3OH]
    pd(35,37) =  &
        +k(657)*n(idx_HEj)  &
        +k(658)*n(idx_HEj)

    !d[CH3OH_dot]/d[CH3OH]
    pd(37,37) =  &
        -k(716)*n(idx_Nj)  &
        -k(495)*n(idx_Hj)  &
        -k(809)*n(idx_Oj)  &
        -k(810)*n(idx_Oj)  &
        -k(249)  &
        -k(1224)  &
        -k(494)*n(idx_Hj)  &
        -k(1120)  &
        -k(1122)  &
        -k(398)*n(idx_CHj)  &
        -k(714)*n(idx_Nj)  &
        -k(248)  &
        -k(657)*n(idx_HEj)  &
        -k(397)*n(idx_CHj)  &
        -k(658)*n(idx_HEj)  &
        -k(366)*n(idx_Cj)  &
        -k(821)*n(idx_O2j)  &
        -k(440)*n(idx_CH3j)  &
        -k(713)*n(idx_Nj)  &
        -k(861)*n(idx_SIj)  &
        -k(1121)  &
        -k(367)*n(idx_Cj)  &
        -k(493)*n(idx_Hj)  &
        -k(715)*n(idx_Nj)  &
        -k(580)*n(idx_H3j)

    !d[CH3OH_DUST_dot]/d[CH3OH]
    pd(45,37) =  &
        +k(1224)

    !d[HCO+_dot]/d[CH3OH]
    pd(70,37) =  &
        +k(495)*n(idx_Hj)

    !d[H+_dot]/d[CH3OH]
    pd(71,37) =  &
        -k(495)*n(idx_Hj)  &
        -k(494)*n(idx_Hj)  &
        -k(493)*n(idx_Hj)

    !d[C+_dot]/d[CH3OH]
    pd(73,37) =  &
        -k(366)*n(idx_Cj)  &
        -k(367)*n(idx_Cj)

    !d[CH+_dot]/d[CH3OH]
    pd(75,37) =  &
        -k(397)*n(idx_CHj)  &
        -k(398)*n(idx_CHj)

    !d[H2CO+_dot]/d[CH3OH]
    pd(76,37) =  &
        +k(809)*n(idx_Oj)  &
        +k(713)*n(idx_Nj)

    !d[NO+_dot]/d[CH3OH]
    pd(79,37) =  &
        +k(715)*n(idx_Nj)

    !d[SI+_dot]/d[CH3OH]
    pd(80,37) =  &
        -k(861)*n(idx_SIj)

    !d[O2+_dot]/d[CH3OH]
    pd(89,37) =  &
        -k(821)*n(idx_O2j)

    !d[O+_dot]/d[CH3OH]
    pd(92,37) =  &
        -k(810)*n(idx_Oj)  &
        -k(809)*n(idx_Oj)

    !d[OH+_dot]/d[CH3OH]
    pd(93,37) =  &
        +k(657)*n(idx_HEj)

    !d[CH3+_dot]/d[CH3OH]
    pd(94,37) =  &
        +k(493)*n(idx_Hj)  &
        +k(367)*n(idx_Cj)  &
        -k(440)*n(idx_CH3j)  &
        +k(658)*n(idx_HEj)  &
        +k(580)*n(idx_H3j)  &
        +k(397)*n(idx_CHj)  &
        +k(716)*n(idx_Nj)

    !d[N+_dot]/d[CH3OH]
    pd(96,37) =  &
        -k(713)*n(idx_Nj)  &
        -k(715)*n(idx_Nj)  &
        -k(714)*n(idx_Nj)  &
        -k(716)*n(idx_Nj)

    !d[HE+_dot]/d[CH3OH]
    pd(103,37) =  &
        -k(657)*n(idx_HEj)  &
        -k(658)*n(idx_HEj)

    !d[H3+_dot]/d[CH3OH]
    pd(106,37) =  &
        -k(580)*n(idx_H3j)

    !d[H3CO+_dot]/d[CH3OH]
    pd(107,37) =  &
        +k(714)*n(idx_Nj)  &
        +k(1121)  &
        +k(440)*n(idx_CH3j)  &
        +k(494)*n(idx_Hj)  &
        +k(398)*n(idx_CHj)  &
        +k(366)*n(idx_Cj)  &
        +k(821)*n(idx_O2j)  &
        +k(810)*n(idx_Oj)

    !d[SIOH+_dot]/d[CH3OH]
    pd(115,37) =  &
        +k(861)*n(idx_SIj)

    !d[CH_dot]/d[CO2]
    pd(2,38) =  &
        -k(924)*n(idx_CH)

    !d[O_dot]/d[CO2]
    pd(3,38) =  &
        +k(666)*n(idx_HEj)  &
        -k(1060)*n(idx_O)  &
        +k(1133)  &
        +k(838)*n(idx_OHj)  &
        +k(497)*n(idx_Hj)  &
        +k(253)

    !d[H2_dot]/d[CO2]
    pd(6,38) =  &
        +k(583)*n(idx_H3j)

    !d[C_dot]/d[CO2]
    pd(7,38) =  &
        +k(668)*n(idx_HEj)

    !d[H_dot]/d[CO2]
    pd(8,38) =  &
        -k(972)*n(idx_H)  &
        +k(515)*n(idx_H2j)

    !d[OH_dot]/d[CO2]
    pd(10,38) =  &
        +k(972)*n(idx_H)

    !d[O2_dot]/d[CO2]
    pd(11,38) =  &
        +k(669)*n(idx_HEj)  &
        +k(1060)*n(idx_O)  &
        +k(822)*n(idx_O2Hj)

    !d[HCO_dot]/d[CO2]
    pd(14,38) =  &
        +k(924)*n(idx_CH)  &
        +k(751)*n(idx_NHj)

    !d[NO_dot]/d[CO2]
    pd(17,38) =  &
        +k(653)*n(idx_HNOj)  &
        +k(1013)*n(idx_N)  &
        +k(720)*n(idx_Nj)

    !d[SI_dot]/d[CO2]
    pd(18,38) =  &
        -k(1104)*n(idx_SI)

    !d[CN_dot]/d[CO2]
    pd(24,38) =  &
        +k(621)*n(idx_HCNj)

    !d[CO_dot]/d[CO2]
    pd(25,38) =  &
        +k(1060)*n(idx_O)  &
        +k(667)*n(idx_HEj)  &
        +k(750)*n(idx_NHj)  &
        +k(924)*n(idx_CH)  &
        +k(813)*n(idx_Oj)  &
        +k(399)*n(idx_CHj)  &
        +k(1104)*n(idx_SI)  &
        +k(972)*n(idx_H)  &
        +k(1133)  &
        +k(253)  &
        +k(1013)*n(idx_N)  &
        +k(368)*n(idx_Cj)  &
        +k(417)*n(idx_CH2j)

    !d[N2_dot]/d[CO2]
    pd(26,38) =  &
        +k(735)*n(idx_N2Hj)

    !d[CH3_dot]/d[CO2]
    pd(28,38) =  &
        +k(448)*n(idx_CH4j)

    !d[N_dot]/d[CO2]
    pd(30,38) =  &
        +k(749)*n(idx_NHj)  &
        -k(1013)*n(idx_N)

    !d[SIO_dot]/d[CO2]
    pd(34,38) =  &
        +k(1104)*n(idx_SI)

    !d[HE_dot]/d[CO2]
    pd(35,38) =  &
        +k(669)*n(idx_HEj)  &
        +k(668)*n(idx_HEj)  &
        +k(667)*n(idx_HEj)  &
        +k(666)*n(idx_HEj)

    !d[CO2_dot]/d[CO2]
    pd(38,38) =  &
        -k(621)*n(idx_HCNj)  &
        -k(253)  &
        -k(720)*n(idx_Nj)  &
        -k(669)*n(idx_HEj)  &
        -k(666)*n(idx_HEj)  &
        -k(399)*n(idx_CHj)  &
        -k(924)*n(idx_CH)  &
        -k(583)*n(idx_H3j)  &
        -k(1104)*n(idx_SI)  &
        -k(813)*n(idx_Oj)  &
        -k(497)*n(idx_Hj)  &
        -k(751)*n(idx_NHj)  &
        -k(822)*n(idx_O2Hj)  &
        -k(417)*n(idx_CH2j)  &
        -k(515)*n(idx_H2j)  &
        -k(368)*n(idx_Cj)  &
        -k(448)*n(idx_CH4j)  &
        -k(735)*n(idx_N2Hj)  &
        -k(838)*n(idx_OHj)  &
        -k(749)*n(idx_NHj)  &
        -k(667)*n(idx_HEj)  &
        -k(972)*n(idx_H)  &
        -k(1259)  &
        -k(653)*n(idx_HNOj)  &
        -k(1060)*n(idx_O)  &
        -k(750)*n(idx_NHj)  &
        -k(1133)  &
        -k(1013)*n(idx_N)  &
        -k(668)*n(idx_HEj)

    !d[CO2_DUST_dot]/d[CO2]
    pd(57,38) =  &
        +k(1259)

    !d[HCO+_dot]/d[CO2]
    pd(70,38) =  &
        +k(497)*n(idx_Hj)  &
        +k(399)*n(idx_CHj)

    !d[H+_dot]/d[CO2]
    pd(71,38) =  &
        -k(497)*n(idx_Hj)

    !d[C+_dot]/d[CO2]
    pd(73,38) =  &
        +k(669)*n(idx_HEj)  &
        -k(368)*n(idx_Cj)

    !d[CH2+_dot]/d[CO2]
    pd(74,38) =  &
        -k(417)*n(idx_CH2j)

    !d[CH+_dot]/d[CO2]
    pd(75,38) =  &
        -k(399)*n(idx_CHj)

    !d[H2CO+_dot]/d[CO2]
    pd(76,38) =  &
        +k(417)*n(idx_CH2j)

    !d[NO+_dot]/d[CO2]
    pd(79,38) =  &
        +k(751)*n(idx_NHj)

    !d[CO+_dot]/d[CO2]
    pd(87,38) =  &
        +k(720)*n(idx_Nj)  &
        +k(368)*n(idx_Cj)  &
        +k(666)*n(idx_HEj)

    !d[O2+_dot]/d[CO2]
    pd(89,38) =  &
        +k(668)*n(idx_HEj)  &
        +k(813)*n(idx_Oj)

    !d[O+_dot]/d[CO2]
    pd(92,38) =  &
        -k(813)*n(idx_Oj)  &
        +k(667)*n(idx_HEj)

    !d[OH+_dot]/d[CO2]
    pd(93,38) =  &
        -k(838)*n(idx_OHj)

    !d[CH4+_dot]/d[CO2]
    pd(95,38) =  &
        -k(448)*n(idx_CH4j)

    !d[N+_dot]/d[CO2]
    pd(96,38) =  &
        -k(720)*n(idx_Nj)

    !d[HCN+_dot]/d[CO2]
    pd(97,38) =  &
        -k(621)*n(idx_HCNj)

    !d[NH+_dot]/d[CO2]
    pd(98,38) =  &
        -k(751)*n(idx_NHj)  &
        -k(750)*n(idx_NHj)  &
        -k(749)*n(idx_NHj)

    !d[H2+_dot]/d[CO2]
    pd(102,38) =  &
        -k(515)*n(idx_H2j)

    !d[HE+_dot]/d[CO2]
    pd(103,38) =  &
        -k(667)*n(idx_HEj)  &
        -k(666)*n(idx_HEj)  &
        -k(669)*n(idx_HEj)  &
        -k(668)*n(idx_HEj)

    !d[HNO+_dot]/d[CO2]
    pd(104,38) =  &
        -k(653)*n(idx_HNOj)  &
        +k(750)*n(idx_NHj)

    !d[H3+_dot]/d[CO2]
    pd(106,38) =  &
        -k(583)*n(idx_H3j)

    !d[HCO2+_dot]/d[CO2]
    pd(110,38) =  &
        +k(735)*n(idx_N2Hj)  &
        +k(822)*n(idx_O2Hj)  &
        +k(653)*n(idx_HNOj)  &
        +k(749)*n(idx_NHj)  &
        +k(515)*n(idx_H2j)  &
        +k(838)*n(idx_OHj)  &
        +k(448)*n(idx_CH4j)  &
        +k(583)*n(idx_H3j)  &
        +k(621)*n(idx_HCNj)

    !d[N2H+_dot]/d[CO2]
    pd(112,38) =  &
        -k(735)*n(idx_N2Hj)

    !d[O2H+_dot]/d[CO2]
    pd(113,38) =  &
        -k(822)*n(idx_O2Hj)

    !d[O_dot]/d[H2CN]
    pd(3,39) =  &
        -k(1061)*n(idx_O)

    !d[HCN_dot]/d[H2CN]
    pd(5,39) =  &
        +k(255)  &
        +k(974)*n(idx_H)  &
        +k(1014)*n(idx_N)  &
        +k(1136)

    !d[H2_dot]/d[H2CN]
    pd(6,39) =  &
        +k(974)*n(idx_H)  &
        +k(1061)*n(idx_O)

    !d[H_dot]/d[H2CN]
    pd(8,39) =  &
        +k(255)  &
        -k(974)*n(idx_H)  &
        +k(1136)

    !d[N_dot]/d[H2CN]
    pd(30,39) =  &
        -k(1014)*n(idx_N)

    !d[NH_dot]/d[H2CN]
    pd(31,39) =  &
        +k(1014)*n(idx_N)

    !d[H2CN_dot]/d[H2CN]
    pd(39,39) =  &
        -k(974)*n(idx_H)  &
        -k(1299)  &
        -k(1014)*n(idx_N)  &
        -k(255)  &
        -k(1136)  &
        -k(1061)*n(idx_O)

    !d[OCN_dot]/d[H2CN]
    pd(44,39) =  &
        +k(1061)*n(idx_O)

    !d[H2CN_DUST_dot]/d[H2CN]
    pd(65,39) =  &
        +k(1299)

    !d[H2_dot]/d[H2SIO]
    pd(6,40) =  &
        +k(1144)  &
        +k(500)*n(idx_Hj)  &
        +k(258)

    !d[H_dot]/d[H2SIO]
    pd(8,40) =  &
        +k(676)*n(idx_HEj)  &
        +2.d0*k(1145)

    !d[SIO_dot]/d[H2SIO]
    pd(34,40) =  &
        +k(1144)  &
        +k(1145)  &
        +k(258)

    !d[HE_dot]/d[H2SIO]
    pd(35,40) =  &
        +k(676)*n(idx_HEj)

    !d[H2SIO_dot]/d[H2SIO]
    pd(40,40) =  &
        -k(1248)  &
        -k(1144)  &
        -k(1145)  &
        -k(500)*n(idx_Hj)  &
        -k(676)*n(idx_HEj)  &
        -k(258)

    !d[H2SIO_DUST_dot]/d[H2SIO]
    pd(49,40) =  &
        +k(1248)

    !d[H+_dot]/d[H2SIO]
    pd(71,40) =  &
        -k(500)*n(idx_Hj)

    !d[HE+_dot]/d[H2SIO]
    pd(103,40) =  &
        -k(676)*n(idx_HEj)

    !d[SIOH+_dot]/d[H2SIO]
    pd(115,40) =  &
        +k(676)*n(idx_HEj)  &
        +k(500)*n(idx_Hj)

    !d[HNC_dot]/d[HNCO]
    pd(4,41) =  &
        +k(1005)*n(idx_C)

    !d[C_dot]/d[HNCO]
    pd(7,41) =  &
        -k(1005)*n(idx_C)

    !d[CO_dot]/d[HNCO]
    pd(25,41) =  &
        +k(503)*n(idx_Hj)  &
        +k(264)  &
        +k(1153)  &
        +k(1005)*n(idx_C)

    !d[NH_dot]/d[HNCO]
    pd(31,41) =  &
        +k(264)  &
        +k(1153)

    !d[HNCO_dot]/d[HNCO]
    pd(41,41) =  &
        -k(1153)  &
        -k(1005)*n(idx_C)  &
        -k(503)*n(idx_Hj)  &
        -k(264)  &
        -k(1225)

    !d[HNCO_DUST_dot]/d[HNCO]
    pd(46,41) =  &
        +k(1225)

    !d[H+_dot]/d[HNCO]
    pd(71,41) =  &
        -k(503)*n(idx_Hj)

    !d[NH2+_dot]/d[HNCO]
    pd(91,41) =  &
        +k(503)*n(idx_Hj)

    !d[O_dot]/d[NO2]
    pd(3,42) =  &
        +k(277)  &
        +2.d0*k(1020)*n(idx_N)  &
        +k(1165)  &
        -k(1076)*n(idx_O)

    !d[H2_dot]/d[NO2]
    pd(6,42) =  &
        +k(596)*n(idx_H3j)

    !d[H_dot]/d[NO2]
    pd(8,42) =  &
        -k(987)*n(idx_H)

    !d[OH_dot]/d[NO2]
    pd(10,42) =  &
        +k(987)*n(idx_H)  &
        +k(596)*n(idx_H3j)  &
        +k(505)*n(idx_Hj)

    !d[O2_dot]/d[NO2]
    pd(11,42) =  &
        +k(819)*n(idx_Oj)  &
        +k(1022)*n(idx_N)  &
        +k(1076)*n(idx_O)

    !d[CH2_dot]/d[NO2]
    pd(12,42) =  &
        -k(886)*n(idx_CH2)

    !d[H2CO_dot]/d[NO2]
    pd(13,42) =  &
        +k(886)*n(idx_CH2)  &
        +k(910)*n(idx_CH3)

    !d[NO_dot]/d[NO2]
    pd(17,42) =  &
        +k(946)*n(idx_CN)  &
        +k(277)  &
        +k(987)*n(idx_H)  &
        +k(1076)*n(idx_O)  &
        +2.d0*k(1021)*n(idx_N)  &
        +k(953)*n(idx_CO)  &
        +k(886)*n(idx_CH2)  &
        +k(1042)*n(idx_NH)  &
        +k(1165)

    !d[CN_dot]/d[NO2]
    pd(24,42) =  &
        -k(946)*n(idx_CN)

    !d[CO_dot]/d[NO2]
    pd(25,42) =  &
        -k(953)*n(idx_CO)

    !d[N2_dot]/d[NO2]
    pd(26,42) =  &
        +k(1022)*n(idx_N)  &
        +k(1020)*n(idx_N)

    !d[CH3_dot]/d[NO2]
    pd(28,42) =  &
        -k(910)*n(idx_CH3)

    !d[N_dot]/d[NO2]
    pd(30,42) =  &
        -k(1020)*n(idx_N)  &
        -k(1021)*n(idx_N)  &
        -k(1022)*n(idx_N)

    !d[NH_dot]/d[NO2]
    pd(31,42) =  &
        -k(1042)*n(idx_NH)

    !d[HNO_dot]/d[NO2]
    pd(36,42) =  &
        +k(1042)*n(idx_NH)  &
        +k(910)*n(idx_CH3)

    !d[CO2_dot]/d[NO2]
    pd(38,42) =  &
        +k(953)*n(idx_CO)

    !d[NO2_dot]/d[NO2]
    pd(42,42) =  &
        -k(819)*n(idx_Oj)  &
        -k(505)*n(idx_Hj)  &
        -k(953)*n(idx_CO)  &
        -k(1022)*n(idx_N)  &
        -k(1076)*n(idx_O)  &
        -k(596)*n(idx_H3j)  &
        -k(1294)  &
        -k(946)*n(idx_CN)  &
        -k(1165)  &
        -k(1042)*n(idx_NH)  &
        -k(987)*n(idx_H)  &
        -k(277)  &
        -k(886)*n(idx_CH2)  &
        -k(1020)*n(idx_N)  &
        -k(1021)*n(idx_N)  &
        -k(910)*n(idx_CH3)

    !d[OCN_dot]/d[NO2]
    pd(44,42) =  &
        +k(946)*n(idx_CN)

    !d[NO2_DUST_dot]/d[NO2]
    pd(62,42) =  &
        +k(1294)

    !d[H+_dot]/d[NO2]
    pd(71,42) =  &
        -k(505)*n(idx_Hj)

    !d[NO+_dot]/d[NO2]
    pd(79,42) =  &
        +k(819)*n(idx_Oj)  &
        +k(596)*n(idx_H3j)  &
        +k(505)*n(idx_Hj)

    !d[O+_dot]/d[NO2]
    pd(92,42) =  &
        -k(819)*n(idx_Oj)

    !d[H3+_dot]/d[NO2]
    pd(106,42) =  &
        -k(596)*n(idx_H3j)

    !d[CH_dot]/d[O2H]
    pd(2,43) =  &
        -k(939)*n(idx_CH)  &
        -k(938)*n(idx_CH)

    !d[O_dot]/d[O2H]
    pd(3,43) =  &
        +k(1172)  &
        +k(991)*n(idx_H)  &
        -k(1078)*n(idx_O)

    !d[H2_dot]/d[O2H]
    pd(6,43) =  &
        +k(992)*n(idx_H)

    !d[H_dot]/d[O2H]
    pd(8,43) =  &
        +k(282)  &
        -k(991)*n(idx_H)  &
        -k(993)*n(idx_H)  &
        +k(1171)  &
        -k(992)*n(idx_H)

    !d[H2O_dot]/d[O2H]
    pd(9,43) =  &
        +k(1101)*n(idx_OH)  &
        +k(991)*n(idx_H)

    !d[OH_dot]/d[O2H]
    pd(10,43) =  &
        +k(955)*n(idx_CO)  &
        -k(1101)*n(idx_OH)  &
        +k(1078)*n(idx_O)  &
        +2.d0*k(993)*n(idx_H)  &
        +k(938)*n(idx_CH)  &
        +k(1172)

    !d[O2_dot]/d[O2H]
    pd(11,43) =  &
        +k(992)*n(idx_H)  &
        +k(915)*n(idx_CH3)  &
        +k(1078)*n(idx_O)  &
        +k(939)*n(idx_CH)  &
        +k(282)  &
        +k(1171)  &
        +k(1025)*n(idx_N)  &
        +k(1004)*n(idx_HCO)  &
        +k(1101)*n(idx_OH)

    !d[CH2_dot]/d[O2H]
    pd(12,43) =  &
        +k(939)*n(idx_CH)

    !d[H2CO_dot]/d[O2H]
    pd(13,43) =  &
        +k(1004)*n(idx_HCO)

    !d[HCO_dot]/d[O2H]
    pd(14,43) =  &
        +k(938)*n(idx_CH)  &
        -k(1004)*n(idx_HCO)

    !d[CO_dot]/d[O2H]
    pd(25,43) =  &
        -k(955)*n(idx_CO)

    !d[CH3_dot]/d[O2H]
    pd(28,43) =  &
        -k(915)*n(idx_CH3)

    !d[CH4_dot]/d[O2H]
    pd(29,43) =  &
        +k(915)*n(idx_CH3)

    !d[N_dot]/d[O2H]
    pd(30,43) =  &
        -k(1025)*n(idx_N)

    !d[NH_dot]/d[O2H]
    pd(31,43) =  &
        +k(1025)*n(idx_N)

    !d[CO2_dot]/d[O2H]
    pd(38,43) =  &
        +k(955)*n(idx_CO)

    !d[O2H_dot]/d[O2H]
    pd(43,43) =  &
        -k(1101)*n(idx_OH)  &
        -k(939)*n(idx_CH)  &
        -k(955)*n(idx_CO)  &
        -k(991)*n(idx_H)  &
        -k(1004)*n(idx_HCO)  &
        -k(1025)*n(idx_N)  &
        -k(993)*n(idx_H)  &
        -k(1078)*n(idx_O)  &
        -k(1172)  &
        -k(915)*n(idx_CH3)  &
        -k(938)*n(idx_CH)  &
        -k(1171)  &
        -k(1304)  &
        -k(282)  &
        -k(992)*n(idx_H)

    !d[O2H_DUST_dot]/d[O2H]
    pd(64,43) =  &
        +k(1304)

    !d[O_dot]/d[OCN]
    pd(3,44) =  &
        -k(1079)*n(idx_O)  &
        -k(1080)*n(idx_O)  &
        +k(284)  &
        +k(698)*n(idx_HEj)  &
        +k(1173)  &
        +k(994)*n(idx_H)

    !d[HCN_dot]/d[OCN]
    pd(5,44) =  &
        +k(994)*n(idx_H)

    !d[C_dot]/d[OCN]
    pd(7,44) =  &
        -k(875)*n(idx_C)

    !d[H_dot]/d[OCN]
    pd(8,44) =  &
        -k(994)*n(idx_H)  &
        -k(996)*n(idx_H)  &
        -k(995)*n(idx_H)

    !d[OH_dot]/d[OCN]
    pd(10,44) =  &
        +k(996)*n(idx_H)

    !d[O2_dot]/d[OCN]
    pd(11,44) =  &
        -k(1056)*n(idx_O2)  &
        +k(1080)*n(idx_O)  &
        -k(1055)*n(idx_O2)

    !d[NO_dot]/d[OCN]
    pd(17,44) =  &
        +k(1055)*n(idx_O2)  &
        -k(1054)*n(idx_NO)  &
        +k(1079)*n(idx_O)

    !d[CN_dot]/d[OCN]
    pd(24,44) =  &
        +k(1173)  &
        +k(284)  &
        +k(379)*n(idx_Cj)  &
        +k(996)*n(idx_H)  &
        +k(875)*n(idx_C)  &
        +k(699)*n(idx_HEj)  &
        +k(1080)*n(idx_O)

    !d[CO_dot]/d[OCN]
    pd(25,44) =  &
        +k(1079)*n(idx_O)  &
        +k(875)*n(idx_C)  &
        +k(1056)*n(idx_O2)  &
        +k(995)*n(idx_H)

    !d[N2_dot]/d[OCN]
    pd(26,44) =  &
        +k(1054)*n(idx_NO)

    !d[NH_dot]/d[OCN]
    pd(31,44) =  &
        +k(995)*n(idx_H)

    !d[HE_dot]/d[OCN]
    pd(35,44) =  &
        +k(698)*n(idx_HEj)  &
        +k(699)*n(idx_HEj)

    !d[CO2_dot]/d[OCN]
    pd(38,44) =  &
        +k(1055)*n(idx_O2)  &
        +k(1054)*n(idx_NO)

    !d[NO2_dot]/d[OCN]
    pd(42,44) =  &
        +k(1056)*n(idx_O2)

    !d[OCN_dot]/d[OCN]
    pd(44,44) =  &
        -k(284)  &
        -k(994)*n(idx_H)  &
        -k(1079)*n(idx_O)  &
        -k(875)*n(idx_C)  &
        -k(1226)  &
        -k(1080)*n(idx_O)  &
        -k(699)*n(idx_HEj)  &
        -k(1173)  &
        -k(1056)*n(idx_O2)  &
        -k(698)*n(idx_HEj)  &
        -k(1055)*n(idx_O2)  &
        -k(379)*n(idx_Cj)  &
        -k(996)*n(idx_H)  &
        -k(995)*n(idx_H)  &
        -k(1054)*n(idx_NO)

    !d[HNCO_DUST_dot]/d[OCN]
    pd(46,44) =  &
        +k(1226)

    !d[C+_dot]/d[OCN]
    pd(73,44) =  &
        -k(379)*n(idx_Cj)

    !d[CN+_dot]/d[OCN]
    pd(86,44) =  &
        +k(698)*n(idx_HEj)

    !d[CO+_dot]/d[OCN]
    pd(87,44) =  &
        +k(379)*n(idx_Cj)

    !d[O+_dot]/d[OCN]
    pd(92,44) =  &
        +k(699)*n(idx_HEj)

    !d[HE+_dot]/d[OCN]
    pd(103,44) =  &
        -k(698)*n(idx_HEj)  &
        -k(699)*n(idx_HEj)

    !d[CH3OH_dot]/d[CH3OH_DUST]
    pd(37,45) =  &
        +k(1321)

    !d[CH3OH_DUST_dot]/d[CH3OH_DUST]
    pd(45,45) =  &
        -k(1321)

    !d[HNCO_dot]/d[HNCO_DUST]
    pd(41,46) =  &
        +k(1323)

    !d[HNCO_DUST_dot]/d[HNCO_DUST]
    pd(46,46) =  &
        -k(1323)

    !d[H2CO_dot]/d[H2CO_DUST]
    pd(13,47) =  &
        +k(1318)

    !d[H2CO_DUST_dot]/d[H2CO_DUST]
    pd(47,47) =  &
        -k(1318)

    !d[SIH4_dot]/d[SIH4_DUST]
    pd(32,48) =  &
        +k(1326)

    !d[SIH4_DUST_dot]/d[SIH4_DUST]
    pd(48,48) =  &
        -k(1326)

    !d[H2SIO_dot]/d[H2SIO_DUST]
    pd(40,49) =  &
        +k(1329)

    !d[H2SIO_DUST_dot]/d[H2SIO_DUST]
    pd(49,49) =  &
        -k(1329)

    !d[SIC_dot]/d[SIC_DUST]
    pd(21,50) =  &
        +k(1327)

    !d[SIC_DUST_dot]/d[SIC_DUST]
    pd(50,50) =  &
        -k(1327)

    !d[SIC2_dot]/d[SIC2_DUST]
    pd(19,51) =  &
        +k(1330)

    !d[SIC2_DUST_dot]/d[SIC2_DUST]
    pd(51,51) =  &
        -k(1330)

    !d[SIC3_dot]/d[SIC3_DUST]
    pd(20,52) =  &
        +k(1331)

    !d[SIC3_DUST_dot]/d[SIC3_DUST]
    pd(52,52) =  &
        -k(1331)

    !d[CH4_dot]/d[CH4_DUST]
    pd(29,53) =  &
        +k(1308)

    !d[CH4_DUST_dot]/d[CH4_DUST]
    pd(53,53) =  &
        -k(1308)

    !d[CO_dot]/d[CO_DUST]
    pd(25,54) =  &
        +k(1314)

    !d[CO_DUST_dot]/d[CO_DUST]
    pd(54,54) =  &
        -k(1314)

    !d[H2O_dot]/d[H2O_DUST]
    pd(9,55) =  &
        +k(1310)

    !d[H2O_DUST_dot]/d[H2O_DUST]
    pd(55,55) =  &
        -k(1310)

    !d[NO_dot]/d[NO_DUST]
    pd(17,56) =  &
        +k(1317)

    !d[NO_DUST_dot]/d[NO_DUST]
    pd(56,56) =  &
        -k(1317)

    !d[CO2_dot]/d[CO2_DUST]
    pd(38,57) =  &
        +k(1324)

    !d[CO2_DUST_dot]/d[CO2_DUST]
    pd(57,57) =  &
        -k(1324)

    !d[N2_dot]/d[N2_DUST]
    pd(26,58) =  &
        +k(1315)

    !d[N2_DUST_dot]/d[N2_DUST]
    pd(58,58) =  &
        -k(1315)

    !d[HCN_dot]/d[HCN_DUST]
    pd(5,59) =  &
        +k(1312)

    !d[HCN_DUST_dot]/d[HCN_DUST]
    pd(59,59) =  &
        -k(1312)

    !d[NH3_dot]/d[NH3_DUST]
    pd(16,60) =  &
        +k(1309)

    !d[NH3_DUST_dot]/d[NH3_DUST]
    pd(60,60) =  &
        -k(1309)

    !d[O2_dot]/d[O2_DUST]
    pd(11,61) =  &
        +k(1320)

    !d[O2_DUST_dot]/d[O2_DUST]
    pd(61,61) =  &
        -k(1320)

    !d[NO2_dot]/d[NO2_DUST]
    pd(42,62) =  &
        +k(1325)

    !d[NO2_DUST_dot]/d[NO2_DUST]
    pd(62,62) =  &
        -k(1325)

    !d[HNO_dot]/d[HNO_DUST]
    pd(36,63) =  &
        +k(1319)

    !d[HNO_DUST_dot]/d[HNO_DUST]
    pd(63,63) =  &
        -k(1319)

    !d[O2H_dot]/d[O2H_DUST]
    pd(43,64) =  &
        +k(1322)

    !d[O2H_DUST_dot]/d[O2H_DUST]
    pd(64,64) =  &
        -k(1322)

    !d[H2CN_dot]/d[H2CN_DUST]
    pd(39,65) =  &
        +k(1316)

    !d[H2CN_DUST_dot]/d[H2CN_DUST]
    pd(65,65) =  &
        -k(1316)

    !d[MG_dot]/d[MG_DUST]
    pd(15,66) =  &
        +k(1311)

    !d[MG_DUST_dot]/d[MG_DUST]
    pd(66,66) =  &
        -k(1311)

    !d[HNC_dot]/d[HNC_DUST]
    pd(4,67) =  &
        +k(1313)

    !d[HNC_DUST_dot]/d[HNC_DUST]
    pd(67,67) =  &
        -k(1313)

    !d[SIO_dot]/d[SIO_DUST]
    pd(34,69) =  &
        +k(1328)

    !d[SIO_DUST_dot]/d[SIO_DUST]
    pd(69,69) =  &
        -k(1328)

    !d[E_dot]/d[HCO+]
    pd(1,70) =  &
        -k(331)*n(idx_E)

    !d[CH_dot]/d[HCO+]
    pd(2,70) =  &
        -k(467)*n(idx_CH)

    !d[HNC_dot]/d[HCO+]
    pd(4,70) =  &
        -k(649)*n(idx_HNC)

    !d[HCN_dot]/d[HCO+]
    pd(5,70) =  &
        -k(630)*n(idx_HCN)

    !d[C_dot]/d[HCO+]
    pd(7,70) =  &
        -k(387)*n(idx_C)

    !d[H_dot]/d[HCO+]
    pd(8,70) =  &
        +k(331)*n(idx_E)  &
        +k(856)*n(idx_OH)  &
        +k(1149)

    !d[H2O_dot]/d[HCO+]
    pd(9,70) =  &
        -k(567)*n(idx_H2O)

    !d[OH_dot]/d[HCO+]
    pd(10,70) =  &
        -k(855)*n(idx_OH)  &
        -k(856)*n(idx_OH)

    !d[CH2_dot]/d[HCO+]
    pd(12,70) =  &
        -k(430)*n(idx_CH2)

    !d[H2CO_dot]/d[HCO+]
    pd(13,70) =  &
        -k(636)*n(idx_H2CO)

    !d[HCO_dot]/d[HCO+]
    pd(14,70) =  &
        +k(150)*n(idx_MG)  &
        -k(637)*n(idx_HCO)

    !d[MG_dot]/d[HCO+]
    pd(15,70) =  &
        -k(150)*n(idx_MG)

    !d[SI_dot]/d[HCO+]
    pd(18,70) =  &
        -k(862)*n(idx_SI)

    !d[SIH2_dot]/d[HCO+]
    pd(22,70) =  &
        -k(638)*n(idx_SIH2)

    !d[CO_dot]/d[HCO+]
    pd(25,70) =  &
        +k(430)*n(idx_CH2)  &
        +k(630)*n(idx_HCN)  &
        +k(862)*n(idx_SI)  &
        +k(649)*n(idx_HNC)  &
        +k(567)*n(idx_H2O)  &
        +k(387)*n(idx_C)  &
        +k(467)*n(idx_CH)  &
        +k(638)*n(idx_SIH2)  &
        +k(855)*n(idx_OH)  &
        +k(331)*n(idx_E)  &
        +k(800)*n(idx_NH)  &
        +k(636)*n(idx_H2CO)  &
        +k(637)*n(idx_HCO)  &
        +k(788)*n(idx_NH2)  &
        +k(640)*n(idx_SIH)  &
        +k(641)*n(idx_SIO)  &
        +k(639)*n(idx_SIH4)

    !d[NH2_dot]/d[HCO+]
    pd(27,70) =  &
        -k(788)*n(idx_NH2)

    !d[NH_dot]/d[HCO+]
    pd(31,70) =  &
        -k(800)*n(idx_NH)

    !d[SIH4_dot]/d[HCO+]
    pd(32,70) =  &
        -k(639)*n(idx_SIH4)

    !d[SIH_dot]/d[HCO+]
    pd(33,70) =  &
        -k(640)*n(idx_SIH)

    !d[SIO_dot]/d[HCO+]
    pd(34,70) =  &
        -k(641)*n(idx_SIO)

    !d[H2CO_DUST_dot]/d[HCO+]
    pd(47,70) =  &
        +k(1282)

    !d[HCO+_dot]/d[HCO+]
    pd(70,70) =  &
        -k(638)*n(idx_SIH2)  &
        -k(640)*n(idx_SIH)  &
        -k(862)*n(idx_SI)  &
        -k(430)*n(idx_CH2)  &
        -k(150)*n(idx_MG)  &
        -k(636)*n(idx_H2CO)  &
        -k(855)*n(idx_OH)  &
        -k(630)*n(idx_HCN)  &
        -k(567)*n(idx_H2O)  &
        -k(649)*n(idx_HNC)  &
        -k(641)*n(idx_SIO)  &
        -k(331)*n(idx_E)  &
        -k(800)*n(idx_NH)  &
        -k(856)*n(idx_OH)  &
        -k(467)*n(idx_CH)  &
        -k(637)*n(idx_HCO)  &
        -k(387)*n(idx_C)  &
        -k(1149)  &
        -k(788)*n(idx_NH2)  &
        -k(639)*n(idx_SIH4)  &
        -k(1282)

    !d[CH2+_dot]/d[HCO+]
    pd(74,70) =  &
        +k(467)*n(idx_CH)

    !d[CH+_dot]/d[HCO+]
    pd(75,70) =  &
        +k(387)*n(idx_C)

    !d[H2CO+_dot]/d[HCO+]
    pd(76,70) =  &
        +k(637)*n(idx_HCO)

    !d[MG+_dot]/d[HCO+]
    pd(77,70) =  &
        +k(150)*n(idx_MG)

    !d[NH3+_dot]/d[HCO+]
    pd(78,70) =  &
        +k(788)*n(idx_NH2)

    !d[SIH2+_dot]/d[HCO+]
    pd(84,70) =  &
        +k(640)*n(idx_SIH)

    !d[SIH3+_dot]/d[HCO+]
    pd(85,70) =  &
        +k(638)*n(idx_SIH2)

    !d[CO+_dot]/d[HCO+]
    pd(87,70) =  &
        +k(1149)

    !d[H2O+_dot]/d[HCO+]
    pd(90,70) =  &
        +k(855)*n(idx_OH)

    !d[NH2+_dot]/d[HCO+]
    pd(91,70) =  &
        +k(800)*n(idx_NH)

    !d[CH3+_dot]/d[HCO+]
    pd(94,70) =  &
        +k(430)*n(idx_CH2)

    !d[SIH+_dot]/d[HCO+]
    pd(100,70) =  &
        +k(862)*n(idx_SI)

    !d[H3CO+_dot]/d[HCO+]
    pd(107,70) =  &
        +k(636)*n(idx_H2CO)

    !d[H3O+_dot]/d[HCO+]
    pd(108,70) =  &
        +k(567)*n(idx_H2O)

    !d[HCNH+_dot]/d[HCO+]
    pd(109,70) =  &
        +k(630)*n(idx_HCN)  &
        +k(649)*n(idx_HNC)

    !d[HCO2+_dot]/d[HCO+]
    pd(110,70) =  &
        +k(856)*n(idx_OH)

    !d[SIH5+_dot]/d[HCO+]
    pd(114,70) =  &
        +k(639)*n(idx_SIH4)

    !d[SIOH+_dot]/d[HCO+]
    pd(115,70) =  &
        +k(641)*n(idx_SIO)

    !d[E_dot]/d[H+]
    pd(1,71) =  &
        -k(1217)*n(idx_E)

    !d[CH_dot]/d[H+]
    pd(2,71) =  &
        -k(79)*n(idx_CH)

    !d[O_dot]/d[H+]
    pd(3,71) =  &
        +k(497)*n(idx_CO2)  &
        -k(90)*n(idx_O)

    !d[HNC_dot]/d[H+]
    pd(4,71) =  &
        -k(2)*n(idx_HNC)

    !d[HCN_dot]/d[H+]
    pd(5,71) =  &
        -k(82)*n(idx_HCN)  &
        +k(2)*n(idx_HNC)

    !d[H2_dot]/d[H+]
    pd(6,71) =  &
        +k(492)*n(idx_CH2)  &
        +2.d0*k(495)*n(idx_CH3OH)  &
        +k(508)*n(idx_SIH4)  &
        +k(506)*n(idx_SIH2)  &
        +k(501)*n(idx_HCO)  &
        +k(496)*n(idx_CH4)  &
        +k(504)*n(idx_HNO)  &
        +k(509)*n(idx_SIH)  &
        +k(500)*n(idx_H2SIO)  &
        +k(507)*n(idx_SIH3)  &
        +k(498)*n(idx_H2CO)  &
        +k(499)*n(idx_H2CO)  &
        +k(494)*n(idx_CH3OH)

    !d[H_dot]/d[H+]
    pd(8,71) =  &
        +k(83)*n(idx_HCO)  &
        +k(79)*n(idx_CH)  &
        -k(1198)*n(idx_H)  &
        +k(78)*n(idx_CH4)  &
        +k(98)*n(idx_SIH4)  &
        +k(88)*n(idx_NO)  &
        +k(95)*n(idx_SIC)  &
        +k(100)*n(idx_SIO)  &
        +k(498)*n(idx_H2CO)  &
        +k(80)*n(idx_H2CO)  &
        +k(96)*n(idx_SIH2)  &
        +k(84)*n(idx_MG)  &
        +k(89)*n(idx_O2)  &
        +k(97)*n(idx_SIH3)  &
        +k(86)*n(idx_NH3)  &
        +k(87)*n(idx_NH)  &
        +k(1217)*n(idx_E)  &
        +k(85)*n(idx_NH2)  &
        +k(93)*n(idx_SIC2)  &
        +k(82)*n(idx_HCN)  &
        +k(99)*n(idx_SIH)  &
        +k(94)*n(idx_SIC3)  &
        +k(92)*n(idx_SI)  &
        +k(90)*n(idx_O)  &
        +k(91)*n(idx_OH)  &
        +k(76)*n(idx_CH2)  &
        +k(77)*n(idx_CH3)  &
        +k(81)*n(idx_H2O)

    !d[H2O_dot]/d[H+]
    pd(9,71) =  &
        -k(81)*n(idx_H2O)  &
        +k(493)*n(idx_CH3OH)

    !d[OH_dot]/d[H+]
    pd(10,71) =  &
        +k(505)*n(idx_NO2)  &
        -k(91)*n(idx_OH)

    !d[O2_dot]/d[H+]
    pd(11,71) =  &
        -k(89)*n(idx_O2)

    !d[CH2_dot]/d[H+]
    pd(12,71) =  &
        -k(492)*n(idx_CH2)  &
        -k(76)*n(idx_CH2)

    !d[H2CO_dot]/d[H+]
    pd(13,71) =  &
        -k(498)*n(idx_H2CO)  &
        -k(499)*n(idx_H2CO)  &
        -k(80)*n(idx_H2CO)

    !d[HCO_dot]/d[H+]
    pd(14,71) =  &
        -k(502)*n(idx_HCO)  &
        -k(501)*n(idx_HCO)  &
        -k(83)*n(idx_HCO)

    !d[MG_dot]/d[H+]
    pd(15,71) =  &
        -k(84)*n(idx_MG)

    !d[NH3_dot]/d[H+]
    pd(16,71) =  &
        -k(86)*n(idx_NH3)

    !d[NO_dot]/d[H+]
    pd(17,71) =  &
        -k(88)*n(idx_NO)

    !d[SI_dot]/d[H+]
    pd(18,71) =  &
        -k(92)*n(idx_SI)

    !d[SIC2_dot]/d[H+]
    pd(19,71) =  &
        -k(93)*n(idx_SIC2)

    !d[SIC3_dot]/d[H+]
    pd(20,71) =  &
        -k(94)*n(idx_SIC3)

    !d[SIC_dot]/d[H+]
    pd(21,71) =  &
        -k(95)*n(idx_SIC)

    !d[SIH2_dot]/d[H+]
    pd(22,71) =  &
        -k(506)*n(idx_SIH2)  &
        -k(96)*n(idx_SIH2)

    !d[SIH3_dot]/d[H+]
    pd(23,71) =  &
        -k(97)*n(idx_SIH3)  &
        -k(507)*n(idx_SIH3)

    !d[CO_dot]/d[H+]
    pd(25,71) =  &
        +k(502)*n(idx_HCO)  &
        +k(503)*n(idx_HNCO)

    !d[NH2_dot]/d[H+]
    pd(27,71) =  &
        -k(85)*n(idx_NH2)

    !d[CH3_dot]/d[H+]
    pd(28,71) =  &
        -k(77)*n(idx_CH3)

    !d[CH4_dot]/d[H+]
    pd(29,71) =  &
        -k(78)*n(idx_CH4)  &
        -k(496)*n(idx_CH4)

    !d[NH_dot]/d[H+]
    pd(31,71) =  &
        -k(87)*n(idx_NH)

    !d[SIH4_dot]/d[H+]
    pd(32,71) =  &
        -k(98)*n(idx_SIH4)  &
        -k(508)*n(idx_SIH4)

    !d[SIH_dot]/d[H+]
    pd(33,71) =  &
        -k(509)*n(idx_SIH)  &
        -k(99)*n(idx_SIH)

    !d[SIO_dot]/d[H+]
    pd(34,71) =  &
        -k(100)*n(idx_SIO)

    !d[HE_dot]/d[H+]
    pd(35,71) =  &
        -k(1199)*n(idx_HE)

    !d[HNO_dot]/d[H+]
    pd(36,71) =  &
        -k(504)*n(idx_HNO)

    !d[CH3OH_dot]/d[H+]
    pd(37,71) =  &
        -k(495)*n(idx_CH3OH)  &
        -k(493)*n(idx_CH3OH)  &
        -k(494)*n(idx_CH3OH)

    !d[CO2_dot]/d[H+]
    pd(38,71) =  &
        -k(497)*n(idx_CO2)

    !d[H2SIO_dot]/d[H+]
    pd(40,71) =  &
        -k(500)*n(idx_H2SIO)

    !d[HNCO_dot]/d[H+]
    pd(41,71) =  &
        -k(503)*n(idx_HNCO)

    !d[NO2_dot]/d[H+]
    pd(42,71) =  &
        -k(505)*n(idx_NO2)

    !d[HCO+_dot]/d[H+]
    pd(70,71) =  &
        +k(83)*n(idx_HCO)  &
        +k(495)*n(idx_CH3OH)  &
        +k(499)*n(idx_H2CO)  &
        +k(497)*n(idx_CO2)

    !d[H+_dot]/d[H+]
    pd(71,71) =  &
        -k(1198)*n(idx_H)  &
        -k(498)*n(idx_H2CO)  &
        -k(84)*n(idx_MG)  &
        -k(86)*n(idx_NH3)  &
        -k(493)*n(idx_CH3OH)  &
        -k(504)*n(idx_HNO)  &
        -k(87)*n(idx_NH)  &
        -k(509)*n(idx_SIH)  &
        -k(99)*n(idx_SIH)  &
        -k(92)*n(idx_SI)  &
        -k(79)*n(idx_CH)  &
        -k(500)*n(idx_H2SIO)  &
        -k(2)*n(idx_HNC)  &
        -k(89)*n(idx_O2)  &
        -k(497)*n(idx_CO2)  &
        -k(503)*n(idx_HNCO)  &
        -k(78)*n(idx_CH4)  &
        -k(85)*n(idx_NH2)  &
        -k(501)*n(idx_HCO)  &
        -k(495)*n(idx_CH3OH)  &
        -k(1217)*n(idx_E)  &
        -k(506)*n(idx_SIH2)  &
        -k(91)*n(idx_OH)  &
        -k(80)*n(idx_H2CO)  &
        -k(82)*n(idx_HCN)  &
        -k(95)*n(idx_SIC)  &
        -k(76)*n(idx_CH2)  &
        -k(96)*n(idx_SIH2)  &
        -k(81)*n(idx_H2O)  &
        -k(492)*n(idx_CH2)  &
        -k(508)*n(idx_SIH4)  &
        -k(505)*n(idx_NO2)  &
        -k(507)*n(idx_SIH3)  &
        -k(83)*n(idx_HCO)  &
        -k(502)*n(idx_HCO)  &
        -k(98)*n(idx_SIH4)  &
        -k(496)*n(idx_CH4)  &
        -k(1199)*n(idx_HE)  &
        -k(77)*n(idx_CH3)  &
        -k(100)*n(idx_SIO)  &
        +k(2)*n(idx_HNC)  &
        -k(88)*n(idx_NO)  &
        -k(94)*n(idx_SIC3)  &
        -k(494)*n(idx_CH3OH)  &
        -k(90)*n(idx_O)  &
        -k(93)*n(idx_SIC2)  &
        -k(97)*n(idx_SIH3)  &
        -k(499)*n(idx_H2CO)

    !d[CH2+_dot]/d[H+]
    pd(74,71) =  &
        +k(76)*n(idx_CH2)

    !d[CH+_dot]/d[H+]
    pd(75,71) =  &
        +k(492)*n(idx_CH2)  &
        +k(79)*n(idx_CH)

    !d[H2CO+_dot]/d[H+]
    pd(76,71) =  &
        +k(80)*n(idx_H2CO)

    !d[MG+_dot]/d[H+]
    pd(77,71) =  &
        +k(84)*n(idx_MG)

    !d[NH3+_dot]/d[H+]
    pd(78,71) =  &
        +k(86)*n(idx_NH3)

    !d[NO+_dot]/d[H+]
    pd(79,71) =  &
        +k(88)*n(idx_NO)  &
        +k(504)*n(idx_HNO)  &
        +k(505)*n(idx_NO2)

    !d[SI+_dot]/d[H+]
    pd(80,71) =  &
        +k(92)*n(idx_SI)  &
        +k(509)*n(idx_SIH)

    !d[SIC2+_dot]/d[H+]
    pd(81,71) =  &
        +k(93)*n(idx_SIC2)

    !d[SIC3+_dot]/d[H+]
    pd(82,71) =  &
        +k(94)*n(idx_SIC3)

    !d[SIC+_dot]/d[H+]
    pd(83,71) =  &
        +k(95)*n(idx_SIC)

    !d[SIH2+_dot]/d[H+]
    pd(84,71) =  &
        +k(96)*n(idx_SIH2)  &
        +k(507)*n(idx_SIH3)

    !d[SIH3+_dot]/d[H+]
    pd(85,71) =  &
        +k(508)*n(idx_SIH4)  &
        +k(97)*n(idx_SIH3)

    !d[CO+_dot]/d[H+]
    pd(87,71) =  &
        +k(498)*n(idx_H2CO)  &
        +k(501)*n(idx_HCO)

    !d[O2+_dot]/d[H+]
    pd(89,71) =  &
        +k(89)*n(idx_O2)

    !d[H2O+_dot]/d[H+]
    pd(90,71) =  &
        +k(81)*n(idx_H2O)

    !d[NH2+_dot]/d[H+]
    pd(91,71) =  &
        +k(503)*n(idx_HNCO)  &
        +k(85)*n(idx_NH2)

    !d[O+_dot]/d[H+]
    pd(92,71) =  &
        +k(90)*n(idx_O)

    !d[OH+_dot]/d[H+]
    pd(93,71) =  &
        +k(91)*n(idx_OH)

    !d[CH3+_dot]/d[H+]
    pd(94,71) =  &
        +k(496)*n(idx_CH4)  &
        +k(493)*n(idx_CH3OH)  &
        +k(77)*n(idx_CH3)

    !d[CH4+_dot]/d[H+]
    pd(95,71) =  &
        +k(78)*n(idx_CH4)

    !d[HCN+_dot]/d[H+]
    pd(97,71) =  &
        +k(82)*n(idx_HCN)

    !d[NH+_dot]/d[H+]
    pd(98,71) =  &
        +k(87)*n(idx_NH)

    !d[SIH4+_dot]/d[H+]
    pd(99,71) =  &
        +k(98)*n(idx_SIH4)

    !d[SIH+_dot]/d[H+]
    pd(100,71) =  &
        +k(99)*n(idx_SIH)  &
        +k(506)*n(idx_SIH2)

    !d[SIO+_dot]/d[H+]
    pd(101,71) =  &
        +k(100)*n(idx_SIO)

    !d[H2+_dot]/d[H+]
    pd(102,71) =  &
        +k(1198)*n(idx_H)  &
        +k(502)*n(idx_HCO)

    !d[H3CO+_dot]/d[H+]
    pd(107,71) =  &
        +k(494)*n(idx_CH3OH)

    !d[HEH+_dot]/d[H+]
    pd(111,71) =  &
        +k(1199)*n(idx_HE)

    !d[SIOH+_dot]/d[H+]
    pd(115,71) =  &
        +k(500)*n(idx_H2SIO)

    !d[E_dot]/d[HOC+]
    pd(1,72) =  &
        -k(336)*n(idx_E)

    !d[H2_dot]/d[HOC+]
    pd(6,72) =  &
        -k(6)*n(idx_H2)  &
        +k(6)*n(idx_H2)

    !d[H_dot]/d[HOC+]
    pd(8,72) =  &
        +k(336)*n(idx_E)

    !d[CO_dot]/d[HOC+]
    pd(25,72) =  &
        +k(336)*n(idx_E)

    !d[H2CO_DUST_dot]/d[HOC+]
    pd(47,72) =  &
        +k(1227)

    !d[HCO+_dot]/d[HOC+]
    pd(70,72) =  &
        +k(6)*n(idx_H2)

    !d[HOC+_dot]/d[HOC+]
    pd(72,72) =  &
        -k(1227)  &
        -k(6)*n(idx_H2)  &
        -k(336)*n(idx_E)

    !d[E_dot]/d[C+]
    pd(1,73) =  &
        -k(1215)*n(idx_E)

    !d[CH_dot]/d[C+]
    pd(2,73) =  &
        +k(366)*n(idx_CH3OH)  &
        -k(16)*n(idx_CH)  &
        +k(370)*n(idx_H2CO)

    !d[O_dot]/d[C+]
    pd(3,73) =  &
        +k(377)*n(idx_O2)  &
        -k(1194)*n(idx_O)

    !d[H2_dot]/d[C+]
    pd(6,73) =  &
        +k(381)*n(idx_SIH2)  &
        +k(375)*n(idx_NH3)  &
        -k(529)*n(idx_H2)  &
        -k(1200)*n(idx_H2)

    !d[C_dot]/d[C+]
    pd(7,73) =  &
        +k(24)*n(idx_SIC3)  &
        +k(23)*n(idx_SIC2)  &
        +k(1215)*n(idx_E)  &
        +k(25)*n(idx_SIC)  &
        +k(20)*n(idx_NH3)  &
        +k(21)*n(idx_NO)  &
        +k(26)*n(idx_SIH2)  &
        +k(19)*n(idx_MG)  &
        +k(22)*n(idx_SI)  &
        +k(17)*n(idx_H2CO)  &
        +k(18)*n(idx_HCO)  &
        +k(27)*n(idx_SIH3)  &
        +k(15)*n(idx_CH2)  &
        +k(16)*n(idx_CH)

    !d[H_dot]/d[C+]
    pd(8,73) =  &
        +k(374)*n(idx_NH2)  &
        +k(380)*n(idx_OH)  &
        +k(529)*n(idx_H2)  &
        +k(376)*n(idx_NH)  &
        -k(1206)*n(idx_H)  &
        +k(371)*n(idx_H2O)  &
        +k(382)*n(idx_SIH)  &
        +k(372)*n(idx_H2O)

    !d[H2O_dot]/d[C+]
    pd(9,73) =  &
        -k(371)*n(idx_H2O)  &
        -k(372)*n(idx_H2O)

    !d[OH_dot]/d[C+]
    pd(10,73) =  &
        -k(380)*n(idx_OH)

    !d[O2_dot]/d[C+]
    pd(11,73) =  &
        -k(378)*n(idx_O2)  &
        -k(377)*n(idx_O2)

    !d[CH2_dot]/d[C+]
    pd(12,73) =  &
        -k(15)*n(idx_CH2)

    !d[H2CO_dot]/d[C+]
    pd(13,73) =  &
        -k(370)*n(idx_H2CO)  &
        -k(17)*n(idx_H2CO)  &
        -k(369)*n(idx_H2CO)

    !d[HCO_dot]/d[C+]
    pd(14,73) =  &
        -k(373)*n(idx_HCO)  &
        +k(367)*n(idx_CH3OH)  &
        -k(18)*n(idx_HCO)

    !d[MG_dot]/d[C+]
    pd(15,73) =  &
        -k(19)*n(idx_MG)

    !d[NH3_dot]/d[C+]
    pd(16,73) =  &
        -k(375)*n(idx_NH3)  &
        -k(20)*n(idx_NH3)

    !d[NO_dot]/d[C+]
    pd(17,73) =  &
        -k(21)*n(idx_NO)

    !d[SI_dot]/d[C+]
    pd(18,73) =  &
        -k(22)*n(idx_SI)

    !d[SIC2_dot]/d[C+]
    pd(19,73) =  &
        -k(23)*n(idx_SIC2)

    !d[SIC3_dot]/d[C+]
    pd(20,73) =  &
        -k(24)*n(idx_SIC3)

    !d[SIC_dot]/d[C+]
    pd(21,73) =  &
        -k(25)*n(idx_SIC)

    !d[SIH2_dot]/d[C+]
    pd(22,73) =  &
        -k(381)*n(idx_SIH2)  &
        -k(26)*n(idx_SIH2)

    !d[SIH3_dot]/d[C+]
    pd(23,73) =  &
        -k(27)*n(idx_SIH3)

    !d[CN_dot]/d[C+]
    pd(24,73) =  &
        +k(379)*n(idx_OCN)

    !d[CO_dot]/d[C+]
    pd(25,73) =  &
        +k(378)*n(idx_O2)  &
        +k(369)*n(idx_H2CO)  &
        +k(383)*n(idx_SIO)  &
        +k(373)*n(idx_HCO)  &
        +k(368)*n(idx_CO2)

    !d[NH2_dot]/d[C+]
    pd(27,73) =  &
        -k(374)*n(idx_NH2)

    !d[N_dot]/d[C+]
    pd(30,73) =  &
        -k(1193)*n(idx_N)

    !d[NH_dot]/d[C+]
    pd(31,73) =  &
        -k(376)*n(idx_NH)

    !d[SIH_dot]/d[C+]
    pd(33,73) =  &
        -k(382)*n(idx_SIH)

    !d[SIO_dot]/d[C+]
    pd(34,73) =  &
        -k(383)*n(idx_SIO)

    !d[CH3OH_dot]/d[C+]
    pd(37,73) =  &
        -k(367)*n(idx_CH3OH)  &
        -k(366)*n(idx_CH3OH)

    !d[CO2_dot]/d[C+]
    pd(38,73) =  &
        -k(368)*n(idx_CO2)

    !d[OCN_dot]/d[C+]
    pd(44,73) =  &
        -k(379)*n(idx_OCN)

    !d[CH4_DUST_dot]/d[C+]
    pd(53,73) =  &
        +k(1265)

    !d[HCO+_dot]/d[C+]
    pd(70,73) =  &
        +k(371)*n(idx_H2O)  &
        +k(18)*n(idx_HCO)  &
        +k(370)*n(idx_H2CO)

    !d[HOC+_dot]/d[C+]
    pd(72,73) =  &
        +k(372)*n(idx_H2O)

    !d[C+_dot]/d[C+]
    pd(73,73) =  &
        -k(382)*n(idx_SIH)  &
        -k(374)*n(idx_NH2)  &
        -k(23)*n(idx_SIC2)  &
        -k(27)*n(idx_SIH3)  &
        -k(15)*n(idx_CH2)  &
        -k(377)*n(idx_O2)  &
        -k(25)*n(idx_SIC)  &
        -k(370)*n(idx_H2CO)  &
        -k(17)*n(idx_H2CO)  &
        -k(376)*n(idx_NH)  &
        -k(21)*n(idx_NO)  &
        -k(1194)*n(idx_O)  &
        -k(16)*n(idx_CH)  &
        -k(383)*n(idx_SIO)  &
        -k(1206)*n(idx_H)  &
        -k(379)*n(idx_OCN)  &
        -k(19)*n(idx_MG)  &
        -k(372)*n(idx_H2O)  &
        -k(18)*n(idx_HCO)  &
        -k(378)*n(idx_O2)  &
        -k(1200)*n(idx_H2)  &
        -k(381)*n(idx_SIH2)  &
        -k(1193)*n(idx_N)  &
        -k(371)*n(idx_H2O)  &
        -k(380)*n(idx_OH)  &
        -k(26)*n(idx_SIH2)  &
        -k(368)*n(idx_CO2)  &
        -k(20)*n(idx_NH3)  &
        -k(22)*n(idx_SI)  &
        -k(369)*n(idx_H2CO)  &
        -k(529)*n(idx_H2)  &
        -k(1265)  &
        -k(24)*n(idx_SIC3)  &
        -k(373)*n(idx_HCO)  &
        -k(375)*n(idx_NH3)  &
        -k(366)*n(idx_CH3OH)  &
        -k(367)*n(idx_CH3OH)  &
        -k(1215)*n(idx_E)

    !d[CH2+_dot]/d[C+]
    pd(74,73) =  &
        +k(15)*n(idx_CH2)  &
        +k(369)*n(idx_H2CO)  &
        +k(1200)*n(idx_H2)

    !d[CH+_dot]/d[C+]
    pd(75,73) =  &
        +k(529)*n(idx_H2)  &
        +k(1206)*n(idx_H)  &
        +k(373)*n(idx_HCO)  &
        +k(16)*n(idx_CH)

    !d[H2CO+_dot]/d[C+]
    pd(76,73) =  &
        +k(17)*n(idx_H2CO)

    !d[MG+_dot]/d[C+]
    pd(77,73) =  &
        +k(19)*n(idx_MG)

    !d[NH3+_dot]/d[C+]
    pd(78,73) =  &
        +k(20)*n(idx_NH3)

    !d[NO+_dot]/d[C+]
    pd(79,73) =  &
        +k(21)*n(idx_NO)

    !d[SI+_dot]/d[C+]
    pd(80,73) =  &
        +k(383)*n(idx_SIO)  &
        +k(22)*n(idx_SI)

    !d[SIC2+_dot]/d[C+]
    pd(81,73) =  &
        +k(23)*n(idx_SIC2)

    !d[SIC3+_dot]/d[C+]
    pd(82,73) =  &
        +k(24)*n(idx_SIC3)

    !d[SIC+_dot]/d[C+]
    pd(83,73) =  &
        +k(381)*n(idx_SIH2)  &
        +k(382)*n(idx_SIH)  &
        +k(25)*n(idx_SIC)

    !d[SIH2+_dot]/d[C+]
    pd(84,73) =  &
        +k(26)*n(idx_SIH2)

    !d[SIH3+_dot]/d[C+]
    pd(85,73) =  &
        +k(27)*n(idx_SIH3)

    !d[CN+_dot]/d[C+]
    pd(86,73) =  &
        +k(376)*n(idx_NH)  &
        +k(1193)*n(idx_N)

    !d[CO+_dot]/d[C+]
    pd(87,73) =  &
        +k(379)*n(idx_OCN)  &
        +k(380)*n(idx_OH)  &
        +k(1194)*n(idx_O)  &
        +k(377)*n(idx_O2)  &
        +k(368)*n(idx_CO2)

    !d[O+_dot]/d[C+]
    pd(92,73) =  &
        +k(378)*n(idx_O2)

    !d[CH3+_dot]/d[C+]
    pd(94,73) =  &
        +k(367)*n(idx_CH3OH)

    !d[HCN+_dot]/d[C+]
    pd(97,73) =  &
        +k(374)*n(idx_NH2)  &
        +k(375)*n(idx_NH3)

    !d[H3CO+_dot]/d[C+]
    pd(107,73) =  &
        +k(366)*n(idx_CH3OH)

    !d[E_dot]/d[CH2+]
    pd(1,74) =  &
        -k(298)*n(idx_E)  &
        -k(296)*n(idx_E)  &
        -k(297)*n(idx_E)

    !d[CH_dot]/d[CH2+]
    pd(2,74) =  &
        +k(1112)  &
        +k(298)*n(idx_E)

    !d[O_dot]/d[CH2+]
    pd(3,74) =  &
        -k(422)*n(idx_O)

    !d[H2_dot]/d[CH2+]
    pd(6,74) =  &
        -k(531)*n(idx_H2)  &
        +k(1110)  &
        +k(616)*n(idx_H)  &
        +k(296)*n(idx_E)

    !d[C_dot]/d[CH2+]
    pd(7,74) =  &
        +k(296)*n(idx_E)  &
        +k(297)*n(idx_E)

    !d[H_dot]/d[CH2+]
    pd(8,74) =  &
        +2.d0*k(297)*n(idx_E)  &
        +k(419)*n(idx_H2O)  &
        +k(298)*n(idx_E)  &
        -k(616)*n(idx_H)  &
        +k(422)*n(idx_O)  &
        +k(737)*n(idx_N)  &
        +k(531)*n(idx_H2)  &
        +k(1111)

    !d[H2O_dot]/d[CH2+]
    pd(9,74) =  &
        -k(419)*n(idx_H2O)

    !d[OH_dot]/d[CH2+]
    pd(10,74) =  &
        +k(421)*n(idx_O2)

    !d[O2_dot]/d[CH2+]
    pd(11,74) =  &
        -k(421)*n(idx_O2)

    !d[CH2_dot]/d[CH2+]
    pd(12,74) =  &
        +k(37)*n(idx_NO)

    !d[H2CO_dot]/d[CH2+]
    pd(13,74) =  &
        -k(418)*n(idx_H2CO)

    !d[HCO_dot]/d[CH2+]
    pd(14,74) =  &
        -k(420)*n(idx_HCO)

    !d[NO_dot]/d[CH2+]
    pd(17,74) =  &
        -k(37)*n(idx_NO)

    !d[CO_dot]/d[CH2+]
    pd(25,74) =  &
        +k(420)*n(idx_HCO)  &
        +k(417)*n(idx_CO2)

    !d[CH3_dot]/d[CH2+]
    pd(28,74) =  &
        +k(418)*n(idx_H2CO)

    !d[N_dot]/d[CH2+]
    pd(30,74) =  &
        -k(737)*n(idx_N)

    !d[CO2_dot]/d[CH2+]
    pd(38,74) =  &
        -k(417)*n(idx_CO2)

    !d[CH4_DUST_dot]/d[CH2+]
    pd(53,74) =  &
        +k(1279)

    !d[HCO+_dot]/d[CH2+]
    pd(70,74) =  &
        +k(418)*n(idx_H2CO)  &
        +k(422)*n(idx_O)  &
        +k(421)*n(idx_O2)

    !d[H+_dot]/d[CH2+]
    pd(71,74) =  &
        +k(1112)

    !d[C+_dot]/d[CH2+]
    pd(73,74) =  &
        +k(1110)

    !d[CH2+_dot]/d[CH2+]
    pd(74,74) =  &
        -k(417)*n(idx_CO2)  &
        -k(531)*n(idx_H2)  &
        -k(296)*n(idx_E)  &
        -k(1111)  &
        -k(298)*n(idx_E)  &
        -k(737)*n(idx_N)  &
        -k(1110)  &
        -k(297)*n(idx_E)  &
        -k(37)*n(idx_NO)  &
        -k(419)*n(idx_H2O)  &
        -k(420)*n(idx_HCO)  &
        -k(1112)  &
        -k(418)*n(idx_H2CO)  &
        -k(422)*n(idx_O)  &
        -k(616)*n(idx_H)  &
        -k(421)*n(idx_O2)  &
        -k(1279)

    !d[CH+_dot]/d[CH2+]
    pd(75,74) =  &
        +k(616)*n(idx_H)  &
        +k(1111)

    !d[H2CO+_dot]/d[CH2+]
    pd(76,74) =  &
        +k(417)*n(idx_CO2)

    !d[NO+_dot]/d[CH2+]
    pd(79,74) =  &
        +k(37)*n(idx_NO)

    !d[CH3+_dot]/d[CH2+]
    pd(94,74) =  &
        +k(420)*n(idx_HCO)  &
        +k(531)*n(idx_H2)

    !d[HCN+_dot]/d[CH2+]
    pd(97,74) =  &
        +k(737)*n(idx_N)

    !d[H3CO+_dot]/d[CH2+]
    pd(107,74) =  &
        +k(419)*n(idx_H2O)

    !d[E_dot]/d[CH+]
    pd(1,75) =  &
        -k(295)*n(idx_E)

    !d[CH_dot]/d[CH+]
    pd(2,75) =  &
        +k(33)*n(idx_MG)  &
        +k(36)*n(idx_SI)  &
        +k(34)*n(idx_NH3)  &
        +k(32)*n(idx_HCO)  &
        +k(35)*n(idx_NO)

    !d[O_dot]/d[CH+]
    pd(3,75) =  &
        -k(415)*n(idx_O)  &
        +k(413)*n(idx_O2)

    !d[HNC_dot]/d[CH+]
    pd(4,75) =  &
        -k(408)*n(idx_HNC)

    !d[HCN_dot]/d[CH+]
    pd(5,75) =  &
        -k(406)*n(idx_HCN)

    !d[H2_dot]/d[CH+]
    pd(6,75) =  &
        +k(405)*n(idx_H2O)  &
        +k(411)*n(idx_NH)  &
        +k(410)*n(idx_NH2)  &
        +k(416)*n(idx_OH)  &
        +k(615)*n(idx_H)  &
        -k(530)*n(idx_H2)

    !d[C_dot]/d[CH+]
    pd(7,75) =  &
        +k(1109)  &
        +k(406)*n(idx_HCN)  &
        +k(295)*n(idx_E)  &
        +k(404)*n(idx_H2O)  &
        +k(401)*n(idx_H2CO)  &
        +k(408)*n(idx_HNC)

    !d[H_dot]/d[CH+]
    pd(8,75) =  &
        +k(403)*n(idx_H2O)  &
        +k(295)*n(idx_E)  &
        +k(242)  &
        -k(615)*n(idx_H)  &
        +k(409)*n(idx_N)  &
        +k(415)*n(idx_O)  &
        +k(530)*n(idx_H2)

    !d[H2O_dot]/d[CH+]
    pd(9,75) =  &
        -k(403)*n(idx_H2O)  &
        -k(405)*n(idx_H2O)  &
        -k(404)*n(idx_H2O)

    !d[OH_dot]/d[CH+]
    pd(10,75) =  &
        +k(412)*n(idx_O2)  &
        -k(416)*n(idx_OH)

    !d[O2_dot]/d[CH+]
    pd(11,75) =  &
        -k(414)*n(idx_O2)  &
        -k(412)*n(idx_O2)  &
        -k(413)*n(idx_O2)

    !d[CH2_dot]/d[CH+]
    pd(12,75) =  &
        +k(398)*n(idx_CH3OH)  &
        +k(402)*n(idx_H2CO)

    !d[H2CO_dot]/d[CH+]
    pd(13,75) =  &
        +k(397)*n(idx_CH3OH)  &
        -k(401)*n(idx_H2CO)  &
        -k(400)*n(idx_H2CO)  &
        -k(402)*n(idx_H2CO)

    !d[HCO_dot]/d[CH+]
    pd(14,75) =  &
        -k(407)*n(idx_HCO)  &
        -k(32)*n(idx_HCO)  &
        +k(414)*n(idx_O2)

    !d[MG_dot]/d[CH+]
    pd(15,75) =  &
        -k(33)*n(idx_MG)

    !d[NH3_dot]/d[CH+]
    pd(16,75) =  &
        -k(34)*n(idx_NH3)

    !d[NO_dot]/d[CH+]
    pd(17,75) =  &
        -k(35)*n(idx_NO)

    !d[SI_dot]/d[CH+]
    pd(18,75) =  &
        -k(36)*n(idx_SI)

    !d[CO_dot]/d[CH+]
    pd(25,75) =  &
        +k(400)*n(idx_H2CO)  &
        +k(407)*n(idx_HCO)  &
        +k(399)*n(idx_CO2)

    !d[NH2_dot]/d[CH+]
    pd(27,75) =  &
        -k(410)*n(idx_NH2)

    !d[N_dot]/d[CH+]
    pd(30,75) =  &
        -k(409)*n(idx_N)

    !d[NH_dot]/d[CH+]
    pd(31,75) =  &
        -k(411)*n(idx_NH)

    !d[CH3OH_dot]/d[CH+]
    pd(37,75) =  &
        -k(398)*n(idx_CH3OH)  &
        -k(397)*n(idx_CH3OH)

    !d[CO2_dot]/d[CH+]
    pd(38,75) =  &
        -k(399)*n(idx_CO2)

    !d[CH4_DUST_dot]/d[CH+]
    pd(53,75) =  &
        +k(1273)

    !d[HCO+_dot]/d[CH+]
    pd(70,75) =  &
        +k(405)*n(idx_H2O)  &
        +k(402)*n(idx_H2CO)  &
        +k(32)*n(idx_HCO)  &
        +k(413)*n(idx_O2)  &
        +k(399)*n(idx_CO2)

    !d[H+_dot]/d[CH+]
    pd(71,75) =  &
        +k(1109)

    !d[C+_dot]/d[CH+]
    pd(73,75) =  &
        +k(615)*n(idx_H)  &
        +k(242)

    !d[CH2+_dot]/d[CH+]
    pd(74,75) =  &
        +k(530)*n(idx_H2)  &
        +k(407)*n(idx_HCO)

    !d[CH+_dot]/d[CH+]
    pd(75,75) =  &
        -k(295)*n(idx_E)  &
        -k(404)*n(idx_H2O)  &
        -k(1273)  &
        -k(415)*n(idx_O)  &
        -k(615)*n(idx_H)  &
        -k(409)*n(idx_N)  &
        -k(397)*n(idx_CH3OH)  &
        -k(35)*n(idx_NO)  &
        -k(403)*n(idx_H2O)  &
        -k(411)*n(idx_NH)  &
        -k(402)*n(idx_H2CO)  &
        -k(408)*n(idx_HNC)  &
        -k(399)*n(idx_CO2)  &
        -k(36)*n(idx_SI)  &
        -k(412)*n(idx_O2)  &
        -k(407)*n(idx_HCO)  &
        -k(33)*n(idx_MG)  &
        -k(416)*n(idx_OH)  &
        -k(400)*n(idx_H2CO)  &
        -k(32)*n(idx_HCO)  &
        -k(414)*n(idx_O2)  &
        -k(530)*n(idx_H2)  &
        -k(405)*n(idx_H2O)  &
        -k(1109)  &
        -k(401)*n(idx_H2CO)  &
        -k(34)*n(idx_NH3)  &
        -k(406)*n(idx_HCN)  &
        -k(413)*n(idx_O2)  &
        -k(398)*n(idx_CH3OH)  &
        -k(242)  &
        -k(410)*n(idx_NH2)

    !d[H2CO+_dot]/d[CH+]
    pd(76,75) =  &
        +k(403)*n(idx_H2O)

    !d[MG+_dot]/d[CH+]
    pd(77,75) =  &
        +k(33)*n(idx_MG)

    !d[NH3+_dot]/d[CH+]
    pd(78,75) =  &
        +k(34)*n(idx_NH3)

    !d[NO+_dot]/d[CH+]
    pd(79,75) =  &
        +k(35)*n(idx_NO)

    !d[SI+_dot]/d[CH+]
    pd(80,75) =  &
        +k(36)*n(idx_SI)

    !d[CN+_dot]/d[CH+]
    pd(86,75) =  &
        +k(411)*n(idx_NH)  &
        +k(409)*n(idx_N)

    !d[CO+_dot]/d[CH+]
    pd(87,75) =  &
        +k(415)*n(idx_O)  &
        +k(416)*n(idx_OH)  &
        +k(412)*n(idx_O2)

    !d[O+_dot]/d[CH+]
    pd(92,75) =  &
        +k(414)*n(idx_O2)

    !d[CH3+_dot]/d[CH+]
    pd(94,75) =  &
        +k(397)*n(idx_CH3OH)  &
        +k(400)*n(idx_H2CO)

    !d[HCN+_dot]/d[CH+]
    pd(97,75) =  &
        +k(410)*n(idx_NH2)

    !d[H3CO+_dot]/d[CH+]
    pd(107,75) =  &
        +k(398)*n(idx_CH3OH)  &
        +k(401)*n(idx_H2CO)

    !d[H3O+_dot]/d[CH+]
    pd(108,75) =  &
        +k(404)*n(idx_H2O)

    !d[HCNH+_dot]/d[CH+]
    pd(109,75) =  &
        +k(408)*n(idx_HNC)  &
        +k(406)*n(idx_HCN)

    !d[E_dot]/d[H2CO+]
    pd(1,76) =  &
        -k(307)*n(idx_E)  &
        -k(308)*n(idx_E)  &
        -k(309)*n(idx_E)  &
        -k(310)*n(idx_E)  &
        -k(1218)*n(idx_E)

    !d[CH_dot]/d[H2CO+]
    pd(2,76) =  &
        -k(56)*n(idx_CH)  &
        -k(460)*n(idx_CH)

    !d[O_dot]/d[H2CO+]
    pd(3,76) =  &
        +k(307)*n(idx_E)

    !d[HNC_dot]/d[H2CO+]
    pd(4,76) =  &
        -k(647)*n(idx_HNC)

    !d[HCN_dot]/d[H2CO+]
    pd(5,76) =  &
        -k(628)*n(idx_HCN)

    !d[H2_dot]/d[H2CO+]
    pd(6,76) =  &
        +k(308)*n(idx_E)

    !d[H_dot]/d[H2CO+]
    pd(8,76) =  &
        +2.d0*k(309)*n(idx_E)  &
        +k(310)*n(idx_E)

    !d[H2O_dot]/d[H2CO+]
    pd(9,76) =  &
        -k(564)*n(idx_H2O)

    !d[O2_dot]/d[H2CO+]
    pd(11,76) =  &
        -k(550)*n(idx_O2)

    !d[CH2_dot]/d[H2CO+]
    pd(12,76) =  &
        -k(424)*n(idx_CH2)  &
        +k(307)*n(idx_E)  &
        -k(40)*n(idx_CH2)

    !d[H2CO_dot]/d[H2CO+]
    pd(13,76) =  &
        +k(195)*n(idx_NH3)  &
        +k(1218)*n(idx_E)  &
        +k(40)*n(idx_CH2)  &
        -k(549)*n(idx_H2CO)  &
        +k(56)*n(idx_CH)  &
        +k(149)*n(idx_MG)  &
        +k(229)*n(idx_SI)  &
        +k(204)*n(idx_NO)  &
        +k(137)*n(idx_HCO)

    !d[HCO_dot]/d[H2CO+]
    pd(14,76) =  &
        +k(424)*n(idx_CH2)  &
        +k(781)*n(idx_NH2)  &
        +k(460)*n(idx_CH)  &
        -k(642)*n(idx_HCO)  &
        +k(628)*n(idx_HCN)  &
        +k(549)*n(idx_H2CO)  &
        -k(137)*n(idx_HCO)  &
        +k(564)*n(idx_H2O)  &
        +k(310)*n(idx_E)  &
        +k(647)*n(idx_HNC)

    !d[MG_dot]/d[H2CO+]
    pd(15,76) =  &
        -k(149)*n(idx_MG)

    !d[NH3_dot]/d[H2CO+]
    pd(16,76) =  &
        -k(195)*n(idx_NH3)

    !d[NO_dot]/d[H2CO+]
    pd(17,76) =  &
        -k(204)*n(idx_NO)

    !d[SI_dot]/d[H2CO+]
    pd(18,76) =  &
        -k(229)*n(idx_SI)

    !d[CO_dot]/d[H2CO+]
    pd(25,76) =  &
        +k(308)*n(idx_E)  &
        +k(309)*n(idx_E)  &
        +k(642)*n(idx_HCO)

    !d[NH2_dot]/d[H2CO+]
    pd(27,76) =  &
        -k(781)*n(idx_NH2)

    !d[CH3_dot]/d[H2CO+]
    pd(28,76) =  &
        +k(453)*n(idx_CH4)

    !d[CH4_dot]/d[H2CO+]
    pd(29,76) =  &
        -k(453)*n(idx_CH4)

    !d[N_dot]/d[H2CO+]
    pd(30,76) =  &
        +k(797)*n(idx_NH)

    !d[NH_dot]/d[H2CO+]
    pd(31,76) =  &
        -k(797)*n(idx_NH)

    !d[O2H_dot]/d[H2CO+]
    pd(43,76) =  &
        +k(550)*n(idx_O2)

    !d[H2CO_DUST_dot]/d[H2CO+]
    pd(47,76) =  &
        +k(1285)

    !d[HCO+_dot]/d[H2CO+]
    pd(70,76) =  &
        +k(550)*n(idx_O2)  &
        +k(137)*n(idx_HCO)

    !d[CH2+_dot]/d[H2CO+]
    pd(74,76) =  &
        +k(40)*n(idx_CH2)  &
        +k(460)*n(idx_CH)

    !d[CH+_dot]/d[H2CO+]
    pd(75,76) =  &
        +k(56)*n(idx_CH)

    !d[H2CO+_dot]/d[H2CO+]
    pd(76,76) =  &
        -k(642)*n(idx_HCO)  &
        -k(549)*n(idx_H2CO)  &
        -k(453)*n(idx_CH4)  &
        -k(308)*n(idx_E)  &
        -k(424)*n(idx_CH2)  &
        -k(309)*n(idx_E)  &
        -k(310)*n(idx_E)  &
        -k(195)*n(idx_NH3)  &
        -k(40)*n(idx_CH2)  &
        -k(628)*n(idx_HCN)  &
        -k(797)*n(idx_NH)  &
        -k(781)*n(idx_NH2)  &
        -k(137)*n(idx_HCO)  &
        -k(204)*n(idx_NO)  &
        -k(149)*n(idx_MG)  &
        -k(564)*n(idx_H2O)  &
        -k(307)*n(idx_E)  &
        -k(229)*n(idx_SI)  &
        -k(550)*n(idx_O2)  &
        -k(56)*n(idx_CH)  &
        -k(647)*n(idx_HNC)  &
        -k(1285)  &
        -k(460)*n(idx_CH)  &
        -k(1218)*n(idx_E)

    !d[MG+_dot]/d[H2CO+]
    pd(77,76) =  &
        +k(149)*n(idx_MG)

    !d[NH3+_dot]/d[H2CO+]
    pd(78,76) =  &
        +k(781)*n(idx_NH2)  &
        +k(195)*n(idx_NH3)

    !d[NO+_dot]/d[H2CO+]
    pd(79,76) =  &
        +k(204)*n(idx_NO)

    !d[SI+_dot]/d[H2CO+]
    pd(80,76) =  &
        +k(229)*n(idx_SI)

    !d[CH3+_dot]/d[H2CO+]
    pd(94,76) =  &
        +k(424)*n(idx_CH2)

    !d[H3CO+_dot]/d[H2CO+]
    pd(107,76) =  &
        +k(797)*n(idx_NH)  &
        +k(642)*n(idx_HCO)  &
        +k(549)*n(idx_H2CO)  &
        +k(453)*n(idx_CH4)

    !d[H3O+_dot]/d[H2CO+]
    pd(108,76) =  &
        +k(564)*n(idx_H2O)

    !d[HCNH+_dot]/d[H2CO+]
    pd(109,76) =  &
        +k(647)*n(idx_HNC)  &
        +k(628)*n(idx_HCN)

    !d[E_dot]/d[MG+]
    pd(1,77) =  &
        -k(1220)*n(idx_E)

    !d[MG_dot]/d[MG+]
    pd(15,77) =  &
        +k(1220)*n(idx_E)

    !d[MG_DUST_dot]/d[MG+]
    pd(66,77) =  &
        +k(1300)

    !d[MG+_dot]/d[MG+]
    pd(77,77) =  &
        -k(1220)*n(idx_E)  &
        -k(1300)

    !d[E_dot]/d[NH3+]
    pd(1,78) =  &
        -k(344)*n(idx_E)  &
        -k(345)*n(idx_E)

    !d[O_dot]/d[NH3+]
    pd(3,78) =  &
        -k(829)*n(idx_O)

    !d[H2_dot]/d[NH3+]
    pd(6,78) =  &
        +k(829)*n(idx_O)

    !d[H_dot]/d[NH3+]
    pd(8,78) =  &
        +k(344)*n(idx_E)  &
        +2.d0*k(345)*n(idx_E)

    !d[CH2_dot]/d[NH3+]
    pd(12,78) =  &
        -k(435)*n(idx_CH2)

    !d[HCO_dot]/d[NH3+]
    pd(14,78) =  &
        -k(190)*n(idx_HCO)

    !d[MG_dot]/d[NH3+]
    pd(15,78) =  &
        -k(191)*n(idx_MG)

    !d[NH3_dot]/d[NH3+]
    pd(16,78) =  &
        +k(193)*n(idx_SI)  &
        +k(190)*n(idx_HCO)  &
        +k(192)*n(idx_NO)  &
        +k(191)*n(idx_MG)

    !d[NO_dot]/d[NH3+]
    pd(17,78) =  &
        -k(192)*n(idx_NO)

    !d[SI_dot]/d[NH3+]
    pd(18,78) =  &
        -k(193)*n(idx_SI)

    !d[NH2_dot]/d[NH3+]
    pd(27,78) =  &
        +k(344)*n(idx_E)  &
        +k(435)*n(idx_CH2)

    !d[NH_dot]/d[NH3+]
    pd(31,78) =  &
        +k(345)*n(idx_E)

    !d[NH3_DUST_dot]/d[NH3+]
    pd(60,78) =  &
        +k(1284)

    !d[HCO+_dot]/d[NH3+]
    pd(70,78) =  &
        +k(190)*n(idx_HCO)

    !d[MG+_dot]/d[NH3+]
    pd(77,78) =  &
        +k(191)*n(idx_MG)

    !d[NH3+_dot]/d[NH3+]
    pd(78,78) =  &
        -k(191)*n(idx_MG)  &
        -k(344)*n(idx_E)  &
        -k(435)*n(idx_CH2)  &
        -k(829)*n(idx_O)  &
        -k(192)*n(idx_NO)  &
        -k(345)*n(idx_E)  &
        -k(1284)  &
        -k(190)*n(idx_HCO)  &
        -k(193)*n(idx_SI)

    !d[NO+_dot]/d[NH3+]
    pd(79,78) =  &
        +k(192)*n(idx_NO)

    !d[SI+_dot]/d[NH3+]
    pd(80,78) =  &
        +k(193)*n(idx_SI)

    !d[CH3+_dot]/d[NH3+]
    pd(94,78) =  &
        +k(435)*n(idx_CH2)

    !d[HNO+_dot]/d[NH3+]
    pd(104,78) =  &
        +k(829)*n(idx_O)

    !d[E_dot]/d[NO+]
    pd(1,79) =  &
        -k(346)*n(idx_E)

    !d[O_dot]/d[NO+]
    pd(3,79) =  &
        +k(346)*n(idx_E)

    !d[MG_dot]/d[NO+]
    pd(15,79) =  &
        -k(152)*n(idx_MG)

    !d[NO_dot]/d[NO+]
    pd(17,79) =  &
        +k(230)*n(idx_SI)  &
        +k(152)*n(idx_MG)

    !d[SI_dot]/d[NO+]
    pd(18,79) =  &
        -k(230)*n(idx_SI)

    !d[N_dot]/d[NO+]
    pd(30,79) =  &
        +k(346)*n(idx_E)

    !d[NO_DUST_dot]/d[NO+]
    pd(56,79) =  &
        +k(1278)

    !d[MG+_dot]/d[NO+]
    pd(77,79) =  &
        +k(152)*n(idx_MG)

    !d[NO+_dot]/d[NO+]
    pd(79,79) =  &
        -k(230)*n(idx_SI)  &
        -k(152)*n(idx_MG)  &
        -k(346)*n(idx_E)  &
        -k(1278)

    !d[SI+_dot]/d[NO+]
    pd(80,79) =  &
        +k(230)*n(idx_SI)

    !d[E_dot]/d[SI+]
    pd(1,80) =  &
        -k(1223)*n(idx_E)

    !d[CH_dot]/d[SI+]
    pd(2,80) =  &
        -k(477)*n(idx_CH)

    !d[O_dot]/d[SI+]
    pd(3,80) =  &
        -k(1213)*n(idx_O)

    !d[H2_dot]/d[SI+]
    pd(6,80) =  &
        -k(1203)*n(idx_H2)

    !d[H_dot]/d[SI+]
    pd(8,80) =  &
        +k(573)*n(idx_H2O)  &
        +k(860)*n(idx_OH)  &
        +k(477)*n(idx_CH)  &
        -k(1210)*n(idx_H)

    !d[H2O_dot]/d[SI+]
    pd(9,80) =  &
        -k(573)*n(idx_H2O)

    !d[OH_dot]/d[SI+]
    pd(10,80) =  &
        -k(860)*n(idx_OH)

    !d[MG_dot]/d[SI+]
    pd(15,80) =  &
        -k(154)*n(idx_MG)

    !d[SI_dot]/d[SI+]
    pd(18,80) =  &
        +k(154)*n(idx_MG)  &
        +k(1223)*n(idx_E)

    !d[CH3_dot]/d[SI+]
    pd(28,80) =  &
        +k(861)*n(idx_CH3OH)

    !d[CH3OH_dot]/d[SI+]
    pd(37,80) =  &
        -k(861)*n(idx_CH3OH)

    !d[SIH4_DUST_dot]/d[SI+]
    pd(48,80) =  &
        +k(1232)

    !d[MG+_dot]/d[SI+]
    pd(77,80) =  &
        +k(154)*n(idx_MG)

    !d[SI+_dot]/d[SI+]
    pd(80,80) =  &
        -k(154)*n(idx_MG)  &
        -k(861)*n(idx_CH3OH)  &
        -k(477)*n(idx_CH)  &
        -k(1232)  &
        -k(860)*n(idx_OH)  &
        -k(1223)*n(idx_E)  &
        -k(1203)*n(idx_H2)  &
        -k(1213)*n(idx_O)  &
        -k(573)*n(idx_H2O)  &
        -k(1210)*n(idx_H)

    !d[SIC+_dot]/d[SI+]
    pd(83,80) =  &
        +k(477)*n(idx_CH)

    !d[SIH2+_dot]/d[SI+]
    pd(84,80) =  &
        +k(1203)*n(idx_H2)

    !d[SIH+_dot]/d[SI+]
    pd(100,80) =  &
        +k(1210)*n(idx_H)

    !d[SIO+_dot]/d[SI+]
    pd(101,80) =  &
        +k(860)*n(idx_OH)  &
        +k(1213)*n(idx_O)

    !d[SIOH+_dot]/d[SI+]
    pd(115,80) =  &
        +k(861)*n(idx_CH3OH)  &
        +k(573)*n(idx_H2O)

    !d[E_dot]/d[SIC2+]
    pd(1,81) =  &
        -k(351)*n(idx_E)

    !d[C_dot]/d[SIC2+]
    pd(7,81) =  &
        +k(351)*n(idx_E)

    !d[SIC_dot]/d[SIC2+]
    pd(21,81) =  &
        +k(351)*n(idx_E)

    !d[SIC2_DUST_dot]/d[SIC2+]
    pd(51,81) =  &
        +k(1243)

    !d[SIC2+_dot]/d[SIC2+]
    pd(81,81) =  &
        -k(351)*n(idx_E)  &
        -k(1243)

    !d[E_dot]/d[SIC3+]
    pd(1,82) =  &
        -k(352)*n(idx_E)

    !d[C_dot]/d[SIC3+]
    pd(7,82) =  &
        +k(352)*n(idx_E)

    !d[SIC2_dot]/d[SIC3+]
    pd(19,82) =  &
        +k(352)*n(idx_E)

    !d[SIC3_DUST_dot]/d[SIC3+]
    pd(52,82) =  &
        +k(1245)

    !d[SIC3+_dot]/d[SIC3+]
    pd(82,82) =  &
        -k(1245)  &
        -k(352)*n(idx_E)

    !d[E_dot]/d[SIC+]
    pd(1,83) =  &
        -k(350)*n(idx_E)

    !d[O_dot]/d[SIC+]
    pd(3,83) =  &
        -k(832)*n(idx_O)

    !d[C_dot]/d[SIC+]
    pd(7,83) =  &
        +k(832)*n(idx_O)  &
        +k(350)*n(idx_E)

    !d[SI_dot]/d[SIC+]
    pd(18,83) =  &
        +k(350)*n(idx_E)

    !d[CN_dot]/d[SIC+]
    pd(24,83) =  &
        +k(745)*n(idx_N)

    !d[N_dot]/d[SIC+]
    pd(30,83) =  &
        -k(745)*n(idx_N)

    !d[SIC_DUST_dot]/d[SIC+]
    pd(50,83) =  &
        +k(1242)

    !d[SI+_dot]/d[SIC+]
    pd(80,83) =  &
        +k(745)*n(idx_N)

    !d[SIC+_dot]/d[SIC+]
    pd(83,83) =  &
        -k(832)*n(idx_O)  &
        -k(350)*n(idx_E)  &
        -k(745)*n(idx_N)  &
        -k(1242)

    !d[SIO+_dot]/d[SIC+]
    pd(101,83) =  &
        +k(832)*n(idx_O)

    !d[E_dot]/d[SIH2+]
    pd(1,84) =  &
        -k(356)*n(idx_E)  &
        -k(355)*n(idx_E)  &
        -k(354)*n(idx_E)

    !d[O_dot]/d[SIH2+]
    pd(3,84) =  &
        -k(834)*n(idx_O)

    !d[H2_dot]/d[SIH2+]
    pd(6,84) =  &
        +k(354)*n(idx_E)

    !d[H_dot]/d[SIH2+]
    pd(8,84) =  &
        +2.d0*k(355)*n(idx_E)  &
        +k(834)*n(idx_O)  &
        +k(356)*n(idx_E)

    !d[OH_dot]/d[SIH2+]
    pd(10,84) =  &
        +k(863)*n(idx_O2)

    !d[O2_dot]/d[SIH2+]
    pd(11,84) =  &
        -k(863)*n(idx_O2)

    !d[SI_dot]/d[SIH2+]
    pd(18,84) =  &
        +k(355)*n(idx_E)  &
        +k(354)*n(idx_E)

    !d[SIH_dot]/d[SIH2+]
    pd(33,84) =  &
        +k(356)*n(idx_E)

    !d[SIH4_DUST_dot]/d[SIH2+]
    pd(48,84) =  &
        +k(1235)

    !d[SIH2+_dot]/d[SIH2+]
    pd(84,84) =  &
        -k(355)*n(idx_E)  &
        -k(354)*n(idx_E)  &
        -k(356)*n(idx_E)  &
        -k(863)*n(idx_O2)  &
        -k(834)*n(idx_O)  &
        -k(1235)

    !d[SIOH+_dot]/d[SIH2+]
    pd(115,84) =  &
        +k(834)*n(idx_O)  &
        +k(863)*n(idx_O2)

    !d[E_dot]/d[SIH3+]
    pd(1,85) =  &
        -k(358)*n(idx_E)  &
        -k(357)*n(idx_E)

    !d[O_dot]/d[SIH3+]
    pd(3,85) =  &
        -k(835)*n(idx_O)

    !d[H2_dot]/d[SIH3+]
    pd(6,85) =  &
        +k(835)*n(idx_O)  &
        -k(1205)*n(idx_H2)  &
        +k(358)*n(idx_E)

    !d[H_dot]/d[SIH3+]
    pd(8,85) =  &
        +k(357)*n(idx_E)

    !d[SIH2_dot]/d[SIH3+]
    pd(22,85) =  &
        +k(357)*n(idx_E)

    !d[SIH_dot]/d[SIH3+]
    pd(33,85) =  &
        +k(358)*n(idx_E)

    !d[SIH4_DUST_dot]/d[SIH3+]
    pd(48,85) =  &
        +k(1237)

    !d[SIH3+_dot]/d[SIH3+]
    pd(85,85) =  &
        -k(1237)  &
        -k(358)*n(idx_E)  &
        -k(357)*n(idx_E)  &
        -k(1205)*n(idx_H2)  &
        -k(835)*n(idx_O)

    !d[SIH5+_dot]/d[SIH3+]
    pd(114,85) =  &
        +k(1205)*n(idx_H2)

    !d[SIOH+_dot]/d[SIH3+]
    pd(115,85) =  &
        +k(835)*n(idx_O)

    !d[E_dot]/d[CN+]
    pd(1,86) =  &
        -k(304)*n(idx_E)

    !d[CH_dot]/d[CN+]
    pd(2,86) =  &
        -k(54)*n(idx_CH)

    !d[O_dot]/d[CN+]
    pd(3,86) =  &
        -k(217)*n(idx_O)

    !d[HCN_dot]/d[CN+]
    pd(5,86) =  &
        +k(480)*n(idx_H2CO)  &
        -k(66)*n(idx_HCN)

    !d[H2_dot]/d[CN+]
    pd(6,86) =  &
        -k(532)*n(idx_H2)

    !d[C_dot]/d[CN+]
    pd(7,86) =  &
        +k(738)*n(idx_N)  &
        -k(28)*n(idx_C)  &
        +k(304)*n(idx_E)

    !d[H_dot]/d[CN+]
    pd(8,86) =  &
        -k(127)*n(idx_H)  &
        +k(532)*n(idx_H2)

    !d[H2O_dot]/d[CN+]
    pd(9,86) =  &
        -k(561)*n(idx_H2O)  &
        -k(562)*n(idx_H2O)

    !d[OH_dot]/d[CN+]
    pd(10,86) =  &
        +k(561)*n(idx_H2O)  &
        -k(226)*n(idx_OH)

    !d[O2_dot]/d[CN+]
    pd(11,86) =  &
        -k(69)*n(idx_O2)  &
        -k(482)*n(idx_O2)

    !d[CH2_dot]/d[CN+]
    pd(12,86) =  &
        -k(38)*n(idx_CH2)

    !d[H2CO_dot]/d[CN+]
    pd(13,86) =  &
        -k(65)*n(idx_H2CO)  &
        -k(480)*n(idx_H2CO)

    !d[HCO_dot]/d[CN+]
    pd(14,86) =  &
        -k(67)*n(idx_HCO)  &
        -k(481)*n(idx_HCO)

    !d[NO_dot]/d[CN+]
    pd(17,86) =  &
        -k(68)*n(idx_NO)

    !d[CN_dot]/d[CN+]
    pd(24,86) =  &
        +k(226)*n(idx_OH)  &
        +k(184)*n(idx_NH2)  &
        +k(28)*n(idx_C)  &
        +k(127)*n(idx_H)  &
        +k(200)*n(idx_NH)  &
        +k(54)*n(idx_CH)  &
        +k(68)*n(idx_NO)  &
        +k(67)*n(idx_HCO)  &
        +k(217)*n(idx_O)  &
        +k(69)*n(idx_O2)  &
        +k(65)*n(idx_H2CO)  &
        +k(66)*n(idx_HCN)  &
        +k(38)*n(idx_CH2)  &
        +k(64)*n(idx_CO)

    !d[CO_dot]/d[CN+]
    pd(25,86) =  &
        -k(64)*n(idx_CO)  &
        +k(482)*n(idx_O2)  &
        +k(481)*n(idx_HCO)

    !d[NH2_dot]/d[CN+]
    pd(27,86) =  &
        -k(184)*n(idx_NH2)

    !d[N_dot]/d[CN+]
    pd(30,86) =  &
        +k(304)*n(idx_E)  &
        -k(738)*n(idx_N)

    !d[NH_dot]/d[CN+]
    pd(31,86) =  &
        +k(562)*n(idx_H2O)  &
        -k(200)*n(idx_NH)

    !d[HCN_DUST_dot]/d[CN+]
    pd(59,86) =  &
        +k(1277)

    !d[HCO+_dot]/d[CN+]
    pd(70,86) =  &
        +k(562)*n(idx_H2O)  &
        +k(480)*n(idx_H2CO)  &
        +k(67)*n(idx_HCO)

    !d[H+_dot]/d[CN+]
    pd(71,86) =  &
        +k(127)*n(idx_H)

    !d[C+_dot]/d[CN+]
    pd(73,86) =  &
        +k(28)*n(idx_C)

    !d[CH2+_dot]/d[CN+]
    pd(74,86) =  &
        +k(38)*n(idx_CH2)

    !d[CH+_dot]/d[CN+]
    pd(75,86) =  &
        +k(54)*n(idx_CH)

    !d[H2CO+_dot]/d[CN+]
    pd(76,86) =  &
        +k(65)*n(idx_H2CO)

    !d[NO+_dot]/d[CN+]
    pd(79,86) =  &
        +k(482)*n(idx_O2)  &
        +k(68)*n(idx_NO)

    !d[CN+_dot]/d[CN+]
    pd(86,86) =  &
        -k(304)*n(idx_E)  &
        -k(67)*n(idx_HCO)  &
        -k(482)*n(idx_O2)  &
        -k(68)*n(idx_NO)  &
        -k(127)*n(idx_H)  &
        -k(65)*n(idx_H2CO)  &
        -k(38)*n(idx_CH2)  &
        -k(481)*n(idx_HCO)  &
        -k(66)*n(idx_HCN)  &
        -k(562)*n(idx_H2O)  &
        -k(69)*n(idx_O2)  &
        -k(561)*n(idx_H2O)  &
        -k(217)*n(idx_O)  &
        -k(64)*n(idx_CO)  &
        -k(738)*n(idx_N)  &
        -k(28)*n(idx_C)  &
        -k(184)*n(idx_NH2)  &
        -k(200)*n(idx_NH)  &
        -k(54)*n(idx_CH)  &
        -k(1277)  &
        -k(480)*n(idx_H2CO)  &
        -k(226)*n(idx_OH)  &
        -k(532)*n(idx_H2)

    !d[CO+_dot]/d[CN+]
    pd(87,86) =  &
        +k(64)*n(idx_CO)

    !d[N2+_dot]/d[CN+]
    pd(88,86) =  &
        +k(738)*n(idx_N)

    !d[O2+_dot]/d[CN+]
    pd(89,86) =  &
        +k(69)*n(idx_O2)

    !d[NH2+_dot]/d[CN+]
    pd(91,86) =  &
        +k(184)*n(idx_NH2)

    !d[O+_dot]/d[CN+]
    pd(92,86) =  &
        +k(217)*n(idx_O)

    !d[OH+_dot]/d[CN+]
    pd(93,86) =  &
        +k(226)*n(idx_OH)

    !d[HCN+_dot]/d[CN+]
    pd(97,86) =  &
        +k(561)*n(idx_H2O)  &
        +k(66)*n(idx_HCN)  &
        +k(481)*n(idx_HCO)  &
        +k(532)*n(idx_H2)

    !d[NH+_dot]/d[CN+]
    pd(98,86) =  &
        +k(200)*n(idx_NH)

    !d[E_dot]/d[CO+]
    pd(1,87) =  &
        -k(305)*n(idx_E)

    !d[CH_dot]/d[CO+]
    pd(2,87) =  &
        +k(423)*n(idx_CH2)  &
        -k(55)*n(idx_CH)  &
        -k(459)*n(idx_CH)

    !d[O_dot]/d[CO+]
    pd(3,87) =  &
        +k(1132)  &
        -k(218)*n(idx_O)  &
        +k(305)*n(idx_E)  &
        +k(852)*n(idx_OH)

    !d[HCN_dot]/d[CO+]
    pd(5,87) =  &
        -k(135)*n(idx_HCN)

    !d[H2_dot]/d[CO+]
    pd(6,87) =  &
        -k(534)*n(idx_H2)  &
        -k(533)*n(idx_H2)

    !d[C_dot]/d[CO+]
    pd(7,87) =  &
        +k(459)*n(idx_CH)  &
        +k(305)*n(idx_E)  &
        -k(29)*n(idx_C)

    !d[H_dot]/d[CO+]
    pd(8,87) =  &
        +k(533)*n(idx_H2)  &
        +k(534)*n(idx_H2)  &
        -k(128)*n(idx_H)

    !d[H2O_dot]/d[CO+]
    pd(9,87) =  &
        -k(124)*n(idx_H2O)  &
        -k(563)*n(idx_H2O)

    !d[OH_dot]/d[CO+]
    pd(10,87) =  &
        -k(852)*n(idx_OH)  &
        -k(227)*n(idx_OH)  &
        +k(563)*n(idx_H2O)

    !d[O2_dot]/d[CO+]
    pd(11,87) =  &
        -k(74)*n(idx_O2)

    !d[CH2_dot]/d[CO+]
    pd(12,87) =  &
        -k(39)*n(idx_CH2)  &
        -k(423)*n(idx_CH2)

    !d[H2CO_dot]/d[CO+]
    pd(13,87) =  &
        -k(71)*n(idx_H2CO)  &
        -k(485)*n(idx_H2CO)

    !d[HCO_dot]/d[CO+]
    pd(14,87) =  &
        -k(72)*n(idx_HCO)  &
        +k(485)*n(idx_H2CO)

    !d[NH3_dot]/d[CO+]
    pd(16,87) =  &
        -k(194)*n(idx_NH3)  &
        -k(793)*n(idx_NH3)

    !d[NO_dot]/d[CO+]
    pd(17,87) =  &
        -k(73)*n(idx_NO)

    !d[CO_dot]/d[CO+]
    pd(25,87) =  &
        +k(72)*n(idx_HCO)  &
        +k(124)*n(idx_H2O)  &
        +k(185)*n(idx_NH2)  &
        +k(73)*n(idx_NO)  &
        +k(201)*n(idx_NH)  &
        +k(128)*n(idx_H)  &
        +k(53)*n(idx_CH4)  &
        +k(227)*n(idx_OH)  &
        +k(55)*n(idx_CH)  &
        +k(71)*n(idx_H2CO)  &
        +k(29)*n(idx_C)  &
        +k(39)*n(idx_CH2)  &
        +k(74)*n(idx_O2)  &
        +k(135)*n(idx_HCN)  &
        +k(218)*n(idx_O)  &
        +k(194)*n(idx_NH3)

    !d[NH2_dot]/d[CO+]
    pd(27,87) =  &
        -k(185)*n(idx_NH2)  &
        +k(793)*n(idx_NH3)  &
        -k(780)*n(idx_NH2)

    !d[CH3_dot]/d[CO+]
    pd(28,87) =  &
        +k(452)*n(idx_CH4)

    !d[CH4_dot]/d[CO+]
    pd(29,87) =  &
        -k(452)*n(idx_CH4)  &
        -k(53)*n(idx_CH4)

    !d[N_dot]/d[CO+]
    pd(30,87) =  &
        +k(796)*n(idx_NH)

    !d[NH_dot]/d[CO+]
    pd(31,87) =  &
        +k(780)*n(idx_NH2)  &
        -k(201)*n(idx_NH)  &
        -k(796)*n(idx_NH)

    !d[CO_DUST_dot]/d[CO+]
    pd(54,87) =  &
        +k(1276)

    !d[HCO+_dot]/d[CO+]
    pd(70,87) =  &
        +k(423)*n(idx_CH2)  &
        +k(533)*n(idx_H2)  &
        +k(793)*n(idx_NH3)  &
        +k(796)*n(idx_NH)  &
        +k(563)*n(idx_H2O)  &
        +k(72)*n(idx_HCO)  &
        +k(780)*n(idx_NH2)  &
        +k(485)*n(idx_H2CO)  &
        +k(452)*n(idx_CH4)  &
        +k(852)*n(idx_OH)  &
        +k(459)*n(idx_CH)

    !d[H+_dot]/d[CO+]
    pd(71,87) =  &
        +k(128)*n(idx_H)

    !d[HOC+_dot]/d[CO+]
    pd(72,87) =  &
        +k(534)*n(idx_H2)

    !d[C+_dot]/d[CO+]
    pd(73,87) =  &
        +k(1132)  &
        +k(29)*n(idx_C)

    !d[CH2+_dot]/d[CO+]
    pd(74,87) =  &
        +k(39)*n(idx_CH2)

    !d[CH+_dot]/d[CO+]
    pd(75,87) =  &
        +k(55)*n(idx_CH)

    !d[H2CO+_dot]/d[CO+]
    pd(76,87) =  &
        +k(71)*n(idx_H2CO)

    !d[NH3+_dot]/d[CO+]
    pd(78,87) =  &
        +k(194)*n(idx_NH3)

    !d[NO+_dot]/d[CO+]
    pd(79,87) =  &
        +k(73)*n(idx_NO)

    !d[CO+_dot]/d[CO+]
    pd(87,87) =  &
        -k(55)*n(idx_CH)  &
        -k(74)*n(idx_O2)  &
        -k(53)*n(idx_CH4)  &
        -k(305)*n(idx_E)  &
        -k(201)*n(idx_NH)  &
        -k(485)*n(idx_H2CO)  &
        -k(73)*n(idx_NO)  &
        -k(71)*n(idx_H2CO)  &
        -k(128)*n(idx_H)  &
        -k(534)*n(idx_H2)  &
        -k(1132)  &
        -k(29)*n(idx_C)  &
        -k(459)*n(idx_CH)  &
        -k(227)*n(idx_OH)  &
        -k(793)*n(idx_NH3)  &
        -k(780)*n(idx_NH2)  &
        -k(124)*n(idx_H2O)  &
        -k(194)*n(idx_NH3)  &
        -k(135)*n(idx_HCN)  &
        -k(796)*n(idx_NH)  &
        -k(185)*n(idx_NH2)  &
        -k(72)*n(idx_HCO)  &
        -k(1276)  &
        -k(852)*n(idx_OH)  &
        -k(39)*n(idx_CH2)  &
        -k(563)*n(idx_H2O)  &
        -k(218)*n(idx_O)  &
        -k(452)*n(idx_CH4)  &
        -k(423)*n(idx_CH2)  &
        -k(533)*n(idx_H2)

    !d[O2+_dot]/d[CO+]
    pd(89,87) =  &
        +k(74)*n(idx_O2)

    !d[H2O+_dot]/d[CO+]
    pd(90,87) =  &
        +k(124)*n(idx_H2O)

    !d[NH2+_dot]/d[CO+]
    pd(91,87) =  &
        +k(185)*n(idx_NH2)

    !d[O+_dot]/d[CO+]
    pd(92,87) =  &
        +k(218)*n(idx_O)

    !d[OH+_dot]/d[CO+]
    pd(93,87) =  &
        +k(227)*n(idx_OH)

    !d[CH4+_dot]/d[CO+]
    pd(95,87) =  &
        +k(53)*n(idx_CH4)

    !d[HCN+_dot]/d[CO+]
    pd(97,87) =  &
        +k(135)*n(idx_HCN)

    !d[NH+_dot]/d[CO+]
    pd(98,87) =  &
        +k(201)*n(idx_NH)

    !d[E_dot]/d[N2+]
    pd(1,88) =  &
        -k(338)*n(idx_E)

    !d[CH_dot]/d[N2+]
    pd(2,88) =  &
        -k(59)*n(idx_CH)

    !d[O_dot]/d[N2+]
    pd(3,88) =  &
        -k(826)*n(idx_O)  &
        -k(219)*n(idx_O)

    !d[HCN_dot]/d[N2+]
    pd(5,88) =  &
        -k(136)*n(idx_HCN)

    !d[H2_dot]/d[N2+]
    pd(6,88) =  &
        -k(540)*n(idx_H2)  &
        +k(456)*n(idx_CH4)

    !d[C_dot]/d[N2+]
    pd(7,88) =  &
        -k(30)*n(idx_C)

    !d[H_dot]/d[N2+]
    pd(8,88) =  &
        +k(457)*n(idx_CH4)  &
        +k(731)*n(idx_H2CO)  &
        +k(540)*n(idx_H2)

    !d[H2O_dot]/d[N2+]
    pd(9,88) =  &
        -k(126)*n(idx_H2O)  &
        -k(570)*n(idx_H2O)

    !d[OH_dot]/d[N2+]
    pd(10,88) =  &
        +k(570)*n(idx_H2O)  &
        -k(228)*n(idx_OH)

    !d[O2_dot]/d[N2+]
    pd(11,88) =  &
        -k(174)*n(idx_O2)

    !d[CH2_dot]/d[N2+]
    pd(12,88) =  &
        -k(42)*n(idx_CH2)

    !d[H2CO_dot]/d[N2+]
    pd(13,88) =  &
        -k(171)*n(idx_H2CO)  &
        -k(731)*n(idx_H2CO)

    !d[HCO_dot]/d[N2+]
    pd(14,88) =  &
        -k(172)*n(idx_HCO)  &
        -k(732)*n(idx_HCO)

    !d[MG_dot]/d[N2+]
    pd(15,88) =  &
        -k(151)*n(idx_MG)

    !d[NH3_dot]/d[N2+]
    pd(16,88) =  &
        -k(198)*n(idx_NH3)

    !d[NO_dot]/d[N2+]
    pd(17,88) =  &
        -k(173)*n(idx_NO)

    !d[CN_dot]/d[N2+]
    pd(24,88) =  &
        -k(70)*n(idx_CN)

    !d[CO_dot]/d[N2+]
    pd(25,88) =  &
        +k(732)*n(idx_HCO)  &
        -k(75)*n(idx_CO)

    !d[N2_dot]/d[N2+]
    pd(26,88) =  &
        +k(75)*n(idx_CO)  &
        +k(219)*n(idx_O)  &
        +k(174)*n(idx_O2)  &
        +k(228)*n(idx_OH)  &
        +k(456)*n(idx_CH4)  &
        +k(30)*n(idx_C)  &
        +k(187)*n(idx_NH2)  &
        +k(175)*n(idx_N)  &
        +k(172)*n(idx_HCO)  &
        +k(136)*n(idx_HCN)  &
        +k(151)*n(idx_MG)  &
        +k(171)*n(idx_H2CO)  &
        +k(198)*n(idx_NH3)  &
        +k(42)*n(idx_CH2)  &
        +k(173)*n(idx_NO)  &
        +k(202)*n(idx_NH)  &
        +k(70)*n(idx_CN)  &
        +k(126)*n(idx_H2O)  &
        +k(731)*n(idx_H2CO)  &
        +k(59)*n(idx_CH)  &
        +k(457)*n(idx_CH4)

    !d[NH2_dot]/d[N2+]
    pd(27,88) =  &
        -k(187)*n(idx_NH2)

    !d[CH4_dot]/d[N2+]
    pd(29,88) =  &
        -k(456)*n(idx_CH4)  &
        -k(457)*n(idx_CH4)

    !d[N_dot]/d[N2+]
    pd(30,88) =  &
        +2.d0*k(338)*n(idx_E)  &
        -k(175)*n(idx_N)  &
        +k(826)*n(idx_O)

    !d[NH_dot]/d[N2+]
    pd(31,88) =  &
        -k(202)*n(idx_NH)

    !d[N2_DUST_dot]/d[N2+]
    pd(58,88) =  &
        +k(1272)

    !d[HCO+_dot]/d[N2+]
    pd(70,88) =  &
        +k(172)*n(idx_HCO)  &
        +k(731)*n(idx_H2CO)

    !d[C+_dot]/d[N2+]
    pd(73,88) =  &
        +k(30)*n(idx_C)

    !d[CH2+_dot]/d[N2+]
    pd(74,88) =  &
        +k(42)*n(idx_CH2)  &
        +k(456)*n(idx_CH4)

    !d[CH+_dot]/d[N2+]
    pd(75,88) =  &
        +k(59)*n(idx_CH)

    !d[H2CO+_dot]/d[N2+]
    pd(76,88) =  &
        +k(171)*n(idx_H2CO)

    !d[MG+_dot]/d[N2+]
    pd(77,88) =  &
        +k(151)*n(idx_MG)

    !d[NH3+_dot]/d[N2+]
    pd(78,88) =  &
        +k(198)*n(idx_NH3)

    !d[NO+_dot]/d[N2+]
    pd(79,88) =  &
        +k(826)*n(idx_O)  &
        +k(173)*n(idx_NO)

    !d[CN+_dot]/d[N2+]
    pd(86,88) =  &
        +k(70)*n(idx_CN)

    !d[CO+_dot]/d[N2+]
    pd(87,88) =  &
        +k(75)*n(idx_CO)

    !d[N2+_dot]/d[N2+]
    pd(88,88) =  &
        -k(187)*n(idx_NH2)  &
        -k(172)*n(idx_HCO)  &
        -k(59)*n(idx_CH)  &
        -k(173)*n(idx_NO)  &
        -k(198)*n(idx_NH3)  &
        -k(75)*n(idx_CO)  &
        -k(30)*n(idx_C)  &
        -k(338)*n(idx_E)  &
        -k(202)*n(idx_NH)  &
        -k(570)*n(idx_H2O)  &
        -k(70)*n(idx_CN)  &
        -k(457)*n(idx_CH4)  &
        -k(219)*n(idx_O)  &
        -k(151)*n(idx_MG)  &
        -k(175)*n(idx_N)  &
        -k(1272)  &
        -k(228)*n(idx_OH)  &
        -k(126)*n(idx_H2O)  &
        -k(540)*n(idx_H2)  &
        -k(731)*n(idx_H2CO)  &
        -k(826)*n(idx_O)  &
        -k(732)*n(idx_HCO)  &
        -k(174)*n(idx_O2)  &
        -k(136)*n(idx_HCN)  &
        -k(456)*n(idx_CH4)  &
        -k(171)*n(idx_H2CO)  &
        -k(42)*n(idx_CH2)

    !d[O2+_dot]/d[N2+]
    pd(89,88) =  &
        +k(174)*n(idx_O2)

    !d[H2O+_dot]/d[N2+]
    pd(90,88) =  &
        +k(126)*n(idx_H2O)

    !d[NH2+_dot]/d[N2+]
    pd(91,88) =  &
        +k(187)*n(idx_NH2)

    !d[O+_dot]/d[N2+]
    pd(92,88) =  &
        +k(219)*n(idx_O)

    !d[OH+_dot]/d[N2+]
    pd(93,88) =  &
        +k(228)*n(idx_OH)

    !d[CH3+_dot]/d[N2+]
    pd(94,88) =  &
        +k(457)*n(idx_CH4)

    !d[N+_dot]/d[N2+]
    pd(96,88) =  &
        +k(175)*n(idx_N)

    !d[HCN+_dot]/d[N2+]
    pd(97,88) =  &
        +k(136)*n(idx_HCN)

    !d[NH+_dot]/d[N2+]
    pd(98,88) =  &
        +k(202)*n(idx_NH)

    !d[N2H+_dot]/d[N2+]
    pd(112,88) =  &
        +k(732)*n(idx_HCO)  &
        +k(570)*n(idx_H2O)  &
        +k(540)*n(idx_H2)

    !d[E_dot]/d[O2+]
    pd(1,89) =  &
        -k(347)*n(idx_E)

    !d[CH_dot]/d[O2+]
    pd(2,89) =  &
        -k(62)*n(idx_CH)  &
        -k(474)*n(idx_CH)

    !d[O_dot]/d[O2+]
    pd(3,89) =  &
        +2.d0*k(347)*n(idx_E)  &
        +k(436)*n(idx_CH2)  &
        +k(474)*n(idx_CH)  &
        +k(1168)  &
        +k(392)*n(idx_C)  &
        +k(743)*n(idx_N)  &
        +k(805)*n(idx_NH)

    !d[C_dot]/d[O2+]
    pd(7,89) =  &
        -k(31)*n(idx_C)  &
        -k(392)*n(idx_C)

    !d[H_dot]/d[O2+]
    pd(8,89) =  &
        +k(821)*n(idx_CH3OH)  &
        +k(552)*n(idx_H2CO)

    !d[O2_dot]/d[O2+]
    pd(11,89) =  &
        +k(138)*n(idx_HCO)  &
        +k(31)*n(idx_C)  &
        +k(153)*n(idx_MG)  &
        +k(188)*n(idx_NH2)  &
        +k(62)*n(idx_CH)  &
        +k(117)*n(idx_H2CO)  &
        +k(206)*n(idx_NO)  &
        +k(231)*n(idx_SI)  &
        +k(199)*n(idx_NH3)  &
        +k(821)*n(idx_CH3OH)  &
        +k(552)*n(idx_H2CO)  &
        +k(45)*n(idx_CH2)

    !d[CH2_dot]/d[O2+]
    pd(12,89) =  &
        -k(45)*n(idx_CH2)  &
        -k(436)*n(idx_CH2)

    !d[H2CO_dot]/d[O2+]
    pd(13,89) =  &
        -k(552)*n(idx_H2CO)  &
        -k(117)*n(idx_H2CO)

    !d[HCO_dot]/d[O2+]
    pd(14,89) =  &
        -k(645)*n(idx_HCO)  &
        -k(138)*n(idx_HCO)

    !d[MG_dot]/d[O2+]
    pd(15,89) =  &
        -k(153)*n(idx_MG)

    !d[NH3_dot]/d[O2+]
    pd(16,89) =  &
        -k(199)*n(idx_NH3)

    !d[NO_dot]/d[O2+]
    pd(17,89) =  &
        -k(206)*n(idx_NO)

    !d[SI_dot]/d[O2+]
    pd(18,89) =  &
        -k(231)*n(idx_SI)

    !d[CO_dot]/d[O2+]
    pd(25,89) =  &
        +k(645)*n(idx_HCO)

    !d[NH2_dot]/d[O2+]
    pd(27,89) =  &
        -k(188)*n(idx_NH2)

    !d[N_dot]/d[O2+]
    pd(30,89) =  &
        -k(743)*n(idx_N)

    !d[NH_dot]/d[O2+]
    pd(31,89) =  &
        -k(805)*n(idx_NH)

    !d[CH3OH_dot]/d[O2+]
    pd(37,89) =  &
        -k(821)*n(idx_CH3OH)

    !d[O2_DUST_dot]/d[O2+]
    pd(61,89) =  &
        +k(1271)

    !d[HCO+_dot]/d[O2+]
    pd(70,89) =  &
        +k(138)*n(idx_HCO)  &
        +k(474)*n(idx_CH)  &
        +k(552)*n(idx_H2CO)

    !d[C+_dot]/d[O2+]
    pd(73,89) =  &
        +k(31)*n(idx_C)

    !d[CH2+_dot]/d[O2+]
    pd(74,89) =  &
        +k(45)*n(idx_CH2)

    !d[CH+_dot]/d[O2+]
    pd(75,89) =  &
        +k(62)*n(idx_CH)

    !d[H2CO+_dot]/d[O2+]
    pd(76,89) =  &
        +k(436)*n(idx_CH2)  &
        +k(117)*n(idx_H2CO)

    !d[MG+_dot]/d[O2+]
    pd(77,89) =  &
        +k(153)*n(idx_MG)

    !d[NH3+_dot]/d[O2+]
    pd(78,89) =  &
        +k(199)*n(idx_NH3)

    !d[NO+_dot]/d[O2+]
    pd(79,89) =  &
        +k(206)*n(idx_NO)  &
        +k(743)*n(idx_N)

    !d[SI+_dot]/d[O2+]
    pd(80,89) =  &
        +k(231)*n(idx_SI)

    !d[CO+_dot]/d[O2+]
    pd(87,89) =  &
        +k(392)*n(idx_C)

    !d[O2+_dot]/d[O2+]
    pd(89,89) =  &
        -k(31)*n(idx_C)  &
        -k(805)*n(idx_NH)  &
        -k(1271)  &
        -k(62)*n(idx_CH)  &
        -k(117)*n(idx_H2CO)  &
        -k(347)*n(idx_E)  &
        -k(153)*n(idx_MG)  &
        -k(436)*n(idx_CH2)  &
        -k(45)*n(idx_CH2)  &
        -k(206)*n(idx_NO)  &
        -k(645)*n(idx_HCO)  &
        -k(743)*n(idx_N)  &
        -k(231)*n(idx_SI)  &
        -k(552)*n(idx_H2CO)  &
        -k(199)*n(idx_NH3)  &
        -k(138)*n(idx_HCO)  &
        -k(392)*n(idx_C)  &
        -k(188)*n(idx_NH2)  &
        -k(821)*n(idx_CH3OH)  &
        -k(474)*n(idx_CH)  &
        -k(1168)

    !d[NH2+_dot]/d[O2+]
    pd(91,89) =  &
        +k(188)*n(idx_NH2)

    !d[O+_dot]/d[O2+]
    pd(92,89) =  &
        +k(1168)

    !d[HNO+_dot]/d[O2+]
    pd(104,89) =  &
        +k(805)*n(idx_NH)

    !d[H3CO+_dot]/d[O2+]
    pd(107,89) =  &
        +k(821)*n(idx_CH3OH)

    !d[O2H+_dot]/d[O2+]
    pd(113,89) =  &
        +k(645)*n(idx_HCO)

    !d[E_dot]/d[H2O+]
    pd(1,90) =  &
        -k(315)*n(idx_E)  &
        -k(314)*n(idx_E)  &
        -k(313)*n(idx_E)

    !d[CH_dot]/d[H2O+]
    pd(2,90) =  &
        -k(57)*n(idx_CH)  &
        -k(461)*n(idx_CH)

    !d[O_dot]/d[H2O+]
    pd(3,90) =  &
        +k(853)*n(idx_OH)  &
        +k(313)*n(idx_E)  &
        -k(824)*n(idx_O)  &
        +k(314)*n(idx_E)

    !d[HNC_dot]/d[H2O+]
    pd(4,90) =  &
        -k(560)*n(idx_HNC)

    !d[HCN_dot]/d[H2O+]
    pd(5,90) =  &
        -k(557)*n(idx_HCN)

    !d[H2_dot]/d[H2O+]
    pd(6,90) =  &
        -k(535)*n(idx_H2)  &
        +k(824)*n(idx_O)  &
        +k(313)*n(idx_E)  &
        +k(740)*n(idx_N)

    !d[C_dot]/d[H2O+]
    pd(7,90) =  &
        -k(384)*n(idx_C)

    !d[H_dot]/d[H2O+]
    pd(8,90) =  &
        +k(1141)  &
        +k(535)*n(idx_H2)  &
        +k(739)*n(idx_N)  &
        +k(315)*n(idx_E)  &
        +2.d0*k(314)*n(idx_E)

    !d[H2O_dot]/d[H2O+]
    pd(9,90) =  &
        +k(186)*n(idx_NH2)  &
        +k(118)*n(idx_H2CO)  &
        +k(120)*n(idx_MG)  &
        -k(556)*n(idx_H2O)  &
        +k(57)*n(idx_CH)  &
        +k(196)*n(idx_NH3)  &
        +k(41)*n(idx_CH2)  &
        +k(119)*n(idx_HCO)  &
        +k(121)*n(idx_NO)  &
        +k(122)*n(idx_O2)  &
        +k(123)*n(idx_SI)

    !d[OH_dot]/d[H2O+]
    pd(10,90) =  &
        +k(384)*n(idx_C)  &
        +k(782)*n(idx_NH2)  &
        +k(425)*n(idx_CH2)  &
        -k(853)*n(idx_OH)  &
        +k(315)*n(idx_E)  &
        +k(557)*n(idx_HCN)  &
        +k(554)*n(idx_CO)  &
        +k(560)*n(idx_HNC)  &
        +k(556)*n(idx_H2O)  &
        +k(461)*n(idx_CH)  &
        +k(555)*n(idx_H2CO)  &
        +k(559)*n(idx_HCO)

    !d[O2_dot]/d[H2O+]
    pd(11,90) =  &
        -k(122)*n(idx_O2)

    !d[CH2_dot]/d[H2O+]
    pd(12,90) =  &
        -k(425)*n(idx_CH2)  &
        -k(41)*n(idx_CH2)

    !d[H2CO_dot]/d[H2O+]
    pd(13,90) =  &
        -k(118)*n(idx_H2CO)  &
        -k(555)*n(idx_H2CO)

    !d[HCO_dot]/d[H2O+]
    pd(14,90) =  &
        -k(119)*n(idx_HCO)  &
        -k(558)*n(idx_HCO)  &
        -k(559)*n(idx_HCO)

    !d[MG_dot]/d[H2O+]
    pd(15,90) =  &
        -k(120)*n(idx_MG)

    !d[NH3_dot]/d[H2O+]
    pd(16,90) =  &
        -k(196)*n(idx_NH3)

    !d[NO_dot]/d[H2O+]
    pd(17,90) =  &
        -k(121)*n(idx_NO)

    !d[SI_dot]/d[H2O+]
    pd(18,90) =  &
        -k(123)*n(idx_SI)

    !d[CO_dot]/d[H2O+]
    pd(25,90) =  &
        +k(558)*n(idx_HCO)  &
        -k(554)*n(idx_CO)

    !d[NH2_dot]/d[H2O+]
    pd(27,90) =  &
        -k(782)*n(idx_NH2)  &
        -k(186)*n(idx_NH2)

    !d[CH3_dot]/d[H2O+]
    pd(28,90) =  &
        +k(454)*n(idx_CH4)

    !d[CH4_dot]/d[H2O+]
    pd(29,90) =  &
        -k(454)*n(idx_CH4)

    !d[N_dot]/d[H2O+]
    pd(30,90) =  &
        +k(798)*n(idx_NH)  &
        -k(740)*n(idx_N)  &
        -k(739)*n(idx_N)

    !d[NH_dot]/d[H2O+]
    pd(31,90) =  &
        -k(798)*n(idx_NH)

    !d[H2O_DUST_dot]/d[H2O+]
    pd(55,90) =  &
        +k(1281)

    !d[HCO+_dot]/d[H2O+]
    pd(70,90) =  &
        +k(554)*n(idx_CO)  &
        +k(119)*n(idx_HCO)

    !d[CH2+_dot]/d[H2O+]
    pd(74,90) =  &
        +k(461)*n(idx_CH)  &
        +k(41)*n(idx_CH2)

    !d[CH+_dot]/d[H2O+]
    pd(75,90) =  &
        +k(384)*n(idx_C)  &
        +k(57)*n(idx_CH)

    !d[H2CO+_dot]/d[H2O+]
    pd(76,90) =  &
        +k(118)*n(idx_H2CO)  &
        +k(559)*n(idx_HCO)

    !d[MG+_dot]/d[H2O+]
    pd(77,90) =  &
        +k(120)*n(idx_MG)

    !d[NH3+_dot]/d[H2O+]
    pd(78,90) =  &
        +k(196)*n(idx_NH3)  &
        +k(782)*n(idx_NH2)

    !d[NO+_dot]/d[H2O+]
    pd(79,90) =  &
        +k(740)*n(idx_N)  &
        +k(121)*n(idx_NO)

    !d[SI+_dot]/d[H2O+]
    pd(80,90) =  &
        +k(123)*n(idx_SI)

    !d[O2+_dot]/d[H2O+]
    pd(89,90) =  &
        +k(122)*n(idx_O2)  &
        +k(824)*n(idx_O)

    !d[H2O+_dot]/d[H2O+]
    pd(90,90) =  &
        -k(1141)  &
        -k(57)*n(idx_CH)  &
        -k(535)*n(idx_H2)  &
        -k(123)*n(idx_SI)  &
        -k(122)*n(idx_O2)  &
        -k(425)*n(idx_CH2)  &
        -k(554)*n(idx_CO)  &
        -k(558)*n(idx_HCO)  &
        -k(853)*n(idx_OH)  &
        -k(461)*n(idx_CH)  &
        -k(314)*n(idx_E)  &
        -k(313)*n(idx_E)  &
        -k(560)*n(idx_HNC)  &
        -k(186)*n(idx_NH2)  &
        -k(782)*n(idx_NH2)  &
        -k(196)*n(idx_NH3)  &
        -k(1281)  &
        -k(315)*n(idx_E)  &
        -k(118)*n(idx_H2CO)  &
        -k(740)*n(idx_N)  &
        -k(798)*n(idx_NH)  &
        -k(555)*n(idx_H2CO)  &
        -k(384)*n(idx_C)  &
        -k(119)*n(idx_HCO)  &
        -k(41)*n(idx_CH2)  &
        -k(454)*n(idx_CH4)  &
        -k(120)*n(idx_MG)  &
        -k(557)*n(idx_HCN)  &
        -k(824)*n(idx_O)  &
        -k(121)*n(idx_NO)  &
        -k(559)*n(idx_HCO)  &
        -k(556)*n(idx_H2O)  &
        -k(739)*n(idx_N)

    !d[NH2+_dot]/d[H2O+]
    pd(91,90) =  &
        +k(186)*n(idx_NH2)

    !d[OH+_dot]/d[H2O+]
    pd(93,90) =  &
        +k(1141)

    !d[CH3+_dot]/d[H2O+]
    pd(94,90) =  &
        +k(425)*n(idx_CH2)

    !d[HNO+_dot]/d[H2O+]
    pd(104,90) =  &
        +k(739)*n(idx_N)

    !d[H3CO+_dot]/d[H2O+]
    pd(107,90) =  &
        +k(555)*n(idx_H2CO)

    !d[H3O+_dot]/d[H2O+]
    pd(108,90) =  &
        +k(853)*n(idx_OH)  &
        +k(558)*n(idx_HCO)  &
        +k(798)*n(idx_NH)  &
        +k(535)*n(idx_H2)  &
        +k(454)*n(idx_CH4)  &
        +k(556)*n(idx_H2O)

    !d[HCNH+_dot]/d[H2O+]
    pd(109,90) =  &
        +k(557)*n(idx_HCN)  &
        +k(560)*n(idx_HNC)

    !d[E_dot]/d[NH2+]
    pd(1,91) =  &
        -k(342)*n(idx_E)  &
        -k(343)*n(idx_E)

    !d[CH_dot]/d[NH2+]
    pd(2,91) =  &
        -k(60)*n(idx_CH)  &
        -k(472)*n(idx_CH)

    !d[O_dot]/d[NH2+]
    pd(3,91) =  &
        -k(828)*n(idx_O)  &
        +k(778)*n(idx_O2)

    !d[HNC_dot]/d[NH2+]
    pd(4,91) =  &
        -k(776)*n(idx_HNC)

    !d[HCN_dot]/d[NH2+]
    pd(5,91) =  &
        -k(774)*n(idx_HCN)

    !d[H2_dot]/d[NH2+]
    pd(6,91) =  &
        -k(543)*n(idx_H2)

    !d[H_dot]/d[NH2+]
    pd(8,91) =  &
        +2.d0*k(342)*n(idx_E)  &
        +k(343)*n(idx_E)  &
        +k(828)*n(idx_O)  &
        +k(742)*n(idx_N)  &
        +k(543)*n(idx_H2)

    !d[H2O_dot]/d[NH2+]
    pd(9,91) =  &
        -k(772)*n(idx_H2O)  &
        -k(773)*n(idx_H2O)

    !d[OH_dot]/d[NH2+]
    pd(10,91) =  &
        +k(773)*n(idx_H2O)  &
        +k(779)*n(idx_O2)

    !d[O2_dot]/d[NH2+]
    pd(11,91) =  &
        -k(779)*n(idx_O2)  &
        -k(778)*n(idx_O2)

    !d[CH2_dot]/d[NH2+]
    pd(12,91) =  &
        -k(434)*n(idx_CH2)  &
        -k(43)*n(idx_CH2)

    !d[H2CO_dot]/d[NH2+]
    pd(13,91) =  &
        -k(771)*n(idx_H2CO)  &
        -k(770)*n(idx_H2CO)

    !d[HCO_dot]/d[NH2+]
    pd(14,91) =  &
        +k(771)*n(idx_H2CO)  &
        -k(775)*n(idx_HCO)  &
        -k(181)*n(idx_HCO)

    !d[NH3_dot]/d[NH2+]
    pd(16,91) =  &
        -k(182)*n(idx_NH3)

    !d[NO_dot]/d[NH2+]
    pd(17,91) =  &
        -k(183)*n(idx_NO)

    !d[NH2_dot]/d[NH2+]
    pd(27,91) =  &
        -k(777)*n(idx_NH2)  &
        +k(60)*n(idx_CH)  &
        +k(181)*n(idx_HCO)  &
        +k(43)*n(idx_CH2)  &
        +k(182)*n(idx_NH3)  &
        +k(183)*n(idx_NO)

    !d[N_dot]/d[NH2+]
    pd(30,91) =  &
        +k(803)*n(idx_NH)  &
        +k(342)*n(idx_E)  &
        -k(742)*n(idx_N)

    !d[NH_dot]/d[NH2+]
    pd(31,91) =  &
        +k(774)*n(idx_HCN)  &
        +k(777)*n(idx_NH2)  &
        +k(343)*n(idx_E)  &
        +k(434)*n(idx_CH2)  &
        +k(775)*n(idx_HCO)  &
        +k(776)*n(idx_HNC)  &
        +k(472)*n(idx_CH)  &
        -k(803)*n(idx_NH)  &
        +k(770)*n(idx_H2CO)  &
        +k(772)*n(idx_H2O)

    !d[NH3_DUST_dot]/d[NH2+]
    pd(60,91) =  &
        +k(1280)

    !d[HCO+_dot]/d[NH2+]
    pd(70,91) =  &
        +k(181)*n(idx_HCO)

    !d[CH2+_dot]/d[NH2+]
    pd(74,91) =  &
        +k(472)*n(idx_CH)  &
        +k(43)*n(idx_CH2)

    !d[CH+_dot]/d[NH2+]
    pd(75,91) =  &
        +k(60)*n(idx_CH)

    !d[H2CO+_dot]/d[NH2+]
    pd(76,91) =  &
        +k(775)*n(idx_HCO)

    !d[NH3+_dot]/d[NH2+]
    pd(78,91) =  &
        +k(543)*n(idx_H2)  &
        +k(803)*n(idx_NH)  &
        +k(773)*n(idx_H2O)  &
        +k(182)*n(idx_NH3)  &
        +k(771)*n(idx_H2CO)  &
        +k(777)*n(idx_NH2)

    !d[NO+_dot]/d[NH2+]
    pd(79,91) =  &
        +k(183)*n(idx_NO)

    !d[NH2+_dot]/d[NH2+]
    pd(91,91) =  &
        -k(183)*n(idx_NO)  &
        -k(1280)  &
        -k(343)*n(idx_E)  &
        -k(828)*n(idx_O)  &
        -k(43)*n(idx_CH2)  &
        -k(772)*n(idx_H2O)  &
        -k(771)*n(idx_H2CO)  &
        -k(742)*n(idx_N)  &
        -k(776)*n(idx_HNC)  &
        -k(803)*n(idx_NH)  &
        -k(778)*n(idx_O2)  &
        -k(342)*n(idx_E)  &
        -k(60)*n(idx_CH)  &
        -k(775)*n(idx_HCO)  &
        -k(472)*n(idx_CH)  &
        -k(774)*n(idx_HCN)  &
        -k(181)*n(idx_HCO)  &
        -k(770)*n(idx_H2CO)  &
        -k(777)*n(idx_NH2)  &
        -k(182)*n(idx_NH3)  &
        -k(434)*n(idx_CH2)  &
        -k(773)*n(idx_H2O)  &
        -k(543)*n(idx_H2)  &
        -k(779)*n(idx_O2)

    !d[CH3+_dot]/d[NH2+]
    pd(94,91) =  &
        +k(434)*n(idx_CH2)

    !d[HNO+_dot]/d[NH2+]
    pd(104,91) =  &
        +k(828)*n(idx_O)  &
        +k(779)*n(idx_O2)

    !d[H2NO+_dot]/d[NH2+]
    pd(105,91) =  &
        +k(778)*n(idx_O2)

    !d[H3CO+_dot]/d[NH2+]
    pd(107,91) =  &
        +k(770)*n(idx_H2CO)

    !d[H3O+_dot]/d[NH2+]
    pd(108,91) =  &
        +k(772)*n(idx_H2O)

    !d[HCNH+_dot]/d[NH2+]
    pd(109,91) =  &
        +k(774)*n(idx_HCN)  &
        +k(776)*n(idx_HNC)

    !d[N2H+_dot]/d[NH2+]
    pd(112,91) =  &
        +k(742)*n(idx_N)

    !d[E_dot]/d[O+]
    pd(1,92) =  &
        -k(1222)*n(idx_E)

    !d[CH_dot]/d[O+]
    pd(2,92) =  &
        -k(473)*n(idx_CH)  &
        +k(816)*n(idx_HCN)  &
        -k(61)*n(idx_CH)

    !d[O_dot]/d[O+]
    pd(3,92) =  &
        +k(209)*n(idx_CO)  &
        +k(214)*n(idx_NH3)  &
        +k(211)*n(idx_H2O)  &
        +k(203)*n(idx_NH)  &
        +k(212)*n(idx_HCO)  &
        +k(132)*n(idx_H)  &
        +k(1222)*n(idx_E)  &
        +k(216)*n(idx_OH)  &
        +k(213)*n(idx_NH2)  &
        +k(208)*n(idx_CH4)  &
        +k(44)*n(idx_CH2)  &
        +k(61)*n(idx_CH)  &
        +k(210)*n(idx_H2CO)  &
        +k(215)*n(idx_O2)

    !d[HCN_dot]/d[O+]
    pd(5,92) =  &
        -k(816)*n(idx_HCN)  &
        -k(815)*n(idx_HCN)

    !d[H2_dot]/d[O+]
    pd(6,92) =  &
        -k(544)*n(idx_H2)

    !d[C_dot]/d[O+]
    pd(7,92) =  &
        +k(812)*n(idx_CN)  &
        -k(1196)*n(idx_C)

    !d[H_dot]/d[O+]
    pd(8,92) =  &
        -k(132)*n(idx_H)  &
        +k(473)*n(idx_CH)  &
        +k(544)*n(idx_H2)  &
        +k(820)*n(idx_OH)  &
        +k(804)*n(idx_NH)

    !d[H2O_dot]/d[O+]
    pd(9,92) =  &
        +k(809)*n(idx_CH3OH)  &
        -k(211)*n(idx_H2O)

    !d[OH_dot]/d[O+]
    pd(10,92) =  &
        +k(814)*n(idx_H2CO)  &
        +k(811)*n(idx_CH4)  &
        -k(820)*n(idx_OH)  &
        -k(216)*n(idx_OH)  &
        +k(810)*n(idx_CH3OH)

    !d[O2_dot]/d[O+]
    pd(11,92) =  &
        -k(215)*n(idx_O2)  &
        +k(819)*n(idx_NO2)

    !d[CH2_dot]/d[O+]
    pd(12,92) =  &
        -k(44)*n(idx_CH2)

    !d[H2CO_dot]/d[O+]
    pd(13,92) =  &
        -k(210)*n(idx_H2CO)  &
        -k(814)*n(idx_H2CO)

    !d[HCO_dot]/d[O+]
    pd(14,92) =  &
        -k(212)*n(idx_HCO)  &
        -k(817)*n(idx_HCO)

    !d[NH3_dot]/d[O+]
    pd(16,92) =  &
        -k(214)*n(idx_NH3)

    !d[CN_dot]/d[O+]
    pd(24,92) =  &
        -k(812)*n(idx_CN)

    !d[CO_dot]/d[O+]
    pd(25,92) =  &
        +k(817)*n(idx_HCO)  &
        +k(813)*n(idx_CO2)  &
        -k(209)*n(idx_CO)

    !d[N2_dot]/d[O+]
    pd(26,92) =  &
        -k(818)*n(idx_N2)

    !d[NH2_dot]/d[O+]
    pd(27,92) =  &
        -k(213)*n(idx_NH2)

    !d[CH4_dot]/d[O+]
    pd(29,92) =  &
        -k(811)*n(idx_CH4)  &
        -k(208)*n(idx_CH4)

    !d[N_dot]/d[O+]
    pd(30,92) =  &
        +k(818)*n(idx_N2)  &
        +k(815)*n(idx_HCN)

    !d[NH_dot]/d[O+]
    pd(31,92) =  &
        -k(804)*n(idx_NH)  &
        -k(203)*n(idx_NH)

    !d[CH3OH_dot]/d[O+]
    pd(37,92) =  &
        -k(809)*n(idx_CH3OH)  &
        -k(810)*n(idx_CH3OH)

    !d[CO2_dot]/d[O+]
    pd(38,92) =  &
        -k(813)*n(idx_CO2)

    !d[NO2_dot]/d[O+]
    pd(42,92) =  &
        -k(819)*n(idx_NO2)

    !d[H2O_DUST_dot]/d[O+]
    pd(55,92) =  &
        +k(1270)

    !d[HCO+_dot]/d[O+]
    pd(70,92) =  &
        +k(814)*n(idx_H2CO)  &
        +k(815)*n(idx_HCN)  &
        +k(212)*n(idx_HCO)

    !d[H+_dot]/d[O+]
    pd(71,92) =  &
        +k(132)*n(idx_H)

    !d[CH2+_dot]/d[O+]
    pd(74,92) =  &
        +k(44)*n(idx_CH2)

    !d[CH+_dot]/d[O+]
    pd(75,92) =  &
        +k(61)*n(idx_CH)

    !d[H2CO+_dot]/d[O+]
    pd(76,92) =  &
        +k(210)*n(idx_H2CO)  &
        +k(809)*n(idx_CH3OH)

    !d[NH3+_dot]/d[O+]
    pd(78,92) =  &
        +k(214)*n(idx_NH3)

    !d[NO+_dot]/d[O+]
    pd(79,92) =  &
        +k(812)*n(idx_CN)  &
        +k(818)*n(idx_N2)  &
        +k(804)*n(idx_NH)  &
        +k(816)*n(idx_HCN)  &
        +k(819)*n(idx_NO2)

    !d[CO+_dot]/d[O+]
    pd(87,92) =  &
        +k(209)*n(idx_CO)  &
        +k(473)*n(idx_CH)  &
        +k(1196)*n(idx_C)

    !d[O2+_dot]/d[O+]
    pd(89,92) =  &
        +k(215)*n(idx_O2)  &
        +k(820)*n(idx_OH)  &
        +k(813)*n(idx_CO2)

    !d[H2O+_dot]/d[O+]
    pd(90,92) =  &
        +k(211)*n(idx_H2O)

    !d[NH2+_dot]/d[O+]
    pd(91,92) =  &
        +k(213)*n(idx_NH2)

    !d[O+_dot]/d[O+]
    pd(92,92) =  &
        -k(820)*n(idx_OH)  &
        -k(804)*n(idx_NH)  &
        -k(203)*n(idx_NH)  &
        -k(1270)  &
        -k(213)*n(idx_NH2)  &
        -k(816)*n(idx_HCN)  &
        -k(216)*n(idx_OH)  &
        -k(819)*n(idx_NO2)  &
        -k(815)*n(idx_HCN)  &
        -k(813)*n(idx_CO2)  &
        -k(810)*n(idx_CH3OH)  &
        -k(132)*n(idx_H)  &
        -k(812)*n(idx_CN)  &
        -k(212)*n(idx_HCO)  &
        -k(814)*n(idx_H2CO)  &
        -k(818)*n(idx_N2)  &
        -k(61)*n(idx_CH)  &
        -k(210)*n(idx_H2CO)  &
        -k(214)*n(idx_NH3)  &
        -k(544)*n(idx_H2)  &
        -k(473)*n(idx_CH)  &
        -k(209)*n(idx_CO)  &
        -k(809)*n(idx_CH3OH)  &
        -k(817)*n(idx_HCO)  &
        -k(44)*n(idx_CH2)  &
        -k(811)*n(idx_CH4)  &
        -k(1222)*n(idx_E)  &
        -k(215)*n(idx_O2)  &
        -k(208)*n(idx_CH4)  &
        -k(211)*n(idx_H2O)  &
        -k(1196)*n(idx_C)

    !d[OH+_dot]/d[O+]
    pd(93,92) =  &
        +k(544)*n(idx_H2)  &
        +k(817)*n(idx_HCO)  &
        +k(216)*n(idx_OH)

    !d[CH3+_dot]/d[O+]
    pd(94,92) =  &
        +k(811)*n(idx_CH4)

    !d[CH4+_dot]/d[O+]
    pd(95,92) =  &
        +k(208)*n(idx_CH4)

    !d[NH+_dot]/d[O+]
    pd(98,92) =  &
        +k(203)*n(idx_NH)

    !d[H3CO+_dot]/d[O+]
    pd(107,92) =  &
        +k(810)*n(idx_CH3OH)

    !d[E_dot]/d[OH+]
    pd(1,93) =  &
        -k(349)*n(idx_E)

    !d[CH_dot]/d[OH+]
    pd(2,93) =  &
        -k(63)*n(idx_CH)  &
        -k(476)*n(idx_CH)

    !d[O_dot]/d[OH+]
    pd(3,93) =  &
        +k(845)*n(idx_HNC)  &
        +k(840)*n(idx_H2CO)  &
        -k(831)*n(idx_O)  &
        +k(476)*n(idx_CH)  &
        +k(846)*n(idx_N2)  &
        +k(792)*n(idx_NH2)  &
        +k(841)*n(idx_H2O)  &
        +k(848)*n(idx_OH)  &
        +k(349)*n(idx_E)  &
        +k(849)*n(idx_SI)  &
        +k(842)*n(idx_HCN)  &
        +k(837)*n(idx_CN)  &
        +k(847)*n(idx_NO)  &
        +k(851)*n(idx_SIO)  &
        +k(438)*n(idx_CH2)  &
        +k(844)*n(idx_HCO)  &
        +k(850)*n(idx_SIH)  &
        +k(838)*n(idx_CO2)  &
        +k(807)*n(idx_NH)  &
        +k(839)*n(idx_CO)  &
        +k(394)*n(idx_C)

    !d[HNC_dot]/d[OH+]
    pd(4,93) =  &
        -k(845)*n(idx_HNC)

    !d[HCN_dot]/d[OH+]
    pd(5,93) =  &
        -k(842)*n(idx_HCN)

    !d[H2_dot]/d[OH+]
    pd(6,93) =  &
        -k(546)*n(idx_H2)

    !d[C_dot]/d[OH+]
    pd(7,93) =  &
        -k(394)*n(idx_C)

    !d[H_dot]/d[OH+]
    pd(8,93) =  &
        +k(831)*n(idx_O)  &
        +k(1174)  &
        +k(546)*n(idx_H2)  &
        +k(349)*n(idx_E)  &
        +k(744)*n(idx_N)

    !d[H2O_dot]/d[OH+]
    pd(9,93) =  &
        -k(841)*n(idx_H2O)  &
        -k(221)*n(idx_H2O)

    !d[OH_dot]/d[OH+]
    pd(10,93) =  &
        +k(223)*n(idx_NH3)  &
        +k(220)*n(idx_H2CO)  &
        +k(224)*n(idx_NO)  &
        +k(222)*n(idx_HCO)  &
        -k(848)*n(idx_OH)  &
        +k(46)*n(idx_CH2)  &
        +k(221)*n(idx_H2O)  &
        +k(189)*n(idx_NH2)  &
        +k(225)*n(idx_O2)  &
        +k(63)*n(idx_CH)

    !d[O2_dot]/d[OH+]
    pd(11,93) =  &
        -k(225)*n(idx_O2)

    !d[CH2_dot]/d[OH+]
    pd(12,93) =  &
        +k(458)*n(idx_CH4)  &
        -k(438)*n(idx_CH2)  &
        -k(46)*n(idx_CH2)

    !d[H2CO_dot]/d[OH+]
    pd(13,93) =  &
        -k(220)*n(idx_H2CO)  &
        -k(840)*n(idx_H2CO)

    !d[HCO_dot]/d[OH+]
    pd(14,93) =  &
        -k(843)*n(idx_HCO)  &
        -k(222)*n(idx_HCO)  &
        -k(844)*n(idx_HCO)

    !d[NH3_dot]/d[OH+]
    pd(16,93) =  &
        -k(223)*n(idx_NH3)

    !d[NO_dot]/d[OH+]
    pd(17,93) =  &
        -k(847)*n(idx_NO)  &
        -k(224)*n(idx_NO)

    !d[SI_dot]/d[OH+]
    pd(18,93) =  &
        -k(849)*n(idx_SI)

    !d[CN_dot]/d[OH+]
    pd(24,93) =  &
        -k(837)*n(idx_CN)

    !d[CO_dot]/d[OH+]
    pd(25,93) =  &
        -k(839)*n(idx_CO)  &
        +k(843)*n(idx_HCO)

    !d[N2_dot]/d[OH+]
    pd(26,93) =  &
        -k(846)*n(idx_N2)

    !d[NH2_dot]/d[OH+]
    pd(27,93) =  &
        -k(792)*n(idx_NH2)  &
        -k(189)*n(idx_NH2)

    !d[CH4_dot]/d[OH+]
    pd(29,93) =  &
        -k(458)*n(idx_CH4)

    !d[N_dot]/d[OH+]
    pd(30,93) =  &
        -k(744)*n(idx_N)

    !d[NH_dot]/d[OH+]
    pd(31,93) =  &
        -k(807)*n(idx_NH)

    !d[SIH_dot]/d[OH+]
    pd(33,93) =  &
        -k(850)*n(idx_SIH)

    !d[SIO_dot]/d[OH+]
    pd(34,93) =  &
        -k(851)*n(idx_SIO)

    !d[CO2_dot]/d[OH+]
    pd(38,93) =  &
        -k(838)*n(idx_CO2)

    !d[H2O_DUST_dot]/d[OH+]
    pd(55,93) =  &
        +k(1275)

    !d[HCO+_dot]/d[OH+]
    pd(70,93) =  &
        +k(839)*n(idx_CO)  &
        +k(222)*n(idx_HCO)

    !d[CH2+_dot]/d[OH+]
    pd(74,93) =  &
        +k(46)*n(idx_CH2)  &
        +k(476)*n(idx_CH)

    !d[CH+_dot]/d[OH+]
    pd(75,93) =  &
        +k(394)*n(idx_C)  &
        +k(63)*n(idx_CH)

    !d[H2CO+_dot]/d[OH+]
    pd(76,93) =  &
        +k(220)*n(idx_H2CO)  &
        +k(844)*n(idx_HCO)

    !d[NH3+_dot]/d[OH+]
    pd(78,93) =  &
        +k(792)*n(idx_NH2)  &
        +k(223)*n(idx_NH3)

    !d[NO+_dot]/d[OH+]
    pd(79,93) =  &
        +k(224)*n(idx_NO)  &
        +k(744)*n(idx_N)

    !d[SIH2+_dot]/d[OH+]
    pd(84,93) =  &
        +k(850)*n(idx_SIH)

    !d[O2+_dot]/d[OH+]
    pd(89,93) =  &
        +k(225)*n(idx_O2)  &
        +k(831)*n(idx_O)

    !d[H2O+_dot]/d[OH+]
    pd(90,93) =  &
        +k(546)*n(idx_H2)  &
        +k(843)*n(idx_HCO)  &
        +k(221)*n(idx_H2O)  &
        +k(848)*n(idx_OH)

    !d[NH2+_dot]/d[OH+]
    pd(91,93) =  &
        +k(189)*n(idx_NH2)  &
        +k(807)*n(idx_NH)

    !d[O+_dot]/d[OH+]
    pd(92,93) =  &
        +k(1174)

    !d[OH+_dot]/d[OH+]
    pd(93,93) =  &
        -k(189)*n(idx_NH2)  &
        -k(851)*n(idx_SIO)  &
        -k(546)*n(idx_H2)  &
        -k(349)*n(idx_E)  &
        -k(792)*n(idx_NH2)  &
        -k(844)*n(idx_HCO)  &
        -k(476)*n(idx_CH)  &
        -k(843)*n(idx_HCO)  &
        -k(839)*n(idx_CO)  &
        -k(222)*n(idx_HCO)  &
        -k(221)*n(idx_H2O)  &
        -k(831)*n(idx_O)  &
        -k(845)*n(idx_HNC)  &
        -k(223)*n(idx_NH3)  &
        -k(842)*n(idx_HCN)  &
        -k(63)*n(idx_CH)  &
        -k(840)*n(idx_H2CO)  &
        -k(458)*n(idx_CH4)  &
        -k(807)*n(idx_NH)  &
        -k(850)*n(idx_SIH)  &
        -k(438)*n(idx_CH2)  &
        -k(394)*n(idx_C)  &
        -k(841)*n(idx_H2O)  &
        -k(846)*n(idx_N2)  &
        -k(744)*n(idx_N)  &
        -k(848)*n(idx_OH)  &
        -k(849)*n(idx_SI)  &
        -k(220)*n(idx_H2CO)  &
        -k(46)*n(idx_CH2)  &
        -k(837)*n(idx_CN)  &
        -k(1174)  &
        -k(224)*n(idx_NO)  &
        -k(847)*n(idx_NO)  &
        -k(838)*n(idx_CO2)  &
        -k(1275)  &
        -k(225)*n(idx_O2)

    !d[CH3+_dot]/d[OH+]
    pd(94,93) =  &
        +k(438)*n(idx_CH2)

    !d[HCN+_dot]/d[OH+]
    pd(97,93) =  &
        +k(837)*n(idx_CN)

    !d[SIH+_dot]/d[OH+]
    pd(100,93) =  &
        +k(849)*n(idx_SI)

    !d[HNO+_dot]/d[OH+]
    pd(104,93) =  &
        +k(847)*n(idx_NO)

    !d[H3CO+_dot]/d[OH+]
    pd(107,93) =  &
        +k(840)*n(idx_H2CO)

    !d[H3O+_dot]/d[OH+]
    pd(108,93) =  &
        +k(841)*n(idx_H2O)  &
        +k(458)*n(idx_CH4)

    !d[HCNH+_dot]/d[OH+]
    pd(109,93) =  &
        +k(845)*n(idx_HNC)  &
        +k(842)*n(idx_HCN)

    !d[HCO2+_dot]/d[OH+]
    pd(110,93) =  &
        +k(838)*n(idx_CO2)

    !d[N2H+_dot]/d[OH+]
    pd(112,93) =  &
        +k(846)*n(idx_N2)

    !d[SIOH+_dot]/d[OH+]
    pd(115,93) =  &
        +k(851)*n(idx_SIO)

    !d[E_dot]/d[CH3+]
    pd(1,94) =  &
        -k(299)*n(idx_E)  &
        -k(300)*n(idx_E)  &
        -k(301)*n(idx_E)  &
        -k(1216)*n(idx_E)

    !d[CH_dot]/d[CH3+]
    pd(2,94) =  &
        +k(300)*n(idx_E)  &
        +k(301)*n(idx_E)

    !d[O_dot]/d[CH3+]
    pd(3,94) =  &
        -k(445)*n(idx_O)  &
        +k(443)*n(idx_O2)  &
        -k(444)*n(idx_O)

    !d[H2_dot]/d[CH3+]
    pd(6,94) =  &
        +k(445)*n(idx_O)  &
        +k(617)*n(idx_H)  &
        +k(1115)  &
        +k(795)*n(idx_NH)  &
        +k(300)*n(idx_E)  &
        +k(446)*n(idx_OH)

    !d[H_dot]/d[CH3+]
    pd(8,94) =  &
        +2.d0*k(301)*n(idx_E)  &
        -k(617)*n(idx_H)  &
        +k(444)*n(idx_O)  &
        +k(1116)  &
        +k(299)*n(idx_E)

    !d[OH_dot]/d[CH3+]
    pd(10,94) =  &
        -k(446)*n(idx_OH)

    !d[O2_dot]/d[CH3+]
    pd(11,94) =  &
        -k(443)*n(idx_O2)

    !d[CH2_dot]/d[CH3+]
    pd(12,94) =  &
        +k(299)*n(idx_E)

    !d[H2CO_dot]/d[CH3+]
    pd(13,94) =  &
        -k(441)*n(idx_H2CO)

    !d[HCO_dot]/d[CH3+]
    pd(14,94) =  &
        -k(442)*n(idx_HCO)  &
        -k(47)*n(idx_HCO)

    !d[MG_dot]/d[CH3+]
    pd(15,94) =  &
        -k(48)*n(idx_MG)

    !d[NO_dot]/d[CH3+]
    pd(17,94) =  &
        -k(49)*n(idx_NO)

    !d[CO_dot]/d[CH3+]
    pd(25,94) =  &
        +k(442)*n(idx_HCO)

    !d[CH3_dot]/d[CH3+]
    pd(28,94) =  &
        +k(47)*n(idx_HCO)  &
        +k(48)*n(idx_MG)  &
        +k(1216)*n(idx_E)  &
        +k(49)*n(idx_NO)

    !d[CH4_dot]/d[CH3+]
    pd(29,94) =  &
        +k(447)*n(idx_SIH4)  &
        +k(441)*n(idx_H2CO)  &
        +k(440)*n(idx_CH3OH)

    !d[NH_dot]/d[CH3+]
    pd(31,94) =  &
        -k(795)*n(idx_NH)

    !d[SIH4_dot]/d[CH3+]
    pd(32,94) =  &
        -k(447)*n(idx_SIH4)

    !d[CH3OH_dot]/d[CH3+]
    pd(37,94) =  &
        -k(440)*n(idx_CH3OH)

    !d[CH4_DUST_dot]/d[CH3+]
    pd(53,94) =  &
        +k(1286)

    !d[HCO+_dot]/d[CH3+]
    pd(70,94) =  &
        +k(47)*n(idx_HCO)  &
        +k(445)*n(idx_O)  &
        +k(441)*n(idx_H2CO)

    !d[CH2+_dot]/d[CH3+]
    pd(74,94) =  &
        +k(1116)  &
        +k(617)*n(idx_H)

    !d[CH+_dot]/d[CH3+]
    pd(75,94) =  &
        +k(1115)

    !d[H2CO+_dot]/d[CH3+]
    pd(76,94) =  &
        +k(444)*n(idx_O)  &
        +k(446)*n(idx_OH)

    !d[MG+_dot]/d[CH3+]
    pd(77,94) =  &
        +k(48)*n(idx_MG)

    !d[NO+_dot]/d[CH3+]
    pd(79,94) =  &
        +k(49)*n(idx_NO)

    !d[SIH3+_dot]/d[CH3+]
    pd(85,94) =  &
        +k(447)*n(idx_SIH4)

    !d[CH3+_dot]/d[CH3+]
    pd(94,94) =  &
        -k(1216)*n(idx_E)  &
        -k(443)*n(idx_O2)  &
        -k(795)*n(idx_NH)  &
        -k(49)*n(idx_NO)  &
        -k(444)*n(idx_O)  &
        -k(446)*n(idx_OH)  &
        -k(299)*n(idx_E)  &
        -k(300)*n(idx_E)  &
        -k(617)*n(idx_H)  &
        -k(1115)  &
        -k(447)*n(idx_SIH4)  &
        -k(301)*n(idx_E)  &
        -k(48)*n(idx_MG)  &
        -k(1116)  &
        -k(440)*n(idx_CH3OH)  &
        -k(442)*n(idx_HCO)  &
        -k(441)*n(idx_H2CO)  &
        -k(445)*n(idx_O)  &
        -k(1286)  &
        -k(47)*n(idx_HCO)

    !d[CH4+_dot]/d[CH3+]
    pd(95,94) =  &
        +k(442)*n(idx_HCO)

    !d[H3CO+_dot]/d[CH3+]
    pd(107,94) =  &
        +k(443)*n(idx_O2)  &
        +k(440)*n(idx_CH3OH)

    !d[HCNH+_dot]/d[CH3+]
    pd(109,94) =  &
        +k(795)*n(idx_NH)

    !d[E_dot]/d[CH4+]
    pd(1,95) =  &
        -k(302)*n(idx_E)  &
        -k(303)*n(idx_E)

    !d[O_dot]/d[CH4+]
    pd(3,95) =  &
        -k(823)*n(idx_O)

    !d[H2_dot]/d[CH4+]
    pd(6,95) =  &
        +k(618)*n(idx_H)  &
        +k(1123)

    !d[H_dot]/d[CH4+]
    pd(8,95) =  &
        -k(618)*n(idx_H)  &
        +k(1124)  &
        +2.d0*k(302)*n(idx_E)  &
        +k(303)*n(idx_E)

    !d[H2O_dot]/d[CH4+]
    pd(9,95) =  &
        -k(451)*n(idx_H2O)

    !d[OH_dot]/d[CH4+]
    pd(10,95) =  &
        +k(823)*n(idx_O)

    !d[O2_dot]/d[CH4+]
    pd(11,95) =  &
        -k(52)*n(idx_O2)

    !d[CH2_dot]/d[CH4+]
    pd(12,95) =  &
        +k(302)*n(idx_E)

    !d[H2CO_dot]/d[CH4+]
    pd(13,95) =  &
        -k(50)*n(idx_H2CO)  &
        -k(450)*n(idx_H2CO)

    !d[NH3_dot]/d[CH4+]
    pd(16,95) =  &
        -k(51)*n(idx_NH3)

    !d[CO_dot]/d[CH4+]
    pd(25,95) =  &
        -k(449)*n(idx_CO)

    !d[CH3_dot]/d[CH4+]
    pd(28,95) =  &
        +k(448)*n(idx_CO2)  &
        +k(451)*n(idx_H2O)  &
        +k(450)*n(idx_H2CO)  &
        +k(449)*n(idx_CO)  &
        +k(303)*n(idx_E)

    !d[CH4_dot]/d[CH4+]
    pd(29,95) =  &
        +k(52)*n(idx_O2)  &
        +k(51)*n(idx_NH3)  &
        +k(50)*n(idx_H2CO)

    !d[CO2_dot]/d[CH4+]
    pd(38,95) =  &
        -k(448)*n(idx_CO2)

    !d[CH4_DUST_dot]/d[CH4+]
    pd(53,95) =  &
        +k(1289)

    !d[HCO+_dot]/d[CH4+]
    pd(70,95) =  &
        +k(449)*n(idx_CO)

    !d[CH2+_dot]/d[CH4+]
    pd(74,95) =  &
        +k(1123)

    !d[H2CO+_dot]/d[CH4+]
    pd(76,95) =  &
        +k(50)*n(idx_H2CO)

    !d[NH3+_dot]/d[CH4+]
    pd(78,95) =  &
        +k(51)*n(idx_NH3)

    !d[O2+_dot]/d[CH4+]
    pd(89,95) =  &
        +k(52)*n(idx_O2)

    !d[CH3+_dot]/d[CH4+]
    pd(94,95) =  &
        +k(823)*n(idx_O)  &
        +k(618)*n(idx_H)  &
        +k(1124)

    !d[CH4+_dot]/d[CH4+]
    pd(95,95) =  &
        -k(448)*n(idx_CO2)  &
        -k(823)*n(idx_O)  &
        -k(449)*n(idx_CO)  &
        -k(618)*n(idx_H)  &
        -k(50)*n(idx_H2CO)  &
        -k(302)*n(idx_E)  &
        -k(303)*n(idx_E)  &
        -k(51)*n(idx_NH3)  &
        -k(1289)  &
        -k(1123)  &
        -k(451)*n(idx_H2O)  &
        -k(52)*n(idx_O2)  &
        -k(1124)  &
        -k(450)*n(idx_H2CO)

    !d[H3CO+_dot]/d[CH4+]
    pd(107,95) =  &
        +k(450)*n(idx_H2CO)

    !d[H3O+_dot]/d[CH4+]
    pd(108,95) =  &
        +k(451)*n(idx_H2O)

    !d[HCO2+_dot]/d[CH4+]
    pd(110,95) =  &
        +k(448)*n(idx_CO2)

    !d[E_dot]/d[N+]
    pd(1,96) =  &
        -k(1221)*n(idx_E)

    !d[CH_dot]/d[N+]
    pd(2,96) =  &
        -k(58)*n(idx_CH)  &
        -k(469)*n(idx_CH)

    !d[O_dot]/d[N+]
    pd(3,96) =  &
        +k(729)*n(idx_O2)  &
        +k(728)*n(idx_NO)

    !d[HCN_dot]/d[N+]
    pd(5,96) =  &
        -k(162)*n(idx_HCN)

    !d[H2_dot]/d[N+]
    pd(6,96) =  &
        -k(539)*n(idx_H2)  &
        +k(718)*n(idx_CH4)  &
        +k(725)*n(idx_NH3)

    !d[C_dot]/d[N+]
    pd(7,96) =  &
        +k(721)*n(idx_CO)

    !d[H_dot]/d[N+]
    pd(8,96) =  &
        +2.d0*k(719)*n(idx_CH4)  &
        +k(539)*n(idx_H2)  &
        +k(716)*n(idx_CH3OH)  &
        +k(727)*n(idx_NH)  &
        +k(469)*n(idx_CH)  &
        +k(717)*n(idx_CH4)  &
        +k(713)*n(idx_CH3OH)  &
        +k(718)*n(idx_CH4)  &
        +k(715)*n(idx_CH3OH)

    !d[H2O_dot]/d[N+]
    pd(9,96) =  &
        -k(161)*n(idx_H2O)

    !d[OH_dot]/d[N+]
    pd(10,96) =  &
        -k(170)*n(idx_OH)

    !d[O2_dot]/d[N+]
    pd(11,96) =  &
        -k(730)*n(idx_O2)  &
        -k(729)*n(idx_O2)  &
        -k(169)*n(idx_O2)

    !d[CH2_dot]/d[N+]
    pd(12,96) =  &
        +k(723)*n(idx_H2CO)  &
        -k(156)*n(idx_CH2)

    !d[H2CO_dot]/d[N+]
    pd(13,96) =  &
        -k(160)*n(idx_H2CO)  &
        -k(722)*n(idx_H2CO)  &
        -k(723)*n(idx_H2CO)

    !d[HCO_dot]/d[N+]
    pd(14,96) =  &
        -k(163)*n(idx_HCO)  &
        -k(724)*n(idx_HCO)

    !d[MG_dot]/d[N+]
    pd(15,96) =  &
        -k(164)*n(idx_MG)

    !d[NH3_dot]/d[N+]
    pd(16,96) =  &
        -k(726)*n(idx_NH3)  &
        -k(166)*n(idx_NH3)  &
        -k(725)*n(idx_NH3)

    !d[NO_dot]/d[N+]
    pd(17,96) =  &
        +k(730)*n(idx_O2)  &
        +k(716)*n(idx_CH3OH)  &
        +k(720)*n(idx_CO2)  &
        -k(728)*n(idx_NO)  &
        -k(168)*n(idx_NO)

    !d[CN_dot]/d[N+]
    pd(24,96) =  &
        -k(158)*n(idx_CN)

    !d[CO_dot]/d[N+]
    pd(25,96) =  &
        +k(724)*n(idx_HCO)  &
        -k(159)*n(idx_CO)  &
        -k(721)*n(idx_CO)

    !d[NH2_dot]/d[N+]
    pd(27,96) =  &
        -k(165)*n(idx_NH2)

    !d[CH3_dot]/d[N+]
    pd(28,96) =  &
        +k(715)*n(idx_CH3OH)

    !d[CH4_dot]/d[N+]
    pd(29,96) =  &
        -k(157)*n(idx_CH4)  &
        -k(719)*n(idx_CH4)  &
        -k(717)*n(idx_CH4)  &
        -k(718)*n(idx_CH4)

    !d[N_dot]/d[N+]
    pd(30,96) =  &
        +k(165)*n(idx_NH2)  &
        +k(167)*n(idx_NH)  &
        +k(170)*n(idx_OH)  &
        +k(160)*n(idx_H2CO)  &
        +k(157)*n(idx_CH4)  &
        +k(168)*n(idx_NO)  &
        +k(58)*n(idx_CH)  &
        +k(156)*n(idx_CH2)  &
        +k(162)*n(idx_HCN)  &
        +k(159)*n(idx_CO)  &
        +k(1221)*n(idx_E)  &
        +k(161)*n(idx_H2O)  &
        +k(169)*n(idx_O2)  &
        +k(158)*n(idx_CN)  &
        +k(164)*n(idx_MG)  &
        +k(717)*n(idx_CH4)  &
        +k(163)*n(idx_HCO)  &
        +k(166)*n(idx_NH3)  &
        -k(1211)*n(idx_N)

    !d[NH_dot]/d[N+]
    pd(31,96) =  &
        +k(714)*n(idx_CH3OH)  &
        -k(727)*n(idx_NH)  &
        +k(722)*n(idx_H2CO)  &
        +k(726)*n(idx_NH3)  &
        +k(713)*n(idx_CH3OH)  &
        -k(167)*n(idx_NH)

    !d[CH3OH_dot]/d[N+]
    pd(37,96) =  &
        -k(714)*n(idx_CH3OH)  &
        -k(713)*n(idx_CH3OH)  &
        -k(716)*n(idx_CH3OH)  &
        -k(715)*n(idx_CH3OH)

    !d[CO2_dot]/d[N+]
    pd(38,96) =  &
        -k(720)*n(idx_CO2)

    !d[NH3_DUST_dot]/d[N+]
    pd(60,96) =  &
        +k(1269)

    !d[HCO+_dot]/d[N+]
    pd(70,96) =  &
        +k(163)*n(idx_HCO)  &
        +k(722)*n(idx_H2CO)

    !d[CH2+_dot]/d[N+]
    pd(74,96) =  &
        +k(156)*n(idx_CH2)

    !d[CH+_dot]/d[N+]
    pd(75,96) =  &
        +k(58)*n(idx_CH)

    !d[H2CO+_dot]/d[N+]
    pd(76,96) =  &
        +k(713)*n(idx_CH3OH)  &
        +k(160)*n(idx_H2CO)

    !d[MG+_dot]/d[N+]
    pd(77,96) =  &
        +k(164)*n(idx_MG)

    !d[NH3+_dot]/d[N+]
    pd(78,96) =  &
        +k(166)*n(idx_NH3)

    !d[NO+_dot]/d[N+]
    pd(79,96) =  &
        +k(168)*n(idx_NO)  &
        +k(729)*n(idx_O2)  &
        +k(721)*n(idx_CO)  &
        +k(715)*n(idx_CH3OH)  &
        +k(723)*n(idx_H2CO)

    !d[CN+_dot]/d[N+]
    pd(86,96) =  &
        +k(158)*n(idx_CN)  &
        +k(469)*n(idx_CH)

    !d[CO+_dot]/d[N+]
    pd(87,96) =  &
        +k(159)*n(idx_CO)  &
        +k(720)*n(idx_CO2)

    !d[N2+_dot]/d[N+]
    pd(88,96) =  &
        +k(1211)*n(idx_N)  &
        +k(727)*n(idx_NH)  &
        +k(728)*n(idx_NO)

    !d[O2+_dot]/d[N+]
    pd(89,96) =  &
        +k(169)*n(idx_O2)

    !d[H2O+_dot]/d[N+]
    pd(90,96) =  &
        +k(161)*n(idx_H2O)

    !d[NH2+_dot]/d[N+]
    pd(91,96) =  &
        +k(165)*n(idx_NH2)  &
        +k(726)*n(idx_NH3)

    !d[O+_dot]/d[N+]
    pd(92,96) =  &
        +k(730)*n(idx_O2)

    !d[OH+_dot]/d[N+]
    pd(93,96) =  &
        +k(170)*n(idx_OH)

    !d[CH3+_dot]/d[N+]
    pd(94,96) =  &
        +k(716)*n(idx_CH3OH)  &
        +k(717)*n(idx_CH4)

    !d[CH4+_dot]/d[N+]
    pd(95,96) =  &
        +k(157)*n(idx_CH4)

    !d[N+_dot]/d[N+]
    pd(96,96) =  &
        -k(163)*n(idx_HCO)  &
        -k(719)*n(idx_CH4)  &
        -k(720)*n(idx_CO2)  &
        -k(165)*n(idx_NH2)  &
        -k(168)*n(idx_NO)  &
        -k(730)*n(idx_O2)  &
        -k(714)*n(idx_CH3OH)  &
        -k(159)*n(idx_CO)  &
        -k(1269)  &
        -k(160)*n(idx_H2CO)  &
        -k(727)*n(idx_NH)  &
        -k(1221)*n(idx_E)  &
        -k(169)*n(idx_O2)  &
        -k(1211)*n(idx_N)  &
        -k(156)*n(idx_CH2)  &
        -k(729)*n(idx_O2)  &
        -k(170)*n(idx_OH)  &
        -k(166)*n(idx_NH3)  &
        -k(726)*n(idx_NH3)  &
        -k(58)*n(idx_CH)  &
        -k(717)*n(idx_CH4)  &
        -k(723)*n(idx_H2CO)  &
        -k(158)*n(idx_CN)  &
        -k(728)*n(idx_NO)  &
        -k(722)*n(idx_H2CO)  &
        -k(167)*n(idx_NH)  &
        -k(157)*n(idx_CH4)  &
        -k(721)*n(idx_CO)  &
        -k(162)*n(idx_HCN)  &
        -k(469)*n(idx_CH)  &
        -k(718)*n(idx_CH4)  &
        -k(725)*n(idx_NH3)  &
        -k(713)*n(idx_CH3OH)  &
        -k(716)*n(idx_CH3OH)  &
        -k(715)*n(idx_CH3OH)  &
        -k(539)*n(idx_H2)  &
        -k(724)*n(idx_HCO)  &
        -k(164)*n(idx_MG)  &
        -k(161)*n(idx_H2O)

    !d[HCN+_dot]/d[N+]
    pd(97,96) =  &
        +k(162)*n(idx_HCN)  &
        +k(718)*n(idx_CH4)

    !d[NH+_dot]/d[N+]
    pd(98,96) =  &
        +k(167)*n(idx_NH)  &
        +k(539)*n(idx_H2)  &
        +k(724)*n(idx_HCO)

    !d[H3CO+_dot]/d[N+]
    pd(107,96) =  &
        +k(714)*n(idx_CH3OH)

    !d[HCNH+_dot]/d[N+]
    pd(109,96) =  &
        +k(719)*n(idx_CH4)

    !d[N2H+_dot]/d[N+]
    pd(112,96) =  &
        +k(725)*n(idx_NH3)

    !d[E_dot]/d[HCN+]
    pd(1,97) =  &
        -k(327)*n(idx_E)

    !d[CH_dot]/d[HCN+]
    pd(2,97) =  &
        -k(464)*n(idx_CH)

    !d[HNC_dot]/d[HCN+]
    pd(4,97) =  &
        -k(627)*n(idx_HNC)

    !d[HCN_dot]/d[HCN+]
    pd(5,97) =  &
        +k(197)*n(idx_NH3)  &
        +k(125)*n(idx_H2O)  &
        +k(130)*n(idx_H)  &
        +k(134)*n(idx_O2)  &
        +k(133)*n(idx_NO)  &
        -k(624)*n(idx_HCN)

    !d[H2_dot]/d[HCN+]
    pd(6,97) =  &
        -k(536)*n(idx_H2)

    !d[C_dot]/d[HCN+]
    pd(7,97) =  &
        -k(386)*n(idx_C)

    !d[H_dot]/d[HCN+]
    pd(8,97) =  &
        -k(130)*n(idx_H)  &
        +k(536)*n(idx_H2)  &
        +k(327)*n(idx_E)

    !d[H2O_dot]/d[HCN+]
    pd(9,97) =  &
        -k(125)*n(idx_H2O)  &
        -k(566)*n(idx_H2O)

    !d[OH_dot]/d[HCN+]
    pd(10,97) =  &
        -k(854)*n(idx_OH)

    !d[O2_dot]/d[HCN+]
    pd(11,97) =  &
        -k(134)*n(idx_O2)

    !d[CH2_dot]/d[HCN+]
    pd(12,97) =  &
        -k(427)*n(idx_CH2)

    !d[H2CO_dot]/d[HCN+]
    pd(13,97) =  &
        -k(623)*n(idx_H2CO)

    !d[HCO_dot]/d[HCN+]
    pd(14,97) =  &
        -k(626)*n(idx_HCO)  &
        -k(625)*n(idx_HCO)

    !d[NH3_dot]/d[HCN+]
    pd(16,97) =  &
        -k(197)*n(idx_NH3)  &
        -k(794)*n(idx_NH3)

    !d[NO_dot]/d[HCN+]
    pd(17,97) =  &
        -k(133)*n(idx_NO)

    !d[CN_dot]/d[HCN+]
    pd(24,97) =  &
        +k(624)*n(idx_HCN)  &
        +k(386)*n(idx_C)  &
        +k(464)*n(idx_CH)  &
        +k(623)*n(idx_H2CO)  &
        +k(566)*n(idx_H2O)  &
        +k(427)*n(idx_CH2)  &
        +k(799)*n(idx_NH)  &
        +k(622)*n(idx_CO)  &
        +k(785)*n(idx_NH2)  &
        +k(627)*n(idx_HNC)  &
        +k(327)*n(idx_E)  &
        +k(621)*n(idx_CO2)  &
        +k(854)*n(idx_OH)  &
        +k(625)*n(idx_HCO)

    !d[CO_dot]/d[HCN+]
    pd(25,97) =  &
        -k(622)*n(idx_CO)  &
        +k(626)*n(idx_HCO)

    !d[NH2_dot]/d[HCN+]
    pd(27,97) =  &
        +k(794)*n(idx_NH3)  &
        -k(785)*n(idx_NH2)

    !d[CH3_dot]/d[HCN+]
    pd(28,97) =  &
        +k(455)*n(idx_CH4)

    !d[CH4_dot]/d[HCN+]
    pd(29,97) =  &
        -k(455)*n(idx_CH4)

    !d[NH_dot]/d[HCN+]
    pd(31,97) =  &
        -k(799)*n(idx_NH)

    !d[CO2_dot]/d[HCN+]
    pd(38,97) =  &
        -k(621)*n(idx_CO2)

    !d[HCN_DUST_dot]/d[HCN+]
    pd(59,97) =  &
        +k(1283)

    !d[HCO+_dot]/d[HCN+]
    pd(70,97) =  &
        +k(622)*n(idx_CO)

    !d[H+_dot]/d[HCN+]
    pd(71,97) =  &
        +k(130)*n(idx_H)

    !d[CH2+_dot]/d[HCN+]
    pd(74,97) =  &
        +k(464)*n(idx_CH)

    !d[CH+_dot]/d[HCN+]
    pd(75,97) =  &
        +k(386)*n(idx_C)

    !d[H2CO+_dot]/d[HCN+]
    pd(76,97) =  &
        +k(625)*n(idx_HCO)

    !d[NH3+_dot]/d[HCN+]
    pd(78,97) =  &
        +k(197)*n(idx_NH3)  &
        +k(785)*n(idx_NH2)

    !d[NO+_dot]/d[HCN+]
    pd(79,97) =  &
        +k(133)*n(idx_NO)

    !d[O2+_dot]/d[HCN+]
    pd(89,97) =  &
        +k(134)*n(idx_O2)

    !d[H2O+_dot]/d[HCN+]
    pd(90,97) =  &
        +k(854)*n(idx_OH)  &
        +k(125)*n(idx_H2O)

    !d[NH2+_dot]/d[HCN+]
    pd(91,97) =  &
        +k(799)*n(idx_NH)

    !d[CH3+_dot]/d[HCN+]
    pd(94,97) =  &
        +k(427)*n(idx_CH2)

    !d[HCN+_dot]/d[HCN+]
    pd(97,97) =  &
        -k(125)*n(idx_H2O)  &
        -k(197)*n(idx_NH3)  &
        -k(427)*n(idx_CH2)  &
        -k(794)*n(idx_NH3)  &
        -k(455)*n(idx_CH4)  &
        -k(566)*n(idx_H2O)  &
        -k(464)*n(idx_CH)  &
        -k(130)*n(idx_H)  &
        -k(386)*n(idx_C)  &
        -k(621)*n(idx_CO2)  &
        -k(625)*n(idx_HCO)  &
        -k(536)*n(idx_H2)  &
        -k(626)*n(idx_HCO)  &
        -k(627)*n(idx_HNC)  &
        -k(623)*n(idx_H2CO)  &
        -k(799)*n(idx_NH)  &
        -k(1283)  &
        -k(622)*n(idx_CO)  &
        -k(134)*n(idx_O2)  &
        -k(785)*n(idx_NH2)  &
        -k(854)*n(idx_OH)  &
        -k(133)*n(idx_NO)  &
        -k(327)*n(idx_E)  &
        -k(624)*n(idx_HCN)

    !d[H3CO+_dot]/d[HCN+]
    pd(107,97) =  &
        +k(623)*n(idx_H2CO)

    !d[H3O+_dot]/d[HCN+]
    pd(108,97) =  &
        +k(566)*n(idx_H2O)

    !d[HCNH+_dot]/d[HCN+]
    pd(109,97) =  &
        +k(794)*n(idx_NH3)  &
        +k(536)*n(idx_H2)  &
        +k(627)*n(idx_HNC)  &
        +k(626)*n(idx_HCO)  &
        +k(624)*n(idx_HCN)  &
        +k(455)*n(idx_CH4)

    !d[HCO2+_dot]/d[HCN+]
    pd(110,97) =  &
        +k(621)*n(idx_CO2)

    !d[E_dot]/d[NH+]
    pd(1,98) =  &
        -k(341)*n(idx_E)

    !d[CH_dot]/d[NH+]
    pd(2,98) =  &
        -k(471)*n(idx_CH)

    !d[O_dot]/d[NH+]
    pd(3,98) =  &
        +k(765)*n(idx_NO)  &
        -k(768)*n(idx_O)  &
        +k(757)*n(idx_H2O)

    !d[HNC_dot]/d[NH+]
    pd(4,98) =  &
        -k(761)*n(idx_HNC)

    !d[HCN_dot]/d[NH+]
    pd(5,98) =  &
        -k(759)*n(idx_HCN)

    !d[H2_dot]/d[NH+]
    pd(6,98) =  &
        +k(756)*n(idx_H2O)  &
        -k(542)*n(idx_H2)  &
        -k(541)*n(idx_H2)

    !d[C_dot]/d[NH+]
    pd(7,98) =  &
        -k(391)*n(idx_C)

    !d[H_dot]/d[NH+]
    pd(8,98) =  &
        +k(542)*n(idx_H2)  &
        +k(741)*n(idx_N)  &
        +k(341)*n(idx_E)

    !d[H2O_dot]/d[NH+]
    pd(9,98) =  &
        -k(757)*n(idx_H2O)  &
        -k(177)*n(idx_H2O)  &
        -k(756)*n(idx_H2O)  &
        -k(755)*n(idx_H2O)  &
        -k(758)*n(idx_H2O)

    !d[OH_dot]/d[NH+]
    pd(10,98) =  &
        +k(758)*n(idx_H2O)  &
        -k(769)*n(idx_OH)  &
        +k(766)*n(idx_O2)

    !d[O2_dot]/d[NH+]
    pd(11,98) =  &
        -k(766)*n(idx_O2)  &
        -k(767)*n(idx_O2)  &
        -k(180)*n(idx_O2)

    !d[CH2_dot]/d[NH+]
    pd(12,98) =  &
        -k(433)*n(idx_CH2)

    !d[H2CO_dot]/d[NH+]
    pd(13,98) =  &
        -k(176)*n(idx_H2CO)  &
        -k(754)*n(idx_H2CO)  &
        -k(753)*n(idx_H2CO)

    !d[HCO_dot]/d[NH+]
    pd(14,98) =  &
        -k(760)*n(idx_HCO)  &
        +k(751)*n(idx_CO2)

    !d[NH3_dot]/d[NH+]
    pd(16,98) =  &
        -k(178)*n(idx_NH3)

    !d[NO_dot]/d[NH+]
    pd(17,98) =  &
        -k(765)*n(idx_NO)  &
        -k(179)*n(idx_NO)

    !d[CN_dot]/d[NH+]
    pd(24,98) =  &
        -k(748)*n(idx_CN)

    !d[CO_dot]/d[NH+]
    pd(25,98) =  &
        -k(752)*n(idx_CO)  &
        +k(750)*n(idx_CO2)

    !d[N2_dot]/d[NH+]
    pd(26,98) =  &
        -k(762)*n(idx_N2)

    !d[NH2_dot]/d[NH+]
    pd(27,98) =  &
        +k(754)*n(idx_H2CO)  &
        -k(763)*n(idx_NH2)

    !d[N_dot]/d[NH+]
    pd(30,98) =  &
        +k(769)*n(idx_OH)  &
        +k(748)*n(idx_CN)  &
        +k(341)*n(idx_E)  &
        +k(1157)  &
        +k(433)*n(idx_CH2)  &
        +k(471)*n(idx_CH)  &
        +k(764)*n(idx_NH)  &
        +k(768)*n(idx_O)  &
        +k(541)*n(idx_H2)  &
        -k(741)*n(idx_N)  &
        +k(763)*n(idx_NH2)  &
        +k(760)*n(idx_HCO)  &
        +k(759)*n(idx_HCN)  &
        +k(749)*n(idx_CO2)  &
        +k(761)*n(idx_HNC)  &
        +k(755)*n(idx_H2O)  &
        +k(753)*n(idx_H2CO)  &
        +k(752)*n(idx_CO)  &
        +k(391)*n(idx_C)  &
        +k(762)*n(idx_N2)  &
        +k(767)*n(idx_O2)

    !d[NH_dot]/d[NH+]
    pd(31,98) =  &
        -k(764)*n(idx_NH)  &
        +k(176)*n(idx_H2CO)  &
        +k(177)*n(idx_H2O)  &
        +k(180)*n(idx_O2)  &
        +k(178)*n(idx_NH3)  &
        +k(179)*n(idx_NO)

    !d[CO2_dot]/d[NH+]
    pd(38,98) =  &
        -k(750)*n(idx_CO2)  &
        -k(751)*n(idx_CO2)  &
        -k(749)*n(idx_CO2)

    !d[NH3_DUST_dot]/d[NH+]
    pd(60,98) =  &
        +k(1274)

    !d[HCO+_dot]/d[NH+]
    pd(70,98) =  &
        +k(752)*n(idx_CO)  &
        +k(754)*n(idx_H2CO)

    !d[H+_dot]/d[NH+]
    pd(71,98) =  &
        +k(1157)

    !d[CH2+_dot]/d[NH+]
    pd(74,98) =  &
        +k(471)*n(idx_CH)

    !d[CH+_dot]/d[NH+]
    pd(75,98) =  &
        +k(391)*n(idx_C)

    !d[H2CO+_dot]/d[NH+]
    pd(76,98) =  &
        +k(760)*n(idx_HCO)  &
        +k(176)*n(idx_H2CO)

    !d[NH3+_dot]/d[NH+]
    pd(78,98) =  &
        +k(763)*n(idx_NH2)  &
        +k(178)*n(idx_NH3)  &
        +k(757)*n(idx_H2O)

    !d[NO+_dot]/d[NH+]
    pd(79,98) =  &
        +k(766)*n(idx_O2)  &
        +k(751)*n(idx_CO2)  &
        +k(179)*n(idx_NO)

    !d[N2+_dot]/d[NH+]
    pd(88,98) =  &
        +k(741)*n(idx_N)

    !d[O2+_dot]/d[NH+]
    pd(89,98) =  &
        +k(180)*n(idx_O2)

    !d[H2O+_dot]/d[NH+]
    pd(90,98) =  &
        +k(769)*n(idx_OH)  &
        +k(177)*n(idx_H2O)

    !d[NH2+_dot]/d[NH+]
    pd(91,98) =  &
        +k(542)*n(idx_H2)  &
        +k(758)*n(idx_H2O)  &
        +k(764)*n(idx_NH)

    !d[OH+_dot]/d[NH+]
    pd(93,98) =  &
        +k(768)*n(idx_O)

    !d[CH3+_dot]/d[NH+]
    pd(94,98) =  &
        +k(433)*n(idx_CH2)

    !d[HCN+_dot]/d[NH+]
    pd(97,98) =  &
        +k(748)*n(idx_CN)

    !d[NH+_dot]/d[NH+]
    pd(98,98) =  &
        -k(757)*n(idx_H2O)  &
        -k(433)*n(idx_CH2)  &
        -k(542)*n(idx_H2)  &
        -k(752)*n(idx_CO)  &
        -k(755)*n(idx_H2O)  &
        -k(769)*n(idx_OH)  &
        -k(766)*n(idx_O2)  &
        -k(759)*n(idx_HCN)  &
        -k(760)*n(idx_HCO)  &
        -k(753)*n(idx_H2CO)  &
        -k(748)*n(idx_CN)  &
        -k(765)*n(idx_NO)  &
        -k(176)*n(idx_H2CO)  &
        -k(391)*n(idx_C)  &
        -k(754)*n(idx_H2CO)  &
        -k(762)*n(idx_N2)  &
        -k(180)*n(idx_O2)  &
        -k(750)*n(idx_CO2)  &
        -k(764)*n(idx_NH)  &
        -k(756)*n(idx_H2O)  &
        -k(178)*n(idx_NH3)  &
        -k(1274)  &
        -k(749)*n(idx_CO2)  &
        -k(758)*n(idx_H2O)  &
        -k(751)*n(idx_CO2)  &
        -k(471)*n(idx_CH)  &
        -k(767)*n(idx_O2)  &
        -k(341)*n(idx_E)  &
        -k(177)*n(idx_H2O)  &
        -k(179)*n(idx_NO)  &
        -k(741)*n(idx_N)  &
        -k(541)*n(idx_H2)  &
        -k(768)*n(idx_O)  &
        -k(761)*n(idx_HNC)  &
        -k(1157)  &
        -k(763)*n(idx_NH2)

    !d[HNO+_dot]/d[NH+]
    pd(104,98) =  &
        +k(756)*n(idx_H2O)  &
        +k(750)*n(idx_CO2)

    !d[H3+_dot]/d[NH+]
    pd(106,98) =  &
        +k(541)*n(idx_H2)

    !d[H3CO+_dot]/d[NH+]
    pd(107,98) =  &
        +k(753)*n(idx_H2CO)

    !d[H3O+_dot]/d[NH+]
    pd(108,98) =  &
        +k(755)*n(idx_H2O)

    !d[HCNH+_dot]/d[NH+]
    pd(109,98) =  &
        +k(759)*n(idx_HCN)  &
        +k(761)*n(idx_HNC)

    !d[HCO2+_dot]/d[NH+]
    pd(110,98) =  &
        +k(749)*n(idx_CO2)

    !d[N2H+_dot]/d[NH+]
    pd(112,98) =  &
        +k(765)*n(idx_NO)  &
        +k(762)*n(idx_N2)

    !d[O2H+_dot]/d[NH+]
    pd(113,98) =  &
        +k(767)*n(idx_O2)

    !d[E_dot]/d[SIH4+]
    pd(1,99) =  &
        -k(359)*n(idx_E)  &
        -k(360)*n(idx_E)

    !d[H2_dot]/d[SIH4+]
    pd(6,99) =  &
        -k(547)*n(idx_H2)  &
        +k(359)*n(idx_E)

    !d[H_dot]/d[SIH4+]
    pd(8,99) =  &
        +k(360)*n(idx_E)  &
        +k(547)*n(idx_H2)

    !d[H2O_dot]/d[SIH4+]
    pd(9,99) =  &
        -k(575)*n(idx_H2O)

    !d[SIH2_dot]/d[SIH4+]
    pd(22,99) =  &
        +k(359)*n(idx_E)

    !d[SIH3_dot]/d[SIH4+]
    pd(23,99) =  &
        +k(490)*n(idx_CO)  &
        +k(575)*n(idx_H2O)  &
        +k(360)*n(idx_E)

    !d[CO_dot]/d[SIH4+]
    pd(25,99) =  &
        -k(490)*n(idx_CO)

    !d[SIH4_DUST_dot]/d[SIH4+]
    pd(48,99) =  &
        +k(1239)

    !d[HCO+_dot]/d[SIH4+]
    pd(70,99) =  &
        +k(490)*n(idx_CO)

    !d[SIH4+_dot]/d[SIH4+]
    pd(99,99) =  &
        -k(359)*n(idx_E)  &
        -k(547)*n(idx_H2)  &
        -k(1239)  &
        -k(360)*n(idx_E)  &
        -k(575)*n(idx_H2O)  &
        -k(490)*n(idx_CO)

    !d[H3O+_dot]/d[SIH4+]
    pd(108,99) =  &
        +k(575)*n(idx_H2O)

    !d[SIH5+_dot]/d[SIH4+]
    pd(114,99) =  &
        +k(547)*n(idx_H2)

    !d[E_dot]/d[SIH+]
    pd(1,100) =  &
        -k(353)*n(idx_E)

    !d[CH_dot]/d[SIH+]
    pd(2,100) =  &
        -k(478)*n(idx_CH)

    !d[O_dot]/d[SIH+]
    pd(3,100) =  &
        -k(833)*n(idx_O)

    !d[H2_dot]/d[SIH+]
    pd(6,100) =  &
        -k(1204)*n(idx_H2)  &
        +k(620)*n(idx_H)

    !d[C_dot]/d[SIH+]
    pd(7,100) =  &
        -k(395)*n(idx_C)

    !d[H_dot]/d[SIH+]
    pd(8,100) =  &
        +k(353)*n(idx_E)  &
        +k(395)*n(idx_C)  &
        -k(620)*n(idx_H)  &
        +k(833)*n(idx_O)  &
        +k(1180)

    !d[H2O_dot]/d[SIH+]
    pd(9,100) =  &
        -k(574)*n(idx_H2O)

    !d[SI_dot]/d[SIH+]
    pd(18,100) =  &
        +k(353)*n(idx_E)  &
        +k(574)*n(idx_H2O)  &
        +k(478)*n(idx_CH)

    !d[SIH4_DUST_dot]/d[SIH+]
    pd(48,100) =  &
        +k(1233)

    !d[CH2+_dot]/d[SIH+]
    pd(74,100) =  &
        +k(478)*n(idx_CH)

    !d[SI+_dot]/d[SIH+]
    pd(80,100) =  &
        +k(1180)  &
        +k(620)*n(idx_H)

    !d[SIC+_dot]/d[SIH+]
    pd(83,100) =  &
        +k(395)*n(idx_C)

    !d[SIH3+_dot]/d[SIH+]
    pd(85,100) =  &
        +k(1204)*n(idx_H2)

    !d[SIH+_dot]/d[SIH+]
    pd(100,100) =  &
        -k(395)*n(idx_C)  &
        -k(1204)*n(idx_H2)  &
        -k(574)*n(idx_H2O)  &
        -k(1233)  &
        -k(1180)  &
        -k(478)*n(idx_CH)  &
        -k(620)*n(idx_H)  &
        -k(833)*n(idx_O)  &
        -k(353)*n(idx_E)

    !d[SIO+_dot]/d[SIH+]
    pd(101,100) =  &
        +k(833)*n(idx_O)

    !d[H3O+_dot]/d[SIH+]
    pd(108,100) =  &
        +k(574)*n(idx_H2O)

    !d[E_dot]/d[SIO+]
    pd(1,101) =  &
        -k(363)*n(idx_E)

    !d[CH_dot]/d[SIO+]
    pd(2,101) =  &
        -k(479)*n(idx_CH)

    !d[O_dot]/d[SIO+]
    pd(3,101) =  &
        +k(1190)  &
        -k(836)*n(idx_O)  &
        +k(363)*n(idx_E)

    !d[H2_dot]/d[SIO+]
    pd(6,101) =  &
        -k(548)*n(idx_H2)

    !d[C_dot]/d[SIO+]
    pd(7,101) =  &
        -k(396)*n(idx_C)

    !d[H_dot]/d[SIO+]
    pd(8,101) =  &
        +k(548)*n(idx_H2)

    !d[O2_dot]/d[SIO+]
    pd(11,101) =  &
        +k(836)*n(idx_O)

    !d[CH2_dot]/d[SIO+]
    pd(12,101) =  &
        -k(439)*n(idx_CH2)

    !d[H2CO_dot]/d[SIO+]
    pd(13,101) =  &
        +k(439)*n(idx_CH2)

    !d[HCO_dot]/d[SIO+]
    pd(14,101) =  &
        -k(139)*n(idx_HCO)

    !d[MG_dot]/d[SIO+]
    pd(15,101) =  &
        -k(155)*n(idx_MG)

    !d[NO_dot]/d[SIO+]
    pd(17,101) =  &
        -k(207)*n(idx_NO)  &
        +k(747)*n(idx_N)

    !d[SI_dot]/d[SIO+]
    pd(18,101) =  &
        +k(479)*n(idx_CH)  &
        +k(363)*n(idx_E)  &
        +k(746)*n(idx_N)

    !d[CO_dot]/d[SIO+]
    pd(25,101) =  &
        -k(491)*n(idx_CO)  &
        +k(396)*n(idx_C)

    !d[N_dot]/d[SIO+]
    pd(30,101) =  &
        -k(747)*n(idx_N)  &
        -k(746)*n(idx_N)

    !d[SIO_dot]/d[SIO+]
    pd(34,101) =  &
        +k(155)*n(idx_MG)  &
        +k(207)*n(idx_NO)  &
        +k(139)*n(idx_HCO)

    !d[CO2_dot]/d[SIO+]
    pd(38,101) =  &
        +k(491)*n(idx_CO)

    !d[H2SIO_DUST_dot]/d[SIO+]
    pd(49,101) =  &
        +k(1246)

    !d[HCO+_dot]/d[SIO+]
    pd(70,101) =  &
        +k(479)*n(idx_CH)  &
        +k(139)*n(idx_HCO)

    !d[MG+_dot]/d[SIO+]
    pd(77,101) =  &
        +k(155)*n(idx_MG)

    !d[NO+_dot]/d[SIO+]
    pd(79,101) =  &
        +k(207)*n(idx_NO)  &
        +k(746)*n(idx_N)

    !d[SI+_dot]/d[SIO+]
    pd(80,101) =  &
        +k(747)*n(idx_N)  &
        +k(439)*n(idx_CH2)  &
        +k(836)*n(idx_O)  &
        +k(491)*n(idx_CO)  &
        +k(1190)  &
        +k(396)*n(idx_C)

    !d[SIO+_dot]/d[SIO+]
    pd(101,101) =  &
        -k(747)*n(idx_N)  &
        -k(491)*n(idx_CO)  &
        -k(139)*n(idx_HCO)  &
        -k(746)*n(idx_N)  &
        -k(155)*n(idx_MG)  &
        -k(1246)  &
        -k(479)*n(idx_CH)  &
        -k(836)*n(idx_O)  &
        -k(207)*n(idx_NO)  &
        -k(363)*n(idx_E)  &
        -k(439)*n(idx_CH2)  &
        -k(1190)  &
        -k(548)*n(idx_H2)  &
        -k(396)*n(idx_C)

    !d[SIOH+_dot]/d[SIO+]
    pd(115,101) =  &
        +k(548)*n(idx_H2)

    !d[E_dot]/d[H2+]
    pd(1,102) =  &
        -k(306)*n(idx_E)

    !d[CH_dot]/d[H2+]
    pd(2,102) =  &
        -k(103)*n(idx_CH)  &
        -k(513)*n(idx_CH)

    !d[O_dot]/d[H2+]
    pd(3,102) =  &
        -k(527)*n(idx_O)

    !d[HCN_dot]/d[H2+]
    pd(5,102) =  &
        -k(108)*n(idx_HCN)

    !d[H2_dot]/d[H2+]
    pd(6,102) =  &
        +k(102)*n(idx_CH4)  &
        -k(517)*n(idx_H2)  &
        +k(103)*n(idx_CH)  &
        +k(112)*n(idx_NH)  &
        +k(114)*n(idx_O2)  &
        +k(111)*n(idx_NH3)  &
        +k(113)*n(idx_NO)  &
        +k(108)*n(idx_HCN)  &
        +k(110)*n(idx_NH2)  &
        +k(101)*n(idx_CH2)  &
        +k(106)*n(idx_H2CO)  &
        +k(115)*n(idx_OH)  &
        +k(512)*n(idx_CH4)  &
        +k(104)*n(idx_CN)  &
        +k(518)*n(idx_H2CO)  &
        +k(109)*n(idx_HCO)  &
        +k(105)*n(idx_CO)  &
        +k(129)*n(idx_H)  &
        +k(107)*n(idx_H2O)

    !d[C_dot]/d[H2+]
    pd(7,102) =  &
        -k(510)*n(idx_C)

    !d[H_dot]/d[H2+]
    pd(8,102) =  &
        -k(129)*n(idx_H)  &
        +k(516)*n(idx_CO)  &
        +k(511)*n(idx_CH2)  &
        +k(528)*n(idx_OH)  &
        +k(524)*n(idx_NH)  &
        +k(527)*n(idx_O)  &
        +k(1135)  &
        +k(510)*n(idx_C)  &
        +k(521)*n(idx_HE)  &
        +k(523)*n(idx_N)  &
        +k(526)*n(idx_O2)  &
        +k(519)*n(idx_H2O)  &
        +k(522)*n(idx_N2)  &
        +k(515)*n(idx_CO2)  &
        +k(517)*n(idx_H2)  &
        +k(512)*n(idx_CH4)  &
        +2.d0*k(306)*n(idx_E)  &
        +k(518)*n(idx_H2CO)  &
        +k(514)*n(idx_CN)  &
        +k(525)*n(idx_NO)  &
        +k(513)*n(idx_CH)

    !d[H2O_dot]/d[H2+]
    pd(9,102) =  &
        -k(519)*n(idx_H2O)  &
        -k(107)*n(idx_H2O)

    !d[OH_dot]/d[H2+]
    pd(10,102) =  &
        -k(115)*n(idx_OH)  &
        -k(528)*n(idx_OH)

    !d[O2_dot]/d[H2+]
    pd(11,102) =  &
        -k(114)*n(idx_O2)  &
        -k(526)*n(idx_O2)

    !d[CH2_dot]/d[H2+]
    pd(12,102) =  &
        -k(101)*n(idx_CH2)  &
        -k(511)*n(idx_CH2)

    !d[H2CO_dot]/d[H2+]
    pd(13,102) =  &
        -k(518)*n(idx_H2CO)  &
        -k(106)*n(idx_H2CO)

    !d[HCO_dot]/d[H2+]
    pd(14,102) =  &
        -k(520)*n(idx_HCO)  &
        -k(109)*n(idx_HCO)

    !d[NH3_dot]/d[H2+]
    pd(16,102) =  &
        -k(111)*n(idx_NH3)

    !d[NO_dot]/d[H2+]
    pd(17,102) =  &
        -k(525)*n(idx_NO)  &
        -k(113)*n(idx_NO)

    !d[CN_dot]/d[H2+]
    pd(24,102) =  &
        -k(104)*n(idx_CN)  &
        -k(514)*n(idx_CN)

    !d[CO_dot]/d[H2+]
    pd(25,102) =  &
        +k(520)*n(idx_HCO)  &
        -k(516)*n(idx_CO)  &
        -k(105)*n(idx_CO)

    !d[N2_dot]/d[H2+]
    pd(26,102) =  &
        -k(522)*n(idx_N2)

    !d[NH2_dot]/d[H2+]
    pd(27,102) =  &
        -k(110)*n(idx_NH2)

    !d[CH4_dot]/d[H2+]
    pd(29,102) =  &
        -k(512)*n(idx_CH4)  &
        -k(102)*n(idx_CH4)

    !d[N_dot]/d[H2+]
    pd(30,102) =  &
        -k(523)*n(idx_N)

    !d[NH_dot]/d[H2+]
    pd(31,102) =  &
        -k(524)*n(idx_NH)  &
        -k(112)*n(idx_NH)

    !d[HE_dot]/d[H2+]
    pd(35,102) =  &
        -k(521)*n(idx_HE)

    !d[CO2_dot]/d[H2+]
    pd(38,102) =  &
        -k(515)*n(idx_CO2)

    !d[HCO+_dot]/d[H2+]
    pd(70,102) =  &
        +k(516)*n(idx_CO)  &
        +k(109)*n(idx_HCO)  &
        +k(518)*n(idx_H2CO)

    !d[H+_dot]/d[H2+]
    pd(71,102) =  &
        +k(1135)  &
        +k(129)*n(idx_H)

    !d[CH2+_dot]/d[H2+]
    pd(74,102) =  &
        +k(101)*n(idx_CH2)  &
        +k(513)*n(idx_CH)

    !d[CH+_dot]/d[H2+]
    pd(75,102) =  &
        +k(103)*n(idx_CH)  &
        +k(510)*n(idx_C)

    !d[H2CO+_dot]/d[H2+]
    pd(76,102) =  &
        +k(106)*n(idx_H2CO)

    !d[NH3+_dot]/d[H2+]
    pd(78,102) =  &
        +k(111)*n(idx_NH3)

    !d[NO+_dot]/d[H2+]
    pd(79,102) =  &
        +k(113)*n(idx_NO)

    !d[CN+_dot]/d[H2+]
    pd(86,102) =  &
        +k(104)*n(idx_CN)

    !d[CO+_dot]/d[H2+]
    pd(87,102) =  &
        +k(105)*n(idx_CO)

    !d[O2+_dot]/d[H2+]
    pd(89,102) =  &
        +k(114)*n(idx_O2)

    !d[H2O+_dot]/d[H2+]
    pd(90,102) =  &
        +k(528)*n(idx_OH)  &
        +k(107)*n(idx_H2O)

    !d[NH2+_dot]/d[H2+]
    pd(91,102) =  &
        +k(110)*n(idx_NH2)  &
        +k(524)*n(idx_NH)

    !d[OH+_dot]/d[H2+]
    pd(93,102) =  &
        +k(527)*n(idx_O)  &
        +k(115)*n(idx_OH)

    !d[CH3+_dot]/d[H2+]
    pd(94,102) =  &
        +k(511)*n(idx_CH2)  &
        +k(512)*n(idx_CH4)

    !d[CH4+_dot]/d[H2+]
    pd(95,102) =  &
        +k(102)*n(idx_CH4)

    !d[HCN+_dot]/d[H2+]
    pd(97,102) =  &
        +k(514)*n(idx_CN)  &
        +k(108)*n(idx_HCN)

    !d[NH+_dot]/d[H2+]
    pd(98,102) =  &
        +k(112)*n(idx_NH)  &
        +k(523)*n(idx_N)

    !d[H2+_dot]/d[H2+]
    pd(102,102) =  &
        -k(523)*n(idx_N)  &
        -k(515)*n(idx_CO2)  &
        -k(522)*n(idx_N2)  &
        -k(516)*n(idx_CO)  &
        -k(306)*n(idx_E)  &
        -k(101)*n(idx_CH2)  &
        -k(102)*n(idx_CH4)  &
        -k(112)*n(idx_NH)  &
        -k(113)*n(idx_NO)  &
        -k(528)*n(idx_OH)  &
        -k(521)*n(idx_HE)  &
        -k(518)*n(idx_H2CO)  &
        -k(111)*n(idx_NH3)  &
        -k(510)*n(idx_C)  &
        -k(107)*n(idx_H2O)  &
        -k(526)*n(idx_O2)  &
        -k(527)*n(idx_O)  &
        -k(115)*n(idx_OH)  &
        -k(108)*n(idx_HCN)  &
        -k(520)*n(idx_HCO)  &
        -k(511)*n(idx_CH2)  &
        -k(129)*n(idx_H)  &
        -k(512)*n(idx_CH4)  &
        -k(104)*n(idx_CN)  &
        -k(525)*n(idx_NO)  &
        -k(103)*n(idx_CH)  &
        -k(114)*n(idx_O2)  &
        -k(514)*n(idx_CN)  &
        -k(513)*n(idx_CH)  &
        -k(524)*n(idx_NH)  &
        -k(109)*n(idx_HCO)  &
        -k(105)*n(idx_CO)  &
        -k(517)*n(idx_H2)  &
        -k(106)*n(idx_H2CO)  &
        -k(110)*n(idx_NH2)  &
        -k(1135)  &
        -k(519)*n(idx_H2O)

    !d[HNO+_dot]/d[H2+]
    pd(104,102) =  &
        +k(525)*n(idx_NO)

    !d[H3+_dot]/d[H2+]
    pd(106,102) =  &
        +k(520)*n(idx_HCO)  &
        +k(517)*n(idx_H2)

    !d[H3O+_dot]/d[H2+]
    pd(108,102) =  &
        +k(519)*n(idx_H2O)

    !d[HCO2+_dot]/d[H2+]
    pd(110,102) =  &
        +k(515)*n(idx_CO2)

    !d[HEH+_dot]/d[H2+]
    pd(111,102) =  &
        +k(521)*n(idx_HE)

    !d[N2H+_dot]/d[H2+]
    pd(112,102) =  &
        +k(522)*n(idx_N2)

    !d[O2H+_dot]/d[H2+]
    pd(113,102) =  &
        +k(526)*n(idx_O2)

    !d[E_dot]/d[HE+]
    pd(1,103) =  &
        -k(1219)*n(idx_E)

    !d[CH_dot]/d[HE+]
    pd(2,103) =  &
        -k(663)*n(idx_CH)  &
        -k(142)*n(idx_CH)  &
        +k(678)*n(idx_HCN)

    !d[O_dot]/d[HE+]
    pd(3,103) =  &
        +k(666)*n(idx_CO2)  &
        +k(711)*n(idx_SIO)  &
        +k(673)*n(idx_H2CO)  &
        +k(696)*n(idx_NO)  &
        +k(670)*n(idx_CO)  &
        +k(698)*n(idx_OCN)  &
        +k(697)*n(idx_O2)  &
        +k(683)*n(idx_HCO)

    !d[HNC_dot]/d[HE+]
    pd(4,103) =  &
        -k(684)*n(idx_HNC)  &
        -k(685)*n(idx_HNC)  &
        -k(686)*n(idx_HNC)

    !d[HCN_dot]/d[HE+]
    pd(5,103) =  &
        -k(677)*n(idx_HCN)  &
        -k(680)*n(idx_HCN)  &
        -k(678)*n(idx_HCN)  &
        -k(679)*n(idx_HCN)

    !d[H2_dot]/d[HE+]
    pd(6,103) =  &
        +2.d0*k(708)*n(idx_SIH4)  &
        -k(116)*n(idx_H2)  &
        +k(654)*n(idx_CH2)  &
        +k(659)*n(idx_CH4)  &
        +k(656)*n(idx_CH3)  &
        +k(690)*n(idx_NH2)  &
        +k(706)*n(idx_SIH3)  &
        -k(537)*n(idx_H2)  &
        +k(671)*n(idx_H2CO)  &
        +k(692)*n(idx_NH3)  &
        +k(704)*n(idx_SIH2)  &
        +k(660)*n(idx_CH4)  &
        +k(709)*n(idx_SIH4)

    !d[C_dot]/d[HE+]
    pd(7,103) =  &
        +k(686)*n(idx_HNC)  &
        +k(668)*n(idx_CO2)  &
        +k(701)*n(idx_SIC3)  &
        +k(702)*n(idx_SIC)  &
        +k(664)*n(idx_CN)  &
        -k(140)*n(idx_C)

    !d[H_dot]/d[HE+]
    pd(8,103) =  &
        +k(700)*n(idx_OH)  &
        +k(677)*n(idx_HCN)  &
        +k(672)*n(idx_H2CO)  &
        +k(687)*n(idx_HNO)  &
        +k(679)*n(idx_HCN)  &
        +k(691)*n(idx_NH2)  &
        +k(537)*n(idx_H2)  &
        +k(659)*n(idx_CH4)  &
        +k(685)*n(idx_HNC)  &
        +k(710)*n(idx_SIH)  &
        +k(707)*n(idx_SIH3)  &
        +k(674)*n(idx_H2O)  &
        +k(681)*n(idx_HCO)  &
        +k(676)*n(idx_H2SIO)  &
        +k(693)*n(idx_NH3)  &
        +k(694)*n(idx_NH)  &
        +k(705)*n(idx_SIH2)  &
        +k(661)*n(idx_CH4)  &
        +k(684)*n(idx_HNC)  &
        +k(663)*n(idx_CH)  &
        +k(709)*n(idx_SIH4)  &
        +k(655)*n(idx_CH2)  &
        -k(131)*n(idx_H)

    !d[H2O_dot]/d[HE+]
    pd(9,103) =  &
        -k(675)*n(idx_H2O)  &
        -k(144)*n(idx_H2O)  &
        -k(674)*n(idx_H2O)

    !d[OH_dot]/d[HE+]
    pd(10,103) =  &
        +k(658)*n(idx_CH3OH)  &
        -k(700)*n(idx_OH)  &
        +k(675)*n(idx_H2O)

    !d[O2_dot]/d[HE+]
    pd(11,103) =  &
        -k(147)*n(idx_O2)  &
        -k(697)*n(idx_O2)  &
        +k(669)*n(idx_CO2)

    !d[CH2_dot]/d[HE+]
    pd(12,103) =  &
        -k(654)*n(idx_CH2)  &
        -k(655)*n(idx_CH2)

    !d[H2CO_dot]/d[HE+]
    pd(13,103) =  &
        -k(672)*n(idx_H2CO)  &
        -k(673)*n(idx_H2CO)  &
        -k(143)*n(idx_H2CO)  &
        -k(671)*n(idx_H2CO)

    !d[HCO_dot]/d[HE+]
    pd(14,103) =  &
        -k(681)*n(idx_HCO)  &
        -k(682)*n(idx_HCO)  &
        -k(683)*n(idx_HCO)

    !d[NH3_dot]/d[HE+]
    pd(16,103) =  &
        -k(146)*n(idx_NH3)  &
        -k(692)*n(idx_NH3)  &
        -k(693)*n(idx_NH3)

    !d[NO_dot]/d[HE+]
    pd(17,103) =  &
        -k(695)*n(idx_NO)  &
        +k(688)*n(idx_HNO)  &
        -k(696)*n(idx_NO)

    !d[SI_dot]/d[HE+]
    pd(18,103) =  &
        +k(712)*n(idx_SIO)  &
        -k(148)*n(idx_SI)  &
        +k(703)*n(idx_SIC)

    !d[SIC3_dot]/d[HE+]
    pd(20,103) =  &
        -k(701)*n(idx_SIC3)

    !d[SIC_dot]/d[HE+]
    pd(21,103) =  &
        -k(702)*n(idx_SIC)  &
        -k(703)*n(idx_SIC)

    !d[SIH2_dot]/d[HE+]
    pd(22,103) =  &
        -k(704)*n(idx_SIH2)  &
        -k(705)*n(idx_SIH2)

    !d[SIH3_dot]/d[HE+]
    pd(23,103) =  &
        -k(706)*n(idx_SIH3)  &
        -k(707)*n(idx_SIH3)

    !d[CN_dot]/d[HE+]
    pd(24,103) =  &
        -k(664)*n(idx_CN)  &
        +k(699)*n(idx_OCN)  &
        -k(665)*n(idx_CN)

    !d[CO_dot]/d[HE+]
    pd(25,103) =  &
        +k(682)*n(idx_HCO)  &
        -k(670)*n(idx_CO)  &
        +k(667)*n(idx_CO2)

    !d[N2_dot]/d[HE+]
    pd(26,103) =  &
        -k(689)*n(idx_N2)  &
        -k(145)*n(idx_N2)

    !d[NH2_dot]/d[HE+]
    pd(27,103) =  &
        -k(690)*n(idx_NH2)  &
        -k(691)*n(idx_NH2)

    !d[CH3_dot]/d[HE+]
    pd(28,103) =  &
        +k(662)*n(idx_CH4)  &
        +k(657)*n(idx_CH3OH)  &
        -k(656)*n(idx_CH3)

    !d[CH4_dot]/d[HE+]
    pd(29,103) =  &
        -k(661)*n(idx_CH4)  &
        -k(660)*n(idx_CH4)  &
        -k(662)*n(idx_CH4)  &
        -k(141)*n(idx_CH4)  &
        -k(659)*n(idx_CH4)

    !d[N_dot]/d[HE+]
    pd(30,103) =  &
        +k(689)*n(idx_N2)  &
        +k(695)*n(idx_NO)  &
        +k(679)*n(idx_HCN)  &
        +k(680)*n(idx_HCN)  &
        +k(665)*n(idx_CN)  &
        +k(685)*n(idx_HNC)

    !d[NH_dot]/d[HE+]
    pd(31,103) =  &
        -k(694)*n(idx_NH)

    !d[SIH4_dot]/d[HE+]
    pd(32,103) =  &
        -k(709)*n(idx_SIH4)  &
        -k(708)*n(idx_SIH4)

    !d[SIH_dot]/d[HE+]
    pd(33,103) =  &
        -k(710)*n(idx_SIH)

    !d[SIO_dot]/d[HE+]
    pd(34,103) =  &
        -k(712)*n(idx_SIO)  &
        -k(711)*n(idx_SIO)

    !d[HE_dot]/d[HE+]
    pd(35,103) =  &
        +k(710)*n(idx_SIH)  &
        +k(688)*n(idx_HNO)  &
        +k(706)*n(idx_SIH3)  &
        +k(686)*n(idx_HNC)  &
        +k(711)*n(idx_SIO)  &
        +k(694)*n(idx_NH)  &
        +k(662)*n(idx_CH4)  &
        +k(678)*n(idx_HCN)  &
        +k(131)*n(idx_H)  &
        +k(692)*n(idx_NH3)  &
        +k(680)*n(idx_HCN)  &
        +k(146)*n(idx_NH3)  &
        +k(677)*n(idx_HCN)  &
        +k(147)*n(idx_O2)  &
        +k(657)*n(idx_CH3OH)  &
        +k(672)*n(idx_H2CO)  &
        +k(687)*n(idx_HNO)  &
        +k(679)*n(idx_HCN)  &
        +k(691)*n(idx_NH2)  &
        +k(700)*n(idx_OH)  &
        +k(709)*n(idx_SIH4)  &
        +k(689)*n(idx_N2)  &
        +k(537)*n(idx_H2)  &
        +k(697)*n(idx_O2)  &
        +k(659)*n(idx_CH4)  &
        +k(656)*n(idx_CH3)  &
        +k(664)*n(idx_CN)  &
        +k(684)*n(idx_HNC)  &
        +k(683)*n(idx_HCO)  &
        +k(704)*n(idx_SIH2)  &
        +k(116)*n(idx_H2)  &
        +k(707)*n(idx_SIH3)  &
        +k(658)*n(idx_CH3OH)  &
        +k(698)*n(idx_OCN)  &
        +k(696)*n(idx_NO)  &
        +k(144)*n(idx_H2O)  &
        +k(142)*n(idx_CH)  &
        +k(674)*n(idx_H2O)  &
        +k(681)*n(idx_HCO)  &
        +k(695)*n(idx_NO)  &
        +k(676)*n(idx_H2SIO)  &
        +k(669)*n(idx_CO2)  &
        +k(701)*n(idx_SIC3)  &
        +k(660)*n(idx_CH4)  &
        +k(673)*n(idx_H2CO)  &
        +k(690)*n(idx_NH2)  &
        +k(143)*n(idx_H2CO)  &
        +k(666)*n(idx_CO2)  &
        +k(693)*n(idx_NH3)  &
        +k(670)*n(idx_CO)  &
        +k(654)*n(idx_CH2)  &
        +k(702)*n(idx_SIC)  &
        +k(148)*n(idx_SI)  &
        +k(665)*n(idx_CN)  &
        +k(685)*n(idx_HNC)  &
        +k(708)*n(idx_SIH4)  &
        +k(705)*n(idx_SIH2)  &
        +k(661)*n(idx_CH4)  &
        +k(699)*n(idx_OCN)  &
        +k(140)*n(idx_C)  &
        +k(668)*n(idx_CO2)  &
        +k(675)*n(idx_H2O)  &
        +k(703)*n(idx_SIC)  &
        +k(141)*n(idx_CH4)  &
        +k(712)*n(idx_SIO)  &
        +k(671)*n(idx_H2CO)  &
        +k(667)*n(idx_CO2)  &
        +k(1219)*n(idx_E)  &
        +k(663)*n(idx_CH)  &
        +k(145)*n(idx_N2)  &
        +k(655)*n(idx_CH2)

    !d[HNO_dot]/d[HE+]
    pd(36,103) =  &
        -k(688)*n(idx_HNO)  &
        -k(687)*n(idx_HNO)

    !d[CH3OH_dot]/d[HE+]
    pd(37,103) =  &
        -k(658)*n(idx_CH3OH)  &
        -k(657)*n(idx_CH3OH)

    !d[CO2_dot]/d[HE+]
    pd(38,103) =  &
        -k(668)*n(idx_CO2)  &
        -k(669)*n(idx_CO2)  &
        -k(666)*n(idx_CO2)  &
        -k(667)*n(idx_CO2)

    !d[H2SIO_dot]/d[HE+]
    pd(40,103) =  &
        -k(676)*n(idx_H2SIO)

    !d[OCN_dot]/d[HE+]
    pd(44,103) =  &
        -k(699)*n(idx_OCN)  &
        -k(698)*n(idx_OCN)

    !d[HCO+_dot]/d[HE+]
    pd(70,103) =  &
        +k(672)*n(idx_H2CO)

    !d[H+_dot]/d[HE+]
    pd(71,103) =  &
        +k(662)*n(idx_CH4)  &
        +k(688)*n(idx_HNO)  &
        +k(537)*n(idx_H2)  &
        +k(131)*n(idx_H)  &
        +k(675)*n(idx_H2O)

    !d[C+_dot]/d[HE+]
    pd(73,103) =  &
        +k(654)*n(idx_CH2)  &
        +k(669)*n(idx_CO2)  &
        +k(703)*n(idx_SIC)  &
        +k(679)*n(idx_HCN)  &
        +k(663)*n(idx_CH)  &
        +k(670)*n(idx_CO)  &
        +k(140)*n(idx_C)  &
        +k(665)*n(idx_CN)  &
        +k(685)*n(idx_HNC)

    !d[CH2+_dot]/d[HE+]
    pd(74,103) =  &
        +k(673)*n(idx_H2CO)  &
        +k(660)*n(idx_CH4)

    !d[CH+_dot]/d[HE+]
    pd(75,103) =  &
        +k(142)*n(idx_CH)  &
        +k(659)*n(idx_CH4)  &
        +k(656)*n(idx_CH3)  &
        +k(680)*n(idx_HCN)  &
        +k(683)*n(idx_HCO)  &
        +k(655)*n(idx_CH2)

    !d[H2CO+_dot]/d[HE+]
    pd(76,103) =  &
        +k(143)*n(idx_H2CO)

    !d[NH3+_dot]/d[HE+]
    pd(78,103) =  &
        +k(146)*n(idx_NH3)

    !d[NO+_dot]/d[HE+]
    pd(79,103) =  &
        +k(687)*n(idx_HNO)

    !d[SI+_dot]/d[HE+]
    pd(80,103) =  &
        +k(708)*n(idx_SIH4)  &
        +k(711)*n(idx_SIO)  &
        +k(702)*n(idx_SIC)  &
        +k(704)*n(idx_SIH2)  &
        +k(148)*n(idx_SI)  &
        +k(710)*n(idx_SIH)

    !d[SIC2+_dot]/d[HE+]
    pd(81,103) =  &
        +k(701)*n(idx_SIC3)

    !d[SIH2+_dot]/d[HE+]
    pd(84,103) =  &
        +k(707)*n(idx_SIH3)

    !d[CN+_dot]/d[HE+]
    pd(86,103) =  &
        +k(698)*n(idx_OCN)  &
        +k(677)*n(idx_HCN)  &
        +k(684)*n(idx_HNC)

    !d[CO+_dot]/d[HE+]
    pd(87,103) =  &
        +k(681)*n(idx_HCO)  &
        +k(671)*n(idx_H2CO)  &
        +k(666)*n(idx_CO2)

    !d[N2+_dot]/d[HE+]
    pd(88,103) =  &
        +k(145)*n(idx_N2)

    !d[O2+_dot]/d[HE+]
    pd(89,103) =  &
        +k(147)*n(idx_O2)  &
        +k(668)*n(idx_CO2)

    !d[H2O+_dot]/d[HE+]
    pd(90,103) =  &
        +k(144)*n(idx_H2O)

    !d[NH2+_dot]/d[HE+]
    pd(91,103) =  &
        +k(693)*n(idx_NH3)

    !d[O+_dot]/d[HE+]
    pd(92,103) =  &
        +k(695)*n(idx_NO)  &
        +k(700)*n(idx_OH)  &
        +k(712)*n(idx_SIO)  &
        +k(667)*n(idx_CO2)  &
        +k(699)*n(idx_OCN)  &
        +k(697)*n(idx_O2)

    !d[OH+_dot]/d[HE+]
    pd(93,103) =  &
        +k(657)*n(idx_CH3OH)  &
        +k(674)*n(idx_H2O)

    !d[CH3+_dot]/d[HE+]
    pd(94,103) =  &
        +k(658)*n(idx_CH3OH)  &
        +k(661)*n(idx_CH4)

    !d[CH4+_dot]/d[HE+]
    pd(95,103) =  &
        +k(141)*n(idx_CH4)

    !d[N+_dot]/d[HE+]
    pd(96,103) =  &
        +k(689)*n(idx_N2)  &
        +k(690)*n(idx_NH2)  &
        +k(678)*n(idx_HCN)  &
        +k(664)*n(idx_CN)  &
        +k(696)*n(idx_NO)  &
        +k(694)*n(idx_NH)

    !d[NH+_dot]/d[HE+]
    pd(98,103) =  &
        +k(692)*n(idx_NH3)  &
        +k(686)*n(idx_HNC)  &
        +k(691)*n(idx_NH2)

    !d[SIH+_dot]/d[HE+]
    pd(100,103) =  &
        +k(705)*n(idx_SIH2)  &
        +k(706)*n(idx_SIH3)  &
        +k(709)*n(idx_SIH4)

    !d[H2+_dot]/d[HE+]
    pd(102,103) =  &
        +k(116)*n(idx_H2)

    !d[HE+_dot]/d[HE+]
    pd(103,103) =  &
        -k(658)*n(idx_CH3OH)  &
        -k(697)*n(idx_O2)  &
        -k(144)*n(idx_H2O)  &
        -k(669)*n(idx_CO2)  &
        -k(707)*n(idx_SIH3)  &
        -k(700)*n(idx_OH)  &
        -k(677)*n(idx_HCN)  &
        -k(116)*n(idx_H2)  &
        -k(685)*n(idx_HNC)  &
        -k(143)*n(idx_H2CO)  &
        -k(683)*n(idx_HCO)  &
        -k(710)*n(idx_SIH)  &
        -k(696)*n(idx_NO)  &
        -k(682)*n(idx_HCO)  &
        -k(142)*n(idx_CH)  &
        -k(673)*n(idx_H2CO)  &
        -k(140)*n(idx_C)  &
        -k(667)*n(idx_CO2)  &
        -k(661)*n(idx_CH4)  &
        -k(694)*n(idx_NH)  &
        -k(698)*n(idx_OCN)  &
        -k(703)*n(idx_SIC)  &
        -k(699)*n(idx_OCN)  &
        -k(665)*n(idx_CN)  &
        -k(692)*n(idx_NH3)  &
        -k(655)*n(idx_CH2)  &
        -k(695)*n(idx_NO)  &
        -k(705)*n(idx_SIH2)  &
        -k(672)*n(idx_H2CO)  &
        -k(688)*n(idx_HNO)  &
        -k(537)*n(idx_H2)  &
        -k(668)*n(idx_CO2)  &
        -k(679)*n(idx_HCN)  &
        -k(674)*n(idx_H2O)  &
        -k(662)*n(idx_CH4)  &
        -k(660)*n(idx_CH4)  &
        -k(701)*n(idx_SIC3)  &
        -k(664)*n(idx_CN)  &
        -k(681)*n(idx_HCO)  &
        -k(657)*n(idx_CH3OH)  &
        -k(687)*n(idx_HNO)  &
        -k(708)*n(idx_SIH4)  &
        -k(689)*n(idx_N2)  &
        -k(684)*n(idx_HNC)  &
        -k(147)*n(idx_O2)  &
        -k(686)*n(idx_HNC)  &
        -k(670)*n(idx_CO)  &
        -k(711)*n(idx_SIO)  &
        -k(680)*n(idx_HCN)  &
        -k(706)*n(idx_SIH3)  &
        -k(141)*n(idx_CH4)  &
        -k(146)*n(idx_NH3)  &
        -k(691)*n(idx_NH2)  &
        -k(148)*n(idx_SI)  &
        -k(145)*n(idx_N2)  &
        -k(678)*n(idx_HCN)  &
        -k(704)*n(idx_SIH2)  &
        -k(675)*n(idx_H2O)  &
        -k(671)*n(idx_H2CO)  &
        -k(690)*n(idx_NH2)  &
        -k(654)*n(idx_CH2)  &
        -k(663)*n(idx_CH)  &
        -k(676)*n(idx_H2SIO)  &
        -k(1219)*n(idx_E)  &
        -k(666)*n(idx_CO2)  &
        -k(656)*n(idx_CH3)  &
        -k(702)*n(idx_SIC)  &
        -k(693)*n(idx_NH3)  &
        -k(712)*n(idx_SIO)  &
        -k(659)*n(idx_CH4)  &
        -k(709)*n(idx_SIH4)  &
        -k(131)*n(idx_H)

    !d[HEH+_dot]/d[HE+]
    pd(111,103) =  &
        +k(682)*n(idx_HCO)

    !d[SIOH+_dot]/d[HE+]
    pd(115,103) =  &
        +k(676)*n(idx_H2SIO)

    !d[E_dot]/d[HNO+]
    pd(1,104) =  &
        -k(335)*n(idx_E)

    !d[CH_dot]/d[HNO+]
    pd(2,104) =  &
        -k(468)*n(idx_CH)

    !d[HNC_dot]/d[HNO+]
    pd(4,104) =  &
        -k(650)*n(idx_HNC)

    !d[HCN_dot]/d[HNO+]
    pd(5,104) =  &
        -k(631)*n(idx_HCN)

    !d[C_dot]/d[HNO+]
    pd(7,104) =  &
        -k(389)*n(idx_C)

    !d[H_dot]/d[HNO+]
    pd(8,104) =  &
        +k(335)*n(idx_E)

    !d[H2O_dot]/d[HNO+]
    pd(9,104) =  &
        -k(569)*n(idx_H2O)

    !d[OH_dot]/d[HNO+]
    pd(10,104) =  &
        -k(857)*n(idx_OH)

    !d[CH2_dot]/d[HNO+]
    pd(12,104) =  &
        -k(431)*n(idx_CH2)

    !d[H2CO_dot]/d[HNO+]
    pd(13,104) =  &
        -k(551)*n(idx_H2CO)

    !d[HCO_dot]/d[HNO+]
    pd(14,104) =  &
        -k(643)*n(idx_HCO)

    !d[NO_dot]/d[HNO+]
    pd(17,104) =  &
        +k(487)*n(idx_CO)  &
        +k(551)*n(idx_H2CO)  &
        +k(468)*n(idx_CH)  &
        +k(653)*n(idx_CO2)  &
        +k(643)*n(idx_HCO)  &
        +k(650)*n(idx_HNC)  &
        +k(389)*n(idx_C)  &
        +k(789)*n(idx_NH2)  &
        +k(483)*n(idx_CN)  &
        +k(631)*n(idx_HCN)  &
        +k(569)*n(idx_H2O)  &
        +k(335)*n(idx_E)  &
        +k(801)*n(idx_NH)  &
        +k(431)*n(idx_CH2)  &
        -k(205)*n(idx_NO)  &
        +k(857)*n(idx_OH)  &
        +k(733)*n(idx_N2)

    !d[CN_dot]/d[HNO+]
    pd(24,104) =  &
        -k(483)*n(idx_CN)

    !d[CO_dot]/d[HNO+]
    pd(25,104) =  &
        -k(487)*n(idx_CO)

    !d[N2_dot]/d[HNO+]
    pd(26,104) =  &
        -k(733)*n(idx_N2)

    !d[NH2_dot]/d[HNO+]
    pd(27,104) =  &
        -k(789)*n(idx_NH2)

    !d[NH_dot]/d[HNO+]
    pd(31,104) =  &
        -k(801)*n(idx_NH)

    !d[HNO_dot]/d[HNO+]
    pd(36,104) =  &
        +k(205)*n(idx_NO)

    !d[CO2_dot]/d[HNO+]
    pd(38,104) =  &
        -k(653)*n(idx_CO2)

    !d[HNO_DUST_dot]/d[HNO+]
    pd(63,104) =  &
        +k(1296)

    !d[HCO+_dot]/d[HNO+]
    pd(70,104) =  &
        +k(487)*n(idx_CO)

    !d[CH2+_dot]/d[HNO+]
    pd(74,104) =  &
        +k(468)*n(idx_CH)

    !d[CH+_dot]/d[HNO+]
    pd(75,104) =  &
        +k(389)*n(idx_C)

    !d[H2CO+_dot]/d[HNO+]
    pd(76,104) =  &
        +k(643)*n(idx_HCO)

    !d[NH3+_dot]/d[HNO+]
    pd(78,104) =  &
        +k(789)*n(idx_NH2)

    !d[NO+_dot]/d[HNO+]
    pd(79,104) =  &
        +k(205)*n(idx_NO)

    !d[H2O+_dot]/d[HNO+]
    pd(90,104) =  &
        +k(857)*n(idx_OH)

    !d[NH2+_dot]/d[HNO+]
    pd(91,104) =  &
        +k(801)*n(idx_NH)

    !d[CH3+_dot]/d[HNO+]
    pd(94,104) =  &
        +k(431)*n(idx_CH2)

    !d[HCN+_dot]/d[HNO+]
    pd(97,104) =  &
        +k(483)*n(idx_CN)

    !d[HNO+_dot]/d[HNO+]
    pd(104,104) =  &
        -k(389)*n(idx_C)  &
        -k(733)*n(idx_N2)  &
        -k(631)*n(idx_HCN)  &
        -k(468)*n(idx_CH)  &
        -k(1296)  &
        -k(487)*n(idx_CO)  &
        -k(569)*n(idx_H2O)  &
        -k(551)*n(idx_H2CO)  &
        -k(801)*n(idx_NH)  &
        -k(431)*n(idx_CH2)  &
        -k(643)*n(idx_HCO)  &
        -k(335)*n(idx_E)  &
        -k(789)*n(idx_NH2)  &
        -k(857)*n(idx_OH)  &
        -k(653)*n(idx_CO2)  &
        -k(205)*n(idx_NO)  &
        -k(483)*n(idx_CN)  &
        -k(650)*n(idx_HNC)

    !d[H3CO+_dot]/d[HNO+]
    pd(107,104) =  &
        +k(551)*n(idx_H2CO)

    !d[H3O+_dot]/d[HNO+]
    pd(108,104) =  &
        +k(569)*n(idx_H2O)

    !d[HCNH+_dot]/d[HNO+]
    pd(109,104) =  &
        +k(650)*n(idx_HNC)  &
        +k(631)*n(idx_HCN)

    !d[HCO2+_dot]/d[HNO+]
    pd(110,104) =  &
        +k(653)*n(idx_CO2)

    !d[N2H+_dot]/d[HNO+]
    pd(112,104) =  &
        +k(733)*n(idx_N2)

    !d[E_dot]/d[H2NO+]
    pd(1,105) =  &
        -k(312)*n(idx_E)  &
        -k(311)*n(idx_E)

    !d[H2_dot]/d[H2NO+]
    pd(6,105) =  &
        +k(312)*n(idx_E)

    !d[H_dot]/d[H2NO+]
    pd(8,105) =  &
        +k(311)*n(idx_E)  &
        +k(1297)

    !d[NO_dot]/d[H2NO+]
    pd(17,105) =  &
        +k(312)*n(idx_E)

    !d[HNO_dot]/d[H2NO+]
    pd(36,105) =  &
        +k(311)*n(idx_E)

    !d[HNO_DUST_dot]/d[H2NO+]
    pd(63,105) =  &
        +k(1297)

    !d[H2NO+_dot]/d[H2NO+]
    pd(105,105) =  &
        -k(312)*n(idx_E)  &
        -k(311)*n(idx_E)  &
        -k(1297)

    !d[E_dot]/d[H3+]
    pd(1,106) =  &
        -k(317)*n(idx_E)  &
        -k(316)*n(idx_E)

    !d[CH_dot]/d[H3+]
    pd(2,106) =  &
        -k(581)*n(idx_CH)

    !d[O_dot]/d[H3+]
    pd(3,106) =  &
        -k(599)*n(idx_O)  &
        -k(600)*n(idx_O)

    !d[HNC_dot]/d[H3+]
    pd(4,106) =  &
        -k(590)*n(idx_HNC)

    !d[HCN_dot]/d[H3+]
    pd(5,106) =  &
        -k(588)*n(idx_HCN)

    !d[H2_dot]/d[H3+]
    pd(6,106) =  &
        +k(581)*n(idx_CH)  &
        +k(593)*n(idx_N2)  &
        +k(606)*n(idx_SIH)  &
        +k(577)*n(idx_C)  &
        +k(598)*n(idx_O2)  &
        +k(591)*n(idx_HNO)  &
        +k(600)*n(idx_O)  &
        +k(585)*n(idx_CO)  &
        +k(595)*n(idx_NH)  &
        +k(596)*n(idx_NO2)  &
        +k(588)*n(idx_HCN)  &
        +k(592)*n(idx_MG)  &
        +k(586)*n(idx_H2CO)  &
        +k(579)*n(idx_CH3)  &
        +k(597)*n(idx_NO)  &
        +k(316)*n(idx_E)  &
        +k(584)*n(idx_CO)  &
        +k(578)*n(idx_CH2)  &
        +k(607)*n(idx_SIO)  &
        +k(605)*n(idx_SIH4)  &
        +k(587)*n(idx_H2O)  &
        +k(580)*n(idx_CH3OH)  &
        +k(1147)  &
        +k(590)*n(idx_HNC)  &
        +k(603)*n(idx_SIH2)  &
        +k(602)*n(idx_SI)  &
        +k(589)*n(idx_HCO)  &
        +k(604)*n(idx_SIH3)  &
        +k(601)*n(idx_OH)  &
        +k(583)*n(idx_CO2)  &
        +k(594)*n(idx_NH2)  &
        +k(582)*n(idx_CN)

    !d[C_dot]/d[H3+]
    pd(7,106) =  &
        -k(577)*n(idx_C)

    !d[H_dot]/d[H3+]
    pd(8,106) =  &
        +3.d0*k(317)*n(idx_E)  &
        +k(316)*n(idx_E)  &
        +k(1146)  &
        +k(599)*n(idx_O)  &
        +k(592)*n(idx_MG)

    !d[H2O_dot]/d[H3+]
    pd(9,106) =  &
        -k(587)*n(idx_H2O)  &
        +k(580)*n(idx_CH3OH)

    !d[OH_dot]/d[H3+]
    pd(10,106) =  &
        +k(596)*n(idx_NO2)  &
        -k(601)*n(idx_OH)

    !d[O2_dot]/d[H3+]
    pd(11,106) =  &
        -k(598)*n(idx_O2)

    !d[CH2_dot]/d[H3+]
    pd(12,106) =  &
        -k(578)*n(idx_CH2)

    !d[H2CO_dot]/d[H3+]
    pd(13,106) =  &
        -k(586)*n(idx_H2CO)

    !d[HCO_dot]/d[H3+]
    pd(14,106) =  &
        -k(589)*n(idx_HCO)

    !d[MG_dot]/d[H3+]
    pd(15,106) =  &
        -k(592)*n(idx_MG)

    !d[NO_dot]/d[H3+]
    pd(17,106) =  &
        -k(597)*n(idx_NO)

    !d[SI_dot]/d[H3+]
    pd(18,106) =  &
        -k(602)*n(idx_SI)

    !d[SIH2_dot]/d[H3+]
    pd(22,106) =  &
        -k(603)*n(idx_SIH2)

    !d[SIH3_dot]/d[H3+]
    pd(23,106) =  &
        -k(604)*n(idx_SIH3)

    !d[CN_dot]/d[H3+]
    pd(24,106) =  &
        -k(582)*n(idx_CN)

    !d[CO_dot]/d[H3+]
    pd(25,106) =  &
        -k(585)*n(idx_CO)  &
        -k(584)*n(idx_CO)

    !d[N2_dot]/d[H3+]
    pd(26,106) =  &
        -k(593)*n(idx_N2)

    !d[NH2_dot]/d[H3+]
    pd(27,106) =  &
        -k(594)*n(idx_NH2)

    !d[CH3_dot]/d[H3+]
    pd(28,106) =  &
        -k(579)*n(idx_CH3)

    !d[NH_dot]/d[H3+]
    pd(31,106) =  &
        -k(595)*n(idx_NH)

    !d[SIH4_dot]/d[H3+]
    pd(32,106) =  &
        -k(605)*n(idx_SIH4)

    !d[SIH_dot]/d[H3+]
    pd(33,106) =  &
        -k(606)*n(idx_SIH)

    !d[SIO_dot]/d[H3+]
    pd(34,106) =  &
        -k(607)*n(idx_SIO)

    !d[HNO_dot]/d[H3+]
    pd(36,106) =  &
        -k(591)*n(idx_HNO)

    !d[CH3OH_dot]/d[H3+]
    pd(37,106) =  &
        -k(580)*n(idx_CH3OH)

    !d[CO2_dot]/d[H3+]
    pd(38,106) =  &
        -k(583)*n(idx_CO2)

    !d[NO2_dot]/d[H3+]
    pd(42,106) =  &
        -k(596)*n(idx_NO2)

    !d[HCO+_dot]/d[H3+]
    pd(70,106) =  &
        +k(584)*n(idx_CO)

    !d[H+_dot]/d[H3+]
    pd(71,106) =  &
        +k(1147)

    !d[HOC+_dot]/d[H3+]
    pd(72,106) =  &
        +k(585)*n(idx_CO)

    !d[CH2+_dot]/d[H3+]
    pd(74,106) =  &
        +k(581)*n(idx_CH)

    !d[CH+_dot]/d[H3+]
    pd(75,106) =  &
        +k(577)*n(idx_C)

    !d[H2CO+_dot]/d[H3+]
    pd(76,106) =  &
        +k(589)*n(idx_HCO)

    !d[MG+_dot]/d[H3+]
    pd(77,106) =  &
        +k(592)*n(idx_MG)

    !d[NH3+_dot]/d[H3+]
    pd(78,106) =  &
        +k(594)*n(idx_NH2)

    !d[NO+_dot]/d[H3+]
    pd(79,106) =  &
        +k(596)*n(idx_NO2)

    !d[SIH2+_dot]/d[H3+]
    pd(84,106) =  &
        +k(606)*n(idx_SIH)

    !d[SIH3+_dot]/d[H3+]
    pd(85,106) =  &
        +k(603)*n(idx_SIH2)

    !d[H2O+_dot]/d[H3+]
    pd(90,106) =  &
        +k(601)*n(idx_OH)  &
        +k(599)*n(idx_O)

    !d[NH2+_dot]/d[H3+]
    pd(91,106) =  &
        +k(595)*n(idx_NH)

    !d[OH+_dot]/d[H3+]
    pd(93,106) =  &
        +k(600)*n(idx_O)

    !d[CH3+_dot]/d[H3+]
    pd(94,106) =  &
        +k(578)*n(idx_CH2)  &
        +k(580)*n(idx_CH3OH)

    !d[CH4+_dot]/d[H3+]
    pd(95,106) =  &
        +k(579)*n(idx_CH3)

    !d[HCN+_dot]/d[H3+]
    pd(97,106) =  &
        +k(582)*n(idx_CN)

    !d[SIH4+_dot]/d[H3+]
    pd(99,106) =  &
        +k(604)*n(idx_SIH3)

    !d[SIH+_dot]/d[H3+]
    pd(100,106) =  &
        +k(602)*n(idx_SI)

    !d[H2+_dot]/d[H3+]
    pd(102,106) =  &
        +k(1146)

    !d[HNO+_dot]/d[H3+]
    pd(104,106) =  &
        +k(597)*n(idx_NO)

    !d[H2NO+_dot]/d[H3+]
    pd(105,106) =  &
        +k(591)*n(idx_HNO)

    !d[H3+_dot]/d[H3+]
    pd(106,106) =  &
        -k(584)*n(idx_CO)  &
        -k(579)*n(idx_CH3)  &
        -k(1147)  &
        -k(592)*n(idx_MG)  &
        -k(589)*n(idx_HCO)  &
        -k(591)*n(idx_HNO)  &
        -k(593)*n(idx_N2)  &
        -k(582)*n(idx_CN)  &
        -k(599)*n(idx_O)  &
        -k(585)*n(idx_CO)  &
        -k(597)*n(idx_NO)  &
        -k(607)*n(idx_SIO)  &
        -k(581)*n(idx_CH)  &
        -k(577)*n(idx_C)  &
        -k(600)*n(idx_O)  &
        -k(605)*n(idx_SIH4)  &
        -k(1146)  &
        -k(317)*n(idx_E)  &
        -k(587)*n(idx_H2O)  &
        -k(594)*n(idx_NH2)  &
        -k(588)*n(idx_HCN)  &
        -k(595)*n(idx_NH)  &
        -k(316)*n(idx_E)  &
        -k(583)*n(idx_CO2)  &
        -k(578)*n(idx_CH2)  &
        -k(603)*n(idx_SIH2)  &
        -k(601)*n(idx_OH)  &
        -k(580)*n(idx_CH3OH)  &
        -k(602)*n(idx_SI)  &
        -k(596)*n(idx_NO2)  &
        -k(586)*n(idx_H2CO)  &
        -k(604)*n(idx_SIH3)  &
        -k(590)*n(idx_HNC)  &
        -k(598)*n(idx_O2)  &
        -k(606)*n(idx_SIH)

    !d[H3CO+_dot]/d[H3+]
    pd(107,106) =  &
        +k(586)*n(idx_H2CO)

    !d[H3O+_dot]/d[H3+]
    pd(108,106) =  &
        +k(587)*n(idx_H2O)

    !d[HCNH+_dot]/d[H3+]
    pd(109,106) =  &
        +k(588)*n(idx_HCN)  &
        +k(590)*n(idx_HNC)

    !d[HCO2+_dot]/d[H3+]
    pd(110,106) =  &
        +k(583)*n(idx_CO2)

    !d[N2H+_dot]/d[H3+]
    pd(112,106) =  &
        +k(593)*n(idx_N2)

    !d[O2H+_dot]/d[H3+]
    pd(113,106) =  &
        +k(598)*n(idx_O2)

    !d[SIH5+_dot]/d[H3+]
    pd(114,106) =  &
        +k(605)*n(idx_SIH4)

    !d[SIOH+_dot]/d[H3+]
    pd(115,106) =  &
        +k(607)*n(idx_SIO)

    !d[E_dot]/d[H3CO+]
    pd(1,107) =  &
        -k(320)*n(idx_E)  &
        -k(322)*n(idx_E)  &
        -k(319)*n(idx_E)  &
        -k(318)*n(idx_E)  &
        -k(321)*n(idx_E)

    !d[CH_dot]/d[H3CO+]
    pd(2,107) =  &
        +k(319)*n(idx_E)  &
        -k(462)*n(idx_CH)

    !d[HNC_dot]/d[H3CO+]
    pd(4,107) =  &
        -k(648)*n(idx_HNC)

    !d[HCN_dot]/d[H3CO+]
    pd(5,107) =  &
        -k(629)*n(idx_HCN)

    !d[H2_dot]/d[H3CO+]
    pd(6,107) =  &
        +k(320)*n(idx_E)

    !d[H_dot]/d[H3CO+]
    pd(8,107) =  &
        +2.d0*k(322)*n(idx_E)  &
        +k(320)*n(idx_E)  &
        +k(321)*n(idx_E)

    !d[H2O_dot]/d[H3CO+]
    pd(9,107) =  &
        +k(319)*n(idx_E)  &
        -k(565)*n(idx_H2O)

    !d[OH_dot]/d[H3CO+]
    pd(10,107) =  &
        +k(318)*n(idx_E)

    !d[CH2_dot]/d[H3CO+]
    pd(12,107) =  &
        +k(318)*n(idx_E)

    !d[H2CO_dot]/d[H3CO+]
    pd(13,107) =  &
        +k(783)*n(idx_NH2)  &
        +k(462)*n(idx_CH)  &
        +k(321)*n(idx_E)  &
        +k(648)*n(idx_HNC)  &
        +k(629)*n(idx_HCN)  &
        +k(565)*n(idx_H2O)

    !d[HCO_dot]/d[H3CO+]
    pd(14,107) =  &
        +k(322)*n(idx_E)

    !d[CO_dot]/d[H3CO+]
    pd(25,107) =  &
        +k(320)*n(idx_E)

    !d[NH2_dot]/d[H3CO+]
    pd(27,107) =  &
        -k(783)*n(idx_NH2)

    !d[CH3OH_DUST_dot]/d[H3CO+]
    pd(45,107) =  &
        +k(1250)

    !d[CH2+_dot]/d[H3CO+]
    pd(74,107) =  &
        +k(462)*n(idx_CH)

    !d[NH3+_dot]/d[H3CO+]
    pd(78,107) =  &
        +k(783)*n(idx_NH2)

    !d[H3CO+_dot]/d[H3CO+]
    pd(107,107) =  &
        -k(319)*n(idx_E)  &
        -k(318)*n(idx_E)  &
        -k(648)*n(idx_HNC)  &
        -k(565)*n(idx_H2O)  &
        -k(783)*n(idx_NH2)  &
        -k(322)*n(idx_E)  &
        -k(1250)  &
        -k(320)*n(idx_E)  &
        -k(321)*n(idx_E)  &
        -k(629)*n(idx_HCN)  &
        -k(462)*n(idx_CH)

    !d[H3O+_dot]/d[H3CO+]
    pd(108,107) =  &
        +k(565)*n(idx_H2O)

    !d[HCNH+_dot]/d[H3CO+]
    pd(109,107) =  &
        +k(648)*n(idx_HNC)  &
        +k(629)*n(idx_HCN)

    !d[E_dot]/d[H3O+]
    pd(1,108) =  &
        -k(324)*n(idx_E)  &
        -k(325)*n(idx_E)  &
        -k(323)*n(idx_E)  &
        -k(326)*n(idx_E)

    !d[CH_dot]/d[H3O+]
    pd(2,108) =  &
        -k(463)*n(idx_CH)

    !d[O_dot]/d[H3O+]
    pd(3,108) =  &
        +k(324)*n(idx_E)

    !d[HNC_dot]/d[H3O+]
    pd(4,108) =  &
        -k(610)*n(idx_HNC)

    !d[HCN_dot]/d[H3O+]
    pd(5,108) =  &
        -k(609)*n(idx_HCN)

    !d[H2_dot]/d[H3O+]
    pd(6,108) =  &
        +k(324)*n(idx_E)  &
        +k(325)*n(idx_E)  &
        +k(385)*n(idx_C)

    !d[C_dot]/d[H3O+]
    pd(7,108) =  &
        -k(385)*n(idx_C)

    !d[H_dot]/d[H3O+]
    pd(8,108) =  &
        +k(323)*n(idx_E)  &
        +k(324)*n(idx_E)  &
        +2.d0*k(326)*n(idx_E)  &
        +k(1287)

    !d[H2O_dot]/d[H3O+]
    pd(9,108) =  &
        +k(609)*n(idx_HCN)  &
        +k(611)*n(idx_SI)  &
        +k(610)*n(idx_HNC)  &
        +k(426)*n(idx_CH2)  &
        +k(608)*n(idx_H2CO)  &
        +k(784)*n(idx_NH2)  &
        +k(323)*n(idx_E)  &
        +k(612)*n(idx_SIH2)  &
        +k(613)*n(idx_SIH)  &
        +k(614)*n(idx_SIO)  &
        +k(463)*n(idx_CH)

    !d[OH_dot]/d[H3O+]
    pd(10,108) =  &
        +k(325)*n(idx_E)  &
        +k(326)*n(idx_E)

    !d[CH2_dot]/d[H3O+]
    pd(12,108) =  &
        -k(426)*n(idx_CH2)

    !d[H2CO_dot]/d[H3O+]
    pd(13,108) =  &
        -k(608)*n(idx_H2CO)

    !d[SI_dot]/d[H3O+]
    pd(18,108) =  &
        -k(611)*n(idx_SI)

    !d[SIH2_dot]/d[H3O+]
    pd(22,108) =  &
        -k(612)*n(idx_SIH2)

    !d[NH2_dot]/d[H3O+]
    pd(27,108) =  &
        -k(784)*n(idx_NH2)

    !d[SIH_dot]/d[H3O+]
    pd(33,108) =  &
        -k(613)*n(idx_SIH)

    !d[SIO_dot]/d[H3O+]
    pd(34,108) =  &
        -k(614)*n(idx_SIO)

    !d[H2O_DUST_dot]/d[H3O+]
    pd(55,108) =  &
        +k(1287)

    !d[HCO+_dot]/d[H3O+]
    pd(70,108) =  &
        +k(385)*n(idx_C)

    !d[CH2+_dot]/d[H3O+]
    pd(74,108) =  &
        +k(463)*n(idx_CH)

    !d[NH3+_dot]/d[H3O+]
    pd(78,108) =  &
        +k(784)*n(idx_NH2)

    !d[SIH2+_dot]/d[H3O+]
    pd(84,108) =  &
        +k(613)*n(idx_SIH)

    !d[SIH3+_dot]/d[H3O+]
    pd(85,108) =  &
        +k(612)*n(idx_SIH2)

    !d[CH3+_dot]/d[H3O+]
    pd(94,108) =  &
        +k(426)*n(idx_CH2)

    !d[SIH+_dot]/d[H3O+]
    pd(100,108) =  &
        +k(611)*n(idx_SI)

    !d[H3CO+_dot]/d[H3O+]
    pd(107,108) =  &
        +k(608)*n(idx_H2CO)

    !d[H3O+_dot]/d[H3O+]
    pd(108,108) =  &
        -k(612)*n(idx_SIH2)  &
        -k(784)*n(idx_NH2)  &
        -k(1287)  &
        -k(614)*n(idx_SIO)  &
        -k(326)*n(idx_E)  &
        -k(463)*n(idx_CH)  &
        -k(324)*n(idx_E)  &
        -k(325)*n(idx_E)  &
        -k(613)*n(idx_SIH)  &
        -k(608)*n(idx_H2CO)  &
        -k(609)*n(idx_HCN)  &
        -k(426)*n(idx_CH2)  &
        -k(610)*n(idx_HNC)  &
        -k(323)*n(idx_E)  &
        -k(385)*n(idx_C)  &
        -k(611)*n(idx_SI)

    !d[HCNH+_dot]/d[H3O+]
    pd(109,108) =  &
        +k(609)*n(idx_HCN)  &
        +k(610)*n(idx_HNC)

    !d[SIOH+_dot]/d[H3O+]
    pd(115,108) =  &
        +k(614)*n(idx_SIO)

    !d[E_dot]/d[HCNH+]
    pd(1,109) =  &
        -k(328)*n(idx_E)  &
        -k(329)*n(idx_E)  &
        -k(330)*n(idx_E)

    !d[CH_dot]/d[HCNH+]
    pd(2,109) =  &
        -k(466)*n(idx_CH)  &
        -k(465)*n(idx_CH)

    !d[HNC_dot]/d[HCNH+]
    pd(4,109) =  &
        +k(330)*n(idx_E)  &
        +k(787)*n(idx_NH2)  &
        +k(635)*n(idx_H2CO)  &
        +k(466)*n(idx_CH)  &
        +k(429)*n(idx_CH2)

    !d[HCN_dot]/d[HCNH+]
    pd(5,109) =  &
        +k(634)*n(idx_H2CO)  &
        +k(786)*n(idx_NH2)  &
        +k(428)*n(idx_CH2)  &
        +k(465)*n(idx_CH)  &
        +k(329)*n(idx_E)

    !d[H_dot]/d[HCNH+]
    pd(8,109) =  &
        +k(330)*n(idx_E)  &
        +k(1302)  &
        +2.d0*k(328)*n(idx_E)  &
        +k(329)*n(idx_E)

    !d[CH2_dot]/d[HCNH+]
    pd(12,109) =  &
        -k(429)*n(idx_CH2)  &
        -k(428)*n(idx_CH2)

    !d[H2CO_dot]/d[HCNH+]
    pd(13,109) =  &
        -k(634)*n(idx_H2CO)  &
        -k(635)*n(idx_H2CO)

    !d[CN_dot]/d[HCNH+]
    pd(24,109) =  &
        +k(328)*n(idx_E)

    !d[NH2_dot]/d[HCNH+]
    pd(27,109) =  &
        -k(786)*n(idx_NH2)  &
        -k(787)*n(idx_NH2)

    !d[HCN_DUST_dot]/d[HCNH+]
    pd(59,109) =  &
        +k(1302)

    !d[CH2+_dot]/d[HCNH+]
    pd(74,109) =  &
        +k(466)*n(idx_CH)  &
        +k(465)*n(idx_CH)

    !d[NH3+_dot]/d[HCNH+]
    pd(78,109) =  &
        +k(787)*n(idx_NH2)  &
        +k(786)*n(idx_NH2)

    !d[CH3+_dot]/d[HCNH+]
    pd(94,109) =  &
        +k(429)*n(idx_CH2)  &
        +k(428)*n(idx_CH2)

    !d[H3CO+_dot]/d[HCNH+]
    pd(107,109) =  &
        +k(634)*n(idx_H2CO)  &
        +k(635)*n(idx_H2CO)

    !d[HCNH+_dot]/d[HCNH+]
    pd(109,109) =  &
        -k(786)*n(idx_NH2)  &
        -k(466)*n(idx_CH)  &
        -k(1302)  &
        -k(635)*n(idx_H2CO)  &
        -k(330)*n(idx_E)  &
        -k(787)*n(idx_NH2)  &
        -k(428)*n(idx_CH2)  &
        -k(634)*n(idx_H2CO)  &
        -k(429)*n(idx_CH2)  &
        -k(465)*n(idx_CH)  &
        -k(328)*n(idx_E)  &
        -k(329)*n(idx_E)

    !d[E_dot]/d[HCO2+]
    pd(1,110) =  &
        -k(333)*n(idx_E)  &
        -k(332)*n(idx_E)  &
        -k(334)*n(idx_E)

    !d[O_dot]/d[HCO2+]
    pd(3,110) =  &
        -k(825)*n(idx_O)  &
        +k(333)*n(idx_E)

    !d[C_dot]/d[HCO2+]
    pd(7,110) =  &
        -k(388)*n(idx_C)

    !d[H_dot]/d[HCO2+]
    pd(8,110) =  &
        +k(1288)  &
        +k(333)*n(idx_E)  &
        +k(332)*n(idx_E)

    !d[H2O_dot]/d[HCO2+]
    pd(9,110) =  &
        -k(568)*n(idx_H2O)

    !d[OH_dot]/d[HCO2+]
    pd(10,110) =  &
        +k(334)*n(idx_E)

    !d[O2_dot]/d[HCO2+]
    pd(11,110) =  &
        +k(825)*n(idx_O)

    !d[CO_dot]/d[HCO2+]
    pd(25,110) =  &
        +k(334)*n(idx_E)  &
        +k(333)*n(idx_E)  &
        -k(486)*n(idx_CO)

    !d[CO2_dot]/d[HCO2+]
    pd(38,110) =  &
        +k(486)*n(idx_CO)  &
        +k(568)*n(idx_H2O)  &
        +k(388)*n(idx_C)  &
        +k(332)*n(idx_E)

    !d[CO2_DUST_dot]/d[HCO2+]
    pd(57,110) =  &
        +k(1288)

    !d[HCO+_dot]/d[HCO2+]
    pd(70,110) =  &
        +k(825)*n(idx_O)  &
        +k(486)*n(idx_CO)

    !d[CH+_dot]/d[HCO2+]
    pd(75,110) =  &
        +k(388)*n(idx_C)

    !d[H3O+_dot]/d[HCO2+]
    pd(108,110) =  &
        +k(568)*n(idx_H2O)

    !d[HCO2+_dot]/d[HCO2+]
    pd(110,110) =  &
        -k(825)*n(idx_O)  &
        -k(1288)  &
        -k(333)*n(idx_E)  &
        -k(332)*n(idx_E)  &
        -k(568)*n(idx_H2O)  &
        -k(334)*n(idx_E)  &
        -k(388)*n(idx_C)  &
        -k(486)*n(idx_CO)

    !d[E_dot]/d[HEH+]
    pd(1,111) =  &
        -k(337)*n(idx_E)

    !d[H2_dot]/d[HEH+]
    pd(6,111) =  &
        -k(538)*n(idx_H2)

    !d[H_dot]/d[HEH+]
    pd(8,111) =  &
        +k(337)*n(idx_E)  &
        -k(619)*n(idx_H)

    !d[HE_dot]/d[HEH+]
    pd(35,111) =  &
        +k(619)*n(idx_H)  &
        +k(337)*n(idx_E)  &
        +k(538)*n(idx_H2)

    !d[H2+_dot]/d[HEH+]
    pd(102,111) =  &
        +k(619)*n(idx_H)

    !d[H3+_dot]/d[HEH+]
    pd(106,111) =  &
        +k(538)*n(idx_H2)

    !d[HEH+_dot]/d[HEH+]
    pd(111,111) =  &
        -k(619)*n(idx_H)  &
        -k(337)*n(idx_E)  &
        -k(538)*n(idx_H2)

    !d[E_dot]/d[N2H+]
    pd(1,112) =  &
        -k(340)*n(idx_E)  &
        -k(339)*n(idx_E)

    !d[CH_dot]/d[N2H+]
    pd(2,112) =  &
        -k(470)*n(idx_CH)

    !d[O_dot]/d[N2H+]
    pd(3,112) =  &
        -k(827)*n(idx_O)

    !d[HNC_dot]/d[N2H+]
    pd(4,112) =  &
        -k(651)*n(idx_HNC)

    !d[HCN_dot]/d[N2H+]
    pd(5,112) =  &
        -k(632)*n(idx_HCN)

    !d[C_dot]/d[N2H+]
    pd(7,112) =  &
        -k(390)*n(idx_C)

    !d[H_dot]/d[N2H+]
    pd(8,112) =  &
        +k(339)*n(idx_E)  &
        +k(1303)

    !d[H2O_dot]/d[N2H+]
    pd(9,112) =  &
        -k(571)*n(idx_H2O)

    !d[OH_dot]/d[N2H+]
    pd(10,112) =  &
        -k(858)*n(idx_OH)

    !d[CH2_dot]/d[N2H+]
    pd(12,112) =  &
        -k(432)*n(idx_CH2)

    !d[H2CO_dot]/d[N2H+]
    pd(13,112) =  &
        -k(736)*n(idx_H2CO)

    !d[HCO_dot]/d[N2H+]
    pd(14,112) =  &
        -k(644)*n(idx_HCO)

    !d[CO_dot]/d[N2H+]
    pd(25,112) =  &
        -k(488)*n(idx_CO)

    !d[N2_dot]/d[N2H+]
    pd(26,112) =  &
        +k(390)*n(idx_C)  &
        +k(790)*n(idx_NH2)  &
        +k(571)*n(idx_H2O)  &
        +k(488)*n(idx_CO)  &
        +k(644)*n(idx_HCO)  &
        +k(858)*n(idx_OH)  &
        +k(632)*n(idx_HCN)  &
        +k(827)*n(idx_O)  &
        +k(735)*n(idx_CO2)  &
        +k(736)*n(idx_H2CO)  &
        +k(802)*n(idx_NH)  &
        +k(651)*n(idx_HNC)  &
        +k(339)*n(idx_E)  &
        +k(432)*n(idx_CH2)  &
        +k(470)*n(idx_CH)

    !d[NH2_dot]/d[N2H+]
    pd(27,112) =  &
        -k(790)*n(idx_NH2)

    !d[N_dot]/d[N2H+]
    pd(30,112) =  &
        +k(340)*n(idx_E)

    !d[NH_dot]/d[N2H+]
    pd(31,112) =  &
        +k(340)*n(idx_E)  &
        -k(802)*n(idx_NH)

    !d[CO2_dot]/d[N2H+]
    pd(38,112) =  &
        -k(735)*n(idx_CO2)

    !d[N2_DUST_dot]/d[N2H+]
    pd(58,112) =  &
        +k(1303)

    !d[HCO+_dot]/d[N2H+]
    pd(70,112) =  &
        +k(488)*n(idx_CO)

    !d[CH2+_dot]/d[N2H+]
    pd(74,112) =  &
        +k(470)*n(idx_CH)

    !d[CH+_dot]/d[N2H+]
    pd(75,112) =  &
        +k(390)*n(idx_C)

    !d[H2CO+_dot]/d[N2H+]
    pd(76,112) =  &
        +k(644)*n(idx_HCO)

    !d[NH3+_dot]/d[N2H+]
    pd(78,112) =  &
        +k(790)*n(idx_NH2)

    !d[H2O+_dot]/d[N2H+]
    pd(90,112) =  &
        +k(858)*n(idx_OH)

    !d[NH2+_dot]/d[N2H+]
    pd(91,112) =  &
        +k(802)*n(idx_NH)

    !d[OH+_dot]/d[N2H+]
    pd(93,112) =  &
        +k(827)*n(idx_O)

    !d[CH3+_dot]/d[N2H+]
    pd(94,112) =  &
        +k(432)*n(idx_CH2)

    !d[H3CO+_dot]/d[N2H+]
    pd(107,112) =  &
        +k(736)*n(idx_H2CO)

    !d[H3O+_dot]/d[N2H+]
    pd(108,112) =  &
        +k(571)*n(idx_H2O)

    !d[HCNH+_dot]/d[N2H+]
    pd(109,112) =  &
        +k(632)*n(idx_HCN)  &
        +k(651)*n(idx_HNC)

    !d[HCO2+_dot]/d[N2H+]
    pd(110,112) =  &
        +k(735)*n(idx_CO2)

    !d[N2H+_dot]/d[N2H+]
    pd(112,112) =  &
        -k(632)*n(idx_HCN)  &
        -k(858)*n(idx_OH)  &
        -k(827)*n(idx_O)  &
        -k(432)*n(idx_CH2)  &
        -k(1303)  &
        -k(571)*n(idx_H2O)  &
        -k(735)*n(idx_CO2)  &
        -k(644)*n(idx_HCO)  &
        -k(390)*n(idx_C)  &
        -k(470)*n(idx_CH)  &
        -k(651)*n(idx_HNC)  &
        -k(339)*n(idx_E)  &
        -k(488)*n(idx_CO)  &
        -k(790)*n(idx_NH2)  &
        -k(736)*n(idx_H2CO)  &
        -k(802)*n(idx_NH)  &
        -k(340)*n(idx_E)

    !d[E_dot]/d[O2H+]
    pd(1,113) =  &
        -k(348)*n(idx_E)

    !d[CH_dot]/d[O2H+]
    pd(2,113) =  &
        -k(475)*n(idx_CH)

    !d[O_dot]/d[O2H+]
    pd(3,113) =  &
        -k(830)*n(idx_O)

    !d[HNC_dot]/d[O2H+]
    pd(4,113) =  &
        -k(652)*n(idx_HNC)

    !d[HCN_dot]/d[O2H+]
    pd(5,113) =  &
        -k(633)*n(idx_HCN)

    !d[H2_dot]/d[O2H+]
    pd(6,113) =  &
        -k(545)*n(idx_H2)

    !d[C_dot]/d[O2H+]
    pd(7,113) =  &
        -k(393)*n(idx_C)

    !d[H_dot]/d[O2H+]
    pd(8,113) =  &
        +k(348)*n(idx_E)

    !d[H2O_dot]/d[O2H+]
    pd(9,113) =  &
        -k(572)*n(idx_H2O)

    !d[OH_dot]/d[O2H+]
    pd(10,113) =  &
        -k(859)*n(idx_OH)

    !d[O2_dot]/d[O2H+]
    pd(11,113) =  &
        +k(545)*n(idx_H2)  &
        +k(437)*n(idx_CH2)  &
        +k(475)*n(idx_CH)  &
        +k(822)*n(idx_CO2)  &
        +k(646)*n(idx_HCO)  &
        +k(553)*n(idx_H2CO)  &
        +k(348)*n(idx_E)  &
        +k(489)*n(idx_CO)  &
        +k(791)*n(idx_NH2)  &
        +k(830)*n(idx_O)  &
        +k(652)*n(idx_HNC)  &
        +k(572)*n(idx_H2O)  &
        +k(734)*n(idx_N2)  &
        +k(393)*n(idx_C)  &
        +k(633)*n(idx_HCN)  &
        +k(808)*n(idx_NO)  &
        +k(859)*n(idx_OH)  &
        +k(484)*n(idx_CN)  &
        +k(806)*n(idx_NH)

    !d[CH2_dot]/d[O2H+]
    pd(12,113) =  &
        -k(437)*n(idx_CH2)

    !d[H2CO_dot]/d[O2H+]
    pd(13,113) =  &
        -k(553)*n(idx_H2CO)

    !d[HCO_dot]/d[O2H+]
    pd(14,113) =  &
        -k(646)*n(idx_HCO)

    !d[NO_dot]/d[O2H+]
    pd(17,113) =  &
        -k(808)*n(idx_NO)

    !d[CN_dot]/d[O2H+]
    pd(24,113) =  &
        -k(484)*n(idx_CN)

    !d[CO_dot]/d[O2H+]
    pd(25,113) =  &
        -k(489)*n(idx_CO)

    !d[N2_dot]/d[O2H+]
    pd(26,113) =  &
        -k(734)*n(idx_N2)

    !d[NH2_dot]/d[O2H+]
    pd(27,113) =  &
        -k(791)*n(idx_NH2)

    !d[NH_dot]/d[O2H+]
    pd(31,113) =  &
        -k(806)*n(idx_NH)

    !d[CO2_dot]/d[O2H+]
    pd(38,113) =  &
        -k(822)*n(idx_CO2)

    !d[O2H_DUST_dot]/d[O2H+]
    pd(64,113) =  &
        +k(1298)

    !d[HCO+_dot]/d[O2H+]
    pd(70,113) =  &
        +k(489)*n(idx_CO)

    !d[CH2+_dot]/d[O2H+]
    pd(74,113) =  &
        +k(475)*n(idx_CH)

    !d[CH+_dot]/d[O2H+]
    pd(75,113) =  &
        +k(393)*n(idx_C)

    !d[H2CO+_dot]/d[O2H+]
    pd(76,113) =  &
        +k(646)*n(idx_HCO)

    !d[NH3+_dot]/d[O2H+]
    pd(78,113) =  &
        +k(791)*n(idx_NH2)

    !d[H2O+_dot]/d[O2H+]
    pd(90,113) =  &
        +k(859)*n(idx_OH)

    !d[NH2+_dot]/d[O2H+]
    pd(91,113) =  &
        +k(806)*n(idx_NH)

    !d[OH+_dot]/d[O2H+]
    pd(93,113) =  &
        +k(830)*n(idx_O)

    !d[CH3+_dot]/d[O2H+]
    pd(94,113) =  &
        +k(437)*n(idx_CH2)

    !d[HCN+_dot]/d[O2H+]
    pd(97,113) =  &
        +k(484)*n(idx_CN)

    !d[HNO+_dot]/d[O2H+]
    pd(104,113) =  &
        +k(808)*n(idx_NO)

    !d[H3+_dot]/d[O2H+]
    pd(106,113) =  &
        +k(545)*n(idx_H2)

    !d[H3CO+_dot]/d[O2H+]
    pd(107,113) =  &
        +k(553)*n(idx_H2CO)

    !d[H3O+_dot]/d[O2H+]
    pd(108,113) =  &
        +k(572)*n(idx_H2O)

    !d[HCNH+_dot]/d[O2H+]
    pd(109,113) =  &
        +k(633)*n(idx_HCN)  &
        +k(652)*n(idx_HNC)

    !d[HCO2+_dot]/d[O2H+]
    pd(110,113) =  &
        +k(822)*n(idx_CO2)

    !d[N2H+_dot]/d[O2H+]
    pd(112,113) =  &
        +k(734)*n(idx_N2)

    !d[O2H+_dot]/d[O2H+]
    pd(113,113) =  &
        -k(553)*n(idx_H2CO)  &
        -k(806)*n(idx_NH)  &
        -k(652)*n(idx_HNC)  &
        -k(348)*n(idx_E)  &
        -k(633)*n(idx_HCN)  &
        -k(489)*n(idx_CO)  &
        -k(646)*n(idx_HCO)  &
        -k(1298)  &
        -k(822)*n(idx_CO2)  &
        -k(437)*n(idx_CH2)  &
        -k(545)*n(idx_H2)  &
        -k(572)*n(idx_H2O)  &
        -k(734)*n(idx_N2)  &
        -k(393)*n(idx_C)  &
        -k(808)*n(idx_NO)  &
        -k(484)*n(idx_CN)  &
        -k(830)*n(idx_O)  &
        -k(859)*n(idx_OH)  &
        -k(791)*n(idx_NH2)  &
        -k(475)*n(idx_CH)

    !d[E_dot]/d[SIH5+]
    pd(1,114) =  &
        -k(361)*n(idx_E)  &
        -k(362)*n(idx_E)

    !d[H2_dot]/d[SIH5+]
    pd(6,114) =  &
        +k(361)*n(idx_E)

    !d[H_dot]/d[SIH5+]
    pd(8,114) =  &
        +k(1249)  &
        +k(362)*n(idx_E)

    !d[H2O_dot]/d[SIH5+]
    pd(9,114) =  &
        -k(576)*n(idx_H2O)

    !d[SIH3_dot]/d[SIH5+]
    pd(23,114) =  &
        +k(361)*n(idx_E)

    !d[SIH4_dot]/d[SIH5+]
    pd(32,114) =  &
        +k(576)*n(idx_H2O)  &
        +k(362)*n(idx_E)

    !d[SIH4_DUST_dot]/d[SIH5+]
    pd(48,114) =  &
        +k(1249)

    !d[H3O+_dot]/d[SIH5+]
    pd(108,114) =  &
        +k(576)*n(idx_H2O)

    !d[SIH5+_dot]/d[SIH5+]
    pd(114,114) =  &
        -k(361)*n(idx_E)  &
        -k(576)*n(idx_H2O)  &
        -k(362)*n(idx_E)  &
        -k(1249)

    !d[E_dot]/d[SIOH+]
    pd(1,115) =  &
        -k(364)*n(idx_E)  &
        -k(365)*n(idx_E)

    !d[H_dot]/d[SIOH+]
    pd(8,115) =  &
        +k(365)*n(idx_E)

    !d[OH_dot]/d[SIOH+]
    pd(10,115) =  &
        +k(364)*n(idx_E)

    !d[SI_dot]/d[SIOH+]
    pd(18,115) =  &
        +k(364)*n(idx_E)

    !d[SIO_dot]/d[SIOH+]
    pd(34,115) =  &
        +k(365)*n(idx_E)

    !d[H2SIO_DUST_dot]/d[SIOH+]
    pd(49,115) =  &
        +k(1247)

    !d[SIOH+_dot]/d[SIOH+]
    pd(115,115) =  &
        -k(1247)  &
        -k(364)*n(idx_E)  &
        -k(365)*n(idx_E)

  end subroutine jex

end module krome_ode
