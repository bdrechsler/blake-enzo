
!############### MODULE ##############
module krome_commons
  implicit none

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
  integer,parameter::idx_E=1
  integer,parameter::idx_CH=2
  integer,parameter::idx_O=3
  integer,parameter::idx_HNC=4
  integer,parameter::idx_HCN=5
  integer,parameter::idx_H2=6
  integer,parameter::idx_C=7
  integer,parameter::idx_H=8
  integer,parameter::idx_H2O=9
  integer,parameter::idx_OH=10
  integer,parameter::idx_O2=11
  integer,parameter::idx_CH2=12
  integer,parameter::idx_H2CO=13
  integer,parameter::idx_HCO=14
  integer,parameter::idx_MG=15
  integer,parameter::idx_NH3=16
  integer,parameter::idx_NO=17
  integer,parameter::idx_SI=18
  integer,parameter::idx_SIC2=19
  integer,parameter::idx_SIC3=20
  integer,parameter::idx_SIC=21
  integer,parameter::idx_SIH2=22
  integer,parameter::idx_SIH3=23
  integer,parameter::idx_CN=24
  integer,parameter::idx_CO=25
  integer,parameter::idx_N2=26
  integer,parameter::idx_NH2=27
  integer,parameter::idx_CH3=28
  integer,parameter::idx_CH4=29
  integer,parameter::idx_N=30
  integer,parameter::idx_NH=31
  integer,parameter::idx_SIH4=32
  integer,parameter::idx_SIH=33
  integer,parameter::idx_SIO=34
  integer,parameter::idx_HE=35
  integer,parameter::idx_HNO=36
  integer,parameter::idx_CH3OH=37
  integer,parameter::idx_CO2=38
  integer,parameter::idx_H2CN=39
  integer,parameter::idx_H2SIO=40
  integer,parameter::idx_HNCO=41
  integer,parameter::idx_NO2=42
  integer,parameter::idx_O2H=43
  integer,parameter::idx_OCN=44
  integer,parameter::idx_CH3OH_DUST=45
  integer,parameter::idx_HNCO_DUST=46
  integer,parameter::idx_H2CO_DUST=47
  integer,parameter::idx_SIH4_DUST=48
  integer,parameter::idx_H2SIO_DUST=49
  integer,parameter::idx_SIC_DUST=50
  integer,parameter::idx_SIC2_DUST=51
  integer,parameter::idx_SIC3_DUST=52
  integer,parameter::idx_CH4_DUST=53
  integer,parameter::idx_CO_DUST=54
  integer,parameter::idx_H2O_DUST=55
  integer,parameter::idx_NO_DUST=56
  integer,parameter::idx_CO2_DUST=57
  integer,parameter::idx_N2_DUST=58
  integer,parameter::idx_HCN_DUST=59
  integer,parameter::idx_NH3_DUST=60
  integer,parameter::idx_O2_DUST=61
  integer,parameter::idx_NO2_DUST=62
  integer,parameter::idx_HNO_DUST=63
  integer,parameter::idx_O2H_DUST=64
  integer,parameter::idx_H2CN_DUST=65
  integer,parameter::idx_MG_DUST=66
  integer,parameter::idx_HNC_DUST=67
  integer,parameter::idx_E_DUST=68
  integer,parameter::idx_SIO_DUST=69
  integer,parameter::idx_HCOj=70
  integer,parameter::idx_Hj=71
  integer,parameter::idx_HOCj=72
  integer,parameter::idx_Cj=73
  integer,parameter::idx_CH2j=74
  integer,parameter::idx_CHj=75
  integer,parameter::idx_H2COj=76
  integer,parameter::idx_MGj=77
  integer,parameter::idx_NH3j=78
  integer,parameter::idx_NOj=79
  integer,parameter::idx_SIj=80
  integer,parameter::idx_SIC2j=81
  integer,parameter::idx_SIC3j=82
  integer,parameter::idx_SICj=83
  integer,parameter::idx_SIH2j=84
  integer,parameter::idx_SIH3j=85
  integer,parameter::idx_CNj=86
  integer,parameter::idx_COj=87
  integer,parameter::idx_N2j=88
  integer,parameter::idx_O2j=89
  integer,parameter::idx_H2Oj=90
  integer,parameter::idx_NH2j=91
  integer,parameter::idx_Oj=92
  integer,parameter::idx_OHj=93
  integer,parameter::idx_CH3j=94
  integer,parameter::idx_CH4j=95
  integer,parameter::idx_Nj=96
  integer,parameter::idx_HCNj=97
  integer,parameter::idx_NHj=98
  integer,parameter::idx_SIH4j=99
  integer,parameter::idx_SIHj=100
  integer,parameter::idx_SIOj=101
  integer,parameter::idx_H2j=102
  integer,parameter::idx_HEj=103
  integer,parameter::idx_HNOj=104
  integer,parameter::idx_H2NOj=105
  integer,parameter::idx_H3j=106
  integer,parameter::idx_H3COj=107
  integer,parameter::idx_H3Oj=108
  integer,parameter::idx_HCNHj=109
  integer,parameter::idx_HCO2j=110
  integer,parameter::idx_HEHj=111
  integer,parameter::idx_N2Hj=112
  integer,parameter::idx_O2Hj=113
  integer,parameter::idx_SIH5j=114
  integer,parameter::idx_SIOHj=115
  integer,parameter::idx_CR=116
  integer,parameter::idx_g=117
  integer,parameter::idx_Tgas=118
  integer,parameter::idx_dummy=119
  integer,parameter::nrea=1331
  integer,parameter::nmols=115
  integer,parameter::nspec=119
  integer,parameter::natoms=7
  integer,parameter::ndust=0
  integer,parameter::ndustTypes=0
  integer,parameter::nPhotoBins=0
  integer,parameter::nPhotoRea=0

  integer,parameter::idx_atom_C=1
  integer,parameter::idx_atom_H=2
  integer,parameter::idx_atom_Mg=3
  integer,parameter::idx_atom_O=4
  integer,parameter::idx_atom_N=5
  integer,parameter::idx_atom_Si=6
  integer,parameter::idx_atom_He=7

  !cooling index
  integer,parameter::idx_cool_h2 = 1
  integer,parameter::idx_cool_h2gp = 2
  integer,parameter::idx_cool_atomic = 3
  integer,parameter::idx_cool_cen = 3
  integer,parameter::idx_cool_hd = 4
  integer,parameter::idx_cool_metal = 5
  integer,parameter::idx_cool_z = 5
  integer,parameter::idx_cool_dh = 6
  integer,parameter::idx_cool_enthalpic = 6
  integer,parameter::idx_cool_dust = 7
  integer,parameter::idx_cool_compton = 8
  integer,parameter::idx_cool_cie = 9
  integer,parameter::idx_cool_cont = 10
  integer,parameter::idx_cool_continuum = 10
  integer,parameter::idx_cool_expansion = 11
  integer,parameter::idx_cool_exp = 11
  integer,parameter::idx_cool_ff = 12
  integer,parameter::idx_cool_bss = 12
  integer,parameter::idx_cool_custom = 13
  integer,parameter::idx_cool_co = 14
  integer,parameter::idx_cool_zcie = 15
  integer,parameter::idx_cool_zcienouv = 16
  integer,parameter::idx_cool_zextend = 17
  integer,parameter::idx_cool_gh = 18
  integer,parameter::ncools = 18

  !heating index
  integer,parameter::idx_heat_chem = 1
  integer,parameter::idx_heat_compress = 2
  integer,parameter::idx_heat_compr = 2
  integer,parameter::idx_heat_photo = 3
  integer,parameter::idx_heat_dh = 4
  integer,parameter::idx_heat_enthalpic = 4
  integer,parameter::idx_heat_av = 5
  integer,parameter::idx_heat_photoav = 5
  integer,parameter::idx_heat_cr = 6
  integer,parameter::idx_heat_dust = 7
  integer,parameter::idx_heat_xray = 8
  integer,parameter::idx_heat_viscous = 9
  integer,parameter::idx_heat_visc = 9
  integer,parameter::idx_heat_custom = 10
  integer,parameter::idx_heat_zcie = 11
  integer,parameter::nheats = 11

  real*8::arr_k(nrea)

  !commons for rate tables
  !modify ktab_n according to the required precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating)
  integer::krome_thermo_toggle
  !$omp threadprivate(krome_thermo_toggle)

  !debug bit flag, print and array with fallback values for extreme environments
  integer:: red_flag
  real*8::n_global(nspec)
  integer, save :: nprint_negative=10
  !$omp threadprivate(n_global,nprint_negative,red_flag)

  !commons for implicit RHS
  integer::arr_r1(nrea)
  integer::arr_r2(nrea)
  integer::arr_p1(nrea)
  integer::arr_p2(nrea)
  integer::arr_p3(nrea)
  integer::arr_p4(nrea)

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !synchronize the gamma and mu from hydro code
  real*8::hydro_gamma, hydro_mu
  !$omp threadprivate(hydro_gamma,hydro_mu)

  !commons for frequency bins

  ! Draine dust absorption data loaded from file, via load_kabs
  ! in krome_photo module
  real*8::find_Av_draine_kabs(nPhotoBins)

  !commons for H2 photodissociation (Solomon)
  ! note: paramters here are set depending on the data
  ! but if you have a different file you should modify them
  integer,parameter::H2pdData_nvibX=15
  integer,parameter::H2pdData_nvibB=37
  real*8::H2pdData_dE(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_pre(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_EX(H2pdData_nvibX)
  integer::H2pdData_binMap(H2pdData_nvibX,H2pdData_nvibB)

  !commons for dust optical properties

  !square of turbulence velocity for broadening
  real*8::broadeningVturb2

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
  real*8::user_crflux,user_dgomega,user_zeta,user_rad,user_gRad,user_gArea,user_fr,user_therm,user_desorb,user_h2desorb,user_crdesorb,user_uvcr
  !$omp threadprivate(user_crflux,user_dgomega,user_zeta,user_rad,user_gRad,user_gArea,user_fr,user_therm,user_desorb,user_h2desorb,user_crdesorb,user_uvcr)

  !commons for anytab

  !physical commons
  real*8::phys_Tcmb
  real*8::phys_zredshift
  real*8::phys_orthoParaRatio
  real*8::phys_metallicity
  real*8::phys_Tfloor
  !$omp threadprivate(phys_Tcmb)
  !$omp threadprivate(phys_zredshift)
  !$omp threadprivate(phys_orthoParaRatio)
  !$omp threadprivate(phys_metallicity)
  !$omp threadprivate(phys_Tfloor)

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !total metallicity relative to solar Z/Z_solar
  real*8::total_Z
  real*8::dust2gas_ratio

  !commons for dust tabs (cool,H2,Tdust)
  integer,parameter::dust_tab_imax=50, dust_tab_jmax=50
  real*8::dust_tab_ngas(dust_tab_imax)
  real*8::dust_tab_Tgas(dust_tab_jmax)
  real*8::dust_mult_Tgas,dust_mult_ngas
  real*8::dust_table_AvVariable_log

  real*8::dust_tab_cool(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_heat(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_Tdust(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_H2(dust_tab_imax, dust_tab_jmax)

  !commons for exp(-a) table
  integer,parameter::exp_table_na=int(1d5)
  real*8,parameter::exp_table_aMax=1d4,exp_table_aMin=0d0
  real*8,parameter::exp_table_multa=(exp_table_na-1) &
      / (exp_table_aMax-exp_table_aMin)
  real*8,parameter::exp_table_da=1d0/exp_table_multa
  real*8::exp_table(exp_table_na)

  !stores the last evaluation of the rates in the fex
  real*8::last_coe(nrea)
  !$omp threadprivate(last_coe)

  !xsecs from file variables

  ! Gibbs free energy data from file variables

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !Habing flux for the photoelectric heating by dust
  ! and clumping factor for H2 formation
  ! on dust by Jura/Gnedin
  real*8::GHabing,Ghabing_thin,clump_factor
  !$omp threadprivate(GHabing,GHabing_thin)

  !partition functions common vars

  !verbatim reactions
  character*50::reactionNames(nrea)

end module krome_commons
