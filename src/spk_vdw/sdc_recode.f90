!!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>
!!

!!
!       This implementation is a recode based on the semp_disp_corr.f90
!       written by Erik McNellis and Joerg Meyer in 2008 at the FHI in
!       Berlin. The Input and Output routines were stripped down to 
!       the minimum required in order to integrate with the Atomic
!       Simulation Environment via f2py.
!       Michelitsch Georg, 2013 @ Technische Universitaet Muenchen
!!     

module types
implicit none
integer,parameter,public                                 :: sedc_rk = kind(1D0)
end module

module sdc_recode

use types, only: sedc_rk

implicit none


! TODO: [ ] Error Handling, mutually exclusive statements, etc -> design a check_parameters subroutine
! TODO: [ ] Have stdout-messages of initial module be redirected to an logfile
! TODO: [ ] Set the public/private routines
! TODO: [ ] If there are problems with the C6-parameters,etc it might be caused by the initialization
!           using '= Jnm6_p_mol_to_eVA6 * array' ; this might be circumvented to do explicit multiplication
!           at every access to the parameters or by a subroutine which checks for 'nnnn' and returns '.false.'
!           if so.

! -------------- Variable declaration, PUBLIC ------------------------------------------------------
! -------------- Data type declaration 15pos a.c., max_exponent=300 --------------------------------
!integer,parameter,public                        :: sedc_rk = selected_real_kind(15,300)
!integer,parameter,public                                 :: sedc_rk = kind(1D0)

! -------------- Calculation parameters ------------------------------------------------------------
character(len=50),public,SAVE                            :: logfile_name='sdc_recode.log'
character(len=8),public,SAVE                             :: sedc_scheme='OBS',sedc_xc=''
logical,public,SAVE                                      :: sedc_pbc_img_fixed_nshells=.false.,&
                                                          & sedc_do_num_F=.false.,&
                                                          & sedc_pbc_backfold_coord=.false.,&
                                                          & sedc_variable_cell=.false.,&
                                                          & sedc_do_PBC=.false.
real (kind=sedc_rk),public,SAVE                          :: sedc_pbc_energy_tol= 1e-6_sedc_rk,&
                                                          & sedc_pbc_force_tol= 1e-7_sedc_rk
real (kind=sedc_rk),public,SAVE                          :: sedc_atom_dist_tol= 1e-2_sedc_rk
integer,public,SAVE                                      :: sedc_pbc_n_shells = 0,&
                                                          & sedc_n_groups = 0,&
                                                          & sedc_print_level = 1
logical,public,save                                      :: sedc_do_num_stress=.false.
logical                                                  :: sedc_pbc_file_read=.false.
character(len=3),public                                  :: sedc_pbc_switch='ABC'
logical,dimension(:),allocatable,public,save             :: sedc_pbc_g_only_intra
real (kind=sedc_rk),dimension(:),allocatable,public,save :: sedc_skip_atom
logical,dimension(:),allocatable,public,save             :: sedc_pbc_g_skip
logical,dimension(:),allocatable,public,save             :: sedc_pbc_g_fold
integer,dimension(:,:),allocatable,public,save           :: sedc_pbc_g_switches
logical,public,save                                      :: sedc_do_standalone = .false.

! -------------- Geometry data ---------------------------------------------------------------------
integer,public,SAVE                                             :: sedc_n_ions = 0
integer,dimension(:),allocatable,public,save                    :: sedc_species,&
                                                                 & sedc_groups
real (kind=sedc_rk),dimension(3,3),public,save                  :: sedc_cell_vectors = 0.0_sedc_rk
real (kind=sedc_rk),dimension(:,:),allocatable,public,save      :: sedc_cart_coord,&
                                                                 & sedc_pbc_g_cells

! -------------- I/O routines ----------------------------------------------------------------------
integer,public,SAVE                                :: sedc_stdout_unit=16,sedc_stderr_unit=0
integer,parameter                                  :: max_char_len = 80

! -------------- Results ---------------------------------------------------------------------------
real (kind=sedc_rk),public,SAVE                                 :: sedc_energy = 0.0_sedc_rk
real (kind=sedc_rk),dimension(:),allocatable,public,save        :: sedc_stress
real (kind=sedc_rk),dimension(:,:),allocatable,public,save      :: sedc_forces 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ------------- Variable declaration, PRIVATE ------------------------------------------------------
! -------------- I/O routines ----------------------------------------------------------------------
character(len=4),parameter                         :: lstart='%  ',&
                                                    & tab = '   '
integer                                            :: n_ions_to_print = 0
integer,dimension(:),allocatable                   :: species_to_print
real (kind=sedc_rk),dimension(:,:),allocatable     :: ion_coordinates_to_print
character(len=80)                                  :: xyz_comment='',&
                                                    & xyz_F_comment=''

! -------------- Geometry data ---------------------------------------------------------------------
real(kind=sedc_rk),dimension(:,:),allocatable,save :: c6_ij,R0_ij,internal_cart_coord,&
                                                    & gcell_x_inv_cellvecs

! -------------- Program status --------------------------------------------------------------------
logical,SAVE                                       :: sedc_is_initialized=.false.,&
                                                    & sedc_print_xyz=.false.
integer                                            :: sedc_exec_status = 0

! -------------- Internal Mapping ------------------------------------------------------------------
real (kind=sedc_rk)                                :: central_cell_energy = 0.0_sedc_rk,&
                                                    & sedc_energy_store = 0.0_sedc_rk
real (kind=sedc_rk)                                :: central_cell_volume = 0.0_sedc_rk
real (kind=sedc_rk),SAVE                           :: F_abs_max = 0.0_sedc_rk ! from function-conversion
integer,SAVE                                       :: n_species = 0,&
                                                    & single_point_number = 0
real (kind=sedc_rk),parameter                      :: num_F_coord_shift = 5e-7_sedc_rk
real (kind=sedc_rk),SAVE                           :: E_num_F                 ! from function-conversion
real (kind=sedc_rk),dimension(3,3),SAVE            :: inv_matrix              ! from function-conversion
real (kind=sedc_rk),SAVE                           :: det                     ! from function-conversion
real (kind=sedc_rk),SAVE                           :: sedc_skip_group         ! from function-conversion
real (kind=sedc_rk),SAVE                           :: ions_ij_to_c6_ij        ! from function-conversion  
real (kind=sedc_rk),SAVE                           :: ions_ij_to_R0_ij        ! from function-conversion  
integer,dimension(:),allocatable,save              :: unique_species,ion_species_indices
real (kind=sedc_rk),parameter                      :: num_stress_eps = 1E-3_sedc_rk

! -------------- Verbosity -------------------------------------------------------------------------              
integer,parameter                                  :: xyz_prn_lvl = 2            ! Level for XYZ file printing
integer,parameter                                  :: high_output_prn_lvl = 3    ! High level output
integer,parameter                                  :: high_and_xyz_prn_lvl = 4   ! High level + XYZ file printing
integer,parameter                                  :: debug_prn_lvl = 5          ! Debug (max)

! ------------- Periodic boundary conditions internal mapping --------------------------------------
real (kind=sedc_rk),dimension(:,:),allocatable     :: group_E_decomposed,&
                                                    & group_E_decomposed_store,&
                                                    & sedc_forces_store
real (kind=sedc_rk),dimension(:),allocatable       :: sedc_stress_store
integer,SAVE                                       :: pbc_shells_counted=0
integer,SAVE                                       :: n_cells_in_this_shell     ! from function-conversion
integer,SAVE                                       :: n_ions_in_this_shell      ! from function-conversion
integer,dimension(3),save                          :: pbc_ion_in_s=0,pbc_cell_in_s=0
real (kind=sedc_rk),SAVE                           :: E_num_F_pbc               ! from function-conversion
integer,dimension(:),allocatable                   :: pbc_g_int_scaling

! ------------- Array properties -------------------------------------------------------------------
logical                                            :: scheme_dft_density_indep = .true.
integer,parameter                                  :: param_array_size = 112,&
                                                    & param_label_width = 25
real,parameter                                     :: nnnn = -1.0

! ------------- NIST 2002 constants, N_A Avogadros number ------------------------------------------
real,parameter                                     :: eV_per_J = 6.24150965E18,&
                                                    & N_A = 6.02214179E23
real (kind=sedc_rk),parameter                      :: Jnm6_p_mol_to_eVA6 = (1E6 * eV_per_J) / N_A

! ------------- Periodic table of the elements -----------------------------------------------------
character(len=3),parameter,dimension(param_array_size),private :: periodic_table = (/ &
& 'H  ',                                                                                                'He ', &
& 'Li ','Be ',                                                            'B  ','C  ','N  ','O  ','F  ','Ne ', &
& 'Na ','Mg ',                                                            'Al ','Si ','P  ','S  ','Cl ','Ar ', &
& 'K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ', &
& 'Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ', &
& 'Cs ','Ba ', &
& 'La ','Ce ','Pr ','Nd ','Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ', &
& 'Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ', &
& 'Fr ','Ra ', &
& 'Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ', &
& 'Rf ','Db ','Sg ','Bh ','Hs ','Mt ','Ds ','Rg ','Uub' /)

! TODO: [ ] Move these huge data tables somewhere else (i.e. external library)
!-----------------        Begin OBS scheme parameters and variables        -----------------!

integer,parameter   :: n_param_ar_obs = 3
real (kind=sedc_rk) :: lambda_obs = 0.0_sedc_rk, n_obs = 0.0_sedc_rk

! atomic polarizabilities in 1E-24cm^3 = 1A taken from
! CRC Handbook of Chemistry and Physics, 88th Edition 2007-2008, p.10-194f
! (when several values from different sources are given the one with 
! the best estimated accuracy was chosen)
real (kind=sedc_rk),dimension(param_array_size),parameter :: alpha_obs = (/ &
& 0.666793,                                                                             0.2050522, &
& 24.33,  5.6,                                                  3.03,1.76, 1.1,0.802,0.557,0.3956, &
& 24.11,10.06,                                                   6.8,5.38,3.63,  2.9, 2.18,1.6411, &
&  43.4, 22.8,17.8,14.6,12.4,11.6, 9.4, 8.4, 7.5, 6.8, 6.2,5.75,8.12,6.07,4.31, 3.77, 3.05,2.4844, &
&  47.3, 27.6,22.7,17.9,15.7,12.8,11.4, 9.6, 8.6, 4.8, 7.2,7.36,10.2, 7.7, 6.6,  5.5, 5.35, 4.044, &
& 59.42, 39.7,&
&  31.1, 29.6,28.2,31.4,30.1,28.8,27.7,23.5,25.5,24.5,23.6,22.7,21.8,21.0,21.9, &
&  16.2, 13.1,11.1, 9.7, 8.5, 7.6, 6.5, 5.8,5.02, 7.6, 6.8, 7.4, 6.8, 6.0, 5.3, &
&  48.6, 38.3,&
&  32.1, 32.1,25.4,24.9,24.8,24.5,23.3,23.0,22.7,20.5,19.7,23.8,18.2,17.5,nnnn, &
&  nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

! atomic ionization energies in eV taken from
! http://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page) (2008/03/19)
! (seems to be identical to CRC Handbook of Chemistry and Physics, 88th Edition 2007-2008, p.10-203f)
real (kind=sedc_rk),dimension(param_array_size),parameter :: I_obs = (/ & 
& 13.59844,24.58741, &
& 5.39172, 9.3227,                                         8.29803,11.2603,14.53414,13.61806,17.42282, 21.5646, &
& 5.13908,7.64624,                                         5.98577,8.15169,10.48669,10.36001,12.96764,15.75962, &
& 4.34066,6.11316,6.5615,6.8281,6.7462,6.7665,7.43402,7.9024,7.8810,7.6398,7.72638,9.3942,5.9993,7.8994,9.7886, &
& 9.75238,11.81381,13.99961, &
&  4.17713, 5.6949,6.2171, 6.6339,6.75885,7.09243, 7.2800,7.3605, 7.4589,8.3369, 7.5762,8.9938,5.78636,7.3439,&
& 8.6084,9.0096,10.45126, 12.1298, &
& 3.8939, 5.2117, &
& 5.5769, 5.5387,5.4730, 5.5250, 5.5820, 5.6436, 5.6704,6.1501, 5.8638,5.9389, 6.0215,6.1077,6.18431,6.25416, 5.4259, &
& 6.82507, 7.5496,7.8640, 7.8335, 8.4382, 8.9670, 8.9587,9.2255,10.4375,6.1082,7.41666,7.2856, 8.4170, 9.2242,10.7485, &
& 4.0727, 5.2784, &
& 5.17, 6.3067,5.8900,6.19405, 6.2657, 6.0262, 5.9738,5.9915, 6.1979,6.2817,   6.42,   6.5,   6.58,   6.65,    4.9, &
&  6.0,   nnnn,  nnnn,   nnnn,   nnnn,   nnnn,   nnnn,  nnnn,   nnnn /)

! covalent radii in Angstrom taken from
!	http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)	(2008/03/18)
! (seem to be identical to data at http://www.webelements.com)
real (kind=sedc_rk),dimension(param_array_size),parameter      :: R0_obs = (/ &
& 0.37,                                                                                0.32, &
& 1.34,0.90,                                                  0.82,0.77,0.75,0.73,0.71,0.69, &
& 1.54,1.30,                                                  1.18,1.11,1.06,1.02,0.99,0.97, &
& 1.96,1.74,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.19,1.16,1.14,1.10, &
& 2.11,1.92,1.62,1.48,1.37,1.45,1.56,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33,1.30, &
& 2.25,1.98, &
& 1.69,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,1.60, &
& 1.50,1.38,1.46,1.59,1.28,1.37,1.28,1.44,1.49,1.48,1.47,1.46,nnnn,nnnn,1.45, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

!-----------------        End OBS scheme parameters and variables          -----------------!
!-----------------       Begin G06 scheme parameters and variables         -----------------!

integer,parameter   :: n_param_ar_g06 = 2
real (kind=sedc_rk) :: s6_g06 = 0.0_sedc_rk, d_g06 = 0.0_sedc_rk

real (kind=sedc_rk),dimension(param_array_size),parameter :: C6i_g06 = Jnm6_p_mol_to_eVA6 * (/ &
&  0.14,                                                                                                 0.08, &
&  1.61, 1.61,                                                             3.13, 1.75, 1.23, 0.70, 0.75, 0.63, &
&  5.71, 5.71,                                                            10.79, 9.23, 7.84, 5.57, 5.07, 4.61, &
& 10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,16.99,17.10,16.37,12.64,12.47,12.01, &
& 24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,37.32,38.71,38.44,31.74,31.50,29.99, &
&  nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, &
&  nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn /)

real (kind=sedc_rk),dimension(param_array_size),parameter :: R0_g06 = (/ &
& 1.001,                                                                                                1.012, &
& 0.825,1.408,                                                            1.485,1.452,1.397,1.342,1.287,1.243, &
& 1.144,1.364,                                                            1.639,1.716,1.705,1.683,1.639,1.595, &
& 1.485,1.474,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.650,1.727,1.760,1.771,1.749,1.727, &
& 1.628,1.606,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.672,1.804,1.881,1.892,1.892,1.881, &
&  nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, &
&  nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, &
&  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn /)

!-----------------         End G06 scheme parameters and variables         -----------------!
!-----------------      Begin JCHS scheme parameters and variables         -----------------!

integer,parameter   :: n_param_ar_jchs = 3
real (kind=sedc_rk) :: SR_jchs = 0.0_sedc_rk, d_jchs = 0.0_sedc_rk, &
& s6_jchs = 0.0_sedc_rk

real (kind=sedc_rk),dimension(param_array_size),parameter :: C6i_jchs = Jnm6_p_mol_to_eVA6 * (/ &
& 0.16,                                                                                nnnn, &
& nnnn,nnnn,                                                  nnnn,1.65,1.11,0.70,0.57,0.45, &
& nnnn,nnnn,                                                  nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

! Parameters taken from J . Am. Chem. Soc.,Vol. 114, No. 20, (1992)
! as indicated in JCHS reference
real (kind=sedc_rk),dimension(param_array_size),parameter :: N_eff_jchs = (/ &
& 0.80,                                                                                nnnn, &
& nnnn,nnnn,                                                  nnnn,2.49,2.82,3.15,3.48,nnnn, &
& nnnn,nnnn,                                                  nnnn,nnnn,nnnn,nnnn,5.10,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,6.00,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

real (kind=sedc_rk),dimension(param_array_size),parameter :: Bondi_vdW_radii = (/ &
& 1.20,                                                                                1.40, &
& 1.82,nnnn,                                                  nnnn,1.70,1.55,1.52,1.47,1.54, &
& 2.27,1.73,                                                  nnnn,2.10,1.80,1.80,1.75,1.88, &
& 2.75,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,1.63,1.40,1.39,1.87,nnnn,1.85,1.90,1.85,2.02, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,1.63,1.72,1.58,1.93,2.17,nnnn,2.06,1.98,2.16, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,1.75,1.66,1.55,1.96,2.02,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,1.86,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

!-----------------        End JCHS scheme parameters and variables         -----------------!
!-----------------        Begin TS scheme parameters and variables         -----------------!

integer,parameter   :: n_param_ar_ts = 3
real (kind=sedc_rk) :: SR_ts = 0.0_sedc_rk, d_ts = 0.0_sedc_rk
real,parameter      :: Ha_to_eV = 27.2113838668, Bohr_to_Ang = 0.52917720859 ! NIST 2002
real (kind=sedc_rk),dimension(:),allocatable,public,save :: sedc_ts_veff_div_vfree

!~ ! These two arrays are the database (in Ha * Bohr ^ 6 and Bohr ^ 3) published in 
!~ ! (X. Chu and A. Dalgarno, J. Chem. Phys. 121, 4083 (2004))
!~ real (kind=sedc_rk),dimension(param_array_size),parameter :: C6i_ts = Ha_to_eV * (Bohr_to_Ang ** 6) * (/ &
!~ & nnnn,                                                                                1.42, &
!~ & 1392, 227,                                                  99.5,46.6,24.2,15.6,9.52,6.20, &
!~ & 1518, 626,                                                   528, 305, 185, 134,94.6,64.2, &
!~ & 3923,2163,1383,1044, 832, 602, 552, 482, 408, 373, 253, 284, 498, 354, 246, 210, 162, 130, &
!~ & 4769,3175,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, 779, 659, 492, 445, 385,nnnn, &
!~ & nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!~ & nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

!~ real (kind=sedc_rk),dimension(param_array_size),parameter :: alpha_ts = (Bohr_to_Ang ** 3) * (/ &
!~ & nnnn,                                                                                1.38, &
!~ &  164,  38,                                                    21,  12, 7.4, 5.4, 3.8,2.67, &
!~ &  163,  71,                                                    60,  37,  25,19.6,  15,11.1, &
!~ &  294, 160, 120,  98,  84,  78,  63,  56,  51,  48,  42,  40,  60,  41,  29,  25,  20,16.7, &
!~ &  320, 199,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,  75,  60,  44,  40,  35,nnnn, &
!~ & nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!~ & nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!~ & nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

!!! included within CASTEP up to 5.5.2c:
! These two arrays are the database (in Ha * Bohr ^ 6 and Bohr ^ 3) published in 
! (X. Chu and A. Dalgarno, J. Chem. Phys. 121, 4083 (2004)),
! with missing values replaced by those computed by AT and MS
!real (kind=sedc_rk),dimension(param_array_size),parameter :: C6i_ts = Ha_to_eV * (Bohr_to_Ang ** 6) * (/ &
!&   6.50,                                                                                                    1.42, &
!& 1392.0, 227.0,                                                               99.5, 46.6, 24.2, 15.6, 9.52, 6.20, &
!& 1518.0, 626.0,                                                              528.0,305.0,185.0,134.0, 94.6, 64.2, &
!& 3923.0,2163.0,1383.0,1044.0,832.0,602.0,552.0,482.0,408.0,373.0,253.0,284.0,498.0,354.0,246.0,210.0,162.0,130.0, &
!& 4769.0,3175.0,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn,339.0, nnnn,779.0,659.0,492.0,445.0,385.0,286.0, &
!&   nnnn,  nnnn, &
!&   nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn,nnnn, &
!&   nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn,298.0, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn,nnnn, &
!&   nnnn,  nnnn, &
!&   nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn,nnnn, &
!&   nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn /)
!
!real (kind=sedc_rk),dimension(param_array_size),parameter :: alpha_ts = (Bohr_to_Ang ** 3) * (/ &
!& 4.50,                                                                                  1.38, &
!&164.0, 38.0,                                                   21.0,12.0, 7.4, 5.4, 3.8,2.67, &
!&163.0, 71.0,                                                   60.0,37.0,25.0,19.6,15.0,11.1, &
!&294.0,160.0,120.0,98.0,84.0,78.0,63.0,56.0,51.0,48.0,42.0,40.0,60.0,41.0,29.0,25.0,20.0,16.7, &
!&320.0,199.0, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,50.6,nnnn,75.0,60.0,44.0,40.0,35.0,27.3, &
!& nnnn, nnnn, &
!& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,36.5,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!& nnnn, nnnn, &
!& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)
!
! Van der Waals radii used by A.T. and M.S. in Bohr
!real (kind=sedc_rk),dimension(param_array_size),parameter :: R0i_ts = Bohr_to_Ang * (/ &
!& 3.10,                                                                                2.65, &
!& 4.16,4.17,                                                  3.89,3.59,3.34,3.19,3.04,2.91, &
!& 3.73,4.27,                                                  4.33,4.20,4.01,3.86,3.71,3.55, &
!& 3.71,4.65,4.59,4.51,4.44,3.99,3.97,4.23,4.18,3.82,3.76,4.02,4.19,4.20,4.11,4.04,3.93,3.82, &
!& 3.72,4.54,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,3.82,nnnn,nnnn,nnnn,nnnn,nnnn,4.17,4.08, &
!& nnnn,nnnn, &
!& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,3.86,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!& nnnn,nnnn, &
!& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
!& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

! latest values for vdW-TS, PRL 102, 073005 (2009)
!
! Mitch: Latest values for Lanthanoids from Alexander Tkatchenko / Vivekanand
!        Gobre, private communication (TS-scheme!)
real (kind=sedc_rk),dimension(param_array_size),parameter :: alpha_ts = (Bohr_to_Ang ** 3) * (/ &
& 4.50,                                                                                  1.38, &
&164.2, 38.0,                                                   21.0,12.0, 7.4, 5.4, 3.8,2.67, &
&162.7, 71.0,                                                   60.0,37.0,25.0,19.6,15.0,11.1, &
&292.9,160.0,120.0,98.0,84.0,78.0,63.0,56.0,50.0,48.0,42.0,40.0,60.0,41.0,29.0,25.0,20.0,16.8, &
&319.2,199.0, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,56.1,23.68,50.6,39.7,75.0,60.0,44.0,37.65,35.0,27.3, &
& nnnn, 275.0, &
& 213.7,204.7,215.8,208.4,200.2,192.1,184.2,158.3,169.5,164.64,156.3,150.2,144.3,138.9,137.2, &
& nnnn, nnnn, nnnn,nnnn,nnnn,42.51,39.68,36.5,33.9,nnnn,61.8,49.02,nnnn,nnnn,nnnn, &
& nnnn, nnnn, &
& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)
!
real (kind=sedc_rk),dimension(param_array_size),parameter :: C6i_ts = Ha_to_eV * (Bohr_to_Ang ** 6) * (/ &
& 6.50,                                                                                           1.46, &
& 1387.0, 214.0,                                                               99.5, 46.6, 24.2, 15.6, 9.52, 6.38, &
& 1556.0, 627.0,                                                              528.0,305.0,185.0,134.0, 94.6, 64.3, &
& 3897.0,2221.0,1383.0,1044.0,832.0,602.0,552.0,482.0,408.0,373.0,253.0,284.0,498.0,354.0,246.0,210.0,162.0,129.6, &
& 4691.0,3170.0,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn,469.0,157.5,339.0, 452.0,779.0,659.0,492.0,396.0,385.0,285.9, &
&   nnnn,5727.0, &
& 3884.5,3708.33,3911.84,3908.75,3847.68,3708.69,3511.71,2781.53,3124.41,2984.29,2839.95,2724.12,2576.78,2387.53,2371.80, &
&   nnnn,  nnnn,  nnnn,  nnnn, nnnn,359.1,347.1,298.0,392.0, nnnn,697.0,571.0, nnnn, nnnn,nnnn, &
&   nnnn,  nnnn, &
&   nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn,nnnn, &
&   nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn /)
!
real (kind=sedc_rk),dimension(param_array_size),parameter :: R0i_ts = Bohr_to_Ang * (/ &
& 3.10,                                                                                2.65, &
& 4.16,4.17,                                                  3.89,3.59,3.34,3.19,3.04,2.91, &
& 3.73,4.27,                                                  4.33,4.20,4.01,3.86,3.71,3.55, &
& 3.71,4.65,4.59,4.51,4.44,3.99,3.97,4.23,4.18,3.82,3.76,4.02,4.19,4.20,4.11,4.04,3.93,3.82, &
& 3.72,4.54,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,3.95,3.66,3.82,3.99,nnnn,nnnn,nnnn,4.22,4.17,4.08, &
& nnnn,4.77, &
& 3.16,3.26,3.28,3.30,3.27,3.32,3.40,3.62,3.42,3.26,3.24,3.30,3.26,3.22,3.20, &
& nnnn,nnnn,nnnn,nnnn,nnnn,4.00,3.92,3.86,3.98,nnnn,4.31,4.32,nnnn,nnnn,nnnn, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

!-----------------         End TS scheme parameters and variables           -----------------!
!-----------------       Begin TSsurf scheme parameters and variables       -----------------!
!In this scheme all parameters have to be put in manually, C6, alpha and R_vdW 
! in Ha*bohr**6, bohr**3 and in bohr
!Also the Hirshfeld volume has to be put in manually in bohr**3
!The periodic table data of TS is still put, but only to be overwritten!

!
! Mitch: Latest values for Lanthanoids from Alexander Tkatchenko / Vivekanand
!        Gobre, private communication (TS-scheme - as adsorbate only!)

  integer,parameter   :: n_param_ar_tssurf = 3
  real (kind=sedc_rk) :: SR_tssurf = 0.0_sedc_rk, d_tssurf = 0.0_sedc_rk
!REINI
  real (kind=sedc_rk),dimension(:),allocatable,public,save :: sedc_tssurf_vfree_div_vbulk
  logical, public,save:: sedc_tssurf_coeffs_not_set=.true.

real (kind=sedc_rk),dimension(param_array_size),parameter :: alpha_tssurf = (Bohr_to_Ang ** 3) * (/ &
& 4.50,                                                                                    1.38, &
& nnnn, nnnn,                                                     21.0,12.0, 7.4, 5.4, 3.8,2.67, &
& nnnn, nnnn,                                                     60.0,37.0,25.0,19.6,15.0,11.1, &
& nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,10.22,10.88,13.77,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,13.90,15.36, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn, nnnn, &
& 213.7,204.7,215.8,208.4,200.2,192.1,184.2,158.3,169.5,164.64,156.3,150.2,144.3,138.9,137.2, &
& nnnn, nnnn, nnnn,nnnn,nnnn,13.2,14.45,15.62,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn, nnnn, &
& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn, nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)
!
real (kind=sedc_rk),dimension(param_array_size),parameter :: C6i_tssurf = Ha_to_eV * (Bohr_to_Ang ** 6) * (/ &
& 6.50,                                                                                                   1.46, &
& nnnn,  nnnn,                                                              99.5, 46.6, 24.2, 15.6, 9.52, 6.38, &
& nnnn,  nnnn,                                                             528.0,305.0,185.0,134.0, 94.6, 64.3, &
& nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, 59.2, 58.9,46.0,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn,102.0,122.0,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,  nnnn, &
& 3884.5,3708.33,3911.84,3908.75,3847.68,3708.69,3511.71,2781.53,3124.41,2984.29,2839.95,2724.12,2576.78,2387.53,2371.80, &
& nnnn,  nnnn,  nnnn,  nnnn, nnnn, 98.0, 120.5, 133.9, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, &
& nnnn,  nnnn, &
& nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn, nnnn,nnnn, &
& nnnn,  nnnn,  nnnn,  nnnn, nnnn, nnnn, nnnn, nnnn, nnnn /)
!
real (kind=sedc_rk),dimension(param_array_size),parameter :: R0i_tssurf = Bohr_to_Ang * (/ &
& 3.10,                                                                                   2.65, &
& nnnn, nnnn,                                                    3.89,3.59,3.34,3.19,3.04,2.91, &
& nnnn, nnnn,                                                    4.33,4.20,4.01,3.86,3.71,3.55, &
& nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,2.28, 2.40, 2.82,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,3.06, 2.57, nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn, nnnn, &
& 3.16,3.26,3.28,3.30,3.27,3.32,3.40,3.62,3.42,3.26,3.24,3.30,3.26,3.22,3.20, &
& nnnn,nnnn,nnnn,nnnn,nnnn,2.71,2.80,2.91,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn, &
& nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn,nnnn /)

!  real (kind=sedc_rk),dimension(:),allocatable  :: C6i_tssurf
!  real (kind=sedc_rk),dimension(:),allocatable  :: alpha_tssurf
!  real (kind=sedc_rk),dimension(:),allocatable  :: R0i_tssurf

!parameters are being calculated as in all other schemes and are additionally modified 
! by sedc_tssurf_vfree_div_vbulk

!-----------------       End TS-SURF scheme parameters and variables         -----------------!


contains      

subroutine sedc(exec_status)
implicit none
integer,intent(out)         :: exec_status

exec_status = 0

! Print eveything to sdc_recode.log sedc_stdout_unit was = 6 !
! TODO: [ ] Maybe think of better way?
open(unit=sedc_stdout_unit,file=logfile_name,status='unknown',action='write',position='append')


if (sedc_exec_status == 0) then
    select case(sedc_scheme)
        case("OBS")
            call initialize_main('init_obs','Phys. Rev. B 73, 205101, (2006)')
            call sedc_main('OBS')
        case("G06")
            call initialize_main('init_g06','J. Comput. Chem. 27, 1787, (2006)')
            call sedc_main('G06')
        case("JCHS")
            call initialize_main('init_jchs','J. Comput. Chem. 28, 555, (2007)')
            call sedc_main('JCHS')
        case("TS")
            call initialize_main('init_ts','Phys. Rev. Lett. 102, 073005 (2009)')
            call sedc_main('TS')
        case("TS-SURF")
            call initialize_main('init_tssurf','Phys. Rev. Lett (2012), in press')
            call sedc_main('TS-SURF')
        case ("")
            call sedc_errquit('SEDC Initialization ERROR!',&
            & 'No dispersion correction scheme given.')
        case default
            call sedc_errquit('SEDC Initialization ERROR!',&
            & "Requested '" // trim(sedc_scheme) // "' scheme is not implemented.")
    end select
endif

close(unit=sedc_stdout_unit)

exec_status = sedc_exec_status

end subroutine

! 'OBS' scheme functional initialization routine
subroutine init_obs(xc)
implicit none
character(len=*),intent(in)                                   :: xc
real (kind=sedc_rk),dimension(n_param_ar_obs)              :: param_errvals_obs = 0.0_sedc_rk
character(len=param_label_width),dimension(n_param_ar_obs) :: param_ar_labels_obs = ''
character(len=param_label_width),dimension(n_param_ar_obs) :: param_ar_names_obs = ''

select case(xc)
    case("PZ81","PW91")
        lambda_obs = 7.5E-4_sedc_rk; n_obs = 8.0_sedc_rk
    case default
        call sedc_errquit("SEDC Initialization ERROR!","Requested '" &
        & // xc // "' functional is not implemented in the 'OBS' scheme.")
        return
end select

if (sedc_print_level == debug_prn_lvl) &
& write (sedc_stdout_unit,"(5(A),/,8(A),ES7.1,A,F4.1)") &
& lstart,tab,"'OBS' functionals available",tab,": PZ81,PW91", &
& lstart,tab,"'OBS'/",xc," parameters",tab,tab,": 'Lambda' = ",lambda_obs, &
& "; 'n' = ",n_obs

! OBS scheme has an overall factor -1 in the exponential argument;
! multiply in here
lambda_obs = -1 * lambda_obs

param_errvals_obs = (/ nnnn,nnnn,nnnn /)

param_ar_labels_obs(1) = 'Polarizability [A^3]'
param_ar_labels_obs(2) = 'Ionization potential [eV]'
param_ar_labels_obs(3) = 'Covalent radius [A]'

param_ar_names_obs(1)  = 'alpha_obs'
param_ar_names_obs(2)  = 'I_obs'
param_ar_names_obs(3)  = 'R0_obs'

call put_c6_and_r0_ij('OBS',param_errvals_obs,param_ar_labels_obs,&
& param_ar_names_obs,n_param_ar_obs,alpha_obs,I_obs,R0_obs)

end subroutine


! 'OBS' dispersion coefficient / vdW radius mean computation routine
subroutine c6_r0_ij_obs(ion_i,ion_j,species_i,species_j,c6_ij_obs,r0_ij_obs)
integer,intent(in)              :: ion_i,ion_j,species_i,species_j
real (kind=sedc_rk),intent(out) :: c6_ij_obs,r0_ij_obs

c6_ij_obs = &
& 1.5_sedc_rk * (alpha_obs(species_i) * alpha_obs(species_j) &
&             *  I_obs(species_i) * I_obs(species_j)) &
&             / (I_obs(species_i) + I_obs(species_j))

r0_ij_obs = R0_obs(species_i) + R0_obs(species_j)
end subroutine


! 'OBS' scheme damping function * dispersion coefficients
subroutine c_f_obs(ion_i,ion_j,r,cf,cdf_dr)
implicit none
integer,intent(in)              :: ion_i,ion_j
real (kind=sedc_rk),intent(in)  :: r
real (kind=sedc_rk),intent(out) :: cf,cdf_dr
real (kind=sedc_rk)             :: x_ij=0,arg=0,c_exp_of_arg=0,c=0

call ions_ij_to_R0_ij_ff(ion_i,ion_j)
x_ij = r / ions_ij_to_R0_ij

call ions_ij_to_c6_ij_ff(ion_i,ion_j)
c = ions_ij_to_c6_ij

arg = lambda_obs * (x_ij ** n_obs)
c_exp_of_arg = exp(arg) * c
cf = c_exp_of_arg - c
cdf_dr = n_obs * (arg/r) * c_exp_of_arg
end subroutine


! 'G06' scheme initialization routine
subroutine init_g06(xc)
implicit none
character(len=*),intent(in)                                   :: xc
real (kind=sedc_rk)                                        :: n = 0.0_sedc_rk
real (kind=sedc_rk),dimension(n_param_ar_g06)              :: param_errvals_g06 = 0.0_sedc_rk
character(len=param_label_width),dimension(n_param_ar_g06) :: param_ar_labels_g06 = ''
character(len=param_label_width),dimension(n_param_ar_g06) :: param_ar_names_g06 = ''

select case(xc)
    case("PBE")
        s6_g06 = 0.75_sedc_rk
    case("BLYP")
        s6_g06 = 1.20_sedc_rk
    case("BP86","B3LYP")
        s6_g06 = 1.05_sedc_rk
    case("TPSS")
        s6_g06 = 1.00_sedc_rk
    case default
        call sedc_errquit("SEDC Initialization ERROR!","Requested '" &
        & // xc // "' functional is not implemented in the 'G06' scheme.")
        return
end select

! Grimme does not optimize 'd'; it is kept fixed at 20
d_g06 = 20.0_sedc_rk

if (sedc_print_level == debug_prn_lvl) &
& write (sedc_stdout_unit,"(5(A),/,7(A),A,F5.2,A,F4.1)") &
& lstart,tab,"'G06' functionals available",tab,": PBE,BLYP,BP86,B3LYP,TPSS", &
& lstart,tab,"'G06'/",xc," parameters",tab,tab,": 's6' = ",s6_g06,",'d' = ",d_g06

n = 1.0 * nnnn
param_errvals_g06 = (/ Jnm6_p_mol_to_eVA6 * n,n /)

param_ar_labels_g06(1) = 'C_6 coefficient [eV*A^6]'
param_ar_labels_g06(2) = 'van der Waals radius [A]'

param_ar_names_g06(1)  = 'C6i_g06'
param_ar_names_g06(2)  = 'R0_g06'

call put_c6_and_r0_ij('G06',param_errvals_g06,param_ar_labels_g06,&
& param_ar_names_g06,n_param_ar_g06,C6i_g06,R0_g06)
if (sedc_exec_status /= 0) return

c6_ij= s6_g06 * c6_ij
R0_ij = (-1.0_sedc_rk * d_g06)/R0_ij

end subroutine


! 'G06' dispersion coefficient / vdW radius mean computation routine
subroutine c6_r0_ij_g06(ion_i,ion_j,species_i,species_j,c6_ij_g06,r0_ij_g06)
integer,intent(in)              :: ion_i,ion_j,species_i,species_j
real (kind=sedc_rk),intent(out) :: c6_ij_g06,r0_ij_g06

c6_ij_g06 = sqrt( C6i_g06(species_i) * C6i_g06(species_j) )
r0_ij_g06 = R0_g06(species_i) + R0_g06(species_j)

end subroutine


! 'G06' scheme damping function * dispersion coefficients
subroutine c_f_g06(ion_i,ion_j,r,cf,cdf_dr)
implicit none
integer,intent(in)              :: ion_i,ion_j
real (kind=sedc_rk),intent(in)  :: r
real (kind=sedc_rk),intent(out) :: cf,cdf_dr
real (kind=sedc_rk)             :: exp_of_arg=0,f=0,m=0

call ions_ij_to_R0_ij_ff(ion_i,ion_j)
m = ions_ij_to_R0_ij
exp_of_arg = exp( r * m + d_g06 )
f          = -1 / (1 + exp_of_arg)             ! -f
call ions_ij_to_c6_ij_ff(ion_i,ion_j)
cf         = ions_ij_to_c6_ij * f ! E
cdf_dr     = cf * f * exp_of_arg * m           ! nabla -f 
end subroutine


! 'JCHS' scheme initialization routine
subroutine init_jchs(xc)
implicit none
character(len=*),intent(in)                                    :: xc
real (kind=sedc_rk)                                         :: n = 0.0_sedc_rk
real (kind=sedc_rk),dimension(n_param_ar_jchs)              :: param_errvals_jchs = 0.0_sedc_rk
character(len=param_label_width),dimension(n_param_ar_jchs) :: param_ar_labels_jchs = ''
character(len=param_label_width),dimension(n_param_ar_jchs) :: param_ar_names_jchs = ''

! Parameters are chosen as the published combination with lowest average error
! of the mixed test sets, of the 'LP', counterpoise corrected 'LP' or 'aQZ' basis sets
select case(xc)
    case("PBE")
        SR_jchs = 1.02_sedc_rk; d_jchs = 33.0_sedc_rk; s6_jchs = 1.00_sedc_rk
    case("BLYP")
        SR_jchs = 1.25_sedc_rk; d_jchs = 23.0_sedc_rk; s6_jchs = 0.90_sedc_rk
    case("B3LYP","TPSS")
        SR_jchs = 0.93_sedc_rk; d_jchs = 35.0_sedc_rk; s6_jchs = 1.00_sedc_rk
    case default
        call sedc_errquit("SEDC Initialization ERROR!","Requested '" &
        & // xc // "' functional is not implemented in the 'JCHS' scheme.")
        return
end select

if (sedc_print_level == debug_prn_lvl) &
& write (sedc_stdout_unit,"(5(A),/,7(A),A,F4.1,2(A,F5.2))") &
& lstart,tab,"'JCHS' functionals available",tab,": PBE,BLYP,B3LYP,TPSS", &
& lstart,tab,"'JCHS'/",xc," parameters",tab,tab,": 'd' = ", &
& d_jchs,"; 'S_R' = ",SR_jchs,"; 's6' = ",s6_jchs

n = 1.0 * nnnn
param_errvals_jchs = (/ Jnm6_p_mol_to_eVA6 * n,n,n /)

param_ar_labels_jchs(1) = 'C_6 coefficient [eV*A^6]'
param_ar_labels_jchs(2) = 'Effective electron number'
param_ar_labels_jchs(3) = 'van der Waals radius [A]'

param_ar_names_jchs(1)  = 'C6i_jchs'
param_ar_names_jchs(2)  = 'N_eff_jchs'
param_ar_names_jchs(3)  = 'R0_jchs'

call put_c6_and_r0_ij('JCHS',param_errvals_jchs,param_ar_labels_jchs,&
& param_ar_names_jchs,n_param_ar_jchs,C6i_jchs,N_eff_jchs,Bondi_vdw_radii)
if (sedc_exec_status /= 0) return

c6_ij = s6_jchs * c6_ij
R0_ij = (-1.0_sedc_rk * (d_jchs/SR_jchs))/R0_ij

end subroutine


! 'JCHS' dispersion coefficient / vdW radius mean computation routine
subroutine c6_r0_ij_jchs(ion_i,ion_j,species_i,species_j,c6_ij_jchs,r0_ij_jchs)
integer,intent(in)              :: ion_i,ion_j,species_i,species_j
real (kind=sedc_rk),intent(out) :: c6_ij_jchs,r0_ij_jchs
real (kind=sedc_rk)             :: c_1=0,c_2=0,n_1=0,n_2=0
c_1 = C6i_jchs(species_i);   c_2 = C6i_jchs(species_j)
n_1 = N_eff_jchs(species_i); n_2 = N_eff_jchs(species_j)

c6_ij_jchs = (2 * (c_1 * c_1 * c_2 * c_2 * n_1 * n_2) ** (1.0_sedc_rk/3.0_sedc_rk)) / &
& ((c_1 * n_2 * n_2) ** (1.0_sedc_rk/3.0_sedc_rk) + (c_2 * n_1 * n_1) ** (1.0_sedc_rk/3.0_sedc_rk))

r0_ij_jchs = (Bondi_vdw_radii(species_i) ** 3 + Bondi_vdw_radii(species_j) ** 3) / &
& (Bondi_vdw_radii(species_i) ** 2 + Bondi_vdw_radii(species_j) ** 2)

end subroutine


! 'JCHS' scheme damping function * dispersion coefficients
subroutine c_f_jchs(ion_i,ion_j,r,cf,cdf_dr)
implicit none
integer,intent(in)              :: ion_i,ion_j
real (kind=sedc_rk),intent(in)  :: r
real (kind=sedc_rk),intent(out) :: cf,cdf_dr
real (kind=sedc_rk)             :: exp_of_arg=0,m=0,f=0

call ions_ij_to_R0_ij_ff(ion_i,ion_j)
m = ions_ij_to_R0_ij
! exp(-d(r_ij/(S_R*R^0_ij)-1))
exp_of_arg = exp( r * m + d_jchs )
f          = -1 / (1 + exp_of_arg)             ! -f	
call ions_ij_to_c6_ij_ff(ion_i,ion_j)
cf         = ions_ij_to_c6_ij * f ! E
cdf_dr     = cf * f * m * exp_of_arg           ! nabla -f
end subroutine


! 'TS' scheme initialization routine
subroutine init_ts(xc)
implicit none
character(len=*),intent(in)                               :: xc
real (kind=sedc_rk)                                       :: n = 0.0_sedc_rk
real (kind=sedc_rk),dimension(n_param_ar_ts)              :: param_errvals_ts = 0.0_sedc_rk
character(len=param_label_width),dimension(n_param_ar_ts) :: param_ar_labels_ts = ''
character(len=param_label_width),dimension(n_param_ar_ts) :: param_ar_names_ts = ''
integer                                                   :: i

scheme_dft_density_indep = .false.

!mijp: cannot do this as don't know order of execution of clauses in 'if'
!if((.not.(allocated(sedc_ts_veff_div_vfree))) .or. (size(sedc_ts_veff_div_vfree) /= sedc_n_ions)) then
!     call sedc_errquit("SEDC Initialization ERROR!",&
!     & "'TS' scheme requires atomic V_eff / V_free ratios to function.")
!     return
if (allocated(sedc_ts_veff_div_vfree)) then
      if (size(sedc_ts_veff_div_vfree) /= sedc_n_ions) then
         call sedc_errquit("SEDC Initialization ERROR!",&
         & "'TS' scheme requires atomic V_eff / V_free ratios to function.")
         return
      end if
else
         call sedc_errquit("SEDC Initialization ERROR!",&
         & "'TS' scheme requires atomic V_eff / V_free ratios to function.")
         return
endif


select case(xc)
    case("PBE")
        SR_ts = 0.94_sedc_rk
    case("RPBE")
        SR_ts = 0.77_sedc_rk
    case("PBE0","PBE1PBE")
        SR_ts = 0.96_sedc_rk
    case("BLYP")
        SR_ts = 0.62_sedc_rk
    case("B3LYP","AM05")
        SR_ts = 0.84_sedc_rk
    case default
        call sedc_errquit("SEDC Initialization ERROR!","Requested '" &
        & // xc // "' functional is not implemented in the 'TS' scheme.")
        return
end select

d_ts = 20.0_sedc_rk

if (sedc_print_level == debug_prn_lvl) &
& write (sedc_stdout_unit,"(5(A),/,7(A),A,F4.1,A,F5.2)") &
& lstart,tab,"'TS' functionals available",tab,": PBE,RPBE,PBE0,BLYP,B3LYP,AM05", &
& lstart,tab,"'TS'/",xc," parameters",tab,tab,": 'd' = ", &
& d_ts,"; 'S_R' = ",SR_ts

n = 1.0 * nnnn
param_errvals_ts = (/ Ha_to_eV * (Bohr_to_Ang ** 6) * n,(Bohr_to_Ang ** 3) * n, Bohr_to_Ang * n /)

param_ar_labels_ts(1) = 'C_6 coefficient [eV*A^6]'
param_ar_labels_ts(2) = 'Polarizability [A^3]'
param_ar_labels_ts(3) = 'van der Waals radius [A]'

param_ar_names_ts(1)  = 'C6i_ts'
param_ar_names_ts(2)  = 'alpha_ts'
param_ar_names_ts(3)  = 'R0i_ts'

call put_c6_and_r0_ij('TS',param_errvals_ts,param_ar_labels_ts,param_ar_names_ts,n_param_ar_ts, &
& C6i_ts,alpha_ts,R0i_ts)

if (sedc_exec_status /= 0) return

! TODO: [ ] This is purely an stdout-debug print statement, removed
if (sedc_print_level == debug_prn_lvl) &
&  call print_vector('     Atom ',sedc_ts_veff_div_vfree, &
&  (/ (periodic_table(sedc_species(i)),i = 1,sedc_n_ions) /),'V_eff / V_free','F',12,8)

R0_ij = (-1.0_sedc_rk * (d_ts/SR_ts))/R0_ij

end subroutine

subroutine c6_r0_ij_ts(ion_i,ion_j,species_i,species_j,c6_ij_ts,r0_ij_ts)
integer,intent(in)              :: ion_i,ion_j,species_i,species_j
real (kind=sedc_rk),intent(out) :: c6_ij_ts,r0_ij_ts

! V_ratio = U
! C6_AA = C6_A * U_A ^ 2
! alpha_AA = alpha_A * U_A
! R0_AA = R0_A * U_A ^ (1/3)
! 
! C6_AB = 2 * C6_AA * C6_BB / (( C6_AA * alpha_BB) / alpha_AA + ( C6_BB * alpha_AA) / alpha_BB))
! R0_AB = R0_AA + R0_BB
! 
! Following expression is simplified (one U_A * U_B term cancels between numerator / denominator,
! and one U_{A,B} term cancels in each fraction in the denominator)
 c6_ij_ts = (sedc_ts_veff_div_vfree(ion_i) * C6i_ts(species_i) * &
&           sedc_ts_veff_div_vfree(ion_j) * C6i_ts(species_j) * 2.0_sedc_rk) / &
&          ((alpha_ts(species_j) * C6i_ts(species_i)) / alpha_ts(species_i) + &
&           (alpha_ts(species_i) * C6i_ts(species_j)) / alpha_ts(species_j))

 r0_ij_ts = (R0i_ts(species_i) * sedc_ts_veff_div_vfree(ion_i) ** (0.3333333333333333333_sedc_rk) + &
&           R0i_ts(species_j) * sedc_ts_veff_div_vfree(ion_j) ** (0.3333333333333333333_sedc_rk))
end subroutine


! 'TS' scheme damping function * dispersion coefficients
subroutine c_f_ts(ion_i,ion_j,r,cf,cdf_dr)
implicit none
integer,intent(in)              :: ion_i,ion_j
real (kind=sedc_rk),intent(in)  :: r
real (kind=sedc_rk),intent(out) :: cf,cdf_dr
real (kind=sedc_rk)             :: exp_of_arg=0,f=0
exp_of_arg = exp( r * R0_ij(ion_i,ion_j) + d_ts )
f      = -1 / (1 + exp_of_arg)                    ! -f		
cf     = c6_ij(ion_i,ion_j) * f                   ! E
cdf_dr = cf * f * R0_ij(ion_i,ion_j) * exp_of_arg ! nabla -f
end subroutine

! 'TSsurf' scheme initialization routine
subroutine init_tssurf(xc)
  implicit none
  character(len=*),intent(in)                                   :: xc
  real (kind=sedc_rk)                                           :: n = 0.0_sedc_rk
  real (kind=sedc_rk),dimension(n_param_ar_tssurf)              :: param_errvals_tssurf = 0.0_sedc_rk
  character(len=param_label_width),dimension(n_param_ar_tssurf) :: param_ar_labels_tssurf = ''
  character(len=param_label_width),dimension(n_param_ar_tssurf) :: param_ar_names_tssurf = ''
  integer                                                       :: i

  scheme_dft_density_indep = .false.

  !mijp: cannot do this as don't know order of execution of clauses in 'if'
  !if((.not.(allocated(sedc_ts_veff_div_vfree))) .or. (size(sedc_ts_veff_div_vfree) /= sedc_n_ions)) then
  !     call sedc_errquit("SEDC Initialization ERROR!",&
  !     & "'TS1103' scheme requires atomic V_eff / V_free ratios to function.")
  !     return


  if (allocated(sedc_ts_veff_div_vfree)) then
        if (size(sedc_ts_veff_div_vfree) /= sedc_n_ions) then
          call sedc_errquit("SEDC Initialization ERROR!",&
          & "'TS-SURF' scheme requires atomic V_eff / V_free ratios to function.")
          return
        end if
  else
          call sedc_errquit("SEDC Initialization ERROR!",&
          & "'TS-SURF' scheme requires atomic V_eff / V_free ratios to function.")
          return
  endif

!REINI 
  if (allocated(sedc_tssurf_vfree_div_vbulk)) then
     if (size(sedc_tssurf_vfree_div_vbulk) /= sedc_n_ions) then
        call sedc_errquit("SEDC Initialization ERROR!", &
&       "'TS-SURF' scheme requires sedc_tssurf_vfree_div_vbulk to function.")
        return
     end if
  else
     call sedc_errquit("SEDC Initialization ERROR!", &
&    "'TS-SURF' scheme requires sedc_tssurf_vfree_div_vbulk to function.")
     return
  endif

  select case(xc)
    case("PBE")
        SR_tssurf = 0.94_sedc_rk
    case("RPBE")
        SR_tssurf = 0.77_sedc_rk
    case("PBE0","PBE1PBE")
        SR_tssurf = 0.96_sedc_rk
    case("BLYP")
        SR_tssurf = 0.62_sedc_rk
    case("B3LYP","AM05")
        SR_tssurf = 0.84_sedc_rk
    case default
        call sedc_errquit("SEDC Initialization ERROR!","Requested '" &
        & // xc // "' functional is not implemented in the 'TS-SURF' scheme.")
        return
  end select

  d_tssurf= 20.0_sedc_rk

  if (sedc_print_level == debug_prn_lvl) &
&    write (sedc_stdout_unit,"(5(A),/,7(A),A,F4.1,A,F5.2)") &
&           lstart,tab,"'TS-SURF' functionals available",tab,": PBE,RPBE,PBE0,BLYP,B3LYP,AM05", &
&           lstart,tab,"'TS-SURF'/",xc," parameters",tab,tab,": 'd' = ", &
&           d_tssurf,"; 'S_R' = ",SR_tssurf

  n = 1.0 * nnnn
  param_errvals_tssurf = (/ Ha_to_eV * (Bohr_to_Ang ** 6) * n,(Bohr_to_Ang ** 3) * n, Bohr_to_Ang * n /)

  param_ar_labels_tssurf(1) = 'C_6 coefficient [eV*A^6]'
  param_ar_labels_tssurf(2) = 'Polarizability [A^3]'
  param_ar_labels_tssurf(3) = 'van der Waals radius [A]'

  param_ar_names_tssurf(1)  = 'C6i_tssurf'
  param_ar_names_tssurf(2)  = 'alpha_tssurf'
  param_ar_names_tssurf(3)  = 'R0i_tssurf'

  call put_c6_and_r0_ij('TS-SURF', &
&   param_errvals_tssurf,param_ar_labels_tssurf,param_ar_names_tssurf,n_param_ar_tssurf, &
&   C6i_tssurf,alpha_tssurf,R0i_tssurf)


! TODO: [ ] This is purely a debug-stdout statement, removed.  
   if (sedc_print_level == debug_prn_lvl) &
&    call print_vector('     Atom ',sedc_ts_veff_div_vfree, &
&    (/ (periodic_table(sedc_species(i)),i = 1,sedc_n_ions) /),'V_eff / V_free','F',12,8)
!REINI
  if (sedc_print_level == debug_prn_lvl) &
&    call print_vector('     Atom ',sedc_tssurf_vfree_div_vbulk, &
&    (/ (periodic_table(sedc_species(i)),i = 1,sedc_n_ions) /),'V_free / V_bulk','F',12,8)

  R0_ij = (-1.0_sedc_rk * (d_tssurf/SR_tssurf))/R0_ij

  if (sedc_exec_status /= 0) return

end subroutine

subroutine c6_r0_ij_tssurf(ion_i,ion_j,species_i,species_j,c6_ij_tssurf,r0_ij_tssurf)
  integer,intent(in)              :: ion_i,ion_j,species_i,species_j
  real (kind=sedc_rk),intent(out) :: c6_ij_tssurf,r0_ij_tssurf
  ! see comments above in 
  !   subroutine c6_r0_ij_ts
  
  real (kind=sedc_rk)             :: temp_vratio_i, temp_vratio_j
 
  
  temp_vratio_i= sedc_ts_veff_div_vfree(ion_i) * sedc_tssurf_vfree_div_vbulk(ion_i)
  temp_vratio_j= sedc_ts_veff_div_vfree(ion_j) * sedc_tssurf_vfree_div_vbulk(ion_j)
  !write(*,*) species_i, ' ', species_j
  !write(*,*) sedc_tssurf_vfree_div_vbulk(species_i), ' ', sedc_tssurf_vfree_div_vbulk(species_j)

  c6_ij_tssurf = (temp_vratio_i * C6i_tssurf(species_i) * temp_vratio_j * C6i_tssurf(species_j) *  &
&                 2.0_sedc_rk) / &
&                ((alpha_tssurf(species_j) * C6i_tssurf(species_i)) / alpha_tssurf(species_i) + &
&                 (alpha_tssurf(species_i) * C6i_tssurf(species_j)) / alpha_tssurf(species_j))
  r0_ij_tssurf = (R0i_tssurf(species_i) * temp_vratio_i ** (0.3333333333333333333_sedc_rk) + &
&                 R0i_tssurf(species_j) * temp_vratio_j ** (0.3333333333333333333_sedc_rk))


!   c6_ij_ts1103 = (vbulk_vfree_ratio_ts1103(species_i) * sedc_ts_veff_div_vfree(ion_i) * C6i_ts1103(species_i) * &
! &                 vbulk_vfree_ratio_ts1103(species_j) * sedc_ts_veff_div_vfree(ion_j) * C6i_ts1103(species_j) *  &
! &                 2.0_sedc_rk) / &
! &                ((alpha_ts1103(species_j) * C6i_ts1103(species_i)) / alpha_ts1103(species_i) + &
! &                 (alpha_ts1103(species_i) * C6i_ts1103(species_j)) / alpha_ts1103(species_j))
!   r0_ij_ts1103 = (R0i_ts1103(species_i) * vbulk_vfree_ratio_ts1103(species_i) * sedc_ts_veff_div_vfree(ion_i) ** &
! &                 (0.3333333333333333333_sedc_rk) + &
! &                 R0i_ts1103(species_j) * vbulk_vfree_ratio_ts1103(species_j) * sedc_ts_veff_div_vfree(ion_j) ** &
! &                 (0.3333333333333333333_sedc_rk))

end subroutine

! 'TSsurf' -> 'TS' scheme damping function * dispersion coefficients
subroutine c_f_tssurf(ion_i,ion_j,r,cf,cdf_dr)
  implicit none
  integer,intent(in)              :: ion_i,ion_j
  real (kind=sedc_rk),intent(in)  :: r
  real (kind=sedc_rk),intent(out) :: cf,cdf_dr
  real (kind=sedc_rk)             :: exp_of_arg=0,f=0

  !if (sedc_do_standalone) then
  !  exp_of_arg = exp( r * R0_ij(ion_i,ion_j) * 0.30_sedc_rk + d_tssurf )
  !else
  exp_of_arg = exp( r * R0_ij(ion_i,ion_j) + d_tssurf )
  !endif
  f      = -1 / (1 + exp_of_arg) 

  if (sedc_do_standalone) then
    f = -1.0_sedc_rk
    exp_of_arg = 0.0_sedc_rk
  end if

  cf     = c6_ij(ion_i,ion_j) * f                   ! E
  cdf_dr = cf * f * R0_ij(ion_i,ion_j) * exp_of_arg ! nabla -f
end subroutine 

!!! Quick'n'dirty-f2py fix mitch circumventing the usage of external function declaration
subroutine scheme_selector_c6(ion_i,ion_j,species_i,species_j,c6_ij_obs,r0_ij_obs)
integer,intent(in)              :: ion_i,ion_j,species_i,species_j
real (kind=sedc_rk),intent(out) :: c6_ij_obs,r0_ij_obs

select case(sedc_scheme)
case("OBS")
    call c6_r0_ij_obs(ion_i,ion_j,species_i,species_j,c6_ij_obs,r0_ij_obs)
case("G06")
    call c6_r0_ij_g06(ion_i,ion_j,species_i,species_j,c6_ij_obs,r0_ij_obs)
case("JCHS")
    call c6_r0_ij_jchs(ion_i,ion_j,species_i,species_j,c6_ij_obs,r0_ij_obs)
case("TS")
    call c6_r0_ij_ts(ion_i,ion_j,species_i,species_j,c6_ij_obs,r0_ij_obs)
case("TS-SURF")
    call c6_r0_ij_tssurf(ion_i,ion_j,species_i,species_j,c6_ij_obs,r0_ij_obs)
case ("")
    call sedc_errquit('SEDC Initialization ERROR!',&
    & "No C6 parameter scheme given.")
case default
    call sedc_errquit('SEDC Initialization ERROR!',&
    & "Requested '" // trim(sedc_scheme) // "' scheme is not implemented.")
end select
end subroutine


! TODO: [x] dummy function
subroutine sedc_main(c_f)
implicit none
character(len=*),intent(in)                                  :: c_f

if (sedc_exec_status /= 0) return

! begin: stress tensor
call sedc_central_E_F_stress(c_f)
! end: stress tensor

central_cell_energy = sedc_energy

if (sedc_do_PBC) call sedc_pbc(c_f,sedc_print_level >= high_output_prn_lvl,sedc_pbc_energy_tol)

! begin: stress tensor
  ! stress gets divided by volume (i.e. is in pressure units afterwards)
  ! obviously, central_cell_volume must have been properly set before (currently done in initialize_main)
  sedc_stress = sedc_stress / central_cell_volume
! end: stress tensor

if (sedc_print_xyz) call sedc_prn_xyz()

if (sedc_do_num_F.or.sedc_do_num_stress) then
    sedc_energy_store = sedc_energy
    allocate(sedc_forces_store(3,sedc_n_ions))
    sedc_forces_store = sedc_forces
    allocate(sedc_stress_store(6))
    sedc_stress_store = sedc_stress
    allocate(group_E_decomposed_store(3,max(1,sedc_n_groups)))
    group_E_decomposed_store = group_E_decomposed
    
    if (sedc_do_num_F) call sedc_central_num_F(c_f,sedc_do_PBC)

    ! begin: stress tensor
    if (sedc_do_num_stress) call sedc_num_stress(c_f)
    ! end: stress tensor
    
    sedc_energy = sedc_energy_store
    sedc_forces = sedc_forces_store
    deallocate(sedc_forces_store)
    sedc_stress = sedc_stress_store
    deallocate(sedc_stress_store)
    group_E_decomposed = group_E_decomposed_store
    deallocate(group_E_decomposed_store)
endif
! TODO: [ ] This function call was removed, keep for later preliminary testing
        call sedc_print_results(trim(sedc_scheme),trim(sedc_xc),sedc_n_ions)

end subroutine

! TODO: dummy function
subroutine initialize_main(init_c,reference)
implicit none
!external                                       :: init_c
! Mitch: circumventing usage of EXTERNAL (f2py doesnt like it)
character(len=*)                                  :: init_c
character(len=*)                                  :: reference
integer                                        :: species_i=0,i=0,g=0,n_in_g=0,jcount=0
integer,dimension(:),allocatable               :: tmp_unique_species
character(len=max_char_len)                    :: a = '',b = ''
real (kind=sedc_rk),dimension(:,:),allocatable :: folded_coords
logical                                        :: do_exit = .false.

if (.not.(allocated(sedc_species))) then
    do_exit = .true.
    a = trim(a) // "Species array ('sedc_species')"
endif

if (.not.(allocated(sedc_cart_coord))) then
    if (do_exit) a = trim(a) // ' and'
    do_exit = .true.
    b = trim(b) // "Coordinate array ('sedc_cart_coord')"
endif

if (do_exit) then
    a = 'SEDC ERROR! ' // trim(a)
    b = trim(b) // ' not allocated.'
    call sedc_errquit(a,b)
    return
endif

sedc_is_initialized=.false. !mijp

if (.not.sedc_is_initialized) then

    if (sedc_print_level > 0) then
    
        call print_header(1,3,'DFT-SEDC: DFT SemiEmpirical Dispersion interaction Correction module')

        write (sedc_stdout_unit,"(A,4(/,(3A)))") lstart,lstart, &
        & tab,'Copyright (c) 2008 Erik McNellis, Joerg Meyer, Fritz-Haber-Institut', &
        & lstart,tab,'der MPG. Distributed under the GNU Lesser General Public License (LGPL).', &
        & lstart,tab,'Please cite "Erik R. McNellis, Joerg Meyer, and Karsten Reuter,', &
        & lstart,tab,'Phys. Rev. B 80, 205414 (2009)" in all publications using this module.'

        if (sedc_print_level >= high_output_prn_lvl) then 
            write(a,"(3(A))") adjustl("Using scheme (see: "),trim(reference),")"
            write(b,"(A)") "Parameters optimized for DFT x-c functional"
            write (sedc_stdout_unit,"(A,/,2(A),A56,3(A),/,2(A),A56,2(A))") lstart, &
            & lstart,tab,a,": '",trim(adjustl(sedc_scheme)),"'", &
            & lstart,tab,adjustl(b),":  ",trim(adjustl(sedc_xc))
            if (sedc_print_level == debug_prn_lvl) write (sedc_stdout_unit,"(A)") lstart
        endif
    endif
    
    call group_atoms()
    if (sedc_exec_status /= 0) return
    
    if (allocated(group_E_decomposed)) deallocate(group_E_decomposed) !mijp
    allocate(group_E_decomposed(3,max(1,sedc_n_groups)))
    if (allocated(internal_cart_coord)) deallocate(internal_cart_coord) !mijp
    allocate(internal_cart_coord(3,sedc_n_ions))
    if (allocated(sedc_forces)) deallocate(sedc_forces) !mijp
    allocate(sedc_forces(3,sedc_n_ions))
!!! begin: stress tensor
    if (allocated(sedc_stress)) deallocate(sedc_stress) !mijp
    allocate(sedc_stress(6))
!!! end: stress tensor
    
    ! numerical stress calculation currently requires variable_cell to be set
    ! already at initialization stage
    if (sedc_do_num_stress) sedc_variable_cell  = .true.
endif

call reinitialize_main()

if (.not.sedc_is_initialized) then
    
    if ((sedc_print_level == xyz_prn_lvl) .or. &
    & (sedc_print_level == high_and_xyz_prn_lvl)) sedc_print_xyz = .true.

    if (allocated(ion_species_indices)) deallocate(ion_species_indices) !mijp
    allocate(ion_species_indices(sedc_n_ions))
    ion_species_indices = 0
        
    allocate(tmp_unique_species(sedc_n_ions))
    tmp_unique_species = 0

    n_species=0 !mijp
    do i = 1,sedc_n_ions
        species_i = sedc_species(i)
        if (.not.any(mask= tmp_unique_species == species_i)) then
            n_species = n_species + 1
            tmp_unique_species(n_species) = species_i
            where(sedc_species == species_i) &
            & ion_species_indices = n_species
        endif
    enddo
    
    if (allocated(unique_species)) deallocate(unique_species) !mijp
    allocate(unique_species(n_species))
    unique_species = 0
        
    unique_species = tmp_unique_species(1:n_species)
        
    deallocate(tmp_unique_species)

!    call init_c(trim(sedc_xc))
    call scheme_selector_init(trim(init_c))
    if (sedc_exec_status /= 0) return
        
    sedc_pbc_backfold_coord = .false.

!REINI

    if (sedc_do_PBC) then
        if (any(sedc_pbc_g_fold)) sedc_pbc_backfold_coord = .true.
        allocate(folded_coords(3,sedc_n_ions)); jcount = 0
        do g = 1,sedc_n_groups
            n_in_g = sedc_groups(g)
            call fold_group_coords( g, jcount + 1, jcount + n_in_g, .true., &
            & folded_coords(:,(jcount + 1):(jcount + n_in_g)) )
            if (sedc_exec_status /= 0) return
            jcount = jcount + n_in_g
        enddo
        jcount=0
        if (any(sedc_pbc_g_skip)) then
            do g = 1,sedc_n_groups
                n_in_g = sedc_groups(g)
                if (sedc_pbc_g_skip(g)) then
                    do i = 1, n_in_g
                        sedc_skip_atom(i+jcount) = 0.0_sedc_rk
                    end do
                end if
                jcount = jcount + n_in_g
             enddo
         endif
    endif

    call check_distance_matrix(internal_cart_coord,internal_cart_coord,'atom','atom',2)
    if (sedc_exec_status /= 0) return
        
    if (sedc_do_PBC) then
        call check_distance_matrix(internal_cart_coord,folded_coords,'atom','image of atom',1)
        if (sedc_exec_status /= 0) return
    
        call check_distance_matrix(folded_coords,internal_cart_coord,'image of atom','atom',2)
        if (sedc_exec_status /= 0) return
        deallocate(folded_coords)
    endif
        
    if (sedc_print_level >= high_output_prn_lvl) then
        call setting_summary()
    endif
    
    sedc_is_initialized = .true.
    
else

    if (.not.scheme_dft_density_indep) then
        !call init_c(trim(sedc_xc))
        call scheme_selector_init(trim(init_c))
        if (sedc_exec_status /= 0) return
    endif
    
endif

end subroutine 


!!! Quick'n'dirty-f2py fix mitch circumventing the usage of external function declaration
subroutine scheme_selector_init(xc)
character(len=*),intent(in)                                   :: xc

select case(sedc_scheme)
case("OBS")
    call init_obs(trim(sedc_xc)) 
case("G06")
    call init_g06(trim(sedc_xc)) 
case("JCHS")
    call init_jchs(trim(sedc_xc)) 
case("TS")
    call init_ts(trim(sedc_xc)) 
case("TS-SURF")
    call init_tssurf(trim(sedc_xc)) 
case ("")
    call sedc_errquit('SEDC Initialization ERROR!',&
    & "No SEDC_XC parameter scheme given.")
case default
    call sedc_errquit('SEDC Initialization ERROR!',&
    & "Requested '" // trim(sedc_xc) // "' functional is not implemented.")
end select
end subroutine

subroutine reinitialize_main()
implicit none
integer :: jcount,g,n_in_g

!!! begin: stress tensor
  ! calculate volume of central cell (parallelepipedial product of sedc_cell_vectors!)
    call det_ff( sedc_cell_vectors(:,1), sedc_cell_vectors(:,2), sedc_cell_vectors(:,3) ) 
    central_cell_volume = abs(det)
!!! end: stress tensor

if (sedc_do_PBC) then
    call sedc_init_pbc()
    if (sedc_exec_status /= 0) return
endif

single_point_number = single_point_number + 1
internal_cart_coord = 0.0_sedc_rk

if (sedc_pbc_backfold_coord) then
    jcount = 0
    do g = 1,sedc_n_groups
        n_in_g = sedc_groups(g)
        if (sedc_pbc_g_fold(g)) then
            call fold_group_coords( g, jcount + 1, jcount + n_in_g, .false., &
            & internal_cart_coord(:,(jcount + 1):(jcount + n_in_g)))
            if (sedc_exec_status /= 0) return
        else
            internal_cart_coord(:,(jcount + 1):(jcount + n_in_g)) = &
            & sedc_cart_coord(:,(jcount + 1):(jcount + n_in_g))
        endif
        jcount = jcount + n_in_g
    enddo
else
    internal_cart_coord = sedc_cart_coord
endif

sedc_energy = 0.0_sedc_rk; sedc_forces = 0.0_sedc_rk
group_E_decomposed = 0.0_sedc_rk; sedc_stress = 0.0_sedc_rk

end subroutine


subroutine check_distance_matrix(coord_ar_one,coord_ar_two,msg_one,msg_two,i_start)
implicit none
real (kind=sedc_rk),dimension(:,:),intent(in) :: coord_ar_one,coord_ar_two
character(len=*),intent(in)                      :: msg_one,msg_two
integer,intent(in)                            :: i_start
character(len=max_char_len)                   :: a = '', b = '', err_msg = ''
real (kind=sedc_rk),dimension(3)              :: v = 0.0_sedc_rk
real (kind=sedc_rk)                           :: r = 0.0_sedc_rk
integer                                       :: s = 0, i = 0, j = 0

s = size(coord_ar_one) / 3


do i = i_start, s
    do j = 1, i-1
        v = coord_ar_two(:,j) - coord_ar_one(:,i)
        r = sqrt(dot_product(v,v))
        ! TODO: [ ] remove this debug statement!!!
        ! be careful about order of matrix elements, atoms object has to be transposed before
        ! handed to sdc_recode - routines!
            !write (a,"(I5)") i; write (b,"(I5)") j
            !write (*,*) coord_ar_one(0,i),' i=',i,'/',j
            !write (*,*) coord_ar_one(1,i),' i=',i,'/',j
            !write (*,*) coord_ar_one(2,i),' i=',i,'/',j
            !write (*,*) '/---------------------------/'
            !write (a,"(I5)") i; write (b,"(I5)") j
            !write (*,*) coord_ar_one(0,j),' ',i,'/j=',j
            !write (*,*) coord_ar_one(1,j),' ',i,'/j=',j
            !write (*,*) coord_ar_one(2,j),' ',i,'/j=',j
            !write (sedc_stdout_unit,"(7(A,1X),ES7.1,A)") 'Distance between',trim(msg_one),&
            !& trim(adjustl(a)),'and',trim(msg_two),trim(adjustl(b)), &
            !& '=',r,' [A]'
        if (r < sedc_atom_dist_tol) then
            write (a,"(I5)") i; write (b,"(I5)") j
            write (err_msg,"(7(A,1X),ES7.1,A)") 'Distance between',trim(msg_one),&
            & trim(adjustl(a)),'and',trim(msg_two),trim(adjustl(b)), &
            & '<',sedc_atom_dist_tol,' [A]'
            call sedc_errquit('SEDC Initialization ERROR!',err_msg)
            return
        endif
    enddo
enddo

end subroutine


subroutine fold_group_coords(g,a,b,check,folded_coords)
implicit none
integer,intent(in)                                       :: g,a,b
logical,intent(in)                                       :: check
real (kind=sedc_rk),dimension(3,(b - a + 1)),intent(out) :: folded_coords
real (kind=sedc_rk),dimension(3,3)                       :: cell,inv_cell
real (kind=sedc_rk),dimension(3)                         :: fr,fo
integer                                                  :: i,j

folded_coords = 0.0_sedc_rk; cell = sedc_pbc_g_cells(:,(3*g-2):(3*g))
call inv_matrix_ff(cell)
inv_cell = inv_matrix; if (sedc_exec_status /= 0) return; j = 0

if (any(mask= sedc_pbc_g_switches(:,g) == 1)) then
    do i = a,b
        fr = matmul(inv_cell,sedc_cart_coord(:,i))
        fo = modulo(fr,1.0_sedc_rk)
        if (check) where ( &
        & (fr <  0.0_sedc_rk .and. sedc_pbc_g_switches(1:3,g) == 0) .or. & 
        & (fr >= 1.0_sedc_rk .and. sedc_pbc_g_switches(4:6,g) == 0)) fo = fr
        j = j + 1; folded_coords(:,j) = matmul(cell,fo)
    enddo
else
    folded_coords(:,1:(b - a + 1)) = sedc_cart_coord(:,a:b)
endif

end subroutine

subroutine put_c6_and_r0_ij(c6_r0_ij,ar_errvals,ar_labels,ar_names,n_ar, &
& ar_1,ar_2,ar_3,ar_4,ar_5,ar_6)
implicit none
character(len=*),intent(in)                             :: c6_r0_ij
integer,intent(in)                                   :: n_ar
real (kind=sedc_rk),dimension(:),intent(in)          :: ar_errvals
character(len=*),dimension(:),intent(in)                :: ar_labels,ar_names
real (kind=sedc_rk),dimension(:),intent(in),optional :: ar_1,ar_2,ar_3,ar_4,ar_5,ar_6
real (kind=sedc_rk),dimension(:,:),allocatable       :: tmp_parray
character(len=10),dimension(:),allocatable           :: labels
character(len=max_char_len)                          :: a,b
integer                                              :: i=0,j=0

if (.not.sedc_is_initialized) then
    a = ''

    allocate(tmp_parray(n_species,n_ar))

    do i = 1, n_species
        j = unique_species(i)
        if (present(ar_1)) tmp_parray(i,1) = ar_1(j)
        if (present(ar_2)) tmp_parray(i,2) = ar_2(j)
        if (present(ar_3)) tmp_parray(i,3) = ar_3(j)
        if (present(ar_4)) tmp_parray(i,4) = ar_4(j)
        if (present(ar_5)) tmp_parray(i,5) = ar_5(j)
        if (present(ar_6)) tmp_parray(i,6) = ar_6(j)
    enddo

    do i = 1, n_species
        do j = 1, n_ar
            if (tmp_parray(i,j) == ar_errvals(j)) then
                if (a /= '') a = trim(a) // ","
                a = trim(a) // " " // "'" // trim(ar_names(j)) // "'" 
            endif
        enddo

        if (a /= '') then
            call sedc_errquit( &
                & "SEDC ERROR! Parameter Not Available in Scheme '" // trim(sedc_scheme) &
                & // "':", "Element '" // trim(periodic_table(unique_species(i))) &
                & // "' not in array;" // trim(a))
            return
        endif
    enddo
endif

if (scheme_dft_density_indep) then
    call fill_r0_c6_arrays(c6_r0_ij,unique_species,n_species)
else
    call fill_r0_c6_arrays(c6_r0_ij,sedc_species,sedc_n_ions)
endif

!REINI
if (sedc_scheme .eq. 'TS-SURF') &
    & call print_matrix('Ion dep. parameters:',1,tmp_parray,ar_labels,&
    & (/ (periodic_table(unique_species(i)), i = 1, n_species) /),'ES',10,3)

        ! TODO: [ ] This is only a stdout-debug statement, removed.
if (sedc_print_level == debug_prn_lvl) then
      
    if (scheme_dft_density_indep) then
        allocate(labels(n_species))
        labels = (/ (periodic_table(unique_species(i)), i = 1, n_species) /)
    else
        allocate(labels(sedc_n_ions))
        do i = 1,sedc_n_ions
            write(a,"(I4)") i; write(b,"(4(A))") &
            & trim(adjustl(periodic_table(sedc_species(i)))),'(',trim(adjustl(a)),')'
            labels(i) = trim(b)
        enddo
    endif    

    if (sedc_scheme /= 'TS-SURF' .and. .not.sedc_is_initialized) & 
        & call print_matrix('Species dep. parameters:',4,tmp_parray,ar_labels,&
        & (/ (periodic_table(unique_species(i)), i = 1, n_species) /),'ES',10,3)
        
    call print_matrix('C_6ij [eV*A^6]',5,c6_ij,labels,labels,'ES',10,3)
    call print_matrix('     R0_ij [A]',5,R0_ij,labels,labels,'ES',10,3)
    
    deallocate(labels)
endif

if (allocated(tmp_parray)) deallocate(tmp_parray)

end subroutine


subroutine fill_r0_c6_arrays(c6_r0_ij,species,s)
implicit none
character(len=*),intent(in)        :: c6_r0_ij
integer,dimension(:),intent(in) :: species
integer,intent(in)              :: s
integer                         :: i,j
real (kind=sedc_rk)             :: c_ij_value=0,r0_ij_value=0

if (.not.sedc_is_initialized) then
    if (allocated(c6_ij)) deallocate(c6_ij) !mijp
    if (allocated(R0_ij)) deallocate(R0_ij) !mijp
    allocate(c6_ij(s,s)); allocate(R0_ij(s,s))
endif

c6_ij = 0.0_sedc_rk; R0_ij = 0.0_sedc_rk
do i = 1,s
    do j = 1,i
        !!! Mitch
        !call c6_r0_ij(i,j,species(i),species(j),c_ij_value,r0_ij_value)
        call scheme_selector_c6(i,j,species(i),species(j),c_ij_value,r0_ij_value)
        c6_ij(i,j) = c_ij_value;  c6_ij(j,i) = c_ij_value
        R0_ij(i,j) = r0_ij_value; R0_ij(j,i) = r0_ij_value
    enddo
enddo
end subroutine


subroutine sedc_central_E_F_stress(c_f)
        implicit none
        character(len=*),intent(in)                       :: c_f
        integer                                        :: i = 0, j = 0
        real (kind=sedc_rk)                            :: tmp_E
        real (kind=sedc_rk),dimension(3)               :: F_i = 0.0_sedc_rk
        real (kind=sedc_rk),dimension(:,:),allocatable :: pair_terms

        ! TODO: [ ] It seems this block does output handling, keep until final implementation is finished
        ! TODO: [ ] Delete associated global variables
        if (sedc_print_xyz) then
            n_ions_to_print = sedc_n_ions

            if (allocated(species_to_print)) deallocate(species_to_print) !mijp
            allocate(species_to_print(n_ions_to_print))
            species_to_print = sedc_species

            if (allocated(ion_coordinates_to_print)) deallocate(ion_coordinates_to_print) !mijp
            allocate(ion_coordinates_to_print(3,n_ions_to_print))
            ion_coordinates_to_print = internal_cart_coord

            xyz_comment   = '! Coordinates in [A] of central cell'
            xyz_F_comment = '! Coords [A], SEDC force corr. [eV/A] of central cell'
        endif
        ! END TODO

        allocate(pair_terms(sedc_n_ions,sedc_n_ions))
        pair_terms = 0.0_sedc_rk



        do i = 2,sedc_n_ions
            do j = 1,i-1
                tmp_E = sedc_energy
                call pair_E_F_stress(c_f,i,j,internal_cart_coord(:,j),F_i)
                sedc_forces(:,j) = sedc_forces(:,j) - F_i
                tmp_E = sedc_energy - tmp_E
                pair_terms(i,j) = tmp_E
                pair_terms(j,i) = tmp_E
            enddo
        enddo

        if (sedc_n_groups > 1) then

            i = 0
            do j = 1,sedc_n_groups
                group_E_decomposed(1,j) = sum(pair_terms((i+1):(i+sedc_groups(j)),(i+1):(i+sedc_groups(j))))
                group_E_decomposed(2,j) = sum(pair_terms((i+1):(i+sedc_groups(j)),:)) - group_E_decomposed(1,j)
                i = i + sedc_groups(j)
            enddo
            
            group_E_decomposed = group_E_decomposed / 2.0_sedc_rk
        endif

        deallocate(pair_terms)

end subroutine


!!! Mitch-f2py-quick'n'dirty-fix
!!! begin case-selecting code
subroutine scheme_selector(ion_i, ion_j, r, cf, cdf_dr)
integer,intent(in)              :: ion_i,ion_j
real (kind=sedc_rk),intent(in)  :: r
real (kind=sedc_rk),intent(out) :: cf,cdf_dr

select case(sedc_scheme)
case("OBS")
    call c_f_obs(ion_i, ion_j, r, cf, cdf_dr)
case("G06")
    call c_f_g06(ion_i, ion_j, r, cf, cdf_dr)
case("JCHS")
    call c_f_jchs(ion_i, ion_j, r, cf, cdf_dr)
case("TS")
    call c_f_ts(ion_i, ion_j, r, cf, cdf_dr)
case("TS-SURF")
    call c_f_tssurf(ion_i, ion_j, r, cf, cdf_dr)
case ("")
    call sedc_errquit('SEDC Initialization ERROR!',&
    & "No dispersion correction scheme given.")
case default
    call sedc_errquit('SEDC Initialization ERROR!',&
    & "Requested '" // trim(sedc_scheme) // "' scheme is not implemented.")
end select
end subroutine


!!! begin: stress tensor
subroutine pair_E_F_stress(c_f,ion_i,ion_j,j_coords,F_i)
implicit none
character(len=*),intent(in)                     :: c_f
integer,intent(in)                           :: ion_i,ion_j
real (kind=sedc_rk),dimension(3),intent(in)  :: j_coords
real (kind=sedc_rk),dimension(3),intent(out) :: F_i
real (kind=sedc_rk),dimension(3)             :: e_ij = 0.0_sedc_rk, r_ij = 0.0_sedc_rk
real (kind=sedc_rk)                          :: r=0,cf=0,cdf_dr=0,r_6=0,r_12=0

! Form the vector r_ij pointing at ion j from ion i
r_ij = j_coords - internal_cart_coord(:,ion_i)

! The norm of this vector
r = sqrt(dot_product(r_ij,r_ij))

! Form unit vector e_ij from r_ij
e_ij = r_ij / r

! Compute c6_ij * f(r) and c6_ij * df(r)/dr
!!! Mitch
!call c_f(ion_i,ion_j,r,cf,cdf_dr)
call scheme_selector(ion_i,ion_j,r,cf,cdf_dr)

! Energy:
call sedc_skip_group_ff(ion_i, ion_j)
! TODO: [ ] This is a debug statement to check for correct group-skipping
!if (sedc_skip_group == 0) then 
!write(*,*) "debug: ",sedc_skip_group," ion_i: ",ion_i," ion_j: ",ion_j
!endif

r_6 = r ** (-6); sedc_energy = sedc_energy + cf * r_6 * sedc_skip_group

! F = -nabla E ==>
! F_j = (-( c_6_ij * (df(r)/dr * r^-6 - 6 * f(r) * r^-7)) * unit vector e_ij (i --> j)
! F_i = -F_j
call sedc_skip_group_ff(ion_i, ion_j)
F_i = r_6 * (cdf_dr - (6.0_sedc_rk/r) * cf) * e_ij * sedc_skip_group

if (sedc_do_standalone) then
  r_12 = r ** (-12); sedc_energy = sedc_energy + 8.0_sedc_rk * cf * cf * r_12 * sedc_skip_group
  F_i = F_i + 8.0_sedc_rk * r_12 * (cdf_dr - (12.0_sedc_rk/r) * cf * cf) * e_ij * sedc_skip_group
endif

! Update analytical forces for ion i only, 
! return F_i as output for use elsewhere
sedc_forces(:,ion_i) = sedc_forces(:,ion_i) + F_i 

  sedc_stress(1) = sedc_stress(1) + F_i(1) * r_ij(1)
  sedc_stress(2) = sedc_stress(2) + F_i(2) * r_ij(2)
  sedc_stress(3) = sedc_stress(3) + F_i(3) * r_ij(3)
!   sedc_stress(4) = sedc_stress(4) + F_i(2) * r_ij(3)
!   sedc_stress(5) = sedc_stress(5) + F_i(3) * r_ij(1)
!   sedc_stress(6) = sedc_stress(6) + F_i(1) * r_ij(2)
  ! this should average out possible summation errors...
  sedc_stress(4) = sedc_stress(4) + 0.5_sedc_rk * ( F_i(2) * r_ij(3) + F_i(3) * r_ij(2) )
  sedc_stress(5) = sedc_stress(5) + 0.5_sedc_rk * ( F_i(3) * r_ij(1) + F_i(1) * r_ij(3) )
  sedc_stress(6) = sedc_stress(6) + 0.5_sedc_rk * ( F_i(1) * r_ij(2) + F_i(2) * r_ij(1) )

end subroutine
!!! end: stress tensor


subroutine sedc_central_num_F(c_f,dopbc)
implicit none
character(len=*),intent(in)                       :: c_f
logical,intent(in)                             :: dopbc
real (kind=sedc_rk),dimension(:,:),allocatable :: num_forces
integer                                        :: i = 0,j = 0
real (kind=sedc_rk)                            :: tmp_coord,dE

allocate(num_forces(3,sedc_n_ions)); num_forces = 0.0_sedc_rk

if (dopbc) then
    do i = 1,3
        do j = 1,sedc_n_ions
        
            dE = 0.0_sedc_rk
            tmp_coord = internal_cart_coord(i,j)
           
           ! TODO: [ ] from fct-conv, check 
            call E_num_F_pbc_ff(c_f)

            internal_cart_coord(i,j) = tmp_coord + &
            & num_F_coord_shift; dE  = E_num_F_pbc
            
            internal_cart_coord(i,j) = tmp_coord - &
            & num_F_coord_shift; dE  = dE - E_num_F_pbc
            
            internal_cart_coord(i,j) = tmp_coord

            ! _Subtract_ the gradient contribution to form the force
            num_forces(i,j) = num_forces(i,j) &
            & - dE / (2.0_sedc_rk * num_F_coord_shift)

        enddo
    enddo
else
    do i = 1,3
        do j = 1,sedc_n_ions
        
            dE = 0.0_sedc_rk
            tmp_coord = internal_cart_coord(i,j)
           
            ! TODO: [ ] from fct-conv, check 
            call E_num_F_ff(c_f)
            
            internal_cart_coord(i,j) = tmp_coord + &
            & num_F_coord_shift; dE  = E_num_F
            
            internal_cart_coord(i,j) = tmp_coord - &
            & num_F_coord_shift; dE  = dE - E_num_F            
            internal_cart_coord(i,j) = tmp_coord

            ! _Subtract_ the gradient contribution to form the force
            num_forces(i,j) = num_forces(i,j) &
            & - dE / (2.0_sedc_rk * num_F_coord_shift)

        enddo
    enddo
endif

! TODO: [ ] This function does indeed do some final calculations, therefore we will need it!
!           This is not a very tidy solution, maybe change the subroutine's behaviour
call sedc_num_F_prettyprint(sedc_n_ions,num_forces)

deallocate(num_forces)

end subroutine

subroutine sedc_num_F_prettyprint(array_size,num_forces)
implicit none
integer,intent(in)                            :: array_size
real (kind=sedc_rk),dimension(:,:),intent(in) :: num_forces
integer                                       :: i=0
real (kind=sedc_rk),dimension(3,array_size)   :: F_minus_num_F
real (kind=sedc_rk),dimension(array_size)     :: delta_F_abs_vector,F_num_abs_vector,ratio_vector
real (kind=sedc_rk)                           :: delta_F_max_abs=0,F_num_max_abs=0,max_ratio=0
character(len=5),dimension(array_size)        :: labels_1
character(len=16),dimension(3),parameter      :: labels_2 = (/ 'delta Fx','delta Fy','delta Fz' /)
!~ character(len=15),dimension(3),parameter      :: F_labels = (/'Num Fx','Num Fy','Num Fz'/)

F_minus_num_F    = 0.0_sedc_rk; delta_F_abs_vector = 0.0_sedc_rk
F_num_abs_vector = 0.0_sedc_rk;       ratio_vector = 0.0_sedc_rk

! Row labels for matrix printing (species of ions)
labels_1 = (/ (periodic_table(sedc_species(i)),i = 1,sedc_n_ions) /)

! Form the analytical - numerical forces difference
F_minus_num_F   = sedc_forces_store - num_forces
call F_abs_max_ff(num_forces,sedc_n_ions)
F_num_max_abs   = F_abs_max
call F_abs_max_ff(F_minus_num_F,sedc_n_ions)
delta_F_max_abs = F_abs_max

! Also compute the max. ratio between delta_F and num_F
delta_F_abs_vector = sqrt(sum(F_minus_num_F * F_minus_num_F,dim=1))
F_num_abs_vector   = sqrt(sum(num_forces * num_forces,dim=1))
! Protect 'max_ratio' from 'Infty' values - the numerical gradient _can_ be 0
where (F_num_abs_vector /= 0.0_sedc_rk) &
& ratio_vector     = delta_F_abs_vector / F_num_abs_vector
max_ratio          = maxval(ratio_vector)

!~ call print_header(3,1,'SEDC Numerical Forces')
!~ call print_matrix('Atom',sedc_n_ions,3,num_forces,labels_1,F_labels,'ES',17,5)

! TODO: [ ] these print statements should most likely be removed after final implementation,
!           however we should keep them for debugging.
!       [ ] Make the printed variables global so that they can be accessed from ASE if needed
call print_header(3,1,'SEDC Numerical Force Check')

write (sedc_stdout_unit,"(A,/,5(A),/,6(A),ES9.3,A,/,A3)" ) lstart, &
& lstart,tab,'Analytical - numerical F difference',tab,': delta F = F - F_num', &
& lstart,tab,'Coordinate displacement used',tab,tab,': ',num_F_coord_shift,' [A]', &
& lstart

write (sedc_stdout_unit,"(2(7(A),ES9.3,A,/),6(A),ES9.3,A)") & 
& lstart,tab,tab,'Per atom |F_num|_max',  tab,tab,'= ',F_num_max_abs,' [eV/A] ', &
& lstart,tab,tab,'Per atom |delta F|_max',tab,tab,'= ',delta_F_max_abs,' [eV/A] ', &
& lstart,tab,tab,'(|delta F| / |F_num|)_max', tab,'= ',max_ratio,' [eV/A] '

! TODO: [ ] This can be safely omitted from the initial port since it only contains printing statements!
call print_matrix('Atom',3,F_minus_num_F,labels_1,labels_2,'ES',13,5)

end subroutine


subroutine sedc_num_stress(c_f)
implicit none
character(len=*),intent(in)                                  :: c_f
! quantities which are changed when calculating the lattice distorted cells with the 'main' routines
real (kind=sedc_rk), dimension(:,:), allocatable :: sedc_cart_coord0
real (kind=sedc_rk), dimension(:,:), allocatable :: cell_vectors0
real (kind=sedc_rk) :: central_cell_volume0

real (kind=sedc_rk), parameter :: zero = 0.0_sedc_rk
real (kind=sedc_rk), parameter :: one = 1.0_sedc_rk
real (kind=sedc_rk), dimension(3,3), parameter :: &
&    unity = reshape( (/ one, zero, zero,  zero, one, zero,  zero, zero, one /), (/ 3,3 /) )
! real (kind=sedc_rk), dimension(3,3), parameter :: &
! &    unity = reshape( &
! &    (/ 1.0_sedc_rk, 0.0_sedc_rk, 0.0_sedc_rk, &
! &       0.0_sedc_rk, 1.0_sedc_rk, 0.0_sedc_rk, &
! &       0.0_sedc_rk, 1.0_sedc_rk, 1.0_sedc_rk /), (/ 3,3 /) )
real (kind=sedc_rk), parameter :: eps = num_stress_eps
! mapping of 6D stress array 'stress_6' to 3x3 stress matrix 'stress_33':
!	stress_6(1) = stress_33(1,1) i.e. Sxx
!	stress_6(2) = stress_33(2,2) i.e. Syy
!	stress_6(3) = stress_33(3,3) i.e. Szz
!	stress_6(4) = stress_33(2,3) i.e. Syz and Szy
!	stress_6(5) = stress_33(1,3) i.e. Sxz and Szx
!	stress_6(6) = stress_33(1,2) i.e. Sxy and Syx
integer, dimension(2,6), parameter :: &
&    stress_6to33 = reshape( (/ 1,1, 2,2, 3,3, 2,3, 1,3, 1,2 /), (/ 2,6 /) )

real (kind=sedc_rk),dimension(:,:), allocatable :: inv_cell
real (kind=sedc_rk),dimension(:,:), allocatable :: fract_coord
real (kind=sedc_rk),dimension(:,:), allocatable :: cell_dist
real (kind=sedc_rk) :: dE
real (kind=sedc_rk),dimension(:), allocatable :: num_stress
integer :: i,j,n

! store quantities which are changed when calculating the lattice distorted cells with the 'main' routines
allocate(sedc_cart_coord0(3,sedc_n_ions)); sedc_cart_coord0 = sedc_cart_coord
allocate(cell_vectors0(3,3)); cell_vectors0 = sedc_cell_vectors
central_cell_volume0 = central_cell_volume


! allocate work arrays
allocate(inv_cell(3,3))
call inv_matrix_ff(cell_vectors0)
inv_cell = inv_matrix
allocate(fract_coord(3,sedc_n_ions))
do i = 1,sedc_n_ions
   fract_coord(:,i) = matmul(inv_cell,sedc_cart_coord(:,i))
end do
allocate(cell_dist(3,3))
allocate(num_stress(6)); num_stress = 0.0_sedc_rk

call print_header(3,1,'SEDC Numerical Stress Check')

do n = 1,6

    dE = 0.0_sedc_rk

    write (sedc_stdout_unit,"(A,/,3(A),I1,A)") &
    &        lstart,lstart,tab,'Doing stress displacement ', n, ' (+)'
    cell_dist = unity
    i = stress_6to33(1,n); j = stress_6to33(2,n)
    cell_dist(i,j) = cell_dist(i,j) + eps
    sedc_cell_vectors = matmul(cell_dist,cell_vectors0)
    do i = 1,sedc_n_ions
      sedc_cart_coord(:,i) = matmul(sedc_cell_vectors,fract_coord(:,i))
    end do
    ! using dummy arguments (not required anyway - since sedc_initialized.eq.true should always hold here)
    call reinitialize_main()
    call sedc_central_E_F_stress(c_f)
    if (sedc_do_PBC) call sedc_pbc(c_f,.false.,num_F_coord_shift)
    dE = dE + (sedc_energy - sedc_energy_store)

    write (sedc_stdout_unit,"(A,/,3(A),I1,A)") &
    &        lstart,lstart,tab,'Doing stress displacement ', n, ' (-)'
    cell_dist = unity
    i = stress_6to33(1,n); j = stress_6to33(2,n)
    cell_dist(i,j) = cell_dist(i,j) - eps
    sedc_cell_vectors = matmul(cell_dist,cell_vectors0)
    do i = 1,sedc_n_ions
        sedc_cart_coord(:,i) = matmul(sedc_cell_vectors,fract_coord(:,i))
    end do
    ! using dummy arguments (not required anyway - since sedc_initialized.eq.true should always hold here)
    call reinitialize_main()
    call sedc_central_E_F_stress(c_f)
    if (sedc_do_PBC) call sedc_pbc(c_f,.false.,num_F_coord_shift)
    dE = dE - (sedc_energy - sedc_energy_store)

    num_stress(n) = dE / (2.0_sedc_rk * eps)

    ! abuse global stress array - since it contains bogus anyway at this point...
    sedc_stress = 0.0_sedc_rk
    sedc_stress(n) = num_stress(n) / central_cell_volume0

end do

num_stress = num_stress / central_cell_volume

! output final comparison
write (sedc_stdout_unit,"(A,/,3(A))") lstart,lstart,tab,'Numerical stress'
call sedc_print_stress(num_stress)
write (sedc_stdout_unit,"(A,/,3(A))") lstart,lstart,tab,'Analytical - Numerical stress'
call sedc_print_stress(sedc_stress_store - num_stress)

! deallocate work arrays
deallocate(num_stress)
deallocate(cell_dist)
deallocate(fract_coord)
deallocate(inv_cell)

! restore quantities from above
sedc_cart_coord = sedc_cart_coord0; deallocate(sedc_cart_coord0)
sedc_cell_vectors = cell_vectors0; deallocate(cell_vectors0)

end subroutine

! print routine for stress tensor only 
! (-> could also be called from sedc_print_results to avoid code duplication?!)
subroutine sedc_print_stress(stress)
implicit none
real (kind=sedc_rk),dimension(:),intent(in) :: stress
character(len=12),dimension(3),parameter    :: &
& s_labels_1 = (/ '[eV/A^3]  Sx','          Sy','          Sz' /),&
& s_labels_2 = (/ '          Sx','          Sy','          Sz' /)
real (kind=sedc_rk),dimension(3,3)          :: s_tmp

! Let's print the stress tensor only if the cell matrix is invertible
! (Stress tensor is a NaN otherwise)
if (central_cell_volume > 0) then
   ! Stress tensor printed as:
   ! S{xx,yy,zz} = sedc_stress(1:3)
   ! S{xy,yx} = sedc_stress(6)
   ! S{xz,zx} = sedc_stress(5)
   ! S{yz,zy} = sedc_stress(4)
   s_tmp = 0.0_sedc_rk;         s_tmp(1,1) = stress(1)
   s_tmp(2,2) = stress(2); s_tmp(3,3) = stress(3)
   s_tmp(2,3) = stress(4); s_tmp(3,2) = stress(4)
   s_tmp(3,1) = stress(5); s_tmp(1,3) = stress(5)
   s_tmp(2,1) = stress(6); s_tmp(1,2) = stress(6)
   call print_matrix('    Stress tensor',3,s_tmp,s_labels_2,s_labels_1,'ES',12,4)
endif

end subroutine

!!! end: stress tensor

subroutine group_atoms()
implicit none
character(len=max_char_len)                    :: err_msg = '',a = '', b = ''
integer                                        :: tot_n_grouped_ions = 0
integer,dimension(:),allocatable               :: tmp_pbc_groups
logical,dimension(:),allocatable               :: tmp_pbc_g_fold
integer,dimension(:,:),allocatable             :: tmp_pbc_g_switches
real (kind=sedc_rk),dimension(:,:),allocatable :: tmp_pbc_g_cells
logical                                        :: pbc_settings,sw_alloc

pbc_settings = .false.; sw_alloc = allocated(sedc_pbc_g_switches)

!mijp: cannot do this as don't know order of execution of clauses!
!pbc_settings = &
!& allocated(sedc_groups) .and. &
!& allocated(sedc_pbc_g_fold) .and. sw_alloc .and. &
!& allocated(sedc_pbc_g_cells) .and. &
!& size(sedc_groups) == sedc_n_groups .and. &
!& size(sedc_pbc_g_fold) == sedc_n_groups .and. &
!& size(sedc_pbc_g_switches) == 6 * sedc_n_groups .and. &
!& size(sedc_pbc_g_cells) == 9 * sedc_n_groups
pbc_settings = &
& allocated(sedc_groups) .and. &
& allocated(sedc_pbc_g_fold) .and. sw_alloc .and. &
& allocated(sedc_pbc_g_cells)
if (pbc_settings) pbc_settings=pbc_settings .and. &
& size(sedc_groups) == sedc_n_groups .and. &
& size(sedc_pbc_g_fold) == sedc_n_groups .and. &
& size(sedc_pbc_g_switches) == 6 * sedc_n_groups .and. &
& size(sedc_pbc_g_cells) == 9 * sedc_n_groups

if (.not.pbc_settings) then
    if (sedc_n_groups > 0) &
        & call sedc_errquit('SEDC PBC Initialization Warning:', &
        & 'One or more PBC setting arrays not properly allocated - PBC reset.',.true.)
    if (sw_alloc)                    deallocate(sedc_pbc_g_switches)
    if (allocated(sedc_groups))      deallocate(sedc_groups)
    if (allocated(sedc_pbc_g_fold))  deallocate(sedc_pbc_g_fold)
    if (allocated(sedc_pbc_g_cells)) deallocate(sedc_pbc_g_cells)
    
    if (allocated(sedc_pbc_g_skip))  deallocate(sedc_pbc_g_skip)
    if (allocated(sedc_pbc_g_only_intra) ) deallocate(sedc_pbc_g_only_intra)
    
    sedc_n_groups = 0
    allocate(sedc_groups(sedc_n_groups)); sedc_groups = 0
    allocate(sedc_pbc_g_fold(sedc_n_groups))
    allocate(sedc_pbc_g_skip(sedc_n_groups))
    allocate(sedc_pbc_g_only_intra(sedc_n_groups))
    allocate(sedc_pbc_g_switches(6,sedc_n_groups))
    allocate(sedc_pbc_g_cells(3,sedc_n_groups))
endif

tot_n_grouped_ions = sum(sedc_groups)

!If we've somehow grouped more ions than exist, that's not good
if (tot_n_grouped_ions > sedc_n_ions) then
    write (a,"(I5)") tot_n_grouped_ions; write (b,"(I5)") sedc_n_ions
    write (err_msg,"(5(A))") 'Number of atoms in PBC control groups (= '&
    & ,trim(adjustl(a)),') > total number of atoms (= ',trim(adjustl(b)),')'
    call sedc_errquit('SEDC PBC Initialization ERROR!',err_msg)
    return

else if (tot_n_grouped_ions < sedc_n_ions) then

    allocate(tmp_pbc_groups(sedc_n_groups))
    tmp_pbc_groups = sedc_groups
    deallocate(sedc_groups)
    allocate(sedc_groups(sedc_n_groups + 1)); sedc_groups = 0
    sedc_groups(1:sedc_n_groups) = tmp_pbc_groups
    deallocate(tmp_pbc_groups)

    sedc_groups(sedc_n_groups + 1) = sedc_n_ions - tot_n_grouped_ions

    allocate(tmp_pbc_g_cells(3,3*sedc_n_groups))
    tmp_pbc_g_cells = sedc_pbc_g_cells
    deallocate(sedc_pbc_g_cells)
    allocate(sedc_pbc_g_cells(3,3*(sedc_n_groups+1)))
    sedc_pbc_g_cells = 0.0_sedc_rk
    sedc_pbc_g_cells(:,1:(3*sedc_n_groups)) = tmp_pbc_g_cells
    deallocate(tmp_pbc_g_cells)

    sedc_pbc_g_cells(:,(3*sedc_n_groups + 1):(3*sedc_n_groups + 3)) = &
    & sedc_cell_vectors
        
    allocate(tmp_pbc_g_fold(sedc_n_groups))
    tmp_pbc_g_fold = sedc_pbc_g_fold
    deallocate(sedc_pbc_g_fold)
    allocate(sedc_pbc_g_fold(sedc_n_groups + 1))
    sedc_pbc_g_fold = sedc_pbc_backfold_coord
    sedc_pbc_g_fold(1:sedc_n_groups) = tmp_pbc_g_fold
    deallocate(tmp_pbc_g_fold)
!REINI
    allocate(tmp_pbc_g_fold(sedc_n_groups))
    tmp_pbc_g_fold = sedc_pbc_g_skip
    deallocate(sedc_pbc_g_skip)
    allocate(sedc_pbc_g_skip(sedc_n_groups + 1))
    sedc_pbc_g_skip = .false.
    sedc_pbc_g_skip(1:sedc_n_groups) = tmp_pbc_g_fold
    deallocate(tmp_pbc_g_fold)

    allocate(tmp_pbc_g_fold(sedc_n_groups))
    tmp_pbc_g_fold = sedc_pbc_g_only_intra
    deallocate(sedc_pbc_g_only_intra)
    allocate(sedc_pbc_g_only_intra(sedc_n_groups + 1))
    sedc_pbc_g_only_intra = .false.
    sedc_pbc_g_only_intra(1:sedc_n_groups) = tmp_pbc_g_fold
    deallocate(tmp_pbc_g_fold)

    allocate(tmp_pbc_g_switches(6,sedc_n_groups))
    tmp_pbc_g_switches = sedc_pbc_g_switches
    deallocate(sedc_pbc_g_switches)
    allocate(sedc_pbc_g_switches(6,(sedc_n_groups + 1)))
    sedc_pbc_g_switches = -1
    sedc_pbc_g_switches(:,1:sedc_n_groups) = tmp_pbc_g_switches
    deallocate(tmp_pbc_g_switches)
        
    sedc_n_groups = sedc_n_groups + 1
endif

end subroutine


subroutine sedc_init_pbc()
implicit none
integer                                        :: i = 0,k = 0,n = 0,lc = 0,cc = 0,rc = 0
character(len=1),dimension(3)                  :: switch_array
integer,dimension(3)                           :: indices = 0
integer,dimension(:,:),allocatable             :: tmp_pbc_g_switches
real (kind=sedc_rk),dimension(3,3)             :: inv_cell = 0.0_sedc_rk,tmp_matrix = 0.0_sedc_rk

if (.not.sedc_is_initialized) then
    
    if (sedc_pbc_img_fixed_nshells .and. sedc_pbc_n_shells < 1) then
        sedc_do_PBC = .false.
        return
    endif
    
    allocate(tmp_pbc_g_switches(6,sedc_n_groups))
    tmp_pbc_g_switches = 0; switch_array = ''
    read (sedc_pbc_switch,"(3(A1))") switch_array(1:3)
    do i = 1,3
        select case (switch_array(i))
            case("A")
                tmp_pbc_g_switches(1,:) = 1
                tmp_pbc_g_switches(4,:) = 1
            case("B")
                tmp_pbc_g_switches(2,:) = 1
                tmp_pbc_g_switches(5,:) = 1
            case("C")
                tmp_pbc_g_switches(3,:) = 1
                tmp_pbc_g_switches(6,:) = 1
            case("")
            case default
                call sedc_errquit('SEDC PBC Initialization ERROR!',&
                & "PBC switch parameter value '" &
                & // trim(sedc_pbc_switch) // "' not understood")
                return
        end select
    enddo
    
    if (all(mask= sedc_pbc_g_switches == 0)) then
        sedc_do_PBC = .false.
        return
    endif
    
    where (sedc_pbc_g_switches * sedc_pbc_g_switches /= sedc_pbc_g_switches) &
    & sedc_pbc_g_switches = tmp_pbc_g_switches
    deallocate(tmp_pbc_g_switches)
    
    do i = 1,sedc_n_groups
        do k = 1,3
!            if (all(mask= sedc_pbc_g_cells(:,(3*(i-1) + k)) == 0)) &
!            & sedc_pbc_g_cells(:,(3*(i-1) + k)) = sedc_cell_vectors(:,k)
! jm 20101108: 
! Hack (aka hot fix) to cope with sedc_is_initialized effectively being disabled.
! Almost destroys grouping functionality as originally designed
! (-> e.g. for different PBCs of adsorbate und substrate),
! but due to sedc_parse_input currently commented out in dftd.F90,
! no user way to input groups anyway (i.e. functionality not exposed).
            sedc_pbc_g_cells(:,(3*(i-1) + k)) = sedc_cell_vectors(:,k)
        enddo
    enddo

    if (sedc_variable_cell) then
        if (allocated(gcell_x_inv_cellvecs)) deallocate(gcell_x_inv_cellvecs) !mijp
        allocate(gcell_x_inv_cellvecs(3,3 * sedc_n_groups))
        call inv_matrix_ff(sedc_cell_vectors)
        inv_cell = inv_matrix
        if (sedc_exec_status /= 0) return
        
        do i = 1,sedc_n_groups
            tmp_matrix = sedc_pbc_g_cells(:,(3 * i - 2):(3 * i))
            gcell_x_inv_cellvecs(:,(3 * i - 2):(3 * i)) = &
            & matmul(inv_cell,tmp_matrix)
        enddo
    endif

    if (allocated(pbc_g_int_scaling)) deallocate(pbc_g_int_scaling) !mijp
    allocate(pbc_g_int_scaling(sedc_n_groups)); pbc_g_int_scaling = 2

    do i = 1,sedc_n_groups
        if (all(mask= sedc_pbc_g_switches(:,i) == 0)) then
            pbc_g_int_scaling(i) = 1
            sedc_pbc_g_fold(i) = .false.
        endif
    enddo

    sedc_pbc_energy_tol = abs(sedc_pbc_energy_tol)
    sedc_pbc_force_tol  = abs(sedc_pbc_force_tol)

    ! The following algorithm helps to determine the exact number of 
    ! atoms / cells per image shell
    pbc_ion_in_s = 0; pbc_cell_in_s = 0

    do k = 1,sedc_n_groups

        indices = (/ 3,1,2 /)
        do i = 1,3
            ! left column, center column and right column indices
            lc = indices(1); cc = indices(2); rc = indices(3)

            n = (sedc_pbc_g_switches(cc,k) + sedc_pbc_g_switches((cc+3),k))

            pbc_cell_in_s(1) = pbc_cell_in_s(1) + n
            pbc_ion_in_s(1)  = pbc_ion_in_s(1)  + n * sedc_groups(k)

            n = ( sedc_pbc_g_switches(rc,k)     * sedc_pbc_g_switches(lc,k) + &
                & sedc_pbc_g_switches(rc,k)     * sedc_pbc_g_switches((lc+3),k) + &
                & sedc_pbc_g_switches((rc+3),k) * sedc_pbc_g_switches(lc,k) + &
                & sedc_pbc_g_switches((rc+3),k) * sedc_pbc_g_switches((lc+3),k))

            pbc_cell_in_s(2) = pbc_cell_in_s(2) + n
            pbc_ion_in_s(2)  = pbc_ion_in_s(2)  + n * sedc_groups(k)

            n = ( sedc_pbc_g_switches((lc+3),k) * sedc_pbc_g_switches(cc,k) * &
                & sedc_pbc_g_switches(rc,k)     + sedc_pbc_g_switches(lc,k) * &
                & sedc_pbc_g_switches((cc+3),k) * sedc_pbc_g_switches((rc+3),k))

            pbc_cell_in_s(3) = pbc_cell_in_s(3) + n
            pbc_ion_in_s(3) = pbc_ion_in_s(3)   + n * sedc_groups(k)

            indices = cshift(indices,1)
        enddo

        n = ( sedc_pbc_g_switches(1,k) * sedc_pbc_g_switches(2,k) * &
            & sedc_pbc_g_switches(3,k) + sedc_pbc_g_switches(4,k) * &
            & sedc_pbc_g_switches(5,k) * sedc_pbc_g_switches(6,k))

        pbc_cell_in_s(3) = pbc_cell_in_s(3) + n
        pbc_ion_in_s(3)  = pbc_ion_in_s(3)  + n * sedc_groups(k)
    enddo
    
else if (sedc_variable_cell) then

    do i = 1,sedc_n_groups
        tmp_matrix = gcell_x_inv_cellvecs(:,(3 * i - 2):(3 * i))
        sedc_pbc_g_cells(:,(3 * i - 2):(3 * i)) = &
        & matmul(tmp_matrix,sedc_cell_vectors)
    enddo
        
    ! TODO: [ ] This is only a stdout-debug statement, removed.        
    if (sedc_print_level == debug_prn_lvl) &
    & call print_pbc_group_cells(sedc_n_groups)
    
endif

end subroutine

! TODO: [ ] purely for debugging!! remove later
subroutine print_pbc_group_cells(n_pbc_groups)
implicit none
integer,intent(in)                            :: n_pbc_groups
integer                                       :: i,j,k
character(len=max_char_len)                   :: a,b
character(len=6),dimension(3),parameter       :: abc = (/ 'A vec.','B vec.','C vec.' /)
character(len=1),dimension(3),parameter       :: g_cell_col_labels = (/ 'X','Y','Z' /)
character(len=16),dimension(3 * n_pbc_groups) :: g_cell_row_labels

if (sedc_variable_cell) write (sedc_stdout_unit,"(A,/,3(A))") &
& lstart,lstart,tab,'Variable cell calculation requested.'

k = 0; a = ''
do i = 1,n_pbc_groups
    do j = 1,3
        k = k + 1; write (b,"(I4)") i
        write (a,"(A,1X,2(A),2X,A)") 'G',trim(adjustl(b)),';',abc(j)
        g_cell_row_labels(k) = a
    enddo
enddo

call print_matrix(' PBC Group Cells:',3,sedc_pbc_g_cells,g_cell_row_labels,g_cell_col_labels,'F',12,6)

end subroutine

subroutine sedc_pbc(c_f,prn,etol)
implicit none
character(len=*),intent(in)                       :: c_f
logical,intent(in)                             :: prn
real (kind=sedc_rk),intent(in)                 :: etol
real (kind=sedc_rk)                            :: shell_dE = 0, shell_dF = 0
real (kind=sedc_rk),dimension(:,:),allocatable :: tmp_forces
character(len=max_char_len)                    :: a = ''

pbc_shells_counted = 0; shell_dE = 1e2; shell_dF = 1e2

! TODO: [ ] This seems to be printing options, delete after final implementation
! TODO: [ ] Delete associated variables
if (sedc_print_xyz) then
    xyz_comment   = trim(xyz_comment)   // ' + per. image'
    xyz_F_comment = trim(xyz_F_comment) // ' + per. image'
endif

if (prn) then
    write (a,"(3(A,5X),A)") &
    & 'Shell','Total E corr.','|dE| / atom','|dF|max / atom'
    call print_header(3,1,'SEDC PBC Interaction Energy')
    call print_header(1,2,trim(a))
endif
! END TODO

allocate(tmp_forces(3,sedc_n_ions)); tmp_forces = 0.0_sedc_rk

do while ((sedc_pbc_img_fixed_nshells  .and. &
    & pbc_shells_counted < sedc_pbc_n_shells) .or. &
    & (.not.sedc_pbc_img_fixed_nshells .and. &
    & (shell_dE > etol .or. shell_dF > sedc_pbc_force_tol)))

    pbc_shells_counted = pbc_shells_counted + 1
    
    shell_dE = sedc_energy; tmp_forces = sedc_forces
    
    ! TODO: n_cells_in_this_shell: func->subrout conversion!
    call n_cells_in_this_shell_ff(pbc_shells_counted)
    call sedc_pbc_compute_shell(c_f,pbc_shells_counted,&
    & n_cells_in_this_shell)
    
    tmp_forces = sedc_forces - tmp_forces

    ! TODO: F_abs_max: func->subrout conversion!
    call F_abs_max_ff(tmp_forces,sedc_n_ions)
    shell_dF = F_abs_max
    shell_dE = abs(sedc_energy - shell_dE) / sedc_n_ions
        
    ! TODO: [ ] This seems to be a printing statement, remove with final commit
    if (prn) write (sedc_stdout_unit,"(A3,A1,I5,F17.8,2(ES17.6))") &
        & lstart,tab,pbc_shells_counted,sedc_energy,shell_dE,shell_dF
enddo

deallocate(tmp_forces)

end subroutine

! Computes interaction with one single shell
subroutine sedc_pbc_compute_shell(c_f,shell,cells_in_shell)
implicit none
character(len=*),intent(in)                       :: c_f
integer,intent(in)                             :: shell,cells_in_shell
integer,dimension(:),allocatable               :: tmp_species_to_print,image_indices
real (kind=sedc_rk),dimension(:,:),allocatable :: tmp_print_coord,coord_array
real (kind=sedc_rk),dimension(6)               :: tmp_stress
integer,dimension(3,cells_in_shell)            :: cell_index_list
real (kind=sedc_rk),dimension(3)               :: offset_v
real (kind=sedc_rk)                            :: tmp_E,delta
integer                                        :: i,j,g,start_ion,n_in_g,n_cell_in_g,k
integer                                        :: n_counted,array_start,array_stop,ions_in_shell

n_counted = 0; array_start = 0; array_stop = 0

call n_ions_in_this_shell_ff(shell)
ions_in_shell = n_ions_in_this_shell

allocate(image_indices(ions_in_shell)); image_indices = 0
allocate(coord_array(3,ions_in_shell)); coord_array   = 0.0_sedc_rk

do g = 1,sedc_n_groups

    call image_list(shell,cells_in_shell,g,n_cell_in_g,cell_index_list)
   
    n_in_g    = sedc_groups(g)
    start_ion = n_counted + 1
    n_counted = n_counted + n_in_g
    
    do i = 1,n_cell_in_g
        array_start = array_stop + 1
        array_stop  = array_start + n_in_g - 1
        
        offset_v = matmul(sedc_pbc_g_cells(:,(3*g-2):(3*g)),cell_index_list(:,i))
        
        do j = 1,3; coord_array(j,array_start:array_stop) = offset_v(j); enddo
        
        coord_array(:,array_start:array_stop) = &
        & coord_array(:,array_start:array_stop) + internal_cart_coord(:,start_ion:n_counted)
        
        k = start_ion
        do j = array_start,array_stop
            image_indices(j) = k
            k = k + 1
        enddo

    enddo

enddo

n_counted = 0
do g = 1,sedc_n_groups

!REINI Adaptation
    if (sedc_pbc_g_skip(g)) then
        n_in_g = sedc_groups(g)
        n_counted = n_counted + n_in_g
        cycle
    end if

    n_in_g = sedc_groups(g)
    tmp_E  = sedc_energy
!!! begin: stress tensor
    tmp_stress = sedc_stress
!!! end: stress tensor

    !!!!$OMP PARALLEL
    call sedc_pbc_shell_group_int(c_f,n_counted + 1,n_counted + n_in_g,&
    & ions_in_shell,image_indices,coord_array)
    !!!!$OMP END PARALLEL
    
    ! Now rescale the PBC interaction energy to avoid double counting
    ! (and keep it where desirable)
    delta = (sedc_energy - tmp_E) / pbc_g_int_scaling(g)
    sedc_energy = tmp_E + delta
    group_E_decomposed(3,g) = group_E_decomposed(3,g) + delta

!!! begin: stress tensor
    ! keep in sync with energy expression
    sedc_stress = tmp_stress + (sedc_stress - tmp_stress) / pbc_g_int_scaling(g)
!!! end: stress tensor

    n_counted = n_counted + n_in_g
enddo

! If the image ion coordinates are to be printed, the print
! coordinate arrays need to be reallocated, extended and have shell image
! ion coordinates added
! TODO: [ ] another printing statement, keep for debug, remove for final commit
if (sedc_print_xyz) then
    allocate(tmp_print_coord(3,n_ions_to_print)); tmp_print_coord = 0.0_sedc_rk
    tmp_print_coord(:,1:n_ions_to_print) = &
    & ion_coordinates_to_print(:,1:n_ions_to_print)
    deallocate(ion_coordinates_to_print)
        
    allocate(ion_coordinates_to_print(3,(n_ions_to_print + ions_in_shell)))
    ion_coordinates_to_print = 0.0_sedc_rk
    ion_coordinates_to_print(:,1:n_ions_to_print) = &
    & tmp_print_coord(:,1:n_ions_to_print)
    deallocate(tmp_print_coord)
        
    ion_coordinates_to_print(:,(n_ions_to_print + 1):(n_ions_to_print + ions_in_shell)) = &
    & coord_array(:,1:ions_in_shell)
    
    allocate(tmp_species_to_print(n_ions_to_print)); tmp_species_to_print = 0
    tmp_species_to_print(1:n_ions_to_print) = &
    & species_to_print(1:n_ions_to_print)
    deallocate(species_to_print)
    
    allocate(species_to_print(n_ions_to_print + ions_in_shell)); species_to_print = 0
    species_to_print(1:n_ions_to_print) = &
    & tmp_species_to_print(1:n_ions_to_print)
    deallocate(tmp_species_to_print)
    
    species_to_print((n_ions_to_print + 1):(n_ions_to_print + ions_in_shell)) = &
    & (/ (sedc_species(image_indices(i)), i = 1,ions_in_shell) /)
    
    n_ions_to_print = n_ions_to_print + ions_in_shell
endif

deallocate(image_indices)
deallocate(coord_array)

end subroutine

subroutine sedc_pbc_shell_group_int(c_f,start_ion,stop_ion,n_in_s,image_indices,coords)
implicit none
character(len=*),intent(in)                            :: c_f
integer,intent(in)                                  :: start_ion,stop_ion,n_in_s
integer,dimension(n_in_s),intent(in)                :: image_indices
real (kind=sedc_rk),dimension(3,n_in_s),intent(in)  :: coords
real (kind=sedc_rk),dimension(3)                    :: F_i=0
integer                                             :: i=0,j=0

do i = start_ion,stop_ion
    do j = 1,n_in_s
        call pair_E_F_stress(c_f,i,image_indices(j),coords(:,j),F_i)
    enddo
enddo

end subroutine


subroutine image_list(shell,list_length,group,no_cells_counted,list)
implicit none
integer,intent(in)                           :: shell,list_length,group
integer,intent(out)                          :: no_cells_counted
integer,dimension(3,list_length),intent(out) :: list
integer,dimension(3)                         :: column_indices = 0
integer                                      :: p=0,n=0,cc=0,rc=0,lc=0,i=0,j=0,k=0,l=0
integer                                      :: j_sgn=0,edge_minus_one=0,j_sgn_sh=0,q=0
integer                                      :: j_stop=0,k_start=0,k_stop=0,l_start=0,l_stop=0

no_cells_counted = 0; list = 0

column_indices = (/ 3,1,2 /)

do i = 1,3

    lc = column_indices(1); cc = column_indices(2); rc = column_indices(3)

    if ((sedc_pbc_g_switches(cc,group) + &
    &    sedc_pbc_g_switches((cc+3),group)) /= 0) then

! Outer j loop only has two steps; +shell and -shell. If the positive
! periodicity switch P is 0 ('turned off'), the j loop must start at 
! -shell to construct the inner k,l loop boundaries correctly. Hence 'q' construct.
        q      = 1 - sedc_pbc_g_switches(cc,group) ! P = 0 ==> q = 1 and v.v.
        j      = 0 + shell * (sedc_pbc_g_switches(cc,group) - q)
        j_stop = 0 - shell * sedc_pbc_g_switches((cc+3),group)

        j_sgn  = sign(1,j)

        do while (j_sgn * j >= j_sgn * j_stop)
        
            j_sgn = sign(1,j)
            j_sgn_sh = j_sgn * shell

            p = (j_sgn * (j_sgn + 1))/2; n = 1 - p
            edge_minus_one = j_sgn_sh  - j_sgn

! k index loop boundaries: Includes edge on neg side, but not on pos side.
            k_start = 0 - j_sgn_sh * &
            & (sedc_pbc_g_switches(rc,group) * n + &
            &  sedc_pbc_g_switches((rc+3),group) * p)
            
            k_stop  = 0 + edge_minus_one * &
            & (sedc_pbc_g_switches(rc,group) * p + &
            &  sedc_pbc_g_switches((rc+3),group) * n)

! l index loop boundaries: Identical to k loop, but negative
            l_start = 0 + j_sgn_sh * &
            & (sedc_pbc_g_switches(lc,group) * p + &
            &  sedc_pbc_g_switches((lc+3),group) * n)

            l_stop  = 0 - edge_minus_one * &
            & (sedc_pbc_g_switches(lc,group) * n + &
            &  sedc_pbc_g_switches((lc+3),group) * p)
                    
            do k = k_start,k_stop,j_sgn
                do l = l_start,l_stop,-j_sgn
                    no_cells_counted = no_cells_counted + 1
                    list(rc,no_cells_counted) = k
                    list(cc,no_cells_counted) = j
                    list(lc,no_cells_counted) = l
                enddo
            enddo
            j = j - 2 * j_sgn_sh
        enddo
    endif
    column_indices = cshift(column_indices,1)
enddo

! The loop above misses two corners of the hypercube,
! namely (shell,shell,shell) and (-shell,-shell,-shell)
if (all(mask= sedc_pbc_g_switches(1:3,group) == 1)) then
    no_cells_counted = no_cells_counted + 1
    list(:,no_cells_counted) = shell
endif

if (all(mask= sedc_pbc_g_switches(4:6,group) == 1)) then
    no_cells_counted = no_cells_counted + 1
    list(:,no_cells_counted) = -shell
endif

end subroutine

! TODO: [ ] purely for debugging, remove later!!
subroutine print_vector(title,vector,row_labels,col_label,elem_fmt,w,d)
implicit none
real (kind=sedc_rk),dimension(:),intent(in) :: vector
integer,intent(in)                          :: w,d
character(len=*),dimension(:),intent(in)       :: row_labels
character(len=*),intent(in)                    :: title,col_label,elem_fmt
character(len=max_char_len)                 :: a
character(len=20)                           :: t_fmt,e_fmt
integer                                     :: l,m

a = ''; t_fmt = ''; e_fmt = ''

l = max(len(title),len(row_labels))
m = max(w,len(col_label)) + 1

write (t_fmt,"(2(A,I3.3),A)") '(A',l,',A',m,')'
write (a,t_fmt) title,col_label
write (e_fmt,"(A,I3.3,2(A),2(I3.3,A))") '(2(A),A',l,',',elem_fmt,m,'.',d,')'

m = size(vector)
call print_header(1,2, trim(a))
do l = 1,m
    write (sedc_stdout_unit,e_fmt) lstart,tab,adjustr(row_labels(l)),vector(l)
enddo

end subroutine

subroutine print_matrix(title,ncol,matrix,row_labels,col_labels,elem_fmt,w,d)
implicit none
real (kind=sedc_rk),dimension(:,:),intent(in) :: matrix
integer,intent(in)                            :: ncol,w,d
character(len=*),dimension(:),intent(in)         :: row_labels
character(len=*),dimension(:),intent(in)         :: col_labels
character(len=*),intent(in)                      :: title,elem_fmt
character(len=max_char_len)                   :: a,b
character(len=11)                             :: t_fmt,l_fmt,e_fmt
integer                                       :: i,j,k,l,m,nrow

a = ''; b = ''; t_fmt = ''; l_fmt = ''; e_fmt = ''
i = 0; j = 0; k = 0; l = 0
nrow = size(matrix,2)

write (t_fmt,"(A,I3.3,A)") 'A',max(len(title),len(row_labels)),')'
b = '(' // trim(t_fmt)
write (a,b) title
t_fmt = '(2(A),'// trim(t_fmt)

m = max(w,len(col_labels)) + 1
write (l_fmt,"(A,I3.3,A)") '(A,A',m,')'
write (e_fmt,"(2(A),2(I3.3,A))") '(',elem_fmt,m,'.',d,')'

m = size(matrix,1)

do while (l < m)
    k = min(m-l,ncol); b = ''

    do i = 1,k
        write (b,l_fmt) trim(b),adjustr(col_labels(l+i))
    enddo

    call print_header(1,2, trim(a) // trim(b))

    do i = 1,nrow
        write (sedc_stdout_unit,t_fmt,advance="no") &
        & lstart,tab,adjustr(row_labels(i))
        do j = 1,k
            write (sedc_stdout_unit,e_fmt,advance="no") matrix(l+j,i)
        enddo
        write (sedc_stdout_unit,"()")
    enddo

    l = l + ncol
enddo

end subroutine

subroutine print_header(ntabs,htype,header)
implicit none
integer,intent(in)       :: ntabs,htype
character(len=*),intent(in) :: header
integer                  :: i=0,l=0

if (htype /= 3 ) write(sedc_stdout_unit,"(A)") lstart

l = len_trim(header)
write (sedc_stdout_unit,"(A3)",advance="no") lstart
do i = 1,ntabs
    write (sedc_stdout_unit,"(A1)",advance="no") tab
enddo

if (htype == 2) then
    write (sedc_stdout_unit,"(A,/,A3)",advance="no") header,lstart
else
    write (sedc_stdout_unit,"(2X,A,/,A3)",advance="no") header,lstart
endif

do i = 1,ntabs
    write (sedc_stdout_unit,"(A1)",advance="no") tab
enddo

if (htype == 2) then
    do i = 1,l-1
        write (sedc_stdout_unit,"(A1)",advance="no") '-'
    enddo
    write (sedc_stdout_unit,"(A1)") '-'
else
    write (sedc_stdout_unit,"(A1)",advance="no") '+'
    do i = 1,l+2
        write (sedc_stdout_unit,"(A1)",advance="no") '-'
    enddo
    write (sedc_stdout_unit,"(A1)") '+'
endif

end subroutine


subroutine setting_summary()
implicit none
character(len=max_char_len) :: a,b,c,d
character(len=3)            :: e
integer                     :: i = 0, j = 0, k = 0, n = 0

if (sedc_print_level >= high_output_prn_lvl) &
    & call print_header(3,1,'SEDC Parameter / Setting Summary')

a = 'High'; if (sedc_print_level == debug_prn_lvl) a = 'Debug'
write (sedc_stdout_unit,"(A,/,9(A))") lstart, &
& lstart,tab,'Output print level',tab,tab,tab,tab,': ',trim(a)

a = 'No'; if (sedc_print_xyz) a = 'Yes'
write (sedc_stdout_unit,"(7(A))") lstart,tab, &
& 'Printing XYZ file of geometry / forces ',tab,tab,': ',trim(a)

if (sedc_do_PBC) then 
    write (sedc_stdout_unit,"(3(A))") &
    & lstart,tab,'Periodic Boundary Condition (PBC) interaction'

    if (sedc_pbc_img_fixed_nshells) then
        write (sedc_stdout_unit,"(6(A),I3)") &
        & lstart,tab,tab,'geometry fixed, no. image shells',&
        & tab,': ',sedc_pbc_n_shells
    else
        write (sedc_stdout_unit,"(5(A),ES7.1,A,/,5(A),ES7.1,A)") &
        & lstart,tab,'energy / atom conv. threshold; shell |dE|',&
        & tab,'   : ',sedc_pbc_energy_tol,' [eV]', &
        & lstart,tab,' force / atom conv. threshold; shell |dF|max',&
        & tab,': ',sedc_pbc_force_tol,' [eV/A]'
    endif

    write (a,"(A,5X,A,3X,3(A4,1X),3(2X,A))") &
    & "  PBC Groups:","Atoms in Group",'A','B','C',adjustr("Fold"), adjustr("Skip"), &
    & adjustr("Skip intra")  
    call print_header(1,2,a)
    
    j = 0
    do i = 1,sedc_n_groups
        n = sedc_groups(i)
        write (sedc_stdout_unit,"(2(A),I13)",advance="no") lstart,tab,i
        write (c,"(I5)") j + 1; write (d,"(I5)") n + j
        write (a,"(4A)") trim(adjustl(c)), &
        & "(",trim(adjustl(periodic_table(sedc_species(j+1)))),")"
        write (b,"(4A)") trim(adjustl(d)), &
        & "(",trim(adjustl(periodic_table(sedc_species(n+j)))),")"
        write (sedc_stdout_unit,"(1X,A8,A3,A8,2X)",advance="no") &
        & adjustr(trim(a))," - ",adjustr(trim(b))
            do k = 1,3
                e = ' 0 '
                if ( sedc_pbc_g_switches(k,i)   == 1 .and. &
                   & sedc_pbc_g_switches(k+3,i) == 1) then
                    e = '+/-'
                else
                    if (sedc_pbc_g_switches(k,i)   == 1) e = ' + '
                    if (sedc_pbc_g_switches(k+3,i) == 1) e = ' - '
                endif
                write (sedc_stdout_unit,"(A5)",advance="no") e
            enddo
        write (sedc_stdout_unit,"(2X,L2,4X,L2,6X,L2)") &
        &    sedc_pbc_g_fold(i), sedc_pbc_g_skip(i), sedc_pbc_g_only_intra(i) !REINI
        j = j + n
    enddo
    ! TODO: [ ] This is purely a debug statement, removed.
    if (sedc_print_level == debug_prn_lvl) &
        & call print_pbc_group_cells(sedc_n_groups)
endif

end subroutine

! TODO: [ ] This is used temporarily for debugging and might be adapted to write to a log later
subroutine sedc_print_results(scheme,xc,n_ions)
implicit none
character(len=*),intent(in)                   :: scheme,xc
integer,intent(in)                         :: n_ions
character(len=8),dimension(3),parameter    :: &
                                              & F_labels   = (/ '      Fx','      Fy','      Fz' /)
character(len=5),dimension(n_ions)         :: elem_labels
character(len=15),dimension(:),allocatable :: group_labels
integer                                    :: i = 0

if (sedc_print_level > 0) then
    if (sedc_print_level >= high_output_prn_lvl) then
    
        if (sedc_n_groups > 1) then
            allocate(group_labels(sedc_n_groups))
            
            do i = 1,sedc_n_groups
                write (group_labels(i),"(A,I3)") 'Gr. ',i
            enddo
            
            call print_matrix('Corr. E decompos.',3,group_E_decomposed,group_labels,&
            & (/'Intra Gr.','Inter Gr.','   PBC   '/),'F',12,8)
            
            deallocate(group_labels)
        endif
    
        if (sedc_print_level == debug_prn_lvl) then

            elem_labels = (/ (periodic_table(sedc_species(i)),i = 1,sedc_n_ions) /)
            call print_matrix('Force Corr [eV/A]',3,&
            & sedc_forces,elem_labels,F_labels,'ES',12,4)
            
            call sedc_print_stress(sedc_stress)
            
        endif
    endif
    
    write (sedc_stdout_unit,"(A)") lstart

    if (sedc_print_level >= high_output_prn_lvl) then
        write (sedc_stdout_unit,"(7(A),3X,A,F15.8,A,/,7(A),1X,A,F15.8,A)") &
        & lstart,tab,"'",scheme,"'/",xc," structure energy corr.   ", &
        & "=",central_cell_energy," [eV]", &
        & lstart,tab,"'",scheme,"'/",xc," PBC image interaction corr.", &
        & "=",sedc_energy - central_cell_energy," [eV]"
    else if (sedc_do_PBC) then
        write (sedc_stdout_unit,"(7(A),4X,A,I15)") &
        & lstart,tab,"'",scheme,"'/",xc," PBC image shells counted", &
        & "=",pbc_shells_counted
    endif        

    call F_abs_max_ff(sedc_forces,sedc_n_ions)
    write (sedc_stdout_unit,"(7(A),5X,A,F15.8,A,/,7(A),10X,A,F15.8,A)") &
    & lstart,tab,"'",scheme,"'/",xc," total energy correction","=",sedc_energy," [eV]", &
    & lstart,tab,"'",scheme,"'/",xc," correction |F|max ", &
    & "=",F_abs_max," [eV/A]"
endif

end subroutine

! TODO: [ ] keep this as a fallback in case we need to do testing or some quirks with
!           external programs which need an xyz-input, etc..
subroutine sedc_prn_xyz()
implicit none
integer                         :: i = 0,species_i = 0
character(len=130)              :: xyz_filename='',xyz_F_filename=''
character(len=30)               :: char_ions_to_print = '',a
character(len=5)                :: symbol=''
character(len=5)                :: sedc_file_bname='debug'
character(len=*),parameter         :: xyz_geom_file_ext = '_sedc_geom.xyz',&
                                 & xyz_F_file_ext = '_sedc_forces.xyz'
integer,parameter               :: xyz_file_unit = 315,&
                                 & xyz_F_file_unit = 416

! Output filenames from basename and extensions
xyz_filename = ''; xyz_F_filename = ''
if (single_point_number > 1 ) then
    write (a,"(I10)") single_point_number
    xyz_filename   = trim(sedc_file_bname) // '_step_' // &
                   & trim(adjustl(a)) // trim(xyz_geom_file_ext)
    xyz_filename   = trim(xyz_filename)

    xyz_F_filename = trim(sedc_file_bname) // '_step_' // &
                   & trim(adjustl(a)) // trim(xyz_F_file_ext)
    xyz_F_filename = trim(xyz_F_filename)
else
    xyz_filename   = trim(sedc_file_bname) // trim(xyz_geom_file_ext)
    xyz_filename   = trim(xyz_filename)

    xyz_F_filename = trim(sedc_file_bname) // trim(xyz_F_file_ext)
    xyz_F_filename = trim(xyz_F_filename)
endif

! Print the xyz file from the stored coordinates selected to print
write (char_ions_to_print, "(I30)") n_ions_to_print
xyz_comment = trim(xyz_comment) // ' atoms'

open(UNIT=xyz_file_unit, FILE=xyz_filename, STATUS="REPLACE",position="REWIND")

write (xyz_file_unit,"(A)") adjustl(trim(char_ions_to_print))
write (xyz_file_unit,"(A)") xyz_comment

do i = 1,n_ions_to_print
    species_i = species_to_print(i)
    symbol = trim(periodic_table(species_i))
    write (xyz_file_unit, "(A5,3(F18.9))") adjustl(symbol), &
    & ion_coordinates_to_print(:,i)
enddo

close(xyz_file_unit)

deallocate(species_to_print)
deallocate(ion_coordinates_to_print)
n_ions_to_print = 0

! Now print the coordinates and force vectors of the central cell
write (char_ions_to_print, "(I30)") sedc_n_ions
xyz_F_comment = trim(xyz_F_comment) // ' atoms'

open(UNIT=xyz_F_file_unit, FILE=xyz_F_filename, STATUS="REPLACE",position="REWIND")

write (xyz_F_file_unit,"(A)") adjustl(trim(char_ions_to_print))
write (xyz_F_file_unit,"(A)") xyz_F_comment

do i = 1,sedc_n_ions
    species_i = sedc_species(i)
    symbol = trim(periodic_table(species_i))
    write (xyz_F_file_unit, "(A5,6(F12.6))") adjustl(symbol), &
    & internal_cart_coord(:,i),sedc_forces(:,i)
enddo

close(xyz_F_file_unit)

xyz_comment = ''; xyz_F_comment = ''

end subroutine






! TODO: [ ] func->subrout ; check if works, remove after final commit
!function F_abs_max(force_matrix,n_force_vectors)
!implicit none
!real (kind=sedc_rk)                              :: F_abs_max
!real (kind=sedc_rk),dimension(:,:),intent(in)    :: force_matrix
!integer,intent(in)                               :: n_force_vectors
!real (kind=sedc_rk),dimension(3,n_force_vectors) :: tmp
!real (kind=sedc_rk),dimension(n_force_vectors)   :: F
!tmp = force_matrix * force_matrix
!F = sqrt(sum(tmp,dim=1)); F_abs_max = maxval(F)
!end function
subroutine F_abs_max_ff(force_matrix,n_force_vectors)
implicit none
real (kind=sedc_rk),dimension(:,:),intent(in)    :: force_matrix
integer,intent(in)                               :: n_force_vectors
real (kind=sedc_rk),dimension(3,n_force_vectors) :: tmp
real (kind=sedc_rk),dimension(n_force_vectors)   :: F
tmp = force_matrix * force_matrix
F = sqrt(sum(tmp,dim=1)); F_abs_max = maxval(F)
end subroutine 


! TODO: [ ] made a subroutine, check if its working!
!integer function n_cells_in_this_shell(shell)
!implicit none
!integer,intent(in) :: shell
!n_cells_in_this_shell = pbc_cell_in_s(1) + &
!& pbc_cell_in_s(2) * (2 * shell - 1) + &
!& pbc_cell_in_s(3) * (3 * shell * (shell - 1) + 1)
!end function



subroutine n_cells_in_this_shell_ff(shell)
implicit none
integer,intent(in) :: shell
n_cells_in_this_shell = pbc_cell_in_s(1) + &
& pbc_cell_in_s(2) * (2 * shell - 1) + &
& pbc_cell_in_s(3) * (3 * shell * (shell - 1) + 1)
end subroutine

! TODO: [ ] fct->subrout, check if its working!
subroutine n_ions_in_this_shell_ff(shell)
implicit none
integer,intent(in) :: shell
n_ions_in_this_shell = pbc_ion_in_s(1) + &
& pbc_ion_in_s(2) * (2 * shell - 1) + &
& pbc_ion_in_s(3) * (3 * shell * (shell - 1) + 1)
end subroutine 

! TODO: [ ] func->subrout conversion, check if it works
!function E_num_F(c_f)
!implicit none
!real (kind=sedc_rk)              :: E_num_F
!!external                         :: c_f
  !character(len=*),intent(in)                                  :: c_f
!real (kind=sedc_rk)              :: r,cf = 0.0_sedc_rk,cdf_dr = 0.0_sedc_rk
!real (kind=sedc_rk),dimension(3) :: r_ij
!integer                          :: j,k
!E_num_F = 0.0_sedc_rk
!do j = 2,sedc_n_ions
    !do k = 1,j-1
        !r_ij = internal_cart_coord(:,k) - internal_cart_coord(:,j)
        !r = sqrt(dot_product(r_ij,r_ij))
        !!!! Mitch
        !!call c_f(j,k,r,cf,cdf_dr)
        !call scheme_selector(j,k,r,cf,cdf_dr)
        !E_num_F = E_num_F + cf * (r ** (-6))
    !enddo
!enddo
!end function
subroutine E_num_F_ff(c_f)
implicit none
character(len=*),intent(in)         :: c_f
real (kind=sedc_rk)              :: r,cf = 0.0_sedc_rk,cdf_dr = 0.0_sedc_rk
real (kind=sedc_rk),dimension(3) :: r_ij
integer                          :: j,k
E_num_F = 0.0_sedc_rk
do j = 2,sedc_n_ions
    do k = 1,j-1
        r_ij = internal_cart_coord(:,k) - internal_cart_coord(:,j)
        r = sqrt(dot_product(r_ij,r_ij))
        !!! Mitch
        !call c_f(j,k,r,cf,cdf_dr)
        call scheme_selector(j,k,r,cf,cdf_dr)
        E_num_F = E_num_F + cf * (r ** (-6))
    enddo
enddo
end subroutine

! TODO: [ ]  func->subrout conversion, check if it works
!function E_num_F_pbc(c_f)
!implicit none
!real (kind=sedc_rk)              :: E_num_F_pbc
!!external                         :: c_f
  !character(len=*),intent(in)                                  :: c_f
!real (kind=sedc_rk)              :: r,cf = 0.0_sedc_rk,cdf_dr = 0.0_sedc_rk
!real (kind=sedc_rk),dimension(3) :: r_ij
!integer                          :: j,k

!E_num_F_pbc = 0.0_sedc_rk
!do j = 2,sedc_n_ions
    !do k = 1,j-1
        !r_ij = internal_cart_coord(:,k) - internal_cart_coord(:,j)
        !r = sqrt(dot_product(r_ij,r_ij))
        !!!! Mitch
        !!call c_f(j,k,r,cf,cdf_dr)
        !call scheme_selector(j,k,r,cf,cdf_dr)
        !E_num_F_pbc = E_num_F_pbc + cf * (r ** (-6))
    !enddo
!enddo

!sedc_energy = 0.0_sedc_rk
!call sedc_pbc(c_f,.false.,num_F_coord_shift)
!E_num_F_pbc = E_num_F_pbc + sedc_energy
!end function
subroutine E_num_F_pbc_ff(c_f)
implicit none
character(len=*),intent(in)         :: c_f
real (kind=sedc_rk)              :: r,cf = 0.0_sedc_rk,cdf_dr = 0.0_sedc_rk
real (kind=sedc_rk),dimension(3) :: r_ij
integer                          :: j,k

E_num_F_pbc = 0.0_sedc_rk
do j = 2,sedc_n_ions
    do k = 1,j-1
        r_ij = internal_cart_coord(:,k) - internal_cart_coord(:,j)
        r = sqrt(dot_product(r_ij,r_ij))
        !!! Mitch
        !call c_f(j,k,r,cf,cdf_dr)
        call scheme_selector(j,k,r,cf,cdf_dr)
        E_num_F_pbc = E_num_F_pbc + cf * (r ** (-6))
    enddo
enddo

sedc_energy = 0.0_sedc_rk
call sedc_pbc(c_f,.false.,num_F_coord_shift)
E_num_F_pbc = E_num_F_pbc + sedc_energy
end subroutine







! TODO: [ ] fct->subrout, check if it works :)
subroutine ions_ij_to_c6_ij_ff(ion_i,ion_j)
implicit none
integer,intent(in)  :: ion_i,ion_j
integer             :: spec_ind_i = 0,spec_ind_j = 0
spec_ind_i = ion_species_indices(ion_i)
spec_ind_j = ion_species_indices(ion_j)
ions_ij_to_c6_ij = c6_ij(spec_ind_i,spec_ind_j)
end subroutine 

! TODO: [ ] fct->subrout, check if it works :)
subroutine ions_ij_to_R0_ij_ff(ion_i,ion_j)
implicit none
integer,intent(in)  :: ion_i,ion_j
integer             :: spec_ind_i = 0,spec_ind_j = 0
spec_ind_i = ion_species_indices(ion_i)
spec_ind_j = ion_species_indices(ion_j)
ions_ij_to_R0_ij = R0_ij(spec_ind_i,spec_ind_j)
end subroutine               

! TODO: [ ] func->subrout conversion, check if it works!
!function inv_matrix(M)
!implicit none
!real (kind=sedc_rk),dimension(3,3)            :: inv_matrix,I
!real (kind=sedc_rk),dimension(3,3),intent(in) :: M
!real (kind=sedc_rk)                           :: det_M
!integer                                       :: j
!I = 0.0_sedc_rk; I(1,1) = 1.0_sedc_rk; I(2,2) = 1.0_sedc_rk; I(3,3) = 1.0_sedc_rk
!det_M = 0.0_sedc_rk; det_M = det( M(:,1), M(:,2), M(:,3) )

!if (det_M == 0.0_sedc_rk) then
    !call sedc_errquit("SEDC PBC Initialization ERROR!", &
    !& "Cell matrix not invertible - exiting.")
    !return
!endif

!do j = 1, 3
    !inv_matrix(1,j) = det(  I(:,j),  M(:,2),  M(:,3) )
    !inv_matrix(2,j) = det(  M(:,1),  I(:,j),  M(:,3) )
    !inv_matrix(3,j) = det(  M(:,1),  M(:,2),  I(:,j) )
    !inv_matrix(:,j) = inv_matrix(:,j) / det_M
!enddo
!end function inv_matrix



subroutine inv_matrix_ff(M)
implicit none
real (kind=sedc_rk),dimension(3,3)            :: I
real (kind=sedc_rk),dimension(3,3),intent(in) :: M
real (kind=sedc_rk)                           :: det_M
integer                                       :: j
I = 0.0_sedc_rk; I(1,1) = 1.0_sedc_rk; I(2,2) = 1.0_sedc_rk; I(3,3) = 1.0_sedc_rk
call det_ff( M(:,1), M(:,2), M(:,3) )
det_M = 0.0_sedc_rk; det_M = det

if (det_M == 0.0_sedc_rk) then
    call sedc_errquit("SEDC PBC Initialization ERROR!", &
    & "Cell matrix not invertible - exiting.")
    return
endif

do j = 1, 3
    call det_ff(  I(:,j),  M(:,2),  M(:,3) )
    inv_matrix(1,j) = det
    call det_ff(  M(:,1),  I(:,j),  M(:,3) )
    inv_matrix(2,j) = det
    call det_ff(  M(:,1),  M(:,2),  I(:,j) )
    inv_matrix(3,j) = det
    inv_matrix(:,j) = inv_matrix(:,j) / det_M
enddo
end subroutine

! TODO; [ ] func->subrout: check if it works!
subroutine det_ff(a,b,c)
implicit none
real (kind=sedc_rk),dimension(3),intent(in) :: a,b,c
real (kind=sedc_rk),dimension(3)            :: b_cross_c
b_cross_c(1) = b(2) * c(3) - b(3) * c(2)
b_cross_c(2) = b(3) * c(1) - b(1) * c(3)
b_cross_c(3) = b(1) * c(2) - b(2) * c(1)
det = dot_product(a, b_cross_c )
end subroutine



! TODO: [ ]  fct->subrout conversion; check for consistency
subroutine sedc_skip_group_ff(i, j)
implicit none
integer, intent(in)       :: i, j
integer                   :: g, ion_count

! TODO: [x] we needed to change that since it'd have been changed in input routines! -> chgd back
!if (.not.sedc_pbc_file_read) then
if (.not.sedc_pbc_file_read) then
    sedc_skip_group = 1.0_sedc_rk
    return
end if

if (any(sedc_pbc_g_only_intra)) then
    ion_count = 0
    sedc_skip_group = 1.0_sedc_rk * sedc_skip_atom(i) * sedc_skip_atom(j)
    do g = 1, sedc_n_groups
        if (i<=(sedc_groups(g)+ion_count).and.i>ion_count .and. &
        &   j<=(sedc_groups(g)+ion_count).and.j>ion_count.and.sedc_pbc_g_only_intra(g) ) then
            sedc_skip_group = 0.0_sedc_rk * sedc_skip_atom(i) * sedc_skip_atom(j)
        end if
        ion_count = ion_count + sedc_groups(g)
    end do
    
else
    sedc_skip_group = 1.0_sedc_rk * sedc_skip_atom(i) * sedc_skip_atom(j)
end if

end subroutine 



subroutine sedc_errquit(header,msg,warning)
implicit none
character(len=*),intent(in)    :: header,msg
!f2py intent(in)            :: header,msg
logical,optional,intent(in) :: warning
!f2py intent(in)            :: warning
integer                     :: out_unit = 0

if (present(warning)) then
   if (warning) then
      out_unit = sedc_stdout_unit
   else
      call sedc_exit_cleanup()
      sedc_exec_status = 1
      out_unit = sedc_stderr_unit
   end if
else
   call sedc_exit_cleanup()
   sedc_exec_status = 1
   out_unit = sedc_stderr_unit
endif

write (out_unit,"(A,/,2(3(A),/),A)") lstart, &
& lstart,tab,header, &
& lstart,tab,msg, &
& lstart

end subroutine

! TODO: still a dummy function
! TODO: [ ] be careful with deallocation, since it would impede access from python;
!           python automatically allocates when we hand over variables
subroutine sedc_exit_cleanup()
implicit none

if (allocated(unique_species))           deallocate(unique_species)
if (allocated(ion_species_indices))      deallocate(ion_species_indices)
if (allocated(c6_ij))                    deallocate(c6_ij)
if (allocated(R0_ij))                    deallocate(R0_ij)
if (allocated(sedc_pbc_g_switches))      deallocate(sedc_pbc_g_switches)
if (allocated(sedc_pbc_g_cells))         deallocate(sedc_pbc_g_cells)
if (allocated(sedc_groups))              deallocate(sedc_groups)
if (allocated(sedc_pbc_g_fold))          deallocate(sedc_pbc_g_fold)
if (allocated(pbc_g_int_scaling))        deallocate(pbc_g_int_scaling)
if (allocated(sedc_species))             deallocate(sedc_species)
if (allocated(sedc_cart_coord))          deallocate(sedc_cart_coord)
if (allocated(sedc_forces))              deallocate(sedc_forces)
if (allocated(internal_cart_coord))      deallocate(internal_cart_coord)
if (allocated(species_to_print))         deallocate(species_to_print)
if (allocated(ion_coordinates_to_print)) deallocate(ion_coordinates_to_print)
if (allocated(sedc_stress))              deallocate(sedc_stress)
if (allocated(gcell_x_inv_cellvecs))     deallocate(gcell_x_inv_cellvecs)
if (allocated(sedc_ts_veff_div_vfree))   deallocate(sedc_ts_veff_div_vfree)
if (allocated(group_E_decomposed))       deallocate(group_E_decomposed)

sedc_n_ions = 0
sedc_is_initialized = .false.

end subroutine

! TODO: [ ] This is simply for debugging purpose, remove after final implementation ----------------
! TODO: [!] Obsolete, we can already hand over arrays from ASE
!subroutine debug_coords()
!implicit none

!sedc_n_ions = 12
!allocate(sedc_species(sedc_n_ions))
!allocate(internal_cart_coord(3,sedc_n_ions))

!! Structure is a PBE optimized benzene molecule
!sedc_species(1:12) = (/ 1,1,1,1,1,1,6,6,6,6,6,6 /)

!internal_cart_coord(:,1) = (/ -6.001018828d0, 5.422821473d0, 0.0d0/)
!internal_cart_coord(:,2) = (/ -3.848586176d0, 4.180242787d0, 0.0d0/)
!internal_cart_coord(:,3) = (/ -1.695807130d0, 5.423341420d0, 0.0d0/)
!internal_cart_coord(:,4) = (/ -1.695987996d0, 7.908905417d0, 0.0d0/)
!internal_cart_coord(:,5) = (/ -3.848874784d0, 9.151764075d0, 0.0d0/)
!internal_cart_coord(:,6) = (/ -6.001461399d0, 7.908625446d0, 0.0d0/)

!internal_cart_coord(:,7) = (/  -5.056698050d0, 5.968478606d0, 0.0d0/)
!internal_cart_coord(:,8) = (/  -3.848405283d0, 5.270870457d0, 0.0d0/)
!internal_cart_coord(:,9) = (/  -2.640443485d0, 5.968438610d0, 0.0d0/)
!internal_cart_coord(:,10) = (/ -2.640570479d0, 7.363661575d0, 0.0d0/)
!internal_cart_coord(:,11) = (/ -3.848589995d0, 8.061116407d0, 0.0d0/)
!internal_cart_coord(:,12) = (/ -5.056809650d0, 7.363488260d0, 0.0d0/)

!sedc_cell_vectors(:,1) = (/ 7.6972109, 0.0, 0.0 /)		! A vec
!sedc_cell_vectors(:,2) = (/ -3.8486054, 6.6659801, 0.0 /)	! B vec
!sedc_cell_vectors(:,3) = (/ 0.0, 0.0, 10.0 /)			! C vec
!!sedc_cell_vectors(:,3) = (/ 0.0, 0.0, 20.7593236 /)

!!sedc_cell_vectors(:,1) = (/ 5.0, 0.0, 0.0 /)
!!sedc_cell_vectors(:,2) = (/ 0.0, 5.0, 0.0 /)
!!sedc_cell_vectors(:,3) = (/ 0.0, 0.0, 5.0 /)

!end subroutine
! --------------------------------------------------------------------------------------------------


end module sdc_recode           
