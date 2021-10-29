module write_mod
  use define
  use fitsio_mod
  use output_sum
  use utility
  implicit none
  !---- this should be accesible within this module.
  character(len=128) :: fname_backup
  private :: fname_backup
contains
!------------------
!-- 2021-05-20, bug-fixed
!-- 2020-10-17, divided the routine into three smaller routines.
!               now the routine saves FITS data for each observer.
!-- 2020-08-30, radial profiles for Stokes parameters are saved in a table format.
!               PEELOFF routines can merge the previous outputs with the current ones.
!-- 2017-07-12, Jout, Jabs, Jin versus xfreq are now saved in a table format.
!-- 2017-06-28, Kwang-il Seon
!------------------
  subroutine write_output_car(filename,grid)
  use define
  implicit none
  character(len=*),    intent(in)    :: filename
  type(grid_type),     intent(inout) :: grid
  integer          :: k
  character(len=4) :: filename_end

  call write_output_basic(trim(filename),grid)

  if (par%save_peeloff) then
     do k = 1, par%nobs
        if (par%nobs == 1) then
           filename_end = ''
        else
           write(filename_end,'(a,i3.3)') '_',k
        endif
        if (par%save_peeloff_2D) call write_output_peeling_2D(trim(filename),grid,observer(k),suffix=trim(filename_end))
        if (par%save_peeloff_3D) call write_output_peeling_3D(trim(filename),grid,observer(k),suffix=trim(filename_end))
     enddo
  endif
  if (par%save_all_photons) call write_output_allph(trim(filename))
  end subroutine write_output_car
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine write_output_basic(filename,grid)
  use define
  implicit none
  character(len=*), intent(in)    :: filename
  type(grid_type),  intent(inout) :: grid
  !--------------------
  integer            :: unit0,unit,status=0
  character(len=128) :: filename1, filename2
  logical            :: file_exist, merge_ok
  integer            :: colnum
  real(real64)       :: nph1, nph2, nph_tot
  real(real64)       :: exetime1, exetime
  real(real64)       :: nscatt_HI1, nscatt_dust1, nscatt_tot1
  real(real64)       :: raccept1, fluxfac1
  real(real64), allocatable :: arr_1D(:), arr_2D(:,:), arr_3D(:,:,:), arr_4D(:,:,:,:)

  !--- check the previous FITS output.
  merge_ok = par%out_merge
  if (merge_ok) then
     inquire(file=trim(filename),exist=file_exist)
     if (.not. file_exist) merge_ok = .false.
  endif

  exetime = par%exetime
  nph_tot = par%no_photons
  if (merge_ok) then
     !--- make a backup file.
     if (par%save_backup) then
        fname_backup = name_for_backup(get_base_name(trim(filename)))
        filename2    = trim(fname_backup)//'.fits.gz'
        !--- copy_file modifies the time stamp of the file.
        !--- call copy_file(trim(par%out_file), trim(filename2), status)
        call execute_command_line("cp -p "//trim(par%out_file)//" "//trim(filename2), exitstat=status)
     endif

     call fits_open_old(unit0,trim(filename),status)
     call fits_move_to_next_hdu(unit0,status)

     call fits_get_keyword(unit0,'exetime',  exetime1,    status)
     call fits_get_keyword(unit0,'nphotons', nph1,        status)
     call fits_get_keyword(unit0,'Nsc_HI',   nscatt_HI1,  status)
     call fits_get_keyword(unit0,'Nsc_dust', nscatt_dust1,status)
     call fits_get_keyword(unit0,'Nsc_tot',  nscatt_tot1, status)
     !--- do not update par%exetime, which will give information for each run, not for all accumulated runs.
     !--- update the number of scatterings to improve statistics (2020.09.06).
     exetime         = par%exetime + exetime1
     nph2            = par%no_photons
     nph_tot         = nph1 + nph2
     par%nscatt_HI   = (nscatt_HI1  *nph1 + par%nscatt_HI  *nph2)/nph_tot
     par%nscatt_dust = (nscatt_dust1*nph1 + par%nscatt_dust*nph2)/nph_tot
     par%nscatt_tot  = (nscatt_tot1 *nph1 + par%nscatt_tot *nph2)/nph_tot

     if (trim(par%source_geometry) == 'stellar_illumination') then
        call fits_get_keyword(unit0,'Raccept',  raccept1, status)
        call fits_get_keyword(unit0,'fluxfac',  fluxfac1, status)
        par%acceptance_rate  = (raccept1 *nph1 + par%acceptance_rate *nph2)/nph_tot
        par%flux_factor      = (fluxfac1 *nph1 + par%flux_factor     *nph2)/nph_tot
     endif

     if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%xfreq)
     call fits_get_column_number(unit0,'Jout',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D,status)
     grid%Jout = (arr_1D * nph1 + grid%Jout * nph2)/nph_tot

     if (par%save_Jabs) then
        call fits_get_column_number(unit0,'Jabs',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D,status)
        grid%Jabs = (arr_1D * nph1 + grid%Jabs * nph2)/nph_tot
     endif
     if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
        call fits_get_column_number(unit0,'Jabs2',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D,status)
        grid%Jabs2 = (arr_1D * nph1 + grid%Jabs2 * nph2)/nph_tot
     endif
     if (par%save_Jin) then
        call fits_get_column_number(unit0,'Jin',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D,status)
        grid%Jin = (arr_1D * nph1 + grid%Jin * nph2)/nph_tot
     endif
     if (allocated(arr_1D)) deallocate(arr_1D)

#ifdef CALCP
     if (associated(grid%Pa)) then
        !--- save P_alapha(x,y,z) : scattering number per atom per photon
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%Pa)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        grid%Pa = (arr_3D * nph1 + grid%Pa * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     if (associated(grid%P1)) then
        !--- save radial or z profile, P_alpha
        if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%P1)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_1D,status)
        grid%P1 = (arr_1D * nph1 + grid%P1 * nph2)/nph_tot
        if (allocated(arr_1D)) deallocate(arr_1D)
     endif

     if (associated(grid%P2)) then
        !--- save cylindrical (r, z) profile, P_alpha
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%P2)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        grid%P2 = (arr_2D * nph1 + grid%P2 * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif
#endif

#ifdef CALCPnew
     if (associated(grid%Pa_new)) then
        !--- save P_alapha(x,y,z) : scattering number per atom per photon
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%Pa_new)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        grid%Pa_new = (arr_3D * nph1 + grid%Pa_new * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     if (associated(grid%P1_new)) then
        !--- save radial or z profile, P_alpha
        if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%P1_new)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_1D,status)
        grid%P1_new = (arr_1D * nph1 + grid%P1_new * nph2)/nph_tot
        if (allocated(arr_1D)) deallocate(arr_1D)
     endif

     if (associated(grid%P2_new)) then
        !--- save cylidrical (r, z)  profile, P_alpha
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%P2_new)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        grid%P2_new = (arr_2D * nph1 + grid%P2_new * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif
#endif

#ifdef CALCJ
     if (associated(grid%J)) then
        !--- save J(nu,x,y,z) : mean intensity spectrum at (x,y,z)
        if (.not. allocated(arr_4D)) allocate(arr_4D, source=grid%J)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_4D,status)
        grid%J = (arr_4D * nph1 + grid%J * nph2)/nph_tot
        if (allocated(arr_4D)) deallocate(arr_4D)
     endif

     if (associated(grid%J1)) then
        !--- save radial or z profile of mean intensity spectrum J(nu)
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%J1)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        grid%J1 = (arr_2D * nph1 + grid%J1 * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif

     if (associated(grid%J2)) then
        !--- save cylindrical (r, z) profile of mean intensity spectrum J(nu)
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%J2)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        grid%J2 = (arr_3D * nph1 + grid%J2 * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif
#endif
     call fits_close(unit0,status)
  endif

  !--- open main FITS file.
  call fits_open_new(unit,trim(filename),status)

  !--- xfreq, (relative) frequency, velocity, and dlambda
  call fits_append_table_column(unit,'Xfreq',   grid%xfreq,   status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'velocity',grid%velocity,status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'lambda',  grid%lambda,  status,bitpix=par%out_bitpix)

  !--- Jout, emerging spectrum
  call fits_append_table_column(unit,'Jout',grid%Jout,status,bitpix=par%out_bitpix)

  !--- Jabs, absorbed spectrum
  if (par%save_Jabs) then
     call fits_append_table_column(unit,'Jabs',grid%Jabs,status,bitpix=par%out_bitpix)
  endif
  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
     call fits_append_table_column(unit,'Jabs2',grid%Jabs2,status,bitpix=par%out_bitpix)
  endif

  !--- Jin, input spectrum
  if (par%save_Jin) then
     call fits_append_table_column(unit,'Jin',grid%Jin,status,bitpix=par%out_bitpix)
  endif

  call fits_put_keyword(unit,'ExeTime',  exetime,            'Excution Time (min)',status)
  call fits_put_keyword(unit,'Nprocs',   mpar%num_procs,     'No. of Threads',status)
  call fits_put_keyword(unit,'recoil',   par%recoil,         'recoil',status)
  call fits_put_keyword(unit,'coreskip', par%core_skip,      'coreskip algorithm',status)
  call fits_put_keyword(unit,'xyz_sym',  par%xyz_symmetry,   'xyz_symmetry',status)
  call fits_put_keyword(unit,'xy_per',   par%xy_periodic,    'xy_periodic',status)
  call fits_put_keyword(unit,'save_all', par%save_all,       'save all 3D output?',status)
  call fits_put_keyword(unit,'save_Jin', par%save_Jin,       'save input J?',status)
  call fits_put_keyword(unit,'save_Jab', par%save_Jabs,      'save absorbed J?',status)
  call fits_put_keyword(unit,'nphotons', nph_tot,            'number of photons',status)
  call fits_put_keyword(unit,'taumax',   par%taumax,         'tau_max',status)
  call fits_put_keyword(unit,'tauhomo',  par%tauhomo,        'tau_homo',status)
  call fits_put_keyword(unit,'N_HImax',  par%N_HImax,        'N(HI)_max cm^-2',status)
  call fits_put_keyword(unit,'N_HIhomo', par%N_HIhomo,       'N(HI)_homo cm^-2',status)
  call fits_put_keyword(unit,'temp',     par%temperature,    'temperature (K)',status)
  call fits_put_keyword(unit,'Vexp',     par%Vexp,           'Velocity (km/s)',status)
  call fits_put_keyword(unit,'DGR',      par%DGR,            'Dust-to-Gas Ratio',status)
  call fits_put_keyword(unit,'atau3',    par%atau3,          '(a*tau_max)^(1/3)',status)
  call fits_put_keyword(unit,'voigta',   grid%voigt_amean,   'voigt_a',status)
  call fits_put_keyword(unit,'Xfreq1',   grid%xfreq_min,     'Xfreq_min',status)
  call fits_put_keyword(unit,'Xfreq2',   grid%xfreq_max,     'Xfreq_max',status)
  call fits_put_keyword(unit,'Dxfreq',   grid%dxfreq,        'Dxfreq', status)
  call fits_put_keyword(unit,'Dlambda',  grid%dlambda,       'Dlambda (angstrom)', status)
  call fits_put_keyword(unit,'I_unit',   par%intensity_unit, 'Intensity Unit (0:no dimension, 1:cm^-2 A^-1)', status)
  call fits_put_keyword(unit,'Dfreq',    grid%Dfreq_ref,     'Doppler Freq. (Hz)',status)
  call fits_put_keyword(unit,'Nsc_dust', par%nscatt_dust,    'Nscatt_dust/photon',status)
  call fits_put_keyword(unit,'Nsc_HI',   par%nscatt_HI,      'Nscatt_HI/photon',status)
  call fits_put_keyword(unit,'Nsc_tot',  par%nscatt_tot,     'Nscatt_tot/photon',status)
  call fits_put_keyword(unit,'nx',       grid%nx,            'No. of x cells',status)
  call fits_put_keyword(unit,'ny',       grid%ny,            'No. of y cells',status)
  call fits_put_keyword(unit,'nz',       grid%nz,            'No. of z cells',status)
  call fits_put_keyword(unit,'xmax',     par%xmax,           'xmax',status)
  call fits_put_keyword(unit,'ymax',     par%ymax,           'ymax',status)
  call fits_put_keyword(unit,'zmax',     par%zmax,           'zmax',status)
  call fits_put_keyword(unit,'EXTNAME',  'Spectrum',         'Spectrum',status)
#ifdef CALCP
  call fits_put_keyword(unit,'calc_P',   .true.,             'calculated Pa within media?',status)
#else
  call fits_put_keyword(unit,'calc_P',   .false.,            'calculated Pa within media?',status)
#endif
#ifdef CALCPnew
  call fits_put_keyword(unit,'calc_Pnew',   .true.,          'calculated Pa within media?',status)
#else
  call fits_put_keyword(unit,'calc_Pnew',   .false.,         'calculated Pa within media?',status)
#endif
#ifdef CALCJ
  call fits_put_keyword(unit,'calc_J',   .true.,             'calculated J within media?',status)
#else
  call fits_put_keyword(unit,'calc_J',   .false.,            'calculated J within media?',status)
#endif
  if (trim(par%source_geometry) == 'stellar_illumination') then
     call fits_put_keyword(unit,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
     call fits_put_keyword(unit,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call fits_put_keyword(unit,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     call fits_put_keyword(unit,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
     call fits_put_keyword(unit,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
     call fits_put_keyword(unit,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
  endif
!---------------------------------

#ifdef CALCP
  if (associated(grid%Pa)) then
     !--- save P_alapha(x,y,z) : scattering number per atom per photon
     call fits_append_image(unit,grid%Pa,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_3D',  'P_alpha (number of scattering)',status)
  endif

  if (associated(grid%P1)) then
     !--- save radial or z profile, P_alpha
     call fits_append_image(unit,grid%P1,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_1D',  'P_alpha (number of scattering)',status)
  endif

  if (associated(grid%P2)) then
     !--- save cylindrical (r, z) profile, P_alpha
     call fits_append_image(unit,grid%P2,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_2D',  'P_alpha (number of scattering)',status)
  endif
#endif

#ifdef CALCPnew
  if (associated(grid%Pa_new)) then
     !--- save P_alapha(x,y,z) : scattering number per atom per photon
     call fits_append_image(unit,grid%Pa_new,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_3D',  'P_alpha (new_method, number of scattering)',status)
  endif

  if (associated(grid%P1_new)) then
     !--- save radial or z profile, P_alpha
     call fits_append_image(unit,grid%P1_new,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_1D',  'P_alpha (new_method, number of scattering)',status)
  endif

  if (associated(grid%P2_new)) then
     !--- save cylindrical (r, z) profile, P_alpha
     call fits_append_image(unit,grid%P2_new,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_2D',  'P_alpha (new_method, number of scattering)',status)
  endif
#endif

#ifdef CALCJ
  if (associated(grid%J)) then
     !--- save J(nu,x,y,z) : mean intensity spectrum at (x,y,z)
     call fits_append_image(unit,grid%J,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Jx_3D','J(x) (mean intensity)',status)
  endif

  if (associated(grid%J1)) then
     !--- save radial or z profile of mean intensity spectrum J(nu)
     call fits_append_image(unit,grid%J1,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Jx_1D','J(x) (mean intensity)',status)
  endif

  if (associated(grid%J2)) then
     !--- save cylindrical (r, z) profile of mean intensity spectrum J(nu)
     call fits_append_image(unit,grid%J2,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Jx_2D','J(x) (mean intensity)',status)
  endif
#endif

  !--- close the FITS file.
  call fits_close(unit,status)
  end subroutine write_output_basic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine write_output_peeling_2D(filename,grid,obs,suffix)
  use define
  implicit none
  character(len=*),    intent(in)    :: filename
  type(grid_type),     intent(inout) :: grid
  type(observer_type), intent(inout) :: obs
  character(len=*), optional, intent(in) :: suffix
  !--------------------
  integer            :: unit0,unit,status=0
  character(len=128) :: filename1, filename2, filename_end
  logical            :: file_exist, merge_ok
  real(real64)       :: nph1, nph2, nph_tot
  real(real64)       :: cd1_1, cd1_2, cd2_1, cd2_2
  real(real64)       :: crpix1, crpix2, crval1, crval2
  integer            :: equinox = 2000
  real(real64), allocatable :: arr_2D(:,:)

  if (present(suffix)) then
     filename_end = trim(suffix)
  else
     filename_end = ''
  endif

  !--- Initialize FITS file name.
  filename1 = trim(get_base_name(filename))//'_obs2D'//trim(filename_end)//'.fits.gz'

  !--- check the previous FITS output.
  merge_ok = par%out_merge
  if (merge_ok) then
     inquire(file=trim(filename1),exist=file_exist)
     if (.not. file_exist) merge_ok = .false.
  endif

  !--- number of photons
  nph_tot = par%no_photons

  !--- read the previous outputs.
  if (merge_ok) then
     !--- make a backup file.
     if (par%save_backup) then
        filename2    = trim(fname_backup)//'_obs2D'//trim(filename_end)//'.fits.gz'
        !call copy_file(trim(filename1), trim(filename2), status)
        call execute_command_line("cp -p "//trim(filename1)//" "//trim(filename2), exitstat=status)
     endif

     if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%scatt_2D)
     call fits_open_old(unit0,trim(filename1),status)
     call fits_get_keyword(unit0,'nphotons', nph1, status)
     nph2    = par%no_photons
     nph_tot = nph1 + nph2

     call fits_read_image(unit0,arr_2D,status)
     obs%scatt_2D = (arr_2D * nph1 + obs%scatt_2D * nph2)/nph_tot
     call fits_move_to_next_hdu(unit0,status)
     call fits_read_image(unit0,arr_2D,status)
     obs%direc_2D = (arr_2D * nph1 + obs%direc_2D * nph2)/nph_tot
     if (allocated(arr_2D)) deallocate(arr_2D)

     if (par%save_direc0) then
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%direc0_2D)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        obs%direc0_2D = (arr_2D * nph1 + obs%direc0_2D * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif

     call fits_close(unit0,status)
  endif

  !--- We will make the radial spectral or intensity profiles.
  if (par%save_radial_profile) call make_radial_intensity(grid,obs,use_2D_data=.true.)

  !--- open FITS file for peel-off.
  call fits_open_new(unit,trim(filename1),status)

  !--- write scattered data
  call fits_append_image(unit,obs%scatt_2D,status,bitpix=par%out_bitpix)

  !--- header keyword for 2D image
  cd1_1  = par%dxim
  cd1_2  = 0.0_wp
  cd2_1  = 0.0_wp
  cd2_2  = par%dyim
  crpix1 = (par%nxim+1)/2.0_wp
  crpix2 = (par%nyim+1)/2.0_wp
  crval1 = 0.0_wp
  crval2 = 0.0_wp
  call fits_put_keyword(unit,'EXTNAME','Scattered','J(freq,x,y) (intensity)',status)
  call fits_put_keyword(unit,'EQUINOX' ,equinox,   'Equinox of Ref. Coord.' ,status)
  call fits_put_keyword(unit,'CD1_1'   ,cd1_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_1'   ,cd2_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD1_2'   ,cd1_2  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_2'   ,cd2_2  ,   'Degree / Pixel' ,status)
  call fits_put_keyword(unit,'CRPIX1'  ,crpix1 ,   'Reference Pixel in X',status)
  call fits_put_keyword(unit,'CRPIX2'  ,crpix2 ,   'Reference Pixel in Y',status)
  call fits_put_keyword(unit,'CRVAL1'  ,crval1 ,   'R.A. (Degree)',status)
  call fits_put_keyword(unit,'CRVAL2'  ,crval2 ,   'Dec  (Degree)',status)
  call fits_put_keyword(unit,'CTYPE1'  ,'RA--TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'CTYPE2'  ,'DEC-TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
  call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)
  call fits_put_keyword(unit,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
  call fits_put_keyword(unit,'nphotons', nph_tot,          'number of photons',status)
  call fits_put_keyword(unit,'alpha',    obs%alpha, 'alpha (degree)',status)
  call fits_put_keyword(unit,'beta',     obs%beta,  'beta (degree)',status)
  call fits_put_keyword(unit,'gamma',    obs%gamma, 'gamma (degree)',status)
  call fits_put_keyword(unit,'obsx',     obs%x,     'Observer X coordinate',status)
  call fits_put_keyword(unit,'obsy',     obs%y,     'Observer Y coordinate',status)
  call fits_put_keyword(unit,'obsz',     obs%z,     'Observer Z coordinate',status)
  call fits_put_keyword(unit,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
  call fits_put_keyword(unit,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
  call fits_put_keyword(unit,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
  if (trim(par%source_geometry) == 'stellar_illumination') then
     call fits_put_keyword(unit,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
     call fits_put_keyword(unit,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call fits_put_keyword(unit,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     call fits_put_keyword(unit,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
     call fits_put_keyword(unit,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
     call fits_put_keyword(unit,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
  endif

  !--- write direct data
  call fits_append_image(unit,obs%direc_2D,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','Direct','J(x,y) (intensity)',status)

  if (par%save_direc0) then
     call fits_append_image(unit,obs%direc0_2D,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Direct0','J(x,y) (intensity)',status)
  endif

  if (associated(obs%radial_I)) then
     call fits_append_table_column(unit,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','RadialI','I(r) (radial profile)',status)
  endif

  !--- close FITS file
  call fits_close(unit,status)

  if (par%use_stokes) then
     !--- open new FITS file
     filename1 = trim(get_base_name(filename))//'_stokes_2D'//trim(filename_end)//'.fits.gz'

     !--- check the previous FITS output.
     merge_ok = par%out_merge
     if (merge_ok) then
        inquire(file=trim(filename1),exist=file_exist)
        if (.not. file_exist) merge_ok = .false.
     endif

     !--- read the previous outputs.
     if (merge_ok) then
        !--- make a backup file.
        if (par%save_backup) then
           filename2 = trim(fname_backup)//'_stokes_2D'//trim(filename_end)//'.fits.gz'
           !call copy_file(trim(filename1), trim(filename2), status)
           call execute_command_line("cp -p "//trim(filename1)//" "//trim(filename2), exitstat=status)
        endif

        call fits_open_old(unit0,trim(filename1),status)
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%I_2D)
        call fits_read_image(unit0,arr_2D,status)
        obs%I_2D = (arr_2D * nph1 + obs%I_2D * nph2)/nph_tot
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        obs%Q_2D = (arr_2D * nph1 + obs%Q_2D * nph2)/nph_tot
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        obs%U_2D = (arr_2D * nph1 + obs%U_2D * nph2)/nph_tot
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        obs%V_2D = (arr_2D * nph1 + obs%V_2D * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
        call fits_close(unit0,status)
     endif

     !--- Making an average of radial_pol is a non-sense, because pol = sqrt(Q^2 + U^2)/I.
     !--- Hence, we will make the radial profiles for Stokes parameters from the 3D Stokes parameters,
     !---        instead of averaging the pre-existing two radial profiles.
     if (par%save_radial_profile) call make_radial_stokes(grid,obs,use_2D_data=.true.)

     !--- open FITS file for Stokes parameters.
     call fits_open_new(unit,trim(filename1),status)

     !--- write Stokes I data
     call fits_append_image(unit,obs%I_2D,status,bitpix=par%out_bitpix)

     !--- header keyword for 2D image
     cd1_1  = par%dxim
     cd1_2  = 0.0_wp
     cd2_1  = 0.0_wp
     cd2_2  = par%dyim
     crpix1 = (par%nxim+1)/2.0_wp
     crpix2 = (par%nyim+1)/2.0_wp
     crval1 = 0.0_wp
     crval2 = 0.0_wp
     call fits_put_keyword(unit,'EXTNAME','Stokes_I','Stokes I image',status)
     call fits_put_keyword(unit,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
     call fits_put_keyword(unit,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
     call fits_put_keyword(unit,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
     call fits_put_keyword(unit,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
     call fits_put_keyword(unit,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
     call fits_put_keyword(unit,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
     call fits_put_keyword(unit,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
     call fits_put_keyword(unit,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
     call fits_put_keyword(unit,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
     call fits_put_keyword(unit,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
     call fits_put_keyword(unit,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
     call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
     call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)
     call fits_put_keyword(unit,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
     call fits_put_keyword(unit,'nphotons', nph_tot,          'number of photons',status)
     call fits_put_keyword(unit,'alpha',    obs%alpha, 'alpha (degree)',status)
     call fits_put_keyword(unit,'beta',     obs%beta,  'beta (degree)',status)
     call fits_put_keyword(unit,'gamma',    obs%gamma, 'gamma (degree)',status)
     call fits_put_keyword(unit,'obsx',     obs%x,     'Observer X coordinate',status)
     call fits_put_keyword(unit,'obsy',     obs%y,     'Observer Y coordinate',status)
     call fits_put_keyword(unit,'obsz',     obs%z,     'Observer Z coordinate',status)
     call fits_put_keyword(unit,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
     call fits_put_keyword(unit,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
     call fits_put_keyword(unit,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
     if (trim(par%source_geometry) == 'stellar_illumination') then
        call fits_put_keyword(unit,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
        call fits_put_keyword(unit,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
        call fits_put_keyword(unit,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
        call fits_put_keyword(unit,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
        call fits_put_keyword(unit,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
        call fits_put_keyword(unit,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
     endif

     !--- write Stokes Q, U, V data
     call fits_append_image(unit,obs%Q_2D,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Stokes_Q','Stokes Q image',status)
     call fits_append_image(unit,obs%U_2D,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Stokes_U','Stokes U image',status)
     call fits_append_image(unit,obs%V_2D,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Stokes_V','Stokes V image',status)

     !--- radial Stokes profiles
     if (associated(obs%radial_r)) then
        call fits_append_table_column(unit,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'Q',  obs%radial_Q,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'U',  obs%radial_U,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'V',  obs%radial_V,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'pol',obs%radial_pol,status,bitpix=par%out_bitpix)
        call fits_put_keyword(unit,'EXTNAME','Stokes_radial','Stokes radial profile',status)
     endif

     !--- close the FITS file
     call fits_close(unit,status)
  endif
  end subroutine write_output_peeling_2D
  !-----------------------------------------------------------------
  subroutine write_output_peeling_3D(filename,grid,obs,suffix)
  use define
  implicit none
  character(len=*),    intent(in)    :: filename
  type(grid_type),     intent(inout) :: grid
  type(observer_type), intent(inout) :: obs
  character(len=*), optional, intent(in) :: suffix
  !--------------------
  integer            :: unit0,unit,status=0
  character(len=128) :: filename1, filename2, filename_end
  logical            :: file_exist, merge_ok
  real(real64)       :: nph1, nph2, nph_tot
  real(real64)       :: cd1_1, cd1_2, cd2_1, cd2_2
  real(real64)       :: crpix1, crpix2, crval1, crval2
  integer            :: equinox = 2000
  real(real64), allocatable :: arr_2D(:,:), arr_3D(:,:,:)

  if (present(suffix)) then
     filename_end = trim(suffix)
  else
     filename_end = ''
  endif

  !--- Initialize FITS file name.
  filename1 = trim(get_base_name(filename))//'_obs'//trim(filename_end)//'.fits.gz'

  !--- check the previous FITS output.
  merge_ok = par%out_merge
  if (merge_ok) then
     inquire(file=trim(filename1),exist=file_exist)
     if (.not. file_exist) merge_ok = .false.
  endif

  !--- number of photons
  nph_tot = par%no_photons

  !--- read the previous outputs.
  if (merge_ok) then
     !--- make a backup file.
     if (par%save_backup) then
        filename2    = trim(fname_backup)//'_obs'//trim(filename_end)//'.fits.gz'
        !call copy_file(trim(filename1), trim(filename2), status)
        call execute_command_line("cp -p "//trim(filename1)//" "//trim(filename2), exitstat=status)
     endif

     if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%scatt)
     call fits_open_old(unit0,trim(filename1),status)
     call fits_get_keyword(unit0,'nphotons', nph1, status)
     nph2    = par%no_photons
     nph_tot = nph1 + nph2

     call fits_read_image(unit0,arr_3D,status)
     obs%scatt = (arr_3D * nph1 + obs%scatt * nph2)/nph_tot
     call fits_move_to_next_hdu(unit0,status)
     call fits_read_image(unit0,arr_3D,status)
     obs%direc = (arr_3D * nph1 + obs%direc * nph2)/nph_tot
     if (allocated(arr_3D)) deallocate(arr_3D)

     if (par%save_direc0) then
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%direc0)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        obs%direc0 = (arr_3D * nph1 + obs%direc0 * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     call fits_close(unit0,status)
  endif

  !--- We will make the radial spectral or intensity profiles.
  if (par%save_radial_profile) call make_radial_intensity(grid,obs)

  !--- open FITS file for peel-off.
  call fits_open_new(unit,trim(filename1),status)

  !--- write scattered data
  call fits_append_image(unit,obs%scatt,status,bitpix=par%out_bitpix)

  !--- header keyword for 3D spectral image
  cd1_1  = par%dxim
  cd1_2  = 0.0_wp
  cd2_1  = 0.0_wp
  cd2_2  = par%dyim
  crpix1 = (par%nxim+1)/2.0_wp
  crpix2 = (par%nyim+1)/2.0_wp
  crval1 = 0.0_wp
  crval2 = 0.0_wp
  call fits_put_keyword(unit,'EXTNAME','Scattered','J(freq,x,y) (intensity)',status)
  call fits_put_keyword(unit,'EQUINOX' ,equinox,   'Equinox of Ref. Coord.' ,status)
  call fits_put_keyword(unit,'CD1_1'   ,cd1_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_1'   ,cd2_1  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD1_2'   ,cd1_2  ,   'Degree / Pixel',status)
  call fits_put_keyword(unit,'CD2_2'   ,cd2_2  ,   'Degree / Pixel' ,status)
  call fits_put_keyword(unit,'CRPIX1'  ,crpix1 ,   'Reference Pixel in X',status)
  call fits_put_keyword(unit,'CRPIX2'  ,crpix2 ,   'Reference Pixel in Y',status)
  call fits_put_keyword(unit,'CRVAL1'  ,crval1 ,   'R.A. (Degree)',status)
  call fits_put_keyword(unit,'CRVAL2'  ,crval2 ,   'Dec  (Degree)',status)
  call fits_put_keyword(unit,'CTYPE1'  ,'RA--TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'CTYPE2'  ,'DEC-TAN', 'Coordinate Type',status)
  call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
  call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)
  call fits_put_keyword(unit,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
  call fits_put_keyword(unit,'Xfreq1'  , grid%xfreq_min,   'Xfreq_min',status)
  call fits_put_keyword(unit,'Xfreq2'  , grid%xfreq_max,   'Xfreq_max',status)
  call fits_put_keyword(unit,'Dxfreq'  , grid%dxfreq,      'Dxfreq', status)
  call fits_put_keyword(unit,'Dlambda',  grid%dlambda,     'Dlambda (angstrom)', status)
  call fits_put_keyword(unit,'I_unit',   par%intensity_unit, 'Intensity Unit (0:no dimension, 1:cm^-2 A^-1)', status)
  call fits_put_keyword(unit,'Dfreq'   , grid%Dfreq_ref,   'Doppler Freq.  (Hz)',status)
  call fits_put_keyword(unit,'nphotons', nph_tot,          'number of photons',status)
  call fits_put_keyword(unit,'alpha',    obs%alpha, 'alpha (degree)',status)
  call fits_put_keyword(unit,'beta',     obs%beta,  'beta (degree)',status)
  call fits_put_keyword(unit,'gamma',    obs%gamma, 'gamma (degree)',status)
  call fits_put_keyword(unit,'obsx',     obs%x,     'Observer X coordinate',status)
  call fits_put_keyword(unit,'obsy',     obs%y,     'Observer Y coordinate',status)
  call fits_put_keyword(unit,'obsz',     obs%z,     'Observer Z coordinate',status)
  call fits_put_keyword(unit,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
  call fits_put_keyword(unit,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
  call fits_put_keyword(unit,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
  if (trim(par%source_geometry) == 'stellar_illumination') then
     call fits_put_keyword(unit,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
     call fits_put_keyword(unit,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call fits_put_keyword(unit,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     call fits_put_keyword(unit,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
     call fits_put_keyword(unit,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
     call fits_put_keyword(unit,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
  endif

  !--- write direct data
  call fits_append_image(unit,obs%direc,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','Direct','J(freq,x,y) (intensity)',status)

  if (par%save_direc0) then
     call fits_append_image(unit,obs%direc0,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Direct0','J(freq,x,y) (intensity)',status)
  endif

  if (associated(obs%radial_I)) then
     call fits_append_table_column(unit,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','RadialI','I(r) (radial profile)',status)
  endif

  !--- close FITS file
  call fits_close(unit,status)

  if (par%use_stokes) then
     !--- open new FITS file
     filename1 = trim(get_base_name(filename))//'_stokes'//trim(filename_end)//'.fits.gz'

     !--- check the previous FITS output.
     merge_ok = par%out_merge
     if (merge_ok) then
        inquire(file=trim(filename1),exist=file_exist)
        if (.not. file_exist) merge_ok = .false.
     endif

     !--- read the previous outputs.
     if (merge_ok) then
        !--- make a backup file.
        if (par%save_backup) then
           filename2 = trim(fname_backup)//'_stokes'//trim(filename_end)//'.fits.gz'
           !call copy_file(trim(filename1), trim(filename2), status)
           call execute_command_line("cp -p "//trim(filename1)//" "//trim(filename2), exitstat=status)
        endif

        call fits_open_old(unit0,trim(filename1),status)
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%I)
        call fits_read_image(unit0,arr_3D,status)
        obs%I = (arr_3D * nph1 + obs%I * nph2)/nph_tot
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        obs%Q = (arr_3D * nph1 + obs%Q * nph2)/nph_tot
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        obs%U = (arr_3D * nph1 + obs%U * nph2)/nph_tot
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        obs%V = (arr_3D * nph1 + obs%V * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
        call fits_close(unit0,status)
     endif

     !--- Making an average of radial_pol is a non-sense, because pol = sqrt(Q^2 + U^2)/I.
     !--- Hence, we will make the radial profiles for Stokes parameters from the 3D Stokes parameters,
     !---        instead of averaging the pre-existing two radial profiles.
     if (par%save_radial_profile) call make_radial_stokes(grid,obs)

     !--- open FITS file for Stokes parameters.
     call fits_open_new(unit,trim(filename1),status)

     !--- write Stokes I data
     call fits_append_image(unit,obs%I,status,bitpix=par%out_bitpix)

     !--- header keyword for 3D spectral image
     cd1_1  = par%dxim
     cd1_2  = 0.0_wp
     cd2_1  = 0.0_wp
     cd2_2  = par%dyim
     crpix1 = (par%nxim+1)/2.0_wp
     crpix2 = (par%nyim+1)/2.0_wp
     crval1 = 0.0_wp
     crval2 = 0.0_wp
     call fits_put_keyword(unit,'EXTNAME','Stokes_I','Stokes I image',status)
     call fits_put_keyword(unit,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
     call fits_put_keyword(unit,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
     call fits_put_keyword(unit,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
     call fits_put_keyword(unit,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
     call fits_put_keyword(unit,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
     call fits_put_keyword(unit,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
     call fits_put_keyword(unit,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
     call fits_put_keyword(unit,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
     call fits_put_keyword(unit,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
     call fits_put_keyword(unit,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
     call fits_put_keyword(unit,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
     call fits_put_keyword(unit,'DISTANCE', par%distance,     'Distance',status)
     call fits_put_keyword(unit,'DISTUNIT', par%distance_unit,'Distance Unit',status)
     call fits_put_keyword(unit,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
     call fits_put_keyword(unit,'Xfreq1' ,grid%xfreq_min, 'Xfreq_min',status)
     call fits_put_keyword(unit,'Xfreq2' ,grid%xfreq_max, 'Xfreq_max',status)
     call fits_put_keyword(unit,'Dxfreq' ,grid%dxfreq,    'Dxfreq', status)
     call fits_put_keyword(unit,'Dlambda',grid%dlambda,   'Dlambda (angstrom)', status)
     call fits_put_keyword(unit,'I_unit' ,par%intensity_unit, 'Intensity Unit (0:no dimension, 1:cm^-2 A^-1)', status)
     call fits_put_keyword(unit,'Dfreq'  ,grid%Dfreq_ref, 'Doppler Freq. (Hz)',status)
     call fits_put_keyword(unit,'nphotons', nph_tot,      'number of photons',status)
     call fits_put_keyword(unit,'alpha',    obs%alpha, 'alpha (degree)',status)
     call fits_put_keyword(unit,'beta',     obs%beta,  'beta (degree)',status)
     call fits_put_keyword(unit,'gamma',    obs%gamma, 'gamma (degree)',status)
     call fits_put_keyword(unit,'obsx',     obs%x,     'Observer X coordinate',status)
     call fits_put_keyword(unit,'obsy',     obs%y,     'Observer Y coordinate',status)
     call fits_put_keyword(unit,'obsz',     obs%z,     'Observer Z coordinate',status)
     call fits_put_keyword(unit,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
     call fits_put_keyword(unit,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
     call fits_put_keyword(unit,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
     if (trim(par%source_geometry) == 'stellar_illumination') then
        call fits_put_keyword(unit,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
        call fits_put_keyword(unit,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
        call fits_put_keyword(unit,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
        call fits_put_keyword(unit,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
        call fits_put_keyword(unit,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
        call fits_put_keyword(unit,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
     endif

     !--- write Stokes Q, U, V data
     call fits_append_image(unit,obs%Q,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Stokes_Q','Stokes Q image',status)
     call fits_append_image(unit,obs%U,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Stokes_U','Stokes U image',status)
     call fits_append_image(unit,obs%V,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Stokes_V','Stokes V image',status)

     !--- radial Stokes profiles
     if (associated(obs%radial_r)) then
        call fits_append_table_column(unit,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'Q',  obs%radial_Q,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'U',  obs%radial_U,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'V',  obs%radial_V,  status,bitpix=par%out_bitpix)
        call fits_append_table_column(unit,'pol',obs%radial_pol,status,bitpix=par%out_bitpix)
        call fits_put_keyword(unit,'EXTNAME','Stokes_radial','Stokes radial profile',status)
     endif

     !--- close the FITS file
     call fits_close(unit,status)
  endif
  end subroutine write_output_peeling_3D
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine write_output_allph(filename)
  use define
  implicit none
  character(len=*),  intent(in) :: filename

  !--------------------
  integer            :: unit0,unit,status=0
  character(len=128) :: filename1, filename2
  logical            :: file_exist, merge_ok
  integer            :: colnum
  integer            :: nph1, nph
  real(real64), allocatable :: arr_1D_old(:), arr_1D_new(:)

  !--- Initialize FITS file name.
  filename1 = trim(get_base_name(filename))//'_allph.fits.gz'

  !--- check the previous FITS output.
  merge_ok = par%out_merge
  if (merge_ok) then
     inquire(file=trim(filename1),exist=file_exist)
     if (.not. file_exist) merge_ok = .false.
  endif

  if (merge_ok) then
     if (par%save_backup) then
        filename2    = trim(fname_backup)//'_allph.fits.gz'
        !call copy_file(trim(filename1), trim(filename2), status)
        call execute_command_line("cp -p "//trim(filename1)//" "//trim(filename2), exitstat=status)
     endif

     call fits_open_old(unit0,trim(filename1),status)
     call fits_move_to_next_hdu(unit0,status)
     call fits_get_keyword(unit0,'NAXIS2', nph1, status)

     if (.not. allocated(arr_1D_old)) allocate(arr_1D_old(nph1))
     if (.not. allocated(arr_1D_new)) allocate(arr_1D_new(par%nphotons))
     nph = nph1 + par%nphotons

     !--- column : r
     call fits_get_column_number(unit0,'r',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%rp
     if (associated(allph%rp))       deallocate(allph%rp) 
     if (.not.associated(allph%rp))  allocate(allph%rp(nph))
     allph%rp(1:nph1)     = arr_1D_old(:)
     allph%rp(nph1+1:nph) = arr_1D_new(:)

     !--- column : xfreq1
     call fits_get_column_number(unit0,'xfreq1',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%xfreq1
     if (associated(allph%xfreq1))       deallocate(allph%xfreq1)
     if (.not.associated(allph%xfreq1))  allocate(allph%xfreq1(nph))
     allph%xfreq1(1:nph1)     = arr_1D_old(:)
     allph%xfreq1(nph1+1:nph) = arr_1D_new(:)

     !--- column : xfreq2
     call fits_get_column_number(unit0,'xfreq2',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%xfreq2
     if (associated(allph%xfreq2))       deallocate(allph%xfreq2)
     if (.not.associated(allph%xfreq2))  allocate(allph%xfreq2(nph))
     allph%xfreq2(1:nph1)     = arr_1D_old(:)
     allph%xfreq2(nph1+1:nph) = arr_1D_new(:)

     !--- column : NSCATT_HI
     call fits_get_column_number(unit0,'NSCATT_HI',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%nscatt_HI
     if (associated(allph%nscatt_HI))       deallocate(allph%nscatt_HI)
     if (.not.associated(allph%nscatt_HI))  allocate(allph%nscatt_HI(nph))
     allph%nscatt_HI(1:nph1)     = arr_1D_old(:)
     allph%nscatt_HI(nph1+1:nph) = arr_1D_new(:)

     !--- column : NSCATT_dust
     call fits_get_column_number(unit0,'NSCATT_dust',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%nscatt_dust
     if (associated(allph%nscatt_dust))       deallocate(allph%nscatt_dust)
     if (.not.associated(allph%nscatt_dust))  allocate(allph%nscatt_dust(nph))
     allph%nscatt_dust(1:nph1)     = arr_1D_old(:)
     allph%nscatt_dust(nph1+1:nph) = arr_1D_new(:)

     if (par%use_stokes) then
        !--- column : I
        call fits_get_column_number(unit0,'I',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%I
        if (associated(allph%I))       deallocate(allph%I)
        if (.not.associated(allph%I))  allocate(allph%I(nph))
        allph%I(1:nph1)     = arr_1D_old(:)
        allph%I(nph1+1:nph) = arr_1D_new(:)

        !--- column : Q
        call fits_get_column_number(unit0,'Q',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%Q
        if (associated(allph%Q))       deallocate(allph%Q)
        if (.not.associated(allph%Q))  allocate(allph%Q(nph))
        allph%Q(1:nph1)     = arr_1D_old(:)
        allph%Q(nph1+1:nph) = arr_1D_new(:)

        !--- column : U
        call fits_get_column_number(unit0,'U',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%U
        if (associated(allph%U))       deallocate(allph%U)
        if (.not.associated(allph%U))  allocate(allph%U(nph))
        allph%U(1:nph1)     = arr_1D_old(:)
        allph%U(nph1+1:nph) = arr_1D_new(:)

        !--- column : V
        call fits_get_column_number(unit0,'V',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%V
        if (associated(allph%V))       deallocate(allph%V)
        if (.not.associated(allph%V))  allocate(allph%V(nph))
        allph%V(1:nph1)     = arr_1D_old(:)
        allph%V(nph1+1:nph) = arr_1D_new(:)
     endif

     if (trim(par%source_geometry) /= 'point') then
        !--- column : r0
        call fits_get_column_number(unit0,'r0',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%rp0
        if (associated(allph%rp0))       deallocate(allph%rp0)
        if (.not.associated(allph%rp0))  allocate(allph%rp0(nph))
        allph%rp0(1:nph1)     = arr_1D_old(:)
        allph%rp0(nph1+1:nph) = arr_1D_new(:)
     endif
     call fits_close(unit0,status)

     if (allocated(arr_1D_old)) deallocate(arr_1D_old)
     if (allocated(arr_1D_new)) deallocate(arr_1D_new)
  endif

  call fits_open_new(unit,trim(filename1),status)
  call fits_append_table_column(unit,'r',          allph%rp,         status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'xfreq1',     allph%xfreq1,     status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'xfreq2',     allph%xfreq2,     status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'NSCATT_HI',  allph%nscatt_HI,  status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'NSCATT_dust',allph%nscatt_dust,status,bitpix=par%out_bitpix)
  if (par%use_stokes) then
     call fits_append_table_column(unit,'I',       allph%I,          status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'Q',       allph%Q,          status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'U',       allph%U,          status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'V',       allph%V,          status,bitpix=par%out_bitpix)
  endif
  if (trim(par%source_geometry) /= 'point') then
     call fits_append_table_column(unit,'r0',      allph%rp0,        status,bitpix=par%out_bitpix)
  endif
  call fits_close(unit,status)
  end subroutine write_output_allph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module write_mod
