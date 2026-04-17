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
!-- 2020-11-02, allph arrays are defined to MPI-3 shared memory.
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
  integer          :: nobs_save, k
  character(len=2) :: filename_end

  call write_output_basic(trim(filename),grid)

#ifdef PEELINGOFF
  if (par%peel_average) then
     nobs_save = 1
  else
     nobs_save = par%nobs
  endif
  do k = 1, nobs_save
     if (nobs_save == 1) then
        filename_end = ''
     else
        write(filename_end,'(i2.2)') k
     endif
     call write_output_peeling(trim(filename),grid,observer(k), suffix=trim(filename_end))
  enddo
#endif
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
  real(real64)       :: nphotons1, nphotons2, nphotons
  real(real64)       :: exetime1, exetime
  real(real64)       :: nscatt_HI1, nscatt_dust1, nscatt_tot1
  real(real64), allocatable :: arr_1D(:), arr_2D(:,:), arr_3D(:,:,:), arr_4D(:,:,:,:)

  !--- check the previous FITS output.
  merge_ok = par%out_merge
  if (merge_ok) then
     inquire(file=trim(filename),exist=file_exist)
     if (.not. file_exist) merge_ok = .false.
  endif

  exetime  = par%exetime
  nphotons = par%no_photons
  if (merge_ok) then
     !--- make a backup file.
     if (par%save_backup) then
        fname_backup = name_for_backup(get_base_name(trim(filename)))
        filename2    = trim(fname_backup)//'.fits.gz'
        !--- copy_file modifies the time stamp of the file.
        !--- call copy_file(trim(par%out_file), trim(fname_backup), status)
        call execute_command_line("cp -p "//trim(par%out_file)//" "//trim(filename2), exitstat=status)
     endif

     call fits_open_old(unit0,trim(filename),status)
     call fits_move_to_next_hdu(unit0,status)

     call fits_get_keyword(unit0,'exetime',  exetime1,    status)
     call fits_get_keyword(unit0,'nphotons', nphotons1,   status)
     call fits_get_keyword(unit0,'Nsc_HI',   nscatt_HI1,  status)
     call fits_get_keyword(unit0,'Nsc_dust', nscatt_dust1,status)
     call fits_get_keyword(unit0,'Nsc_tot',  nscatt_tot1, status)
     !--- do not update par%exetime, which will give information for each run, not for all accumulated runs.
     !--- update the number of scatterings to improve statistics (2020.09.06).
     exetime         = par%exetime + exetime1
     nphotons2       = par%no_photons
     nphotons        = nphotons1 + nphotons2
     par%nscatt_HI   = (nscatt_HI1  *nphotons1 + par%nscatt_HI  *nphotons2)/nphotons
     par%nscatt_dust = (nscatt_dust1*nphotons1 + par%nscatt_dust*nphotons2)/nphotons
     par%nscatt_tot  = (nscatt_tot1 *nphotons1 + par%nscatt_tot *nphotons2)/nphotons

     colnum = 2
     if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%xfreq)
     call fits_read_table_column(unit0,colnum,arr_1D,status)
     grid%Jout = (arr_1D * nphotons1 + grid%Jout * nphotons2)/nphotons

     if (par%save_Jabs) then
        colnum = colnum + 1
        call fits_read_table_column(unit0,colnum,arr_1D,status)
        grid%Jabs = (arr_1D * nphotons1 + grid%Jabs * nphotons2)/nphotons
     endif
     if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
        colnum = colnum + 1
        call fits_read_table_column(unit0,colnum,arr_1D,status)
        grid%Jabs2 = (arr_1D * nphotons1 + grid%Jabs2 * nphotons2)/nphotons
     endif
     if (par%save_Jin) then
        colnum = colnum + 1
        call fits_read_table_column(unit0,colnum,arr_1D,status)
        grid%Jin = (arr_1D * nphotons1 + grid%Jin * nphotons2)/nphotons
     endif
     if (allocated(arr_1D)) deallocate(arr_1D)

#ifdef CALCP
     !if (par%save_all) then
     if (associated(grid%Pa)) then
        !--- save P_alapha(x,y,z) : scattering number per atom per photon
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%Pa)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        grid%Pa = (arr_3D * nphotons1 + grid%Pa * nphotons2)/nphotons
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     if (associated(grid%P1)) then
        !--- save radial or z profile, P_alpha
        if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%P1)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_1D,status)
        grid%P1 = (arr_1D * nphotons1 + grid%P1 * nphotons2)/nphotons
        if (allocated(arr_1D)) deallocate(arr_1D)
     endif

     if (associated(grid%P2)) then
        !--- save cylindrical (r, z) profile, P_alpha
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%P2)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        grid%P2 = (arr_2D * nphotons1 + grid%P2 * nphotons2)/nphotons
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif
#endif

#ifdef CALCPnew
     !if (par%save_all) then
     if (associated(grid%Pa_new)) then
        !--- save P_alapha(x,y,z) : scattering number per atom per photon
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%Pa_new)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        grid%Pa_new = (arr_3D * nphotons1 + grid%Pa_new * nphotons2)/nphotons
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     if (associated(grid%P1_new)) then
        !--- save radial or z profile, P_alpha
        if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%P1_new)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_1D,status)
        grid%P1_new = (arr_1D * nphotons1 + grid%P1_new * nphotons2)/nphotons
        if (allocated(arr_1D)) deallocate(arr_1D)
     endif

     if (associated(grid%P2_new)) then
        !--- save cylidrical (r, z)  profile, P_alpha
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%P2_new)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        grid%P2_new = (arr_2D * nphotons1 + grid%P2_new * nphotons2)/nphotons
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif
#endif

#ifdef CALCJ
     !if (par%save_all) then
     if (associated(grid%J)) then
        !--- save J(nu,x,y,z) : mean intensity spectrum at (x,y,z)
        if (.not. allocated(arr_4D)) allocate(arr_4D, source=grid%J)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_4D,status)
        grid%J = (arr_4D * nphotons1 + grid%J * nphotons2)/nphotons
        if (allocated(arr_4D)) deallocate(arr_4D)
     endif

     if (associated(grid%J1)) then
        !--- save radial or z profile of mean intensity spectrum J(nu)
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%J1)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        grid%J1 = (arr_2D * nphotons1 + grid%J1 * nphotons2)/nphotons
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif

     if (associated(grid%J2)) then
        !--- save cylindrical (r, z) profile of mean intensity spectrum J(nu)
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%J2)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_3D,status)
        grid%J2 = (arr_3D * nphotons1 + grid%J2 * nphotons2)/nphotons
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
  call fits_append_table_column(unit,'dlambda', grid%dlambda, status,bitpix=par%out_bitpix)

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
  call fits_put_keyword(unit,'Nprocs',   mpar%nproc,         'No. of Threads',status)
  call fits_put_keyword(unit,'recoil',   par%recoil,         'recoil',status)
  call fits_put_keyword(unit,'coreskip', par%core_skip,      'coreskip algorithm',status)
  call fits_put_keyword(unit,'xyz_sym',  par%xyz_symmetry,   'xyz_symmetry',status)
  call fits_put_keyword(unit,'xy_per',   par%xy_periodic,    'xy_periodic',status)
  call fits_put_keyword(unit,'save_all', par%save_all,       'save all 3D output?',status)
  call fits_put_keyword(unit,'save_Jin', par%save_Jin,       'save input J?',status)
  call fits_put_keyword(unit,'save_Jab', par%save_Jabs,      'save absorbed J?',status)
  call fits_put_keyword(unit,'nphotons', nphotons,           'number of photons',status)
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
!---------------------------------

#ifdef CALCP
  !if (par%save_all) then
  if (associated(grid%Pa)) then
     !--- save P_alapha(x,y,z) : scattering number per atom per photon
     call fits_append_image(unit,grid%Pa,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_3D',  'P_alpha (number of scattering)',status)
  endif

  !if ((grid%nx == grid%ny .and. grid%nx == grid%nz) .or. par%xy_periodic) then
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
  !if (par%save_all) then
  if (associated(grid%Pa_new)) then
     !--- save P_alapha(x,y,z) : scattering number per atom per photon
     call fits_append_image(unit,grid%Pa_new,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Pa_3D',  'P_alpha (new_method, number of scattering)',status)
  endif

  !if ((grid%nx == grid%ny .and. grid%nx == grid%nz) .or. par%xy_periodic) then
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
  !if (par%save_all) then
  if (associated(grid%J)) then
     !--- save J(nu,x,y,z) : mean intensity spectrum at (x,y,z)
     call fits_append_image(unit,grid%J,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Jx_3D','J(x) (mean intensity)',status)
  endif

  !if ((grid%nx == grid%ny .and. grid%nx == grid%nz) .or. par%xy_periodic) then
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
#ifdef PEELINGOFF
  subroutine write_output_peeling(filename,grid,obs,suffix)
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
  real(real64)       :: nphotons1, nphotons2, nphotons
  real(real64)       :: nph1, nph2, nph_tot
  real(real64)       :: fobs = 1.0, fobs1 = 1.0, fobs2 = 1.0
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
  if (par%peel_average) fobs = par%nobs

  !--- check the previous FITS output.
  merge_ok = par%out_merge
  if (merge_ok) then
     inquire(file=trim(filename1),exist=file_exist)
     if (.not. file_exist) merge_ok = .false.
  endif

  !--- number of photons
  nphotons = par%no_photons

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
     call fits_get_keyword(unit0,'nphotons', nphotons1, status)
     if (status == 202) then
        nphotons1 = par%no_photons
        status    = 0
     endif
     nphotons2 = par%no_photons
     nphotons  = nphotons1 + nphotons2

     !--- average number of observers (2020.10.03)
     !--- FITSIO Error Code 202 means that specified keyword name was not found.
     call fits_get_keyword(unit0,'fobs', fobs1, status)
     if (status == 202) then
        fobs1  = 1.0_wp
        status = 0
     endif
     if (par%peel_average) fobs2 = par%nobs
     fobs    = (fobs1 * nphotons1 + fobs2 * nphotons2) / (nphotons1 + nphotons2)
     nph1    = fobs1 * nphotons1
     nph2    = fobs2 * nphotons2
     nph_tot = nph1  + nph2

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

     if (associated(obs%radial_spec)) then
        call fits_move_to_next_hdu(unit0,status)
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%radial_spec)
        call fits_read_image(unit0,arr_2D,status)
        obs%radial_spec = (arr_2D * nph1 + obs%radial_spec * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif

     call fits_close(unit0,status)
  endif

  !--- open FITS file for peel-off.
  call fits_open_new(unit,trim(filename1),status)

  !--- write scattered data
  call fits_append_image(unit,obs%scatt,status,bitpix=par%out_bitpix)

  !--- header keyword for spectral image
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
  call fits_put_keyword(unit,'Xfreq1'  , grid%xfreq_min,   'Xfreq_min',status)
  call fits_put_keyword(unit,'Xfreq2'  , grid%xfreq_max,   'Xfreq_max',status)
  call fits_put_keyword(unit,'Dxfreq'  , grid%dxfreq,      'Dxfreq', status)
  call fits_put_keyword(unit,'Dfreq'   , grid%Dfreq_ref,   'Doppler Freq.  (Hz)',status)
  call fits_put_keyword(unit,'nphotons', nphotons,         'number of photons',status)
  call fits_put_keyword(unit,'fobs'    , fobs,             'Average # of Observers',status)

  !--- write direct data
  call fits_append_image(unit,obs%direc,status,bitpix=par%out_bitpix)
  call fits_put_keyword(unit,'EXTNAME','Direct','J(freq,x,y) (intensity)',status)

  if (par%save_direc0) then
     call fits_append_image(unit,obs%direc0,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','Direct0','J(freq,x,y) (intensity)',status)
  endif

  if (associated(obs%radial_spec)) then
     call fits_append_image(unit,obs%radial_spec,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','RadialSpec','J(freq,r) (radial spectrum)',status)
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
        if (associated(obs%radial_r)) then
           !--- Making an average of radial_pol is a non-sense, because pol = sqrt(Q^2 + U^2)/I.
           !--- Hence, we will make the radial profiles for Stokes parameters from the 3D Stokes parameters,
           !---        instead of averaging the pre-existing two radial profiles.
           call make_radial_stokes(grid,obs)
        endif
        call fits_close(unit0,status)
     endif

     !--- open FITS file for Stokes parameters.
     call fits_open_new(unit,trim(filename1),status)

     !--- write Stokes I data
     call fits_append_image(unit,obs%I,status,bitpix=par%out_bitpix)

     !--- header keyword for spectral image
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
     call fits_put_keyword(unit,'Xfreq1' ,grid%xfreq_min, 'Xfreq_min',status)
     call fits_put_keyword(unit,'Xfreq2' ,grid%xfreq_max, 'Xfreq_max',status)
     call fits_put_keyword(unit,'Dxfreq' ,grid%dxfreq,    'Dxfreq', status)
     call fits_put_keyword(unit,'Dfreq'  ,grid%Dfreq_ref, 'Doppler Freq. (Hz)',status)
     call fits_put_keyword(unit,'nphotons',nphotons,       'number of photons',status)
     call fits_put_keyword(unit,'fobs'    ,fobs,           'Average # of Observers',status)

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

  if (par%save_sightline_tau) then
     !--- Initialize FITS file name.
     filename1 = trim(get_base_name(filename))//'_tau'//trim(filename_end)//'.fits.gz'

     !--- check the previous FITS output.
     merge_ok = par%out_merge
     if (merge_ok) then
        inquire(file=trim(filename1),exist=file_exist)
        if (.not. file_exist) merge_ok = .false.
     endif

     if (merge_ok) then
        !--- make a backup file.
        if (par%save_backup) then
           filename2    = trim(fname_backup)//'_tau'//trim(filename_end)//'.fits.gz'
           !call copy_file(trim(filename1), trim(filename2), status)
           call execute_command_line("cp -p "//trim(filename1)//" "//trim(filename2), exitstat=status)
        endif
 
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%tau_HI)
        call fits_open_old(unit0,trim(filename1),status)
        call fits_read_image(unit0,arr_3D,status)
        obs%tau_HI = (arr_3D * nph1 + obs%tau_HI * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)

        if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%N_HI)
        call fits_move_to_next_hdu(unit0,status)
        call fits_read_image(unit0,arr_2D,status)
        obs%N_HI = (arr_2D * nph1 + obs%N_HI * nph2)/nph_tot

        if (par%DGR > 0.0_wp) then
           call fits_move_to_next_hdu(unit0,status)
           call fits_read_image(unit0,arr_2D,status)
           obs%tau_dust = (arr_2D * nph1 + obs%tau_dust * nph2)/nph_tot
        endif

        if (allocated(arr_2D)) deallocate(arr_2D)
        call fits_close(unit0,status)
     endif

     call fits_open_new(unit,trim(filename1),status)

     !--- write Images for optical depths
     call fits_append_image(unit,obs%tau_HI,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','TAU_HI','HI optical depth',status)

     !--- write keywords
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

     !--- write Images for optical depths
     call fits_append_image(unit,obs%N_HI,status,bitpix=par%out_bitpix)
     call fits_put_keyword(unit,'EXTNAME','N_HI','HI column density',status)

     !--- write Images for optical depths
     if (par%DGR > 0.0_wp) then
        call fits_append_image(unit,obs%tau_dust,status,bitpix=par%out_bitpix)
        call fits_put_keyword(unit,'EXTNAME','TAU_dust','dust column density',status)
     endif

     !--- close fits file
     call fits_close(unit,status)
  endif
  end subroutine write_output_peeling
#endif
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
  type(all_photons_type)    :: all_ph

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

     !--- Note that allph arrays are created as MPI-3 shared memories and thus cannot be destroyed with "deallocate."
     !--- We need to create a new allph structure to expand the arrays (2020-11-02).
     !--- destroy_mem routine in memory_mod should be called by all threads, but this module is called by p_rank = 0 (2020-11-08).
     !--- column : r
     call fits_get_column_number(unit0,'r',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%rp
     call create_mem(all_ph%rp, [nph])
     all_ph%rp(1:nph1)     = arr_1D_old(:)
     all_ph%rp(nph1+1:nph) = arr_1D_new(:)

     !--- column : xfreq1
     call fits_get_column_number(unit0,'xfreq1',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%xfreq1
     call create_mem(all_ph%xfreq1, [nph])
     all_ph%xfreq1(1:nph1)     = arr_1D_old(:)
     all_ph%xfreq1(nph1+1:nph) = arr_1D_new(:)

     !--- column : xfreq2
     call fits_get_column_number(unit0,'xfreq2',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%xfreq2
     call create_mem(all_ph%xfreq2, [nph])
     all_ph%xfreq2(1:nph1)     = arr_1D_old(:)
     all_ph%xfreq2(nph1+1:nph) = arr_1D_new(:)

     !--- column : NSCATT_HI
     call fits_get_column_number(unit0,'NSCATT_HI',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%nscatt_HI
     call create_mem(all_ph%nscatt_HI, [nph])
     all_ph%nscatt_HI(1:nph1)     = arr_1D_old(:)
     all_ph%nscatt_HI(nph1+1:nph) = arr_1D_new(:)

     !--- column : NSCATT_dust
     call fits_get_column_number(unit0,'NSCATT_dust',colnum,status)
     call fits_read_table_column(unit0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%nscatt_dust
     call create_mem(all_ph%nscatt_dust, [nph])
     all_ph%nscatt_dust(1:nph1)     = arr_1D_old(:)
     all_ph%nscatt_dust(nph1+1:nph) = arr_1D_new(:)

     if (par%use_stokes) then
        !--- column : I
        call fits_get_column_number(unit0,'I',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%I
        call create_mem(all_ph%I, [nph])
        all_ph%I(1:nph1)     = arr_1D_old(:)
        all_ph%I(nph1+1:nph) = arr_1D_new(:)

        !--- column : Q
        call fits_get_column_number(unit0,'Q',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%Q
        call create_mem(all_ph%Q, [nph])
        all_ph%Q(1:nph1)     = arr_1D_old(:)
        all_ph%Q(nph1+1:nph) = arr_1D_new(:)

        !--- column : U
        call fits_get_column_number(unit0,'U',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%U
        call create_mem(all_ph%U, [nph])
        all_ph%U(1:nph1)     = arr_1D_old(:)
        all_ph%U(nph1+1:nph) = arr_1D_new(:)

        !--- column : V
        call fits_get_column_number(unit0,'V',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%V
        call create_mem(all_ph%V, [nph])
        all_ph%V(1:nph1)     = arr_1D_old(:)
        all_ph%V(nph1+1:nph) = arr_1D_new(:)
     endif

     if (trim(par%source_geometry) /= 'point') then
        !--- column : r0
        call fits_get_column_number(unit0,'r0',colnum,status)
        call fits_read_table_column(unit0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%rp0
        call create_mem(all_ph%rp0, [nph])
        all_ph%rp0(1:nph1)     = arr_1D_old(:)
        all_ph%rp0(nph1+1:nph) = arr_1D_new(:)
     endif
     call fits_close(unit0,status)

     if (allocated(arr_1D_old)) deallocate(arr_1D_old)
     if (allocated(arr_1D_new)) deallocate(arr_1D_new)
  else
     all_ph%rp          => allph%rp
     all_ph%xfreq1      => allph%xfreq1
     all_ph%xfreq2      => allph%xfreq2
     all_ph%nscatt_HI   => allph%nscatt_HI
     all_ph%nscatt_dust => allph%nscatt_dust
     if (par%use_stokes) then
        all_ph%I => allph%I
        all_ph%Q => allph%Q
        all_ph%U => allph%U
        all_ph%V => allph%V
     endif
     if (trim(par%source_geometry) /= 'point') then
        all_ph%rp0 => allph%rp0
     endif
  endif

  call fits_open_new(unit,trim(filename1),status)
  call fits_append_table_column(unit,'r',          all_ph%rp,         status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'xfreq1',     all_ph%xfreq1,     status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'xfreq2',     all_ph%xfreq2,     status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'NSCATT_HI',  all_ph%nscatt_HI,  status,bitpix=par%out_bitpix)
  call fits_append_table_column(unit,'NSCATT_dust',all_ph%nscatt_dust,status,bitpix=par%out_bitpix)
  if (par%use_stokes) then
     call fits_append_table_column(unit,'I',       all_ph%I,          status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'Q',       all_ph%Q,          status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'U',       all_ph%U,          status,bitpix=par%out_bitpix)
     call fits_append_table_column(unit,'V',       all_ph%V,          status,bitpix=par%out_bitpix)
  endif
  if (trim(par%source_geometry) /= 'point') then
     call fits_append_table_column(unit,'r0',      all_ph%rp0,        status,bitpix=par%out_bitpix)
  endif
  call fits_close(unit,status)
  end subroutine write_output_allph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module write_mod
