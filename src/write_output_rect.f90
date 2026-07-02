module write_output_rect
  use define
  use iofile_mod
  use output_sum_rect
  use utility
  implicit none
  !---- this should be accesible within this module.
  character(len=128) :: fname_backup
  private :: fname_backup
  private :: write_output_basic, write_output_peeling_2D, write_output_peeling_3D, write_output_allph
  public  :: write_output_outside
contains
!------------------
!-- 2021-05-20, bug-fixed
!-- 2020-11-02, allph arrays are defined to MPI-3 shared memory.
!-- 2020-10-17, divided the routine into three smaller routines.
!               now the routine saves FITS data for each observer.
!-- 2020-08-30, radial profiles for Stokes parameters are saved in a table format.
!               PEELOFF routines can merge the previous outputs with the current ones.
!-- 2017-07-12, Jout, Jabs, Jin versus xfreq are now saved in a table format.
!-- 2017-06-28, Kwang-il Seon
!------------------
  subroutine write_output_outside(filename,grid)
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
  end subroutine write_output_outside
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine write_output_basic(filename,grid)
  use define
#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
  use octree_mod, only: amr_grid
#endif
  implicit none
  character(len=*), intent(in)    :: filename
  type(grid_type),  intent(inout) :: grid
  !--------------------
  type(io_file_type) :: iofh0, iofh
  integer            :: status=0
  character(len=128) :: filename1, filename2
  logical            :: file_exist, merge_ok
  integer            :: colnum
  real(real64)       :: nph1, nph2, nph_tot
  real(real64)       :: exetime1, exetime
  real(real64)       :: nscatt_gas1, nscatt_dust1, nscatt_tot1
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
        filename2    = trim(fname_backup)//trim(io_file_extension(par%file_format))
        call copy_file(trim(par%out_file), trim(filename2), status)
     endif

     call io_open_old(iofh0,trim(filename),status)
     call io_move_to_next_section(iofh0,status)

     call io_get_keyword(iofh0,'exetime',  exetime1,    status)
     call io_get_keyword(iofh0,'nphotons', nph1,        status)
     call io_get_keyword(iofh0,'Nsc_gas',  nscatt_gas1, status)
     call io_get_keyword(iofh0,'Nsc_dust', nscatt_dust1,status)
     call io_get_keyword(iofh0,'Nsc_tot',  nscatt_tot1, status)
     !--- do not update par%exetime, which will give information for each run, not for all accumulated runs.
     !--- update the number of scatterings to improve statistics (2020.09.06).
     exetime         = par%exetime + exetime1
     nph2            = par%no_photons
     nph_tot         = nph1 + nph2
     par%nscatt_gas  = (nscatt_gas1 *nph1 + par%nscatt_gas *nph2)/nph_tot
     par%nscatt_dust = (nscatt_dust1*nph1 + par%nscatt_dust*nph2)/nph_tot
     par%nscatt_tot  = (nscatt_tot1 *nph1 + par%nscatt_tot *nph2)/nph_tot

     if (trim(par%source_geometry) == 'stellar_illumination' .or. trim(par%source_geometry) == 'point_illumination') then
        call io_get_keyword(iofh0,'Raccept',  raccept1, status)
        call io_get_keyword(iofh0,'fluxfac',  fluxfac1, status)
        par%acceptance_rate  = (raccept1 *nph1 + par%acceptance_rate *nph2)/nph_tot
        par%flux_factor      = (fluxfac1 *nph1 + par%flux_factor     *nph2)/nph_tot
     endif

     if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%xfreq)
     call io_get_column_number(iofh0,'Jout',colnum,status)
     call io_read_table_column(iofh0,colnum,arr_1D,status)
     grid%Jout = (arr_1D * nph1 + grid%Jout * nph2)/nph_tot

     if (par%save_Jabs) then
        call io_get_column_number(iofh0,'Jabs',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D,status)
        grid%Jabs = (arr_1D * nph1 + grid%Jabs * nph2)/nph_tot
     endif
     if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
        call io_get_column_number(iofh0,'Jabs2',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D,status)
        grid%Jabs2 = (arr_1D * nph1 + grid%Jabs2 * nph2)/nph_tot
     endif
     if (par%save_Jin) then
        call io_get_column_number(iofh0,'Jin',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D,status)
        grid%Jin = (arr_1D * nph1 + grid%Jin * nph2)/nph_tot
     endif
     if (allocated(arr_1D)) deallocate(arr_1D)

     !--- Jmu image (escaped spectrum binned by mu)
     if (associated(grid%Jmu)) then
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%Jmu)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        grid%Jmu = (arr_2D * nph1 + grid%Jmu * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif

#ifdef CALCP
     if (associated(grid%Pa)) then
        !--- save P_alapha(x,y,z) : scattering number per atom per photon
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%Pa)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        grid%Pa = (arr_3D * nph1 + grid%Pa * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     if (associated(grid%P1)) then
        !--- save radial or z profile, P_alpha
        if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%P1)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_1D,status)
        grid%P1 = (arr_1D * nph1 + grid%P1 * nph2)/nph_tot
        if (allocated(arr_1D)) deallocate(arr_1D)
     endif

     if (associated(grid%P2)) then
        !--- save cylindrical (r, z) profile, P_alpha
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%P2)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        grid%P2 = (arr_2D * nph1 + grid%P2 * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif
     if (par%use_amr_grid) then
        if (associated(amr_grid%Pa)) then
           if (.not. allocated(arr_1D)) allocate(arr_1D, source=amr_grid%Pa)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_1D,status)
           amr_grid%Pa = (arr_1D * nph1 + amr_grid%Pa * nph2)/nph_tot
           if (allocated(arr_1D)) deallocate(arr_1D)
        endif
        if (associated(amr_grid%P1)) then
           if (.not. allocated(arr_1D)) allocate(arr_1D, source=amr_grid%P1)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_1D,status)
           amr_grid%P1 = (arr_1D * nph1 + amr_grid%P1 * nph2)/nph_tot
           if (allocated(arr_1D)) deallocate(arr_1D)
        endif
        if (associated(amr_grid%P2)) then
           if (.not. allocated(arr_2D)) allocate(arr_2D, source=amr_grid%P2)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_2D,status)
           amr_grid%P2 = (arr_2D * nph1 + amr_grid%P2 * nph2)/nph_tot
           if (allocated(arr_2D)) deallocate(arr_2D)
        endif
     endif
#endif

#ifdef CALCPnew
     if (associated(grid%Pa_new)) then
        !--- save P_alapha(x,y,z) : scattering number per atom per photon
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%Pa_new)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        grid%Pa_new = (arr_3D * nph1 + grid%Pa_new * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     if (associated(grid%P1_new)) then
        !--- save radial or z profile, P_alpha
        if (.not. allocated(arr_1D)) allocate(arr_1D, source=grid%P1_new)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_1D,status)
        grid%P1_new = (arr_1D * nph1 + grid%P1_new * nph2)/nph_tot
        if (allocated(arr_1D)) deallocate(arr_1D)
     endif

     if (associated(grid%P2_new)) then
        !--- save cylidrical (r, z)  profile, P_alpha
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%P2_new)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        grid%P2_new = (arr_2D * nph1 + grid%P2_new * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif
     if (par%use_amr_grid) then
        if (associated(amr_grid%Pa_new)) then
           if (.not. allocated(arr_1D)) allocate(arr_1D, source=amr_grid%Pa_new)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_1D,status)
           amr_grid%Pa_new = (arr_1D * nph1 + amr_grid%Pa_new * nph2)/nph_tot
           if (allocated(arr_1D)) deallocate(arr_1D)
        endif
        if (associated(amr_grid%P1_new)) then
           if (.not. allocated(arr_1D)) allocate(arr_1D, source=amr_grid%P1_new)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_1D,status)
           amr_grid%P1_new = (arr_1D * nph1 + amr_grid%P1_new * nph2)/nph_tot
           if (allocated(arr_1D)) deallocate(arr_1D)
        endif
        if (associated(amr_grid%P2_new)) then
           if (.not. allocated(arr_2D)) allocate(arr_2D, source=amr_grid%P2_new)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_2D,status)
           amr_grid%P2_new = (arr_2D * nph1 + amr_grid%P2_new * nph2)/nph_tot
           if (allocated(arr_2D)) deallocate(arr_2D)
        endif
     endif
#endif

#ifdef CALCJ
     if (associated(grid%J)) then
        !--- save J(nu,x,y,z) : mean intensity spectrum at (x,y,z)
        if (.not. allocated(arr_4D)) allocate(arr_4D, source=grid%J)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_4D,status)
        grid%J = (arr_4D * nph1 + grid%J * nph2)/nph_tot
        if (allocated(arr_4D)) deallocate(arr_4D)
     endif

     if (associated(grid%J1)) then
        !--- save radial or z profile of mean intensity spectrum J(nu)
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=grid%J1)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        grid%J1 = (arr_2D * nph1 + grid%J1 * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif

     if (associated(grid%J2)) then
        !--- save cylindrical (r, z) profile of mean intensity spectrum J(nu)
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=grid%J2)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        grid%J2 = (arr_3D * nph1 + grid%J2 * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif
     if (par%use_amr_grid) then
        if (associated(amr_grid%J)) then
           if (.not. allocated(arr_2D)) allocate(arr_2D, source=amr_grid%J)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_2D,status)
           amr_grid%J = (arr_2D * nph1 + amr_grid%J * nph2)/nph_tot
           if (allocated(arr_2D)) deallocate(arr_2D)
        endif
        if (associated(amr_grid%J1)) then
           if (.not. allocated(arr_2D)) allocate(arr_2D, source=amr_grid%J1)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_2D,status)
           amr_grid%J1 = (arr_2D * nph1 + amr_grid%J1 * nph2)/nph_tot
           if (allocated(arr_2D)) deallocate(arr_2D)
        endif
        if (associated(amr_grid%J2)) then
           if (.not. allocated(arr_3D)) allocate(arr_3D, source=amr_grid%J2)
           call io_move_to_next_section(iofh0,status)
           call io_read_image(iofh0,arr_3D,status)
           amr_grid%J2 = (arr_3D * nph1 + amr_grid%J2 * nph2)/nph_tot
           if (allocated(arr_3D)) deallocate(arr_3D)
        endif
     endif
#endif
     call io_close(iofh0,status)
  endif

  !--- open main FITS file.
  call io_open_new(iofh,trim(filename),status)

  !--- xfreq, (relative) frequency, velocity, and wavelength
  call io_append_table_column(iofh,'Xfreq',   grid%xfreq,   status,bitpix=par%out_bitpix)
  call io_append_table_column(iofh,'velocity',grid%velocity,status,bitpix=par%out_bitpix)
  !call io_append_table_column(iofh,'wavelength', grid%wavelength, status,bitpix=par%out_bitpix)
  !--- wavelength should be saved in 64 bits. (2025.10.08)
  call io_append_table_column(iofh,'wavelength', grid%wavelength, status,bitpix=-64)

  !--- Jout, emerging spectrum
  call io_append_table_column(iofh,'Jout',grid%Jout,status,bitpix=par%out_bitpix)

  !--- Jabs, absorbed spectrum
  if (par%save_Jabs) then
     call io_append_table_column(iofh,'Jabs',grid%Jabs,status,bitpix=par%out_bitpix)
  endif
  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
     call io_append_table_column(iofh,'Jabs2',grid%Jabs2,status,bitpix=par%out_bitpix)
  endif

  !--- Jin, input spectrum
  if (par%save_Jin) then
     call io_append_table_column(iofh,'Jin',grid%Jin,status,bitpix=par%out_bitpix)
  endif

  call io_put_keyword(iofh,'ExeTime',  exetime,            'Excution Time (min)',status)
  call io_put_keyword(iofh,'Nprocs',   mpar%nproc,         'No. of Threads',status)
  call io_put_keyword(iofh,'recoil',   par%recoil,         'recoil',status)
  call io_put_keyword(iofh,'coreskip', par%core_skip,      'coreskip algorithm',status)
  call io_put_keyword(iofh,'xyz_sym',  par%xyz_symmetry,   'xyz_symmetry',status)
  call io_put_keyword(iofh,'xy_per',   par%xy_periodic,    'xy_periodic',status)
  call io_put_keyword(iofh,'save_all', par%save_all,       'save all 3D output?',status)
  call io_put_keyword(iofh,'save_Jin', par%save_Jin,       'save input J?',status)
  call io_put_keyword(iofh,'save_Jab', par%save_Jabs,      'save absorbed J?',status)
  call io_put_keyword(iofh,'nphotons', nph_tot,            'number of photons',status)
  call io_put_keyword(iofh,'taumax',   par%taumax,         'tau_max',status)
  call io_put_keyword(iofh,'tauhomo',  par%tauhomo,        'tau_homo',status)
  call io_put_keyword(iofh,'Ngasmax',  par%N_gasmax,       'N(gas)_max cm^-2',status)
  call io_put_keyword(iofh,'Ngashomo', par%N_gashomo,      'N(gas)_homo cm^-2',status)
  call io_put_keyword(iofh,'temp',     par%temperature,    'temperature (K)',status)
  call io_put_keyword(iofh,'Vexp',     par%Vexp,           'Velocity (km/s)',status)
  call io_put_keyword(iofh,'DGR',      par%DGR,            'Dust-to-Gas Ratio',status)
  call io_put_keyword(iofh,'atau3',    par%atau3,          '(a*tau_max)^(1/3)',status)
  call io_put_keyword(iofh,'voigta',   grid%voigt_amean,   'voigt_a',status)
  call io_put_keyword(iofh,'Xfreq1',   grid%xfreq_min,     'Xfreq_min',status)
  call io_put_keyword(iofh,'Xfreq2',   grid%xfreq_max,     'Xfreq_max',status)
  call io_put_keyword(iofh,'Dxfreq',   grid%dxfreq,        'Dxfreq', status)
  call io_put_keyword(iofh,'Dwave',    grid%dwave,         'Dwavelength (angstrom)', status)
  call io_put_keyword(iofh,'I_unit',   par%intensity_unit, 'Intensity Unit (0:no dimension, 1:cm^-2 A^-1)', status)
  call io_put_keyword(iofh,'Dfreq',    grid%Dfreq_ref,     'Doppler Freq. (Hz)',status)
  call io_put_keyword(iofh,'Nsc_dust', par%nscatt_dust,    'Nscatt_dust/photon',status)
  call io_put_keyword(iofh,'Nsc_gas',  par%nscatt_gas,     'Nscatt_gas/photon',status)
  call io_put_keyword(iofh,'Nsc_tot',  par%nscatt_tot,     'Nscatt_tot/photon',status)
  call io_put_keyword(iofh,'nx',       grid%nx,            'No. of x cells',status)
  call io_put_keyword(iofh,'ny',       grid%ny,            'No. of y cells',status)
  call io_put_keyword(iofh,'nz',       grid%nz,            'No. of z cells',status)
  call io_put_keyword(iofh,'xmax',     par%xmax,           'xmax',status)
  call io_put_keyword(iofh,'ymax',     par%ymax,           'ymax',status)
  call io_put_keyword(iofh,'zmax',     par%zmax,           'zmax',status)
  call io_put_keyword(iofh,'EXTNAME',  'Spectrum',         'Spectrum',status)
#ifdef CALCP
  call io_put_keyword(iofh,'calc_P',   .true.,             'calculated Pa within media?',status)
#else
  call io_put_keyword(iofh,'calc_P',   .false.,            'calculated Pa within media?',status)
#endif
#ifdef CALCPnew
  call io_put_keyword(iofh,'calc_Pnew',   .true.,          'calculated Pa within media?',status)
#else
  call io_put_keyword(iofh,'calc_Pnew',   .false.,         'calculated Pa within media?',status)
#endif
#ifdef CALCJ
  call io_put_keyword(iofh,'calc_J',   .true.,             'calculated J within media?',status)
#else
  call io_put_keyword(iofh,'calc_J',   .false.,            'calculated J within media?',status)
#endif
  if (trim(par%source_geometry) == 'stellar_illumination') then
     call io_put_keyword(iofh,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
     call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     call io_put_keyword(iofh,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
     call io_put_keyword(iofh,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
     call io_put_keyword(iofh,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
  endif
  if (trim(par%source_geometry) == 'point_illumination') then
     call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
  endif
!---------------------------------

  !--- save Jmu(xfreq, mu): escaped spectrum binned by mu = cos(theta_z)
  if (associated(grid%Jmu)) then
     call io_append_image(iofh,grid%Jmu,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Jmu',         'J(xfreq, mu) escaped spectrum vs cos(theta_z)',status)
     call io_put_keyword(iofh,'CTYPE1', 'XFREQ',       'axis 1 = xfreq',status)
     call io_put_keyword(iofh,'CRPIX1', 1.0_wp,        'reference pixel',status)
     call io_put_keyword(iofh,'CRVAL1', grid%xfreq_min + 0.5_wp*grid%dxfreq, 'value at CRPIX1',status)
     call io_put_keyword(iofh,'CDELT1', grid%dxfreq,   'xfreq step per pixel',status)
     call io_put_keyword(iofh,'CTYPE2', 'MU',          'axis 2 = cos(theta_z)',status)
     call io_put_keyword(iofh,'CRPIX2', 1.0_wp,        'reference pixel',status)
     call io_put_keyword(iofh,'CRVAL2', par%mu_min + 0.5_wp*par%dmu, 'mu value at CRPIX2 (bin center)',status)
     call io_put_keyword(iofh,'CDELT2', par%dmu,       'mu step per pixel',status)
     call io_put_keyword(iofh,'nmu',    par%nmu,       'number of mu bins',status)
     call io_put_keyword(iofh,'mu_min', par%mu_min,    'mu range lower edge',status)
     call io_put_keyword(iofh,'dmu',    par%dmu,       'mu bin width',status)
  endif

#ifdef CALCP
  if (associated(grid%Pa)) then
     !--- save P_alapha(x,y,z) : scattering number per atom per photon
     call io_append_image(iofh,grid%Pa,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Pa_3D',  'P_alpha (number of scattering)',status)
  endif

  if (associated(grid%P1)) then
     !--- save radial or z profile, P_alpha
     call io_append_image(iofh,grid%P1,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Pa_1D',  'P_alpha (number of scattering)',status)
  endif

  if (associated(grid%P2)) then
     !--- save cylindrical (r, z) profile, P_alpha
     call io_append_image(iofh,grid%P2,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Pa_2D',  'P_alpha (number of scattering)',status)
  endif
  if (par%use_amr_grid) then
     if (associated(amr_grid%Pa)) then
        call io_append_image(iofh,amr_grid%Pa,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Pa_AMR', 'P_alpha per leaf (number of scattering)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
     if (associated(amr_grid%P1)) then
        call io_append_image(iofh,amr_grid%P1,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Pa_1D',  'P_alpha (number of scattering)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
     if (associated(amr_grid%P2)) then
        call io_append_image(iofh,amr_grid%P2,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Pa_2D',  'P_alpha (number of scattering)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
  endif
#endif

#ifdef CALCPnew
  if (associated(grid%Pa_new)) then
     !--- save P_alapha(x,y,z) : scattering number per atom per photon
     call io_append_image(iofh,grid%Pa_new,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Pa_3D_new',  'P_alpha (new_method, number of scattering)',status)
  endif

  if (associated(grid%P1_new)) then
     !--- save radial or z profile, P_alpha
     call io_append_image(iofh,grid%P1_new,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Pa_1D_new',  'P_alpha (new_method, number of scattering)',status)
  endif

  if (associated(grid%P2_new)) then
     !--- save cylindrical (r, z) profile, P_alpha
     call io_append_image(iofh,grid%P2_new,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Pa_2D_new',  'P_alpha (new_method, number of scattering)',status)
  endif
  if (par%use_amr_grid) then
     if (associated(amr_grid%Pa_new)) then
        call io_append_image(iofh,amr_grid%Pa_new,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Pa_AMR_new', 'P_alpha per leaf (new_method)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
     if (associated(amr_grid%P1_new)) then
        call io_append_image(iofh,amr_grid%P1_new,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Pa_1D_new',  'P_alpha (new_method, number of scattering)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
     if (associated(amr_grid%P2_new)) then
        call io_append_image(iofh,amr_grid%P2_new,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Pa_2D_new',  'P_alpha (new_method, number of scattering)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
  endif
#endif

#ifdef CALCJ
  if (associated(grid%J)) then
     !--- save J(nu,x,y,z) : mean intensity spectrum at (x,y,z)
     call io_append_image(iofh,grid%J,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Jx_3D','J(x) (mean intensity)',status)
  endif

  if (associated(grid%J1)) then
     !--- save radial or z profile of mean intensity spectrum J(nu)
     call io_append_image(iofh,grid%J1,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Jx_1D','J(x) (mean intensity)',status)
  endif

  if (associated(grid%J2)) then
     !--- save cylindrical (r, z) profile of mean intensity spectrum J(nu)
     call io_append_image(iofh,grid%J2,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Jx_2D','J(x) (mean intensity)',status)
  endif
  if (par%use_amr_grid) then
     if (associated(amr_grid%J)) then
        call io_append_image(iofh,amr_grid%J,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Jx_AMR','J(x) per leaf (mean intensity)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
     if (associated(amr_grid%J1)) then
        call io_append_image(iofh,amr_grid%J1,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Jx_1D','J(x) (mean intensity)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
     if (associated(amr_grid%J2)) then
        call io_append_image(iofh,amr_grid%J2,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Jx_2D','J(x) (mean intensity)',status)
        call put_amr_JPa_axes(iofh,status)
     endif
  endif
#endif

  !--- close the FITS file.
  call io_close(iofh,status)
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
  type(io_file_type) :: iofh0, iofh
  integer            :: status=0
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
  filename1 = trim(get_base_name(filename))//'_obs2D'//trim(filename_end)//trim(io_file_extension(par%file_format))

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
        filename2    = trim(fname_backup)//'_obs2D'//trim(filename_end)//trim(io_file_extension(par%file_format))
        call copy_file(trim(filename1), trim(filename2), status)
     endif

     if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%scatt_2D)
     call io_open_old(iofh0,trim(filename1),status)
     call io_get_keyword(iofh0,'nphotons', nph1, status)
     nph2    = par%no_photons
     nph_tot = nph1 + nph2

     call io_read_image(iofh0,arr_2D,status)
     obs%scatt_2D = (arr_2D * nph1 + obs%scatt_2D * nph2)/nph_tot
     call io_move_to_next_section(iofh0,status)
     call io_read_image(iofh0,arr_2D,status)
     obs%direc_2D = (arr_2D * nph1 + obs%direc_2D * nph2)/nph_tot
     if (allocated(arr_2D)) deallocate(arr_2D)

     if (par%save_direc0) then
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%direc0_2D)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        obs%direc0_2D = (arr_2D * nph1 + obs%direc0_2D * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
     endif

     call io_close(iofh0,status)
  endif

  !--- We will make the radial intensity profiles.
  if (par%save_radial_profile) call make_radial_intensity(grid,obs,use_2D_data=.true.)

  !--- open FITS file for peel-off.
  call io_open_new(iofh,trim(filename1),status)

  !--- write scattered data
  call io_append_image(iofh,obs%scatt_2D,status,bitpix=par%out_bitpix)

  !--- header keyword for 2D image
  cd1_1  = par%dxim
  cd1_2  = 0.0_wp
  cd2_1  = 0.0_wp
  cd2_2  = par%dyim
  crpix1 = (par%nxim+1)/2.0_wp
  crpix2 = (par%nyim+1)/2.0_wp
  crval1 = 0.0_wp
  crval2 = 0.0_wp
  call io_put_keyword(iofh,'EXTNAME','Scattered','J(x,y) (intensity)',status)
  call io_put_keyword(iofh,'EQUINOX' ,equinox,   'Equinox of Ref. Coord.' ,status)
  call io_put_keyword(iofh,'CD1_1'   ,cd1_1  ,   'Degree / Pixel',status)
  call io_put_keyword(iofh,'CD2_1'   ,cd2_1  ,   'Degree / Pixel',status)
  call io_put_keyword(iofh,'CD1_2'   ,cd1_2  ,   'Degree / Pixel',status)
  call io_put_keyword(iofh,'CD2_2'   ,cd2_2  ,   'Degree / Pixel' ,status)
  call io_put_keyword(iofh,'CRPIX1'  ,crpix1 ,   'Reference Pixel in X',status)
  call io_put_keyword(iofh,'CRPIX2'  ,crpix2 ,   'Reference Pixel in Y',status)
  call io_put_keyword(iofh,'CRVAL1'  ,crval1 ,   'R.A. (Degree)',status)
  call io_put_keyword(iofh,'CRVAL2'  ,crval2 ,   'Dec  (Degree)',status)
  call io_put_keyword(iofh,'CTYPE1'  ,'RA--TAN', 'Coordinate Type',status)
  call io_put_keyword(iofh,'CTYPE2'  ,'DEC-TAN', 'Coordinate Type',status)
  call io_put_keyword(iofh,'DISTANCE', par%distance,     'Distance',status)
  call io_put_keyword(iofh,'DISTUNIT', par%distance_unit,'Distance Unit',status)
  call io_put_keyword(iofh,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
  call io_put_keyword(iofh,'nphotons', nph_tot,          'number of photons',status)
  call io_put_keyword(iofh,'alpha',    obs%alpha, 'alpha (degree)',status)
  call io_put_keyword(iofh,'beta',     obs%beta,  'beta (degree)',status)
  call io_put_keyword(iofh,'gamma',    obs%gamma, 'gamma (degree)',status)
  call io_put_keyword(iofh,'obsx',     obs%x,     'Observer X coordinate',status)
  call io_put_keyword(iofh,'obsy',     obs%y,     'Observer Y coordinate',status)
  call io_put_keyword(iofh,'obsz',     obs%z,     'Observer Z coordinate',status)
  call io_put_keyword(iofh,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
  call io_put_keyword(iofh,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
  call io_put_keyword(iofh,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
  if (trim(par%source_geometry) == 'stellar_illumination') then
     call io_put_keyword(iofh,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
     call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     call io_put_keyword(iofh,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
     call io_put_keyword(iofh,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
     call io_put_keyword(iofh,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
  endif
  if (trim(par%source_geometry) == 'point_illumination') then
     call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
  endif

  !--- write direct data
  call io_append_image(iofh,obs%direc_2D,status,bitpix=par%out_bitpix)
  call io_put_keyword(iofh,'EXTNAME','Direct','J(x,y) (intensity)',status)

  if (par%save_direc0) then
     call io_append_image(iofh,obs%direc0_2D,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Direct0','J(x,y) (intensity)',status)
  endif

  if (associated(obs%radial_I)) then
     call io_append_table_column(iofh,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
     call io_append_table_column(iofh,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','RadialI','I(r) (radial profile)',status)
  endif

  !--- close FITS file
  call io_close(iofh,status)

  if (par%use_stokes) then
     !--- open new FITS file
     filename1 = trim(get_base_name(filename))//'_stokes_2D'//trim(filename_end)//trim(io_file_extension(par%file_format))

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
           filename2 = trim(fname_backup)//'_stokes_2D'//trim(filename_end)//trim(io_file_extension(par%file_format))
           call copy_file(trim(filename1), trim(filename2), status)
        endif

        call io_open_old(iofh0,trim(filename1),status)
        if (.not. allocated(arr_2D)) allocate(arr_2D, source=obs%I_2D)
        call io_read_image(iofh0,arr_2D,status)
        obs%I_2D = (arr_2D * nph1 + obs%I_2D * nph2)/nph_tot
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        obs%Q_2D = (arr_2D * nph1 + obs%Q_2D * nph2)/nph_tot
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        obs%U_2D = (arr_2D * nph1 + obs%U_2D * nph2)/nph_tot
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_2D,status)
        obs%V_2D = (arr_2D * nph1 + obs%V_2D * nph2)/nph_tot
        if (allocated(arr_2D)) deallocate(arr_2D)
        call io_close(iofh0,status)
     endif

     !--- Making an average of radial_pol is a non-sense, because pol = sqrt(Q^2 + U^2)/I.
     !--- Hence, we will make the radial profiles for Stokes parameters from the 3D Stokes parameters,
     !---        instead of averaging the pre-existing two radial profiles.
     if (par%save_radial_profile) call make_radial_stokes(grid,obs,use_2D_data=.true.)

     !--- open FITS file for Stokes parameters.
     call io_open_new(iofh,trim(filename1),status)

     !--- write Stokes I data
     call io_append_image(iofh,obs%I_2D,status,bitpix=par%out_bitpix)

     !--- header keyword for 2D image
     cd1_1  = par%dxim
     cd1_2  = 0.0_wp
     cd2_1  = 0.0_wp
     cd2_2  = par%dyim
     crpix1 = (par%nxim+1)/2.0_wp
     crpix2 = (par%nyim+1)/2.0_wp
     crval1 = 0.0_wp
     crval2 = 0.0_wp
     call io_put_keyword(iofh,'EXTNAME','Stokes_I','Stokes I image',status)
     call io_put_keyword(iofh,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
     call io_put_keyword(iofh,'CD1_1'  ,cd1_1  ,   'Degree / Pixel',status)
     call io_put_keyword(iofh,'CD2_1'  ,cd2_1  ,   'Degree / Pixel',status)
     call io_put_keyword(iofh,'CD1_2'  ,cd1_2  ,   'Degree / Pixel',status)
     call io_put_keyword(iofh,'CD2_2'  ,cd2_2  ,   'Degree / Pixel' ,status)
     call io_put_keyword(iofh,'CRPIX1' ,crpix1 ,   'Reference Pixel in X',status)
     call io_put_keyword(iofh,'CRPIX2' ,crpix2 ,   'Reference Pixel in Y',status)
     call io_put_keyword(iofh,'CRVAL1' ,crval1 ,   'R.A. (Degree)',status)
     call io_put_keyword(iofh,'CRVAL2' ,crval2 ,   'Dec  (Degree)',status)
     call io_put_keyword(iofh,'CTYPE1' ,'RA--TAN', 'Coordinate Type',status)
     call io_put_keyword(iofh,'CTYPE2' ,'DEC-TAN', 'Coordinate Type',status)
     call io_put_keyword(iofh,'DISTANCE', par%distance,     'Distance',status)
     call io_put_keyword(iofh,'DISTUNIT', par%distance_unit,'Distance Unit',status)
     call io_put_keyword(iofh,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
     call io_put_keyword(iofh,'nphotons', nph_tot,          'number of photons',status)
     call io_put_keyword(iofh,'alpha',    obs%alpha, 'alpha (degree)',status)
     call io_put_keyword(iofh,'beta',     obs%beta,  'beta (degree)',status)
     call io_put_keyword(iofh,'gamma',    obs%gamma, 'gamma (degree)',status)
     call io_put_keyword(iofh,'obsx',     obs%x,     'Observer X coordinate',status)
     call io_put_keyword(iofh,'obsy',     obs%y,     'Observer Y coordinate',status)
     call io_put_keyword(iofh,'obsz',     obs%z,     'Observer Z coordinate',status)
     call io_put_keyword(iofh,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
     call io_put_keyword(iofh,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
     call io_put_keyword(iofh,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
     if (trim(par%source_geometry) == 'stellar_illumination') then
        call io_put_keyword(iofh,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
        call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
        call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
        call io_put_keyword(iofh,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
        call io_put_keyword(iofh,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
        call io_put_keyword(iofh,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
     endif
     if (trim(par%source_geometry) == 'point_illumination') then
        call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
        call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     endif

     !--- write Stokes Q, U, V data
     call io_append_image(iofh,obs%Q_2D,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Stokes_Q','Stokes Q image',status)
     call io_append_image(iofh,obs%U_2D,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Stokes_U','Stokes U image',status)
     call io_append_image(iofh,obs%V_2D,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Stokes_V','Stokes V image',status)

     !--- radial Stokes profiles
     if (associated(obs%radial_r)) then
        call io_append_table_column(iofh,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'Q',  obs%radial_Q,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'U',  obs%radial_U,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'V',  obs%radial_V,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'pol',obs%radial_pol,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Stokes_radial','Stokes radial profile',status)
     endif

     !--- close the FITS file
     call io_close(iofh,status)
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
  type(io_file_type) :: iofh0, iofh
  integer            :: status=0
  character(len=128) :: filename1, filename2, filename_end
  logical            :: file_exist, merge_ok
  real(real64)       :: nph1, nph2, nph_tot
  real(real64)       :: cd1_1, cd1_2, cd2_1, cd2_2, cd3_3
  real(real64)       :: crpix1, crpix2, crval1, crval2, crpix3, crval3
  integer            :: equinox = 2000
  real(real64), allocatable :: arr_3D(:,:,:)

  if (present(suffix)) then
     filename_end = trim(suffix)
  else
     filename_end = ''
  endif

  !--- Initialize FITS file name.
  filename1 = trim(get_base_name(filename))//'_obs'//trim(filename_end)//trim(io_file_extension(par%file_format))

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
        filename2    = trim(fname_backup)//'_obs'//trim(filename_end)//trim(io_file_extension(par%file_format))
        call copy_file(trim(filename1), trim(filename2), status)
     endif

     if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%scatt)
     call io_open_old(iofh0,trim(filename1),status)
     call io_get_keyword(iofh0,'nphotons', nph1, status)
     nph2    = par%no_photons
     nph_tot = nph1 + nph2

     call io_read_image(iofh0,arr_3D,status)
     obs%scatt = (arr_3D * nph1 + obs%scatt * nph2)/nph_tot
     call io_move_to_next_section(iofh0,status)
     call io_read_image(iofh0,arr_3D,status)
     obs%direc = (arr_3D * nph1 + obs%direc * nph2)/nph_tot
     if (allocated(arr_3D)) deallocate(arr_3D)

     if (par%save_direc0) then
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%direc0)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        obs%direc0 = (arr_3D * nph1 + obs%direc0 * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif
     if (par%save_dust_scattered) then
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%scatt_dust)
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        obs%scatt_dust = (arr_3D * nph1 + obs%scatt_dust * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
     endif

     call io_close(iofh0,status)
  endif

  !--- We will make the radial intensity profile.
  if (par%save_radial_profile) call make_radial_intensity(grid,obs)

  !--- open FITS file for peel-off.
  call io_open_new(iofh,trim(filename1),status)

  !--- write scattered data
  call io_append_image(iofh,obs%scatt,status,bitpix=par%out_bitpix)

  !--- header keywords for 3D spectral cube.
  !    Data layout: obs%scatt(par%nxfreq, par%nxim, par%nyim) ⇒
  !      NAXIS1 = par%nxfreq  (wavelength axis; index 1 holds the largest
  !                            wavelength because xfreq increases with i and
  !                            wavelength = (1 - vtherm*xfreq/c) * lambda0)
  !      NAXIS2 = par%nxim    (RA pixel)
  !      NAXIS3 = par%nyim    (DEC pixel)
  cd1_1  = -grid%dwave             ! Angstrom / pixel (negative: lambda decreases with i)
  cd2_2  =  par%dxim               ! degree / pixel (RA)
  cd3_3  =  par%dyim               ! degree / pixel (DEC)
  crpix1 = 1.0_wp                  ! reference pixel for spectral axis
  crpix2 = (par%nxim+1)/2.0_wp     ! center of RA axis
  crpix3 = (par%nyim+1)/2.0_wp     ! center of DEC axis
  crval1 = grid%wavelength(1)      ! wavelength at pixel 1 (Angstrom)
  crval2 = 0.0_wp
  crval3 = 0.0_wp
  call io_put_keyword(iofh,'EXTNAME','Scattered','J(freq,x,y) (intensity)',status)
  call io_put_keyword(iofh,'EQUINOX' ,equinox,   'Equinox of Ref. Coord.' ,status)
  !--- Spectral axis (NAXIS1)
  call io_put_keyword(iofh,'CTYPE1'  ,'WAVE'    ,'vacuum wavelength (FITS WCS)', status)
  call io_put_keyword(iofh,'CUNIT1'  ,'Angstrom','wavelength unit', status)
  call io_put_keyword(iofh,'CRPIX1'  ,crpix1    ,'reference pixel for spectral axis', status)
  call io_put_keyword(iofh,'CRVAL1'  ,crval1    ,'wavelength at CRPIX1 (Angstrom)', status)
  call io_put_keyword(iofh,'CD1_1'   ,cd1_1     ,'wavelength step (Angstrom/pixel)', status)
  !--- Spatial axes (NAXIS2 = RA, NAXIS3 = DEC)
  call io_put_keyword(iofh,'CTYPE2'  ,'RA--TAN' ,'Coordinate Type',status)
  call io_put_keyword(iofh,'CUNIT2'  ,'deg'     ,'RA unit',status)
  call io_put_keyword(iofh,'CRPIX2'  ,crpix2    ,'Reference Pixel in X',status)
  call io_put_keyword(iofh,'CRVAL2'  ,crval2    ,'R.A. (Degree)',status)
  call io_put_keyword(iofh,'CD2_2'   ,cd2_2     ,'Degree / Pixel',status)
  call io_put_keyword(iofh,'CTYPE3'  ,'DEC-TAN' ,'Coordinate Type',status)
  call io_put_keyword(iofh,'CUNIT3'  ,'deg'     ,'DEC unit',status)
  call io_put_keyword(iofh,'CRPIX3'  ,crpix3    ,'Reference Pixel in Y',status)
  call io_put_keyword(iofh,'CRVAL3'  ,crval3    ,'Dec  (Degree)',status)
  call io_put_keyword(iofh,'CD3_3'   ,cd3_3     ,'Degree / Pixel',status)
  call io_put_keyword(iofh,'DISTANCE', par%distance,     'Distance',status)
  call io_put_keyword(iofh,'DISTUNIT', par%distance_unit,'Distance Unit',status)
  call io_put_keyword(iofh,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
  call io_put_keyword(iofh,'Xfreq1'  , grid%xfreq_min,   'Xfreq_min',status)
  call io_put_keyword(iofh,'Xfreq2'  , grid%xfreq_max,   'Xfreq_max',status)
  call io_put_keyword(iofh,'Dxfreq'  , grid%dxfreq,      'Dxfreq', status)
  call io_put_keyword(iofh,'Dwave',    grid%dwave,       'Dwavelength (angstrom)', status)
  call io_put_keyword(iofh,'I_unit',   par%intensity_unit, 'Intensity Unit (0:no dimension, 1:cm^-2 A^-1)', status)
  call io_put_keyword(iofh,'Dfreq'   , grid%Dfreq_ref,   'Doppler Freq.  (Hz)',status)
  call io_put_keyword(iofh,'nphotons', nph_tot,          'number of photons',status)
  call io_put_keyword(iofh,'alpha',    obs%alpha, 'alpha (degree)',status)
  call io_put_keyword(iofh,'beta',     obs%beta,  'beta (degree)',status)
  call io_put_keyword(iofh,'gamma',    obs%gamma, 'gamma (degree)',status)
  call io_put_keyword(iofh,'obsx',     obs%x,     'Observer X coordinate',status)
  call io_put_keyword(iofh,'obsy',     obs%y,     'Observer Y coordinate',status)
  call io_put_keyword(iofh,'obsz',     obs%z,     'Observer Z coordinate',status)
  call io_put_keyword(iofh,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
  call io_put_keyword(iofh,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
  call io_put_keyword(iofh,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
  if (trim(par%source_geometry) == 'stellar_illumination') then
     call io_put_keyword(iofh,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
     call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     call io_put_keyword(iofh,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
     call io_put_keyword(iofh,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
     call io_put_keyword(iofh,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
  endif
  if (trim(par%source_geometry) == 'point_illumination') then
     call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
     call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
  endif

  call io_put_keyword(iofh,'ScatDust', par%save_dust_scattered,'Save Dust-Scattered Light', status)

  !--- write direct data
  call io_append_image(iofh,obs%direc,status,bitpix=par%out_bitpix)
  call io_put_keyword(iofh,'EXTNAME','Direct','J(freq,x,y) (intensity)',status)

  if (par%save_direc0) then
     call io_append_image(iofh,obs%direc0,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Direct0','J(freq,x,y) (intensity)',status)
  endif

  if (par%save_dust_scattered) then
     call io_append_image(iofh,obs%scatt_dust,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','ScatDust','J(freq,x,y) (intensity)',status)
  endif

  if (associated(obs%radial_I)) then
     call io_append_table_column(iofh,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
     call io_append_table_column(iofh,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','RadialI','I(r) (radial profile)',status)
  endif

  !--- close FITS file
  call io_close(iofh,status)

  if (par%use_stokes) then
     !--- open new FITS file
     filename1 = trim(get_base_name(filename))//'_stokes'//trim(filename_end)//trim(io_file_extension(par%file_format))

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
           filename2 = trim(fname_backup)//'_stokes'//trim(filename_end)//trim(io_file_extension(par%file_format))
           call copy_file(trim(filename1), trim(filename2), status)
        endif

        call io_open_old(iofh0,trim(filename1),status)
        if (.not. allocated(arr_3D)) allocate(arr_3D, source=obs%I)
        call io_read_image(iofh0,arr_3D,status)
        obs%I = (arr_3D * nph1 + obs%I * nph2)/nph_tot
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        obs%Q = (arr_3D * nph1 + obs%Q * nph2)/nph_tot
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        obs%U = (arr_3D * nph1 + obs%U * nph2)/nph_tot
        call io_move_to_next_section(iofh0,status)
        call io_read_image(iofh0,arr_3D,status)
        obs%V = (arr_3D * nph1 + obs%V * nph2)/nph_tot
        if (allocated(arr_3D)) deallocate(arr_3D)
        call io_close(iofh0,status)
     endif

     !--- Making an average of radial_pol is a non-sense, because pol = sqrt(Q^2 + U^2)/I.
     !--- Hence, we will make the radial profiles for Stokes parameters from the 3D Stokes parameters,
     !---        instead of averaging the pre-existing two radial profiles.
     if (par%save_radial_profile) call make_radial_stokes(grid,obs)

     !--- open FITS file for Stokes parameters.
     call io_open_new(iofh,trim(filename1),status)

     !--- write Stokes I data
     call io_append_image(iofh,obs%I,status,bitpix=par%out_bitpix)

     !--- header keywords for 3D Stokes spectral cube. Same data layout as
     !    'Scattered' above: NAXIS1=nxfreq (wavelength), NAXIS2=nxim (RA), NAXIS3=nyim (DEC).
     cd1_1  = -grid%dwave
     cd2_2  =  par%dxim
     cd3_3  =  par%dyim
     crpix1 = 1.0_wp
     crpix2 = (par%nxim+1)/2.0_wp
     crpix3 = (par%nyim+1)/2.0_wp
     crval1 = grid%wavelength(1)
     crval2 = 0.0_wp
     crval3 = 0.0_wp
     call io_put_keyword(iofh,'EXTNAME','Stokes_I','Stokes I image',status)
     call io_put_keyword(iofh,'EQUINOX',equinox,   'Equinox of Ref. Coord.' ,status)
     call io_put_keyword(iofh,'CTYPE1' ,'WAVE'    ,'vacuum wavelength (FITS WCS)', status)
     call io_put_keyword(iofh,'CUNIT1' ,'Angstrom','wavelength unit', status)
     call io_put_keyword(iofh,'CRPIX1' ,crpix1    ,'reference pixel for spectral axis', status)
     call io_put_keyword(iofh,'CRVAL1' ,crval1    ,'wavelength at CRPIX1 (Angstrom)', status)
     call io_put_keyword(iofh,'CD1_1'  ,cd1_1     ,'wavelength step (Angstrom/pixel)', status)
     call io_put_keyword(iofh,'CTYPE2' ,'RA--TAN' ,'Coordinate Type',status)
     call io_put_keyword(iofh,'CUNIT2' ,'deg'     ,'RA unit',status)
     call io_put_keyword(iofh,'CRPIX2' ,crpix2    ,'Reference Pixel in X',status)
     call io_put_keyword(iofh,'CRVAL2' ,crval2    ,'R.A. (Degree)',status)
     call io_put_keyword(iofh,'CD2_2'  ,cd2_2     ,'Degree / Pixel',status)
     call io_put_keyword(iofh,'CTYPE3' ,'DEC-TAN' ,'Coordinate Type',status)
     call io_put_keyword(iofh,'CUNIT3' ,'deg'     ,'DEC unit',status)
     call io_put_keyword(iofh,'CRPIX3' ,crpix3    ,'Reference Pixel in Y',status)
     call io_put_keyword(iofh,'CRVAL3' ,crval3    ,'Dec  (Degree)',status)
     call io_put_keyword(iofh,'CD3_3'  ,cd3_3     ,'Degree / Pixel',status)
     call io_put_keyword(iofh,'DISTANCE', par%distance,     'Distance',status)
     call io_put_keyword(iofh,'DISTUNIT', par%distance_unit,'Distance Unit',status)
     call io_put_keyword(iofh,'DIST_CM',  par%distance2cm,  'Distance Unit (cm)',status)
     call io_put_keyword(iofh,'Xfreq1' ,grid%xfreq_min, 'Xfreq_min',status)
     call io_put_keyword(iofh,'Xfreq2' ,grid%xfreq_max, 'Xfreq_max',status)
     call io_put_keyword(iofh,'Dxfreq' ,grid%dxfreq,    'Dxfreq', status)
     call io_put_keyword(iofh,'Dwave',  grid%dwave,     'Dwavelength (angstrom)', status)
     call io_put_keyword(iofh,'I_unit' ,par%intensity_unit, 'Intensity Unit (0:no dimension, 1:cm^-2 A^-1)', status)
     call io_put_keyword(iofh,'Dfreq'  ,grid%Dfreq_ref, 'Doppler Freq. (Hz)',status)
     call io_put_keyword(iofh,'nphotons', nph_tot,      'number of photons',status)
     call io_put_keyword(iofh,'alpha',    obs%alpha, 'alpha (degree)',status)
     call io_put_keyword(iofh,'beta',     obs%beta,  'beta (degree)',status)
     call io_put_keyword(iofh,'gamma',    obs%gamma, 'gamma (degree)',status)
     call io_put_keyword(iofh,'obsx',     obs%x,     'Observer X coordinate',status)
     call io_put_keyword(iofh,'obsy',     obs%y,     'Observer Y coordinate',status)
     call io_put_keyword(iofh,'obsz',     obs%z,     'Observer Z coordinate',status)
     call io_put_keyword(iofh,'rot_cenx', par%rotation_center_x,  'Rotation Center X coordinate',status)
     call io_put_keyword(iofh,'rot_ceny', par%rotation_center_y,  'Rotation Center Y coordinate',status)
     call io_put_keyword(iofh,'rot_cenz', par%rotation_center_z,  'Rotation Center Z coordinate',status)
     if (trim(par%source_geometry) == 'stellar_illumination') then
        call io_put_keyword(iofh,'LimbDark', par%stellar_limb_darkening,'Limb Darkening Function', status)
        call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
        call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
        call io_put_keyword(iofh,'RSTAR',    par%stellar_radius         ,'stellar radius', status)
        call io_put_keyword(iofh,'REXOSPH',  par%rmax                   ,'exosphere radius', status)
        call io_put_keyword(iofh,'DIST_S2P', par%distance_star_to_planet,'distance between star and planet', status)
     endif
     if (trim(par%source_geometry) == 'point_illumination') then
        call io_put_keyword(iofh,'Raccept',  par%acceptance_rate,'Acceptance Rate', status)
        call io_put_keyword(iofh,'fluxfac',  par%flux_factor,    'Flux (or luminosity) factor', status)
     endif

     !--- write Stokes Q, U, V data
     call io_append_image(iofh,obs%Q,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Stokes_Q','Stokes Q image',status)
     call io_append_image(iofh,obs%U,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Stokes_U','Stokes U image',status)
     call io_append_image(iofh,obs%V,status,bitpix=par%out_bitpix)
     call io_put_keyword(iofh,'EXTNAME','Stokes_V','Stokes V image',status)

     !--- radial Stokes profiles
     if (associated(obs%radial_r)) then
        call io_append_table_column(iofh,'r',  obs%radial_r,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'I',  obs%radial_I,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'Q',  obs%radial_Q,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'U',  obs%radial_U,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'V',  obs%radial_V,  status,bitpix=par%out_bitpix)
        call io_append_table_column(iofh,'pol',obs%radial_pol,status,bitpix=par%out_bitpix)
        call io_put_keyword(iofh,'EXTNAME','Stokes_radial','Stokes radial profile',status)
     endif

     !--- close the FITS file
     call io_close(iofh,status)
  endif
  end subroutine write_output_peeling_3D
  !-------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine write_output_allph(filename)
  use define
  implicit none
  character(len=*),  intent(in) :: filename

  !--------------------
  type(io_file_type) :: iofh0, iofh
  integer            :: status=0
  character(len=128) :: filename1, filename2
  logical            :: file_exist, merge_ok
  integer            :: colnum
  integer            :: nph1, nph
  real(real64), allocatable :: arr_1D_old(:), arr_1D_new(:)
  type(all_photons_type)    :: all_ph

  !--- Initialize FITS file name.
  filename1 = trim(get_base_name(filename))//'_allph'//trim(io_file_extension(par%file_format))

  !--- check the previous FITS output.
  merge_ok = par%out_merge
  if (merge_ok) then
     inquire(file=trim(filename1),exist=file_exist)
     if (.not. file_exist) merge_ok = .false.
  endif

  if (merge_ok) then
     if (par%save_backup) then
        filename2    = trim(fname_backup)//'_allph'//trim(io_file_extension(par%file_format))
        call copy_file(trim(filename1), trim(filename2), status)
     endif

     call io_open_old(iofh0,trim(filename1),status)
     call io_move_to_next_section(iofh0,status)
     call io_get_keyword(iofh0,'NAXIS2', nph1, status)

     if (.not. allocated(arr_1D_old)) allocate(arr_1D_old(nph1))
     if (.not. allocated(arr_1D_new)) allocate(arr_1D_new(par%nphotons))
     nph = nph1 + par%nphotons

     !--- Note that allph arrays are created as MPI-3 shared memories and thus cannot be destroyed with "deallocate."
     !--- We need to create a new allph structure to expand the arrays (2020-11-02).
     !--- destroy_mem routine in memory_mod should be called by all threads, but this module is called by p_rank = 0 (2020-11-08).
     !--- column : r
     call io_get_column_number(iofh0,'r',colnum,status)
     call io_read_table_column(iofh0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%rp
     call create_mem(all_ph%rp, [nph])
     all_ph%rp(1:nph1)     = arr_1D_old(:)
     all_ph%rp(nph1+1:nph) = arr_1D_new(:)

     !--- column : xfreq1
     call io_get_column_number(iofh0,'xfreq1',colnum,status)
     call io_read_table_column(iofh0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%xfreq1
     call create_mem(all_ph%xfreq1, [nph])
     all_ph%xfreq1(1:nph1)     = arr_1D_old(:)
     all_ph%xfreq1(nph1+1:nph) = arr_1D_new(:)

     !--- column : xfreq2
     call io_get_column_number(iofh0,'xfreq2',colnum,status)
     call io_read_table_column(iofh0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%xfreq2
     call create_mem(all_ph%xfreq2, [nph])
     all_ph%xfreq2(1:nph1)     = arr_1D_old(:)
     all_ph%xfreq2(nph1+1:nph) = arr_1D_new(:)

     !--- column : NSCATT_gas
     call io_get_column_number(iofh0,'NSCATT_gas',colnum,status)
     call io_read_table_column(iofh0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%nscatt_gas
     call create_mem(all_ph%nscatt_gas, [nph])
     all_ph%nscatt_gas(1:nph1)     = arr_1D_old(:)
     all_ph%nscatt_gas(nph1+1:nph) = arr_1D_new(:)

     !--- column : NSCATT_dust
     call io_get_column_number(iofh0,'NSCATT_dust',colnum,status)
     call io_read_table_column(iofh0,colnum,arr_1D_old,status)
     arr_1D_new(:) = allph%nscatt_dust
     call create_mem(all_ph%nscatt_dust, [nph])
     all_ph%nscatt_dust(1:nph1)     = arr_1D_old(:)
     all_ph%nscatt_dust(nph1+1:nph) = arr_1D_new(:)

     if (par%use_stokes) then
        !--- column : I
        call io_get_column_number(iofh0,'I',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%I
        call create_mem(all_ph%I, [nph])
        all_ph%I(1:nph1)     = arr_1D_old(:)
        all_ph%I(nph1+1:nph) = arr_1D_new(:)

        !--- column : Q
        call io_get_column_number(iofh0,'Q',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%Q
        call create_mem(all_ph%Q, [nph])
        all_ph%Q(1:nph1)     = arr_1D_old(:)
        all_ph%Q(nph1+1:nph) = arr_1D_new(:)

        !--- column : U
        call io_get_column_number(iofh0,'U',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%U
        call create_mem(all_ph%U, [nph])
        all_ph%U(1:nph1)     = arr_1D_old(:)
        all_ph%U(nph1+1:nph) = arr_1D_new(:)

        !--- column : V
        call io_get_column_number(iofh0,'V',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%V
        call create_mem(all_ph%V, [nph])
        all_ph%V(1:nph1)     = arr_1D_old(:)
        all_ph%V(nph1+1:nph) = arr_1D_new(:)
     endif

     if (trim(par%source_geometry) /= 'point') then
        !--- column : r0
        call io_get_column_number(iofh0,'r0',colnum,status)
        call io_read_table_column(iofh0,colnum,arr_1D_old,status)
        arr_1D_new(:) = allph%rp0
        call create_mem(all_ph%rp0, [nph])
        all_ph%rp0(1:nph1)     = arr_1D_old(:)
        all_ph%rp0(nph1+1:nph) = arr_1D_new(:)
     endif
     call io_close(iofh0,status)

     if (allocated(arr_1D_old)) deallocate(arr_1D_old)
     if (allocated(arr_1D_new)) deallocate(arr_1D_new)
  else
     all_ph%rp          => allph%rp
     all_ph%xfreq1      => allph%xfreq1
     all_ph%xfreq2      => allph%xfreq2
     all_ph%nscatt_gas  => allph%nscatt_gas
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

  call io_open_new(iofh,trim(filename1),status)
  call io_append_table_column(iofh,'r',          all_ph%rp,         status,bitpix=par%out_bitpix)
  call io_append_table_column(iofh,'xfreq1',     all_ph%xfreq1,     status,bitpix=par%out_bitpix)
  call io_append_table_column(iofh,'xfreq2',     all_ph%xfreq2,     status,bitpix=par%out_bitpix)
  call io_append_table_column(iofh,'NSCATT_gas', all_ph%nscatt_gas, status,bitpix=par%out_bitpix)
  call io_append_table_column(iofh,'NSCATT_dust',all_ph%nscatt_dust,status,bitpix=par%out_bitpix)
  if (par%use_stokes) then
     call io_append_table_column(iofh,'I',       all_ph%I,          status,bitpix=par%out_bitpix)
     call io_append_table_column(iofh,'Q',       all_ph%Q,          status,bitpix=par%out_bitpix)
     call io_append_table_column(iofh,'U',       all_ph%U,          status,bitpix=par%out_bitpix)
     call io_append_table_column(iofh,'V',       all_ph%V,          status,bitpix=par%out_bitpix)
  endif
  if (trim(par%source_geometry) /= 'point') then
     call io_append_table_column(iofh,'r0',      all_ph%rp0,        status,bitpix=par%out_bitpix)
  endif
  call io_close(iofh,status)
  end subroutine write_output_allph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--- Write the bin-axis keywords for the AMR CALC* sections so that readers
!--- can reconstruct the radial / z bin grids without the input AMR file.
  subroutine put_amr_JPa_axes(iofh,status)
  use octree_mod, only: amr_grid
  type(io_file_type), intent(inout) :: iofh
  integer,            intent(inout) :: status
  call io_put_keyword(iofh,'geom_JPa', amr_grid%geometry_JPa, 'JPa geometry (3/2/1/-1)',status)
  select case (amr_grid%geometry_JPa)
  case (3)
     call io_put_keyword(iofh,'nleaf', amr_grid%nleaf, 'number of AMR leaves',status)
  case (2)
     call io_put_keyword(iofh,'nr',   amr_grid%nr_JPa,   'number of radial bins',status)
     call io_put_keyword(iofh,'rmax', amr_grid%rmax_JPa, 'outer radius of radial bins',status)
     call io_put_keyword(iofh,'dr',   amr_grid%dr_JPa,   'radial bin width',status)
     call io_put_keyword(iofh,'nz',   amr_grid%nz_JPa,   'number of z bins',status)
     call io_put_keyword(iofh,'zmin', amr_grid%zmin,     'lower edge of z bins',status)
     call io_put_keyword(iofh,'dz',   amr_grid%dz_JPa,   'z bin width',status)
  case (-1)
     call io_put_keyword(iofh,'nz',   amr_grid%nz_JPa, 'number of z bins',status)
     call io_put_keyword(iofh,'zmin', amr_grid%zmin,   'lower edge of z bins',status)
     call io_put_keyword(iofh,'dz',   amr_grid%dz_JPa, 'z bin width',status)
  case default
     call io_put_keyword(iofh,'nr',   amr_grid%nr_JPa,   'number of radial bins',status)
     call io_put_keyword(iofh,'rmax', amr_grid%rmax_JPa, 'outer radius of radial bins',status)
     call io_put_keyword(iofh,'dr',   amr_grid%dr_JPa,   'radial bin width',status)
  end select
  end subroutine put_amr_JPa_axes
#endif

end module write_output_rect
