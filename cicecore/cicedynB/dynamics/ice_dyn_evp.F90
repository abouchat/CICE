!=======================================================================
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (submitted 2013).  The 
! revised elastic-viscous-plastic method.  Ocean Modelling.
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module ice_dyn_evp

      use ice_kinds_mod
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, c4, p027, p055, p111, p166, &
          p2, p222, p25, p333, p5, c1, c6, c2
      use ice_dyn_shared, only: stepu, dyn_prep1, dyn_prep2, dyn_finish, &
          ndte, yield_curve, ecci, denom1, arlx1i, fcor_blk, uvel_init,  &
          vvel_init, basal_stress_coeff, basalstress, Ktens, revp
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength, icepack_query_parameters

      implicit none
      private
      public :: evp

!=======================================================================

      contains

!=======================================================================

! Elastic-viscous-plastic dynamics driver
!
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied 
! rotated on u grid and multiplied by aice.  strairxT = 0 in this 
! case so operations in dyn_prep1 are pointless but carried out to 
! minimise code changes.
#endif
!
! author: Elizabeth C. Hunke, LANL

      subroutine evp (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks, ncat, nx_global, ny_global
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT, &
          strairx, strairy, uocn, vocn, ss_tltx, ss_tlty, iceumask, fm, &
          strtltx, strtlty, strocnx, strocny, strintx, strinty, taubx, tauby, &
          strocnxT, strocnyT, strax, stray, &
          Tbu, hwater, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, uarear, tinyarea, to_ugrid, t2ugrid_vector, u2tgrid_vector, &
          grid_type, HTE, HTN
      use ice_state, only: aice, vice, vsno, uvel, vvel, divu, shear, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop, timer_evp_1d, timer_evp_2d
      use ice_dyn_evp_1d
      use ice_dyn_shared, only: evp_kernel_ver, yield_curve

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: & 
         ksub           , & ! subcycle step
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      integer (kind=int_kind), dimension(max_blocks) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         tmass    , & ! total mass of ice and snow (kg/m^2)
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         forcex   , & ! work array: combined atm stress and ocn tilt, x
         forcey   , & ! work array: combined atm stress and ocn tilt, y
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         strtmp       ! stress combinations for momentum equation

      logical (kind=log_kind) :: calc_strair

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask, &  ! ice extent mask (T-cell)
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block
      
      character(len=*), parameter :: subname = '(evp)'

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(fld2(nx_block,ny_block,2,max_blocks))

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !-----------------------------------------------------------------
      ! boundary updates
      ! commented out because the ghost cells are freshly 
      ! updated after cleanup_itd
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (aice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vsno,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         do j = 1, ny_block 
         do i = 1, nx_block 
            rdg_conv (i,j,iblk) = c0 
            rdg_shear(i,j,iblk) = c0 
            divu (i,j,iblk) = c0 
            shear(i,j,iblk) = c0 
         enddo
         enddo

      !-----------------------------------------------------------------
      ! preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep1 (nx_block,           ny_block,           & 
                         ilo, ihi,           jlo, jhi,           &
                         aice    (:,:,iblk), vice    (:,:,iblk), & 
                         vsno    (:,:,iblk), tmask   (:,:,iblk), & 
                         strairxT(:,:,iblk), strairyT(:,:,iblk), & 
                         strairx (:,:,iblk), strairy (:,:,iblk), & 
                         tmass   (:,:,iblk), icetmask(:,:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (icetmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call to_ugrid(tmass,umass)
      call to_ugrid(aice_init, aiu)

      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO or other forcing
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------
      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. calc_strair) then       
         strairx(:,:,:) = strax(:,:,:)
         strairy(:,:,:) = stray(:,:,:)
      else
         call t2ugrid_vector(strairx)
         call t2ugrid_vector(strairy)
      endif      

! tcraig, tcx, threading here leads to some non-reproducbile results and failures in icepack_ice_strength
! need to do more debugging
      !$TCXOMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep2 (nx_block,             ny_block,             & 
                         ilo, ihi,             jlo, jhi,             &
                         icellt(iblk),         icellu(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         indxui      (:,iblk), indxuj      (:,iblk), & 
                         aiu       (:,:,iblk), umass     (:,:,iblk), & 
                         umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), & 
                         umask     (:,:,iblk),                       & 
                         uocn      (:,:,iblk), vocn      (:,:,iblk), & 
                         strairx   (:,:,iblk), strairy   (:,:,iblk), & 
                         ss_tltx   (:,:,iblk), ss_tlty   (:,:,iblk), &  
                         icetmask  (:,:,iblk), iceumask  (:,:,iblk), & 
                         fm        (:,:,iblk), dt,                   & 
                         strtltx   (:,:,iblk), strtlty   (:,:,iblk), & 
                         strocnx   (:,:,iblk), strocny   (:,:,iblk), & 
                         strintx   (:,:,iblk), strinty   (:,:,iblk), & 
                         taubx     (:,:,iblk), tauby     (:,:,iblk), & 
                         waterx    (:,:,iblk), watery    (:,:,iblk), & 
                         forcex    (:,:,iblk), forcey    (:,:,iblk), & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         uvel_init (:,:,iblk), vvel_init (:,:,iblk), &
                         uvel      (:,:,iblk), vvel      (:,:,iblk), &
                         Tbu       (:,:,iblk))

      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------

         strength(:,:,iblk) = c0  ! initialize
         do ij = 1, icellt(iblk)
            i = indxti(ij, iblk)
            j = indxtj(ij, iblk)
            call icepack_ice_strength (ncat,                 &
                                      aice    (i,j,  iblk), & 
                                      vice    (i,j,  iblk), & 
                                      aice0   (i,j,  iblk), & 
                                      aicen   (i,j,:,iblk), &  
                                      vicen   (i,j,:,iblk), & 
                                      strength(i,j,  iblk) )
         enddo  ! ij

         ! load velocity into array for boundary updates
         fld2(:,:,1,iblk) = uvel(:,:,iblk)
         fld2(:,:,2,iblk) = vvel(:,:,iblk)

      enddo  ! iblk
      !$TCXOMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in dyn_prep2
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call ice_timer_stop(timer_bound)

      ! unload
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         uvel(:,:,iblk) = fld2(:,:,1,iblk)
         vvel(:,:,iblk) = fld2(:,:,2,iblk)
      enddo
      !$OMP END PARALLEL DO

      if (maskhalo_dyn) then
         call ice_timer_start(timer_bound)
         halomask = 0
         where (iceumask) halomask = 1
         call ice_HaloUpdate (halomask,          halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

      !-----------------------------------------------------------------
      ! basal stress coefficients (landfast ice)
      !-----------------------------------------------------------------
      
      if (basalstress) then
       !$OMP PARALLEL DO PRIVATE(iblk)
       do iblk = 1, nblocks
         call basal_stress_coeff (nx_block,         ny_block,       &
                                  icellu  (iblk),                   &
                                  indxui(:,iblk),   indxuj(:,iblk), &
                                  vice(:,:,iblk),   aice(:,:,iblk), &
                                  hwater(:,:,iblk), Tbu(:,:,iblk))
       enddo
       !$OMP END PARALLEL DO
      endif
      call ice_timer_start(timer_evp_2d)
      if (evp_kernel_ver > 0) then
        !write(*,*)'Entering evp_kernel version ',evp_kernel_ver
        if (trim(grid_type) == 'tripole') then
          call abort_ice('(ice_dyn_evp): &
             & Kernel not tested on tripole grid. Set evp_kernel_ver=0')
        endif
        call evp_copyin(                                                &
          nx_block,ny_block,nblocks,nx_global+2*nghost,ny_global+2*nghost,&
          HTE,HTN,                                                      &
!v1          dxhy,dyhx,cyp,cxp,cym,cxm,tinyarea,                           &
!v1          waterx,watery,                                                &
          icetmask, iceumask,                                           &
          cdn_ocn,aiu,uocn,vocn,forcex,forcey,Tbu,        &
          umassdti,fm,uarear,tarear,strintx,strinty,uvel_init,vvel_init,&
          strength,uvel,vvel,dxt,dyt,                                   &
          stressp_1 ,stressp_2, stressp_3, stressp_4,                   &
          stressm_1 ,stressm_2, stressm_3, stressm_4,                   &
          stress12_1,stress12_2,stress12_3,stress12_4                   )
        if (evp_kernel_ver == 2) then
          call ice_timer_start(timer_evp_1d)
          call evp_kernel_v2()
          call ice_timer_stop(timer_evp_1d)
!v1        else if (evp_kernel_ver == 1) then
!v1          call evp_kernel_v1()
        else
          write(*,*)'Kernel: evp_kernel_ver = ',evp_kernel_ver
          call abort_ice('(ice_dyn_evp): Kernel not implemented.')
        endif
        call evp_copyout(                                               &
          nx_block,ny_block,nblocks,nx_global+2*nghost,ny_global+2*nghost,&
!strocn          uvel,vvel, strocnx,strocny, strintx,strinty,                  &
          uvel,vvel, strintx,strinty,                                   &
          stressp_1, stressp_2, stressp_3, stressp_4,                   &
          stressm_1, stressm_2, stressm_3, stressm_4,                   &
          stress12_1,stress12_2,stress12_3,stress12_4,                  &
          divu,rdg_conv,rdg_shear,shear,taubx,tauby                     )

      else ! evp_kernel_ver == 0 (Standard CICE)

      do ksub = 1,ndte        ! subcycling

      !-----------------------------------------------------------------
      ! stress tensor equation, total surface stress
      !-----------------------------------------------------------------

         !$TCXOMP PARALLEL DO PRIVATE(iblk,strtmp)
         do iblk = 1, nblocks

!            if (trim(yield_curve) == 'ellipse') then
               call stress (nx_block,             ny_block,             & 
                            ksub,                 icellt(iblk),         & 
                            indxti      (:,iblk), indxtj      (:,iblk), & 
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &     
                            dxt       (:,:,iblk), dyt       (:,:,iblk), & 
                            dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
                            cxp       (:,:,iblk), cyp       (:,:,iblk), & 
                            cxm       (:,:,iblk), cym       (:,:,iblk), & 
                            tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
                            strength  (:,:,iblk),                       & 
                            stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                            stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                            stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                            stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                            stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                            stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                            shear     (:,:,iblk), divu      (:,:,iblk), & 
                            rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
                            strtmp    (:,:,:) )
!            endif               ! yield_curve

      !-----------------------------------------------------------------
      ! momentum equation
      !-----------------------------------------------------------------

            call stepu (nx_block,            ny_block,           &
                        icellu       (iblk), Cdn_ocn (:,:,iblk), & 
                        indxui     (:,iblk), indxuj    (:,iblk), &
                        ksub,                                    &
                        aiu      (:,:,iblk), strtmp  (:,:,:),    & 
                        uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                        waterx   (:,:,iblk), watery  (:,:,iblk), & 
                        forcex   (:,:,iblk), forcey  (:,:,iblk), & 
                        umassdti (:,:,iblk), fm      (:,:,iblk), & 
                        uarear   (:,:,iblk),                     & 
                        strintx  (:,:,iblk), strinty (:,:,iblk), &
                        taubx    (:,:,iblk), tauby   (:,:,iblk), & 
                        uvel_init(:,:,iblk), vvel_init(:,:,iblk),&
                        uvel     (:,:,iblk), vvel    (:,:,iblk), &
                        Tbu      (:,:,iblk))

            ! load velocity into array for boundary updates
            fld2(:,:,1,iblk) = uvel(:,:,iblk)
            fld2(:,:,2,iblk) = vvel(:,:,iblk)
         enddo
         !$TCXOMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld2,               halo_info_mask, &
                                 field_loc_NEcorner, field_type_vector)
         else
            call ice_HaloUpdate (fld2,               halo_info, &
                                 field_loc_NEcorner, field_type_vector)
         endif
         call ice_timer_stop(timer_bound)

         ! unload
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            uvel(:,:,iblk) = fld2(:,:,1,iblk)
            vvel(:,:,iblk) = fld2(:,:,2,iblk)
         enddo
         !$OMP END PARALLEL DO

      enddo                     ! subcycling

!       write (nu_diag,*) 'divune,', divune(37,110)    !#ab#
!       write (nu_diag,*) 'divucapne,', divucapne(37,110)    !#ab#
!       write (nu_diag,*) 'dmin,', dmin(37,110)    !#ab#
!       write (nu_diag,*) 'tinyarea,', tinyarea(37,110)    !#ab#
!       write (nu_diag,*) 'dmin*tinyarea,', dmin*tinyarea(37,110)!#ab#
!
!       write (nu_diag,*) 'divu_out,', divu(37,110)    !#ab#
!       write (nu_diag,*) '1/tarea,', tarear(37,110)    !#ab#
!
!       write (nu_diag,*) 'mshearne,', mshearne(37,110)    !#ab#
!       write (nu_diag,*) 'mshearcapne,', mshearcapne(37,110)    !#ab#
!       write (nu_diag,*) 'smin,', smin*tinyarea(37,110)    !#ab#
!       write (nu_diag,*) 'mshear_out,', shear(37,110)    !#ab#
!
!       write (nu_diag,*) 'Mne,', Mne(37,110)    !#ab#
!       write (nu_diag,*) 'Deltane,', Deltane(37,110)    !#ab#



      endif  ! evp_kernel_ver
      call ice_timer_stop(timer_evp_2d)

      deallocate(fld2)
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      ! Force symmetry across the tripole seam
      if (trim(grid_type) == 'tripole') then
      if (maskhalo_dyn) then
         !-------------------------------------------------------
         ! set halomask to zero because ice_HaloMask always keeps
         ! local copies AND tripole zipper communication
         !-------------------------------------------------------
         halomask = 0
         call ice_HaloMask(halo_info_mask, halo_info, halomask)

         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloDestroy(halo_info_mask)
      else
         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                              field_loc_center,  field_type_scalar)
      endif   ! maskhalo
      endif   ! tripole

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call dyn_finish                               & 
              (nx_block,           ny_block,           & 
               icellu      (iblk), Cdn_ocn (:,:,iblk), & 
               indxui    (:,iblk), indxuj    (:,iblk), & 
               uvel    (:,:,iblk), vvel    (:,:,iblk), & 
               uocn    (:,:,iblk), vocn    (:,:,iblk), & 
               aiu     (:,:,iblk), fm      (:,:,iblk), & 
               strintx (:,:,iblk), strinty (:,:,iblk), &
               strairx (:,:,iblk), strairy (:,:,iblk), &
               strocnx (:,:,iblk), strocny (:,:,iblk), & 
               strocnxT(:,:,iblk), strocnyT(:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

      call u2tgrid_vector(strocnxT)    ! shift
      call u2tgrid_vector(strocnyT)

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine evp

!=======================================================================

! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! author: Elizabeth C. Hunke, LANL

      subroutine stress (nx_block,   ny_block,   & 
                         ksub,       icellt,     & 
                         indxti,     indxtj,     & 
                         uvel,       vvel,       & 
                         dxt,        dyt,        & 
                         dxhy,       dyhx,       & 
                         cxp,        cyp,        & 
                         cxm,        cym,        & 
                         tarear,     tinyarea,   & 
                         strength,               & 
                         stressp_1,  stressp_2,  & 
                         stressp_3,  stressp_4,  & 
                         stressm_1,  stressm_2,  & 
                         stressm_3,  stressm_4,  & 
                         stress12_1, stress12_2, & 
                         stress12_3, stress12_4, & 
                         shear,      divu,       & 
                         rdg_conv,   rdg_shear,  & 
                         str )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(out) :: &
         str          ! stress combinations

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        mshearne ,mshearnw, mshearse, mshearsw    , & ! max. shear rate (eps_II)
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delta
        zetane, zetanw, zetase, zetasw            , & ! bulk viscosity
        divucapne, divucapnw, divucapse, divucapsw        , & ! divergence viscous capping
        mshearcapne ,mshearcapnw, mshearcapse, mshearcapsw, & ! max. shear rate viscous capping
        Deltacapne, Deltacapnw, Deltacapse, Deltacapsw    , & ! Delta viscous capping
        puny                                      , & ! puny
        c0ne, c0nw, c0se, c0sw                    , & ! useful combinations
        c1ne, c1nw, c1se, c1sw                    , &
        Mne, Mnw, Mse, Msw                        , &
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp, tmp, divu_tmp         , &
        tmpne, tmpnw, tmpsw, tmpse, shear_tmp     , &
        P_ne, P_nw, P_se, P_sw                        ! replacement pressure or ice strength, given the regularization

      character(len=*), parameter :: subname = '(stress)'

      real (kind=dbl_kind) :: &
         smin   = 1.0_dbl_kind , &
         dmin   = 1.0_dbl_kind , &
         pi                    , &
         phi                   , &
         sinphi                , &
         delta

      character(len=char_len):: & ! #ab#
         capping  ! 'max' or 'add_min' or 'tanh', capping of the viscous coefficients

      character(len=char_len):: & ! #ab#
         regularization  ! 'RP' or 'noRP', choose if using replacement pressure or not.

      call icepack_query_parameters(pi_out=pi)
      phi    = pi/c6  ! angle of friction = 30 degrees
      sinphi = sin(phi)
      delta  = p166

      capping        = 'max' ! 'max', 'add_min', or 'tanh'
      regularization = 'RP'      ! 'RP', or 'noRP'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
         divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
         divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
         divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
         tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
         tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
         tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
                 -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
                 -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
                 -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
                 -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )
         
         !-----------------------------------------------------------------
         ! stresses for ellipse #ab#
         !-----------------------------------------------------------------           
         if (trim(yield_curve) == 'ellipse') then
         
            ! Delta (in the denominator of zeta, eta)
            Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
            Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
            Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))
            Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))

            ! Viscous capping on Delta
            if (trim(capping) == 'max') then
               Deltacapne = max(Deltane,dmin*tinyarea(i,j))
               Deltacapnw = max(Deltanw,dmin*tinyarea(i,j))
               Deltacapsw = max(Deltasw,dmin*tinyarea(i,j))
               Deltacapse = max(Deltase,dmin*tinyarea(i,j))
            elseif (trim(capping) == 'add_min') then
               Deltacapne = (Deltane+(dmin*tinyarea(i,j)))
               Deltacapnw = (Deltanw+(dmin*tinyarea(i,j)))
               Deltacapsw = (Deltasw+(dmin*tinyarea(i,j)))
               Deltacapse = (Deltase+(dmin*tinyarea(i,j)))
            elseif (trim(capping) == 'tanh') then
               Deltacapne = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/Deltane))
               Deltacapnw = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/Deltanw))
               Deltacapsw = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/Deltasw))
               Deltacapse = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/Deltase))
            endif

         !-----------------------------------------------------------------
         ! on last subcycle, save quantities for mechanical redistribution
         !-----------------------------------------------------------------
            if (ksub == ndte) then
               if (trim(regularization) == 'RP') then
                  divu_tmp = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
                  tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
                  rdg_conv(i,j)  = -min(divu_tmp,c0)
                  rdg_shear(i,j) = p5*(tmp-abs(divu_tmp))
                  divu(i,j) = divune * tarear(i,j)

                  ! diagnostic only
                  ! shear = sqrt(tension**2 + shearing**2)
                  shear_tmp = p25*tarear(i,j)*sqrt( &
                     (tensionne + tensionnw + tensionse + tensionsw)**2 &
                     +  (shearne +   shearnw +   shearse +   shearsw)**2)
                  shear(i,j) = p25*tarear(i,j)*sqrt((tensionne)**2 + (shearne)**2)

               elseif (trim(regularization) == 'noRP') then
                  divu_tmp = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
                  tmp = p25*(((Deltane**2)/Deltacapne) + ((Deltanw**2)/Deltacapnw) + ((Deltase**2)/Deltacapse) + ((Deltasw**2)/Deltacapsw))   * tarear(i,j)
                  rdg_conv(i,j)  = -min(divu_tmp,c0)
                  rdg_shear(i,j) = p5*(tmp-abs(divu_tmp))
                  divu(i,j) = divune * tarear(i,j)

                  ! diagnostic only
                  ! shear = sqrt(tension**2 + shearing**2)
                  shear_tmp = p25*tarear(i,j)*sqrt( &
                     (tensionne + tensionnw + tensionse + tensionsw)**2 &
                     +  (shearne +   shearnw +   shearse +   shearsw)**2)
                  shear(i,j) = p25*tarear(i,j)*sqrt((tensionne)**2 + (shearne)**2)

               endif
            endif

         !-----------------------------------------------------------------
         ! strength/Delta                   ! kg/s
         !-----------------------------------------------------------------
            c0ne = strength(i,j)/Deltacapne
            c0nw = strength(i,j)/Deltacapnw
            c0sw = strength(i,j)/Deltacapsw
            c0se = strength(i,j)/Deltacapse

            c1ne = c0ne*arlx1i
            c1nw = c0nw*arlx1i
            c1sw = c0sw*arlx1i
            c1se = c0se*arlx1i

            c0ne = c1ne*ecci
            c0nw = c1nw*ecci
            c0sw = c1sw*ecci
            c0se = c1se*ecci

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
         !-----------------------------------------------------------------
            if (trim(regularization) == 'RP') then
               stressp_1(i,j) = (stressp_1(i,j)*(c1-arlx1i*revp) + c1ne*(divune*(c1+Ktens) - Deltane*(c1-Ktens))) &
                                 * denom1
               stressp_2(i,j) = (stressp_2(i,j)*(c1-arlx1i*revp) + c1nw*(divunw*(c1+Ktens) - Deltanw*(c1-Ktens))) &
                                 * denom1
               stressp_3(i,j) = (stressp_3(i,j)*(c1-arlx1i*revp) + c1sw*(divusw*(c1+Ktens) - Deltasw*(c1-Ktens))) &
                                 * denom1
               stressp_4(i,j) = (stressp_4(i,j)*(c1-arlx1i*revp) + c1se*(divuse*(c1+Ktens) - Deltase*(c1-Ktens))) &
                                 * denom1

               stressm_1(i,j) = (stressm_1(i,j)*(c1-arlx1i*revp) + c0ne*tensionne*(c1+Ktens)) * denom1
               stressm_2(i,j) = (stressm_2(i,j)*(c1-arlx1i*revp) + c0nw*tensionnw*(c1+Ktens)) * denom1
               stressm_3(i,j) = (stressm_3(i,j)*(c1-arlx1i*revp) + c0sw*tensionsw*(c1+Ktens)) * denom1
               stressm_4(i,j) = (stressm_4(i,j)*(c1-arlx1i*revp) + c0se*tensionse*(c1+Ktens)) * denom1

               stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) + c0ne*shearne*p5*(c1+Ktens)) * denom1
               stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) + c0nw*shearnw*p5*(c1+Ktens)) * denom1
               stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) + c0sw*shearsw*p5*(c1+Ktens)) * denom1
               stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) + c0se*shearse*p5*(c1+Ktens)) * denom1

            elseif (trim(regularization) == 'noRP') then
               stressp_1(i,j) = (stressp_1(i,j)*(c1-arlx1i*revp) + c1ne*(divune*(c1+Ktens) - Deltacapne*(c1-Ktens))) &
                                 * denom1
               stressp_2(i,j) = (stressp_2(i,j)*(c1-arlx1i*revp) + c1nw*(divunw*(c1+Ktens) - Deltacapnw*(c1-Ktens))) &
                                 * denom1
               stressp_3(i,j) = (stressp_3(i,j)*(c1-arlx1i*revp) + c1sw*(divusw*(c1+Ktens) - Deltacapsw*(c1-Ktens))) &
                                 * denom1
               stressp_4(i,j) = (stressp_4(i,j)*(c1-arlx1i*revp) + c1se*(divuse*(c1+Ktens) - Deltacapse*(c1-Ktens))) &
                                 * denom1

               stressm_1(i,j) = (stressm_1(i,j)*(c1-arlx1i*revp) + c0ne*tensionne*(c1+Ktens)) * denom1
               stressm_2(i,j) = (stressm_2(i,j)*(c1-arlx1i*revp) + c0nw*tensionnw*(c1+Ktens)) * denom1
               stressm_3(i,j) = (stressm_3(i,j)*(c1-arlx1i*revp) + c0sw*tensionsw*(c1+Ktens)) * denom1
               stressm_4(i,j) = (stressm_4(i,j)*(c1-arlx1i*revp) + c0se*tensionse*(c1+Ktens)) * denom1

               stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) + c0ne*shearne*p5*(c1+Ktens)) * denom1
               stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) + c0nw*shearnw*p5*(c1+Ktens)) * denom1
               stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) + c0sw*shearsw*p5*(c1+Ktens)) * denom1
               stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) + c0se*shearse*p5*(c1+Ktens)) * denom1

            endif

            if (i == 39) then
               if (j == 112) then
                  write (nu_diag,*) 'stress_I,', p5*stressp_1(i,j)/strength(i,j)    !#ab#
                  write (nu_diag,*) 'stress_II,', (p5*sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))/ strength(i,j) !#ab#
                  write (nu_diag,*) 'divu,', divune * tarear(i,j)*86400.0    !#ab#
                  write (nu_diag,*) 'shear,', tarear(i,j)*sqrt((tensionne)**2 + (shearne)**2)*86400.0    !#ab#
               endif
            endif

         !-----------------------------------------------------------------
         ! stresses for Mohr-Coulomb #ab#
         !-----------------------------------------------------------------                
         elseif (trim(yield_curve) == 'MC') then

            ! Max shear, eps_II = [(e_11 - e_22)^2 + (e_12)^2]^0.5
            mshearne = sqrt(tensionne**2 + shearne**2)
            mshearnw = sqrt(tensionnw**2 + shearnw**2)
            mshearsw = sqrt(tensionsw**2 + shearsw**2)
            mshearse = sqrt(tensionse**2 + shearse**2)

            ! Viscous capping on eps_I and eps_II
            if (trim(capping) == 'max') then
               mshearcapne = max(mshearne,smin*tinyarea(i,j))
               mshearcapnw = max(mshearnw,smin*tinyarea(i,j))
               mshearcapsw = max(mshearsw,smin*tinyarea(i,j))
               mshearcapse = max(mshearse,smin*tinyarea(i,j))

               divucapne = max(abs(divune),dmin*tinyarea(i,j))
               divucapnw = max(abs(divunw),dmin*tinyarea(i,j))
               divucapsw = max(abs(divusw),dmin*tinyarea(i,j))
               divucapse = max(abs(divuse),dmin*tinyarea(i,j))

            elseif (trim(capping) == 'add_min') then
               mshearcapne = mshearne+(smin*tinyarea(i,j))
               mshearcapnw = mshearnw+(smin*tinyarea(i,j))
               mshearcapsw = mshearsw+(smin*tinyarea(i,j))
               mshearcapse = mshearse+(smin*tinyarea(i,j))

               divucapne = abs(divune)+(dmin*tinyarea(i,j))
               divucapnw = abs(divunw)+(dmin*tinyarea(i,j))
               divucapsw = abs(divusw)+(dmin*tinyarea(i,j))
               divucapse = abs(divuse)+(dmin*tinyarea(i,j))

            elseif (trim(capping) == 'tanh') then
               mshearcapne = (smin*tinyarea(i,j))/tanh(smin*tinyarea(i,j)*(c1/mshearne))
               mshearcapnw = (smin*tinyarea(i,j))/tanh(smin*tinyarea(i,j)*(c1/mshearnw))
               mshearcapsw = (smin*tinyarea(i,j))/tanh(smin*tinyarea(i,j)*(c1/mshearsw))
               mshearcapse = (smin*tinyarea(i,j))/tanh(smin*tinyarea(i,j)*(c1/mshearse))

               divucapne = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/abs(divune)))
               divucapnw = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/abs(divunw)))
               divucapsw = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/abs(divusw)))
               divucapse = (dmin*tinyarea(i,j))/tanh(dmin*tinyarea(i,j)*(c1/abs(divuse)))

            endif

! !             ! Bulk viscosity
! !             zetane = p5*strength(i,j)/divucapne
! !             zetanw = p5*strength(i,j)/divucapnw
! !             zetasw = p5*strength(i,j)/divucapsw
! !             zetase = p5*strength(i,j)/divucapse

            if (trim(regularization) == 'RP') then
! !                ! Replacement pressure
! !                P_ne  = c2*zetane*abs(divune) ! = c0ne*abs(divune)
! !                P_nw  = c2*zetanw*abs(divunw) ! = c0nw*abs(divunw)
! !                P_sw  = c2*zetasw*abs(divusw) ! = c0sw*abs(divusw)
! !                P_se  = c2*zetase*abs(divuse) ! = c0se*abs(divuse)

               Mne = (sinphi*(abs(divune)-divune))/mshearcapne
               Mnw = (sinphi*(abs(divunw)-divunw))/mshearcapnw
               Msw = (sinphi*(abs(divusw)-divusw))/mshearcapsw
               Mse = (sinphi*(abs(divuse)-divuse))/mshearcapse
            elseif (trim(regularization) == 'noRP') then
! !                ! Pressure is ice strength
! !                P_ne  = strength(i,j) ! = c0ne*divucapne
! !                P_nw  = strength(i,j) ! = c0nw*divucapnw
! !                P_sw  = strength(i,j) ! = c0sw*divucapsw
! !                P_se  = strength(i,j) ! = c0se*divucapse

               Mne = (sinphi*(divucapne-divune))/mshearcapne
               Mnw = (sinphi*(divucapnw-divunw))/mshearcapnw
               Msw = (sinphi*(divucapsw-divusw))/mshearcapsw
               Mse = (sinphi*(divucapse-divuse))/mshearcapse
            endif

            Deltane = sqrt(divune**2 + Mne*(tensionne**2 + shearne**2))
            Deltanw = sqrt(divunw**2 + Mnw*(tensionnw**2 + shearnw**2))
            Deltasw = sqrt(divusw**2 + Msw*(tensionsw**2 + shearsw**2))
            Deltase = sqrt(divuse**2 + Mse*(tensionse**2 + shearse**2))

         !-----------------------------------------------------------------
         ! on last subcycle, save quantities for mechanical redistribution #ab#
         !-----------------------------------------------------------------
            if (ksub == ndte) then
               if (trim(regularization) == 'RP') then
                  tmpne = (Deltane**2) / abs(divune)
                  tmpnw = (Deltanw**2) / abs(divunw)
                  tmpsw = (Deltasw**2) / abs(divusw)
                  tmpse = (Deltase**2) / abs(divuse)

                  if (divune == 0.0) then
                     tmpne = 0.0
                  endif
                  if  (divunw == 0.0) then
                     tmpnw = 0.0
                  endif
                  if  (divuse == 0.0) then
                     tmpse = 0.0
                  endif
                  if  (divusw == 0.0) then
                     tmpsw = 0.0
                  endif

               elseif (trim(regularization) == 'noRP') then
                  tmpne = (Deltane**2) / divucapne
                  tmpnw = (Deltanw**2) / divucapnw
                  tmpsw = (Deltasw**2) / divucapsw
                  tmpse = (Deltase**2) / divucapse
               endif

               divu_tmp = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
               tmp   = p25*(tmpne + tmpnw + tmpse + tmpsw)   * tarear(i,j)
               rdg_conv(i,j)  = -min(divu_tmp,c0)
               rdg_shear(i,j) = p5*(tmp-abs(divu_tmp))
               divu(i,j) = divune * tarear(i,j)

               ! diagnostic only
               ! shear = sqrt(tension**2 + shearing**2)
               shear_tmp = p25*tarear(i,j)*sqrt( &
                  (tensionne + tensionnw + tensionse + tensionsw)**2 &
                  +  (shearne +   shearnw +   shearse +   shearsw)**2)
               shear(i,j) = tarear(i,j)*sqrt((tensionne)**2 + (shearne)**2)

            endif

         !-----------------------------------------------------------------
         ! constants for stress formulation                   ! kg/s #ab#
         !-----------------------------------------------------------------
            c0ne = strength(i,j)/divucapne ! = 2*zetane
            c0nw = strength(i,j)/divucapnw ! = 2*zetanw
            c0sw = strength(i,j)/divucapsw ! = 2*zetasw
            c0se = strength(i,j)/divucapse ! = 2*zetase

            c1ne = c0ne*arlx1i ! = 2*zetane/alpha
            c1nw = c0nw*arlx1i ! = 2*zetanw/alpha
            c1sw = c0sw*arlx1i ! = 2*zetasw/alpha
            c1se = c0se*arlx1i ! = 2*zetase/alpha

            c0ne = c1ne*Mne ! = (2*zetane/alpha)*Mne
            c0nw = c1nw*Mnw ! = (2*zetanw/alpha)*Mnw
            c0sw = c1sw*Msw ! = (2*zetasw/alpha)*Msw
            c0se = c1se*Mse ! = (2*zetase/alpha)*Mse

         !-----------------------------------------------------------------
         ! the stresses for Mohr-Coulomb                     ! kg/s^2 #ab#
         ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
         !-----------------------------------------------------------------
            if (trim(regularization) == 'RP') then
               stressp_1(i,j) = (stressp_1(i,j)*(c1-arlx1i*revp) + c1ne*(divune*(c1+Ktens) - abs(divune)*(c1-Ktens))) &
                           * denom1
               stressp_2(i,j) = (stressp_2(i,j)*(c1-arlx1i*revp) + c1nw*(divunw*(c1+Ktens) - abs(divunw)*(c1-Ktens))) &
                           * denom1
               stressp_3(i,j) = (stressp_3(i,j)*(c1-arlx1i*revp) + c1sw*(divusw*(c1+Ktens) - abs(divusw)*(c1-Ktens))) &
                           * denom1
               stressp_4(i,j) = (stressp_4(i,j)*(c1-arlx1i*revp) + c1se*(divuse*(c1+Ktens) - abs(divuse)*(c1-Ktens))) &
                           * denom1

               stressm_1(i,j) = (stressm_1(i,j)*(c1-arlx1i*revp) + c0ne*tensionne*(c1+Ktens)) * denom1
               stressm_2(i,j) = (stressm_2(i,j)*(c1-arlx1i*revp) + c0nw*tensionnw*(c1+Ktens)) * denom1
               stressm_3(i,j) = (stressm_3(i,j)*(c1-arlx1i*revp) + c0sw*tensionsw*(c1+Ktens)) * denom1
               stressm_4(i,j) = (stressm_4(i,j)*(c1-arlx1i*revp) + c0se*tensionse*(c1+Ktens)) * denom1

               stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) + c0ne*shearne*p5*(c1+Ktens)) * denom1
               stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) + c0nw*shearnw*p5*(c1+Ktens)) * denom1
               stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) + c0sw*shearsw*p5*(c1+Ktens)) * denom1
               stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) + c0se*shearse*p5*(c1+Ktens)) * denom1

            elseif (trim(regularization) == 'noRP') then
               stressp_1(i,j) = (stressp_1(i,j)*(c1-arlx1i*revp) + c1ne*(divune*(c1+Ktens) - divucapne*(c1-Ktens))) &
                           * denom1
               stressp_2(i,j) = (stressp_2(i,j)*(c1-arlx1i*revp) + c1nw*(divunw*(c1+Ktens) - divucapnw*(c1-Ktens))) &
                           * denom1
               stressp_3(i,j) = (stressp_3(i,j)*(c1-arlx1i*revp) + c1sw*(divusw*(c1+Ktens) - divucapsw*(c1-Ktens))) &
                           * denom1
               stressp_4(i,j) = (stressp_4(i,j)*(c1-arlx1i*revp) + c1se*(divuse*(c1+Ktens) - divucapse*(c1-Ktens))) &
                           * denom1

               stressm_1(i,j) = (stressm_1(i,j)*(c1-arlx1i*revp) + c0ne*tensionne*(c1+Ktens)) * denom1
               stressm_2(i,j) = (stressm_2(i,j)*(c1-arlx1i*revp) + c0nw*tensionnw*(c1+Ktens)) * denom1
               stressm_3(i,j) = (stressm_3(i,j)*(c1-arlx1i*revp) + c0sw*tensionsw*(c1+Ktens)) * denom1
               stressm_4(i,j) = (stressm_4(i,j)*(c1-arlx1i*revp) + c0se*tensionse*(c1+Ktens)) * denom1

               stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) + c0ne*shearne*p5*(c1+Ktens)) * denom1
               stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) + c0nw*shearnw*p5*(c1+Ktens)) * denom1
               stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) + c0sw*shearsw*p5*(c1+Ktens)) * denom1
               stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) + c0se*shearse*p5*(c1+Ktens)) * denom1
            endif

            if (i == 39) then
               if (j == 112) then
                  write (nu_diag,*) 'stress_I,', p5*stressp_1(i,j)/strength(i,j)    !#ab#
                  write (nu_diag,*) 'stress_II,', (p5*sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))/ strength(i,j) !#ab#
                  write (nu_diag,*) 'divu,', divune * tarear(i,j)*86400.0    !#ab#
                  write (nu_diag,*) 'shear,', tarear(i,j)*sqrt((tensionne)**2 + (shearne)**2)*86400.0    !#ab#
               endif
            endif

! !          !-----------------------------------------------------------------
! !          ! on last subcycle, save quantities for mechanical redistribution #ab#
! !          !-----------------------------------------------------------------
! !             if (ksub == ndte) then
! !                tmpne = (c2*zetane*Deltane**2) / P_ne
! !                tmpnw = (c2*zetanw*Deltanw**2) / P_nw
! !                tmpsw = (c2*zetasw*Deltasw**2) / P_sw
! !                tmpse = (c2*zetase*Deltase**2) / P_se
! !
! !                if (P_ne == 0.0) then
! !                   tmpne = 0.0
! !                endif
! !                if  (P_nw == 0.0) then
! !                   tmpnw = 0.0
! !                endif
! !                if  (P_se == 0.0) then
! !                   tmpse = 0.0
! !                endif
! !                if  (P_sw == 0.0) then
! !                   tmpsw = 0.0
! !                endif
! !
! !                divu_tmp = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
! !                tmp   = p25*(tmpne + tmpnw + tmpse + tmpsw)   * tarear(i,j)
! !                rdg_conv(i,j)  = -min(divu_tmp,c0)
! !                rdg_shear(i,j) = p5*(tmp-abs(divu_tmp))
! !              divu(i,j) = divune * tarear(i,j)
! !
! !                ! diagnostic only
! !                ! shear = sqrt(tension**2 + shearing**2)
! !                shear_tmp = p25*tarear(i,j)*sqrt( &
! !                     (tensionne + tensionnw + tensionse + tensionsw)**2 &
! !                    +  (shearne +   shearnw +   shearse +   shearsw)**2)
! !                shear(i,j) = p25*tarear(i,j)*sqrt((tensionne)**2 + (shearne)**2)
! !
! !
! !             endif
! !
! !          !-----------------------------------------------------------------
! !          ! constants for stress formulation                   ! kg/s #ab#
! !          !-----------------------------------------------------------------
! !             c0ne = c2*zetane*arlx1i
! !             c0nw = c2*zetanw*arlx1i
! !             c0sw = c2*zetasw*arlx1i
! !             c0se = c2*zetase*arlx1i
! !
! !             c1ne = c0ne*Mne
! !             c1nw = c0nw*Mnw
! !             c1sw = c0sw*Msw
! !             c1se = c0se*Mse
! !
! !          !-----------------------------------------------------------------
! !          ! the stresses for Mohr-Coulomb                     ! kg/s^2 #ab#
! !          ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
! !          !-----------------------------------------------------------------
! !             stressp_1(i,j) = (stressp_1(i,j)*(c1-arlx1i*revp) + c0ne*(divune*(c1+Ktens)) - (P_ne*arlx1i)*(c1-Ktens)) &
! !                           * denom1
! !             stressp_2(i,j) = (stressp_2(i,j)*(c1-arlx1i*revp) + c0nw*(divunw*(c1+Ktens)) - (P_nw*arlx1i)*(c1-Ktens)) &
! !                           * denom1
! !             stressp_3(i,j) = (stressp_3(i,j)*(c1-arlx1i*revp) + c0sw*(divusw*(c1+Ktens)) - (P_sw*arlx1i)*(c1-Ktens)) &
! !                           * denom1
! !             stressp_4(i,j) = (stressp_4(i,j)*(c1-arlx1i*revp) + c0se*(divuse*(c1+Ktens)) - (P_se*arlx1i)*(c1-Ktens)) &
! !                           * denom1
! !
! !             stressm_1(i,j) = (stressm_1(i,j)*(c1-arlx1i*revp) + c1ne*tensionne*(c1+Ktens)) * denom1
! !             stressm_2(i,j) = (stressm_2(i,j)*(c1-arlx1i*revp) + c1nw*tensionnw*(c1+Ktens)) * denom1
! !             stressm_3(i,j) = (stressm_3(i,j)*(c1-arlx1i*revp) + c1sw*tensionsw*(c1+Ktens)) * denom1
! !             stressm_4(i,j) = (stressm_4(i,j)*(c1-arlx1i*revp) + c1se*tensionse*(c1+Ktens)) * denom1
! !
! !             stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) + c1ne*shearne*p5*(c1+Ktens)) * denom1
! !             stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) + c1nw*shearnw*p5*(c1+Ktens)) * denom1
! !             stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) + c1sw*shearsw*p5*(c1+Ktens)) * denom1
! !             stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) + c1se*shearse*p5*(c1+Ktens)) * denom1
! !
! !             if (i == 39) then
! !                if (j == 112) then
! !                   write (nu_diag,*) 'stress_I,', p5*stressp_1(i,j)/strength(i,j)    !#ab#
! !                   write (nu_diag,*) 'stress_II,', (p5*sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))/ strength(i,j) !#ab#
! !                   write (nu_diag,*) 'divu,', divune * tarear(i,j)*86400.0    !#ab#
! !                   write (nu_diag,*) 'shear,', tarear(i,j)*sqrt((tensionne)**2 + (shearne)**2)*86400.0    !#ab#
! !                endif
! !             endif
! !
         endif

      !-----------------------------------------------------------------
      ! Eliminate underflows.
      ! The following code is commented out because it is relatively 
      ! expensive and most compilers include a flag that accomplishes
      ! the same thing more efficiently.  This code is cheaper than
      ! handling underflows if the compiler lacks a flag; uncomment
      ! it in that case.  The compiler flag is often described with the 
      ! phrase "flush to zero".
      !-----------------------------------------------------------------

!      call icepack_query_parameters(puny_out=puny)
!      call icepack_warnings_flush(nu_diag)
!      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
!         file=__FILE__, line=__LINE__)

!      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1(i,j) + stressp_2(i,j)
         ssigps  = stressp_3(i,j) + stressp_4(i,j)
         ssigpe  = stressp_1(i,j) + stressp_4(i,j)
         ssigpw  = stressp_2(i,j) + stressp_3(i,j)
         ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
         ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

         ssigmn  = stressm_1(i,j) + stressm_2(i,j)
         ssigms  = stressm_3(i,j) + stressm_4(i,j)
         ssigme  = stressm_1(i,j) + stressm_4(i,j)
         ssigmw  = stressm_2(i,j) + stressm_3(i,j)
         ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
         ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

         ssig12n = stress12_1(i,j) + stress12_2(i,j)
         ssig12s = stress12_3(i,j) + stress12_4(i,j)
         ssig12e = stress12_1(i,j) + stress12_4(i,j)
         ssig12w = stress12_2(i,j) + stress12_3(i,j)
         ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
         ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

         csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
         csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
         csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
         csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
         
         csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
         csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
         csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
         csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
         
         csig12ne = p222*stress12_1(i,j) + ssig122 &
                  + p055*stress12_3(i,j)
         csig12nw = p222*stress12_2(i,j) + ssig121 &
                  + p055*stress12_4(i,j)
         csig12sw = p222*stress12_3(i,j) + ssig122 &
                  + p055*stress12_1(i,j)
         csig12se = p222*stress12_4(i,j) + ssig121 &
                  + p055*stress12_2(i,j)

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij
 1005    format (a30,2x,f12.10) ! float
 1015    format (a30,2x,e15.8) ! exponential
      end subroutine stress

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
