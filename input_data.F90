module input_data

use esmf
use netcdf
use mpas_netcdf_interface, only : open_netcdf, close_netcdf, get_netcdf_dims, netcdf_err, get_netcdf_var
use utils_mod, only      : error_handler
use program_setup, only  : nvars_to_blend, variables_to_blend, interp_cartesian_winds
use model_grid, only     : mpas_mesh_type

implicit none

private

! Public variables accessible to other modules
integer, allocatable, public :: nVertLevelsPerVariable(:)
integer, public              :: nz_input = -1    ! number of input grid atm layers
integer, public              :: nzp1_input = -1  ! number of input grid atm layer interfaces
integer, public              :: nsoil_input = -1 ! number of input soil levels

! Public subroutines
public :: read_input_data, input_data_error_checks, set_nVertLevelsPerVariable, project_ull_to_uxyz, project_uxyz_to_ull

contains

subroutine set_nVertLevelsPerVariable(localpet,the_file)

   integer, intent(in)                   :: localpet
   character(len=500), intent(in)        :: the_file

   character(len=50)            :: vname
   integer                      :: error, ncid, rc
   integer                      :: id_var, i, j, k, nz
   integer                      :: xtype, ndims_var, natts, dimids(NF90_MAX_VAR_DIMS), dims(4)
   character(len=500)           :: varname
   character(len=NF90_MAX_NAME) :: dim_names(4)

   ! open the netCDF file with the data we want to read
   call open_netcdf(trim(the_file),ncid)

   if (localpet==0) write(*,*) "GETTING INPUT DATA FROM"//trim(adjustl(the_file))

   do i = 1,nvars_to_blend
       
      vname = trim(adjustl(variables_to_blend(i)))
      if (localpet == 0 ) write(*,*) 'in set_nVertLevels, ', vname

      if (interp_cartesian_winds .and. i==nvars_to_blend) then

         call get_netcdf_dims(ncid,'nVertLevels',nz_input) ! define public variable
         nz = nz_input
         nVertLevelsPerVariable(i) = nz

      else

         ! get variable ID--use direct nf90 call here because we just want the id
         error=nf90_inq_varid(ncid, trim(vname), id_var)
         call netcdf_err(error, 'reading field id' )
         
         ! get dimension information about the variable
         dims(:) = 1
         dim_names(:) = ''
         error = nf90_Inquire_Variable(ncid, id_var, varname, xtype, ndims_var, dimids, natts) ! Output is varname, xtype,ndims_var, dimids
         call netcdf_err(error, 'inquiring variable')
         do j = 1,ndims_var
            error = nf90_inquire_dimension( ncid, dimids(j), name=dim_names(j), len=dims(j) )
            dim_names(j) = trim(adjustl(dim_names(j))) ! can't do trim on an array, so do it here
         enddo
         if ( localpet == 0 ) write(*,fmt='(a28, a30,i2,4i8)')'    variable, ndims, dims = ',trim(adjustl(varname)),ndims_var,dims(1:ndims_var)
   
         if ( any(dim_names == 'nEdges') ) then
            write(*,*)'Variable '//trim(adjustl(vname))//' is defined on cell edges'
            write(*,*)'This code can only deal with variables on cell centers'
            call error_handler("Not a cell-centered variable", -23)
         endif
   
         if ( any(dim_names.eq."nVertLevels")) then
            call get_netcdf_dims(ncid,'nVertLevels',nz_input) ! define public variable
            nz = nz_input
         else if ( any(dim_names.eq."nSoilLevels")) then
            call get_netcdf_dims(ncid,'nSoilLevels',nsoil_input) ! define public variable
            nz = nsoil_input
         else if ( any(dim_names.eq."nVertLevelsP1")) then
            call get_netcdf_dims(ncid,'nVertLevelsP1',nzp1_input) ! define public variable
            nz = nzp1_input
         else ! 2D variable
            nz = 1
         endif
   
         nVertLevelsPerVariable(i) = nz ! public from this module, allocated in mpas_blending.f90

      endif

   end do

   call close_netcdf(trim(the_file),ncid)

end subroutine set_nVertLevelsPerVariable

subroutine read_input_data(localpet,input_mesh,mpas_mesh,the_file,bundle)
   integer, intent(in)                   :: localpet
   type(esmf_mesh), intent(in)           :: input_mesh
   type(mpas_mesh_type), intent(in)      :: mpas_mesh
   character(len=500), intent(in)        :: the_file
   type(esmf_fieldbundle), intent(inout) :: bundle

   character(len=50)               :: vname
   integer                         :: nCells_input, nCellsPerPET
   integer                         :: lowBound, nz, i, j
   integer                         :: error, rc, ncid
   type(esmf_field),allocatable    :: fields(:)
   real(esmf_kind_r8), allocatable :: dummy2(:,:), dummy3(:,:,:)
   real(esmf_kind_r8), allocatable :: dummy_uzonal(:,:,:), dummy_umerid(:,:,:), dummy_ux(:,:,:), dummy_uy(:,:,:), dummy_uz(:,:,:) 
   real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:) 
   real(esmf_kind_r8), pointer     :: varptrx(:,:), varptry(:,:), varptrz(:,:)

   integer :: save_zonal, save_meridional

   nCells_input = mpas_mesh%nCells
   nCellsPerPET = mpas_mesh%nCellsPerPET

   allocate(fields(nvars_to_blend))
   call ESMF_FieldBundleGet(bundle, fieldList=fields, &
                            itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN EMSF_FieldBundleGet", rc)

   !-------------------------------------
   ! Get the data from the input file
   !-------------------------------------

   ! open the netCDF file with the data we want to read
   call open_netcdf(trim(the_file),ncid)

   if ( interp_cartesian_winds ) then
      save_zonal = -1
      save_meridional = -1
   endif

   do i = 1,nvars_to_blend

      vname = trim(adjustl(variables_to_blend(i)))
      if (localpet==0) write(*,*) "- INIT FIELD ", vname

      if ( interp_cartesian_winds ) then
          if ( vname .eq. 'uReconstructZonal' ) then
              save_zonal = i
          elseif ( vname .eq. 'uReconstructMeridional' ) then 
              save_meridional = i
          endif
      endif


      lowBound = 1
      nz = nVertLevelsPerVariable(i)

      if ( interp_cartesian_winds .and. vname .eq. 'wCartesian' ) then 

             if ( save_zonal .ne. -1 .and. save_meridional .ne. -1 ) then 

                 ! this converts uzonal, umeridional to ux, uy, uz 
                 ! -  replaces the variables varptrx,varptry with ux,uy values
		 call project_ull_to_uxyz(mpas_mesh,fields(save_zonal),fields(save_meridional),fields(i))

		 if (localpet == 0 ) write(*,*) "[interp_cartesian_winds] Converted uReconstructZonal, uReconstructMeridional to u_x, u_y, u_z."

             else
                 if (localpet == 0) write(*,*) "ERROR [interp_cartesian_winds]: Could not find uReconstructZonal or uReconstructMeridional to do wind conversion."
                 stop
             endif

      else ! all other fields

          if ( nz.eq.1 ) then
    
             allocate(dummy2(nCells_input,1))
             call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
             if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
                call error_handler("IN FieldGet", rc)
             call get_netcdf_var(ncid,trim(vname),(/1/),(/nCells_input/),dummy2(:,1))
             do j = 1, nCellsPerPET
                varptr(j) = dummy2(mpas_mesh%elemIDs(j),1)
             enddo
             nullify(varptr)
             deallocate(dummy2)
          else
    
             allocate(dummy3(nz,nCells_input,1))
             call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
             if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
                call error_handler("IN FieldGet", rc)
             call get_netcdf_var(ncid,trim(vname),(/1,1/),(/nz,nCells_input/),dummy3(:,:,1))
             do j = 1, nCellsPerPET
                varptr2(j,:) = dummy3(:,mpas_mesh%elemIDs(j),1)
             enddo
             nullify(varptr2)
             deallocate(dummy3)
          endif

      endif
           
      if (localpet==0) write(*,*) "- SET ON MESH ", trim(vname)

   enddo ! end loop over number of variables

   ! update the fields in target_bundle
   call ESMF_FieldBundleAddReplace(bundle, fields, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN ESMF_FieldBundleAddReplace", rc)

   deallocate(fields)

   call close_netcdf(trim(the_file),ncid)

   if (localpet==0) write(*,*) "DONE GETTING DATA FROM ",trim(adjustl(the_file))
   if (localpet==0) write(*,*)''
   
end subroutine read_input_data

subroutine input_data_error_checks(file1, file2)
   character(len=*), intent(in) :: file1, file2

   integer :: ncid1, ncid2, nCells1, nCells2

   ! Open file1 and get info
   call open_netcdf(trim(file1),ncid1)
   call get_netcdf_dims(ncid1,'nCells',nCells1)

   ! Open file2 and get info
   call open_netcdf(trim(file2),ncid2)
   call get_netcdf_dims(ncid2,'nCells',nCells2)

   ! Close files--no longer needed
   call close_netcdf(trim(file1),ncid1)
   call close_netcdf(trim(file2),ncid2)

   if ( nCells1 .ne. nCells2 ) then
      call error_handler("Mismatched horizontal dimensions in input data", -111)
   endif

end subroutine input_data_error_checks



subroutine project_ull_to_uxyz(mpas_mesh,field_zonal,field_meridional,field_wCartesian)
   ! This code projects zonal/meridional winds to Cartesian coordinates (x,y,z)
   ! field_zonal, field_meridional are input zonal/meridional winds, but are output
   ! as ux, uy; field_wCartesian is output as uz
   type(mpas_mesh_type), intent(in)  :: mpas_mesh
   type(esmf_field), intent(in) :: field_zonal, field_meridional, field_wCartesian

   integer :: i,j,k
   real(esmf_kind_r8) :: sinlat, sinlon, coslat, coslon, u, v
   real(esmf_kind_r8), pointer :: varptrx(:,:), varptry(:,:), varptrz(:,:) 
   integer :: rc, nCellsPerPET

   nCellsPerPET = mpas_mesh%nCellsPerPET

   ! read in uzonal and umeridional
   call ESMF_FieldGet(field_zonal, farrayPtr=varptrx, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldGet", rc)

   call ESMF_FieldGet(field_meridional, farrayPtr=varptry, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldGet", rc)

   ! index i = wCartesian
   call ESMF_FieldGet(field_wCartesian, farrayPtr=varptrz, rc=rc) 
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldGet", rc)

   do j = 1, nCellsPerPET
         i = mpas_mesh%elemIDs(j)
	 sinlat = sin(mpas_mesh%latCell(i))
	 sinlon = sin(mpas_mesh%lonCell(i))
	 coslat = cos(mpas_mesh%latCell(i))
	 coslon = cos(mpas_mesh%lonCell(i))
         do k = 1, nz_input
           u = varptrx(j,k)
           v = varptry(j,k)
           varptrx(j,k) = - u*sinlon - v*sinlat*coslon
           varptry(j,k) = u*coslon - v*sinlat*sinlon
           varptrz(j,k) = v*coslat 
         enddo
   enddo

   nullify(varptrx)
   nullify(varptry)
   nullify(varptrz)

end subroutine project_ull_to_uxyz 


subroutine project_uxyz_to_ull(mpas_mesh,field_ux,field_uy,field_uz,field_zonal,field_meridional)

   type(mpas_mesh_type), intent(in)  :: mpas_mesh
   type(esmf_field), intent(in) :: field_ux, field_uy, field_uz
   type(esmf_field), intent(out) :: field_zonal, field_meridional

   integer :: i,j,k
   real(esmf_kind_r8) :: sinlat, sinlon, coslat, coslon, u, v
   real(esmf_kind_r8), pointer :: varptrx(:,:), varptry(:,:), varptrz(:,:), ptru(:,:), ptrv(:,:)
   integer :: rc, nCellsPerPET

   type(esmf_mesh) :: mesh
   integer :: nz


   ! This code projects Cartesian-coordinate winds back to zonal/meridional winds 
   nCellsPerPET = mpas_mesh%nCellsPerPET

   ! read in ux, uy, uz winds
   call ESMF_FieldGet(field_ux, farrayPtr=varptrx, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldGet", rc)

   call ESMF_FieldGet(field_uy, farrayPtr=varptry, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldGet", rc)

   call ESMF_FieldGet(field_uz, farrayPtr=varptrz, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldGet", rc)

   ! get the mesh from field_ux
   call ESMF_FieldGet(field_ux, mesh=mesh, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldGet", rc)
   nz = size(varptrx,2)

   ! create temporary ESMF fields to hold the converted zonal/meridional winds for output
   field_zonal = ESMF_FieldCreate(mesh, &
                    typekind=ESMF_TYPEKIND_R8, &
                    meshloc=ESMF_MESHLOC_ELEMENT, &
                    name='uReconstructZonal', &
                    ungriddedLBound=(/1/), &
                    ungriddedUBound=(/nz/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldCreate", rc)


   field_meridional = ESMF_FieldCreate(mesh, &
                    typekind=ESMF_TYPEKIND_R8, &
                    meshloc=ESMF_MESHLOC_ELEMENT, &
                    name='uReconstructMeridional', &
                    ungriddedLBound=(/1/), &
                    ungriddedUBound=(/nz/), rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
       call error_handler("IN FieldCreate", rc)

   ! Now get the ESMF fields to store the converted winds
   call ESMF_FieldGet(field_zonal, farrayPtr=ptru, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN FieldGet", rc)

   call ESMF_FieldGet(field_meridional, farrayPtr=ptrv, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__line__,file=__file__)) &
      call error_handler("IN FieldGet", rc)

   do j = 1, nCellsPerPET
         i = mpas_mesh%elemIDs(j)
	 sinlat = sin(mpas_mesh%latCell(i))
	 sinlon = sin(mpas_mesh%lonCell(i))
	 coslat = cos(mpas_mesh%latCell(i))
	 coslon = cos(mpas_mesh%lonCell(i))
         do k = 1, nz_input
           ptru(j,k) = -varptrx(j,k)*sinlon + varptry(j,k)*coslon
           ptrv(j,k) = -varptrx(j,k)*sinlat*coslon - varptry(j,k)*sinlon*sinlat + varptrz(j,k)*coslat
         enddo
   enddo

   nullify(ptru)
   nullify(ptrv)
   nullify(varptrx)
   nullify(varptry)
   nullify(varptrz)

end subroutine project_uxyz_to_ull


 
end module input_data
