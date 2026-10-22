module write_data

use esmf
use netcdf
use mpi
use utils_mod,     only : error_handler
use program_setup, only : variables_to_blend, &
                          nlat, nlon, lat_ll, lon_ll, &
                          nvars_to_blend, dx_in_degrees, dy_in_degrees, &
                          large_scale_file, output_blended_edge_normal_wind !, valid_time
use model_grid, only  : mpas_mesh_type, &
                        lons_output_grid, lats_output_grid, &
                        grid_files_heirarchy
use input_data, only  : nVertLevelsPerVariable, &
                        nz_input, nzp1_input, nsoil_input, project_ull_to_uxyz, project_uxyz_to_ull
use mpas_netcdf_interface, only : netcdf_err, open_netcdf, close_netcdf, get_netcdf_var, get_netcdf_dims

implicit none

private

! Public subroutines visible to other modules
public :: write_to_file, write_to_file_latlon

contains

subroutine write_to_file(localpet,mpas_mesh,input_bundle,get_metadata_from_template,file_template,output_file)

   integer, intent(in)                 :: localpet
   type(mpas_mesh_type),        intent(in)  :: mpas_mesh
   type(esmf_fieldbundle), intent(in)  :: input_bundle
   logical, intent(in)                 :: get_metadata_from_template
   character(len=*), intent(in)        :: file_template ! netcdf template with dimensions/vars/attrs we want to carry to output file
   character(len=*), intent(in)        :: output_file ! name of file you want to write

   character(len=500)                :: varname
   character(len=20)                 :: tempstr(1,19)

   logical                          :: found_it
   integer, parameter               :: Datestrlen=19
   integer                          :: error, ncidin, ncidout, rc, i, nz, f
   integer                          :: header_buffer_val = 16384
   integer                          :: dim_time, dim_z, dim_zp1, dim_soil, my_dim_z
   integer                          :: dim_ncells, nCells, nEdges, ncid, dim_nedges
   integer                          :: nCells_test, nEdges_test
   integer                          :: id_lat, id_lon, id_times, id_var, id_var2, id_mask
   integer                          :: id_dim, ndims, nvars, ngatts,unlimdimid
   integer                          :: idims,dimsval,ivars,idims2,igatts
   integer                          :: varstype,varsndims,varsdimids(4),varsnatts, ivarsnatts
   character(len=100)               :: DIMSNAME,VARSNAME,ATTSNAME
   real(esmf_kind_r8), allocatable  :: dum1d(:), dum1dt(:,:), dum2d(:,:), dum2dt(:,:,:)
   type(esmf_field), allocatable    :: fields(:)
   integer, allocatable             :: edgesOnCell(:,:),cellsOnEdge(:,:)
   real, allocatable                :: zonal(:,:), meridional(:,:), edge_normal_wind(:,:)
   real, allocatable                :: edgeNormalVectors(:,:)
   character(len=500)               :: fnames(3)

   nCells = mpas_mesh%nCells
   nEdges = mpas_mesh%nEdges

   if (localpet == 0) then

      !--- open the file that we are outputting; nf90_create opens in define mode
      error = nf90_create(trim(adjustl(output_file)), NF90_NETCDF4, ncidout) ! make the type match the input file???
      call netcdf_err(error, 'CREATING FILE '//trim(adjustl(output_file)) )
      write(*,*)'Writing '//trim(adjustl(output_file))

     !error = nf90_redef(ncidout) ! enter define mode

      !--- open the template file and copy over stuff
      if ( get_metadata_from_template ) then
         error=nf90_open(trim(file_template),nf90_nowrite,ncidin)
         call netcdf_err(error, 'opening template file')

         ! copy dimensions from template file
         error=nf90_inquire(ncidin,ndims,nvars,ngatts,unlimdimid)
         do idims=1,ndims
            error=nf90_inquire_dimension(ncidin,idims,DIMSNAME,dimsval)
            error=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),dimsval,idims2)
            write(*,*)'copied ',trim(adjustl(DIMSNAME))
         enddo

         ! specifically get the dimension IDs for nCells and Time 
         !  incase we need them later
         ! these really should be there, so no error checking....
         error = nf90_inq_dimid(ncidin,'nCells', dim_ncells)
         error = nf90_inq_dimid(ncidin,'nEdges', dim_nedges)
         error = nf90_inq_dimid(ncidin,'Time', dim_time)
         error = nf90_inq_dimid(ncidin,'nVertLevels', dim_z)

         ! copy variables and their attributes from template file
         do i = 1,nvars_to_blend
            varname = variables_to_blend(i)
            error=nf90_inq_varid(ncidin, trim(adjustl(varname)), id_var)
            error=nf90_inquire_variable(ncidin,id_var,VARSNAME,varstype,varsndims,varsdimids,varsnatts)
            error = nf90_def_var(ncidout,trim(adjustl(varname)), varstype, varsdimids(1:varsndims), id_var2)
            do ivarsnatts=1,varsnatts
               error=nf90_inq_attname(ncidin,id_var,ivarsnatts,ATTSNAME)
               error=nf90_copy_att(ncidin,id_var,ATTSNAME,ncidout,id_var2)
            enddo
         enddo

         ! copy global attributes from template file input_file_template
         do igatts=1,ngatts
            error=nf90_inq_attname(ncidin,nf90_global,igatts,ATTSNAME)
            error=nf90_copy_att(ncidin,nf90_global,ATTSNAME,ncidout,nf90_global)
         enddo

        ! we no longer need this the template file
        error = nf90_close(ncidin)
        call netcdf_err(error, 'CLOSING TEMPLATE FILE')

      else ! not get_metadata_from_template

         !--- define nCells and Time
         error = nf90_def_dim(ncidout, 'nCells', nCells, dim_ncells)
         call netcdf_err(error, 'Defining nCells')
         error = nf90_def_dim(ncidout, 'Time', NF90_UNLIMITED, dim_time)
         call netcdf_err(error, 'Defining Time')

         !--- define vertical dimensions if needed
         ! These variables are defaulted to -1, and will be -1 unless at
         !   least one of the variables to blend is on the given 3D levels
         if ( nz_input > 0 ) then 
            error = nf90_def_dim(ncidout, 'nVertLevels', nz_input, dim_z)
            call netcdf_err(error, 'Defining nVertLevels')
         endif
         if ( nzp1_input > 0 ) then 
            error = nf90_def_dim(ncidout, 'nVertLevelsP1', nzp1_input, dim_zp1)
            call netcdf_err(error, 'Defining nVertLevelsP1')
         endif
         if ( nsoil_input > 0 ) then 
            error = nf90_def_dim(ncidout, 'nSoilLevels', nsoil_input, dim_soil)
            call netcdf_err(error, 'Defining nSoilLevels')
         endif

         do i = 1,nvars_to_blend
            varname = variables_to_blend(i)
            if ( nVertLevelsPerVariable(i) == 1 ) then
               error = nf90_def_var(ncidout, trim(adjustl(varname)), NF90_FLOAT, (/dim_ncells, dim_time/), id_var2)
            else
               nz = nVertLevelsPerVariable(i)
               if ( nz == nz_input ) then
                  my_dim_z = dim_z
               else if ( nz == nsoil_input ) then
                  my_dim_z = dim_soil
               else if ( nz == nzp1_input ) then
                  my_dim_z = dim_zp1
               endif
               error = nf90_def_var(ncidout, trim(adjustl(varname)), NF90_FLOAT, (/my_dim_z, dim_ncells, dim_time/), id_var2)
            endif
           !do ivarsnatts=1,varsnatts
           !   error=nf90_inq_attname(ncidin,id_var,ivarsnatts,ATTSNAME)
           !   error=nf90_copy_att(ncidin,id_var,ATTSNAME,ncidout,id_var2)
           !enddo
         enddo 

      endif !get_metadata_from_template

      ! Define latitude/longitude. Make sure it's not already there first
      !  If it's there, error == nf90_noerr
      ! Either way, if it's already there, or we define it, we also
      !  define id_lat, id_lon
      error = nf90_inq_varid(ncidout, 'latCell', id_lat)
      if (error /= nf90_noerr) then
        !error = nf90_def_var(ncidout, 'latCell', NF90_FLOAT, (/dim_ncells, dim_time/), id_lat)
         error = nf90_def_var(ncidout, 'latCell', NF90_FLOAT, (/dim_ncells/), id_lat)
         call netcdf_err(error, 'Defining latCell')
      endif

      error = nf90_inq_varid(ncidout, 'lonCell', id_lon)
      if (error /= nf90_noerr) then
        !error = nf90_def_var(ncidout, 'lonCell', NF90_FLOAT, (/dim_ncells, dim_time/), id_lon)
         error = nf90_def_var(ncidout, 'lonCell', NF90_FLOAT, (/dim_ncells /), id_lon)
         call netcdf_err(error, 'Defining lonCell')
      endif

      ! All done with dimensions and variables so exit define mode
      error = nf90_enddef(ncidout)

      allocate(dum1d(nCells))
      allocate(dum1dt(nCells,1))
   else
      allocate(dum1d(0))
      allocate(dum1dt(0,0))
   endif ! localpet = 0

   allocate(fields(nvars_to_blend))
   call ESMF_FieldBundleGet(input_bundle, fieldList=fields, &
                    itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=error) 
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleGet", error)

   do i = 1,nvars_to_blend
      varname = variables_to_blend(i)
      !call ESMF_FieldGet(fields(i),name=varname,rc=error)
      !call ESMF_FieldGet(fields(i), dimCount=ndims, rc=rc)

      !if (ndims==2) then
      if ( nVertLevelsPerVariable(i) == 1 ) then

        !if (localpet==0) print*,"- CALL FieldGather for 2d field ", trim(varname)
         call ESMF_FieldGather(fields(i), dum1d, rootPet=0, rc=error)
         if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", error)

         if (localpet==0) then
           !print*, "- WRITE TO FILE ", trim(varname)
            dum1dt(:,1) = dum1d
            error = nf90_inq_varid(ncidout, trim(adjustl(varname)), id_var2)
            call netcdf_err(error, 'Getting ID for '//trim(adjustl(varname)) )
            error = nf90_put_var( ncidout, id_var2, dum1dt, count=(/nCells,1/) )
            call netcdf_err(error, 'WRITING RECORD')
         endif

      else ! 3d variable

         nz = nVertLevelsPerVariable(i)

         if (localpet==0) then
            allocate(dum2d(nCells,nz)) ! this is the dimension order things are assigned to the ESMF field
            allocate(dum2dt(nz,nCells,1))
         else
            allocate(dum2d(0,0))
            allocate(dum2dt(0,0,0))
         endif

        !if (localpet==0) print*,"- CALL FieldGather for 3d field ", trim(varname)
         call ESMF_FieldGather(fields(i), dum2d, rootPet=0, rc=error)
         if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldGather", error)

         if (localpet==0) then
           !print*, "- WRITE TO FILE ", trim(varname)
            dum2dt(:,:,1) = transpose(dum2d)
            error = nf90_inq_varid( ncidout, trim(adjustl(varname)), id_var2)
            call netcdf_err(error, 'Getting ID for '//trim(adjustl(varname)) )
            error = nf90_put_var( ncidout, id_var2, dum2dt, count=(/nz,nCells,1/) )
            call netcdf_err(error, 'WRITING RECORD' )
         endif
         
         deallocate(dum2d, dum2dt)

      endif ! end if 3d vars
   enddo ! loop over fields

   !--- write latitude, longitude

   ! latitude
   if (localpet ==0) then
      dum1dt(:,1) = mpas_mesh%latCell
      error = nf90_inq_varid(ncidout, 'latCell', id_lat)
     !error = nf90_put_var( ncidout, id_lat, dum1dt, start=(/1,1/), count=(/nCells,1/) )
      error = nf90_put_var( ncidout, id_lat, dum1dt(:,1), start=(/1/), count=(/nCells/) )
      call netcdf_err(error, 'WRITING LATITUDE RECORD' )
   endif

   ! longitude
   if (localpet ==0) then
      dum1dt(:,1) = mpas_mesh%lonCell
      error = nf90_inq_varid(ncidout, 'lonCell', id_lon)
     !error = nf90_put_var( ncidout, id_lon, dum1dt, start=(/1,1/), count=(/nCells,1/) )
      error = nf90_put_var( ncidout, id_lon, dum1dt(:,1), start=(/1/), count=(/nCells/) )
      call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
   endif

   ! bdyCell-need to create the variable
   if (localpet ==0) then
      dum1dt(:,1) = mpas_mesh%bdyMaskCellOrig
      error = nf90_def_var(ncidout, 'bdyMaskCell', NF90_INT, (/dim_ncells /), id_mask)
      error = nf90_put_var( ncidout, id_mask, dum1dt(:,1), start=(/1/), count=(/nCells/) )
      call netcdf_err(error, 'WRITING bdyMaskCell RECORD' )

     !dum1dt(:,1) = mpas_mesh%bdyMaskCell
     !error = nf90_def_var(ncidout, 'bdyMaskCell', NF90_INT, (/dim_ncells /), id_mask)
     !error = nf90_put_var( ncidout, id_mask, dum1dt(:,1), start=(/1/), count=(/nCells/) )
     !call netcdf_err(error, 'WRITING bdyMaskCell RECORD' )
   endif

   !times
!  if (localpet==0)  print*,"- WRITE TO FILE TARGET GRID Times"
!  if (localpet ==0) then
!     tempstr(1,:) = valid_time
!     error = nf90_put_var( ncidout, id_times, tempstr, start = (/1,1/), count=(/Datestrlen,1/))
!     call netcdf_err(error, 'WRITING TIMES RECORD' )
!  endif

   if (localpet==0) then
      !error = nf90_enddef(ncidout, header_buffer_val,4,0,4)
      !call netcdf_err(error, 'DEFINING STUFF' )
      error = nf90_close(ncidout)
      call netcdf_err(error, 'CLOSING OUTPUT FILE' )
   endif

   deallocate(fields, dum1d, dum1dt)

   ! We're not quite done yet.  If we need to output MPAS's "u" variable,
   !   velocity normal to the cell edges, then we need to derive it.
   ! The subroutine call is simple, but there are a variety of ways
   !   to go about getting the needed input data.
   ! We will use the file we just wrote to get the input data, and
   !   then overwrite "u" in the blended field
   ! We really only want to call this at the very end when outputting final blended file
   !   ...the conditional should maybe change.
   if ( output_blended_edge_normal_wind .and. get_metadata_from_template .and. localpet == 0 ) then

      ! Make sure we blended uReconstructZonal and uReconstructMeridional
      !  If not, we can't get a meaningful blended u, so exit
      if ( .not. any(variables_to_blend .eq. 'uReconstructZonal') .or. &
           .not. any(variables_to_blend .eq. 'uReconstructMeridional')) then
         call error_handler("Missing variables for u", -76)
      endif

      allocate(zonal(nz_input,nCells))
      allocate(meridional(nz_input,nCells))
      allocate(edge_normal_wind(nz_input,mpas_mesh%nEdges))
      allocate(edgeNormalVectors(3,mpas_mesh%nEdges))

      ! We need a few more variables about the geometry for the mesh we are outputting on
      allocate(edgesOnCell(mpas_mesh%maxEdges,nCells))
      allocate(cellsOnEdge(2,mpas_mesh%nEdges))
      call open_netcdf(trim(grid_files_heirarchy(1)),ncid)
      call get_netcdf_var(ncid,'edgesOnCell',(/1,1/),(/mpas_mesh%maxEdges,nCells/),edgesOnCell)
      call get_netcdf_var(ncid,'cellsOnEdge',(/1,1/),(/2,mpas_mesh%nEdges/),cellsOnEdge)
      call close_netcdf(trim(grid_files_heirarchy(1)),ncid)

      ! open the file that we just wrote (probably has blended data)
      !  for reading only. get data, then close it.
      call open_netcdf(trim(output_file),ncid)
      call get_netcdf_var(ncid,'uReconstructZonal',(/1,1/),(/nz_input,nCells/),zonal)
      call get_netcdf_var(ncid,'uReconstructMeridional',(/1,1/),(/nz_input,nCells/),meridional)
      call close_netcdf(trim(output_file),ncid)

      ! project the blended cell-centered winds to the edges
      ! We need edgeNormalVectors to do this, which is not in a "grid.nc" file, but only
      !  in static.nc and init.nc file types.  We can try to get it from a file. Otherwise,
      !  we need to derive it

      ! try file_template, large_scale_file, and grid_files_heirarchy(1)
      found_it = .false.
      fnames = (/ file_template, large_scale_file, grid_files_heirarchy(1) /)
      do f = 1,3
         call open_netcdf(trim(fnames(f)),ncid)
         error = nf90_inq_varid(ncid, 'edgeNormalVectors', id_var)
         if ( error == 0 ) then
            call get_netcdf_var(ncid,'edgeNormalVectors',(/1,1/),(/3,mpas_mesh%nEdges/),edgeNormalVectors)
            ! make sure the values aren't 0 everywhere
            if ( all(edgeNormalVectors .eq. 0)) then
               write(*,*)'found edgeNormalVectors from '//trim(fnames(f))//' but it is everywhere 0 '
            else
               found_it = .true. ! we found valid data!
               write(*,*)'got edgeNormalVectors from '//trim(fnames(f))
            endif
         endif
         call close_netcdf(trim(fnames(f)),ncid)
         if ( found_it ) exit ! get out of loop
      enddo
      if ( .not. found_it ) then ! we have to derive it :(
         write(*,*)'deriving edgeNormalVectors'
         call derive_edgeNormalVectors(grid_files_heirarchy(1), mpas_mesh, cellsOnEdge, edgeNormalVectors)
      endif

      call project_to_edges(mpas_mesh,zonal,meridional,edgesOnCell,edgeNormalVectors,edge_normal_wind)

      ! now open up the blended file for writing, write to it, and close
      error = nf90_open(path=trim(adjustl(output_file)), mode=nf90_write,ncid=ncid)
      error = nf90_def_var(ncid, 'u', NF90_FLOAT, (/dim_z, dim_nedges, dim_time/), id_var)
      error = nf90_put_var(ncid, id_var, edge_normal_wind)
      call netcdf_err(error, 'WRITING U')
      error = nf90_put_att(ncid, id_var, "units", "m s^{-1}")
      error = nf90_put_att(ncid, id_var, "long_name", "Horizontal normal velocity at edges")
      error = nf90_close(ncid)
      call netcdf_err(error, 'CLOSING OUTPUT FILE' )

      deallocate(zonal,meridional,edge_normal_wind)
      deallocate(edgeNormalVectors)
      deallocate(edgesOnCell,cellsOnEdge)
   endif

end subroutine write_to_file

subroutine write_to_file_latlon(localpet,nmeshes,input_bundle,output_file)

   integer, intent(in)                 :: localpet, nmeshes
   type(esmf_fieldbundle), intent(in)  :: input_bundle(nmeshes)
   character(len=*), intent(in)        :: output_file ! name of file you want to write

   character(len=500)                :: varname
   character(len=20)                 :: tempstr(1,19)

   integer, parameter               :: Datestrlen=19
   integer                          :: error, ncidout, rc, i, k, nz
   integer                          :: header_buffer_val = 16384
   integer                          :: dim_time, dim_lon, dim_lat, dim_z, dim_zp1, dim_soil, my_dim_z
   integer                          :: dim_str, ndims
   integer                          :: id_lat, id_lon, id_times, id_var, id_var2
   integer                          :: i_target, j_target
   real(esmf_kind_r8), allocatable  :: dum2d(:,:), dum2dt(:,:,:), &
                                       dum3d(:,:,:), dum3dt(:,:,:,:)
   type(esmf_field), allocatable    :: fields(:)

   i_target = nlon ! from program_setup
   j_target = nlat
   if (localpet ==0 ) then
      allocate(dum2d(i_target,j_target))
      allocate(dum2dt(i_target,j_target,1))
   else
      allocate(dum2d(0,0))
      allocate(dum2dt(0,0,0))
   endif

   if (localpet == 0) then

      !--- open the file
      error = nf90_create(output_file, NF90_NETCDF4, ncidout) ! make the type match the input file???
      call netcdf_err(error, 'CREATING FILE '//trim(output_file) )
      write(*,*)'Writing '//trim(adjustl(output_file))

      !--- define dimensions
      error = nf90_def_dim(ncidout, 'nmeshes', NF90_UNLIMITED , dim_time) ! dim_time is a misnomer...
      call netcdf_err(error, 'DEFINING Time DIMENSION' )
      error = nf90_def_dim(ncidout, 'longitude', i_target, dim_lon)
      call netcdf_err(error, 'DEFINING STAGGERED LON DIMENSION' )
      error = nf90_def_dim(ncidout, 'latitude', j_target, dim_lat)
      call netcdf_err(error, 'DEFINING LAT DIMENSION' )
      if ( nz_input > 0 ) then 
         error = nf90_def_dim(ncidout, 'nVertLevels', nz_input, dim_z)
         call netcdf_err(error, 'DEFINING VERTICAL DIMENSION' )
      endif
      if ( nzp1_input > 0 ) then 
         error = nf90_def_dim(ncidout, 'nVertLevelsP1', nzp1_input, dim_zp1)
         call netcdf_err(error, 'DEFINING VERTICALP1 DIMENSION' )
      endif
      if ( nsoil_input > 0 ) then 
         error = nf90_def_dim(ncidout, 'nSoilLevels', nsoil_input, dim_soil)
         call netcdf_err(error, 'DEFINING SOIL DIMENSION' )
      endif

      error = nf90_put_att(ncidout, NF90_GLOBAL, 'DX', dx_in_degrees)
      call netcdf_err(error, 'DEFINING DX GLOBAL ATTRIBUTE')
      error = nf90_put_att(ncidout, NF90_GLOBAL, 'DY', dy_in_degrees)
      call netcdf_err(error, 'DEFINING DY GLOBAL ATTRIBUTE')

      !--- define variables
      do i = 1,nvars_to_blend
         varname = variables_to_blend(i)
         if ( nVertLevelsPerVariable(i) == 1 ) then
            error = nf90_def_var(ncidout, trim(adjustl(varname)), NF90_FLOAT, (/dim_lon, dim_lat, dim_time/), id_var)
         else
            nz = nVertLevelsPerVariable(i)
            if ( nz == nz_input ) then
               my_dim_z = dim_z
            else if ( nz == nsoil_input ) then
               my_dim_z = dim_soil
            else if ( nz == nzp1_input ) then
               my_dim_z = dim_zp1
            endif
            error = nf90_def_var(ncidout, trim(adjustl(varname)), NF90_FLOAT, (/dim_lon, dim_lat, my_dim_z, dim_time/), id_var)
         endif
         !do ivarsnatts=1,varsnatts
         !   error=nf90_inq_attname(ncidin,id_var,ivarsnatts,ATTSNAME)
         !   error=nf90_copy_att(ncidin,id_var,ATTSNAME,ncidout,id_var2)
         !enddo
         !error = nf90_put_att(ncidout, id_var, "units", target_units(i))
         !call netcdf_err(error, 'DEFINING UNITS' )
         !error = nf90_put_att(ncidout, id_var, "description", target_longname(i))
         !call netcdf_err(error, 'DEFINING LONG_NAME' )
      enddo 

      !--- define lats/lons
      error = nf90_def_var(ncidout, 'XLONG', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lon)
      call netcdf_err(error, 'DEFINING GEOLON FIELD' )
      error = nf90_put_att(ncidout, id_lon, "description", "LONGITUDE, WEST IS NEGATIVE")
      call netcdf_err(error, 'DEFINING GEOLON NAME' )
      error = nf90_put_att(ncidout, id_lon, "units", "degree_east")
      call netcdf_err(error, 'DEFINING GEOLON UNITS' )

      error = nf90_def_var(ncidout, 'XLAT', NF90_FLOAT, (/dim_lon,dim_lat, dim_time/), id_lat)
      call netcdf_err(error, 'DEFINING GEOLAT FIELD' )
      error = nf90_put_att(ncidout, id_lat, "description", "LATITUDE, SOUTH IS NEGATIVE")
      call netcdf_err(error, 'DEFINING GEOLAT NAME' )
      error = nf90_put_att(ncidout, id_lat, "units", "degree_north")
      call netcdf_err(error, 'DEFINING GEOLAT UNITS' )

     !error = nf90_def_var(ncidout, 'Times', NF90_CHAR, (/dim_str, dim_time/), id_times)
     !call netcdf_err(error, 'DEFINING Times FIELD' )
     !error = nf90_put_att(ncidout, id_times, "description", "Times")
     !call netcdf_err(error, 'DEFINING Times NAME' )

      ! All done with dimensions and variables so exit define mode
      error = nf90_enddef(ncidout)

   endif ! localpet = 0

   allocate(fields(nvars_to_blend))

   do k = 1,nmeshes
      call ESMF_FieldBundleGet(input_bundle(k), fieldList=fields, &
                       itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                       rc=error)  ! can probably move outside the loop
      if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleGet", error)

      do i = 1,nvars_to_blend
         varname = variables_to_blend(i)
         !call ESMF_FieldGet(fields(i),name=varname,rc=error)
         !call ESMF_FieldGet(fields(i), dimCount=ndims, rc=rc)

         !if (ndims==2) then
         if ( nVertLevelsPerVariable(i) == 1 ) then

           !if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname)
            call ESMF_FieldGather(fields(i), dum2d, rootPet=0, rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
               call error_handler("IN FieldGather", error)

            if (localpet==0) then
              !print*, "- WRITE TO FILE ", trim(varname)
               dum2dt(:,:,1) = dum2d
               error = nf90_inq_varid( ncidout, trim(adjustl(varname)), id_var2)
               error = nf90_put_var( ncidout, id_var2, dum2dt, start=(/1,1,k/), count=(/i_target,j_target,1/))
               call netcdf_err(error, 'WRITING RECORD')
            endif

         else ! 3d variable

            nz = nVertLevelsPerVariable(i)

            if ( nz == nz_input ) then
               my_dim_z = dim_z
            else if ( nz == nsoil_input ) then
               my_dim_z = dim_soil
            else if ( nz == nzp1_input ) then
               my_dim_z = dim_zp1
            endif

            if (localpet==0) then
               allocate(dum3d(i_target,j_target,nz))
               allocate(dum3dt(i_target,j_target,nz,1))
            else
               allocate(dum3d(0,0,0))
               allocate(dum3dt(0,0,0,0))
            endif

           !if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID ", trim(varname)

            call ESMF_FieldGather(fields(i), dum3d, rootPet=0, rc=error)
            if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
               call error_handler("IN FieldGather", error)

            if (localpet==0) then
               dum3dt(:,:,:,1) = dum3d
               error = nf90_inq_varid( ncidout, trim(adjustl(varname)), id_var2)
               error = nf90_put_var( ncidout, id_var2, dum3dt, start=(/1,1,1,k/), count=(/i_target,j_target,nz,1/))
               call netcdf_err(error, 'WRITING RECORD' )
           endif

           deallocate(dum3d, dum3dt)

         endif ! end if 3d vars
      enddo ! loop over fields (i)
   enddo ! loop over grids (k)

   deallocate(fields)

!--- write latitude, longitude

   ! longitude
!  if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID LONGITUDE"
!  call ESMF_FieldGather(longitude_target_grid, dum2d, rootPet=0, rc=error)
!  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!     call error_handler("IN FieldGather", error)
   if (localpet ==0) then
      do k = 1,nmeshes
         dum2dt(:,:,1) = lons_output_grid
         error = nf90_put_var( ncidout, id_lon, dum2dt, start=(/1,1,k/), count=(/i_target,j_target,1/))
         call netcdf_err(error, 'WRITING LONGITUDE RECORD' )
      enddo
   endif

   ! latitude
!  if (localpet==0) print*,"- CALL FieldGather FOR TARGET GRID LATITUDE"
!  call ESMF_FieldGather(latitude_target_grid, dum2d, rootPet=0, rc=error)
!  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
!     call error_handler("IN FieldGather", error)
   if (localpet ==0) then
      do k = 1,nmeshes
         dum2dt(:,:,1) = lats_output_grid
         error = nf90_put_var( ncidout, id_lat, dum2dt, start=(/1,1,k/), count=(/i_target,j_target,1/))
         call netcdf_err(error, 'WRITING LATITUDE RECORD' )
      enddo
   endif

   !times
!  if (localpet==0)  print*,"- WRITE TO FILE TARGET GRID Times"
!     if (localpet ==0) then
!     tempstr(1,:) = valid_time
!     error = nf90_put_var( ncidout, id_times, tempstr, start = (/1,1/), count=(/Datestrlen,1/))
!     call netcdf_err(error, 'WRITING TIMES RECORD' )
!  endif

   if (localpet==0) then
     !error = nf90_enddef(ncidout, header_buffer_val,4,0,4)
     !call netcdf_err(error, 'DEFINING STUFF' )
      error = nf90_close(ncidout)
      call netcdf_err(error, 'CLOSING FILE' )
   endif

   deallocate(dum2d, dum2dt)

end subroutine write_to_file_latlon

subroutine project_to_edges(mpas_mesh,zonal,meridional,edgesOnCell,edgeNormalVectors,edge_normal_wind)

   type(mpas_mesh_type), intent(in)  :: mpas_mesh
   real,  dimension(nz_input,mpas_mesh%nCells), intent(in) :: zonal, meridional
   integer,  dimension(mpas_mesh%maxEdges,mpas_mesh%nCells), intent(in) :: edgesOnCell
   real,  dimension(3,mpas_mesh%nEdges),        intent(in) :: edgeNormalVectors
   real,  dimension(nz_input,mpas_mesh%nEdges), intent(inout) :: edge_normal_wind

   integer :: jEdge, iEdge, i, j, k
   real,  dimension(:,:), allocatable :: east, north

   ! This code taken from DART models/mpas_atm/model_mod, in turn taken from MPAS

   allocate( east(3,mpas_mesh%nCells))
   allocate( north(3,mpas_mesh%nCells))

   do i = 1, mpas_mesh%nCells
      east(1,i) = -sin(mpas_mesh%lonCell(i))
      east(2,i) =  cos(mpas_mesh%lonCell(i))
      east(3,i) =  0.0 !dble(0.0)
      call r3_normalize(east(1,i), east(2,i), east(3,i))

      north(1,i) = -cos(mpas_mesh%lonCell(i))*sin(mpas_mesh%latCell(i))
      north(2,i) = -sin(mpas_mesh%lonCell(i))*sin(mpas_mesh%latCell(i))
      north(3,i) =  cos(mpas_mesh%latCell(i))
      call r3_normalize(north(1,i), north(2,i), north(3,i))
   enddo

 ! Project data from the cell centers to the edges
   edge_normal_wind = 0.0 ! dble(0.0)
   do i = 1, mpas_mesh%nCells
      do jEdge = 1, mpas_mesh%nEdgesOnCell(i)
         iEdge = edgesOnCell(jEdge, i)
         do k = 1, nz_input
           edge_normal_wind(k,iEdge) = edge_normal_wind(k,iEdge) + 0.5 * zonal(k,i)   &
                     * (edgeNormalVectors(1,iEdge) * east(1,i)  &
                     +  edgeNormalVectors(2,iEdge) * east(2,i)  &
                     +  edgeNormalVectors(3,iEdge) * east(3,i)) &
                     + 0.5 * meridional(k,i)            &
                     * (edgeNormalVectors(1,iEdge) * north(1,i) &
                     +  edgeNormalVectors(2,iEdge) * north(2,i) &
                     +  edgeNormalVectors(3,iEdge) * north(3,i))
         enddo
      enddo
   enddo
   deallocate(east,north)

end subroutine project_to_edges


subroutine r3_normalize(ax, ay, az) !CSS added from DART models/mpas_atm/model_mod, in turn taken from MPAS

   !normalizes the vector (ax, ay, az)

   real, intent(inout) :: ax, ay, az
   real :: mi

   mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
   ax = ax * mi
   ay = ay * mi
   az = az * mi

end subroutine r3_normalize

subroutine derive_edgeNormalVectors(fname, mpas_mesh, cellsOnEdge, edgeNormalVectors)

   character(len=*), intent(in) :: fname
   type(mpas_mesh_type), intent(in)  :: mpas_mesh
   integer, intent(in), dimension(2,mpas_mesh%nEdges) :: cellsOnEdge
   real, intent(inout), dimension(3,mpas_mesh%nEdges) :: edgeNormalVectors

   integer :: iEdge,cell1, cell2, nCells, nEdges
   integer :: ncid
   real, allocatable :: xCell(:), yCell(:), zCell(:)
   real, allocatable :: xEdge(:), yEdge(:), zEdge(:)

   nCells = mpas_mesh%nCells
   nEdges = mpas_mesh%nEdges

   allocate(xCell(nCells),yCell(nCells),zCell(nCells))
   allocate(xEdge(nEdges),yEdge(nEdges),zEdge(nEdges))

   call open_netcdf(trim(fname),ncid)
   call get_netcdf_var(ncid,'xCell',(/1/),(/nCells/),xCell)
   call get_netcdf_var(ncid,'yCell',(/1/),(/nCells/),yCell)
   call get_netcdf_var(ncid,'zCell',(/1/),(/nCells/),zCell)
   call get_netcdf_var(ncid,'xEdge',(/1/),(/nEdges/),xEdge)
   call get_netcdf_var(ncid,'yEdge',(/1/),(/nEdges/),yEdge)
   call get_netcdf_var(ncid,'zEdge',(/1/),(/nEdges/),zEdge)
   call close_netcdf(trim(fname),ncid)

   ! Algorithm is from ...MPAS/src/operators/mpas_vector_operations.F (mpas_initialize_vectors)

   edgeNormalVectors = 0

    ! Initialize normal unit vectors at each edge
    ! These vectors point from cell to cell.
    ! At boundaries, one cell does not exist, so it points from cell to edge or from edge to cell.
   do iEdge = 1,nEdges
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)

      if ( (cell1 == nCells+1) .or. (cell1 == 0) ) then ! this is a boundary edge ! CSS added cell1 == 0 condition
        ! the normal points from the edge location to the cell location
         edgeNormalVectors(1,iEdge) = xCell(cell2) - xEdge(iEdge)
         edgeNormalVectors(2,iEdge) = yCell(cell2) - yEdge(iEdge)
         edgeNormalVectors(3,iEdge) = zCell(cell2) - zEdge(iEdge)
      else if ( (cell2 == nCells+1) .or. (cell2 == 0) ) then ! this is a boundary edge ! CSS added cell2 == 0 condition
        ! the normal points from the cell location to the edge location
         edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(cell1)
         edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(cell1)
         edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell1)
      else ! this is not a boundary cell
        ! the normal points from the cell 1 to cell2
        ! mrp problem: on periodic domains, vectors on edges of domain point the wrong way.
         edgeNormalVectors(1,iEdge) = xCell(cell2) - xCell(cell1)
         edgeNormalVectors(2,iEdge) = yCell(cell2) - yCell(cell1)
         edgeNormalVectors(3,iEdge) = zCell(cell2) - zCell(cell1)
      end if
      call mpas_unit_vec_in_r3(edgeNormalVectors(:,iEdge))
   end do

   deallocate(xCell,yCell,zCell)
   deallocate(xEdge,yEdge,zEdge)

end subroutine derive_edgeNormalVectors

subroutine mpas_unit_vec_in_r3(xin)
    real, dimension(3), intent(inout) :: xin !< Input/Output: Vector and unit vector
    real :: mag
    mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
    xin(:) = xin(:) / mag
end subroutine mpas_unit_vec_in_r3

end module write_data
