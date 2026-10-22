! Code to read the namelist file
module program_setup

use esmf
use ESMF_LogPublicMod
use utils_mod, only : error_handler

implicit none

private

! Private variables only visible to this module
integer, parameter :: nvars_max  = 200

! Public variables accessible to other modules
type(ESMF_LogKind_Flag), public :: LogType
integer, public :: nvars_to_blend
 
! Public namelist &share variables
character(len=500), public  :: large_scale_file = "NULL" ! Full path of MPAS file with large-scale information
character(len=500), public  :: small_scale_file = "NULL" ! Full path of MPAS file with small-scale information
character(len=500), public  :: output_blended_filename = "./blended.nc" ! Full path of output file to be created with blended fields
character(len=500), public  :: grid_info_file = "NULL"   ! Full path of text file with information about meshes
character(len=500), public  :: interp_method = "bilinear" ! Interpolation method
character(len=500), public  :: interp_method_going_up = "NULL" ! Interpolation method when coarsening (going upscale)
character(len=500), public  :: interp_method_going_down = "NULL" ! Interpolation method when going downscale (coarse to fine mesh)
character(len=500), public  :: extrap_method = "mpas_nearest"     ! Extrapolation method; default is mpas_nearest
integer, public             :: extrap_num_levels_creep = 64       ! extrap_method = 'esmf_creep' or 'esmf_creep_nrst_d' ; number of levels for creep fill
character(len=500), public  :: variables_to_blend(nvars_max) = "NULL" ! netCDF variable names to blend
logical, public             :: average_upscale_before_interp = .false. ! If true, average to the next coarsest mesh before interpolating to that mesh
logical, public             :: output_intermediate_files_up = .false. ! If true, output a file for each intermediate mesh when progressively upscaling
logical, public             :: output_intermediate_files_down = .false. ! If true, output a file for each intermediate mesh when progressively downscaling
logical, public             :: output_upscaled_data_on_native_mesh = .false. ! If true, output a file for each intermediate mesh interpolated onto the original native mesh
logical, public             :: output_blended_edge_normal_wind = .false. ! If true, output blended edge normal wind (u)
logical, public             :: esmf_log = .false. ! logging information
logical, public             :: smooth_going_downscale = .false. ! If true, run a smoother after interpolating downscale when reconstructing fields
real,    public             :: smoother_dimensionless_coefficient = 0.0 ! Not used if smooth_going_downscale = .false.; parameter for smoothing
logical, public             :: interp_cartesian_winds = .false. ! If true, convert zonal/meridional winds (if provided) to cartesian winds before interpolation, and then converted back to zonal/meridional winds for output 

! Public namelist &latlon_output variables
logical, public            :: output_latlon_grid = .false. ! Switch to turn on lat-lon output (for diagnostics purposes)
logical, public            :: is_regional = .false. ! If lat-lon mesh is regional, set to .true.  If global, set to .false.
character(len=500), public  :: extrap_method_latlon = "none"     ! Extrapolation method for outputting lat/lon grid. Defaults to none. If using creep, uses extrap_num_levels_creep.
integer, public            :: nlat = -1 ! Number of latitude points on output lat-lon grid
integer, public            :: nlon = -1 ! Number of longitude points on output lat-lon grid
real, public               :: lat_ll = -1.0 ! Latitude (degrees) of lower left corner of output lat-lon grid
real, public               :: lon_ll = -1.0 ! Longitude (degrees) of lower left corner of output lat-lon grid
real, public               :: dx_in_degrees = -9999. ! Horizontal grid spacing (degrees) for longitude in the output lat-lon grid
real, public               :: dy_in_degrees = -9999. ! Horizontal grid spacing (degrees) for latitude  in output lat-lon grid

! Public subroutines
public :: read_setup_namelist

! Namelists
namelist /share/ large_scale_file, small_scale_file, output_blended_filename, &
         variables_to_blend, grid_info_file, interp_method, interp_method_going_up, interp_method_going_down,extrap_method, &
         extrap_num_levels_creep, average_upscale_before_interp, output_intermediate_files_up, &
         output_intermediate_files_down, output_upscaled_data_on_native_mesh, output_blended_edge_normal_wind, esmf_log, &
         smooth_going_downscale, smoother_dimensionless_coefficient, interp_cartesian_winds 

namelist /latlon_output/ output_latlon_grid, is_regional, extrap_method_latlon, &
                         nlat, nlon, lat_ll, lon_ll, dx_in_degrees, dy_in_degrees

contains

subroutine read_setup_namelist(filename,localpet)

   character(len=*), intent(in) :: filename
   integer, intent(in) :: localpet

   integer :: i, ierr, unum
 
   ! Read the namelist
   unum = 41
   open(unum, file=filename, iostat=ierr)
   if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
   read(unum, nml=share, iostat=ierr)
   if (ierr /= 0) call error_handler("READING SETUP NAMELIST SHARE.", ierr)
   read(unum, nml=latlon_output, iostat=ierr)
   if (ierr /= 0) call error_handler("READING SETUP NAMELIST LATLON_OUTPUT.", ierr)
   close (unum)

   ! Figure out how many variables there are to blend
   nvars_to_blend = 0
   do i = 1,nvars_max
     if ( variables_to_blend(i) == "NULL" ) exit
     nvars_to_blend = nvars_to_blend + 1
   enddo
   if ( nvars_to_blend == 0 ) call error_handler("NO VARIABLES TO BLEND.", 1)

   ! Adding an extra variable to blend - wCartesian
   if ( interp_cartesian_winds ) then
       nvars_to_blend = nvars_to_blend + 1
       variables_to_blend(nvars_to_blend) = 'wCartesian'
   endif

   ! Define the ESMF log type
   if (esmf_log) then
     LogType = ESMF_LOGKIND_MULTI_ON_ERROR
   else
     LogType = ESMF_LOGKIND_NONE
   endif 

   ! Set dy to dx if dx is set but dy isn't
   if ( dy_in_degrees == -9999. ) then
      dy_in_degrees = dx_in_degrees
      if (localpet == 0 ) write(*,*)'setting dy_in_degrees to dx_in_degrees = ',dy_in_degrees
   endif

   ! If interp_method_going_up isn't set in the namelist, set it to interp_method (for backwards compatibility)
   if ( interp_method_going_up == "NULL" ) then
      interp_method_going_up = interp_method
      if (localpet == 0 ) write(*,*)'interp_method_going_up unset; setting to interp_method, which is ',interp_method
   endif

   ! If interp_method_going_down isn't set in the namelist, set it to interp_method (for backwards compatibility)
   if ( interp_method_going_down == "NULL" ) then
      interp_method_going_down = interp_method
      if (localpet == 0 ) write(*,*)'interp_method_going_down unset; setting to interp_method, which is ',interp_method
   endif

   if ( interp_cartesian_winds ) then 
      if (localpet == 0 ) write(*,*) 'Conversion from zonal/meridional winds to Cartesian wind vectors prior to interpolation is turned on.' 
   endif



end subroutine read_setup_namelist

end module program_setup
