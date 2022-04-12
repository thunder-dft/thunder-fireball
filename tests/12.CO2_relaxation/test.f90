program test
  implicit none
  
  type Tspecies
     integer :: nZ
  end type Tspecies

  logical :: file_exists
  integer :: nspecies
!  Tspecies, pointer :: species
  integer :: ispecies
  integer :: nZ
  CHARACTER (len=128) :: Fdata_location

  write (*,*) ' Reading: Fdata.inp  '

  INQUIRE(FILE="Fdata.inp", EXIST=file_exists)   ! file_exists will be TRUE if the file                                                                    
  ! exists and FALSE otherwise                                                                              
  if ( file_exists ) then
     open (unit = 11, file = 'Fdata.inp', status = 'old')
  else
     write(*,*) 'ERROR: Could not open: "Fdata.inp"'
     call exit(1)
  end if

  read (11,*) nspecies
!  allocate (species (nspecies))
  do ispecies = 1, nspecies
!     read (11,*) species(ispecies)%nZ
     read (11,*) nZ
     write(*,*) nZ
  end do

  read (11,*) Fdata_location
  write (*,*) Fdata_location
  close (11)
end program test

