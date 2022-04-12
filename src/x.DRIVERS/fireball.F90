#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

! copyright info:
!
!                                    @Copyright 2016
!                                  Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! Program Description
! ===========================================================================
!> This is the main driver for FIREBALL.
! ==========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        program fireball
        use M_precision

! /SYSTEM
        use M_species
        use M_configuraciones
        use M_neighbors
        use M_neighbors_PP
        use M_atom_functions

! /FDATA
        use M_Fdata_1c
        use M_Fdata_2c

! /ASSEMBLERS
        use M_assemble_2c
        use M_assemble_3c
        use M_assemble_ewald
        use M_assemble_usr
        use M_assemble_vxc
        use M_assemble_PP_2c
        use M_assemble_PP_3c

! /DASSEMBLERS
        use M_Dassemble_2c
        use M_Dassemble_PP_2c
        use M_Dassemble_3c
        use M_Dassemble_vxc
        use M_Dassemble_vxc_3c
        use M_Dassemble_usr
        use M_Dassemble_ewald
        use M_build_forces
! /MD
        use M_dynamics

! /SOLVESH
        use M_kspace
        use M_density_matrix

! /SCF
        use M_charges

        implicit none

        interface
           subroutine Qmixer (t, iscf_iteration, sigma)
             use M_configuraciones
             implicit none
             integer, intent (in) :: iscf_iteration
             type(T_structure), target :: t
             real, intent (inout) :: sigma
           end subroutine Qmixer

           subroutine absorption (t)
             use M_species
             use M_configuraciones
             use M_atom_functions
             implicit none
             type(T_structure), target :: t         !< the structure to be used
           end subroutine absorption

           subroutine dos (t)
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t         !< the structure to be used
           end subroutine dos

           subroutine writeout_xyz (t, ebs, uii_uee, uxcdcc)
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t         ! the structure to be used
             real, intent (in) :: ebs               ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc   ! short-range energies
           end subroutine writeout_xyz

           subroutine writeout_energies (t, ebs, uii_uee, uxcdcc)
             use M_assemble_blocks
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t         ! the structure to be used
             real, intent (in) :: ebs               ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc   ! short-range energies
           end subroutine writeout_energies

          end interface

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                        !< counter over atoms
        integer in1                          !< defines species number
        integer iscf_iteration
        integer istructure, iseparate
        integer itime_step
        real sigma
        character (len = 25) :: slogfile

        logical :: file_exists

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        real timei, timef

! Energies
        real ebs                                     ! band-structure energy
        real uii_uee, uxcdcc                         ! short-range energies
        real etot                                    ! total energy

! Temporary for debugging forces
        real etot_previous

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Temporary assignment
        etot_previous = 0.0d0

! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call cpu_time (time_begin)
        iseparate = 1
        open (unit = ilogfile, file = 'output.log', status = 'replace')
        call welcome_fireball

! ===========================================================================
! ---------------------------------------------------------------------------
!             R E A D   I N   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
        write (ilogfile,'(A)') 'Fdata Setup '
        write (ilogfile,'(A)') '=========== '
        write (ilogfile,*)
        call read_Fdata_location
        call read_info

        write (ilogfile,'(A)') 'Hamiltonian Interactions (Fdata) '
        write (ilogfile,'(A)') '================================ '
        write (ilogfile,*)
        call read_Fdata_1c
        call read_Fdata_2c
        call read_Fdata_3c

! Read in the wavefunctions
        write (ilogfile,*)
        write (ilogfile,'(A)') 'Sankey-Niklewski wave-functions (Fdata) '
        write (ilogfile,'(A)') '======================================= '
        call read_wavefunctions

! Read parameters from structures.inp file
        write (ilogfile,'(A)') 'Structures '
        write (ilogfile,'(A)') '========== '
        write (ilogfile,*)

        ! file_exists will be TRUE if the file
        INQUIRE(FILE="structures.inp", EXIST=file_exists)   
                                                            
        ! exists and FALSE otherwise
        if ( file_exists ) then
           open (unit = 1, file = 'structures.inp', status = 'old')
        else
           write(*,*) 'ERROR: Could not open: "structures.inp"'
           call exit(1)
        end if

        read (1, *) nstructures
        if (nstructures .gt. 999) then
          stop ' Cannot calculate more than 999 structures!  '
        end if
        write (ilogfile,'(2x, A, I4, A2)') 'Number of structures to calculate: ', nstructures, '  '
        write (ilogfile, *)
        allocate (structures (nstructures))

        write (ilogfile,*) ' Reading parameters from the structure input:  '
        write (ilogfile,*) '  structures.inp'
        write (ilogfile,*)
        call read_parameters

        write (ilogfile,*)
        write (ilogfile,'(A)') 'Execution '
        write (ilogfile,'(A)') '========= '
        write (ilogfile,*)

! Loop over all structures
! This loop can be made parallel if each subroutine in lightning
! is modified to take s as a parameter, rather than reference s directly
!$omp parallel do private (s, slogfile, sigma, iscf_iteration, timei, timef) &
!$omp             private (ebs, uii_uee, uxcdcc, etot)
        do istructure = 1, nstructures
          s => structures(istructure)
          read (1, *) s%basisfile
          if (iseparate .eq. 1) then
            s%logfile = istructure + 1000
            s%inpfile = istructure + 2000
            slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
            slogfile = trim(slogfile)//'.log'
            open (unit = s%logfile, file = slogfile, status = 'replace')
          end if
          write (ilogfile, 100) s%basisfile

          write (s%logfile,'(A)') 'Structure'
          write (s%logfile,'(A)') '========='
          write (s%logfile,*)

          write (s%logfile, *) ' Structure = ', istructure
          write (s%logfile, *)

          ! Read in the coordinates and parameters
          call read_positions (s)

          ! Set the charges for istructure
          call read_charges (s)
          call set_constraints (s)

! Molecular-dynamics loop
! ---------------------------------------------------------------------------
! ===========================================================================
! FIXME - We set nstepi and nstepf here!
          !nstepi = 1
          !nstepf = 100
          call set_gear ()
          do itime_step = nstepi, nstepf
            write (s%logfile, *)
            write (s%logfile, '(A, I5, A1, I5, A1)') 'Molecular-Dynamics Loop  Step: (', itime_step, '/', nstepf, ')'
            write (s%logfile, '(A)') '==============================================='
            write (s%logfile, *)

! ===========================================================================
! ---------------------------------------------------------------------------
!           N E I G H B O R S   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
            call driver_neighbors (s)
            call driver_neighbors_PP (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!              A S S E M B L E   T H E   H A M I L T O N I A N
! ---------------------------------------------------------------------------
! ===========================================================================
! Assemble the Hamiltonian matrix:
            write (s%logfile, *)
            write (s%logfile, *) ' Two-center non-charge dependent assemblers. '
            call assemble_S (s)
            call assemble_T (s)
            call assemble_dipole_z (s)
            call assemble_svnl (s)
            call assemble_vnl_2c (s)

            write (s%logfile,*) ' Three-center non-charge dependent assemblers.'
            call assemble_vnl_3c (s)

! Put scf loop here
            sigma = 999.0d0
            iscf_iteration = 1
            do while (sigma .gt. scf_tolerance_set .and.                     &
      &               iscf_iteration .le. max_scf_iterations_set)
              write (s%logfile, *)
              write (s%logfile, '(A, I5, A7, I5, A1)') 'Self-Consistent Field step: ', &
                   & iscf_iteration, ' (max: ', max_scf_iterations_set, ')'
              write (s%logfile, '(A)') '----------------------------------------------------'
              write (s%logfile, *)
              write (s%logfile, *) ' Two-center charge dependent assemblers. '
              call assemble_vna_2c (s)
              call assemble_ewaldsr (s)
              call assemble_ewaldlr (s)
              write (s%logfile, *) ' Three-center charge dependent assemblers. '
              call assemble_vna_3c (s)
              call assemble_vxc (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                         D I A G O N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculating the overlap matrix in K space
              write (s%logfile, *) ' Kspace '
              call cpu_time(timei)
              call driver_kspace (s, iscf_iteration)
              call cpu_time(timef)
              write (s%logfile, *)
              write (s%logfile, *) ' kspace time: ', timef - timei

              call density_matrix (s)
              if (iwriteout_density .eq. 1) call writeout_density (s)

              call calculate_charges (s)
              if (iwriteout_charges .eq. 1) call writeout_charges (s)
              call Qmixer (s, iscf_iteration, sigma)

! ===========================================================================
! ---------------------------------------------------------------------------
!                       T O T A L   E N E R G I E S
! ---------------------------------------------------------------------------
! ===========================================================================
! short-range interactions (double-counting interactions)
              call calculate_ebs (s, ebs)
              uii_uee = 0.0d0; uxcdcc = 0.0d0
              call assemble_uee (s, uii_uee)
              call assemble_uxc (s, uxcdcc)

              ! Evaluate total energy
              etot = ebs + uii_uee + uxcdcc

              call writeout_energies (s, ebs, uii_uee, uxcdcc)
              if (iwriteout_xyz .eq. 1) call writeout_xyz (s, ebs, uii_uee, uxcdcc)

! End scf loop
              if (sigma .gt. 0.0d0) then
                iscf_iteration = iscf_iteration + 1
              else
                exit
              end if
            end do ! End scf loop
            etot_previous = etot

! ===========================================================================
! ---------------------------------------------------------------------------
!                               F O R C E S
! ---------------------------------------------------------------------------
! ===========================================================================
! Initialize forces arrays
            call initialize_forces (s)
            call densityPP_matrix (s)
            call cape_matrix (s)

            write (s%logfile, *)
            write (s%logfile,'(A)') 'Forces '
            write (s%logfile,'(A)') '------ '
            write (s%logfile, *)

! Assemble the derivative blocks needed for forces
            write (s%logfile, *) ' Two-center non-charge dependent Dassemblers.'
            call Dassemble_S (s)
            call Dassemble_T (s)
            call Dassemble_svnl (s)
            call Dassemble_vnl_2c (s)
            call Dassemble_vna_2c (s)
            call Dassemble_rho_2c (s)
            call Dassemble_vxc (s)
            call Dassemble_ewaldsr (s)

            write (s%logfile,*) ' Three-center non-charge dependent Dassemblers.'
            call Dassemble_vna_3c (s)
            call Dassemble_vxc_3c (s)

! short-range interactions (double-counting interactions)
            call Dassemble_uee (s)
            call Dassemble_uxc (s)

            call build_forces (s)
            if (iwriteout_forces .eq. 1) call writeout_forces (s)

            write (s%logfile,*)
            write (s%logfile,*) ' Total Forces: '
            do iatom = 1, s%natoms
              write (s%logfile, 512)  iatom, s%forces(iatom)%ftot
            end do
            call md (s, itime_step)
            
! Output the coordinates to a .xyz file
            slogfile = s%basisfile(:len(trim(s%basisfile))-4)
            slogfile = trim(slogfile)//'.xyz'
            open (unit = s%inpfile, file = slogfile, status = 'unknown',     &
     &            position = 'append')
            write (s%inpfile, *) s%natoms
            write (s%inpfile, *) etot, T_instantaneous
            do iatom = 1, s%natoms
              in1 = s%atom(iatom)%imass
              write (s%inpfile,*) species(in1)%symbol, s%atom(iatom)%ratom
            end do
            close (unit = s%inpfile)

! ===========================================================================
! ---------------------------------------------------------------------------
!           Free memory in the molecular dynamics loop
! ---------------------------------------------------------------------------
! ===========================================================================
            call destroy_kspace (s)
            call destroy_denmat (s)
            call destroy_neighbors (s)
            call destroy_neighbors_PP (s)
          end do ! end molecular dynamics loop

! ===========================================================================
! ---------------------------------------------------------------------------
!    P O S T   O P T I M I Z A T I O N   C H A R A C T E R I Z A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculate the absorption spectrum.
          if (iwriteout_abs .eq. 1) call absorption (s)

! Destroy final arrays
          !call destroy_kspace (s)
          !call destroy_denmat (s)
          call destroy_charges (s)

          !call destroy_neighbors (s)
          !call destroy_neighbors_PP (s)

          ! destroy neighbors last
          deallocate (s%xl) ! where to put this?

! Calculate the electronic density of states.
! We do this after destroying some arrays so that we can optimize the
! memory usage.
          if (iwriteout_dos .eq. 1) call dos (s)

          write (s%logfile,*)
          write (s%logfile,'(A)') 'FIREBALL EXECUTION COMPLETE'
          if (iseparate .eq. 1) close (s%logfile)

        end do ! end loop over all structures
        close (1)
        write (ilogfile,*)

! ===========================================================================
! ---------------------------------------------------------------------------
!               D E S T R O Y   A R R A Y S - F I N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Destroy datafile storage
        call destroy_Fdata_1C
        call destroy_Fdata_2C
        call destroy_Fdata_3c

! Destroy SYSTEM information.
        call destroy_positions
        call destroy_species

        call cpu_time (time_end)

        write (ilogfile,'(A, F9.2, A)') 'FIREBALL RUNTIME : ',               &
     &    time_end-time_begin, ' [sec]  '
        write (*,'(A, F9.2, A)') 'FIREBALL RUNTIME : ',                      &
     &    time_end-time_begin, ' [sec]  '
        write (ilogfile,'(A)') 'FIREBALL EXECUTION COMPLETE'
        close (ilogfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Structure: ', a25)

512     format (2x, 'ftot =',i6 ,3(2x,f15.6))

! End Program
! ===========================================================================
        stop
        end program fireball
