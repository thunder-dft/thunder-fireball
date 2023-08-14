! copyright info:
!
!                             @Copyright 2022
!                           Fireball Committee
! Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek
! Arizona State University - Otto F. Sankey

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! California Institute of Technology - Brandon Keith
! Czech Institute of Physics - Prokop Hapala
! Czech Institute of Physics - Vladimír Zobač
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Synfuels China Technology Co., Ltd. - Pengju Ren
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

! /GLOBAL
        use M_precision
        use M_welcome

! /SYSTEM
        use M_species
        use M_configuraciones
        use M_kpoints
        use M_neighbors
        use M_neighbors_PP
        use M_atom_functions
        use M_vdW

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
        use M_Dassemble_PP_3c
        use M_Dassemble_vxc
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

! /MPI
        use M_mpi

! /SOCKETS
        use M_sockets

        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
        include '../include/constants.h'

! Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over atoms and neighbors
        integer in1

        integer iscf_iteration
        integer istructure, iseparate
        integer itime_step

        integer nssh                        !< number of shells

        real sigma                          !< difference for SCF
        real rms                            !< RMS of the forces

        character (len = 25) :: slogfile
        character (len = 25) :: sjsonfile

! --------------------------------------------------------------------------
! Socket communication
! --------------------------------------------------------------------------
        ! length of the headers of the driver/wrapper communication protocol
        integer, parameter :: msglen = 12
        integer socket                   ! socket ID & address of the server

        ! degrees of freedom - should be 3*natoms
        ! this is for sanity checking the communications
        integer idof

        ! need to pass double precision energy and forces
        double precision energy

        ! cell information received - currently not used
        double precision, dimension (9) :: xcell
        double precision, dimension (9) :: xcell_inverse

        character (len = 12) :: header

        ! message buffer for sending forces and other information
        double precision, allocatable :: msgbuffer (:)

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin, time_end
        real time_scf_begin, time_scf_end
        real time_forces_begin, time_forces_end
        real time_initial, time_final

! Energies
        real ebs                                 ! band-structure energy
        real efermi                              ! Fermi energy
        real uii_uee, uxcdcc                     ! short-range energies
        real etot                                ! total energy
        real vdW                                 ! van der Waal's energy

        real getot                               ! grand total energy
        real getot_initial                       ! initial grand total

! Temporary for debugging forces
        real etot_previous
        real ftot_previous

        interface
           subroutine Qmixer (t, iscf_iteration, sigma)
             use M_configuraciones
             use M_charges
             implicit none
             integer, intent (in) :: iscf_iteration
             type(T_structure), target :: t
             real, intent (inout) :: sigma
           end subroutine Qmixer

           subroutine writeout_energies (t, ebs, uii_uee, uxcdcc)
             use M_assemble_blocks
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t         ! the structure to be used
             real, intent (in) :: ebs               ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc   ! short-range energies
           end subroutine writeout_energies

           subroutine writeout_xyz (t, ebs, uii_uee, uxcdcc)
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t         ! the structure to be used
             real, intent (in) :: ebs               ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc   ! short-range energies
           end subroutine writeout_xyz

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
          end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Temporary assignment
        etot_previous = 0.0d0
        ftot_previous = 0.0d0

! ===========================================================================
! ---------------------------------------------------------------------------
!                            M P I   B E G I N
! ---------------------------------------------------------------------------
! ===========================================================================
        call initialize_mpi

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

        write (ilogfile,*) ' Reading parameters from the structure input:  '
        write (ilogfile,*) '  structures.inp'
        write (ilogfile,*)
        call read_parameters

        write (ilogfile, *)
        open (unit = 1, file = 'structures.inp', status = 'old')
        read (1, *) nstructures
        if (nstructures .gt. 999) then
          stop ' Cannot calculate more than 999 structures!  '
        end if
        write (ilogfile,'(2x, A, I4, A2)') 'Number of structures to calculate: ', nstructures, '  '
        write (ilogfile, *)
        allocate (structures (nstructures))

        write (ilogfile,*)
        write (ilogfile,'(A)') 'Execution '
        write (ilogfile,'(A)') '========= '
        call initialize_vdW  ! only if structures.vdW exists
        write (ilogfile,*)

! ===========================================================================
! ---------------------------------------------------------------------------
!                          S O C K E T   B E G I N
! ---------------------------------------------------------------------------
! ===========================================================================
! Open socket communiction for interactive dynamics using ipi
! Calls the interface to the POSIX sockets library to open a communication
! channel.
        if (ipi .eq. 1) then
          write (ilogfile, '(A)') 'Open socket for i-pi with ASE communcation'
          call open_socket (socket, inet, port, host)
        end if

! Loop over all structures
! This loop can be made parallel if each subroutine in lightning
! is modified to take s as a parameter, rather than reference s directly
        do istructure = 1, nstructures
          s => structures(istructure)
          read (1, *) s%basisfile
          if (iseparate .eq. 1) then
            s%inpfile = istructure + 1000
            s%logfile = istructure + 2000
            s%jsonfile = istructure + 3000
            slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
            sjsonfile = trim(slogfile)//'.json'
            slogfile = trim(slogfile)//'.log'
            open (unit = s%logfile, file = slogfile, status = 'replace')
            open (unit = s%jsonfile, file = sjsonfile, status = 'replace')
          end if
          write (s%jsonfile,'(A)') '{"fireball":['
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
          call read_vdW_parameters (s)

!         ! SOCKET - GET POSITIONS
          if (ipi .eq. 1) then
            do while (.true.)
              call readbuffer (socket, header, msglen)

              if (trim(header) .eq. "STATUS") then
                call writebuffer (socket, "READY       ", msglen)
                call readbuffer (socket, header, msglen)
                if (trim(header) .eq. "POSDATA") then
                  ! we will read this cell information; however, FIREBALL will
                  ! not use them - useful later when we do constant pressure.
                  call readbuffer (socket, xcell, 9)
                  call readbuffer (socket, xcell_inverse, 9)

                  ! read the actual positions
                  call readbuffer (socket, idof)
                  if (idof .ne. s%natoms) then
                    stop ' socket is communicating wrong number of atoms! '
                  end if
                  allocate (msgbuffer(3*s%natoms))
                  call readbuffer (socket, msgbuffer, 3*s%natoms)
                  deallocate (msgbuffer)
                  exit
                end if
              end if
            end do
          end if

! Molecular-dynamics loop
! ---------------------------------------------------------------------------
! All molecular-dynamics is done by ASE.
! ===========================================================================
          do itime_step = nstepi, nstepf

            ! we write out another step so rewind two lines of the json file and put '},'
            if (itime_step .ne. nstepi) then
              backspace (s%jsonfile)
              backspace (s%jsonfile)
              write (s%jsonfile,'(A)') '},'
            end if

            write (s%jsonfile,'(A)') '{'
            write (s%jsonfile,'(A, I5, A)') '      "nstep":', itime_step, ','
            write (s%jsonfile,'(A)') '      "cell":['
            write (s%jsonfile,'(A, 2x, 3(F15.6, A), A)')                      &
     &        '      [', s%lattice(1)%a(1), ',', s%lattice(1)%a(2), ',',      &
     &                   s%lattice(1)%a(3),'],'
            write (s%jsonfile,'(A, 2x, 3(F15.6, A), A)')                      &
     &        '      [', s%lattice(2)%a(1), ',', s%lattice(2)%a(2), ',',      &
     &                   s%lattice(2)%a(3),'],'
            write (s%jsonfile,'(A, 2x, 3(F15.6, A), A)')                      &
     &        '      [', s%lattice(3)%a(1), ',', s%lattice(3)%a(2), ',',      &
     &                   s%lattice(3)%a(3),']],'

            write (s%jsonfile,'(A)') '      "numbers":['
            do iatom = 1, s%natoms - 1
              in1 = s%atom(iatom)%imass
              write (s%jsonfile,'(16x, i3, A)') species(in1)%nZ, ','
            end do
            in1 = s%atom(s%natoms)%imass
            write (s%jsonfile,'(16x, i3, A)') species(in1)%nZ, '],'

            write (s%jsonfile,'(A)') '      "positions":['
            do iatom = 1, s%natoms - 1
              write (s%jsonfile,'(A, 6x, 3(F15.6, A), A)')                    &
     &          '      [', s%atom(iatom)%ratom(1), ',',                       &
     &                     s%atom(iatom)%ratom(2), ',',                       &
     &                     s%atom(iatom)%ratom(3),'],'
            end do
            write (s%jsonfile,'(A, 6x, 3(F15.6, A), A)')                      &
     &        '      [', s%atom(s%natoms)%ratom(1), ',',                      &
     &                   s%atom(s%natoms)%ratom(2), ',',                      &
     &                   s%atom(s%natoms)%ratom(3),']],'

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
            call driver_neighbors_vdW (s)

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

! ===========================================================================
! ---------------------------------------------------------------------------
!              S C F   L O O P
! ---------------------------------------------------------------------------
! Note that self-consistency is preformed regardless of method used.
! But, in Harris, we just do a single scf cycle.
! ===========================================================================
            sigma = 999.0d0
            iscf_iteration = 1
            do while (sigma .gt. scf_tolerance_set .and.                      &
      &               iscf_iteration .le. max_scf_iterations_set - 1)
              call cpu_time (time_scf_begin)
              write (s%logfile, *)
              write (s%logfile, '(A, I5, A7, I5, A1)') 'Self-Consistent Field step: ', &
                   & iscf_iteration, ' (max: ', max_scf_iterations_set, ')'
              write (s%logfile, '(A)') '----------------------------------------------------'
              write (s%logfile, *)

              write (s%logfile, *) ' Two-center charge dependent assemblers. '
              call assemble_vna_2c (s)
              call assemble_ewaldsr (s)
              call assemble_ewaldlr (s)

              call cpu_time (time_initial)
              write (s%logfile, *)
              write (s%logfile, *) ' Three-center charge dependent assemblers. '
              call assemble_vna_3c (s)

              write (s%logfile, *)
              write (s%logfile, *) ' Exchange-correlation assemblers. '
              call assemble_vxc (s)

              call cpu_time (time_final)
              write (s%logfile, *) ' vna_3c, vxc time: ', time_final - time_initial

! ===========================================================================
! ---------------------------------------------------------------------------
!                         D I A G O N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculating the overlap matrix in K space
              write (s%logfile, *)
              write (s%logfile, *) ' Kspace '
              call cpu_time (time_initial)
              call driver_kspace (s, iscf_iteration)
              call density_matrix (s, efermi)
              if (iwriteout_density .eq. 1) call writeout_density (s)
              call cpu_time (time_final)
              write (s%logfile, *) ' kspace time: ', time_final - time_initial

              if (ifix_CHARGES .ne. 1) then
                call calculate_charges (s)
                call Qmixer (s, iscf_iteration, sigma)
              end if
              if (iwriteout_charges .eq. 1) call writeout_charges (s)

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

              if (sigma .gt. scf_tolerance_set .and.                      &
      &           iscf_iteration .le. max_scf_iterations_set - 1 .and.    &
      &           ifix_CHARGES .ne. 1) then
                write (s%logfile, *) ' Destroy some SCF arrays... '
                call destroy_denmat (s)
                call destroy_assemble_ewald (s)
                call destroy_assemble_vxc (s)
                call destroy_assemble_vna (s)
              end if
              call cpu_time (time_scf_end)
              write (s%logfile, *)
              write (s%logfile, *) ' SCF ENERGY time: ', time_scf_end - time_scf_begin

! After building the density matrix, then we can free up ewald and denmat arrays
! - we reallocate these during the next SCF cycle anyways.
! We also free up the vna and vxc arrays if this is not converged.
              if (sigma .gt. 0.0d0) then
                iscf_iteration = iscf_iteration + 1
              else
                exit
              end if
              if (ifix_CHARGES .eq. 1) exit
            end do

! Write out the charges to .json file
            write (s%jsonfile,'(A)') '      "charges":['
            do iatom = 1, s%natoms - 1
              in1 = s%atom(iatom)%imass
              nssh = species(in1)%nssh
              if (nssh .eq. 1) then
                write (s%jsonfile,'(A, 6x, (F15.6, A), A)')                   &
      &            '      [', s%atom(iatom)%shell(1)%Qin,'],'
              else if (nssh .eq. 2) then
                write (s%jsonfile,'(A, 6x, 2(F15.6, A), A)')                  &
      &            '      [', s%atom(iatom)%shell(1)%Qin, ',',                &
      &                       s%atom(iatom)%shell(2)%Qin,'],'
              else if (nssh .eq. 3) then
                write (s%jsonfile,'(A, 6x, 3(F15.6, A), A)')                  &
      &            '      [', s%atom(iatom)%shell(1)%Qin, ',',                &
      &                       s%atom(iatom)%shell(2)%Qin, ',',                &
      &                       s%atom(iatom)%shell(3)%Qin,'],'
              else if (nssh .eq. 4) then
                write (s%jsonfile,'(A, 6x, 4(F15.6, A), A)')                  &
      &            '      [', s%atom(iatom)%shell(1)%Qin, ',',                &
      &                       s%atom(iatom)%shell(2)%Qin, ',',                &
      &                       s%atom(iatom)%shell(3)%Qin, ',',                &
      &                       s%atom(iatom)%shell(4)%Qin,'],'
              else if (nssh .eq. 5) then
                write (s%jsonfile,'(A, 6x, 5(F15.6, A), A)')                  &
      &            '      [', s%atom(iatom)%shell(1)%Qin, ',',                &
      &                       s%atom(iatom)%shell(2)%Qin, ',',                &
      &                       s%atom(iatom)%shell(3)%Qin, ',',                &
      &                       s%atom(iatom)%shell(4)%Qin, ',',                &
      &                       s%atom(iatom)%shell(5)%Qin,'],'
              end if
            end do
            in1 = s%atom(s%natoms)%imass
            nssh = species(in1)%nssh
            if (nssh .eq. 1) then
              write (s%jsonfile,'(A, 6x, (F15.6, A), A)')                     &
      &          '      [', s%atom(s%natoms)%shell(1)%Qin,']],'
            else if (nssh .eq. 2) then
              write (s%jsonfile,'(A, 6x, 2(F15.6, A), A)')                    &
      &          '      [', s%atom(s%natoms)%shell(1)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(2)%Qin,']],'
            else if (nssh .eq. 3) then
              write (s%jsonfile,'(A, 6x, 3(F15.6, A), A)')                    &
      &          '      [', s%atom(s%natoms)%shell(1)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(2)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(3)%Qin,']],'
            else if (nssh .eq. 4) then
              write (s%jsonfile,'(A, 6x, 4(F15.6, A), A)')                    &
      &          '      [', s%atom(s%natoms)%shell(1)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(2)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(3)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(4)%Qin,']],'
            else if (nssh .eq. 5) then
              write (s%jsonfile,'(A, 6x, 5(F15.6, A), A)')                    &
      &          '      [', s%atom(s%natoms)%shell(1)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(2)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(3)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(4)%Qin, ',',               &
      &                     s%atom(s%natoms)%shell(5)%Qin,']],'
            end if
            call writeout_energies (s, ebs, uii_uee, uxcdcc)

            ! json output for Fermi energy
            write (s%jsonfile,'(A, F15.6, A)') '      "fermi":', efermi, ','

            ! json output for energy
            write (s%jsonfile,'(A, F15.6, A)') '      "energy":', etot, ','

! ===========================================================================
! ---------------------------------------------------------------------------
!                               F O R C E S
! ---------------------------------------------------------------------------
! ===========================================================================
            call cpu_time (time_forces_begin)
            call initialize_forces (s)
            call densityPP_matrix (s)
            call cape_matrix (s)

! After building the density matrix, then we can free up the kspace memory
            call destroy_kspace (s)

            write (s%logfile, *)
            write (s%logfile,'(A)') 'Forces '
            write (s%logfile,'(A)') '------ '

! Assemble the derivative blocks needed for forces
            write (s%logfile, *) ' Two-center non-charge dependent Dassemblers.'
            call Dassemble_S (s)
            call Dassemble_T (s)
            call Dassemble_dipole_z (s)
            call Dassemble_svnl (s)
            call Dassemble_vnl_2c (s)

            write (s%logfile, *)
            write (s%logfile, *) ' Two-center charge dependent Dassemblers.'
            call Dassemble_vna_2c (s)
            call Dassemble_ewaldsr (s)
            call Dassemble_ewaldlr (s)

            write (s%logfile, *)
            write (s%logfile,*) ' Three-center non-charge dependent assemblers.'
            call Dassemble_vnl_3c (s)

            call cpu_time (time_initial)
            write (s%logfile,*) ' Three-center charge dependent Dassemblers.'
            call Dassemble_vna_3c (s)

            write (s%logfile, *)
            write (s%logfile, *) ' Exchange-correlation Dassemblers. '
            call Dassemble_vxc (s)
            call cpu_time (time_final)
            write (s%logfile, *) ' vxc, vna_3c forces time: ', time_final - time_initial

! short-range interactions (double-counting interactions)
            call Dassemble_uee (s)
            call Dassemble_uxc (s)

            call build_forces (s, rms)

! Add in the van der Waals energy and forces if there is an input file
            call calculate_vdW (s, vdW)

            if (iwriteout_forces .eq. 1) call writeout_forces (s)
            write (s%logfile,*)
            write (s%logfile,*) ' Total Forces:'
            do iatom = 1, s%natoms
              write (s%logfile, 512)  iatom, s%forces(iatom)%ftot
            end do

            ! json output for forces
            write (s%jsonfile,'(A)') '      "forces":['
            do iatom = 1, s%natoms - 1
              write (s%jsonfile,'(A, 3x, 3(F15.6, A), A)')                    &
     &          '      [', s%forces(iatom)%ftot(1), ',',                      &
     &                     s%forces(iatom)%ftot(2), ',',                      &
     &                     s%forces(iatom)%ftot(3),'],'
            end do
            write (s%jsonfile,'(A, 3x, 3(F15.6, A), A)')                      &
     &        '      [', s%forces(s%natoms)%ftot(1), ',',                     &
     &                   s%forces(s%natoms)%ftot(2), ',',                     &
     &                   s%forces(s%natoms)%ftot(3),']],'
            write (s%jsonfile,'(A, F15.6, A)') '      "RMS":', rms
            write (s%jsonfile,'(A)') '}'
            write (s%jsonfile,'(A)') ']}'

            call cpu_time (time_forces_end)
            write (s%logfile, *)
            write (s%logfile, *) ' FORCES time: ', time_forces_end - time_forces_begin

            write (s%logfile, *)
            write (s%logfile, '(A)') ' Grand Total Energy '
            write (s%logfile, '(A)') ' ------------------ '
            write (s%logfile,601) tkinetic
            write (s%logfile,602) vdW
            getot = etot + tkinetic + vdW
            write (s%logfile,603) getot
            write (s%logfile,604) getot/s%natoms
            if (itime_step .eq. nstepi) getot_initial = getot/s%natoms
! Check energy conservation
            write (s%logfile,605) 1000.d0*(getot/s%natoms - getot_initial)

! ===========================================================================
! ---------------------------------------------------------------------------
!                            S O C K E T S
! ---------------------------------------------------------------------------
! ===========================================================================
            if (ipi .eq. 1) then
              allocate (msgbuffer(3*s%natoms))
              do while (.true.)
                call readbuffer (socket, header, msglen)

                ! SOCKET - SEND FORCES
                if (trim(header) .eq. "STATUS") then
                  call writebuffer (socket, "HAVEDATA    ", msglen)
                  call readbuffer (socket, header, msglen)
                  if (trim(header) .eq. "GETFORCE") then
                    call writebuffer (socket, "FORCEREADY  ", msglen)
                    energy = etot/P_Hartree
                    call writebuffer (socket, energy)
                    call writebuffer (socket, s%natoms)
                    do iatom = 1, s%natoms
                      msgbuffer(3*(iatom-1)+1:3*iatom) = s%forces(iatom)%ftot*(P_abohr/P_Hartree)
                    end do
                    call writebuffer (socket, msgbuffer, 3*s%natoms)

                    ! sends the cell and fantasy friction information
                    ! for constant pressure - currently not used
                    call writebuffer (socket, xcell, 9)
                    call writebuffer (socket, 1)
                    call writebuffer (socket, ' ', 1)
                    exit
                  end if
                end if
              end do

              ! SOCKET - GET POSITIONS
              do while (.true.)
                call readbuffer (socket, header, msglen)
                if (trim(header) .eq. "STATUS") then
                  call writebuffer (socket, "READY       ", msglen)
                  call readbuffer (socket, header, msglen)
                  if (trim(header) .eq. "POSDATA") then
                    ! we will read this cell information; however, FIREBALL will
                    ! not use them - useful later when we do constant pressure.
                    call readbuffer (socket, xcell, 9)
                    call readbuffer (socket, xcell_inverse, 9)

                    ! read the actual positions
                    call readbuffer (socket, idof)
                    if (idof .ne. s%natoms) then
                      stop ' socket is communicating wrong number of atoms! '
                    end if
                    call readbuffer (socket, msgbuffer, 3*s%natoms)
                    do iatom = 1, s%natoms
                     s%atom(iatom)%ratom = msgbuffer(3*(iatom-1)+1:3*iatom)*P_abohr
                    end do
                    exit
                  end if
                end if
              end do
              deallocate (msgbuffer)
            end if

            ! because we are using ASE, then ratom never really gets shifted
            if (ishiftO .eq. 1) then
              do iatom = 1, s%natoms
                s%atom(iatom)%ratom = s%atom(iatom)%ratom + shifter
              end do
            end if
            call writeout_xyz (s, ebs, uii_uee, uxcdcc)
            if (ishiftO .eq. 1) then
              do iatom = 1, s%natoms
                s%atom(iatom)%ratom = s%atom(iatom)%ratom - shifter
              end do
            end if

! ===========================================================================
! ---------------------------------------------------------------------------
!           Free memory in the molecular dynamics loop
! ---------------------------------------------------------------------------
! ===========================================================================
            call destroy_denmat (s)
            call destroy_denmat_PP (s)
            call destroy_Dassemble_2c (s)
            call destroy_assemble_2c (s)
            call destroy_assemble_vna (s)
            call destroy_Dassemble_vxc (s)
            call destroy_assemble_vxc (s)
            call destroy_Dassemble_PP_2c (s)
            call destroy_assemble_PP_2c (s)
            call destroy_assemble_ewald (s)
            call destroy_Dassemble_ewald (s)
            call destroy_forces (s)
            call destroy_neighbors (s)
            call destroy_neighbors_PP (s)
            call destroy_neighbors_vdW (s)
          end do ! end molecular dynamics loop

! ===========================================================================
! ---------------------------------------------------------------------------
!    P O S T   O P T I M I Z A T I O N   C H A R A C T E R I Z A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculate the absorption spectrum.
!         if (iwriteout_abs .eq. 1) call absorption (s)

! Destroy final arrays
          call destroy_charges (s)
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

! ===========================================================================
! ---------------------------------------------------------------------------
!                             M P I   E N D
! ---------------------------------------------------------------------------
! ===========================================================================
        call finalize_mpi

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Structure: ', a25)

512     format (2x, 'f_total =',i6 ,3(2x,f15.6))

601     format (2x, '                           Nuclear Kinetic Energy = ', f18.8)
602     format (2x, '                            van der Waal''s Energy = ', f18.8)
603     format (2x, ' Grand Total Energy (Nuclear Kinetic + Potential) = ', f18.8)
604     format (2x, '                      Grand Total Energy per Atom = ', f18.8)
605     format (2x, '                               deltaE/atom  (meV) = ', f18.8)

! End Program
! ===========================================================================
        stop
        end program fireball
