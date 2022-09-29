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

        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                           !< counter over atoms
        integer ineigh                          !< counter over structures
        integer ikpoint                         !< counter over kpoints

        integer iscf_iteration
        integer istructure, iseparate
        integer itime_step
        real sigma
        character (len = 25) :: slogfile

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        real timei, timef

! Energies
        real ebs                                 ! band-structure energy
        real uii_uee, uxcdcc                     ! short-range energies
        real etot                                ! total energy

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

              write (s%logfile, *) ' Exchange-correlation assemblers. '
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

! After building the density matrix, then we can free up ewald and denmat arrays
! - we reallocate these during the next SCF cycle anyways.
! We also free up the vna and vxc arrays if this is not converged.
              call destroy_assemble_ewald (s)
              if (sigma .gt. scf_tolerance_set .and.                          &
      &           iscf_iteration .le. max_scf_iterations_set - 1) then
                do iatom = 1, s%natoms
                  do ineigh = 1, s%neighbors(iatom)%neighn
                    deallocate (s%denmat(iatom)%neighbors(ineigh)%block)
                    deallocate (s%vna(iatom)%neighbors(ineigh)%block)
                    deallocate (s%vna(iatom)%neighbors(ineigh)%blocko)
                    deallocate (s%vxc(iatom)%neighbors(ineigh)%block)
                  end do
                  deallocate (s%denmat(iatom)%neighbors)
                  deallocate (s%vna(iatom)%neighbors)
                  deallocate (s%vxc(iatom)%neighbors)
                end do
                deallocate (s%denmat)
                deallocate (s%vna)
                deallocate (s%vxc)

                ! destroy the coefficients, but only if not converged
                do ikpoint = 1, s%nkpoints
                  deallocate (s%kpoints(ikpoint)%c)
                end do
              end if

! End scf loop
              if (sigma .gt. 0.0d0) then
                iscf_iteration = iscf_iteration + 1
              else
                exit
              end if
              if (ifix_CHARGES .eq. 1) exit
            end do

            call writeout_energies (s, ebs, uii_uee, uxcdcc)
            if (iwriteout_xyz .eq. 1) call writeout_xyz (s, ebs, uii_uee, uxcdcc)

! ===========================================================================
! ---------------------------------------------------------------------------
!                               F O R C E S
! ---------------------------------------------------------------------------
! ===========================================================================
! Initialize forces arrays
            call initialize_forces (s)
            call densityPP_matrix (s)
            call cape_matrix (s)

! After building the density matrix, then we can free up the kspace memory
            call destroy_kspace (s)

            write (s%logfile, *)
            write (s%logfile,'(A)') 'Forces '
            write (s%logfile,'(A)') '------ '
            write (s%logfile, *)

! Assemble the derivative blocks needed for forces
            write (s%logfile, *) ' Two-center non-charge dependent Dassemblers.'
            call Dassemble_S (s)
            call Dassemble_T (s)
            call Dassemble_dipole_z (s)
            call Dassemble_svnl (s)
            call Dassemble_vnl_2c (s)

            write (s%logfile, *) ' Two-center charge dependent Dassemblers.'
            call Dassemble_vna_2c (s)
            call Dassemble_ewaldsr (s)
            call Dassemble_ewaldlr (s)

            write (s%logfile,*) ' Three-center non-charge dependent assemblers.'
            call Dassemble_vnl_3c (s)

            write (s%logfile,*) ' Three-center charge dependent Dassemblers.'
            call Dassemble_vna_3c (s)

            write (s%logfile, *) ' Exchange-correlation Dassemblers. '
            call Dassemble_vxc (s)

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

            write (s%logfile, *)
            write (s%logfile, '(A)') ' Grand Total Energy '
            write (s%logfile, '(A)') ' ------------------ '
            write (s%logfile,601) tkinetic
            getot = etot + tkinetic
            write (s%logfile,602) getot
            write (s%logfile,603) getot/s%natoms
            if (itime_step .eq. nstepi) getot_initial = getot/s%natoms
! Check energy conservation
            write (s%logfile,604) 1000.d0*(getot/s%natoms - getot_initial)

! ===========================================================================
! ---------------------------------------------------------------------------
!           Free memory in the molecular dynamics loop
! ---------------------------------------------------------------------------
! ===========================================================================
            call destroy_denmat (s)
            call destroy_Dassemble_2c (s)
            call destroy_assemble_2c (s)
            call destroy_Dassemble_vxc (s)
            call destroy_assemble_vxc (s)
            call destroy_Dassemble_vnl (s)
            call destroy_assemble_PP_2c (s)
            call destroy_Dassemble_ewald (s)
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

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Structure: ', a25)

512     format (2x, 'f_total =',i6 ,3(2x,f15.6))

601     format (2x, '                           Nuclear Kinetic Energy = ', f18.8)
602     format (2x, ' Grand Total Energy (Nuclear Kinetic + Potential) = ', f18.8)
603     format (2x, '                      Grand Total Energy per Atom = ', f18.8)
604     format (2x, '                               deltaE/atom  (meV) = ', f18.8)

! End Program
! ===========================================================================
        stop
        end program fireball
