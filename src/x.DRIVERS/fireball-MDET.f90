! copyright info:
!
!                             @Copyright 2012
!                           Fireball Committee
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
!> This is the main driver for the LIGHTNING version of FIREBALL.
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
        program fireball_MDET

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

! /SOLVESH
        use M_kspace
        use M_density_matrix

! /SCF
        use M_charges

! /DASSEMBLERS
        use M_Dassemble_2c
        use M_Dassemble_PP_2c
        use M_Dassemble_usr
        use M_build_forces

! /MD
        use M_dynamics

!/ NAC
        use M_non_adiabatic

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
             type(T_structure), target :: t           !< the structure to be used
           end subroutine absorption

           subroutine dos (t)
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t           !< the structure to be used
           end subroutine dos

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
        integer iscf_iteration               !< iteration step
        integer istep                        !< time step
        integer istructure, iseparate
        integer itime_step
        integer icurrent_state         !< current state of excited electron

        real sigma
        real time                       !< simulation time

        character (len = 25) :: slogfile

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        real timei, timef

! Energies
        real atomic_energy                   ! total atomic energy
        real ebs                             ! band-structure energy
        real uii_uee, uxcdcc                 ! short-range energies
        real etot                            ! total energy
        real etot_per_atom                   ! total energy per atom
        real getot                           ! grand total energy
        real getot_per_atom                  ! grand total energy per atom
        real getot_initial                   ! initial grand total
        real deltaE                          ! energy difference (conservation)

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call cpu_time (time_begin)
        iseparate = 1
        open (unit = ilogfile, file = 'output.log', status = 'replace')
        call welcome

! ===========================================================================
! ---------------------------------------------------------------------------
!             R E A D   I N   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
        call read_Fdata_location
        call read_info
        call read_Fdata_1c
        call read_Fdata_2c
        call read_Fdata_3c

        write (ilogfile, *)
        open (unit = 1, file = 'structures.inp', status = 'old')
        read (1, *) nstructures
        if (nstructures .gt. 999) then
          stop ' Cannot calculate more than 999 structures! '
        end if
        write (ilogfile, *) ' Number of structure calculating = ', nstructures
        allocate (structures (nstructures))

! Read parameters from structures.inp file
        write (ilogfile,*)
        write (ilogfile,*) ' Reading parameters from the structure input. '
        call read_parameters

! Loop over all structures
! This loop can be made parallel if each subroutine in lightning
! is modified to take s as a parameter, rather than reference s directly
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
          write (ilogfile, 100) slogfile
          write (s%logfile, *) ' Structure = ', istructure

          ! read positions for istructure
          call read_positions (s)

          ! Set the charges for istructure
          call read_charges (s)

          call set_constraints (s)
          call initialize_mdet (s, icurrent_state)

! Molecular-dynamics loop
! ---------------------------------------------------------------------------
! ===========================================================================
          itime_step = 0
          time = 0.0d0
          do istep = nstepi, nstepf
             itime_step = itime_step + 1
             write (s%logfile,*)
             write (s%logfile,101)
             write (s%logfile,*) ' M O L E C U L A R   D Y N A M I C S   S T E P = ', itime_step
             write (s%logfile,101)

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
            write (s%logfile, *) ' Calling two-center non-charge dependent assemblers. '
            call assemble_S (s)
            call assemble_T (s)
            call assemble_dipole_z (s)
            call assemble_svnl (s)
            call assemble_vnl_2c (s)

            write (s%logfile,*) ' Calling three-center non-charge dependent assemblers. '
            call assemble_vnl_3c (s)

! Put scf loop here
            sigma = 999.0d0
            iscf_iteration = 1
            do while (sigma .gt. scf_tolerance_set .and.                  &
      &                 iscf_iteration .le. max_scf_iterations_set)
              write (s%logfile, *)
              write (s%logfile, *) ' Begin scf step = ', iscf_iteration
              write (s%logfile, *) ' Calling two-center charge dependent assemblers. '
              call assemble_vna_2c (s)
              call assemble_ewaldsr (s)
              call assemble_ewaldlr (s)

              write (s%logfile, *) ' Calling three-center charge dependent assemblers. '
              call assemble_vna_3c (s)
              call assemble_vxc (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                         D I A G O N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculating the overlap matrix in K space
              write (s%logfile, *) ' Calling kspace: '
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
              call calculate_ebs (s, ebs)

! Add in short-range (double-counting) contributions
              call assemble_uee (s, uii_uee)
              call assemble_uxc (s, uxcdcc)

              ! Writing out the energy pieces
              write (s%logfile, *)
              write (s%logfile, *) ' ---------- T H E  T O T A L  E N E R G Y ----------- '
              write (s%logfile, 202) ebs
              write (s%logfile, 203) uii_uee
              write (s%logfile, 205) uxcdcc

              ! Evaluate total energy
              etot = ebs + uii_uee + uxcdcc
              write (s%logfile, 207) etot

              ! Total energy per atom
              etot_per_atom = etot/s%natoms
              write (s%logfile, 208) etot_per_atom

              ! Cohesive Energy
              atomic_energy = 0.0d0
              do iatom = 1, s%natoms
                in1 = s%atom(iatom)%imass
                atomic_energy = atomic_energy + species(in1)%atomicE
              end do
              write (s%logfile, 209) atomic_energy
              write (s%logfile, 210) etot - atomic_energy
              write (s%logfile, *)
              write (s%logfile, 211) (etot - atomic_energy)/s%natoms
              write (s%logfile, *) ' ----------------------------------------------------- '

! End scf loop
              if (sigma .gt. 0.0d0) then
                iscf_iteration = iscf_iteration + 1
              else
                exit
              end if
            end do

! ===========================================================================
! ---------------------------------------------------------------------------
!                             F O R C E S
! ---------------------------------------------------------------------------
! ===========================================================================
            call initialize_forces (s)
            call densityPP_matrix (s)

! Assemble the derivative blocks needed for forces - two-center
            write (s%logfile, *)
            write (s%logfile, *) ' Calling two-center non-charge dependent Dassemblers. '
            call Dassemble_T (s)
            call Dassemble_svnl (s)
            call Dassemble_vnl_2c (s)

! Assemble the derivative blocks needed for forces - three-center
            write (s%logfile,*) ' Calling three-center non-charge dependent Dassemblers. '
!           call Dassemble_vnl_3c (s)

            if (sigma .lt. scf_tolerance_set .or.                            &
     &          iscf_iteration .le. max_scf_iterations_set) then
! Assemble the derivative blocks needed for forces - two-center
              write (s%logfile, *) ' Calling two-center charge dependent Dassemblers. '
              call Dassemble_vna_2c (s)
!             call Dassemble_ewaldsr (s)
!             call Dassemble_ewaldlr (s)

! Assemble the derivative blocks needed for forces - three-center
              write (s%logfile, *) ' Calling three-center charge dependent Dassemblers. '
!             call Dassemble_vna_3c (s)
!             call Dassemble_vxc (s)
            end if

! Add in short-range (double-counting) contributions
            call Dassemble_uee (s)
!           call Dassemble_uxc (s)

            ! Build forces
            call build_forces (s)
            if (iwriteout_forces .eq. 1) call writeout_forces (s)

            write (s%logfile,*)
            write (s%logfile,*) ' Total Forces: '
            write (s%logfile,101)
            do iatom = 1, s%natoms
              write (s%logfile,302) iatom, s%forces(iatom)%ftot
            end do
            write (s%logfile,101)

            call md (s, itime_step)
            time = time + dt

! ===========================================================================
! ---------------------------------------------------------------------------
!                 G R A N D   T O T A L   E N E R G I E S
! ---------------------------------------------------------------------------
! ===========================================================================
            write (s%logfile,*)
            write (s%logfile,*) ' ------ T H E  G R A N D  T O T A L  E N E R G Y ------ '
            write (s%logfile,400) T_average
            write (s%logfile,*)
            write (s%logfile,401) tkinetic
            getot = etot + tkinetic
            write (s%logfile,402) getot
            getot_per_atom = getot/s%natoms
            write (s%logfile,403) etot_per_atom
            write (s%logfile,404) getot_per_atom
            write (s%logfile, *) ' ----------------------------------------------------- '
            if (itime_step .eq. nstepi) getot_initial = getot_per_atom

! Check energy conservation
            deltaE = 1000.d0*(getot_per_atom - getot_initial)
            write (s%logfile,405) deltaE

! Output the coordinates to a .xyz file
            slogfile = s%basisfile(:len(trim(s%basisfile))-4)
            slogfile = trim(slogfile)//'.xyz'
            open (unit = s%inpfile, file = slogfile, status = 'unknown', position = 'append')
            write (s%inpfile, *) s%natoms
            write (s%inpfile, *) etot, T_instantaneous
            do iatom = 1, s%natoms
              in1 = s%atom(iatom)%imass
              write (s%inpfile,*) species(in1)%symbol, s%atom(iatom)%ratom
            end do
            close (unit = s%inpfile)

! ===========================================================================
! ---------------------------------------------------------------------------
!              E L E C T R O N I C    T R A N S I T I O N S
! ---------------------------------------------------------------------------
! ===========================================================================
!           call calculate_nac (s)
!           call evolve_ks_states (s, itime_step)

! Destroy Hamiltonian matrix elements storage
            call destroy_assemble_2c (s)
            call destroy_assemble_PP_2c (s)
            call destroy_assemble_vxc_McWEDA (s)

            call destroy_kspace (s)
            call destroy_denmat (s)

            ! destroy neighbors last
            call destroy_neighbors (s)
            call destroy_neighbors_PP (s)
          end do ! end molecular dynamics loop

          call destroy_charges (s)
          deallocate (s%xl)

          if (iseparate .eq. 1) close (s%logfile)
        end do ! end loop over all structures
        close (1)

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

        write (ilogfile,*) ' LIGHTNING RUNTIME : ', time_end-time_begin, '[sec] '
        write (*,*) ' LIGHTNING RUNTIME : ', time_end-time_begin, '[sec] '
        close (ilogfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Working on structure - ', a25)
101     format (75 ('='))

202     format (2x, '           ebs = ', f15.6)
203     format (2x, '     uii - uee = ', f15.6)
205     format (2x, '        uxcdcc = ', f15.6)
207     format (2x, '          ETOT = ', f15.6)
208     format (2x, '     Etot/atom = ', f15.6)
209     format (2x, ' Atomic Energy = ', f15.6)
210     format (2x, '     CohesiveE = ', f15.6)
211     format (2x, ' Cohesive Energy per atom  = ', f15.6)

302     format (2x, ' iatom = ', i4, ' f = ', 3e14.6)

400     format (2x, '                       average temperature = ', f15.6)
401     format (2x, '                    nuclear kinetic energy = ', f15.6)
402     format (2x, ' Grand Total = nuclear kinetic + potential = ', f15.6)
403     format (2x, '                     total energy per atom = ', f15.6)
404     format (2x, '               grand total energy per atom = ', f15.6)
405     format (2x, '                        deltaE/atom  (meV) = ', f15.6)

! End Program
! ===========================================================================
        stop
        end program fireball_MDET
