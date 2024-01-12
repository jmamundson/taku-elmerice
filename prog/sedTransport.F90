!*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Solver to calculate rate of sediment excavation/deposition. Based on
! *  Brinkerhoff et al. (2017).
! *  
! ******************************************************************************
! *
! *  Authors: Jason Amundson
! *  Email:   jmamundson@alaska.edu
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 11.1.2024
! *
! ****************************************************************************/

SUBROUTINE calc_Qb( Model, Solver, dt, TransientSimulation )
  ! calculates the balance flux (per unit width) as a function of position
  
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   !USE quadpack ! I think this can do the numerical integration?
   
   IMPLICIT NONE ! don't treat variables that start with i, j, k, l, m, and n as integers
!------------------------------------------------------------------------------
   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
   
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
   REAL(KIND=dp), ALLOCATABLE :: Xcoord(:), balanceRate(:), balanceFlux(:)

   TYPE(Nodes_t) :: ElementNodes
   TYPE(Element_t) :: Element
   INTEGER :: i

   TYPE(Variable_t), POINTER :: bVar ! pointer to the Balance Rate

   
   SAVE :: Xcoord, balanceRate
   
   CALL GetElementNodes( ElementNodes, Element, Solver )

   Xcoord = ElementNodes % x(i) ! x-coordinate of nodes (?)
   
   bVar => VariableGet( Solver % Mesh % Variables, 'Balance Rate', UnfoundFatal=.TRUE.)
   balanceRate = bVar % Values(i) ! balance rate at the nodes (?)


   !balanceFlux = integrate(Xcoord,balanceRate)
   

   

END SUBROUTINE calc_Qb
