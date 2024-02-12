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

SUBROUTINE integrate( Model, Solver, dt, TransientSimulation )
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
   REAL(KIND=dp), ALLOCATABLE ::NodalBVar(:), STIFF(:,:), LOAD(:), FORCE(:)

   TYPE(Nodes_t) :: ElementNodes
   TYPE(Element_t), POINTER :: Element, BoundaryElement
   INTEGER :: i, t, n, m, istat, nd

   TYPE(Variable_t), POINTER :: bVar, PointerToVariable ! pointer to the Balance Rate
   REAL(KIND=dp), POINTER :: bVarVals(:), VariableValues(:)
   REAL(KIND=dp) :: Norm
   INTEGER, POINTER :: Permutation(:), bVarPerm(:), NodeIndexes(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

   LOGICAL :: AllocationsDone
   
   SAVE :: NodalBVar, AllocationsDone, FORCE, LOAD, STIFF
   
   WRITE(SolverName, '(A)') 'Integrate'
   CALL INFO(SolverName,"Start",Level=1)
   
   IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalBVar)
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), NodalBVar(N), &
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
   END IF
   
   ! pick up solver and forcing variable
   PointerToVariable => Solver % Variable
   Permutation  => PointerToVariable % Perm
   VariableValues => PointerToVariable % Values
   bVar => VariableGet( Solver % Mesh % Variables, 'ToBeIntegrated', UnfoundFatal=.TRUE.)
   bVarVals => bVar % Values
   bVarPerm => bVar % Perm
   
   CALL DefaultInitialize()
   ! bulk assembly 
   Solver % variable % values = 0
   DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()

     NodeIndexes => Element % NodeIndexes

     DO i=1,n
       NodalBVar(i) = bVarVals(bVarPerm(NodeIndexes(i))) !-1.0_dp*MIN(bVarVals(bVarPerm(NodeIndexes(i))), 0.0_dp)
     END DO
   
     CALL LocalMatrix (  STIFF, FORCE, Element, n , NodalBVar )
     CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO
   
   ! Neumann conditions 
   DO t=1,Solver % Mesh % NUmberOfBoundaryElements
     BoundaryElement => GetBoundaryElement(t)
     NodeIndexes => BoundaryElement % NodeIndexes
     IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE
     n = GetElementNOFNodes()

     NodalBVar(1:n) = bVarVals(bVarPerm(NodeIndexes(1:n)))

     STIFF = 0.0D00
     FORCE = 0.0D00
!    CALL LocalMatrixBC(  STIFF, FORCE, BoundaryElement, n, NodalBVar)
!    CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO
   CALL DefaultFinishAssembly()
   CALL DefaultDirichletBCs()
   Norm = DefaultSolve()
   
 CONTAINS
   !------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n, var)
     REAL(KIND=dp) :: STIFF(:,:), FORCE(:) , var(:)
     INTEGER :: n, nd
     TYPE(Element_t), POINTER :: Element
     !------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n+1), dBasisdx(n+1,3), ddBasisddx(n,3,3), detJ,grad, IPvalue, s
     LOGICAL :: Stat
     INTEGER :: t, p,q ,dim
     TYPE(GaussIntegrationPoints_t) :: IP

     TYPE(Nodes_t) :: Nodes
     SAVE Nodes
     !------------------------------------------------------------------------------
     CALL GetElementNodes( Nodes )


     STIFF = 0.0d0
     FORCE = 0.0d0
     !PRINT *, var(1:n), grad
     
     IP = GaussPoints( Element )
     s = 0._dp
     DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )      

       grad  = SUM( var(1:n) * dBasisdx(1:n,1) )
       IPvalue = SUM( var(1:n) * Basis(1:n) )

       !PRINT *, "grad=", grad, "IPvalue", IPvalue
BLOCK
       real(kind=dp) :: h, coeff

       h = ElementDiameter(Element,nodes)
       coeff =h/2
       
       !FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ  * Basis(1:n)
       FORCE(1:n) = FORCE(1:n) + IPvalue * IP % s(t) * DetJ  * (Basis(1:n))

       s = s + IP % s(t) * detJ * ipvalue


       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % s(t) * detJ * coeff * dBasisdx(q,1)*dBasisdx(p,1)
           STIFF(p,q) = STIFF(p,q) + IP % s(t) * detJ * Basis(p)*dBasisdx(q,1) 
         END DO
       END DO
END BLOCK
     END DO

#if 0
BLOCK
     INTEGER :: n1,n2
     REAL(KIND=dP) :: x1,x2

     n1 = nodeindexes(1)
     n2 = nodeindexes(2)
     x1 = nodes % x(1)
     x2 = nodes % x(2)

     IF(x1>x2) THEN
       i = n1; n1 = n2; n2 = i
     END IF
    
     Solver  % Variable % Values(Solver % variable % perm(n2)) =  &
       Solver  % Variable % Values(Solver % Variable % perm(n1))  + s
END BLOCK
#endif

     !------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix
   !------------------------------------------------------------------------------


   
 END SUBROUTINE integrate
