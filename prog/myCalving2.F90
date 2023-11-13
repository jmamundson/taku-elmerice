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
! *  Solver to force calving events in 2D
! *
! *  Relies on TwoMeshes.F90 for mesh migration and interpolation 
! *  following calving, and FrontDisplacement.F90 for mesh update computation
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
! *  Original Date: 3.5.2014
! *
! ****************************************************************************/
SUBROUTINE Find_Calving (Model, Solver, dt, TransientSimulation )

   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE
!------------------------------------------------------------------------------
   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
   Type(Nodes_t) :: CornerElementNodes, CurrentElementNodes, &
        TargetNodes, CalvedNodes
   Type(Nodes_t), POINTER :: Nodes0
   Type(Nodes_t), TARGET :: EvalNodes
   TYPE(Matrix_t), POINTER :: Vel1Matrix !For finding neighbours
   TYPE(ValueList_t), POINTER :: Material, BC, SolverParams
   TYPE(Element_t), POINTER :: CurrentElement, CornerElement, Element
   TYPE(Mesh_t), POINTER :: Mesh, EvalMesh, Mesh0 => Null()
   TYPE(Variable_t), POINTER :: Vel1Sol, TimeVar, Var
   REAL (KIND=dp) :: t, CalvingSize, xcoord, ycoord, NodeDepth, Length, &
        NodeLength, Normal(3), CornerNormal(3), BedSecond, BedSecondDiff, &
        beddiff, BedToler, dx, dy, LocalDist, LocalDistNode, PropDistNode, &
        normalcond, work(8), &
#ifdef USE_ISO_C_BINDINGS
        rt0, rt
#else
        rt0, rt, RealTime
#endif

   REAL (KIND=DP), POINTER :: Calving1Values(:), Calving2Values(:), &
        FrontValues(:), WorkReal(:)
   REAL (KIND=dp), ALLOCATABLE :: CumDist(:), PropCumDist(:), &
        TargetCumDist(:), TargetPropDist(:), CalvingSizeNodes(:)
   INTEGER :: DIM, i, j, n, NoNodes, MaxN, FrontNodes, BotNodes, TopNodes, &
        BotCornerIndex, BotSecondIndex, county, GoToNode, PrevNode, &
        NextNode, NoNeighbours, MaxNeighbours, DOFs
   INTEGER, POINTER :: Vel1Perm(:), Vel1InvPerm(:), OrderPerm(:), &
        FrontPerm(:)=>NULL(), InvFrontPerm(:), TopPerm(:)=>NULL(), &
        BotPerm(:)=>NULL(), Permutation(:), NodeNeighbours(:,:), &
        NumNeighbours(:)
   INTEGER, ALLOCATABLE :: ThisNodeNeighbours(:)
   LOGICAL :: FirstTime = .TRUE., CalvingOccurs, RemeshOccurs, &
        CornerCalving, KeepLooking, Found, Debug = .FALSE., &
        BasalFS, CornerBadBed, CornerBadSlope
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FrontMaskName, BotMaskName, &
        TopMaskName, BasalFSVarName

   SAVE :: FirstTime, DIM, Permutation, Calving1Values, Calving2Values, &
        NoNodes, Mesh, FrontMaskName, FrontPerm, FrontValues, FrontNodes, &
        BotPerm, BotNodes, TopPerm, TopNodes, EvalMesh, Material, EvalNodes, &
        NodeNeighbours, NumNeighbours, Mesh0, Nodes0, BasalFS, BasalFSVarName, &
        Vel1InvPerm,NodeNeighbours,NumNeighbours

   ! Only solver parameter should be the calving event size   
   SolverParams => GetSolverParams()
   
   rt0 = RealTime() !time it!

   SolverName = "Calving"

   Debug = .FALSE.

   IF(ParEnv % Pes > 1) CALL Fatal(SolverName, "Solver doesnt work in parallel!")

   Timevar => VariableGet( Model % Variables,'Time', UnfoundFatal=.TRUE.)
   t = TimeVar % Values(1)
   dt = Model % Solver % dt
   
   RemeshOccurs = .FALSE. !For undercutting

   !*****************************Get Variables*******************************

   !This is the main solver variable, which only has a value on the calving
   !face of the glacier, and corresponds to the magnitude of retreat.
   IF(.NOT. ASSOCIATED(Solver % Variable)) CALL Fatal('Calving', 'Variable not associated!')

   IF(FirstTime) THEN

      DIM = CoordinateSystemDimension()
      IF(DIM /= 2) CALL Fatal('Calving','Solver only works in 2D!')
      MaxN = Model % Mesh % MaxElementNodes
      
      Mesh => Solver % Mesh
      NoNodes = SIZE( Mesh % Nodes % x )

      FrontMaskName = 'Calving Front Mask'
      TopMaskName = 'Top Surface Mask'
      BotMaskName = 'Bottom Surface Mask'
      ALLOCATE( FrontPerm(NoNodes), TopPerm(NoNodes), BotPerm(NoNodes))

      CALL MakePermUsingMask( Model,Solver,Mesh,FrontMaskName, &
           .FALSE., FrontPerm, FrontNodes )
      CALL MakePermUsingMask( Model,Solver,Mesh,TopMaskName, &
           .FALSE., TopPerm, TopNodes )
      CALL MakePermUsingMask( Model,Solver,Mesh,BotMaskName, &
           .FALSE., BotPerm, BotNodes )

      !Holds the variable values
      ALLOCATE( FrontValues(FrontNodes * 2) )

      Mesh0 => AllocateMesh()
      Mesh0 = Mesh
      Mesh0 % Name = TRIM(Mesh % Name)//'_initial'
      CALL Info('Calving','Created initial reference mesh to remap front, maintaining quality')
      ALLOCATE( Nodes0 )
      ALLOCATE( Nodes0 % x(NoNodes), Nodes0 % y(NoNodes), Nodes0 % z(NoNodes) )
      Nodes0 % x = Mesh % Nodes % x
      Nodes0 % y = Mesh % Nodes % y
      Nodes0 % z = Mesh % Nodes % z
      Mesh0 % Nodes => Nodes0

      Permutation => FrontPerm
      Calving1Values => FrontValues(1::2)
      Calving2Values => FrontValues(2::2)

      !Initialize
      Calving1Values = 0.0_dp
      Calving2Values = 0.0_dp

   ENDIF

   Solver % Variable % Values => FrontValues
   Solver % Variable % Perm => FrontPerm

   Var => VariableGet(Solver % Mesh % Variables, ComponentName(Solver % Variable % Name, 1), .TRUE.)
   Var % Values => Calving1Values
   Var % Perm => Permutation
   Var => VariableGet(Solver % Mesh % Variables, ComponentName(Solver % Variable % Name, 2), .TRUE.)
   Var % Values => Calving2Values
   Var % Perm => Permutation

   IF(FirstTime .OR. Solver % Mesh % Changed) THEN
      FirstTime = .FALSE.

      !STRATEGY: Finding neighbours on the fly works fine UNLESS you are in a recursive subroutine
      !Then it messes up, because at each level, ThisNodeNeighbours is deallocate and reallocated,
      !meaning that when you jump back up, info is already overwritten. SO:
      !Keep the structure as it was with CalvingNeighbours, cycle as below to fill it, and ALSO
      !create an array to hold the number of neighbours for each node

      !Get the Matrix of the N-S Solver
      Vel1Sol => VariableGet( Solver % Mesh % Variables, 'Velocity 1', UnfoundFatal=.TRUE.)
      Vel1Perm => Vel1Sol % Perm
      Vel1Matrix => Vel1Sol % Solver % Matrix

      !Vel * DIM + Pressure...
      DOFs = DIM + 1
      MaxNeighbours = DIM * 10  !totally arbitrary...

      !Create inverse perm to lookup Matrix later
      ALLOCATE(Vel1InvPerm(MAXVAL(Vel1Perm)*DOFs)) !TODO DEALLOCATE
      !2D array to hold each nodes neighbours
      ALLOCATE(NodeNeighbours(NoNodes,MaxNeighbours))
      !1D array to hold number of neighbours for each node
      ALLOCATE(NumNeighbours(NoNodes))
      NodeNeighbours = 0
      NumNeighbours = 0
      Vel1InvPerm = 0

      j = 0
      DO i=1,SIZE(Vel1Perm)
         IF(Vel1Perm(i) == 0) CYCLE
         j = j + 1
         Vel1InvPerm( (Vel1Perm(i)*DOFs-2) : (Vel1Perm(i)*DOFs) ) = j !The 2 here is suspect...
      END DO

      DO i = 1,NoNodes
         CALL FindNodeNeighbours(i) !Updates the allocatable array 'ThisNodeNeighbours'
         NumNeighbours(i) = SIZE(ThisNodeNeighbours)
         NodeNeighbours(i,1:NumNeighbours(i)) = ThisNodeNeighbours
      END DO
   END IF

   !MaxN = Model % Mesh % MaxElementNodes

   ! Determine the glacier length
   Length = 0.0
   DO i = 1, NoNodes
     IF(Permutation(i) == 0) CYCLE
     NodeLength = Model % Nodes % x(i)
     IF(NodeLength > Length) Length = NodeLength
   END DO
   PRINT *, '**** Glacier Length [m] = ',Length
   
   ! Get Calving Size from solver parameters
   CalvingSize = GetCReal(SolverParams,'Calving Size',Found )
   PRINT *, '**** Calving Event Size [m] = ', CalvingSize

   ! Determine horizontal component of calving, so that glacier
   ! terminus is vertical post-calving
   DO j = 1, NoNodes
      IF(Permutation(j) == 0) CYCLE
      Calving1Values(Permutation(j)) = Length - CalvingSize - &
            Mesh % Nodes % x(j)
      IF(Calving1Values(Permutation(j)) .GE. 0.0_dp) &
           Calving1Values(Permutation(j)) = 0.0_dp
   END DO

   IF (CalvingSize .GT. 0.0_dp) CalvingOccurs = .TRUE.

   
   !---------------------------------------
   !   CALVING DONE
   !---------------------------------------

   !At this point, the 'calving' solution is done.  However...
   !Here we solve a problem to do with undercutting:
   !Progressive undercutting by melting of the calving front
   !can lead to a situation where 'Front' nodes start to look like
   !basal nodes.  However, they are 'officially' front nodes and so
   !don't have a friction law, grounding dynamics OR (most importantly
   !I think), a bed constraint.  Thus, it is necessary to check for this
   !occurring, and shift the bed nodes appropriately.
   !
   !Strategy:
   !-Identify the corner node by BotPerm and FrontPerm
   !
   !-Use BCelement connections to find the second to bottom
   !    node on the calving front
   !
   !-Check a condition: either
   !    second node is 'near' bed OR
   !    BCelement slope is below some critical level
   !
   !-If condition is met (i.e. we need to take action)
   !    Calving1Values @CornerNode = (X@2nd - X@corner)
   !    CalvingOccurs = .TRUE.  <-- will this have any unforeseen consequences?
   !
   !NOTE: This works in tandem with a section of TwoMeshes.f90 which does the actual
   !deformation

   !Get the node index of the bottom corner
   !NOTE: this could be 'FirstTime'd if it was also 'SAVE'd
   DO i=1,NoNodes
      IF(BotPerm(i) > 0 .AND. FrontPerm(i) > 0) THEN
         BotCornerIndex = i
      END IF
   END DO

   !Loop boundary elements, we're looking for the BCelement
   !containing BotCornerIndex and ANOTHER FrontPerm node...
   DO i=Mesh % NumberOfBulkElements+1,&
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      CurrentElement => Mesh % Elements(i)
      IF(.NOT.(ANY(CurrentElement % nodeindexes == BotCornerIndex))) CYCLE
      IF(ALL(FrontPerm(CurrentElement % nodeindexes) .GT. 0)) THEN
         !We have a winner
         CornerElement => Mesh % Elements(i)
         DO j=1,2
            IF(CurrentElement % NodeIndexes(j) .NE. BotCornerIndex) THEN
               BotSecondIndex = CurrentElement % NodeIndexes(j)
            END IF
         END DO
      END IF
   END DO

   !Check corner node isn't already calving
   IF(Calving1Values(Permutation(BotCornerIndex)) .LT. 0.0_dp) THEN
      CornerCalving = .TRUE.
   ELSE
      CornerCalving = .FALSE.
   END IF

   !Get normal vector:
   CALL GetElementNodes(CornerElementNodes, CornerElement)
   CornerNormal = NormalVector(CornerElement, CornerElementNodes)

   IF(Debug) PRINT *, 'Debug Calving, corner normal is: ' , &
        CornerNormal(1), CornerNormal(2), CornerNormal(3)

   IF(BasalFS) THEN
     BedSecond = ListGetRealAtNode( Material, "Min "//BasalFSVarName, &
          BotSecondIndex, UnfoundFatal=.TRUE. )

     IF(Debug) PRINT *, 'Debug Calving, second node bed is: ',&
          BedSecond,' and y coord is: ', Model % Nodes % y(BotSecondIndex)

     PRINT *, 'Debug Calving, second node bed is: ',&
          BedSecond,' and y coord is: ', Model % Nodes % y(BotSecondIndex)

     BedSecondDiff = Model % Nodes % y(BotSecondIndex) - BedSecond
   END IF

   !TODO - unhardcode these
   BedToler = 2.0_dp
   normalcond = 0.95_dp

   CornerBadBed = BasalFS .AND. (BedSecondDiff < BedToler)
   CornerBadSlope = ABS(CornerNormal(2)) > normalcond

   !If the slope normal is above threshold, or the second node is too close to the bed,
   !move the corner node to the second, via 'calving'
   IF((CornerBadSlope .OR. CornerBadBed) .AND. (.NOT. CornerCalving)) THEN

      IF(Debug) PRINT *,'Debug Calving, migrating mesh'
      county = 1
      GoToNode = BotSecondIndex
      PrevNode = BotCornerIndex
      KeepLooking = .TRUE.
      DO WHILE (KeepLooking)
         !Check if we should shift more than one node forward...
         IF(Debug) PRINT *, 'Debug calving: looking!'
         KeepLooking = .FALSE.
         DO i=Mesh % NumberOfBulkElements+1,&
              Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
            CurrentElement => Mesh % Elements(i)
            IF(.NOT.(ALL(FrontPerm(CurrentElement % nodeindexes) .GT. 0))) CYCLE
            !This point reached by Front BC elements, stupid way of doing it, but whatever
            IF(.NOT.(ANY(CurrentElement % nodeindexes == GoToNode))) CYCLE
            IF(ANY(CurrentElement % nodeindexes == PrevNode)) CYCLE
            !We only get here if element is the next one along from previous
            CALL GetElementNodes(CurrentElementNodes, Currentelement)
            Normal = NormalVector(CurrentElement, CurrentElementNodes)
            DO j=1,2
              IF(CurrentElement % NodeIndexes(j) .NE. GoToNode) &
                   NextNode = CurrentElement % NodeIndexes(j)
            END DO

            IF(BasalFS) THEN
              beddiff = Model % Nodes % y(NextNode) - ListGetRealAtNode( Material, &
                   "Min "//BasalFSVarName, NextNode, UnfoundFatal=.TRUE. )
            END IF

            CornerBadBed = BasalFS .AND. (beddiff < BedToler)
            CornerBadSlope = ABS(CornerNormal(2)) > normalcond

            IF(CornerBadBed .OR. CornerBadSlope) THEN
               PrevNode = GoToNode
               GoToNode = NextNode
               county = county + 1
               IF(Debug) PRINT *, 'Debug calving: Found another shift'
               KeepLooking = .TRUE.
               EXIT
            END IF
         END DO
      END DO

      RemeshOccurs = .TRUE.
   ELSE
      county = 0
   END IF


   IF(CalvingOccurs .OR. RemeshOccurs) THEN

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Find the FrontDisplacement (= Calving 2 <-sif) for each frontal node
      !resulting from the shift in the corner node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Use FrontPerm to construct ordered node list
      !I think MakeMaskUsingPerm already ordered the nodes, from the comments:
      !! The bandwidth optimization for lines results to perfectly ordered
      !! permutations. If there is only one line the 1st node should be the
      !! lower left corner.
      ALLOCATE(InvFrontPerm(FrontNodes),&
           CumDist(FrontNodes),&
           PropCumDist(FrontNodes),&
           TargetCumDist(FrontNodes),&
           TargetPropDist(FrontNodes),&
           TargetNodes % x(FrontNodes),&
           TargetNodes % y(FrontNodes),&
           TargetNodes % z(FrontNodes),&
           CalvedNodes % x(FrontNodes))

      !InvFrontPerm(FrontNodes) points to nodeindexes in order they appear...
      !InvFrontPerm(1) points to 'lower left' corner, according to MakePermUsingMask
      DO i=1,NoNodes
         IF(FrontPerm(i) .GT. 0) THEN
            InvFrontPerm(FrontPerm(i)) = i
         END IF
      END DO

      ALLOCATE(OrderPerm(FrontNodes), WorkReal(FrontNodes))
      OrderPerm = [(i,i=1,FrontNodes)]
      DO i=1,FrontNodes
         WorkReal(i) = Mesh0 % Nodes % y(InvFrontPerm(i))
      END DO
      CALL SortD( FrontNodes, WorkReal, OrderPerm )
      DEALLOCATE(WorkReal)

      DO i=1,FrontNodes
         j = InvFrontPerm(OrderPerm(i))
         CalvedNodes % x(i) = Mesh % Nodes % x(j) + Calving1Values(Permutation(j))
         IF(Debug) PRINT *,'Debug Calving, CalvedNodes node: ',j,' is ',&
              CalvedNodes % x(i),' init coord: ',Mesh % Nodes % x(j)
      END DO

      !cycle through in order
      !First get target distribution from Mesh0
      TargetCumDist(1) = 0.0_dp
      DO i=2,FrontNodes
        dx = Mesh0 % Nodes % x(InvFrontPerm(OrderPerm(i))) - Mesh0 % Nodes % x(InvFrontPerm(OrderPerm(i-1)))
        dy = Mesh0 % Nodes % y(InvFrontPerm(OrderPerm(i))) - Mesh0 % Nodes % y(InvFrontPerm(OrderPerm(i-1)))

        TargetCumDist(i) = TargetCumDist(i-1) + (((dx**2) + (dy**2)) ** 0.5_dp)
      END DO
      TargetPropDist = TargetCumDist / MAXVAL(TargetCumDist)

      !Now find the length segments of our current line
      !If RemeshOccurs (because of bad corner node), county dictates
      !the offset from the previous bottom node of the new calving front
      CumDist(1:county+1) = 0.0_dp
      DO i=county+2,FrontNodes
        !sum coord magnitude from base upwards to give front 'length'
        !keep cumulative total
        !allocate proporitional y (and x) distances (i.e. out of 1)
        dx = CalvedNodes % x(i) - CalvedNodes % x(i-1)
        dy = Mesh % Nodes % y(InvFrontPerm(OrderPerm(i))) - Mesh % Nodes % y(InvFrontPerm(OrderPerm(i-1)))

        CumDist(i) = CumDist(i-1) + (((dx**2) + (dy**2)) ** 0.5_dp)
        !Remember first one is corner node...
        IF(Debug) PRINT *, 'Debug Calving: CumDist at node: ',&
             InvFrontPerm(OrderPerm(i)),' is ',CumDist(i)
        IF(Debug) PRINT *, 'Debug Calving: TargetDist at node: ',&
             InvFrontPerm(OrderPerm(i)),' is ',TargetCumDist(i)
      END DO
      PropCumDist = CumDist / MAXVAL(CumDist)

      !Loop each front node
      TargetNodes % x(1) = CalvedNodes % x(county+1)
      TargetNodes % y(1) = Mesh % Nodes % y(InvFrontPerm(OrderPerm(county+1)))
      TargetNodes % x(FrontNodes) = CalvedNodes % x(FrontNodes)
      TargetNodes % y(FrontNodes) = Mesh % Nodes % y(InvFrontPerm(OrderPerm(FrontNodes)))

      DO i=2,FrontNodes-1
        !and find nearest two nodes to interpolate
        DO j=county+2,FrontNodes
          IF(PropCumDist(j) .GT. TargetPropDist(i)) THEN
            !lin interp between j and j-1
            LocalDist = PropCumDist(j) - PropCumDist(j-1)
            LocalDistNode = TargetPropDist(i) - PropCumDist(j-1)

            PropDistNode = LocalDistNode / LocalDist
            IF(Debug) PRINT *, 'Debug Calving: PropDist at node: ',&
                 InvFrontPerm(OrderPerm(i)),' is ',PropDistNode

            TargetNodes % x(i) = ((1 - PropDistNode) * CalvedNodes % x(j-1))  + &
                 (PropDistNode * CalvedNodes % x(j))
            TargetNodes % y(i) = ((1 - PropDistNode) * Mesh % Nodes % y(InvFrontPerm(OrderPerm(j-1))))  + &
                 (PropDistNode * Mesh % Nodes % y(InvFrontPerm(OrderPerm(j))))
            EXIT
          END IF
        END DO
      END DO

      !At this point, we have obtained, for each FrontNode, a TargetNode % x and y
      !Thus, it simply remains to compute the two components of the displacement

      !Calving 1 = Diff X  (New % x - Old % x)
      !Calving 2 = Diff Y  (New % y - Old % y)

      DO i=1,FrontNodes
         Calving1Values(Permutation(InvFrontPerm(OrderPerm(i)))) = TargetNodes % x(i) &
              - Mesh % Nodes % x(InvFrontPerm(OrderPerm(i)))

         Calving2Values(Permutation(InvFrontPerm(OrderPerm(i)))) = TargetNodes % y(i) &
              - Mesh % Nodes % y(InvFrontPerm(OrderPerm(i)))

         IF(Debug) THEN
            PRINT *,'Debug Calving: Node: ',InvFrontPerm(OrderPerm(i)),' pos x: ',&
                 Mesh % nodes % x(InvFrontPerm(OrderPerm(i))),&
                 ' pos y: ',Mesh % nodes % y(InvFrontPerm(OrderPerm(i)))
            PRINT *,'Moving to: x: ',TargetNodes % x(i),' y: ',TargetNodes % y(i)
            PRINT *,'Displacement 1: ',Calving1Values(Permutation(InvFrontPerm(OrderPerm(i)))),&
                 'Displacement 2: ',Calving2Values(Permutation(InvFrontPerm(OrderPerm(i))))
         END IF
       END DO
     END IF

   CALL ListAddLogical( Model % Simulation, 'CalvingOccurs', CalvingOccurs )
   CALL ListAddLogical( Model % Simulation, 'RemeshOccurs', RemeshOccurs )


 CONTAINS

   SUBROUTINE FindNodeNeighbours(NodeNumber)
     INTEGER :: NodeNumber, i, count

     NoNeighbours = Vel1Matrix % Rows((Vel1Perm(NodeNumber)*DOFs)+1) - Vel1Matrix % Rows(Vel1Perm(NodeNumber)*DOFs)
     IF(MOD(NoNeighbours, DOFs).NE. 0) CALL Fatal(SolverName,"This shouldn't have happened...")
     !Each neighbour appears once per DOF, and there's also the current node thus: (x/DOFS) - 1...
     NoNeighbours = (NoNeighbours / DOFs) - 1
     IF(NoNeighbours .GT. MaxNeighbours) CALL Fatal(SolverName,"Need more neighbour slots!")

     IF(ALLOCATED(ThisNodeNeighbours)) DEALLOCATE(ThisNodeNeighbours)
     ALLOCATE(ThisNodeNeighbours(NoNeighbours))
     ThisNodeNeighbours = 0

     count = 0
     DO i=Vel1Matrix % Rows(Vel1Perm(NodeNumber)*DOFs),&
          (Vel1Matrix % Rows((Vel1Perm(NodeNumber)*DOFs)+1)-1)
        IF(MOD(i,DOFs) .NE. 0) CYCLE !Stored DOF1, DOF2, DOF3, only need every 3rd
        IF(Vel1InvPerm(Vel1Matrix % Cols(i)) == NodeNumber) CYCLE !Not our own neighbour
        count = count + 1
        ThisNodeNeighbours(count) = &
             Vel1InvPerm(Vel1Matrix % Cols(i))
     END DO

   END SUBROUTINE FindNodeNeighbours

  
END SUBROUTINE Find_Calving
