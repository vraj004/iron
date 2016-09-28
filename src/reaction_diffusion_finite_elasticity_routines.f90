!> \file
!> \authors Vijay Rajagopal
!> \brief This module handles all routines pertaining to reaction diffusion coupled with finite elasticity.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): Vijay Rajagopal
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module handles all routines pertaining to reaction-diffusion coupled with finite elasticity.


MODULE REACTION_DIFFUSION_FINITE_ELASTICITY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS
  USE PROBLEM_CONSTANTS
  USE REACTION_DIFFUSION_EQUATION_ROUTINES
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TYPES

#include "macros.h"

  IMPLICIT NONE

  PUBLIC ReacDiffFiniteElastic_EquationsSetSetup
  PUBLIC ReacDiffFiniteElastic_EquationsSetSolutionMethodSet

  PUBLIC REAC_DIFF_FINITE_ELASTIC_PROBLEM_SETUP
  PUBLIC ReacDiffFiniteElastic_ProblemSpecificationSet
 
  PUBLIC ReacDiffFiniteElastic_FiniteElementCalculate

  PUBLIC REAC_DIFF_FINITE_ELASTIC_PRE_SOLVE
  PUBLIC REAC_DIFF_FINITE_ELASTIC_POST_SOLVE

  PUBLIC ReacDiffFiniteElastic_ControlLoopPreLoop
  PUBLIC ReacDiffFiniteElastic_ControlLoopPostLoop
  
  PUBLIC ReacDiffFiniteElastic_UpdateGeometricField

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a ReactionDiffusions finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE ReacDiffFiniteElastic_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("ReacDiffFiniteElastic_EquationsSetSolutionMethodSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a "// &
          & "ReactionDiffusion-finite elasticity type equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        !\todo what are problem constants doing here???
      CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a ReactionDiffusion finite elasticity equation type of a multi physics equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ReacDiffFiniteElastic_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_EquationsSetSolutionMethodSet",ERR,ERROR)
    EXITS("ReacDiffFiniteElastic_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE ReacDiffFiniteElastic_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the ReactionDiffusion finite elasticity equation.
  SUBROUTINE ReacDiffFiniteElastic_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("ReacDiffFiniteElastic_EquationsSetSetup",ERR,ERROR,*999)

    CALL FlagError("ReacDiffFiniteElastic_EquationsSetSetup is not implemented.",ERR,ERROR,*999)

    EXITS("ReacDiffFiniteElastic_EquationsSetSetup")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_EquationsSetSetup",ERR,ERROR)
    EXITS("ReacDiffFiniteElastic_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE ReacDiffFiniteElastic_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a ReactionDiffusions finite elasticity equation finite element equations set.
  SUBROUTINE ReacDiffFiniteElastic_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("ReacDiffFiniteElastic_FiniteElementCalculate",ERR,ERROR,*999)

    CALL FlagError("ReacDiffFiniteElastic_FiniteElementCalculate is not implemented.",ERR,ERROR,*999)

    EXITS("ReacDiffFiniteElastic_FiniteElementCalculate")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_FiniteElementCalculate",ERR,ERROR)
    EXITS("ReacDiffFiniteElastic_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE ReacDiffFiniteElastic_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a ReactionDiffusion finite elasticity problem type .
  SUBROUTINE ReacDiffFiniteElastic_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("ReacDiffFiniteElastic_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a ReactionDiffusion finite elasticity type of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_REACTION_DIFFUSION_FINITE_ELASTICITY_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("ReactionDiffusion finite elasticity problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("ReacDiffFiniteElastic_ProblemSpecificationSet")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_ProblemSpecificationSet",err,error)
    EXITS("ReacDiffFiniteElastic_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE ReacDiffFiniteElastic_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the ReactionDiffusion finite elasticity problem.
  SUBROUTINE REAC_DIFF_FINITE_ELASTIC_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: REACTION_DIFFUSION_SUB_LOOP,ELASTICITY_SUB_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS,REACTION_DIFFUSION_SOLVERS,ELASTICITY_SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("REAC_DIFF_FINITE_ELASTIC_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(REACTION_DIFFUSION_SUB_LOOP)
    NULLIFY(ELASTICITY_SUB_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(REACTION_DIFFUSION_SOLVERS)
    NULLIFY(ELASTICITY_SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(CELLML_EQUATIONS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a ReactionDiffusion-finite elasticity problem.", &
        & err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   Transient Strang-split reaction diffusion, simple finite elasticity  
      !--------------------------------------------------------------------
      CASE(PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a ReactionDiffusions finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up the control sub loop for reaction diffusion
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(REACTION_DIFFUSION_SUB_LOOP,'REACTION_DIFFUSION_TIME_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(REACTION_DIFFUSION_SUB_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(REACTION_DIFFUSION_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up the control sub loop for finite elasicity
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(ELASTICITY_SUB_LOOP,'ELASTICITY_LOAD_INCREMENT_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(ELASTICITY_SUB_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(ELASTICITY_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
          
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
            !Sub-loops are finished when parent is finished
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a ReactionDiffusions finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the reaction diffusion sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(REACTION_DIFFUSION_SUB_LOOP,REACTION_DIFFUSION_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(REACTION_DIFFUSION_SOLVERS,3,ERR,ERROR,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"First ODE Solver",ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !Set the third solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,3,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Second ODE Solver",ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)

            !Get the finite elasticity sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(ELASTICITY_SOLVERS,1,ERR,ERROR,*999)
            !Set the finite elasticity solver to be a nonlinear solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the reaction diffusion solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(REACTION_DIFFUSION_SUB_LOOP,REACTION_DIFFUSION_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(REACTION_DIFFUSION_SOLVERS,ERR,ERROR,*999)

            !Get the finite elasticity solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(ELASTICITY_SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a ReactionDiffusions finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)

            !Get the reaction diffusion sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(REACTION_DIFFUSION_SUB_LOOP,REACTION_DIFFUSION_SOLVERS,ERR,ERROR,*999)
            !Create the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)

            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Get the finite elasticity solver and create the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            
            !Get the reaction diffusion sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(REACTION_DIFFUSION_SUB_LOOP,REACTION_DIFFUSION_SOLVERS,ERR,ERROR,*999)
            !Get the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             

            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the creation of the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a ReactionDiffusions finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(REACTION_DIFFUSION_SUB_LOOP,REACTION_DIFFUSION_SOLVERS,ERR,ERROR,*999)
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)            
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Create the CellML equations for the second DAE solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)            
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,3,SOLVER,ERR,ERROR,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)

          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(REACTION_DIFFUSION_SUB_LOOP,REACTION_DIFFUSION_SOLVERS,ERR,ERROR,*999)
            !Get the CellML equations for the first DAE solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
            !Get the CellML equations for the second DAE solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(REACTION_DIFFUSION_SOLVERS,3,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)

          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a ReactionDiffusions finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a ReactionDiffusions finite elasticity equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a transient reaction diffusion quasistatic finite elasticity equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("REAC_DIFF_FINITE_ELASTIC_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("REAC_DIFF_FINITE_ELASTIC_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE REAC_DIFF_FINITE_ELASTIC_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the ReactionDiffusions finite elasticity problem pre-solve.
  SUBROUTINE REAC_DIFF_FINITE_ELASTIC_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("REAC_DIFF_FINITE_ELASTIC_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a ReactionDiffusion-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              CALL REACTION_DIFFUSION_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a ReactionDiffusions finite elasticity type of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a reaction diffusion finite elasticity type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REAC_DIFF_FINITE_ELASTIC_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("REAC_DIFF_FINITE_ELASTIC_PRE_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE REAC_DIFF_FINITE_ELASTIC_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the reaction-diffusions finite elasticity problem post solve.
  SUBROUTINE REAC_DIFF_FINITE_ELASTIC_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("REAC_DIFF_FINITE_ELASTIC_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a reaction-diffusion-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE)
            SELECT CASE(SOLVER%SOLVE_TYPE)
            CASE(SOLVER_DAE_TYPE)
              CALL REACTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_TYPE)
              CALL REACTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_NONLINEAR_TYPE)
              CALL FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Solver solve type "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a reaction diffusion finite elasticity type of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a reaction diffusion finite elasticity type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("REAC_DIFF_FINITE_ELASTIC_POST_SOLVE")
    RETURN
999 ERRORSEXITS("REAC_DIFF_FINITE_ELASTIC_POST_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE REAC_DIFF_FINITE_ELASTIC_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the reaction-diffusions finite elasticity problem pre-control loop.
  SUBROUTINE ReacDiffFiniteElastic_ControlLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("ReacDiffFiniteElastic_ControlLoopPreLoop",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a reaction-diffusion-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_REACTION_DIFFUSION_FINITE_ELASTICITY_TYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              !do nothing ???
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE)
                !CALL ReacDiffFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The third problem specification of "// &
                  & TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for reaction-diffusion-finite elasticity problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL FiniteElasticity_ControlTimeLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for reaction-diffusion finite elasticity problem type."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The second problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("ReacDiffFiniteElastic_ControlLoopPreLoop")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_ControlLoopPreLoop",ERR,ERROR)
    EXITS("ReacDiffFiniteElastic_ControlLoopPreLoop")
    RETURN 1
    
  END SUBROUTINE ReacDiffFiniteElastic_ControlLoopPreLoop

  !
  !================================================================================================================================
  !

  !>Sets up the reaction-diffusion finite elasticity problem post-control loop.
  SUBROUTINE ReacDiffFiniteElastic_ControlLoopPostLoop(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG) :: equations_set_idx
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION   
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,LOCAL_ERROR,METHOD
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ELASTICITY_SUB_LOOP,REACTION_DIFFUSION_SUB_LOOP

    ENTERS("ReacDiffFiniteElastic_ControlLoopPostLoop",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a reaction diffusion-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_REACTION_DIFFUSION_FINITE_ELASTICITY_TYPE)
              !the reaction diffusion time loop - output of the reaction diffusion fields
              CALL REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
                & " is not valid for a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            CALL ReacDiffFiniteElastic_UpdateGeometricField(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
          CASE DEFAULT
            !do nothing
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - output the finite elasticity fields 
        IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          !Export the dependent field for this time step
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            PROBLEM=>CONTROL_LOOP%PROBLEM
            IF(ASSOCIATED(PROBLEM)) THEN
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(ELASTICITY_SUB_LOOP)
              !Get the solver. The first solver of the second sub loop will contain the finite elasticity dependent field equation set
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              !Loop over the equations sets associated with the solver
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      NULLIFY(DEPENDENT_REGION)
                      CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                      FILENAME="MainTime_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                        & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                      METHOD="FORTRAN"
                      CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for equations set index "// &
                        & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                        & " in the solver mapping."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                ELSE
                  CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
              ENDIF
              IF((PROBLEM%SPECIFICATION(3)==PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE)) THEN
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(REACTION_DIFFUSION_SUB_LOOP)
                !Get the solver. The second solver of the first sub loop will contain the reaction diffusion equation set
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,REACTION_DIFFUSION_SUB_LOOP,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(REACTION_DIFFUSION_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                !Loop over the equations sets associated with the solver
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        NULLIFY(DEPENDENT_REGION)
                        CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                        FILENAME="MainTime_RD_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                          & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                        METHOD="FORTRAN"
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                        
                        WRITE(*,*) TIME_LOOP%ITERATION_NUMBER
                        
                      ELSE
                        LOCAL_ERROR="Equations set is not associated for equations set index "// &
                          & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                          & " in the solver mapping."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !equations_set_idx
                  ELSE
                    CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF !PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE
            ELSE
              CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("ReactinoDiffusionFiniteElastic_ControlLoopPostLoop")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_ControlLoopPostLoop",ERR,ERROR)
    EXITS("ReacDiffFiniteElastic_ControlLoopPostLoop")
    RETURN 1
    
  END SUBROUTINE ReacDiffFiniteElastic_ControlLoopPostLoop


  !
  !================================================================================================================================
  !


  !>Update the the reaction diffusion equation geometric field from the finite elasticity dependent field (deformed geometry)
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE ReacDiffFiniteElastic_UpdateGeometricField(CONTROL_LOOP,CALC_CLOSEST_GAUSS_POINT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    LOGICAL, INTENT(IN) :: CALC_CLOSEST_GAUSS_POINT !<If true then the closest finite elasticity Gauss point for each bioelectrics node is calculated
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_ELASTICITY, &
      & CONTROL_LOOP_REACTION_DIFFUSION
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_ELASTICITY,GEOMETRIC_FIELD_REACTION_DIFFUSION,GEOMETRIC_FIELD_ELASTICITY
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_REACTION_DIFFUSION,INDEPENDENT_FIELD_REACTION_DIFFUSION,DEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_DEP_M,FIELD_VAR_GEO_M,FIELD_VAR_IND_FE,FIELD_VAR_IND_M,FIELD_VAR_IND_M_2
    INTEGER(INTG) :: component_idx,element_idx,ne,start_elem,START_ELEMENT,start_element_idx
    INTEGER(INTG) :: DEPENDENT_FIELD_INTERPOLATION,GEOMETRIC_FIELD_INTERPOLATION
    INTEGER(INTG) :: node_idx,node_idx_2,GAUSS_POINT,gauss_idx,fibre_idx
    INTEGER(INTG) :: nodes_in_Xi_1,nodes_in_Xi_2,nodes_in_Xi_3,n3,n2,n1,dof_idx,dof_idx2,idx,my_element_idx
    REAL(DP) :: VALUE,VALUE_LEFT,DISTANCE,VELOCITY,VELOCITY_MAX,OLD_DIST
    REAL(DP) :: XI(3),PREVIOUS_NODE(3),DIST_INIT,SARCO_LENGTH_INIT,TIME_STEP,DIST
    REAL(DP), POINTER :: GAUSS_POSITIONS(:,:)

    ENTERS("ReacDiffFiniteElastic_UpdateGeometricField",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(CONTROL_LOOP_REACTION_DIFFUSION)
    NULLIFY(PROBLEM)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(INDEPENDENT_FIELD_ELASTICITY)
    NULLIFY(DEPENDENT_FIELD_REACTION_DIFFUSION)
    NULLIFY(INDEPENDENT_FIELD_REACTION_DIFFUSION)
    NULLIFY(DEPENDENT_FIELD_ELASTICITY)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(GEOMETRIC_FIELD_REACTION_DIFFUSION)
    NULLIFY(GEOMETRIC_FIELD_ELASTICITY)
    NULLIFY(ELEMENTS_TOPOLOGY)
    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(FIELD_VAR_DEP_M)
    NULLIFY(FIELD_VAR_GEO_M)
    NULLIFY(FIELD_VAR_IND_FE)
    NULLIFY(FIELD_VAR_IND_M)
    NULLIFY(FIELD_VAR_IND_M_2)
    NULLIFY(GAUSS_POSITIONS)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a reaction-diffusion-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_REACTION_DIFFUSION_FINITE_ELASTICITY_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(3))

            CASE(PROBLEM_STRANG_REACTION_DIFFUSION_ELASTICITY_SUBTYPE)

              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
              !get the REACTION_DIFFUSION sub loop, solvers, solver, and finally geometric and field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_REACTION_DIFFUSION,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_REACTION_DIFFUSION,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    GEOMETRIC_FIELD_REACTION_DIFFUSION=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_REACTION_DIFFUSION)) THEN
                      CALL FlagError("Geometric field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_MAPPING)
              NULLIFY(EQUATIONS_SET)
              NULLIFY(SOLVER_EQUATIONS)
              !get the finite elasticity sub loop, solvers, solver, and finally the dependent field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    DEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(DEPENDENT_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              DO component_idx=1,GEOMETRIC_FIELD_REACTION_DIFFUSION%VARIABLES(1)%NUMBER_OF_COMPONENTS
                !check for identical interpolation of the fields
                GEOMETRIC_FIELD_INTERPOLATION= &
                & GEOMETRIC_FIELD_REACTION_DIFFUSION%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
                DEPENDENT_FIELD_INTERPOLATION=DEPENDENT_FIELD_ELASTICITY%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
                IF(GEOMETRIC_FIELD_INTERPOLATION==DEPENDENT_FIELD_INTERPOLATION) THEN
                  !copy the dependent field components to the geometric field
                  CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,GEOMETRIC_FIELD_REACTION_DIFFUSION, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & component_idx,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="The interpolation type of component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR, &
                    & ERROR))//" of field number "// &
                    & TRIM(NUMBER_TO_VSTRING(GEOMETRIC_FIELD_REACTION_DIFFUSION%USER_NUMBER,"*",ERR, &
                    & ERROR))//" does not coincide with the interpolation type of field number " &
                    & //TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_ELASTICITY%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO
            CASE DEFAULT
              LOCAL_ERROR="The third problem specification of "// &
                & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(2),"*",ERR,ERROR))// &
                & " is not valid for a reaction diffusion finite elasticity of a multi physics problem."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The second problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(2),"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("ReacDiffFiniteElastic_UpdateGeometricField")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_UpdateGeometricField",ERR,ERROR)
    EXITS("ReacDiffFiniteElastic_UpdateGeometricField")
    RETURN 1
    
  END SUBROUTINE ReacDiffFiniteElastic_UpdateGeometricField

  !
  !================================================================================================================================
  !

  !>Interpolates the finite elasticity independent field from the biolectrics independent field.
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE ReacDiffFiniteElastic_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ELASTICITY,CONTROL_LOOP_REACTION_DIFFUSION
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_REACTION_DIFFUSION,INDEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING,NODES_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE_U,FIELD_VARIABLE_V,FIELD_VARIABLE_FE
    INTEGER(INTG) :: node_idx,element_idx,gauss_idx,ne
    INTEGER(INTG) :: nearestGP,inElement,dof_idx
    INTEGER(INTG) :: NUMBER_OF_GAUSS_POINTS
    REAL(DP) :: ACTIVE_STRESS
    REAL(DP) :: TITIN_STRESS_UNBOUND,TITIN_STRESS_BOUND,TITIN_STRESS_CROSS_FIBRE_UNBOUND,TITIN_STRESS_CROSS_FIBRE_BOUND,ACTIVATION
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_GAUSS_POINTS=64
    INTEGER(INTG) :: NUMBER_OF_NODES(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: ACTIVE_STRESS_VALUES(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: TITIN_STRESS_VALUES_UNBOUND(MAX_NUMBER_OF_GAUSS_POINTS),TITIN_STRESS_VALUES_BOUND(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP):: ACTIVATION_VALUES(MAX_NUMBER_OF_GAUSS_POINTS)

    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_REACTION_DIFFUSION)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(FIELD_VARIABLE_U)
    NULLIFY(FIELD_VARIABLE_V)
    NULLIFY(FIELD_VARIABLE_FE)
    
    ENTERS("ReacDiffFiniteElastic_IndependentFieldInterpolate",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        IF(.NOT.ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(problem%specification,1)<3) THEN
          CALL FlagError("Problem specification must have three entries for a reaction-diffusion-finite elasticity problem.", &
            & err,error,*999)
        END IF
        SELECT CASE(PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
          IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
            !--- REACTION_DIFFUSION ---
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_REACTION_DIFFUSION,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_REACTION_DIFFUSION,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_REACTION_DIFFUSION=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_REACTION_DIFFUSION)) CALL FlagError("Independent field is not associated.", &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF

            !--- FINITE ELASTICITY ---
            NULLIFY(SOLVERS)
            NULLIFY(SOLVER)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) CALL FlagError("Independent field is not associated.",ERR, &
                    & ERROR,*999)
                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF

            !--- NOW INTERPOLATE ---
            ELEMENTS_MAPPING=>INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
            NODES_MAPPING=>INDEPENDENT_FIELD_REACTION_DIFFUSION%DECOMPOSITION% &
              & DOMAIN(INDEPENDENT_FIELD_REACTION_DIFFUSION%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES

            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_REACTION_DIFFUSION,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_U,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_REACTION_DIFFUSION,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE_V,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_FE,ERR,ERROR,*999)

            !loop over the finite elasticity elements
            !first process the internal and boundary elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%BOUNDARY_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              
              NUMBER_OF_GAUSS_POINTS=INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY% &
                & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP &
                & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS

              IF(NUMBER_OF_GAUSS_POINTS>MAX_NUMBER_OF_GAUSS_POINTS) CALL FlagError( & 
                & "NUMBER_OF_GAUSS_POINTS is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)
              NUMBER_OF_NODES=0
              ACTIVE_STRESS_VALUES=0.0_DP
              TITIN_STRESS_VALUES_UNBOUND=0.0_DP
              TITIN_STRESS_VALUES_BOUND=0.0_DP
              TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND=0.0_DP
              TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND=0.0_DP
              ACTIVATION_VALUES=0.0_DP
              
              !loop over the reaction-diffusions nodes
              DO node_idx=1,NODES_MAPPING%NUMBER_OF_LOCAL
                dof_idx=FIELD_VARIABLE_V%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION, & 
                  & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & dof_idx,inElement,ERR,ERROR,*999) !component 5 of variable V contains inElem info (LOCAL NUMBERING!!!)

                !check if the reaction-diffusions node is located within the finite elasticity element
                IF(inElement==ne) THEN
                  !component 4 of variable V contains Nearest Gauss Point info
                  dof_idx=FIELD_VARIABLE_V%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                    & VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & dof_idx,nearestGP,ERR,ERROR,*999)
                  IF(nearestGP>MAX_NUMBER_OF_GAUSS_POINTS) CALL FlagError( &
                    & "Nearest Gauss Point is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)
                  !component 1 of variable U contains the active stress
                  dof_idx=FIELD_VARIABLE_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                    & VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)
                  
                  !count the number of reaction-diffusions nodes that are closest to each finite elasticity Gauss point
                  NUMBER_OF_NODES(nearestGP)=NUMBER_OF_NODES(nearestGP)+1
                  !add up the active stress value
                  ACTIVE_STRESS_VALUES(nearestGP)=ACTIVE_STRESS_VALUES(nearestGP)+ACTIVE_STRESS

                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                    !component 2 of variable U contains the titin stress unbound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_UNBOUND,ERR,ERROR,*999)
                    !component 3 of variable U contains the titin stress bound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_BOUND,ERR,ERROR,*999)
                    !component 4 of variable U contains the titin XF-stress (cross-fibre directions) unbound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_UNBOUND,ERR,ERROR,*999)
                    !component 5 of variable U contains the titin XF-stress (cross-fibre directions) bound
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_BOUND,ERR,ERROR,*999)
                    !component 6 of variable U contains the titin activation
                    dof_idx=FIELD_VARIABLE_U%COMPONENTS(6)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_REACTION_DIFFUSION,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,dof_idx,ACTIVATION,ERR,ERROR,*999)

                    TITIN_STRESS_VALUES_UNBOUND(nearestGP)=TITIN_STRESS_VALUES_UNBOUND(nearestGP)+TITIN_STRESS_UNBOUND
                    TITIN_STRESS_VALUES_BOUND(nearestGP)=TITIN_STRESS_VALUES_BOUND(nearestGP)+TITIN_STRESS_BOUND
                    TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(nearestGP)=TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(nearestGP) + &
                      & TITIN_STRESS_CROSS_FIBRE_UNBOUND
                    TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(nearestGP)=TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(nearestGP) + &
                      & TITIN_STRESS_CROSS_FIBRE_BOUND
                    ACTIVATION_VALUES(nearestGP)=ACTIVATION_VALUES(nearestGP)+ACTIVATION
                  ENDIF

                ENDIF
              ENDDO

              !loop over the finite elasticity Gauss points
              DO gauss_idx=1,NUMBER_OF_GAUSS_POINTS
                !make sure we don't divide by zero
                IF(NUMBER_OF_NODES(gauss_idx)<=0) THEN
                  ACTIVE_STRESS=0.0_DP
                  TITIN_STRESS_UNBOUND=0.0_DP
                  TITIN_STRESS_BOUND=0.0_DP
                  TITIN_STRESS_CROSS_FIBRE_UNBOUND=0.0_DP
                  TITIN_STRESS_CROSS_FIBRE_BOUND=0.0_DP
                  ACTIVATION=0.0_DP
                ELSE
                  ACTIVE_STRESS=ACTIVE_STRESS_VALUES(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_UNBOUND=TITIN_STRESS_VALUES_UNBOUND(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_BOUND=TITIN_STRESS_VALUES_BOUND(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_CROSS_FIBRE_UNBOUND=TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                  TITIN_STRESS_CROSS_FIBRE_BOUND=TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(gauss_idx)/ &
                    & NUMBER_OF_NODES(gauss_idx)
                  ACTIVATION=ACTIVATION_VALUES(gauss_idx)/NUMBER_OF_NODES(gauss_idx)
                ENDIF

                dof_idx=FIELD_VARIABLE_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,dof_idx,ACTIVE_STRESS,ERR,ERROR,*999)

                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_UNBOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_BOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_UNBOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(5)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_STRESS_CROSS_FIBRE_BOUND,ERR,ERROR,*999)
                  dof_idx=FIELD_VARIABLE_FE%COMPONENTS(6)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,ne)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,dof_idx,ACTIVATION,ERR,ERROR,*999)
                ENDIF

              ENDDO !gauss_idx
            ENDDO !element_idx

            !now the ghost elements -- get the relevant info from the other computational nodes
            CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Independent field interpolation is not implemented for problem subtype " &
            & //TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("ReacDiffFiniteElastic_IndependentFieldInterpolate")
    RETURN
999 ERRORS("ReacDiffFiniteElastic_IndependentFieldInterpolate",ERR,ERROR)
    EXITS("ReacDiffFiniteElastic_IndependentFieldInterpolate")
    RETURN 1

  END SUBROUTINE ReacDiffFiniteElastic_IndependentFieldInterpolate

  !
  !================================================================================================================================
  !


END MODULE REACTION_DIFFUSION_FINITE_ELASTICITY_ROUTINES

