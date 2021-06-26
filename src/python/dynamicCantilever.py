#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script for a dynamically loaded cantilever using OpenCMISS calls in python.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): Chris Bradley
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> Main script
# Add Python bindings directory to PATH
import sys, os
import math

#import fpectl

# Intialise OpenCMISS-Iron
from opencmiss.iron import iron

#fpectl.turnoff_sigfpe()

iron.PetscOptionsSetValue("-ksp_view_rhs_binary","")
iron.PetscOptionsSetValue("-ksp_view_mat_binary","")

ELEMENT_CONSTANT = 0
LINEAR_LAGRANGE = 1
QUADRATIC_LAGRANGE = 2
CUBIC_LAGRANGE = 3
CUBIC_HERMITE = 4

INCOMPRESSIBLE_MATERIAL = 1
COMPRESSIBLE_MATERIAL = 2

ST_VENANT_KIRCHOFF_MATERIAL = 1
MOONEY_RIVLIN_MATERIAL = 2

STATIC = 1
DYNAMIC = 2

DIRICHLET_BCS = 1
NEUMANN_BCS = 2

LUMPED = 1
UNLUMPED = 2

# Set problem parameters

# Cantilever dimensions
length = 50.0 # mm
height = 10.0 # mm
width = 10.0 # mm

#timeDependence = STATIC
timeDependence = DYNAMIC

# Time
startTime = 0.0 # ms
stopTime = 2000.1 # ms
timeStep = 50 # ms
#timeStep = 5 # ms

# Loading

boundaryConditionType = DIRICHLET_BCS

if (boundaryConditionType == DIRICHLET_BCS):
    maxDisplacement = -0.33*height;
elif (boundaryConditionType == NEUMANN_BCS):
    #maxForce = 0.0 # N.mm^-2
    maxForce = -50.0 # N.mm^-2
else:
    print('Invalid boundary condition type')
    exit()
   
frequency = 0.001 # kHz or ms^-1
phase = 0.0

numberOfLoadIncrements = 1

# Material properties
materialLaw = ST_VENANT_KIRCHOFF_MATERIAL
#materialLaw = MOONEY_RIVLIN_MATERIAL

#materialCompressibilty = INCOMPRESSIBLE_MATERIAL
materialCompressibility = COMPRESSIBLE_MATERIAL

# Look at silicon rubber
density = 1.1 # mg.mm^-3
#density = 0.0 # mg.mm^-3

#massLumping = LUMPED
massLumping = UNLUMPED

if (materialLaw == ST_VENANT_KIRCHOFF_MATERIAL):
    poissonsRatio = 0.47
    youngsModulus = 1000.0 # mg.mm^-1.ms^-2
    lameLambda = poissonsRatio*youngsModulus/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio))
    lameMu = youngsModulus/(2.0*(1.0+poissonsRatio))
elif (materialLaw == MOONEY_RIVLIN_MATERIAL):
    mooneyRivlin1 = 800.0 # mg.mm^-1.ms^-2
    mooneyRivlin2 = -7000.0 # mg.mm^-1.ms^-2
else:
    print('Invalid material law')
    exit()


if (materialLaw == MOONEY_RIVLIN_MATERIAL):
    pInit = -mooneyRivlin1 # Initial hydrostatic pressure
else:
    pInit = 0.0 # Initial hydrostatic pressure

pRef = pInit # Reference hydrostatic pressure

numberOfLoadIncrements = 1

# Interpolation
#geometricInterpolation = CUBIC_HERMITE
geometricInterpolation = LINEAR_LAGRANGE
#geometricInterpolation = QUADRATIC_LAGRANGE
fibreInterpolation = geometricInterpolation
pressureInterpolation = ELEMENT_CONSTANT
#pressureInterpolation = LINEAR_LAGRANGE

# Should not need to change anything below here

coordinateSystemUserNumber = 1
regionUserNumber = 1
linearLagrangeBasisUserNumber = 1
quadraticLagrangeBasisUserNumber = 2
cubicLagrangeBasisUserNumber = 3
cubicHermiteBasisUserNumber = 4
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
decomposerUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
elasticityEquationsSetUserNumber = 1
elasticityEquationsSetFieldUserNumber = 3
elasticityDependentFieldUserNumber = 4
elasticityMaterialsFieldUserNumber = 5
bcCellMLUserNumber = 1
bcCellMLModelsFieldUserNumber = 6
bcCellMLStateFieldUserNumber = 7
bcCellMLParametersFieldUserNumber = 8
bcCellMLIntermediateFieldUserNumber = 9
elasticityProblemUserNumber = 1

if len(sys.argv) > 1:
    if len(sys.argv) > 6:
        sys.exit('Error: too many arguments- currently only accepting 5 options: numberXElements numberYElements numberZElements')
    if len(sys.argv) >= 4:
        numberOfGlobalXElements = int(sys.argv[1])
        numberOfGlobalYElements = int(sys.argv[2])
        numberOfGlobalZElements = int(sys.argv[3])
    else:
        numberOfGlobalXElements = 3
        numberOfGlobalYElements = 1
        numberOfGlobalZElements = 1
else:
    numberOfGlobalXElements = 3
    numberOfGlobalYElements = 1
    numberOfGlobalZElements = 1

haveLinearLagrange = False
haveQuadraticLagrange = False
haveCubicLagrange = False
haveCubicHermite = False
if (geometricInterpolation == LINEAR_LAGRANGE):
    haveLinearLagrange = True
    numberOfNodesXi = 2
    numberOfGaussXi = 2
    if(pressureInterpolation != ELEMENT_CONSTANT):
        sys.exit('Invalid pressure interpolation for linear Lagrange geometric interpolation')
elif (geometricInterpolation == QUADRATIC_LAGRANGE):
    haveQuadraticLagrange = True
    numberOfNodesXi = 3
    numberOfGaussXi = 3
    if(pressureInterpolation == ELEMENT_CONSTANT):
        a = 1
    elif (pressureInterpolation == LINEAR_LAGRANGE):
        haveLinearLagrange = True
    else:
        sys.exit('Invalid pressure interpolation for quadratic Lagrange geometric interpolation')
elif (geometricInterpolation == CUBIC_LAGRANGE):
    haveCubicLagrange = True
    numberOfNodesXi = 4
    numberOfGaussXi = 4
    if(pressureInterpolation == ELEMENT_CONSTANT):
        a
    elif (pressureInterpolation == LINEAR_LAGRANGE):
        haveLinearLagrange = True        
    elif (pressureInterpolation == QUADRATIC_LAGRANGE):
        haveQuadraticLagrange = True
    else:
        sys.exit('Invalid pressure interpolation for cubic Lagrange geometric interpolation')
elif (geometricInterpolation == CUBIC_HERMITE):
    haveCubicHermite = True
    numberOfGaussXi = 4
    if(pressureInterpolation == ELEMENT_CONSTANT):
        numberOfNodesXi = 2
    elif (pressureInterpolation == LINEAR_LAGRANGE):
        haveLinearLagrange = True        
        numberOfNodesXi = 2
    elif (pressureInterpolation == QUADRATIC_LAGRANGE):
        haveQuadraticLagrange = True
        numberOfNodesXi = 3
    elif (pressureInterpolation == CUBIC_LAGRANGE):
        haveCubicLagrange = True
        numberOfNodesXi = 4
    else:
        sys.exit('Invalid pressure interpolation for cubic Hermite geometric interpolation')
else:
    sys.exit('Invalid geometric interpolation')

if (numberOfGlobalZElements == 0):
    numberOfDimensions = 2
    numberOfElements = numberOfGlobalXElements*numberOfGlobalYElements
    numberOfXNodes = numberOfGlobalXElements*(numberOfNodesXi-1)+1
    numberOfYNodes = numberOfGlobalYElements*(numberOfNodesXi-1)+1
    numberOfNodes = numberOfXNodes*numberOfYNodes
else:
    numberOfDimensions = 3
    numberOfElements = numberOfGlobalXElements*numberOfGlobalYElements*numberOfGlobalZElements
    numberOfXNodes = numberOfGlobalXElements*(numberOfNodesXi-1)+1
    numberOfYNodes = numberOfGlobalYElements*(numberOfNodesXi-1)+1
    numberOfZNodes = numberOfGlobalZElements*(numberOfNodesXi-1)+1
    numberOfNodes = numberOfXNodes*numberOfYNodes*numberOfZNodes

numberOfXi = numberOfDimensions    
numberOfGauss = pow(numberOfGaussXi,numberOfXi)

worldRegion = iron.Region()
iron.Context.WorldRegionGet(worldRegion)

iron.OutputSetOn("DynamicCantilever")
   
# Get the computational nodes info
computationEnvironment = iron.ComputationEnvironment()
iron.Context.ComputationEnvironmentGet(computationEnvironment)
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()

worldWorkGroup = iron.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)

# Create a rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber,iron.Context)
coordinateSystem.DimensionSet(numberOfDimensions)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,worldRegion)
region.LabelSet("CantileverRegion")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

numberOfMeshComponents = 0
linearLagrangeMeshComponent = 0
quadraticLagrangeMeshComponent = 0
cubicLagrangeMeshComponent = 0
cubicHermiteMeshComponent = 0

if (haveCubicLagrange):
    # Define cubic Lagrange basis
    cubicLagrangeBasis = iron.Basis()
    cubicLagrangeBasis.CreateStart(cubicLagrangeBasisUserNumber,iron.Context)
    cubicLagrangeBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    cubicLagrangeBasis.NumberOfXiSet(numberOfXi)
    cubicLagrangeBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*numberOfXi)
    cubicLagrangeBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfXi)
    cubicLagrangeBasis.CreateFinish()
    numberOfMeshComponents = numberOfMeshComponents + 1
    cubicLagrangeMeshComponent = numberOfMeshComponents
    
if (haveQuadraticLagrange):
    # Define quadratic basis
    quadraticLagrangeBasis = iron.Basis()
    quadraticLagrangeBasis.CreateStart(quadraticLagrangeBasisUserNumber,iron.Context)
    quadraticLagrangeBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    quadraticLagrangeBasis.NumberOfXiSet(numberOfXi)
    quadraticLagrangeBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfXi)
    quadraticLagrangeBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfXi)
    quadraticLagrangeBasis.CreateFinish()
    numberOfMeshComponents = numberOfMeshComponents + 1
    quadraticLagrangeMeshComponent = numberOfMeshComponents
    
if (haveLinearLagrange):
    # Define linear Lagrange basis
    linearLagrangeBasis = iron.Basis()
    linearLagrangeBasis.CreateStart(linearLagrangeBasisUserNumber,iron.Context)
    linearLagrangeBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    linearLagrangeBasis.NumberOfXiSet(numberOfXi)
    linearLagrangeBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi)
    linearLagrangeBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfXi)
    linearLagrangeBasis.CreateFinish()
    numberOfMeshComponents = numberOfMeshComponents + 1
    linearLagrangeMeshComponent = numberOfMeshComponents

if (haveCubicHermite):    
    # Define cubic Hermite basis
    cubicHermiteBasis = iron.Basis()
    cubicHermiteBasis.CreateStart(cubicHermiteBasisUserNumber,iron.Context)
    cubicHermiteBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    cubicHermiteBasis.NumberOfXiSet(numberOfXi)
    cubicHermiteBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*numberOfXi)
    cubicHermiteBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfXi)
    cubicHermiteBasis.CreateFinish()
    numberOfMeshComponents = numberOfMeshComponents + 1
    cubicHermiteMeshComponent = numberOfMeshComponents

geometricMeshComponent = 0    
if (geometricInterpolation == LINEAR_LAGRANGE):
    geometricMeshComponent = linearLagrangeMeshComponent
elif (geometricInterpolation == QUADRATIC_LAGRANGE):
    geometricMeshComponent = quadraticLagrangeMeshComponent
elif (geometricInterpolation == CUBIC_LAGRANGE):
    geometricMeshComponent = cubicLagrangeMeshComponent
elif (geometricInterpolation == CUBIC_HERMITE):
    geometricMeshComponent = cubicHermiteMeshComponent

fibreMeshComponent = geometricMeshComponent

pressureMeshComponent = 0    
if (pressureInterpolation == ELEMENT_CONSTANT):
    pressureMeshComponent = geometricMeshComponent
if (pressureInterpolation == LINEAR_LAGRANGE):
    pressureMeshComponent = linearLagrangeMeshComponent
elif (pressureInterpolation == QUADRATIC_LAGRANGE):
    pressureMeshComponent = quadraticLagrangeMeshComponent

# Start the creation of a generated mesh in the region
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
                      
if (haveLinearLagrange):
    if (haveQuadraticLagrange):
        if (haveCubicLagrange):
            if (haveCubicHermite):
                generatedMesh.BasisSet([cubicLagrangeBasis,quadraticLagrangeBasis,linearLagrangeBasis,cubicHermiteBasis])
            else:
                generatedMesh.BasisSet([cubicLagrangeBasis,quadraticLagrangeBasis,linearLagrangeBasis])
        else:
            if (haveCubicHermite):
                generatedMesh.BasisSet([quadraticLagrangeBasis,linearLagrangeBasis,cubicHermiteBasis])
            else:
                generatedMesh.BasisSet([quadraticLagrangeBasis,linearLagrangeBasis])
    else:
        if (haveCubicLagrange):
            if (haveCubicHermite):
                generatedMesh.BasisSet([cubicLagrangeBasis,linearLagrangeBasis,cubicHermiteBasis])
            else:
                generatedMesh.BasisSet([cubicLagrangeBasis,linearLagrangeBasis])
        else:
            if (haveCubicHermite):
                generatedMesh.BasisSet([linearLagrangeBasis,cubicHermiteBasis])
            else:
                generatedMesh.BasisSet([linearLagrangeBasis])
else:
    if (haveQuadraticLagrange):
        if (haveCubicLagrange):
            if (haveCubicHermite):
                generatedMesh.BasisSet([cubicLagrangeBasis,quadraticLagrangeBasis,cubicHermiteBasis])
            else:
                generatedMesh.BasisSet([cubicLagrangeBasis,quadraticLagrangeBasis])
        else:
            if (haveCubicHermite):
                generatedMesh.BasisSet([quadraticLagrangeBasis,cubicHermiteBasis])
            else:
                generatedMesh.BasisSet([quadraticLagrangeBasis])
    else:
        if (haveCubicLagrange):
            if (haveCubicHermiteBasi):
                generatedMesh.BasisSet([cubicLagrangeBasis,cubicHermiteBasis])
            else:
                generatedMesh.BasisSet([cubicLagrangeBasis])
        else:
            if (haveCubicHermite):
                generatedMesh.BasisSet([cubicHermiteBasis])
            else:
                sys.exit('No basis functions have been used')

if (numberOfDimensions == 2):
    generatedMesh.ExtentSet([length,height])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements])
else:
    generatedMesh.ExtentSet([length,height,width])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements])
    
# Finish the creation of a generated mesh in the region
mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.CreateFinish()

# Create a decomposer
decomposer = iron.Decomposer()
decomposer.CreateStart(decomposerUserNumber,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

print("Geometric mesh component = %d" % geometricMeshComponent)

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,geometricMeshComponent)
if (numberOfDimensions == 3):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,geometricMeshComponent)
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.DecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,fibreMeshComponent)
fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,fibreMeshComponent)
if (numberOfDimensions == 3):
    fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,fibreMeshComponent)
fibreField.CreateFinish()

# Create the elasticity equations_set
elasticityEquationsSetField = iron.Field()
elasticityEquationsSet = iron.EquationsSet()
if (timeDependence == STATIC):
    if (materialLaw == ST_VENANT_KIRCHOFF_MATERIAL):
        print("Not implemented")
        exit()
    elif (materialLaw == MOONEY_RIVLIN_MATERIAL):
        elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                               iron.EquationsSetTypes.FINITE_ELASTICITY,
                                               iron.EquationsSetSubtypes.MOONEY_RIVLIN]
elif (timeDependence == DYNAMIC):
    if (materialCompressibility == INCOMPRESSIBLE_MATERIAL):
        if (materialLaw == ST_VENANT_KIRCHOFF_MATERIAL):
            elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                                   iron.EquationsSetTypes.FINITE_ELASTICITY,
                                                   iron.EquationsSetSubtypes.DYNAMIC_ST_VENANT_KIRCHOFF]
        elif (materialLaw == MOONEY_RIVLIN_MATERIAL):
            elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                                   iron.EquationsSetTypes.FINITE_ELASTICITY,
                                                   iron.EquationsSetSubtypes.DYNAMIC_MOONEY_RIVLIN]
    elif (materialCompressibility == COMPRESSIBLE_MATERIAL):
        if (materialLaw == ST_VENANT_KIRCHOFF_MATERIAL):
            elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                                   iron.EquationsSetTypes.FINITE_ELASTICITY,
                                                   iron.EquationsSetSubtypes.DYNAMIC_COMP_ST_VENANT_KIRCHOFF]
        elif (materialLaw == MOONEY_RIVLIN_MATERIAL):
            elasticityEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                                   iron.EquationsSetTypes.FINITE_ELASTICITY,
                                                   iron.EquationsSetSubtypes.DYNAMIC_COMP_MOONEY_RIVLIN]
    else:
        print("Invalid material compressibility")
        exit()
else:
    print("Invalid time dependence")
    exit()
    
elasticityEquationsSet.CreateStart(elasticityEquationsSetUserNumber,region,fibreField, \
                                   elasticityEquationsSetSpecification,elasticityEquationsSetFieldUserNumber, \
                                   elasticityEquationsSetField)
elasticityEquationsSet.CreateFinish()

# Create the dependent field
elasticityDependentField = iron.Field()
elasticityEquationsSet.DependentCreateStart(elasticityDependentFieldUserNumber,elasticityDependentField)
elasticityDependentField.VariableLabelSet(iron.FieldVariableTypes.U,"ElasticityDependent")
elasticityDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"ElasticityTraction")
elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,geometricMeshComponent)
elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,geometricMeshComponent)
elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,geometricMeshComponent)
if (numberOfDimensions == 3):
    elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,geometricMeshComponent)
    elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,geometricMeshComponent)
if (materialCompressibility == INCOMPRESSIBLE_MATERIAL):
    # Set the pressure to be nodally based and use the second mesh component
    if (pressureInterpolation == ELEMENT_CONSTANT):
        elasticityDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,numberOfDimensions+1,\
                                                           iron.FieldInterpolationTypes.ELEMENT_BASED)
        elasticityDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,\
                                                           iron.FieldInterpolationTypes.ELEMENT_BASED)                
    else:
        elasticityDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,numberOfDimensions+1,\
                                                           iron.FieldInterpolationTypes.NODE_BASED)
        elasticityDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,\
                                                           iron.FieldInterpolationTypes.NODE_BASED)
    elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,numberOfDimensions+1,pressureMeshComponent)
    elasticityDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,pressureMeshComponent)
elasticityEquationsSet.DependentCreateFinish()

# Initialise elasticity dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
    elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
if (numberOfDimensions == 3):
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
        elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
if (materialCompressibility == INCOMPRESSIBLE_MATERIAL):
    iron.Field.ComponentValuesInitialiseDP(
        elasticityDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,numberOfDimensions+1,pInit)

# Create the material field
elasticityMaterialsField = iron.Field()
elasticityEquationsSet.MaterialsCreateStart(elasticityMaterialsFieldUserNumber,elasticityMaterialsField)
elasticityMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
elasticityEquationsSet.MaterialsCreateFinish()

# Set materials parameters
if (timeDependence == STATIC):
    if (materialLaw == ST_VENANT_KIRCHOFF_MATERIAL):
        print("Not implemented")
        exit()
    elif (materialLaw == MOONEY_RIVLIN_MATERIAL):    
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             1,mooneyRivlin1)
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             2,mooneyRivlin2)
elif (timeDependence == DYNAMIC):
    if (materialLaw == ST_VENANT_KIRCHOFF_MATERIAL):
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             1,density)
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             2,lameLambda)
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             3,lameMu)
    elif (materialLaw == MOONEY_RIVLIN_MATERIAL):    
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             1,density)
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             2,mooneyRivlin1)
        elasticityMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             3,mooneyRivlin2)
else:
    print("Invalid time dependence")
    exit()
    
# Create elasticity equations
elasticityEquations = iron.Equations()
elasticityEquationsSet.EquationsCreateStart(elasticityEquations)
elasticityEquations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
elasticityEquations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
elasticityEquations.OutputTypeSet(iron.EquationsOutputTypes.MATRIX)
if (massLumping == LUMPED):
    elasticityEquations.LumpingTypeSet(iron.EquationsLumpingTypes.LUMPED)
elasticityEquationsSet.EquationsCreateFinish()

# Use analytic Jacobian
elasticityEquations.JacobianCalculationTypeSet(1,iron.FieldVariableTypes.U,iron.EquationsJacobianCalculated.ANALYTIC)
#elasticityEquations.JacobianCalculationTypeSet(1,iron.FieldVariableTypes.U,iron.EquationsJacobianCalculated.FINITE_DIFFERENCE)
#elasticityEquations.JacobianFiniteDifferenceStepSizeSet(1,iron.FieldVariableTypes.U,1.0e-8)

if(timeDependence == DYNAMIC):
    # Create CellML equations for the temporal solid boundary conditions
    bcCellML = iron.CellML()
    bcCellML.CreateStart(bcCellMLUserNumber,region)
    if (boundaryConditionType == DIRICHLET_BCS):
        bcCellMLIdx = bcCellML.ModelImport("input/sinusoiddispbc.cellml")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/maxDisplacement")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/frequency")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/phase")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/geometricX")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/geometricY")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/geometricZ")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/deformedX")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/deformedY")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/deformedZ")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/velocityX")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/velocityY")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/velocityZ")
    else:
        bcCellMLIdx = bcCellML.ModelImport("input/sinusoidforcebc.cellml")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/maxForce")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/frequency")
        bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/phase")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/forceX")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/forceY")
        bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/forceZ")
    bcCellML.CreateFinish()

    # Create CellML <--> OpenCMISS field maps
    bcCellML.FieldMapsCreateStart()
    if (boundaryConditionType == DIRICHLET_BCS):
        # Map geometric positions 
        bcCellML.CreateFieldToCellMLMap(geometricField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                                bcCellMLIdx,"main/geometricX",iron.FieldParameterSetTypes.VALUES)
        bcCellML.CreateFieldToCellMLMap(geometricField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                                bcCellMLIdx,"main/geometricY",iron.FieldParameterSetTypes.VALUES)
        if (numberOfDimensions == 3):
            bcCellML.CreateFieldToCellMLMap(geometricField,iron.FieldVariableTypes.U,3,iron.FieldParameterSetTypes.VALUES,
	                                    bcCellMLIdx,"main/geometricZ",iron.FieldParameterSetTypes.VALUES)
        # Map deformed position to ensure dependent field isn't cleared when the forces are copied back
        bcCellML.CreateFieldToCellMLMap(elasticityDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                                bcCellMLIdx,"main/deformedX",iron.FieldParameterSetTypes.VALUES)
        bcCellML.CreateFieldToCellMLMap(elasticityDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                                bcCellMLIdx,"main/deformedY",iron.FieldParameterSetTypes.VALUES)
        if (numberOfDimensions == 3):
            bcCellML.CreateFieldToCellMLMap(elasticityDependentField,iron.FieldVariableTypes.U,3,iron.FieldParameterSetTypes.VALUES,
	                                    bcCellMLIdx,"main/deformedZ",iron.FieldParameterSetTypes.VALUES)
        # Map displacement to dependent field
        bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/deformedX",iron.FieldParameterSetTypes.VALUES,
	                                elasticityDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
        bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/deformedY",iron.FieldParameterSetTypes.VALUES,
	                                elasticityDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES)
        if (numberOfDimensions == 3):
            bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/deformedZ",iron.FieldParameterSetTypes.VALUES,
	                                    elasticityDependentField,iron.FieldVariableTypes.U,3,iron.FieldParameterSetTypes.VALUES)
    else:
        # Map forces to ensure dependent field isn't cleared when the forces are copied back
        bcCellML.CreateFieldToCellMLMap(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,iron.FieldParameterSetTypes.VALUES,
	                                bcCellMLIdx,"main/forceX",iron.FieldParameterSetTypes.VALUES)
        bcCellML.CreateFieldToCellMLMap(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,2,iron.FieldParameterSetTypes.VALUES,
	                                bcCellMLIdx,"main/forceY",iron.FieldParameterSetTypes.VALUES)
        if (numberOfDimensions == 3):
            bcCellML.CreateFieldToCellMLMap(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,3,iron.FieldParameterSetTypes.VALUES,
	                                    bcCellMLIdx,"main/forceZ",iron.FieldParameterSetTypes.VALUES)
        # Map forces to dependent traction field
        bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/forceX",iron.FieldParameterSetTypes.VALUES,
	                                elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,iron.FieldParameterSetTypes.VALUES)
        bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/forceY",iron.FieldParameterSetTypes.VALUES,
	                                elasticityDependentField,iron.FieldVariableTypes.DELUDELN,2,iron.FieldParameterSetTypes.VALUES)
        if (numberOfDimensions == 3):
            bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/forceZ",iron.FieldParameterSetTypes.VALUES,
	                                    elasticityDependentField,iron.FieldVariableTypes.DELUDELN,3,iron.FieldParameterSetTypes.VALUES)
    bcCellML.FieldMapsCreateFinish()

    # Create the CellML models field
    bcCellMLModelsField = iron.Field()
    bcCellML.ModelsFieldCreateStart(bcCellMLModelsFieldUserNumber,bcCellMLModelsField)
    bcCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"BCModelMap")
    bcCellML.ModelsFieldCreateFinish()
    
    # Only evaluate BC at the end of the cantilever
    bcCellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0)
    if (numberOfDimensions == 2):
        nodeNumber = numberOfNodes
        nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
        if (nodeDomain == computationalNodeNumber):
            bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                           1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,bcCellMLIdx)
    else:
        for zNodeIdx in range(1,numberOfZNodes+1):
            nodeNumber=zNodeIdx*numberOfXNodes*numberOfYNodes
            nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
            if (nodeDomain == computationalNodeNumber):
                bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                               1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,bcCellMLIdx)

    # Create the CellML state field
    bcCellMLStateField = iron.Field()
    bcCellML.StateFieldCreateStart(bcCellMLStateFieldUserNumber,bcCellMLStateField)
    bcCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCState")
    bcCellML.StateFieldCreateFinish()
    
    # Create the CellML parameters field
    bcCellMLParametersField = iron.Field()
    bcCellML.ParametersFieldCreateStart(bcCellMLParametersFieldUserNumber,bcCellMLParametersField)
    bcCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"BCParameters")
    bcCellML.ParametersFieldCreateFinish()
    
    # Get the component numbers
    if (boundaryConditionType == DIRICHLET_BCS):
        maxDisplacementComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/maxDisplacement")
    else:
        maxForceComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/maxForce")
    frequencyComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/frequency")
    phaseComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/phase")
    # Set up the parameters field
    if (boundaryConditionType == DIRICHLET_BCS):
        bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,maxDisplacementComponentNumber,maxDisplacement)
    else:
        bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,maxForceComponentNumber,maxForce)
    bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,frequencyComponentNumber,frequency)
    bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,phaseComponentNumber,phase)

    # Create the CELL intermediate field
    bcCellMLIntermediateField = iron.Field()
    bcCellML.IntermediateFieldCreateStart(bcCellMLIntermediateFieldUserNumber,bcCellMLIntermediateField)
    bcCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCIntermediate")
    bcCellML.IntermediateFieldCreateFinish()
    
# Define the elasticity problem
elasticityProblem = iron.Problem()
if (timeDependence == STATIC):
    elasticityProblemSpecification = [iron.ProblemClasses.ELASTICITY,
                                      iron.ProblemTypes.FINITE_ELASTICITY,
                                      iron.ProblemSubtypes.STATIC_FINITE_ELASTICITY]
elif (timeDependence == DYNAMIC):
    elasticityProblemSpecification = [iron.ProblemClasses.ELASTICITY,
                                      iron.ProblemTypes.FINITE_ELASTICITY,
                                      iron.ProblemSubtypes.DYNAMIC_FINITE_ELASTICITY]
else:
    print("Invalid time dependence")
    exit()
    
elasticityProblem.CreateStart(elasticityProblemUserNumber,iron.Context,elasticityProblemSpecification)
elasticityProblem.CreateFinish()

# Create the elasticity problem control loop
elasticityProblem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
elasticityProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
if (timeDependence == STATIC):
    controlLoop.LabelSet('LoadIncrementLoop')
    controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
elif (timeDependence == DYNAMIC):   
    controlLoop.LabelSet('TimeLoop')
    controlLoop.TimesSet(startTime,stopTime,timeStep)
    controlLoop.TimeOutputSet(1)
else:
    print("Invalid time dependence")
    exit()
elasticityProblem.ControlLoopCreateFinish()

# Create elasticity problem solvers
bcCellMLEvaluationSolver = iron.Solver()
elasticityDynamicSolver = iron.Solver()
elasticityNonlinearSolver = iron.Solver()
elasticityLinearSolver = iron.Solver()
elasticityProblem.SolversCreateStart()
if (timeDependence == STATIC):
    elasticityProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,elasticityNonlinearSolver)
elif (timeDependence == DYNAMIC):
    # Get the BC CellML solver
    elasticityProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
    bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
    # Get the dynamic solver
    elasticityProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,elasticityDynamicSolver)
    #elasticityDynamicSolver.OutputTypeSet(iron.SolverOutputTypes.MONITOR)
    elasticityDynamicSolver.OutputTypeSet(iron.SolverOutputTypes.MATRIX)
    elasticityDynamicSolver.DynamicSchemeSet(iron.DynamicSchemeTypes.NEWMARK1)
    #elasticityDynamicSolver.DynamicSchemeSet(iron.DynamicSchemeTypes.NEWMARK2)
    elasticityDynamicSolver.DynamicNonlinearSolverGet(elasticityNonlinearSolver)
else:
    print("Invalid time dependence")
    exit()    
elasticityNonlinearSolver.OutputTypeSet(iron.SolverOutputTypes.MONITOR)
#elasticityNonlinearSolver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
elasticityNonlinearSolver.OutputTypeSet(iron.SolverOutputTypes.MATRIX)
#elasticityNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
elasticityNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS)
elasticityNonlinearSolver.NewtonAbsoluteToleranceSet(1e-8)
elasticityNonlinearSolver.NewtonSolutionToleranceSet(1e-8)
elasticityNonlinearSolver.NewtonRelativeToleranceSet(1e-8)
elasticityNonlinearSolver.NewtonLinearSolverGet(elasticityLinearSolver)
elasticityLinearSolver.linearType = iron.LinearSolverTypes.DIRECT
elasticityProblem.SolversCreateFinish()

# Create elasticity solver equations and add elasticity equations set to solver equations
elasticitySolverEquations = iron.SolverEquations()
elasticityProblem.SolverEquationsCreateStart()
if (timeDependence == STATIC):
    elasticityNonlinearSolver.SolverEquationsGet(elasticitySolverEquations)
elif (timeDependence == DYNAMIC):
    elasticityDynamicSolver.SolverEquationsGet(elasticitySolverEquations)
else:
    print("Invalid time dependence")
    exit()    
elasticitySolverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
elasticityEquationsSetIndex = elasticitySolverEquations.EquationsSetAdd(elasticityEquationsSet)
elasticityProblem.SolverEquationsCreateFinish()

if (timeDependence != STATIC):
    bcEquations = iron.CellMLEquations()
    elasticityProblem.CellMLEquationsCreateStart()
    bcCellMLEvaluationSolver.CellMLEquationsGet(bcEquations)
    bcEquationsIndex = bcEquations.CellMLAdd(bcCellML)
    elasticityProblem.CellMLEquationsCreateFinish()
    
# Prescribe boundary conditions (absolute nodal parameters)
elasticityBoundaryConditions = iron.BoundaryConditions()
elasticitySolverEquations.BoundaryConditionsCreateStart(elasticityBoundaryConditions)

if (numberOfDimensions == 2):
    #Set left hand edge to be built in.
    for yNodeIdx in range(0,numberOfYNodes):
        nodeNumber=yNodeIdx*numberOfXNodes+1
        nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
        if (nodeDomain == computationalNodeNumber):
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,
                                                 iron.BoundaryConditionsTypes.FIXED,0.0)
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,
                                                 iron.BoundaryConditionsTypes.FIXED,0.0)
else:
    for zNodeIdx in range(1,numberOfZNodes+1):
        for yNodeIdx in range(1,numberOfYNodes+1):
            # Set left hand build in nodes ot no displacement
            nodeNumber=(zNodeIdx-1)*numberOfXNodes*numberOfYNodes+(yNodeIdx-1)*numberOfXNodes+1
            nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
            if (nodeDomain == computationalNodeNumber):
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,3,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)

if (boundaryConditionType == DIRICHLET_BCS):
    if (timeDependence == STATIC):
        initialDisplacement = maxDisplacement
    else:
        initialDisplacement = 0.0
    if (numberOfDimensions == 2):
        #Set downward displacement on the right hand edge (value will come from CellML)
        nodeNumber = numberOfNodes
        nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
        if (nodeDomain == computationalNodeNumber):
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,
                                                 iron.BoundaryConditionsTypes.FIXED,0.0)
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,
                                                 iron.BoundaryConditionsTypes.FIXED,initialDisplacement)
    else:
        # Set downward displacement on right-hand edge (value will come from CellML)
        for zNodeIdx in range(1,numberOfZNodes+1):
            nodeNumber=zNodeIdx*numberOfXNodes*numberOfYNodes
            nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
            if (nodeDomain == computationalNodeNumber):
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,
                                                     iron.BoundaryConditionsTypes.FIXED,initialDisplacement)
                elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,3,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
else:
    if (timeDependence == STATIC):
        initialForce = maxForce
    else:
        initialForce = 0.0
    if (numberOfDimensions == 2):
        #Set downward force on the right hand edge (value will come from CellML)
        nodeNumber = numberOfNodes
        nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
        if (nodeDomain == computationalNodeNumber):
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeNumber,1,
                                                 iron.BoundaryConditionsTypes.FIXED,0.0)
            elasticityBoundaryConditions.AddNode(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeNumber,2,
                                                 iron.BoundaryConditionsTypes.FIXED,initialForce)
    else:
        # Set downward force on right-hand edge (value will come from CellML)
        for zNodeIdx in range(1,numberOfZNodes+1):
            nodeNumber=zNodeIdx*numberOfXNodes*numberOfYNodes
            nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
            if (nodeDomain == computationalNodeNumber):
                elasticityBoundaryConditions.SetNode(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeNumber,1,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
                elasticityBoundaryConditions.SetNode(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeNumber,2,
                                                     iron.BoundaryConditionsTypes.FIXED,initialForce)
                elasticityBoundaryConditions.SetNode(elasticityDependentField,iron.FieldVariableTypes.DELUDELN,1,1,nodeNumber,3,
                                                     iron.BoundaryConditionsTypes.FIXED,0.0)
if (materialCompressibility == INCOMPRESSIBLE_MATERIAL):         
    # Set reference pressure
    if (pressureInterpolation == ELEMENT_CONSTANT):
        if (numberOfElements > 1):
            elementNumber = 1
            elementDomain = decomposition.ElementDomainGet(elementNumber)
            if (elementDomain == computationalNodeNumber):
                elasticityBoundaryConditions.SetElement(elasticityDependentField,iron.FieldVariableTypes.U,elementNumber,numberOfDimensions+1,
                                                        iron.BoundaryConditionsTypes.FIXED,pRef)
    else:
        nodeNumber = 1
        nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
        if (nodeDomain == computationalNodeNumber):
            elasticityBoundaryConditions.SetNode(elasticityDependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,numberOfDimensions+1,
                                                 iron.BoundaryConditionsTypes.FIXED,pRef)
                
elasticitySolverEquations.BoundaryConditionsCreateFinish()

if (timeDependence == DYNAMIC):
    #Set initial velocities
    if (numberOfDimensions == 2):
        #Set downward force on the right hand edge (value will come from CellML)
        for yNodeIdx in range(1,numberOfYNodes+1):
            for xNodeIdx in range(1,numberOfXNodes+1):
                nodeNumber=xNodeIdx+(yNodeIdx-1)*numberOfXNodes
                nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
                if (nodeDomain == computationalNodeNumber):
                    initialVelocity = 2.0*math.pi*frequency*maxDisplacement*float(xNodeIdx-1)/float(numberOfXNodes-1)
                    elasticityDependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.INITIAL_VELOCITY,
                                                                      1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,0.0)
                    elasticityDependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.INITIAL_VELOCITY,
                                                                      1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,initialVelocity)
    else:
        for zNodeIdx in range(1,numberOfZNodes+1):
            for yNodeIdx in range(1,numberOfYNodes+1):
                for xNodeIdx in range(1,numberOfXNodes+1):
                    nodeNumber=xNodeIdx+(yNodeIdx-1)*numberOfXNodes+(zNodeIdx-1)*numberOfXNodes*numberOfYNodes
                    nodeDomain = decomposition.NodeDomainGet(nodeNumber,geometricMeshComponent)
                    if (nodeDomain == computationalNodeNumber):
                        initialVelocity = 2.0*math.pi*frequency*maxDisplacement*float(xNodeIdx-1)/float(numberOfXNodes-1)
                        elasticityDependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.INITIAL_VELOCITY,
                                                                          1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,0.0)
                        elasticityDependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.INITIAL_VELOCITY,
                                                                          1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,initialVelocity)
                        elasticityDependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.INITIAL_VELOCITY,
                                                                          1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,0.0)
                        
# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("Cantilever","FORTRAN")
fields.ElementsExport("Cantilever","FORTRAN")
fields.Finalise()

# Solve the elasticity problem
elasticityProblem.Solve()
                    
# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("Cantilever","FORTRAN")
fields.ElementsExport("Cantilever","FORTRAN")
fields.Finalise()

