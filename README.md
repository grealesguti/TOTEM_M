# README

## Introduction

This code provides a FEM code for the simulation and optimization of thermoeletromechanical devices.



## Input file

This file describes all boundary conditions, other input files, material properties and topology optimization operators.

Multiple examples can be found in the Benchmarks and TECTO folders.

### Material properties

### Nodal and boundary conditions

### TO operators supported

## Mesh file
This code is expected to work with .inp mesh files with the Abaqus format.

All elements in the mesh must be of the same type. Allowed elements include:

Surfaces
* 
* 

Volumes
* 
*

### Add new element types

## Using gmsh

* Create the mesh
* Export using .inp format
* Export all nodes and elements!

## Running the code