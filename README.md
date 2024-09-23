# PGABLE V0.2.1

PGABLE, Copyright (c) 2024, University of Waterloo  
Copying, use and development for non-commercial purposes permitted. All rights for commercial use reserved; for more information contact Stephen Mann (smann@uwaterloo.ca).  
This software is unsupported.

## Introduction

PGABLE is a Matlab toolkit for teaching and learning about geometric algebra. It is the successor of GABLE.

## Installation

To install PGABLE, run `addpath(genpath('.../PGABLE'))` where `.../PGABLE` is the path to the PGABLE folder.

## Basics

Geometric algebra has various models. These various models differ in two ways:
1.  They may differ in the underlying Clifford algebra.
2.  They may differ in the geometric interpretation of their elements.

PGABLE currently implements Projective (or Plane-based) Geometric Algebra, called PGA, and Ordinary Geometric Algebra, called OGA. It is assumed that the user will be using one model at a time.

To receive help with this plugin, run `help PGABLE`.