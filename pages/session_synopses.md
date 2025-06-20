---
layout: page
show_meta: false
title: "Session Synopses"
subheadline: "More details about each session"
permalink: "/session_synopses/"
---

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Introduction to Numerical Software
### Introduction to Numerical Software
[Slides](../files/ATPESC-2024-Track-5-Talk-1-Yang-IntroToNumericalSoftware.pdf)

In the first hour of our day with you we will give a brief introduction
to our agenda, logistics and a very brief overview of the various numerical packages
you will learn about. We will also discuss important components of computational science and
engineering and exascale computing. We will give an overview of the FASTMath Institute and the 
Exascale Computing Project and its achievements.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Structured Meshes (with AMReX)
### Structured Meshes (with AMReX)
[Slides](../files/ATPESC-2024-Track-5-MyersZhang-AMREX.pdf)

[Lesson](../lessons/amrex/)

Block-structured adaptive mesh refinement (AMR) provides a natural framework
in which to focus computing power on the most critical parts of the problem in
the most computationally efficient way possible.  AMReX supports the development
of block-structured AMR algorithms for solving systems of partial differential
equations (PDE's) and other algorithms that require structured mesh and/or
particle discretizations.   We will begin with an overview of block-structured
AMR, and then discuss the
features of AMReX we might want to use to solve a multiphysics problem on
machines from laptops to supercomputers.  Hands-on exercises will include passive
scalar advection with time-dependent adaptivity, and the use of native linear
solvers to impose incompressibility on a flow with particles around an obstacle.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Unstructured Meshes (with MFEM/PUMI)
### Unstructured Meshes (with MFEM/PUMI)
[Slides](../files/ATPESC-2024-Track-5-TomovShephard-unstructured-mesh-tools.pdf)

[Lesson](../lessons/mfem_convergence/)

Unstructured meshes can yield required levels of accuracy using fewer degrees of
freedom at the cost of more complex parallel data structures and algorithms. To
support the ability of application code developers to take advantage of unstructured
meshes, FASTMath develops core tools to support the development of unstructured
mesh simulation capabilities. This lecture will first introduce the highly extendible
MFEM high order finite element solver library and then overview the PUMI unstructured
mesh tools developed to support mesh adaptation, load balancing and PIC calculations.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Nonlinear Solvers (with PETSc)
### Nonlinear Solvers (with PETSc)
[Slides](../files/ATPESC-2024-Track-5-TobyIsaac-PETSc.pdf)

[Lesson](../lessons/nonlinear_solvers_petsc/)

We will begin with a quick overview of iterative solvers for nonlinear systems,
and then take a deeper look into Newton-Krylov methods and how to use them via
the PETSc Scalable Nonlinear Equation Solvers (SNES) component. We will do
some hands-on exploration with a classic computational fluid dynamics benchmark,
the lid-driven cavity problem. We will end by looking at how nonlinear composition
and preconditioning can be used to construct a wide array of nonlinear solvers from
the algorithmic building blocks in SNES, and demonstrate how these techniques can
handle particularly difficult nonlinearities. Time permitting, we will also present
some bonus or "extra credit" material briefly looking at how SNES solvers can be
run using PETSc's GPU back-ends.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Time Integration (with SUNDIALS)
### Time Integration (with SUNDIALS)
[Slides](../files/ATPESC-2024-Track-5-TimeIntegration.pdf)

[Lesson](../lessons/time_integrators_sundials/)

In this lecture we will discuss the role and impact of high order, adaptive, and
flexible time integration libraries in solution accuracy and computational
efficiency of large-scale simulations.  Due to the wide variety of
backgrounds among ATPESC participants, we will briefly discuss

* the location of time integrators in the HPC landscape, and their reliance on scalable nonlinear and linear solver libraries,
* the different categories of time integration methods (explicit/implicit/IMEX),
* the basic theoretical properties of time integration methods (order of accuracy, linear stability),
* the role of temporal adaptivity for improving accuracy and efficiency,
* an overview of DOE time integration packages.

We will spend approximately half of the time period in lecture, followed by
hands-on exercises that examine stability, accuracy, temporal adaptivity,
and the role of problem-specific preconditioning.  All of the hands-on
exercises focus on time-dependent PDEs, and use the SUNDIALS' ARKODE
library for time integration, along with the AMReX library for spatial
semi-discretization.  These exercises utilize both MPI and GPUs, allowing
participants to explore both on-node and multi-node performance.


<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Iterative Solvers & Algebraic Multigrid (with Trilinos, Belos & MueLu)
### Iterative Solvers & Algebraic Multigrid (with Trilinos, Belos & MueLu)
[Slides](../files/ATPESC-2024-Track-5-GlusaHarper-IterativeLinearSolversAndAlgebraicMultigridWithTrilinosBelosMueLu.pdf)

[Lesson](../lessons/krylov_amg_muelu/)

In this session, attendees will learn about linear solvers and preconditioners
available in the Trilinos project.  We will focus on Krylov solvers such
as conjugate gradients (CG) and generalized minimum residual (GMRES); simple
preconditioners like Jacobi, Gauss-Seidel, and Chebyshev polynomials; and
scalable aggregation-based algebraic multigrid preconditioning.  The two
hands-on lessons will provide an opportunity to run a variety of stand-alone
examples that demonstrate some of the many Trilinos solver capabilities on a
model linear problem.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Krylov Solvers & Algebraic Multigrid (with HYPRE)
### Krylov Solvers & Algebraic Multigrid (with HYPRE)
[Slides](../files/ATPESC-2024-Track-5-Iterative-Solvers-hypre-DOK.pdf)

[Lesson](../lessons/krylov_amg_hypre/)

This session will present the basic concepts of iterative linear solvers with focus on
Krylov solvers, including the generalized minimum residual method (GMRES),
preconditioning and algebraic multigrid (AMG) methods. We will provide a brief
description of the high performance linear solvers library HYPRE, its
interfaces, and its most used multigrid solvers, BoomerAMG and PFMG, including
a brief discussion of the effect of their data structures on performance.
The lesson includes hands-on examples with structured and unstructured solvers
from the HYPRE library applied to several test problems on CPUs and GPUs.


<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Optimization (with TAO)
### Optimization (with TAO)
[Slides](../files/ATPESC-2024-Track-5-ToddMunson-TAO.pdf)

[Lesson](../lessons/numerical_optimization_tao/)

This lecture will provide an introduction to numerical optimization with a
theoretical focus on simulation-based problems. We will introduce the user
interfaces for the Toolkit for Advanced Optimization (TAO) package within
the PETSc library and exercise several gradient-based algorithms on a scalable
synthetic test problem. We will observe and discuss the relative convergence
of different classes of algorithms and sensitivity analysis methods in a
parallel environment. Finally, we will also introduce constraints and study
their effect on the problem solution.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Direct Solvers (with SuperLU/STRUMPACK)
### Direct Solvers (with SuperLU/STRUMPACK)
[Slides](../files/ATPESC-2024-Track-5-Li-Ghysels-DirectSolvers.pdf)

[Lesson](../lessons/superlu_dist/)

In this lecture we discuss several aspects of direct solvers for large
sparse linear systems. Sparse direct solvers, based on Gaussian
elimination, are popular for a variety of applications mainly because
of their numerical robustness and ease of use. Two software packages
are introduced: SuperLU(_Dist) and STRUMPACK. We give a brief
introduction to the underlying algorithms used in these solvers and
illustrate their use.

In the second part of the lecture we introduce the concept of
data-sparse matrix approximation, where a matrix is hierarchically
partitioned in sub-blocks, which are compressed using low-rank
approximations. Such a hierarchical matrix approximation allows for
linear algebra operations that are much more efficient than the
traditional O(n^3) dense linear algebra operations, such as LU
decomposition or matrix-matrix multiplication, while also reducing
memory usage. We then show how STRUMPACK uses such hierarchical matrix
approximation to construct efficient preconditioners based on
approximate sparse LU factorization.

* Session 1 (11:15am-11:45am) includes a hands-on demo session with
SuperLU.
* Session 2 (11:45am-12:15pm) has a hands-on demo session with STRUMPACK.
* Session 3 (12:15pm-12:30pm) includes a Q&A session (no hands-on).


<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Working with Numerical Packages in Practice
### Working with Numerical Packages in Practice
[Slides](../files/ATPESC-2024-Track-5-Talk-9-AnnAlmgren-wrapup.pdf)

Developing high-quality, large-scale scientific computing applications in science and
engineering involves expertise in many areas. Typically, no one person or group has
all of the essential expertise and skills. Numerical software libraries and packages
are a key way we share capability and know-how. Learning to leverage numerical
packages to address new scientific computing challenges is part of becoming a member of
the scientific computing community. In this wrap-up session of the day, we
briefly discuss key trade-offs in using numerical packages in practice.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### Extreme-Scale Numerical Algorithms and Software (Panel)
### Extreme-Scale Numerical Algorithms and Software (Panel)

This will be a 45 minute panel question and answer period for ATPESC learners to ask
questions about working with numerical package and the community of numerical package
developers. If you have question(s) you know you would like to ask, we encourage attendees
to submit questions *ahead* of time via the
[submission form](https://forms.gle/SmVroQBAq145g4F79).
However, we expect there will also be ample opportunity to indicate your desire to ask
questions via the main (Amphitheater) [slack channel](https://app.slack.com/client/TMW2FLNCQ/C029YJ95ULQ),
and then you may be called upon to un-mute and ask your question.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

### SME Speed Dating OPTIONAL ACTIVITY
### SME Speed Dating OPTIONAL ACTIVITY

This is an opportunity for you to meet and talk (1-on-1 or sometimes 2-on-1) with subject matter experts
(SMEs) about their work and numerical packages they support. Each *speed date* will be a 20 minute
meeting. Attendees may select up to
3 SMEs they would like to meet using
[this form](https://docs.google.com/document/d/1-uCcBU5zOBmVkGOqOT2dg4eXDtLl6TA3wfHfBgcRfi0/edit?usp=sharing)
Requests will be accommodated on a first-come, first-served basis until all available slots are taken.

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
