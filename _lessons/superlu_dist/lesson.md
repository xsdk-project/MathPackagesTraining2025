---
layout: page-fullwidth
order: 7
title: "Sparse, Direct Solvers with SuperLU"
subheadline: "Linear Solvers"
teaser: "Role and Use of Direct Solvers in Ill-Conditioned Problems"
permalink: "lessons/superlu_dist/"
use_math: true
lesson: true
header:
 image_fullwidth: "matrices.png"
---

## To begin this lesson

- Enter the lesson directory
```
cd {{site.handson_root}}/superlu_dist
```

There are two folders:
```
EXAMPLE/ : example drivers

```

## The Example Source Code

README file explains various example drivers.

## Running the Example

The current directory contains relatively small sparse matrices: g4.rua, g20.rua, and big.rua.

Set the following directory that contains larger matrix files:
```
export matdir=/eagle/ATPESC2025/usr/MathPackages/datafiles
```

Set to use one OpenMP thread:
```
export OMP_NUM_THREADS=1
```

### Run 1 (test 2D algorithm driver, with GPU)
```
mpiexec -n 2 pddrive -r 1 -c 2 ${matdir}/torso3.mtx 2>&1 | tee run1.txt
```

#### Expected Behavior/Output
```
**************************************************
.. options:
**    Fact                      :    0
**    Equil                     :    1
**    DiagInv                   :    1
**    ParSymbFact               :    0
**    ColPerm                   :    4
**    RowPerm                   :    1
**    ReplaceTinyPivot          :    0
**    IterRefine                :    2
**    Trans                     :    0
**    SymPattern                :    0
**    lookahead_etree           :    0
**    Use_TensorCore            :    0
**    Use 3D algorithm          :    0
**    num_lookaheads            :   10
** parameters that can be altered by environment variables:
**    superlu_relax             :   60
**    superlu_maxsup            :  256
**    min GEMM m*k*n to use GPU : 5000
**    GPU buffer size           :  256000000
**    GPU streams               :    4
**    estimated fill ratio      :    5
**************************************************
       Matrix size min_mn    259156
       Nonzeros in L       122213768
       Nonzeros in U       122247101
       nonzeros in L+U     244201713
       nonzeros in LSUB    16559556

GPU Driver version:   v 11040
GPU Devices:

0 : NVIDIA A100-SXM4-40GB 8 0
  Global memory:   40536 mb
  Shared memory:   48 kb
  Constant memory: 64 kb
  Block registers: 65536
	
** Memory Usage **********************************
** Total highmark (MB):
    Sum-of-all :  2548.54 | Avg :  1274.27  | Max :  1289.68
    Max at rank 1, different stages (MB):
        . symbfact          153.10
        . distribution      735.40
      	. numfact          1258.87
** NUMfact space (MB): (sum-of-all-processes)
    L\U :         2078.13 |  Total :  2548.54
    	. max at rank 1, max L+U memory (MB):  1052.52
        . max at rank 1, peak buffer (MB):      237.15
**************************************************

** number of Tiny Pivots:        0

   Sol  0: ||X-Xtrue||/||X|| = 3.972251e-15
**************************************************
**** Time (seconds) ****
     EQUIL time            0.017
     ROWPERM time          0.086
     COLPERM time          2.839
     SYMBFACT time         0.549
     DISTRIBUTE time       1.693
     FACTOR time          11.032
     Factor flops  3.479445e+11	Mflops	31540.87
     SOLVE time            0.217
     Solve flops	   4.846323e+08	Mflops	 2230.42
     REFINEMENT time       0.465	 Steps       2
```

### Run 2 (test 3D algorithm driver, with GPU)
```
mpiexec -n 2 pddrive3d -r 1 -c 1 -d 2  ${matdir}/torso3.mtx 2>&1 | tee run2.txt
```

#### Expected Behavior/Output
```
  ....

     FACTOR time           5.213
  ....
```

### Run 3 (test 3D algorithm driver, without GPU)
```
export SUPERLU_ACC_OFFLOAD=0
mpiexec -n 2 pddrive3d -r 1 -c 1 -d 2  ${matdir}/torso3.mtx 2>&1 | tee run3.txt
```

#### Expected Behavior/Output
```
 ....
     FACTOR time          14.373
 ....
```

### Run 4 (try different sparsity ordering options)
```
mpiexec -n 1 pddrive -q 0 big.rua
mpiexec -n 1 pddrive -q 1 big.rua
mpiexec -n 1 pddrive -q 2 big.rua
mpiexec -n 1 pddrive -q 3 big.rua
mpiexec -n 1 pddrive -q 4 big.rua
mpiexec -n 1 pddrive -q 5 big.rua
```
#### Expected Behavior/Output
```
0: NATURAL,         nonozeros in L+U: 15497854
1: MMD_ATA,         nonozeros in L+U:   248854
2: MMD_AT_PLUS_A,   nonozeros in L+U:   262794
3: COLAMD,          nonozeros in L+U:   210766
4: METIS_AT_PLUS_A, nonozeros in L+U:   245710
5: PARMETIS,        nonozeros in L+U:   245710
```

### Further Reading

https://github.com/xiaoyeli/superlu_dist
