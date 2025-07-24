# MathPackagesTraining
A gh-pages site to host SWC style training materials for various HPC math packages


# Website

The public site for this repo is https://xsdk-project.github.io/MathPackagesTraining2025/.
After pushing to the Repo, changes should be visible within minutes.

## To Render Locally

install Jekyll, see http://jekyllrb.com/

install Ruby dependencies:
```
bundle install
```

(See [this stackoverflow question](https://stackoverflow.com/questions/40385493/how-to-run-bundle-install-as-normal-user) if `bundle` complains about write permissions.)

Clone or move to the MathPackagesTraining2025 directory and start the Jekyll server:

```
git clone https://github.com/xsdk-project/MathPackagesTraining2025.git
bundle exec jekyll serve
```

Then point your web broswer at http://localhost:4000/MathPackagesTraining2025/


# Polaris

## Accounts

If you have an active ACLF account and are a member of ATPESC_Instructors, then
you can access polaris node now. **Anyone who will be participating as an
instructor and does not have an ALCF account or does not have access to the
ATPESC_Instructors project should request access to those ASAP**.

You need to be a part of two groups on polaris:

- `ATPESC_Instructors`: to submit jobs on the account and for write access the installation directory for libraries, `/eagle/ATPESC2025/usr/MathPackages`
- `ATPESC2025`: for write access to the examples directory `/eagle/ATPESC2025/EXAMPLES/track-5-numerical`

## Quick Start

See also [Getting Started on Polaris](https://docs.alcf.anl.gov/polaris/getting-started/)

To connect:

```
ssh polaris.alcf.anl.gov
```

All compute nodes [are the same](https://docs.alcf.anl.gov/polaris/hardware-overview/machine-overview/): a 32 core EPYC Milan 7543P node with 4 NVIDIA A100 GPUs connected via NVLink.

To request an interactive session:

```
qsub -I -l select=1 -l filesystems=home:eagle -l walltime=1:00:00 -q debug -A ATPESC_Instructors
```

More control over running non-interactive jobs is described in [Running jobs on Polaris](https://docs.alcf.anl.gov/polaris/running-jobs/)

### Modules

The following module commands have been tested and found to work when building and installing both Trilinos and PETSc/TAO:

```
module use /soft/modulefiles
module use /eagle/ATPESC2025/usr/modulefiles
module load track-5-numerical 
```
For reference:

```
polaris-login-01:/eagle/ATPESC2025/usr/modulefiles/track-5-numerical> cat 2025.lua
family("atpesct_track_5")
load("PrgEnv-gnu")
load("nvhpc-mixed")
load("craype-accel-nvidia80")
load("spack-pe-base")
load("cmake")
load("ninja")
load("cray-libsci")
```

## CUDA options

The A100 GPU has `CUDA Capability: 8.0`  i.e the corresponding compile options are:

```
nvcc -gencode arch=compute_80,code=sm_80
```

With cmake - the likely option is: `-DCMAKE_CUDA_ARCHITECTURES=80`

## Install software

Install software at `/eagle/projects/ATPESC2025/usr/MathPackages` - for ex: `/eagle/projects/ATPESC2025/usr/MathPackages/petsc-3.19.4`

And then copy over needed tutorial binaries, datafiles etc. over to `/eagle/projects/ATPESC2025/EXAMPLES/track-5-numerical` into appropriate folders - for ex: (from last year)

```
balay@thetagpu06:~$ ls -l /eagle/projects/ATPESC2025/EXAMPLES/track-5-numerical
total 40
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 amrex
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 hand_coded_heat
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 krylov_amg_hypre
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 krylov_amg_muelu
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 mfem-pumi-lesson
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:33 nonlinear_solvers_petsc
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:21 numerical_optimization_tao
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 rank_structured_strumpack
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 superlu
drwxrwsr-x 2 balay ATPESC_Instructors 4096 Aug  2 12:13 time_integrators_sundials
```

It is recommended to use a compute node when building and installing, because parallel make is limited on login nodes, and because the login nodes do not have GPUs (though they do have the GPU compilers and environments: building on the host is possible as long as running code on the GPU is not required).


## Internet proxy

If you need internet access from a node (for instance, to download packages) add the proxy commands to your environment given [here](https://docs.alcf.anl.gov/polaris/getting-started/?h=proxy).

