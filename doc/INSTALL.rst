Setting up your account to use Watermelon
=========================================
Watermelon is installed on comp3/5.
You can run it by adding this line to your ~/.bashrc:

module use /nfs/med-bfx-common/software/bfx_modules

Once added to bashrc, you can load watermelon so:

::

  $ module load watermelon
   loaded watermelon/0.1

  $ watermelon --help
   ...

FWIW, your experience with "screen" will be better if you add a file ~/.screenrc 
containing these two lines:

caption     always        "%{+b rk}%H%{gk} |%c %{yk}%d.%m.%Y | %72=Load: %l %{wk}"
hardstatus alwayslastline "%?%{yk}%-Lw%?%{wb}%n*%f %t%?(%u)%?%?%{yk}%+Lw%?"


Installing on a new system (e.g. comp9)
=======================================

 * TBD (modulefiles, Snakemake & other requirements)
