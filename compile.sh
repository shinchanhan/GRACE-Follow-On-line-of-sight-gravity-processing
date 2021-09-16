src=$HOME/GRACEFO/src
bin=$HOME/GRACEFO/bin

gfortran -o $bin/Read_GRFO_L1B.e   $src/Read_GRFO_L1B.f90    $HOME/ftools/src/SBRT.interp1.f90  $HOME/ftools/src/SBRT.SPLINE.f90
rm *.mod

gfortran -o $bin/Compute_LGD.e  $src/Compute_LGD.f90

gfortran -o $bin/rr2lgd.e $src/rr2lgd.f90 $HOME/ftools/src/FFTCC.f $HOME/ftools/src/SBRT.GetFrequency.f90 $HOME/ftools/src/SBRT.gradient.f90 $HOME/ftools/src/SBRT.filter_butterworth.f90

