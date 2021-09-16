src=$HOME/GRACEFO/src
bin=$HOME/GRACEFO/bin

#gfortran -o $bin/Read_GRFO_L1B.e   $src/Read_GRFO_L1B.f90    $HOME/ftools/src/SBRT.interp1.f90  $HOME/ftools/src/SBRT.SPLINE.f90
#rm *.mod

#gfortran -o $bin/Compute_LGD_AOD1B.e  $src/Compute_LGD_AOD1B.f90
#gfortran -o $bin/Compute_LGD.e     $src/Compute_LGD.f90
#gfortran -o $bin/Compute_LOSACC.e  $src/Compute_LOSACC.f90   $HOME/ftools/src/SBRT.gradient.f90 $HOME/ftools/src/SBRT.filter_butterworth.f90


#gfortran -o $bin/rr2lgd.e $src/rr2lgd.f90 $HOME/ftools/src/FFTCC.f $HOME/ftools/src/SBRT.GetFrequency.f90 $HOME/ftools/src/SBRT.gradient.f90 $HOME/ftools/src/SBRT.filter_butterworth.f90


#gfortran -o $bin/NEQ_Gnr_cmat_sph.e  $src/NEQ_Gnr_cmat_sph.f90
#gfortran -o $bin/NEQ_Gnr_cvec_sph.e  $src/NEQ_Gnr_cvec_sph.f90

gfortran -o $bin/NEQ_Gnr_cmat_blk.e  $src/NEQ_Gnr_cmat_blk.f90
gfortran -o $bin/NEQ_Gnr_cvec_blk.e  $src/NEQ_Gnr_cvec_blk.f90


#gfortran -o $bin/NEQ_Comb.e  $src/NEQ_Comb.f90
#gfortran -o $bin/NEQ_Solv.e        $src/NEQ_Solv.f90       $HOME/sht/src/SBRT.rotsph.f90 -llapack -lblas
#gfortran -o $bin/NEQ_Solv_glb.e   $src/NEQ_Solv_glb.f90  -llapack -lblas -fopenmp


#gfortran -o $bin/Compute_AliasCorr.e $src/Compute_AliasCorr.f90   $HOME/ftools/src/SBRT.filter_butterworth.f90 $HOME/ftools/src/SBRT.interp1.f90 $HOME/ftools/src/SBRT.SPLINE.f90


#gfortran -o $bin/Slp_Inv.e $src/Slp_Inv.f90  -llapack -lblas

#gfortran -o $bin/Anal_Grtk.e $src/Anal_Grtk.f90 

#gfortran -o $bin/scap.e $src/scap.f90 $HOME/sht/src/SBRT.rotsph.f90


#gfortran -o $bin/Compute_LGD_otide.e $src/Compute_LGD_otide.f90 $src/AA.f


#gfortran -o $bin/NEQ_Gnr_cmat_otide.e  $src/NEQ_Gnr_cmat_otide.f90  $src/AA.f
#gfortran -o $bin/NEQ_Gnr_cvec_otide.e  $src/NEQ_Gnr_cvec_otide.f90  $src/AA.f
#gfortran -o $bin/NEQ_Comb_otide.e      $src/NEQ_Comb_otide.f90
#gfortran -o $bin/NEQ_Solv_otide.e      $src/NEQ_Solv_otide.f90      -llapack -lblas


#gfortran -o $bin/lgd_orbcmp.e $src/lgd_orbcmp.f90 $HOME/ftools/misc/FFTCC.f $HOME/ftools/misc/SBRT.GetFrequency.f90 $HOME/ftools/misc/SBRT.gradient.f90 $HOME/ftools/misc/SBRT.filter_butterworth.f90



# gfortran -o $bin/Find_Scaps.e             $src/Find_Scaps.f90  $HOME/sht/src/SBRT.rotsph.f90
# gfortran -o $bin/NEQ_Gnr_cmat_sph_ext.e   $src/NEQ_Gnr_cmat_sph_ext.f90      -llapack -lblas
