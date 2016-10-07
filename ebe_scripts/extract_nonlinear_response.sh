#!/usr/bin/env bash

collision_system="PbPb2760"
# first print header
echo "# cen rho_422 rho_523 rho_6222 rho_633 chi_422 chi_523 chi_6222 chi_633" > $collision_system\_nonlinear_response_coefficients_ALICE.dat
echo "# cen v4_L v4(Psi2) v5_L v5(Psi23) v6_L v6(Psi2) v6(Psi3)" > $collision_system\_nonlinear_response_vn_decomposition_ALICE.dat
echo "# cen rho_422 rho_523 rho_6222 rho_633 chi_422 chi_523 chi_6222 chi_633" > $collision_system\_nonlinear_response_coefficients_CMS.dat
echo "# cen v4_L v4(Psi2) v5_L v5(Psi23) v6_L v6(Psi2) v6(Psi3)" > $collision_system\_nonlinear_response_vn_decomposition_CMS.dat
echo "# cen rho_422 rho_523 rho_6222 rho_633 chi_422 chi_523 chi_6222 chi_633" > $collision_system\_nonlinear_response_coefficients_ATLAS.dat
echo "# cen v4_L v4(Psi2) v5_L v5(Psi23) v6_L v6(Psi2) v6(Psi3)" > $collision_system\_nonlinear_response_vn_decomposition_ATLAS.dat

idx=0
centrality=(2.5 7.5 15 25 35 45)
for icen in 0-5 5-10 10-20 20-30 30-40 40-50
do
        foldername=`echo $collision_system\_IPG_C$icen\_Tsw_145`
        # ALICE
        # retrive v4, v5, and v6
        v4L=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "v4_L" | cut -f 3-5 -d " "`
        v4Psi2=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "v4(Psi2)" | cut -f 3-5 -d " "`
        v5L=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "v5_L" | cut -f 3-5 -d " "`
        v5Psi23=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "v5(Psi23)" | cut -f 3-5 -d " "`
        v6L=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "v6_L" | cut -f 3-5 -d " "`
        v6Psi2=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "v6(Psi2)" | cut -f 3-5 -d " "`
        v6Psi3=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "v6(Psi3)" | cut -f 3-5 -d " "`
        # retrive rho_xxx and chi_xxx
        rho422=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "rho_422" | cut -f 3-5 -d " "`
        rho523=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "rho_523" | cut -f 3-5 -d " "`
        rho6222=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "rho_6222" | cut -f 3-5 -d " "`
        rho633=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "rho_633" | cut -f 3-5 -d " "`
        chi422=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "chi_422" | cut -f 3-5 -d " "`
        chi523=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "chi_523" | cut -f 3-5 -d " "`
        chi6222=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "chi_6222" | cut -f 3-5 -d " "`
        chi633=`cat $foldername/non_linear_response_coefficients_ALICE.dat | grep "chi_633" | cut -f 3-5 -d " "`
        # print out data
        echo ${centrality[idx]} $v4L $v4Psi2 $v5L $v5Psi23 $v6_L $v6Psi2 $v6Psi3 >> $collision_system\_nonlinear_response_vn_decomposition_ALICE.dat
        echo ${centrality[idx]} $rho422 $rho523 $rho6222 $rho633 $chi422 $chi523 $chi6222 $chi633 >> $collision_system\_nonlinear_response_coefficients_ALICE.dat
        # CMS
        # retrive v4, v5, and v6
        v4L=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "v4_L" | cut -f 3-5 -d " "`
        v4Psi2=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "v4(Psi2)" | cut -f 3-5 -d " "`
        v5L=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "v5_L" | cut -f 3-5 -d " "`
        v5Psi23=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "v5(Psi23)" | cut -f 3-5 -d " "`
        v6L=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "v6_L" | cut -f 3-5 -d " "`
        v6Psi2=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "v6(Psi2)" | cut -f 3-5 -d " "`
        v6Psi3=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "v6(Psi3)" | cut -f 3-5 -d " "`
        # retrive rho_xxx and chi_xxx
        rho422=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "rho_422" | cut -f 3-5 -d " "`
        rho523=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "rho_523" | cut -f 3-5 -d " "`
        rho6222=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "rho_6222" | cut -f 3-5 -d " "`
        rho633=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "rho_633" | cut -f 3-5 -d " "`
        chi422=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "chi_422" | cut -f 3-5 -d " "`
        chi523=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "chi_523" | cut -f 3-5 -d " "`
        chi6222=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "chi_6222" | cut -f 3-5 -d " "`
        chi633=`cat $foldername/non_linear_response_coefficients_CMS.dat | grep "chi_633" | cut -f 3-5 -d " "`
        # print out data
        echo ${centrality[idx]} $v4L $v4Psi2 $v5L $v5Psi23 $v6L $v6Psi2 $v6Psi3 >> $collision_system\_nonlinear_response_vn_decomposition_CMS.dat
        echo ${centrality[idx]} $rho422 $rho523 $rho6222 $rho633 $chi422 $chi523 $chi6222 $chi633 >> $collision_system\_nonlinear_response_coefficients_CMS.dat
        # ATLAS
        # retrive v4, v5, and v6
        v4L=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "v4_L" | cut -f 3-5 -d " "`
        v4Psi2=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "v4(Psi2)" | cut -f 3-5 -d " "`
        v5L=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "v5_L" | cut -f 3-5 -d " "`
        v5Psi23=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "v5(Psi23)" | cut -f 3-5 -d " "`
        v6L=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "v6_L" | cut -f 3-5 -d " "`
        v6Psi2=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "v6(Psi2)" | cut -f 3-5 -d " "`
        v6Psi3=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "v6(Psi3)" | cut -f 3-5 -d " "`
        # retrive rho_xxx and chi_xxx
        rho422=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "rho_422" | cut -f 3-5 -d " "`
        rho523=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "rho_523" | cut -f 3-5 -d " "`
        rho6222=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "rho_6222" | cut -f 3-5 -d " "`
        rho633=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "rho_633" | cut -f 3-5 -d " "`
        chi422=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "chi_422" | cut -f 3-5 -d " "`
        chi523=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "chi_523" | cut -f 3-5 -d " "`
        chi6222=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "chi_6222" | cut -f 3-5 -d " "`
        chi633=`cat $foldername/non_linear_response_coefficients_ATLAS.dat | grep "chi_633" | cut -f 3-5 -d " "`
        # print out data
        echo ${centrality[idx]} $v4L $v4Psi2 $v5L $v5Psi23 $v6L $v6Psi2 $v6Psi3 >> $collision_system\_nonlinear_response_vn_decomposition_ATLAS.dat
        echo ${centrality[idx]} $rho422 $rho523 $rho6222 $rho633 $chi422 $chi523 $chi6222 $chi633 >> $collision_system\_nonlinear_response_coefficients_ATLAS.dat

        idx=$((idx+1))
done
