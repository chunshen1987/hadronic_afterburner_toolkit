#!/usr/bin/env bash

collision_system="PbPb2760"
# first print header
echo "# cen dN/deta(ch)  dN/deta_err(ch)" > $collision_system\_charged_hadron_dNdy.dat
echo "# cen v2{2}(ch)  v2{2}_err(ch)" > $collision_system\_charged_hadron_v2\{2\}.dat
echo "# cen v3{2}(ch)  v3{2}_err(ch)" > $collision_system\_charged_hadron_v3\{2\}.dat
echo "# cen v4{2}(ch)  v4{2}_err(ch)" > $collision_system\_charged_hadron_v4\{2\}.dat
echo "# cen v5{2}(ch)  v5{2}_err(ch)" > $collision_system\_charged_hadron_v5\{2\}.dat
echo "# cen v6{2}(ch)  v6{2}_err(ch)" > $collision_system\_charged_hadron_v6\{2\}.dat
echo "# cen dN/dy dN/dy_err (pi+ K+ p Lambda Xi_m Omega phi)" > $collision_system\_pid_dNdy.dat
echo "# cen <pT> <pT>_err (pi+ K+ p Lambda Xi_m Omega phi)" > $collision_system\_pid_meanpT.dat

idx=0
centrality=(2.5 7.5 15 25 35 45)
for icen in 0-5 5-10 10-20 20-30 30-40 40-50
do
        # charged hadron multiplicity
        foldername=`echo $collision_system\_IPG_C$icen\_Tsw_145`
        dNdy=`cat $foldername/charged_hadron_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $dNdy >> $collision_system\_charged_hadron_dNdy.dat

        # charged hadron v2{2}
        v2=`cat $foldername/charged_hadron_integrated_observables.dat | grep "v_2{2}(ALICE)" | cut -f 2,4 -d " "` 
        echo ${centrality[idx]} $v2 >> $collision_system\_charged_hadron_v2\{2\}.dat
        
        # charged hadron v3{2}
        v3=`cat $foldername/charged_hadron_integrated_observables.dat | grep "v_3{2}(ALICE)" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $v3 >> $collision_system\_charged_hadron_v3\{2\}.dat
        
        # charged hadron v4{2}
        v4=`cat $foldername/charged_hadron_integrated_observables.dat | grep "v_4{2}(ALICE)" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $v4 >> $collision_system\_charged_hadron_v4\{2\}.dat
        
        # charged hadron v5{2}
        v5=`cat $foldername/charged_hadron_integrated_observables.dat | grep "v_5{2}(ALICE)" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $v5 >> $collision_system\_charged_hadron_v5\{2\}.dat
        
        # charged hadron v6{2}
        v6=`cat $foldername/charged_hadron_integrated_observables.dat | grep "v_6{2}(ALICE)" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $v6 >> $collision_system\_charged_hadron_v6\{2\}.dat

        # pid dN/dy
        pion=`cat $foldername/pion_p_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        kaon=`cat $foldername/kaon_p_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        proton=`cat $foldername/proton_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        Lambda=`cat $foldername/Lambda_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        Xi_m=`cat $foldername/Xi_m_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        Omega=`cat $foldername/Omega_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        phi=`cat $foldername/phi_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $pion $kaon $proton $Lambda $Xi_m $Omega $phi >> $collision_system\_pid_dNdy.dat
        
        # pid <pT>
        pion=`cat $foldername/pion_p_integrated_observables.dat | grep "<pT>" | cut -f 2,4 -d " "`
        kaon=`cat $foldername/kaon_p_integrated_observables.dat | grep "<pT>" | cut -f 2,4 -d " "`
        proton=`cat $foldername/proton_integrated_observables.dat | grep "<pT>" | cut -f 2,4 -d " "`
        Lambda=`cat $foldername/Lambda_integrated_observables.dat | grep "<pT>" | cut -f 2,4 -d " "`
        Xi_m=`cat $foldername/Xi_m_integrated_observables.dat | grep "<pT>" | cut -f 2,4 -d " "`
        Omega=`cat $foldername/Omega_integrated_observables.dat | grep "<pT>" | cut -f 2,4 -d " "`
        phi=`cat $foldername/phi_integrated_observables.dat | grep "<pT>" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $pion $kaon $proton $Lambda $Xi_m $Omega $phi >> $collision_system\_pid_meanpT.dat

        idx=$((idx+1))
done
