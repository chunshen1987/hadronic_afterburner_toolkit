#!/usr/bin/env bash

folder=$1
exp=$2

filename="charged_hadron_C2{4}_$exp.dat"
filename1="charged_hadron_v2{2}_$exp.dat"
filename13="charged_hadron_v3{2}_$exp.dat"
filename14="charged_hadron_v4{2}_$exp.dat"
filename15="charged_hadron_v5{2}_$exp.dat"
filename2="charged_hadron_v2{4}_$exp.dat"
filename3="pid_dNdy.dat"
filename4="pid_meanpT.dat"
filename5="charged_hadron_v24_over_v22_$exp.dat"
filename51="charged_hadron_v34_over_v32_$exp.dat"
filename52="charged_hadron_v26_over_v24_$exp.dat"
filename53="symmetric_cumulants_$exp.dat"

idx=0
centrality=(2.5 7.5 15 25 35 45 55 65 75 85 95)
(
    cd $folder
    # first print header
    echo "# cen dN/deta(ch)  C2{4} C2{4}_err" > $filename
    echo "# cen dN/deta(ch)  v2{2} v2{2}_err" > $filename1
    echo "# cen dN/deta(ch)  v3{2} v3{2}_err" > $filename13
    echo "# cen dN/deta(ch)  v4{2} v4{2}_err" > $filename14
    echo "# cen dN/deta(ch)  v5{2} v5{2}_err" > $filename15
    echo "# cen dN/deta(ch)  v2{4} v2{4}_err" > $filename2
    echo "# cen dN/deta(ch)  dN/deta(ch)_err  dN/dy dN/dy_err(pi^+, pi^-, K^+, K^-, p, pbar)" > $filename3
    echo "# cen dN/deta(ch)  <pT>  <pT>_err (charged, pi^+, K^+, p)" > $filename4
    echo "# cen dN/deta(ch)  v2{4}/v2{2}(ch) v2{4}/v2{2}_err(ch) F(v2) F(v2)_err" > $filename5
    echo "# cen dN/deta(ch)  v3{4}/v3{2}(ch) v3{4}/v3{2}_err(ch) F(v3) F(v3)_err" > $filename51
    echo "# cen dN/deta(ch)  v2{6}/v2{4}(ch) v2{6}/v2{4}_err(ch) gamma1 gamma1_err" > $filename52
    if [ $exp = "ALICE" ]
    then
        echo "# cen dN/deta(ch)  SC32  SC32_err  SC42  SC42_err" > $filename53
    fi
    for icen in 00-05 05-10 10-20 20-30 30-40 40-50 50-60 60-70 70-80 80-90 90-100
    do
        echo $icen
        # charged hadron multiplicity
        dNdy=`cat ./$icen/charged_hadron_integrated_observables.dat | grep "dN/dy=" | cut -f 2,4 -d " "`
        if [ $exp = "ALICE" ];
        then
            dNdyCut=`cat ./$icen/charged_hadron_integrated_observables.dat | grep "dN/dy(pT>0.2,|eta|<0.8)=" | cut -f 2,4 -d " "`
        elif [ $exp = "ATLAS" ];
        then
             dNdyCut=`cat ./$icen/charged_hadron_integrated_observables.dat | grep "dN/dy(pT>0.4,|eta|<2.5)=" | cut -f 2,4 -d " "`  
        else
            dNdyCut=`echo $dNdy | cut -f 1 -d " "`
        fi
        C24=`cat ./$icen/charged_hadron_vn4_ALICE.dat | head -n 3 | tail -n 1 | awk {'print $4, $5'}`
        v22=`cat ./$icen/charged_hadron_integrated_observables.dat | sed 's/phenix/PHENIX/g' | grep "$exp" | grep "v_2" | awk {'print $2, $4'}`
        v32=`cat ./$icen/charged_hadron_integrated_observables.dat | sed 's/phenix/PHENIX/g' | grep "$exp" | grep "v_3" | awk {'print $2, $4'}`
        v42=`cat ./$icen/charged_hadron_integrated_observables.dat | sed 's/phenix/PHENIX/g' | grep "$exp" | grep "v_4" | awk {'print $2, $4'}`
        v52=`cat ./$icen/charged_hadron_integrated_observables.dat | sed 's/phenix/PHENIX/g' | grep "$exp" | grep "v_5" | awk {'print $2, $4'}`
        v24=`cat ./$icen/charged_hadron_vn4_ALICE.dat | head -n 3 | tail -n 1 | awk {'print $2, $3'}`
        echo ${centrality[idx]} $dNdyCut $C24 >> $filename
        echo ${centrality[idx]} $dNdyCut $v22 >> $filename1
        echo ${centrality[idx]} $dNdyCut $v32 >> $filename13
        echo ${centrality[idx]} $dNdyCut $v42 >> $filename14
        echo ${centrality[idx]} $dNdyCut $v52 >> $filename15
        echo ${centrality[idx]} $dNdyCut $v24 >> $filename2

        # v2{4}/v2{2}
        ratio=`head -n 4 ./$icen/charged_hadron_vn4_over_vn2_$exp.dat | tail -n 1 | awk {'print $2, $3, $4, $5'}`
        echo ${centrality[idx]} $dNdyCut $ratio >> $filename5

        # v3{4}/v3{2}
        ratio=`head -n 5 ./$icen/charged_hadron_vn4_over_vn2_$exp.dat | tail -n 1 | awk {'print $2, $3, $4, $5'}`
        echo ${centrality[idx]} $dNdyCut $ratio >> $filename51

        # v2{6}/v2{4}
        ratio=`head -n 2 ./$icen/charged_hadron_vn6_over_vn4_$exp.dat | tail -n 1 | awk {'print $2, $3, $4, $5'}`
        echo ${centrality[idx]} $dNdyCut $ratio >> $filename52

        # symmetric cumulants
        if [ $exp = "ALICE" ]
        then
            SC32=`head -n 2 ./$icen/symmetric_cumulant_$exp.dat | tail -n 1 | awk {'print $2, $3'}`
            SC42=`head -n 3 ./$icen/symmetric_cumulant_$exp.dat | tail -n 1 | awk {'print $2, $3'}`
            echo ${centrality[idx]} $dNdyCut $SC32 $SC42 >> $filename53
        fi
        # pid <pT>
        charged=`cat ./$icen/charged_hadron_integrated_observables.dat | grep -m 1 "<pT>" | cut -f 2,4 -d " "`
        pion_p=`cat ./$icen/pion_p_integrated_observables.dat | grep -m 1 "<pT>" | cut -f 2,4 -d " "`
        kaon_p=`cat ./$icen/kaon_p_integrated_observables.dat | grep -m 1 "<pT>" | cut -f 2,4 -d " "`
        pion_m=`cat ./$icen/pion_m_integrated_observables.dat | grep -m 1 "<pT>" | cut -f 2,4 -d " "`
        kaon_m=`cat ./$icen/kaon_m_integrated_observables.dat | grep -m 1 "<pT>" | cut -f 2,4 -d " "`
        proton=`cat ./$icen/proton_integrated_observables.dat | grep -m 1 "<pT>" | cut -f 2,4 -d " "`
        pbar=`cat ./$icen/anti_proton_integrated_observables.dat | grep -m 1 "<pT>" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $dNdyCut $charged $pion_p $pion_m $kaon_p $kaon_m $proton $pbar >> $filename4

        # pid dN/dy
        pion_p=`cat ./$icen/pion_p_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        kaon_p=`cat ./$icen/kaon_p_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        pion_m=`cat ./$icen/pion_m_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        kaon_m=`cat ./$icen/kaon_m_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        proton=`cat ./$icen/proton_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        pbar=`cat ./$icen/anti_proton_integrated_observables.dat | grep "dN/dy" | cut -f 2,4 -d " "`
        echo ${centrality[idx]} $dNdy $pion_p $pion_m $kaon_p $kaon_m $proton $pbar >> $filename3

        idx=$((idx+1))
    done
)
