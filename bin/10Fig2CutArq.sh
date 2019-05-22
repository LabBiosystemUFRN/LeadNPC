#!/bin/bash

arqName=$1
workDir="/home/clovis/Dropbox/Chumbo/figuras/Figura2/redes/net/"
curDir=$(pwd)

cd $workDir

#for arqName in $(ls *.txt); do
    echo $arqName
    refName=$(basename $arqName ".txt")
    
    sed  '/N.\tROUNDED_RECTANGLE/d' $arqName > tmp

    let line=$(grep -En "node_a.+node_b.+weight" tmp|cut -f 1 -d ":")-1
    sed -r 's/#//g' tmp > tmp1
    head -$line tmp1 >"../NodesEdges/"${refName}"Nodes.txt"
    let line=${line}+1
    tail -n +$line tmp1 >"../NodesEdges/"${refName}"Edges.txt"

    rm tmp tmp1
    
#done;

cd ${curDir}





