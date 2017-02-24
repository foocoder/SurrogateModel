#!/bin/bash


# -*- coding: utf-8 -*-
## ---- Program Info Start----
#FileName:     batch.sh
#
#Author:       Fuchen Duan
#
#Email:        slow295185031@gmail.com
#
#CreatedAt:    2015-11-03 13:19:25
## ---- Program Info End  ----

#Process Dofia Algorithm
for run in {0..0}
do
    num=0;
    savepath=$1$run
    if [ ! -d "./data/$savepath" ];then
        mkdir -p ./data/$savepath;
    fi
    if [ -f "./data/$savepath/runtime.log" ];then
        rm ./data/$savepath/runtime.log;
    fi
    if [ -f "./data/$savepath/HV.dat" ]; then
        rm ./data/$savepath/HV.dat;
    fi
    if [ -f "./data/logs" ];then
        rm ./data/logs;
    fi

    for i in `ls ../databaseAllCat/itemID:*:bin`
    do
        make clean> /dev/null 2>&1
        filename=`echo $i|grep -o -E "[0-9]+"`
        fileinfo=`cat ../databaseAllCat/itemID:$filename:fileinfo`
        filerow=`echo $fileinfo|awk '{print $1}'`
        filecol=`echo $fileinfo|awk '{print $2}'`
        let num=num+1
        echo $run,$num,$filename,$filerow
        make CFLAGS="-DN=$filerow -DM=$filecol"> /dev/null 2>&1
        ./a.out ../databaseAllCat ./data/$savepath itemID:$filename:bin $filename | tee -a ./data/$savepath/runtime.log
        echo -n $filename >> ./data/$savepath/HV.dat
        echo -n -e "\t" >> ./data/$savepath/HV.dat
        ../QHV/project < ./data/$savepath/${filename}.fit >> ./data/$savepath/HV.dat;
        if [ $? -ne 0 ]; then
            echo $run : $filename >> ./data/logs
        fi
    done
done
