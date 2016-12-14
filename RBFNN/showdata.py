#!/usr/bin/env python
# -*- coding: utf-8 -*-
## ---- Program Info Start----
#FileName:  showdata.py
#
#Author:    Fuchen Duan
#
#Email:     slow295185031@gmail.com
#
#CreatedAt: 2016-12-02 19:16:35
## ---- Program Info End  ----

from matplotlib import pyplot as plt

def drawPic( filename ):

    with open( './testdata', 'rb' ) as fdata:
        x         = [];
        y_real    = [];
        y_predict = [];
        for eachline in fdata.readlines():
            elems = eachline.strip().split();
            x.append( float(elems[0]) );
            y_predict.append( float(elems[1]) );
            y_real.append( float(elems[2]) );
        plt.figure(figsize=(12, 8))
        plt.plot(x, y_real, 'k-')

        # plot learned model
        plt.plot(x, y_predict, 'r-', linewidth=2)

    with open( './traindata', 'rb' ) as fdata:
        x = [];
        y = [];
        for eachline in fdata.readlines():
            elems = eachline.strip().split();
            x.append( float(elems[0]) );
            y.append( float(elems[1]) );
        plt.plot( x, y, 'gs' );

    with open( './paradata', 'rb' ) as fdata:
        x = [];
        for eachline in fdata.readlines():
            elems = eachline.strip().split();
            x.append( float(elems[1]) );
        plt.plot( x, [0 for i in range(len(x))], 'r^' );
    # with open( './errdata', 'rb' ) as fdata:
        # x = [];
        # for eachline in fdata.readlines():
            # elems = eachline.strip();
            # x.append( float(elems) );
        # plt.figure(figsize=(12, 8))
        # plt.plot( range(len(x)), x, 'g-' );

    print "Reading Data Finish..."

    plt.show();

if __name__ == "__main__":
    drawPic('./traindata');
