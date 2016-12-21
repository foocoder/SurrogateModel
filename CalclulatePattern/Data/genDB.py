#!/usr/bin/env python
# -*- coding: utf-8 -*-
## ---- Program Info Start----
#FileName:  genDB.py
#
#Author:    Fuchen Duan
#
#Email:     slow295185031@gmail.com
#
#CreatedAt: 2016-12-09 15:36:23
## ---- Program Info End  ----

import random
import sys

def genRandomDB( iRow, iColumn, dDensity, strFileName ):
    random.seed();
    DB = '\n'.join([ ' '.join([ '0' if random.uniform(0.0,1.0) > dDensity else '1' for j in xrange(iColumn) ]) for i in xrange(iRow) ]);
    with open( strFileName, 'wb' ) as f:
        f.write( DB );

if __name__ == '__main__':
    genRandomDB( int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4] );

