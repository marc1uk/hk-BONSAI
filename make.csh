#!/bin/csh
#setenv ROOTSYS /usr/local/sklib_g77/root_v5.28.00b
setenv ROOTSYS /cern/root/

if ($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH  `$ROOTSYS/bin/root-config --libdir`:$LD_LIBRARY_PATH
else
setenv LD_LIBRARY_PATH  `$ROOTSYS/bin/root-config --libdir`
endif

set path = ( $ROOTSYS/bin $path )

make  $1
