#!/bin/sh

# impact_atmosphere_posted.sh
#  purpose is to loop over Mi or over pCO2
#  source impact_atmosphere_posted.sh

#  these are initial inventories expressed in bars
   pCO=0.2
   pCO2=10
   pH2=5
   pH2O=500
   pCH4=0
   pN2=2
   pNH3=0
   minpH2O=10   #  smallest amount of H2O for a global event

#   planet=Super_Earth   # 12 characters max
#   planet=Mars
    planet=Earth
#   Mplanet=1   # multiples Earth
#   Aplanet=1  # Astronomical units

#   buffer=1     # no buffer, decrement atm Oxygen with Fe to FeO
#   buffer=5    # IW
#   buffer=2    # QFM
#   buffer=3     # QFM-1
#   buffer=4     # QFM-2
#   buffer=6    # QFI

   eta=0.5   # fraction of impact energy spent evaporating water and heating atmosphere
   buff="buffer"

   mkdir $planet            #
   mkdir $planet/"pCO2=${pCO2}"   # and endure taunts if it already exists


# I use this here when looping over Mi
   for Mi in 22 22.2 22.4 22.6 22.8 23 23.2 23.4 23.6 23.8 24 24.2 24.4 24.6 24.8 25 25.2
   do
    mkdir $planet/"pCO2=${pCO2}"/"Mi=${Mi}"   # and endure taunts if it already exists
# I use this here when looping over buffer or pCO2 or...
    for buffer in  1  # 2 3 4 5 6
    do

# this next bit is to put nice labels on files
     if [ $buffer -eq 5 ]
     then
       fluff="IW"
     elif [ $buffer -eq 6 ]
     then
       fluff="QFI"
     elif [ $buffer -eq 2 ]
    then
       fluff="QFM"
     elif [ $buffer -eq 3 ]
     then
       fluff="QFM1"
     elif [ $buffer -eq 4 ]
     then
       fluff="QFM2"
     elif [ $buffer -eq 1 ]
     then
       fluff="nobuffer"
     else
       fluff="zero"
     fi

#  put parameters into file ocean.para that is read in by fortran program
     echo ${buffer} ${Mi} ${pCO} ${pCO2} ${pH2} ${pH2O} ${pCH4} ${pN2} ${pNH3} ${planet} ${minpH2O} ${eta} > IW.input
     
     gfortran IW_posted.f -o pomo
     
#  the stuff on the screen goes into a file
     ./pomo    > "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.screen"

# outputs are stored in labeled files
     cp IW.out "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.out"
     cp IW.outcolumns "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.outcolumns"
     cp IW.pressure "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.pressure"
     cp IW.drypressure "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.drypressure"
     cp IW.column "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.column"
     cp IW.650columns "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.650columns"  # this is for equilibrating at 650 K

# put the files in subdirectories
#
     mv "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.screen" $planet/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.screen"
     mv "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.out" $planet/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.out"
     mv "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.outcolumns" $planet/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.outcolumns"
     mv "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.650columns" $planet/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.650columns"
     mv "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.pressure" $planet/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.pressure"
     mv "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.drypressure" $planet/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.drypressure"
     mv "IW_Mi=${Mi}_pCO2${pCO2}_pH2O${pH2O}_${fluff}.column" $planet/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.column"

#  create plot file from many terses
     cat < IW.pressure >> terses  # this is the way to append lines I think?
     cat < IW.drypressure >> curses  # this is the way to append lines I think?
     cat < IW.column >> verses  # this is the way to append lines I think?
     cat < IW.650columns >> hearses  # this is the way to append lines I think?

    done  # buffer loop
   done  # Mi loop

# create prettier file
# it takes the first line - the headers - and puts it in new file
# this strips redundant headers
   head -1 terses > plot_file_p
   fgrep + terses >> plot_file_p   # what the fuck is this?

   head -1 verses > plot_file_N
   fgrep + verses >> plot_file_N

   head -1 curses > plot_file_dryP
   fgrep + curses >> plot_file_dryP

   head -1 hearses > plot_file_N650
   fgrep + hearses >> plot_file_N650

   cp plot_file_p "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_p"
   cp plot_file_N "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_N"
   cp plot_file_dryP "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_dryP"
   cp plot_file_N650 "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_N650"
   mv "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_p" $planet/"pCO2=${pCO2}"/"pCO2=${pCO2}_pH2O=${pH2O}.plot_file_p"
   mv "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_N" $planet/"pCO2=${pCO2}"/"pCO2=${pCO2}_pH2O=${pH2O}.plot_file_N"
   mv "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_dryP" $planet/"pCO2=${pCO2}"/"pCO2=${pCO2}_pH2O=${pH2O}.plot_file_dryP"
   mv "IW_pCO2${pCO2}_pH2O${pH2O}.plot_file_N650" $planet/"pCO2=${pCO2}"/"pCO2=${pCO2}_pH2O=${pH2O}.plot_file_N650"

# clean up a little
   rm terses
   rm verses
   rm curses
   rm hearses

   mv IW.pressure IW.pre
   mv IW.column IW.col
   mv IW.drypressure IW.dry
   mv IW.650columns IW.650columns

#
# Created by Kevin Zahnle on 6/5/12.
# Copyright 2012 NASA.


