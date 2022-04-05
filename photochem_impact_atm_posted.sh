#!/bin/sh

# photochem_impact_atm_posted.sh
#  source photochem_impact_atm_posted.sh
#  this runs a single case of photochem_implicit_posted.f using IW_posted.f as the input
   eta=0.5

# total initial inventories expressed in bars
   pCO=0
   pCO2=1
   pH2=0
   pH2O=500
   pCH4=0
   pN2=1
   pNH3=0
   minpH2O=10  # smallest amount of H2O for a global event

    planet=Earth

   buffer=1
#   buffer=1 # Fe titration by CO2, H2O
#   buffer=5    # IW
#   buffer=2    # QFM
#   buffer=3     # QFM-1
#   buffer=4     # QFM-2
#   buffer=6    # QFI

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
     fluff="Fe"
   else
     fluff="zero"
   fi

#  stratospheric mixing ratio is a free parameter
#  1e-6 seems dry, 1e-5 seems moist, 1e-4 seems et
   fH2O=1.e-7

   mkdir Chem            #
   mkdir Chem/"pCO2=${pCO2}"   # and endure taunts if it already exists

   Mi=23.4       # impact mass
                  # Vesta is 23.4

   mkdir Chem/"pCO2=${pCO2}"/"Mi=${Mi}"   # and endure taunts if it already exists
   mkdir Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"fH2O=${fH2O}"   # and endure taunts if it already exists

#  put parameters into file IW.input that is read in by fortran program
     echo ${buffer} ${Mi} ${pCO} ${pCO2} ${pH2} ${pH2O} ${pCH4} ${pN2} ${pNH3} ${planet} ${minpH2O} ${eta} > IW.input

     gfortran IW_posted.f -o domo
     
#  the file IW.parameters is written by IW_posted.f for photochem_implicit_posted.f

#  the stuff on the screen goes into a file
     ./domo > "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.screen"
     
# outputs are stored in labeled files
     cp IW.out "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.out"
     cp IW.outcolumns "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.outcolumns"
     cp IW.pressure "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.pressure"
     cp IW.drypressure "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.drypressure"
     cp IW.column "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.column"
     cp IW.650columns "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.650columns"

# put the files in subdirectories
#
     mv "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.screen" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.screen"
     mv "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.out" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.out"
     mv "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.outcolumns" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.columnsout"
     mv "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.650columns" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.650columns"
     mv "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.pressure" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.pressure"
     mv "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.drypressure" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.drypressure"
     mv "IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.column" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"IW_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}.column"

     echo ${fH2O} ${buffer} ${planet} ${Mi} > photochem.input

     gfortran photochem_implicit_posted.f -o tomo

#  the stuff on the screen goes into a file
     ./tomo > "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.screen"

# outputs are stored in labeled files
     cp evolve.out "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.out"
     cp evolve.p "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.p"
     cp evolve.mix "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.mix"

     mv "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.screen" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"fH2O=${fH2O}"/"evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.screen"
     mv "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.out" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"fH2O=${fH2O}"/"evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.out"
     mv "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.p" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"fH2O=${fH2O}"/"evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.p"
     mv "evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.mix" Chem/"pCO2=${pCO2}"/"Mi=${Mi}"/"fH2O=${fH2O}"/"evolve_Mi=${Mi}_pCO2=${pCO2}_pH2O=${pH2O}_${fluff}_f=${fH2O}.mix"

#
# Created by Kevin Zahnle on 6/5/12.
# Copyright 2012 NASA.


