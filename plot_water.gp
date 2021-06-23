set terminal pdf
set output 'qsat_water.pdf'
set title "Comparision of different qsat approximations over water"
set xlabel 'Temperature [degC]'
set ylabel 'difference wrt. Wagner Pruss [%]'
plot 'water.dat' u 1:(100 - 100*$3/$2) w l t 'Imp.Mag.', \
     'water.dat' u 1:(100 - 100*$4/$2) w l t 'Huang', \
     'water.dat' u 1:(100 - 100*$5/$2) w l t 'DALES'
