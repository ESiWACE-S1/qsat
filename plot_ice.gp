set terminal pdf
set output 'qsat_ice.pdf'
set title "Comparision of different qsat approximations over ice"
set xlabel 'Temperature [degC]'
set ylabel 'difference wrt. Wagner Pruss [%]'
plot 'ice.dat' u 1:(100 - 100*$3/$2) w l t 'Imp.Mag.', \
     'ice.dat' u 1:(100 - 100*$4/$2) w l t 'Huang', \
     'ice.dat' u 1:(100 - 100*$5/$2) w l t 'DALES'
