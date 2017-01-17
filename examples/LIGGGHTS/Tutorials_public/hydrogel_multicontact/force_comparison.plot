#!/usr/bin/gnuplot
set terminal pngcairo dashed

set output 'force.png'
set ytics 5
set key left
set xlabel 'Compression level {/Symbol D} [mm]'
set ylabel 'Force F on the top plate [N]'

set style line 1 lt 1 pt 4 lw 2 lc rgb 'red'
set style line 2 lt 3 pt 0 lw 2 lc rgb 'orange'
set style line 3 lt 1 pt 1 lw 2 lc rgb 'blue'
set style line 4 lt 2 pt 2 lw 2 lc rgb 'cyan'

set xrange [0:25.001]
set yrange [0:16]
plot 'post/piston_force_multicontact.txt' u (( $1- 60000)/30000*25):4 w l ls 1 t 'multicontact contraction', \
     'post/piston_force_multicontact.txt' u ((-$1+120000)/30000*25):4 w l ls 2 t 'multicontact expansion', \
     'post/piston_force_default.txt'      u (( $1- 60000)/30000*25):4 w l ls 3 t 'default contraction', \
     'post/piston_force_default.txt'      u ((-$1+120000)/30000*25):4 w l ls 4 t 'default expansion'
