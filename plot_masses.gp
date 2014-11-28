set term postscript color solid enhanced eps "Roman" 28
set output "masses_C4.eps"

set pointsize 1.5 
set style line 1 lw 5
set ylabel 'C4(t)' offset 1.5,0.0
set xlabel 't' offset 0.0,0.3

set grid front
set style fill transparent solid 0.2 noborder
set style line 1 lt -1 lw 1
set key box linestyle 1 t r Left spacing 2.5 width 0.5 samplen 0.9 reverse
set bar 1.5
set size ratio 0.75
set sample 500 

set xrange [10:23]
set yrange [0.01:2]
#set ytics 0.01
#set format y '%.2f'
set xtics 2.0
set logscale y
plot "./plots/C4_p00.dat"  u 1:2:3 t 'p=0'  w yerrorbars pt 4  lw 1.5 lt 1 lc rgb 'red', \
     "./plots/C4_p11.dat"  u 1:2:3 t 'p=1'  w yerrorbars pt 4  lw 1.5 lt 1 lc rgb 'blue', \
     "./plots/C4_p22.dat"  u 1:2:3 t 'p=2'  w yerrorbars pt 4  lw 1.5 lt 1 lc rgb 'green', \
     "./plots/C4_p33.dat"  u 1:2:3 t 'p=3'  w yerrorbars pt 4  lw 1.5 lt 1 lc rgb 'orange', \
     "./plots/C4_p44.dat"  u 1:2:3 t 'p=4'  w yerrorbars pt 4  lw 1.5 lt 1 lc rgb 'black'

