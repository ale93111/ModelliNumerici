#!/gnuplot

FILE_IN = 'out.txt'
set terminal png enhanced
#pngcairo
#gnuplotting.org

set output "henon_energy_histogram.png"

set style data histogram

set xlabel 'Energy'
set ylabel 'Occorrenze'
#set logscale y
#set yrange [1e-5:2]
#set format y "10^{%L}"
plot FILE_IN using 1:2 with points pointtype 7 linecolor rgb 'blue' title 'Inizio' , '' using 1:3 with points pointtype 7 linecolor rgb 'red' title 'Fine'
#plot for [ROW=1:2] FILE_IN using ROW