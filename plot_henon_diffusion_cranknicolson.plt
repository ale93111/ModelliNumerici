#!/gnuplot

FILE_IN = 'henon_diffusion_cranknicolson.txt'
set terminal png enhanced

set output "henon_diffusion_cranknicolson.png"

set xlabel 'Action'
set ylabel 'ro %'
#set logscale y
#set yrange [1e-5:2]
#set format y "10^{%L}"
plot FILE_IN using 1:2 with points pointtype 7 linecolor rgb 'blue' title 'initial' , '' using 1:3 with points pointtype 7 linecolor rgb 'red' title 'final'
