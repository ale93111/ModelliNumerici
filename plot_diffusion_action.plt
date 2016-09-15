#!/gnuplot

FILE_IN = 'diffusion_action.txt'
set terminal png enhanced

#pngcairo
#gnuplotting.org

set output "henon_diffusion_action.png"


set xlabel 'Action'
set ylabel 'Diffusion coeff.'

#set logscale y
#set logscale x
#set yrange [1e-5:2]
#set format y "10^{%L}"

plot FILE_IN using 1:2 with points pointtype 7 linecolor rgb 'blue' title 'Diffusion coeff. VS Action' 
