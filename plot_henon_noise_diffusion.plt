#!/gnuplot

FILE_IN = 'henon_noise_diffusion.txt'
set terminal png enhanced

set output "henon_noise_diffusion.png"

set xlabel 'Action'
set ylabel 'Occurrences'
set xrange [0:0.18]

binwidth=0.001
set boxwidth binwidth

bin(x,width)=width*floor(x/width) + width/2.0

plot 'henon_noise_diffusion.txt' using (bin($1,binwidth)):(1.0) smooth freq with points pointtype 7 linecolor rgb 'blue' title 'initial',\
							  '' using (bin($2,binwidth)):(1.0) smooth freq with points pointtype 7 linecolor rgb 'red' title 'final'
