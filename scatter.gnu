set terminal png \
    size 4096,2048 \
    font ",36" \
    linewidth 2
set output "plots/scatter.png"
set grid mxtics mytics
set grid xtics ytics

set multiplot layout 1,2

L = 1 # we plot (x/L, y/L, z/L)

set title sprintf("Particle grid, L = %f", L)

set xlabel "x"
set ylabel "y"
set xrange [-L:L]
set yrange [-L:L]
unset key
plot "data/scatter.dat" using 1:2 title "Particles" pt 7 ps 3 

set xlabel "x"
set ylabel "z"
set xrange [-L:L]
set yrange [-L:L]
unset key
plot "data/scatter.dat" using 1:3 title "Particles" pt 7 ps 3 