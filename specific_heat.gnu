set terminal png \
    size 2048,1536 \
    font ",40" \
    linewidth 2
set output "plots/specific_heat.png"
set grid mxtics mytics
set grid xtics ytics
set xtics 0.1
set ytics 2
set yrange [0:32]
set pointsize 3

set xlabel "Densità"
set ylabel "Calore specifico isobaro (c_P/k_B)"
set title "Calore specifico isobaro in funzione della densità"

plot "data/averages_heat_NPT2.dat" using 6:8 with linespoints pt 7 lc "red"  title "NPT, T=2.0", \
     "data/averages_heat_NPT1.dat" using 6:8 with linespoints pt 5 lc "blue" title "NPT, T=1.0", \