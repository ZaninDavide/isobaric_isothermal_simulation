set terminal png \
    size 2048,1536 \
    font ",30" \
    linewidth 2
set output "plots/plateau.png"
set grid mxtics mytics
set grid xtics ytics

set xlabel "Punti per blocco (B)"
set ylabel "Deviazione standard delle medie sui blocchi"

plot "data/plateau_00.dat" using 1:2 title "Std Avg Energy", \
     "data/plateau_00.dat" using 1:3 title "Std Avg Compressibility"