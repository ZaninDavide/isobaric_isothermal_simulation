set terminal png \
    size 2048,1536 \
    font ",36" \
    linewidth 2
set output "plots/density_histogram.png"
set grid mxtics mytics
set grid xtics ytics
set xtics 0.1

set title "Distribuzione di probabilità della densità"
set xlabel "Densità ρ"
set ylabel "Frequenza"

# Histogram
MAX_DENSITY_BIN = 1.0
DENSITY_BINS = 300.0
width = MAX_DENSITY_BIN / DENSITY_BINS
set xrange [0:MAX_DENSITY_BIN]
set boxwidth width # 1 / BINS
set style fill solid 0.20 #fillstyle
set tics out nomirror
plot "data/density_histogram.dat" using ($1 + width/2.0):2 with boxes notitle