set terminal png \
    size 2048,1536 \
    font ",36" \
    linewidth 2
set output "plots/cumulative_density_histogram.png"
set grid mxtics mytics
set grid xtics ytics
set xtics 0.1

set title "Distribuzione di probabilità della densità in funzione della pressione"
set xlabel "Densità ρ"
set ylabel "Frequenza / frequenza massima"


stats "data/density_histogram0,0550.dat" using 2 nooutput
max00550 = STATS_max
stats "data/density_histogram0,0500.dat" using 2 nooutput
max00500 = STATS_max
stats "data/density_histogram0,0425.dat" using 2 nooutput
max00425 = STATS_max
stats "data/density_histogram0,0350.dat" using 2 nooutput
max00350 = STATS_max
stats "data/density_histogram0,0150.dat" using 2 nooutput
max00150 = STATS_max
stats "data/density_histogram0,0050.dat" using 2 nooutput
max00050 = STATS_max
stats "data/density_histogram0,0010.dat" using 2 nooutput
max00010 = STATS_max

# Histogram
MAX_DENSITY_BIN = 1.0
DENSITY_BINS = 300.0
width = MAX_DENSITY_BIN / DENSITY_BINS
set xrange [0:MAX_DENSITY_BIN]
set yrange [:1.05]
set boxwidth width # 1 / BINS
set style fill solid 0.20 #fillstyle
set tics out nomirror
set ytics 0.1


plot "data/density_histogram0,0500.dat" using 1:($2/max00500) with lines title "T = 1.0, P = 0.0500", \
     "data/density_histogram0,0425.dat" using 1:($2/max00425) with lines title "T = 1.0, P = 0.0425", \
     "data/density_histogram0,0350.dat" using 1:($2/max00350) with lines title "T = 1.0, P = 0.0350", \
     "data/density_histogram0,0010.dat" using 1:($2/max00010) with lines title "T = 1.0, P = 0.0010"

#    "data/density_histogram0,0150.dat" using 1:($2/max00150) with lines title "T = 1.0, P = 0.0150", \
#    "data/density_histogram0,0050.dat" using 1:($2/max00050) with lines title "T = 1.0, P = 0.0050", \