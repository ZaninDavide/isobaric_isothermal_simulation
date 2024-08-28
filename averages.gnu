set terminal png \
    size 2048*2,1536*2 \
    font ",40" \
    linewidth 2
set output "plots/averages.png"
set grid mxtics mytics
set grid xtics ytics



set pointsize 3
#set xrange [0:0.7] 

set multiplot layout 2,2

unset xlabel
set ylabel "Pressione P"
set xrange [0:1.2]
set yrange [0.001: 100]
set xtics 0.1
set key top left Left reverse
set logscale y
plot "data/averages_NPT1nolog.dat" using 6:7 with points pt 2 lc "cyan" title "NPT, T=1.0, nolog", \
     "data/averages_NPT2nolog.dat" using 6:7 with points pt 2 lc "orange" title "NPT, T=2.0, nolog", \
     "data/averages_NPT1.dat" using 6:7 with points pt 6 lc "blue" title "NPT, T=1.0", \
     "data/averages_NVT1.dat" using 6:($3*$6*$2) with linespoints pt 1 lc "black" title "NVT, T=1.0", \
     "data/averages_NPT2.dat" using 6:7 with points pt 6 lc "red"  title "NPT, T=2.0", \
     "data/averages_NVT2.dat" using 6:($3*$6*$2) with linespoints pt 1 lc "black"  title "NVT, T=2.0", \
     "data/averages.dat"      using 6:7 with linespoints pt 1 lc "green" notitle, \

unset xlabel
set xrange [0.82:40]
set yrange [0.001: 100]
unset logscale
unset ylabel
set key top right Right
set logscale y
set xtics 5.0
plot "data/averages_NPT1nolog.dat" using (1.0/$6):7 with points pt 2 lc "cyan" title "NPT, T=1.0, nolog", \
     "data/averages_NPT2nolog.dat" using (1.0/$6):7 with points pt 2 lc "orange" title "NPT, T=2.0, nolog", \
     "data/averages_NPT1.dat" using (1.0/$6):7 with points pt 6 lc "blue" title "NPT, T=1.0", \
     "data/averages_NVT1.dat" using (1.0/$6):($3*$6*$2) with linespoints pt 1 lc "black" title "NVT, T=1.0", \
     "data/averages_NPT2.dat" using (1.0/$6):7 with points pt 6 lc "red"  title "NPT, T=2.0", \
     "data/averages_NVT2.dat" using (1.0/$6):($3*$6*$2) with linespoints pt 1 lc "black"  title "NVT, T=2.0", \
     "data/averages.dat"      using (1.0/$6):7 with linespoints pt 1 lc "green" notitle, \

set xlabel "Densità ρ"
set ylabel "Compressibilità C=P/ρk_B"
set xrange [0:1.2]
set yrange [:13]
set xtics 0.1
set ytics 0.1
set key top left Left reverse
unset title
plot "data/averages_NPT1nolog.dat" using 6:3 with points pt 2 lc "cyan" title "NPT, T=1.0, nolog", \
     "data/averages_NPT2nolog.dat" using 6:3 with points pt 2 lc "orange" title "NPT, T=2.0, nolog", \
     "data/averages_NPT1.dat" using 6:3 with points pt 6 lc "blue" title "NPT, T=1.0", \
     "data/averages_NVT1.dat" using 6:3 with linespoints pt 1 lc "black" title "NVT, T=1.0", \
     "data/averages_NPT2.dat" using 6:3 with points pt 6 lc "red"  title "NPT, T=2.0", \
     "data/averages_NVT2.dat" using 6:3 with linespoints pt 1 lc "black"  title "NVT, T=2.0", \
     "data/averages.dat"      using 6:3 with linespoints pt 1 lc "green" notitle, \

set xlabel "Volume per particella V/N"
set key top right Right
set xrange [0.82:40]
set yrange [:13]
set ytics 0.1
set xtics 5
unset title
unset ylabel
plot "data/averages_NPT1nolog.dat" using (1.0/$6):3 with points pt 2 lc "cyan" title "NPT, T=1.0, nolog", \
     "data/averages_NPT2nolog.dat" using (1.0/$6):3 with points pt 2 lc "orange" title "NPT, T=2.0, nolog", \
     "data/averages_NPT1.dat" using (1.0/$6):3 with points pt 6 lc "blue" title "NPT, T=1.0", \
     "data/averages_NVT1.dat" using (1.0/$6):3 with linespoints pt 1 lc "black" title "NVT, T=1.0", \
     "data/averages_NPT2.dat" using (1.0/$6):3 with points pt 6 lc "red"  title "NPT, T=2.0", \
     "data/averages_NVT2.dat" using (1.0/$6):3 with linespoints pt 1 lc "black"  title "NVT, T=2.0", \
     "data/averages.dat"      using (1.0/$6):3 with linespoints pt 1 lc "green" notitle