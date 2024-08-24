set terminal png \
    size 2048*2,1536*3 \
    enhanced font ",40" \
    linewidth 4
set output "plots/multiplot_plateau.png"
set grid xtics ytics

set grid

# set multiplot layout 5,3
set multiplot layout 3,2

array Titles[6] = [ "P = 0.001", "P = 0.050", "P = 0.500", "P = 5.0001", "P = 50", "P = 80" ]

array indices[6] = [ 0, 1, 2, 3, 4, 5 ]
# array indices[6] = [7,8,9,11,12,13]

do for [j=0:5] {
    i = indices[j+1]

    points = 20000
    pointsE = 20000
    pointsC = 20000

    set title Titles[i+1]
    set xlabel "Punti per blocco (B)"
    set ylabel "s_B/N_B"
    set xrange [0:points]
    set yrange [0:]
    # set logscale y
    set key horizontal top left Left reverse
    set format y '%.1t×10^{%T}';

    file = sprintf("data/plateau_%02d.dat", i)

    errorE(B) = 1.0 / (  exp(B/pointsE) - 1 + 1.0    )
    errorC(B) = 1.0 / (  exp(2*B/pointsC) - 1 + 1.0  )
    
    f(x) = c*(1 - exp(-abs(a)*x))
    a = 1.0/5000.0
    stats file using 2 nooutput
    c = STATS_max
    fit f(x) file every ::0::(pointsE - 1) using 1:2 via a, c
    ff(x) = c - d*exp(-a*x)
    d = c
    fit ff(x) file every ::0::(pointsE - 1) using 1:2:(0):(errorE($1)) xyerrors via a
    fit ff(x) file every ::0::(pointsE - 1) using 1:2:(0):(errorE($1)) xyerrors via c, d
    # title_f = sprintf('Fit Energia: %f - %f⋅exp(-B/%f))', c, d, 1/a)

    g(x) = cc*(1 - exp(-aa*x))
    aa = 1.0/5000
    stats file using 3 nooutput
    cc = STATS_max
    fit g(x) file every ::0::(pointsC - 1)  using 1:3 via aa, cc 
    gg(x) = cc - dd*exp(-aa*x)
    dd = cc
    fit gg(x) file every ::0::(pointsC - 1) using 1:3:(0):(errorC($1)) xyerrors via aa 
    fit gg(x) file every ::0::(pointsC - 1) using 1:3:(0):(errorC($1)) xyerrors via cc, dd 
    # title_g = sprintf('Fit Compressibilità: %f - %f⋅exp(-B/%f))', cc, dd, 1/aa)

    plot file using 1:2 title "Std Media Energia" lc "blue", \
         file using 1:3 title "Std Media Compressibilità" lc "red", \
         ff(x) with lines title sprintf("Energia: τ = %d, σ = %.1e", floor(1/a) + 1, c) lc "cyan", \
         gg(x) with lines title sprintf("Compressibilità: τ = %d, σ = %.1e", floor(1/aa) + 1, cc) lc "magenta"
}

unset multiplot
