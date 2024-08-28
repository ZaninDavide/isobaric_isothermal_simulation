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

    points = 80000
    pointsE = 80000
    pointsC = 80000
    pointsV = 80000

    set title Titles[i+1]
    set xlabel "Punti per blocco (B)"
    set ylabel "s_B/N_B"
    set xrange [0:points]
    set yrange [0:]
    # set logscale y
    # set key horizontal top left Left reverse
    set key horizontal bottom right Right
    set format y '%.1t×10^{%T}';

    file = sprintf("data/plateau_%02d.dat", i)

    errorE(B) = 0.001 # exp(B/pointsE) - 1 + 1.0  
    errorC(B) = 0.001 # exp(2*B/pointsC) - 1 + 1.0
    errorV(B) = 0.001 # exp(5*B/pointsC) - 1 + 1.0
    
#   f(x) = c*(1 - exp(-abs(a)*x))
#   a = 1.0/5000.0
#   stats file using 2 nooutput
#   c = STATS_max
#   fit f(x) file every ::0::(pointsE - 1) using 1:2 via a, c
#   ff(x) = c - d*exp(-a*x)
#   d = c
#   fit ff(x) file every ::0::(pointsE - 1) using 1:2:(0):(errorE($1)) xyerrors via a
#   fit ff(x) file every ::0::(pointsE - 1) using 1:2:(0):(errorE($1)) xyerrors via c, d
#   # title_f = sprintf('Fit Energia: %f - %f⋅exp(-B/%f))', c, d, 1/a)
#
#   g(x) = cc*(1 - exp(-aa*x))
#   aa = 1.0/5000
#   stats file using 3 nooutput
#   cc = STATS_max
#   fit g(x) file every ::0::(pointsC - 1)  using 1:3 via aa, cc 
#   gg(x) = cc - dd*exp(-aa*x)
#   dd = cc
#   fit gg(x) file every ::0::(pointsC - 1) using 1:3:(0):(errorC($1)) xyerrors via aa 
#   fit gg(x) file every ::0::(pointsC - 1) using 1:3:(0):(errorC($1)) xyerrors via cc, dd 
#   # title_g = sprintf('Fit Compressibilità: %f - %f⋅exp(-B/%f))', cc, dd, 1/aa)

    h(x) = ccc*(1.0 - exp(-aaa*x))
    aaa = 1.0/4000.0
    stats file using 4 nooutput
    ccc = STATS_max
    fit h(x) file every ::0::(pointsV - 1)  using 1:4 via aaa, ccc
    hh(x) = cccc - dddd*exp(-aaaa*x)
    aaaa = aaa
    cccc = ccc
    dddd = ccc
    fit hh(x) file every ::0::(pointsV - 1) using 1:4:(0):(errorV($1)) xyerrors via aaaa 
    fit hh(x) file every ::0::(pointsV - 1) using 1:4:(0):(errorV($1)) xyerrors via aaaa, cccc, dddd 
    # title_g = sprintf('Fit Volume: %f - %f⋅exp(-B/%f))', cc, dd, 1/aa)

#   plot file using 1:2 title "Std Media Energia" lc "blue", \
#        file using 1:3 title "Std Media Compressibilità" lc "red", \
#        file using 1:4 title "Std Media Volume" lc "green", \
#        ff(x) with lines title sprintf("Energia: τ = %d, σ = %.1e", floor(1/a) + 1, c) lc "cyan", \
#        gg(x) with lines title sprintf("Compressibilità: τ = %d, σ = %.1e", floor(1/aa) + 1, cc) lc "magenta", \
#        hh(x) with lines title sprintf("Volume: τ = %d, σ = %.1e", floor(1/aaa) + 1, ccc) lc "black"

    plot file using 1:4 title "Std Media Volume" lc "green", \
         h(x) with lines title sprintf("Volume: τ = %d, σ = %.1e", floor(1/aaa) + 1, ccc) lc "black"
#        hh(x) with lines title sprintf("Volume: τ = %d, σ = %.1e", floor(1/aaa) + 1, ccc) lc "black"

}

unset multiplot
