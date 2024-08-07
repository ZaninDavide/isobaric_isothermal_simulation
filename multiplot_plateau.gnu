set terminal png \
    size 1600*2,1200*3 \
    enhanced font ",28" \
    linewidth 4
set output "plots/multiplot_plateau.png"
set grid mxtics mytics
set grid xtics ytics

set grid

# set multiplot layout 5,3
set multiplot layout 3,2

array Titles[14] = [ "ρ=0.001", "ρ=0.010", "ρ=0.100", "ρ=0.200", "ρ=0.250", "ρ=0.300", "ρ=0.400", "ρ=0.500", "ρ=0.600", "ρ=0.700", "ρ=0.750", "ρ=0.800", "ρ=1.000", "ρ=1.200" ]

# array indices[6] = [0,1,2,3,5,6]
array indices[6] = [7,8,9,11,12,13]

do for [j=0:5] {
    i = indices[j+1]

    points = 1000
    pointsE = 1000
    pointsC = 1250

    set title Titles[i+1]
    set xlabel "Punti per blocco (B)"
    set ylabel "s_B/N_B"
    set xrange [0:points]
    set key horizontal top

    file = sprintf("data/plateau_%02d.dat", i)

    errorE(B) = exp(B/pointsE) - 1 + 0.1
    errorC(B) = exp(2*B/pointsC) - 1 + 0.1
    
    f(x) = c*(1 - exp(-a*x))
    a = 0.05
    c = 0.01
    fit f(x) file every ::0::(pointsE - 1) using 1:2 via a, c
    ff(x) = c - d*exp(-a*x)
    d = c
    fit ff(x) file every ::0::(pointsE - 1) using 1:2:(0):(errorE($1)) xyerrors via a, c, d
    # title_f = sprintf('Fit Energia: %f - %f⋅exp(-B/%f))', c, d, 1/a)

    g(x) = cc*(1 - exp(-aa*x))
    aa = 0.05
    cc = 0.01
    fit g(x) file every ::0::(pointsC - 1)  using 1:3 via aa, cc 
    gg(x) = cc - dd*exp(-aa*x)
    dd = cc
    fit gg(x) file every ::0::(pointsC - 1) using 1:3:(0):(errorC($1)) xyerrors via aa, cc, dd 
    # title_g = sprintf('Fit Compressibilità: %f - %f⋅exp(-B/%f))', cc, dd, 1/aa)

    plot file using 1:2 title "Std Media Energia" lc "blue", \
         file using 1:3 title "Std Media Compressibilità" lc "red", \
         ff(x) with lines title sprintf("Energia: τ = %d, σ = %.4f", floor(1/a) + 1, c) lc "blue", \
         gg(x) with lines title sprintf("Compressibilità: τ = %d, σ = %.4f", floor(1/aa) + 1, cc) lc "red"
}

unset multiplot