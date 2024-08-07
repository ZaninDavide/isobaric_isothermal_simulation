set terminal png \
    size 2048*3,1536*5 \
    font ",45" \
    linewidth 2
set output "plots/multiplot_measures.png"
set grid mxtics mytics
set grid xtics ytics

set multiplot layout 5,3
# set multiplot layout 3,2

array Titles[14] = [ "ρ=0.001", "ρ=0.010", "ρ=0.100", "ρ=0.200", "ρ=0.250", "ρ=0.300", "ρ=0.400", "ρ=0.500", "ρ=0.600", "ρ=0.700", "ρ=0.750", "ρ=0.800", "ρ=1.000", "ρ=1.200" ]

do for [i=0:13] {
    set title titles(i)
    set xlabel "Step"
    set ylabel "Osservabile"

    file = sprintf("data/measures_%02d.dat", i)
    set key horizontal top

    plot file using 1:2 title "Energia" lc "blue", \
         file using 1:3 title "Compressibilità" lc "red"   
}

unset multiplot