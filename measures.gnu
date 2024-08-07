set terminal png \
    size 2048,1536 \
    font ",30" \
    linewidth 2
set output "plots/measures.png"
set grid mxtics mytics
set grid xtics ytics

set xrange [50000:100000]


#  1. Step     2. PotentialEnergyPerParticle     3. Compressibility
plot "data/measures_00.dat" using 1:2 title "Energia potenziale media", \
     "data/measures_00.dat" using 1:3 title "Compressibilità", \
     "data/measures_00.dat" using 1:(($4)**(1.0/3.0)) title "Lato"