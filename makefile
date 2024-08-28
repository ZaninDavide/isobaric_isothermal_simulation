main.x: main.c
	gcc -o main.x main.c -lm -O2

build: main.x

run: 
	./main.x

# Define the variable for all .dat files in the data directory
DAT_FILES = $(wildcard data/*.dat)

# Define the rule for generating PNG files from Gnuplot scripts
plots/%.png: %.gnu $(DAT_FILES)
	gnuplot $<

# Define the plots target that depends on all the specific PNG files
plots: plots/histogram.png plots/measures.png plots/plateau.png plots/scatter.png plots/density_histogram.png