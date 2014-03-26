#********************************************************************
#*                   Spectra
#********************************************************************
#--------> Time
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Taps=8, Channels=1024"
set output "GPU [t][S].eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set key top left
set logscale x
set grid
plot '../k40-simple-spectra.dat' using 1:3 w lines lt -1 lw 3 lc rgb 'black' title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-spectra.dat' using 1:2 w lines lt 6 lw 3 lc rgb 'black' title 'Simple [FIR]'\
,'../k40-ldg-spectra.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'blue' title 'ldg [FIR]'\
,'../k40-ldg-spectra.dat' using 1:3 w lines lt 6 lw 3 lc rgb 'blue' title 'ldg [PPF=FIR+FFT]'
reset;

#--------> Flops
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Taps=8, Channels=1024"
set output "GPU [f][S].eps"
set ylabel "Flops" offset 0,0 font "Times-Italic,22"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set key top left
set logscale x
set grid
plot '../k40-simple-spectra.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'black' title 'Gflops-Simple [FIR]'\
,'../k40-ldg-spectra.dat' using 1:5 w lines lt 6 lw 3 lc rgb 'green' title 'Gflops-ldg [FIR]'\
,'../k40-restrict-spectra.dat' using 1:5 w lines lt 6 lw 3 lc rgb 'red' title 'Gflops-restrict [FIR]'\
,'../k40-fmaf-spectra.dat' using 1:5 w lines lt 6 lw 3 lc rgb 'blue' title 'Gflops-fmaf [FIR]'\
,'../k40-fmaf-ldg-spectra.dat' using 1:5 w lines lt 6 lw 3 lc rgb 'yellow' title 'Gflops-fmaf-ldg [FIR]'
#,'../k40-simple-spectra.dat' using 1:7 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
reset;

#--------> bandwidth
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Taps=8, Channels=1024"
set output "GPU [b][S].eps"
set ylabel "Bandwidth" offset 0,0 font "Times-Italic,22"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set key top left
set logscale x
set grid
plot '../k40-simple-spectra.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'black' title 'Bandwidth-Simple [FIR]'\
,'../k40-ldg-spectra.dat' using 1:4 w lines lt 6 lw 3 lc rgb 'green' title 'ldg [FIR]'\
,'../k40-restrict-spectra.dat' using 1:4 w lines lt 6 lw 3 lc rgb 'red' title 'restrict [FIR]'\
,'../k40-fmaf-spectra.dat' using 1:4 w lines lt 6 lw 3 lc rgb 'blue' title 'fmaf [FIR]'\
,'../k40-fmaf-ldg-spectra.dat' using 1:4 w lines lt 6 lw 3 lc rgb 'yellow' title 'fmaf-ldg [FIR]'
#,'../k40-simple-spectra.dat' using 1:7 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
#,'../k40-simple-spectra.dat' using 1:4 w lines lt 6 lw 3 lc rgb 'black' title 'Gflops [FFT]'\
#,'../k40-simple-spectra.dat' using 1:4 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
reset;


#********************************************************************
#*                   Taps
#********************************************************************
#--------> Time
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Channels=1024, Spectra=240000"
set output "GPU [t][T].eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("5" 5, "8" 8, "16" 16, "32" 32, "64" 64)
set key top left
set logscale x
set grid
plot '../k40-simple-taps.dat' using 8:3 w lines lt -1 lw 3 lc rgb 'black' title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-taps.dat' using 8:2 w lines lt 6 lw 3 lc rgb 'black' title 'Simple [FIR]'\
,'../k40-ldg-taps.dat' using 8:2 w lines lt -1 lw 3 lc rgb 'blue' title 'ldg [FIR]'\
,'../k40-ldg-taps.dat' using 8:3 w lines lt 6 lw 3 lc rgb 'blue' title 'ldg [PPF=FIR+FFT]'
reset;

#--------> Flops
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Channels=1024, Spectra=240000"
set output "GPU [f][T].eps"
set ylabel "Flops" offset 0,0 font "Times-Italic,22"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("5" 5, "8" 8, "16" 16, "32" 32, "64" 64)
set key top left
set logscale x
set grid
plot '../k40-simple-taps.dat' using 8:5 w lines lt -1 lw 3 lc rgb 'black' title 'Gflops-Simple [FIR]'\
,'../k40-ldg-taps.dat' using 8:5 w lines lt 6 lw 3 lc rgb 'green' title 'Gflops-ldg [FIR]'\
,'../k40-restrict-taps.dat' using 8:5 w lines lt 6 lw 3 lc rgb 'red' title 'Gflops-restrict [FIR]'\
,'../k40-fmaf-taps.dat' using 8:5 w lines lt 6 lw 3 lc rgb 'blue' title 'Gflops-fmaf [FIR]'\
,'../k40-fmaf-ldg-taps.dat' using 8:5 w lines lt 6 lw 3 lc rgb 'yellow' title 'Gflops-fmaf-ldg [FIR]'
#,'../k40-simple-spectra.dat' using 8:7 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
reset;

#--------> bandwidth
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Channels=1024, Spectra=240000"
set output "GPU [b][T].eps"
set ylabel "Bandwidth" offset 0,0 font "Times-Italic,22"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("5" 5, "8" 8, "16" 16, "32" 32, "64" 64)
set key top left
set logscale x
set grid
plot '../k40-simple-taps.dat' using 8:4 w lines lt -1 lw 3 lc rgb 'black' title 'Simple [FIR]'\
,'../k40-ldg-taps.dat' using 8:4 w lines lt 6 lw 3 lc rgb 'green' title 'ldg [FIR]'\
,'../k40-restrict-taps.dat' using 8:4 w lines lt 6 lw 3 lc rgb 'red' title 'restrict [FIR]'\
,'../k40-fmaf-taps.dat' using 8:4 w lines lt 6 lw 3 lc rgb 'blue' title 'fmaf [FIR]'\
,'../k40-fmaf-ldg-taps.dat' using 8:4 w lines lt 6 lw 3 lc rgb 'yellow' title 'fmaf-ldg [FIR]'
#,'../k40-simple-spectra.dat' using 1:7 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
#,'../k40-simple-spectra.dat' using 1:4 w lines lt 6 lw 3 lc rgb 'black' title 'Gflops [FFT]'\
#,'../k40-simple-spectra.dat' using 1:4 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
reset;


#********************************************************************
#*                   Channels
#********************************************************************
#--------> Time
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Taps=8, Spectra=120000"
set output "GPU [t][C].eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("256" 256, "512" 512, "1024" 1024, "2048" 2048)
set key top left
set logscale x
set grid
plot '../k40-simple-channels.dat' using 11:3 w lines lt -1 lw 3 lc rgb 'black' title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-channels.dat' using 11:2 w lines lt 6 lw 3 lc rgb 'black' title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:2 w lines lt -1 lw 3 lc rgb 'blue' title 'ldg [FIR]'\
,'../k40-ldg-channels.dat' using 11:3 w lines lt 6 lw 3 lc rgb 'blue' title 'ldg [PPF=FIR+FFT]'
reset;

#--------> Flops
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) Taps=8, Spectra=120000"
set output "GPU [f][C].eps"
set ylabel "Flops" offset 0,0 font "Times-Italic,22"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("256" 256, "512" 512, "1024" 1024, "2048" 2048)
set key top left
set logscale x
set grid
plot '../k40-simple-channels.dat' using 11:5 w lines lt -1 lw 3 lc rgb 'black' title 'Gflops-Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:5 w lines lt 6 lw 3 lc rgb 'green' title 'Gflops-ldg [FIR]'\
,'../k40-restrict-channels.dat' using 11:5 w lines lt 6 lw 3 lc rgb 'red' title 'Gflops-restrict [FIR]'\
,'../k40-fmaf-channels.dat' using 11:5 w lines lt 6 lw 3 lc rgb 'blue' title 'Gflops-fmaf [FIR]'\
,'../k40-fmaf-ldg-channels.dat' using 11:5 w lines lt 6 lw 3 lc rgb 'yellow' title 'Gflops-fmaf-ldg [FIR]'
#,'../k40-simple-spectra.dat' using 8:7 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
reset;

#--------> bandwidth
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) taps=8, Spectra=120000"
set output "GPU [b][C].eps"
set ylabel "Bandwidth" offset 0,0 font "Times-Italic,22"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("256" 256, "512" 512, "1024" 1024, "2048" 2048)
set key top left
set logscale x
set grid
plot '../k40-simple-channels.dat' using 11:4 w lines lt -1 lw 3 lc rgb 'black' title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:4 w lines lt 6 lw 3 lc rgb 'green' title 'ldg [FIR]'\
,'../k40-restrict-channels.dat' using 11:4 w lines lt 6 lw 3 lc rgb 'red' title 'restrict [FIR]'\
,'../k40-fmaf-channels.dat' using 11:4 w lines lt 6 lw 3 lc rgb 'blue' title 'fmaf [FIR]'\
,'../k40-fmaf-ldg-channels.dat' using 11:4 w lines lt 6 lw 3 lc rgb 'yellow' title 'fmaf-ldg [FIR]'
#,'../k40-simple-spectra.dat' using 1:7 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
#,'../k40-simple-spectra.dat' using 1:4 w lines lt 6 lw 3 lc rgb 'black' title 'Gflops [FFT]'\
#,'../k40-simple-spectra.dat' using 1:4 w lines lt 3 lw 3 lc rgb 'black' title 'Gflops [PPF=FIR+FFT]'
reset;


#********************************************************************
#*                   ILP - blocks
#********************************************************************
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (k40) taps=8, Channels = 1024"
set output "GPU [t][Bl].eps"
set ylabel "Time" offset 0,0 font "Times-Italic,22"
set xlabel "nBlocks" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("1" 1, "2" 2, "3" 3, "4" 4,  "5" 5,  "6" 6, "7" 7,  "8" 8)
set key top left
#set logscale x
set grid
plot '../k40-ldg-blocks.dat' using 12:3 w lines lt -1 lw 3 lc rgb 'black' title 'ldg-240 [FIR]'\
, '../k40-ldg-blocks-60.dat' using 12:3 w lines lt -1 lw 3 lc rgb 'black' title 'ldg-60 [FIR]'\

reset;
