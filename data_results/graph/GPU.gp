#********************************************************************
#*	default line styles
#********************************************************************
set style line 10 linetype -1 linewidth 2 linecolor rgb 'black'
set style line 11 linetype  6 linewidth 2 linecolor rgb 'black'
set style line 30 linetype  3 linewidth 2 linecolor rgb 'green'
set style line 20 linetype -1 linewidth 2 linecolor rgb 'blue'
set style line 21 linetype  6 linewidth 2 linecolor rgb 'blue'

#********************************************************************
#*                   Spectra
#********************************************************************
#--------> Time
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (Tesla k40m) Taps=8, Channels=1024"
set output "GPU-t-S.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set key top left
set logscale x
set grid
plot '../k40-simple-spectra.dat' using 1:3 w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-spectra.dat' using 1:2 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-spectra.dat' using 1:3 w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-spectra.dat' using 1:2 w lines ls 21 title 'ldg [FIR]'

#--------> Flops
set output "GPU-f-S.eps"
set ylabel "Flops" offset 0,0 font "Times-Italic,22"
set key top left
set grid
plot '../k40-simple-spectra.dat' using 1:(($5*$2+25600*$1)/($3)) w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-spectra.dat' using 1:5 w lines ls 11 title 'Simple [FIR]'\
,'../k40-simple-spectra.dat' using 1:(25600*$1/($3-$2)) w lines ls 30 title 'cuFFT'\
,'../k40-ldg-spectra.dat' using 1:(($5*$2+25600*$1)/($3)) w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-spectra.dat' using 1:5 w lines ls 21 title 'ldg [FIR]'

#--------> bandwidth
#set format y "%.1e";
set output "GPU-b-S.eps"
set ylabel "Bandwidth" offset 0,0 font "Times-Italic,22"
set key top left
set logscale x
set grid
plot '../k40-simple-spectra.dat' using 1:4 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-spectra.dat' using 1:4 w lines ls 21 title 'ldg [FIR]'


#********************************************************************
#*                   Taps
#********************************************************************
#--------> Time
set output "GPU-t-T.eps"
set title "GPU (Tesla k40m) Channels=1024, Spectra=240000"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("5" 5, "8" 8, "16" 16, "32" 32, "64" 64)
set key top left
set logscale x
set grid
plot '../k40-simple-taps.dat' using 8:3 w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-taps.dat' using 8:2 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-taps.dat' using 8:3 w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-taps.dat' using 8:2 w lines ls 21 title 'ldg [FIR]'

#--------> Flops
set output "GPU-f-T.eps"
set ylabel "Flops" offset 0,0 font "Times-Italic,22"
set key top left
set logscale x
set grid
plot '../k40-simple-taps.dat' using 8:(6.144e+9/($3-$2)) w lines ls 30 title 'cuFFT'\
,'../k40-simple-taps.dat' using 8:((6.144e+9*($3-$2) + $5*$2)/($3)) w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-taps.dat' using 8:5 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-taps.dat' using 8:((6.144e+9*($3-$2) + $5*$2)/($3)) w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-taps.dat' using 8:5 w lines ls 21 title 'ldg [FIR]'

#--------> bandwidth
set output "GPU-b-T.eps"
set ylabel "Bandwidth" offset 0,0 font "Times-Italic,22"
set key top left
set logscale x
set grid
plot '../k40-simple-taps.dat' using 8:4 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-taps.dat' using 8:4 w lines ls 21 title 'ldg [FIR]'

#********************************************************************
#*                   Channels
#********************************************************************
#--------> Time
set title "GPU (Tesla k40m) Taps=8, Spectra=120000"
set output "GPU-t-C.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("256" 256, "512" 512, "1024" 1024, "2048" 2048)
set key top left
set logscale x
set grid
plot '../k40-simple-channels.dat' using 11:3 w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-channels.dat' using 11:2 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:3 w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-channels.dat' using 11:2 w lines ls 21 title 'ldg [FIR]'

#--------> Flops
set output "GPU-f-C.eps"
set ylabel "Flops" offset 0,0 font "Times-Italic,22"
set key top left
set logscale x
set grid
plot '../k40-simple-channels.dat' using 11:($11*3e5*log($11)/($3-$2)/log(2)) w lines ls 30 title 'cuFFT'\
,'../k40-simple-channels.dat' using 11:($5*$2 + $11*2.5*120000*0.3*log($11)*($3-$2))/($3) w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-channels.dat' using 11:5 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:($5*$2 + $11*2.5*120000*0.3*log($11)*($3-$2))/($3) w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-channels.dat' using 11:5 w lines ls 21 title 'ldg [FIR]'

#--------> bandwidth
set output "GPU-b-C.eps"
set ylabel "Bandwidth" offset 0,0 font "Times-Italic,22"
set key top left
set logscale x
set grid
plot '../k40-simple-channels.dat' using 11:4 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:4 w lines ls 21 title 'ldg [FIR]'


#********************************************************************
#*                   ILP - blocks
#********************************************************************
set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (Tesla k40m) Taps=8, Channels=1024, Spectra=240000"
set output "GPU-t-Bl.eps"
set ylabel "Time" offset 0,0 font "Times-Italic,22"
set xlabel "nBlocks" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("1" 1, "2" 2, "3" 3, "4" 4,  "5" 5,  "6" 6, "7" 7,  "8" 8)
set key top left
#set logscale x
unset logscale x
set grid
plot '../k40-ldg-blocks.dat' using 12:3 w lines ls 21 title 'ldg [FIR]'
#, '../k40-ldg-blocks-60.dat' using 12:3 w lines lt -1 lw 3 lc rgb 'black' title 'ldg-60 [FIR]'\

reset;
