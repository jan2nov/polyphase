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
,'../k40-simple-channels.dat' using 11:($5*$2 + $11*3e5*log($11)/log(2))/($3) w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-channels.dat' using 11:5 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:($5*$2 + $11*3e5*log($11)/log(2))/($3) w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
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

set term postscript eps enhanced color 'Times-Italic' 18
set title "GPU (Tesla k40m) Taps=8, Channels=1024"
set output "GPU-streams.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set logscale x
set key top left
set grid
plot for [IDX=0:4] '../k40-stream-spectra.dat' i IDX u 1:2 w lines lw 3 title columnheader(1)
reset;



#------------------------------------------------------------------
#-------------------> Serial spectra
set term postscript eps enhanced 'Times-Italic' 16
set output "CPU Serial S.eps"

set size 1.0, 1.0
set origin 0.0, 0.0
set multiplot

#----------------------------- Left -------------------------------------------
set size 0.48,0.93
set origin 0.0,0.0
set ylabel "Time[s]" offset 0,0 font "Times-Italic,20"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set grid
set logscale x
unset key
#set yrange[0:5]
plot '../k40-simple-spectra.dat' using 1:3 w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-spectra.dat' using 1:2 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-spectra.dat' using 1:3 w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-spectra.dat' using 1:2 w lines ls 21 title 'ldg [FIR]'
#------------------------------------------------------------------------

#----------------------------- Right -------------------------------------------
set size 0.48,0.93
set origin 0.50,0.0
set ylabel "GFLOPS" offset 0,0 font "Times-Italic,20"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set ytics ("150" 1.5e10, "200" 2e11, "250" 2.5e11, "300" 3e11, "350" 3.5e11, "400" 4e11)
set logscale x
set grid
unset key
set yrange[1.5e11:4e11]
plot '../k40-simple-spectra.dat' using 1:(($5*$2+25600*$1)/($3)) w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-spectra.dat' using 1:5 w lines ls 11 title 'Simple [FIR]'\
,'../k40-simple-spectra.dat' using 1:(25600*$1/($3-$2)) w lines ls 30 title 'cuFFT'\
,'../k40-ldg-spectra.dat' using 1:(($5*$2+25600*$1)/($3)) w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
,'../k40-ldg-spectra.dat' using 1:5 w lines ls 21 title 'ldg [FIR]'
#------------------------------------------------------------------------


#--------------------------------- Legend -----------------------------------------
set size 1.07,0.165
set origin 0.05,0.875


unset key
unset xlabel
unset ylabel
unset xtics
unset ytics
unset title
unset logscale x
unset zeroaxis
unset border
set xrange[-1:10]
set yrange[-3:1]
set label 1 "Naive" at -0.3,0.1 font "Times-Italic,18"
set label 2 "Sequential" at 2.7,0.1 font "Times-Italic,18"
set label 3 "Intrinsic" at 5.7,0.1 font "Times-Italic,18"
set label 4 "Fastest" at -0.3,-1.9 font "Times-Italic,18"
set label 5 "F. cache" at 2.7,-1.9 font "Times-Italic,18"
#set label 6 "p_{m}(a=1.1)" at 5.7,-1.9 font "Times-Italic,18"
plot (x>1.0 && x<2.0 ? 0:1/0) lw 5 lt 1 lc rgb 'black',\
(x>4.0 && x<5.0 ? 0:1/0) lw 5 lt 2 lc rgb 'black',\
(x>7.0 && x<8.0 ? 0:1/0) lw 5 lt 4 lc rgb 'black',\
(x>1.0 && x<2.0 ? -2:1/0) lw 5 lt 5 lc rgb 'black',\
(x>4.0 && x<5.0 ? -2:1/0) lw 5 lt 7 lc rgb 'black'
#(x>7.0 && x<8.0 ? -2:1/0) lw 5 lt 6 lc rgb 'black'

unset multiplot
reset;
