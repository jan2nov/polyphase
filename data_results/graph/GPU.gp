#********************************************************************
#*	default line styles
#********************************************************************
set style line 10 linetype -1 linewidth 2 linecolor rgb 'black'
set style line 11 linetype  6 linewidth 2 linecolor rgb 'black'
set style line 30 linetype -1 linewidth 2 linecolor rgb 'green'
set style line 31 linetype  6 linewidth 2 linecolor rgb 'green'
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
,'../k40-simple-taps.dat' using 8:((6.144e+9 + $5*$2)/($3)) w lines ls 10 title 'Simple [PPF=FIR+FFT]'\
,'../k40-simple-taps.dat' using 8:5 w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-taps.dat' using 8:((6.144e+9 + $5*$2)/($3)) w lines ls 20 title 'ldg [PPF=FIR+FFT]'\
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

#------------------------------------------------------------------
#-------------------> Serial spectra
set term postscript eps enhanced 'Times-Italic' 16
set output "GPU-spectra-tf.eps"

set size 1.0, 1.0
set origin 0.0, 0.0
set multiplot

#----------------------------- Left -------------------------------------------
set size 0.48,0.93
set origin 0.0,0.0
unset title
set ylabel "Time[s]" offset 0,0 font "Times-Italic,20"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set xrange[15000:120000]
set grid
set logscale x
unset key
#set yrange[0:5]
plot '../gtx580-spectra.dat' using 1:2 w lines ls 10 title 'FIR'\
,'../gtx580-spectra.dat' using 1:($3-$2) w lines ls 11 title 'FFT'\
,'../k40-ldg-spectra.dat' using 1:2 w lines ls 20 title 'FIR'\
,'../k40-ldg-spectra.dat' using 1:($3-$2) w lines ls 21 title 'FFT'\
,'../gtx750ti-spectra.dat' using 1:2 w lines ls 30 title 'FIR'\
,'../gtx750ti-spectra.dat' using 1:($3-$2) w lines ls 31 title 'FFT'

#------------------------------------------------------------------------

#----------------------------- Right -------------------------------------------
set size 0.48,0.93
set origin 0.50,0.0
set ylabel "GFLOPS" offset 0,0 font "Times-Italic,20"
set xlabel "Spectra" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("15000" 15000, "30000" 30000, "60000" 60000, "120000" 120000, "240000" 240000)
set ytics ("100" 1e11, "150" 1.5e11, "200" 2e11, "250" 2.5e11, "300" 3e11, "350" 3.5e11, "400" 4e11)
set logscale x
set grid
unset key
set yrange[0.8e11:3.2e11]
plot '../gtx580-spectra.dat' using 1:5 w lines ls 10 title 'FIR'\
,'../k40-ldg-spectra.dat' using 1:5 w lines ls 20 title 'FIR'\
,'../gtx750ti-spectra.dat' using 1:5 w lines ls 30 title 'FIR'
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
set label 1 "Fermi [FIR]" at -0.6,0.1 font "Times-Italic,18"
set label 2 "Kepler [FIR]" at 2.4,0.1 font "Times-Italic,18"
set label 3 "Maxwell [FIR]" at 5.4,0.1 font "Times-Italic,18"
set label 4 "Fermi [FFT]" at -0.6,-1.9 font "Times-Italic,18"
set label 5 "Kepler [FFT]" at 2.4,-1.9 font "Times-Italic,18"
set label 6 "Maxwell [FFT]" at 5.3,-1.9 font "Times-Italic,18"
plot (x>1.0 && x<2.0 ? 0:1/0) ls 10,\
(x>4.0 && x<5.0 ? 0:1/0) ls 20,\
(x>7.0 && x<8.0 ? 0:1/0) ls 30,\
(x>1.0 && x<2.0 ? -2:1/0) ls 11,\
(x>4.0 && x<5.0 ? -2:1/0) ls 21,\
(x>7.0 && x<8.0 ? -2:1/0) ls 31

unset multiplot

reset
#********************************************************************
#*	default line styles
#********************************************************************
set style line 10 linetype -1 linewidth 2 linecolor rgb 'black'
set style line 11 linetype  6 linewidth 2 linecolor rgb 'black'
set style line 20 linetype -1 linewidth 2 linecolor rgb 'blue'
set style line 21 linetype  6 linewidth 2 linecolor rgb 'blue'
set style line 30 linetype -1 linewidth 2 linecolor rgb 'green'
set style line 31 linetype  6 linewidth 2 linecolor rgb 'green'

#--------------------- TAPS --------------------------------------
#------------------------------------------------------------------
#-------------------> Serial spectra
set term postscript eps enhanced 'Times-Italic' 16
set output "GPU-taps-tf.eps"

set size 1.0, 1.0
set origin 0.0, 0.0
set multiplot
#----------------------------- Left -------------------------------------------
set size 0.48,0.93
set origin 0.0,0.0
set ylabel "Time[s]" offset 0,0 font "Times-Italic,20"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("5" 5, "8" 8, "16" 16, "32" 32, "64" 64)
set grid
set logscale x
unset key
plot '../gtx580-taps.dat' using 8:2 w lines ls 10 title 'FIR'\
,'../gtx580-taps.dat' using 8:($3-$2) w lines ls 11 title 'FFT'\
,'../k40-ldg-taps.dat' using 8:2 w lines ls 20 title 'FIR'\
,'../k40-ldg-taps.dat' using 8:($3-$2) w lines ls 21 title 'FFT'\
,'../gtx750ti-taps.dat' using 8:2 w lines ls 30 title 'FIR'\
,'../gtx750ti-taps.dat' using 8:($3-$2) w lines ls 31 title 'FFT'
#------------------------------------------------------------------------

#----------------------------- Right -------------------------------------------
set size 0.48,0.93
set origin 0.50,0.0
set ylabel "GFLOPS" offset 0,0 font "Times-Italic,20"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("5" 5, "8" 8, "16" 16, "32" 32, "64" 64)
set ytics ("100" 1e11, "150" 1.5e11,"200" 2e11, "250" 2.5e11,"300" 3e11, "350" 3.5e11,"400" 4e11, "450" 4.5e11,"500" 5e11, "550" 5.5e11, "600" 6e11, "650" 6.5e11, "700" 7e11, "750" 7.5e11, "800" 8e11, "850" 8.5e11)
set logscale x
set yrange[0.8e11:3.8e11]
set grid
unset key
plot '../gtx580-taps.dat' using 8:5 w lines ls 10 title '[FIR]'\
,'../k40-ldg-taps.dat' using 8:5 w lines ls 20 title '[FIR]'\
,'../gtx750ti-taps.dat' using 8:5 w lines ls 30 title '[FIR]'
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
set label 1 "Fermi [FIR]" at -0.6,0.1 font "Times-Italic,18"
set label 2 "Kepler [FIR]" at 2.4,0.1 font "Times-Italic,18"
set label 3 "Maxwell [FIR]" at 5.4,0.1 font "Times-Italic,18"
set label 4 "Fermi [FFT]" at -0.6,-1.9 font "Times-Italic,18"
set label 5 "Kepler [FFT]" at 2.4,-1.9 font "Times-Italic,18"
set label 6 "Maxwell [FFT]" at 5.3,-1.9 font "Times-Italic,18"
plot (x>1.0 && x<2.0 ? 0:1/0) ls 10,\
(x>4.0 && x<5.0 ? 0:1/0) ls 20,\
(x>7.0 && x<8.0 ? 0:1/0) ls 30,\
(x>1.0 && x<2.0 ? -2:1/0) ls 11,\
(x>4.0 && x<5.0 ? -2:1/0) ls 21,\
(x>7.0 && x<8.0 ? -2:1/0) ls 31
unset multiplot
reset;

#********************************************************************
#*	default line styles
#********************************************************************
set style line 10 linetype -1 linewidth 2 linecolor rgb 'black'
set style line 11 linetype  6 linewidth 2 linecolor rgb 'black'
set style line 20 linetype -1 linewidth 2 linecolor rgb 'blue'
set style line 21 linetype  6 linewidth 2 linecolor rgb 'blue'
set style line 30 linetype -1 linewidth 2 linecolor rgb 'green'
set style line 31 linetype  6 linewidth 2 linecolor rgb 'green'
#------------------------ Channels --------------------------------
#-------------------> Serial spectra
set term postscript eps enhanced 'Times-Italic' 16
set output "GPU-channels-tf.eps"

set size 1.0, 1.0
set origin 0.0, 0.0
set multiplot

#----------------------------- Left -------------------------------------------
set size 0.48,0.93
set origin 0.0,0.0
set ylabel "Time[s]" offset 0,0 font "Times-Italic,20"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("256" 256, "512" 512, "1024" 1024, "2048" 2048, "4096" 4096)
set grid
set logscale x
unset key
plot '../gtx580-channels.dat' using 11:2 w lines ls 10 title 'Simple [FIR]'\
,'../gtx580-channels.dat' using 11:($3-$2) w lines ls 11 title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:2 w lines ls 20 title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:($3-$2) w lines ls 21 title 'Simple [FIR]'\
,'../gtx750ti-channels.dat' using 11:2 w lines ls 30 title 'Simple [FIR]'\
,'../gtx750ti-channels.dat' using 11:($3-$2) w lines ls 31 title 'Simple [FIR]'\
#------------------------------------------------------------------------

#----------------------------- Right -------------------------------------------
set size 0.48,0.93
set origin 0.50,0.0
set ylabel "GFLOPS" offset 0,0 font "Times-Italic,20"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,20"
set xtics ("256" 256, "512" 512, "1024" 1024, "2048" 2048, "4096" 4096)
set ytics ("100" 1e11, "150" 1.5e11,"200" 2e11, "250" 2.5e11,"300" 3e11, "350" 3.5e11,"400" 4e11, "450" 4.5e11,"500" 5e11, "550" 5.5e11, "600" 6e11, "650" 6.5e11, "700" 7e11, "750" 7.5e11, "800" 8e11, "850" 8.5e11)
set logscale x
set yrange[0.8e11:3.2e+11]
set grid
unset key
plot '../gtx580-channels.dat' using 11:5 w lines ls 10 title 'Simple [FIR]'\
,'../k40-ldg-channels.dat' using 11:5 w lines ls 20 title 'Simple [FIR]'\
,'../gtx750ti-channels.dat' using 11:5 w lines ls 30 title 'ldg [FIR]'

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
set label 1 "Fermi [PPF]" at -0.6,0.1 font "Times-Italic,18"
set label 2 "Fermi [FIR]" at 2.4,0.1 font "Times-Italic,18"
set label 3 "Kepler [PPF]" at 5.4,0.1 font "Times-Italic,18"
set label 4 "Kepler [FIR]" at -0.6,-1.9 font "Times-Italic,18"
set label 5 "cuFFT" at 2.4,-1.9 font "Times-Italic,18"
#set label 6 "p_{m}(a=1.1)" at 5.7,-1.9 font "Times-Italic,18"
plot (x>1.0 && x<2.0 ? 0:1/0) ls 10,\
(x>4.0 && x<5.0 ? 0:1/0) ls 11,\
(x>7.0 && x<8.0 ? 0:1/0) ls 20,\
(x>1.0 && x<2.0 ? -2:1/0) ls 21,\
(x>4.0 && x<5.0 ? -2:1/0) ls 30
#(x>7.0 && x<8.0 ? -2:1/0) lw 5 lt 6 lc rgb 'black'

unset multiplot
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

reset
set term postscript eps enhanced color 'Times-Italic' 18
set title ""
set output "GPU-threads.eps"
set grid
set logscale x
set xtics ("32" 32, "64" 64,"128" 128, "192" 192, "256" 256, "384" 384, "512" 512, "768" 768, "1024" 1024)
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Threads" rotate by 90  offset 1,0 font "Times-Italic,22"
plot '../gtx580-nthreads.dat' using 10:2 with lines lw 3,\
'../gtx750ti-nthreads.dat' using 10:2 with lines lw 3,\
'../k40-nthreads.dat' using 10:2 with lines lw 3
unset output
