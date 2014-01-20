set term postscript eps enhanced color 'Times-Italic' 18
set output "FIR_time_Channels.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("512" 512,"768" 768,"1024" 1024,"2048" 2048,"4096" 4096)
set key top left
set logscale x
set grid
plot 'dataZen/pwi3_C_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'gray' title 'Intrinsic kernel'\
,'dataZen/pwi4_C_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'black' title 'Align kernel'\
,'dataZen/pwi6_C_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'green' title 'Cache kernel'
reset;



set term postscript eps enhanced color 'Times-Italic' 18
set output "bandwidth_Channels.eps"
set ylabel "Bandwidth[GB]" offset 0,0 font "Times-Italic,22"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("80" 8.0e+11,"90" 9.0e+10,"100" 1.0e+11,"110" 1.1e+11, "120" 1.2e+11, "130" 1.3e+11, "140" 1.4e+11, "150" 1.5e+11, "160" 1.6e+11, "170" 1.7e+11, "180" 1.8e+11)
set xtics ("512" 512,"768" 768,"1024" 1024,"2048" 2048,"4096" 4096)
set key top right
set yrange[8.8e+10:1.7e+11]
set logscale x
set grid
plot 'dataZen/pwi3_C_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'gray' title 'Intrinsic kernel'\
,'dataZen/pwi4_C_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'black' title 'Align kernel'\
,'dataZen/pwi6_C_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'green' title 'Cache kernel'
reset;



set term postscript eps enhanced color 'Times-Italic' 18
set output "flops_Channels.eps"
set ylabel "Gflops" offset 0,0 font "Times-Italic,22"
set xlabel "Channels" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("0.1" 1e8, "1" 1e9, "10" 1e10,"25" 2.5e10,"35" 3.5e10,"40" 4.0e10,"45" 4.5e10,"50" 5.0e10,"75" 7.5e10, "100" 1e11)
set xtics ("512" 512,"768" 768,"1024" 1024,"2048" 2048,"4096" 4096)
set yrange[2.0e+10:4.0e+10]
set key top right
set logscale x
set grid
plot 'dataZen/pwi3_C_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'gray' title 'Intrinsic kernel'\
,'dataZen/pwi4_C_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'black' title 'Align kernel'\
,'dataZen/pwi6_C_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'green' title 'Cache kernel'
reset;
