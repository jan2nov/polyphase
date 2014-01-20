set term postscript eps enhanced color 'Times-Italic' 18
set output "FIR_time_Taps.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("8" 8,"12" 12,"16" 16,"20" 20,"24" 24)
set key top left
set grid
plot 'dataZen/pwi3_T_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'gray' title 'Intrinsic kernel'\
,'dataZen/pwi4_T_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'black' title 'Aling kernel'\
,'dataZen/pwi6_T_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'green' title 'Cache kernel'
reset;



set term postscript eps enhanced color 'Times-Italic' 18
set output "bandwidth_Taps.eps"
set ylabel "Bandwidth[GB]" offset 0,0 font "Times-Italic,22"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("120" 1.2e+11,"130" 1.3e+11,"140" 1.4e+11, "150" 1.5e+11, "160" 1.6e+11, "170" 1.7e+11, "180" 1.8e+11, "190" 1.9e+11, "200" 2.0e+11, "210" 2.1e+11)
set xtics ("8" 8,"12" 12,"16" 16,"20" 20,"24" 24)
set key top left
set yrange[1.2e+11:2.1e+11]
set grid
plot 'dataZen/pwi3_T_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'gray' title 'Intrinsic kernel'\
,'dataZen/pwi4_T_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'black' title 'Aling kernel'\
,'dataZen/pwi6_T_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'green' title 'Cache kernel'
reset;





set term postscript eps enhanced color 'Times-Italic' 18
set output "flops_Taps.eps"
set ylabel "Gflops" offset 0,0 font "Times-Italic,22"
set xlabel "Taps" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("0.1" 1e8, "1" 1e9, "25" 2.5e10,"30" 3.0e10,"35" 3.5e10,"40" 4.0e10,"45" 4.5e10,"50" 5.0e10,"75" 7.5e10, "100" 1e11)
set xtics ("8" 8,"12" 12,"16" 16,"20" 20,"24" 24)
set yrange[2.5e10:5.5e+10]
set key top left
set grid
plot 'dataZen/pwi3_T_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'gray' title 'Intrinsic kernel'\
,'dataZen/pwi4_T_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'black' title 'Aling kernel'\
,'dataZen/pwi6_T_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'green' title 'Cache kernel'
reset;
