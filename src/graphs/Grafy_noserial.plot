set term postscript eps enhanced color 'Times-Italic' 18
set output "FIR_time.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Size of datafile[MB] Log" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000,"30000" 30000,"60000" 60000,"120000" 120000,"240000" 240000)
set key top left
set logscale x
set grid
plot 'dataZen/SSE42_pwi4_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'black' title 'Parallel with intrinsic XEON 5650 (SSE4.2)'\
,'dataZen/AVX_pwi4_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'
reset;



set term postscript eps enhanced color 'Times-Italic' 18
set output "Whole_time.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Number of spectra calculated" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000,"30000" 30000,"60000" 60000,"120000" 120000,"240000" 240000)
set key top left
set logscale x
set grid
plot 'dataZen/SSE42_pwi4_perf.dat' using 1:3 w lines lt -1 lw 3 lc rgb 'black' title 'Parallel with intrinsic XEON 5650 (SSE4.2)'\
,'dataZen/AVX_pwi4_perf.dat' using 1:3 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'
reset;




set term postscript eps enhanced color 'Times-Italic' 18
set output "bandwidth.eps"
set ylabel "Bandwidth[GB]" offset 0,0 font "Times-Italic,22"
set xlabel "Number of spectra calculated" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("140" 1.4e+11, "150" 1.5e+11, "160" 1.6e+11, "170" 1.7e+11, "180" 1.8e+11, "190" 1.9e+11, "200" 2.0e+11, "210" 2.1e+11)
set xtics ("15000" 15000,"30000" 30000,"60000" 60000,"120000" 120000,"240000" 240000)
set key top right horizontal
set yrange[1.4e+11:2.1e+11]
set logscale x
set grid
plot 'dataZen/SSE42_pwi4_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'black' title 'Parallel with intrinsic XEON 5650 (SSE4.2)'\
,'dataZen/AVX_pwi4_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'
reset;





set term postscript eps enhanced color 'Times-Italic' 18
set output "flops.eps"
set ylabel "Gflops" offset 0,0 font "Times-Italic,22"
set xlabel "Number of spectra calculated" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("0.1" 1e8, "1" 1e9, "10" 1e10,"25" 2.5e10,"35" 3.5e10,"40" 4.0e10,"45" 4.5e10,"50" 5.0e10,"75" 7.5e10, "100" 1e11)
set xtics ("15000" 15000,"30000" 30000,"60000" 60000,"120000" 120000,"240000" 240000)
set yrange[3.0e+10:5.0e+10]
set key top right horizontal
set logscale x
set grid
plot 'dataZen/SSE42_pwi4_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'black' title 'Parallel with intrinsic XEON 5650 (SSE4.2)'\
,'dataZen/AVX_pwi4_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'
reset;



