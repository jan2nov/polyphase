set term postscript eps enhanced color 'Times-Italic' 18
set output "FIR_time.eps"
set ylabel "Time[s]" offset 0,0 font "Times-Italic,22"
set xlabel "Size of datafile[MB] Log" rotate by 90  offset 1,0 font "Times-Italic,22"
set xtics ("15000" 15000,"30000" 30000,"60000" 60000,"120000" 120000,"240000" 240000)
set key top left
set logscale x
set grid
plot 'dataZen/SSE42_pwi4_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'black' title 'Parallel with intrinsic XEON 5650 (SSE4.2)'\
,'dataZen/AVX_pwi4_perf.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'\
,'dataZen/gpu_perf-k20m.dat' using 1:2 w lines lt -1 lw 3 lc rgb '#32CD32' title 'Cache kernel without stream on Tesla k20m'\
,'dataZen/gpu_perf-480.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'navyblue' title 'Cache kernel without stream on GTX 480'
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
,'dataZen/AVX_pwi4_perf.dat' using 1:3 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'\
,'dataZen/gpu_perf-gtx480-all.dat' using 1:2 w lines lt -1 lw 3 lc rgb 'navyblue' title 'Cache kernel without stream on GTX 480'
reset;




set term postscript eps enhanced color 'Times-Italic' 18
set output "bandwidth.eps"
set ylabel "Bandwidth[GB]" offset 0,0 font "Times-Italic,22"
set xlabel "Number of spectra calculated" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("140" 1.4e+11, "160" 1.60e+11, "180" 1.8e+11, "200" 2.0e+11, "220" 2.2e+11, "240" 2.4e+11, "260" 2.6e+11, "280" 2.8e+11, "300" 3.0e+11, "320" 3.2e+11, "340" 3.4e+11, "360" 3.6e+11\
,"380" 3.8e+11, "400" 4.0e+11, "420" 4.2e+11)
set xtics ("15000" 15000,"30000" 30000,"60000" 60000,"120000" 120000,"240000" 240000)
set key top right horizontal
set yrange[1.4e+11:4.2e+11]
set logscale x
set grid
plot 'dataZen/SSE42_pwi4_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'black' title 'Parallel with intrinsic XEON 5650 (SSE4.2)'\
,'dataZen/AVX_pwi4_perf.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'\
,'dataZen/gpu_perf-k20m.dat' using 1:4 w lines lt -1 lw 3 lc rgb '#32CD32' title 'Cache kernel without stream on Tesla k20m'\
,'dataZen/gpu_perf-480.dat' using 1:4 w lines lt -1 lw 3 lc rgb 'navyblue' title 'Cache kernel without stream on GTX 480'
reset;





set term postscript eps enhanced color 'Times-Italic' 18
set output "flops.eps"
set ylabel "Gflops" offset 0,0 font "Times-Italic,22"
set xlabel "Number of spectra calculated" rotate by 90  offset 1,0 font "Times-Italic,22"
set ytics ("0.1" 1e8, "1" 1e9, "10" 1e10,"25" 2.5e10,"35" 3.5e10,"40" 4.0e10,"45" 4.5e10,"50" 5.0e10,"75" 7.5e10, "100" 1e11, "115" 1.15e+11, "125" 1.25e+11, "140" 1.4e+11)
set xtics ("15000" 15000,"30000" 30000,"60000" 60000,"120000" 120000,"240000" 240000)
set yrange[3.0e+10:1.4e+11]
set key top right horizontal
set logscale x
set grid
plot 'dataZen/SSE42_pwi4_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'black' title 'Parallel with intrinsic XEON 5650 (SSE4.2)'\
,'dataZen/AVX_pwi4_perf.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'red' title 'Parallel with intrinsic XEON E5-2650 (AVX)'\
,'dataZen/gpu_perf-k20m.dat' using 1:5 w lines lt -1 lw 3 lc rgb '#32CD32' title 'Cache kernel without stream on Tesla k20m'\
,'dataZen/gpu_perf-480.dat' using 1:5 w lines lt -1 lw 3 lc rgb 'navyblue' title 'Cache kernel without stream on GTX 480'
reset;



