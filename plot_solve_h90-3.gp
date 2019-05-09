set terminal post eps color enhanced font "Times, 15" size 5cm, 4cm
set out "plot_solve_h90.eps"
set lmargin 2.8
set rmargin 0.5 
set palette rgb 33,13,10 
set cbrange [0:90]
set cbtics (10, 30, 50, 70, 90)
unset colorbox
set datafile separator ","
set xtics (50, 100, 150, 200) font ", 12" offset -1, 0
set xlabel "Degree of the algebra" offset 0, 0.5
set ylabel "Time (ms)" offset 5.5, -1 
set ytics font ", 12" offset 0.9, 0
unset y2label
set yrange [0.06:3000]
set logscale y
plot 'benchmarks/solve_h90-3.txt' u 2:($3*1000):1 notitle with points pt 7 palette
