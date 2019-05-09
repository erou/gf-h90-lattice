set terminal post eps color enhanced font "Times, 15" size 5cm, 4cm
set out "plot_embed.eps"
set rmargin 5
set lmargin 0.8
set datafile separator ","
set palette rgb 33,13,10 
set cbrange [0:90]
set cbtics (10, 30, 50, 70, 90) font ", 12" offset -1, 0
set colorbox user origin 0.9, 0.176 size 0.04, 0.76
set xtics (50, 100, 150, 200) font ", 12"
set xlabel "Degree of the destination algebra" offset 1.5, 0.5
set y2label "Level of the destination algebra" offset -1.3, 0
unset ylabel
set format y ''
set yrange [0.06:3000]
set logscale y

plot "<(sed -n '116,179p' benchmarks/embed-3.txt)" u 4:($5*1000):3 notitle with points pt 7 palette
