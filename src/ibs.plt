ps = 0;

if (ps == 0) \
  set terminal x11; \
else if (ps == 1) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 18; \
else if (ps == 2) \
  set terminal pdf enhanced color solid linewidth 2 font "Times-Roman, 18";

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";
set style line 4 lt 1 lw 1 lc rgb "dark-orange";

if (ps == 1) \
  set output "ibs.eps"; \
else if (ps == 2) \
  set output "ibs.pdf"

set title "Impact of IBS: {/Symbol e}_{x,y} = 16 pm.rad" \
          . ", {/Symbol s_d} = 1{/Symbol \264}10^{-3}, I_b = 5 nC";
set xlabel "Natural {/Symbol s}_s [cm]"; set ylabel "[pm.rad]";
set y2label "[10^{-3}]";
set ytics nomirror; set y2tics;
plot "ibs.out" using 1:2 title "{/Symbol e}_{x,y}" with lines ls 1, \
     "ibs.out" using 1:(1e3*$7) axis x1y2 title "{/Symbol s_d}" \
     with lines ls 2;

if (ps == 0) pause mouse "click on graph to cont.\n";
