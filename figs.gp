
######################################
###Figure with examples
######################################
set term pdfcairo enhanced size 7.5,3 font "Arial-Bold,12"
set out "Model_examples.pdf"




################################################################

################################################################
set multi lay 2,2 title '' font "Arial-Bold,12


#set style line 101 lc rgb 'black' lt 1 lw 2.5
#set border 15 front ls 101
#set tics nomirror out scale 1.0



################################################################
#Example traces of bursting neurons with increasing gspk or gAHP
################################################################


unset border
unset xrange
unset yrange
set yrange[-62:-5]
unset xtic
unset ytic
skip=1
gap = 2
shift = 10
y1=-5


set label 1 "Control" at 0*(shift+gap)+5,y1 front center tc rgb "black"
set label 2 "gSPK=2nS" at 1*(shift+gap)+5,y1 front center tc rgb "blue"
set label 3 "gSPK=4nS" at 2*(shift+gap)+5,y1 front center tc rgb "blue"
set label 4 "gSPK=6nS" at 3*(shift+gap)+5,y1 front center tc rgb "blue"


set label 5 "1s" at 3*(shift+gap)+9+.5,-62 front center
set arrow from 3*(shift+gap)+9,-59 to 3*(shift+gap)+10,-59 front nohead lw 4
set label 6 "20mV" at 3.5,-45 front left
set arrow from 3,-55 to 3,-35 front nohead lw 4

plot "data/2_0_0.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+0*gap):0/0):($2-0*50) w l lw 2 lc rgb 'black' t"",\
     "data/2_0_1.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+1*(gap+shift)):0/0):($2-0*50) w l lw 2 lc rgb 'blue' t"",\
     "data/2_0_2.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+2*(gap+shift)):0/0):($2-0*50) w l lw 2 lc rgb 'blue' t"",\
     "data/2_0_3.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+3*(gap+shift)):0/0):($2-0*50) w l lw 2 lc rgb 'blue' t"",\



set label 1 "Control" at 0*(shift+gap)+5,y1 front center tc rgb "black"
set label 2 "gAHP=6nS" at 1*(shift+gap)+5,y1 front center tc rgb "green"
set label 3 "gAHP=12nS" at 2*(shift+gap)+5,y1 front center tc rgb "green"
set label 4 "gAHP=18nS" at 3*(shift+gap)+5,y1 front center tc rgb "green"
plot "data/2_1_0.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+0*gap):0/0):($2-0*50) w l lw 2 lc rgb 'black' t"",\
     "data/2_1_1.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+1*(gap+shift)):0/0):($2-0*50) w l lw 2 lc rgb 'green' t"",\
     "data/2_1_2.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+2*(gap+shift)):0/0):($2-0*50) w l lw 2 lc rgb 'green' t"",\
     "data/2_1_3.sp" every skip u ($1>=(30) & $1<=(40) ?($1-30+3*(gap+shift)):0/0):($2-0*50) w l lw 2 lc rgb 'green' t"",\



unset arrow
unset label
unset yrange



#skip=1
set size 1,0.925*1.0/2
set origin 0.00, 0.02
set arrow 1 from 2*(shift+gap)+9,graph -0.05 to first 2*(shift+gap)+10,graph -.05 front nohead lw 4
set label 1 "1s" at 2*(shift+gap)+9+.5,graph -.095 front center
set arrow 2 from graph .035,first 5 to graph .035, first 25  front nohead lw 4
set label 2 "20Hz" at graph 0.0385, first 15 left front


set arrow 3 from 0*(shift+gap),graph 1.025 to first 0*(shift+gap)+10, graph 1.025 front nohead lw 3 lc rgb "black"
set label 3 "Control (100% Burst Capable)"at first 0*(shift+gap)+5, graph 1.085 front center tc rgb"black"
set arrow 4 from 1*(shift+gap),graph 1.025 to first 1*(shift+gap)+10, graph 1.025 front nohead lw 3 lc rgb "blue"
set label 4 "gSPK=15nS (0% Burst Capable)"at first 1*(shift+gap)+5, graph 1.085 front center tc rgb"blue"
set arrow 5 from 2*(shift+gap),graph 1.025 to first 2*(shift+gap)+10, graph 1.025 front nohead lw 3 lc rgb "green"
set label 5 "gAHP=30nS (0% Burst Capable)"at first 2*(shift+gap)+5, graph 1.085 front center tc rgb"green"

t0=65.25
dt=6
#set xrange[t0:t0+dt]
#set xrange[70:80]
set xrange[0:2*(shift+gap)+10]
#set yrange[-10:*]
#unset xrange
#set yrange[0:4]
#set arrow 3 from graph 0,first 1.8335 to graph 1,first 1.8335 front nohead lw 1 lc rgb "red"
#set arrow 4 from graph 0,second 74 to graph 1,second 74 front nohead lw 3 dt"-" lc rgb "orange"
plot "data/3_0_0.hst" every skip u ($1>=(69.5) & $1<=(79.5) ?($1-69.5+0*(gap+shift)):0/0):2  w p pt 5 ps .1 lt rgb "#10a0a0a0" axis x1y2 title '',\
     "data/3_0_0.sp"             u ($1>=(69.5) & $1<=(79.5) ?($1-69.5+0*(gap+shift)):0/0):2 w l lw 3 lc rgb "black" title"",\
     "data/3_1_0.hst" every skip u ($1>=(69.5) & $1<=(79.5) ?($1-69.5+1*(gap+shift)):0/0):2  w p pt 5 ps .1 lt rgb "#10a0a0a0" axis x1y2 title '',\
     "data/3_1_0.sp"             u ($1>=(69.5) & $1<=(79.5) ?($1-69.5+1*(gap+shift)):0/0):2 w l lw 3 lc rgb "black" title"",\
     "data/3_2_0.hst" every skip u ($1>=(69.5) & $1<=(79.5) ?($1-69.5+2*(gap+shift)):0/0):2  w p pt 5 ps .1 lt rgb "#10a0a0a0" axis x1y2 title '',\
     "data/3_2_0.sp"             u ($1>=(69.5) & $1<=(79.5) ?($1-69.5+2*(gap+shift)):0/0):2 w l lw 3 lc rgb "black" title""





