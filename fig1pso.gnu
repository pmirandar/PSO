set ylabel "Ad/mb"
set xlabel "Au/mt"
set key right
#set key box
mb = 2839 
mt = 168260
set pointsize 0.5
#set palette color
#set colorbox vertical
#set samples 2
plot [0:1][0:1] "run.dat" using (($2)/mt):(($3)/mb) title '4<Chi2' pointtype 30 lc "red"
