set xr[0.5:2]
set yr[1.0:2.5]
set xlabel 'O1-H3'
set ylabel 'H2-H3'
appV(x, y) = 754.942260742188 * ( ( (1 - exp(-0.5 * x) ) ) ** 4 ) - 2605.0244140625 * ( ( (1 - exp(-0.5 * x) ) ) ** 3 ) + 4188.0419921875  * ( ( (1 - exp(-0.5 * x) ) ) ** 2 ) + 394.569519042969 * ( ( (1 - exp(-0.5 * x) ) ) ** 1 ) - 167.891464233398 * ( ( (1 - exp(-0.5 * x) ) ) ** 1 ) * ( ( (1 - exp(-0.5 * y) ) ) ** 1 ) + 1550.93994140625 - 4187.25 * ( ( (1 - exp(-0.5 * y) ) ) ** 1 )

set isosamples 50 
set palette gray 

splot 'h2o.ref' using 3:4:1 w l
replot appV(x, y)
# replot appV(x, y) with pm3d
