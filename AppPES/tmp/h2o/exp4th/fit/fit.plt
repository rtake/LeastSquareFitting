set xr[0:3]
set yr[0:3]
set xlabel 'O1-H2'
set ylabel 'O1-H3'
appV(x, y) = 344.127 * ( ( (1 - exp(-0.5 * x) ) ) ** 2 ) - 140.528 * ( ( (1 - exp(-0.5 * x) ) ) ** 1 ) * ( ( (1 - exp(-0.5 * y) ) ) ** 2 ) + ( -9.50408 - 242.609 ) * ( ( (1 - exp(-0.5 * x) ) ) ** 1 )

set palette gray 
splot appV(x, y) with pm3d 
