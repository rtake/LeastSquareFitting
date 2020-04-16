set xr[0.0:5.0]
set yr[0.0:5.0]
set xlabel 'x1 (O1 - H3)'
set ylabel 'x2 (H2 - H3)'

app( x, y ) = (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 4 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 3 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + ( -119.142417907715 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 3 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 3 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 2 ) * ( (1 - exp(-0.5 * x) ) ** 2 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 2 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 2 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 2 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (  298.010620117188 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 2 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 2 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 2 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 3 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 2 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 2 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 2 ) + (  -77.517105102539 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 2 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 1 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 3 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 4 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 3 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 3 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 2 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 2 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 2 ) * ( (1 - exp(-0.5 * y) ) ** 2 ) + (   69.074951171875 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + ( -213.078674316406 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 2 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 1 ) * ( (1 - exp(-0.5 * y) ) ** 3 ) + (  -69.178894042969 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 0 ) + (  -26.502862930298 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 1 ) + (   55.280437469482 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 2 ) + (    0.000000000000 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 3 ) + (   17.691976547241 ) * ( (1 - exp(-0.5 * 0.96533) ) ** 0 ) * ( (1 - exp(-0.5 * x) ) ** 0 ) * ( (1 - exp(-0.5 * y) ) ** 4 )

set isosamples 50
set palette gray

splot app(x, y) with pm3d

