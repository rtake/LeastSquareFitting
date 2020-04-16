# MIN

0 1
O          0.215688960000          0.000000000000         -1.343610340000
H          0.215688960000          0.759337000000         -0.747567340000
H          0.215688960000         -2.359828796000          0.508741692000
Options
fitorder=4
reference=h2o_lup_LUPOUTt.log.xyz
opt=internal
maxoptitr=10000
stepsize=0.02
fit=leastsquare
