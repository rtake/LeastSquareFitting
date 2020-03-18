#/bin/csh

set i=0
while ( $i < 180 )
	echo $i 
	# cp ./h2o_h2h3.xyz ./h2o_h2h3_${i}.xyz
	# cat ./header ./h2o_h2h3_${i}.xyz > ./h2o_h2h3_${i}.gjf
	g09 ./h2o_h2h3_${i}.gjf
	@ i = $i + 10
end

