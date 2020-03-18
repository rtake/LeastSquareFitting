#/bin/csh

set i=0
while ( $i < 10 )
	echo $i 
	# cat ./header ./h2o_o1h3_0${i}.xyz > ./h2o_o1h3_0${i}.gjf
	g09 ./h2o_o1h3_0${i}.gjf
	@ i = $i + 1
end

while ( $i < 20 )
	echo $i
	# cat ./header ./h2o_o1h3_${i}.xyz > ./h2o_o1h3_${i}.gjf
	g09 ./h2o_o1h3_${i}.gjf
	@ i = $i + 1
end

