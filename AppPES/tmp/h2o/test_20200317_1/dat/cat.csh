#/bin/csh

set i=0
while ( $i < 10 )
	echo $i 
	cat ./h2o_o1h3_0${i}.xyz >> list.xyz
	@ i = $i + 1
end

while ( $i < 20 )
	echo $i
	cat ./h2o_o1h3_${i}.xyz >> list.xyz
	@ i = $i + 1
end

