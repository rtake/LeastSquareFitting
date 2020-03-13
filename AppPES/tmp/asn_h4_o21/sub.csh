#/bin/csh

set i=5
while ( $i < 10 ) 
	#cat ./header ./h4_o21_0${i}.xyz > ./h4_o21_0${i}.gjf
	g09 ./h4_o21_0${i}.gjf 
	@ i = $i + 1
end

while ( $i < 21 )
	#cat ./header ./h4_o21_${i}.xyz > ./h4_o21_${i}.gjf
	g09 ./h4_o21_${i}.gjf 
	@ i = $i + 1
end

