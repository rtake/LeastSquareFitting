#/bin/csh

set i=5
while ( $i < 10 ) 
	#cat ./header ./n20_h23_0${i}.xyz > ./n20_h23_0${i}.gjf
	g09 ./n20_h23_0${i}.gjf 
	@ i = $i + 1
end

while ( $i < 21 )
	#cat ./header ./n20_h23_${i}.xyz > ./n20_h23_${i}.gjf
	g09 ./n20_h23_${i}.gjf 
	@ i = $i + 1
end

