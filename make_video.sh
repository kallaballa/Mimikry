cd result
cp 0000000000_00099_result.png lastResult.png

ls -1 *_00099_result.png | while read f; do 
	PREFIX="`echo $f | cut -d"_" -f1-2`"; 

	if [ ! -f ${PREFIX}_frame.png ]; then
  	compare -metric AE -fuzz 10% $f lastResult.png null: &> /dev/null
  	if [ $? -ne 0 ]; then
			convert $f ${PREFIX}_dft.png ${PREFIX}_bdft.png +append ${PREFIX}_append.png; 
			convert ${PREFIX}_append.png ../footer.png -append ${PREFIX}_montage.png;
			convert ${PREFIX}_montage.png -font Monospace -background none -stroke black -strokewidth 1 -fill white -gravity south -pointsize 28 label:"`echo ${PREFIX} | cut -d"_" -f1`" -composite ${PREFIX}_frame.png
			cp $f lastResult.png
			echo "processed $f"
		else
			echo "skipped $f"
		fi
	fi;
done

mencoder mf://*_00099_frame.png -mf w=300:h=225:fps=1:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
