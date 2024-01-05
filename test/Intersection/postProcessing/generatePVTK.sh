head -n 3 template.pvtk > out.pvtk
find ../DUMP/ -maxdepth 1 -type f -iname "*.vtk" -print0 | xargs -0 -I {} echo "{}" >> out.pvtk
tail -n 1 template.pvtk >> out.pvtk
NN=$(grep -ri "../DUMP" out.pvtk | wc -l)
echo $NN
sed -i "s/numberOfPieces=\"1\"/numberOfPieces=\"$NN\"/" out.pvtk
sed -i "s@../DUMP@    <Piece fileName=\"../DUMP@" out.pvtk
sed -i -e '3,$s@.vtk@.vtk\" />@' out.pvtk
