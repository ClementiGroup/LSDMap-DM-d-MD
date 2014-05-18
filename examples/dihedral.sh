# example of how to use LSDMap with dihedral metric using .xvg files

set -e

sed -i 's/metric=.*/metric=dihedral/g' config.ini
sed -i 's/status=.*/status=constant/g' config.ini
sed -i 's/epsilon=.*/epsilon=0.15/g' config.ini

rm -rf xvgfiles
mkdir xvgfiles

# create .xvg files (dihedral angles should be given in radians, so use flag -rad)
g_chi -psi -all -xvg none -rad -s ala12_1000.gro -f ala12_1000.gro &> /dev/null
mv *.xvg xvgfiles/
rm -f chi.log

xvglist=
for file in xvgfiles/psi*.xvg; do xvglist="$xvglist $file"; done

# run LSDMap with all .xvg files
lsdmap -f config.ini -c $xvglist -o psiALA
