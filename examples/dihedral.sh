# example of how to use LSDMap with dihedral metric using .xvg files

set -e

sed -i 's/metric=.*/metric=dihedral/g' config.ini
sed -i 's/status=.*/status=constant/g' config.ini
sed -i 's/epsilon=.*/epsilon=0.15/g' config.ini

# run LSDMap with all .xvg files
lsdmap -f config.ini -c xvgfiles/psiALA* -o psiALA
