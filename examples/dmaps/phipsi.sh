#!/bin/sh

cd lsdmap
g_chi -phi -psi -all -xvg none -s lsdmap.gro -f lsdmap.gro 1>/dev/null

cd ..

cd fit
g_chi -phi -psi -all -xvg none -s fit.gro -f fit.gro 1>/dev/null
cd ..


g_chi -phi -psi -all -xvg none -s confall.gro -f confall.gro 1>/dev/null
