#!/bin/bash

sed -i '' 's/iter=.*/iter=0/g' settings
rm -rf iter*
