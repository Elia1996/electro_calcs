#!/bin/bash

cd ..
git clone https://github.com/Elia1996/myplot.git
myplot_dir=$(pwd)/myplot
cd ..
cd electro_calcs
echo "export PYTHONPATH=$PYTHONPATH:$(pwd)/lib:$myplot_dir" >> ~/.bashrc
echo "To complete installation close and reopen terminal
you should have now a new path added to PYTHONPATH, and 
you will be able to import cross_talk and plotter module
from whatever script!!!"
