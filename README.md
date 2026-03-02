# riptide_optimization

## To build and execute debug mode
cmake -S . -B build/ -G"Ninja Multi-Config"
cmake --build build/ --config Debug

## To start program in visualization mode
./build/Debug/simulate -g geometry/main.gdml -m macros/run1.mac -v

## To start 2D image visualization program
./build/analysis/Debug/plot2D

x1=83.5 x2=153.5
x1min=-25.1 x2min=45.4 

step = 3mm
  x1=94.90 mm
  x2=186.40 mm

  x1=1.9 mm
  x2=180.4 mm