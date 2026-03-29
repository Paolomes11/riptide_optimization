#!/bin/bash
set -e

# Colori per output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo "=== RIPTIDE Test Suite ==="

# Directory build
BUILD_DIR="build"
PROJECT_ROOT=$(pwd)

if [ ! -d "$BUILD_DIR" ]; then
    echo -e "${RED}[ERROR] Build directory not found. Please build the project first.${NC}"
    exit 1
fi

# 1. Test lens_cutter_main
echo -e "\n[TEST 1] lens_cutter_main"
./$BUILD_DIR/Debug/lens_cutter_main --id LB4553 > /dev/null
if [ $? -eq 0 ]; then
    echo -e "${GREEN}[PASS] lens_cutter_main basic execution${NC}"
else
    echo -e "${RED}[FAIL] lens_cutter_main execution failed${NC}"
    exit 1
fi

# 2. Test lens_simulation_main (dry-run/batch)
echo -e "\n[TEST 2] lens_simulation_main (batch, few pairs)"
# Creiamo un config temporaneo per il test rapido
TEST_CONFIG="tests/test_config.json"
cat > $TEST_CONFIG <<EOF
{
  "x_min": 70.0,
  "x_max": 80.0,
  "dx": 5.0,
  "n_photons": 10,
  "lower_percentile": 0.05,
  "upper_percentile": 0.05,
  "pairs": [[75.0, 160.0]]
}
EOF

./$BUILD_DIR/Debug/lens_simulation_main -g geometry/main.gdml -l -b --config $TEST_CONFIG --output tests/test_lens.root > /dev/null
if [ -f "tests/test_lens.root" ]; then
    echo -e "${GREEN}[PASS] lens_simulation_main generated ROOT file${NC}"
else
    echo -e "${RED}[FAIL] lens_simulation_main failed to generate output${NC}"
    exit 1
fi

# 3. Test optimization_main (batch, few pairs)
echo -e "\n[TEST 3] optimization_main (batch, few pairs)"
./$BUILD_DIR/Debug/optimization_main -g geometry/main.gdml -o -b --config $TEST_CONFIG --output tests/test_opt.root > /dev/null
if [ -f "tests/test_opt.root" ]; then
    echo -e "${GREEN}[PASS] optimization_main generated ROOT file${NC}"
else
    echo -e "${RED}[FAIL] optimization_main failed to generate output${NC}"
    exit 1
fi

# 4. Test plot2D
echo -e "\n[TEST 4] plot2D"
./$BUILD_DIR/analysis/Debug/plot2D -i tests/test_opt.root -c $TEST_CONFIG -o tests/test_2d.png > /dev/null
if [ -f "tests/test_2d.png" ]; then
    echo -e "${GREEN}[PASS] plot2D generated PNG${NC}"
else
    echo -e "${RED}[FAIL] plot2D failed to generate PNG${NC}"
    exit 1
fi

# 5. Test plot3D
echo -e "\n[TEST 5] plot3D"
./$BUILD_DIR/analysis/Debug/plot3D -i tests/test_opt.root -c $TEST_CONFIG -o tests/test_3d.png > /dev/null
if [ -f "tests/test_3d.png" ]; then
    echo -e "${GREEN}[PASS] plot3D generated PNG${NC}"
else
    echo -e "${RED}[FAIL] plot3D failed to generate PNG${NC}"
    exit 1
fi

echo -e "\n=== All tests completed successfully! ==="
rm -f tests/test_lens.root tests/test_opt.root tests/test_2d.png tests/test_3d.png tests/test_config.json
