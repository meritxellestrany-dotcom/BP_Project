NACCESS_BIN="../NACCESS/naccess"
SCAN_DIR="ala_scanning"

cd "$SCAN_DIR" || exit

for pdb in *.pdb; do
    "$NACCESS_BIN" "$pdb"
done