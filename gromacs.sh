#!/bin/bash

# GROMACS MD Simulation Setup Script
# Usage: ./MDsimscript.sh -i PDBfile -o outdir -t time (e.g., 100ns)

# Default values
INPUT_PDB=""
OUTPUT_DIR=""
SIM_TIME="10ns"
FORCE_FIELD="amber99sb-ildn"
WATER_MODEL="tip3p"

# Parse command line arguments
while getopts "i:o:t:" opt; do
  case $opt in
    i) INPUT_PDB="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    t) SIM_TIME="$OPTARG" ;;
    *) echo "Usage: $0 -i input.pdb -o output_dir [-t simulation_time (e.g., 100ns)]"
       exit 1 ;;
  esac
done

# Check required arguments
if [ -z "$INPUT_PDB" ] || [ -z "$OUTPUT_DIR" ]; then
  echo "Error: Missing required arguments"
  echo "Usage: $0 -i input.pdb -o output_dir [-t simulation_time (e.g., 100ns)]"
  exit 1
fi

# Extract base filename
BASENAME=$(basename "$INPUT_PDB" .pdb)

# Extract time and convert
TIME_VALUE=$(echo "$SIM_TIME" | grep -oE '[0-9]+')
TIME_UNIT=$(echo "$SIM_TIME" | grep -oE '[a-zA-Z]+')

case $TIME_UNIT in
  ns) SIM_TIME_NS=$TIME_VALUE ;;
  ps) SIM_TIME_NS=$(echo "$TIME_VALUE / 1000" | bc -l) ;;
  us) SIM_TIME_NS=$(echo "$TIME_VALUE * 1000" | bc -l) ;;
  ms) SIM_TIME_NS=$(echo "$TIME_VALUE * 1000000" | bc -l) ;;
  *)  echo "Error: Unknown time unit '$TIME_UNIT'. Use ns, ps, us, or ms."
      exit 1 ;;
esac

TOTAL_STEPS=$(awk "BEGIN { printf \"%d\", $SIM_TIME_NS * 1000000 / 2 }")

echo "Preparing to run MD simulation for $SIM_TIME ($SIM_TIME_NS ns, $TOTAL_STEPS steps)"

# Prep output dir
mkdir -p "$OUTPUT_DIR"
cp "$INPUT_PDB" "$OUTPUT_DIR/"
cd "$OUTPUT_DIR" || exit 1

# 1. Generate topology
echo "1. Preparing protein structure..."
gmx pdb2gmx -f "$BASENAME.pdb" -o "${BASENAME}_processed.gro" -p "${BASENAME}_topol.top" -water "$WATER_MODEL" -ff "$FORCE_FIELD" -ignh -ter || exit 1

# 2. Define box
echo "2. Defining simulation box..."
gmx editconf -f "${BASENAME}_processed.gro" -o "${BASENAME}_box.gro" -c -d 1.0 -bt cubic || exit 1

# 3. Solvate
echo "3. Adding water..."
gmx solvate -cp "${BASENAME}_box.gro" -cs spc216.gro -o "${BASENAME}_solv.gro" -p "${BASENAME}_topol.top" || exit 1

# 4. Add ions
echo "4. Adding ions..."
cat > ions.mdp << EOF
integrator = steep
emtol = 1000.0
nsteps = 500
EOF

gmx grompp -f ions.mdp -c "${BASENAME}_solv.gro" -p "${BASENAME}_topol.top" -o "${BASENAME}_ions.tpr" || exit 1
echo 13 | gmx genion -s "${BASENAME}_ions.tpr" -o "${BASENAME}_solv_ions.gro" -p "${BASENAME}_topol.top" -pname NA -nname CL -neutral || exit 1

# 5. Energy minimization
echo "5. Energy minimization..."
cat > em.mdp << EOF
integrator  = steep
nsteps      = 50000
emtol       = 1000.0
emstep      = 0.01
nstlist     = 1
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF

gmx grompp -f em.mdp -c "${BASENAME}_solv_ions.gro" -p "${BASENAME}_topol.top" -o "${BASENAME}_em.tpr" || exit 1
gmx mdrun -v -deffnm "${BASENAME}_em" || exit 1

# 6. NVT Equilibration
echo "6. NVT equilibration..."
cat > nvt.mdp << EOF
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
continuation= no
constraint_algorithm = lincs
constraints = all-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
DispCorr    = EnerPres
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300
gen_vel     = yes
EOF

gmx grompp -f nvt.mdp -c "${BASENAME}_em.gro" -r "${BASENAME}_em.gro" -p "${BASENAME}_topol.top" -o "${BASENAME}_nvt.tpr" || exit 1
gmx mdrun -v -deffnm "${BASENAME}_nvt" || exit 1

# 7. NPT Equilibration
echo "7. NPT equilibration..."
cat > npt.mdp << EOF
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
continuation= yes
constraint_algorithm = lincs
constraints = all-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
DispCorr    = EnerPres
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
EOF

gmx grompp -f npt.mdp -c "${BASENAME}_nvt.gro" -r "${BASENAME}_nvt.gro" -t "${BASENAME}_nvt.cpt" -p "${BASENAME}_topol.top" -o "${BASENAME}_npt.tpr" || exit 1
gmx mdrun -v -deffnm "${BASENAME}_npt" || exit 1

# 8. Production MD
echo "8. Running production MD..."
cat > md.mdp << EOF
integrator  = md
nsteps      = $TOTAL_STEPS
dt          = 0.002
nstxout     = 5000
nstvout     = 5000
nstenergy   = 5000
nstlog      = 5000
nstxout-compressed = 5000
continuation= yes
constraint_algorithm = lincs
constraints = all-bonds
lincs_iter  = 1
lincs_order = 4
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
DispCorr    = EnerPres
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
gen_vel     = no
EOF

gmx grompp -f md.mdp -c "${BASENAME}_npt.gro" -t "${BASENAME}_npt.cpt" -p "${BASENAME}_topol.top" -o "${BASENAME}_md.tpr" || exit 1
gmx mdrun -v -deffnm "${BASENAME}_md" || exit 1

echo "✅ MD simulation completed. Results saved in $OUTPUT_DIR"

