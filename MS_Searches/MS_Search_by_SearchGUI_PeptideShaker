#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 2  # Only need 2 cores for job submission
#SBATCH -t 1-0:00:00  # Enough time for job submission
#SBATCH --mem=50G  # Memory for job submission
#SBATCH --output=Job_Commands/run_all.log  # Store master log in Job_Commands

source activate MS_Searches

# Enable nullglob to handle cases where no .mgf files exist
shopt -s nullglob  

### **Step 1: Set Up Environment**
echo "Setting up FASTA environment..."

FASTA_FILES_DIR="FASTA_FILES"
FASTA_DB_DIR="FASTA_DB"
SOFTWARE_MASTER="Software_Master"

# Create necessary directories
mkdir -p "$FASTA_DB_DIR" "Job_Commands" "$SOFTWARE_MASTER"

# Get the list of FASTA files
FASTA_FILES=("$FASTA_FILES_DIR"/*.fasta)

# Download and extract SearchGUI into Software_Master (if not already downloaded)
if [ ! -f "$SOFTWARE_MASTER/SearchGUI/SearchGUI-4.3.15/SearchGUI-4.3.15.jar" ]; then
    echo "Downloading SearchGUI..."
    mkdir -p "$SOFTWARE_MASTER/SearchGUI"
    wget -q -O "$SOFTWARE_MASTER/SearchGUI.tar.gz" "https://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/4.3.15/SearchGUI-4.3.15-mac_and_linux.tar.gz"
    tar -xzf "$SOFTWARE_MASTER/SearchGUI.tar.gz" -C "$SOFTWARE_MASTER/SearchGUI"
fi

# Download and extract PeptideShaker into Software_Master (if not already downloaded)
if [ ! -f "$SOFTWARE_MASTER/PeptideShaker/PeptideShaker-3.0.11/PeptideShaker-3.0.11.jar" ]; then
    echo "Downloading PeptideShaker..."
    mkdir -p "$SOFTWARE_MASTER/PeptideShaker"
    wget -q -O "$SOFTWARE_MASTER/PeptideShaker.zip" "https://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/3.0.11/PeptideShaker-3.0.11.zip"
    unzip -q "$SOFTWARE_MASTER/PeptideShaker.zip" -d "$SOFTWARE_MASTER/PeptideShaker"
fi

echo "FASTA database and software setup complete."

### **Step 2: Create Decoy Databases for Each FASTA File**
for fasta_file in "${FASTA_FILES[@]}"; do
    fasta_basename=$(basename "$fasta_file" .fasta)
    expected_decoy_fasta="${FASTA_DB_DIR}/${fasta_basename}_concatenated_target_decoy.fasta"

    # Generate decoy FASTA if it doesn't exist
    if [ ! -f "$expected_decoy_fasta" ]; then
        echo "Generating decoy FASTA for $fasta_basename..."
        java -Xmx128G -Djava.awt.headless=true -cp "$SOFTWARE_MASTER/SearchGUI/SearchGUI-4.3.15/SearchGUI-4.3.15.jar" \
            eu.isas.searchgui.cmd.FastaCLI -in "$fasta_file" -decoy

        # Rename the generated decoy FASTA to match the expected name
        mv "$FASTA_FILES_DIR/${fasta_basename}_concatenated_target_decoy.fasta" "$expected_decoy_fasta"
    fi
done

### **Step 3: Submit Parallel Jobs for Each FASTA and Each Dataset**
dataset_folders=(roots stems 10DayNodules 14DayNodules 28DayNodules flowers leaves buds seeds WholePlant FC FM RC RM SC SM)

for fasta_file in "${FASTA_FILES[@]}"; do
    fasta_basename=$(basename "$fasta_file" .fasta)
    decoy_fasta="${FASTA_DB_DIR}/${fasta_basename}_concatenated_target_decoy.fasta"

    for folder in "${dataset_folders[@]}"; do
        job_name="MS_${fasta_basename}_${folder}"
        job_script="Job_Commands/job_${fasta_basename}_${folder}.sh"
        log_file="Job_Commands/job_${fasta_basename}_${folder}.log"

        # Write job script
        cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 14-0:00:00
#SBATCH --mem=550G
#SBATCH --output=${log_file}
#SBATCH --job-name=$job_name

sleep $(( RANDOM % 5 + 1 ))

source activate MS_Searches

cd .. # This is necessary if you are not sending jobs automatically 

WORKDIR=$(pwd)/

# Define job-specific variables
folder="$folder"
fasta_basename="$fasta_basename"
FASTA_COPY="\$folder/\${fasta_basename}_concatenated_target_decoy.fasta"
LOG_FILE="${log_file}"

echo "Processing dataset folder: \$folder with FASTA: \$fasta_basename" | tee -a "\$LOG_FILE"

# Define paths
SOFTWARE_DIR="\$folder/Software"
SEARCHGUI_DIR="\$SOFTWARE_DIR/SearchGUI"
PEPTIDESHAKER_DIR="\$SOFTWARE_DIR/PeptideShaker"

# Copy SearchGUI and PeptideShaker from Software_Master to the dataset folder
mkdir -p "\$SOFTWARE_DIR"
cp -r "Software_Master/SearchGUI" "\$SEARCHGUI_DIR"
cp -r "Software_Master/PeptideShaker" "\$PEPTIDESHAKER_DIR"

# Define JAR paths
SEARCHGUI_JAR="\$SEARCHGUI_DIR/SearchGUI-4.3.15/SearchGUI-4.3.15.jar"
PEPTIDESHAKER_JAR="\$PEPTIDESHAKER_DIR/PeptideShaker-3.0.11/PeptideShaker-3.0.11.jar"

# Define output folders
SEARCHGUI_OUTPUT="\$folder/SearchGUI"
PEPTIDESHAKER_OUTPUT="\$folder/PeptideShaker"

# Ensure output directories exist
mkdir -p "\$SEARCHGUI_OUTPUT"
mkdir -p "\$PEPTIDESHAKER_OUTPUT"

# Copy FASTA file to dataset folder
cp "$decoy_fasta" "\$FASTA_COPY"

# Unique output name
OUTPUT_NAME="\${fasta_basename}_\${folder}"

# Full SearchGUI command
SEARCHGUI_CMD="java -Xmx400G -Djava.awt.headless=true -cp \"\$SEARCHGUI_JAR\" eu.isas.searchgui.cmd.SearchCLI \
    -spectrum_files \"../RAW_files/\$folder/\" \
    -fasta_file \"\$WORKDIR\$FASTA_COPY\" \
    -output_folder \"\$SEARCHGUI_OUTPUT/\" \
    -id_params \"FASTA_FILES/020225.par\" \
    -xtandem 1 -msgf 1 -omssa 1 -comet 1 -output_data 1 -output_default_name \"\$OUTPUT_NAME\""

echo "Running SearchGUI command:" | tee -a "\$LOG_FILE"
echo "\$SEARCHGUI_CMD" | tee -a "\$LOG_FILE"
eval \$SEARCHGUI_CMD | tee -a "\$LOG_FILE"

#mkdir "$PEPTIDESHAKER_OUTPUT/Reports/" 

# Full PeptideShaker command
PEPTIDESHAKER_CMD="java -Xmx400G -cp \"\$PEPTIDESHAKER_JAR\" eu.isas.peptideshaker.cmd.PeptideShakerCLI \
    -reference \"\$OUTPUT_NAME\" \
    -identification_files \"\$SEARCHGUI_OUTPUT/\$OUTPUT_NAME.zip\" \
    -zip \"\$PEPTIDESHAKER_OUTPUT/\${OUTPUT_NAME}-PeptideShaker.zip\" \
    -out \"\$PEPTIDESHAKER_OUTPUT/\${OUTPUT_NAME}-PeptideShaker.psdb\" \
    -reports \"0,6,9\""

echo "Running PeptideShaker command:" | tee -a "\$LOG_FILE"
echo "\$PEPTIDESHAKER_CMD" | tee -a "\$LOG_FILE"
eval \$PEPTIDESHAKER_CMD | tee -a "\$LOG_FILE"
unzip -j -o "\$PEPTIDESHAKER_OUTPUT/\${OUTPUT_NAME}-PeptideShaker.zip" "reports/*" -d "\$PEPTIDESHAKER_OUTPUT/Reports" | tee -a "\$LOG_FILE"
sleep 180

EOF

        chmod +x "$job_script"
        #sbatch "$job_script"
    done
done

echo "Done"

#    -xtandem 1 -msgf 1 -omssa 1 -comet 1 -output_data 1 -output_default_name \"\$OUTPUT_NAME\""
 
