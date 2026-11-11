#!/bin/bash
#SBATCH -p scc-cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 2-0:00:00
#SBATCH --mem=350G
#SBATCH --output=%x.log

module load java

# Define dataset folders
dataset_folders=(roots stems 10DayNodules 14DayNodules 28DayNodules flowers leaves buds seeds WholePlant FC FM RC RM SC SM)

# Define output directory (one level above current directory)
OUTPUT_DIR="../MS_searches"
mkdir -p "$OUTPUT_DIR"

# Define FASTA database URL
FASTA_URL="https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.9/20220708_MtrunA17r5.0-ANR-EGN-r1.9_Annotation.zip"
FASTA_ZIP="fasta_database.zip"
FASTA_FILE="MtrunA17r5.0-ANR-EGN-r1.9.prot.fasta"

# Download and extract FASTA database
if [ ! -f "$FASTA_FILE" ]; then
    wget -O "$FASTA_ZIP" "$FASTA_URL"
    unzip "$FASTA_ZIP" "$FASTA_FILE"
    rm "$FASTA_ZIP"
fi

# Define SearchGUI and PeptideShaker download URLs
SEARCHGUI_URL="https://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/4.3.11/SearchGUI-4.3.11-mac_and_linux.tar.gz"
PEPTIDESHAKER_URL="https://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/3.0.11/PeptideShaker-3.0.11.zip"

# Define installation directories
SOFTWARE_DIR="SOFTWARE"
SEARCHGUI_DIR="$SOFTWARE_DIR/SearchGUI"
PEPTIDESHAKER_DIR="$SOFTWARE_DIR/PeptideShaker"
mkdir -p "$SEARCHGUI_DIR" "$PEPTIDESHAKER_DIR"

# Download and extract SearchGUI
if [ ! -f "$SEARCHGUI_DIR/SearchGUI.jar" ]; then
    wget -O "$SOFTWARE_DIR/SearchGUI.tar.gz" "$SEARCHGUI_URL"
    tar -xzf "$SOFTWARE_DIR/SearchGUI.tar.gz" -C "$SEARCHGUI_DIR" --strip-components=1
    rm "$SOFTWARE_DIR/SearchGUI.tar.gz"
fi

# Download and extract PeptideShaker
if [ ! -f "$PEPTIDESHAKER_DIR/PeptideShaker.jar" ]; then
    wget -O "$SOFTWARE_DIR/PeptideShaker.zip" "$PEPTIDESHAKER_URL"
    unzip "$SOFTWARE_DIR/PeptideShaker.zip" -d "$PEPTIDESHAKER_DIR"
    rm "$SOFTWARE_DIR/PeptideShaker.zip"
fi

# Define paths to SearchGUI and PeptideShaker jars
SEARCHGUI_JAR="$SEARCHGUI_DIR/SearchGUI.jar"
PEPTIDESHAKER_JAR="$PEPTIDESHAKER_DIR/PeptideShaker.jar"

# Loop through dataset folders
for folder in "${dataset_folders[@]}"; do
    if [ -d "$folder" ]; then
        MGF_FILES=("$folder"/*.mgf)
        if [ -e "${MGF_FILES[0]}" ]; then
            echo "Processing $folder..."
            
            # Define output directories for SearchGUI and PeptideShaker
            SEARCHGUI_OUTPUT="$OUTPUT_DIR/$folder/SearchGUI"
            PEPTIDESHAKER_OUTPUT="$OUTPUT_DIR/$folder/PeptideShaker"
            mkdir -p "$SEARCHGUI_OUTPUT" "$PEPTIDESHAKER_OUTPUT"

            # Run SearchGUI
            java -Xms64G -Xmx300G -jar "$SEARCHGUI_JAR" -spectrum_files "$folder"/ -output_folder "$SEARCHGUI_OUTPUT" -id_params search_parameters.par -db "$FASTA_FILE"

            # Run PeptideShaker
            java -Xms64G -Xmx300G -jar "$PEPTIDESHAKER_JAR" -experiment "$folder" -sample "$folder" -replicate 1 -identification_files "$SEARCHGUI_OUTPUT" -out "$PEPTIDESHAKER_OUTPUT/$folder.cpsx"
        fi
    fi
done

echo "Done"
