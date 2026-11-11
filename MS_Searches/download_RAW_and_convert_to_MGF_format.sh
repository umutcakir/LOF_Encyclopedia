#!/bin/bash
#SBATCH -p scc-cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 2-0:00:00
#SBATCH --mem=350G
#SBATCH --output=%x.log

module load mono

# Define the base URLs
BASE_URLS=(
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2016/10/PXD002692/"
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/11/PXD013606/"
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/07/PXD022278/"
)

CHECKSUM_URLS=(
    "https://www.ebi.ac.uk/pride/ws/archive/v3/files/checksum/PXD002692"
    "https://www.ebi.ac.uk/pride/ws/archive/v3/files/checksum/PXD013606"
    "https://www.ebi.ac.uk/pride/ws/archive/v3/files/checksum/PXD022278"
)

FILTER_PATTERN="_CEM_.*\.raw|_ALR_.*\.raw|MedicagoCompleteProteome_WholePlant\.raw|FC[0-9]*\.raw|FM[0-9]*\.raw|RC[0-9]*\.raw|RM[0-9]*\.raw|SC[0-9]*\.raw|SM[0-9]*\.raw"

# Create directories
mkdir -p RAW MGF SOFTWARE roots stems 10DayNodules 14DayNodules 28DayNodules flowers leaves buds seeds WholePlant FC FM RC RM SC SM

# Download and unzip ThermoRawFileParser
wget -O SOFTWARE/ThermoRawFileParser.zip "https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.5/ThermoRawFileParser1.4.5.zip"
unzip SOFTWARE/ThermoRawFileParser.zip -d SOFTWARE/

# Process each dataset
for i in "${!BASE_URLS[@]}"; do
    BASE_URL=${BASE_URLS[$i]}
    CHECKSUM_URL=${CHECKSUM_URLS[$i]}
    
    # Get the list of all files and filter for the desired patterns
    wget -q -O - "$BASE_URL" | grep -Eo 'href="[^"]+"' | cut -d'"' -f2 | grep -E "$FILTER_PATTERN" > raw_files.txt
    
    echo "Downloading, $filename" 
    # Download filtered files
    while read -r filename; do
        wget -q -P RAW "$BASE_URL$filename"
    done < raw_files.txt

    # Verify downloaded files
    if [ -z "$(ls -A RAW/*.raw 2>/dev/null)" ]; then
        echo "No matching files were downloaded from $BASE_URL."
        continue
    fi
    
    echo "Download from $BASE_URL completed successfully."

    # Download official checksum file
    wget -O RAW/checksums_raw.txt "$CHECKSUM_URL"

    # Extract only filename and MD5 checksum for files containing the filter pattern
    grep -E "$FILTER_PATTERN" RAW/checksums_raw.txt | awk -F'\t' '{print $2, $1}' > RAW/checksums.md5

    # Verify checksums against official list
    cd RAW || exit
    md5sum -c checksums.md5 >> checksum_verification.log 2>&1

    if grep -q "FAILED" checksum_verification.log; then
        echo "Warning: Some files failed checksum verification! Check checksum_verification.log for details."
    else
        echo "All files from $BASE_URL passed checksum verification."
    fi
    cd ..

    # Move files into appropriate folders
    for file in RAW/*.raw; do
        if [[ "$file" == *"roots"* ]]; then
            mv "$file" roots/
        elif [[ "$file" == *"stems"* ]]; then
            mv "$file" stems/
        elif [[ "$file" == *"10DayNodules"* ]]; then
            mv "$file" 10DayNodules/
        elif [[ "$file" == *"14DayNodules"* ]]; then
            mv "$file" 14DayNodules/
        elif [[ "$file" == *"28DayNodules"* ]]; then
            mv "$file" 28DayNodules/
        elif [[ "$file" == *"flowers"* ]]; then
            mv "$file" flowers/
        elif [[ "$file" == *"leaves"* ]]; then
            mv "$file" leaves/
        elif [[ "$file" == *"buds"* ]]; then
            mv "$file" buds/
        elif [[ "$file" == *"seeds"* ]]; then
            mv "$file" seeds/
        elif [[ "$file" == *"MedicagoCompleteProteome_WholePlant"* ]]; then
            mv "$file" WholePlant/
        elif [[ "$file" =~ FC[0-9]*\.raw ]]; then
            mv "$file" FC/
        elif [[ "$file" =~ FM[0-9]*\.raw ]]; then
            mv "$file" FM/
        elif [[ "$file" =~ RC[0-9]*\.raw ]]; then
            mv "$file" RC/
        elif [[ "$file" =~ RM[0-9]*\.raw ]]; then
            mv "$file" RM/
        elif [[ "$file" =~ SC[0-9]*\.raw ]]; then
            mv "$file" SC/
        elif [[ "$file" =~ SM[0-9]*\.raw ]]; then
            mv "$file" SM/
        fi
    done

done

echo "Download, integrity check, and sorting complete."

# Convert RAW files to MGF format using ThermoRawFileParser
for folder in roots stems 10DayNodules 14DayNodules 28DayNodules flowers leaves buds seeds WholePlant FC FM RC RM SC SM; do
    mono SOFTWARE/ThermoRawFileParser.exe -d="$folder/" -f=0
done

echo "RAW to MGF conversion complete."

# Move all .raw files back to RAW folder after conversion
for folder in roots stems 10DayNodules 14DayNodules 28DayNodules flowers leaves buds seeds WholePlant FC FM RC RM SC SM; do
    mv "$folder"/*.raw RAW/
done

echo "All .raw files moved back to RAW folder."

