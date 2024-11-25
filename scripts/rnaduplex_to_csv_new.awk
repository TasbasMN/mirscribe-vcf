# Script to process RNA structure predictions from FASTA format
# Input format: FASTA header (>ENSG...-MIMAT...-type) followed by sequence and structure prediction
# Output format: ENSG,MIMAT,structure,pos1,pos2,range1,range2,energy,type

BEGIN {
    # Define field separators
    FS = "-"    # For splitting header line
    OFS = ","   # For output format
}

# Process header lines (starting with ">")
/^>/ {
    # Remove ">" from the header
    sub(/^>/, "", $1)
    
    # Store components from header
    metadata = $1 OFS $2
    thirdColumn = $3
    
    # Get the next line (structure prediction line)
    getline
    
    # Replace "&" with "," in structure notation
    gsub(/&/, ",")
    
    # Remove all spaces around plus signs first
    gsub(/\+[ ]+/, "+", $0)
    gsub(/[ ]+\+/, "+", $0)
    
    # Match the line into components using regex
    match($0, /([.()]+[.,][.()]+)[ ]+([0-9,]+)[ ]*:[ ]*([0-9,]+)[ ]*\([ ]*([+-]?[0-9.]+)[ ]*\)/, parts)
    
    if (length(parts) > 0) {
        # Check if energy value is positive (no minus sign or has plus sign)
        if (parts[4] !~ /^-/) {
            parts[4] = "0.00"
        }
        # Construct output with matched components
        $0 = parts[1] OFS parts[2] OFS parts[3] OFS parts[4]
    }
    
    # Print final output
    print metadata OFS $0 OFS thirdColumn
}
