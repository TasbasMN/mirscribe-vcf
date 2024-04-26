BEGIN {
    # Set the input field separator (FS) to "-"
    # Set the output field separator (OFS) to ","
    FS = "-"
    OFS = ","
}

# Process lines that start with ">"
/^>/ {
    # Remove the leading ">" from the first field
    sub(/^>/, "", $1)
    
    # Store the first and second columns in the "metadata" variable
    metadata = $1 OFS $2
    
    # Store the third column in the "thirdColumn" variable
    thirdColumn = $3
    
    # Read the next line (details line) into the current record
    getline
    
    # Replace all occurrences of "&" with "," in the current record
    gsub(/&/, ",")
    
    # Remove all occurrences of ":" in the current record
    gsub(/:/, "")
    
    # Replace multiple consecutive spaces with a single comma in the current record
    gsub(/  +/, ",")
    
    # Split the current record into an array called "arr" using "," as the delimiter
    # Store the number of elements in the array in the variable "n"
    n = split($0, arr, ",")
    
    # Remove any parentheses from the last entry in the "arr" array
    gsub(/[()]/, "", arr[n])
   

    
    # Reconstruct the current record by concatenating the elements of the "arr" array
    # Start with the first element
    $0 = arr[1]
    
    # Loop through the remaining elements (from index 2 to n) and append them to the current record
    # Use the output field separator (OFS) to separate the elements
    for (i = 2; i <= n; i++) {
        $0 = $0 OFS arr[i]
    }

    # Remove all occurrences of " " in the current record (edge case for positive energy values)
    gsub(/ /, "")

    # Print the combined line:
    # - First, print the "metadata" variable
    # - Then, print the reconstructed current record
    # - Finally, print the "thirdColumn" variable
    # Use the output field separator (OFS) to separate the fields
    print metadata OFS $0 OFS thirdColumn
}

