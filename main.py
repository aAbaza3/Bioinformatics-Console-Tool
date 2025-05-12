# Import necessary modules
from Bio import SeqIO  # Used to read FASTA or GenBank sequence files
from collections import Counter  # Used to count nucleotide frequency

def load_sequences(file_path, file_format):
    """
    Load sequences from a file (FASTA or GenBank).
    Args:
        file_path (str): Path to the input file.
        file_format (str): Format of the file ('fasta' or 'genbank').
    Returns:
        list: A list of SeqRecord objects.
    """
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        if not sequences:
            print(" No sequences found in the file.")
        return sequences
    except Exception as e:
        print(f"Error while reading the file: {e}")
        return []

def validate_sequence(seq, seq_type):
    """
    Validate if a sequence contains only valid DNA or RNA nucleotides.
    Args:
        seq (str): The sequence to validate.
        seq_type (str): 'DNA' or 'RNA'.
    Returns:
        bool: True if valid, False otherwise.
    """
    valid_dna = {'A', 'T', 'C', 'G'}
    valid_rna = {'A', 'U', 'C', 'G'}
    seq = seq.upper()
    if seq_type.upper() == 'DNA':
        return set(seq).issubset(valid_dna)
    elif seq_type.upper() == 'RNA':
        return set(seq).issubset(valid_rna)
    else:
        print("Unknown sequence type. Please enter 'DNA' or 'RNA'.")
        return False

def count_nucleotide_frequency(seq):
    """
    Count the frequency of each nucleotide in the sequence.
    Args:
        seq (str): The nucleotide sequence.
    Returns:
        dict: A dictionary with nucleotides as keys and counts as values.
    """
    seq = seq.upper()
    return dict(Counter(seq))

def transcribe_dna_to_rna(seq):
    """
    Transcribe a DNA sequence into RNA by replacing T with U.
    Args:
        seq (str): DNA sequence.
    Returns:
        str: RNA sequence.
    """
    return seq.upper().replace('T', 'U')

def main():
    """
    Main function providing a simple console-based user interface
    for loading sequences and performing bioinformatics operations.
    """
    print("=== Bioinformatics Console Tool ===")
    file_path = input("Enter the file path (e.g., data/sample.fasta): ").strip()
    file_format = input(" Enter the file format (fasta or genbank): ").strip().lower()

    sequences = load_sequences(file_path, file_format)
    if not sequences:
        return

    while True:
        print("\n Options:")
        print("1. Validate sequence")
        print("2. Count nucleotide frequency")
        print("3. Transcribe DNA to RNA")
        print("4.  Exit")

        choice = input("Enter your choice: ").strip()

        if choice == '1':
            seq_type = input("Enter sequence type (DNA or RNA): ").strip()
            for record in sequences:
                is_valid = validate_sequence(str(record.seq), seq_type)
                print(f" Sequence '{record.id}' is {'valid' if is_valid else 'invalid'}.")

        elif choice == '2':
            for record in sequences:
                freq = count_nucleotide_frequency(str(record.seq))
                print(f" Nucleotide frequency in '{record.id}': {freq}")

        elif choice == '3':
            for record in sequences:
                rna_seq = transcribe_dna_to_rna(str(record.seq))
                print(f" RNA sequence for '{record.id}': {rna_seq}")

        elif choice == '4':
            print("Thank you for using the program. Goodbye!")
            break

        else:
            print(" Invalid choice. Please try again.")

# Run the program
if __name__ == "__main__":
    main()
