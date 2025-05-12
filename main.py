# Import necessary modules
from Bio import SeqIO
from collections import Counter

def load_sequences(file_path, file_format):
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        if not sequences:
            print("No sequences found in the file.")
        return sequences
    except Exception as e:
        print(f"Error while reading the file: {e}")
        return []

def validate_sequence(seq, seq_type):
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
    return dict(Counter(seq.upper()))

def transcribe_dna_to_rna(seq):
    return seq.upper().replace('T', 'U')

def calculate_hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        print("Sequences must be of equal length to compute Hamming Distance.")
        return None
    return sum(base1 != base2 for base1, base2 in zip(seq1.upper(), seq2.upper()))

def main():
    print("=== Bioinformatics Console Tool ===")
    file_path = input("Enter the file path (e.g., data/sample.fasta): ").strip()
    file_format = input("Enter the file format (fasta or genbank): ").strip().lower()

    sequences = load_sequences(file_path, file_format)
    if not sequences:
        return

    while True:
        print("\nOptions:")
        print("1. Validate sequence")
        print("2. Count nucleotide frequency")
        print("3. Transcribe DNA to RNA")
        print("4. Calculate Hamming Distance")
        print("5. Exit")

        choice = input("Enter your choice: ").strip()

        if choice == '1':
            seq_type = input("Enter sequence type (DNA or RNA): ").strip()
            for record in sequences:
                is_valid = validate_sequence(str(record.seq), seq_type)
                print(f"Sequence '{record.id}' is {'valid' if is_valid else 'invalid'}.")

        elif choice == '2':
            for record in sequences:
                freq = count_nucleotide_frequency(str(record.seq))
                print(f"Nucleotide frequency in '{record.id}': {freq}")

        elif choice == '3':
            for record in sequences:
                rna_seq = transcribe_dna_to_rna(str(record.seq))
                print(f"RNA sequence for '{record.id}': {rna_seq}")

        elif choice == '4':
            file1 = input("Enter path to first sequence file: ").strip()
            file2 = input("Enter path to second sequence file: ").strip()
            file_format1 = input("Enter format for first file (fasta or genbank): ").strip().lower()
            file_format2 = input("Enter format for second file (fasta or genbank): ").strip().lower()

            sequences1 = load_sequences(file1, file_format1)
            sequences2 = load_sequences(file2, file_format2)

            if not sequences1 or not sequences2:
                print("Failed to load sequences from one or both files.")
                continue

            seq1 = str(sequences1[0].seq)
            seq2 = str(sequences2[0].seq)

            distance = calculate_hamming_distance(seq1, seq2)
            if distance is not None:
                print(f"Hamming Distance between '{sequences1[0].id}' and '{sequences2[0].id}': {distance}")

        elif choice == '5':
            print("Thank you for using the program. Goodbye!")
            break

        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main()
