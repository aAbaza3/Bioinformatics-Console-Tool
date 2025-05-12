
````markdown
# 🧬 Bioinformatics Console Tool

This is a simple Python console application designed to help users analyze DNA and RNA sequences. The tool reads biological sequence data from common file formats such as FASTA and GenBank, then offers several analysis options via a user-friendly command-line interface.

## ✅ Features

- Read sequence data from **FASTA** or **GenBank** files
- Validate DNA or RNA sequences
- Count nucleotide frequency
- Transcribe DNA into RNA

## 📁 Supported File Formats

- `.fasta`  
- `.fa`  
- `.gb` or `.gbk` (GenBank)

## 🛠 Requirements

- Python 3.x
- [Biopython](https://biopython.org/) library

Install Biopython using pip:
```bash
pip install biopython
````

## 🚀 How to Run

1. Clone or download this repository
2. Make sure you have a valid `.fasta` or `.gbk` file in your project directory
3. Run the Python script:

```bash
python bio_tool.py
```

4. Follow the interactive menu to select analysis options

## 📋 Example Menu

```
Options:
1. Validate sequence
2. Count nucleotide frequency
3. Transcribe DNA to RNA
4. Exit
```

## 📌 Notes

* DNA is considered valid if it only contains A, T, C, G.
* RNA is considered valid if it only contains A, U, C, G.
* For transcription, T is replaced by U to simulate the process of forming RNA from a DNA strand.

## 📧 Author

* Project created by a student learning bioinformatics and Python 💻🧬


