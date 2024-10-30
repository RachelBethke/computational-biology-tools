import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import random
from . import core

# Rest of your code stays the same
def load_haplotype_data(file_path):
    """
    Loads haplotype data from a specified text file.
    The file is assumed to contain blocks of haplotypes separated by double newlines.

    Args: 
        file_path (str): The path to the haplotype data file.
    
    Returns: pd.DataFrame: A DataFrame containing the processed haplotype data.
    """
    with open(file_path, 'r') as file:
        blocks = file.read().split('\n\n')  # Split file into blocks
        block = random.choice(blocks)  # Select one block at random
        processed_lines = []
        for line in block.strip().split('\n'):
            stripped_line = line.strip()
            if stripped_line.replace(",", "").isdigit():
                processed_lines.append(stripped_line.split(','))
        return pd.DataFrame(processed_lines).astype(int) 

# Process the file and display results
def process_file(haplotype_data):
    """
    Processes the haplotype data by calculating diversity metrics and updating the UI with the results.
    
    Args: 
        haplotype_data (pd.DataFrame): The processed haplotype data.
    """
    try:
        pi = core.calc_pi(haplotype_data)
        watts_theta = core.calc_watterson(haplotype_data)
        tajima_d = core.tajimas_d(haplotype_data)
        
        # Update result display
        result_text = (
            f"Nucleotide Diversity (π): {pi:.4f}\n"
            f"Watterson's Theta: {watts_theta:.4f}\n"
            f"Tajima's D: {tajima_d:.4f}"
        )
        result_label.config(text=result_text)
        
        # Create plot
        create_plot(haplotype_data)
        
    except Exception as e:
        handle_file_error(e)

def create_plot(haplotype_data):
    """
    Creates and displays a plot of the haplotype data.
    
    Args: 
        haplotype_data (pd.DataFrame): The haplotype data to be plotted.
    """
    # Clear any existing plot
    for widget in plot_frame.winfo_children():
        widget.destroy()
    
    # Create new plot
    fig, ax = plt.subplots(figsize=(8, 4))
    frequencies = haplotype_data.mean()
    ax.bar(range(len(frequencies)), frequencies)
    ax.set_title('Allele Frequencies')
    ax.set_xlabel('Position')
    ax.set_ylabel('Frequency')
    
    # Add plot to GUI
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill='both', expand=True)

# Function to handle file loading and error handling
def load_file():
    """
    Opens a file dialog to select a haplotype file, processes the file, and displays the results.
    Calls error handling if the file cannot be processed.
    """
    try:
        file_path = filedialog.askopenfilename(
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if file_path:
            haplotype_data = load_haplotype_data(file_path)
            process_file(haplotype_data)
    except Exception as e:
        handle_file_error(e)

# Function to handle errors during file processing
def handle_file_error(e):
    """
    Handles errors that occur during the file loading process and displays a message box with the error information.
    
    Args:
        e (Exception): The exception object that contains error details.
    """
    messagebox.showerror("Error", f"Failed to load file: {e}")

# Tkinter UI setup
root = tk.Tk()  # Create the main window
root.title("Conservation Genomics Tool")  # Set the window title
root.geometry("800x600")  # Set window size

# Create main container frame
main_frame = tk.Frame(root)
main_frame.pack(expand=True, fill='both', padx=10, pady=10)

# Add a button to upload a haplotype file
upload_button = tk.Button(main_frame, text="Upload Haplotype File", command=load_file)
upload_button.pack(pady=10)  # Add padding to the button

# Add a label to display the results
result_label = tk.Label(main_frame, text="Results will be displayed here")
result_label.pack(pady=20)  # Add padding to the label

# Add frame for plotting
plot_frame = tk.Frame(main_frame)
plot_frame.pack(expand=True, fill='both', pady=10)

def main():
    root.mainloop()

if __name__ == "__main__":
    main()