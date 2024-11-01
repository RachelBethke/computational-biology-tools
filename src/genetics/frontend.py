import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import random
from . import core

def load_haplotype_data(file_path):
    """
    Loads haplotype data from a specified text file.
    The file is assumed to contain blocks of haplotypes separated by double newlines.

    Args: 
        file_path (str): The path to the haplotype data file.
    
    Returns: tuple: (pd.DataFrame, list) - The processed haplotype data and list of all blocks
    """
    with open(file_path, 'r') as file:
        blocks = file.read().split('\n\n')  # Split file into blocks
        return blocks
    
def process_block(block):
    """
    Process a single block of haplotype data.

    Args:
        block (str): String containing the block data
    
    Returns:
        pd.DataFrame: Processed haplotype data
    """
    processed_lines = []
    for line in block.strip().split('\n'):
        stripped_line = line.strip()
        if stripped_line.replace(",", "").isdigit():
            processed_lines.append(stripped_line.split(','))
    
    if not processed_lines:  #check if we have any valid data
        raise ValueError("No valid haplotype data was found in selected block")
        
    df = pd.DataFrame(processed_lines).astype(int)
    if df.empty:
        raise ValueError("Empty haplotype data after processing")
    
    return df

def block_choice(blocks):
    """
    Handles the selection of a block from the dropdown menu
    """
    selection = block_var.get()
    if not hasattr(block_choice, 'blocks'):
        return
        
    if selection == "Random":
        block = random.choice(block_choice.blocks)
    else:
        block_idx = int(selection.split()[1]) - 1
        block = block_choice.blocks[block_idx]
    
    haplotype_data = process_block(block)
    process_file(haplotype_data)

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

    plt.close('all') # Deal with memory leakage
    
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
            # Load and store blocks
            blocks = load_haplotype_data(file_path)
            block_choice.blocks = blocks
            
            # Update dropdown menu
            block_options = ["Random"] + [f"Block {i+1}" for i in range(len(blocks)) if blocks[i].strip()]
            block_dropdown['values'] = block_options
            block_dropdown['state'] = 'readonly'
            block_var.set("Random")
            
            # Process initial block - use random block to start
            initial_block = random.choice(blocks)
            haplotype_data = process_block(initial_block)
            process_file(haplotype_data)
            
    except Exception as e:
        handle_file_error(e)

# Function to handle errors during file processing
# TODO: Make the message more specific, when you upload wrong files it says division by zero
def handle_file_error(e):
    """
    Handles errors that occur during the file loading process and displays 
    a message box with the error information.
    
    Args:
        e (Exception): The exception object that contains error details.
    """
    if "division by zero" in str(e):
        error_msg = "Invalid block selection: No valid haplotype data found"
    else:
        error_msg = f"Failed to process file: {str(e)}"
    messagebox.showerror("Error", error_msg)

#Tkinter UI setup
root = tk.Tk()  # makes the main window
root.title("Conservation Genomics Tool")  #window title
root.geometry("800x600")  # window size

def show(): 
    tk.label.config( text = tk.licked.get() ) 

# Create main container frame
main_frame = tk.Frame(root)
main_frame.pack(expand=True, fill='both', padx=10, pady=10)

# Add a button to upload a haplotype file
upload_button = tk.Button(main_frame, text="Upload Haplotype File", command=load_file)
upload_button.pack(pady=10)  # Add padding to the button

block_frame = tk.Frame(main_frame)
block_frame.pack(pady=5)

block_var = tk.StringVar(value="Random")
block_dropdown = ttk.Combobox(
    block_frame, 
    textvariable=block_var,
    values=["Random"],
    state='disabled'
)
block_dropdown.bind('<<ComboboxSelected>>', block_choice)
block_dropdown.pack(side=tk.LEFT, padx=5)

result_label = tk.Label(main_frame, text="Results will be displayed here")
result_label.pack(pady=20)  # Add padding to the label

# Add frame for plotting
plot_frame = tk.Frame(main_frame)
plot_frame.pack(expand=True, fill='both', pady=10)

def main():
    root.mainloop()

if __name__ == "__main__":
    main()