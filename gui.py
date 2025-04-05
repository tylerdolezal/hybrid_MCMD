import tkinter as tk
from tkinter import ttk, filedialog, messagebox

surface_check_var = False

def update_randomize_visibility(*args):
    if grain_boundary_var.get():
        # Randomization
        randomize_label = tk.Label(frame, text="Randomize:").grid(row=2, column=2, sticky="w")
        tk.Checkbutton(frame, variable=randomize_var).grid(row=2, column=3, sticky="w")


def update_hybrid_visibility(*args):
    if hybrid_run_var.get():
        md_params_label.grid(row=12, column=0, columnspan=2, pady=5)
        tk.Label(root, text="MD Steps:").grid(row=13, column=0, sticky="w", padx=10)
        md_steps_entry.grid(row=13, column=1, sticky="w")
        md_interval_label.grid(row=14, column=0, sticky="w", padx=10)
        md_interval_entry.grid(row=14, column=1, sticky="w")

    """Enable or disable MD input fields based on Hybrid Run selection."""
    state = "normal" if hybrid_run_var.get() else "disabled"
    
    # Update the entry fields to be editable or not
    md_steps_entry.config(state=state)
    md_interval_entry.config(state=state)


def clean_value(value):
    """Convert empty lists, tuples, or specific values into 'None'."""
    if isinstance(value, (list, tuple)) and not any(value):  # Check for empty list/tuple
        return "None"
    if value in ["[]", "()", "[()]"]:  # Handle string representations
        return "None"
    return value  # Return value if it's valid

import os
def generate_input_file():
    """Generate or update the input file based on user selections."""
    script_directory = os.getcwd()  # Get the directory where the script is executed
    file_path = os.path.join(script_directory, "input_file")  # Set output file name
    if continue_run_var.get():
        # If "Continue Run" is selected, modify existing input_file
        if os.path.exists(file_path):
            with open(file_path, "r") as file:
                lines = file.readlines()
            
            updated_lines = []
            for line in lines:
                if line.startswith("continue_run:"):
                    updated_lines.append("continue_run: True\n")
                else:
                    updated_lines.append(line)
            
            with open(file_path, "w") as file:
                file.writelines(updated_lines)
            
            messagebox.showinfo("Info", f"Reading inputs from existing input_file!")
            root.destroy()  # Close GUI
        else:
            messagebox.showerror("Error", "input_file does not exist! Run without 'Continue Run' checked first.")
        return
    
    if adsorbed_var.get():
        messagebox.showinfo("Warning", f"Make sure to provide surface as: [(adsorbate, None)] for already adsorbed species.")
    # Define input file content
    md_params_value = f"0, {tsim_var.get()}, {md_potential_var.get()}" if not hybrid_run_var.get() else f"{md_steps_var.get()}, {md_temp_var.get()}, {md_potential_var.get()}"
    
    input_data = {
        "continue_run": continue_run_var.get(),
        "composition": composition_var.get(),
        "crystal_shape": crystal_shape_var.get(),
        "size": size_var.get(),
        "local_swap": local_swap_var.get(),
        "grain_boundary": grain_boundary_var.get(),
        "randomize": randomize_var.get() if grain_boundary_var.get() else False,
        "md_params": md_params_value,
        "md_interval": md_interval_var.get() if hybrid_run_var.get() else str(int(num_mc_steps_var.get())+1),
        "num_mc_steps": num_mc_steps_var.get(),
        "batch_mode": batch_mode_var.get(),
        "additives": additives_var.get(),
        "vacancies": vacancies_var.get(),
        "metal_library": metal_library_var.get(),
        "surface": surface_var.get(),
        "freeeze_threshold": freeze_var.get(),
        "supcomp_command": supcomp_command_var.get()
    }
    
    # Clean up values before writing to file
    input_data = {key: clean_value(value) for key, value in input_data.items()}

    # Write to input_file in the same directory as the script
    with open(file_path, "w") as file:
        file.writelines(f"{key}: {value}\n" for key, value in input_data.items())
    
    # if GB true, check for POSCAR-gb file
    if grain_boundary_var.get():
        if not os.path.exists("POSCAR-gb"):
            messagebox.showerror("Error", "POSCAR-gb file not found! When using GB mode, please ensure the file exists in the working directory.")
            root.destroy()  # Close GUI
        
    root.destroy()  # Close GUI
    messagebox.showinfo("Success", f"input_file generated successfully!")


# Create the main application window
root = tk.Tk()
root.title("Hybrid MCMD Input File Generator")
root.geometry("1200x600")

# Define variables
composition_var = tk.StringVar(value="Ni=0.7, Cr=0.3")
crystal_shape_var = tk.StringVar(value="fcc")
grain_boundary_var = tk.BooleanVar(value=False)
randomize_var = tk.BooleanVar(value=False)
hybrid_run_var = tk.BooleanVar(value=False)
md_steps_var = tk.StringVar(value="1000")
md_temp_var = tk.StringVar(value="1073")
tsim_var = tk.StringVar(value="1073")
md_potential_var = tk.StringVar(value="pfp")
num_mc_steps_var = tk.StringVar(value="6000")
md_interval_var = tk.StringVar(value="100")
size_var = tk.StringVar(value="6")
supcomp_command_var = tk.StringVar(value="lmp_serial")
continue_run_var = tk.BooleanVar(value=False)
additives_var = tk.StringVar(value="[()]")
vacancies_var = tk.BooleanVar(value=False)
metal_library_var = tk.StringVar(value="[]")
surface_var = tk.StringVar(value="[()]")
freeze_var = tk.StringVar(value="0.0")
batch_mode_var = tk.StringVar(value="[]")
local_swap_var = tk.BooleanVar(value=True)
adsorbed_var = tk.BooleanVar(value=False)

# Layout frame
frame = tk.Frame(root)
frame.grid(row=0, column=0, padx=10, pady=10)

# Checkboxes (side-by-side)
tk.Label(frame, text="Continue Run:").grid(row=0, column=0, sticky="w")
tk.Checkbutton(frame, variable=continue_run_var, command=generate_input_file).grid(row=0, column=1, sticky="w")

tk.Label(frame, text="Hybrid Run:").grid(row=1, column=0, sticky="w")
tk.Checkbutton(frame, variable=hybrid_run_var, command=update_hybrid_visibility).grid(row=1, column=1, sticky="w")

tk.Label(frame, text="Vacancies:").grid(row=0, column=2, sticky="w")
tk.Checkbutton(frame, variable=vacancies_var).grid(row=0, column=3, sticky="w")

tk.Label(frame, text="Grain Boundary:").grid(row=2, column=0, sticky="w")
tk.Checkbutton(frame, variable=grain_boundary_var, command=update_randomize_visibility).grid(row=2, column=1, sticky="w")

tk.Label(frame, text="Local Swap:").grid(row=1, column=2, sticky="w")
tk.Checkbutton(frame, variable=local_swap_var).grid(row=1, column=3, sticky="w")

tk.Label(frame, text="Already Adsorbed (surface only)?").grid(row=2, column=2, sticky="w")
tk.Checkbutton(frame, variable=adsorbed_var).grid(row=2, column=3, sticky="w")

preamble_fields = [
    ("Continue Run", continue_run_var, "Picking up where we left off?", "True or False"),
    ("Hybrid Run", hybrid_run_var, "Use MD in MCMD?", "True or False"),
    ("Vacancies", vacancies_var, "Place vacancies in the lattice?", "True or False"),
    ("Grain Boundary", grain_boundary_var, "Using a predefined grain boundary cell?", "True or False"),
    ("Randomize", randomize_var, "Shuffle chemical symbols based on composition?", "True or False"),
    ("Local Swap", local_swap_var, "Swap metals locally?", "True or False"),
    ("Adsorbed?", adsorbed_var, "Adsrobate already relaxed?", "True or False")]

# MD Parameters
md_params_label = tk.Label(root, text="MD Parameters:")
md_steps_entry = tk.Entry(root, textvariable=md_steps_var, width=10)
md_temp_entry = tk.Entry(root, textvariable=tsim_var, width=10)
md_potential_entry = tk.Entry(root, textvariable=md_potential_var, width=10)

md_interval_label = tk.Label(root, text="MD Interval:")
md_interval_entry = tk.Entry(root, textvariable=md_interval_var, width=10)

# Additional input fields
fields = [
    ("Composition", composition_var, "Composition of the undoped simulation lattice", "Ni=0.7, Cr=0.3"),
    ("Crystal Shape", crystal_shape_var, "Crystal type of the undoped simulation lattice", "fcc or bcc or hcp"),
    ("Superlattice Multiplier", size_var, "Multiply lattice by (m, m, m) to build superlattice", "6"),
    ("Tsim (K)", tsim_var, "Simulation temperature in Kelvin", "1073"),
    ("MD Potential", md_potential_var, "Potential type used with LAMMPS (or CHGnet)", "pfp or meam or eam or chgnet"),
    ("MC Steps", num_mc_steps_var, "Number of Monte Carlo steps", "6000"),
    ("Additives", additives_var, "Dopants or interstitials added to the system", "[(Cr, B, 9), (None, B, 9), etc] or None"),
    ("Batch Mode", batch_mode_var, "Iterate compositions for multiple simulations", "[0.01, 0.02, etc] or None"),
    ("Surface", surface_var, "Define adsorbate(s) and surface geometry", "[(O, hexagonal), ...] or None"),
    ("Freeze Threshold", freeze_var, "Define z-coordinate freeze threshold (Å)", "0.0 or float value"),
    ("Metal Library", metal_library_var, "Available elements for swapping/selection", "[Cr, Fe, Mo, etc] or None"),
    ("Execution Command", supcomp_command_var, "Command to execute LAMMPS", "lmp_serial"),
]

# Layout input fields in the left column
for i, (label, var, _, _) in enumerate(fields, start=1):
    tk.Label(root, text=label + ":").grid(row=i, column=0, sticky="w", padx=10)
    
    if isinstance(var, tk.StringVar):
        tk.Entry(root, textvariable=var, width=15).grid(row=i, column=1, sticky="w")

# Generate File Button
tk.Button(root, text="Generate Input File", command=generate_input_file).grid(row=15, column=0, columnspan=2, pady=10)

md_fields = [
    ("MD Steps", md_steps_var, "Number of MD steps per MD Simulation", "1000"),
    ("MD Interval", md_interval_var, "Interval, in MC steps, between MD simulations", "100")]

# **Static Table on the Right**
table_frame = tk.Frame(root)
table_frame.grid(row=1, column=2, rowspan=12, padx=20, pady=10, sticky="n")
tk.Label(root, text="See TUTORIAL.md for complete description").grid(row=0, column=2, columnspan=2, pady=5)

columns = ("Parameter", "Description", "Example")
table = ttk.Treeview(table_frame, columns=columns, show="headings", height=len(preamble_fields+fields+md_fields))

# Define column headings
table.heading("Parameter", text="Parameter", anchor="w")
table.heading("Description", text="Description", anchor="w")
table.heading("Example", text="Example", anchor="w")

# Define column width
table.column("Parameter", width=150, anchor="w")
table.column("Description", width=350, anchor="w")
table.column("Example", width=250, anchor="w")

for label, _, description, example in preamble_fields:
    table.insert("", "end", values=(label, description, example))  # Add row
    table.insert("", "end", values=("", "", ""))  # Add

table.insert("", "end", values=("~"*150, "~"*350, "~"*250))  # Add
table.insert("", "end", values=("", "", ""))  # Add
# Insert **static** descriptions with spacing
for label, _, description, example in fields:
    table.insert("", "end", values=(label, description, example))  # Add row
    table.insert("", "end", values=("", "", ""))  # Add empty row for spacing

table.insert("", "end", values=("~"*150, "~"*350, "~"*250))  # Add
table.insert("", "end", values=("", "", ""))  # Add

for label, _, description, example in md_fields:
    table.insert("", "end", values=(label, description, example))  # Add row
    table.insert("", "end", values=("", "", ""))  # Add

# Add horizontal scrollbar
scroll_x = ttk.Scrollbar(table_frame, orient="horizontal", command=table.xview)
table.configure(xscrollcommand=scroll_x.set)
scroll_x.pack(side="bottom", fill="x")

# Pack the Treeview widget
table.pack(fill="both", expand=True)

root.mainloop()
