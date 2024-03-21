import tkinter as tk

GRID_SIZE = 10

root = tk.Tk()
root.title("Grid with Values")

# Create a frame to hold the grid
frame = tk.Frame(root)
frame.pack()

# Generate values for the grid
value_layers = []
for i in range(GRID_SIZE):
    row_values = []
    for j in range(GRID_SIZE):
        tgf_beta = round(0.72 - (i * 0.08), 2)
        il_1_beta = il_6 = tnfa = round(j * 0.08, 2)
        row_values.append([tgf_beta, il_1_beta, il_6, tnfa])

    value_layers.append(row_values)

# Display the grid with values
for i in range(GRID_SIZE):
    for j in range(GRID_SIZE):
        cell_values = value_layers[i][j]
        cell_text = f"TGF-Beta: {cell_values[0]}\nIL-1 Beta: {cell_values[1]}\nIL-6: {cell_values[2]}\nTNF Alpha: {cell_values[3]}"

        cell_frame = tk.Frame(frame, relief=tk.RIDGE, borderwidth=2)
        cell_frame.grid(row=i, column=j, padx=5, pady=5)

        label = tk.Label(cell_frame, text=cell_text, justify=tk.LEFT)
        label.pack()

root.mainloop()
