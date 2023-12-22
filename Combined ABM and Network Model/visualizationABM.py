'''
Visual Improvements:
 - Add grid lines to the left and top
 - Correct spacing between elements in the key
 - Center middle the key and other info on the bottom
'''

import tkinter as tk
import time
import random
from NetworkModel import NetworkModel
import numpy as np


class GridSimulation:

    def __init__(self, master):
        self.master = master
        self.master.title("Grid Simulation")

        self.titles = ["Weight of TGFB: ", "Weight of TNFA: ", "Weight of IL-6: ", "Weight of IL-1B: ", "Collagen: "]
        self.results_dict = None

        self.grid_size = 10
        self.current_step = 0
        self.value_layers = []

        # Initialize grid
        self.grid = [[0] * self.grid_size for _ in range(self.grid_size)]

        # Create canvas for visualization
        self.canvas = tk.Canvas(master, width=800, height=500)
        self.canvas.pack()

        # Create grid lines with each cell having lines on all four sides
        cell_size = 500 // self.grid_size
        for i in range(self.grid_size):
            x = i * cell_size
            y = 0
            self.canvas.create_line(x, y, x, y + 500, fill="black")  # Vertical lines
            self.canvas.create_line(y, x, y + 500, x, fill="black")  # Horizontal lines

            # Create grid squares and store item IDs in the grid list
            row = []
            for j in range(self.grid_size):
                rect = self.canvas.create_rectangle(i * cell_size, j * cell_size, (i + 1) * cell_size,
                                                    (j + 1) * cell_size)
                row.append(rect)
            self.grid.append(row)

        # Create movable object (isosceles triangle pointing to the right)
        x0, y0 = 0, 0
        x1, y1 = cell_size, cell_size / 2
        x2, y2 = 0, cell_size
        self.object = self.canvas.create_polygon(x0, y0, x1, y1, x2, y2, fill="#3498db")

        # Create step label
        self.step_label = tk.Label(master, text="Step: 0")
        self.step_label.pack()

        # Create key frame
        self.key_frame = tk.Frame(master)
        self.key_frame.pack()

        # Create key label
        self.key_label = tk.Label(self.key_frame, text="Fibroblast Position: (0, 0)")
        self.key_label.pack()

        # Create key canvas
        self.key_canvas = tk.Canvas(self.key_frame, width=420, height=50)
        self.key_canvas.pack(side=tk.BOTTOM)

        # Add shapes and labels to the key
        self.add_key()

        # Create start/stop button
        self.start_button_text = tk.StringVar()
        self.start_button_text.set("Start Simulation")
        self.start_button = tk.Button(master, textvariable=self.start_button_text, command=self.toggle_simulation)
        self.start_button.pack(side=tk.LEFT)

        # Create reset button
        self.reset_button = tk.Button(master, text="Reset Simulation", command=self.reset_simulation)
        self.reset_button.pack(side=tk.LEFT)

        # Create toggle button for showing the previous position blue box
        self.show_prev_box_var = tk.BooleanVar()
        self.show_prev_box_var.set(True)  # Set the initial state to True
        self.show_prev_box_button = tk.Checkbutton(master, text="Show Previous Location",
                                                   variable=self.show_prev_box_var,
                                                   command=self.update_show_prev_box)
        self.show_prev_box_button.pack(side=tk.LEFT)

        # Simulation state
        self.simulation_running = False
        self.object_position = (0, 0)  # Initial position of the fibroblast

        # Keep track of the current yellow rectangle
        self.current_yellow_rectangle = None

        # Keep track of the current blue rectangle
        self.current_blue_rectangle = None

        self.info_canvas = tk.Canvas(master, width=300, height=300)
        self.info_canvas.place(x=501, y=1)

        # Create initial changing values
        self.changing_values = [0.25] * 5

        # Create a graph canvas
        self.graph_canvas = tk.Canvas(master, width=300, height=200)
        self.graph_canvas.place(x=501, y=300)

        # Display the changing values on the info canvas
        self.update_info_canvas()

        self.SPECIFIC_TIME = 10
        self.simulation_counter = 0

        # Initialize Network Model
        self.network_model_sim = NetworkModel(specific_time=self.SPECIFIC_TIME)

        for i in range(self.grid_size):
            row = []
            for j in range(self.grid_size):
                tgf_beta = (0.72 - (i * 0.08))
                il_1_beta = il_6 = tnfa = j * 0.08
                row.append([tgf_beta, il_1_beta, il_6, tnfa])
            self.value_layers.append(row)

        # Initialize default input parameters for all Network Model simulations
        self.initialize_default_params()

        self.abm_array = [0.25] * 5

    def initialize_default_params(self):



        # Define default parameters
        default_tau_values = np.array(
            [1, 1.000000e-01, 10, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01,
             1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1,
             1.000000e-01, 1.000000e-01, 10, 10, 1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01,
             1.000000e-01,
             1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01,
             1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01,
             1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01,
             1.000000e-01,
             1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1, 10,
             1.000000e-01, 1.000000e-01, 10, 1.000000e-01, 10, 10, 1.000000e-01, 10, 10, 10, 1, 1, 1, 1, 10, 10, 10, 10,
             10,
             10, 1, 1, 10, 10])

        default_y_max_values = np.array(
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5.000000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1,
             1, 1, 1, 1, 5.000000e-01, 5.000000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        # index 1: TGFB, index 3: IL-6, index 4: IL-1B, index 5: TNFA
        w_values = np.array(
            [2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01,
             2.500000e-01,
             2.500000e-01, 2.500000e-01, 2.500000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        default_n_values = np.array(
            [1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00,
             1.400000e+00,
             1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00])
        default_ec50_values = np.array(
            [6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01,
             6.000000e-01,
             6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01])

        # # Set default parameters for all ABM simulations
        # for i in range(self.grid_size):
        #     for j in range(self.grid_size):
        #         # index 1: TGFB, index 3: IL-6, index 4: IL-1B, index 5: TNFA
        #         w_values[1] = self.value_layers[i][j][0]
        #         w_values[3] = self.value_layers[i][j][1]
        #         w_values[4] = self.value_layers[i][j][1]
        #         w_values[5] = self.value_layers[i][j][1]
        #         self.network_model_sim[i][j].set_params(default_tau_values, default_y_max_values, w_values,
        #                                                 default_n_values, default_ec50_values)

    def toggle_simulation(self):
        if self.simulation_running:
            self.stop_simulation()
        else:
            self.start_simulation()

    def start_simulation(self):
        self.simulation_running = True
        self.start_button_text.set("Stop Simulation")
        self.animate_button()
        self.run_simulation()

    def stop_simulation(self):
        self.simulation_running = False
        self.start_button_text.set("Start Simulation")
        self.animate_button()

    def reset_simulation(self):
        self.stop_simulation()
        self.current_step = 0
        cell_size = 500 // self.grid_size
        self.object_position = (0, 0)  # Reset position of the fibroblast
        self.canvas.coords(self.object, 0, 0, cell_size, cell_size / 2, 0, cell_size)
        self.step_label.config(text="Step: 0")

        # Reset key label
        self.key_label.config(text="Fibroblast Position: (0, 0)")

        # Clear the previous yellow and blue rectangles
        if self.current_yellow_rectangle is not None:
            self.canvas.delete(self.current_yellow_rectangle)

        if self.current_blue_rectangle is not None:
            self.canvas.delete(self.current_blue_rectangle)

        # Reset the yellow and blue rectangles
        self.current_yellow_rectangle = None
        self.current_blue_rectangle = None

        # Reset changing values to 0.25
        self.changing_values = [0.25] * 5

        # Clear info canvas
        self.info_canvas.delete("changing_values")
        self.update_info_canvas()

    def animate_button(self):
        # Animate button press
        self.start_button.config(relief=tk.SUNKEN)
        self.master.after(100, self.reset_button_state)

    def reset_button_state(self):
        # Reset button state after animation
        self.start_button.config(relief=tk.RAISED)

    def run_simulation(self):
        while self.simulation_running:
            # Update simulation state here (move object, update grid, etc.)
            self.move_object()

            # # Run ABM simulations for each cell in the grid
            # for i in range(self.grid_size):
            #     for j in range(self.grid_size):
            #         y0 = np.array(
            #             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            #              0, 0, 0, 0, 0,
            #              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            #              0, 0, 0, 0, 0,
            #              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            #         self.ode_simulations[i][j].specific_time = (10 * self.simulation_counter) + 10
            #         self.results_dict = self.ode_simulations[i][j].run_simulation_and_return_dict(
            #             [(10 * self.simulation_counter), (10 * self.simulation_counter) + 10], y0)
            #
            # self.abm_array = ABMSimulation(
            #     self.ode_simulations[self.object_position[0]][self.object_position[1]]).run_abm_simulation(
            #     self.results_dict)

            # Update key label
            self.key_label.config(text=f"Fibroblast Position: {self.object_position}")

            # Increment step
            self.current_step += 1
            self.step_label.config(text=f"Step: {self.current_step}")

            # Update changing values
            self.update_changing_values()

            # Increment simulation counter
            self.simulation_counter += 1

            # Add a delay to make the simulation visually observable
            time.sleep(0.5)
            self.master.update()

    def move_object(self):
        # Example: Move object to a random direction
        prev_x, prev_y = self.object_position
        possible_moves = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # Right, Down, Left, Up

        # Keep selecting a new direction until the fibroblast moves to a different position
        while True:
            move = random.choice(possible_moves)
            new_x = max(0, min(prev_x + move[0], self.grid_size - 1))
            new_y = max(0, min(prev_y + move[1], self.grid_size - 1))

            if (new_x, new_y) != (prev_x, prev_y) and self.grid[new_x][new_y] != 1:
                break

        # Check if the new position is within the grid boundaries
        if 0 <= new_x < self.grid_size and 0 <= new_y < self.grid_size:
            self.object_position = (new_x, new_y)
            x, y = self.object_position

            cell_size = 500 // self.grid_size
            self.canvas.coords(self.object, x * cell_size, y * cell_size, (x + 1) * cell_size,
                               y * cell_size + cell_size / 2, x * cell_size, (y + 1) * cell_size)

            # Delete the current yellow rectangle if it exists
            if self.current_yellow_rectangle is not None:
                self.canvas.delete(self.current_yellow_rectangle)

            # Delete the current blue rectangle if it exists and the toggle button is checked
            if self.current_blue_rectangle is not None and self.show_prev_box_var.get():
                self.canvas.delete(self.current_blue_rectangle)

            # Create a new yellow rectangle at the current position
            self.current_yellow_rectangle = self.canvas.create_rectangle(x * cell_size, y * cell_size,
                                                                         (x + 1) * cell_size,
                                                                         (y + 1) * cell_size, fill="yellow")

            # Create a new blue rectangle at the previous position if the toggle button is checked
            if self.show_prev_box_var.get():
                self.current_blue_rectangle = self.canvas.create_rectangle(prev_x * cell_size, prev_y * cell_size,
                                                                           (prev_x + 1) * cell_size,
                                                                           (prev_y + 1) * cell_size,
                                                                           fill="blue")

            # Raise the triangle to the front
            self.canvas.tag_raise(self.object)

    def update_info_canvas(self):
        # Display changing values on the info canvas
        title_text = "Important Cytokines and Collagen"
        self.info_canvas.create_text(150, 10, text=title_text, font=("Arial", 12, "bold underline"),
                                     tags="changing_values")

        for i, value in enumerate(self.changing_values):
            text = self.titles[i] + str(self.changing_values[i])
            self.info_canvas.create_text(150, 40 + i * 30, text=text, font=("Arial", 10), tags="changing_values")

    def update_changing_values(self):
        # Update changing values randomly
        for i in range(len(self.changing_values)):
            self.changing_values[i] = self.abm_array[i]

        # Clear and update the info canvas
        self.info_canvas.delete("changing_values")
        self.update_info_canvas()

    def update_show_prev_box(self):
        # Update the display based on the state of the show_prev_box toggle button
        if not self.show_prev_box_var.get() and self.current_blue_rectangle is not None:
            self.canvas.delete(self.current_blue_rectangle)

    def add_key(self):
        cell_size = 30
        text_padding = 20
        shape_text_space = 80
        key_space = 30

        # Add blue background shape to the key
        self.key_canvas.create_rectangle(10, 10, 10 + cell_size, 10 + cell_size, fill="blue")
        # Add blue background label
        self.key_canvas.create_text(10 + cell_size + text_padding, 10 + cell_size // 2,
                                    text="Previous Location                           ",
                                    anchor=tk.W)

        # Add space between shapes and text
        offset = cell_size + 2 * text_padding + shape_text_space

        # Add yellow background shape to the key
        self.key_canvas.create_rectangle(10 + offset, 10, 10 + offset + cell_size, 10 + cell_size, fill="yellow")
        # Add yellow background label
        self.key_canvas.create_text(10 + offset + cell_size + text_padding, 10 + cell_size // 2, text="New Location",
                                    anchor=tk.W)

        # Add space between shapes and text
        offset += cell_size + 2 * text_padding + shape_text_space

        # Add fibroblast (triangle) shape to the key
        x0, y0 = 10 + offset, 10
        x1, y1 = x0 + cell_size, y0 + cell_size / 2
        x2, y2 = x0, y0 + cell_size
        self.key_canvas.create_polygon(x0, y0, x1, y1, x2, y2, fill="#3498db")
        # Add fibroblast label
        self.key_canvas.create_text(x1 + text_padding, y0 + cell_size // 2, text="Fibroblast", anchor=tk.W)


if __name__ == "__main__":
    root = tk.Tk()

    simulation = GridSimulation(root)

    root.mainloop()
