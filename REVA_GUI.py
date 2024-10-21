from image_processor import *
from tkinter import filedialog, messagebox
import customtkinter as ctk
import zarr
import os
import numpy as np
from PIL import Image, ImageTk
import threading
import tkinter as tk
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
matplotlib.use('TkAgg')
from tkinter import ttk



class ImageViewer:
	last_crop_coords = None
	def __init__(self):
		"""
		Initialize the ImageViewer and display the image.
		"""
		self.original_image = None  # Convert the NumPy array to a PIL image
		self.processed_image = None  # Image after applying adjustments

		# Variables to track zoom level
		self.zoom_level = 0.1  # Start at 50%
		self.max_zoom = 4.0  # Maximum zoom level (400%)
		self.min_zoom = 0.1  # Minimum zoom level (10%)

		# Initialize 'hist_canvas' to None
		self.hist_canvas = None

		# Create a new Toplevel window for the image viewer
		self.top = tk.Toplevel()
		self.top.title("Image Viewer")
		self.top.geometry("1000x600")  # Reduced window size

		# Create a main frame to hold the control panel and canvas
		self.main_frame = tk.Frame(self.top)
		self.main_frame.pack(fill=tk.BOTH, expand=True)

		# Create a canvas for the controls to enable scrolling
		self.controls_canvas = tk.Canvas(self.main_frame, width=350)
		self.controls_canvas.pack(side=tk.LEFT, fill=tk.Y)

		# Add scrollbar to the controls canvas
		self.controls_scrollbar = tk.Scrollbar(self.main_frame, orient=tk.VERTICAL, command=self.controls_canvas.yview)
		self.controls_scrollbar.pack(side=tk.LEFT, fill=tk.Y)

		self.controls_canvas.configure(yscrollcommand=self.controls_scrollbar.set)

		# Create a frame inside the canvas to hold the controls
		self.controls_frame = tk.LabelFrame(self.controls_canvas)

		# Create a window in the canvas to hold the controls frame
		self.controls_canvas.create_window((0, 0), window=self.controls_frame, anchor='nw')

		# Bind the configure event to update the scroll region
		self.controls_frame.bind("<Configure>", lambda event: self.controls_canvas.configure(
			scrollregion=self.controls_canvas.bbox("all")))

		# Create a frame for the canvas (right side)
		self.image_frame = tk.LabelFrame(self.main_frame)
		self.image_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

		# Set up the canvas to display the image
		self.canvas = tk.Canvas(self.image_frame, background='black')
		self.canvas.pack(fill=tk.BOTH, expand=True)

		# Zoom in and zoom out buttons
		self.zoom_in_button = tk.Button(self.image_frame, text="+", command=self.zoom_in)
		self.zoom_in_button.place(x=10, y=10)  # Top-left corner with some padding
		self.zoom_out_button = tk.Button(self.image_frame, text="-", command=self.zoom_out)
		self.zoom_out_button.place(x=10, y=50)  # Next to zoom in button

		# Initialize variables for cropping
		self.rect_id = None
		self.coord_text_id = None
		self.crop_coords = None

		# Initialize variables for controls
		self.init_control_variables()

		# Add image adjustment controls
		self.add_controls()

		# Bind mouse wheel and button events for zooming and panning
		self.canvas.bind("<MouseWheel>", self.zoom)
		self.canvas.bind("<Button-4>", self.zoom)  # For Linux systems
		self.canvas.bind("<Button-5>", self.zoom)
		self.canvas.bind("<ButtonPress-1>", self.pan_start)
		self.canvas.bind("<B1-Motion>", self.pan_move)


	def init_control_variables(self):
		"""Initialize variables for image adjustment controls."""
		# Windowing (-w)
		if not hasattr(self, 'window_min'):
			self.window_min = tk.DoubleVar(value=0)
		else:
			self.window_min.set(0)
		if not hasattr(self, 'window_max'):
			self.window_max = tk.DoubleVar(value=255)
		else:
			self.window_max.set(255)
		# BreakPoint (-bk)
		if not hasattr(self, 'bk_active'):
			self.bk_active = tk.BooleanVar(value=False)
		else:
			self.bk_active.set(False)
		if not hasattr(self, 'bk_index'):
			self.bk_index = tk.IntVar(value=0)
		else:
			self.bk_index.set(0)
		# Flythrough (-f)
		if not hasattr(self, 'f_active'):
			self.f_active = tk.BooleanVar(value=False)
		else:
			self.f_active.set(False)
		# PNG Stack (-png)
		if not hasattr(self, 'png_active'):
			self.png_active = tk.BooleanVar(value=False)
		else:
			self.png_active.set(False)
		# Downgrade Sample (-8b)
		if not hasattr(self, 'b8_active'):
			self.b8_active = tk.BooleanVar(value=False)
		else:
			self.b8_active.set(False)
		# Cropping
		if not hasattr(self, 'cropping_active'):
			self.cropping_active = tk.BooleanVar(value=False)
		else:
			self.cropping_active.set(False)

	def add_controls(self):
		"""Add image adjustment controls to the controls frame."""
		# Controls Frame Title
		controls_title = tk.Label(self.controls_frame, text="Adjustments", font=("Helvetica", 12, "bold"))
		controls_title.pack(pady=(5, 5))

		self.add_file_selection_controls()

		# Controls that require user interaction (without activation checkboxes)
		self.add_interactive_controls()

		# Controls that need checkboxes (grouped together)
		self.add_checkbox_controls()

		# Reset Button
		reset_button = tk.Button(self.controls_frame, text="Reset", command=self.reset_adjustments)
		reset_button.pack(pady=5)

	def add_file_selection_controls(self):
		"""Add controls for selecting the folder, zarr file, and slices."""
		selection_frame = tk.LabelFrame(self.controls_frame, text='Select Folder', padx=5, pady=5)
		selection_frame.pack(fill='x', padx=5, pady=5)

		# Folder Selection
		folder_frame = tk.Frame(selection_frame)
		folder_frame.pack(padx=0, pady=2)

		folder_label = tk.Label(folder_frame, text='Folder')
		folder_label.pack(side=tk.LEFT)

		# Entry to show folder path, fill X to expand to the remaining space
		self.folder_entry = tk.Entry(folder_frame, width=20)
		self.folder_entry.pack(side=tk.LEFT, fill='x', expand=True, padx=5)

		browse_btn = tk.Button(folder_frame, text='...', command=self.select_folder)
		browse_btn.pack(side=tk.LEFT, padx=5)

		# Zarr file selection
		zarr_frame = tk.Frame(selection_frame)
		zarr_frame.pack(fill='x', padx=(0,10), pady=2)

		zarr_label = tk.Label(zarr_frame, text='Zarr File:')
		zarr_label.pack(side=tk.LEFT)

		self.zarr_var = tk.StringVar()
		self.zarr_option_menu = ttk.Combobox(zarr_frame, textvariable=self.zarr_var, state="readonly", width=20)

		# Bind the mouse wheel to the custom handler
		self.zarr_option_menu.bind("<Enter>", self.bind_mousewheel)
		self.zarr_option_menu.bind("<Leave>", self.unbind_mousewheel)
		self.zarr_option_menu.pack(side=tk.LEFT, padx=5)

		# Slice selection (-bk)
		slice_frame = tk.Frame(selection_frame)
		slice_frame.pack(fill='x', padx=0, pady=2)
		slice_label = tk.Label(slice_frame, text='Slices (-bk):')
		slice_label.pack(side=tk.LEFT)

		self.slice_entry = tk.Entry(slice_frame, width=10)
		self.slice_entry.pack(side=tk.LEFT, padx=5)

		# Load button
		load_button = tk.Button(selection_frame, text='Load', command=self.load_selected_zarr)
		load_button.pack(pady=5)

	def select_folder(self):
		"""Open a dialog to select the folder containing zarr files."""
		folder_selected = filedialog.askdirectory()
		if folder_selected:
			self.handle_post_selection(folder_selected)
			self.folder_entry.delete(0, tk.END)
			self.folder_entry.insert(0, folder_selected)

	def handle_post_selection(self, folder_path):
		"""Handle actions after a folder is selected."""
		self.folder_path = folder_path
		zarr_folders = [f for f in os.listdir(folder_path) if
						f.endswith('.zarr') and os.path.isdir(os.path.join(folder_path, f))]
		if zarr_folders:
			# If there are .zarr folders, display them in the option menu
			self.update_zarr_option_menu(zarr_folders)
		else:
			# No .zarr folders found, display a message
			messagebox.showerror('Error', 'No .zarr files found in the selected folder. Please select another folder.')

	def update_zarr_option_menu(self, zarr_folders):
		"""Update the zarr option menu with the list of zarr files."""
		self.zarr_files = zarr_folders
		self.zarr_option_menu['values'] = self.zarr_files
		if self.zarr_files:
			self.zarr_option_menu.current(0)
		else:
			self.zarr_var.set('')

	def disable_mouse_wheel(self, event):
		"""Prevent mouse wheel events from propagating beyond the Combobox."""
		return "break"

	def bind_mousewheel(self, event):
		"""Bind the mouse wheel to the Combobox when the mouse enters it."""
		self.zarr_option_menu.bind_all("<MouseWheel>", self.disable_mouse_wheel)

	def unbind_mousewheel(self, event):
		"""Unbind the mouse wheel from the Combobox when the mouse leaves it."""
		self.zarr_option_menu.unbind_all("<MouseWheel>")


	def load_selected_zarr(self):
		"""Load the selected zarr file and slices."""
		zarr_file = self.zarr_var.get()
		if not zarr_file:
			messagebox.showwarning("No Zarr File Selected", "Please select a Zarr file.")
			return
		# Load the zarr file
		zarr_path = os.path.join(self.folder_path, zarr_file)
		try:
			zarr_group = zarr.open(zarr_path, mode='r')
			stitched_array = zarr_group['/muse/stitched']
		except Exception as e:
			messagebox.showerror("Error Loading Zarr File", f"An error occurred while loading the Zarr file:\n{e}")
			return
		# Get slices from slice_entry
		slice_str = self.slice_entry.get()
		if slice_str:
			# Assume slices are specified as indices, e.g., "0", "1", "2"
			try:
				slices = [int(s.strip()) for s in slice_str.split(',')]
			except ValueError:
				messagebox.showerror("Invalid Slice Indices", "Please enter valid slice indices separated by commas.")
				return
			# Load the specified slices
		else:
			slices = [20]
		try:
			zarr_slice = np.array(stitched_array[slices, :, :])
			# Normalize the image data
			if zarr_slice.max() == zarr_slice.min():
				image_data = np.full(zarr_slice.shape, 128, dtype=np.uint8)
			else:
				image_data = ((zarr_slice - zarr_slice.min()) / (zarr_slice.max() - zarr_slice.min()) * 255).astype(
					np.uint8)
			# Update the image
			self.load_image_data(image_data)
		except Exception as e:
			messagebox.showerror("Error Loading Slices", f"An error occurred while loading the slices:\n{e}")
			return

	def load_image_data(self, image_data):
		"""Load the image data into the viewer."""
		# Check if image_data is 3D (multiple slices) or 2D (single slice)
		if image_data.ndim == 3:
			# Multiple slices; display the first slice by default
			self.image_data = image_data
			self.current_slice_index = 0
			self.total_slices = image_data.shape[0]
			image_to_display = image_data[self.current_slice_index]
		else:
			# Single slice
			self.image_data = image_data
			self.current_slice_index = 0
			self.total_slices = 1
			image_to_display = image_data

		# Convert the NumPy array to a PIL image
		self.original_image = Image.fromarray(image_to_display)
		self.processed_image = self.original_image.copy()
		# Initialize or update the histogram
		self.initialize_histogram()
		self.adjust_image()
		# Update slice navigation controls
		self.update_slice_navigation()

		# Slice navigation methods
	def update_slice_navigation(self):
		"""Update or create slice navigation controls if there are multiple slices."""
		if self.total_slices > 1:
			# Create navigation frame if it doesn't exist
			if not hasattr(self, 'nav_frame'):
				self.nav_frame = tk.Frame(self.controls_frame)
				self.nav_frame.pack(fill='x', padx=5, pady=5)

				prev_button = tk.Button(self.nav_frame, text='Previous Slice', command=self.show_previous_slice)
				prev_button.pack(side='left', padx=5)

				next_button = tk.Button(self.nav_frame, text='Next Slice', command=self.show_next_slice)
				next_button.pack(side='left', padx=5)

				self.slice_label = tk.Label(self.nav_frame,
												text=f"Slice {self.current_slice_index + 1} of {self.total_slices}")
				self.slice_label.pack(side='left', padx=5)
			else:
				# Update the slice label
				self.slice_label.config(text=f"Slice {self.current_slice_index + 1} of {self.total_slices}")
		else:
			# Hide navigation controls if only one slice
			if hasattr(self, 'nav_frame'):
				self.nav_frame.pack_forget()

	def show_previous_slice(self):
		"""Display the previous slice."""
		if self.current_slice_index > 0:
			self.current_slice_index -= 1
			self.update_displayed_slice()

	def show_next_slice(self):
		"""Display the next slice."""
		if self.current_slice_index < self.total_slices - 1:
			self.current_slice_index += 1
			self.update_displayed_slice()

	def update_displayed_slice(self):
		"""Update the displayed image to the current slice."""
		image_to_display = self.image_data[self.current_slice_index]

		# Convert the NumPy array to a PIL image
		self.original_image = Image.fromarray(image_to_display)
		self.processed_image = self.original_image.copy()
		self.adjust_image()
		self.initialize_histogram()
		self.update_slice_navigation()

	def add_interactive_controls(self):
		"""Add controls that require user interaction and remove activation checkboxes."""
		interactive_frame = tk.LabelFrame(self.controls_frame, text="Interactive Controls", padx=5, pady=5)
		interactive_frame.pack(fill="x", padx=5, pady=5)

		# Windowing (-w) without activation checkbox
		self.add_windowing_controls(interactive_frame)

	def add_windowing_controls(self, parent_frame):
		"""Add windowing controls to the parent frame."""
		w_frame = tk.Frame(parent_frame)
		w_frame.pack(fill="x", padx=0, pady=0)

		# Histogram with draggable markers
		self.add_histogram_with_draggable_markers(w_frame)

	def add_histogram_with_draggable_markers(self, parent_frame):
		"""Add a histogram with draggable markers for windowing control."""
		# Create a figure for the histogram
		self.hist_fig = Figure(figsize=(2.5, 1.5), dpi=100)
		self.hist_ax = self.hist_fig.add_subplot(111)

		self.hist_ax.set_facecolor('black')
		self.hist_fig.set_facecolor('black')

		# Remove axis numbers and ticks
		self.hist_ax.set_xticks([])
		self.hist_ax.set_yticks([])
		self.hist_ax.set_xlabel('')
		self.hist_ax.set_ylabel('')

		# Create a canvas to display the histogram in Tkinter
		self.hist_canvas = FigureCanvasTkAgg(self.hist_fig, master=parent_frame)
		self.hist_canvas.get_tk_widget().pack(fill='x')

		# Initialize the draggable markers
		self.min_marker = None
		self.max_marker = None
		self.min_line = None
		self.max_line = None
		self.dragging_marker = None

		# Plot the histogram and initialize the markers
		self.initialize_histogram()

		# Connect the event handlers
		self.cid_press = self.hist_fig.canvas.mpl_connect('button_press_event', self.on_press)
		self.cid_release = self.hist_fig.canvas.mpl_connect('button_release_event', self.on_release)
		self.cid_motion = self.hist_fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

	def initialize_histogram(self):
		"""Initial plot of the histogram with markers and lines."""
		# Clear the histogram
		self.hist_ax.clear()

		# Remove axis numbers and ticks
		self.hist_ax.set_xticks([])
		self.hist_ax.set_yticks([])
		self.hist_ax.set_xlabel('')
		self.hist_ax.set_ylabel('')

		if not self.original_image:
			return


		# Get the pixel intensities from the original image
		img_np = np.array(self.original_image.convert('L'))  # Convert to grayscale
		flat_img = img_np.flatten()


		# Plot the histogram
		counts, bins, patches = self.hist_ax.hist(flat_img, bins=256, range=(0, 255), color='gray')

		# Initialize markers and lines
		min_intensity = self.window_min.get()
		max_intensity = self.window_max.get()
		max_y = max(counts)  # Get the maximum y-value for placing the markers

		# Draw lines
		self.min_line = self.hist_ax.axvline(min_intensity, color='red', linewidth=3)
		self.max_line = self.hist_ax.axvline(max_intensity, color='green', linewidth=3)

		# Draw markers (triangles pointing down)
		self.min_marker, = self.hist_ax.plot(min_intensity, max_y, marker='o', color='red', markersize=6, picker=8)
		self.max_marker, = self.hist_ax.plot(max_intensity, max_y, marker='o', color='green', markersize=6, picker=8)

		# Set limits
		self.hist_ax.set_xlim(0, 255)

		# Redraw the canvas
		self.hist_canvas.draw()

	def on_press(self, event):
		"""Handle mouse press events on the histogram."""
		if event.inaxes != self.hist_ax:
			return
		contains_min, _ = self.min_marker.contains(event)
		contains_max, _ = self.max_marker.contains(event)
		if contains_min:
			self.dragging_marker = 'min'
			self.press_event = event
		elif contains_max:
			self.dragging_marker = 'max'
			self.press_event = event
		else:
			self.dragging_marker = None

	def on_motion(self, event):
		"""Handle mouse motion events on the histogram."""
		if self.dragging_marker is None:
			return
		if event.inaxes != self.hist_ax:
			return
		x = event.xdata
		if x is None:
			return
		x = max(0, min(255, x))  # Clamp x between 0 and 255
		if self.dragging_marker == 'min':
			if x >= self.window_max.get():
				x = self.window_max.get() - 1
			self.window_min.set(x)
			self.update_markers()
		elif self.dragging_marker == 'max':
			if x <= self.window_min.get():
				x = self.window_min.get() + 1
			self.window_max.set(x)
			self.update_markers()
		self.hist_canvas.draw()
		self.adjust_image()

	def on_release(self, event):
		"""Handle mouse release events on the histogram."""
		self.dragging_marker = None
		self.press_event = None

	def update_markers(self):
		"""Update the positions of the draggable markers."""
		min_intensity = self.window_min.get()
		max_intensity = self.window_max.get()

		# Get current y-axis limits
		y_min, y_max = self.hist_ax.get_ylim()

		# Calculate marker_y at 70% of the y-axis
		marker_fraction = 0.7
		marker_y = y_min + (y_max - y_min) * marker_fraction

		# Update markers
		self.min_marker.set_data([min_intensity], [marker_y])
		self.max_marker.set_data([max_intensity], [marker_y])

		# Update lines
		self.min_line.set_xdata([min_intensity, min_intensity])
		self.max_line.set_xdata([max_intensity, max_intensity])

	def update_histogram(self):
		"""Update the markers and lines on the histogram."""
		# Update markers and lines
		self.update_markers()

		# Redraw the canvas
		self.hist_canvas.draw()

	def add_checkbox_controls(self):
		"""Add controls that need checkboxes and keep textboxes."""
		checkbox_frame = tk.LabelFrame(self.controls_frame, text="Options", padx=5, pady=5)
		checkbox_frame.pack(fill="x", padx=5, pady=5)

		# Flythrough (-f)
		self.add_flythrough_controls(checkbox_frame)

		# PNG Stack (-png)
		self.add_png_controls(checkbox_frame)

		# Downgrade Sample (-8b)
		self.add_downgrade_controls(checkbox_frame)


		self.add_cropping_controls(checkbox_frame)


	def add_cropping_controls(self, parent_frame):
		crop_frame = tk.Frame(parent_frame)
		crop_frame.pack(fill="x", padx=0, pady=2)
		crop_check = tk.Checkbutton(crop_frame, text="Enable Cropping", variable=self.cropping_active,
								command=self.toggle_cropping_mode)
		crop_check.pack(anchor="w")

	def add_breakpoint_controls(self, parent_frame):
		bk_frame = tk.Frame(parent_frame)
		bk_frame.pack(fill="x", padx=0, pady=2)
		bk_label = tk.Label(bk_frame, text="Break Index:")
		bk_label.pack(side=tk.LEFT, padx=(5, 0))
		bk_entry = tk.Entry(bk_frame, textvariable=self.bk_index, width=5)
		bk_entry.pack(side=tk.LEFT)

	def add_flythrough_controls(self, parent_frame):
		f_frame = tk.Frame(parent_frame)
		f_frame.pack(fill="x", padx=0, pady=2)
		f_check = tk.Checkbutton(f_frame, text="Flythrough (-f)", variable=self.f_active)
		f_check.pack(anchor="w")

	def add_png_controls(self, parent_frame):
		png_frame = tk.Frame(parent_frame)
		png_frame.pack(fill="x", padx=0, pady=2)
		png_check = tk.Checkbutton(png_frame, text="PNG Stack (-png)", variable=self.png_active)
		png_check.pack(anchor="w")

	def add_downgrade_controls(self, parent_frame):
		b8_frame = tk.Frame(parent_frame)
		b8_frame.pack(fill="x", padx=0, pady=2)
		b8_check = tk.Checkbutton(b8_frame, text="Downgrade Sample (-8b)", variable=self.b8_active,
								  command=self.adjust_image)
		b8_check.pack(anchor="w")

	# Image adjustment methods
	def adjust_image(self, event=None):
		"""Adjust the image based on the current settings."""
		# Start with the original image
		img = self.original_image.copy()

		# Apply crop if any
		if self.crop_coords is not None:
			img = img.crop(self.crop_coords)

		# Convert PIL image to numpy array
		img_np = np.array(img)

		# Apply windowing
		min_intensity = self.window_min.get()
		max_intensity = self.window_max.get()
		# Ensure min < max
		if min_intensity >= max_intensity:
			min_intensity = max_intensity - 1
			self.window_min.set(min_intensity)
		# Apply windowing
		img_np = np.clip(img_np, min_intensity, max_intensity)
		# Normalize to 0-255
		if max_intensity - min_intensity > 0:
			img_np = ((img_np - min_intensity) * (255 / (max_intensity - min_intensity))).astype(np.uint8)
		else:
			img_np = np.zeros_like(img_np)

		# If Downgrade Sample (-8b) is active
		if self.b8_active.get():
			# Assuming the original image is of higher bit depth, scale it to 8-bit
			img_np = ((img_np / img_np.max()) * 255).astype(np.uint8)

		# Convert back to PIL image
		img = Image.fromarray(img_np)

		# Update the processed image
		self.processed_image = img

		# Update the display
		self.update_image()

		# Update the markers and lines on the histogram
		if self.hist_canvas is not None:
			self.update_histogram()


	def toggle_cropping_mode(self):
		"""Enable or disable cropping mode based on the checkbox state."""
		if self.cropping_active.get():
			# Disable panning
			self.canvas.unbind("<ButtonPress-1>")
			self.canvas.unbind("<B1-Motion>")
			# Bind cropping events to left mouse button
			self.canvas.bind("<ButtonPress-1>", self.on_button_press)
			self.canvas.bind("<B1-Motion>", self.on_move_press)
			self.canvas.bind("<ButtonRelease-1>", self.on_button_release)
		else:
			# Unbind cropping events
			self.canvas.unbind("<ButtonPress-1>")
			self.canvas.unbind("<B1-Motion>")
			self.canvas.unbind("<ButtonRelease-1>")
			# Re-bind panning
			self.canvas.bind("<ButtonPress-1>", self.pan_start)
			self.canvas.bind("<B1-Motion>", self.pan_move)
			# Reset any existing rectangle and coordinates
			if self.rect_id:
				self.canvas.delete(self.rect_id)
				self.rect_id = None
			if self.coord_text_id:
				self.canvas.delete(self.coord_text_id)
				self.coord_text_id = None


	def on_button_press(self, event):
		"""Start drawing the rectangle."""
		if not self.cropping_active.get():
			return
		# Translate mouse coordinates to canvas coordinates
		self.start_x = self.canvas.canvasx(event.x)
		self.start_y = self.canvas.canvasy(event.y)
		# Create rectangle if not yet existing
		if not self.rect_id:
			self.rect_id = self.canvas.create_rectangle(self.start_x, self.start_y, self.start_x, self.start_y,
													outline='red')
		# Create text label for coordinates
		if not self.coord_text_id:
			self.coord_text_id = self.canvas.create_text(self.start_x, self.start_y - 10, anchor="sw", text="")


	def on_move_press(self, event):
		"""Update the rectangle as the mouse moves."""
		if not self.cropping_active.get():
			return
		# Translate mouse coordinates to canvas coordinates
		curX = self.canvas.canvasx(event.x)
		curY = self.canvas.canvasy(event.y)
		# Update rectangle coordinates
		self.canvas.coords(self.rect_id, self.start_x, self.start_y, curX, curY)
		# Update the coordinate text
		x0 = min(self.start_x, curX)
		y0 = min(self.start_y, curY)
		x1 = max(self.start_x, curX)
		y1 = max(self.start_y, curY)
		# Transform to image coordinates considering zoom level
		x0_img = int(x0 / self.zoom_level)
		y0_img = int(y0 / self.zoom_level)
		x1_img = int(x1 / self.zoom_level)
		y1_img = int(y1 / self.zoom_level)
		# Update the text
		coord_text = f"({x0_img}, {y0_img}) - ({x1_img}, {y1_img})"
		self.canvas.itemconfig(self.coord_text_id, text=coord_text)
		# Position the text near the rectangle
		self.canvas.coords(self.coord_text_id, x0, y0 - 10)


	def on_button_release(self, event):
		"""Finalize the rectangle and perform the crop."""
		if not self.cropping_active.get():
			return
		# Translate mouse coordinates to canvas coordinates
		self.end_x = self.canvas.canvasx(event.x)
		self.end_y = self.canvas.canvasy(event.y)
		# Get rectangle coordinates
		x0 = min(self.start_x, self.end_x)
		y0 = min(self.start_y, self.end_y)
		x1 = max(self.start_x, self.end_x)
		y1 = max(self.start_y, self.end_y)
		# Transform to image coordinates considering zoom level
		x0_img = int(x0 / self.zoom_level)
		y0_img = int(y0 / self.zoom_level)
		x1_img = int(x1 / self.zoom_level)
		y1_img = int(y1 / self.zoom_level)
		# Store crop coordinates
		self.crop_coords = (x0_img, y0_img, x1_img, y1_img)
		# Update the class variable with the crop coordinates
		ImageViewer.last_crop_coords = self.crop_coords
		# Adjust the image to apply the crop
		self.adjust_image()
		# Remove the rectangle and text from the canvas
		self.canvas.delete(self.rect_id)
		self.rect_id = None
		self.canvas.delete(self.coord_text_id)
		self.coord_text_id = None
		# Deactivate cropping mode
		self.cropping_active.set(False)
		self.toggle_cropping_mode()

	def reset_adjustments(self):
		"""Reset image adjustments to default values."""
		# Reset variables without re-initializing them
		self.window_min.set(0)
		self.window_max.set(255)
		self.bk_active.set(False)
		self.bk_index.set(0)
		self.f_active.set(False)
		self.png_active.set(False)
		self.b8_active.set(False)
		self.cropping_active.set(False)

		# Reset crop coordinates
		self.crop_coords = None

		# Reset zoom level
		self.zoom_level = 0.1

		# Reset processed image
		self.processed_image = self.original_image.copy()
		self.update_image()

		# Update histogram
		if self.hist_canvas is not None:
			self.initialize_histogram()
			self.update_histogram()

		# Ensure event bindings are correct
		self.toggle_cropping_mode()

	def zoom_in(self):
		"""Zoom in the image."""
		if self.zoom_level < self.max_zoom:
			self.zoom_level *= 1.1
			self.update_image()

	def zoom_out(self):
		"""Zoom out the image."""
		if self.zoom_level > self.min_zoom:
			self.zoom_level /= 1.1
			self.update_image()

	def zoom(self, event):
		"""Handle zooming with mouse wheel."""
		if event.delta > 0 or event.num == 4:
			self.zoom_in()
		elif event.delta < 0 or event.num == 5:
			self.zoom_out()

	def pan_start(self, event):
		"""Begin panning the image."""
		self.canvas.scan_mark(event.x, event.y)

	def pan_move(self, event):
		"""Handle panning motion."""
		self.canvas.scan_dragto(event.x, event.y, gain=1)

	def update_image(self):
		"""Update the canvas with the current zoom level and adjustments."""
		if self.processed_image is None:
			return
		# Calculate new size based on zoom level
		new_width = int(self.processed_image.width * self.zoom_level)
		new_height = int(self.processed_image.height * self.zoom_level)

		# Resize the image using high-quality resampling
		resized_image = self.processed_image.resize((new_width, new_height), Image.LANCZOS)
		self.img_tk = ImageTk.PhotoImage(resized_image)

		# Clear the canvas and display the new image
		self.canvas.delete("all")
		self.canvas.create_image(0, 0, anchor="nw", image=self.img_tk)

		# Configure the scroll region to the size of the image (needed for panning)
		self.canvas.config(scrollregion=(0, 0, new_width, new_height))

		# Bring zoom buttons to the front
		self.zoom_in_button.lift()
		self.zoom_out_button.lift()


class PhaseHandler:
	def __init__(self, phase_name, parent_frame, mandatory_keys, optional_keys, entry_widgets):
		self.phase_name = phase_name
		self.parent_frame = parent_frame
		self.mandatory_keys = mandatory_keys
		self.optional_keys = optional_keys
		self.entry_widgets = entry_widgets
		self.checkbox_vars = {}
		# Create the phase frame
		self.create_phase_frame()

	def create_phase_frame(self):
		"""Create the phase frame and add mandatory and optional inputs"""
		self.frame = ctk.CTkFrame(self.parent_frame)

		# Add mandatory inputs label and content
		self._add_section_label("Mandatory Inputs")
		self.add_mandatory_inputs()

		# Add optional inputs label and content
		if self.optional_keys is not None:
			self._add_separator()
			self._add_section_label("Optional Inputs")
			self.add_optional_inputs()
			self.frame.pack(fill="both", expand=True)

	def _add_separator(self):
		"""Add a horizontal line separator."""
		separator = ctk.CTkFrame(self.frame, height=2, fg_color="gray")  # Create a thin gray frame
		separator.pack(fill="x", padx=10, pady=(40,20))  # Add some padding for spacing

	def _add_section_label(self, text):
		"""Helper method to add section labels"""
		label = ctk.CTkLabel(self.frame, text=text, font=("Calibri", 14, "bold"))
		label.pack(pady=10)

	def add_mandatory_inputs(self):
		"""Add mandatory inputs with entry and browse functionality for input_folder"""
		for key in self.mandatory_keys:
			input_row_frame = ctk.CTkFrame(self.frame, height=30)
			input_row_frame.pack(fill='x', pady=5)
			# Add Label
			self._add_input_label(input_row_frame, cmdInputs[key]['name'])

			# If the key is input_folder, add an entry widget with a browse button
			if key == '-if' or key == '-o':
				self._add_entry_with_browse_button(input_row_frame, f"{self.phase_name}_{key}")
				cmdInputs[key]['active'] = True
			else:
				if cmdInputs[key]['variable'] != []:
					# For other keys, add a regular entry widget
					self._add_entry_widget(input_row_frame, f"{self.phase_name}_{key}")

	def create_container_frame(self):
		"""Create the container frame to hold the scrollable frames side by side."""
		container_frame = ctk.CTkFrame(self.frame, height=50)
		container_frame.pack(fill='x')
		return container_frame

	def create_types_frame(self, container_frame):
		"""Create the scrollable frame for keys with types."""
		types_frame = ctk.CTkScrollableFrame(container_frame, width=400, height=50)
		types_frame.pack(side='left', padx=10, pady=5, fill='both', expand=True)
		return types_frame

	def create_no_types_frame(self, container_frame):
		"""Create the scrollable frame for keys without types."""
		no_types_frame = ctk.CTkScrollableFrame(container_frame, width=400, height=50)
		no_types_frame.pack(side='left', padx=10, pady=5, fill='both', expand=True)
		return no_types_frame

	def add_input_row(self, key, frame):
		"""Add a single input row to the given frame, including the checkbox and the label."""
		input_row_frame = ctk.CTkFrame(frame)
		input_row_frame.pack(fill='x', pady=5)

		# Add the checkbox for the key
		checkbox = self.add_checkbox(key, input_row_frame)

		# Add the label for the key
		self._add_input_label(input_row_frame, cmdInputs[key]['name'])

		# Add the entry widget for the output folder or other optional keys
		if key == '-o':
			self._add_entry_with_browse_button(input_row_frame, f"{self.phase_name}_{key}")
		elif cmdInputs[key]['types'] != []:
			self._add_entry_widget(input_row_frame, f"{self.phase_name}_{key}")

	def add_checkbox(self, key, frame):
		"""Add a checkbox to the frame for the given key."""
		checkbox_var = ctk.BooleanVar()
		checkbox = ctk.CTkCheckBox(
			frame,
			fg_color='#333333',
			checkbox_height=30,
			checkbox_width=30,
			corner_radius=10,
			border_width=2,
			text="",  # No text on the checkbox itself, label will be separate
			variable=checkbox_var,
			command=lambda k=key: self.toggle_optional_input(k)
		)

		# Adjust padding depending on whether the key has types or not
		if cmdInputs[key]['types'] == []:
			checkbox.pack(side="right")  # Closer padding for checkboxes without widgets
		else:
			checkbox.pack(side="right")  # Default padding for checkboxes with widgets
		self.checkbox_vars[key] = checkbox_var
		return checkbox

	def add_optional_inputs(self):
		"""Add optional inputs, separating keys with 'types' into one frame and those without 'types' into another frame."""
		# Create the container frame for holding the two frames side by side
		container_frame = self.create_container_frame()
		# Initialize frames for keys with and without types
		types_frame, no_types_frame = None, None
		# Check if any keys have types, and create the appropriate frames
		if any(cmdInputs[key]['types'] != [] for key in self.optional_keys):
			types_frame = self.create_types_frame(container_frame)
		if any(cmdInputs[key]['types'] == [] for key in self.optional_keys):
			no_types_frame = self.create_no_types_frame(container_frame)
		# Iterate over optional keys and add input rows to the appropriate frame
		for key in self.optional_keys:
			if self.requires_multiple_input(key):
				if cmdInputs[key]['types'] != [] and types_frame:
					self.add_multiple_optional_inputs(key, types_frame)
				elif cmdInputs[key]['types'] == [] and no_types_frame:
					self.add_multiple_optional_inputs(key, no_types_frame)
			else:
				if cmdInputs[key]['types'] != [] and types_frame:
					self.add_input_row(key, types_frame)
				elif cmdInputs[key]['types'] == [] and no_types_frame:
					self.add_input_row(key, no_types_frame)

	def requires_multiple_input(self, key):
		"""Determine if the key requires multiple input widgets"""
		variable = cmdInputs[key]['variable']
		if not variable:
			return False
		return len(variable) > 1

	def add_multiple_optional_inputs(self, key, frame):
		"""Add multiple widgets for optional input keys that require multiple input"""
		# Create a subframe inside the given frame to hold multiple inputs
		input_row_frame = ctk.CTkFrame(frame)
		input_row_frame.pack(fill='x')

		labels = cmdInputs[key].get('names', [])
		variable= cmdInputs[key]['variable']

		# Adjust values assignment
		if isinstance(variable, list):
			if len(variable) > 0 and isinstance(variable[0], list):
				# If variable is a list of lists, extract the first list
				values = variable[0]
			else:
				# If variable is a list of values, use it directly
				values = variable
		else:
			# If variable is a single value, convert it to a list
			values = [variable]

		self.add_checkbox(key, input_row_frame)

		self._add_input_label(input_row_frame, cmdInputs[key]['name'])

		# Now add entry widgets for each value
		for i, label_text in enumerate(labels):
			sub_frame = ctk.CTkFrame(input_row_frame)
			sub_frame.pack(side='left')

			# Create a label for each input
			point_label = ctk.CTkLabel(sub_frame, text=label_text, width=20, anchor='w', font=('Calibri', 10))
			point_label.pack(side='left', padx=(5, 5))

			# Ensure that values[i] exists
			if i < len(values):
				value = values[i]
			else:
				value = ''  # Default value if not provided

			# Create an entry widget for each value and pre-fill with current value
			entry = ctk.CTkEntry(sub_frame, width=50, height=10, corner_radius=5)
			entry.insert(0, str(value))
			entry.pack(side='left', padx = (0,5))

			# Bind each entry to update cmdInputs when the user types
			entry.bind("<KeyRelease>", self.make_update_callback(key, i))
			# Store the entry widget in self.entry_widgets
			phase_key = f"{self.phase_name}_{key}_{i}"
			self.entry_widgets[phase_key] = entry

	def make_update_callback(self, key, idx):
		"""Create a callback function that captures key and idx."""
		def callback(event):
			self.update_multiple_optional_input(key, idx, event.widget.get())
		return callback

	def update_multiple_optional_input(self, key, idx, value):
		"""Update the value in cmdInputs for keys with multiple inputs."""
		if value.strip() == '' or value.strip() == '-':
			return
		# Convert the value to the appropriate type
		try:
			if 'float' in cmdInputs[key]['types']:
				value_converted = float(value)
			elif 'int' in cmdInputs[key]['types']:
				value_converted = int(value)
			else:
				value_converted = value  # Keep as string if no type specified

			variable = cmdInputs[key]['variable']
			# Ensure the variable is structured correctly
			if isinstance(variable, list):
				if len(variable) > 0 and isinstance(variable[0], list):
					# variable is a list of lists
					if idx < len(variable[0]):
						variable[0][idx] = value_converted
					else:
						# Extend the list if necessary
						variable[0].extend([''] * (idx - len(variable[0]) + 1))
						variable[0][idx] = value_converted
				else:
					# variable is a list of values
					if idx < len(variable):
						variable[idx] = value_converted
					else:
						# Extend the list if necessary
						variable.extend([''] * (idx - len(variable) + 1))
						variable[idx] = value_converted
			else:
				# variable is a single value
				cmdInputs[key]['variable'] = [value_converted]

			print(f"Updated {cmdInputs[key]['name']}[{idx}] to {value_converted}")
		except ValueError:
			# Handle the case where the value cannot be converted
			print(f"Invalid input for {cmdInputs[key]['name']}[{idx}]")

	def _add_input_label(self, frame, text):
		"""Helper method to add input labels"""
		label = ctk.CTkLabel(frame, text=text, width=150, font=("Calibri", 12, "bold"), anchor='w')
		label.pack(side="left", padx=5)

	def _add_entry_widget(self, frame, phase_key):
		"""Helper method to add a regular entry widget"""
		entry = ctk.CTkEntry(frame, placeholder_text='Select...' ,width=200, text_color="white", justify='left', corner_radius=10, height=20)
		entry.pack(side="left")
		entry.bind("<KeyRelease>", lambda event, k=phase_key: self.update_cmd_input(k))
		self.entry_widgets[phase_key] = entry

	def _add_entry_with_browse_button(self, frame, phase_key):
		"""Helper method to add an entry widget with a small browse button inside"""
		# Create a frame to hold both the entry and the browse button
		entry_browse_frame = ctk.CTkFrame(frame)
		entry_browse_frame.pack(side="left")

		# Create the entry widget
		entry = ctk.CTkEntry(entry_browse_frame, placeholder_text='Browse...',width=500, text_color="white", corner_radius=10, height=20)
		entry.pack(side="left")
		entry.bind("<KeyRelease>", lambda event, k=phase_key: self.update_cmd_input(k))
		# Create the small browse button inside the entry frame
		browse_button = ctk.CTkButton(
			entry_browse_frame,
			text="...",
			width=10,  # Small button width
			command=lambda k=phase_key: self.browse_folder(k),
			font=("Calibri", 14, "bold"),
			fg_color="#2C2C2C",  # Button color
			text_color='white',
			corner_radius=10,
			border_width=2,  # Thickness of the border (adjustable)
		)
		browse_button.pack(side="left")

		# Store entry widget for later access
		self.entry_widgets[phase_key] = entry

	def update_cmd_input(self, phase_key):
		"""Update cmdInputs variable when the entry value is modified"""
		key=phase_key.split("_")[-1]
		entry = self.entry_widgets.get(phase_key)
		if entry:
			input_value = entry.get()
			if input_value:
				cmdInputs[key]["variable"][0] = input_value
			else:
				cmdInputs[key]['variable'][0] = ""

			print(f"Live updated {cmdInputs[key]['name']} to {cmdInputs[key]['variable'][0]}")

	def toggle_optional_input(self, key):
		"""Toggle the activation state of optional input based on the checkbox state."""
		phase_key = f"{self.phase_name}_{key}"

		checkbox_state = self.checkbox_vars.get(key)
		if checkbox_state is None:
			return

		# Get the checkbox state
		is_active = checkbox_state.get()

		# If there is an entry associated with the key, toggle the entry state
		entry = self.entry_widgets.get(phase_key)
		if entry:
			if is_active:
				self._activate_input(entry, key)
			else:
				self._deactivate_input(entry, key)
		else:
			# If there's no entry, just update the cmdInputs for this key
			cmdInputs[key]['active'] = is_active  # Update the 'active' state
			if cmdInputs[key]['types'] == []:  # If no specific input, set the variable as active/inactive
				cmdInputs[key]['variable'] = [is_active]  # Example, this could be True/False or other logic
			print(f"{cmdInputs[key]['name']} is now {'active' if is_active else 'inactive'}")  # Debug output


	def _activate_input(self, entry, key):
		"""Activate the input by enabling the entry"""
		cmdInputs[key]['active'] = True
		entry.configure(state="normal")  # Enable the entry widget
		print(f"{cmdInputs[key]['name']} is now active")  # Debug output


	def _deactivate_input(self, entry, key):
		"""Deactivate the input by disabling the entry"""
		cmdInputs[key]['active'] = False
		entry.configure(state="disabled")  # Disable the entry widget
		print(f"{cmdInputs[key]['name']} is now inactive")  # Debug output

	def browse_folder(self, phase_key):
		"""Browse for a folder and set the path to the corresponding entry widget"""
		folder_path = filedialog.askdirectory()
		if folder_path:
			if not folder_path.endswith(os.path.sep):
				folder_path += os.path.sep
			# Update the corresponding entry widget with the selected folder path
			entry = self.entry_widgets.get(phase_key)
			if entry:
				entry.delete(0, ctk.END)
				entry.insert(0, folder_path)

			# Automatically check the checkbox if the folder is browsed for optional keys
			key = phase_key.split("_")[-1]
			cmdInputs[key]['variable'][0] = folder_path  # Update cmdInputs
			print(f"Live updated {cmdInputs[key]['name']} to {cmdInputs[key]['variable'][0]}")

			if key in self.optional_keys:
				self.checkbox_vars[key].set(True)  # Activate the checkbox
				self._activate_input(entry, key)

class Stage1Handler(PhaseHandler):
	def __init__(self, parent_frame, entry_widgets):
		# Define mandatory keys:
		mandatory_keys = ['-if', '-o']
		optional_keys = ['-r', '-f', '-i']
		phase_name = 'Stage_1'
		super().__init__(phase_name, parent_frame, mandatory_keys, optional_keys, entry_widgets)


class Stage2Handler(PhaseHandler):
	def __init__(self, parent_frame, entry_widgets):
		# Define mandatory keys and call parent class
		self.scroll_frame = None
		mandatory_keys = ['-if', '-o']
		optional_keys = ['-tr', '-of', '-os', '-fg', '-f']
		phase_name = 'Stage_2'
		super().__init__(phase_name, parent_frame, mandatory_keys, optional_keys, entry_widgets)

	def browse_folder(self, phase_key):
		"""Browse for a folder and set the path to the corresponding entry widget."""
		folder_path = filedialog.askdirectory()  # Open folder dialog
		if folder_path:
			# Ensure the folder path has the appropriate separator at the end
			if not folder_path.endswith(os.path.sep):
				folder_path += os.path.sep

			# Call the method to update the entry widget with the selected folder path
			self.update_entry_with_path(phase_key, folder_path)

			# Call the method to handle post-selection logic (popup, checkboxes, etc.)
			self.handle_post_selection(phase_key, folder_path)

	def update_entry_with_path(self, phase_key, folder_path):
		"""Update the entry widget with the selected folder path and force a UI refresh."""
		entry = self.entry_widgets.get(phase_key)
		if entry:
			entry.delete(0, ctk.END)  # Clear any existing text
			entry.insert(0, folder_path)  # Insert the selected folder path
			self.parent_frame.update()  # Force UI update to ensure the path is displayed

	def handle_post_selection(self, phase_key, folder_path):
		"""Handle actions after the folder path is selected (scrape .zarr folders, display in popup)."""

		# Scrape all subfolders that end with ".zarr"
		zarr_folders = [f for f in os.listdir(folder_path) if
						f.endswith('.zarr') and os.path.isdir(os.path.join(folder_path, f))]

		# Keep the cmdInputs[key]['variable'] storing the folder path
		key = phase_key.split("_")[-1]
		cmdInputs[key]['variable'][0] = folder_path  # Keep storing the folder path in 'variable'

		# Display the .zarr folders in a popup window
		if key == '-if':
			self.show_zarr_in_scrollable_container(zarr_folders)

	def show_zarr_in_scrollable_container(self, zarr_folders):
		"""Display the .zarr folders in a scrollable frame, each with a rating dropdown next to it."""
		# Clear existing scroll frame if it exists
		if self.scroll_frame is not None:
			self.scroll_frame.pack_forget()  # Remove the scroll frame from the display
			self.scroll_frame.destroy()  # Destroy all the widgets inside it

		# Create a new scrollable frame inside the main window
		self.scroll_frame = ctk.CTkScrollableFrame(self.frame, width=400, height=300)
		self.scroll_frame.pack(padx=20, pady=20, fill="both", expand=True)

		# Initialize cmdInputs['-c']['variable'][0] as an empty list to store ratings
		cmdInputs['-c']['variable'] = [[]]  # Reset the ratings list to handle the new selection

		# Add each .zarr folder with a dropdown menu for rating
		for index, folder in enumerate(zarr_folders):
			folder_frame = ctk.CTkFrame(self.scroll_frame)
			folder_frame.pack(fill='x', padx=10, pady=5)

			# Label for the folder
			folder_label = ctk.CTkLabel(folder_frame, text=folder)
			folder_label.pack(side='left', padx=10)

			# Dropdown menu for rating selection with a pre-text
			rating_var = ctk.StringVar(value="Select Rating")  # Use StringVar with pre-text
			rating_dropdown = ctk.CTkComboBox(folder_frame, variable=rating_var, values=["1", "2", "3", "4"],
											  fg_color='#0093E9', border_color='#FBAB7E', dropdown_fg_color='#0093E9')
			rating_dropdown.pack(side='right', padx=10)

			# Bind the rating change event to update the cmdInputs directly
			rating_var.trace_add("write", lambda *args, idx=index, var=rating_var: self.update_rating(idx, var.get()))

		# Add a Close or Done button at the bottom of the scrollable frame
		# Add some padding to separate it from the folder listings
		close_button_frame = ctk.CTkFrame(self.scroll_frame)
		close_button_frame.pack(fill='both', pady=(10, 10))  # Padding at top to separate from previous content

		close_button = ctk.CTkButton(close_button_frame, text='Finish', command=self.on_done)
		close_button.pack(pady=5, padx=5, anchor='w')  # Centering the button horizontally

	def update_rating(self, index, rating):
		"""Update the rating in cmdInputs['-c']['variable'][0] when a valid rating is selected."""
		# Skip if the rating is still the pre-text ("Select Rating")
		if rating == "Select Rating" or rating == '':
			return

		# Convert the rating from string to integer before storing it
		rating_int = int(rating)

		# Ensure the list has enough elements to update by index
		ratings_list = cmdInputs['-c']['variable'][0]

		# If we're updating an existing index, replace it; otherwise, extend the list
		if len(ratings_list) > index:
			ratings_list[index] = rating_int  # Replace existing rating
		else:
			# Add 0s to fill in gaps if necessary, then append the new rating
			ratings_list.extend([0] * (index - len(ratings_list)))
			ratings_list.append(rating_int)

		# Update the list in cmdInputs['-c']['variable'][0]
		cmdInputs['-c']['variable'][0] = ratings_list
		print(f"Updated Ratings: {cmdInputs['-c']['variable'][0]}")

	def on_done(self):
		"""Handle the Done button action."""
		print(f"Final Ratings: {cmdInputs['-c']['variable'][0]}")
		self.scroll_frame.pack_forget()


class Stage3Handler(PhaseHandler):
	def __init__(self, parent_frame, entry_widgets):
		# Initialize variables
		self.crop_points = None
		self.crop_point_entries = []
		self.checkbox_vars = {}
		# Define mandatory and optional keys
		mandatory_keys = ['-if', '-o', '-fi', '-f']  # '-fi' for Crop and Histogram is mandatory
		optional_keys = None
		phase_name = 'Stage_3'
		# Initialize the parent class
		super().__init__(phase_name, parent_frame, mandatory_keys, optional_keys, entry_widgets)

	def add_mandatory_inputs(self):
		"""Add mandatory inputs with entry and browse functionality."""
		for key in self.mandatory_keys:
			input_row_frame = ctk.CTkFrame(self.frame)
			input_row_frame.pack(fill='x', pady=5)

			# Add label for the key
			self._add_input_label(input_row_frame, cmdInputs[key]['name'])

			if key == '-if' or key == '-o':
				# Add entry widget with browse button
				self._add_entry_with_browse_button(input_row_frame, f"{self.phase_name}_{key}")
				cmdInputs[key]['active'] = True

			elif key == '-fi':
				# Add crop points input
				self.add_crop_points_input(key, input_row_frame)
				cmdInputs[key]['active'] = True

			elif key == '-f':
				# Add checkbox for '-f'
				self.add_checkbox('-f', input_row_frame)
				cmdInputs[key]['active'] = True

	def add_crop_points_input(self, key, frame):
		"""Add four entry widgets for crop points."""
		# Create a subframe for crop point inputs
		input_row_frame = ctk.CTkFrame(frame)
		input_row_frame.pack(fill='x', pady=5)

		# Define labels and get current crop points
		labels = ["Height Min:", "Height Max:", "Width Min:", "Width Max:"]
		self.crop_points = cmdInputs[key]['variable'][0]  # Should be a list of four integers
		self.crop_point_entries = []  # Store entry widget references

		# Create entry widgets for each crop point
		for i, label_text in enumerate(labels):
			sub_frame = ctk.CTkFrame(input_row_frame)
			sub_frame.pack(side='left', padx=5, pady=5)

			# Create label for each crop point
			point_label = ctk.CTkLabel(sub_frame, text=label_text, width=80, anchor='w')
			point_label.pack(side='left')

			# Create entry widget with current value
			entry = ctk.CTkEntry(sub_frame, width=60)
			if i < len(self.crop_points):
				entry.insert(0, str(self.crop_points[i]))
			else:
				entry.insert(0, '0')  # Default value if not enough elements
			entry.pack(side='left')
			self.crop_point_entries.append(entry)

			# Bind entry to update cmdInputs on user input
			entry.bind("<KeyRelease>", lambda event, idx=i: self.update_crop_point(key, idx, event.widget.get()))

	def update_crop_point(self, key, idx, value):
		"""Update cmdInputs when crop point value changes."""
		value = value.strip()
		if value == '' or value == '-':
			return
		try:
			value_int = int(value)
			# Ensure the list has enough elements
			while idx >= len(cmdInputs[key]['variable'][0]):
				cmdInputs[key]['variable'][0].append(0)
			cmdInputs[key]['variable'][0][idx] = value_int
		except ValueError:
			pass  # Ignore invalid input

	def add_checkbox(self, key, frame):
		"""Add a checkbox for the given key."""
		checkbox_var = tk.BooleanVar(value=False)
		checkbox = ctk.CTkCheckBox(
			frame,
			text="",
			variable=checkbox_var,
			command=lambda k=key: self.toggle_optional_input(k)
		)
		checkbox.pack(side='left')
		self.checkbox_vars[key] = checkbox_var

	def browse_folder(self, phase_key):
		"""Browse for a folder and update the corresponding entry widget."""
		folder_path = filedialog.askdirectory()  # Open folder dialog
		if folder_path:
			if not folder_path.endswith(os.path.sep):
				folder_path += os.path.sep

			self.update_entry_with_path(phase_key, folder_path)
			self.handle_post_selection(phase_key, folder_path)

	def update_entry_with_path(self, phase_key, folder_path):
		"""Update the entry widget with the selected folder path."""
		entry = self.entry_widgets.get(phase_key)

		if entry:
			entry.delete(0, ctk.END)
			entry.insert(0, folder_path)
			self.parent_frame.update()

	def handle_post_selection(self, phase_key, folder_path):
		"""Handle actions after a folder is selected."""
		key = phase_key.split("_")[-1]
		cmdInputs[key]['variable'][0] = folder_path
		if key == '-if':
			# Open ImageViewer and wait until it closes
			image_viewer = ImageViewer()
			self.parent_frame.wait_window(image_viewer.top)

			# Get crop coordinates from ImageViewer
			crop_coords = ImageViewer.last_crop_coords
			if crop_coords is not None:
				# Map coordinates to Height and Width
				x0_img, y0_img, x1_img, y1_img = crop_coords
				crop_points = [y0_img, y1_img, x0_img, x1_img]

				# Update cmdInputs
				cmdInputs['-fi']['variable'][0] = crop_points

				# Update existing entry widgets
				for idx, value in enumerate(crop_points):
					if idx < len(self.crop_point_entries):
						entry_widget = self.crop_point_entries[idx]
						entry_widget.delete(0, ctk.END)
						entry_widget.insert(0, str(value))
					else:
						print(f"Index {idx} out of bounds for crop_point_entries")
			else:
				# No crop performed; handle if necessary
				pass  # Optional: Handle cases where no cropping is done

	def toggle_optional_input(self, key):
		"""Toggle the active state of an optional input."""
		cmdInputs[key]['active'] = self.checkbox_vars[key].get()

class Stage4Handler(PhaseHandler):
	def __init__(self, parent_frame, entry_widgets):
		# Define mandatory keys and call parent class
		mandatory_keys = ['-if', '-o']
		optional_keys = ['-e', '-d', '-ms', '-ma', '-bt', '-bk', '-br', '-f', '-png', '-w']
		phase_name = 'Stage_4'
		super().__init__(phase_name, parent_frame, mandatory_keys, optional_keys, entry_widgets)


class REVA_GUI(ctk.CTk):
	def __init__(self):
		super().__init__()
		self.title("REVA_GUI")
		self.geometry("1200x700")

		# Set the theme to dark for a consistent black background
		ctk.set_appearance_mode("dark-blue")
		self.entry_widgets = {}

		# Create a frame for the buttons on the left
		self.buttons_frame = ctk.CTkFrame(self, fg_color='#2C2C2C')
		self.buttons_frame.pack(side="left", padx=10, pady=(10, 0), fill="y", expand=False)

		# Create buttons for each phase
		self.survey_btn = ctk.CTkButton(self.buttons_frame, text="Stage 1", font=("Calibri", 14, "bold"),
										command=lambda: self.show_phase_frame("Stage_1"), width=200)
		self.survey_btn.pack(padx=5, pady=5)

		self.compile_btn = ctk.CTkButton(self.buttons_frame, text="Stage 2", font=("Calibri", 14, "bold"),
										 command=lambda: self.show_phase_frame("Stage_2"), width=200)
		self.compile_btn.pack(padx=5, pady=5)

		self.finish_btn = ctk.CTkButton(self.buttons_frame, text="Stage 3", font=("Calibri", 14, "bold"),
										command=lambda: self.show_phase_frame("Stage_3"), width=200)
		self.finish_btn.pack(padx=5, pady=5)

		# Create a frame for content on the right (where survey/compile/finish frames will appear)
		self.content_frame = ctk.CTkFrame(self)
		self.content_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

		# Add Process button
		self.run_button = ctk.CTkButton(self.content_frame, text="Process", command=self.run_process,
										font=("Calibri", 14, "bold"), width=200)
		self.run_button.pack(padx=5, pady=5, side='bottom', fill='x')

		# Add a spacer that expands to take up available space
		self.spacer = ctk.CTkFrame(self.buttons_frame)
		self.spacer.pack(padx=5, pady=5, fill='both', expand=True)

		# Dictionary to hold PhaseHandler objects for different phases
		self.frames = {}

	def show_phase_frame(self, phase_name):
		"""Show the relevant phase frame based on the button clicked"""
		self.hide_frames()
		if phase_name not in self.frames:
			self.frames[phase_name] = self.create_phase_handler(phase_name)
		self.frames[phase_name].frame.pack(fill="both", expand=True)
		self.activate_stage(phase_name)

	def activate_stage(self, phase_name):
		# Deactivate all stages first
		self.deactivate_all_stages()
		# Activate only the selected stage
		if phase_name == "Stage_2":
			cmdInputs['-c']['active'] = True  # Activate Stage 2
			print("Stage 2 activated")
		elif phase_name == 'Stage_1':
			cmdInputs['-s']['active'] = True  # Activate Stage 1
			print("Stage 1 activated")
		elif phase_name == 'Stage_3':
			cmdInputs['-fi']['active'] = True  # Activate Stage 3
			print("Stage 3 activated")


	def deactivate_all_stages(self):
		"""Deactivate all stages by setting their 'active' state to False."""
		cmdInputs['-c']['active'] = False  # Deactivate Stage 2
		cmdInputs['-s']['active'] = False  # Deactivate Stage 1
		cmdInputs['-fi']['active'] = False  # Deactivate Stage 3
		cmdInputs['-p']['active'] = False  # Deactivate Stage 4

	def hide_frames(self):
		"""Hide all phase frames"""
		for frame in self.frames.values():
			frame.frame.pack_forget()

	def create_phase_handler(self, phase_name):
		"""Create PhaseHandler object for the selected phase"""
		if phase_name == "Stage_1":
			return Stage1Handler(self.content_frame, self.entry_widgets)
		elif phase_name == "Stage_2":
			return Stage2Handler(self.content_frame, self.entry_widgets)
		elif phase_name == "Stage_3":
			return Stage3Handler(self.content_frame, self.entry_widgets)
		elif phase_name == "Stage_4":
			return Stage4Handler(self.content_frame, self.entry_widgets)

	def run_process(self):
		"""Run the overall process when clicking the Process button"""
		# Disable the button to prevent multiple clicks
		self.run_button.configure(state="disabled")
		self.run_button.configure(text="Processing...")

		# Start the processing in a new thread
		processing_thread = threading.Thread(target=self.process_task)
		processing_thread.start()

	def process_task(self):
		"""The heavy processing task"""
		if cmdInputs['-r']['active']:
			recursive_active_path = cmdInputs['-if']['variable'][0]
			flist = findAllDir(recursive_active_path)
			for f in flist:
				# Perform the processing task
				globals(f, recursive_active_path, cmdInputs)
		else:
			try:
				selected_path = cmdInputs['-if']['variable'][0]  # Get the input path from cmdInputs
				split_path = selected_path.split("/")  # Split the path by '/'
				file_name = split_path[-2]  # Extract the second-to-last element as file name
				input_path = "/".join(
					split_path[:-2]) + "/"  # Join the elements before the second-to-last as the input path
				# Perform the processing task
				globals(file_name, input_path, cmdInputs)
			except Exception:
				# Example usage with a paraphrased message
				messagebox.showerror("Error", "Please select the specific sample you wish to process.")
		# Schedule the finish_process method to run on the main thread
		self.after(0, self.finish_process)

	def finish_process(self):
		print("Processing completed!")
		# Re-enable the button and restore text
		self.run_button.configure(state="normal", text="Process")



# Create and run the app
if __name__ == "__main__":
	app = REVA_GUI()
	app.mainloop()