import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import math
import ctypes

try:
    # Windows 8.1 and later
    ctypes.windll.shcore.SetProcessDpiAwareness(1) 
except:
    # Windows 8.0 or earlier (optional fallback)
    try:
        ctypes.windll.user32.SetProcessDPIAware()
    except:
        pass # Not all versions of Windows have this function

# --- Styling Constants ---
FONT_MAIN = ("Space Grotesk", 9)
FONT_BOLD = ("Space Grotesk", 9, "bold")
FONT_HEADER = ("Space Grotesk", 12, "bold")
BG_COLOR = "#e2ffff"
ACCENT_COLOR = "#0055ff"

class PIDSimulatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("PIDLab10")
        self.root.geometry("1300x900")
        self.root.configure(bg=BG_COLOR)

        # Variables for UI State
        self.mode_var = tk.StringVar(value='Free PID Mode')
        self.motor_type_var = tk.StringVar(value='NEO')
        self.pid_type_var = tk.StringVar(value='Speed PID')
        
        # Data storage for inputs
        self.entries = {}
        
        self._setup_ui()
        
        # Trigger initial visibility and simulation
        self.update_load_inputs_visibility()
        self.update_motor_inputs_visibility()
        self.run_simulation()

    def _setup_ui(self):
        # --- Main Layout ---
        # --- Left Panel Container ---
        self.left_container = tk.Frame(self.root, bg="#ffffff", width=580)
        self.left_container.pack(side=tk.LEFT, fill=tk.Y)

        # Canvas (scrollable area)
        self.left_canvas = tk.Canvas(
            self.left_container,
            bg=BG_COLOR,
            highlightthickness=0,
            width=560
        )
        self.left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        style = ttk.Style()
        style.theme_use('clam')
        style.configure("Vertical.TScrollbar", troughcolor=BG_COLOR, background=BG_COLOR, foreground=ACCENT_COLOR, handlecolor=ACCENT_COLOR, bordercolor=ACCENT_COLOR, borderwidth=1, width="20", relief="solid", arrowcolor=ACCENT_COLOR)
        style.map(
            "Vertical.TScrollbar",
            background=[('disabled', BG_COLOR)], # Thumb color when disabled
            troughcolor=[('disabled', BG_COLOR)],     # Trough color when disabled (optional)
            bordercolor=[('disabled', ACCENT_COLOR)],
            width=[('disabled', 20)]
        )
        # Scrollbar (outside canvas)
        self.left_scrollbar = ttk.Scrollbar(
            self.left_container,
            orient="vertical",
            command=self.left_canvas.yview,
            style="Vertical.TScrollbar"
        )
        self.left_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.left_canvas.configure(yscrollcommand=self.left_scrollbar.set)
        self.left_frame = tk.Frame(self.left_canvas, bg=BG_COLOR, padx=10, pady=10)
        self.left_canvas.create_window(
            (0, 0),
            window=self.left_frame,
            anchor="nw"
        )
        self.left_frame.bind(
            "<Configure>",
            lambda e: self.left_canvas.configure(
                scrollregion=self.left_canvas.bbox("all")
            )
        )

        # Right Panel: Graphs
        self.right_frame = tk.Frame(self.root, bg="white")
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.appheading = tk.Text(self.left_frame, wrap="word", bg=BG_COLOR, width=14, height=4, highlightthickness=0, state=tk.NORMAL, borderwidth=0)
        self.appheading.pack(expand=True)
        self.appheading.tag_configure("Main", font=("Space Grotesk", 18, "bold"), foreground=ACCENT_COLOR)
        self.appheading.tag_configure("Georgia", font=("Georgia", 18), foreground=ACCENT_COLOR)
        self.appheading.tag_configure("GeorgiaB", font=("Georgia", 18, "bold"), foreground="#ff0000", border=1)
        self.appheading.tag_configure("MainS", font=("Space Grotesk", 8), foreground=ACCENT_COLOR, justify='center')
        self.appheading.insert(tk.END, "PID", "Main")
        self.appheading.insert(tk.END, "Lab ", "Georgia")
        self.appheading.insert(tk.END, "10", "GeorgiaB")
        self.appheading.insert(tk.END, "\nby glerpstudios", "MainS")
        self.appheading.configure(state="disabled")
        self.appheading.bindtags((str(self.appheading), str(root), "all"))
        # B1: System Mode
        self._create_header(self.left_frame, "1. System Mode")
        self.mode_dropdown = ttk.Combobox(self.left_frame, textvariable=self.mode_var, state="readonly", 
                                          values=['Free PID Mode', 'Horizontal Load-Bearing PID Mode', 'Vertical Load-Bearing PID Mode'])
        self.mode_dropdown.pack(fill=tk.X, pady=5)
        self.mode_dropdown.bind("<<ComboboxSelected>>", self.on_layout_change)

        # --- B2: Load Parameters (Conditional) ---
        self.load_frame = tk.Frame(self.left_frame, bg=BG_COLOR, bd=1, relief=tk.SOLID)
        
        tk.Label(self.load_frame, text="2. Load Parameters", font=FONT_BOLD, bg=BG_COLOR).grid(row=0, column=0, columnspan=3, sticky="w", pady=5)
        
        # B2i, B2ii
        tk.Label(self.load_frame, text="Distance (m)", font=FONT_MAIN, bg=BG_COLOR).grid(row=1, column=0)
        self.entries['load_dist'] = self._create_entry(self.load_frame, row=2, col=0)
        
        tk.Label(self.load_frame, text="Mass (kg)", font=FONT_MAIN, bg=BG_COLOR).grid(row=1, column=1)
        self.entries['load_mass'] = self._create_entry(self.load_frame, row=2, col=1)

        # B2iii
        tk.Label(self.load_frame, text="--Or--", font=("Space Grotesk", 8, "italic"), bg=BG_COLOR).grid(row=2, column=2, padx=5)

        # B2iv
        tk.Label(self.load_frame, text="Moment of Inertia (kg·m²)", font=FONT_MAIN, bg=BG_COLOR).grid(row=1, column=3)
        self.entries['load_inertia_direct'] = self._create_entry(self.load_frame, row=2, col=3)

        # B3: Motor Type
        self._create_header(self.left_frame, "3. Motor Type")
        self.motor_dropdown = ttk.Combobox(self.left_frame, textvariable=self.motor_type_var, state="readonly",
                                           values=["NEO", "Kraken X60", "Custom..."])
        self.motor_dropdown.pack(fill=tk.X, pady=5)
        self.motor_dropdown.bind("<<ComboboxSelected>>", self.on_layout_change)

        # --- B4: Custom Motor Parameters (Conditional) ---
        self.custom_motor_frame = tk.Frame(self.left_frame, bg=BG_COLOR, bd=1, relief=tk.SOLID, padx=5, pady=5)

        tk.Label(self.custom_motor_frame, text="4. Motor Constants", font=FONT_BOLD, bg=BG_COLOR).pack(anchor="w")

        # B4i: Motor Inertia
        self._create_labeled_entry(self.custom_motor_frame, "Motor inertia (kg·m²)", 'motor_inertia')
        
        # B4ii: Empirical Free Running Current
        self._create_labeled_entry(self.custom_motor_frame, "Empirical Free Running Current (A)", 'free_current')

        # B4iii: Variable Area (Torque Constants)
        self.b4iii_frame = tk.Frame(self.custom_motor_frame, bg=BG_COLOR)
        self.b4iii_frame.pack(fill=tk.X, pady=5)

        # B4iv: Continuous Current
        self._create_labeled_entry(self.custom_motor_frame, "Continuous Current (A)", 'imax')

        # B4v: Motor Resistance
        self._create_labeled_entry(
            self.custom_motor_frame,
            "Motor Resistance (Ohms)",
            "motor_resistance"
        )

        # B4vi: Viscous Damping Area
        self.b4vi_frame = tk.Frame(self.custom_motor_frame, bg=BG_COLOR)
        self.b4vi_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(self.b4vi_frame, text="Empirical Free Speed (RPM)", font=FONT_MAIN, bg=BG_COLOR).grid(row=0, column=0)
        self.entries['free_speed'] = self._create_entry(self.b4vi_frame, row=1, col=0)
        
        tk.Label(self.b4vi_frame, text="--Or--", bg=BG_COLOR).grid(row=1, column=1)
        
        tk.Label(self.b4vi_frame, text="Viscous Damping Factor (Nm·s/rad)", font=FONT_MAIN, bg=BG_COLOR).grid(row=0, column=2)
        self.entries['viscous_damping_direct'] = self._create_entry(self.b4vi_frame, row=1, col=2)

        # --- B5: PID Coefficients ---
        self._create_header(self.left_frame, "5. PID Coefficients")
        pid_frame = tk.Frame(self.left_frame, bg=BG_COLOR)
        pid_frame.pack(fill=tk.X)
        
        tk.Label(pid_frame, text="K_p", bg=BG_COLOR).pack(side=tk.LEFT, expand=True)
        self.entries['kp'] = self._create_entry(pid_frame, side=tk.LEFT, default="0.1")
        
        tk.Label(pid_frame, text="K_i", bg=BG_COLOR).pack(side=tk.LEFT, expand=True)
        self.entries['ki'] = self._create_entry(pid_frame, side=tk.LEFT, default="0.0")
        
        tk.Label(pid_frame, text="K_d", bg=BG_COLOR).pack(side=tk.LEFT, expand=True)
        self.entries['kd'] = self._create_entry(pid_frame, side=tk.LEFT, default="0.0")

        # --- B6: PID Type ---
        self._create_header(self.left_frame, "6. Controller Type")
        self.pid_type_dropdown = ttk.Combobox(self.left_frame, textvariable=self.pid_type_var, state="readonly",
                                              values=["Speed PID", "Position PID"])
        self.pid_type_dropdown.pack(fill=tk.X)
        self.pid_type_dropdown.bind("<<ComboboxSelected>>", self.on_ui_change)

        # --- B7: Targets & Error ---
        self._create_header(self.left_frame, "7. Targets & Limits")
        
        tk.Label(self.left_frame, text="Acceptable Error Bound:", bg=BG_COLOR).pack(anchor="w")
        self.entries['error_bound'] = self._create_entry(self.left_frame, default="0.1")

        tk.Label(self.left_frame, text="Targets (Value, Time) - One per line:", bg=BG_COLOR).pack(anchor="w", pady=(5,0))
        self.target_text = tk.Text(self.left_frame, height=5, width=30, font=FONT_MAIN)
        self.target_text.pack(fill=tk.X, pady=2)
        self.target_text.insert("1.0", "10, 0\n0, 5") 
        self.target_text.bind("<KeyRelease>", self.on_ui_change)

        # --- B8: Gearbox (NEW) ---
        self._create_header(self.left_frame, "8. Gearbox")
        tk.Label(self.left_frame, text="Gear Ratio (X:1)", bg=BG_COLOR).pack(anchor="w")
        self.entries['gear_ratio'] = self._create_entry(self.left_frame, default="1.0")

        # --- B9: Feedforward ---
        self._create_header(self.left_frame, "9. Feedforward (Optional)")

        ff_frame = tk.Frame(self.left_frame, bg=BG_COLOR)
        ff_frame.pack(fill=tk.X)

        tk.Label(ff_frame, text="kS (V)", bg=BG_COLOR).pack(side=tk.LEFT, expand=True)
        self.entries['ks'] = self._create_entry(ff_frame, side=tk.LEFT, default="0.0")

        tk.Label(ff_frame, text="kV (V⋅s/rad)", bg=BG_COLOR).pack(side=tk.LEFT, expand=True)
        self.entries['kv_ff'] = self._create_entry(ff_frame, side=tk.LEFT, default="0.0")

        tk.Label(ff_frame, text="kA (V⋅s²/rad)", bg=BG_COLOR).pack(side=tk.LEFT, expand=True)
        self.entries['ka'] = self._create_entry(ff_frame, side=tk.LEFT, default="0.0")

        tk.Label(self.left_frame, text="kG (V, vertical only)", bg=BG_COLOR).pack(anchor="w")
        self.entries['kg'] = self._create_entry(self.left_frame, default="0.0")

        # --- E: Resources Button ---
        tk.Button(self.left_frame, text="Resources & Help", command=self.open_resources, bg="#e6f7ff", fg=ACCENT_COLOR, font=FONT_MAIN).pack(fill=tk.X, pady=20)

        # --- Matplotlib Setup ---
        self.fig = Figure(figsize=(5, 8), dpi=100)
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        self.fig.tight_layout(pad=3.0)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _create_header(self, parent, text):
        tk.Label(parent, text=text, font=FONT_BOLD, fg=ACCENT_COLOR, bg=BG_COLOR).pack(anchor="w", pady=(10, 2))

    def _create_entry(self, parent, row=None, col=None, side=None, default=""):
        entry = tk.Entry(parent, font=FONT_MAIN, width=10)
        if default: entry.insert(0, default)
        entry.bind("<KeyRelease>", self.on_ui_change)
        
        if row is not None and col is not None:
            entry.grid(row=row, column=col, padx=2, pady=2)
        elif side is not None:
            entry.pack(side=side, padx=2)
        else:
            entry.pack(fill=tk.X, pady=2)
        return entry

    def _create_labeled_entry(self, parent, label_text, key):
        tk.Label(parent, text=label_text, font=FONT_MAIN, bg=BG_COLOR).pack(anchor="w")
        self.entries[key] = self._create_entry(parent)

    def on_layout_change(self, event=None):
        self.update_load_inputs_visibility()
        self.update_motor_inputs_visibility()
        self.run_simulation()

    def on_ui_change(self, event=None):
        self.run_simulation()

    def update_load_inputs_visibility(self):
        mode = self.mode_var.get()
        if "Load-Bearing" in mode:
            self.load_frame.pack(fill=tk.X, pady=5, after=self.mode_dropdown)
        else:
            self.load_frame.pack_forget()
        self.update_motor_b4iii_layout()

    def update_motor_inputs_visibility(self):
        m_type = self.motor_type_var.get()
        if m_type == "Custom...":
            self.custom_motor_frame.pack(fill=tk.X, pady=5, after=self.motor_dropdown)
        else:
            self.custom_motor_frame.pack_forget()

    def update_motor_b4iii_layout(self):
        for widget in self.b4iii_frame.winfo_children():
            widget.destroy()

        is_load = "Load-Bearing" in self.mode_var.get()

        if not is_load:
            tk.Label(self.b4iii_frame, text="Velocity Constant K_v (RPM/V)", font=FONT_MAIN, bg=BG_COLOR).grid(row=0, column=0)
            self.entries['kv'] = self._create_entry(self.b4iii_frame, row=1, col=0)
            tk.Label(self.b4iii_frame, text="--Or--", bg=BG_COLOR).grid(row=1, column=1)
            tk.Label(self.b4iii_frame, text="Torque Constant (Nm/A)", font=FONT_MAIN, bg=BG_COLOR).grid(row=0, column=2)
            self.entries['kt_direct'] = self._create_entry(self.b4iii_frame, row=1, col=2)
        else:
            tk.Label(self.b4iii_frame, text="Empirical Stall Torque (Nm)", font=FONT_MAIN, bg=BG_COLOR).grid(row=0, column=0)
            self.entries['stall_torque'] = self._create_entry(self.b4iii_frame, row=1, col=0)
            tk.Label(self.b4iii_frame, text="Empirical Stall Current (A)", font=FONT_MAIN, bg=BG_COLOR).grid(row=2, column=0)
            self.entries['stall_current'] = self._create_entry(self.b4iii_frame, row=3, col=0)
            tk.Label(self.b4iii_frame, text="--Or--", bg=BG_COLOR).grid(row=1, column=1, rowspan=2)
            tk.Label(self.b4iii_frame, text="Torque Constant (Nm/A)", font=FONT_MAIN, bg=BG_COLOR).grid(row=0, column=2)
            self.entries['kt_direct'] = self._create_entry(self.b4iii_frame, row=1, col=2)

    def get_val(self, key, default=0.0):
        if key not in self.entries: return default
        try:
            val = float(self.entries[key].get())
            return val
        except ValueError:
            return default

    # --- Simulation Logic ---
    def run_simulation(self):
        V_supply = 12.0

        # 1. Gather Inputs
        mode = self.mode_var.get()
        motor_type = self.motor_type_var.get()
        pid_mode = self.pid_type_var.get()

        kS = self.get_val('ks')
        kV_ff = self.get_val('kv_ff')
        kA = self.get_val('ka')
        kG = self.get_val('kg')

        
        # --- NEW: Gear Ratio (Defaults to 1 if empty/zero) ---
        gear_ratio = self.get_val('gear_ratio', 1.0)
        if gear_ratio <= 0: gear_ratio = 1.0

        # Load Parameters
        loadDistance = 0.0
        loadMass = 0.0
        loadInertia = 0.0
        
        if "Load-Bearing" in mode:
            direct_inertia = self.get_val('load_inertia_direct', -1.0)
            if direct_inertia >= 0:
                loadInertia = direct_inertia
                loadDistance = self.get_val('load_dist')
                loadMass = self.get_val('load_mass')
            else:
                loadDistance = self.get_val('load_dist')
                loadMass = self.get_val('load_mass')
                loadInertia = loadMass * (loadDistance**2)

        # Motor Parameters (Base values at Motor Shaft)
        motorInertia = 0.0
        iMax = 0.0
        Kt = 0.0
        viscousDamping = 0.0
        coulombFriction = 0.0

        if motor_type == "NEO":
            motorInertia = 0.002      # kg·m² (approx, reflected rotor)
            iMax = 40.0               # A
            R = 0.114                 # Ohms
            Kt = 0.0252               # Nm/A
            Ke = 0.0252               # V·s/rad  (SI: Ke == Kt)
            coulombFriction = 0.0454  # Nm
            viscousDamping = 0.000076

        elif motor_type == "Kraken X60":
            motorInertia = 0.00005
            iMax = 60.0
            R = 0.048
            Kt = 0.035
            Ke = 0.035
            coulombFriction = 0.06
            viscousDamping = 0.0001

        else: # Custom
            motorInertia = self.get_val('motor_inertia')
            iMax = self.get_val('imax')
            free_current = self.get_val('free_current')
            
            # Kt Calculation
            kt_direct = self.get_val('kt_direct', -1.0)
            if "Load-Bearing" not in mode:
                if kt_direct > 0:
                    Kt = kt_direct
                else:
                    kv = self.get_val('kv')
                    if kv > 0:
                        Kt = 60.0 / (kv * 2 * np.pi)
            else:
                if kt_direct > 0:
                    Kt = kt_direct
                else:
                    stall_torque = self.get_val('stall_torque')
                    stall_current = self.get_val('stall_current')
                    denom = stall_current - free_current
                    if denom > 0:
                        Kt = stall_torque / denom
            
            Ke = Kt  # SI Units

            # --- Resistance ---
            R_direct = self.get_val('motor_resistance', -1.0)

            if R_direct > 0:
                R = R_direct
            else:
                stall_current = self.get_val('stall_current', -1.0)
                if stall_current > 0:
                    R = V_supply / stall_current
                else:
                    # Safe fallback (prevents divide-by-zero explosions)
                    R = 0.05
            
            # Friction/Damping Base
            coulombFriction = Kt * free_current
            
            d_direct = self.get_val('viscous_damping_direct', -1.0)
            if d_direct >= 0:
                viscousDamping = d_direct
            else:
                free_speed_rpm = self.get_val('free_speed')
                if free_speed_rpm > 0:
                    rad_s = free_speed_rpm * 2 * np.pi / 60.0
                    viscousDamping = coulombFriction / rad_s

        # --- APPLY GEARBOX TRANSFORMATIONS ---
        G = gear_ratio

        # Reflect motor parameters to output shaft
        Kt = Kt * G
        Ke = Ke * G
        motorInertia = motorInertia * (G**2)
        loadInertia = loadInertia # load inertia is already at output shaft

        coulombFriction = coulombFriction * G
        viscousDamping = viscousDamping * (G**2)


        # PID Coeffs
        Kp = self.get_val('kp')
        Ki = self.get_val('ki')
        Kd = self.get_val('kd')

        # Targets
        raw_targets = self.target_text.get("1.0", tk.END).strip().split('\n')
        targets = [] 
        try:
            for line in raw_targets:
                if ',' in line:
                    val, t = map(float, line.split(','))
                    targets.append((val, t))
            targets.sort(key=lambda x: x[1])
        except ValueError:
            pass 
        
        if not targets:
            targets = [(0.0, 0.0)]

        acceptable_error = self.get_val('error_bound')

        # --- C. Physics Simulation ---
        max_time = targets[-1][1] + 2.0 if targets else 5.0
        dt = 0.001 
        times = np.arange(0, max_time, dt)
        
        theta = 0.0
        omega = 0.0
        alpha = 0.0
        
        integral_error = 0.0
        prev_error = 0.0
        
        pos_log = []
        speed_log = []
        current_log = []
        
        target_idx = 0
        
        for t in times:
            while target_idx < len(targets) - 1 and t >= targets[target_idx+1][1]:
                target_idx += 1
            setpoint = targets[target_idx][0]
            
            # C1. Control Calculation
            current_val = theta if pid_mode == 'Position PID' else omega
            error = setpoint - current_val

            derivative = (error - prev_error) / dt
            prev_error = error

            # --- Feedforward ---
            V_ff = 0.0

            # Static + velocity
            if abs(omega) > 1e-4:
                V_ff += kS * np.sign(omega)
            V_ff += kV_ff * omega
            V_ff += kA * alpha   # uses previous alpha (correct)

            # Gravity feedforward (vertical mode only)
            if mode == 'Vertical Load-Bearing PID Mode':
                V_ff += kG * math.cos(theta)

            # --- PID (anti-windup) ---
            V_pid = (Kp * error) + (Kd * derivative)

            if abs(V_pid + V_ff) < V_supply:
                integral_error += error * dt

            V_pid += Ki * integral_error

            # Total command
            V_cmd = V_ff + V_pid
            V_cmd = max(-V_supply, min(V_supply, V_cmd))


            
            # C2. Torque Calculation
            eps = 1e-3
            friction = coulombFriction * np.tanh(omega / eps) + viscousDamping * omega
            
            net_torque = 0.0
            total_inertia = 0.0

            # Electrical model
            back_emf = Ke * omega
            current = (V_cmd - back_emf) / R

            # Current limiting
            current = max(-iMax, min(iMax, current))

            motor_torque = Kt * current

            
            if mode == 'Free PID Mode':
                total_inertia = motorInertia
                net_torque = motor_torque - friction
            elif mode == 'Horizontal Load-Bearing PID Mode':
                total_inertia = motorInertia + loadInertia
                net_torque = motor_torque - friction
            elif mode == 'Vertical Load-Bearing PID Mode':
                gravity = loadMass * 9.80665 * loadDistance * math.cos(theta)
                total_inertia = motorInertia + loadInertia
                net_torque = motor_torque - gravity - friction

            if total_inertia <= 1e-9: total_inertia = 1e-9

            # C3. Acceleration
            alpha = net_torque / total_inertia
            
            # C4. Integration
            omega = omega + alpha * dt
            theta = theta + omega * dt

            
            pos_log.append(theta)
            speed_log.append(omega)
            current_log.append(current)

        # --- D. Plotting ---
        self.ax1.clear()
        self.ax2.clear()
        
        y_data = pos_log if pid_mode == 'Position PID' else speed_log
        y_label = "Angular Displacement (rad)" if pid_mode == 'Position PID' else "Angular Speed (rad/s)"
        
        self.ax1.plot(times, y_data, color='blue', label='Actual')
        
        for i, (val, start_t) in enumerate(targets):
            end_t = targets[i+1][1] if i < len(targets) - 1 else max_time
            
            self.ax1.hlines(val, start_t, end_t, colors='red', linestyles='dotted')
            self.ax1.text(start_t, val, f"Target {i+1}: {val}", color='red', fontsize=8, verticalalignment='bottom')
            
            if acceptable_error > 0:
                self.ax1.fill_between([start_t, end_t], val - acceptable_error, val + acceptable_error, color='green', alpha=0.2)

            start_idx = int(start_t / dt)
            end_idx = int(end_t / dt)
            start_idx = max(0, min(start_idx, len(y_data)))
            end_idx = max(0, min(end_idx, len(y_data)))
            
            if start_idx < end_idx:
                segment = np.array(y_data[start_idx:end_idx])
                is_ok = np.abs(segment - val) <= acceptable_error
                first_ok_idx = -1
                for k in range(len(is_ok) - 1, -1, -1):
                    if not is_ok[k]:
                        first_ok_idx = k + 1
                        break
                else:
                    first_ok_idx = 0
                
                if first_ok_idx < len(is_ok):
                    ok_time_rel = first_ok_idx * dt
                    ok_time_abs = start_t + ok_time_rel
                    self.ax1.axvline(ok_time_abs, color='magenta', linestyle='dotted')
                    self.ax1.text(ok_time_abs, self.ax1.get_ylim()[0], f"Ok at {ok_time_rel:.2f}s", color='magenta', fontsize=8, rotation=90, verticalalignment='bottom')

        self.ax1.set_ylabel(y_label)
        self.ax1.set_xlabel("Time (s)")
        self.ax1.grid(True, linestyle='--', alpha=0.6)
        self.ax1.margins(y=0.1)
        
        self.ax2.plot(times, current_log, color='orange')
        self.ax2.set_ylabel("Current (A)")
        self.ax2.set_xlabel("Time (s)")
        limit = iMax if iMax > 0 else 1.0
        self.ax2.set_ylim([-limit, limit])
        self.ax2.set_xlim([0, max_time])
        self.ax2.grid(True, linestyle='--', alpha=0.6)

        self.canvas.draw()

    def open_resources(self):
        res_win = tk.Toplevel(self.root)
        res_win.title("Resources")
        res_win.geometry("800x600")
        
        text_area = tk.Text(res_win, wrap=tk.WORD, padx=10, pady=10, font=("Space Grotesk", 10))
        text_area.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = tk.Scrollbar(res_win, command=text_area.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        text_area.config(yscrollcommand=scrollbar.set)

        content = """
PID Control Resources

1. An introduction to PID Control
--------------------------------------------------
The PID controller uses three terms to correct error:

- K_p (Proportional): " The Present." This term produces an output proportional to the current error. If you are far from the target, it pushes hard. If K_p is too high, the motor will oscillate.

- K_i (Integral): "The Past." This term sums up the error over time. If the motor is close to the target but stuck (due to friction or gravity), K_i builds up power to push it over the hump.

- K_d (Derivative): "The Future." This term looks at how fast the error is changing. It acts as a damper. If the motor is rushing toward the target quickly, K_d reduces power to prevent overshooting.

2. How Motor Constants Affect PID
--------------------------------------------------
- Moment of Inertia (kg·m²): Resistance to rotational acceleration. A higher rotational inertia means the motor is harder to speed up or slow down.
  Calculation: J = Mass * (Distance from axis)². This quantity is possessed by both the motor's rotor and any attached load.
- Torque Constant (Nm/A): Relates input current to output torque. A higher torque constant means more torque (turning ability) for the same current.
  Calculation: K_t = 60 / (K_v * 2 * pi)
- Velocity Constant (RPM/V): Relates input voltage to output speed. A lower velocity constant means the motor spins slower for the same voltage, resulting in higher torque. This value is often given on motor spec sheets.
- Coulomb Friction (Nm): A constant torque opposing motion, representing static friction in the motor and load. This value is calculated empirically.
  Calculation: F_c = K_t * Free Current
- Viscous Damping (Nm*s/rad): Represents air resistance and bearing grease friction. This value is calculated empirically.
  Calculation: D = CoulombFriction / FreeSpeed_rad_s

3. Gearboxes and Vertical Compensation with the Feedforward
--------------------------------------------------
Adding a gearbox changes the effective motor constants at the output shaft. A gear ratio is defined as G:1 (e.g., 5:1 means the motor turns 5 times for every 1 turn of the output shaft).

The transformations are:
- Torque Constant: K_t_new = K_t_motor * G
- Velocity Constant: K_v_new = K_v_motor / G
- Motor Inertia: J_new = J_motor * G²

In vertical applications, gravity creates a torque that varies with angle. The feedforward term kG helps counteract this by providing a voltage proportional to the cosine of the angle.
Other parts of the feedforward (kS, kV, kA) help the motor overcome static friction, maintain speed, and accelerate effectively, respectively.

4. Going Deeper: The Control Loop
--------------------------------------------------
The control loop equation is:

r(t) = K_p * e(t) + K_i * ∫ e(τ) dτ + K_d * de(t)/dt

Where:
e(t) = Target - CurrentValue
r(t) = The voltage (Volts) requested by the controller.

5. Getting Real
--------------------------------------------------
To emulate the real world, this simulator uses Differential Equations (specifically Newton's 2nd Law for Rotation):

Sum of Torques = Inertia * Angular Acceleration

We account for:
1. Motor Torque: Provided by electricity.
2. Coulomb Friction: Constant drag (sliding friction).
3. Viscous Damping: Drag that increases with speed.
4. Gravity (in vertical mode): A restoring force based on the angle.
5. Voltage/Current Limits: Motors cannot accept infinite power.

6. Warnings & Disclaimers
--------------------------------------------------
- This is a simulation, not reality. Real motors have inductance, backlash, and noise which are not perfectly modeled here.
- Empirical values (custom inputs) are notoriously difficult to measure accurately without specialized equipment.
- Do not use these calculated PID coefficients on a large industrial machine without starting conservative and tuning manually.
        """
        text_area.insert(tk.END, content)
        text_area.config(state=tk.DISABLED)

if __name__ == "__main__":
    root = tk.Tk()
    app = PIDSimulatorApp(root)
    root.mainloop()