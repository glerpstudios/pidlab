# Python Tkinter + Matplotlib App

This project uses **Tkinter** for the GUI and **Matplotlib** with **NumPy** for plotting.

## Setup Instructions

1. **Clone or download this repository.**

2. **Make sure you have Python installed.**  
   Tkinter is included with most standard Python installations.

   On some Linux distributions :
   ```
   sudo apt-get install python3-tk
   ```
   On Mac OS X :
   ```
   brew install python-tk
    brew install tcl-tk
    env PATH="/usr/local/opt/tcl-tk/bin:$PATH" \
        LDFLAGS="-L/usr/local/opt/tcl-tk/lib" \
        CPPFLAGS="-I/usr/local/opt/tcl-tk/include" \
        PKG_CONFIG_PATH="/usr/local/opt/tcl-tk/lib/pkgconfig" \
        pyenv install 3.12.5
   ```

3. **Create and activate a virtual environment (recommended).**

   On macOS/Linux:
   ```
   python3 -m venv venv
   source venv/bin/activate
   pip install --upgrade pip
   pip install -r requirements.txt
   python -m tkinter
   ```

   On Windows:
   ```
   python -m venv venv
   venv\Scripts\activate
   pip install --upgrade pip
   pip install -r requirements.txt
   python -m tkinter
   ```