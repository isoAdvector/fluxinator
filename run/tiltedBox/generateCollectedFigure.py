import os
import subprocess
from pathlib import Path
from PIL import Image

def generate_latex_figure(output_file="figure.tex"):
    # 1. USER CONFIGURATION - Modify these variables as needed
    # ---------------------------------------------------------
    # Specify the directories to include (in the order you want them)
    solver_dirs = [
        "gravityModelScan_1deltaT",  # Will be first row
        "gravityModelScan_100deltaT",      # Second row
        "gravityModelScan_1000deltaT",     # Third row
        # Add more as needed
    ]
    
    # Column headers (one for each column)
    column_headers = [
        r"interFoam \\ basic",
        r"interIsoFoam \\ basic",
        r"interIsoGoam \\ height",
        r"interIsoGoam \\ geometric",
        r"interIsoGoam \\ gravityRecon",
    ]
    
    # Label for each row
    row_labels = [
        "T = 0.01s",
        "T = 1s",
        "T = 10s",
    ]
    
    # Scaling factors (adjust these to control relative sizes)
    width_scaling = 0.15    # Horizontal scaling factor
    height_scaling = 0.15   # Vertical scaling factor
    # ---------------------------------------------------------

    # Verify directories exist
    solver_dirs = [d for d in solver_dirs if Path(d).exists()]
    if not solver_dirs:
        print("Error: None of the specified directories exist!")
        return

    # Determine number of columns from first directory
    first_dir = Path(solver_dirs[0])
    subdirs = sorted([d for d in first_dir.iterdir() if d.is_dir() and d.name.isdigit() and len(d.name) == 4])
    num_columns = len(subdirs)

    if num_columns == 0:
        print(f"Error: No numbered subdirectories (0000, 0001, etc.) found in {first_dir}!")
        return

    # Get dimensions from first image (assuming all images are same size)
    first_img_path = first_dir / "0000" / "0000.png"
    if not first_img_path.exists():
        print(f"Error: Could not find sample image at {first_img_path}!")
        return

    try:
        with Image.open(first_img_path) as img:
            img_width_px, img_height_px = img.size
            px_to_cm = 0.026458  # Conversion factor (96dpi)
            img_width_cm = img_width_px * px_to_cm * width_scaling
            img_height_cm = img_height_px * px_to_cm * height_scaling
    except Exception as e:
        print(f"Error reading image dimensions: {e}")
        return

    # Adjust column widths
    label_column_width = 0.01
    remaining_width = 0.95 - label_column_width
    column_width = remaining_width / num_columns

    # Calculate paper size based on image dimensions and layout
    # Add margins and account for labels
    paper_width = (img_width_cm * num_columns * 1.1) + 2.0  # 10% spacing + margins
    paper_height = (img_height_cm * len(solver_dirs) * 1.2) + 2.0  # 20% spacing + margins

    # LaTeX template - carefully escaped braces
    latex_template = r"""\documentclass{{article}}
\usepackage{{graphicx}}
\usepackage{{geometry}}
\geometry{{papersize={{{width:.2f}cm,{height:.2f}cm}},total={{{width:.2f}cm,{height:.2f}cm}}}}

\begin{{document}}
\begin{{figure}}
    \centering
{content}
\end{{figure}}
\end{{document}}
"""

    # Generate content
    content = []
    
    for row_idx, solver_dir in enumerate(solver_dirs):
        solver_path = Path(solver_dir)
        
        # Get row label (use directory name if not specified)
        time_label = row_labels[row_idx]
        
        # Row label (rotated)
        row_label = r"\rotatebox[origin=c]{{90}}{{\textbf{{{}}}}}".format(time_label)
        
        # Start row with label - note unescaped \end{minipage}
        row_content = [
            r"\begin{{minipage}}[c]{{{:.3f}\linewidth}}".format(label_column_width),
            r"    \centering",
            f"    {row_label}",
            r"\end{minipage}"  # This one stays unescaped
        ]
        
        # Add columns
        for col_idx in range(num_columns):
            subdir = f"{col_idx:04d}"
            img_path = solver_path / subdir / f"{col_idx:04d}.png"
            
            if not img_path.exists():
                print(f"Warning: Missing image {img_path}")
                continue
            
            # For first row, include column headers
            header = column_headers[col_idx] if row_idx == 0 and col_idx < len(column_headers) else ""
            
            col_block = [
                r"\begin{{minipage}}[c]{{{:.3f}\linewidth}}".format(column_width),
                r"    \centering",
            ]
            
            if header:
                col_block.append(f"    \\textbf{{{header}}} \\\\")
            
            col_block.extend([
                f"    \\includegraphics[width=\\linewidth]{{{img_path}}}",
                r"\end{minipage}"  # This one stays unescaped
            ])
            
            row_content.extend(col_block)
        
        content.append("\n".join(row_content))
    
    # Combine all rows with some vertical space between them
    full_content = "\n\n".join(content)
    
    # Write to file
    with open(output_file, "w") as f:
        f.write(latex_template.format(
            width=paper_width,
            height=paper_height,
            content=full_content
        ))
    
    print(f"Successfully generated LaTeX file: {output_file}")
    print(f"Paper size: {paper_width:.2f}cm x {paper_height:.2f}cm")
    print(f"Based on image size: {img_width_cm:.2f}cm x {img_height_cm:.2f}cm")

    # Compile the LaTeX file
    try:
        print("\nCompiling LaTeX document...")
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", output_file],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if result.returncode == 0:
            print("Successfully compiled PDF document")
            
            # Open the resulting PDF
            pdf_file = output_file.replace('.tex', '.pdf')
            if Path(pdf_file).exists():
                print(f"Opening {pdf_file}...")
                try:
                    subprocess.run(["xdg-open", pdf_file], check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Could not open PDF viewer: {e}")
                except FileNotFoundError:
                    print("Could not find xdg-open - please open the PDF manually")
            else:
                print(f"Error: Expected output file {pdf_file} not found")
        else:
            print("Error compiling LaTeX document:")
            print(result.stderr)
            
    except subprocess.CalledProcessError as e:
        print(f"Error running pdflatex: {e}")
        print(e.stderr)
    except FileNotFoundError:
        print("Error: pdflatex not found. Please ensure LaTeX is installed and in your PATH")

if __name__ == "__main__":
    generate_latex_figure()