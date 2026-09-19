#!/usr/bin/env python3
"""
Script to extract CAR (Commissioning Activity Record) information from RST spreadsheet.
Extracts: CAR number, CAR name, APT program, start time, and duration.
Ignores CARs that are crossed out (strikethrough formatting).
"""

import pandas as pd
import openpyxl
import re
import sys
from pathlib import Path

def is_strikethrough(cell):
    """Check if a cell has strikethrough formatting."""
    if cell.font and cell.font.strike:
        return True
    return False

def extract_apt_from_name(name):
    """Extract APT program number from CAR name if present."""
    # Look for pattern like "(APT 1054)" or "APT 1039"
    match = re.search(r'\(APT\s+(\d+)\)', name)
    if match:
        return match.group(1)
    match = re.search(r'APT\s+(\d+)', name)
    if match:
        return match.group(1)
    return None

def extract_apt_from_description(description):
    """Extract APT program number from step description if present."""
    if pd.isna(description):
        return None
    # Look for pattern like "PID 1054" or "PID 1039"
    match = re.search(r'PID\s+(\d+)', str(description))
    if match:
        return match.group(1)
    return None

def parse_spreadsheet(filepath):
    """Parse the RST commissioning spreadsheet and extract CAR information."""

    # Load workbook with openpyxl to check formatting
    wb = openpyxl.load_workbook(filepath)
    ws = wb.active

    # Also read with pandas for easier data extraction
    df = pd.read_excel(filepath, sheet_name=0)
    
    # The data structure appears to be:
    # Column A (0): CAR ID with name and sometimes duration
    # Column B (1): Step number
    # Column C (2): MET
    # Column D (3): As Run Start Time (UTC)
    # Column E (4): Duration (HH:MM:SS)
    # Column F+ : RP, Step Description, etc.
    
    cars = []
    current_car = None
    
    for idx, row in df.iterrows():
        # Excel rows are 1-indexed, and there may be header rows
        excel_row = idx + 2  # Adjust if your spreadsheet has headers

        car_id_cell = row.iloc[0]  # First column

        # Check if this is a new CAR (starts with "CAR-")
        if pd.notna(car_id_cell) and isinstance(car_id_cell, str) and car_id_cell.startswith('CAR-'):

            # Check if the cell is crossed out (strikethrough)
            openpyxl_cell = ws.cell(row=excel_row, column=1)
            if is_strikethrough(openpyxl_cell):
                print(f"Skipping crossed-out CAR at row {excel_row}: {car_id_cell.split(chr(10))[0]}")
                current_car = None  # Make sure we don't process this CAR
                continue

            # Extract CAR number from cell (first line)
            lines = car_id_cell.split('\n')
            car_number = lines[0].strip()
            car_name = lines[1].strip() if len(lines) > 1 else ""
            
            # Extract APT from name
            apt_program = extract_apt_from_name(car_id_cell)
            
            # Extract duration if present in the cell
            duration_match = re.search(r'Duration:\s*(\d{2}:\d{2}:\d{2})', car_id_cell)
            duration_from_header = duration_match.group(1) if duration_match else None
            
            current_car = {
                'CAR_Number': car_number,
                'CAR_Name': car_name,
                'APT_Program': apt_program,
                'Start_Time': None,
                'Duration': duration_from_header,
                'MET': None
            }
            
        # Check for start time and duration in subsequent row
        if current_car is not None:
            # Get the step number column (column B/1)
            step = row.iloc[1] if len(row) > 1 else None
            
            # Look for the first step (usually xxx.0000)
            if pd.notna(step) and isinstance(step, str) and '.0000' in step:
                # Get MET (column C/2)
                met = row.iloc[2] if len(row) > 2 else None
                if pd.notna(met):
                    current_car['MET'] = str(met).strip()
                
                # Get start time (column D/3)
                start_time = row.iloc[3] if len(row) > 3 else None
                if pd.notna(start_time):
                    current_car['Start_Time'] = str(start_time).strip()
                
                # Get duration from row if not already in header (column E/4)
                if not current_car['Duration']:
                    duration = row.iloc[4] if len(row) > 4 else None
                    if pd.notna(duration):
                        current_car['Duration'] = str(duration).strip()
                
                # Try to extract APT from step description if not found in name
                if not current_car['APT_Program'] and len(row) > 6:
                    description = row.iloc[6]  # Step description column
                    apt_from_desc = extract_apt_from_description(description)
                    if apt_from_desc:
                        current_car['APT_Program'] = apt_from_desc
                
                # Add to list and reset
                cars.append(current_car)
                current_car = None

    wb.close()
    return cars

def format_output_table(cars):
    """Format the CAR data as a clean table."""
    
    # Create DataFrame
    df = pd.DataFrame(cars)
    
    # Reorder columns
    column_order = ['CAR_Number', 'CAR_Name', 'APT_Program', 'MET', 'Start_Time', 'Duration']
    df = df[column_order]
    
    # Fill NaN values with empty string for better display
    df = df.fillna('')
    
    return df

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_car_info.py <spreadsheet_file> [output_csv]")
        print("Example: python extract_car_info.py RST_Commissioning_CAST_AS_RUN.xlsx")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not Path(input_file).exists():
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    
    print(f"Processing: {input_file}\n")

    # Extract CAR data
    cars = parse_spreadsheet(input_file)

    # Format as table
    df = format_output_table(cars)

    # Display results
    print(f"\nFound {len(df)} CAR entries (excluding crossed-out items):\n")
    print(df.to_string(index=False))
    
    # Save to CSV if output file specified
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"\nSaved to: {output_file}")
    
    # Also save to a default output file
    default_output = Path(input_file).stem + "_CAR_Summary.csv"
    df.to_csv(default_output, index=False)
    print(f"Saved to: {default_output}")

if __name__ == "__main__":
    main()