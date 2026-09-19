#!/usr/bin/env python3
"""
Script to extract CAR (Commissioning Activity Record) information from RST spreadsheet.
Extracts: CAR number, CAR name, APT program, UTC start time (ISO format), and duration.
Uses DOY from column E and as-run time from column D.
Ignores CARs that are crossed out (strikethrough formatting) or contain "analysis".
"""

import pandas as pd
import openpyxl
from openpyxl.styles import Font
import re
import sys
from pathlib import Path
from datetime import datetime, timedelta

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

def extract_doy(doy_time_cell):
    """Extract DOY from the DOY\nTime cell."""
    if pd.isna(doy_time_cell):
        return None
    
    doy_time_str = str(doy_time_cell).strip()
    # Split by newline to get DOY
    parts = doy_time_str.split('\n')
    if len(parts) >= 1:
        doy = parts[0].strip()
        try:
            return int(doy)
        except ValueError:
            return None
    return None

def doy_to_iso(doy, time_str, year=2026):
    """
    Convert DOY and time to ISO 8601 format.
    
    Args:
        doy: Day of year (1-366)
        time_str: Time string in HH:MM:SS format
        year: Year (default 2026)
    
    Returns:
        ISO formatted datetime string (YYYY-MM-DDTHH:MM:SS)
    """
    try:
        # Create a datetime object for Jan 1 of the given year
        jan_first = datetime(year, 1, 1)
        # Add the number of days (DOY - 1 because Jan 1 is DOY 1)
        target_date = jan_first + timedelta(days=int(doy) - 1)
        
        # Parse the time
        time_parts = time_str.split(':')
        hour = int(time_parts[0])
        minute = int(time_parts[1])
        second = int(time_parts[2])
        
        # Combine date and time
        full_datetime = target_date.replace(hour=hour, minute=minute, second=second)
        
        # Return ISO format
        return full_datetime.strftime('%Y-%m-%dT%H:%M:%S')
    except Exception as e:
        print(f"Warning: Could not convert DOY {doy} and time {time_str}: {e}")
        return None

def parse_spreadsheet(filepath, year=2026):
    """Parse the RST commissioning spreadsheet and extract CAR information."""
    
    # Load workbook with openpyxl to check formatting
    wb = openpyxl.load_workbook(filepath)
    ws = wb.active
    
    # Also read with pandas for easier data extraction
    df = pd.read_excel(filepath, sheet_name=0)
    
    cars = []
    current_car = None
    
    for idx, row in df.iterrows():
        # Excel rows are 1-indexed, and there may be header rows
        excel_row = idx + 2  # Adjust if your spreadsheet has headers
        
        car_id_cell = row.iloc[0]  # First column (A)
        
        # Check if this is a new CAR (starts with "CAR-")
        if pd.notna(car_id_cell) and isinstance(car_id_cell, str) and car_id_cell.startswith('CAR-'):
            
            # Check if the cell is crossed out (strikethrough)
            openpyxl_cell = ws.cell(row=excel_row, column=1)
            if is_strikethrough(openpyxl_cell):
                print(f"Skipping crossed-out CAR at row {excel_row}: {car_id_cell.split(chr(10))[0]}")
                current_car = None  # Make sure we don't process this CAR
                continue
            
            # Check if the CAR contains "analysis" (case-insensitive)
            if 'analysis' in car_id_cell.lower():
                print(f"Skipping analysis CAR at row {excel_row}: {car_id_cell.split(chr(10))[0]}")
                current_car = None
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
                'Start_Time_UTC': None,
                'Duration': duration_from_header
            }
            
        # Check for start time and duration in subsequent row
        if current_car is not None:
            # Get the step number column (column B/1)
            step = row.iloc[1] if len(row) > 1 else None
            
            # Look for the first step (usually xxx.0000)
            if pd.notna(step) and isinstance(step, str) and '.0000' in step:
                # Column layout:
                # A (0): CAR ID/Name
                # B (1): Step
                # C (2): MET
                # D (3): As Run Start Time (HH:MM:SS)
                # E (4): DOY\nTime format
                # F (5): Duration
                
                # Get DOY from column E (index 4)
                doy_cell = row.iloc[4] if len(row) > 4 else None
                doy = extract_doy(doy_cell)
                
                # Get as-run time from column D (index 3)
                as_run_time = row.iloc[3] if len(row) > 3 else None
                
                # Combine DOY and as-run time to create ISO format
                if doy is not None and pd.notna(as_run_time):
                    as_run_time_str = str(as_run_time).strip()
                    current_car['Start_Time_UTC'] = doy_to_iso(doy, as_run_time_str, year)
                
                # Get duration from column F (index 5) if not already in header
                if not current_car['Duration']:
                    duration = row.iloc[5] if len(row) > 5 else None
                    if pd.notna(duration):
                        current_car['Duration'] = str(duration).strip()
                
                # Try to extract APT from step description if not found in name
                if not current_car['APT_Program'] and len(row) > 7:
                    description = row.iloc[7]  # Step description column
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
    column_order = ['CAR_Number', 'CAR_Name', 'APT_Program', 'Start_Time_UTC', 'Duration']
    df = df[column_order]
    
    # Fill NaN values with empty string for better display
    df = df.fillna('')
    
    return df

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_car_info.py <spreadsheet_file> [output_csv] [--year YYYY]")
        print("Example: python extract_car_info.py RST_Commissioning_CAST_AS_RUN.xlsx")
        print("         python extract_car_info.py RST_Commissioning_CAST_AS_RUN.xlsx output.csv --year 2026")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = None
    year = 2026  # Default year
    
    # Parse command line arguments
    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == '--year' and i + 1 < len(sys.argv):
            year = int(sys.argv[i + 1])
            i += 2
        else:
            output_file = sys.argv[i]
            i += 1
    
    if not Path(input_file).exists():
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    
    print(f"Processing: {input_file}")
    print(f"Using year: {year}\n")
    
    # Extract CAR data
    cars = parse_spreadsheet(input_file, year)
    
    # Format as table
    df = format_output_table(cars)
    
    # Display results
    print(f"\nFound {len(df)} CAR entries (excluding crossed-out and analysis items):\n")
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