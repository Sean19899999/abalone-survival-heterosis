"""Read the revision workbook or the equivalent public CSV snapshot."""
from pathlib import Path
import pandas as pd

def read_sheet(workbook, sheet_name, **kwargs):
    workbook = Path(workbook)
    if workbook.is_file():
        return pd.read_excel(workbook, sheet_name=sheet_name, **kwargs)
    csv_path = workbook.parent / 'data' / f'{sheet_name}.csv'
    if not csv_path.is_file():
        raise FileNotFoundError(f'Missing workbook or exported sheet: {csv_path.name}')
    return pd.read_csv(csv_path, **kwargs)
