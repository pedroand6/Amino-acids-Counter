from pydoc import text
import flet as ft
import pandas as pd

class DropdownID(ft.Dropdown):

    def __init__(self, page, df=pd.DataFrame()):
        super().__init__()
        self.page = page
        self.df = df
        self.updateOptions()

    def before_update(self):
        self.updateOptions()
        return super().before_update()
    
    def updateOptions(self):
        if self.df.empty:
            return
        
        self.options.clear()
        
        for id in self.df['Protein Accession'].unique():
            self.options.append(ft.dropdown.Option(f"{id}"))

        self.value = self.options[1]
