import flet as ft
import pandas as pd

class DropdownID(ft.Dropdown):

    def __init__(self, page, df=pd.DataFrame()):
        super().__init__()
        self.page = page
        self.df = df
        self.width = page.width * 0.5
        self.updateOptions()

    def before_update(self):
        self.updateOptions()
        return super().before_update()
    
    def updateOptions(self):
        if self.df.empty:
            return
        
        self.options.clear()
        
        proteinDic = self.df['Protein Accession'].unique()

        for id in proteinDic:
            self.options.append(ft.dropdown.Option(f"{id}"))

        if self.value == None:
            self.value = proteinDic[0]
