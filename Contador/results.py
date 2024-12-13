import flet as ft
import pandas as pd

class Sequencia(ft.UserControl):

    def __init__(self, page, path):
        super().__init__()
        self.page = page
        self.path = path
        self.df = pd.read_csv(path)
        self.ids = self.df['Protein ID'].unique()

    def build(self):
        self.view = 
        return self.view