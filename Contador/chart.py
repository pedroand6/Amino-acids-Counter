import flet as ft
import pandas as pd

class ProteinChart(ft.UserControl):

    def __init__(self, page, linhas=[], contador=[], aminoacidos=[]):
        super().__init__()
        self.page = page
        self.linhas = linhas
        self.contador = contador
        self.aminoacidos = aminoacidos
        self.barChartGroups = []
    

    def getBarChartGroups(self):
        for i in range(len(self.linhas)):
            self.barChartGroups.append(
                ft.BarChartGroup(
                    x=self.linhas[i], 
                    bar_rods=[
                        ft.BarChartRod(
                            from_y=0,
                            to_y=self.contador[i],
                            width=40,
                            color=ft.Colors.AMBER,
                            tooltip=self.aminoacidos[i],
                            border_radius=1
                        ),
                ])
            )

    def getChartAxis(self):
        labels = []
        for i in range(len(self.linhas)):
            labels.append(
                ft.ChartAxisLabel(value=self.linhas[i], label=ft.Container(ft.Text(f"{self.linhas[i]}"), padding=3))
            )
        
        return ft.ChartAxis(labels=labels, labels_size=40)

    def build(self):
        self.view = ft.BarChart(
            bar_groups=self.barChartGroups,
            border=ft.border.all(1, ft.Colors.GREY_400),
            left_axis=ft.ChartAxis(
                labels_size=40, title=ft.Text("Contagem de Aminoacidos por posiçao"), title_size=40
            ),
            bottom_axis=self.getChartAxis(),
            horizontal_grid_lines=ft.ChartGridLines(
                color=ft.Colors.GREY_300, width=1, dash_pattern=[3, 3]
            ),
            tooltip_bgcolor=ft.Colors.with_opacity(0.5, ft.Colors.GREY_300),
            max_y=80,
            interactive=True,
            expand=True,
        )