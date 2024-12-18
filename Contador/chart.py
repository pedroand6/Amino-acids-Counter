from pydoc import text
import flet as ft
import pandas as pd

class ProteinChart(ft.BarChart):

    def __init__(self, page, linhas=[], contador=[], aminoacidos=[], cdr = [(0,0), (0,0), (0,0)]):
        super().__init__()
        self.page = page
        self.linhas = linhas
        self.contador = contador
        self.aminoacidos = aminoacidos
        self.cdr = cdr

        self.height = 500
        
        self.expand_loose = True
        self.animate = ft.animation.Animation(1000, ft.AnimationCurve.EASE_IN_OUT)
        self.horizontal_grid_lines=ft.ChartGridLines(
            color=ft.Colors.GREY_300, width=1, dash_pattern=[3, 3]
        )
        self.tooltip_bgcolor=ft.Colors.with_opacity(0.5, ft.Colors.GREY_300)
        self.interactive=False
        self.expand=True
        self.border=ft.border.all(1, ft.Colors.GREY_400)
        self.left_axis=ft.ChartAxis(
            labels_size=40, title=ft.Text("Contagem de Aminoacidos"), title_size=40
        )

        self.makeGraph()

    def getBarChartGroups(self):
        barChartGroups = []
        for i in range(len(self.linhas)):
            barChartGroups.append(
                ft.BarChartGroup(
                    x=self.linhas[i], 
                    bar_rods=[
                        ft.BarChartRod(
                            from_y=0,
                            to_y=self.contador[i],
                            width=20,
                            color=ft.Colors.RED,
                            border_radius=1,
                        )],
                )
            )

        for i in barChartGroups[((self.cdr[0][0]-1)//3):(self.cdr[0][0]//3)]:
            i.bar_rods[0].color = ft.Colors.BLUE

        for i in barChartGroups[((self.cdr[1][0]-1)//3):(self.cdr[1][0]//3)]:
            i.bar_rods[0].color = ft.Colors.ORANGE

        for i in barChartGroups[((self.cdr[2][0]-1)//3):(self.cdr[2][0]//3)]:
            i.bar_rods[0].color = ft.Colors.GREEN

        return barChartGroups

    def getChartAxis(self):
        labels = []
        for i in range(len(self.linhas)):
            labels.append(
                ft.ChartAxisLabel(value=self.linhas[i], label=ft.Container(ft.Text(f"{self.linhas[i]}\n{self.aminoacidos[i][0]}", text_align=ft.TextAlign.CENTER)))
            )
        
        return ft.ChartAxis(labels=labels, labels_size=40)
    
    def makeGraph(self):
        if(len(self.contador) > 0):
            self.max_y=max(self.contador)
        self.height = 500
        self.width = 60 * (len(self.linhas)+1)
        self.bar_groups=self.getBarChartGroups()
        self.bottom_axis=self.getChartAxis()

    def before_update(self):
        self.makeGraph()
        return super().before_update()
