import flet as ft
from flet.matplotlib_chart import MatplotlibChart
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import os

def main(page: ft.Page):
    def on_dialog_result(e: ft.FilePickerResultEvent):
        arquivo.value = e.files[0].path
        page.update()

    file_picker = ft.FilePicker(on_result=on_dialog_result)
    page.overlay.append(file_picker)
    page.update()

    btnArquivo = ft.ElevatedButton("Choose file...", on_click=lambda _: file_picker.pick_files())

    arquivo = ft.Text("", size=20)
    ftArquivo = ft.Row([
                arquivo,
                ft.IconButton(icon=ft.icons.SEND, icon_size=20, highlight_color=ft.colors.RED, hover_color=ft.colors.GREY_500, on_click=lambda e: ViewResults(e))
            ])
    
    def ViewResults(e):
        if arquivo.value == "":
            return
        
        page.go("/results")

    def getResults(e, df):
        '''
        (int, dataframe) -> None
        '''
        id = e.control.value
        print(id)
        mask = (df['Protein Accession'].values == id)
        start = (df['Start'].values)[mask]
        end = (df['End'].values)[mask]

        contador = np.full(max(end)-min(start)+1, 0)

        for j in range(len(start)):
            contador[start[j]-1:end[j]] += 1

        linhas = list(range(min(start), max(end)+1))
        data = {'Line': linhas, 'Count': contador}

        result = pd.DataFrame(data=data)

        figure = plt.figure(figsize=(30,8))
        print(plt.style.available)
        #plt.style.use('seaborn-darkgrid')

        plt.bar(linhas, contador, color ='maroon')
        plt.xlabel('Linha')
        plt.ylabel('Contagem')

        plt.title(f'Contagem da proteína {id}')

        plt.plot()

        graph.figure = figure
        page.update()


    idField = ft.TextField("", hint_text="Insira o ID da cobertura", on_submit=lambda e: getResults(e, pd.read_csv(f'{arquivo.value}')))
    graph = MatplotlibChart(plt.figure(figsize=(12,4)), expand=True)

    def route_change(route):
        page.views.clear()
        page.views.append(
            ft.View(
                "/",
                [
                    ft.AppBar(title=ft.Text("Contador"), bgcolor=ft.colors.SURFACE_VARIANT),
                    btnArquivo,
                    ftArquivo
                ],
                bgcolor=ft.Colors.BLUE_200
            )
        )
        if page.route == "/results":
            page.views.append(
                ft.View(
                    "/store",
                    [
                        ft.AppBar(title=ft.Text("Resultados"), bgcolor=ft.colors.SURFACE_VARIANT),
                        idField,
                        ft.Row(
                            [
                                graph
                            ], width=900)
                    ],
                    bgcolor=ft.Colors.BLUE_200
                )
            )
        page.update()

    def view_pop(view):
        page.views.pop()
        top_view = page.views[-1]
        page.go(top_view.route)

    page.on_route_change = route_change
    page.on_view_pop = view_pop
    page.go(page.route)

ft.app(main)
