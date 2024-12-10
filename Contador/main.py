import matplotlib
import flet as ft
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import os

def main(page: ft.Page):
    global arquivo
    arquivo = "Teste"

    def on_dialog_result(e: ft.FilePickerResultEvent):
        print("Selected files:", e.files)
        print("Selected file or directory:", e.path)

        global arquivo
        arquivo = e.files[0].name
        ftArquivo.controls.clear()
        page.update()

    file_picker = ft.FilePicker(on_result=on_dialog_result)
    page.overlay.append(file_picker)
    page.update()

    btnArquivo = ft.ElevatedButton("Choose file...", on_click=lambda _: file_picker.pick_files())

    ftArquivo = ft.Row([
                ft.Text(arquivo, size=40),
                ft.IconButton(icon=ft.icons.SEND, icon_size=20, highlight_color=ft.colors.RED, hover_color=ft.colors.GREY_500, on_click=lambda _: page.go("/results"))
            ])
    
    # df = pd.read_csv()
    # ids = df['Protein ID'].unique()
    # id = input('Insira o ID da proteína: ')

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
            )
        )
        if page.route == "/results":
            page.views.append(
                ft.View(
                    "/store",
                    [
                        ft.AppBar(title=ft.Text("Resultados"), bgcolor=ft.colors.SURFACE_VARIANT),
                        ft.Row(
                            [
                                #ft.MatplotlibChart(figure, expand=True)
                            ], width=900)
                    ],
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


def getResults(id, df):
    '''
    (int, dataframe) -> figure, dataframe
    '''

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
    plt.style.use('seaborn-darkgrid')

    plt.bar(linhas, contador, color ='maroon')
    plt.xlabel('Linha')
    plt.ylabel('Contagem')

    plt.title(f'Contagem da proteína {id}')

    plt.plot()

    return figure, result

ft.app(main)
