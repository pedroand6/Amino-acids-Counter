from math import exp
import flet as ft
from flet.matplotlib_chart import MatplotlibChart
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import chart
import os
import re

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
        peptidios = (df['Peptide'].values)[mask]

        contador = np.full(max(end)-min(start)+1, 0)

        for j in range(len(start)):
            contador[start[j]-1:end[j]] += 1

        linhas = list(range(min(start), max(end)+1))
        data = {'Line': linhas, 'Count': contador}

        aminoacidos = np.full((max(end)-min(start)+1, 1), '')
        for i in range(len(peptidios)):
            thisAminoacidos = re.sub(r'\([^)]*\)', '', peptidios[i]).split(".")
            for j in thisAminoacidos:
                if len(j) > 1:
                    #print(j end[i] - start[i] + 1)
                    aminoacidos[start[i]-1:end[i]] = np.array([list(word) for word in j])


        result = pd.DataFrame(data=data)

        figure = plt.figure(figsize=(30,24))
        print(plt.style.available)

        plt.bar(linhas, contador, color ='maroon')
        plt.xlabel('Linha')
        plt.ylabel('Contagem')

        plt.title(f'Contagem de Aminoacidos por posiçao da proteina {id}')

        plt.plot()

        #graph.figure = figure
        graph.linhas = linhas
        graph.contador = contador
        graph.aminoacidos = list(aminoacidos)
        graph.build()
        page.update()


    idField = ft.TextField("", hint_text="Insira o ID da proteina", on_submit=lambda e: getResults(e, pd.read_csv(f'{arquivo.value}')))
    #graph = MatplotlibChart(plt.figure(), expand=True)
    graph = chart.ProteinChart(page, [1,2,3], [1,2,3], ['A', 'B', 'C'])

    def route_change(route):
        page.views.clear()
        page.views.append(
            ft.View(
                "/",
                [
                    ft.AppBar(title=ft.Text("Contador de Aminoacidos por posiçao"), bgcolor=ft.colors.SURFACE_VARIANT),
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
