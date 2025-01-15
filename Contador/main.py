import time
import flet as ft
from flet.matplotlib_chart import MatplotlibChart
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import chart
from dropdownID import DropdownID
import re
from crowelab_pyir import PyIR
from reverseTranslator import aa2na

def main(page: ft.Page):
    def on_dialog_result(e: ft.FilePickerResultEvent):
        arquivo.value = e.files[0].path
        arqName.value = e.files[0].name
        page.update()

    file_picker = ft.FilePicker(on_result=on_dialog_result)
    page.overlay.append(file_picker)
    page.update()

    btnArquivo = ft.ElevatedButton("Choose file...", on_click=lambda _: file_picker.pick_files())

    arquivo = ft.Text("", size=20)
    arqName = ft.Text("", size=20)
    ftArquivo = ft.Row([
                arqName,
                ft.IconButton(icon=ft.icons.SEND, icon_size=20, highlight_color=ft.Colors.RED, hover_color=ft.Colors.GREY_500, on_click=lambda e: ViewResults(e))
            ])
    
    def ViewResults(e):
        page.go("/results")
        if arquivo.value == "":
            return
        
        page.go("/results")

    def getResults(df):
        '''
        (int, dataframe) -> None
        '''
        page.go("/loading")
        
        id = idField.value
        mask = (df['Protein Accession'].values == id)
        start = (df['Start'].values)[mask]
        end = (df['End'].values)[mask]
        peptidios = (df['Peptide'].values)[mask]

        contador = np.full(max(end)-min(start)+1, 0)

        for j in range(len(start)):
            contador[start[j]-1:end[j]] += 1

        linhas = list(range(min(start), max(end)+1))
        data = {'Line': linhas, 'Count': contador}

        aminoacidos = np.full((max(end), 1), '')
        for i in range(len(peptidios)):
            thisAminoacidos = re.sub(r'\([^)]*\)', '', peptidios[i]).split(".")
            for j in thisAminoacidos:
                if len(j) > 1:
                    aminoacidos[start[i]-1:end[i]] = np.array([list(word) for word in j])


        result = pd.DataFrame(data=data)
        aminStr = ''
        for item in aminoacidos:
            aminStr += item[0]

        #Creating a fasta file to open
        fastaFile = "Sequence.fasta"
        ofile = open(fastaFile, "w")
        ofile.write(">" + id + "\n" + aa2na(aminStr) + "\n")

        #do not forget to close it
        ofile.close()

        #Connecting to IGBlast to get the CDRs
        pyirfile = PyIR(query=fastaFile, args=['--outfmt', 'dict'])
        result = pyirfile.run()

        #Get results
        try:
            cdr1 = (int(list(result.values())[0]['cdr1_start']), int(list(result.values())[0]['cdr1_end']))
        except:
            cdr1 = (0, 0)
        
        try:
            cdr2 = (int(list(result.values())[0]['cdr2_start']), int(list(result.values())[0]['cdr2_end']))
        except:
            cdr2 = (0, 0)
            
        try:
            cdr3 = (int(list(result.values())[0]['cdr3_start']), int(list(result.values())[0]['cdr3_end']))
        except:
            cdr3 = (0, 0)

        figure = plt.figure(figsize=(30,16))

        colors = np.full(len(linhas), 'maroon')
        colors[((cdr1[0]-1)//3):(cdr1[1]//3)] = 'blue'
        colors[((cdr2[0]-1)//3):(cdr2[1]//3)] = 'orange'
        colors[((cdr3[0]-1)//3):(cdr3[1]//3)] = 'green'

        plt.bar(linhas, contador, color = colors)
        plt.xlabel('Linha')
        plt.ylabel('Contagem')

        plt.title(f'Contagem de Aminoacidos por posiçao da proteina {id}')

        plt.plot()

        graphMatPlot.figure = figure

        graph.linhas = linhas
        graph.contador = contador
        graph.aminoacidos = list(aminoacidos)
        graph.cdr = [cdr1, cdr2, cdr3]
        graph.build()
        
        page.go("/results")
        time.sleep(0.05)
        page.go("/loading")
        time.sleep(0.05)
        page.go("/results")
        page.update()

    submitBtn = ft.ElevatedButton(text="Submit", on_click=lambda _: getResults(pd.read_csv(f'{arquivo.value}')))
    idField = DropdownID(page)
    graphMatPlot = MatplotlibChart(plt.figure(), expand=True)
    graph = chart.ProteinChart(page)

    def route_change(route):
        page.views.clear()
        page.views.append(
            ft.View(
                "/",
                [
                    ft.AppBar(title=ft.Text("Contador de Aminoacidos por posiçao"), bgcolor=ft.Colors.ON_SURFACE_VARIANT),
                    btnArquivo,
                    ftArquivo
                ]
            )
        )
        if page.route == "/results":
            idField.df = pd.read_csv(f'{arquivo.value}')
            page.views.append(
                ft.View(
                    "/results",
                    [
                        ft.AppBar(title=ft.Text("Resultados"), bgcolor=ft.Colors.ON_SURFACE_VARIANT),
                        idField,
                        submitBtn,
                        ft.Stack([
                            ft.Row([ft.InteractiveViewer(
                                        boundary_margin=ft.margin.only(0, 50, float(graph.width), 50),
                                        content=graph,
                                        scale_enabled=False,
                                        pan_enabled=False
                                    )], width=page.width, expand=True, scroll=ft.ScrollMode.ALWAYS),
                            ft.Row([
                                ft.Card(
                                    content=ft.Column([
                                        ft.Row([
                                            ft.Text("CDR1"), ft.Icon(ft.icons.SQUARE_ROUNDED, color=ft.Colors.BLUE)
                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                        ft.Row([
                                            ft.Text("CDR2"), ft.Icon(ft.icons.SQUARE_ROUNDED, color=ft.Colors.ORANGE)
                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                        ft.Row([
                                            ft.Text("CDR3"), ft.Icon(ft.icons.SQUARE_ROUNDED, color=ft.Colors.GREEN)
                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                    ]), width=100
                                )], alignment=ft.MainAxisAlignment.END),
                        ], width=page.width, expand=True),
                        ft.Card(
                            content=graphMatPlot,
                            width=page.width
                        ),
                    ],
                    scroll=ft.ScrollMode.AUTO
                )
            )
        if page.route == "/loading":
            page.views.append(
                ft.View(
                    "/loading",
                    [
                        ft.Container(height=page.height//2),
                        ft.ProgressRing(),
                        ft.Container(height=page.height//2),
                    ],
                    horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                    vertical_alignment=ft.MainAxisAlignment.CENTER,
                    scroll=ft.ScrollMode.HIDDEN
                )
            )
        page.update()

    def view_pop(view):
        page.views.pop()
        top_view = page.views[-1]
        page.go(top_view.route)
        
    def updateView(view):
        route = page.route
        page.go("/loading")
        time.sleep(0.05)
        page.go(route)
        page.update()

    page.on_route_change = route_change
    page.on_view_pop = view_pop
    page.on_resized = updateView
    page.go(page.route)

ft.app(main)