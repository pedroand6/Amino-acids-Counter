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
        arqName.value = e.files[0].name
        arquivo.value = e.files[0].path
        page.update()

    file_picker = ft.FilePicker(on_result=on_dialog_result)
    page.overlay.append(file_picker)
    page.update()
    
    arquivo = ft.Text("")
    arqName = ft.Text("Upload a CSV file", size=15, max_lines=1, overflow=ft.TextOverflow.ELLIPSIS, width=page.width*0.12, text_align=ft.TextAlign.CENTER)
    
    def ViewResults(e):
        page.go("/results")
        if arqName.value == "Upload a CSV file":
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

        figure = plt.figure(figsize=(12,6))
        
        plt.style.use("seaborn-v0_8-darkgrid")

        colors = np.full(len(linhas), 'maroon')
        colors[((cdr1[0]-1)//3):(cdr1[1]//3)] = 'blue'
        colors[((cdr2[0]-1)//3):(cdr2[1]//3)] = 'orange'
        colors[((cdr3[0]-1)//3):(cdr3[1]//3)] = 'green'
        
        bars = plt.bar(linhas, contador, color=colors, linewidth=0.7)
        
        plt.title(f'Contagem de aminoácidos da proteina {id}', fontsize=18, weight='bold', color="darkblue", pad=20)
        plt.xlabel("Posição", fontsize=14, labelpad=10)
        plt.ylabel("Contagem", fontsize=14, labelpad=10)
        
        ax = plt.gca()

        # Removendo borda superior e lateral para um visual moderno
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Ajustes finais
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()

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

    submitBtn = ft.ElevatedButton(text="Submit", on_click=lambda _: getResults(pd.read_csv(f'{arquivo.value}')), width=100)
    idField = DropdownID(page)
    graphMatPlot = MatplotlibChart(plt.figure(), expand=True)
    graph = chart.ProteinChart(page)
    
    def changeBg(e):
        btn = e.control
        if btn.bgcolor != ft.colors.GREY_300:
            btn.bgcolor = ft.colors.GREY_300
            btn.update()
        else:
            btn.bgcolor = ft.colors.GREY_200
            btn.update()

    def route_change(route):
        page.views.clear()
        page.views.append(
            ft.View(
                "/",
                [
                    ft.AppBar(title=ft.Text("Amino acids Counter by position"), bgcolor=ft.colors.GREEN_700, color=ft.colors.GREY_200, center_title=True),
                    ft.Column(
                        [
                            ft.Container(
                                width=page.width*0.25,
                                height=page.height*0.5,
                                content=ft.Column(
                                    [
                                        ft.Container(
                                            content=ft.Column(
                                                [
                                                    ft.Icon(name=ft.icons.FILE_UPLOAD_ROUNDED, color=ft.colors.GREEN_400, size=page.width*0.1),
                                                    arqName
                                                ],
                                                alignment=ft.MainAxisAlignment.CENTER,
                                                horizontal_alignment=ft.CrossAxisAlignment.CENTER
                                            ),
                                            opacity=0.8,
                                            shadow=ft.BoxShadow(
                                                spread_radius=2,
                                                blur_radius=15
                                            ),
                                            width=page.width*0.2,
                                            height=page.width*0.2,
                                            bgcolor=ft.colors.GREY_200,
                                            shape=ft.BoxShape.CIRCLE,
                                            on_hover=lambda e: changeBg(e),
                                            on_click=lambda _: file_picker.pick_files()
                                        ),
                                    ],
                                    alignment=ft.MainAxisAlignment.CENTER,
                                    horizontal_alignment=ft.CrossAxisAlignment.CENTER
                                )
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    ft.CupertinoButton(text="Submit", bgcolor=ft.colors.GREEN_700, on_click=lambda e: ViewResults(e), alignment=ft.alignment.bottom_center)
                ],
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                vertical_alignment=ft.MainAxisAlignment.CENTER
            ),
        )
        if page.route == "/results":
            idField.df = pd.read_csv(f'{arquivo.value}')
            page.views.append(
                ft.View(
                    "/results",
                    [
                        ft.AppBar(title=ft.Text("Results"), bgcolor=ft.colors.GREEN_700, color=ft.colors.GREY_200, center_title=True),
                        ft.ResponsiveRow([idField, submitBtn]),
                        ft.Stack([
                            ft.Container(
                                content=ft.Row([ft.InteractiveViewer(
                                        content=graph,
                                        scale_enabled=False,
                                        pan_enabled=False,
                                    )], width=page.width*0.8, scroll=ft.ScrollMode.ALWAYS),
                                border_radius=5,
                                border=ft.Border(None, None, ft.BorderSide(1, ft.colors.GREY_400),
                                                 ft.BorderSide(1, ft.colors.GREY_400)),
                                bgcolor=ft.colors.WHITE,
                                padding=ft.padding.all(10)
                            ),
                            ft.Row([
                                ft.Card(
                                    content=ft.Column([
                                        ft.Row([
                                            ft.Text("CDR1"), ft.Icon(ft.icons.SQUARE_ROUNDED, color=ft.colors.BLUE)
                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                        ft.Row([
                                            ft.Text("CDR2"), ft.Icon(ft.icons.SQUARE_ROUNDED, color=ft.colors.ORANGE)
                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                        ft.Row([
                                            ft.Text("CDR3"), ft.Icon(ft.icons.SQUARE_ROUNDED, color=ft.colors.GREEN)
                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                    ]), width=100
                                )], alignment=ft.MainAxisAlignment.END),
                        ], width=page.width, alignment=ft.alignment.top_center),
                        ft.Card(
                            content=graphMatPlot,
                            width=page.width*0.8,
                        ),
                    ],
                    scroll=ft.ScrollMode.AUTO,
                    horizontal_alignment=ft.CrossAxisAlignment.CENTER
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