from math import exp
import time
import flet as ft
from flet.matplotlib_chart import MatplotlibChart
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import chart
from dropdownID import DropdownID
import re
from crowelab_pyir import PyIR
from reverseTranslator import aa2na
import os
import tempfile

temp_files = []

def main(page: ft.Page):
    page.fonts = {
        "Segoe UI": "/fonts/SegoeUI.ttf",
        "Cinzel": "/fonts/Cinzel-Regular.ttf",
        "Cinzel-Bold": "/fonts/Cinzel-Bold.ttf",
    }

    page.theme_mode = ft.ThemeMode.LIGHT
    page.theme = ft.Theme(font_family="Segoe UI", color_scheme_seed=ft.colors.GREEN)
    os.environ["FLET_SECRET_KEY"] = os.urandom(12).hex()

    file_mapping = {}

    def on_dialog_result(e: ft.FilePickerResultEvent):
        if e.page.web and e.files:
            file_mapping.clear()
            # Obtém a URL para upload do arquivo para o backend
            file_name = e.files[0].name
            file_mapping[file_name] = os.path.basename(tempfile.NamedTemporaryFile(delete=False, suffix=".csv").name)

            upload_url = page.get_upload_url(
                file_mapping[file_name], 
                3600)
            
            # Faz o upload do arquivo
            upload_file = ft.FilePickerUploadFile(file_name, upload_url)
            print(upload_file.upload_url)
            file_picker.upload([upload_file])
        elif e.files:
            arqName.value = e.files[0].name
            arquivo.value = e.files[0].path

    def on_upload_complete(e: ft.FilePickerUploadEvent):
        if file_picker.result != None and file_picker.result.files != None:
            custom_name = file_mapping.get(e.file_name)
            arqName.value = custom_name
            arquivo.value = f"{os.getcwd()+os.sep}Counter{os.sep}src{os.sep}uploads{os.sep+str(custom_name)}"
            page.update()

    def exportGraph(e):
        with tempfile.NamedTemporaryFile(delete=False, dir=os.path.join("Counter", "src", "assets"), suffix=".pdf") as tmp:
            graphMatPlot.figure.savefig(tmp.name, dpi=500)
            page.launch_url(f"{os.path.basename(tmp.name)}")
            temp_files.append(tmp.name)
    
    file_picker = ft.FilePicker(on_result=on_dialog_result, on_upload=on_upload_complete)
    page.overlay.append(file_picker)
    page.update()
    
    arquivo = ft.Text("")
    arqName = ft.Text("Upload a CSV file", width=page.width*0.1, max_lines=1, overflow=ft.TextOverflow.ELLIPSIS, text_align=ft.TextAlign.CENTER)
    
    def ViewResults(e):
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
        try:
            pyirfile = PyIR(query=fastaFile, args=['--outfmt', 'dict'])
            result = pyirfile.run()
        except Exception as e:
            result = {}
            dismissed = False

            def dismissError(e): e = True

            alerta = ft.AlertDialog(modal=False, on_dismiss=lambda _: dismissError(dismissed), title=ft.Text("Error to connect to IGBlast"), content=ft.Text(f"No CDR information can be displayed.\n{e}", selectable=True))
            page.open(alerta)
            while not dismissed: continue

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

        colors = np.full(len(linhas), 'orange')
        colors[((cdr1[0]-1)//3):(cdr1[1]//3)] = 'blue'
        colors[((cdr2[0]-1)//3):(cdr2[1]//3)] = 'maroon'
        colors[((cdr3[0]-1)//3):(cdr3[1]//3)] = 'green'
        
        bars = plt.bar(linhas, contador, color=colors, linewidth=0.7)
        
        plt.title(f'Amino acids Count for protein {id}', fontsize=18, weight='bold', color="darkblue", pad=20)
        plt.xlabel("Position", fontsize=14, labelpad=10)
        plt.ylabel("Spectra Count", fontsize=14, labelpad=10)
        
        ax = plt.gca()

        # Removendo borda superior e lateral para um visual moderno
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Ajustes finais
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        
        legend = {
            'Constant Region': 'orange',
            'CDR1': 'blue',
            'CDR2': 'maroon',
            'CDR3': 'green'
        }
        labels = list(legend.keys())
        handles = [mpatches.Rectangle((0,0),1,1, color=legend[label]) for label in labels]
        plt.legend(handles, labels, fontsize=12)

        plt.plot()

        graphMatPlot.figure = figure

        graph.linhas = linhas
        graph.contador = contador
        graph.aminoacidos = list(aminoacidos)
        graph.cdr = [cdr1, cdr2, cdr3]
        graph.build()

        minValue.data = min(contador)
        maxValue.data = max(contador)

        minValue.value = f"Min count: {minValue.data} \nPosition: {list(linhas)[list(contador).index(minValue.data)]}"
        maxValue.value = f"Max count: {maxValue.data} \nPosition: {list(linhas)[list(contador).index(maxValue.data)]}"
        
        page.go("/results")
        time.sleep(0.05)
        page.go("/loading")
        time.sleep(0.05)
        page.go("/results")
        page.update()

    submitBtn = ft.ElevatedButton(text="Submit", on_click=lambda _: getResults(pd.read_csv(f'{arquivo.value}')), width=100, 
                                  elevation=10)
    
    idField = DropdownID(page)
    
    graphMatPlot = MatplotlibChart(plt.figure(), expand=True)
    graph = chart.ProteinChart(page)

    minValue = ft.Text("Min count: \nPosition:", data=0, text_align=ft.TextAlign.LEFT)
    maxValue = ft.Text("Max count: \nPosition:", data=0, text_align=ft.TextAlign.LEFT)

    def changeBg(e):
        btn = e.control
        if btn.bgcolor != ft.Colors.GREY_300:
            btn.bgcolor = ft.Colors.GREY_300
            btn.update()
        else:
            btn.bgcolor = ft.Colors.GREY_200
            btn.update()

    unifesp_header = ft.AppBar(
        bgcolor="#215a36",
        title=ft.Button(
            content = ft.Row(
                [
                    ft.Image(
                        src="/UnifespLogo.png",
                        height=105,
                        filter_quality=ft.FilterQuality.HIGH,
                    ),
                    ft.Text(
                        "Universidade Federal de São Paulo",
                        font_family="Cinzel",
                        size=35,
                        color=ft.Colors.WHITE,
                        style=ft.TextThemeStyle.TITLE_LARGE,
                        text_align=ft.TextAlign.CENTER,
                    ),
                ],
                expand=1,
                alignment=ft.MainAxisAlignment.CENTER,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
            ),
            url="https://proteomics.bioinfo.unifesp.br/",
            bgcolor=ft.Colors.TRANSPARENT,
            color=ft.Colors.TRANSPARENT,
            style=ft.ButtonStyle(
                bgcolor=ft.colors.TRANSPARENT,
                overlay_color=ft.colors.TRANSPARENT,
                shape=None,
                padding=0,
                elevation=0,
                side=None,
            ),

        ),
        center_title=True,
        toolbar_height=120,
    )

    def route_change(route):
        page.views.clear()
        page.views.append(
            ft.View(
                "/",
                [
                    unifesp_header,
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
                                                    ft.Icon(name=ft.Icons.FILE_UPLOAD_ROUNDED, color=ft.Colors.GREEN_400, size=page.width*page.height*0.00008),
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
                                            width=page.width*0.25,
                                            height=page.height*0.4,
                                            bgcolor=ft.Colors.GREY_200,
                                            shape=ft.BoxShape.CIRCLE,
                                            on_hover=lambda e: changeBg(e),
                                            on_click=lambda _: file_picker.pick_files(allowed_extensions=["csv"])
                                        ),
                                    ],
                                    alignment=ft.MainAxisAlignment.CENTER,
                                    horizontal_alignment=ft.CrossAxisAlignment.CENTER
                                )
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    ft.CupertinoButton(text="Submit", bgcolor=ft.Colors.GREEN_700, color=ft.Colors.WHITE, on_click=lambda e: ViewResults(e), alignment=ft.alignment.bottom_center)
                ],
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                vertical_alignment=ft.MainAxisAlignment.CENTER
            ),
        )
        if page.route == "/results":
            try:
                idField.df = pd.read_csv(f'{arquivo.value}')
                if not set(['Protein Accession','Start', 'End', 'Peptide']).issubset(idField.df.columns):
                    fileWrong = ft.AlertDialog(False, title=ft.Text("Error"),
                                                content=ft.Text("File is incompatible."))
                    page.open(fileWrong)
                    page.go("/")
                    return
            except  Exception as e:
                print(arquivo.value)
                fileWrong = ft.AlertDialog(modal=False, title=ft.Text("Error"),
                                           content=ft.Text(f"Error loading the file.\n{e}"), disabled=False)
                page.open(fileWrong)
                page.go("/")
                return

            page.views.append(
                ft.View(
                    "/results",
                    [
                        ft.AppBar(title=ft.Text("Results"), bgcolor=ft.Colors.GREEN_700, color=ft.Colors.GREY_200, center_title=True),
                        ft.Row([idField, submitBtn], alignment=ft.MainAxisAlignment.CENTER),
                        ft.Tabs(
                            selected_index=0,
                            animation_duration=300,
                            clip_behavior=ft.ClipBehavior.NONE,
                            scrollable=False,
                            height=page.height*0.7,
                            tabs=[
                                ft.Tab(
                                    text="Interactive Graph",
                                    content=ft.Pagelet(
                                        content=ft.Card(
                                            content=ft.Stack([
                                                        ft.Container(
                                                            content=ft.Row([ft.InteractiveViewer(
                                                                    content=graph,
                                                                    scale_enabled=False,
                                                                    pan_enabled=False,
                                                                    width=graph.width*0.8,
                                                                    height=page.height*0.6
                                                                )], width=page.width*0.8, height=page.height*0.6, scroll=ft.ScrollMode.ALWAYS),
                                                        ),
                                                        ft.Row(
                                                            [
                                                               ft.Card(
                                                                    content=ft.Column([
                                                                        ft.Row([
                                                                            ft.Text("CDR1"), ft.Icon(ft.Icons.SQUARE_ROUNDED, color=ft.Colors.BLUE)
                                                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                                                        ft.Row([
                                                                            ft.Text("CDR2"), ft.Icon(ft.Icons.SQUARE_ROUNDED, color=ft.Colors.RED)
                                                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                                                        ft.Row([
                                                                            ft.Text("CDR3"), ft.Icon(ft.Icons.SQUARE_ROUNDED, color=ft.Colors.GREEN)
                                                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                                                    ]), width=100, height=115
                                                                ) 
                                                            ], alignment=ft.MainAxisAlignment.END
                                                        ),
                                                        ft.Row(
                                                            [
                                                               ft.Card(
                                                                    content=ft.Column([
                                                                        ft.Row([
                                                                            maxValue
                                                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                                                        ft.Divider(height=4, color=ft.Colors.BLACK),
                                                                        ft.Row([
                                                                            minValue
                                                                        ], alignment=ft.MainAxisAlignment.CENTER),
                                                                    ]), width=150, height=115
                                                                ) 
                                                            ], alignment=ft.MainAxisAlignment.START
                                                        ),
                                                    ], width=page.width, alignment=ft.alignment.top_center),
                                        ),
                                    )
                                ),
                                ft.Tab(
                                    text="Export Graph",
                                    content=ft.Column(
                                        [
                                            ft.Card(
                                                content=ft.Stack(
                                                    [
                                                        graphMatPlot,
                                                        ft.Row(
                                                            [
                                                                ft.IconButton(icon=ft.Icons.DOWNLOAD_ROUNDED, icon_color=ft.Colors.GREEN_700, padding=20,
                                                                            bgcolor=ft.Colors.WHITE12, tooltip="Export Image",
                                                                            on_click=lambda e: exportGraph(e))
                                                            ], alignment=ft.MainAxisAlignment.END, offset=ft.Offset(-0.01, 0.01)
                                                        )
                                                    ], alignment=ft.alignment.bottom_right
                                                ),
                                                width=page.width*0.8,
                                            ),
                                        ], horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                                        expand=1,
                                    ),
                                )
                            ],
                            expand=1,
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

    def deleteTemp(e):
        for i in temp_files:
            os.remove(i)

        for file in os.listdir(os.path.join("Counter", "src", "uploads")):
            if not os.path.isdir(file) and ".csv" in file:
                os.remove(os.path.join("Counter", "src", "uploads", file))
        
    page.on_close = deleteTemp
    page.on_disconnect = deleteTemp
    page.on_route_change = route_change
    page.on_view_pop = view_pop
    page.on_resized = updateView
    page.go(page.route)


ft.app(
    target=main,
    assets_dir="assets",
    upload_dir="uploads",
    view=None,
    port=8000,
    host="0.0.0.0",
    name ='tool',
    route_url_strategy="path"
)