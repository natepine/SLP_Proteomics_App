import pandas as pd
import numpy as np
import faicons as fa
import plotly.express as px
import plotly.subplots as sp
import plotly.graph_objects as go
import statsmodels.api as smf 

from shiny.types import FileInfo
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_plotly
from sklearn.decomposition import PCA

import requests ## python -m pip install requests
import io

import shinyswatch


def create_plot_list(input_df=None,protein_of_interest=None):
  temp_data_df = input_df.copy()

  timepoints = [0,0,0,0,0,1,1,2,2,4,4,6,6,8,8,8,8,8]

  results_col = ["Nonlinear_R2", "Nonlinear_Asym", "NonLinear_R0", "Linear_R0",
                "Linear_R2", "useLinear", "half_life", "R2", "cv_allTP", "short",
                "Nonlinear_AIC", "Linear_AIC", "half_life_CI_lower", "half_life_CI_upper",'lrc']

  ratio_columns = temp_data_df.filter(regex='scaled_ratio').columns

  # protein_of_interest = 'P25106'

  half_life_df = temp_data_df.copy().sort_values(by=['Protein Id']).set_index('Protein Id|Cell_line')

  half_life_df = half_life_df.loc[half_life_df['Protein Id'].str.contains(protein_of_interest),]


#   protein_data_df = protein_data_df.filter(regex='scaled_ratio')

  list_of_plots = []

  # print(df_len)

  for protein in half_life_df.index.to_list():
    uniprot_id = protein.split('|')[1]
    cell_id = protein.split('_')[-1]
    protein_id = protein.split('|')[2]
    gene_id = protein_id.split('_')[0]
    acession_id = f'{uniprot_id}|{gene_id}'

    sn_val = half_life_df.loc[protein,ratio_columns].values.astype(float)

    resuls_dict = half_life_df.loc[protein,results_col].to_dict()

    # half_life_result = half_life_v4(timepoints, sn_val, R2_cut=0.8,tolerance=1e-6)

    HL_plot = plot_half_life_plotly(timepoints, sn_val, resuls_dict, protein_id=protein,show_all_points=False)

    HL_plot.update_layout(title = dict(text = f"{cell_id}", font = dict(size = 20)))

    list_of_plots.append(HL_plot)
    # HL_plot.show()
    continue
  return(list_of_plots,acession_id)

def combine_plots(plot_list=None,gene_id=None):
  """
  Combines multiple Plotly figures into a single plot with subplots.
  Adapts the layout based on the length of the list, with a maximum of 3 columns.

  Args:
      plot_list (list): A list of Plotly figure objects.

  Returns:
      plotly.graph_objects.Figure: The combined Plotly figure.
  """

  num_plots = len(plot_list)
  ncols = min(3, num_plots)  # Maximum 3 columns
  nrows = (num_plots + ncols - 1) // ncols  # Calculate number of rows

  # Create subplots
  fig = sp.make_subplots(
      rows=nrows,
      cols=ncols,
      subplot_titles=[fig_data.layout.title.text for fig_data in plot_list],
      vertical_spacing=0.05,  # Adjust vertical spacing
      horizontal_spacing=0.1,  # Adjust horizontal spacing
  )

  # Add each figure to a subplot
  for i, fig_data in enumerate(plot_list):
      row = i // ncols + 1
      col = i % ncols + 1
      for trace in fig_data.data:
          fig.add_trace(trace, row=row, col=col)

      # Add annotations from the original figure
      for annotation in fig_data.layout.annotations:
          fig.add_annotation(annotation, row=row, col=col)




  fig.update_xaxes(range=[0, 8])  # Replace xmin and xmax with your desired limits
  fig.update_yaxes(range=[0, 1])  # Replace ymin and ymax with your desired limits


  
  fig.update_layout(plot_bgcolor = 'rgba(0,0,0,0)',
                  font = dict(color = "#909497", size = 18),
                  title = dict(text = f"{gene_id}", font = dict(size = 20),x= 0.5,xanchor= 'center'),
                  margin = dict(t = 100, r = 80, b = 50, l = 120),
                #   height = 500*nrows,
                  # width = 1000,
                  showlegend = False,)
                  #  paper_bgcolor='rgba(0,0,0,0)')

  fig.update_xaxes(linecolor = "gray", range = [-0.01, 8.5],linewidth=2,showgrid=False)
  fig.update_yaxes(range=[0.01, 1.1], linecolor = "gray",linewidth=2,showgrid=False)  

  fig.add_annotation(text = "Timepoint (hr)",
                      xref = "paper",
                      yref = "paper",
                      x = 0.5,
                      y = -0.2,
                      showarrow = False)

  fig.add_annotation(text = "Scaled SN",
                      xref = "paper",
                      yref = "paper",
                      x = -0.08,
                      y = 0.5,
                      showarrow = False,
                      textangle = -90)

  # fig["layout"]["annotations"][0]["font"]["size"] = 20
  # fig["layout"]["annotations"][1]["font"]["size"] = 20
  # fig["layout"]["annotations"][2]["font"]["size"] = 20

  # fig.show()
  return(fig)

def plot_half_life_plotly(tp, sn, half_life_result, protein_id, show_all_points=False):
    """
    Plots the observed data and the fitted curve (linear or nonlinear)
    with dashed lines indicating the half-life using Plotly.  Includes
    confidence interval visualization.

    Args:
        tp (list or numpy.ndarray): Time points.
        sn (list or numpy.ndarray): Signal intensities.
        half_life_result (dict): Output from the half_life_v4 function.
        protein_id (str): Protein ID for the plot title.
        show_all_points (bool): Show all points or just averages.

    Returns:
        plotly.graph_objects.Figure: The Plotly figure.
    """

    # Consistent tp_dense
    tp_dense = np.linspace(min(tp), max(tp), 300)

    # Prepare data
    if show_all_points:
        plot_tp = tp
        plot_sn = sn
        mode = 'markers'  # Use markers for individual points
        marker_size = 7
        alpha = 0.75
    else:
        tpsn = pd.DataFrame({"tp": tp, "sn": sn})
        agg_tpsn = tpsn.groupby("tp")["sn"].mean().reset_index()
        plot_tp = agg_tpsn["tp"]
        plot_sn = agg_tpsn["sn"]
        mode = 'markers'
        marker_size = 10
        alpha = 1

        # Check for empty plot_sn
    if not isinstance(plot_sn, np.ndarray):
        plot_sn = np.array(plot_sn) #Convert to numpy array for .size
    if not plot_sn.size:
        print("Warning: No data to plot after filtering.")
        return go.Figure()  # Return an empty figure

    # Create the base figure
    fig = go.Figure()

    # Add observed data points
    fig.add_trace(go.Scatter(x=plot_tp, y=plot_sn, mode=mode, name="Observed Data",
                             marker=dict(size=marker_size, color='steelblue', opacity=alpha)))

    if half_life_result["useLinear"] == 'true':
        # Linear fit
        slope = np.log(0.5) / half_life_result["half_life"]
        intercept = np.log(half_life_result["Linear_R0"])
        y_fit = np.exp(intercept + slope * tp_dense)
        print(np.isnan(y_fit).sum())

        fig.add_trace(go.Scatter(x=tp_dense, y=y_fit, mode='lines', name="Linear Fit",
                                 line=dict(color='black')))
        half_life_intensity = np.exp(intercept + slope * half_life_result["half_life"])


        # --- Linear Confidence Interval (statsmodels) ---
        df = pd.DataFrame({'tp': tp_dense, 'sn': y_fit})
        lmodel = smf.ols('np.log(sn) ~ tp', data=df).fit()
        predictions = lmodel.get_prediction(df)
        pred_df = predictions.summary_frame(alpha=0.05)

        # fig.add_trace(go.Scatter(x=tp_dense, y=np.exp(pred_df['obs_ci_upper']), mode='lines',
        #                          line=dict(width=0), showlegend=False, name='CI_upper'))
        # fig.add_trace(go.Scatter(x=tp_dense, y=np.exp(pred_df['obs_ci_lower']), mode='lines',
        #                          line=dict(width=0), fill='tonexty', fillcolor='rgba(128,128,128,0.3)',
        #                          showlegend=False, name='CI_lower'))
        # --- End Linear CI ---

    else:
        # Nonlinear fit
        def ssasymp(t, Asym, R0, lrc):
            return Asym + (R0 - Asym) * np.exp(-np.exp(lrc) * t)

        y_fit = ssasymp(tp_dense, half_life_result["Nonlinear_Asym"],
                        half_life_result["NonLinear_R0"], half_life_result["lrc"])
        # print(np.isnan(y_fit).sum())

        if np.isnan(y_fit).sum() == 0:
            

            fig.add_trace(go.Scatter(x=tp_dense, y=y_fit, mode='lines', name="Nonlinear Fit",
                                    line=dict(color='black')))
            half_life_intensity = ssasymp(half_life_result["half_life"], half_life_result["Nonlinear_Asym"],
                                        half_life_result["NonLinear_R0"], half_life_result["lrc"])
        else:
            half_life_intensity = ssasymp(half_life_result["half_life"], half_life_result["Nonlinear_Asym"],
                                        half_life_result["NonLinear_R0"], half_life_result["lrc"])
        # --- Nonlinear Confidence Interval ---
        # if not np.isnan(half_life_result["half_life_CI_lower"]) and not np.isnan(half_life_result["half_life_CI_upper"]):
        #   try:
            # y_lower = ssasymp(half_life_result["half_life_CI_lower"], half_life_result["Nonlinear_Asym"], half_life_result["NonLinear_R0"], half_life_result["lrc"])
            # y_upper = ssasymp(half_life_result["half_life_CI_upper"], half_life_result["Nonlinear_Asym"], half_life_result["NonLinear_R0"], half_life_result["lrc"])

            # #Extend lines over the full range
            # x_lower = [0, half_life_result["half_life_CI_lower"]]
            # y_lower_extend = [y_lower, y_lower]
            # x_upper = [0, half_life_result["half_life_CI_upper"]]
            # y_upper_extend = [y_upper, y_upper]


            # # Create a dense array that can be used
            # dense_x = np.linspace(min(tp_dense), max(tp_dense), 300)
            # y_lower_interp = np.interp(dense_x, x_lower, y_lower_extend, left=np.nan, right=np.nan)
            # y_upper_interp = np.interp(dense_x, x_upper, y_upper_extend, left=np.nan, right=np.nan)

            # # Remove nans so that the plot draws correctly.
            # nan_indices = np.isnan(y_lower_interp) | np.isnan(y_upper_interp)
            # y_lower_interp = y_lower_interp[~nan_indices]
            # y_upper_interp = y_upper_interp[~nan_indices]
            # dense_x_filtered = dense_x[~nan_indices]


            # fig.add_trace(go.Scatter(x=dense_x_filtered, y=y_upper_interp,
            #                             mode='lines', line=dict(width=0),
            #                             showlegend=False, name='CI_upper'))
            # fig.add_trace(go.Scatter(x=dense_x_filtered, y=y_lower_interp,
            #                             mode='lines', line=dict(width=0),
            #                             fill='tonexty', fillcolor='rgba(128,128,128,0.3)',
            #                             showlegend=False, name='CI_lower'))

        #   except:
        #       pass
        # --- End Nonlinear CI ---


    # Add half-life lines (handle inf/nan)
    half_life = half_life_result["half_life"]
    
    if not np.isinf(half_life) and not np.isnan(half_life) and not np.isnan(half_life_intensity):
        # print(half_life,half_life_intensity)
        fig.add_trace(go.Scatter(x=[half_life, half_life], y=[0, half_life_intensity],
                                 mode='lines', line=dict(dash='dash', color='gray'),
                                 name=f"Half-life: {half_life:.2f}"))
        fig.add_trace(go.Scatter(x=[0, half_life], y=[half_life_intensity, half_life_intensity],
                                 mode='lines', line=dict(dash='dash', color='gray'),
                                 showlegend=False))  # Hide duplicate legend entry



    # --- Dynamic Axis Limits ---
    if not isinstance(y_fit, np.ndarray): # Convert y_fit to an array so the max function works
        y_fit = np.array(y_fit) # Convert list to array
    max_y = max(np.max(plot_sn), np.max(y_fit))
    fig.update_layout(xaxis_title="Time Points",
                      yaxis_title="Relative Protein Abundance",
                      title=f"Protein: {protein_id}",
                      yaxis_range=[0, max_y * 1.1],
                      xaxis_range=[0, max(tp) * 1.1],
                      legend=dict(
                        orientation="v",
                        yanchor="top",
                        y=0.99,
                        xanchor="right",
                        x=0.99))
    # ---
    if not np.isinf(half_life) and not np.isnan(half_life) and not np.isnan(half_life_intensity):
        fig.add_annotation(text=f'HL: {half_life:.2f}',
                            x=half_life, 
                            y=half_life_intensity,
                            # ax=50, 
                            # ay=-10,
                            xshift=3,
                            xanchor='left',
                            arrowsize=3,
                            font_size=12)
    else:
        fig.add_annotation(text=f'HL: {half_life:.2f}',
                            x=0.2, 
                            y=4,
                            # ax=50, 
                            # ay=-10,
                            # xshift=3,
                            xanchor='left',
                            arrowsize=3,
                            font_size=14)


    return fig

def get_string_enrichment(input_gene_list, input_bkg_gene_list=[]):

  string_api_url = "https://version-12-0.string-db.org/api"
  output_format = "tsv"
  method = "enrichment"

  ##
  ## Construct the request
  ##

  request_url = "/".join([string_api_url, output_format, method])

  ##
  ## Set parameters
  ##

  my_genes = input_gene_list
  my_bkg_genes = input_bkg_gene_list


  params = {
      "identifiers" : "%0d".join(my_genes), # your protein
      'background_string_identifiers' : "%0d".join(my_bkg_genes), # your background
      "species" : 9606, # NCBI/STRING taxon identifier 
      "caller_identity" : "www.awesome_app.org" # your app name
  }

  ##
  ## Call STRING
  ##

  response = requests.post(request_url, data=params)

  ##
  ## Read and parse the results
  ##

  df = pd.read_csv(io.StringIO(response.text), sep='\t')
  return(df)

def convert_to_stringDB_df(input_accession_list=None):
  string_api_url = "https://version-12-0.string-db.org/api"
  output_format = "tsv"
  method = "get_string_ids"

  ##
  ## Set parameters
  ##

  params = {

      "identifiers" : "\r".join(input_accession_list), # your protein list
      "species" : 9606, # NCBI/STRING taxon identifier 
      "limit" : 1, # only one (best) identifier per input protein
      "echo_query" : 1, # see your input identifiers in the output
      "caller_identity" : "www.awesome_app.org" # your app name

  }

  ##
  ## Construct URL
  ##


  request_url = "/".join([string_api_url, output_format, method])

  ##
  ## Call STRING
  ##

  results = requests.post(request_url, data=params)

  ##
  ## Read and parse the results
  ##
  # string_id =[]
  # for line in results.text.strip().split("\n"):
  #     l = line.split("\t")
  #     input_identifier, string_identifier = l[0], l[2]
  #     string_id.append(string_identifier)
  #     print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
  # string_id
  results.text
  df = pd.read_csv(io.StringIO(results.text), sep='\t')
  stringID_list = df['stringId'].tolist()
  return(stringID_list)

def PCA_data(input_df,PC_x=None,PC_y=None,cell_list=None):
  pivot_data_df = input_df.copy()

  pivot_data_df = pivot_data_df.set_index('Protein Id')

  pivot_data_df = pivot_data_df.filter(regex='scaled_ratio|^Cell_line$')

  list_of_cells = pivot_data_df['Cell_line'].unique().tolist()

  pivoted_df = pd.DataFrame()

  for cell in pivot_data_df['Cell_line'].unique():
      temp = pivot_data_df.loc[pivot_data_df['Cell_line'] == cell,]
      temp = temp.drop(columns=['Cell_line'])

      for column in temp.columns:
          temp = temp.rename(columns={column:column + '|' + cell})

      pivoted_df = pd.concat([pivoted_df,temp],axis=1)

  rgx_str = ''

  for cell in list_of_cells:

      rgx_str = rgx_str + cell[0:3] + '|'

  rgx_str = rgx_str.rstrip('|') 

  pivoted_df = pivoted_df.filter(regex=rgx_str)

  pivoted_df = pivoted_df.dropna()



  pca_out = PCA().fit(pivoted_df)

  np.cumsum(pca_out.explained_variance_ratio_)

  loadings = pca_out.components_

  num_pc = pca_out.n_features_in_

  pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
  loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
  loadings_df['variable'] = pivoted_df.columns.values
  loadings_df = loadings_df.set_index('variable')

  pca_out.explained_variance_

  x = PC_x
  y = PC_y

  pc_len = len(pivoted_df.columns.to_list())

  variance_ratio_dict = {}
  for i in range(pc_len):
    variance_ratio_dict[f'PC{i + 1}'] = round(pca_out.explained_variance_ratio_[i]*100, 2)


  x_var = variance_ratio_dict[PC_x]
  y_var = variance_ratio_dict[PC_y]


  TP_list = ['0h','1h','2h','4h','6h','8h']

  cust_pal = ['#003f5c',
                '#3a5c75',
                '#637b8e',
                '#8c9ba8',
                '#b5bcc3',
                '#dfdfdf']

  loadings_df['timepoint'] = None

  for timepoint in TP_list:
    loadings_df.loc[loadings_df.index.str.contains(timepoint),'timepoint'] = timepoint

  loadings_df['cell_line'] = None

  for cell_line in cell_list:
    loadings_df.loc[loadings_df.index.str.contains(cell_line),'cell_line'] = cell_line


  return(loadings_df,x_var,y_var)

def TMT_PCA_out(input_df,PC_x=None,PC_y=None,Xvariance=None,Yvariance=None):

  loadings_df = input_df.copy()

  marker_symbols  = {'0h':'circle',
                    '1h':'square',
                    '2h':'x',
                    '4h':'star-square',
                    '6h':'diamond',
                    '8h':'hexagon2'}

  fig = px.scatter(loadings_df.filter(regex='PC.'),
                    x=PC_x,
                    y=PC_y,
                    labels={PC_x : f'{PC_x} ({Xvariance}%)',
                            PC_y : f'{PC_y} ({Yvariance}%)' },
                    color=loadings_df['cell_line'].values,
                    symbol=loadings_df['timepoint'].values,
                    color_discrete_sequence=px.colors.qualitative.Prism,
                    symbol_map=marker_symbols)

  fig.update_traces(marker={'size': 15})



  return(fig)



ratio_df = pd.read_csv('data/Ratios_data.csv')
hl_df = pd.read_csv('data/HL_data.csv')
uniprot_df = pd.read_csv('data/idmapping_2025_03_07.tsv',sep='\t').rename(columns={'Entry':'Accession'})

data_df = pd.merge(ratio_df,hl_df,on='Protein Id|Cell_line')
data_df['Accession'] = data_df['Protein Id'].str.split('|').str[1]  
data_df = pd.merge(data_df,uniprot_df,on='Accession',how='left')
data_df = data_df.drop_duplicates(subset=['Protein Id|Cell_line'])
data_df = data_df.sort_values(by='Cell_line')

# data_df['Protein Id'] = str(data_df['Protein Id'] + "|" + data_df['Gene Symbol'])

proteins = sorted(data_df['Protein Id'].unique().tolist()) 
cells = sorted(data_df['Cell_line'].unique().tolist())  # Dynamically get cells and sort them

DSLP_df = pd.read_csv('data/DSLP_analysis_w_uniprot.csv')

PCA_loadings_df,x_var,y_var  = PCA_data(data_df,PC_x='PC1',PC_y='PC2',cell_list=cells)

ICONS = {
    "user": fa.icon_svg('vials'),
    "wallet": fa.icon_svg("list-check"),
    "currency-dollar": fa.icon_svg("stopwatch"),
    "ellipsis": fa.icon_svg("ellipsis"),
}

# Add page title and sidebar
app_ui = ui.page_sidebar(
            ui.sidebar(
                ui.page_fluid(
                    ui.input_checkbox_group("Cell_line",
                                        "Cell Line",
                                        cells,
                                        selected=cells,
                                        inline=True,
                                        width='150px'), 
                    ui.input_action_button("Select_all", "Select All"),
                    ui.input_action_button("Deselect_all", "Deselect All"),
                    ui.hr(),
                    ui.input_selectize("protein_id", "Select Protein", multiple = False, choices=['sp|Q9BXS6|NUSAP_HUMAN'],selected='sp|Q9BXS6|NUSAP_HUMAN',width='800px'),
                    open="desktop"),
                ui.input_dark_mode(mode='light'),
                width = 335),
            ui.page_fluid(
                ui.layout_columns(
                    ui.value_box("Cell Lines", ui.output_ui("cell_line_count"), showcase=ui.img(width="80", 
                     height="80", 
                     src="https://img.icons8.com/arcade/100/body-cells.png",
                     alt="body-cells")),
                    ui.value_box("Unique Proteins", ui.output_ui("protein_count"), showcase=ui.img(width="80", 
                     height="80", 
                     src="https://img.icons8.com/arcade/100/protein.png",
                     alt="protein")),
                    ui.value_box("Short Lived Proteins",ui.output_ui("slp_count"), showcase=ui.img(width="80", 
                     height="80", 
                     src="https://img.icons8.com/arcade/100/timer.png",
                     alt="timer")),
                    fill=False),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Short Lived Proteins"), 
                        ui.output_data_frame("full_slp_table"), 
                        ui.download_button("download_all_SLP_data", "Download CSV"),
                        full_screen=True),
                    ui.card(
                        ui.card_header(
                                "PCA",
                                ui.popover(
                                    ICONS["ellipsis"],
                                    ui.input_radio_buttons(
                                        "cell_for_pca",
                                        'Combined',
                                        cells + ['Combined'],
                                        inline=True,
                                    ),
                                    title="Add a color variable",
                                    placement="top",
                                ),
                                class_="d-flex justify-content-between align-items-center",
                            ),
                        output_widget("pca_plt"),
                        full_screen=True,),
                    ui.card(
                        ui.card_header("SLP Overlap",class_="d-flex justify-content-between align-items-center",),
                        output_widget("slp_overlap_heatmap"),
                        full_screen=True,),
                    ui.card(
                        ui.card_header("Protein Count",class_="d-flex justify-content-between align-items-center",),
                        output_widget("cell_protein_count"),
                        full_screen=True,),
                    col_widths=[12,4,4,4],),
                ui.hr(),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Protein Function"),
                        ui.output_ui("protein_function_text")),
                    ui.card(
                        ui.card_header("Protein SN Values",class_="d-flex justify-content-between align-items-center",),
                        output_widget("SN_ratio_plot"),
                        full_screen=True,),
                    
                    col_widths=[6,6],),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Protein Half-life Plots",class_="d-flex justify-content-between align-items-center",),
                        output_widget("Half_life_plots"),
                        fillable=False,
                        full_screen=True,),),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Half-life Data"), 
                        ui.output_data_frame("HL_table"), 
                        ui.download_button("download_protein_HL_data", "Download CSV"),
                        full_screen=True),),
                ui.hr(),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Short-lived Proteins in All Cell Lines"), 
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.input_selectize(
                                    "Select_cells_w_SLP", "Select Cell Lines:", dict(zip(cells, cells)),multiple=True,),open="desktop"),
                            ui.output_data_frame("create_only_slp_df")), 
                            ui.download_button("download_SLP_in_all_cells", "Download CSV"),
                        full_screen=True),
                        ),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Gene Enrichment for SLP in Selected Cell Lines"), 
                        ui.output_data_frame("gene_enrichment_df"), 
                        ui.download_button("download_gene_enrichment", "Download CSV"),
                        full_screen=True),
                        ),
                ui.hr(),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Short-lived Proteins With Differential Stability"), 
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.input_selectize(
                                    "Select_cells_w_diffSLP", "Select Cell Lines:", dict(zip(cells, cells)),multiple=True,),open="desktop"),
                            ui.output_data_frame("differential_slp")), 
                            ui.download_button("download_DSLP_data", "Download CSV"),
                        full_screen=True),
                        ),
                ),

            
            title="Short Lived Protein Proteomics",
            fillable=True,
            # theme=shinyswatch.theme.darkly,
)


def server(input, output, session):

    @reactive.calc
    def data_df_1():
        return data_df[data_df.Cell_line.isin(input.Cell_line())]
    
    @reactive.Effect
    @reactive.event(input.Cell_line)
    def update_protein_selectize():
        df = data_df_1()
        protein_list = list(set(df['Protein Id'].values.tolist()))
        protein_dict = dict(zip(protein_list, protein_list))
        ui.update_selectize("protein_id",choices=protein_list, selected='sp|Q13772|NCOA4_HUMAN')
        # return protein_dict

    # @render.text
    # def prot_search():
    #     return input.protein_id()

    @reactive.calc
    def filter_PCA_data():
        PCA_loadings_df
        ldngs_df = PCA_loadings_df.loc[PCA_loadings_df['cell_line'].isin(input.Cell_line()),]
        return(ldngs_df)
    
    @render_plotly
    def pca_plt():
        pca_data_df = filter_PCA_data()
        pca_plt = TMT_PCA_out(pca_data_df,PC_x='PC1',PC_y='PC2',Xvariance=x_var,Yvariance=y_var)
        return(pca_plt)

    @render.ui
    def cell_line_count():
        return str(data_df_1()['Cell_line'].nunique())
    
    @render.ui
    def protein_count():
        return str(data_df_1()['Protein Id'].nunique())

    @render.ui
    def slp_count():
        return str(data_df_1().loc[data_df_1()['short'] == True, 'Protein Id'].nunique())
    
    @render.ui
    def protein_function_text():
        df = data_df_1()

        df = df.loc[df['Protein Id']==input.protein_id(),['Protein Id','Function [CC]']]

        if len(df) == 0:
            return('Protein Not In Data. Please search another protein or check spelling')

        function = str(df['Function [CC]'].unique()[0]).replace('FUNCTION: ','')
        protein_ID = str(df['Protein Id'].unique()[0])
        protein = str(df['Protein Id'].unique()[0]).split('|')[1]
        uniprot_url = f'https://www.uniprot.org/uniprotkb/{protein}/entry'

        return ui.HTML(f"<a href='{uniprot_url}' target='_blank'>{protein_ID}</a> {function}")
        

    @render.data_frame
    def full_slp_table():
        df = data_df_1()

        if len(df) == 0:
            return render.DataTable(pd.DataFrame())

        df = df.loc[df['short']== True,
                                 ['Protein Id','Gene Symbol','Number of peptides', 'Cell_line', 'Description', 'half_life', 'half_life_CI_lower',
                                  'half_life_CI_upper']]
        df[['half_life','half_life_CI_lower','half_life_CI_upper']]  = round(df[['half_life','half_life_CI_lower','half_life_CI_upper']],2)
        df["Description"] = df["Description"].str.split(" OS=Homo", expand=True)[0]
        return render.DataTable(df, 
                                filters=True, 
                                styles ={'class' : "display text-center"},
                                width="100%")
    
    @render.download(filename="all_SLP_data.csv")
    def download_all_SLP_data():
            SLP_data = full_slp_table.data_view(selected=False) 
            yield SLP_data.to_csv(index=False)

    
    
    @render_plotly
    def slp_overlap_heatmap():
        df = data_df_1()

        if len(df) == 0:
            fig = go.Figure()
            return(fig)


        df = df.loc[data_df_1()['short'] == True, ['Protein Id', 'Cell_line']]
        df = df.sort_values(by='Cell_line')
        cells = df['Cell_line'].unique()
        overlap_df = pd.DataFrame(index=cells, columns=cells)

        for cell1 in cells:
            for cell2 in cells:
                if cell1 == cell2:
                    overlap_df.loc[cell1, cell2] = df.loc[df['Cell_line'] == cell1, 'Protein Id'].nunique()
                else:
                    shared = df.loc[df['Cell_line'].isin([cell1, cell2]), 'Protein Id'].duplicated(keep=False).sum() // 2 # Correct overlap calculation
                    overlap_df.loc[cell1, cell2] = shared

        overlap_df = overlap_df.sort_index(axis=0)
        hplt = px.imshow(overlap_df, text_auto=True, labels=dict(x="Cell Line", y="Cell Line", color="Protein Overlap"))
        hplt.update_layout(coloraxis_showscale=False)
        return hplt

    @render_plotly
    def cell_protein_count():

        counts_df = data_df_1()

        counts = pd.DataFrame(counts_df.groupby('Cell_line')['Protein Id'].nunique())



        if len(counts) == 0:
            fig = go.Figure()
            return(fig)

        return px.bar(counts, 
                      x=counts.index, 
                      y='Protein Id', 
                      color=counts.index,labels={'x': 'Cell Line', 'y': 'Protein Count'},
                      color_discrete_sequence=px.colors.qualitative.Prism) # More clear labels

    @render_plotly
    def SN_ratio_plot():
        data_df = data_df_1()

        protein_data_df = data_df.loc[data_df['Protein Id']==input.protein_id(),]

        if len(protein_data_df) == 0:
            fig = go.Figure()
            return(fig)

        protein_data_df = protein_data_df.filter(regex='scaled_ratio|^Cell_line$|^Protein.Id$')

        protein_name = protein_data_df['Protein Id'].unique()[0]

        melt_protein_data_df = pd.melt(protein_data_df,id_vars=['Protein Id','Cell_line'],var_name='Time',value_name='Ratio')
        melt_protein_data_df['Time'] = melt_protein_data_df['Time'].str[0]
        melt_protein_data_df['Time'] = melt_protein_data_df['Time'].astype(int)
        melt_protein_data_df
        mean_protein_data_df = melt_protein_data_df.groupby(['Time','Cell_line'])['Ratio'].mean().reset_index()

        fig = px.line(mean_protein_data_df,
                      x='Time',
                      y='Ratio',
                      color='Cell_line',
                      markers=True,
                      title=protein_name,
                      color_discrete_sequence=px.colors.qualitative.Prism)
        
        return fig
    
    @render_plotly
    def Half_life_plots():
        data_df = data_df_1()

        # if data_df['Protein Id|Cell_line'].str.contains(input.protein_id()) == False:
        #     return 

        sum_ratio_columns = data_df.filter(regex='sum')
        data_df = data_df.drop(columns=sum_ratio_columns)
        data_df = data_df.loc[data_df['Protein Id']==input.protein_id(),]

        if len(data_df) == 0:
            fig = go.Figure()
            return(fig)


        # data_df = data_df.loc[data_df['short'] ==True,]
        # data_df = data_df.loc[data_df['R2'] >=0.8,]
        

        plt_list,accsn_id = create_plot_list(input_df=data_df,protein_of_interest=input.protein_id())

        # if len(plt_list) == 0:
        #     fig = go.Figure()
        #     return(fig)

        fig = combine_plots(plot_list = plt_list,gene_id=accsn_id)

        return fig

    @render.data_frame
    def HL_table():
        data_df = data_df_1()

        data_df = data_df.loc[data_df['Protein Id']==input.protein_id(),
                                 [ 'Cell_line','Number of peptides', 'Description', 'half_life', 'half_life_CI_lower',
                                  'half_life_CI_upper']]
        
        # df = data_df_1().loc[data_df_1()['Protein Id'].str.contains(input.protein_id()),
        #                     [ 'Cell_line','Number of peptides', 'Description', 'half_life', 'half_life_CI_lower',
        #                     'half_life_CI_upper']]
        
        if len(data_df) == 0:
            return(render.DataGrid(pd.DataFrame()))
        
        data_df[['half_life','half_life_CI_lower','half_life_CI_upper']]  = round(data_df[['half_life','half_life_CI_lower','half_life_CI_upper']],2)

        data_df["Description"] = data_df["Description"].str.split(" OS=Homo", expand=True)[0]
        return render.DataTable(data_df, 
                                filters=True, 
                                styles ={'class' : "display text-center"},
                                width="100%")

    @render.download(filename="selected_protein_half_life_data.csv")
    def download_protein_HL_data():
            HL_selected_data = HL_table.data_view(selected=False) 
            yield HL_selected_data.to_csv(index=False)



### SLP in all cells 

    @reactive.calc
    def filter_slp_only_df():
        cell_list_of_interest = input.Select_cells_w_SLP()
        list_count = len(input.Select_cells_w_SLP())

        if list_count == 0:
            return (pd.DataFrame(),None,None)


        SLP_df = data_df.copy()
        SLP_df = SLP_df.loc[SLP_df['Cell_line'].isin(cell_list_of_interest)]

        backgroun_list = list(set(SLP_df['Protein Id'].str.split('|').str[1].to_list()))

        SLP_df = SLP_df.loc[SLP_df['short'] == True,]
        SLP_df['SLP_Occurance'] = SLP_df.groupby('Protein Id')['short'].transform('sum')
        SLP_df = SLP_df.loc[SLP_df['SLP_Occurance'] >= list_count,]

        SLP_df = SLP_df[['Protein Id',
                         'Gene Symbol',
                         'Cell_line',
                         'half_life',
                         'R2',
                         'half_life_CI_lower',
                         'half_life_CI_upper',
                         'Description',
                         'Function [CC]',
                         'Gene Ontology (GO)']]
        SLP_df["Description"] = SLP_df["Description"].str.split("OS=Homo", expand=True)[0]

        protein_list = list(set(SLP_df['Protein Id'].str.split('|').str[1].to_list()))

        return(SLP_df,backgroun_list,protein_list)
    
    @render.data_frame
    def create_only_slp_df():
        SLP_df = filter_slp_only_df()[0]

        if len(SLP_df) == 0:
            return render.DataTable(pd.DataFrame(columns=['Protein Id',
                         'Cell_line',
                         'half_life',
                         'R2',
                         'half_life_CI_lower',
                         'half_life_CI_upper',
                         'Description']))
        
        SLP_df[['half_life','half_life_CI_lower','half_life_CI_upper']]  = round(SLP_df[['half_life','half_life_CI_lower','half_life_CI_upper']],2)
        
        return render.DataTable(SLP_df, 
                                filters=True, 
                                styles ={'class' : "display text-center"},
                                width="100%")
    
    
    @render.download(filename="SLPs_in_all_cells.csv")
    def download_SLP_in_all_cells():
            selected_cell_data = create_only_slp_df.data_view(selected=False)
            # file: list[FileInfo] | None = data
            # df = pd.read_csv(file[0]["datapath"])
            yield selected_cell_data.to_csv(index=False)

    @render.data_frame
    def gene_enrichment_df():
        uniprot_bkg_lst,gene_lst = filter_slp_only_df()[1:3]

        if gene_lst == None:
            SDB_results_df = pd.DataFrame(columns=['category','description','number_of_genes','p_value','fdr','number_of_genes_in_background','preferredNames'])
            return render.DataTable(SDB_results_df, filters=True, styles ={'class' : "display text-center"})

        SDB_bk_lst = convert_to_stringDB_df(input_accession_list = uniprot_bkg_lst)

        SDB_results_df = get_string_enrichment(gene_lst, input_bkg_gene_list = SDB_bk_lst)

        if len(SDB_results_df) == 0:
            return(render.DataGrid(pd.DataFrame(columns=['category','description','number_of_genes','p_value','fdr','number_of_genes_in_background','preferredNames'])))

        SDB_results_df = SDB_results_df[['category','description','number_of_genes','p_value','fdr','number_of_genes_in_background','preferredNames']]

        return render.DataTable(SDB_results_df, 
                                filters=True, 
                                styles ={'class' : "display text-center"},
                                width="100%")
    
    @render.download(filename="Gene_enrichment_in_all_cells.csv")
    def download_gene_enrichment():
            GO_selected_data = gene_enrichment_df.data_view(selected=False) 
            yield GO_selected_data.to_csv(index=False)


### Differential SLP output 

    @reactive.calc
    def filter_DSLP():
        cell_list_of_interest = input.Select_cells_w_diffSLP()

        if len(cell_list_of_interest) == 0:
            return pd.DataFrame()
        
        SLP_df = DSLP_df.copy()
        SLP_df[['half_life1','half_life2','CI1_lower','CI1_upper','CI2_lower','CI2_upper']] = round(SLP_df[['half_life1','half_life2','CI1_lower','CI1_upper','CI2_lower','CI2_upper']],2)
        SLP_df = SLP_df.loc[(SLP_df['Cell_line1'].isin(cell_list_of_interest))&(SLP_df['Cell_line2'].isin(cell_list_of_interest))]
        SLP_df = SLP_df.loc[SLP_df['DSLP'] == True,]
        SLP_df = SLP_df.drop(columns=['DSLP'])
        return(SLP_df)
    
    @render.data_frame
    def differential_slp():
        filtered_DSLP_df = filter_DSLP()
        if len(filtered_DSLP_df) == 0:
            return render.DataTable(pd.DataFrame(columns=['Protein Id',
                                        'Cell_line1',
                                        'Cell_line2',
                                        'half_life1',
                                        'half_life2',
                                        'CI1_lower',
                                        'CI1_upper',
                                        'CI2_lower',
                                        'CI2_upper']),filters=True)
        
        return render.DataTable(filtered_DSLP_df, 
                                filters=True, 
                                styles ={'class' : "display text-center"},
                                width="100%")
    
    @render.download(filename="Differential_SLP_data.csv")
    def download_DSLP_data():
            # file: list[FileInfo] | None = data
            # df = pd.read_csv(file[0]["datapath"])
            yield filter_DSLP().to_csv(index=False)
    

    @reactive.effect
    @reactive.event(input.Select_all)
    def _():
        # ui.update_slider("total_bill", value=bill_rng)
        ui.update_checkbox_group("Cell_line", selected=cells)

    @reactive.effect
    @reactive.event(input.Deselect_all)
    def __():
        # ui.update_slider("total_bill", value=bill_rng)
        ui.update_checkbox_group("Cell_line", selected=[])


app = App(app_ui, server)