import streamlit as st

st.title("Create new 1D FDTD Simulation")

st.write("Add new simulation page")

sim_param = {}


sim_param["sim_desc"] = st.text_area(label="Description of the simulation:")
sim_param["isFileUpload"] = st.toggle(label="Toggle to select which mode to load input model/s (Manual Input/File Upload) ")
st.divider()
if sim_param["isFileUpload"]:
    #Create a file upload input 
    st.subheader("Mode: File Upload")
    sim_param["file_upload"] = st.file_uploader(label="Upload *.csv file for the input model/s",
                                                accept_multiple_files=False,
                                                type="csv")
    
    #Read the CSV file and save them to sim_param
    #TODO

else:
    #Create manual inputs
    st.subheader("Mode: Manual Input")
    col1,col2,col3 = st.columns(3)

    with col1:
        sim_param["freq"]= st.number_input("Frequency of interest (Hz):",min_value=0,max_value=None)
    with col2:
        sim_param["src_type"] = st.selectbox(label="Source Type:",options=("Gaussian Pulse Source", "Sinusoidal Source", "Rectangular Pulse Source", "Modulated Sine Pulse Source"))

    with col3:
        sim_param["num_device_layers"] = st.number_input("Number of Device Layers:", min_value=0)

    st.write(f"Number of layers: {sim_param["num_device_layers"]}")

    layer_col1,layer_col2,layer_col3 = st.columns(3)
    for layer in range(1,sim_param["num_device_layers"]+1):
        with layer_col1:
            sim_param[f"layer_{layer}_size"] = st.number_input("Layer Siza (m):", key=f"layer_{layer}_size", min_value=0.0, step=0.00000000000000000001,format="%0.20f")
        with layer_col2:
            sim_param[f"layer_{layer}_mu"] = st.number_input("Relative magnetic permeability (μ_r):", key=f"layer_{layer}_mu", min_value=0.0, step=0.001,format="%0.3f")
        with layer_col3:
            sim_param[f"layer_{layer}_epsilon"] = st.number_input("Relative electric permittivity (ε_r)", key=f"layer_{layer}_epsilon", min_value=0.0, step=0.001,format="%0.3f")
st.divider()
st.header("Simulation Parameters:")
param_col1, param_col2,  = st.columns(2)
with param_col1:
    sim_param["boundary_condition"] = st.selectbox(label="Boundary Condition",options=("Perfectly Absorbing Boundary Condition", "Dirichlet Boundary Condition"))
    sim_param["custom_name"] = st.text_input("Custom Output Filename")
    sim_param["algo"] = st.selectbox(label="Algorithm",options=("Finite Difference Time Domain", "Combined Finite Difference Time Domain and Schwarz's Alternating Method"))
with param_col2:
    sim_param["src_excitation"] = st.selectbox(label="Source Excitation Method",options=("Hard Source Method", "Soft Source Method", "Total Field/Scatter Field Method"))
    sim_param["output_filetype"] = st.selectbox(label="Output File Type",options=("Comma Separated Values (csv) ", "Numpy Array File (npy)", "Hierarchical Data Format version 5 (hdf5)"))
    st.write("")
    st.write("")
    sim_param["is_multithread"] = st.toggle("Enable multithreading in processing",disabled=(sim_param["algo"] == "Finite Difference Time Domain"))

sim_param["num_subdomain"] = st.select_slider("Number of subdomains",options=[1,2,4,8,16,32,64],disabled=not sim_param["is_multithread"])

if st.button("Create new simulation"):
    st.write("Created simulation!1")
    st.write(sim_param)