## Example analysis pipeline for data derived from a SL-SS3xpress scRNA-seq experiment
### Git clone this repository
```
git clone https://github.com/Siegel-Lab/SL-SS3x-protocol.git
git clone https://github.com/<your-username>/<your-repo>.git
cd SL-SS3x-protocol
```

### Create conda environment
```
 conda env create -f SL_SS3xpress_env.yaml
```
### Run the preprocessing script
```
 bash run.sh
```
### After the preprocessing is complete, launch JupyterLab:
```
jupyter-lab
```
In the JupyterLab interface, open:

```
SL_Smart-seq3xpress_protocol_downstream_pipeline.ipynb
```
Follow the notebook cells step-by-step to generate analysis results and plots.

### Example output plot

<img src="./SL_SS3x_example_plot.png">
