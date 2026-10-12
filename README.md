# PEM-Electrolyser
Masters research project


MATLAB/Simulink model of a PEM water electrolyser, built for the MEng thesis *Electrochemical Modelling and Performance Validation of a PEM Hydrogen Electrolyser* (Stellenbosch University). The model predicts cell voltage from current density and temperature, and calibrates its kinetic parameters against your own polarisation data.

## Files

- `PEMWE_Model.slx`: the Simulink model
- `run_model_calibration.m`: calibration script, run automatically by the model
- `PEMVoltage_fit.m`: helper function called by the calibration
- `user_data.xlsx`: **your own** polarisation data (see step 2)
- `all_polarisation_data.xlsx`: literature datasets used for this project only (see below)
- `calibration_only.m`: calibration script used for the thesis project only
- `validation_only.m`: validation script used for the thesis project only

## Using the model with your own data

### 1. Put everything in one folder

The Simulink model (`PEMWE_Model.slx`), `run_model_calibration.m`, any helper function files it calls (`PEMVoltage_fit.m`), and `user_data.xlsx` all need to sit in the same folder. That folder needs to be MATLAB's current working folder when you press Run.

`run_model_calibration.m` looks for `user_data.xlsx` by name in the working folder, not a fixed path. If the file is elsewhere, the calibration either fails outright or silently uses whatever old file happens to still be in the folder.

### 2. Format `user_data.xlsx` correctly

`user_data.xlsx` must contain **your own experimental data**. It is a placeholder name, and the file is not supplied with this repository.

Use one row per experimental point, with three columns and exactly these headers:

- `Temperature_C`
- `CurrentDensity_Acm2`
- `Voltage_V`

Every row recorded at a given temperature becomes one point on that temperature's polarisation curve.

### 3. Set the model's own Constant blocks

These are typed directly onto the canvas. They are not read from the Excel file or computed by the script.

- **Membrane Thickness:** plain micrometres (for example `127`), never pre-converted to centimetres. The script does that conversion internally, so typing an already-converted value makes the membrane appear thousands of times too thin.
- **n_cells**, **Active Area**, **Higher Heating Value**, **Faradaic efficiency**, and **water activity** (0 to 3).
- **Temperature** and **Current**, set for whatever operating point you are simulating. These are separate from the calibration data.

### 4. Wire the calibration into Simulink

Open Model Settings (gear icon), then Callbacks, then **InitFcn**, and paste in:

```matlab
run_model_calibration
```

This makes the calibration run automatically, once, every time you press Run, so it always uses whatever is currently in the folder.

### 5. Press Run and read the Command Window

You should see, in order:

1. How many temperatures were found, and what they are.
2. A warning if there are fewer than 4.
3. The mean held-out R² for Linear and Arrhenius, or a note that there were not enough temperatures to test.
4. Which form was selected, and its coefficients.
5. A final "Done" line stating the form and `i_lim`.

### 6. Calibrating a different dataset

Swap in that dataset's `user_data.xlsx` (same column names), update the Membrane Thickness block to its real thickness, and press Run again.

## About `all_polarisation_data.xlsx`, `calibration_only.m` and `validation_only.m`

These three files were made for the thesis project only.

- `all_polarisation_data.xlsx` holds the literature polarisation datasets used to calibrate and validate the model in the thesis. It is **not** an input to the workflow above, and you should not use it as `user_data.xlsx`. The data belong to the original authors, who are cited in the thesis.
- `calibration_only.m` and `validation_only.m` are the stand-alone calibration and validation scripts used to produce the thesis results from those datasets. They are not needed to run the Simulink model with your own data.

To calibrate the model, supply your own measurements in `user_data.xlsx` and follow the steps above.
