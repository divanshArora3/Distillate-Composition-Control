Here are some steps for fitting an FOPDT model.

1. Run Simulink 2 Tank Model with step tests, doublet tests, etc to the pump rate (0-1).
2. Copy data_manual.txt contents into Excel sheet.
    a. Column 1 = Time
	b. Column 2 = Tank 2 Level Set point (not needed for fitting)
	c. Column 3 = Pump rate (0-1) - Manipulated variable / Controller output
	d. Column 4 = Tank 1 Level (not needed in this study)
	e. Column 5 = Tank 2 Level - Controlled variable / PV / Model and maintain this level with the controller
3. Copy Column 1 (Time), Column 3 (Pump rate), and Column 5 (Tank 2 Level) into the appropriate locations in the Excel sheet
4. Optimize values of Kp, tau_p, theta_p to best fit the data