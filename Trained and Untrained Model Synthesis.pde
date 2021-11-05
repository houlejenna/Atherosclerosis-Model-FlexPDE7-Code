{FlexPDE7 Script, imports data generated by trained and untrained model scripts and creates and exports histories (plots) of oxLDL, chemoattractants, cytokines,
oxLDL consumption, monocyte migration, total macrophages, and total foam cells in a LOW LDL system}

{FlexPDE7 documentation is available at: https://www.pdesolutions.com/help/}

TITLE 'Trained and Untrained Model Synthesis'
COORDINATES
	Cartesian1
VARIABLES
!variable names and tolerances; l=oxLDL, a=chemoattractants, c=cytokines, m=monocytes/macrophages, f=foam cells
l (threshold=0.001), a (threshold=0.001), c (threshold=0.001), m (threshold=0.001), f (threshold=0.001)
SELECT
ERRLIM=1e-3
threads=1
!Assigns a name to the T variable for plot axis formatting purposes
ALIAS(T) = "Time (Days)"
DEFINITIONS

    !Untrained model result imports ***imports files generated by "Untrained Source Model"
    !Naming convention: x#=species and # corresponding to value of sigma m 1-4, "file directory_file name"=import location and import file name
    l1 = table("Untrained Source Model_output/l(0)_USigmam_1")
    l2 = table("Untrained Source Model_output/l(0)_USigmam_2")
    l3 = table("Untrained Source Model_output/l(0)_USigmam_3")
    l4 = table("Untrained Source Model_output/l(0)_USigmam_4")
    a1 = table("Untrained Source Model_output/a(0)_USigmam_1")
    a2 = table("Untrained Source Model_output/a(0)_USigmam_2")
    a3 = table("Untrained Source Model_output/a(0)_USigmam_3")
    a4 = table("Untrained Source Model_output/a(0)_USigmam_4")
    c1 = table("Untrained Source Model_output/c(0)_USigmam_1")
    c2 = table("Untrained Source Model_output/c(0)_USigmam_2")
    c3 = table("Untrained Source Model_output/c(0)_USigmam_3")
    c4 = table("Untrained Source Model_output/c(0)_USigmam_4")
    m1 = table("Untrained Source Model_output/int(m)_USigmam_1")
    m2 = table("Untrained Source Model_output/int(m)_USigmam_2")
    m3 = table("Untrained Source Model_output/int(m)_USigmam_3")
    m4 = table("Untrained Source Model_output/int(m)_USigmam_4")
    f1 = table("Untrained Source Model_output/int(f)_USigmam_1")
    f2 = table("Untrained Source Model_output/int(f)_USigmam_2")
    f3 = table("Untrained Source Model_output/int(f)_USigmam_3")
    f4 = table("Untrained Source Model_output/int(f)_USigmam_4")
    
    !Trained model result imports ***imports files generated by "Trained Source Model"
    !Naming convention: Tx#=species and # corresponding to value of sigma m 1-4, "file directory_file name"=import location and import file name
    Tl1 = table("Trained Source Model_output/l(0)_TSigmam_1")
    Tl2 = table("Trained Source Model_output/l(0)_TSigmam_2")
    Tl3 = table("Trained Source Model_output/l(0)_TSigmam_3")
    Tl4 = table("Trained Source Model_output/l(0)_TSigmam_4")
    Ta1 = table("Trained Source Model_output/a(0)_TSigmam_1")
    Ta2 = table("Trained Source Model_output/a(0)_TSigmam_2")
    Ta3 = table("Trained Source Model_output/a(0)_TSigmam_3")
    Ta4 = table("Trained Source Model_output/a(0)_TSigmam_4")
    Tc1 = table("Trained Source Model_output/c(0)_TSigmam_1")
    Tc2 = table("Trained Source Model_output/c(0)_TSigmam_2")
    Tc3 = table("Trained Source Model_output/c(0)_TSigmam_3")
    Tc4 = table("Trained Source Model_output/c(0)_TSigmam_4")
    Tm1 = table("Trained Source Model_output/int(m)_TSigmam_1")
    Tm2 = table("Trained Source Model_output/int(m)_TSigmam_2")
    Tm3 = table("Trained Source Model_output/int(m)_TSigmam_3")
    Tm4 = table("Trained Source Model_output/int(m)_TSigmam_4")
    Tf1 = table("Trained Source Model_output/int(f)_TSigmam_1")
    Tf2 = table("Trained Source Model_output/int(f)_TSigmam_2")
    Tf3 = table("Trained Source Model_output/int(f)_TSigmam_3")
    Tf4 = table("Trained Source Model_output/int(f)_TSigmam_4")
    
    !Remaining definitions can be from either the trained or untrained model
    !This script doesn't solve PDE's, it combines the results of the two source models to create compound plots, FlexPDE requires the equation section and these definitions for every script
    
    !Diffusion constants; Dl and Da from rescaled experimental values, Dc assumed same as Da, Dm assumed many orders smaller than Dl
    Dl = 1.2e2
    Da = 5.4e3
    Dc = 5.4e3
    Dm = 1.0e-2
    
    !Decay constants; decay taken to be much smaller than other terms, values scaled to time, d_m assumed orders smaller than others
    d_l = 5.0e-1
    d_a = 5.0e-1
    d_c = 5.0e-1
    d_m = 5.0e-3
    
   !Rescalar constants; Mu_a, Mu_c, Mu_f rescaled to time and E_i (expression influx transporters)
    Mu_a = 4.0e1 !chemoattractants per % influx per macrophage
    Mu_c = 4.0e1 !cytokines per % influx per macrophage
    Mu_f = 4.1e-2 !unit LDL per foam cell
    
    !Flux constants (boundary conditions); rescaled to time
    Sigma_a1 = 5.0e2
    Sigma_a2 = 5.0e1
    Sigma_c = 5.0e2
    Xm = 5.0e-3
    
   !Bloodstream concentrations of chemoattractant and cytokines, chemoattractant production scalar, no time dependence
    A0 = 1.0e-1
    C0 = 1.0e-2
    Beta_a = 1.0e0
    
    !Rate of LDL and monocyte entry to the intima
    Sigma_l = 1.0e1
    Sigma_m = 2.0e-5 !2.0e-5 low, 3.0e-5, 4.1e-5, 9.0e-5, 2.0e-4 high
    
    !Epigenetic training modifiers; E_a=chemoattractant production rescalar, E_c=cytokine production rescalar
    E_a = 3 	!Untrained: E_a= 1, trained: E_a=3
    E_c = 4 	!Untrained: E_c=1, trained: E_c=4
    !E_i=influx transporter expression, E_e=efflux transporter expression
    E_i = 116.28 * (ln(l + 7.6) - 1.1)
    E_e = 78.00 * (ln(l + 4.6) - 0.4)
    
    Sigmam1 = 1.0e-6
    Sigmam2 = 4.0e-6
    Sigmam3 = 7.0e-6
    Sigmam4 = 1.0e-5
        
    E_i1 = 116.28 * (ln(l1 + 11.3) - 1.4)
    E_i2 = 116.28 * (ln(l2 + 11.3) - 1.4)
    E_i3 = 116.28 * (ln(l3 + 11.3) - 1.4)
    E_i4 = 116.28 * (ln(l4 + 11.3) - 1.4)
    TE_i1 = 116.28 * (ln(Tl1 + 7.6) - 1.1)
    TE_i2 = 116.28 * (ln(Tl2 + 7.6) - 1.1)
    TE_i3 = 116.28 * (ln(Tl3 + 7.6) - 1.1)
    TE_i4 = 116.28 * (ln(Tl4 + 7.6) - 1.1)
    
    E_e1 = 78.00 * (ln(l1 + 3.2) - 0.2)
    E_e2 = 78.00 * (ln(l2 + 3.2) - 0.2)
    E_e3 = 78.00 * (ln(l3 + 3.2) - 0.2)
    E_e4 = 78.00 * (ln(l4 + 3.2) - 0.2)
    TE_e1 = 78.00 * (ln(Tl1 + 4.6) - 0.4)
    TE_e2 = 78.00 * (ln(Tl2 + 4.6) - 0.4)
    TE_e3 = 78.00 * (ln(Tl3 + 4.6) - 0.4)
    TE_e4 = 78.00 * (ln(Tl4 + 4.6) - 0.4)
    
    Influx = m * l * E_i  / (1.0 + l)
    Efflux = m * l * E_e / (1.0 + l)
    
    ldlc1 = l1*m1 * E_i1 /(1.0+l1)
    ldlc2 = l2*m2* E_i2 /(1.0+l2)
    ldlc3 = l3*m3* E_i3 /(1.0+l3)
    ldlc4 = l4*m4* E_i4 /(1.0+l4)
    Tldlc1 = Tl1*Tm1* TE_i1 /(1.0+Tl1)
    Tldlc2 = Tl2*Tm2* TE_i2 /(1.0+Tl2)
    Tldlc3 = Tl3*Tm3* TE_i3 /(1.0+Tl3)
    Tldlc4 = Tl4*Tm4* TE_i4 /(1.0+Tl4)
    
    mm1 = if (a1 < A0) then 0 else Sigmam1*(1.0+C0*c1)*(a1-A0)
    mm2 = if (a2 < A0) then 0 else Sigmam2*(1.0+C0*c2)*(a2-A0)
    mm3 = if (a3 < A0) then 0 else Sigmam3*(1.0+C0*c3)*(a3-A0)
    mm4 = if (a4 < A0) then 0 else Sigmam4*(1.0+C0*c4)*(a4-A0)
    Tmm1 = if (Ta1 < A0) then 0 else Sigmam1*(1.0+C0*Tc1)*(Ta1-A0)
    Tmm2 = if (Ta2 < A0) then 0 else Sigmam2*(1.0+C0*Tc2)*(Ta2-A0)
    Tmm3 = if (Ta3 < A0) then 0 else Sigmam3*(1.0+C0*Tc3)*(Ta3-A0)
    Tmm4 = if (Ta4 < A0) then 0 else Sigmam4*(1.0+C0*Tc4)*(Ta4-A0)

EQUATIONS
	a:		Da*del2(a) = -E_a * Mu_a * Influx + d_a * a {decay} + dt(a)

	c:		Dc*del2(c) = -E_c * Mu_c * Influx + d_c * c {decay} + dt(c)
	
	l:		Dl*del2(l) = Influx {consumption} + d_l * l {decay} + dt(l)
	
	m:		Dm*del2(m) = Xm * div(m * grad(l)) {chemotaxis} + dt(f) {conversion} + d_m * m {dedifferentiation} + dt(m)
    
    	f:		dt(f) = Mu_f  * (Influx {influx} -  Efflux) {efflux}

BOUNDARIES
	REGION 1
        START(0)
			POINT LOAD(a) = Sigma_a1 * l / (Beta_a + l) + Sigma_a2 * c
			POINT LOAD(c) = -Sigma_c * c
			POINT LOAD(l) = Sigma_l
			POINT LOAD(m) = IF (a < A0) THEN 0 ELSE Sigma_m * (1.0 + C0 * c) * (a - A0)
			POINT LOAD(f) = 0
	LINE TO(1)
			POINT LOAD(a)=0
			POINT LOAD(c)=0
			POINT LOAD(l)=0
			POINT LOAD(m)=0
			POINT LOAD(f)=0

TIME 0 TO 30
plots
!From time 0 to 30 in steps of 30, plot PNGs will be saved at every step (if applicable)
	FOR T=0 BY 30 TO endtime
    
histories
!Compound plots for oxLDL concentration, chemoattractant concentration, cytokine concentration, total intimal macrophage content, total intimal foam cell content 
!history(x1, x2, x3, x4)=plot of x species in untrained macrophages at sigma m 1-4, as " ": plot titles, range(0, #): range y axis, PNG (1024, 2): save plot PNG 1024 pixels
    history(l1, l2, l3, l4) at (0) as "OxLDL at Endothelial Boundary" range(0, 20) PNG (1024, 2)
    history(Tl1, Tl2, Tl3, Tl4) at (0) as "OxLDL at Endothelial Boundary" range (0, 20) PNG (1024, 2)
    history(a1, a2, a3, a4) at (0) as "Chemoattractants at Endothelial Boundary "	range (0, 6000) PNG (1024, 2)
    history(Ta1, Ta2, Ta3, Ta4) at (0) as "Chemoattractants at Endotheilal Boundary" range (0, 6000) PNG (1024, 2)
    history(c1, c2, c3, c4) at (0) as "Cytokines at Endothelial Boundary" range (0,325) PNG (1024, 2)										
    history(Tc1, Tc2, Tc3, Tc4) at (0) as "Cytokines at Endotheial Boundary" range (0,325) PNG (1024, 2)
    history(integral(ldlc1), integral(ldlc2), integral(ldlc3), integral(ldlc4)) as "Rate of OxLDL Consumption" range (0, 17) PNG (1024, 2)
    history(integral(Tldlc1), integral(Tldlc2), integral(Tldlc3), integral(Tldlc4)) as "Rate of OxLDL Consumption" range (0, 17) PNG (1024, 2)
    history(mm1, mm2, mm3, mm4) at (0) as "Monocyte Migration" range(0, 25e-2) PNG (1024, 2)
    history(Tmm1, Tmm2, Tmm3, Tmm4) at (0) as "Monocyte Migration" range(0, 25e-2) PNG (1024, 2)
    history(integral(m1), integral(m2), integral(m3), integral(m4)) as "Total Intimal Macrophages" range (0,3.3) PNG (1024, 2)
    history(integral(Tm1), integral(Tm2), integral(Tm3), integral(Tm4)) as "Total Intimal Macrophages" range (0,3.3) PNG (1024, 2)
    history(f1, f2, f3, f4) as "Total Intimal Foam Cells" range (0,1.8) PNG (1024, 2)
    history(Tf1, Tf2, Tf3, Tf4) at (0) as "Total Intimal Foam Cells" range (0,1.8) PNG (1024, 2)

END
