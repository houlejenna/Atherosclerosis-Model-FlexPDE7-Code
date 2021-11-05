{FlexPDE7 Script, describes solutions to model equations, produces and exports result plots}

{FlexPDE7 documentation is available at: https://www.pdesolutions.com/help/}

TITLE 'Trained Source Model'
COORDINATES
	Cartesian1
VARIABLES
!variable names and tolerances; l=oxLDL, a=chemoattractants, c=cytokines, m=monocytes/macrophages, f=foam cells
l (threshold=0.001), a (threshold=0.001), c (threshold=0.001), m (threshold=0.001), f (threshold=0.001)
SELECT
ERRLIM=1e-3
threads=1
DEFINITIONS
	
    !Diffusion constants; Dl and Da from rescaled experimental values, Dc assumed same as Da, Dm assumed many orders smaller than Dl
    Dl = 5.2e1		!5.2e1, from Dabagh rescaled to 24 hr time and space according to Chalmers
    Da = 5.4e3		!5.4e3, from Paavola rescaled to 24 hr time and space according to Chalmers
    Dc = 5.4e3		!5.4e3, assumed the same as Da 
    Dm = 5.2e-3 	!4 orders of magnitude smaller than Dl, same rescale as Chalmers 
    
    !Decay constants; decay taken to be much smaller than other terms, values scaled to time, d_m assumed orders smaller than others
    d_l = 5.0e-1 
    d_a = 5.0e-1		
    d_c = 5.0e-1
    d_m = 5.0e-3	!2 orders smaller than other decays, same rescale as Chalmers
    
   !Rescalar constants; Mu_a, Mu_c, Mu_f rescaled to time and E_i (expression influx transporters)
    Mu_a = 4.0e1 	!chemoattractants per % influx per macrophage 
    Mu_c = 4.0e0 	!cytokines per % influx per macrophage
    Mu_f = 4.0e-2 	!unit LDL per foam cell
    
    !Flux constants (boundary conditions); rescaled to time
    Sigma_a1 = 5.0e2
    Sigma_a2 = 5.0e0
    Sigma_c = 5.0e-3
    Xm = 5.0e-3
    
   !Bloodstream concentrations of chemoattractant and cytokines, chemoattractant production scalar, no time dependence, same as Chalmers
    A0 = 1.0e-1
    C0 = 1.0e-2 
    Beta_a = 1.0e0
    
    !Rate of LDL and monocyte entry to the intima
    Sigma_l = 1.0e1 	!Low LDL: 1.0e1, high LDL: 4.0e1
    Sigma_m = 1.0e-6 	!Range: (low) 1.0e-6, 4.0e-6, 7.0e-6, 1.0e-5 (high)
    
    !Epigenetic training modifiers; E_a=chemoattractant production rescalar, E_c=cytokine production rescalar
    E_a = 3 	!Untrained: E_a= 1, trained: E_a=3
    E_c = 4 	!Untrained: E_c=1, trained: E_c=4
    !E_i=influx transporter expression, E_e=efflux transporter expression
    E_i = 116.28 * (ln(l + 7.6) - 1.1)
    E_e = 78.00 * (ln(l + 4.6) - 0.4)
    
    !Influx and efflux
    Influx = m * l * E_i  / (1.0 + l)
    Efflux = m * l * E_e / (1.0 + l)

EQUATIONS
!PDEs for the 5 species, rearranged to match FlexPDE requirements
    l:		Dl*del2(l) {diffusion from lumen} = Influx {consumption} + d_l * l {decay} + dt(l) {PDE of l with respect to time}
    
    a:		Da*del2(a) {diffusion from lumen} = -E_a  * Mu_a * Influx {consumption} + d_a * a {decay} + dt(a) {PDE of a with respect to time}

    c:		Dc*del2(c) {diffusion from lumen} = -E_c * Mu_c * Influx {consumption} + d_c * c {decay} + dt(c) {PDE of c with respect to time}
	
    m:		Dm*del2(m) {diffusion from lumen} = Xm * div(m * grad(l)) {chemotaxis} + dt(f) {conversion} + d_m * m {dedifferentiation} + dt(m) {PDE of m with respect to time}
    
    f:		dt(f) {PDE of f with respect to time} = Mu_f  * (Influx - Efflux)

BOUNDARIES
!Boundary conditions
	REGION 1
        START(0)
        !Conditions at the endothelial boundary, x=0
			POINT LOAD(a) = Sigma_a1 * l / (Beta_a + l) + Sigma_a2 * c
			POINT LOAD(c) = -Sigma_c * c
			POINT LOAD(l) = Sigma_l
			POINT LOAD(m) = IF (a < A0) THEN 0 ELSE Sigma_m * (1.0 + C0 * c) * (a - A0)	!Representation of the Heaviside function
			POINT LOAD(f) = 0
	LINE TO(1)
        !Conditions at the media, x=1
			POINT LOAD(a)=0
			POINT LOAD(c)=0
			POINT LOAD(l)=0
			POINT LOAD(m)=0
			POINT LOAD(f)=0

!Time period of 30 days
TIME 0 TO 30
plots
FOR T=0 BY 1 TO 30
   
histories
!Plots for oxLDL concentration, chemoattractant concentration, cytokine concentration, total intimal macrophage content, total intimal foam cell content
!"export file" commands save a png file of each plot with name "name" to the folder where the FlexPDE script is saved
!Naming convention: x(0)=concentration of x species at x=0, int(x)=integral of x species (m or f), T=trained, L/H=low/highLDL, Sigmam_x= sigma m value 1-4
    history(l) at (0) !export file "l(0)_THSigmam_1"			!Modified LDL Endotheilal Boundary
    history(a) at (0) !export file "a(0)_THSigmam_1"			!Chemoattractants on Endothelial Boundary
    history(c) at (0) !export file "c(0)_THSigmam_1"			!ES Cytokines on Endothelial Boundary
    history(integral(m)) !export file "int(m)_THSigmam_1"		!Monocytes / Macrophages Intimal Total
    history(integral(f)) !export file "int(f)_THSigmam_1"		!Foam Cells Intimal Total

END
