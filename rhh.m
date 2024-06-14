function rhh()

	window = figure('Name', "H&H Hydrolase", 'NumberTitle','off', 'MenuBar', 'none', 'Toolbar', 'figure'); %, 'Color', [0.8, 0.8, 0.8]);
	scsz = get(0,'ScreenSize');	%( 'Color','white', ...			% 'fullscreen', 'maximized'
	window.Position = [scsz(3)/3 scsz(4)/3.2 scsz(3)/1.6 scsz(4)/1.6];

	uicontrol('Units', 'normalized', 'Style', 'pushbutton', 'String', 'Reset', 'Pos', [0.03 0.35 0.09 0.04], 'Callback', @ResetAll, 'FontSize', 11);

    % Load data from file
    seq = load('dat.mat');

    time = seq.Traces(:,4); % milliseconds

    time = time - time(1); % Normalize time
	dt = time(2); % Time step
	
	correct = -3; % Bridge, Junction ~?

	% values from: https://neuronaldynamics.epfl.ch/online/Ch2.S2.html#Ch2.T1 
	% https://www.math.mcgill.ca/gantumur/docs/reps/RyanSicilianoHH.pdf new

		I = 0.1;
	% Maximum conductance (mS/cm^2)
    	gNa = [1.2, 1.2, 1.2];
    	gK = [0.36, 0.36, 0.36];
    	gL = [0.003, 0.003, 0.003];
	% Reversal Potential (mV)
    	ENa = [55.17, 55.17, 55.17];
    	EK = [-72.14, -72.14, -72.14];
    	EL = [-49.42, -49.42, -49.42];
	% Membrane Capacitance (uF/cm^2)
		Cm = 0.01;

	GateMult = [1,1,1,1,1,1];
	InitVars = {correct, gNa, gK, gL, ENa, EK, EL, Cm, I};

	Typos = ["KI", "KO", "WT"]; Gpyhos = ["gNa", "gK", "gL"]; Epyhos = ["ENa", "EK", "EL"]; Mpyhos = ["-mV", "I", "Cm"];
	GMypos = ["an", "bn", "am", "bm", "ah", "bh"];

	defaults = struct('Units', 'normalized', 'Style', 'edit', 'Callback', @TakeVars);

	gbuts = cell(3,3);
	ebuts = cell(3,3);
	mbuts = cell(1,3);
	gmuts = cell(1,6);

	for b = 1 : 3
		uicontrol('Units', 'normalized', 'Style', 'text', 'String', Typos(b), 'Pos', [0.05 * b 0.295 0.045 0.03], 'FontSize', 11, 'fontweight', 'bold');
		for c = 1 : 3
			if b == 1
				uicontrol('Units', 'normalized', 'Style', 'text', 'String', Gpyhos(c), 'Pos', [0.016 0.305 - c*0.04 0.045 0.03], 'FontSize', 11, 'fontweight', 'bold');
				uicontrol('Units', 'normalized', 'Style', 'text', 'String', Epyhos(c), 'Pos', [0.016 0.165 - c*0.04 0.045 0.03], 'FontSize', 11, 'fontweight', 'bold');
			end
			gbuts{b,c} = uicontrol(defaults, 'Pos', [0.05 * b 0.31 - c*0.04 0.045 0.03]);
			ebuts{b,c} = uicontrol(defaults, 'Pos', [0.05 * b 0.17 - c*0.04 0.045 0.03]);
		end
		uicontrol('Units', 'normalized', 'Style', 'text', 'String', Mpyhos(b), 'Pos', [0.2 0.305 - b*0.04 0.045 0.03], 'FontSize', 11, 'fontweight', 'bold', 'HorizontalAlignment','right');
		mbuts{1,b} = uicontrol(defaults, 'Pos', [0.25 0.31 - b*0.04 0.045 0.03]);
	end
	for a = 1:6
		uicontrol('Units', 'normalized', 'Style', 'text', 'String', GMypos(a), 'Pos', [0.3 0.31 - a*0.04 0.045 0.03], 'FontSize', 11, 'fontweight', 'bold', 'HorizontalAlignment','right');
		gmuts{1,a} = uicontrol(defaults, 'Pos', [0.35 0.31 - a*0.04 0.04 0.03]);
	end
    % Pre-allocate
		V_data = cell(1,3); phd = zeros(length(time), 1);
    	n = {phd, phd, phd}; m = {phd, phd, phd}; h = {phd, phd, phd};
    	I_Na = {phd, phd, phd}; I_K = {phd, phd, phd}; I_L = {phd, phd, phd};
		V_sim = {phd, phd, phd}; V = 0;

	GiveVars();
	DoCompute();

	function GiveVars()
		for i = 1 : 3
			set(gbuts{i,1}, 'String', num2str(gNa(i)));
			set(gbuts{i,2}, 'String', num2str(gK(i)));
			set(gbuts{i,3}, 'String', num2str(gL(i)));
			set(ebuts{i,1}, 'String', num2str(ENa(i)));
			set(ebuts{i,2}, 'String', num2str(EK(i)));
			set(ebuts{i,3}, 'String', num2str(EL(i)));
		end
		set(mbuts{1}, 'String', num2str(correct));
		set(mbuts{2}, 'String', num2str(I));
		set(mbuts{3}, 'String', num2str(Cm));
		for j = 1:6
			set(gmuts{j}, 'String', num2str(GateMult(j)));
		end
		%ValsnVars = {correct, gNa, gK, gL, ENa, EK, EL, Cm};
	end

	function TakeVars(~, ~)
		for i = 1 : 3
			gNa(i) = str2double(get(gbuts{i,1}, 'String'));
			gK(i) = str2double(get(gbuts{i,2}, 'String'));
			gL(i) = str2double(get(gbuts{i,3}, 'String'));
			ENa(i) = str2double(get(ebuts{i,1}, 'String'));
			EK(i) = str2double(get(ebuts{i,2}, 'String'));
			EL(i) = str2double(get(ebuts{i,3}, 'String'));
		end
		for j = 1:6
			GateMult(j) = str2double(get(gmuts{j}, 'String'));
		end
		correct = str2double(get(mbuts{1}, 'String'));
		I = str2double(get(mbuts{2}, 'String'));
		Cm = str2double(get(mbuts{3}, 'String'));
		DoCompute();
	end

	function ResetAll(~,~)
		correct = InitVars{1}; Cm = InitVars{8}; I = InitVars{9};
		gNa = InitVars{2}; gK = InitVars{3}; gL = InitVars{4};
		ENa = InitVars{5}; EK = InitVars{6}; EL = InitVars{7};
		GateMult = [1,1,1,1,1,1];
		GiveVars();
		DoCompute();
	end

	function DoCompute()
		Initialize();
		ComputeGates();
		ComputePotential();
		PlotUp();
	end

	function Initialize()
		V_data = {seq.Traces(:,1) - correct, seq.Traces(:,2) - correct, seq.Traces(:,3) - correct};
		V = {V_data{1}(1), V_data{2}(1), V_data{3}(1)};
    	[n{1}(1), m{1}(1), h{1}(1)] = GateVars(V{1});
		[n{2}(1), m{2}(1), h{2}(1)] = GateVars(V{2});
		[n{3}(1), m{3}(1), h{3}(1)] = GateVars(V{3});
	end

	function ComputeGates()
		for g = 1:3
			for i = 2 : length(time)
        		V = V_data{g}(i); % Initial
        		[an, bn, am, bm, ah, bh] = RateConst(V);
		
				n{g}(i) = n{g}(i-1) + dt * (an * (1 - n{g}(i-1)) - bn * n{g}(i-1));
				m{g}(i) = m{g}(i-1) + dt * (am * (1 - m{g}(i-1)) - bm * m{g}(i-1));
				h{g}(i) = h{g}(i-1) + dt * (ah * (1 - h{g}(i-1)) - bh * h{g}(i-1));
				%h{g}(i) = h{g}(i-1) + dt * (ah - (ah + bh) * h{g}(i-1));
		
        		% Ionic currents
        		I_Na{g}(i) = gNa(g) * m{g}(i)^3 * h{g}(i) * (V - ENa(g));
        		I_K{g}(i) = gK(g) * n{g}(i)^4 * (V - EK(g));
        		I_L{g}(i) = gL(g) * (V - EL(g));
			end
		end
	end


    function ComputePotential()
		for g = 1:3
        	V_sim{g}(1) = V_data{g}(1); % Initial
        	for i = 2:length(time)
            	I_tot = I - I_Na{g}(i) - I_K{g}(i) - I_L{g}(i);
            	V_sim{g}(i) = V_sim{g}(i-1) + (dt / Cm) * I_tot; % Euler
        	end
		end
	end

	function [n, m, h] = GateVars(V) 
    	[alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h] = RateConst(V);
    	n = alpha_n / (alpha_n + beta_n);
    	m = alpha_m / (alpha_m + beta_m);
    	h = alpha_h / (alpha_h + beta_h);
	end
	
	function [an, bn, am, bm, ah, bh] = RateConst(v)

    	an = GateMult(1) * 0.01*(v+50)/(1-exp(-(v+50)/10));
    	bn = GateMult(2) * 0.125*exp(-(v+60)/80);
    	am = GateMult(3) * 0.1*(v+35)/(1-exp(-(v+35)/10));
    	bm = GateMult(4) * 4.0*exp(-0.0556*(v+60));
    	ah = GateMult(5) * 0.07*exp(-0.05*(v+60));
    	bh = GateMult(6) * 1/(1+exp(-(0.1)*(v+30)));

    	%an = 0.02 * (V - 25) / (1 - exp(-(V - 25) / 9));
    	%bn = -0.002 * (V - 25) / (1 - exp((V - 25) / 9));
    	%am = 0.182 * (V + 35) / (1 - exp(-(V + 35) / 9));
    	%bm = -0.124 * (V + 35) / (1 - exp((V + 35) / 9));
    	%ah = 0.25 * exp(-(V + 90) / 12);
    	%bh = 0.25 * exp((V + 62) / 6) / exp((V + 90) / 12);
	end

	function PlotUp()
    	plotmv = subplot('Position', [0.025 0.44 0.45 0.5]);
    	plot(plotmv, ...
			time, V_data{1}, 'r', ...
			time, V_data{2}, 'g', ...
			time, V_data{3}, 'b', ...
			time, V_sim{1}, 'r--', ...
			time, V_sim{2}, 'g--', ...
			time, V_sim{3}, 'b--');
    	title('Membrane Potential', 'Units', 'normalized', 'Position', [0.88, 0.93, 1], 'Color', [0.3 0.3 0.3]);
    	xlabel('ms', 'Units', 'normalized', 'Position', [0.98, 0.06, 1], 'Color', [0.3 0.3 0.3]);
    	ylabel('mV', 'Units', 'normalized', 'Position', [0.03, 0.96, 1], 'Color', [0.3 0.3 0.3]);
	
    	plotna = subplot('Position', [0.52 0.7 0.45 0.29]);
    	plot(plotna, time, I_Na{1}, 'r', time, I_Na{2}, 'g', time, I_Na{3}, 'b');
    	title('I Na+', 'Units', 'normalized', 'Position', [0.96, 0.91, 1], 'Color', [0.3 0.3 0.3]);
	
    	plotk = subplot('Position', [0.52 0.37 0.45 0.29]);
    	plot(plotk, time, I_K{1}, 'r', time, I_K{2}, 'g', time, I_K{3}, 'b');
    	title('I K+', 'Units', 'normalized', 'Position', [0.96, 0.91, 1], 'Color', [0.3 0.3 0.3]);
	
    	plotl = subplot('Position', [0.52 0.04 0.45 0.29]);
    	plot(plotl, time, I_L{1}, 'r', time, I_L{2}, 'g', time, I_L{3}, 'b');
    	title('I Leak', 'Units', 'normalized', 'Position', [0.96, 0.91, 1], 'Color', [0.3 0.3 0.3]);
	
		curplots = [plotna, plotk, plotl];
		arrayfun(@(x) xlabel(x, 'ms', 'Units', 'normalized', 'Position', [0.98, 0.09, 1], 'Color', [0.3 0.3 0.3]), curplots);
		arrayfun(@(x) ylabel(x, 'uA/cm^2', 'Units', 'normalized', 'Position', [0.04, 0.87, 1], 'Color', [0.3 0.3 0.3]), curplots);
	end

end