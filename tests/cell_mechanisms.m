% cell_mechanisms.m

% mechanism = 'low-threshold Ca';
% mechanism = 'somatic components';
% mechanism = 'amarillo_Ttype'
mechanism = 'Ca activated K';
% mechanism = 'high_thresh_vanilla';
% mechanism = 'h_current_tests';
% mechanism = 'high_thresh_gerstner';
% mechanism = 'schweighofer_gap';
% mechanism = 'somatic mechanisms'
% mechanism = 'axonal_components';
% mechanism = 'calcium_update';

Vrange = [-150:30];

switch mechanism

	case 'low-threshold Ca'
		% NOTE: not GHK! Consider the non linearity

		% [================================================]
		% 		 somatic low threshold calcium
		% [================================================]
		k_inf = @(V_soma) (1 ./ (1 + exp(-1 * (V_soma + 61)   ./ 4.2)));
		l_inf = @(V_soma) (1 ./ (1 + exp((     V_soma + 85.5) ./ 8.5)));

	    tau_k = 1;
	    tau_l = @(V_soma) ((20 * exp((V_soma + 160) ./ 30) ./ (1 + exp((V_soma + 84) ./ 7.3))) +35);
	   
	     figure
	     plot(Vrange, k_inf(Vrange).^3 ,'b'), hold on
	     plot(Vrange, l_inf(Vrange),'g'), hold on
	     plot(Vrange, (k_inf(Vrange).^3) .* l_inf(Vrange),'y'), hold on
	     title('Calcium Low Threshold')
		legend({'k_inf^3' ; 'l_inf' ; })
		title('inward low threshold calcium')





	case 'somatic components'

	 % % Dendrite-soma interaction current
	 %    I_ds  = (g_int / p1) * (V_soma - V_dend); 
	 %    % Inward low-threshold Ca current
	 %    I_CaL = g_CaL * k * k * k * l * (V_soma - V_Ca);
	 %    % Inward Na current
	 %    I_Na_s  = g_Na_s * m * m * m * h * (V_soma - V_Na);
	 %    % Leak current
	 %    I_ls  = g_ls * (V_soma - V_l);
	 %    % H current soma
	 %    I_h_s    =  g_h_s * q_s * (V_soma - V_h); %# change
	 %    % Potassium current
	 %    I_Kdr_s = g_Kdr_s * n * n * n * n * (V_soma - V_K);
	 %    I_K_s   = g_K_s * (x_s ^ 4) * (V_soma - V_K);
	 %    % Axon-soma interaction current
	 %    I_as    = (g_int / (1 - p2)) * (V_soma - V_axon);
	 %    %***** AMPA current *****
	 %    I_amp   = gbar_ampa * g_ampa * (V_dend - V_ampa);
	 %    %***** GABA A current *****
	 %    I_gab_soma   = gbar_gaba_soma * g_gaba_soma * (V_soma - V_gaba_soma);
	        




        
    % somatic Sodium
    m_inf = @(V_soma) 1 ./ (1 + (exp((-30 - V_soma)./ 5.5)));
    h_inf = @(V_soma) 1 ./ (1 + (exp((-70 - V_soma)./-5.8)));
    tau_h = @(V_soma) 3 * exp((-40 - V_soma)./33);

    % Potassium
     n_inf = @(V_soma) 1 ./ (1 + exp( ( -3 - V_soma) ./  10));	
     tau_n = @(V_soma) 5 + (  47 * exp( -(-50 - V_soma) ./  900));
     % n     = delta * dn_dt + Potassium_n;
     
     alpha_x_s = @(V_soma) 0.13 * (V_soma + 25) ./ (1 - exp(-(V_soma + 25) ./ 10));
     beta_x_s  = @(V_soma) 1.69 * exp(-0.0125 * (V_soma + 35));
    
     x_inf_s   = @(V_soma) alpha_x_s(V_soma) ./ (alpha_x_s(V_soma) + beta_x_s(V_soma));
     tau_x_s   = @(V_soma) 1 ./ (alpha_x_s(V_soma) + beta_x_s(V_soma));


     
    

	case 'amarillo_Ttype' % Amarillo 2D CaT

		V1_2mT = -53;
		V_tau_m1 = 2;
		V_tau_m2 = 2;
		% m_T_inf = @(V) 1./(1 + exp((V −- V1_2mT)/ − 6.2);
		% tau_m_T = @(V) (0.612 + 1./(exp((V - −V_tau_m1 )/ −- 16.7) + exp((V −- V_tau_m2 )/18.2)))/3   ;
		% h_T_inf = @(V) 1/(1 + exp[(V −- V1_2mT/4]);
		% tau_h1_T = @(V) (exp[(V −- V_tau_h1 )/66.6])/3; % -- for V < −75 mV
		% tau_2_T = @(V) (28 + exp[(V -− V_tau_h2 )/ − 10.5])/3; % -- for V > −75 mV

		% G(V, Cao, Cai) = z^2F^2V/RT (Cai- −Ca0*exp(−zFV/RT))/ (1 −- exp(−zFV/RT)),

		% dhT/dt = (hT∞(V) − hT)/τhT(V),
		% I = p m
		% h SG(V,Ca,Ca)
		% TTTToi
		% dmT/dt = (mT∞(V) − mT)/τmT(V)


	case 'calcium_update'

		% update Calcium concentration
        dCa_dt = -3 * I_CaH - 0.075 * Ca2Plus;
        Ca2Plus = delta * dCa_dt + Ca2Plus;

		dCa_dt = @(I_CaH, Ca2Plus) -3 * I_CaH - 0.075 * Ca2Plus;

		% I_CaH = r^2 * V_dend;



	case 'high_thresh_vanilla'
		% [================================================]
		% 		 Calcium High Threshold (Ir)
		% [================================================]
		%  activation curves from r

		alpha_r = @(V_dend) 1.7 ./ (1 + exp( -(V_dend - 5) ./ 13.9));
		beta_r  = @(V_dend) 0.02 .* (V_dend + 8.5)./ (exp((V_dend + 8.5) ./ 5) - 1);
		tau_r   = @(V_dend) 5 ./ (alpha_r(V_dend) + beta_r(V_dend));
		r_inf   = @(V_dend) alpha_r(V_dend) ./ (alpha_r(V_dend) + beta_r(V_dend));

		% beta_r_new = @(V_dend)   2.8 ./(1 + exp( (V_dend + 68)  ./  31 )) 
		% tau_r_new   = @(V_dend) 5 ./ (alpha_r(V_dend) + beta_r_new(V_dend));
		% r_inf_new   = @(V_dend) alpha_r(V_dend) ./ (alpha_r(V_dend) + beta_r_new(V_dend));
		

		% dr_dt_new = @(V_dend, r) (r_inf_new(V_dend) - r) ./ tau_r_new(V_dend);
		dr_dt = @(V_dend, r) (r_inf(V_dend) - r) ./ tau_r(V_dend);

		Vrange = -150:20;
		figure
		plot(Vrange, alpha_r(Vrange),'r'), hold on
		plot(Vrange, beta_r(Vrange),'b')
		% plot(Vrange, beta_r_new(Vrange),'b','linewidth',2)
		plot(Vrange, r_inf(Vrange) ,'g')
		% plot(Vrange, r_inf_new(Vrange) ,'g','linewidth',2)
		plot(Vrange, tau_r(Vrange) ,'c')
		% plot(Vrange, tau_r_new(Vrange) ,'c','linewidth',2)

		legend({'alpha' 'beta' 'r inf' 'tau r'})
		title(mechanism)


		figure
		c= 0; rrange = linspace(0,.1,10); lr  = 10;
		for Cr = rrange
			c = c+1;
			plot(Vrange, dr_dt(Vrange,Cr) ,'color' , [0 0 1]*c/lr)
			hold on
			
		end
		legend(num2str( rrange' ))
		title('dr dt') ; xlabel('V dend')

	case 'high_thresh_gerstner'
		% [================================================]
		% 		 Calcium High Threshold (Ir)
		% [================================================]
		%  activation curves from r

		alpha_r = @(V_dend) 1.7 ./ (1 + exp( -(V_dend - 5) ./ 13.9));
		beta_r = @(V_dend)   2.8 ./(1 + exp( (V_dend + 68)  ./  31 )) 
		tau_r   = @(V_dend) 5 ./ (alpha_r(V_dend) + beta_r_new(V_dend));
		r_inf   = @(V_dend) alpha_r(V_dend) ./ (alpha_r(V_dend) + beta_r_new(V_dend));
		

		% dr_dt_new = @(V_dend, r) (r_inf_new(V_dend) - r) ./ tau_r_new(V_dend);
		dr_dt = @(V_dend, r) (r_inf(V_dend) - r) ./ tau_r(V_dend);

		Vrange = -150:20;
		figure
		plot(Vrange, alpha_r(Vrange),'r'), hold on
		plot(Vrange, beta_r(Vrange),'b')
		% plot(Vrange, beta_r_new(Vrange),'b','linewidth',2)
		plot(Vrange, r_inf(Vrange) ,'g')
		% plot(Vrange, r_inf_new(Vrange) ,'g','linewidth',2)
		plot(Vrange, tau_r(Vrange) ,'c')
		% plot(Vrange, tau_r_new(Vrange) ,'c','linewidth',2)

		legend({'alpha' 'beta' 'r inf' 'tau r'})
		title(mechanism)


		figure
		c= 0; rrange = linspace(0,.1,10); lr  = 10;
		for Cr = rrange
			c = c+1;
			plot(Vrange, dr_dt(Vrange,Cr) ,'color' , [0 0 1]*c/lr)
			hold on
			
		end
		legend(num2str( rrange' ))
		title('dr dt') ; xlabel('V dend')

	case 3

		% [================================================]
		% 	substitute beta_r that grows without bounds (?)
		% [================================================]


		beta_r = @(V_dend) 0.02 .* (V_dend + 8.5)./ (exp((V_dend + 8.5) ./ 5) - 1);
		sig =  @(a, X) a(1) ./ (1 + exp( -(X + a(2)) ./ a(3)));

		xdata = [-180:100];
		ydata = beta_r(xdata);
		ydata(1:80) = 1.8;

		a = lsqcurvefit(sig,[eps eps 0],xdata, ydata)

		new_beta_r = @(V_dend)   2.8 ./(1 + exp( (V_dend + 68)  ./  31 )) 
		clf
		plot(xdata, sig(x,xdata),'g')
		hold on
		plot(xdata, new_beta_r(xdata),'y')

		plot(xdata, ydata,'r')
		plot(xdata, alpha_r(xdata),'r')
		title(mechanism)

	case 'schweighofer_gap'


		% [================================================]
		% 		 Gap 
		% [================================================]
		% Schweighofer's gap non linearity (2004)
		fgap = @(DeltaV) (0.8 .* exp(-1.*DeltaV.*DeltaV/100) + 0.2);
		plot(x, fgap(x))
		title(mechanism)


	case 'Ca activated K'


		% [================================================]
		% 		 Ca activated K hyperpolarizing current - s
		% [================================================]     

		%NEW (CaPlus in uM/l)
		alpha_s = @(Ca2Plus, V_dend) .9e-3*log(Ca2Plus)*exp(V_dend/24);
		beta_s  = @(Ca2Plus, V_dend) .75e-3*exp(-V_dend/24);
		s_inf = @(Ca2Plus, V_dend) alpha_s(Ca2Plus, V_dend) ./ (alpha_s(Ca2Plus, V_dend) + beta_s(Ca2Plus, V_dend));
		tau_s = @(Ca2Plus, V_dend) .1 ./ (alpha_s(Ca2Plus,V_dend) + beta_s(Ca2Plus,V_dend));


		% Original : No dependence on V_dend
		alpha_s_orig =@(Ca2Plus) (0.00002 * Ca2Plus) .* (0.00002 * Ca2Plus < 0.01) + 0.01*((0.00002 * Ca2Plus)> 0.01); % nefarious calcium creates instabilities?
		beta_s_orig = 0.015;
		s_inf_orig = @(Ca2Plus) alpha_s_orig(Ca2Plus) ./ (alpha_s_orig(Ca2Plus) + beta_s_orig);
		tau_s_orig = @(Ca2Plus) 1 ./ (alpha_s_orig(Ca2Plus) + beta_s_orig);



			Ca_range = [1 10 100 1000];
			Vd_range = [-130:10];
			clf
			c = 0;
			for ca = Ca_range
				c = c+1;

				figure(2)
				plot(Vd_range, ones(size(Vd_range))*alpha_s_orig( ca),'color', [1 1 0]*c/6), hold on
				plot(Vd_range, alpha_s( ca, Vd_range),'color', [1 0 0]*c/6), hold on
				title('alpha')

				figure(3)
				plot(Vd_range, ones(size(Vd_range))*beta_s_orig,'color', [1 1 0]*c/6), hold on
				plot(Vd_range, beta_s( ca, Vd_range),'color', [1 1 0]*c/6), hold on
				title('beta')

				figure(5)
				plot(Vd_range, ones(size(Vd_range))*s_inf_orig( ca),'color', [0 0 1]*c/6), hold on
				plot(Vd_range, s_inf( ca, Vd_range),'color', [0 0 1]*c/6), hold on
				title('s_inf')

				figure(4)
				plot(Vd_range, tau_s( ca, Vd_range),'color', [0 0 1]*c/6), hold on
				plot(Vd_range, ones(size(Vd_range))*tau_s_orig( ca),'color', [0 0 1]*c/6), hold on
				title('tau_s')

			end

		% plot(x, beta_r(x),'b')
		% plot(x, s_inf(alpha_r(x), beta_r(x)) ,'g')


		ds_dt = @(s_inf, Potassium_s, tau_s	) (s_inf - Potassium_s) / tau_s;

		s = @(delta, ds_dt, Potassium_s) delta * ds_dt + Potassium_s;



	case 'h_current_tests'


		% [================================================]
		% 		 H current
		% [================================================]

	    % #NOTE: new h current has a very steep inactivation.		
		        
		q_inf_old = @(V_dend) 1 ./(1 + exp((V_dend + 80) / 4));
		tau_q_old = @(V_dend) 1 ./(exp(-0.086 * V_dend - 14.6) + exp(0.070 * V_dend - 1.87));


		 % g_h * q * (V_dend + 43)

		 mAlpha	 = @(v) 1e-3*6.43*(v+154.9)./(exp((v+154.9)./11.9)-1); % w V = 154;
		 mBeta	 = @(v) 1e-3*193*exp(v/33.1);

		q_inf_new = @(V_dend)  mAlpha(V_dend) ./(mAlpha(V_dend) + mBeta(V_dend));
		tau_q_new = @(V_dend)       1         ./(mAlpha(V_dend) + mBeta(V_dend));
		% 	if(v == -154.9){
		% 	v = v + 0.000001
		% }
		% From channelpedia
		% mAlpha = 0.001*6.43*(v+154.9)/(exp((v+154.9)/11.9)-1) 
		% mBeta = 0.001*193*exp(v/33.1)
		% mInf = mAlpha/(mAlpha + mBeta)
		% mTau = 1/(mAlpha + mBeta)

		% from Amarillo et al. 2015
		% q_inf_new = @(V) 1./(1 +exp(V+82)/5.49);
		% tau_q_new = @(V) (1./(.8e-3 + 3.5e-6*exp(-0.05787 + V) +exp(-1.87+0.0701*V)))/1.32;


		dq_dt_old = @(V_dend, q) (q_inf_old(V_dend) - q) ./ tau_q_old(V_dend);
		dq_dt_new = @(V_dend, q) (q_inf_new(V_dend) - q) ./ tau_q_new(V_dend);



		plot(Vrange, q_inf_old(Vrange),'r','linewidth', 1.5), hold on, 
		plot(Vrange, q_inf_new(Vrange),'g','linewidth', 1.5)
		title('H current q inf'), xlabel('Vm')
		legend({'q inf old'; 'q inf new'})
		title(mechanism)
		

		figure
			plot(Vrange, tau_q_old(Vrange),'r'), hold on,
			plot(Vrange, tau_q_new(Vrange),'g')
			legend({'tau_old' ;'tau_q_new'}),  xlabel('Vm')
			title('H current')
		
		figure 
			plot(Vrange, mAlpha(Vrange),'r'), hold on
			plot(Vrange, mBeta(Vrange),'g')
			legend({'alpha' ;'beta'}),  xlabel('Vm')

		figure
			plot(Vrange, q_inf_new(Vrange),'r'), hold on
			plot(Vrange, q_inf_old(Vrange),'g')
			legend({'q_inf_new' ;'q_inf_old'}),  xlabel('Vm')




		% figure
		% for qr = linspace(0,1,10)	
		% 	plot(Vrange, dq_dt_old(Vrange, qr),'color', [1 0 0]*qr), hold on
		% 	plot(Vrange, dq_dt_new(Vrange, qr),'color', [0 1 0]*qr)

		% end
		%  xlabel('Vm'); title('dq / dt')
		%  line([-43 -43], [-1 1],'linestyle', '--')


	case 'axonal_components'



		% [================================================]
		% 		 Axonal components
		% [================================================]

		% discontinuity at -25

			m_inf_a   = @(V_axon) 1 / (1 + (exp((-30 - V_axon)/ 5.5)));
		    h_inf_a   = @(V_axon) 1 / (1 + (exp((-60 - V_axon)/-5.8)));
		    tau_h_a   = @(V_axon)    1.5 * exp((-40 - V_axon)/33);

		    % dh_dt_a   = (h_inf_a - Sodium_h_a)/tau_h_a;
     
		     % m_a       = m_inf_a;
		     % h_a       = Sodium_h_a + delta * dh_dt_a;

		% Update potassium components
		     alpha_x_a = @(V_axon) 0.13 * (V_axon + 25) ./ (1 - exp(-(V_axon + 25) / 10));
		     beta_x_a  = @(V_axon) 1.69 * exp(-0.0125 * (V_axon + 35));

		     x_inf_a   = @(alpha_x_a, beta_x_a) alpha_x_a ./ (alpha_x_a + beta_x_a);
		     tau_x_a   = @(alpha_x_a, beta_x_a) 1 ./ (alpha_x_a + beta_x_a);
		    
		     dx_dt_a   = @(A, B, x_a) (x_inf_a(A,B) - x_a) ./ tau_x_a(A,B);
		    
		     % x_a       = delta * dx_dt_a + Potassium_x_a;
		Vrange= [-200:0.1:300];
		plot(Vrange,alpha_x_a(Vrange));hold on
		plot(Vrange,beta_x_a(Vrange),'r')
		plot(Vrange,x_inf_a(alpha_x_a(Vrange), beta_x_a(Vrange)) ,'g')
		plot(Vrange,tau_x_a(alpha_x_a(Vrange), beta_x_a(Vrange)) ,'c')
		legend({'alphaxa' ; 'beta_x_a' ; 'x_inf_a'; 'tau_x_a'})

		for x_a = linspace(0,1,10)
		plot(Vrange, dx_dt_a(alpha_x_a(Vrange), beta_x_a(Vrange),x_a) ,'c')
		hold on
		end


	case 8 % dendritic mechanisms

		
		alpha_r = @(V_dend) 1.7 ./ (1 + exp( -(V_dend - 5) ./ 13.9));
		beta_r  = @(V_dend) 0.02 .* (V_dend + 8.5)./ (exp((V_dend + 8.5) ./ 5) - 1);
		tau_r   = @(V_dend) 5 ./ (alpha_r(V_dend) + beta_r(V_dend));
		r_inf   = @(V_dend) alpha_r(V_dend) ./ (alpha_r(V_dend) + beta_r(V_dend));

		% Ih
		q_inf = @(V_dend) 1 ./(1 + exp((V_dend + 80) / 4));
		tau_q = @(V_dend) 1 ./(exp(-0.086 * V_dend - 14.6) + exp(0.070 * V_dend - 1.87));

		alpha_s = @(Ca2Plus, V_dend) .9e-3*Ca2Plus*exp(V_dend/24);
		beta_s  = @(Ca2Plus, V_dend) .75e-3*exp(-V_dend/24);
		s_inf = @(Ca2Plus, V_dend) alpha_s(Ca2Plus, V_dend) ./ (alpha_s(Ca2Plus, V_dend) + beta_s(Ca2Plus, V_dend));
		tau_s = @(Ca2Plus, V_dend) .1 ./ (alpha_s(Ca2Plus,V_dend) + beta_s(Ca2Plus,V_dend));


		vrange = [-150:20];

		figure
		plot(vrange, 4.5*r_inf(vrange).^2   .*(vrange - 120) ,'r'), hold on
		plot(vrange, 1.2*q_inf(vrange)      .*(vrange + 43 ) ,'g')
		plot(vrange, 55*s_inf(1, vrange)    .*(vrange + 75 ) ,'b')
		plot(vrange, 55*s_inf(10, vrange)   .*(vrange + 75 ) ,'b')
		plot(vrange, 55*s_inf(100, vrange)  .*(vrange + 75 ) ,'b','linewidth',2)

		figure
		plot(vrange, tau_r(vrange).^2   ,'r'), hold on
		plot(vrange, tau_q(vrange)      ,'g')
		plot(vrange, tau_s(1, vrange)   ,'b')
		plot(vrange, tau_s(10, vrange)  ,'b')
		plot(vrange, tau_s(100, vrange) ,'b','linewidth',2)


		figure
		plot(vrange, 4.5*r_inf(vrange).^2  .*(vrange - 120)./tau_r(vrange).^2    ,'r'), hold on
		plot(vrange, 1.2*q_inf(vrange)     .*(vrange + 43 )./tau_q(vrange)       ,'g')
		plot(vrange, 55*s_inf(1, vrange)   .*(vrange + 75 )./tau_s(1, vrange)    ,'b')
		plot(vrange, 55*s_inf(10, vrange)  .*(vrange + 75 )./tau_s(10, vrange)   ,'b')
		plot(vrange, 55*s_inf(100, vrange) .*(vrange +75  )./tau_s(100, vrange)  ,'b','linewidth',2)




		A = [4.5*r_inf(vrange).^2  .*(vrange - 120);
			 1.2*q_inf(vrange)     .*(vrange + 43 );
			 35*s_inf(1, vrange)   .*(vrange + 75 );
			 35*s_inf(10, vrange)   .*(vrange + 75 );
			 35*s_inf(100, vrange)  .*(vrange +75  )];
			 figure
			 area(vrange, A')


		case 'somatic mechanisms'
				% [================================================]
				% 		 somatic mechanisms
				% [================================================]

			%% update somatic components

			k_inf = @(V_soma) (1 / (1 + exp(-1 * (V_soma + 61)   / 4.2)));
			l_inf = @(V_soma)(1 / (1 + exp((     V_soma + 85.5) / 8.5)));

			tau_k = 1;
			tau_l = @(V_soma)((20 * exp((V_soma + 160) / 30) / (1 + exp((V_soma + 84) / 7.3))) +35);
			    
			dk_dt = (k_inf - Calcium_k) / tau_k;
			dl_dt = (l_inf - Calcium_l) / tau_l;

			    % k = delta * dk_dt + Calcium_k;
			    % l = delta * dl_dt + Calcium_l;

			m_inf = 1 / (1 + (exp((-30 - V_soma)/ 5.5)));
			h_inf = 1 / (1 + (exp((-70 - V_soma)/-5.8)));
			tau_h =       3 * exp((-40 - V_soma)/33);

			 m_inf_a   = @(V_axon) 1 / (1 + (exp((-30 - V_axon)/ 5.5)));
		     h_inf_a   = @(V_axon) 1 / (1 + (exp((-60 - V_axon)/-5.8)));
		     tau_h_a   = @(V_axon)     1.5 * exp((-40 - V_axon)/33);

		     % dh_dt_a   =  (h_inf_a  - Sodium_h_a)/tau_h_a;
		     
		     m_a       = m_inf_a;

		     % Update potassium components
		     alpha_x_a = @(V_axon) 0.13 * (V_axon + 25) / (1 - exp(-(V_axon + 25) / 10));
		     beta_x_a  = @(V_axon) 1.69 * exp(-0.0125 * (V_axon + 35));
		    
		     x_inf_a   = alpha_x_a (V_axon)/ (alpha_x_a + beta_x_a);
		     tau_x_a   =         1 / (alpha_x_a + beta_x_a);
		    
		     dx_dt_a   = (x_inf_a - Potassium_x_a) / tau_x_a; 
		    


		figure
		plot(Vrange, k_inf(Vrange),'r'), hold on
		plot(Vrange, l_inf(Vrange),'g')
		plot(Vrange, tau_l(Vrange),'b')
		plot(Vrange, m_inf(Vrange),'c')
		plot(Vrange, s_inf(Vrange),'b','linewidth',2)

		



		A = [4.5*r_inf(Vrange).^2  .*(Vrange - 120);
			 1.2*q_inf(Vrange)     .*(Vrange + 43 );
			 35*s_inf(1, Vrange)   .*(Vrange + 75 );
			 35*s_inf(10, Vrange)   .*(Vrange + 75 );
			 35*s_inf(100, Vrange)  .*(Vrange +75  )];
			 figure
			 area(Vrange, A')





end



 % (1 / (1 + exp(-1 * (V_soma + 61)   / 4.2)))

