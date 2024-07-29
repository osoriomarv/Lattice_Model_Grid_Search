%% 3-D Grid Search Lattice Strain Model
%Fits for Youngs Modulus, Partition Coefficient and Ionic Radius
%Data fit for univariant +3 Ree Cations
%D values taken from Hauri et al, 1994. For long hour expirements 69 Hr @
%1430 C and 2.5 GPa 
%Ionic Radii Taken from Shannon and Prewitt


clear 


%------------------------------------------------------
%Data from Shannon and Prewitt and Hauri et al.
T=1430+273.15; %Temperature of Data from Hauri
Ree_order = cellstr(['La';'Ce';'Nd';'Sm';'Eu';'Dy';'Er';'Yb';'Lu']);%REE +3
ir = [1.16e-10 1.143e-10 1.109e-10 1.079e-10 1.066e-10 1.027e-10 1.004e-10 0.985e-10 0.977e-10]; %Ionic Radius of REE
D_Obs= [0.0515 0.108 0.277 0.462 0.458 0.711 0.66 0.633 0.623]; %D Values from Experimental Data

%------------------------------------------------------
%Strain Lattice Calculation
%Uses equation 2 from Blundy and Wood (1994), Nature

NA = 6.0221409e+23; %Avagadros Number (1/mol)
R = 8.3144598;      %Gas constant (J/(K*mol)), J=(kg*(m^2))/(s^2)

poss_Em=linspace(350e9,400e9,60); %Values of Youngs Modulus to be tested (Pa)
poss_Do=linspace(0.5,0.9,200); %Values of Do to be tested
poss_R=linspace(0.8e-10,1.3e-10,100); %Values of Ionic radius to be tested (m)

Perc_err(numel(poss_Em),numel(poss_Do),numel(poss_R))=0; %Set Matrix

%Number of Values considered
num_E = numel(poss_Em);
num_D = numel(poss_Do);
num_R = numel(poss_R);
num_ir = numel(ir);

%Preallocation
Di= zeros(1,num_ir);

%4 Nested for loop
for i= 1: num_E %Youngs Modulus 
    for j= 1: num_D %D Values
        for k= 1: num_R % Ionic radius to fit
            for e=1: num_ir %9 space vector to fit the REE

            %Equation from Blundy and Wood
            Di(e)= poss_Do(j).* exp((-4.*pi.*poss_Em(i).*NA./(R.*T)).* ...
                ((poss_R(k)./2).*(poss_R(k)-ir(e)).^2-(1/3.*(poss_R(k)-ir(e)).^3)));

          
            end
            %Calculates the Chi Squared from the data 
            %Fit outside the Ir loop to fit data
            Perc_err(i,j,k)=sum(((D_Obs-Di).^2)./D_Obs);
        end
    end
end

%Err is the error and locate holds the locations of the parameters
[err,locate] = min(Perc_err(:));
[x,y,z] = ind2sub(size(Perc_err),locate);

ax = gca; % current axes
ax.FontSize = 12;

%Plots the Chi Squared Values 
figure (1)
plot(log10(Perc_err(:)))
xlabel('Interations')
ylabel('\chi^{2}')
title('\chi^{2} values for Iterations')

%---------------------------------------
%Refit of data from grid search
ir_v=linspace(0.9e-10,1.2e-10,10000); %Ionic Radius vector over which to plot
%Plugs Solved values back into the equation
D_mod= poss_Do(y).* exp((-4.*pi.*poss_Em(x).*NA./(R.*T)).* ...    
    ((poss_R(z)./2).*(poss_R(z)-ir_v).^2-(1/3.*(poss_R(z)-ir_v).^3)));

%Plots the calculate D values and experimental values
figure(2)
xlabel('Ionic Radius (m)')
ylabel('ln (D)')
hold all; box on; grid on
semilogy(ir_v,log(D_mod),'-k',ir,log(D_Obs),'or','markerfacecolor','r')
text(ir,log(D_Obs),Ree_order,'VerticalAlignment','bottom','HorizontalAlignment','right')
title('3-D Lattice Strain Model Ree^{+3}') 
