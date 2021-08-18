% A DISCRETE ELEMENT TOPOLOGY OPTIMIZATION CODE BY CONNOR O'SHAUGHNESSY 2020%
function top(nelx,nely,mass,Diam,rmin,pot,kspr,apot,Mdam)
% nelx, nely = number of particles in x and y directions
% mass = total target mass of the system (optimisation constraint)
% Diam = particle diameter
% rmin = filtering radius
% pot = potential type, can be 'lin_spr', 'e-x', $e+x', tanh(x), sinh(x)
% kspr, apot = potential parameters (for linear type, only kspr is used)
% Mdam = damage matrix with damage x,y positions, radius, and weight factor

disp('DETO STARTED');
tic

dump1 = true;   % if dump=true, the deformation history during iteration 1 is dumped

ptype = 0;  % potential type, lin_spr=1, e-x=2, e+x=3, tanh(x)=4, sinh(x)=5
if (strcmp(pot,'lin_spr')==1)   ptype = 1; end
if (strcmp(pot,'e-x')==1)   ptype = 2; end
if (strcmp(pot,'e+x')==1)   ptype = 3; end
if (strcmp(pot,'tanh(x)')==1)   ptype = 4; end
if (strcmp(pot,'sinh(x)')==1)   ptype = 5; end
if (ptype == 0)   
    msg = sprintf('ERROR: unknown type of interaction: %s',pot);
    disp(msg);
    return
end

% Printing input damage matrix, just to check it is alright
if size(Mdam,1)==0
    disp('Empty damage matrix provided. Optimising structure without damage.')
else
    msg = sprintf('\nNumber of damage scenarios: %d',size(Mdam,1));
    disp(msg);
    disp('X-pos, Y-pos, radius and weight of damage:');
    disp('Xdam  Ydam  Rdam  Wdam');
    for i=1:size(Mdam,1)
        msg = sprintf('%f %f %f %f',Mdam(i,1),Mdam(i,2),Mdam(i,3),Mdam(i,4));
        disp(msg);
    end
end


% INITIALIZATION
m = mass*ones((nely*nelx)-floor(nely/2),1);
cut = 1.01*Diam;
xi = zeros(length(m),1);  yi = xi;
Fxe = zeros(length(m),1);   Fye = Fxe;
Fye(length(m)-floor(nelx/2)) = -0.5;
Fye(length(m)-floor(nelx/2)-1) = -0.5;
Fye(length(m)-floor(nelx/2)+1) = -0.5;
% Initial geometry
for j = 1:nely
    for i = 1:nelx-mod(j+1,2)
        xi((j-1)*nelx+i-floor((j-1)/2)) = Diam*(i-0.5)+(1-mod(j,2))*Diam/2-(nelx*Diam)/2;
        yi((j-1)*nelx+i-floor((j-1)/2)) = Diam*((j-1)*cos(pi/6)+0.0);
    end
end
% Neighbour list and equilibrium distances
neighbour_list
% Variables and plotting
loop = 0;
change = 1.;
x = xi; y = yi;
xy = sprintf('../dump/xy_%i_%i_%1.2f_%1.2f_%1.2f_%i.txt',nelx,nely,mass,Diam,rmin,kspr);
it = sprintf('../dump/it_%i_%i_%1.2f_%1.2f_%1.2f_%i.txt',nelx,nely,mass,Diam,rmin,kspr);
XY = fopen(xy,'w');
Itt = fopen(it,'w');
m_old = rand(length(m),10);
% START ITERATION
while change > 0.004
%     x = xi; y = yi;
    loop = loop+1;
    msg = sprintf('\nITERATION %d',loop);
    disp(msg);
    for i = 1:9
        m_old(:,i) = m_old(:,i+1); % Records prev 10 itteration history
    end
    m_old(:,10) = m;
    % updtaing interaction stiffnesses based on latest m-variables
    for i = 1:length(m)
        for s = 1:nn(i)
            j = N(i,s);
            k(i,s) = m(i)^2 *m(j)^2  * kspr;
        end
    end
    
    % ------------------   damage loop starts here here 
    % before entering loop must define a combined sensitivity vector to sum
    % individual damage weighted sensitivities. 
    
    % Each time a new damage is assessed we must resume the initial confir   x = xi; y = yi;
   
    % For each damage, make list of removed elements ids, record their masses, then set their m =
    % 0, run Qmin, in filtering hard kill by forcing dc/dxi = 0 for removed masses
    
    % Before moving to next damage (or out into optimisation step), restore
    % original masses and add weighted contribution to global sensitivity
    % vector
    
    % After damage loop, record the combined sensitivity vector into dc:
    % this is the dc to which the volume constraint will be applied in
    % OC
    
    % DEM-ANALYSIS assuming harmonic potential
    disp('RUNNING DEM MINIMISATION...');
    Qmin
    %msg = sprintf('\nMax Force on spring is %e',Fmax);
    %disp(msg);
    disp('... DONE');
    
    % Cost function and sensitivities
    %Eer = 0;     % Eer already comes from Qmin - no need to recompute it
    dc = zeros(length(xi),1);
    for i = 1:length(xi)
        for s = 1:nn(i)
            j = N(i,s);
            %Eer = Eer+1/4*m(i)^2*m(j)^2*kspr*(L(i,s)-Li(i,s))^2; % 1/4 cause double counting each interaction
            rij = L(i,s)-Li(i,s);
            if (ptype == 1) potfac = 1./4. * kspr*(rij^2);
                elseif (ptype == 2) potfac = 1./2. * ( kspr/apot * (exp(-apot*rij)/apot+rij ) - kspr/apot/apot );
                elseif (ptype == 3) potfac = 1./2. * ( kspr/apot * (exp(apot*rij)/apot-rij ) - kspr/apot/apot );
                elseif (ptype == 4) potfac = 1./2. * kspr/apot/apot * log(cosh(apot*rij));
                elseif (ptype == 5) potfac = 1./2. * ( kspr/apot/apot * cosh(apot*rij) - kspr/apot/apot );
            end
            %dc(i) = dc(i)-0.5*m(i)*m(j)^2*kspr*(L(i,s)-Li(i,s))^2;
            dc(i) = dc(i) - 2.*m(i)*m(j)^2 * potfac;
        end
       % msg = sprintf('dc(%d) = %e',i,dc(i));
        %disp(msg);
    end

    % Filtering of sensitivities (aka coarse graining)
    if rmin > Diam
        Mesh_filter
    end
    
    % ------------------   damage loop finishes here
    
    
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    OC
    % PER PARTICLE STRESS
    stress_perpart
    change = 0.2;
    for i = 1:10
        if change > max(abs(m-m_old(:,i)))
            change = max(abs(m-m_old(:,i))); % Computes lowest change - Past 10 turns
        end
    end
    % PRINT RESULTS
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%6.4f',Eer) ...
       ' Mass.: ' sprintf('%6.3f ',sum(m)/(length(xi))) ...
       ' ch.:' sprintf('%6.3f',change) ' Time.: ' sprintf('%4.2f',toc)])
    % PLOT DENSITIES
    cc=[1-m,1-m,1-m];
    linkdata on
    scatter(x,y,Diam*50,cc,'filled')
    % WRITE DUMP FILES
    fprintf(XY,[sprintf('%4i\n',length(m)),'Frame.: ',sprintf('%4i\n',loop)]);
    for i=1:length(x)
        fprintf(XY,[sprintf('%4i ',i) sprintf('%6.3f ',x(i)) sprintf('%6.3f ',y(i)) ...
            sprintf('%6.3f ',m(i)) sprintf('%6.3f ',W(i,1)) sprintf('%6.3f ',W(i,2)) ...
            sprintf('%6.3f ',W(i,3)) sprintf('%6.3f ',W(i,4)) sprintf('%6.3f ',W(i,5))...
            sprintf('%6.3f\n',W(i,6))]);
    end
    fprintf(Itt,[' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',Eer) ...
        ' Mass.: ' sprintf('%6.3f',sum(m)/(length(xi))) ...
        ' ch.: ' sprintf('%6.3f',change ) ' Time.: ' sprintf('%6.3f\n',toc)]);
 end
fprintf(Itt,[' Sim time.: ' sprintf( '%6.3f\n', toc)]);
%%%%% POST PROCESSING %%%%%%%
post_processing
disp('\nDETO DONE: SUCCESS');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Enrico Masoero and Connor O'Shaughnessy  %
% School of Engineering Newcastle University                               %
% Please sent your comments to the author: c.o'shaughnessy1@ncl.ac.uk      %
%                                                                          %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but does not guaranty that the code is   %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
