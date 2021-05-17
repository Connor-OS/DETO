vx = zeros(length(x),1); vy = vx;
freeze_step = 0;
tol_min = 5e-10;
tol = 2*tol_min;
nstep = 0;
dmax = 0.01;
%dmax = 0.2;
dt = 0.01;
%dt = 200.;
Eer = 0;
low_m = min(m);
%Fmax = 0;

if (dump1==true && loop==1) df1 = fopen('../dump/it1.dump','w'); end

% Euler integration
while tol > tol_min
    nstep = nstep + 1;
    Fx = Fxe;  Fy = Fye;
    Fmax = 0;
    % Compute force
    for i=1:length(x)
        for s=1:nn(i)
            j = N(i,s);
            dx = x(j)-x(i);     dy = y(j)-y(i);
            L(i,s) = sqrt( dx^2 + dy^2 );
            if j > i
                if (ptype == 1) F = k(i,s) * (L(i,s) - Li(i,s));
                elseif (ptype == 2) F = k(i,s)/apot * (-exp( -apot*(L(i,s) - Li(i,s)))+1 );
                elseif (ptype == 3) F = k(i,s)/apot * (exp( apot*(L(i,s) - Li(i,s)))-1 );
                elseif (ptype == 4) F = k(i,s)/apot * (tanh( apot*(L(i,s) - Li(i,s))) );
                elseif (ptype == 5) F = k(i,s)/apot * (sinh( apot*(L(i,s) - Li(i,s))) );
                end
                    
                %if (abs(F)>Fmax) Fmax=abs(F); end
                Fxx = F * dx/L(i,s);    Fyy = F * dy/L(i,s);
                Fx(i) = Fx(i) + Fxx;    Fx(j) = Fx(j) - Fxx;
                Fy(i) = Fy(i) + Fyy;    Fy(j) = Fy(j) - Fyy;
            end
        end
    end
    % CONSTRAINTS
    %Left support
    Fy(1) = 0;
    %Right support
    Fy(nelx) = 0;
    % Compute scale factor
    vdotf = dot(vx,Fx) + dot(vy,Fy);
    % Overshoot check
    if vdotf < 0
        freeze_step = nstep;
        vx = zeros(length(x),1); vy = vx;
    else
        % Scale velocities
        fdotf = dot(Fx,Fx) + dot(Fy,Fy);
        if fdotf == 0; scale = 0;
        else scale = vdotf/fdotf;
        end
        vx = scale .* Fx;
        vy = scale .* Fy;
    end
    % Limit dmax
    vmax = max(max(abs(vx)),max(abs(vy)));
    if dt*vmax > dmax; dt = dmax/vmax;
    end
    % Euler integration step
    x = x + dt .* vx;
    y = y + dt .* vy;
    for i = 1:length(x)
        if m(i) > 1e-10
            vx(i)= vx(i)+dt/1*Fx(i);
            vy(i)= vy(i)+dt/1*Fy(i);
        else
            vx(i) = 0;
            vy(i) = 0;
        end
    end
    % Tolerance
    Eold = Eer;
    Eer = 0;
    for i = 1:length(xi)
        for s = 1:nn(i)
            j = N(i,s);
            % all potential factors divided by two due to neighbour double counting
            rij = L(i,s)-Li(i,s);
            if (ptype == 1) potfac = 1./4. * kspr*(rij^2);
                elseif (ptype == 2) potfac = 1./2. * ( kspr/apot * (exp(-apot*rij)/apot+rij ) - kspr/apot/apot );
                elseif (ptype == 3) potfac = 1./2. * ( kspr/apot * (exp(apot*rij)/apot-rij ) - kspr/apot/apot );
                elseif (ptype == 4) potfac = 1./2. * kspr/apot/apot * log(cosh(apot*rij));
                elseif (ptype == 5) potfac = 1./2. * ( kspr/apot/apot * cosh(apot*rij) - kspr/apot/apot );
            end
            Eer = Eer + potfac * m(i)^2 *m(j)^2;
            %Eer = Eer+1/4*m(i)^2*m(j)^2*kspr*(L(i,s)-Li(i,s))^2;
        end
    end
    if nstep - freeze_step < 5
        tol = 2*tol_min;
    else
        tol = abs(Eer-Eold)/abs(Eold);
    end
    
    
    % dumping xyz evolution during first iteration, for debugging 
    if (dump1==true && loop==1 && mod(nstep,1000)==1)
        fprintf(df1,[sprintf('%4i\n',length(m)),'Frame.: ',sprintf('%4i\n',nstep)]);
        for i=1:length(x)
        fprintf(df1,[sprintf('%4i ',i) sprintf('%6.3f ',x(i)) sprintf('%6.3f ',y(i)) ...
            sprintf('%6.3f \n',m(i)) ]);
         end
    end
    
end
if (dump1==true && loop==1) fclose(df1); end
