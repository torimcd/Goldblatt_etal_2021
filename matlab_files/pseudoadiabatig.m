function dTdp = pseudoadiabatig(p,T,Rd,Md,cpd,condensablegas,condensedphase)
% Pseudoadiabat for an ideal gas
% Derived from eq 23 in sec9 of Iribarne & Godson, 1981 (p181),
% substituting eq 24-25, simplifiying denominator and setting as
% pseudoadiabatic following text after eq 25
% 
%global Liv    Llv    TLiv   TLlv   Tcpv   Tsati  Tsatl  cpv    psati  psatl
%global t1 t2 t3 t4
switch condensablegas
    case 'h2o'
        Mc = 18.003e-3;
        % get thermodynamic data
        %addpath([stddirpath('runaway'),'general']); %**Pre-run to save LOTS of time
        load satvp_h2o
        load L_h2o
        load cpv_h2o % note: this is on saturation vapour pressure curve
        % interpolate onto T used
        cpv1 = exp(interp1(Tcpv,log(cpv),T,'linear'));
        switch condensedphase
            case 'l'
                %tic
                es = exp(interp1(Tsatl,log(psatl),T,'linear')); % satruation vapour pressure
                %t1=t1+toc;
                L = interp1(TLlv,Llv,T,'linear','extrap'); %bugfix-hack: extrapolate for when using water at low temperature <200K
                %tic
                %t2=t2+toc;
            case 'i'
                es = exp(interp1(Tsati,log(psati),T,'linear')); % satruation vapour pressure
                L = interp1(TLiv,Liv,T,'linear');
            otherwise
                error('Condensed phase should be either l or i');
        end
    otherwise
        error('The condensable gas entered is not recognised');
end
    
epsilon = Mc/Md; % mass ratio. subscript c for condensable
if p==es
    error('(in pseudoadiabatig) Input p=%f and calculated es=%f are identical, which will give rw=Inf and dTdp=NaN',p,es);
end
rw = epsilon*es./(p-es);

dTdp = ((T.*Rd + L.*rw)./(p-es))./...
    (cpd + rw.*cpv1 + L.^2.*rw.*(epsilon+rw)./(Rd.*T.^2));
    