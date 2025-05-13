function [sys,x0,str,ts,simStateCompliance] = distill(t,y,u,flag)

switch flag,
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes();
        
        %%%%%%%%%%%%%%%
        % Derivatives %
        %%%%%%%%%%%%%%%
    case 1,
        sys=mdlDerivatives(t,y,u);
        
        %%%%%%%%%%
        % Update %
        %%%%%%%%%%
    case 2,
        sys=mdlUpdate(t,y,u);
        
        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3,
        sys=mdlOutputs(t,y,u);
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % GetTimeOfNextVarHit %
        %%%%%%%%%%%%%%%%%%%%%%%
    case 4,
        sys=mdlGetTimeOfNextVarHit(t,y,u);
        
        %%%%%%%%%%%%%
        % Terminate %
        %%%%%%%%%%%%%
    case 9,
        sys=mdlTerminate(t,y,u);
        
        %%%%%%%%%%%%%%%%%%%%
        % Unexpected flags %
        %%%%%%%%%%%%%%%%%%%%
    otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
        
end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes()

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 32;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
%
% initialize the initial conditions
x_ss(1:32)=[0.93541941614016
   0.90052553715795
   0.86229645132283
   0.82169940277993
   0.77999079584355
   0.73857168629759
   0.69880490932694
   0.66184253445732
   0.62850777645505
   0.59925269993058
   0.57418567956453
   0.55314422743545
   0.53578454439850
   0.52166550959767
   0.51031495114413
   0.50127509227528
   0.49412891686784
   0.48544992019184
   0.47420248108803
   0.45980349896163
   0.44164297270225
   0.41919109776836
   0.39205549194059
   0.36024592617390
   0.32407993023343
   0.28467681591738
   0.24320921343484
   0.20181568276528
   0.16177269003094
   0.12514970961746
   0.09245832612765
   0.06458317697321];

x0 = x_ss';

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%==========================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
%
% Model source:
% R. Palma and T.F. Edgar, Toward Patient Specific Insulin Therapy: A Novel
%    Insulin Bolus Calculator.  In Proceedings Texas Wisconsin California Control
%    Consortium, Austin, TX, Feb. 7-8, 2011.
%
% Expanded Bergman Minimal model to include meals and insulin
% Parameters for an insulin dependent type-I diabetic

% Inputs (2):
% Reflux ratio is the Manipulated variable
% Reflux Ratio (L/D)
rr = u(1);

% Disturbance variables (DV)
% Feed Flowrate (mol/min)
Feed = u(2);

% Mole Fraction of Feed
x_Feed = u(3);

% States (32):
% x(1) - Reflux Drum Liquid Mole Fraction of Component A
% x(2) - Tray 1 - Liquid Mole Fraction of Component A
% .
% .
% .
% x(17) - Tray 16 - Liquid Mole Fraction of Component A (Feed Location)
% .
% .
% .
% x(31) - Tray 30 - Liquid Mole Fraction of Component A
% x(32) - Reboiler Liquid Mole Fraction of Component A

% Parameters
% Feed Flowrate (mol/min)
%Feed =  24.0/60.0;
% Mole Fraction of Feed
%x_Feed = 0.5;
% Distillate Flowrate (mol/min)
D=0.5*Feed;
% Flowrate of the Liquid in the Rectification Section (mol/min)
L=rr*D;
% Vapor Flowrate in the Column (mol/min)
V=L+D;
% Flowrate of the Liquid in the Stripping Section (mol/min)
FL=Feed+L;
% Relative Volatility = (yA/xA)/(yB/xB) = KA/KB = alpha(A,B)
vol=1.6;
% Total Molar Holdup in the Condenser
atray=0.25;
% Total Molar Holdup on each Tray
acond=0.5;
% Total Molar Holdup in the Reboiler
areb=1.0;
% Vapor Mole Fractions of Component A
% From the equilibrium assumption and mole balances
% 1) vol = (yA/xA) / (yB/xB)
% 2) xA + xB = 1
% 3) yA + yB = 1
y(1)=x(1)*vol/(1+(vol-1)*x(1));
y(2)=x(2)*vol/(1+(vol-1)*x(2));
y(3)=x(3)*vol/(1+(vol-1)*x(3));
y(4)=x(4)*vol/(1+(vol-1)*x(4));
y(5)=x(5)*vol/(1+(vol-1)*x(5));
y(6)=x(6)*vol/(1+(vol-1)*x(6));
y(7)=x(7)*vol/(1+(vol-1)*x(7));
y(8)=x(8)*vol/(1+(vol-1)*x(8));
y(9)=x(9)*vol/(1+(vol-1)*x(9));
y(10)=x(10)*vol/(1+(vol-1)*x(10));
y(11)=x(11)*vol/(1+(vol-1)*x(11));
y(12)=x(12)*vol/(1+(vol-1)*x(12));
y(13)=x(13)*vol/(1+(vol-1)*x(13));
y(14)=x(14)*vol/(1+(vol-1)*x(14));
y(15)=x(15)*vol/(1+(vol-1)*x(15));
y(16)=x(16)*vol/(1+(vol-1)*x(16));
y(17)=x(17)*vol/(1+(vol-1)*x(17));
y(18)=x(18)*vol/(1+(vol-1)*x(18));
y(19)=x(19)*vol/(1+(vol-1)*x(19));
y(20)=x(20)*vol/(1+(vol-1)*x(20));
y(21)=x(21)*vol/(1+(vol-1)*x(21));
y(22)=x(22)*vol/(1+(vol-1)*x(22));
y(23)=x(23)*vol/(1+(vol-1)*x(23));
y(24)=x(24)*vol/(1+(vol-1)*x(24));
y(25)=x(25)*vol/(1+(vol-1)*x(25));
y(26)=x(26)*vol/(1+(vol-1)*x(26));
y(27)=x(27)*vol/(1+(vol-1)*x(27));
y(28)=x(28)*vol/(1+(vol-1)*x(28));
y(29)=x(29)*vol/(1+(vol-1)*x(29));
y(30)=x(30)*vol/(1+(vol-1)*x(30));
y(31)=x(31)*vol/(1+(vol-1)*x(31));
y(32)=x(32)*vol/(1+(vol-1)*x(32));


% Compute xdot
xdot(1) = 1/acond*V*(y(2)-x(1));
xdot(2) = 1/atray*(L*(x(1)-x(2))-V*(y(2)-y(3)));
xdot(3) = 1/atray*(L*(x(2)-x(3))-V*(y(3)-y(4)));
xdot(4) = 1/atray*(L*(x(3)-x(4))-V*(y(4)-y(5)));
xdot(5) = 1/atray*(L*(x(4)-x(5))-V*(y(5)-y(6)));
xdot(6) = 1/atray*(L*(x(5)-x(6))-V*(y(6)-y(7)));
xdot(7) = 1/atray*(L*(x(6)-x(7))-V*(y(7)-y(8)));
xdot(8) = 1/atray*(L*(x(7)-x(8))-V*(y(8)-y(9)));
xdot(9) = 1/atray*(L*(x(8)-x(9))-V*(y(9)-y(10)));
xdot(10) = 1/atray*(L*(x(9)-x(10))-V*(y(10)-y(11)));
xdot(11) = 1/atray*(L*(x(10)-x(11))-V*(y(11)-y(12)));
xdot(12) = 1/atray*(L*(x(11)-x(12))-V*(y(12)-y(13)));
xdot(13) = 1/atray*(L*(x(12)-x(13))-V*(y(13)-y(14)));
xdot(14) = 1/atray*(L*(x(13)-x(14))-V*(y(14)-y(15)));
xdot(15) = 1/atray*(L*(x(14)-x(15))-V*(y(15)-y(16)));
xdot(16) = 1/atray*(L*(x(15)-x(16))-V*(y(16)-y(17)));
xdot(17) = 1/atray*(Feed*x_Feed+L*x(16)-FL*x(17)-V*(y(17)-y(18)));
xdot(18) = 1/atray*(FL*(x(17)-x(18))-V*(y(18)-y(19)));
xdot(19) = 1/atray*(FL*(x(18)-x(19))-V*(y(19)-y(20)));
xdot(20) = 1/atray*(FL*(x(19)-x(20))-V*(y(20)-y(21)));
xdot(21) = 1/atray*(FL*(x(20)-x(21))-V*(y(21)-y(22)));
xdot(22) = 1/atray*(FL*(x(21)-x(22))-V*(y(22)-y(23)));
xdot(23) = 1/atray*(FL*(x(22)-x(23))-V*(y(23)-y(24)));
xdot(24) = 1/atray*(FL*(x(23)-x(24))-V*(y(24)-y(25)));
xdot(25) = 1/atray*(FL*(x(24)-x(25))-V*(y(25)-y(26)));
xdot(26) = 1/atray*(FL*(x(25)-x(26))-V*(y(26)-y(27)));
xdot(27) = 1/atray*(FL*(x(26)-x(27))-V*(y(27)-y(28)));
xdot(28) = 1/atray*(FL*(x(27)-x(28))-V*(y(28)-y(29)));
xdot(29) = 1/atray*(FL*(x(28)-x(29))-V*(y(29)-y(30)));
xdot(30) = 1/atray*(FL*(x(29)-x(30))-V*(y(30)-y(31)));
xdot(31) = 1/atray*(FL*(x(30)-x(31))-V*(y(31)-y(32)));
xdot(32) = 1/areb*(FL*x(31)-(Feed-D)*x(32)-V*y(32));

% derivatives
sys = xdot';
% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,y,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,y,u)

y1 = y(1);

sys = [y1];

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,y,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,y,u)

sys = [];

% end mdlTerminate
